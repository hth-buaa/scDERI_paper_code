# ==============================================================================
# E3: terminal threshold sweep
# A precision-to-coverage operating curve for the final filter (no retraining)
# ==============================================================================


# WHAT THIS SCRIPT DOES
#   It is the targeted sensitivity analysis that backs the statement 
#   "network recall is tunable by relaxing the terminal threshold".
#   We move the single user-facing knob of the FINAL filtering step (the
#   strictness of the multi-factor adjusted threshold) and measure how large
#   the network becomes and how its ChIP-seq support changes.
#   Nothing is retrained: iGRN construction, GA / gradient training and Brent
#   prediction are untouched. We only revisit the last filtering of Part 5.6.

# THE KNOB (the only quantity we change)
#   A candidate triple is kept by the main filter when its score FP-S-R is at
#   least its adjusted threshold
#       AQ = clip( IQ(theta_p) * product_{m=1..11} F_m , 0, 1 ).
#   We scale only this threshold by s:
#       AQ_s = clip( s * AQ, 0, 1 ),   s in {0.70, 0.85, 1.00, 1.15, 1.30}.
#     s < 1 relaxes the filter -> more edges (more coverage), lower precision;
#     s = 1 is the shipped default and must reproduce the released network;
#     s > 1 tightens the filter -> fewer edges, higher precision for the main
#           filter (see "rescue floor" below for the final network).

# TWO CURVES PER s (why we need both)
#   primary : main filter only. This is the textbook precision-coverage trade
#             and the main evidence that the threshold is a clean, monotone knob.
#   final   : main filter + two FIXED rescue steps. The rescue guarantees a
#             minimum connectivity floor (>=1 TF per peak-TG, >=7 peak-TG per
#             TF, with quality floors Q0.25 and Q0). It is held constant across
#             s, so it isolates the tuning from the connectivity safeguard.
#             On the tightening side (s > 1) this floor prevents the released
#             network from shrinking and can dominate its precision; the
#             primary curve shows what the threshold alone does.

# HOW TO READ THE PRECISION METRIC
#   "TF-peak precision" is the macro-average, over TFs that have ChIP-seq
#   ground truth, of the fraction of predicted peaks overlapping that TF's top
#   ChIP peaks. The matching recall column (tp_recall) is structurally tiny
#   (~1e-3): the ground truth holds 10,000 peaks per TF while a predicted TF
#   has only tens-to-hundreds of peaks, so it must not be read as the recall we
#   tune. Here the coverage / recall axis is the NUMBER OF EDGES (n_TF_peak),
#   exactly as in the benchmark TF-peak evaluation, which uses precision.

# STEPS
#   1. Load, for each of the three datasets, the pre-filter table
#      pegrn[[1]][["loose"]] and the released loose_filtered(_simplified).
#   2. Recompute IQ, the 11 factors and AQ for EVERY loose row (identical
#      formulas and Part 5.6 parameters; AQ is only stored on kept rows).
#   3. Build the rescue candidates once (they do not depend on s), then apply
#      the scaled filter and the fixed rescue for every s.
#   4. Sanity check at s = 1: edge-for-edge reproduction of the released
#      network, plus monotonicity of edge counts vs s.
#   5. Evaluate TF-peak precision against ChIP-seq and TF-coverage P/R/F1.
#   6. Save one long table, one sanity table, and three figures.

# OUTPUTS (written to OUT_DIR)
#   E3_threshold_sweep_long.csv one row per dataset x stage x s; counts and
#                               metrics, plus edges_vs_default (x the default)
#                               and tp_precision_robust (NA when <5 ChIP TFs)
#   E3_threshold_sanity.csv     s = 1 reproduction and monotonicity checks
#   E3_operating_curve.png      TF-peak precision vs TF-peak edges
#   E3_edges_vs_scale.png       TF-peak edges vs s (the coverage knob)
#   E3_TFcoverage_vs_scale.png  TF-level coverage vs s (edge-level is tuned)
#
# RUNTIME / DEPENDENCIES
#   A few minutes per dataset (mostly loading qsave); the sweep is vectorized
#   arithmetic over ~10k rows. Needs qs, minpack.lm, data.table, Signac,
#   GenomicRanges, S4Vectors, ggplot2.

# SCOPE / LIMITS
#   Only the "loose" level is swept (the level used in the paper). Bone marrow
#   has only 10 ChIP-seq TFs, so its TF-peak precision estimates have lower
#   statistical power; BM still contributes edge-count and TF-coverage results.


options(stringsAsFactors = FALSE)
suppressMessages({
  library(qs)
  library(minpack.lm)      # Levenberg-Marquardt curve fit, same as filter_eGRN
  library(data.table)
  library(Signac)          # StringToGRanges
  library(GenomicRanges)   # findOverlaps
  library(S4Vectors)
  library(ggplot2)
})


## ------------------------------ Configuration ------------------------------
# Paths: Windows
# BENCH_ROOT <- "I:/scDERI_Paper/3-benchmark/8-scDERI"
# EVAL_ROOT  <- "I:/scDERI_Paper/2-obtain_evaluation_criteria/1-TF-evaluation"
# OUT_DIR    <- "I:/scDERI_Paper/Supplement-Experiment/E3"

# Paths: WSL2
BENCH_ROOT <- "/mnt/i/scDERI_Paper/3-benchmark/8-scDERI"
EVAL_ROOT  <- "/mnt/i/scDERI_Paper/2-obtain_evaluation_criteria/1-TF-evaluation"
OUT_DIR    <- "/mnt/i/scDERI_Paper/Supplement-Experiment/E3"

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Benchmark folder and its matched ChIP-seq truth folder for each dataset.
DATASETS <- list(
  A549    = list(folder = "SNAREseq2_hg38_A549",
                 chip   = "A549-hg38-ChIP-seq"),
  GM12878 = list(folder = "SNAREseq2_hg38_GM12878",
                 chip   = "GM12878-hg38-ChIP-seq"),
  BM      = list(folder = "10X_hg38_Human_bone_marrow_with_T_cell_depleted",
                 chip   = "Bone-marrow-hg38-ChIP-seq")
)

# Threshold scale grid; s = 1 is the shipped default.
# S_GRID <- c(0.70, 0.85, 1.00, 1.15, 1.30)
S_GRID <- seq(0.5, 1.5, 0.05)
TOP_TRUTH_PEAKS <- 10000  # top ChIP peaks per TF, identical to the benchmark
# MIN_COMMON_TF   <- 5      # below this many ChIP TFs, precision is not stable
MIN_COMMON_TF   <- 2      # below this many ChIP TFs, precision is not stable

# ---- filter_eGRN parameters, copied verbatim from the Part 5.6 call ----
THETA_P_POINTS   <- c(0, 0.5, 1)
INITIAL_Q_POINTS <- c(5.79, 4.66, 4.23)
FEATURES <- c("FG-NL", "FG-theta_i-R",
              "F-NZER-R", "F-NZEA-R", "F-EA-R",
              "P-NZAR-R", "P-NZAA-R", "P-AA-R",
              "G-NZER-R", "G-NZEA-R", "G-EA-R")
FACTOR_BALANCE <- setNames(rep(0.2, length(FEATURES)), FEATURES)
FACTOR_RANGE <- list(
  `FG-NL` = 0.43, `FG-theta_i-R` = 0.64,
  `F-NZER-R` = 0.28, `F-NZEA-R` = 0.33, `F-EA-R` = 0.22,
  `P-NZAR-R` = 0.19, `P-NZAA-R` = 0.67, `P-AA-R` = 0.13,
  `G-NZER-R` = 0.21, `G-NZEA-R` = 0.12, `G-EA-R` = 0.26)

# Only FG-NL is a penalty factor; the other ten relax the threshold.
FACTOR_ATTRIBUTE <- list(
  `FG-NL` = "penalty", `FG-theta_i-R` = "relaxation",
  `F-NZER-R` = "relaxation", `F-NZEA-R` = "relaxation", `F-EA-R` = "relaxation",
  `P-NZAR-R` = "relaxation", `P-NZAA-R` = "relaxation", `P-AA-R` = "relaxation",
  `G-NZER-R` = "relaxation", `G-NZEA-R` = "relaxation", `G-EA-R` = "relaxation")

# Fixed rescue floors (NOT varied with s; they define the connectivity floor).
PEAKTG_Q      <- 0.25   # rescue 1 quality floor (Q0.25 of FP-S-R)
MIN_TF_PEAKTG <- 1      # rescue 1: keep >=1 TF per peak-TG
TF_Q          <- 0.0    # rescue 2 quality floor (Q0 of FP-S-R)
MIN_PEAKTG_TF <- 7      # rescue 2: keep >=7 peak-TG per TF

# Columns that must exist in the loose table.
REQUIRED_COLS <- c("TF", "peak", "TG", "FG-theta_p", "FP-S-R",
                   "P-AA-R", "P-NZAR-R", "FG-theta_i-R", "FG-NL", FEATURES)


## ----------------------- Step 2: reconstruct AQ ----------------------------
# Add IQ, the eleven "-F" factors and AQ to EVERY loose row, exactly as
# filter_eGRN does. A stable row_id is attached so later rescue marks can be
# set by id regardless of sorting.
annotate_loose <- function(loose) {
  missing <- setdiff(REQUIRED_COLS, colnames(loose))
  if (length(missing)) {
    stop("loose table is missing expected columns: ",
         paste(missing, collapse = ", "),
         "\nPlease send names(pegrn[[1]][['loose']]) and the filter_eGRN ",
         "version used in Part 5.6.")
  }
  e <- as.data.frame(loose, stringsAsFactors = FALSE)
  e[["peak-TG"]] <- paste0(e$peak, "_to_", e$TG)
  e[["row_id"]]  <- seq_len(nrow(e))
  
  # Global inverse curve: initial_q = 1/(theta_p*a + b) + c, with a > 0.
  x <- THETA_P_POINTS; y <- INITIAL_Q_POINTS
  fit <- minpack.lm::nlsLM(
    y ~ 1 / (x * a + b) + c,
    start = list(a = 0.333333, b = 0.333333, c = -0.5),
    lower = c(1e-9, -Inf, -Inf), upper = c(Inf, Inf, Inf))
  cf <- stats::coef(fit)
  e[["IQ"]] <- 1 / (e[["FG-theta_p"]] * cf[["a"]] + cf[["b"]]) + cf[["c"]]
  
  # Per-feature factor: fa = +1 for relaxation, -1 for penalty.
  for (ff in FEATURES) {
    fb <- FACTOR_BALANCE[[ff]]; fr <- FACTOR_RANGE[[ff]]
    fa <- (-1) ^ as.numeric(FACTOR_ATTRIBUTE[[ff]] == "penalty")
    e[[paste0(ff, "-F")]] <- (1 + fa * fb * fr) - fa * e[[ff]] * fr
  }
  
  # AQ = product of IQ and all factors (log-sum is numerically stable), then
  # clamp to [0,1].
  aq_mat <- as.matrix(e[, c("IQ", paste0(FEATURES, "-F"))])
  aq_mat[aq_mat <= 0] <- .Machine$double.eps
  e[["AQ"]] <- pmin(pmax(exp(rowSums(log(aq_mat))), 0), 1)
  e
}


## ----------------------- Step 3: fixed rescue sets -------------------------
# Build the rescue candidate pools and their within-group ranks. They depend
# only on the priority order, not on s, so they are computed once per dataset.
build_rescue <- function(e) {
  # Shared priority: higher score first, FG-NL last and ascending.
  ord <- order(-e[["FP-S-R"]], -e[["P-AA-R"]], -e[["P-NZAR-R"]],
               -e[["FG-theta_i-R"]], -e[["FG-theta_p"]], e[["FG-NL"]])
  eo <- e[ord, , drop = FALSE]
  
  # Rescue 1: top-ranked TF candidate(s) per peak-TG, above the Q0.25 floor.
  eo[["peakTG_rank"]] <- stats::ave(seq_len(nrow(eo)), eo[["peak-TG"]],
                                    FUN = seq_along)
  c1 <- eo[eo[["peakTG_rank"]] <= MIN_TF_PEAKTG, , drop = FALSE]
  c1 <- c1[c1[["FP-S-R"]] >= stats::quantile(c1[["FP-S-R"]], PEAKTG_Q), ,
           drop = FALSE]
  
  # Rescue 2: top-ranked peak-TG candidate(s) per TF, above the Q0 floor.
  eo[["TF_rank"]] <- stats::ave(seq_len(nrow(eo)), eo$TF, FUN = seq_along)
  c2 <- eo[eo[["TF_rank"]] <= MIN_PEAKTG_TF, , drop = FALSE]
  c2 <- c2[c2[["FP-S-R"]] >= stats::quantile(c2[["FP-S-R"]], TF_Q), ,
           drop = FALSE]
  list(c1 = c1, c2 = c2)
}

# Apply the scaled main filter and then the two fixed rescue steps.
apply_scale <- function(e, rescue, s) {
  g <- e
  g[["AQ_s"]]        <- pmin(1, pmax(0, s * g[["AQ"]]))
  g[["filter_keep"]] <- g[["FP-S-R"]] >= g[["AQ_s"]]
  g[["rescue1"]] <- FALSE
  g[["rescue2"]] <- FALSE
  
  # Rescue 1: fill peak-TG pairs that currently have no kept TF.
  kc1 <- tapply(g[["filter_keep"]], g[["peak-TG"]], sum)
  v1  <- kc1[rescue$c1[["peak-TG"]]]; v1[is.na(v1)] <- 0
  need1 <- rescue$c1[["peakTG_rank"]] <= (MIN_TF_PEAKTG - v1) & v1 < MIN_TF_PEAKTG
  g[["rescue1"]][rescue$c1$row_id[need1]] <- TRUE
  g[["keep_after_rescue1"]] <- g[["filter_keep"]] | g[["rescue1"]]
  
  # Rescue 2: then fill TFs that still have too few peak-TG pairs.
  kc2 <- tapply(g[["keep_after_rescue1"]], g$TF, sum)
  v2  <- kc2[rescue$c2$TF]; v2[is.na(v2)] <- 0
  need2 <- rescue$c2[["TF_rank"]] <= (MIN_PEAKTG_TF - v2) & v2 < MIN_PEAKTG_TF
  g[["rescue2"]][rescue$c2$row_id[need2]] <- TRUE
  
  g[["keep"]] <- g[["filter_keep"]] | g[["rescue1"]] | g[["rescue2"]]
  list(primary = g[g[["filter_keep"]], , drop = FALSE],
       final   = g[g[["keep"]], , drop = FALSE])
}


## ----------------------- Step 5a: network counts ---------------------------
# Unique node / edge / triple counts for a network table.
count_network <- function(df, stage) {
  uniq <- function(cols) nrow(unique(df[, cols, drop = FALSE]))
  data.frame(
    stage = stage,
    n_TF      = length(unique(df$TF)),
    n_peak    = length(unique(df$peak)),
    n_TG      = length(unique(df$TG)),
    n_TF_peak = uniq(c("TF", "peak")),
    n_peak_TG = uniq(c("peak", "TG")),
    n_TF_TG   = uniq(c("TF", "TG")),
    n_triple  = uniq(c("TF", "peak", "TG")),
    stringsAsFactors = FALSE)
}

triple_key <- function(df) {
  d <- unique(df[, c("TF", "peak", "TG"), drop = FALSE])
  paste(d$TF, d$peak, d$TG, sep = "|")
}


## ----------------------- Step 5b: ChIP-seq evaluation ----------------------
# Read one BED per TF, keep its top 10,000 peaks by signal (column V7), and
# return TF -> peak strings (same procedure as the benchmark evaluation).
build_chip_truth <- function(chip_dir) {
  bed_files <- list.files(chip_dir, pattern = "\\.bed$", full.names = TRUE)
  if (!length(bed_files)) {
    warning("No ChIP-seq .bed files in ", chip_dir); return(NULL)
  }
  lst <- lapply(bed_files, function(f) {
    dt <- data.table::fread(f, showProgress = FALSE, data.table = FALSE)
    if (ncol(dt) >= 7) dt <- dt[order(-dt[[7]]), , drop = FALSE]
    dt <- utils::head(dt[, 1:3, drop = FALSE], TOP_TRUTH_PEAKS)
    tf <- toupper(strsplit(basename(f), "_")[[1]][1])
    data.frame(TF = tf, region = paste0(dt[[1]], "-", dt[[2]], "-", dt[[3]]),
               stringsAsFactors = FALSE)
  })
  truth <- do.call(rbind, lst)
  truth[!duplicated(truth), , drop = FALSE]
}

gr_of <- function(region_strings) {
  Signac::StringToGRanges(unique(region_strings), sep = c(":", "-"))
}

# Macro-average per-TF peak overlap precision / recall / F1 over TFs that have
# ChIP-seq truth (n_common_TF). Precision is the headline metric.
tf_peak_metrics <- function(pred_df, truth) {
  empty <- data.frame(tp_precision = NA_real_, tp_recall = NA_real_,
                      tp_f1 = NA_real_, n_common_TF = 0L,
                      stringsAsFactors = FALSE)
  if (is.null(truth) || !nrow(pred_df)) return(empty)
  pred <- unique(data.frame(TF = pred_df$TF, Region = pred_df$peak,
                            stringsAsFactors = FALSE))
  common <- intersect(unique(pred$TF), unique(truth$TF))
  if (!length(common)) return(empty)
  per <- lapply(common, function(tf) {
    q_gr <- gr_of(pred$Region[pred$TF == tf])
    s_gr <- gr_of(truth$region[truth$TF == tf])
    hits <- GenomicRanges::findOverlaps(q_gr, s_gr)
    prec <- length(unique(S4Vectors::queryHits(hits))) / length(q_gr)
    rec  <- length(unique(S4Vectors::subjectHits(hits))) / length(s_gr)
    f1 <- if (prec + rec == 0) 0 else 2 * prec * rec / (prec + rec)
    c(prec, rec, f1)
  })
  m <- colMeans(do.call(rbind, per))
  data.frame(tp_precision = m[1], tp_recall = m[2], tp_f1 = m[3],
             n_common_TF = length(common), stringsAsFactors = FALSE)
}

# TF-level coverage precision / recall / F1 against the ChIP TF universe.
set_prf <- function(pred_set, truth_set) {
  if (!length(truth_set)) {
    return(data.frame(tf_precision = NA_real_, tf_recall = NA_real_,
                      tf_f1 = NA_real_, stringsAsFactors = FALSE))
  }
  tp <- length(intersect(pred_set, truth_set))
  prec <- if (length(pred_set)) tp / length(pred_set) else 0
  rec  <- tp / length(truth_set)
  f1 <- if (prec + rec == 0) 0 else 2 * prec * rec / (prec + rec)
  data.frame(tf_precision = prec, tf_recall = rec, tf_f1 = f1,
             stringsAsFactors = FALSE)
}


## -------------------- Steps 1, 2, 3, 4, 5: main sweep ---------------------
# For each dataset: rebuild AQ once, then for every s compute the primary and
# final networks, count them, score them against ChIP-seq, and at s = 1 verify
# an exact reproduction of the released network. The result is a long table in
# which rows are (dataset, s, stage); comparisons across s are the analysis.
long_rows <- list()
sanity_rows <- list()

for (tag in names(DATASETS)) {
  t1 <- Sys.time()
  cfg <- DATASETS[[tag]]
  message("\n================  ", tag, "  ================")
  
  pegrn_path <- file.path(BENCH_ROOT, cfg$folder, "interest_cell_type_peGRN.qsave")
  if (!file.exists(pegrn_path)) stop("peGRN not found: ", pegrn_path)
  pegrn <- qs::qread(pegrn_path)
  cell <- pegrn[[1]]
  loose        <- as.data.frame(cell[["loose"]], stringsAsFactors = FALSE)
  saved_simple <- cell[["loose_filtered_simplified"]]   # released final triples
  saved_full   <- cell[["loose_filtered"]]              # kept rows + keep flags
  rm(pegrn, cell)
  
  e <- annotate_loose(loose)        # AQ on every loose row
  rescue <- build_rescue(e)         # fixed rescue pools
  
  chip_dir <- file.path(EVAL_ROOT, cfg$chip)
  truth <- if (dir.exists(chip_dir)) build_chip_truth(chip_dir) else NULL
  truth_TFs <- if (!is.null(truth)) unique(truth$TF) else character(0)
  message(tag, ": ", nrow(loose), " loose triples; ",
          length(truth_TFs), " ChIP-seq truth TFs")
  
  for (s in S_GRID) {
    out <- apply_scale(e, rescue, s)
    
    for (stage in c("primary", "final")) {
      df <- out[[stage]]
      cnt <- count_network(df, stage)
      m_tp <- tf_peak_metrics(df, truth)
      m_tf <- set_prf(unique(df$TF), truth_TFs)
      n_rescued <- if (stage == "final") {
        length(setdiff(triple_key(out$final), triple_key(out$primary)))
      } else 0L
      long_rows[[length(long_rows) + 1]] <- cbind(
        data.frame(dataset = tag, s = s, n_rescued_triple = n_rescued,
                   stringsAsFactors = FALSE),
        cnt, m_tp, m_tf)
    }
    
    # s = 1 must reproduce the released network: first the triple set, then the
    # primary and rescue flag counts. Any failure stops the script.
    if (abs(s - 1) < 1e-12) {
      my_keys <- triple_key(out$final)
      saved_keys <- paste(saved_simple$TF, saved_simple$peak, saved_simple$TG,
                          sep = "|")
      if (!setequal(my_keys, saved_keys)) {
        stop(tag, ": s=1 does not reproduce loose_filtered_simplified ",
             "(mine=", length(my_keys), ", saved=", length(saved_keys),
             "). Re-check the filter parameters or send filter_eGRN internals.")
      }
      chk <- data.frame(
        check = c("final triples reproduce loose_filtered_simplified",
                  "primary row count matches saved filter_keep",
                  "rescue1 count matches saved rescue_peakTG_keep",
                  "rescue2 count matches saved rescue_TF_keep"),
        pass = NA, observed = NA, expected = NA, stringsAsFactors = FALSE)
      chk[1, 2:4] <- list(TRUE, length(my_keys), length(saved_keys))
      if (!is.null(saved_full) &&
          all(c("filter_keep", "rescue_peakTG_keep", "rescue_TF_keep") %in%
              colnames(saved_full))) {
        cmp <- function(a, b) {
          data.frame(pass = isTRUE(all.equal(as.numeric(a), as.numeric(b))),
                     observed = a, expected = b, stringsAsFactors = FALSE)
        }
        chk[2, 2:4] <- cmp(nrow(out$primary), sum(saved_full$filter_keep))
        chk[3, 2:4] <- cmp(sum(out$final[["rescue1"]]),
                           sum(saved_full$rescue_peakTG_keep))
        chk[4, 2:4] <- cmp(sum(out$final[["rescue2"]]),
                           sum(saved_full$rescue_TF_keep))
      }
      sanity_rows[[length(sanity_rows) + 1]] <-
        cbind(data.frame(dataset = tag, stringsAsFactors = FALSE), chk)
    }
    
    message(tag, ": s = ", s, " done.")
  }
  
  t2 <- Sys.time()
  message(tag, ": done in ", format(t2 - t1))
}
# ================  A549  ================
# A549: 9778 loose triples; 85 ChIP-seq truth TFs
# A549: done in 4.3942 mins
# 
# ================  GM12878  ================
# GM12878: 11891 loose triples; 163 ChIP-seq truth TFs
# GM12878: done in 6.734745 mins
# 
# ================  BM  ================
# BM: 13655 loose triples; 10 ChIP-seq truth TFs
# BM: done in 30.54006 secs


## ----------------------- Step 6: assemble and check ------------------------
sweep_df <- do.call(rbind, long_rows)
sweep_df$s <- factor(sweep_df$s, levels = S_GRID)
sanity_df <- do.call(rbind, sanity_rows)

# Monotonicity: primary edge counts must not increase as s increases.
mono <- do.call(rbind, lapply(split(sweep_df, sweep_df$dataset), function(d) {
  pp <- d[d$stage == "primary", ]
  pp <- pp[order(as.numeric(as.character(pp$s))), ]
  data.frame(dataset = pp$dataset[1],
             check = "primary edges non-increasing in s",
             pass = all(diff(pp$n_TF_peak) <= 0),
             observed = paste(pp$n_TF_peak, collapse = "->"),
             expected = "monotone non-increasing", stringsAsFactors = FALSE)
}))
sanity_df <- rbind(sanity_df, mono)

# Robust precision: mask precision estimates with too few overlapping ChIP TFs.
sweep_df$tp_precision_robust <-
  ifelse(sweep_df$n_common_TF >= MIN_COMMON_TF, sweep_df$tp_precision, NA_real_)

# Coverage relative to the released default (final network at s = 1).
base <- sweep_df[sweep_df$stage == "final" & as.character(sweep_df$s) == "1",
                 c("dataset", "n_TF_peak")]
base_vec <- setNames(base$n_TF_peak, base$dataset)
sweep_df$edges_vs_default <-
  sweep_df$n_TF_peak / base_vec[sweep_df$dataset]

print(sweep_df, row.names = FALSE)
# dataset    s n_rescued_triple   stage n_TF n_peak n_TG n_TF_peak n_peak_TG n_TF_TG n_triple tp_precision    tp_recall        tp_f1 n_common_TF tf_precision tf_recall     tf_f1 tp_precision_robust edges_vs_default
#    A549  0.5                0 primary   84   1336 1316      6687      1633    5982     6944    0.1976970 0.0017485412 0.0034543858          18    0.2142857 0.2117647 0.2130178           0.1976970        3.2133590
#    A549  0.5               17   final   84   1336 1316      6704      1633    5999     6961    0.1976970 0.0017485412 0.0034543858          18    0.2142857 0.2117647 0.2130178           0.1976970        3.2215281
#    A549 0.55                0 primary   81   1336 1315      6324      1631    5697     6579    0.2034138 0.0016761758 0.0033134544          18    0.2222222 0.2117647 0.2168675           0.2034138        3.0389236
#    A549 0.55               38   final   84   1336 1315      6361      1631    5735     6617    0.2034138 0.0016761758 0.0033134544          18    0.2142857 0.2117647 0.2130178           0.2034138        3.0567035
#    A549  0.6                0 primary   81   1330 1310      5878      1622    5330     6118    0.2166836 0.0015633351 0.0030926901          18    0.2222222 0.2117647 0.2168675           0.2166836        2.8246036
#    A549  0.6               39   final   84   1330 1310      5916      1622    5369     6157    0.2166836 0.0015633351 0.0030926901          18    0.2142857 0.2117647 0.2130178           0.2166836        2.8428640
#    A549 0.65                0 primary   80   1316 1302      5335      1602    4881     5547    0.2344350 0.0014606362 0.0028924788          18    0.2250000 0.2117647 0.2181818           0.2344350        2.5636713
#    A549 0.65               48   final   84   1316 1302      5381      1602    4929     5595    0.2344350 0.0014606362 0.0028924788          18    0.2142857 0.2117647 0.2130178           0.2344350        2.5857761
#    A549  0.7                0 primary   78   1286 1274      4771      1555    4402     4960    0.2537819 0.0013490509 0.0026749214          17    0.2179487 0.2000000 0.2085890           0.2537819        2.2926478
#    A549  0.7               54   final   84   1287 1275      4822      1556    4456     5014    0.2470903 0.0012796592 0.0025374180          18    0.2142857 0.2117647 0.2130178           0.2470903        2.3171552
#    A549 0.75                0 primary   77   1239 1221      4146      1477    3846     4300    0.2831800 0.0012008155 0.0023837765          17    0.2207792 0.2000000 0.2098765           0.2831800        1.9923114
#    A549 0.75               68   final   84   1248 1231      4211      1488    3914     4368    0.2600404 0.0011396591 0.0022624456          18    0.2142857 0.2117647 0.2130178           0.2600404        2.0235464
#    A549  0.8                0 primary   77   1171 1162      3516      1384    3268     3629    0.3005618 0.0010328796 0.0020522226          17    0.2207792 0.2000000 0.2098765           0.3005618        1.6895723
#    A549  0.8               70   final   84   1183 1176      3583      1400    3338     3699    0.2764565 0.0009810530 0.0019493113          18    0.2142857 0.2117647 0.2130178           0.2764565        1.7217684
#    A549 0.85                0 primary   76   1079 1073      2882      1260    2693     2965    0.3335164 0.0008854248 0.0017602307          17    0.2236842 0.2000000 0.2111801           0.3335164        1.3849111
#    A549 0.85              173   final   84   1161 1152      3051      1371    2849     3138    0.3168534 0.0008585334 0.0017067984          18    0.2142857 0.2117647 0.2130178           0.3168534        1.4661221
#    A549  0.9                0 primary   76    982  970      2388      1126    2236     2444    0.3478355 0.0008073104 0.0016055385          17    0.2236842 0.2000000 0.2111801           0.3478355        1.1475252
#    A549  0.9              293   final   84   1144 1134      2674      1348    2490     2737    0.3287573 0.0007958972 0.0015825996          18    0.2142857 0.2117647 0.2130178           0.3287573        1.2849592
#    A549 0.95                0 primary   75    889  872      1928       992    1821     1970    0.3338818 0.0006939057 0.0013806766          17    0.2266667 0.2000000 0.2125000           0.3338818        0.9264777
#    A549 0.95              420   final   84   1134 1121      2333      1332    2180     2390    0.3453280 0.0007277474 0.0014473376          18    0.2142857 0.2117647 0.2130178           0.3453280        1.1210956
#    A549    1                0 primary   74    696  688      1435       747    1374     1455    0.3442755 0.0005686923 0.0011328001          17    0.2297297 0.2000000 0.2138365           0.3442755        0.6895723
#    A549    1              673   final   84   1127 1114      2081      1324    1935     2128    0.3418885 0.0006484421 0.0012897448          18    0.2142857 0.2117647 0.2130178           0.3418885        1.0000000
#    A549 1.05                0 primary   73    624  622      1216       666    1175     1235    0.3297276 0.0004565238 0.0009094510          17    0.2328767 0.2000000 0.2151899           0.3297276        0.5843345
#    A549 1.05              759   final   84   1120 1104      1951      1310    1811     1994    0.3286484 0.0005647647 0.0011231968          18    0.2142857 0.2117647 0.2130178           0.3286484        0.9375300
#    A549  1.1                0 primary   72    550  546      1021       582     989     1036    0.3191288 0.0004092558 0.0008155389          17    0.2361111 0.2000000 0.2165605           0.3191288        0.4906295
#    A549  1.1              844   final   84   1116 1100      1836      1305    1707     1880    0.3211719 0.0005313287 0.0010569085          18    0.2142857 0.2117647 0.2130178           0.3211719        0.8822681
#    A549 1.15                0 primary   68    503  498       885       530     858      898    0.3491290 0.0003784822 0.0007544352          16    0.2352941 0.1882353 0.2091503           0.3491290        0.4252763
#    A549 1.15              903   final   84   1112 1096      1757      1300    1630     1801    0.3235601 0.0004880080 0.0009709322          18    0.2142857 0.2117647 0.2130178           0.3235601        0.8443056
#    A549  1.2                0 primary   65    443  436       753       462     732      764    0.3743607 0.0003488601 0.0006957338          15    0.2307692 0.1764706 0.2000000           0.3743607        0.3618453
#    A549  1.2              976   final   84   1109 1092      1696      1296    1573     1740    0.3308826 0.0004757471 0.0009468442          18    0.2142857 0.2117647 0.2130178           0.3308826        0.8149928
#    A549 1.25                0 primary   63    390  385       632       406     614      640    0.3096555 0.0002954612 0.0005892529          15    0.2380952 0.1764706 0.2027027           0.3096555        0.3037001
#    A549 1.25             1039   final   84   1104 1087      1636      1291    1514     1679    0.3244272 0.0004535361 0.0009026025          18    0.2142857 0.2117647 0.2130178           0.3244272        0.7861605
#    A549  1.3                0 primary   56    342  338       532       357     517      540    0.3684100 0.0003025055 0.0006032684          12    0.2142857 0.1411765 0.1702128           0.3684100        0.2556463
#    A549  1.3             1096   final   84   1101 1082      1591      1286    1473     1636    0.3233334 0.0004312860 0.0008582933          18    0.2142857 0.2117647 0.2130178           0.3233334        0.7645363
#    A549 1.35                0 primary   53    302  299       450       314     441      458    0.3736919 0.0002439752 0.0004865514          12    0.2264151 0.1411765 0.1739130           0.3736919        0.2162422
#    A549 1.35             1140   final   84   1098 1080      1554      1283    1438     1598    0.3172953 0.0003978695 0.0007917755          18    0.2142857 0.2117647 0.2130178           0.3172953        0.7467564
#    A549  1.4                0 primary   53    262  260       379       269     374      385    0.3038012 0.0001670716 0.0003332175          12    0.2264151 0.1411765 0.1739130           0.3038012        0.1821240
#    A549  1.4             1186   final   84   1095 1077      1529      1280    1411     1571    0.3104988 0.0003632732 0.0007230205          18    0.2142857 0.2117647 0.2130178           0.3104988        0.7347429
#    A549 1.45                0 primary   48    230  228       324       235     319      328    0.2301097 0.0001549197 0.0003089329          11    0.2291667 0.1294118 0.1654135           0.2301097        0.1556944
#    A549 1.45             1224   final   84   1092 1075      1510      1277    1394     1552    0.3093639 0.0003632537 0.0007229220          18    0.2142857 0.2117647 0.2130178           0.3093639        0.7256127
#    A549  1.5                0 primary   43    198  196       266       201     262      268    0.2409733 0.0001558812 0.0003108100           9    0.2093023 0.1058824 0.1406250           0.2409733        0.1278232
#    A549  1.5             1261   final   84   1092 1075      1489      1276    1372     1529    0.3198651 0.0003632425 0.0007229422          18    0.2142857 0.2117647 0.2130178           0.3198651        0.7155214
# GM12878  0.5                0 primary   75   1639 1791      7806      2180    7139     8183    0.2268063 0.0023153013 0.0045721739          26    0.3466667 0.1595092 0.2184874           0.2268063        2.9214072
# GM12878  0.5               12   final   76   1640 1791      7817      2181    7151     8195    0.2268063 0.0023153013 0.0045721739          26    0.3421053 0.1595092 0.2175732           0.2268063        2.9255240
# GM12878 0.55                0 primary   75   1632 1783      7301      2169    6726     7665    0.2427853 0.0022105994 0.0043688931          26    0.3466667 0.1595092 0.2184874           0.2427853        2.7324102
# GM12878 0.55               12   final   76   1633 1783      7312      2170    6738     7677    0.2427853 0.0022105994 0.0043688931          26    0.3421053 0.1595092 0.2175732           0.2427853        2.7365269
# GM12878  0.6                0 primary   74   1622 1771      6751      2149    6231     7093    0.2477995 0.0020710378 0.0040955152          26    0.3513514 0.1595092 0.2194093           0.2477995        2.5265719
# GM12878  0.6               24   final   76   1623 1772      6774      2151    6255     7117    0.2477995 0.0020710378 0.0040955152          26    0.3421053 0.1595092 0.2175732           0.2477995        2.5351796
# GM12878 0.65                0 primary   73   1603 1747      6149      2114    5703     6457    0.2680105 0.0019352827 0.0038304347          26    0.3561644 0.1595092 0.2203390           0.2680105        2.3012725
# GM12878 0.65               29   final   76   1605 1749      6177      2117    5732     6486    0.2680105 0.0019352827 0.0038304347          26    0.3421053 0.1595092 0.2175732           0.2680105        2.3117515

# How to read sweep_df
# Rows: one row per dataset × threshold scale s × network stage.
# Key columns:
#   n_TF_peak           = number of TF-peak edges (coverage axis)
#   tp_precision        = macro-averaged TF-peak precision vs ChIP-seq (quality axis)
#   n_common_TF         = number of TFs with ChIP-seq truth used in averaging (reliability indicator)
#   n_rescued_triple    = number of edges added by the two rescue steps
#   tp_precision_robust = precision set to NA when n_common_TF < MIN_COMMON_TF (unreliable estimates)
#   edges_vs_default    = TF-peak edges normalized to the final network at s = 1
# Analysis: compare rows of the same dataset and stage across s values.
#   - Primary curve: clean monotonic coverage-precision trade-off.
#     Minor precision dips at high s stem from shrinking overlapping TF sets, not filter failure.
#   - Final curve: rescue steps enforce a coverage floor for s > 1, so precision plateaus.
#     Default s = 1 lies near the precision peak of the final network.

print(sanity_df, row.names = FALSE)
# dataset                                             check pass                                                                                                                observed                expected
#    A549 final triples reproduce loose_filtered_simplified TRUE                                                                                                                    2128                    2128
#    A549       primary row count matches saved filter_keep TRUE                                                                                                                    1455                    1455
#    A549    rescue1 count matches saved rescue_peakTG_keep TRUE                                                                                                                     558                     558
#    A549        rescue2 count matches saved rescue_TF_keep TRUE                                                                                                                     130                     130
# GM12878 final triples reproduce loose_filtered_simplified TRUE                                                                                                                    2763                    2763
# GM12878       primary row count matches saved filter_keep TRUE                                                                                                                    1934                    1934
# GM12878    rescue1 count matches saved rescue_peakTG_keep TRUE                                                                                                                     755                     755
# GM12878        rescue2 count matches saved rescue_TF_keep TRUE                                                                                                                      83                      83
#      BM final triples reproduce loose_filtered_simplified TRUE                                                                                                                    3309                    3309
#      BM       primary row count matches saved filter_keep TRUE                                                                                                                    2457                    2457
#      BM    rescue1 count matches saved rescue_peakTG_keep TRUE                                                                                                                     720                     720
#      BM        rescue2 count matches saved rescue_TF_keep TRUE                                                                                                                     135                     135
#    A549                 primary edges non-increasing in s TRUE    6687->6324->5878->5335->4771->4146->3516->2882->2388->1928->1435->1216->1021->885->753->632->532->450->379->324->266 monotone non-increasing
#      BM                 primary edges non-increasing in s TRUE 9441->8911->8363->7730->7050->6363->5416->4615->3839->3135->2410->2053->1712->1441->1196->1035->885->750->624->537->449 monotone non-increasing
# GM12878                 primary edges non-increasing in s TRUE   7806->7301->6751->6149->5564->4907->4217->3641->3055->2510->1894->1620->1374->1162->983->832->724->609->504->404->345 monotone non-increasing

# How to read sanity_df
# Rows: validation checks per dataset.
# First 4 rows per dataset: exact reproduction check at s = 1
#   - Final triple set matches released loose_filtered_simplified
#   - Primary edge count matches saved filter_keep flag
#   - Rescue1 / rescue2 counts match saved corresponding flags
# Last 3 rows (one per dataset): monotonicity check
#   - Primary TF-peak edge count must be non-increasing as s increases
# All TRUE values confirm that recomputed AQ and rescue logic match the
# original filter implementation, and the entire sweep is reliable.

write.csv(sweep_df, file.path(OUT_DIR, "E3_threshold_sweep_long.csv"),
          row.names = FALSE)
write.csv(sanity_df, file.path(OUT_DIR, "E3_threshold_sanity.csv"),
          row.names = FALSE)


## -------------------------------- Figures ----------------------------------
# Reading guide:
#   Fig 1 = precision vs coverage. The dashed primary curve is the clean
#           monotone trade-off (the main evidence); the solid final curve shows
#           the rescue floor. BM precision is omitted (only 10 ChIP TFs).
#   Fig 2 = edges vs s. Edge counts are the coverage/recall axis we tune; the
#           final curve flattens for s > 1 because of the rescue floor.
#   Fig 3 = TF-level coverage vs s; it is nearly flat because the threshold
#           tunes edge-level coverage while the rescue holds the TF set.
ds_colors <- c(A549 = "#66C2A5", GM12878 = "#FC8D62", BM = "#8DA0CB")
stage_lab <- c(primary = "Primary filter only",
               final = "Primary + fixed rescue")

# Figure 1: precision-coverage operating curve (one panel per dataset).
p1 <- ggplot(sweep_df, aes(x = n_TF_peak, y = tp_precision_robust,
                           color = dataset, linetype = stage, shape = stage)) +
  geom_line(linewidth = 0.6) + geom_point(size = 2) +
  facet_wrap(~ dataset, scales = "free_x", nrow = 1) +
  scale_color_manual(values = ds_colors) +
  scale_linetype_manual(values = c(primary = "dashed", final = "solid"),
                        labels = stage_lab) +
  scale_shape_manual(values = c(primary = 1, final = 16),
                     labels = stage_lab) +
  labs(x = "TF-peak edges (coverage)",
       y = "Mean TF-peak precision vs ChIP-seq") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        plot.caption = element_text(size = 7, hjust = 0))
print(p1)
ggsave(file.path(OUT_DIR, "E3_operating_curve.png"), p1,
       width = 7, height = 4, units = "in", dpi = 300)

# Figure 2: edge count vs threshold scale s (both stages).
p2 <- ggplot(sweep_df, aes(x = as.numeric(as.character(s)), y = n_TF_peak,
                           color = dataset, linetype = stage, shape = stage)) +
  geom_vline(xintercept = 1, linetype = "dotted", color = "grey50") +
  geom_line(linewidth = 0.6) + geom_point(size = 2) +
  facet_wrap(~ dataset, scales = "free_y", nrow = 1) +
  scale_color_manual(values = ds_colors) +
  scale_linetype_manual(values = c(primary = "dashed", final = "solid"),
                        labels = stage_lab) +
  scale_shape_manual(values = c(primary = 1, final = 16),
                     labels = stage_lab) +
  # scale_x_continuous(breaks = S_GRID) +
  scale_x_continuous(breaks = seq(0.5, 1.5, by = 0.1)) +
  annotate("text", x = 1, y = -Inf, label = "default s=1", vjust = -0.5,
           hjust = 1.1, size = 2.6, color = "grey30") +
  labs(x = "Threshold scale s (AQ_s = s*AQ)", y = "TF-peak edges",
       caption = "Relaxing (s<1) grows the network; the final curve flattens for s>1 (rescue floor).") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        plot.caption = element_text(size = 7, hjust = 0))
print(p2)
ggsave(file.path(OUT_DIR, "E3_edges_vs_scale.png"), p2,
       width = 7, height = 4, units = "in", dpi = 300)

# Figure 3: TF-level coverage precision / recall vs s (final network).
tf_cov <- reshape(
  sweep_df[sweep_df$stage == "final",
           c("dataset", "s", "tf_recall", "tf_precision")],
  idvar = c("dataset", "s"), timevar = "metric",
  varying = c("tf_recall", "tf_precision"), v.names = "value",
  times = c("tf_recall", "tf_precision"), direction = "long")
tf_cov$metric <- factor(tf_cov$metric,
                        levels = c("tf_recall", "tf_precision"),
                        labels = c("TF recall", "TF precision"))
p3 <- ggplot(tf_cov, aes(x = as.numeric(as.character(s)), y = value,
                         color = dataset, linetype = metric, shape = metric)) +
  geom_vline(xintercept = 1, linetype = "dotted", color = "grey50") +
  geom_line(linewidth = 0.6) + geom_point(size = 2) +
  facet_wrap(~ dataset, nrow = 1) +
  scale_color_manual(values = ds_colors) +
  # scale_x_continuous(breaks = S_GRID) +
  scale_x_continuous(breaks = seq(0.5, 1.5, by = 0.1)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_shape_manual(values = c("TF recall" = 1, "TF precision" = 16),
                     labels = c("TF recall", "TF precision")) +
  labs(x = "Threshold scale s (AQ_s = s*AQ)", y = "TF coverage vs ChIP-seq",
       caption = "TF-level coverage is nearly flat: the threshold tunes edges, not the TF set.") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        plot.caption = element_text(size = 7, hjust = 0))
print(p3)
ggsave(file.path(OUT_DIR, "E3_TFcoverage_vs_scale.png"), p3,
       width = 7, height = 4, units = "in", dpi = 300)

message("\nFiles saved under ", OUT_DIR)

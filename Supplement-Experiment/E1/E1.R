# ==============================================================================
# E1: Grouped permutation importance of the six evidence categories for the 
# final peGRN
# ==============================================================================


# Purpose:
#   This script answers "which data category is most decisive for the results" and, in part, parameter/feature sensitivity.
#   It quantifies how much each evidence category contributes at the TERMINAL EDGE SELECTION stage of scDERI: 
#   for each category, the corresponding raw feature(s) are permuted across triples 
#   (preserving their marginal distribution, destroying that information), 
#   the affected rank/normalized features are recomputed exactly as in Section 5.6 of scDERI_benchmark.R, 
#   and only the terminal filter (filter_eGRN) is rerun -- NO model retraining.

# Six evidence groups (permutation target, raw columns):
#   G1 motif   : FP-S               (JASPAR motif match; FP-S-R is directly thresholded)
#   G2 prior   : FG-theta_i         (ChIP-Atlas prior rank)
#   G3 dynamics: FG-theta_p,FG-L    (predicted strength -> IQ threshold; fit loss -> FG-NL)
#   G4 access  : P-NZAR,P-NZAA,P-AA (peak accessibility)
#   G5 TG_expr : G-NZER,G-NZEA,G-EA (target-gene expression)
#   G6 TF_expr : F-NZER,F-NZEA,F-EA (TF expression)

# Two condition types per group (per dataset):
#   permute: B (=20) independent permutations; 
#            edges are re-filtered and compared with the baseline edge set 
#            (Jaccard, kept fraction, added edges, edge count).
#   neutral: deterministic control; the factor range f_r of that group is set to 0,
#            forcing its adjustment factor F_m = 1 for every edge (equivalent to placing every feature at its balance value f_b = 0.2). 
#            G1 has no factor (FP-S-R is the thresholded quantity itself), so permutation IS its control.
#   G3 additionally reports two single-column decompositions (theta_p only / loss only),
#   because FG-theta_p enters the IQ threshold and rescue ordering, 
#   whereas FG-NL is a multiplicative factor.

# Steps:
#   1. Configure datasets, IQ anchor curve (fitted once), 
#      factor balances/ranges and rescue settings (identical to Section 5.6 of scDERI_benchmark.R).
#   2. For each dataset: load interest_cell_type_peGRN.qsave; 
#      rebuild the theta_i rank reference; 
#      rerun the terminal filter on x$loose and verify it reproduces x$loose_filtered_simplified exactly 
#      (sanity check, stops on mismatch).
#   3. Run B permutations per group, the neutral controls and the G3 decomposition.
#   4. Combine everything into a long table and write CSV files.
#   5. Produce a Supplementary-Figure-ready plot (median + IQR of Jaccard/kept_frac,
#      faceted by dataset; neutral controls shown as diamonds).

# Outputs (written to OUT_DIR):
#   E1_permutation_importance.csv      long table: dataset, group, type, b,
#                                      jaccard, kept_frac, added, n_edge
#   E1_summary.csv                     median/IQR per dataset x group x type
#   E1_permutation_importance.pdf/.png pointrange figure (Jaccard and kept_frac)


## ---------------------------- 0. User configuration ---------------------------------
library(qs)
library(minpack.lm)
library(ggplot2)

od <- "/mnt/i/scDERI_Paper/3-benchmark/8-scDERI" # benchmark root
OUT_DIR <- "/mnt/i/scDERI_Paper/Supplement-Experiment/E1" # all outputs go here
# B <- 20 # permutations per group # 43 secs
B <- 1000 # permutations per group # 25.2 mins
SEED_BASE <- 123
DECOMPOSE_DYNAMICS <- TRUE # extra theta_p-only / loss-only runs
VERBOSE <- FALSE # print per-step filter messages

# Dataset registry. The three benchmarks share the same 5.6 settings; 
# set enabled=FALSE to skip. 
# fill its path and settings from its own Section 5.6 run, then set enabled = TRUE.
DATASETS <- list(
  A549 = list(
    path = file.path(od, "SNAREseq2_hg38_A549"),
    igrn = file.path(od, "SNAREseq2_hg38_A549", "3 get iGRN", "interest_cell_type_iGRN.rds"),
    rescue_peakTG = TRUE, peakTG_q = 0.25, min_TF = 1,
    rescue_TF = TRUE, TF_q = 0.0, min_peakTG = 7, enabled = TRUE
  ),
  GM12878 = list(
    path = file.path(od, "SNAREseq2_hg38_GM12878"),
    igrn = file.path(od, "SNAREseq2_hg38_GM12878", "3 get iGRN", "interest_cell_type_iGRN.rds"),
    rescue_peakTG = TRUE, peakTG_q = 0.25, min_TF = 1,
    rescue_TF = TRUE, TF_q = 0.0, min_peakTG = 7, enabled = TRUE
  ),
  BM = list(
    path = file.path(od, "10X_hg38_Human_bone_marrow_with_T_cell_depleted"),
    igrn = file.path(od, "10X_hg38_Human_bone_marrow_with_T_cell_depleted", "3 get iGRN", "interest_cell_type_iGRN.rds"),
    rescue_peakTG = TRUE, peakTG_q = 0.25, min_TF = 1,
    rescue_TF = TRUE, TF_q = 0.0, min_peakTG = 7, enabled = TRUE
  )
)


## ---------------------------- 1. Fixed Section 5.6 settings -------------------------
# Eleven adjustment factors (order irrelevant): balance f_b, range f_r, attribute.
# These values are copied verbatim from Section 5.6 of scDERI_benchmark.R; 
# keep them in sync with the released pipeline.
factor_balance <- list(
  `FG-NL` = 0.2, `FG-theta_i-R` = 0.2,
  `F-NZER-R` = 0.2, `F-NZEA-R` = 0.2, `F-EA-R` = 0.2,
  `P-NZAR-R` = 0.2, `P-NZAA-R` = 0.2, `P-AA-R` = 0.2,
  `G-NZER-R` = 0.2, `G-NZEA-R` = 0.2, `G-EA-R` = 0.2
)
factor_range <- list(
  `FG-NL` = 0.43, `FG-theta_i-R` = 0.64,
  `F-NZER-R` = 0.28, `F-NZEA-R` = 0.33, `F-EA-R` = 0.22,
  `P-NZAR-R` = 0.19, `P-NZAA-R` = 0.67, `P-AA-R` = 0.13,
  `G-NZER-R` = 0.21, `G-NZEA-R` = 0.12, `G-EA-R` = 0.26
)
factor_attribute <- list(
  `FG-NL` = "penalty", `FG-theta_i-R` = "relaxation",
  `F-NZER-R` = "relaxation", `F-NZEA-R` = "relaxation", `F-EA-R` = "relaxation",
  `P-NZAR-R` = "relaxation", `P-NZAA-R` = "relaxation", `P-AA-R` = "relaxation",
  `G-NZER-R` = "relaxation", `G-NZEA-R` = "relaxation", `G-EA-R` = "relaxation"
)

# Inverse quantile curve: IQ = 1 / (theta_p * a + b) + c, fitted once on the anchors.
theta_p_points <- c(0, 0.5, 1)
initial_quantile_points <- c(5.79, 4.66, 4.23)
iq_fit <- minpack.lm::nlsLM(
  initial_quantile_points ~ 1 / (theta_p_points * a + b) + c,
  start = list(a = 0.333333, b = 0.333333, c = -0.5),
  lower = c(1e-9, -Inf, -Inf), upper = c(Inf, Inf, Inf)
)
IQ_A <- as.numeric(coef(iq_fit)["a"])
IQ_B <- as.numeric(coef(iq_fit)["b"])
IQ_C <- as.numeric(coef(iq_fit)["c"])
rm(iq_fit)
message("Fitted IQ curve: 1/(theta_p*", round(IQ_A, 6), "+(", round(IQ_B, 6), "))+", round(IQ_C, 6))
# Fitted IQ curve: 1/(theta_p*0.646435+(0.397096))+3.271714

# Evidence groups: permuted raw columns and the factor columns neutralized in the control.
GROUPS <- list(
  motif    = list(raw =     "FP-S",
                  neutral = character(0)),
  prior    = list(raw =     "FG-theta_i",
                  neutral = "FG-theta_i-R"),
  dynamics = list(raw =     c("FG-theta_p", "FG-L"),
                  neutral = "FG-NL"),
  access   = list(raw =     c("P-NZAR", "P-NZAA", "P-AA"),
                  neutral = c("P-NZAR-R", "P-NZAA-R", "P-AA-R")),
  TG_expr  = list(raw =     c("G-NZER", "G-NZEA", "G-EA"),
                  neutral = c("G-NZER-R", "G-NZEA-R", "G-EA-R")),
  TF_expr  = list(raw =     c("F-NZER", "F-NZEA", "F-EA"),
                  neutral = c("F-NZER-R", "F-NZEA-R", "F-EA-R"))
)
GROUP_ORDER <- names(GROUPS)

PEAK_RANK_RAW <- c("P-NZAR", "P-NZAA", "P-AA")
TG_RANK_RAW   <- c("G-NZER", "G-NZEA", "G-EA")
TF_RANK_RAW   <- c("F-NZER", "F-NZEA", "F-EA")
RANK_RAW <- c(PEAK_RANK_RAW, TG_RANK_RAW, TF_RANK_RAW)


## ---------------------------- 2. Helper functions -----------------------------------
# Rank to (0,1] using minimum ranks, scaled by the largest rank; NAs are preserved.
rank01 <- function(x) {
  r <- rank(x, ties.method = "min", na.last = "keep")
  r / max(r, na.rm = TRUE)
}

# Min-max scale to [0,1]; constant columns return all zeros.
minmax01 <- function(x) {
  lo <- min(x, na.rm = TRUE); hi <- max(x, na.rm = TRUE)
  if (hi - lo <= .Machine$double.eps) rep(0, length(x)) else (x - lo) / (hi - lo)
}

say <- function(verbose, ...) if (isTRUE(verbose)) message(...)

# Build the catalog-wide theta_i -> normalized rank reference (positive theta_i only).
build_theta_i_lookup <- function(igrn_path) {
  ig <- readRDS(igrn_path)[[1]]
  tv <- sort(unlist(ig)); tv <- tv[tv > 0]
  data.frame(value = tv,
             rank = rank(tv, ties.method = "min") / length(tv),
             stringsAsFactors = FALSE)
}

# Independently permute the raw columns of one group, 
# then recompute every derived rank/normalized feature exactly as Sections 5.5-5.6 do. 
# theta_p has no derived column:
# it flows into IQ and the rescue ordering inside the terminal filter.
permute_and_rerank <- function(tri, raw_cols, theta_lookup, verbose = FALSE) {
  tri2 <- tri
  for (cc in raw_cols) {
    say(verbose, "  permuting ", cc)
    tri2[[cc]] <- sample(tri2[[cc]])
  }
  if ("FP-S" %in% raw_cols) {
    say(verbose, "  recomputing FP-S-R")
    tri2[, "FP-S-R"] <- rank01(tri2[, "FP-S"])
  }
  if ("FG-theta_i" %in% raw_cols) {
    say(verbose, "  recomputing FG-theta_i-R from the catalog reference")
    rk <- theta_lookup$rank[match(tri2[, "FG-theta_i"], theta_lookup$value)]
    tri2[, "FG-theta_i-R"] <- ifelse(is.na(rk), 0, rk)
  }
  if ("FG-L" %in% raw_cols) {
    say(verbose, "  recomputing FG-NL")
    tri2[, "FG-NL"] <- minmax01(tri2[, "FG-L"])
  }
  for (cc in intersect(raw_cols, RANK_RAW)) {
    say(verbose, "  recomputing ", cc, "-R")
    tri2[, paste0(cc, "-R")] <- rank01(tri2[, cc])
  }
  tri2
}

# Terminal peGRN filter, reproduced exactly from Section 5.6 of scDERI_benchmark.R.
# cfg carries the per-dataset rescue settings, frange the eleven factor ranges 
# (set to 0 in the neutral control), and verbose toggles step messages.
run_filter_only <- function(tri, cfg, frange = factor_range, verbose = FALSE) {
  e <- tri
  # Remove any columns left by a previous call so the rerun starts from a clean table.
  drop_cols <- intersect(
    c("peak-TG", "IQ", grep("-F$", names(e), value = TRUE),
      "AQ", "filter_keep", "rescue_peakTG_keep", "rescue_TF_keep", "keep"),
    names(e))
  if (length(drop_cols) > 0) e[, drop_cols] <- NULL
  e[, "peak-TG"] <- paste0(e$peak, "_to_", e$TG)
  rownames(e) <- seq_len(nrow(e))
  
  # Per-edge initial quantile from predicted strength.
  e[, "IQ"] <- 1 / (e[, "FG-theta_p"] * IQ_A + IQ_B) + IQ_C
  
  # Feature-specific factors: relaxation F = 1 + f_r*(f_b - x), penalty F = 1 + f_r*(x - f_b).
  for (ff in names(factor_balance)) {
    fb <- factor_balance[[ff]]; fv <- frange[[ff]]
    sign <- (-1) ^ as.numeric(factor_attribute[[ff]] == "penalty")
    e[, paste0(ff, "-F")] <- (1 + sign * fb * fv) - sign * e[, ff] * fv
  }
  
  # Adjusted quantile = IQ * product of factors, clamped to [0,1].
  aq <- as.matrix(e[, c("IQ", paste0(names(factor_balance), "-F"))])
  aq[aq <= 0] <- .Machine$double.eps
  e[, "AQ"] <- pmin(pmax(exp(rowSums(log(aq))), 0), 1)
  
  e$filter_keep <- e[, "FP-S-R"] >= e[, "AQ"]
  e$rescue_peakTG_keep <- FALSE
  e$rescue_TF_keep <- FALSE
  e$keep <- e$filter_keep
  
  # One priority ordering, reused by both rescue steps.
  ord <- order(-e$`FP-S-R`, -e$`P-AA-R`, -e$`P-NZAR-R`,
               -e$`FG-theta_i-R`, -e$`FG-theta_p`, e$`FG-NL`)
  eo <- e[ord, ]
  
  # Rescue step 1: keep at least min_TF TF candidate(s) per peak-TG above a quality floor.
  if (cfg$rescue_peakTG) {
    say(verbose, "  rescue peak-TG (min ", cfg$min_TF, ", q=", cfg$peakTG_q, ")")
    eo$peakTG_rank <- ave(seq_len(nrow(eo)), eo$`peak-TG`, FUN = seq_along)
    cand <- eo[eo$peakTG_rank <= cfg$min_TF, ]
    qt <- quantile(cand$`FP-S-R`, cfg$peakTG_q)
    cand <- cand[cand$`FP-S-R` >= qt, ]
    n_keep <- tapply(e$keep, e$`peak-TG`, sum)
    n_cand <- n_keep[cand$`peak-TG`]
    need <- cand$peakTG_rank <= (cfg$min_TF - n_cand) & n_cand < cfg$min_TF
    e[rownames(cand)[need], "rescue_peakTG_keep"] <- TRUE
    rm(cand, qt, n_keep, n_cand, need)
  }
  e$keep <- e$filter_keep | e$rescue_peakTG_keep
  
  # Rescue step 2: keep at least min_peakTG peak-TG candidate(s) per TF above a quality floor.
  if (cfg$rescue_TF) {
    say(verbose, "  rescue TF (min ", cfg$min_peakTG, ", q=", cfg$TF_q, ")")
    eo$TF_rank <- ave(seq_len(nrow(eo)), eo$TF, FUN = seq_along)
    cand <- eo[eo$TF_rank <= cfg$min_peakTG, ]
    qt <- quantile(cand$`FP-S-R`, cfg$TF_q)
    cand <- cand[cand$`FP-S-R` >= qt, ]
    n_keep <- tapply(e$keep, e$TF, sum)
    n_cand <- n_keep[cand$TF]
    need <- cand$TF_rank <= (cfg$min_peakTG - n_cand) & n_cand < cfg$min_peakTG
    e[rownames(cand)[need], "rescue_TF_keep"] <- TRUE
    rm(cand, qt, n_keep, n_cand, need)
  }
  rm(ord, eo)
  
  e$keep <- e$filter_keep | e$rescue_peakTG_keep | e$rescue_TF_keep
  unique(e[e$keep, c("TF", "peak", "TG")])
}

edge_set <- function(tri) unique(paste(tri$TF, tri$peak, tri$TG, sep = "|"))

# Jaccard, kept fraction, added edges and final edge count relative to the baseline E0.
score_one <- function(E1, E0) {
  inter <- sum(E1 %in% E0)
  data.frame(jaccard = inter / length(union(E0, E1)),
             kept_frac = inter / length(E0),
             added = sum(!E1 %in% E0),
             n_edge = length(E1))
}


## ---------------------------- 3. Main loop ------------------------------------------
dir.create(OUT_DIR, showWarnings = TRUE, recursive = TRUE)
res <- list(); k <- 1
dss <- names(DATASETS)[vapply(DATASETS, function(x) isTRUE(x$enabled), logical(1))]

time1 = Sys.time()
for (ds in dss) {
  cfg <- DATASETS[[ds]]
  message("=== Dataset: ", ds, " ===")
  obj <- qs::qread(file.path(cfg$path, "interest_cell_type_peGRN.qsave"))[[1]]
  tri0 <- obj$loose
  theta_lookup <- build_theta_i_lookup(cfg$igrn)
  
  # Sanity check: rerun the filter on the untouched table and require the saved edges.
  base_keep <- run_filter_only(tri0, cfg, verbose = VERBOSE)
  E0 <- edge_set(base_keep)
  saved <- edge_set(obj$loose_filtered_simplified)
  if (!setequal(E0, saved)) {
    stop("Baseline mismatch for ", ds, ": rerun=", length(E0),
         " saved=", length(saved), ". Check filter parameters before continuing.")
  }
  message("  baseline reproduced: ", length(E0), " edges")
  
  for (g in GROUP_ORDER) {
    message("  group ", g)
    raw <- GROUPS[[g]]$raw
    
    # (a) B permutations of the group's raw columns.
    for (b in seq_len(B)) {
      set.seed(SEED_BASE + (k - 1) * 97 + b)
      tri_p <- permute_and_rerank(tri0, raw, theta_lookup, verbose = VERBOSE)
      E1 <- edge_set(run_filter_only(tri_p, cfg, verbose = VERBOSE))
      s <- score_one(E1, E0)
      res[[k]] <- data.frame(dataset = ds, group = g, type = "permute", b = b, s)
      k <- k + 1
    }
    
    # (b) Deterministic neutral control (f_r = 0), only for groups with a factor.
    if (length(GROUPS[[g]]$neutral) > 0) {
      fr_neutral <- factor_range
      fr_neutral[GROUPS[[g]]$neutral] <- 0
      E1 <- edge_set(run_filter_only(tri0, cfg, frange = fr_neutral, verbose = VERBOSE))
      res[[k]] <- data.frame(dataset = ds, group = g, type = "neutral", b = 0,
                             score_one(E1, E0))
      k <- k + 1
    }
    
    # (c) G3 decomposition: permute theta_p (IQ/rescue ordering) and loss (FG-NL) alone.
    if (DECOMPOSE_DYNAMICS && g == "dynamics") {
      for (single in c("FG-theta_p", "FG-L")) {
        for (b in seq_len(B)) {
          set.seed(SEED_BASE + (k - 1) * 97 + b)
          tri_p <- permute_and_rerank(tri0, single, theta_lookup, verbose = VERBOSE)
          E1 <- edge_set(run_filter_only(tri_p, cfg, verbose = VERBOSE))
          cond <- if (single == "FG-theta_p") "permute_theta_p_only" else "permute_loss_only"
          res[[k]] <- data.frame(dataset = ds, group = g, type = cond, b = b,
                                 score_one(E1, E0))
          k <- k + 1
        }
      }
    }
  }
  # Save an interim copy after each dataset so partial results are never lost.
  write.csv(do.call(rbind, res), file.path(OUT_DIR, "E1_permutation_importance.csv"),
            row.names = FALSE)
}
Res <- do.call(rbind, res)
for (cc in c("jaccard", "kept_frac", "added", "n_edge", "b")) Res[[cc]] <- as.numeric(Res[[cc]])
Res$group <- factor(Res$group, levels = GROUP_ORDER)
Res$dataset <- factor(Res$dataset, levels = dss)
time2 = Sys.time()
message(format(time2 - time1))


## ---------------------------- 4. Summary table --------------------------------------
metric_summary <- function(d) {
  qs <- quantile(d$jaccard, c(0.25, 0.5, 0.75), na.rm = TRUE)
  qk <- quantile(d$kept_frac, c(0.25, 0.5, 0.75), na.rm = TRUE)
  data.frame(n = nrow(d),
             jaccard_median = qs[2], jaccard_q1 = qs[1], jaccard_q3 = qs[3],
             kept_median = qk[2], kept_q1 = qk[1], kept_q3 = qk[3],
             added_median = median(d$added), n_edge_median = median(d$n_edge),
             row.names = NULL)
}
summary_bits <- lapply(
  split(Res, list(Res$dataset, Res$group, Res$type), drop = TRUE),
  function(d) data.frame(dataset = as.character(d$dataset[1]),
                         group = as.character(d$group[1]),
                         type = as.character(d$type[1]),
                         metric_summary(d))
)
Summary <- do.call(rbind, summary_bits)
write.csv(Summary, file.path(OUT_DIR, "E1_summary.csv"), row.names = FALSE)


## ---------------------------- 5. Supplementary figure ------------------------------
# Figure content:
# Publication-ready grouped permutation importance plot for the six evidence categories.
# Results are shown for all three datasets and two similarity metrics (Jaccard, kept fraction).
# The plot visualizes how strongly each evidence category contributes to the final peGRN
# at the terminal edge-selection stage of scDERI (no model retraining).

# How to read this plot:
# - Panel layout: faceted grid. Each row corresponds to one dataset; each column corresponds
# to one evaluation metric. Y-axis scales are free to adapt to each metric's range.
# - X-axis: six evidence groups (motif, prior, dynamics, access, TG_expr, TF_expr).
# - Faint scattered points: raw results from each of the B independent permutations.
# - Vertical pointrange: median value with interquartile range (IQR) across permutations,
# summarizing the central tendency and spread of the perturbation effect.
# - White diamond marker: deterministic neutral control. The adjustment factor of the
# corresponding group is fixed to 1 (factor range f_r = 0), so that feature group
# exerts no influence on the adjusted quantile threshold.
# - Y-axis: similarity score relative to the baseline peGRN. 1 means perfectly identical
# to the baseline result; 0 means no overlap at all.

# How to interpret the results:
# - A lower similarity score (closer to 0) means that shuffling that evidence category
# causes large changes to the final edge set → that category is more important / decisive.
# - A higher similarity score (closer to 1) means the edge set barely changes after
# permutation → that category contributes less to the final selection.
# - Comparing the permutation summary to the neutral control diamond shows whether the
# effect comes from destroying the feature information, or just from removing the
# feature's weighting factor.

# Metrics definition:
# - Jaccard: Jaccard similarity coefficient between the permuted edge set and the baseline.
# Jaccard = |E_perm ∩ E_base| / |E_perm ∪ E_base|, reflecting overall edge-set overlap.
# - kept fraction: proportion of baseline edges that remain after permutation.
# kept_frac = |E_perm ∩ E_base| / |E_base|, reflecting retention of original edges.

# Output: raster PNG and vector PDF of the figure written to OUT_DIR.

# Median with interquartile range for stat_summary (no Hmisc dependency).
median_iqr <- function(x) data.frame(y = median(x), ymin = quantile(x, 0.25),
                                     ymax = quantile(x, 0.75))
plot_df <- Res[Res$type %in% c("permute", "neutral"), ]
p1 <- plot_df; p1$metric <- "Jaccard";        p1$value <- p1$jaccard
p2 <- plot_df; p2$metric <- "kept fraction";  p2$value <- p2$kept_frac
plot_df <- rbind(p1, p2)
plot_df$metric <- factor(plot_df$metric, levels = c("Jaccard", "kept fraction"))

p <- ggplot(plot_df, aes(x = group, y = value, color = group)) +
  # Individual permutations (light jitter).
  geom_jitter(data = subset(plot_df, type == "permute"),
              width = 0.18, alpha = 0.25, size = 0.7, show.legend = FALSE) +
  # Median + IQR across the B permutations.
  stat_summary(data = subset(plot_df, type == "permute"),
               fun.data = median_iqr, geom = "pointrange",
               size = 0.5, show.legend = FALSE) +
  # Deterministic neutral controls (f_r = 0 -> F_m = 1).
  geom_point(data = subset(plot_df, type == "neutral"),
             shape = 23, size = 3, fill = "white", stroke = 0.8, show.legend = FALSE) +
  facet_grid(dataset ~ metric, scales = "free_y") +
  labs(x = NULL, y = "similarity to the baseline peGRN (1 = no change)",
       title = "Grouped permutation importance of the six evidence categories",
       caption = "Pointrange: median with IQR across permutations; point: raw replicate; diamond: neutral control (factor set to 1).") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        strip.text = element_text(face = "bold"),
        plot.caption = element_text(size = 7, hjust = 0))
print(p)
ggsave(file.path(OUT_DIR, "E1_permutation_importance.png"), p, width = 6, height = 4.5, dpi = 300)


## ---------------------------- 6. Console report -------------------------------------
message("\nRows per dataset:"); print(table(Res$dataset))
# A549 GM12878      BM 
# 8005    8005    8005 

message("\nRows per group:"); print(table(Res$group))
# motif    prior dynamics   access  TG_expr  TF_expr 
#  3000     3003     9003     3003     3003     3003 

message("\nRows per type:"); print(table(Res$type))
# neutral              permute    permute_loss_only permute_theta_p_only 
#      15                18000                 3000                 3000

message("\nMedian Jaccard / kept fraction (permute):")
print(Summary[Summary$type == "permute",
              c("dataset", "group", "jaccard_median", "jaccard_q1", "jaccard_q3",
                "kept_median", "kept_q1", "kept_q3")], row.names = FALSE)
# dataset    group jaccard_median jaccard_q1 jaccard_q3 kept_median   kept_q1   kept_q3
#    A549    motif      0.1786214  0.1741525  0.1831183   0.2702068 0.2645677 0.2758459
# GM12878    motif      0.1811308  0.1774987  0.1846343   0.2815780 0.2765110 0.2862830
#      BM    motif      0.1872470  0.1840413  0.1905992   0.2798429 0.2756120 0.2843759
#    A549    prior      0.5303929  0.5252097  0.5356380   0.6715226 0.6663534 0.6762218
# GM12878    prior      0.5342505  0.5296571  0.5386392   0.6706478 0.6663047 0.6749910
#      BM    prior      0.4897038  0.4850815  0.4939493   0.6400725 0.6358416 0.6446812
#    A549 dynamics      0.8570170  0.8531646  0.8610742   0.9149436 0.9125940 0.9182331
# GM12878 dynamics      0.8696970  0.8666385  0.8727688   0.9319580 0.9294245 0.9341296
#      BM dynamics      0.8665073  0.8636876  0.8696906   0.9307948 0.9286038 0.9335147
#    A549   access      0.5641620  0.5590523  0.5693299   0.7415414 0.7363722 0.7462406
# GM12878   access      0.5698369  0.5652854  0.5742777   0.7549765 0.7506334 0.7593196
#      BM   access      0.5550134  0.5504825  0.5586572   0.7588395 0.7546086 0.7624660
#    A549  TG_expr      0.6847755  0.6810383  0.6886886   0.7330827 0.7302632 0.7363722
# GM12878  TG_expr      0.6896201  0.6861815  0.6932041   0.7343467 0.7310894 0.7372421
#      BM  TG_expr      0.6704545  0.6672354  0.6735795   0.7141130 0.7113932 0.7166062
#    A549  TF_expr      0.6195792  0.6156488  0.6238147   0.6785714 0.6752820 0.6823308
# GM12878  TF_expr      0.6152174  0.6119803  0.6184255   0.6652190 0.6623236 0.6681144
#      BM  TF_expr      0.5818081  0.5786090  0.5844968   0.6343306 0.6316108 0.6367483

message("\nNeutral controls (factor set to 1):")
print(Summary[Summary$type == "neutral",
              c("dataset", "group", "jaccard_median", "kept_median",
                "added_median", "n_edge_median")], row.names = FALSE)
 # dataset    group jaccard_median kept_median added_median n_edge_median
 #    A549    prior      0.5432410   0.6080827          254          1548
 # GM12878    prior      0.5510271   0.6116540          304          1994
 #      BM    prior      0.5163117   0.5787247          400          2315
 #    A549 dynamics      0.9129239   0.9360902           54          2046
 # GM12878 dynamics      0.8973730   0.9272530           92          2654
 #      BM dynamics      0.8490840   0.8824418          130          3050
 #    A549   access      0.5635430   0.6188910          209          1526
 # GM12878   access      0.5617345   0.6141875          258          1955
 #      BM   access      0.5602876   0.6122696          307          2333
 #    A549  TG_expr      0.5637527   0.6212406          217          1539
 # GM12878  TG_expr      0.5642550   0.6181687          264          1972
 #      BM  TG_expr      0.5359298   0.5905107          337          2291
 #    A549  TF_expr      0.5262936   0.5878759          249          1500
 # GM12878  TF_expr      0.5229866   0.5805284          304          1908
 #      BM  TF_expr      0.4885599   0.5485041          406          2221

if (DECOMPOSE_DYNAMICS) {
  message("\nDynamics decomposition:")
  print(Summary[Summary$type %in% c("permute_theta_p_only", "permute_loss_only"),
                c("dataset", "group", "type", "jaccard_median",
                  "kept_median", "added_median")], row.names = FALSE)
  # dataset    group                 type jaccard_median kept_median added_median
  #    A549 dynamics    permute_loss_only      0.9078174   0.9445489           86
  # GM12878 dynamics    permute_loss_only      0.9038858   0.9489685          138
  #      BM dynamics    permute_loss_only      0.9039125   0.9483228          162
  #    A549 dynamics permute_theta_p_only      0.9018116   0.9473684          107
  # GM12878 dynamics permute_theta_p_only      0.9266597   0.9667029          119
  #      BM dynamics permute_theta_p_only      0.9441194   0.9752191          109
}

message("\nDone. Outputs written to ", OUT_DIR)

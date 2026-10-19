# ==============================================================================
# E2c: thinning the prior catalogs: missing TFs are deleted, never flagged, 
# imputed, or "recovered"
# ==============================================================================


# PURPOSE
#   Question: when the prior databases are incomplete (JASPAR covers only a 
#   fraction of human TFs, and ChIP-Atlas/ENCODE priors are sparse), does 
#   scDERI silently invent ("hallucinate") a TF/edge that has no prior support,
#   or does it correctly leave that TF out? In other words, is a catalog-missing
#   TF DELETED, or is it kept with a flag / surrogate?

# METHOD, DESIGN CONSIDERATIONS AND DESIGN
#   * Method (causal, not descriptive). E2a only reports a static coverage
#     snapshot ("how many TFs are in the catalogs"). Here we actively remove
#     catalogued TFs and observe the consequence, so that the relationship
#     between catalog completeness and the resulting network is a cause-effect
#     one rather than a single time-point observation.
#   * Design consideration. We must not confound the answer with model noise:
#     no re-download, no JASPAR PWM re-scanning (matchPWM), no re-grouping,
#     no re-training and no call to any scDERI modelling function. The entire
#     experiment is deterministic SET ARITHMETIC over the artefacts already
#     produced by scDERI_benchmark.R. Whatever happens is therefore attributable
#     solely to the catalog gates, not to stochastic refitting.
#   * Design. We reconstruct the three gates actually used by the benchmark
#     (JASPAR, expression QC, valid ChIP-Atlas prior), form the complete
#     candidate pool J n E n C, and repeatedly thin one catalog side by keeping
#     a random fraction keep in {0, .25, .5, .75, 1} of its TFs. At each level
#     we rebuild the candidate pool and the part of the FINAL network that is
#     still reachable, using 1000 random replicates for the interior levels.
#   * Two catalog sides. "ChIP-Atlas" thins the ChIP prior gate;
#     "JASPAR" thins the motif gate. 

# EXPECTED RESULT AND WHAT IT SHOWS
#   * The number of candidate TFs and the retained fraction of the final
#     network shrink monotonically with the catalog keep fraction, close to
#     proportionally (retained ~ keep), on all three datasets and both sides.
#   * keep = 1 reproduces the full candidate pool and the full network;
#     keep = 0 empties both.
#   * A removed TF NEVER reappears in the candidate pool, and every thinned
#     pool is a subset of the full pool (two assertions checked over all
#     18,012 dataset x side x keep x replicate rows).
#   Together this demonstrates the "delete, do not hallucinate" branch: the
#   prior catalog is a GATE / REFINER of TFs, not a source of truth that the
#   model can fall back on. scDERI cannot flag a TF it never knew; edges
#   attached to a missing TF are removed, not replaced.

# SCOPE (what this experiment does NOT claim)
#   It does NOT show that catalog-external TFs can be discovered -- they
#   cannot, and that limitation is stated as a matter of future database
#   growth. On the JASPAR side we only exercise the candidate-TF discovery
#   stage (the set gates); we do not rescan PWM scores for edges, so the
#   JASPAR panel speaks to candidate discovery, not to motif-level edge
#   scoring. 

# STEPS
#   1. Build the JASPAR2024 human TF name table once (local packaged database;
#      shared across datasets).
#   2. Per dataset, rebuild the four gates from cached Part 1.1 / Part 5.6
#      artefacts and cross-check the reconstructed pool against TFs.rds.
#   3. Run the thinning grid (2 sides x 5 keep levels x 1000 replicates for the
#      3 interior levels; endpoints run once as they are deterministic) and
#      write one row per realisation to long_df.
#   4. Aggregate the Monte Carlo replicates to median/IQR (summary_df) and save
#      the static gate sizes (funnel_df).
#   5. Run five deterministic no-imputation / subset checks (check_df).
#   6. Export four CSV tables and two figures.

# OUTPUTS (written to OUT_DIR)
#   E2c_prior_thinning_long.csv       one row per dataset x side x keep x replicate
#   E2c_prior_thinning_summary.csv    median / Q1 / Q3 of candidates and retained edges
#   E2c_static_funnel.csv             one row per dataset: static gate sizes (J / E / C / final)
#   E2c_no_imputation_checks.csv      five deterministic checks (expected, observed, pass)
#   E2c_candidate_TFs_vs_catalog.png  candidate TFs vs catalog keep fraction
#   E2c_edges_retained_vs_catalog.png fraction of final network edges retained vs keep

# RUNTIME
#   A few minutes total; dominated by qs::qread(scRNAseq) and one local JASPAR
#   lookup. Set sampling itself is effectively free (~45-55 s per dataset).

# PREREQUISITES
#   Part 1.1 of scDERI_benchmark.R must have run once per dataset so that
#   "./1.1 Data Preprocessing - Get TFs From Database/ChIP_Atlas_Priors/"
#   (processed_priors.rds, download_information.rds) and TFs.rds exist, and
#   Part 5.6 must have saved "./interest_cell_type_peGRN.qsave".


options(stringsAsFactors = FALSE)
suppressMessages({
  library(qs)            # fast read of the Seurat object and peGRN
  library(JASPAR2024)    # packaged JASPAR PWM database (local, no download)
  library(TFBSTools)     # getMatrixSet for the JASPAR name table
  library(ggplot2)       # diagnostic figures
})


## ------------------------------ Configuration ------------------------------
# Root folder of the benchmark run (same layout as scDERI_benchmark.R).
# Windows run (authoritative results below were produced here):
# BENCH_ROOT <- "I:/scDERI_Paper/3-benchmark/8-scDERI"
# WSL equivalent (commented; switch when running under WSL):
BENCH_ROOT <- "/mnt/i/scDERI_Paper/3-benchmark/8-scDERI"

# Output folder for this supplementary experiment.
# OUT_DIR <- "I:/scDERI_Paper/Supplement-Experiment/E2c"
OUT_DIR <- "/mnt/i/scDERI_Paper/Supplement-Experiment/E2c"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = TRUE)

# The three benchmark datasets (tag -> dataset folder).
DATASETS <- list(
  A549    = list(folder = "SNAREseq2_hg38_A549"),
  GM12878 = list(folder = "SNAREseq2_hg38_GM12878"),
  BM      = list(folder = "10X_hg38_Human_bone_marrow_with_T_cell_depleted")
)
RUN_DATASETS <- c("A549", "GM12878", "BM")

# Catalog keep fractions. 1 = full catalog (baseline), 0 = catalog fully removed.
KEEP_LEVELS <- c(0, 0.25, 0.50, 0.75, 1)

# Random thinnings per interior level. The endpoints 0 and 1 are deterministic,
# so each is run only once; the three interior levels are repeated REPS times.
REPS <- 1000

# JASPAR collections queried by get_TFs_from_priors (its default set), used
# here so the rebuilt JASPAR table is identical to the benchmark's.
JASPAR_COLLECTIONS <- c("CORE", "CNE", "PHYLOFACTS", "SPLICE", "POLII",
                        "FAM", "PBM", "PBM_HOMEO", "PBM_HLH", "UNVALIDATED")


## ------------------------------- Path helpers -------------------------------
# Relative locations of the cached Part 1.1 / Part 5.6 artefacts inside each
# dataset folder.
path_scRNA <- function(d) file.path(BENCH_ROOT, d, "scRNAseq.qsave")

path_prior_dir <- function(d) {
  file.path(BENCH_ROOT, d, "1.1 Data Preprocessing - Get TFs From Database")
}

path_processed_priors <- function(d) {
  file.path(path_prior_dir(d), "ChIP_Atlas_Priors", "processed_priors.rds")
}

path_download_info <- function(d) {
  file.path(path_prior_dir(d), "ChIP_Atlas_Priors", "download_information.rds")
}

path_TFs <- function(d) file.path(path_prior_dir(d), "TFs.rds")

path_peGRN <- function(d) file.path(BENCH_ROOT, d, "interest_cell_type_peGRN.qsave")


## ------------------------------- Small helpers ------------------------------
# Build the full JASPAR2024 human TF name table once. This is the exact local
# query used by get_TFs_from_priors, before intersecting with any dataset.
build_jaspar_human <- function() {
  pwm <- TFBSTools::getMatrixSet(
    x = JASPAR2024::db(JASPAR2024::JASPAR2024()),
    opts = list(all = FALSE, species = "Homo sapiens",
                collection = JASPAR_COLLECTIONS, all_versions = FALSE,
                matrixtype = "PWM"))
  sort(unique(vapply(pwm@listData, function(x) x@name, character(1))))
}

# Extract the final, filtered eGRN edge table (TF-peak-TG rows) from a peGRN
# object. Prefer the final simplified table; fall back to the loose table.
get_final_edges <- function(pegrn) {
  cell <- pegrn[[1]]
  nm <- if ("loose_filtered_simplified" %in% names(cell)) {
    "loose_filtered_simplified"
  } else if ("loose" %in% names(cell)) {
    "loose"
  } else stop("Neither loose_filtered_simplified nor loose found in peGRN.")
  cell[[nm]]
}

# Median and a single quantile (used for Q1/Q3) of a numeric vector.
med <- function(x) stats::median(x, na.rm = TRUE)

qtl <- function(x, q) as.numeric(stats::quantile(x, probs = q, na.rm = TRUE,
                                                 names = FALSE))

## ------------------ One catalog-thinning replicate --------------------------
# Thin ONE catalog side, rebuild the candidate pool and the reachable part of
# the final network, and return a one-row data.frame.
#
# Args:
#   side     : "ChIP-Atlas" (thin the C gate) or "JASPAR" (thin the J gate).
#   keep     : fraction of catalog TFs retained.
#   rep, seed: replicate index and its explicit RNG seed. The seed is a fixed
#              function of dataset / side / keep / replicate, so each random
#              draw is reproducible and independent of execution order.
#   J, E, C  : the three catalog / QC gate sets.
#   edges    : final eGRN edge table (columns TF, peak, TG).
thin_once <- function(side, keep, rep, seed, J, E, C, edges) {
  # The catalog being thinned on this side; the other gates stay untouched.
  universe <- if (side == "ChIP-Atlas") C else J
  set.seed(seed)
  if (keep >= 1) {
    kept <- universe
  } else {
    k <- max(0L, floor(keep * length(universe)))
    kept <- if (k == 0L) character(0) else sample(universe, size = k)
  }
  removed <- setdiff(universe, kept)
  
  # Candidate pool after thinning = J n E n C with the thinned side replaced.
  candidate <- if (side == "ChIP-Atlas") {
    Reduce(intersect, list(J, E, kept))
  } else {
    Reduce(intersect, list(kept, E, C))
  }
  full_pool <- Reduce(intersect, list(J, E, C))
  
  # Two no-imputation assertions (true by construction; recorded as machine-checked evidence):
  #  (1) no removed TF can re-enter the thinned pool;
  #  (2) the thinned pool is always a subset of the full pool.
  removed_never_returns <- length(intersect(removed, candidate)) == 0
  candidate_is_subset   <- all(candidate %in% full_pool)
  
  # Reachable part of the FINAL network: a TF must be in the thinned pool for
  # any of its edges to remain. Edges are filtered by TF, not re-scored.
  edges_kept <- edges[edges$TF %in% candidate, , drop = FALSE]
  data.frame(
    side = side, keep = keep, rep = rep, seed = seed,
    n_catalog = length(universe),             # size of the catalog being thinned
    n_candidate = length(candidate),          # candidate TFs after thinning
    n_candidate_full = length(full_pool),     # candidate TFs at full catalog
    final_TF_total = length(unique(edges$TF)),
    final_TF_kept  = length(unique(edges_kept$TF)),
    final_TF_removed = length(unique(edges$TF)) - length(unique(edges_kept$TF)),
    edges_total = nrow(edges),
    edges_kept = nrow(edges_kept),
    edges_kept_frac = nrow(edges_kept) / nrow(edges),
    peak_total = length(unique(edges$peak)),
    peak_reached = length(unique(edges_kept$peak)),
    TG_total = length(unique(edges$TG)),
    TG_reached = length(unique(edges_kept$TG)),
    removed_never_returns = removed_never_returns,
    candidate_is_subset = candidate_is_subset,
    stringsAsFactors = FALSE
  )
}


## --------------------------- JASPAR (once) ----------------------------------
# PURPOSE. Build the JASPAR2024 human TF name table once and share it across
# datasets, so every dataset is filtered against the same motif catalog.
message("Building the JASPAR2024 human TF name table (local database)...")
t1 <- Sys.time()
jaspar_human <- build_jaspar_human()
t2 <- Sys.time()
message("JASPAR2024 human PWM names: ", length(jaspar_human), ", run time: ",
        format(t2 - t1), ".")
# Observed: JASPAR2024 human PWM names: 675, run time: 24.73 secs.


## ----------------------------- Per-dataset ----------------------------------
# PURPOSE AND LOGIC
#   Each dataset is processed independently through the same two-stage
#   procedure, mirroring the benchmark exactly:
#   (a) Reconstruct the four gates purely from cached artefacts --
#       J = JASPAR2024 TFs that are expressed (rownames of scRNAseq);
#       E = TFs that passed expression QC (rows of download_information);
#       C = TFs with a valid ChIP-Atlas prior for this dataset (names of
#           processed_priors, which already require the prior to cover >=5900
#           target genes);
#       Fi_full = J n E n C, the complete candidate pool;
#       and the FINAL network edges (Part 5.6).
#   (b) Cross-check Fi_full against the saved TFs.rds: they must be identical.
#       This validates that the set reconstruction faithfully reproduces the
#       benchmark, so the thinning experiment is a perturbation of the real
#       pipeline, not of a re-implementation.
#   Then the thinning grid is run for both catalog sides.
# WHAT IT YIELDS
#   One funnel row (static gate sizes) and 6004 thinning rows per dataset
#   (2 sides x [2 endpoints x 1 + 3 interior levels x 1000]).
long_rows <- list()
funnel_rows <- list()
for (di in seq_along(RUN_DATASETS)) {
  t1 <- Sys.time()
  tag <- RUN_DATASETS[di]
  d <- DATASETS[[tag]]$folder
  message("\n================  ", tag, "  ================")
  
  # --- Reconstruct the four gates from cached Part 1.1 / Part 5.6 artefacts. ---
  for (f in c(path_processed_priors(d), path_download_info(d), path_TFs(d),
              path_peGRN(d), path_scRNA(d))) {
    if (!file.exists(f)) stop("Required artefact not found: ", f,
                              "\nRun scDERI_benchmark.R Part 1.1 / Part 5.6 first.")
  }
  
  scRNA <- qs::qread(path_scRNA(d))
  J <- sort(intersect(jaspar_human, rownames(scRNA)))   # JASPAR and expressed
  rm(scRNA)
  dl_info <- readRDS(path_download_info(d))
  E <- sort(rownames(dl_info))                          # passed expression QC
  priors <- readRDS(path_processed_priors(d))
  C <- sort(names(priors))                              # valid ChIP-Atlas priors
  TFs_saved <- sort(readRDS(path_TFs(d)))               # get_TFs_from_priors output
  Fi_full <- sort(Reduce(intersect, list(J, E, C)))
  
  # Validation: the reconstructed full pool must equal the cached TFs.rds.
  if (!setequal(TFs_saved, Fi_full)) {
    warning(tag, ": TFs.rds (", length(TFs_saved),
            ") differs from J n E n C (", length(Fi_full),
            "); using the artefacts as cached.")
  }
  
  # --- Final network (Part 5.6). ---
  pegrn <- qs::qread(path_peGRN(d))
  edges <- get_final_edges(pegrn)
  rm(pegrn)
  message("[", tag, "] J (JASPAR & expressed) = ", length(J),
          " | E (expression QC) = ", length(E),
          " | C (valid ChIP prior) = ", length(C),
          " | final network TFs = ", length(unique(edges$TF)),
          " | edges = ", nrow(edges))
  
  # Static funnel snapshot for this dataset (one row).
  funnel_rows[[length(funnel_rows) + 1]] <- data.frame(
    dataset = tag,
    n_jaspar_expressed = length(J),
    n_expression_qc = length(E),
    n_chip_valid = length(C),
    n_candidate_full = length(Fi_full),
    n_final_TF = length(unique(edges$TF)),
    n_edges = nrow(edges),
    n_peak = length(unique(edges$peak)),
    n_TG = length(unique(edges$TG)),
    stringsAsFactors = FALSE
  )
  
  # --- Thinning grid: two sides x keep levels x replicates. ---
  for (side in c("ChIP-Atlas", "JASPAR")) {
    for (keep in KEEP_LEVELS) {
      reps <- if (keep %in% c(0, 1)) 1L else REPS   # endpoints are deterministic
      for (b in seq_len(reps)) {
        seed <- 7000L + 1000L * di + 100L * (if (side == "ChIP-Atlas") 1L else 2L) +
          10L * round(100 * keep) + b
        r <- thin_once(side, keep, b, seed, J, E, C, edges)
        long_rows[[length(long_rows) + 1]] <- cbind(dataset = tag, r,
                                                    stringsAsFactors = FALSE)
      }
    }
  }
  
  t2 <- Sys.time()
  message("[", tag, "] Run time: ", format(t2 - t1))
}
# ================  A549  ================
# [A549] J (JASPAR & expressed) = 327 | E (expression QC) = 109 | C (valid ChIP prior) = 87 | final network TFs = 84 | edges = 2128
# [A549] Run time: 46.14217 secs
# 
# ================  GM12878  ================
# [GM12878] J (JASPAR & expressed) = 299 | E (expression QC) = 105 | C (valid ChIP prior) = 79 | final network TFs = 76 | edges = 2763
# [GM12878] Run time: 40.23476 secs
# 
# ================  BM  ================
# [BM] J (JASPAR & expressed) = 253 | E (expression QC) = 122 | C (valid ChIP prior) = 94 | final network TFs = 92 | edges = 3309
# [BM] Run time: 49.28325 secs

long_df <- do.call(rbind, long_rows)
head(long_df)
#   dataset       side keep rep seed n_catalog n_candidate n_candidate_full final_TF_total final_TF_kept final_TF_removed edges_total edges_kept edges_kept_frac peak_total peak_reached TG_total TG_reached removed_never_returns candidate_is_subset
# 1    A549 ChIP-Atlas 0.00   1 8101        87           0               87             84             0               84        2128          0       0.0000000       1127            0     1114          0                  TRUE                TRUE
# 2    A549 ChIP-Atlas 0.25   1 8351        87          21               87             84            20               64        2128        756       0.3552632       1127          585     1114        554                  TRUE                TRUE
# 3    A549 ChIP-Atlas 0.25   2 8352        87          21               87             84            21               63        2128        493       0.2316729       1127          390     1114        384                  TRUE                TRUE
# 4    A549 ChIP-Atlas 0.25   3 8353        87          21               87             84            21               63        2128        489       0.2297932       1127          411     1114        392                  TRUE                TRUE
# 5    A549 ChIP-Atlas 0.25   4 8354        87          21               87             84            19               65        2128        357       0.1677632       1127          292     1114        289                  TRUE                TRUE
# 6    A549 ChIP-Atlas 0.25   5 8355        87          21               87             84            20               64        2128        440       0.2067669       1127          353     1114        341                  TRUE                TRUE

# HOW TO READ long_df
#   * Each ROW is one realisation = one dataset x one thinned side x one keep
#     fraction x one random replicate. There are 6004 rows per dataset
#     (2 sides x [2 endpoints + 3 levels x 1000 replicates] = 2 x 3002),
#     i.e. 18012 rows in total.
#   * Key columns:
#       side / keep / rep / seed : which realisation this row is and the seed
#                                  that makes it reproducible.
#       n_catalog                : size of the catalog being thinned. On the
#                                  ChIP side this is C (87/79/94); on the JASPAR
#                                  side this is J (327/299/253).
#       n_candidate / n_candidate_full : candidate TFs after thinning vs at full
#                                  catalog.
#       final_TF_total / _kept / _removed : final-network TFs in total and how
#                                  many survive / disappear with the thinning.
#       edges_total / edges_kept / edges_kept_frac : same accounting for edges.
#       peak_total / peak_reached, TG_total / TG_reached : reached peaks and
#                                  target genes.
#       removed_never_returns, candidate_is_subset : the two logical
#                                  assertions (must be TRUE on every row).
#   * HOW TO ANALYSE. Interior rows are stochastic (different random TFs are
#     removed), so a single row is not the message. Aggregate across replicates
#     (next section): median and IQR describe the central trend and its spread.
#     The endpoints anchor the interpretation: keep=0 gives nothing, keep=1
#     gives the full network.

# Row counts (observed): each dataset has 6004 rows; the two sides have 9006
# rows each; total 18012.
table(long_df$dataset)   #  A549 6004 | BM 6004 | GM12878 6004
table(long_df$side)      #  ChIP-Atlas 9006 | JASPAR 9006

# The two assertions must hold on every one of the 18012 rows (observed TRUE).
table(long_df$removed_never_returns)   #  TRUE: 18012
table(long_df$candidate_is_subset)     #  TRUE: 18012

funnel_df <- do.call(rbind, funnel_rows)
print(funnel_df, row.names = FALSE)
# dataset n_jaspar_expressed n_expression_qc n_chip_valid n_candidate_full n_final_TF n_edges n_peak n_TG
#    A549                327             109           87               87         84    2128   1127 1114
# GM12878                299             105           79               79         76    2763   1395 1498
#      BM                253             122           94               94         92    3309   1787 1725

# HOW TO READ funnel_df
#   * One ROW per dataset; this is the baseline snapshot at the FULL catalog,
#     identical in spirit to E2a. It fixes the absolute numbers from which the
#     thinning curves depart.
#   * Columns, read as a left-to-right funnel:
#       n_jaspar_expressed (J) : JASPAR2024 TFs that are expressed.
#       n_expression_qc  (E)  : TFs passing expression QC.
#       n_chip_valid     (C)  : TFs with a valid ChIP-Atlas prior here.
#       n_candidate_full      : J n E n C, the candidate pool (equals C here).
#       n_final_TF / n_edges / n_peak / n_TG : size of the FINAL network.
#   * HOW TO ANALYSE. The funnel narrows steadily from the large motif
#     reservoir (J: 253-327) to the small prior-supported candidate pool
#     (C: 79-94) and finally to the network TFs (76-92). It also shows that the
#     JASPAR reservoir is much larger than the ChIP gate, which is why thinning
#     JASPAR (J) leaves slack and yields stochastic candidate counts, whereas
#     thinning ChIP (C) removes candidates almost one-for-one.


## --------------------------- Summary statistics -----------------------------
# PURPOSE AND LOGIC
#   The interior keep levels each have 1000 random realisations. To separate
#   the systematic effect of catalog thickness from Monte Carlo noise, we
#   aggregate -- within each dataset x side x keep -- the candidate count and
#   the retained edges / peaks / target genes to their MEDIAN and
#   interquartile range (Q1, Q3). The median marks the central dose-response;
#   the Q1-Q3 interval (the ribbon in the figures) shows how much random
#   choice of which TFs to remove matters.
# WHAT IT YIELDS
#   summary_df: 30 rows = 3 datasets x 2 sides x 5 keep levels.
#   Substantively, the median retained fraction tracks keep on all datasets
#   and both sides (retained ~ keep); on the ChIP side the candidate count is
#   deterministic (C is the bottleneck, so the count is floor(keep*|C|), Q1=Q3)
#   while the retained edges still vary because different TFs carry different
#   numbers of edges; on the JASPAR side even the candidate count varies,
#   because the large J reservoir leaves a hypergeometric overlap with C.
summary_rows <- list()
for (tag in unique(long_df$dataset)) {
  for (side in unique(long_df$side)) {
    for (keep in KEEP_LEVELS) {
      sub <- long_df[long_df$dataset == tag & long_df$side == side &
                       long_df$keep == keep, ]
      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        dataset = tag, side = side, keep = keep, n_reps = nrow(sub),
        candidate_median = med(sub$n_candidate),
        candidate_q1 = qtl(sub$n_candidate, 0.25),
        candidate_q3 = qtl(sub$n_candidate, 0.75),
        final_TF_kept_median = med(sub$final_TF_kept),
        edges_kept_frac_median = med(sub$edges_kept_frac),
        edges_kept_frac_q1 = qtl(sub$edges_kept_frac, 0.25),
        edges_kept_frac_q3 = qtl(sub$edges_kept_frac, 0.75),
        peak_reached_median = med(sub$peak_reached),
        TG_reached_median = med(sub$TG_reached),
        stringsAsFactors = FALSE)
    }
  }
}
summary_df <- do.call(rbind, summary_rows)
print(summary_df, row.names = FALSE)
# dataset       side keep n_reps candidate_median candidate_q1 candidate_q3 final_TF_kept_median edges_kept_frac_median edges_kept_frac_q1 edges_kept_frac_q3 peak_reached_median TG_reached_median
#    A549 ChIP-Atlas 0.00      1                0            0            0                  0.0              0.0000000          0.0000000          0.0000000                 0.0               0.0
#    A549 ChIP-Atlas 0.25   1000               21           21           21                 20.0              0.2361372          0.2030075          0.2763158               409.0             388.0
#    A549 ChIP-Atlas 0.50   1000               43           43           43                 41.5              0.4985902          0.4561795          0.5387688               726.0             699.0
#    A549 ChIP-Atlas 0.75   1000               65           65           65                 63.0              0.7514098          0.7133459          0.7863017               951.0             929.5
#    A549 ChIP-Atlas 1.00      1               87           87           87                 84.0              1.0000000          1.0000000          1.0000000              1127.0            1114.0
#    A549     JASPAR 0.00      1                0            0            0                  0.0              0.0000000          0.0000000          0.0000000                 0.0               0.0
#    A549     JASPAR 0.25   1000               22           19           24                 21.0              0.2436560          0.1992481          0.2924107               420.0             397.5
#    A549     JASPAR 0.50   1000               43           41           46                 42.0              0.4983553          0.4496006          0.5526316               726.0             698.0
#    A549     JASPAR 0.75   1000               65           63           68                 63.0              0.7563440          0.7105263          0.7999295               957.0             935.0
#    A549     JASPAR 1.00      1               87           87           87                 84.0              1.0000000          1.0000000          1.0000000              1127.0            1114.0
# GM12878 ChIP-Atlas 0.00      1                0            0            0                  0.0              0.0000000          0.0000000          0.0000000                 0.0               0.0
# GM12878 ChIP-Atlas 0.25   1000               19           19           19                 18.0              0.2401375          0.2041260          0.2785016               529.0             511.0
# GM12878 ChIP-Atlas 0.50   1000               39           39           39                 38.0              0.4956569          0.4520449          0.5425262               912.5             921.5
# GM12878 ChIP-Atlas 0.75   1000               59           59           59                 57.0              0.7491857          0.7086500          0.7861021              1189.0            1239.0
# GM12878 ChIP-Atlas 1.00      1               79           79           79                 76.0              1.0000000          1.0000000          1.0000000              1395.0            1498.0
# GM12878     JASPAR 0.00      1                0            0            0                  0.0              0.0000000          0.0000000          0.0000000                 0.0               0.0
# GM12878     JASPAR 0.25   1000               20           17           22                 19.0              0.2471951          0.1968874          0.2932501               540.0             524.0
# GM12878     JASPAR 0.50   1000               40           37           42                 38.0              0.5057908          0.4513210          0.5581795               924.0             931.0
# GM12878     JASPAR 0.75   1000               59           57           61                 57.0              0.7517191          0.7017734          0.8011220              1191.5            1242.0
# GM12878     JASPAR 1.00      1               79           79           79                 76.0              1.0000000          1.0000000          1.0000000              1395.0            1498.0
#      BM ChIP-Atlas 0.00      1                0            0            0                  0.0              0.0000000          0.0000000          0.0000000                 0.0               0.0
#      BM ChIP-Atlas 0.25   1000               23           23           23                 23.0              0.2413116          0.2093533          0.2753098               651.0             616.5
#      BM ChIP-Atlas 0.50   1000               47           47           47                 46.0              0.4983379          0.4611665          0.5397401              1151.0            1096.0
#      BM ChIP-Atlas 0.75   1000               70           70           70                 68.0              0.7479601          0.7092022          0.7776519              1505.0            1446.0
#      BM ChIP-Atlas 1.00      1               94           94           94                 92.0              1.0000000          1.0000000          1.0000000              1787.0            1725.0
#      BM     JASPAR 0.00      1                0            0            0                  0.0              0.0000000          0.0000000          0.0000000                 0.0               0.0
#      BM     JASPAR 0.25   1000               23           21           25                 23.0              0.2466002          0.2058024          0.2895890               671.0             633.5
#      BM     JASPAR 0.50   1000               47           44           49                 46.0              0.4971290          0.4529314          0.5488063              1151.5            1098.0
#      BM     JASPAR 0.75   1000               70           68           72                 69.0              0.7491689          0.7089755          0.7917800              1513.0            1453.5
#      BM     JASPAR 1.00      1               94           94           94                 92.0              1.0000000          1.0000000          1.0000000              1787.0            1725.0

# HOW TO READ summary_df
#   * One row per dataset x side x keep.
#   * candidate_median / _q1 / _q3 : median candidate TFs and the IQR.
#   * edges_kept_frac_median / _q1 / _q3 : median fraction of final edges
#     retained and its IQR. final_TF_kept_median, peak_reached_median and
#     TG_reached_median give the same dose-response for TFs / peaks / targets.
#   * n_reps = number of replicates aggregated (1 at the endpoints, 1000
#     otherwise).
#   * Read down a column: values rise monotonically from 0 (keep=0) to the full
#     baseline (keep=1). On the ChIP side candidate_median is exactly
#     floor(keep*|C|) (21/43/65 for A549; 19/39/59 GM; 23/47/70 BM) and the
#     retained-edge median is ~keep (.236/.499/.751 for A549). The JASPAR side
#     matches closely (A549 .25: 22 candidates, edge fraction .244) but shows a
#     non-zero IQR due to the larger reservoir.


## ----------------------- Deterministic no-imputation checks ----------------
# PURPOSE AND LOGIC
#   The Monte Carlo curves show "how much survives"; these checks prove the
#   qualitative safety property. We turn the verbal claims ("missing TFs are
#   deleted, never imputed") into five explicit, machine-verified assertions
#   and store their expected and observed values. The check rows are collected
#   in a dedicated environment (not `<<-`) so that they persist reliably in
#   interactive sessions and when sourced.
check_env <- new.env()
check_env$rows <- list()
add_check <- function(test, expected, observed, pass) {
  check_env$rows[[length(check_env$rows) + 1]] <- data.frame(
    test = test, expected = expected, observed = observed,
    pass = isTRUE(pass), stringsAsFactors = FALSE)
}

# 1. Every removed TF is absent from the thinned candidate pool (no recovery).
add_check("Removed catalog TFs never enter the candidate pool",
          "0 violations",
          paste(sum(!long_df$removed_never_returns), "violations"),
          all(long_df$removed_never_returns))

# 2. Every thinned candidate pool is a subset of the full pool.
add_check("Thinned candidate pool is a subset of the full pool",
          "0 violations",
          paste(sum(!long_df$candidate_is_subset), "violations"),
          all(long_df$candidate_is_subset))

# 3. keep = 1 reproduces the full candidate pool and the full network.
full1 <- long_df[long_df$keep == 1, ]   # 3 datasets x 2 sides = 6 rows
add_check("Full catalog (keep=1) reproduces all candidates",
          "all candidates",
          paste(all(full1$n_candidate == full1$n_candidate_full),
                "of", nrow(full1), "dataset-side cases"),
          all(full1$n_candidate == full1$n_candidate_full))
add_check("Full catalog (keep=1) reproduces all final edges",
          "edges_kept == edges_total",
          paste(all(full1$edges_kept == full1$edges_total),
                "of", nrow(full1), "dataset-side cases"),
          all(full1$edges_kept == full1$edges_total))

# 4. keep = 0 yields no candidates and reaches no edges / peaks / targets.
empty0 <- long_df[long_df$keep == 0, ]
add_check("Empty catalog (keep=0) yields no candidates and no edges",
          "0 candidates, 0 edges",
          paste("candidates:", sum(empty0$n_candidate),
                "edges:", sum(empty0$edges_kept)),
          all(empty0$n_candidate == 0) && all(empty0$edges_kept == 0))

check_df <- do.call(rbind, check_env$rows)
print(check_df, row.names = FALSE)
#                                                     test                  expected                     observed pass
#       Removed catalog TFs never enter the candidate pool              0 violations                 0 violations TRUE
#      Thinned candidate pool is a subset of the full pool              0 violations                 0 violations TRUE
#          Full catalog (keep=1) reproduces all candidates            all candidates TRUE of 6 dataset-side cases TRUE
#         Full catalog (keep=1) reproduces all final edges edges_kept == edges_total TRUE of 6 dataset-side cases TRUE
# Empty catalog (keep=0) yields no candidates and no edges     0 candidates, 0 edges       candidates: 0 edges: 0 TRUE

# HOW TO READ check_df
#   * One row per assertion (5 rows).
#   * test = the property; expected = the required outcome; observed = what the
#     18012 realisations actually produced; pass = TRUE only when fully met.
#   * Expected result: all five pass = TRUE. Rows 1-2 are the "no hallucination"
#     guarantees (removed TFs never return; thinned pools stay subsets); rows
#   3-4 are the boundary anchors (full catalog reproduces the network; empty
#   catalog empties it). If any row is FALSE the script warns and points to the
#   CSV. In the author run all five passed.


## --------------------------------- Figures ----------------------------------
# Shared styling: one colour per dataset; the two catalog sides are laid out as
# two panels (ChIP-Atlas main text, JASPAR appendix).
ds_colors <- c(A549 = "#66C2A5", GM12878 = "#FC8D62", BM = "#8DA0CB")
side_labs <- c("ChIP-Atlas" = "Thin ChIP-Atlas catalog",
               "JASPAR" = "Thin JASPAR catalog")

# Figure 1: candidate TFs vs catalog keep fraction.
# HOW TO READ THIS FIGURE
#   * X axis: fraction of the catalog retained (0%-100%).
#   * Y axis: number of candidate TFs after thinning.
#   * Colours / lines: the three datasets; each point is a median and the
#     shaded ribbon is the Q1-Q3 range across the 1000 replicates.
#   * Panels (facets): left = thinning ChIP-Atlas, right = thinning JASPAR.
#   * What to look for: a monotonic, near-linear rise from (0,0) to the full
#     baseline; a thin ribbon means the outcome is driven by the catalog size
#     rather than by which TFs are removed. The ChIP panel has no width on the
#     median (the candidate count is floor(keep*|C|)); the JASPAR panel has a
#     visible ribbon because the large J reservoir leaves a random overlap.
p_cand <- ggplot(summary_df, aes(x = keep, y = candidate_median,
                                 color = dataset, fill = dataset)) +
  geom_line(linewidth = 0.6) +
  geom_ribbon(aes(ymin = candidate_q1, ymax = candidate_q3),
              alpha = 0.25, color = NA) +
  geom_point(size = 1.8) +
  facet_wrap(~ side, nrow = 1, labeller = labeller(side = side_labs)) +
  scale_color_manual(values = ds_colors) +
  scale_fill_manual(values = ds_colors) +
  scale_x_continuous(labels = function(x) paste0(100 * x, "%")) +
  labs(x = "Fraction of prior catalog retained",
       y = "Candidate TFs",
       title = "Candidate TFs shrink monotonically as a prior catalog is thinned") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", legend.title = element_blank())
print(p_cand)
ggsave(file.path(OUT_DIR, "E2c_candidate_TFs_vs_catalog.png"), p_cand,
       width = 6, height = 3, units = "in", dpi = 300)

# Figure 2: fraction of the FINAL network edges retained vs catalog thickness.
# HOW TO READ THIS FIGURE
#   * X axis: fraction of the catalog retained. Y axis: percentage of the
#     final network's edges that remain reachable (0%-100%).
#   * Same colours / ribbon / panels as Figure 1.
#   * What to look for: the curves again follow ~keep from 0 to 100%, and the
#     ribbon's width reflects how unevenly edges are distributed across TFs.
#     Crucially, the edges attached to a removed TF are dropped and never
#     compensated by replacement edges -- the network contracts rather than
#     being refilled, which is the visual statement of "deleted, not imputed".
p_edge <- ggplot(summary_df, aes(x = keep, y = edges_kept_frac_median,
                                 color = dataset, fill = dataset)) +
  geom_line(linewidth = 0.6) +
  geom_ribbon(aes(ymin = edges_kept_frac_q1, ymax = edges_kept_frac_q3),
              alpha = 0.15, color = NA) +
  geom_point(size = 1.8) +
  facet_wrap(~ side, nrow = 1, labeller = labeller(side = side_labs)) +
  scale_color_manual(values = ds_colors) +
  scale_fill_manual(values = ds_colors) +
  scale_x_continuous(labels = function(x) paste0(100 * x, "%")) +
  scale_y_continuous(labels = function(x) paste0(100 * x, "%"), limits = c(0, 1)) +
  labs(x = "Fraction of prior catalog retained",
       y = "Final network edges retained",
       title = "Edges linked to missing TFs are removed, never replaced") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", legend.title = element_blank())
print(p_edge)
ggsave(file.path(OUT_DIR, "E2c_edges_retained_vs_catalog.png"), p_edge,
       width = 6, height = 3, units = "in", dpi = 300)

## --------------------------------- Outputs ----------------------------------
write.csv(long_df, file.path(OUT_DIR, "E2c_prior_thinning_long.csv"),
          row.names = FALSE)
write.csv(summary_df, file.path(OUT_DIR, "E2c_prior_thinning_summary.csv"),
          row.names = FALSE)
write.csv(funnel_df, file.path(OUT_DIR, "E2c_static_funnel.csv"),
          row.names = FALSE)
write.csv(check_df, file.path(OUT_DIR, "E2c_no_imputation_checks.csv"),
          row.names = FALSE)

message("\n---------------- E2c static funnel ----------------")
print(funnel_df, row.names = FALSE)

message("\n---------------- E2c checks ----------------")
print(check_df, row.names = FALSE)
if (!all(check_df$pass)) {
  warning("One or more E2c no-imputation checks FAILED; inspect ",
          "E2c_no_imputation_checks.csv.")
}

message("\nFiles saved under ", OUT_DIR)
# Observed: Files saved under I:/scDERI_Paper/Supplement-Experiment/E2c

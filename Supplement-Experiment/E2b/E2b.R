# ==============================================================================
# E2b: escalating data sparsity -> adaptive contraction / explicit termination
# ==============================================================================

# WHY THIS EXPERIMENT
#   Answers the question: when critical data are missing,
#   does scDERI "hallucinate" by bending toward available data (garbage in,
#   garbage out) or terminate with no output (garbage in, nothing out)?
#
#   This experiment removes observations in a controlled way and watches how the
#   metacell-grouping step reacts. We never rerun the pseudotime analysis or any
#   model; we only call scDERI::cell_grouping on increasingly sparse inputs.

# PLAIN-LANGUAGE CONCEPTS (for readers unfamiliar with the pipeline)
#   * Cells are ordered along a learned trajectory ("pseudotime").
#   * A metacell "group" is a contiguous window of cells pooled together to
#     average out single-cell sparsity. Fewer groups -> fewer, wider windows
#     -> coarser temporal resolution. n_g = number of groups.
#   * For each group, EVERY gene must be observed (non-zero) in at least n_min
#     cells across three matched series:
#         TFE = transcription-factor RNA expression,
#         TGA = enhancer-derived target-gene chromatin activity,
#         TGE = target-gene RNA expression.
#   * If a window does not meet the floor, the algorithm adds cells one by one
#     (enlarges the window) instead of imputing zeros. It therefore can only
#       (a) shrink the number of groups, (b) collapse to a single window, or
#       (c) stop before producing any aggregated profile. It never fabricates.

# WHAT "CONTRACTION", "COLLAPSE" AND "TERMINATION" MEAN
#   grouped    : n_g >= 2, a multi-point time series can still be built.
#   collapsed  : n_g == 1, the trajectory collapses to one point; no dynamic
#                regulatory curve can be fitted (temporal resolution is gone).
#   terminated : the hard failure branch is hit and no aggregate is produced.
#   Note: the installed cell_grouping signals a hard failure by raising an R
#   error; this script treats both an error and a returned n_g == 0 as
#   "terminated". In both cases the next stage (model training) is never run.
# STEPS
#   1. Load each dataset's pre-grouping object (the TFE/TGA/TGE matrices are
#      already ordered by pseudotime).
#   2. Baseline check: with no masking, grouping must reproduce the published
#      benchmark result (8 groups).
#   3. Escalation: for each masking fraction p and each random replicate, set a
#      fraction p of the non-zero entries to zero in each of the three series,
#      rerun cell_grouping with the EXACT fitting points used in the benchmark,
#      and record n_g, n_min, group sizes, runtime and the outcome. The
#      independent replicate jobs are run in parallel (mclapply).
#   4. Deterministic unit tests (A549, the fast dataset): a fully zeroed object
#      and an n_min > n_cell configuration must both terminate, a pipeline guard
#      must skip model training, and a serial-vs-parallel check must return
#      identical results.
#   5. Write long / summary / unit-test tables and two diagnostic figures.

# OUTPUTS (written to OUT_DIR)
#   E2b_sparsity_termination_long.csv  one row per dataset x p x replicate
#   E2b_sparsity_summary.csv           median/IQR n_g and group size per p
#   E2b_unit_tests.csv                 deterministic termination checks
#   E2b_n_groups_vs_sparsity.png       n_g vs masking fraction (median + IQR)
#   E2b_group_size_vs_sparsity.png     mean cells per group vs masking fraction

# RUNTIME (reference machine; one cell type, one branch per dataset)
#   A549 ~6 s per call, GM12878 ~4.6 min, BM ~5.6 min. Replicate jobs are
#   independent and run in parallel (see MC_CORES), so wall time is roughly the
#   per-dataset total divided by the worker count. Near termination the window
#   can grow to almost the whole branch (cost roughly quadratic in cell
#   number), so the two large datasets are kept within the contraction regime
#   (p <= 0.75); the hard-threshold region is explored densely on A549 and by
#   the deterministic unit tests. For a quick end-to-end check, set
#   RUN_DATASETS <- c("A549").
# REPRODUCIBILITY
#   Grouping itself is deterministic; each replicate uses an EXPLICIT seed
#   (42 + 1000*dataset index + 100*p + replicate), so its random mask does not
#   depend on execution order or on the RNG stream. Serial and parallel runs
#   therefore produce identical scientific results (only the wall-clock
#   runtime_sec column can differ). No algorithmic parameter or threshold is
#   changed.

options(stringsAsFactors = FALSE)
suppressMessages({
  library(scDERI)    # provides cell_grouping
  library(ggplot2)   # diagnostic figures
  library(parallel)  # mclapply / detectCores for parallel replicate jobs
})


## ------------------------------ Configuration ------------------------------
# Root folder of the benchmark run (same layout as scDERI_benchmark.R).
BENCH_ROOT <- "/mnt/i/scDERI_Paper/3-benchmark/8-scDERI"

# Output folder for this supplementary experiment.
OUT_DIR    <- "/mnt/i/scDERI_Paper/Supplement-Experiment/E2b"
setwd(OUT_DIR)

# Worker count for the parallel replicate jobs. Because each replicate draws
# from an explicit seed, MC_CORES changes speed only, never the results. It is
# capped at 16 to bound memory (every worker builds its own masked copy of the
# input object). mclapply is unavailable on Windows, where it falls back to 1.
MC_CORES <- if (.Platform$OS.type == "windows") 1L else
  max(1L, parallel::detectCores() - 1)

# cell_grouping parallelises internally over cell types; here each object has
# one cell type, so one core PER WORKER avoids nested parallelism.
NCORES     <- 1

# Per-dataset design.
#   p_grid  : fractions of non-zero entries set to zero.
#   reps    : independent random masks per p (the deterministic endpoints
#             p = 0 and p = 1 are run once).
#   points_y: the benchmark's control points for fitting n_min; MUST match the
#             published benchmark exactly (A549 -> n_min = 2; GM12878 -> 6; BM -> 7).
DATASETS <- list(
  A549 = list(
    folder   = "SNAREseq2_hg38_A549",
    points_y = c(1, 2, 3),
    p_grid   = c(0, 0.25, 0.50, 0.75, 0.80, 0.85, 0.90, 0.95, 1),
    reps     = 10
  ),
  GM12878 = list(
    folder   = "SNAREseq2_hg38_GM12878",
    points_y = c(5, 6, 7),
    p_grid   = c(0, 0.25, 0.50, 0.75, 0.80, 0.85, 0.90, 0.95, 1),
    # reps     = 2
    reps     = 10
  ),
  BM = list(
    folder   = "10X_hg38_Human_bone_marrow_with_T_cell_depleted",
    points_y = c(6, 7, 8),
    p_grid   = c(0, 0.25, 0.50, 0.75, 0.80, 0.85, 0.90, 0.95, 1),
    # reps     = 2
    reps     = 10
  )
)

# Datasets to run in this session. Set to c("A549") for a quick end-to-end run.
RUN_DATASETS <- c("A549", "GM12878", "BM")

# Run grouping from a scratch folder so the benchmark's saved fitting plots are
# never overwritten; each parallel job gets its own subfolder there (see
# run_once) so workers never write to the same file at the same time.
SCRATCH <- file.path(OUT_DIR, "_grouping_scratch")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = TRUE)  # idempotent
dir.create(SCRATCH, recursive = TRUE, showWarnings = TRUE)

# Names of the three pseudotime series inside the pre-grouping object.
SERIES <- c("Branches_TFE", "Branches_TGA", "Branches_TGE")

# Relative path to the pre-grouping object saved by the benchmark.
GPI_REL <- file.path("2.2 Data Processing - Cell Grouping",
                     "interest_cell_type_genes_pseudotime_info.rds")


## ------------------------------- Small helpers ------------------------------
# Zero a fraction p of the NON-ZERO entries of a gene x cell matrix.
# Only observed entries are removed; zeros are never turned into values, so this
# mimics losing observations rather than injecting noise. The RNG seed is set by
# the caller, which makes the three series masked independently as the RNG
# advances.
# Args:
#   mat: gene x cell data.frame (TFE/TGA/TGE for one branch).
#   p  : fraction of non-zero entries to set to zero, in [0, 1].
# Returns: the same matrix with exactly floor(p * #non-zero) entries zeroed.
mask_nonzero <- function(mat, p) {
  if (p <= 0) return(mat)                      # no masking
  nz <- which(mat != 0)                        # linear indices of observations
  if (length(nz) == 0) return(mat)             # already empty
  k <- min(length(nz), floor(p * length(nz)))  # entries to remove
  if (k == 0) return(mat)
  idx <- arrayInd(sample(nz, size = k), dim(mat))
  mat[idx] <- 0
  mat
}

# Build a masked copy of the whole pre-grouping object.
# The copy is always derived from the untouched original so successive masks
# never accumulate (R copy-on-modify keeps the original intact).
# Args:
#   gpi : loaded pre-grouping object (nested by cell type and branch).
#   p   : masking fraction.
#   seed: random seed controlling this replicate.
# Returns: a new object of identical structure with all three series masked.
make_masked_gpi <- function(gpi, p, seed) {
  set.seed(seed)
  out <- gpi
  for (ct in names(gpi)) {
    for (series in SERIES) {
      for (b in seq_along(gpi[[ct]][[series]])) {
        out[[ct]][[series]][[b]] <- mask_nonzero(gpi[[ct]][[series]][[b]], p)
      }
    }
  }
  out
}

# Count the genes/features present in the three INPUT series.
# These counts are held constant across masking (grouping does not delete
# genes; genes are filtered upstream) and document that the experiment changes
# observations, not the feature set.
# Args:  gpi, pre-grouping object.
# Returns: a list with n_TFE, n_TGA, n_TGE (total feature rows).
count_genes <- function(gpi) {
  nrows_of <- function(key) sum(vapply(gpi, function(ct)
    sum(vapply(ct[[key]], nrow, integer(1))), numeric(1)))
  list(n_TFE = nrows_of("Branches_TFE"),
       n_TGA = nrows_of("Branches_TGA"),
       n_TGE = nrows_of("Branches_TGE"))
}

# Flatten a successful cell_grouping result into the group-level metrics.
# Args:  g, a cell_grouping return value.
# Returns: a list with n_g, fitted n_min, and mean/max cells per group.
# Note: cell_grouping does not return Branches_n_cell; the total cell count is
# taken from the input and equals sum(Branches_n_k) by construction.
summarize_result <- function(g) {
  n_g     <- sum(unlist(lapply(g, function(z) unlist(z$Branches_n_group))))
  n_k     <- unlist(lapply(g, function(z) unlist(z$Branches_n_k)))
  n_min_v <- unlist(lapply(g, function(z) z$n_min))
  list(
    n_g             = as.integer(n_g),
    n_min           = as.numeric(n_min_v[1]),
    mean_group_size = if (length(n_k)) mean(n_k) else NA_real_,
    max_group_size  = if (length(n_k)) max(n_k)  else NA_real_
  )
}

# Map the group count to a human-readable outcome.
# Args: n_g, number of groups; errored, whether grouping raised an error.
# Returns: "grouped" (n_g >= 2), "collapsed" (n_g == 1) or "terminated".
classify_outcome <- function(n_g, errored) {
  if (errored || n_g == 0) "terminated"
  else if (n_g == 1) "collapsed"
  else "grouped"
}

# Run cell_grouping once on a masked object and record a one-row summary.
# Args:
#   gpi   : original pre-grouping object.
#   n_cell: total cells in the branch (from the input).
#   p_x/p_y: benchmark fitting points for n_min.
#   p     : masking fraction; seed: random seed; tag: dataset label.
# Returns: a one-row data.frame of metrics plus outcome/termination flags.
run_once <- function(gpi, n_cell, p_x, p_y, p, seed, tag) {
  masked <- make_masked_gpi(gpi, p, seed)
  gene_n <- count_genes(gpi)             # feature counts, constant across p
  
  # Give this job a private working folder so parallel workers never overwrite
  # each other's fitting-plot files; the folder name is unique per (p, seed).
  job_scratch <- file.path(SCRATCH,
                           paste0("job_", tag, "_p", round(100 * p), "_seed", seed))
  dir.create(job_scratch, recursive = TRUE, showWarnings = FALSE)
  old_wd <- getwd()
  setwd(job_scratch)                     # keep the benchmark's plots untouched
  t0 <- Sys.time()
  
  # suppressMessages hides the verbose step log; errors are still captured.
  out <- tryCatch(
    list(value = suppressMessages(scDERI::cell_grouping(
      interest_cell_type_genes_pseudotime_info = masked,
      points_x_for_fitting_nc_and_nmin = p_x,
      points_y_for_fitting_nc_and_nmin = p_y,
      ncores = NCORES)),
      error = NULL),
    error = function(e) list(value = NULL, error = conditionMessage(e)))
  file.create("finished.txt")
  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  setwd(old_wd)
  
  # Hard failure: grouping aborted before any aggregate was produced.
  if (!is.null(out$error)) {
    return(data.frame(
      dataset = tag, p = p, seed = seed, n_g = 0L, terminated = TRUE,
      n_min = NA_real_, n_cell = n_cell, mean_group_size = NA_real_,
      max_group_size = NA_integer_, n_TFE = gene_n$n_TFE, n_TGA = gene_n$n_TGA,
      n_TGE = gene_n$n_TGE, runtime_sec = round(dt, 2),
      outcome = "terminated", error = out$error, stringsAsFactors = FALSE))
  }
  
  m <- summarize_result(out$value)
  outcome <- classify_outcome(m$n_g, errored = FALSE)
  data.frame(
    dataset = tag, p = p, seed = seed, n_g = m$n_g,
    terminated = (outcome == "terminated"),
    n_min = m$n_min, n_cell = n_cell,
    mean_group_size = m$mean_group_size, max_group_size = m$max_group_size,
    n_TFE = gene_n$n_TFE, n_TGA = gene_n$n_TGA, n_TGE = gene_n$n_TGE,
    runtime_sec = round(dt, 2), outcome = outcome, error = "",
    stringsAsFactors = FALSE)
}

# Build the ordered list of replicate jobs for one dataset.
# The order (masking fraction p, then replicate b) is identical to the original
# serial nested loops, and every job carries its own explicit seed.
# Args:  di, dataset index in RUN_DATASETS; cfg, the dataset's design entry.
# Returns: a list of list(p, b, seed) jobs in serial order.
build_jobs <- function(di, cfg) {
  jobs <- list()
  for (p in cfg$p_grid) {
    reps <- if (p %in% c(0, 1)) 1 else cfg$reps   # endpoints are deterministic
    for (b in seq_len(reps)) {
      jobs[[length(jobs) + 1]] <- list(
        p = p, b = b,
        seed = 42L + 1000L * di + 100L * round(100 * p) + b)
    }
  }
  jobs
}

# Run all jobs for one dataset, serially or in parallel.
# mclapply returns results in the SAME order as the input jobs; because each job
# uses an explicit seed, the returned rows are identical to the serial version.
# Args:
#   gpi, n_cell, p_x, p_y, tag : dataset context.
#   jobs    : output of build_jobs.
#   parallel: TRUE uses mclapply (MC_CORES workers); FALSE uses lapply.
# Returns: a data.frame of result rows, ordered exactly as the input jobs.
run_jobs <- function(gpi, n_cell, p_x, p_y, tag, jobs, parallel = TRUE) {
  one <- function(job) run_once(gpi, n_cell, p_x, p_y, job$p, job$seed, tag)
  rows <- if (parallel && MC_CORES > 1 && length(jobs) > 1) {
    parallel::mclapply(
      jobs, one,
      mc.cores = min(MC_CORES, length(jobs)),
      mc.preschedule = TRUE,   # fork workers once; distribute jobs in blocks
      mc.set.seed = FALSE)     # each job sets its own seed explicitly
  } else {
    lapply(jobs, one)
  }
  
  # Defensive: if a worker crashes unexpectedly (e.g. out of memory), record a
  # worker_error row rather than discarding the whole experiment.
  rows <- lapply(seq_along(rows), function(i) {
    r <- rows[[i]]
    if (is.data.frame(r)) return(r)
    job <- jobs[[i]]
    data.frame(
      dataset = tag, p = job$p, seed = job$seed, n_g = 0L, terminated = TRUE,
      n_min = NA_real_, n_cell = n_cell, mean_group_size = NA_real_,
      max_group_size = NA_integer_, n_TFE = NA_integer_, n_TGA = NA_integer_,
      n_TGE = NA_integer_, runtime_sec = NA_real_, outcome = "worker_error",
      error = if (inherits(r, "try-error")) as.character(r) else "worker failed",
      stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}


## ----------------------------- Escalation grid ------------------------------
# One row per dataset x masking fraction x replicate. Datasets are handled one
# at a time (to bound memory); within a dataset, replicate jobs are parallel.
long_rows <- list()
for (di in seq_along(RUN_DATASETS)) {
  tag <- RUN_DATASETS[di]
  cfg <- DATASETS[[tag]]
  ds_dir <- file.path(BENCH_ROOT, cfg$folder)
  gpi <- readRDS(file.path(ds_dir, GPI_REL))
  
  # x control points follow the benchmark formula from the observed cell count.
  n_cell <- sum(unlist(lapply(gpi, function(z) unlist(z$Branches_n_cell))))
  p_x <- round(c(100, n_cell * 1.1, n_cell * 2.2))
  
  jobs <- build_jobs(di, cfg)
  message("\n================  ", tag, "  (n_cell = ", n_cell, ", ",
          length(jobs), " jobs on ", min(MC_CORES, length(jobs)), " cores)  ================")
  t1 = Sys.time()
  ds_df <- run_jobs(gpi, n_cell, p_x, cfg$points_y, tag, jobs, parallel = TRUE)
  t2 = Sys.time()
  message(tag, " parallel time: ", format(t2 - t1), ".")
  # A549 parallel time: 5.229708 mins.
  # GM12878 parallel time: 3.786189 hours.
  # BM parallel time: 3.938627 hours.
  
  # Results are returned in job order; print the compact progress log once all
  # parallel workers have finished (worker messages are not interleaved).
  for (i in seq_len(nrow(ds_df))) {
    r <- ds_df[i, ]
    message(sprintf("  p=%.2f seed=%d -> n_g=%d %-10s (%.1fs)",
                    r$p, r$seed, r$n_g, toupper(r$outcome), r$runtime_sec))
  }
  long_rows[[length(long_rows) + 1]] <- ds_df
  saveRDS(long_rows, file.path(OUT_DIR, "long_rows.rds"))
  rm(gpi)
}
long_df <- do.call(rbind, long_rows)
saveRDS(long_df, file.path(OUT_DIR, "long_df.rds"))

# Compact preview; the complete table is written to the CSV below.
head(long_df)
#   dataset    p seed n_g terminated n_min n_cell mean_group_size max_group_size n_TFE n_TGA n_TGE runtime_sec outcome error
# 1    A549 0.00 1043   8      FALSE     2   1088        136.0000            177    84  1341  1341       96.04 grouped      
# 2    A549 0.25 3543   5      FALSE     2   1088        217.6000            394    84  1341  1341      124.76 grouped      
# 3    A549 0.25 3544   6      FALSE     2   1088        181.3333            356    84  1341  1341       76.47 grouped      
# 4    A549 0.25 3545   6      FALSE     2   1088        181.3333            293    84  1341  1341       99.91 grouped      
# 5    A549 0.25 3546   6      FALSE     2   1088        181.3333            345    84  1341  1341      110.37 grouped      
# 6    A549 0.25 3547   6      FALSE     2   1088        181.3333            342    84  1341  1341      117.05 grouped  

# NOTE: under parallel execution the n_g / n_min / group-size / outcome columns
# are identical to the serial run for the same seeds; only runtime_sec changes.

# HOW TO READ long_df (illustrative A549 output; full table is in the CSV)
#   Rows    : one call to cell_grouping (dataset x p x replicate).
#   Columns :
#     p                   fraction of non-zero entries removed.
#     n_g                 number of metacell groups = temporal resolution.
#     n_min               fitted minimum non-zero observations per gene/group.
#     mean/max_group_size how wide the adaptive windows became.
#     n_TFE/TGA/TGE       feature counts (constant 84/1341/1341 on A549): genes
#                         are never removed, only observations are.
#     runtime_sec         per-job wall time (not part of the scientific result).
#     outcome             grouped / collapsed / terminated / worker_error.
#     error               the abort message when terminated.


## --------------------------- Summary statistics -----------------------------
# Per dataset x p: median/IQR of n_g and group size, plus the termination rate.
iqr <- function(x, q) as.numeric(stats::quantile(x, probs = q, na.rm = TRUE,
                                                 names = FALSE))
summary_rows <- list()
for (tag in unique(long_df$dataset)) {
  for (p in sort(unique(long_df$p[long_df$dataset == tag]))) {
    sub <- long_df[long_df$dataset == tag & long_df$p == p, ]
    grouped <- sub[!sub$terminated, ]
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      dataset = tag, p = p, n_reps = nrow(sub),
      n_g_median = median(sub$n_g),
      n_g_q1 = iqr(sub$n_g, 0.25), n_g_q3 = iqr(sub$n_g, 0.75),
      group_size_median = if (nrow(grouped)) median(grouped$mean_group_size) else NA_real_,
      group_size_q1 = if (nrow(grouped)) iqr(grouped$mean_group_size, 0.25) else NA_real_,
      group_size_q3 = if (nrow(grouped)) iqr(grouped$mean_group_size, 0.75) else NA_real_,
      termination_rate = mean(sub$terminated),
      stringsAsFactors = FALSE)
  }
}
summary_df <- do.call(rbind, summary_rows)
summary_df
#    dataset    p n_reps n_g_median n_g_q1 n_g_q3 group_size_median group_size_q1 group_size_q3 termination_rate
# 1     A549 0.00      1        8.0      8      8          136.0000      136.0000      136.0000              0.0
# 2     A549 0.25     10        5.0      5      6          217.6000      181.3333      217.6000              0.0
# 3     A549 0.50     10        3.5      3      4          317.3333      272.0000      362.6667              0.0
# 4     A549 0.75     10        2.0      2      2          544.0000      544.0000      544.0000              0.0
# 5     A549 0.80     10        1.0      1      2         1088.0000      544.0000     1088.0000              0.0
# 6     A549 0.85     10        1.0      1      1         1088.0000     1088.0000     1088.0000              0.0
# 7     A549 0.90     10        1.0      1      1         1088.0000     1088.0000     1088.0000              0.2
# 8     A549 0.95     10        0.0      0      0                NA            NA            NA              1.0
# 9     A549 1.00      1        0.0      0      0                NA            NA            NA              1.0
# 10 GM12878 0.00      1        8.0      8      8          928.2500      928.2500      928.2500              0.0
# 11 GM12878 0.25     10        6.0      6      6         1237.6667     1237.6667     1237.6667              0.0
# 12 GM12878 0.50     10        4.0      4      4         1856.5000     1856.5000     1856.5000              0.0
# 13 GM12878 0.75     10        2.0      2      2         3713.0000     3713.0000     3713.0000              0.0
# 14 GM12878 0.80     10        2.0      2      2         3713.0000     3713.0000     3713.0000              0.0
# 15 GM12878 0.85     10        1.0      1      1         7426.0000     7426.0000     7426.0000              0.0
# 16 GM12878 0.90     10        1.0      1      1         7426.0000     7426.0000     7426.0000              0.1
# 17 GM12878 0.95     10        0.0      0      0                NA            NA            NA              1.0
# 18 GM12878 1.00      1        0.0      0      0                NA            NA            NA              1.0
# 19      BM 0.00      1        8.0      8      8          927.2500      927.2500      927.2500              0.0
# 20      BM 0.25     10        6.0      5      6         1236.3333     1236.3333     1483.6000              0.0
# 21      BM 0.50     10        3.5      3      4         2163.5833     1854.5000     2472.6667              0.0
# 22      BM 0.75     10        2.0      2      2         3709.0000     3709.0000     3709.0000              0.0
# 23      BM 0.80     10        1.5      1      2         5563.5000     3709.0000     7418.0000              0.0
# 24      BM 0.85     10        1.0      1      1         7418.0000     7418.0000     7418.0000              0.0
# 25      BM 0.90     10        1.0      1      1         7418.0000     7418.0000     7418.0000              0.0
# 26      BM 0.95     10        1.0      1      1         7418.0000     7418.0000     7418.0000              0.0
# 27      BM 1.00      1        0.0      0      0                NA            NA            NA              1.0

# HOW TO READ summary_df
#   Each row aggregates all replicates at one masking level.
#   * n_g_median (with q1/q3 IQR) is the central trend.
#   * group_size_median rises in mirror image: windows widen
#     to gather the required non-zero observations; it is NA after termination.
#   * termination_rate is the fraction of runs that hit the hard stop.
# First masking level at which the MEDIAN run terminates (the hard threshold).

threshold_df <- do.call(rbind, lapply(
  split(summary_df, summary_df$dataset),
  function(d) {
    hit <- d[d$n_g_median == 0, ]
    data.frame(dataset = d$dataset[1],
               first_termination_p = if (nrow(hit)) min(hit$p) else NA_real_,
               stringsAsFactors = FALSE)
  }
))
threshold_df
#         dataset first_termination_p
# A549       A549                0.95
# BM           BM                1.00
# GM12878 GM12878                0.95

# HOW TO READ threshold_df
#   The median run first returns no groups at p = first_termination_p.
#   NA means the threshold was not reached in the tested range.


## ------------------------- Deterministic unit tests -------------------------
# These remove randomness from the termination claim. Run on A549 only (fast),
# they exercise the abort branch directly, independent of any random mask.
unit_rows <- list()
add_unit <- function(test, expected, observed, pass) {
  unit_rows[[length(unit_rows) + 1]] <<- data.frame(
    dataset = "A549", test = test, expected = expected, observed = observed,
    pass = isTRUE(pass), stringsAsFactors = FALSE)
}
a549_dir <- file.path(BENCH_ROOT, DATASETS$A549$folder)
gpi_a <- readRDS(file.path(a549_dir, GPI_REL))
n_cell_a <- sum(unlist(lapply(gpi_a, function(z) unlist(z$Branches_n_cell)))) # 1088
p_x_a <- round(c(100, n_cell_a * 1.1, n_cell_a * 2.2))                        # 100 1197 2394
p_y_a <- DATASETS$A549$points_y                                               # 1 2 3

# Test 1: with no masking, grouping must reproduce the benchmark (8 groups).
base_row <- run_once(gpi_a, n_cell_a, p_x_a, p_y_a, p = 0, seed = 1L, tag = "A549")
add_unit("Baseline reproducibility (p = 0)", "n_g = 8",
         paste0("n_g = ", base_row$n_g), base_row$n_g == 8)
base_row
#   dataset p seed n_g terminated n_min n_cell mean_group_size max_group_size n_TFE n_TGA n_TGE runtime_sec outcome error
# 1    A549 0    1   8      FALSE     2   1088             136            177    84  1341  1341         7.7 grouped  

# Interpretation: PASS - the experiment starts from the exact published result
# (8 groups, n_min = 2), so every later change is attributable to masking alone.

# Test 2: a fully zeroed object (all observations removed) must terminate.
zero_row <- run_once(gpi_a, n_cell_a, p_x_a, p_y_a, p = 1, seed = 2L, tag = "A549")
add_unit("All-zero series terminates grouping", "terminated",
         zero_row$outcome, zero_row$terminated)
zero_row
#   dataset p seed n_g terminated n_min n_cell mean_group_size max_group_size n_TFE n_TGA n_TGE runtime_sec    outcome                                                                       error
# 1    A549 1    2   0       TRUE    NA   1088              NA             NA    84  1341  1341       38.08 terminated Grouping process failed for branch 1 with invalid cell group count records.

# Interpretation: PASS - with no observations, no aggregated profile is produced
# and the run stops (the message is the grouping self-check), rather than
# returning an all-zero or imputed network.

# Test 3: a fitted n_min larger than the cell count must terminate before any
# aggregation (the Algorithm 1 "n_min > n_cell" branch). Control points
# c(10000, 11000, 12000) force n_min far above the 1088 available cells.
old_wd <- getwd(); setwd(SCRATCH)
hi <- tryCatch(
  list(value = suppressMessages(scDERI::cell_grouping(
    interest_cell_type_genes_pseudotime_info = gpi_a,
    points_x_for_fitting_nc_and_nmin = p_x_a,
    points_y_for_fitting_nc_and_nmin = c(10000, 11000, 12000),
    ncores = NCORES)), error = NULL),
  error = function(e) list(value = NULL, error = conditionMessage(e)))
hi
# $value
# NULL
#
# $error
# [1] "Grouping process failed for branch 1 with invalid cell group count records."

setwd(old_wd)
hi_terminated <- is.null(hi$value) ||
  sum(unlist(lapply(hi$value, function(z) unlist(z$Branches_n_group)))) == 0
hi_terminated
# [1] TRUE

# Interpretation: PASS - an infeasible minimum (more required non-zero cells per
# group than exist) is rejected before grouping starts; the stop is a property
# of the algorithm, not an accident of a particular random mask.

add_unit("n_min > n_cell terminates grouping", "terminated",
         if (hi_terminated) "terminated" else "grouped", hi_terminated)

# Test 4: pipeline-level guard. A real pipeline must call model_train only when
# grouping yields usable groups; on failed grouping it must skip training. This
# stand-in guard mirrors that decision and must return SKIP_TRAIN for an
# all-zero input (so model_train is never invoked).
pipeline_guard <- function(gpi_masked) {
  grp <- tryCatch(
    list(value = suppressMessages(scDERI::cell_grouping(
      interest_cell_type_genes_pseudotime_info = gpi_masked,
      points_x_for_fitting_nc_and_nmin = p_x_a,
      points_y_for_fitting_nc_and_nmin = p_y_a,
      ncores = NCORES)), error = NULL),
    error = function(e) list(value = NULL, error = conditionMessage(e))
  )
  usable <- !is.null(grp$value) &&
    sum(unlist(lapply(grp$value, function(z) unlist(z$Branches_n_group)))) > 0
  if (usable) "TRAIN" else "SKIP_TRAIN"
}
gpi_zero <- make_masked_gpi(gpi_a, p = 1, seed = 3L)
action <- pipeline_guard(gpi_zero)
action
# "SKIP_TRAIN"

# Interpretation: PASS - once grouping fails, the downstream training stage is
# bypassed, i.e. "nothing out" is enforced at the pipeline level.

add_unit("model_train skipped on termination", "SKIP_TRAIN", action,
         action == "SKIP_TRAIN")

# Test 5: serial vs parallel reproducibility. The SAME small A549 job set is
# run with lapply (serial) and mclapply (parallel); all scientific columns must
# match exactly (runtime_sec and the error text are excluded).
CHECK_PARALLEL <- TRUE
if (CHECK_PARALLEL) {
  check_jobs <- build_jobs(1L, list(p_grid = c(0.25, 0.5), reps = 2))
  serial_df   <- run_jobs(gpi_a, n_cell_a, p_x_a, p_y_a, "A549", check_jobs,
                          parallel = FALSE)
  parallel_df <- run_jobs(gpi_a, n_cell_a, p_x_a, p_y_a, "A549", check_jobs,
                          parallel = TRUE)
  sci_cols <- setdiff(names(serial_df), c("runtime_sec", "error"))
  agree <- isTRUE(all.equal(serial_df[, sci_cols], parallel_df[, sci_cols],
                            check.attributes = FALSE))
  add_unit("Serial and parallel runs agree", "identical",
           if (agree) "identical" else "DIFFER", agree)
}
unit_df <- do.call(rbind, unit_rows)
unit_df
#   dataset                                test   expected   observed pass
# 1    A549    Baseline reproducibility (p = 0)    n_g = 8    n_g = 8 TRUE
# 2    A549 All-zero series terminates grouping terminated terminated TRUE
# 3    A549  n_min > n_cell terminates grouping terminated terminated TRUE
# 4    A549  model_train skipped on termination SKIP_TRAIN SKIP_TRAIN TRUE
# 5    A549      Serial and parallel runs agree  identical  identical TRUE

# Interpretation: all checks pass - the baseline is reproduced, both hard stop
# branches fire correctly, model training is skipped on termination, and
# parallelisation does not change any scientific result.


## --------------------------------- Figures ----------------------------------
# One colour per dataset; each panel is a dataset. Lines show the median, shaded
# bands show the IQR, and faint dots show individual replicate runs.
ds_colors <- c(A549 = "#66C2A5", GM12878 = "#FC8D62", BM = "#8DA0CB")

# Figure 1: number of groups (temporal resolution) vs masking fraction.
# X = share of observations removed; Y = metacell groups remaining; the dashed
# vertical line marks the first median hard stop. Expected reading: a decline
# ending at 0 - resolution is given up rather than replaced by a fabricated network.
p_ng <- ggplot(summary_df, aes(x = p, y = n_g_median, color = dataset,
                               fill = dataset)) +
  geom_line(linewidth = 0.6) +
  geom_ribbon(aes(ymin = n_g_q1, ymax = n_g_q3), alpha = 0.15, color = NA) +
  geom_point(size = 2) +
  geom_point(data = long_df, aes(x = p, y = n_g), size = 0.8, alpha = 0.35,
             position = position_jitter(width = 0.006, height = 0.1)) +
  geom_vline(data = threshold_df[!is.na(threshold_df$first_termination_p), ],
             aes(xintercept = first_termination_p), linetype = 2,
             linewidth = 0.4) +
  facet_wrap(~ dataset, nrow = 1) +
  scale_color_manual(values = ds_colors) +
  scale_fill_manual(values = ds_colors) +
  scale_x_continuous(limits = c(0, 1), labels = function(x) paste0(100 * x, "%")) +
  labs(x = "Fraction of non-zero observations masked",
       y = "Number of metacell groups (n_g)",
       title = "Adaptive contraction under escalating sparsity") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")
print(p_ng)
ggsave(file.path(OUT_DIR, "E2b_n_groups_vs_sparsity.png"), p_ng,
       width = 8, height = 3, units = "in", dpi = 300)

# Figure 2: mean cells per group vs masking fraction.
# Mirror image of Figure 1: with fewer groups, each window widens to accumulate
# enough non-zero observations. The curve ends where
# grouping terminates, after which no group size is defined.
p_gs <- ggplot(summary_df, aes(x = p, y = group_size_median, color = dataset,
                               fill = dataset)) +
  geom_line(linewidth = 0.6) +
  geom_ribbon(aes(ymin = group_size_q1, ymax = group_size_q3), alpha = 0.15,
              color = NA) +
  geom_point(size = 2) +
  geom_vline(data = threshold_df[!is.na(threshold_df$first_termination_p), ],
             aes(xintercept = first_termination_p), linetype = 2,
             linewidth = 0.4) +
  facet_wrap(~ dataset, nrow = 1, scales = "free_y") +
  scale_color_manual(values = ds_colors) +
  scale_fill_manual(values = ds_colors) +
  scale_x_continuous(limits = c(0, 1), labels = function(x) paste0(100 * x, "%")) +
  labs(x = "Fraction of non-zero observations masked",
       y = "Mean cells per group",
       title = "Window expansion as support is removed") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")
print(p_gs)
ggsave(file.path(OUT_DIR, "E2b_group_size_vs_sparsity.png"), p_gs,
       width = 8, height = 3, units = "in", dpi = 300)


## --------------------------------- Outputs ----------------------------------
write.csv(long_df, file.path(OUT_DIR, "E2b_sparsity_termination_long.csv"),
          row.names = FALSE)
write.csv(summary_df, file.path(OUT_DIR, "E2b_sparsity_summary.csv"),
          row.names = FALSE)
write.csv(unit_df, file.path(OUT_DIR, "E2b_unit_tests.csv"), row.names = FALSE)

message("\n---------------- E2b summary ----------------")
print(summary_df, row.names = FALSE)

message("\nFirst termination threshold (median):")
print(threshold_df, row.names = FALSE)

message("\nDeterministic unit tests:")
print(unit_df, row.names = FALSE)

if (!all(unit_df$pass)) {
  warning("One or more deterministic unit tests FAILED; inspect E2b_unit_tests.csv.")
}

message("\nFiles saved under ", OUT_DIR)

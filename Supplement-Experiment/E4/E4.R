# ==============================================================================
# E4: GA random-seed stability of the trained regulatory strength
# ==============================================================================

# WHY THIS EXPERIMENT (one paragraph for the reviewer)
#   The genetic algorithm (GA) is the ONLY stochastic component of scDERI.
#   Everything else -- cell grouping, motif scoring, the gradient-ascent local
#   refinement and the final one-dimensional Brent search for the predicted
#   strength theta_p -- is deterministic. The supplement claims that repeated
#   runs with different GA seeds give highly consistent edge sets and
#   predicted strengths. This experiment produces the numbers behind that
#   claim, so reviewers can judge whether any reported edge depends on luck.

# PLAIN-LANGUAGE CONCEPTS (no scDERI background required)
#   * A "training pair" is one TF -> target-gene candidate whose kinetic model
#     is fitted independently.
#   * Fitting one pair = GA global search + gradient-ascent local refinement.
#     It outputs two numbers:
#         theta_p  -- predicted regulatory strength, in [0, 1] (large = strong);
#         fit_loss -- pure prediction error (small = good fit).
#   * The pair is finally called a regulatory edge only when BOTH hold:
#         theta_p >= tau  (the cell-type binarization threshold), and
#         theta_p >= 0.05 (the poor-fit removal used in benchmark Part 4.2).
#   * We retrain the SAME pair under five different GA seeds and ask:
#     (1) does theta_p stay numerically close (coefficient of variation)?
#     (2) does the binary edge call stay the same (agreement)?
#     (3) does the seed-123 run reproduce the published run exactly?

# WHAT IS HELD FIXED (only the seed changes)
#   Training pairs, inputs, time windows, tau and every training hyperparameter
#   are identical to benchmark Part 4.2.2. Seed 123 is the benchmark default;
#   it is the reproduction anchor. The other four seeds are arbitrary but
#   fixed, so every result is reproducible.

# PAIR SELECTION (no cherry-picking)
#   Pairs are ranked by their PRIOR strength theta_i, which is known before
#   training. Ten pairs are then sampled at even positions from the strongest
#   to the weakest prior, so the grid covers both easy and difficult fits.

# PIPELINE MAP (each arrow corresponds to one function below)
#     load_env ------------> read the six cached objects for one dataset
#     even_pairs -----------> pick pairs evenly over the theta_i range
#     subset_to_pair -------> slice the cached objects down to one pair
#     fit_one --------------> model_train(that pair, one seed) -> one result row
#     fit_job_safe ---------> wrapper that turns a failed fit into an error row
#     [run_dataset repeats fit_one over the pair x seed grid, then builds]
#        runs     : raw per-pair x seed results
#        repro    : seed-123 vs published results
#        summary  : per-pair CV / range / agreement
#        ds       : one dataset-level summary line

# STEPS
#   1. Load the cached training set, initial parameters, cell grouping, tau and
#      pseudotime information for each of the three datasets (A549, GM12878,
#      BM). Each dataset has one pseudotime branch.
#   2. Select ten evenly spaced training pairs.
#   3. Fit ONE probe pair to measure the per-fit wall time, then size the
#      pair x seed grid inside the time budget (see the "TIME BUDGET" block).
#   4. Retrain every selected pair under every seed (one model_train call per
#      job; independent jobs run in parallel via mclapply).
#   5. Check that seed 123 reproduces the published model_train_all results.
#   6. Compute per-pair and per-dataset stability summaries and save tables
#      and two figures.

# OUTPUTS (written to OUT_DIR)
#   E4_seed_runs.csv          one row per dataset x pair x seed (raw results)
#   E4_seed_summary.csv       one row per dataset x pair (CV, range, agreement)
#   E4_dataset_summary.csv    one row per dataset (headline stability numbers)
#   E4_seed_reproduction.csv  seed-123 versus the published Part 4.2 results
#   E4_theta_p_by_pair.png    theta_p per pair across seeds, with the tau line
#   E4_cv_vs_prior.png        per-pair CV of theta_p vs prior strength theta_i
#   E4_runs_<dataset>.rds     per-dataset checkpoint (safe against interruption)

# HOW TO RUN
#   * First pass (smoke test): keep SMOKE_TEST <- TRUE; the whole script
#     finishes in about one minute per dataset and only checks that the code
#     runs end to end. Numbers from this pass are NOT reportable.
#   * Real experiment: set SMOKE_TEST <- FALSE and rerun; each dataset then
#     takes roughly 10-20 minutes on a 50+ core machine.

# RUNTIME / PLATFORM
#   One pair is an almost single-core GA job (~7-16 min in the production run).
#   The pair x seed jobs are independent and run in parallel, so the wall time
#   per dataset is roughly one job's time. Windows has no mclapply fork; the
#   script then runs serially (set RUN_DATASETS <- c("A549") for a quick check).

# SCOPE / LIMITS
#   This is a targeted stability test of the random optimizer, NOT a
#   hyperparameter grid: population size, number of generations and mutation
#   settings are fixed to the published values (deferred to future work). Only
#   the selected training pairs are retrained; candidate construction,
#   inference and the terminal filter are untouched.


options(stringsAsFactors = FALSE)
suppressMessages({
  library(scDERI)     # provides model_train, the only algorithmic function called
  library(ggplot2)    # two diagnostic figures
  library(parallel)   # mclapply for parallel pair x seed jobs
})


## ============================= Configuration ===============================
# Folder of the benchmark run (same layout as scDERI_benchmark.R).
BENCH_ROOT <- "/mnt/i/scDERI_Paper/3-benchmark/8-scDERI"

# Folder that receives all E4 tables and figures.
OUT_DIR    <- "/mnt/i/scDERI_Paper/Supplement-Experiment/E4"

# Each job runs in a private scratch folder so that the benchmark's saved
# objects are never overwritten.
SCRATCH    <- file.path(OUT_DIR, "_train_scratch")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = TRUE)
dir.create(SCRATCH, recursive = TRUE, showWarnings = TRUE)

# The three matched datasets and their benchmark folders.
DATASETS <- list(
  A549    = list(folder = "SNAREseq2_hg38_A549"),
  GM12878 = list(folder = "SNAREseq2_hg38_GM12878"),
  BM      = list(folder = "10X_hg38_Human_bone_marrow_with_T_cell_depleted")
)

# Restrict to one dataset, e.g. c("A549"), for a quick check; run all three
# for the real experiment.
RUN_DATASETS <- c("A549", "GM12878", "BM")

# Smoke-test switch.
#   TRUE  -> 4 epochs / early stop 2: verifies the code in about one minute,
#            but the numbers are meaningless and the seed-123 reproduction is
#            EXPECTED to fail (the published run used 100 epochs).
#   FALSE -> 100 epochs / early stop 5: exactly the published Part 4.2.2
#            settings, used for the reportable results.
# SMOKE_TEST <- TRUE
SMOKE_TEST <- FALSE
TRAIN_MAX_EPOCHS <- if (SMOKE_TEST) 4 else 100
TRAIN_EARLY_STOP <- if (SMOKE_TEST) 2 else 5

# Size of the pair x seed grid.
#   N_PAIRS_FULL : number of pairs in the full design (10 = evenly spaced over
#                  the prior range; increase it to cover more pairs, e.g. 20;
#                  set it to a value >= the largest training set (85/76/92) to
#                  train ALL pairs -- even_pairs clips it to the available size).
#   N_SEEDS_FULL : number of seeds in the full design (must not exceed
#                  length(SEEDS), so extend SEEDS first if you want more).
#   MIN_PAIRS / MIN_SEEDS : smallest grid accepted by the time-budget
#                  downgrade (see the "TIME BUDGET" block below).
N_PAIRS_FULL <- 10
N_SEEDS_FULL <- 5
MIN_PAIRS    <- 5
MIN_SEEDS    <- 3

# Explicit GA seeds; the first one (123) is the benchmark default and anchors
# the reproduction check.
SEEDS <- c(123L, 42L, 2024L, 777L, 31337L)

# Optional fixed grid. Leave NULL to let the probe choose the grid
# automatically (recommended). Set it to force a specific grid regardless of
# the time budget, e.g. list(n_pairs = 10, n_seeds = 5).
FORCE_GRID <- NULL

# TIME BUDGET (why the grid can be "downgraded")
#   Step 3 times one fit, then estimates the wall time of a grid as
#       ceil(number_of_jobs / cores) x per-fit-time x safety factor.
#   The candidate grids are tried in order 10x5, 10x3, 5x3, and the LARGEST
#   one whose estimated time fits the remaining budget is used. If even 5x3
#   does not fit, it is still run and flagged as budget_exceeded. The chosen
#   grid is recorded in the output column grid_status, so a downgrade is
#   never hidden. Set TIME_BUDGET_HOURS high (or set FORCE_GRID) on a free
#   machine.
# TIME_BUDGET_HOURS <- 2
TIME_BUDGET_HOURS <- 20
# TIME_SAFETY       <- 1.2   # 20% margin on the probe-based estimate
TIME_SAFETY       <- 12   # 20% margin on the probe-based estimate

# Number of parallel jobs. The default uses every available core up to the
# number of jobs; leave it as is. Set it by hand (e.g. 78L) to force a
# specific core count; set it to 1L for serial execution.
MC_CORES <- if (.Platform$OS.type == "windows") 1L else
  min(parallel::detectCores() - 1, N_PAIRS_FULL * N_SEEDS_FULL)

# Training settings. max_epochs and early_stop follow the SMOKE_TEST switch;
# every other value is identical to benchmark Part 4.2.2 and must not change.
TRAIN_ARGS <- list(
  max_epochs      = TRAIN_MAX_EPOCHS,
  early_stop      = TRAIN_EARLY_STOP,
  eps_theta       = 1e-2,
  eps_loss        = 1e-5,
  popSize_ga      = 50,
  maxiter_ga      = 150,
  pcrossover_ga   = 0.9,
  pmutation_ga    = 0.5,
  parallel_ga     = FALSE,
  iterations_grad = 30,
  ncores_grad     = 1,
  ncores          = 1   # one job per fit; parallelism is handled at the job level
)

# Poor-fit removal used in benchmark Part 4.2: pairs with theta_p < 0.05 are
# dropped from the released network.
QC_THETA_P_MIN <- 0.05

# Relative paths of the six cached benchmark artifacts (folder names are kept
# exactly as in the benchmark, including the original "BUild" spelling).
RELS <- list(
  training_set = "4.1 Build Prediction Model - Select Training Set/interest_cell_type_branch_training_set.rds",
  init_params  = "4.2 BUild Prediction Model/interest_cell_type_branch_init_params.rds",
  group        = "2.2 Data Processing - Cell Grouping/interest_cell_type_group.rds",
  tau          = "3 get iGRN/interest_cell_type_tau.rds",
  gpi          = "2.2 Data Processing - Cell Grouping/interest_cell_type_genes_pseudotime_info.rds",
  ref_all      = "4.2 BUild Prediction Model/interest_cell_type_branch_model_train_all.rds"
)


## ===================== Small statistical helpers (Step 6) ===================
# These three functions only summarize numbers; they never touch the model.

# Coefficient of variation of a numeric vector, sd / mean.
# Used in Step 6 to summarize theta_p and fit_loss across seeds.
# Input: one numeric vector (the five values of one pair across seeds).
# Output: one number; NA when the mean is effectively zero, because CV is then
# not meaningful (the absolute range is reported instead).
safe_cv <- function(x) {
  m <- mean(x, na.rm = TRUE)
  if (is.na(m) || abs(m) < 1e-8) return(NA_real_)
  stats::sd(x, na.rm = TRUE) / m
}

# Fraction of binary calls equal to the majority call.
# Used in Step 6 for the tau threshold call and the final edge call.
# Input: one vector of 0/1 calls (one per seed).
# Output: 1.0 = all seeds agree; 0.6 = 3 of 5 agree; 0.5 = complete disagreement.
call_agreement <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_)
  max(mean(x == 1), mean(x == 0))
}

# Jaccard similarity of two edge sets, |intersection| / |union|.
# Used in Step 6 to compare edge sets produced by different seeds.
# Input: two character vectors of pair IDs.
# Output: 1.0 = identical sets; 0 = no shared edge. Two empty sets count as 1.
jaccard <- function(a, b) {
  u <- union(a, b)
  if (length(u) == 0) 1 else length(intersect(a, b)) / length(u)
}


## ======================= Step 1: data loading ==============================
# Read the six cached objects for one dataset and return them in one list.
# Input: cfg = list(folder = <benchmark subfolder of the dataset>).
# Output: a list with training_set, init_params, group, tau, gpi and ref_all.
# It only reads files; nothing is recomputed.
load_env <- function(cfg) {
  ds_dir <- file.path(BENCH_ROOT, cfg$folder)
  rd <- function(key) readRDS(file.path(ds_dir, RELS[[key]]))
  list(
    training_set = rd("training_set"), # pairs to train + TFE/TGA/TGE features
    init_params  = rd("init_params"),  # starting values / bounds per pair
    group        = rd("group"),        # pseudotime interval lengths
    tau          = rd("tau"),          # cell-type binarization threshold
    gpi          = rd("gpi"),          # cell counts per branch
    ref_all      = rd("ref_all")       # published training results (anchor)
  )
}


# ======================== Step 2: pair selection ============================
# Pick n pairs at even positions over the prior-strength range.
# Input: ts = cached training set; ct = dataset name; n = number of pairs.
# Output: a character vector of pair IDs, strongest prior first.
# If n exceeds the number of available pairs, all pairs are returned.
even_pairs <- function(ts, ct, n) {
  all_df <- ts[[ct]][["all"]]
  n <- min(n, nrow(all_df))
  ord <- order(all_df$theta_i, decreasing = TRUE)
  pick <- unique(round(seq(1, length(ord), length.out = n)))
  rownames(all_df)[ord[pick]]
}

# Slice the cached training set and initial parameters down to ONE pair, while
# preserving the branch1 / "all" layout that model_train expects.
# Input: env (output of load_env); ct; pair_id.
# Output: a two-element list with the subsetted training_set and init_params,
# wrapped so that model_train can be called exactly as in the benchmark.
subset_to_pair <- function(env, ct, pair_id) {
  n_branch <- length(env$training_set[[ct]]) - 1   # one element is "all"
  branches <- paste0("branch", seq_len(n_branch))
  one_ts <- list()
  for (b in branches) {
    one_ts[[b]] <- env$training_set[[ct]][[b]][pair_id, , drop = FALSE]
  }
  one_ts[["all"]] <- env$training_set[[ct]][["all"]][pair_id, , drop = FALSE]
  one_ip <- list()
  for (b in branches) {
    one_ip[[b]] <- env$init_params[[ct]][[b]][pair_id]
  }
  list(
    training_set = setNames(list(one_ts), ct),
    init_params  = setNames(list(one_ip), ct)
  )
}


## ===================== Step 4 unit: fit one pair x seed ====================
# Retrain one pair under one GA seed.
# Input: ct (dataset); pair_id; seed (integer); env (output of load_env).
# Output: a one-row data.frame with theta_p, fit_loss, the binary calls,
# training diagnostics and the error message (empty when successful).
# A failure is caught and returned as the same row with NA results, so one bad
# job never stops the whole grid.
fit_one <- function(ct, pair_id, seed, env) {
  sub <- subset_to_pair(env, ct, pair_id)
  meta <- env$training_set[[ct]][["all"]][pair_id, c("TF", "TG", "theta_i")]
  branches <- paste0("branch", seq_len(length(env$training_set[[ct]]) - 1))
  
  # Work in a private folder; model_train creates "4.2 Build Prediction Model"
  # inside it. on.exit restores the working directory even on error.
  job_dir <- file.path(SCRATCH, paste0("job_", ct, "_", make.names(pair_id),
                                       "_seed", seed))
  dir.create(job_dir, recursive = TRUE, showWarnings = FALSE)
  old_wd <- getwd()
  setwd(job_dir)
  on.exit(setwd(old_wd))
  
  t0 <- Sys.time()
  out <- tryCatch(
    suppressMessages(do.call(
      model_train,
      c(list(
        interest_cell_type_branch_training_set      = sub$training_set,
        interest_cell_type_branch_init_params       = sub$init_params,
        interest_cell_type_group                    = env$group,
        interest_cell_type_tau                      = env$tau,
        interest_cell_type_genes_pseudotime_info    = env$gpi,
        seed_ga                                     = as.integer(seed)
      ), TRAIN_ARGS))),
    error = function(e) e)
  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  
  # Row returned when the fit itself raised an error.
  if (inherits(out, "error")) {
    return(data.frame(
      dataset = ct, pair_id = pair_id, TF = meta$TF, TG = meta$TG,
      theta_i = meta$theta_i, seed = seed, theta_p = NA_real_,
      fit_loss = NA_real_, bin_tau = NA_integer_, qc_pass = NA,
      significant = NA_integer_, best_epoch = NA_integer_, stop_cause = "ERROR",
      runtime_sec = round(dt, 2), error = conditionMessage(out),
      stringsAsFactors = FALSE
    ))
  }
  
  # Extract the integrated one-branch result and the detailed training record.
  r <- out[[ct]]$best_pred_result[pair_id, , drop = FALSE]
  det <- out[[ct]][[branches[1]]][[pair_id]]
  theta_p <- as.numeric(r$theta_p)
  fit_loss <- as.numeric(r$fit_loss)
  bin_tau <- as.integer(r$theta_p_bin)          # 1 if theta_p >= tau
  qc_pass <- theta_p >= QC_THETA_P_MIN          # TRUE if theta_p >= 0.05
  data.frame(
    dataset = ct, pair_id = pair_id, TF = meta$TF, TG = meta$TG,
    theta_i = meta$theta_i, seed = seed, theta_p = theta_p,
    fit_loss = fit_loss, bin_tau = bin_tau, qc_pass = qc_pass,
    significant = as.integer(bin_tau == 1 & qc_pass), # final edge call
    best_epoch = as.integer(det$best_epoch),
    stop_cause = if (is.null(det$stop_cause)) NA_character_ else det$stop_cause,
    runtime_sec = round(dt, 2), error = "",
    stringsAsFactors = FALSE
  )
}

# Defensive wrapper around fit_one for parallel workers.
# Input/output: same as fit_one, but if the worker process itself dies, the
# condition is converted into a WORKER_ERROR row instead of being dropped.
fit_job_safe <- function(ct, pair_id, seed, env) {
  r <- tryCatch(fit_one(ct, pair_id, seed, env), error = function(e) e)
  if (inherits(r, "error")) {
    meta <- env$training_set[[ct]][["all"]][pair_id, c("TF", "TG", "theta_i")]
    return(data.frame(
      dataset = ct, pair_id = pair_id, TF = meta$TF, TG = meta$TG,
      theta_i = meta$theta_i, seed = seed, theta_p = NA_real_,
      fit_loss = NA_real_, bin_tau = NA_integer_, qc_pass = NA,
      significant = NA_integer_, best_epoch = NA_integer_,
      stop_cause = "WORKER_ERROR", runtime_sec = NA_real_,
      error = conditionMessage(r), stringsAsFactors = FALSE
    ))
  }
  r
}


## ===================== Step 3: automatic grid sizing =======================
# Estimate the wall time of a grid from the measured per-fit time.
# Input: n_pairs, n_seeds, t_probe (seconds for one fit).
# Output: estimated seconds, assuming jobs are spread over MC_CORES.
estimate_wall <- function(n_pairs, n_seeds, t_probe) {
  n_jobs <- n_pairs * n_seeds
  cores <- min(MC_CORES, n_jobs)
  ceiling(n_jobs / cores) * t_probe * TIME_SAFETY
}

# Choose the largest feasible grid within the remaining time budget.
# Input: remaining_sec (seconds left before the deadline); t_probe.
# Output: list(n_pairs, n_seeds, status). The grid is forced verbatim when
# FORCE_GRID is set; otherwise the order 10x5 -> 10x3 -> 5x3 is tested.
choose_grid <- function(remaining_sec, t_probe) {
  if (!is.null(FORCE_GRID)) {
    return(list(n_pairs = FORCE_GRID$n_pairs, n_seeds = FORCE_GRID$n_seeds,
                status = "forced"))
  }
  grids <- list(
    c(N_PAIRS_FULL, N_SEEDS_FULL),
    c(N_PAIRS_FULL, MIN_SEEDS),
    c(MIN_PAIRS, MIN_SEEDS)
  )
  for (k in seq_along(grids)) {
    g <- grids[[k]]
    if (estimate_wall(g[1], g[2], t_probe) <= remaining_sec) {
      tag <- if (k == 1) "full" else "downgraded"
      return(list(n_pairs = g[1], n_seeds = g[2],
                  status = paste0(tag, "_", g[1], "x", g[2])))
    }
  }
  # Even the smallest grid is predicted to exceed the budget; run it anyway and
  # flag the overrun so the report stays truthful.
  list(n_pairs = MIN_PAIRS, n_seeds = MIN_SEEDS,
       status = paste0("budget_exceeded_", MIN_PAIRS, "x", MIN_SEEDS))
}


## ============== Steps 2-6 for one dataset (the orchestration) ===============
# Run the complete experiment on one dataset.
# Input: ct (dataset name); cfg (its folder); deadline (global time limit).
# Output: list(runs, summary, repro, dataset_summary), the four tables later
# combined across datasets.
run_dataset <- function(ct, cfg, deadline) {
  message("\n================  ", ct, "  ================")
  # ---- Step 1: load the six cached objects.
  env <- load_env(cfg)
  tau_v <- as.numeric(env$tau[ct])
  message("tau = ", round(tau_v, 4), ", ",
          nrow(env$training_set[[ct]][["all"]]), " training pairs cached.")
  
  # ---- Step 2: evenly spaced pairs; the probe is the middle pair.
  full_pairs <- even_pairs(env$training_set, ct, N_PAIRS_FULL)
  probe_id <- full_pairs[ceiling(length(full_pairs) / 2)]
  message("Probe fit on ", probe_id, " with seed ", SEEDS[1], ".")
  t_probe0 <- Sys.time()
  probe_row <- fit_job_safe(ct, probe_id, SEEDS[1], env)
  t_probe <- as.numeric(difftime(Sys.time(), t_probe0, units = "secs"))
  message("Probe took ", round(t_probe / 60, 2), " min.")
  
  # ---- Step 3: choose the largest grid that fits the remaining time.
  remaining <- as.numeric(difftime(deadline, Sys.time(), units = "secs"))
  grid <- choose_grid(remaining, max(t_probe, 1))
  n_pairs <- grid$n_pairs
  n_seeds <- grid$n_seeds
  # Ensure the probe pair is part of the final grid (reuse its result).
  pairs <- even_pairs(env$training_set, ct, n_pairs)
  if (!(probe_id %in% pairs)) {
    pairs[ceiling(n_pairs / 2)] <- probe_id
    pairs <- unique(pairs)
  }
  # Reorder selected pairs by descending prior strength.
  all_df <- env$training_set[[ct]][["all"]]
  all_ids <- rownames(all_df)[order(all_df$theta_i, decreasing = TRUE)]
  pairs <- unique(all_ids[all_ids %in% pairs])
  use_seeds <- SEEDS[seq_len(n_seeds)]
  message("Grid: ", length(pairs), " pairs x ", n_seeds, " seeds = ",
          length(pairs) * n_seeds, " jobs on ",
          min(MC_CORES, length(pairs) * n_seeds), " cores (", grid$status, ").")
  
  # ---- Step 4: build the pair x seed job list and run it in parallel.
  jobs <- expand.grid(pair_id = pairs, seed = use_seeds,
                      stringsAsFactors = FALSE)
  done_key <- paste(probe_row$pair_id, probe_row$seed, sep = "__")
  jobs <- jobs[paste(jobs$pair_id, jobs$seed, sep = "__") != done_key, ,
               drop = FALSE]   # the probe pair x seed is already fitted
  if (nrow(jobs) > 0) {
    run_one <- function(i) fit_job_safe(ct, jobs$pair_id[i], jobs$seed[i], env)
    raw <- if (MC_CORES > 1) {
      # mc.set.seed = FALSE: every GA uses its own explicit seed, so results
      # are reproducible despite parallel execution.
      parallel::mclapply(seq_len(nrow(jobs)), run_one,
                         mc.cores = min(MC_CORES, nrow(jobs)),
                         mc.preschedule = TRUE, mc.set.seed = FALSE)
    } else lapply(seq_len(nrow(jobs)), run_one)
    # Convert any crashed worker (no data frame returned) into an error row.
    job_rows <- lapply(seq_along(raw), function(i) {
      if (is.data.frame(raw[[i]])) return(raw[[i]])
      pid <- jobs$pair_id[i]
      meta <- env$training_set[[ct]][["all"]][pid, c("TF", "TG", "theta_i")]
      cond <- attr(raw[[i]], "condition")
      msg <- if (is.null(cond)) "worker failed" else conditionMessage(cond)
      data.frame(
        dataset = ct, pair_id = pid, TF = meta$TF, TG = meta$TG,
        theta_i = meta$theta_i, seed = jobs$seed[i], theta_p = NA_real_,
        fit_loss = NA_real_, bin_tau = NA_integer_, qc_pass = NA,
        significant = NA_integer_, best_epoch = NA_integer_,
        stop_cause = "WORKER_ERROR", runtime_sec = NA_real_, error = msg,
        stringsAsFactors = FALSE
      )
    })
  } else job_rows <- list()
  runs <- do.call(rbind, c(list(probe_row), job_rows))
  
  # pair_rank = 1 is the strongest prior; used as the x-axis order in Figure 1.
  pair_rank <- setNames(seq_along(pairs), pairs)
  runs$pair_rank <- as.integer(pair_rank[runs$pair_id])
  runs <- runs[order(runs$pair_rank, runs$seed), ]
  rownames(runs) <- NULL
  saveRDS(runs, file.path(OUT_DIR, paste0("E4_runs_", ct, ".rds")))
  
  # ---- Step 5: seed-123 reproduction of the published model_train_all.
  # Only meaningful under the production settings (SMOKE_TEST = FALSE).
  ref_bpr <- env$ref_all[[ct]]$best_pred_result
  repro <- do.call(rbind, lapply(pairs, function(pid) {
    new <- runs[runs$pair_id == pid & runs$seed == SEEDS[1], ]
    ref_theta <- as.numeric(ref_bpr[pid, "theta_p_branch1"])
    ref_loss  <- as.numeric(ref_bpr[pid, "fit_loss_branch1"])
    d_theta <- abs(new$theta_p - ref_theta)
    rel_loss <- if (is.finite(ref_loss) && abs(ref_loss) > 1e-12) {
      abs(new$fit_loss - ref_loss) / abs(ref_loss)
    } else abs(new$fit_loss - ref_loss)
    data.frame(dataset = ct, pair_id = pid,
               ref_theta_p = ref_theta, seed123_theta_p = new$theta_p,
               abs_err_theta = d_theta,
               ref_fit_loss = ref_loss, seed123_fit_loss = new$fit_loss,
               rel_err_loss = rel_loss,
               pass = isTRUE(d_theta <= 1e-4 && rel_loss <= 1e-3),
               stringsAsFactors = FALSE)
  }))
  
  # ---- Step 6a: per-pair summary (one row per pair across seeds).
  summary <- do.call(rbind, lapply(pairs, function(pid) {
    z <- runs[runs$pair_id == pid, ]
    z <- z[order(z$seed), ]
    th <- z$theta_p
    data.frame(
      dataset = ct, pair_rank = pair_rank[pid], pair_id = pid,
      TF = z$TF[1], TG = z$TG[1], theta_i = z$theta_i[1],
      n_seeds = nrow(z),
      theta_mean = mean(th, na.rm = TRUE),      # mean predicted strength
      theta_sd = stats::sd(th, na.rm = TRUE),
      theta_cv = safe_cv(th),                   # key stability metric
      theta_min = min(th, na.rm = TRUE),
      theta_max = max(th, na.rm = TRUE),
      theta_range = diff(range(th, na.rm = TRUE)), # absolute spread on 0..1
      loss_cv = safe_cv(z$fit_loss),
      bin_tau_agreement = call_agreement(z$bin_tau),  # tau threshold call
      call_agreement = call_agreement(z$significant), # final edge call
      n_significant = sum(z$significant == 1, na.rm = TRUE), # seeds calling edge
      straddles_tau = min(th, na.rm = TRUE) < tau_v &&
        max(th, na.rm = TRUE) >= tau_v,        # TRUE = edge call flips
      median_best_epoch = stats::median(z$best_epoch, na.rm = TRUE),
      stringsAsFactors = FALSE)
  }))
  
  # ---- Step 6b: dataset-level summary.
  # Edge set of each seed = the pairs called significant by that seed.
  edge_sets <- lapply(use_seeds, function(s) {
    z <- runs[runs$seed == s, ]
    z$pair_id[!is.na(z$significant) & z$significant == 1]
  })
  # Mean pairwise Jaccard across all seed pairs.
  pair_jac <- c()
  for (a in seq_len(length(use_seeds) - 1)) {
    for (b in (a + 1):length(use_seeds)) {
      pair_jac <- c(pair_jac, jaccard(edge_sets[[a]], edge_sets[[b]]))
    }
  }
  dataset_summary <- data.frame(
    dataset = ct, tau = tau_v, n_pairs = length(pairs), n_seeds = n_seeds,
    median_cv_theta = median(summary$theta_cv, na.rm = TRUE),  # headline
    max_cv_theta = max(summary$theta_cv, na.rm = TRUE),        # worst pair
    median_cv_loss = median(summary$loss_cv, na.rm = TRUE),
    min_call_agreement = min(summary$call_agreement, na.rm = TRUE),
    n_pairs_call_disagreement = sum(summary$call_agreement < 1, na.rm = TRUE),
    n_pairs_straddle_tau = sum(summary$straddles_tau, na.rm = TRUE),
    mean_edge_set_jaccard = mean(pair_jac),                    # edge stability
    reproduction_all_pass = all(repro$pass),
    reproduction_max_abs_err_theta = max(repro$abs_err_theta),
    reproduction_max_rel_err_loss = max(repro$rel_err_loss),
    wall_time_min = round(sum(runs$runtime_sec, na.rm = TRUE) /
                            max(1, min(MC_CORES, nrow(runs))) / 60, 2),
    grid_status = grid$status,
    stringsAsFactors = FALSE)
  
  message(ct, " done: median CV(theta_p) = ",
          round(dataset_summary$median_cv_theta, 4),
          ", min edge-call agreement = ",
          dataset_summary$min_call_agreement,
          ", seed-123 reproduction all pass = ",
          dataset_summary$reproduction_all_pass, ".")
  list(runs = runs, summary = summary, repro = repro,
       dataset_summary = dataset_summary)
}


## ============================= Run all datasets ============================
deadline <- Sys.time() + TIME_BUDGET_HOURS * 3600
all_res <- list()
for (ct in RUN_DATASETS) {
  t1 <- Sys.time()
  all_res[[ct]] <- run_dataset(ct, DATASETS[[ct]], deadline)
  t2 <- Sys.time()
  saveRDS(all_res, file.path(OUT_DIR, "all_res.rds"))
  message(ct, " run time: ", format(t2 - t1), ".")
}
# ================  A549  ================
# tau = 0.2459, 85 training pairs cached.
# Probe fit on CEBPB_to_ASCC3 with seed 123.
# Probe took 0.6 min.
# Grid: 10 pairs x 5 seeds = 50 jobs on 50 cores (full_10x5).
# A549 done: median CV(theta_p) = 0.0187, min edge-call agreement = 1, seed-123 reproduction all pass = TRUE.
# A549 run time: 5.401873 mins.
# 
# ================  GM12878  ================
# tau = 0.1831, 76 training pairs cached.
# Probe fit on E2F7_to_MCM6 with seed 123.
# Probe took 2.93 min.
# Grid: 10 pairs x 5 seeds = 50 jobs on 50 cores (full_10x5).
# GM12878 done: median CV(theta_p) = 0.0402, min edge-call agreement = 1, seed-123 reproduction all pass = TRUE.
# GM12878 run time: 8.813803 mins.
# 
# ================  BM  ================
# tau = 0.1851, 92 training pairs cached.
# Probe fit on TRPS1_to_TRPS1 with seed 123.
# Probe took 0.61 min.
# Grid: 10 pairs x 5 seeds = 50 jobs on 50 cores (full_10x5).
# BM done: median CV(theta_p) = 0.0107, min edge-call agreement = 1, seed-123 reproduction all pass = TRUE.
# BM run time: 7.184505 mins.


## =================== Combine the four tables across datasets ===============
# runs_all: RAW results. One row per dataset x pair x seed (3 x 10 x 5 = 150
# rows in the full design).
#   dataset      : dataset name
#   pair_id      : "TF_to_TG" label
#   TF, TG       : the two genes of the pair
#   theta_i      : PRIOR strength (input, used to rank pairs)
#   seed         : GA seed
#   theta_p      : PREDICTED strength after training (main continuous output)
#   fit_loss     : pure prediction error (smaller is better)
#   bin_tau      : 1 if theta_p >= tau
#   qc_pass      : TRUE if theta_p >= 0.05 (poor-fit gate)
#   significant  : final edge call = 1 only when bin_tau == 1 and qc_pass
#   best_epoch   : epoch at which the best fit was reached
#   stop_cause   : why training stopped
#   runtime_sec  : fit time in seconds
#   error        : error message when the fit failed ("" otherwise)
#   pair_rank    : rank of the pair by descending theta_i (1 = strongest)
runs_all <- do.call(rbind, lapply(all_res, function(x) x$runs))

# Optional console checks (they only print; remove if not needed).
class(runs_all); dim(runs_all); head(runs_all)
table(runs_all$dataset)     # rows per dataset (50 in the full 10x5 design)
table(runs_all$seed)        # rows per seed
table(runs_all$bin_tau)     # how many seeds cross tau
table(runs_all$qc_pass)     # how many seeds pass the poor-fit gate
table(runs_all$significant) # final edge calls: 1 = edge, 0 = no edge
table(runs_all$stop_cause)  # stop reasons; in the production run "maximum
# training epochs reached" should be rare (it is common in the smoke test)
table(runs_all$pair_rank)   # each pair should have exactly n_seeds rows

# summary_all: PER-PAIR stability. One row per dataset x pair (30 rows).
# Read it as follows:
#   theta_cv       : spread of theta_p across seeds; small = stable,
#                    e.g. 0.02 means about 2% variation around the mean
#   theta_range    : max(theta_p) - min(theta_p); absolute spread on 0..1
#   loss_cv        : spread of fit_loss (typically far smaller than theta_cv)
#   bin_tau_agreement : agreement of the tau threshold call alone
#   call_agreement : 1.0 = all five seeds give the same edge call;
#                    0.8 = 4 of 5 agree; 0.6 = 3 of 5 agree
#   n_significant  : number of seeds (0..5) that call the pair an edge
#   straddles_tau  : TRUE when seeds fall on both sides of tau, i.e. the edge
#                    decision is seed-sensitive for this pair
#   median_best_epoch : training effort, for reference only
summary_all <- do.call(rbind, lapply(all_res, function(x) x$summary))
summary_all
#                          dataset pair_rank          pair_id     TF      TG   theta_i n_seeds theta_mean    theta_sd    theta_cv theta_min theta_max theta_range      loss_cv bin_tau_agreement call_agreement n_significant straddles_tau median_best_epoch
# A549.PBX3_to_TRAF2          A549         1    PBX3_to_TRAF2   PBX3   TRAF2 0.9133666       5  0.9183466 0.009272352 0.010096790 0.9072759 0.9319004 0.024624563 0.0005635432                 1              1             5         FALSE                 3
# A549.CREM_to_APTX           A549         2     CREM_to_APTX   CREM    APTX 0.9126604       5  0.8789985 0.073309684 0.083401373 0.8009537 0.9900000 0.189046271 0.0030492282                 1              1             5         FALSE                18
# A549.FOSL2_to_LMNA          A549         3    FOSL2_to_LMNA  FOSL2    LMNA 0.9121652       5  0.8825380 0.111878422 0.126768963 0.7521395 0.9900000 0.237860503 0.0004407218                 1              1             5         FALSE                 3
# A549.TFE3_to_RBM15          A549         4    TFE3_to_RBM15   TFE3   RBM15 0.9107224       5  0.9079234 0.001569331 0.001728484 0.9056426 0.9094336 0.003791006 0.0117166745                 1              1             5         FALSE                 5
# A549.CEBPB_to_ASCC3         A549         5   CEBPB_to_ASCC3  CEBPB   ASCC3 0.9096866       5  0.9150582 0.013507681 0.014761555 0.8984578 0.9330589 0.034601104 0.0100968293                 1              1             5         FALSE                 6
# A549.SOX13_to_LMNA          A549         6    SOX13_to_LMNA  SOX13    LMNA 0.9075986       5  0.8729955 0.083500893 0.095648706 0.7829566 0.9752630 0.192306438 0.0004305955                 1              1             5         FALSE                 5
# A549.E2F6_to_DLGAP4         A549         7   E2F6_to_DLGAP4   E2F6  DLGAP4 0.9059503       5  0.9214817 0.020928492 0.022711781 0.8966237 0.9435527 0.046928957 0.0065625801                 1              1             5         FALSE                24
# A549.ZNF274_to_AFAP1        A549         8  ZNF274_to_AFAP1 ZNF274   AFAP1 0.9033382       5  0.6704507 0.028897494 0.043101596 0.6447533 0.7120572 0.067303851 0.0004559442                 1              1             5         FALSE                15
# A549.E2F7_to_MCM4           A549         9     E2F7_to_MCM4   E2F7    MCM4 0.8982129       5  0.9175461 0.003642572 0.003969906 0.9131732 0.9207808 0.007607575 0.0028364610                 1              1             5         FALSE                19
# A549.RREB1_to_PTPN2         A549        10   RREB1_to_PTPN2  RREB1   PTPN2 0.8811890       5  0.8811843 0.007519840 0.008533789 0.8712905 0.8895961 0.018305520 0.0017476206                 1              1             5         FALSE                 5
# GM12878.KLF3_to_RAPH1    GM12878         1    KLF3_to_RAPH1   KLF3   RAPH1 0.9135273       5  0.9352720 0.020546765 0.021968760 0.9092694 0.9542254 0.044956003 0.0023009941                 1              1             5         FALSE                 7
# GM12878.HNF4G_to_FTH1    GM12878         2    HNF4G_to_FTH1  HNF4G    FTH1 0.9127288       5  0.8966151 0.067375173 0.075143917 0.7996069 0.9561776 0.156570668 0.0037877213                 1              1             5         FALSE                 9
# GM12878.E2F8_to_PPM1D    GM12878         3    E2F8_to_PPM1D   E2F8   PPM1D 0.9116856       5  0.8902306 0.041421259 0.046528684 0.8508358 0.9499001 0.099064332 0.0044151129                 1              1             5         FALSE                22
# GM12878.XBP1_to_XBP1     GM12878         4     XBP1_to_XBP1   XBP1    XBP1 0.9107738       5  0.7744433 0.038436298 0.049630874 0.7124806 0.8176334 0.105152808 0.0014323369                 1              1             5         FALSE                22
# GM12878.E2F7_to_MCM6     GM12878         5     E2F7_to_MCM6   E2F7    MCM6 0.9100278       5  0.6362923 0.014800033 0.023259804 0.6239217 0.6599970 0.036075308 0.0069549167                 1              1             5         FALSE                30
# GM12878.STAT3_to_HELB    GM12878         6    STAT3_to_HELB  STAT3    HELB 0.9083519       5  0.9172284 0.014199657 0.015481049 0.9057709 0.9413956 0.035624668 0.0408030159                 1              1             5         FALSE                11
# GM12878.ASCL1_to_B3GNTL1 GM12878         7 ASCL1_to_B3GNTL1  ASCL1 B3GNTL1 0.9063763       5  0.9601764 0.016727307 0.017421078 0.9410084 0.9819103 0.040901958 0.0039304277                 1              1             5         FALSE                11
# GM12878.IRF2_to_MX1      GM12878         8      IRF2_to_MX1   IRF2     MX1 0.9038574       5  0.6908934 0.023460836 0.033957244 0.6561389 0.7105815 0.054442598 0.0028279483                 1              1             5         FALSE                17
# GM12878.CEBPG_to_PTPN12  GM12878         9  CEBPG_to_PTPN12  CEBPG  PTPN12 0.8992529       5  0.8447709 0.181387257 0.214717701 0.6146863 0.9900000 0.375313698 0.1487803468                 1              1             5         FALSE                 2
# GM12878.SP4_to_HS2ST1    GM12878        10    SP4_to_HS2ST1    SP4  HS2ST1 0.8858111       5  0.7761761 0.062794066 0.080901825 0.7380243 0.8872953 0.149270913 0.0019710076                 1              1             5         FALSE                26
# BM.NFYC_to_CYCS               BM         1     NFYC_to_CYCS   NFYC    CYCS 0.9140316       5  0.9224142 0.001223737 0.001326667 0.9208113 0.9238826 0.003071277 0.0078762195                 1              1             5         FALSE                11
# BM.TFDP1_to_FANCD2            BM         2  TFDP1_to_FANCD2  TFDP1  FANCD2 0.9135499       5  0.9074665 0.011433567 0.012599437 0.8910898 0.9194111 0.028321312 0.0072473217                 1              1             5         FALSE                 5
# BM.ZNF281_to_GOLGA1           BM         3 ZNF281_to_GOLGA1 ZNF281  GOLGA1 0.9126941       5  0.9208111 0.004577600 0.004971270 0.9131778 0.9255795 0.012401753 0.0012998679                 1              1             5         FALSE                10
# BM.E2F8_to_MCM3               BM         4     E2F8_to_MCM3   E2F8    MCM3 0.9124207       5  0.9057755 0.009821484 0.010843177 0.8944079 0.9202797 0.025871815 0.0007739264                 1              1             5         FALSE                 5
# BM.TRPS1_to_TRPS1             BM         5   TRPS1_to_TRPS1  TRPS1   TRPS1 0.9118202       5  0.9193272 0.002637481 0.002868925 0.9163322 0.9216427 0.005310576 0.0045476709                 1              1             5         FALSE                10
# BM.KLF3_to_FAM107B            BM         6  KLF3_to_FAM107B   KLF3 FAM107B 0.9109270       5  0.9737340 0.036371900 0.037353016 0.9086699 0.9900000 0.081330041 0.1619490439                 1              1             5         FALSE               100
# BM.EGR1_to_STAT1              BM         7    EGR1_to_STAT1   EGR1   STAT1 0.9104529       5  0.9120993 0.009652205 0.010582406 0.9004842 0.9234947 0.023010540 0.0107041660                 1              1             5         FALSE                 8
# BM.SP4_to_NUP205              BM         8    SP4_to_NUP205    SP4  NUP205 0.9093176       5  0.8832936 0.018471692 0.020912291 0.8712519 0.9153562 0.044104357 0.0023772982                 1              1             5         FALSE                 6
# BM.IRF8_to_TARBP1             BM         9   IRF8_to_TARBP1   IRF8  TARBP1 0.9073317       5  0.9162028 0.001298354 0.001417104 0.9140110 0.9171719 0.003160952 0.0601985994                 1              1             5         FALSE                 8
# BM.SATB1_to_ANKRD17           BM        10 SATB1_to_ANKRD17  SATB1 ANKRD17 0.9020393       5  0.8961359 0.010344478 0.011543425 0.8853479 0.9116632 0.026315272 0.0025717035                 1              1             5         FALSE                 5

# repro_all: SEED-123 REPRODUCTION. One row per pair, comparing the new
# seed-123 fit with the published model_train_all result:
#   ref_theta_p / seed123_theta_p : published vs refitted theta_p
#   abs_err_theta  : absolute theta difference (pass when <= 1e-4)
#   ref_fit_loss / seed123_fit_loss : published vs refitted fit_loss
#   rel_err_loss   : relative fit-loss difference (pass when <= 1e-3)
#   pass           : TRUE when both tolerances hold
# Under the production settings (SMOKE_TEST = FALSE), pass should be TRUE for
# essentially every pair; under the smoke test FALSE is expected because the
# published fit used 100 epochs / early stop 5.
repro_all <- do.call(rbind, lapply(all_res, function(x) x$repro))
repro_all
#            dataset          pair_id ref_theta_p seed123_theta_p abs_err_theta ref_fit_loss seed123_fit_loss rel_err_loss pass
# A549.1        A549    PBX3_to_TRAF2   0.9138299       0.9138299             0   0.04345204       0.04345204            0 TRUE
# A549.2        A549     CREM_to_APTX   0.8009537       0.8009537             0   0.03989301       0.03989301            0 TRUE
# A549.3        A549    FOSL2_to_LMNA   0.7521395       0.7521395             0   0.03121006       0.03121006            0 TRUE
# A549.4        A549    TFE3_to_RBM15   0.9094336       0.9094336             0   0.08475388       0.08475388            0 TRUE
# A549.5        A549   CEBPB_to_ASCC3   0.8984578       0.8984578             0   0.06096160       0.06096160            0 TRUE
# A549.6        A549    SOX13_to_LMNA   0.9199545       0.9199545             0   0.05158964       0.05158964            0 TRUE
# A549.7        A549   E2F6_to_DLGAP4   0.8966237       0.8966237             0   0.06401560       0.06401560            0 TRUE
# A549.8        A549  ZNF274_to_AFAP1   0.6447533       0.6447533             0   0.06179992       0.06179992            0 TRUE
# A549.9        A549     E2F7_to_MCM4   0.9131732       0.9131732             0   0.04399578       0.04399578            0 TRUE
# A549.10       A549   RREB1_to_PTPN2   0.8758632       0.8758632             0   0.05501548       0.05501548            0 TRUE
# GM12878.1  GM12878    KLF3_to_RAPH1   0.9092694       0.9092694             0   0.05359946       0.05359946            0 TRUE
# GM12878.2  GM12878    HNF4G_to_FTH1   0.7996069       0.7996069             0   0.02638676       0.02638676            0 TRUE
# GM12878.3  GM12878    E2F8_to_PPM1D   0.8508358       0.8508358             0   0.01908703       0.01908703            0 TRUE
# GM12878.4  GM12878     XBP1_to_XBP1   0.7124806       0.7124806             0   0.02022567       0.02022567            0 TRUE
# GM12878.5  GM12878     E2F7_to_MCM6   0.6246296       0.6246296             0   0.01499179       0.01499179            0 TRUE
# GM12878.6  GM12878    STAT3_to_HELB   0.9057709       0.9057709             0   0.04430699       0.04430699            0 TRUE
# GM12878.7  GM12878 ASCL1_to_B3GNTL1   0.9410084       0.9410084             0   0.01824279       0.01824279            0 TRUE
# GM12878.8  GM12878      IRF2_to_MX1   0.7105815       0.7105815             0   0.02673302       0.02673302            0 TRUE
# GM12878.9  GM12878  CEBPG_to_PTPN12   0.9900000       0.9900000             0   0.03468920       0.03468920            0 TRUE
# GM12878.10 GM12878    SP4_to_HS2ST1   0.7408206       0.7408206             0   0.02638289       0.02638289            0 TRUE
# BM.1            BM     NFYC_to_CYCS   0.9208113       0.9208113             0   0.01968760       0.01968760            0 TRUE
# BM.2            BM  TFDP1_to_FANCD2   0.9194111       0.9194111             0   0.04483241       0.04483241            0 TRUE
# BM.3            BM ZNF281_to_GOLGA1   0.9219445       0.9219445             0   0.04606722       0.04606722            0 TRUE
# BM.4            BM     E2F8_to_MCM3   0.9005993       0.9005993             0   0.05497546       0.05497546            0 TRUE
# BM.5            BM   TRPS1_to_TRPS1   0.9210837       0.9210837             0   0.11773605       0.11773605            0 TRUE
# BM.6            BM  KLF3_to_FAM107B   0.9900000       0.9900000             0   0.06786869       0.06786869            0 TRUE
# BM.7            BM    EGR1_to_STAT1   0.9004842       0.9004842             0   0.03905178       0.03905178            0 TRUE
# BM.8            BM    SP4_to_NUP205   0.8746168       0.8746168             0   0.03099172       0.03099172            0 TRUE
# BM.9            BM   IRF8_to_TARBP1   0.9169240       0.9169240             0   0.02764169       0.02764169            0 TRUE
# BM.10           BM SATB1_to_ANKRD17   0.8880979       0.8880979             0   0.02685235       0.02685235            0 TRUE

# ds_all: DATASET-LEVEL HEADLINE. One row per dataset; these are the numbers
# used in the supplement text:
#   median_cv_theta / max_cv_theta : typical / worst theta_p variability
#   median_cv_loss                 : typical fit-loss variability
#   min_call_agreement             : worst pair's edge-call agreement
#   n_pairs_call_disagreement      : pairs whose edge call varies across seeds
#   n_pairs_straddle_tau           : pairs crossing tau across seeds
#   mean_edge_set_jaccard          : pairwise overlap of seed-wise edge sets
#                                    (1 = identical edge sets)
#   reproduction_all_pass          : TRUE only if every pair is reproduced
#   reproduction_max_abs_err_theta : largest absolute theta deviation
#   reproduction_max_rel_err_loss  : largest relative loss deviation
#   wall_time_min, grid_status     : runtime bookkeeping and grid used
ds_all <- do.call(rbind, lapply(all_res, function(x) x$dataset_summary))
ds_all
#         dataset       tau n_pairs n_seeds median_cv_theta max_cv_theta median_cv_loss min_call_agreement n_pairs_call_disagreement n_pairs_straddle_tau mean_edge_set_jaccard reproduction_all_pass reproduction_max_abs_err_theta reproduction_max_rel_err_loss wall_time_min grid_status
# A549       A549 0.2458855      10       5      0.01873667   0.12676896    0.002292041                  1                         0                    0                     1                  TRUE                              0                             0          1.50   full_10x5
# GM12878 GM12878 0.1831323      10       5      0.04024296   0.21471770    0.003859075                  1                         0                    0                     1                  TRUE                              0                             0          2.45   full_10x5
# BM           BM 0.1851275      10       5      0.01071279   0.03735302    0.005897496                  1                         0                    0                     1                  TRUE                              0                             0          1.63   full_10x5


## =============================== Figures ====================================
ds_colors <- c(A549 = "#66C2A5", GM12878 = "#FC8D62", BM = "#8DA0CB")
tau_df <- ds_all[, c("dataset", "tau")]

# Figure 1 -- how to read it:
#   Each panel is one dataset; x = pairs ranked by prior strength (1 =
#   strongest), y = predicted theta_p; one color per seed; the dashed line is
#   tau. A stable pair shows five points tightly stacked in one column. Points
#   crossing the dashed line within one column are pairs whose edge decision
#   flips across seeds. Tight clusters on the strong-prior side and occasional
#   spread on the weak-prior side are the expected pattern.
p_theta <- ggplot(runs_all, aes(x = factor(pair_rank), y = theta_p,
                                color = factor(seed))) +
  geom_point(size = 1.6, position = position_jitter(width = 0.12, height = 0)) +
  geom_hline(data = tau_df, aes(yintercept = tau), linetype = 2,
             linewidth = 0.4) +
  facet_wrap(~ dataset, scales = "free_y", nrow = 1) +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Training pair (ranked by prior strength theta_i; 1 = strongest)",
       y = "Predicted strength theta_p", color = "GA seed",
       title = "Seed stability of predicted regulatory strength") +
  theme_bw(base_size = 11) +
  theme(legend.position = "right")
p_theta
ggsave(file.path(OUT_DIR, "E4_theta_p_by_pair.png"), p_theta,
       width = 6, height = 3, units = "in", dpi = 300)

# Figure 2 -- how to read it:
#   x = prior strength theta_i, y = CV of theta_p across seeds. Points near the
#   bottom (CV ~ 0) are stable; points high up are seed-sensitive. The desired
#   pattern is that high-CV points, if any, concentrate on the weak-prior side
#   (small theta_i). Such pairs are the most likely to be removed by the
#   downstream QC and terminal filter rather than entering the released network.
p_cv <- ggplot(summary_all, aes(x = theta_i, y = theta_cv, color = dataset)) +
  geom_point(size = 1.8) +
  scale_color_manual(values = ds_colors) +
  labs(x = "Prior strength theta_i", y = "CV of theta_p across seeds",
       color = "Dataset",
       title = "Stability across the prior-strength range") +
  theme_bw(base_size = 11)
p_cv
ggsave(file.path(OUT_DIR, "E4_cv_vs_prior.png"), p_cv,
       width = 5, height = 3, units = "in", dpi = 300)


## ============================== Save tables =================================
write.csv(runs_all,    file.path(OUT_DIR, "E4_seed_runs.csv"), row.names = FALSE)
write.csv(summary_all, file.path(OUT_DIR, "E4_seed_summary.csv"), row.names = FALSE)
write.csv(ds_all,      file.path(OUT_DIR, "E4_dataset_summary.csv"),
          row.names = FALSE)
write.csv(repro_all,   file.path(OUT_DIR, "E4_seed_reproduction.csv"),
          row.names = FALSE)

## ============================== Final checks ===============================
message("\n---------------- E4 dataset summary ----------------")
print(ds_all, row.names = FALSE)
# dataset       tau n_pairs n_seeds median_cv_theta max_cv_theta median_cv_loss min_call_agreement n_pairs_call_disagreement n_pairs_straddle_tau mean_edge_set_jaccard reproduction_all_pass reproduction_max_abs_err_theta reproduction_max_rel_err_loss wall_time_min grid_status
#    A549 0.2458855      10       5      0.01873667   0.12676896    0.002292041                  1                         0                    0                     1                  TRUE                              0                             0          1.50   full_10x5
# GM12878 0.1831323      10       5      0.04024296   0.21471770    0.003859075                  1                         0                    0                     1                  TRUE                              0                             0          2.45   full_10x5
#      BM 0.1851275      10       5      0.01071279   0.03735302    0.005897496                  1                         0                    0                     1                  TRUE                              0                             0          1.63   full_10x5

message("\nSeed-123 reproduction of the published run:")
print(repro_all[, c("dataset", "pair_id", "abs_err_theta",
                    "rel_err_loss", "pass")], row.names = FALSE)
# dataset          pair_id abs_err_theta rel_err_loss pass
#    A549    PBX3_to_TRAF2             0            0 TRUE
#    A549     CREM_to_APTX             0            0 TRUE
#    A549    FOSL2_to_LMNA             0            0 TRUE
#    A549    TFE3_to_RBM15             0            0 TRUE
#    A549   CEBPB_to_ASCC3             0            0 TRUE
#    A549    SOX13_to_LMNA             0            0 TRUE
#    A549   E2F6_to_DLGAP4             0            0 TRUE
#    A549  ZNF274_to_AFAP1             0            0 TRUE
#    A549     E2F7_to_MCM4             0            0 TRUE
#    A549   RREB1_to_PTPN2             0            0 TRUE
# GM12878    KLF3_to_RAPH1             0            0 TRUE
# GM12878    HNF4G_to_FTH1             0            0 TRUE
# GM12878    E2F8_to_PPM1D             0            0 TRUE
# GM12878     XBP1_to_XBP1             0            0 TRUE
# GM12878     E2F7_to_MCM6             0            0 TRUE
# GM12878    STAT3_to_HELB             0            0 TRUE
# GM12878 ASCL1_to_B3GNTL1             0            0 TRUE
# GM12878      IRF2_to_MX1             0            0 TRUE
# GM12878  CEBPG_to_PTPN12             0            0 TRUE
# GM12878    SP4_to_HS2ST1             0            0 TRUE
#      BM     NFYC_to_CYCS             0            0 TRUE
#      BM  TFDP1_to_FANCD2             0            0 TRUE
#      BM ZNF281_to_GOLGA1             0            0 TRUE
#      BM     E2F8_to_MCM3             0            0 TRUE
#      BM   TRPS1_to_TRPS1             0            0 TRUE
#      BM  KLF3_to_FAM107B             0            0 TRUE
#      BM    EGR1_to_STAT1             0            0 TRUE
#      BM    SP4_to_NUP205             0            0 TRUE
#      BM   IRF8_to_TARBP1             0            0 TRUE
#      BM SATB1_to_ANKRD17             0            0 TRUE

# The reproduction check is only meaningful under the production settings.
if (SMOKE_TEST) {
  message("SMOKE_TEST = TRUE: reproduction mismatches are expected because the ",
          "published fit used max_epochs = 100 / early_stop = 5. Set ",
          "SMOKE_TEST <- FALSE for the reportable run.")
} else if (!all(repro_all$pass)) {
  warning("Seed-123 reproduction mismatch; check E4_seed_reproduction.csv.")
}

if (any(runs_all$stop_cause %in% c("ERROR", "WORKER_ERROR"), na.rm = TRUE)) {
  warning("One or more fits failed; see the error column in E4_seed_runs.csv.")
}

message("\nFiles saved under ", OUT_DIR)
# Files saved under /mnt/i/scDERI_Paper/Supplement-Experiment/E4

# ==============================================================================
# E2a: TF coverage funnel
# ==============================================================================

# Purpose
#   Answers questions1: JASPAR/ChIP-Atlas cover only a fraction of all TFs, so
#   what happens to a TF missing from the catalogs? 
#   Answers questions2: which inputs are required before any result is possible? 
# 
#   It quantifies, layer by layer, where transcription factors are lost on the
#   three benchmark datasets (A549, GM12878, BM; KP4 is not included), and proves
#   by set identity that a TF outside the reference catalogs can never enter the 
#   final network and that no TF is imputed or substituted.

# Six nested TF sets (each built exactly as scDERI::get_TFs_from_priors builds it)
#   L1 genes entering the pipeline after QC and HVG selection
#   L2 L1 genes with a JASPAR2024 PWM (package default collections, @name)
#   L3 L2 genes that are also in the built-in ChIP-Atlas TF catalog
#   L4 candidate TFs returned by get_TFs_from_priors (expression QC plus
#      prior download/quality pass; saved as TFs.rds)
#   L5 TFs after matrix construction and interest-cell-type alignment
#   L6 TFs present in the final peGRN
#   A diagnostic set L_expr (expression-QC pass, before any prior download)
#   splits the L3 -> L4 loss into "expression QC" and "prior download/quality".

# Steps
#   1. Load JASPAR2024 with the exact package call and extract the built-in
#      ChIP-Atlas TF vector from the installed function body.
#   2. For each dataset, load the saved pipeline intermediates and build the six
#      nested sets (plus the diagnostic sets).
#   3. Run set-difference and nesting checks (the no-imputation proof).
#   4. Write the summary, long/wide funnel and mechanism tables and raw sets.
#   5. Draw two figures (counts and cumulative retention), both with a broken
#      y-axis so the small L2-L6 counts are readable next to L1 = 10,000.

# Outputs (written to OUT_DIR)
#   E2a_summary.csv            one row per dataset: every layer count and rate
#   E2a_funnel_long.csv        one row per dataset x layer
#   E2a_funnel_wide.csv        wide count matrix (layers x datasets)
#   E2a_mechanism_check.csv    every transition, lost TFs and set checks
#   E2a_funnel_sets.rds        raw named sets for full reproducibility
#   E2a_TF_coverage_counts.png  bar counts with a broken y-axis
#   E2a_cumulative_retention.png retained-vs-L1 lines with a broken y-axis

# Runtime: seconds to a few minutes. Pure set operations, no network access and
# no retraining. Requires the same R environment as scDERI_benchmark.R.

options(stringsAsFactors = FALSE)

suppressMessages({
  library(scDERI)
  library(qs)
  library(JASPAR2024)
  library(TFBSTools)
  library(ggplot2)
  library(ggbreak)
})


## ------------------------------ Configuration ------------------------------
# Root directory of the benchmark run (same as scDERI_benchmark.R).
BENCH_ROOT <- "/mnt/i/scDERI_Paper/3-benchmark/8-scDERI"
# Output directory for the supplementary experiment.
OUT_DIR    <- "/mnt/i/scDERI_Paper/Supplement-Experiment/E2a"

# Dataset registry: analysis folder and the matched cell types excluded from the
# ChIP-Atlas priors (kept here only for provenance; the sets are read from disk).
DATASETS <- list(
  A549 = list(
    folder   = "SNAREseq2_hg38_A549",
    excluded = "A549"
  ),
  GM12878 = list(
    folder   = "SNAREseq2_hg38_GM12878",
    excluded = "GM12878"
  ),
  BM = list(
    folder   = "10X_hg38_Human_bone_marrow_with_T_cell_depleted",
    excluded = c("bone marrow", "NCI-H929", "BLaER1", "OCI-LY1", "OCI-LY3",
                 "erythroblast", "SH-SY5Y", "NB4")
  )
)

# JASPAR2024 retrieval options, identical to the package default and the
# benchmark call (all ten collections, PWM, current versions, human).
JASPAR_SPECIES    <- "Homo sapiens"
JASPAR_COLLECTION <- c("CORE", "CNE", "PHYLOFACTS", "SPLICE", "POLII",
                       "FAM", "PBM", "PBM_HOMEO", "PBM_HLH", "UNVALIDATED")

# Expression-QC thresholds used in the benchmark (identical across datasets).
EXPR_QC <- list(
  min_nonzero_freq          = 0.032,
  min_nonzero_mean          = 1.104,
  min_overall_mean          = 0.046,
  min_nonzero_freq_quantile = 0.387,
  min_nonzero_mean_quantile = 0.382,
  min_overall_mean_quantile = 0.421
)

# Human-readable description of each funnel layer for tables and the figure.
LAYER_STAGE <- c(
  L1 = "Genes after QC and HVG selection",
  L2 = "With a JASPAR2024 PWM",
  L3 = "Also in the ChIP-Atlas prior catalog",
  L4 = "Candidate TFs (expression QC + prior pass)",
  L5 = "After matrix build and cell-type alignment",
  L6 = "TFs in the final peGRN"
)

# Counts printed by the benchmark log, used only as a soft cross-check.
EXPECTED <- list(
  A549    = c(L1 = 10000, L4 = 87, L5 = 85, L6 = 84),
  GM12878 = c(L1 = 10000, L4 = 79, L5 = 77, L6 = 76),
  BM      = c(L1 = 10000, L4 = 94, L5 = 93, L6 = 92)
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = TRUE)


## ------------------------------- File paths --------------------------------
# Each helper resolves one saved pipeline intermediate from a dataset folder.
# Input:  wd, the dataset's working directory.
# Output: a single absolute file path.

path_scrna <- function(wd) file.path(wd, "scRNAseq.qsave")
path_tf <- function(wd) file.path(
  wd, "1.1 Data Preprocessing - Get TFs From Database", "TFs.rds"
)
path_basic <- function(wd) file.path(
  wd, "1.3 Data Preprocessing - Get Expression And Activity Matrix",
  "Basic_Info.qsave"
)
path_pegrn <- function(wd) file.path(wd, "interest_cell_type_peGRN.qsave")
path_gpi <- function(wd) file.path(
  wd, "2.2 Data Processing - Cell Grouping",
  "interest_cell_type_genes_pseudotime_info.rds"
)


## ------------------------------ Small helpers ------------------------------
# Flatten a vector/list of TF names into a sorted unique character vector.
# Input:  x, a vector or list of TF symbols.
# Output: a sorted, unique character vector.
as_char <- function(x) {
  if (is.list(x)) x <- unlist(x)
  sort(unique(as.character(x)))
}

# Percentage n/d, returning NA when the denominator is not a single positive value.
# Input:  n (count), d (denominator).
# Output: percentage value (0-100) or NA_real_.
pct <- function(n, d) ifelse(length(d) == 1 && d > 0, 100 * n / d, NA_real_)

# Collapse the first five names of a vector for compact table reporting.
# Input:  x, a character vector of TF symbols.
# Output: a single character string with up to five names.
example_names <- function(x) paste(utils::head(x, 5), collapse = "; ")


## ----------------------------- Catalog layers -------------------------------
# Fetch all human JASPAR2024 PWM profiles with the exact package call.
# Input:  none (uses the JASPAR_* constants).
# Output: a TFBSTools PWMatrixList.
get_jaspar_pwm_list <- function() {
  db <- JASPAR2024::db(JASPAR2024::JASPAR2024())
  TFBSTools::getMatrixSet(
    x = db,
    opts = list(all = FALSE, species = JASPAR_SPECIES,
                collection = JASPAR_COLLECTION,
                all_versions = FALSE, matrixtype = "PWM")
  )
}

# Restrict JASPAR profiles to measured genes and return unique profile names.
# Input:  pwm, a PWMatrixList; genes, the measured gene symbols.
# Output: a sorted unique character vector of JASPAR profile names (set L2).
jaspar_names_present <- function(pwm, genes) {
  kept <- pwm
  kept@listData <- Filter(function(x) x@name %in% genes, kept@listData)
  sort(unique(vapply(kept@listData, function(x) x@name, character(1))))
}

# Extract the built-in ChIP-Atlas TF vector for one genome assembly directly
# from the installed scDERI function body. This avoids any network download
# while using the exact catalog the package uses.
# Input:  assembly, a genome assembly string (default "hg38").
# Output: a sorted unique character vector of ChIP-Atlas TF symbols.
extract_chip_atlas_names <- function(assembly = "hg38") {
  found <- NULL
  # Recursively walk the parsed function body to find the assignment
  #   TFs_in_ChIP_Atlas = c(...)
  # inside the branch matching the requested assembly.
  walk <- function(node) {
    if (is.call(node)) {
      if (identical(node[[1]], as.name("if")) && length(node) >= 3) {
        cond <- node[[2]]
        branch <- node[[3]]
        if (is.call(cond) && length(cond) == 3 &&
            is.character(cond[[3]]) && cond[[3]] == assembly && is.call(branch)) {
          for (k in seq_along(branch)) {
            el <- branch[[k]]
            if (is.call(el) && length(el) == 3 &&
                (identical(el[[1]], as.name("=")) ||
                 identical(el[[1]], as.name("<-"))) &&
                is.name(el[[2]]) &&
                as.character(el[[2]]) == "TFs_in_ChIP_Atlas") {
              found <<- eval(el[[3]], envir = baseenv())
            }
          }
        }
      }
      for (k in seq_along(node)) walk(node[[k]])
    }
  }
  walk(body(scDERI::get_TFs_from_priors))
  if (is.null(found) || length(found) == 0) {
    stop("Could not extract the built-in ChIP-Atlas TF list for ", assembly,
         ". Please supply the TFs_in_ChIP_Atlas vector for this assembly.")
  }
  sort(unique(as.character(found)))
}

# Replicate the expression-QC stage of get_TFs_from_priors without any download:
# three absolute thresholds and three genome-wide quantile thresholds.
# Input:  scRNA, a Seurat object; universe, the candidate TF set to filter (L3).
# Output: a sorted character vector of TFs passing all six thresholds (L_expr).
expression_qc_passed <- function(scRNA, universe) {
  counts <- Seurat::GetAssayData(scRNA, layer = "counts")
  n_cells <- ncol(counts)
  row_sums  <- Matrix::rowSums(counts)
  n_nonzero <- Matrix::rowSums(counts != 0)
  # Per-gene expression statistics computed on raw counts.
  feat <- data.frame(
    nonzero_freq  = n_nonzero / n_cells,
    nonzero_mean  = ifelse(n_nonzero == 0, 0, row_sums / n_nonzero),
    overall_mean  = row_sums / n_cells,
    row.names     = rownames(counts),
    check.names   = FALSE
  )
  # Genome-wide percentile ranks (0 = lowest, 1 = highest), then restrict to the
  # candidate universe, exactly as the package does.
  feat[paste0(names(feat), "_quantile")] <-
    lapply(feat, function(x) rank(x) / length(x))
  feat <- feat[intersect(universe, rownames(feat)), , drop = FALSE]
  # A TF must satisfy all three absolute and all three quantile thresholds.
  keep <- with(feat,
               nonzero_freq          >= EXPR_QC$min_nonzero_freq &
                 nonzero_mean          >= EXPR_QC$min_nonzero_mean &
                 overall_mean          >= EXPR_QC$min_overall_mean &
                 nonzero_freq_quantile >= EXPR_QC$min_nonzero_freq_quantile &
                 nonzero_mean_quantile >= EXPR_QC$min_nonzero_mean_quantile &
                 overall_mean_quantile >= EXPR_QC$min_overall_mean_quantile)
  sort(rownames(feat)[keep])
}


## --------------------------- Per-dataset funnel ----------------------------
# Build every funnel set for one dataset from its saved intermediates.
# Input:  ds_name, dataset label; cfg, its registry entry; pwm, the shared
#         JASPAR PWMatrixList; chip_all, the built-in ChIP-Atlas vector.
# Output: a list with `sets` (all named TF sets) and `audit` (cross-check notes).
build_funnel <- function(ds_name, cfg, pwm, chip_all) {
  wd <- file.path(BENCH_ROOT, cfg$folder)
  message("\n================  ", ds_name, "  ================")
  audit <- character()
  
  # L1: genes entering the pipeline after QC and feature selection.
  scRNA <- qs::qread(path_scrna(wd))
  L1 <- sort(rownames(scRNA[["RNA"]]))
  message("L1 genes after QC/HVG: ", length(L1))
  
  # L2: JASPAR2024 profile names present in the data.
  L2 <- jaspar_names_present(pwm, L1)
  message("L2 genes with a JASPAR2024 PWM: ", length(L2))
  
  # ChIP-Atlas names present in the data; L3 is the intersection of catalogs.
  chip_in_data <- sort(intersect(chip_all, L1))
  L3 <- sort(intersect(L2, chip_in_data))
  message("L3 TFs also in ChIP-Atlas catalog: ", length(L3),
          " (", length(chip_in_data), " ChIP-Atlas TFs measured in data)")
  
  # Diagnostic: expression-QC pass before any ChIP prior download/quality step.
  L_expr <- expression_qc_passed(scRNA, L3)
  message("L_expr expression-QC pass: ", length(L_expr))
  
  # L4: saved candidate TFs (expression QC and prior download/quality pass).
  L4 <- as_char(readRDS(path_tf(wd)))
  message("L4 candidate TFs (TFs.rds): ", length(L4))
  
  # L5: TFs after matrix construction and cell-type alignment.
  basic <- qs::qread(path_basic(wd))
  L5 <- as_char(basic[["TFs"]])
  message("L5 TFs after matrix/alignment: ", length(L5))
  
  # L6: final peGRN TFs; L6_loose is the same set before terminal filtering.
  pegrn <- qs::qread(path_pegrn(wd))[[1]]
  L6_loose <- as_char(pegrn[["loose"]][["TF"]])
  L6       <- as_char(pegrn[["loose_filtered_simplified"]][["TF"]])
  message("L6 final peGRN TFs: ", length(L6),
          " (", length(L6_loose), " before terminal filtering)")
  
  # Optional audit: TFs retained after pseudotime branch assignment.
  L_branch <- NULL
  if (file.exists(path_gpi(wd))) {
    gpi <- readRDS(path_gpi(wd))
    L_branch <- tryCatch(
      as_char(Reduce(union, gpi[[1]][["Branches_TFs"]])),
      error = function(e) NULL)
  }
  
  # Soft cross-check of the saved layers against the benchmark log.
  for (k in names(EXPECTED[[ds_name]])) {
    got <- length(get(k))
    if (!is.na(got) && got != EXPECTED[[ds_name]][k]) {
      audit <- c(audit, paste0("Count mismatch at ", k, ": got ", got,
                               ", benchmark log ", EXPECTED[[ds_name]][k], "."))
    }
  }
  
  list(sets = list(L1 = L1, L2 = L2, L3 = L3, L_expr = L_expr, L4 = L4,
                   L5 = L5, L6 = L6, L6_loose = L6_loose,
                   L_branch = L_branch, chip_in_data = chip_in_data),
       audit = audit)
}


## --------------------------- Catalog totals (global) ------------------------
# Database-wide catalog sizes, not restricted to any dataset. They define the
# candidate universe; only the intersection with the measured 10,000 genes is
# carried forward. The numbers are computed, never hard-coded.
pwm_all <- get_jaspar_pwm_list()
class(pwm_all); length(pwm_all); head(names(pwm_all))
# [1] "PWMatrixList"
# attr(,"package")
# [1] "TFBSTools"
# [1] 720
# [1] "MA0069.1" "MA0071.1" "MA0074.1" "MA0101.1" "MA0107.1" "MA0115.1"

jaspar_total <- sort(unique(vapply(pwm_all@listData, function(x) x@name,
                                   character(1))))
class(jaspar_total); length(jaspar_total); head(jaspar_total)
# [1] "character"
# [1] 675
# [1] "ALX3" "ARGFX" "ARNT2" "ARNT::HIF1A" "ASCL1" "ATF2"

chip_total <- extract_chip_atlas_names("hg38")
class(chip_total); length(chip_total); head(chip_total)
# [1] "character"
# [1] 1766
# [1] "AATF" "ABI1" "ACSS2" "ACTR5" "ADAR" "ADNP"

message("Catalog totals: ", length(jaspar_total), " human JASPAR2024 PWM names; ",
        length(chip_total), " hg38 ChIP-Atlas TF names.")
# Catalog totals: 675 human JASPAR2024 PWM names; 1766 hg38 ChIP-Atlas TF names.


## --------------------------- Main loop over datasets -----------------------
funnel <- list()
for (ds in names(DATASETS)) {
  funnel[[ds]] <- build_funnel(ds, DATASETS[[ds]], pwm_all, chip_total)
}
# ================  A549  ================
# L1 genes after QC/HVG: 10000
# L2 genes with a JASPAR2024 PWM: 327
# L3 TFs also in ChIP-Atlas catalog: 248 (811 ChIP-Atlas TFs measured in data)
# L_expr expression-QC pass: 109
# L4 candidate TFs (TFs.rds): 87
# L5 TFs after matrix/alignment: 85
# L6 final peGRN TFs: 84 (84 before terminal filtering)
#
# ================  GM12878  ================
# L1 genes after QC/HVG: 10000
# L2 genes with a JASPAR2024 PWM: 299
# L3 TFs also in ChIP-Atlas catalog: 212 (585 ChIP-Atlas TFs measured in data)
# L_expr expression-QC pass: 105
# L4 candidate TFs (TFs.rds): 79
# L5 TFs after matrix/alignment: 77
# L6 final peGRN TFs: 76 (76 before terminal filtering)
#
# ================  BM  ================
# L1 genes after QC/HVG: 10000
# L2 genes with a JASPAR2024 PWM: 253
# L3 TFs also in ChIP-Atlas catalog: 214 (950 ChIP-Atlas TFs measured in data)
# L_expr expression-QC pass: 122
# L4 candidate TFs (TFs.rds): 94
# L5 TFs after matrix/alignment: 93
# L6 final peGRN TFs: 92 (92 before terminal filtering)

# Reading the funnel: the two catalog gates dominate. L1 starts at 10,000 genes;
# requiring a JASPAR PWM is the first and largest reduction (L2 ~= 250-330), and
# intersecting ChIP-Atlas reduces it further (L3 ~= 210-250). Expression QC and
# prior quality then halve the catalog pool (L3 -> L_expr -> L4), whereas the
# downstream matrix, alignment and terminal-filter steps lose almost no TF
# (L5 -> L6). Thus attrition is a catalog/expression effect, not a modeling one.


## ------------------------------- Long table --------------------------------
# Build the long-format funnel: one row per dataset x layer.
# Input:  funnel, the named list returned by build_funnel.
# Output: a data.frame with count, retention vs previous layer and retention
#         vs L1 for every layer of every dataset.
make_long <- function(funnel) {
  rows <- list()
  for (ds in names(funnel)) {
    sets <- funnel[[ds]]$sets
    n1 <- length(sets$L1)
    core <- names(LAYER_STAGE)
    for (i in seq_along(core)) {
      nm <- core[i]
      n <- length(sets[[nm]])
      # No previous-layer value for L1.
      prev <- if (i == 1) NA_real_ else pct(n, length(sets[[core[i - 1]]]))
      rows[[length(rows) + 1]] <- data.frame(
        dataset              = ds,
        layer                = nm,
        stage                = LAYER_STAGE[i],
        n_TFs                = n,
        retained_vs_prev_pct = round(prev, 2),
        retained_vs_L1_pct   = round(pct(n, n1), 2),
        stringsAsFactors     = FALSE)
    }
  }
  do.call(rbind, rows)
}

long_df <- make_long(funnel)
long_df
#     dataset layer                                      stage n_TFs retained_vs_prev_pct retained_vs_L1_pct
# L1     A549    L1           Genes after QC and HVG selection 10000                   NA             100.00
# L2     A549    L2                      With a JASPAR2024 PWM   327                 3.27               3.27
# L3     A549    L3       Also in the ChIP-Atlas prior catalog   248                75.84               2.48
# L4     A549    L4 Candidate TFs (expression QC + prior pass)    87                35.08               0.87
# L5     A549    L5 After matrix build and cell-type alignment    85                97.70               0.85
# L6     A549    L6                     TFs in the final peGRN    84                98.82               0.84
# L11 GM12878    L1           Genes after QC and HVG selection 10000                   NA             100.00
# L21 GM12878    L2                      With a JASPAR2024 PWM   299                 2.99               2.99
# L31 GM12878    L3       Also in the ChIP-Atlas prior catalog   212                70.90               2.12
# L41 GM12878    L4 Candidate TFs (expression QC + prior pass)    79                37.26               0.79
# L51 GM12878    L5 After matrix build and cell-type alignment    77                97.47               0.77
# L61 GM12878    L6                     TFs in the final peGRN    76                98.70               0.76
# L12      BM    L1           Genes after QC and HVG selection 10000                   NA             100.00
# L22      BM    L2                      With a JASPAR2024 PWM   253                 2.53               2.53
# L32      BM    L3       Also in the ChIP-Atlas prior catalog   214                84.58               2.14
# L42      BM    L4 Candidate TFs (expression QC + prior pass)    94                43.93               0.94
# L52      BM    L5 After matrix build and cell-type alignment    93                98.94               0.93
# L62      BM    L6                     TFs in the final peGRN    92                98.92               0.92

# Columns: dataset/layer/stage identify the row; n_TFs is the TF count;
# retained_vs_prev_pct is the stepwise retention (e.g. L3/L2), and
# retained_vs_L1_pct is the cumulative retention from L1. The pattern is
# consistent across datasets: the first step (L1 -> L2) retains only ~2.5-3.3%,
# the catalog intersection and expression QC keep shrinking the pool, while the
# final modeling steps (L5 -> L6) retain ~98-99%.


## ------------------------------- Wide table ---------------------------------
# Reshape counts to a wide matrix: rows are the six layers, columns datasets.
# This is the compact funnel table for the manuscript.
wide_df <- reshape(
  long_df[, c("dataset", "layer", "n_TFs")],
  idvar = "layer", timevar = "dataset", direction = "wide")
names(wide_df) <- sub("n_TFs\\.", "", names(wide_df))
wide_df$layer <- factor(wide_df$layer, levels = names(LAYER_STAGE))
wide_df <- as.data.frame(wide_df[order(wide_df$layer), ], row.names = NULL)
wide_df <- data.frame(
  layer = as.character(wide_df$layer),
  stage = LAYER_STAGE[as.character(wide_df$layer)],
  wide_df[, setdiff(names(wide_df), "layer"), drop = FALSE],
  row.names = NULL, check.names = FALSE)
wide_df

#   layer  stage                                            A549  GM12878   BM
#   L1     Genes after QC and HVG selection                10000    10000 10000
#   L2     With a JASPAR2024 PWM                             327      299   253
#   L3     Also in the ChIP-Atlas prior catalog              248      212   214
#   L4     Candidate TFs (expression QC + prior pass)         87       79    94
#   L5     After matrix build and cell-type alignment         85       77    93
#   L6     TFs in the final peGRN                             84       76    92

# Each row is a funnel layer and the three dataset columns give the counts side
# by side. The three datasets show the same funnel shape: a steep catalog-driven
# decline from L1 to L3/L4, then nearly flat counts through L5 and L6.


## ------------------------------- Summary row --------------------------------
# Build one summary row per dataset, including the diagnostic ChIP-in-data and
# expression-QC counts and the terminal-filter TF loss.
# Input:  funnel, the named list returned by build_funnel.
# Output: a one-row-per-dataset data.frame.
make_summary <- function(funnel) {
  rows <- list()
  for (ds in names(funnel)) {
    s <- funnel[[ds]]$sets
    rows[[length(rows) + 1]] <- data.frame(
      dataset                 = ds,
      n_L1_genes              = length(s$L1),
      n_L2_jaspar             = length(s$L2),
      n_chip_in_data          = length(s$chip_in_data),
      n_L3_catalog_intersect  = length(s$L3),
      n_L_expr_expression     = length(s$L_expr),
      n_L4_candidate          = length(s$L4),
      n_L5_matrix             = length(s$L5),
      n_L6_final              = length(s$L6),
      n_loose_pre_filter      = length(s$L6_loose),
      final_over_L1_pct       = round(pct(length(s$L6), length(s$L1)), 3),
      terminal_filter_TF_loss = length(setdiff(s$L6_loose, s$L6)),
      stringsAsFactors        = FALSE)
  }
  do.call(rbind, rows)
}

summary_df <- make_summary(funnel)
summary_df
#   dataset n_L1_genes n_L2_jaspar n_chip_in_data n_L3_catalog_intersect n_L_expr_expression n_L4_candidate n_L5_matrix n_L6_final n_loose_pre_filter final_over_L1_pct terminal_filter_TF_loss
# 1    A549      10000         327            811                    248                 109             87          85         84                 84              0.84                       0
# 2 GM12878      10000         299            585                    212                 105             79          77         76                 76              0.76                       0
# 3      BM      10000         253            950                    214                 122             94          93         92                 92              0.92                       0

# Key reads per row: n_chip_in_data (811/585/950) is much larger than n_L2_jaspar
# (327/299/253), so JASPAR, not ChIP-Atlas, is the narrower catalog gate on these
# data. final_over_L1_pct shows the network ultimately retains 0.76-0.92% of the
# input genes as TFs. terminal_filter_TF_loss is 0 for every dataset: the final
# edge filter removes edges but never removes a TF.


## ---------------------------- Mechanism checks ------------------------------
# Build the mechanism table: one row per transition or audit, with the counts,
# the lost TF examples and the set-identity result.
# Input:  funnel, the named list returned by build_funnel.
# Output: a data.frame whose status column is PASS/FAIL for every set check.
make_mechanism <- function(funnel) {
  rows <- list()
  # Append one row describing a single transition or audit.
  add <- function(ds, kind, transition, from, to, lost, identity, check) {
    rows[[length(rows) + 1]] <<- data.frame(
      dataset = ds, kind = kind, transition = transition,
      n_from = length(from), n_to = length(to),
      n_lost = length(lost),
      lost_examples = example_names(lost),
      identity_check = identity,
      status = if (is.na(check)) "NOT_CHECKED" else if (check) "PASS" else "FAIL",
      stringsAsFactors = FALSE)
  }
  
  for (ds in names(funnel)) {
    s <- funnel[[ds]]$sets
    
    # Six-layer transitions, with L3 -> L4 split into expression and prior loss.
    add(ds, "transition", "L1->L2 (JASPAR PWM)", s$L1, s$L2,
        setdiff(s$L1, s$L2), "L2 is a subset of L1", all(s$L2 %in% s$L1))
    add(ds, "transition", "L2->L3 (ChIP-Atlas)", s$L2, s$L3,
        setdiff(s$L2, s$L3), "L3 is a subset of L2", all(s$L3 %in% s$L2))
    add(ds, "transition", "L3->L_expr (expression QC)", s$L3, s$L_expr,
        setdiff(s$L3, s$L_expr), "expression-QC pass is a subset of L3",
        all(s$L_expr %in% s$L3))
    add(ds, "transition", "L_expr->L4 (prior download/quality)", s$L_expr, s$L4,
        setdiff(s$L_expr, s$L4), "candidate TFs are a subset of expression pass",
        all(s$L4 %in% s$L_expr))
    add(ds, "transition", "L4->L5 (matrix/alignment)", s$L4, s$L5,
        setdiff(s$L4, s$L5), "L5 is a subset of L4", all(s$L5 %in% s$L4))
    add(ds, "transition", "L5->L6 (branch + terminal)", s$L5, s$L6,
        setdiff(s$L5, s$L6), "L6 is a subset of L5", all(s$L6 %in% s$L5))
    
    # The terminal filter alone must not drop a TF (it removes only edges).
    add(ds, "audit", "loose->L6 (terminal filter)", s$L6_loose, s$L6,
        setdiff(s$L6_loose, s$L6),
        "terminal filtering removes no TF, only edges",
        length(setdiff(s$L6_loose, s$L6)) == 0)
    
    # Branch assignment explains the small L5 -> L6 loss when the object exists.
    if (!is.null(s$L_branch)) {
      add(ds, "audit", "L5->branch (pseudotime)", s$L5, s$L_branch,
          setdiff(s$L5, s$L_branch),
          "TFs lost at pseudotime branch assignment",
          all(s$L_branch %in% s$L5))
    }
    
    # No-imputation guarantees: every final TF must belong to both catalogs and
    # to the input genes. A zero difference means no out-of-catalog TF ever enters.
    add(ds, "audit", "final TFs \\ JASPAR", s$L6,
        intersect(s$L6, s$L2), setdiff(s$L6, s$L2),
        "no final TF lacks a JASPAR2024 PWM",
        length(setdiff(s$L6, s$L2)) == 0)
    add(ds, "audit", "final TFs \\ ChIP-Atlas", s$L6,
        intersect(s$L6, s$L3), setdiff(s$L6, s$L3),
        "no final TF is outside the ChIP-Atlas catalog",
        length(setdiff(s$L6, s$L3)) == 0)
    add(ds, "audit", "final TFs \\ L1", s$L6,
        intersect(s$L6, s$L1), setdiff(s$L6, s$L1),
        "no final TF is absent from the input genes",
        length(setdiff(s$L6, s$L1)) == 0)
  }
  do.call(rbind, rows)
}

mech_df <- make_mechanism(funnel)
mech_df
#    dataset       kind                          transition n_from n_to n_lost                          lost_examples                                identity_check status
# 1     A549 transition                 L1->L2 (JASPAR PWM)  10000  327   9673 A1BG-AS1; A1CF; A2M-AS1; A2ML1; A4GALT                          L2 is a subset of L1   PASS
# 2     A549 transition                 L2->L3 (ChIP-Atlas)    327  248     79     ARGFX; BARHL1; BARHL2; BARX2; CUX2                          L3 is a subset of L2   PASS
# 3     A549 transition          L3->L_expr (expression QC)    248  109    139       BACH2; BATF; BCL11A; BCL6; CEBPA          expression-QC pass is a subset of L3   PASS
# 4     A549 transition L_expr->L4 (prior download/quality)    109   87     22        ARNT2; BNC2; ESRRB; ETS2; FOXN3 candidate TFs are a subset of expression pass   PASS
# 5     A549 transition           L4->L5 (matrix/alignment)     87   85      2                            FOXP1; NFIA                          L5 is a subset of L4   PASS
# 6     A549 transition          L5->L6 (branch + terminal)     85   84      1                                   ETV1                          L6 is a subset of L5   PASS
# 7     A549      audit         loose->L6 (terminal filter)     84   84      0                                         terminal filtering removes no TF, only edges   PASS
# 8     A549      audit             L5->branch (pseudotime)     85   84      1                                   ETV1      TFs lost at pseudotime branch assignment   PASS
# 9     A549      audit                 final TFs \\ JASPAR     84   84      0                                                   no final TF lacks a JASPAR2024 PWM   PASS
# 10    A549      audit             final TFs \\ ChIP-Atlas     84   84      0                                        no final TF is outside the ChIP-Atlas catalog   PASS
# 11    A549      audit                     final TFs \\ L1     84   84      0                                           no final TF is absent from the input genes   PASS
# 12 GM12878 transition                 L1->L2 (JASPAR PWM)  10000  299   9701     A1BG; A1BG-AS1; A1CF; A2M; A2M-AS1                          L2 is a subset of L1   PASS
# 13 GM12878 transition                 L2->L3 (ChIP-Atlas)    299  212     87  ATOH7; BHLHA15; BHLHE22; BHLHE23; BSX                          L3 is a subset of L2   PASS
# 14 GM12878 transition          L3->L_expr (expression QC)    212  105    107    ARNT2; BACH2; BCL11A; BCL6; BHLHE41          expression-QC pass is a subset of L3   PASS
# 15 GM12878 transition L_expr->L4 (prior download/quality)    105   79     26          BNC2; E2F2; ELK3; EMX1; ESRRB candidate TFs are a subset of expression pass   PASS
# 16 GM12878 transition           L4->L5 (matrix/alignment)     79   77      2                           FOXP1; FOXP2                          L5 is a subset of L4   PASS
# 17 GM12878 transition          L5->L6 (branch + terminal)     77   76      1                                 TFAP2B                          L6 is a subset of L5   PASS
# 18 GM12878      audit         loose->L6 (terminal filter)     76   76      0                                         terminal filtering removes no TF, only edges   PASS
# 19 GM12878      audit             L5->branch (pseudotime)     77   76      1                                 TFAP2B      TFs lost at pseudotime branch assignment   PASS
# 20 GM12878      audit                 final TFs \\ JASPAR     76   76      0                                                   no final TF lacks a JASPAR2024 PWM   PASS
# 21 GM12878      audit             final TFs \\ ChIP-Atlas     76   76      0                                        no final TF is outside the ChIP-Atlas catalog   PASS
# 22 GM12878      audit                     final TFs \\ L1     76   76      0                                           no final TF is absent from the input genes   PASS
# 23      BM transition                 L1->L2 (JASPAR PWM)  10000  253   9747    A1BG; A1BG-AS1; A2M; A2M-AS1; A2ML1                          L2 is a subset of L1   PASS
# 24      BM transition                 L2->L3 (ChIP-Atlas)    253  214     39        CREB3L4; CUX2; EGR3; GLI3; HES6                          L3 is a subset of L2   PASS
# 25      BM transition          L3->L_expr (expression QC)    214  122     92   ARNT2; BATF3; BHLHE40; BHLHE41; BNC2          expression-QC pass is a subset of L3   PASS
# 26      BM transition L_expr->L4 (prior download/quality)    122   94     28          BACH2; CTCF; E2F2; ELK3; ELK4 candidate TFs are a subset of expression pass   PASS
# 27      BM transition           L4->L5 (matrix/alignment)     94   93      1                                   NFIA                          L5 is a subset of L4   PASS
# 28      BM transition          L5->L6 (branch + terminal)     93   92      1                                   E2F6                          L6 is a subset of L5   PASS
# 29      BM      audit         loose->L6 (terminal filter)     92   92      0                                         terminal filtering removes no TF, only edges   PASS
# 30      BM      audit             L5->branch (pseudotime)     93   92      1                                   E2F6      TFs lost at pseudotime branch assignment   PASS
# 31      BM      audit                 final TFs \\ JASPAR     92   92      0                                                   no final TF lacks a JASPAR2024 PWM   PASS
# 32      BM      audit             final TFs \\ ChIP-Atlas     92   92      0                                        no final TF is outside the ChIP-Atlas catalog   PASS
# 33      BM      audit                     final TFs \\ L1     92   92      0                                           no final TF is absent from the input genes   PASS

# Columns: kind distinguishes a funnel transition from an audit; n_from/n_to/
# n_lost quantify the step; lost_examples lists the first five dropped TFs; and
# status reports the set-identity result. Every one of the 33 rows is PASS,
# which gives the results:
#   (1) nesting holds at every layer (no TF is added back downstream);
#   (2) the three "final TFs \\ catalog" audits have n_lost = 0, so a TF outside
#       JASPAR or ChIP-Atlas never enters the network and nothing is imputed;
#   (3) the terminal filter has n_lost = 0 (it removes edges, not TFs).


## ----------------------------- Audit messages -------------------------------
# Soft cross-check against the benchmark log. No message means every checked
# layer (L1/L4/L5/L6) matched; a mismatch would be printed per dataset and would
# signal a changed object or parameter setting.
n_mismatch <- sum(vapply(funnel, function(x) length(x$audit), integer(1)))
if (n_mismatch == 0) {
  message("All benchmark cross-checks matched (L1/L4/L5/L6).")
} else {
  for (ds in names(funnel)) {
    for (msg in funnel[[ds]]$audit) message("[", ds, "] ", msg)
  }
}
# Fail the run loudly if any set-identity check did not hold.
if (any(mech_df$status == "FAIL")) {
  warning("One or more set-identity checks FAILED; inspect E2a_mechanism_check.csv.")
}
# All benchmark cross-checks matched (L1/L4/L5/L6).


## ------------------------------- Write outputs ------------------------------
write.csv(summary_df, file.path(OUT_DIR, "E2a_summary.csv"), row.names = FALSE)
write.csv(long_df, file.path(OUT_DIR, "E2a_funnel_long.csv"), row.names = FALSE)
write.csv(wide_df, file.path(OUT_DIR, "E2a_funnel_wide.csv"), row.names = FALSE)
write.csv(mech_df, file.path(OUT_DIR, "E2a_mechanism_check.csv"), row.names = FALSE)
saveRDS(funnel, file.path(OUT_DIR, "E2a_funnel_sets.rds"))

message("\nSummary:")
print(summary_df, row.names = FALSE)

message("\nFunnel:")
print(wide_df, row.names = FALSE)

message("\nMechanism checks:")
print(mech_df[, c("dataset", "transition", "n_from", "n_to", "n_lost", "status")],
      row.names = FALSE)


## --------------------------------- Figure 1 ---------------------------------
# Counts per layer. Bars are colored by dataset and labeled with the count.
# The y-axis is broken between 450 and 9,900 so the tiny L2-L6 bars are visible
# instead of being flattened against the L1 = 10,000 bars; 10,000 is shown as
# 10^4 in the upper panel.
long_df$layer <- factor(long_df$layer, levels = names(LAYER_STAGE))
ds_colors <- c(A549 = "#66C2A5", GM12878 = "#FC8D62", BM = "#8DA0CB")

p_counts <- ggplot(long_df, aes(x = layer, y = n_TFs, fill = dataset)) +
  geom_col(position = position_dodge(width = 0.82), width = 0.72) +
  geom_text(aes(label = ifelse(n_TFs == 10000, "10^4", as.character(n_TFs))),
            position = position_dodge(width = 0.82), vjust = -0.3, size = 3,
            parse = TRUE) +
  scale_fill_manual(values = ds_colors) +
  scale_y_continuous(
    breaks = c(0, 100, 200, 300, 400, 10000),
    limits = c(0, 10050)
  ) +
  ggbreak::scale_y_break(
    breaks = c(450, 9900),
    space = 0.1,
    expand = FALSE
  ) +
  labs(x = NULL, y = "Number of TFs", fill = "Dataset",
       title = "TF coverage funnel") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
print(p_counts)
ggsave(
  filename = file.path(OUT_DIR, "E2a_TF_coverage_counts.png"),
  plot     = p_counts,
  width    = 8,
  height   = 4,
  units    = "in",
  dpi      = 300
)


## --------------------------------- Figure 2 ---------------------------------
# Cumulative retention vs L1. Lines track each dataset along the funnel. The
# y-axis is broken between 4% and 98% because all candidate TFs sit in the
# 0-4% region; the lower panel shows the steady attrition from L2 onward.
p_ret <- ggplot(long_df, aes(x = layer, y = retained_vs_L1_pct,
                             color = dataset, group = dataset)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.8) +
  scale_color_manual(values = ds_colors) +
  scale_y_continuous(
    breaks = c(0, 1, 2, 3, 4, 100),
    limits = c(0, 100.5),
    labels = function(x) paste0(x, "%")
  ) +
  scale_y_break(
    breaks = c(4, 98),
    space = 0.1,
    expand = FALSE
  ) +
  labs(x = NULL, y = "Retained vs L1 (%)", color = "Dataset",
       title = "Cumulative retention") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
print(p_ret)
ggsave(
  filename = file.path(OUT_DIR, "E2a_cumulative_retention.png"),
  plot     = p_ret,
  width    = 4,
  height   = 4,
  units    = "in",
  dpi      = 300
)

message("\nTwo figures saved at 300 DPI:")
message("  - ", file.path(OUT_DIR, "E2a_TF_coverage_counts.png"))
message("  - ", file.path(OUT_DIR, "E2a_cumulative_retention.png"))

message("\nE2a finished. Tables and figures saved under ", OUT_DIR)

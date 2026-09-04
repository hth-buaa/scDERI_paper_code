###############################################################################
#                                                                             #
#   Evaluate the TF-peak relationship prediction performance of each method   #
#                                                                             #
###############################################################################

# Evaluation steps:
# 1. For each method, extract predicted TF-peak relationship
# 2. Compare with ground truth TF-peak relationship from ChIP-seq
# 3. Calculate precision, recall, and F1-score for each method


##### 1. Environment Setup and Parameter Configuration #####
# System and environment configuration
# Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
# rm(list = ls())
options(scipen = 999) # Disable scientific notation

# Load required packages
library(dplyr)
library(Pando)
library(openxlsx)

# Set base paths
base.dir = "/mnt/i/scDERI_Paper/3-benchmark"
parent.dir = "/mnt/i/scDERI_Paper/3-benchmark/9-evaluation/2-TF-peak_relationship/"
eval.dir = "/mnt/i/scDERI_Paper/2-obtain_evaluation_criteria/1-TF-evaluation/"

# Define dataset names
data.ls = c(
  "10X_hg38_Human_bone_marrow_with_T_cell_depleted",
  "SNAREseq2_hg38_GM12878",
  "SNAREseq2_hg38_A549"
)
message("There are ", length(data.ls), " datasets.")


##### 2. Extract Predicted TF-Region Relationships from Each Method #####
# Goal: Extract TF-peak relationships predicted by each GRN inference method
# Each method's results are stored in different formats and locations

# 2.1 scDERI TF-region extraction
message("Extracting TF-region relationships from scDERI...")
work.dir = paste0(base.dir, "/8-scDERI/"); setwd(work.dir)
scDERI_TF_region_links = parallel::mclapply(seq_along(data.ls), mc.cores = parallel::detectCores() - 1, function(i) {
  message("----> Conversion for the ", i, "-th dataset: ", data.ls[i])
  if (!file.exists(paste0(work.dir, data.ls[i], "/interest_cell_type_peGRN.qsave"))) {
    message("----> No eRegulon file on dataset: ", data.ls[i])
    return(NULL)
  }
  res = qs::qread(paste0(work.dir, data.ls[i], "/interest_cell_type_peGRN.qsave"))[[1]][["loose_filtered_simplified"]][, c("TF", "peak")] %>% unique
  colnames(res) = c("TF", "Region")
  return(res)
})
names(scDERI_TF_region_links) = data.ls
scDERI_TF_region_links = scDERI_TF_region_links[!sapply(scDERI_TF_region_links, is.null)]
qs::qsave(scDERI_TF_region_links, paste0(parent.dir, "scDERI_TF_region_links.qsave"))
# scDERI_TF_region_links = qs::qread(paste0(parent.dir, "scDERI_TF_region_links.qsave"))

# 2.2 STREAM TF-region extraction
message("Extracting TF-region relationships from STREAM...")
work.dir = paste0(base.dir, "/1-STREAM/"); setwd(work.dir)
STREAM_TF_region_links = parallel::mclapply(seq_along(data.ls), mc.cores = parallel::detectCores() - 1, function(i) {
  message("----> Conversion for the ", i, "-th dataset: ", data.ls[i])
  if (!file.exists(paste0(work.dir, data.ls[i], "/eRegulons.qsave"))) {
    message("----> No eRegulon file on dataset: ", data.ls[i])
    return(NULL)
  }
  res = qs::qread(paste0(work.dir, data.ls[i], "/eRegulons.qsave")) %>% lapply(., function(x) {
    data.frame(TF = rep(x$TF, length(x$peaks)), Region = x$peaks)
  }) %>% Reduce("rbind", .)
  return(res[!duplicated(res), , drop = FALSE])
})
names(STREAM_TF_region_links) = data.ls
STREAM_TF_region_links = STREAM_TF_region_links[!sapply(STREAM_TF_region_links, is.null)]
qs::qsave(STREAM_TF_region_links, paste0(parent.dir, "STREAM_TF_region_links.qsave"))
# STREAM_TF_region_links = qs::qread(paste0(parent.dir, "STREAM_TF_region_links.qsave"))

# 2.3 Pando TF-region extraction
message("Extracting TF-region relationships from Pando...")
work.dir = paste0(base.dir, "/6-Pando/"); setwd(work.dir)
Pando_TF_region_links = parallel::mclapply(seq_along(data.ls), mc.cores = parallel::detectCores() - 1, function(i) {
  message("----> Conversion for the ", i, "-th dataset: ", data.ls[i])
  if (!file.exists(paste0(work.dir, data.ls[i], ".qs"))) {
    message("----> No eRegulon file on dataset: ", data.ls[i])
    return(NULL)
  }
  grn = qs::qread(paste0(work.dir, data.ls[i], ".qs")) %>% NetworkModules
  dt = grn@meta[, c("tf", "regions")]
  colnames(dt) = c("TF", "Region")
  return(dt[!duplicated(dt),, drop = FALSE])
  message("----> Wrote eRegulons to file: ", "eRegulons_", data.ls[i], ".qsave")
})
names(Pando_TF_region_links) = data.ls
Pando_TF_region_links = Pando_TF_region_links[!sapply(Pando_TF_region_links, is.null)]
qs::qsave(Pando_TF_region_links, paste0(parent.dir, "Pando_TF_region_links.qsave"))
# Pando_TF_region_links = qs::qread(paste0(parent.dir, "Pando_TF_region_links.qsave"))

# 2.4 GLUE TF-region extraction
message("Extracting TF-region relationships from GLUE...")
work.dir = paste0(base.dir, "/4-scGLUE/"); setwd(work.dir)
scGLUE_TF_region_links = parallel::mclapply(seq_along(data.ls), mc.cores = parallel::detectCores() - 1, function(i) {
  message("----> Conversion for the ", i, "-th dataset: ", data.ls[i])
  if (!file.exists(paste0(work.dir, data.ls[i], "/draft_grn.csv"))) {
    message("----> No eRegulon file on dataset: ", data.ls[i])
    return(NULL)
  }
  grn = read.csv(paste0(work.dir, data.ls[i], "/draft_grn.csv"), sep = ",")
  colnames(grn) = c("TF", "Gene", "importance")
  gene2peak = read.csv(paste0(work.dir, data.ls[i], "/gene2peak.csv"), sep = ",")
  colnames(gene2peak) = c("Gene", "Region", "Weight")
  dt = merge(grn, gene2peak, by = "Gene")
  dt = dt[, c("TF", "Region")]
  return(dt[!duplicated(dt),, drop = FALSE])
})
names(scGLUE_TF_region_links) = data.ls
scGLUE_TF_region_links = scGLUE_TF_region_links[!sapply(scGLUE_TF_region_links, is.null)]
qs::qsave(scGLUE_TF_region_links, paste0(parent.dir, "scGLUE_TF_region_links.qsave"))
# scGLUE_TF_region_links = qs::qread(paste0(parent.dir, "scGLUE_TF_region_links.qsave"))

# 2.5 Compile all method results
method.list = list(
  scDERI = scDERI_TF_region_links,
  STREAM = STREAM_TF_region_links,
  scGLUE = scGLUE_TF_region_links,
  Pando = Pando_TF_region_links
)


##### 3. Load Ground Truth TF-Region Relationships from ChIP-seq #####
message("Loading ground truth TF-region relationships from ChIP-seq...")

# Dataset correspondence table
chipseq.hic = tibble::tibble(
  Dataset = c("10X_hg38_Human_bone_marrow_with_T_cell_depleted", "SNAREseq2_hg38_A549", "SNAREseq2_hg38_GM12878"),
  `ChIP.seq` = c("Bone-marrow-hg38-ChIP-seq", "A549-hg38-ChIP-seq", "GM12878-hg38-ChIP-seq"),
  `Hi.C` = c("Bone-marrow-hg38-HiC", "A549-hg38-HiC", "GM12878-hg38-HiC")
)

# Extract TF-peak relationships from ChIP-seq BED files
message("Extracting TF names from ChIP-seq files...")
top.peaks = 10000
ChIP_seq_TF_region_links = pbmcapply::pbmclapply(1:nrow(chipseq.hic), mc.cores = 1, function(i) {
  message("----> Processing ChIP-seq dataset ", i, ": ", chipseq.hic$Dataset[i])
  dir.ll = chipseq.hic$ChIP.seq[i]

  message("Reading BED files in directory: ", dir.ll)
  bed.files = list.files(paste0(eval.dir, dir.ll), pattern = ".bed")

  ChIP_seq_TF_region_links_each = parallel::mclapply(bed.files, mc.cores = parallel::detectCores() - 1, function(bed_file) {
    # Read BED file and extract top peaks by signal value (V7)
    bed.dt = as.data.frame(data.table::fread(paste0(eval.dir, dir.ll, "/", bed_file))) %>%
      arrange(desc(V7)) %>% dplyr::select(V1, V2, V3) %>% head(n = top.peaks)

    # Extract TF name from filename
    if (grepl("hg", dir.ll)) {TF = toupper(strsplit(bed_file, split = "_")[[1]][1])}
    data.frame(TF = rep(TF, nrow(bed.dt)), region = paste0(bed.dt$V1, "-", bed.dt$V2, "-", bed.dt$V3))
  }) %>% Reduce("rbind", .)

  return(ChIP_seq_TF_region_links_each)
})
names(ChIP_seq_TF_region_links) = chipseq.hic$ChIP.seq
qs::qsave(ChIP_seq_TF_region_links, paste0(parent.dir, "ChIP_seq_TF_region_links.qsave"))
# ChIP_seq_TF_region_links = qs::qread(paste0(parent.dir, "ChIP_seq_TF_region_links.qsave"))


##### 4. Define Evaluation Functions #####
# 4.1 Calculate region overlap metrics
cal_overlap_regions = function(
    r1, # Predicted regions in format "chr1:100-200"
    r2  # Ground truth regions in format "chr1:100-200"
) {
  if (is.null(r1) | is.null(r2)) {return(c(precision = 0, recall = 0, f.score = 0))}

  # Convert strings to GRanges and find overlaps
  hits = GenomicRanges::findOverlaps(Signac::StringToGRanges(r1, sep = c(":", "-")),
                                     Signac::StringToGRanges(r2, sep = c(":", "-")))

  # Precision: proportion of predicted regions overlapping ground truth
  prec = length(unique(S4Vectors::queryHits(hits))) / length(r1)

  # Recall: proportion of ground truth regions covered by predictions
  recall = length(unique(S4Vectors::subjectHits(hits))) / length(r2)

  # f1-score: harmonic mean of precision and recall
  if (prec + recall == 0) {f.score = 0} else {f.score = 2 * prec * recall / (prec + recall)}

  return(c(precision = prec, recall = recall, f.score = f.score))
}

# 4.2 Calculate precision, recall, and F-score for TF-region predictions
cal_prec_recall_fscore = function(
    x, # Predicted TF-region relationships
    y, # Ground truth TF-region relationships
    mc.cores = NULL # Number of cores for parallel computation
) {
  if (is.null(mc.cores)) {mc.cores = max(1, parallel::detectCores() - 1)}

  # Identify common TFs between predictions and ground truth
  tfs = intersect(unique(x$TF), unique(y$TF))
  message("--------> ", length(tfs), " common TFs to calculate TF-to-region linkages")
  if (length(tfs) < 1) {return(data.frame(precision = 0, recall = 0, f.score = 0, TF = 0))}

  # Extract regions for each TF from predictions
  r1_list = parallel::mclapply(tfs, function(z) {
    regionsss = unique(x[x$TF == z,, drop = FALSE]$Region)
    if (regionsss[1] == "") {return(NULL)} else {
      unique(x[x$TF == z,, drop = FALSE]$Region) %>% strsplit(., split = ";") %>% sapply(., "[[", 1)
    }
  }, mc.cores = mc.cores)
  names(r1_list) = tfs

  # Extract regions for each TF from ground truth
  r2_list = parallel::mclapply(tfs, function(z) {
    unique(y[y$TF == z,, drop = FALSE]$region)
  }, mc.cores = mc.cores)
  names(r2_list) = tfs

  # Calculate scores for each TF
  tf.scores = parallel::mclapply(tfs, function(z) {
    cal_overlap_regions(r1 = r1_list[[z]], r2 = r2_list[[z]])
  }, mc.cores = mc.cores) %>% Reduce("rbind", .)

  # Format results
  if (is.vector(tf.scores)) {tf.scores = t(as.data.frame(tf.scores))}
  rownames(tf.scores) = tfs
  tf.scores = as.data.frame(tf.scores)
  tf.scores$TF = tfs
  return(tf.scores)
}


##### 5. Evaluate All Methods #####
# 5.1 Calculate evaluation metrics for all methods and datasets
TF_region_benchmark_all = pbmcapply::pbmclapply(data.ls, mc.cores = 1, function(chip.data) {
  message("----> ChIP-seq data: ", chip.data)
  TF_region_benchmark_all_each = lapply(seq_along(method.list), function(i) {
    message("----> Method: ", names(method.list)[i])
    method = method.list[[i]]
    if (!chip.data %in% names(method)) {
      return(c(prec = 0, recall = 0, f.score = 0, TF = 0, method = names(method.list)[i]))
    }
    method.scores = cal_prec_recall_fscore(
      x = method[[chip.data]],
      y = ChIP_seq_TF_region_links[[chipseq.hic$ChIP.seq[chipseq.hic$Dataset == chip.data]]],
      mc.cores = parallel::detectCores() - 1
    )
    method.scores[is.na(method.scores)] = 0.0
    cbind(method.scores, method = rep(names(method.list)[i], nrow(method.scores)))
  }) %>% Reduce("rbind", .)
  TF_region_benchmark_all_each$method = factor(
    TF_region_benchmark_all_each$method,
    levels = c("scDERI", "STREAM", "scGLUE", "Pando")
  )
  return(TF_region_benchmark_all_each)
})
names(TF_region_benchmark_all) = data.ls
for (i in seq_along(TF_region_benchmark_all)) {
  TF_region_benchmark_all[[i]][1:length(TF_region_benchmark_all)] = lapply(TF_region_benchmark_all[[i]][1:length(TF_region_benchmark_all)], as.numeric)
}
qs::qsave(TF_region_benchmark_all, paste0(parent.dir, "TF_region_benchmark_all.qsave"))
# TF_region_benchmark_all = qs::qread(paste0(parent.dir, "TF_region_benchmark_all.qsave"))

wb = createWorkbook()
index_df = data.frame(OriginalName = names(TF_region_benchmark_all), SheetName = paste0("Sheet", 1:length(TF_region_benchmark_all)))
addWorksheet(wb, "Index")
writeData(wb, sheet = "Index", x = index_df)
for (i in 1:length(TF_region_benchmark_all)) {
  sheet_name = index_df$SheetName[i]
  df = TF_region_benchmark_all[[i]]
  df_with_rownames = cbind(rownames(df), df)
  colnames(df_with_rownames)[1] = "RowNames"
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet = sheet_name, x = df_with_rownames, rowNames = FALSE, colNames = TRUE)
}
saveWorkbook(wb, paste0(parent.dir, "TF_region_benchmark_all.xlsx"), overwrite = TRUE)
rm(wb, index_df, i, sheet_name, df, df_with_rownames)

# 5.2 Calculate average metrics per method per dataset
TF_region_benchmark = parallel::mclapply(seq_along(TF_region_benchmark_all), function(i) {
  message("----> Summarizing results for dataset: ", names(TF_region_benchmark_all)[i])
  x = TF_region_benchmark_all[[i]]
  as.data.frame(do.call(rbind, lapply(split(x, f = x$method), function(xx) {c(
    precision = mean(as.numeric(xx$precision)),
    recall = mean(as.numeric(xx$recall)),
    f.score = mean(as.numeric(xx$f.score))
  )})))
}, mc.cores = 1)
names(TF_region_benchmark) = names(TF_region_benchmark_all)

# 5.3 Save and display results
print(TF_region_benchmark)
qs::qsave(TF_region_benchmark, paste0(parent.dir, "TF_region_benchmark.qsave"))
# save.image("2_TF-peak_relationship_benchmark.RData")
# TF_region_benchmark = qs::qread(paste0(parent.dir, "TF_region_benchmark.qsave"))

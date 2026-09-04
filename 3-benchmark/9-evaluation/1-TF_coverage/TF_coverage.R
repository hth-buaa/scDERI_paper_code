######################################################################
#                                                                    #
#   Evaluate the TF coverage prediction performance of each method   #
#                                                                    #
######################################################################

# Evaluation steps:
# 1. For each method, extract predicted TF lists
# 2. Compare with ground truth TFs from ChIP-seq
# 3. Calculate precision, recall, and F1-score for each method


##### 1. Environment Setup and Package Loading #####
# System and environment configuration
# Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
# rm(list = ls())
options(scipen = 999) # Disable scientific notation

# Load required packages
library(dplyr)
library(Pando)


##### 2. Path and Data Configuration #####
base.dir = "/mnt/i/scDERI_Paper/3-benchmark"
parent.dir = "/mnt/i/scDERI_Paper/3-benchmark/9-evaluation/1-TF_coverage/"
chipseq.dir = "/mnt/i/scDERI_Paper/2-obtain_evaluation_criteria/1-TF-evaluation/"

# Define dataset names
data.ls = c(
  "10X_hg38_Human_bone_marrow_with_T_cell_depleted",
  "SNAREseq2_hg38_GM12878",
  "SNAREseq2_hg38_A549"
)
message("There are ", length(data.ls), " datasets")


##### 3. Extract Predicted TFs from Each Method #####
# Goal: Extract TF lists predicted by each GRN inference method
# Each method's results are stored in different formats and locations

# 3.1 scDERI TFs
message("Extracting TFs from scDERI...")
work.dir = paste0(base.dir, "/8-scDERI/"); setwd(work.dir)
scDERI_TFs = pbapply::pblapply(seq_along(data.ls), function(i) {
  message("----> Conversion for the ", i, "-th dataset: ", data.ls[i])
  if (!file.exists(paste0(work.dir, data.ls[i], "/interest_cell_type_peGRN.qsave"))) {
    message("----> No eRegulon file on dataset: ", data.ls[i])
    return(NULL)
  }
  return(qs::qread(paste0(work.dir, data.ls[i], "/interest_cell_type_peGRN.qsave"))[[1]][["loose_filtered_simplified"]][["TF"]] %>% unique)
})
names(scDERI_TFs) = data.ls
scDERI_TFs = scDERI_TFs[!sapply(scDERI_TFs, is.null)]
qs::qsave(scDERI_TFs, paste0(parent.dir, "scDERI_TFs.qsave"))
# scDERI_TFs = qs::qread(paste0(parent.dir, "scDERI_TFs.qsave"))


# 3.2 STREAM TFs
message("Extracting TFs from STREAM...")
work.dir = paste0(base.dir, "/1-STREAM/"); setwd(work.dir)
STREAM_TFs = pbapply::pblapply(seq_along(data.ls), function(i) {
  message("----> Conversion for the ", i, "-th dataset: ", data.ls[i])
  if (!file.exists(paste0(work.dir, data.ls[i], "/eRegulons.qsave"))) {
    message("----> No eRegulon file on dataset: ", data.ls[i])
    return(NULL)
  }
  return(qs::qread(paste0(work.dir, data.ls[i], "/eRegulons.qsave")) %>% lapply(., "[[", "TF") %>% unlist %>% unique)
})
names(STREAM_TFs) = data.ls
STREAM_TFs = STREAM_TFs[!sapply(STREAM_TFs, is.null)]
qs::qsave(STREAM_TFs, paste0(parent.dir, "STREAM_TFs.qsave"))
# STREAM_TFs = qs::qread(paste0(parent.dir, "STREAM_TFs.qsave"))

# 3.3 IReNA TFs
message("Extracting TFs from IReNA...")
work.dir = paste0(base.dir, "/2-IReNA/"); setwd(work.dir)
IReNA_TFs = pbapply::pblapply(seq_along(data.ls), function(i) {
  message("----> Conversion for the ", i, "-th dataset: ", data.ls[i])
  if (!file.exists(paste0(work.dir, data.ls[i], "/IReNA_result.qsave"))) {
    message("----> No eRegulon file on dataset: ", data.ls[i])
    return(NULL)
  }
  return(qs::qread(paste0(work.dir, data.ls[i], "/IReNA_result.qsave"))[["FOSF_RegMTF_Cor_EnTFsTarg"]]$TF %>% as.character %>% unique)
})
names(IReNA_TFs) = data.ls
IReNA_TFs = IReNA_TFs[!sapply(IReNA_TFs, is.null)]
qs::qsave(IReNA_TFs, paste0(parent.dir, "IReNA_TFs.qsave"))
# IReNA_TFs = qs::qread(paste0(parent.dir, "IReNA_TFs.qsave"))

# 3.4 FigR TFs
message("Extracting TFs from FigR...")
work.dir = paste0(base.dir, "/3-FigR/"); setwd(work.dir)
FigR_TFs = pbapply::pblapply(seq_along(data.ls), function(i) {
  message("----> Conversion for the ", i, "-th dataset: ", data.ls[i])
  if (!file.exists(paste0(work.dir, data.ls[i], "/FigR_result.qsave"))) {
    message("----> No eRegulon file on dataset: ", data.ls[i])
    return(NULL)
  }
  return(qs::qread(paste0(work.dir, data.ls[i], "/FigR_result.qsave"))$Motif %>% unique)
})
names(FigR_TFs) = data.ls
FigR_TFs = FigR_TFs[!sapply(FigR_TFs, is.null)]
qs::qsave(FigR_TFs, paste0(parent.dir, "FigR_TFs.qsave"))
# FigR_TFs = qs::qread(paste0(parent.dir, "FigR_TFs.qsave"))

# 3.5 GLUE TFs
message("Extracting TFs from GLUE...")
work.dir = paste0(base.dir, "/4-scGLUE/"); setwd(work.dir)
scGLUE_TFs = pbapply::pblapply(seq_along(data.ls), function(i) {
  message("----> Conversion for the ", i, "-th dataset: ", data.ls[i])
  if (!file.exists(paste0(work.dir, data.ls[i], "/draft_grn.csv"))) {
    message("----> No eRegulon file on dataset: ", data.ls[i])
    return(NULL)
  }
  grn = read.csv(paste0(work.dir, data.ls[i], "/draft_grn.csv"), sep = ",")
  gene2peak = read.csv(paste0(work.dir, data.ls[i], "/gene2peak.csv"), sep = ",")
  grn = grn[grn$target %in% unique(gene2peak$Source), , drop = FALSE]
  dt = grn
  return(unique(split(dt, f = dt$TF) %>% sapply(., nrow) %>% sort(., decreasing = TRUE) %>% names))
})
names(scGLUE_TFs) = data.ls
scGLUE_TFs = scGLUE_TFs[!sapply(scGLUE_TFs, is.null)]
qs::qsave(scGLUE_TFs, paste0(parent.dir, "scGLUE_TFs.qsave"))
# scGLUE_TFs = qs::qread(paste0(parent.dir, "scGLUE_TFs.qsave"))

# 3.6 DIRECT-NET TFs
message("Extracting TFs from DIRECT-NET...")
work.dir = paste0(base.dir, "/5-DIRECT-NET/"); setwd(work.dir)
DIRECT_NET_TFs = pbapply::pblapply(seq_along(data.ls), function(i) {
  message("----> Conversion for the ", i, "-th dataset: ", data.ls[i])
  if (!file.exists(paste0(work.dir, data.ls[i], "/", data.ls[i], "_directnet_out.qs"))) {
    message("----> No eRegulon file on dataset: ", data.ls[i])
    return(NULL)
  }
  df = qs::qread(paste0(work.dir, data.ls[i], "/", data.ls[i], "_directnet_out.qs"))
  dt = df$link
  return(unique(split(dt, f = dt$TF) %>% sapply(., nrow) %>% sort(., decreasing = TRUE) %>% names))
})
names(DIRECT_NET_TFs) = data.ls
DIRECT_NET_TFs = DIRECT_NET_TFs[!sapply(DIRECT_NET_TFs, is.null)]
qs::qsave(DIRECT_NET_TFs, paste0(parent.dir, "DIRECT_NET_TFs.qsave"))
# DIRECT_NET_TFs = qs::qread(paste0(parent.dir, "DIRECT_NET_TFs.qsave"))

# 3.7 Pando TFs
message("Extracting TFs from Pando...")
work.dir = paste0(base.dir, "/6-Pando/"); setwd(work.dir)
Pando_TFs = pbapply::pblapply(seq_along(data.ls), function(i) {
  message("----> Conversion for the ", i, "-th dataset: ", data.ls[i])
  if (!file.exists(paste0(work.dir, data.ls[i], ".qs"))) {
    message("----> No eRegulon file on dataset: ", data.ls[i])
    return(NULL)
  }
  grn = qs::qread(paste0(work.dir, data.ls[i], ".qs")) %>% NetworkModules
  dt = grn@meta
  return(unique(split(dt, f = dt$tf) %>% sapply(., nrow) %>% sort(., decreasing = TRUE) %>% names))
  message("----> Wrote eRegulons to file: ", "eRegulons_", data.ls[i], ".qsave")
})
names(Pando_TFs) = data.ls
Pando_TFs = Pando_TFs[!sapply(Pando_TFs, is.null)]
qs::qsave(Pando_TFs, paste0(parent.dir, "Pando_TFs.qsave"))
# Pando_TFs = qs::qread(paste0(parent.dir, "Pando_TFs.qsave"))

# 3.8 scMEGA TFs
message("Extracting TFs from scMEGA...")
work.dir = paste0(base.dir, "/7-scMEGA/"); setwd(work.dir)
scMEGA_TFs = pbapply::pblapply(seq_along(data.ls), function(i) {
  message("----> Conversion for the ", i, "-th dataset: ", data.ls[i])
  if (!file.exists(paste0(work.dir, data.ls[i], ".RData"))) {
    message("----> No eRegulon file on dataset: ", data.ls[i])
    return(NULL)
  }
  load(paste0(work.dir, data.ls[i], ".RData"))
  dt = tf.gene.cor
  return(unique(split(dt, f = dt$tf) %>% sapply(., nrow) %>% sort(., decreasing = TRUE) %>% names))
  message("----> Wrote eRegulons to file: ", "eRegulons_", data.ls[i], ".qsave")
})
names(scMEGA_TFs) = data.ls
scMEGA_TFs = scMEGA_TFs[!sapply(scMEGA_TFs, is.null)]
qs::qsave(scMEGA_TFs, paste0(parent.dir, "scMEGA_TFs.qsave"))
# scMEGA_TFs = qs::qread(paste0(parent.dir, "scMEGA_TFs.qsave"))


##### 4. Load Ground Truth TFs from ChIP-seq Data #####
message("Extracting TFs from ChIP-seq benchmark data...")
chipseq.hic = tibble::tibble(
  Dataset = c("10X_hg38_Human_bone_marrow_with_T_cell_depleted", "SNAREseq2_hg38_A549", "SNAREseq2_hg38_GM12878"),
  `ChIP.seq` = c("Bone-marrow-hg38-ChIP-seq", "A549-hg38-ChIP-seq", "GM12878-hg38-ChIP-seq"),
  `Hi.C` = c("Bone-marrow-hg38-HiC", "A549-hg38-HiC", "GM12878-hg38-HiC")
)
ChIP_seq_TFs = pbapply::pblapply(1:nrow(chipseq.hic), function(i) {
  message("----> Reading TF ChIP-seq for the ", i, "-th dataset: ", chipseq.hic$Dataset[i])
  if (chipseq.hic$ChIP.seq[i] == "N/A") {return(NULL)}
  tf.lst = list.files(path = paste0(chipseq.dir, chipseq.hic$ChIP.seq[i]), pattern = ".bed") %>%
    strsplit(., split = "_") %>% sapply(., "[", 1) %>% unique
  if (grepl("hg", chipseq.hic$Dataset[i])) {
    message("----> Converting TF symbols to upper case")
    tf.lst = toupper(tf.lst)
  }
  tf.lst
})
names(ChIP_seq_TFs) = chipseq.hic$Dataset
ChIP_seq_TFs = ChIP_seq_TFs[!sapply(ChIP_seq_TFs, is.null)]
ChIP_seq_TFs = ChIP_seq_TFs[sapply(ChIP_seq_TFs, length) > 0]
qs::qsave(ChIP_seq_TFs, paste0(parent.dir, "ChIP_seq_TFs.qsave"))
# ChIP_seq_TFs =  qs::qread(paste0(parent.dir, "ChIP_seq_TFs.qsave"))


##### 5. TF Coverage Evaluation #####
# 5.1 Compile all method results
method.list = list(
  scDERI = scDERI_TFs,
  STREAM = STREAM_TFs,
  scGLUE = scGLUE_TFs,
  Pando = Pando_TFs,
  scMEGA = scMEGA_TFs,
  FigR = FigR_TFs,
  DIRECT_NET = DIRECT_NET_TFs,
  IReNA = IReNA_TFs
)

# 5.2 Define evaluation function
# Calculates precision, recall, and F1-score for predicted vs. ground truth TFs
cal_prec_recall_fscore = function(
    x, # Predicted TF vector
    y # Ground truth TF vector
) {
  ltp = length(intersect(x, y))
  prec = ltp / length(x)
  recall = ltp / length(y)
  sum = prec + recall
  if (sum == 0) {
    return(c(prec = 0, recall = 0, f.score = 0))
  } else {
    return(c(prec = prec, recall = recall, f.score = 2 * prec * recall / sum))
  }
}

# 5.3 Calculate evaluations metrics for all methods and datasets
# Calculate three scores of TF coverage for all methods
TF_coverage_benchmark = pbapply::pblapply(data.ls, function(chip.data) {
  message("----> ChIP-seq data: ", chip.data)
  res = lapply(method.list, function(method) {
    if (!chip.data %in% names(method)) {return(c(prec = 0, recall = 0, f.score = 0))}
    cal_prec_recall_fscore(x = method[[chip.data]], y = ChIP_seq_TFs[[chip.data]])
  }) %>% Reduce("rbind", .)
  rownames(res) = names(method.list)
  res[is.na(res)] = 0
  res
})
names(TF_coverage_benchmark) = data.ls

# 5.4 Save and display results
print(TF_coverage_benchmark)
qs::qsave(TF_coverage_benchmark, paste0(parent.dir, "TF_coverage_benchmark.qsave"))
# save.image("1_TF_coverage_benchmark.RData")
# TF_coverage_benchmark = qs::qread(paste0(parent.dir, "TF_coverage_benchmark.qsave"))

################################################
#                                              #
#       Run FigR on benchmarking dataset       #
#                                              #
################################################


### FigR official website: https://buenrostrolab.github.io/FigR/

# Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
# rm(list = ls())

# sys = "wsl"
sys = "win"

if (sys == "win") {
  setwd("I:/scDERI_Paper/3-benchmark/3-FigR/")
} else {
  setwd("/mnt/i/scDERI_Paper/3-benchmark/3-FigR/")
}

library(FigR)
library(SummarizedExperiment)
library(SeuratObject)
library(sp)
library(IRanges)
library(GenomicRanges)
library(GenomeInfoDb)
library(S4Vectors)
library(BiocGenerics)
library(BSgenome.Hsapiens.UCSC.hg38)
library(pROC)
library(parallel)
library(Signac)
library(qs)

# Parameters
nCores = ceiling(detectCores() / 2)
if (sys == "win") {
  args = c(
    "I:/scDERI_Paper/3-benchmark/0-downsampling/Downsampled-data/10X_hg38_Human_bone_marrow_with_T_cell_depleted",
    "I:/scDERI_Paper/3-benchmark/0-downsampling/Downsampled-data/SNAREseq2_hg38_A549",
    "I:/scDERI_Paper/3-benchmark/0-downsampling/Downsampled-data/SNAREseq2_hg38_GM12878"
  )
} else {
  args = c(
    "/mnt/i/scDERI_Paper/3-benchmark/0-downsampling/Downsampled-data/10X_hg38_Human_bone_marrow_with_T_cell_depleted/",
    "/mnt/i/scDERI_Paper/3-benchmark/0-downsampling/Downsampled-data/SNAREseq2_hg38_A549/",
    "/mnt/i/scDERI_Paper/3-benchmark/0-downsampling/Downsampled-data/SNAREseq2_hg38_GM12878/"
  )
}
dataset = basename(args)

FigR_result_all = list()
for (i in dataset) {FigR_result_all[[i]] = list()}
# FigR_result_all[["10X_hg38_Human_bone_marrow_with_T_cell_depleted"]] = list()
# FigR_result_all[["SNAREseq2_hg38_A549"]] = list()
# FigR_result_all[["SNAREseq2_hg38_GM12878"]] = list()

### Run FigR
for (n in 1:length(args)) {
  # n=1
  message(" ")
  message("Use the ", n, "-th dataset.")
  message("Data is: ", basename(args)[n])

  ### Prepare data
  message("Prepare data...")
  file.path = args[n]
  file.name = strsplit(file.path, split = "\\/") %>% `[[` (1) %>% tail(n = 1)
  if (sys == "win") {
    data.dir = paste0("I:/scDERI_Paper/3-benchmark/3-FigR/", gsub(".qsave", "", file.name))
  } else {
    data.dir = paste0("/mnt/i/scDERI_Paper/3-benchmark/3-FigR/", gsub(".qsave", "", file.name))
  }
  dir.create(data.dir)
  org = strsplit(file.name, split = "_") %>% `[[` (1) %>% `[` (2)

  setwd(data.dir)
  obj = qs::qread(paste0("../../0-downsampling/Downsampled-data/", file.name, "/", file.name, ".qsave"))
  print(dim(obj))
  class(obj)

  rna_counts = GetAssayData(obj, assay = "RNA", layer = "counts")
  scRNAseq = CreateSeuratObject(counts = rna_counts, assay = "RNA")
  atac_counts = GetAssayData(obj, assay = "ATAC", layer = "counts")
  scATACseq = CreateSeuratObject(counts = atac_counts, assay = "ATAC")

  # Sample run code (Demo)
  # set.seed(123)
  # rna_feature_idx = sample(1:nrow(scRNAseq), 500)
  # scRNAseq = scRNAseq[rna_feature_idx, ]

  rnaMat = scRNAseq@assays$RNA@layers$counts
  rownames(rnaMat) = rownames(scRNAseq)
  colnames(rnaMat) = colnames(scRNAseq)

  ATAC.SE_Seurat = scATACseq
  ATAC.SE_Seurat = subset(ATAC.SE_Seurat, cells = colnames(rnaMat))
  ATAC.SE_Seurat = subset(ATAC.SE_Seurat, features = grep("chr", rownames(ATAC.SE_Seurat), value = TRUE))  # Match row names containing "chr" using regular expressions.

  counts_matrix = GetAssayData(ATAC.SE_Seurat, assay = "ATAC", layer = "counts")
  peak_ranges = StringToGRanges(rownames(counts_matrix), sep = c("-", "-"))
  col_data = ATAC.SE_Seurat@meta.data
  ATAC.SE = SummarizedExperiment(
    assays = list(counts = counts_matrix),
    rowRanges = peak_ranges,
    colData = col_data
  ); rm(counts_matrix, peak_ranges, col_data); gc()

  ### Identify DORC (Co-regulatory Open Chromatin Regions)
  message("Do runGenePeakcorr...")
  a = Sys.time(); cisCor = runGenePeakcorr(
    ATAC.se = ATAC.SE,
    RNAmat = rnaMat,
    genome = "hg38",
    nCores = nCores
  ); print(Sys.time() - a); rm(a)

  # Filter peak-gene correlations by p-value
  message("Filter peak-gene correlations by p-value...")
  cisCor.filt = cisCor %>% filter(pvalZ <= 0.05)

  # Determine DORC genes
  message("Determine DORC genes...")
  dorcGenes = cisCor.filt %>% dorcJPlot(cutoff = 2, returnGeneList = TRUE)
  print(length(dorcGenes))
  if (length(dorcGenes) == 0) {stop("No dorcGenes.")}

  # Get DORC scores
  message("Get DORC scores...")
  dorcMat = getDORCScores(ATAC.SE, dorcTab = cisCor.filt, geneList = dorcGenes, nCores = nCores)

  # Run FigR
  message("Run FigR...")
  fig.d = runFigRGRN(
    ATAC.se = ATAC.SE,
    dorcK = min(length(dorcGenes) - 1, 30),
    rnaMat = rnaMat,
    dorcMat = dorcMat,
    dorcTab = cisCor.filt,
    genome = "hg38",
    dorcGenes = dorcGenes,
    nCores = nCores # ceiling(nCores / 20)
  ); gc()

  message("Save result...")
  fig.d_sorted = fig.d[order(-fig.d$Score), ]
  fig.d_sorted = fig.d_sorted[fig.d_sorted$Score > 0.05, ]
  qs::qsave(fig.d_sorted, paste0(data.dir, "/FigR_result.qsave"))

  FigR_result_all[[n]] = fig.d_sorted
  save.image("FigR_result.RData")
}
save.image("../FigR_all_result.RData")

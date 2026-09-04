#################################################
#                                               #
#       Run IReNA on benchmarking dataset       #
#                                               #
#################################################


# Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
# rm(list = ls())

# sys = "wsl"
sys = "win"

if (sys == "win") {
  setwd("I:/scDERI_Paper/3-benchmark/2-IReNA/")
} else {
  setwd("/mnt/i/scDERI_Paper/3-benchmark/2-IReNA/")
}

library(IReNA)
library(Seurat)
library(SeuratObject)
library(Signac)
library(monocle)
library(scCustomize)
library(GenomeInfoDb)
library(S4Vectors)
library(IRanges)
library(BSgenome)
library(GenomicRanges)
library(Biostrings)
library(XVector)
library(BiocIO)
library(rtracklayer)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(igraph)
library(motifmatchr)
library(parallel)

# Parameters
nCores = ceiling(detectCores() / 3)
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

IReNA_result_all = list()
for (i in dataset) {IReNA_result_all[[i]] = list()}
# IReNA_result_all[["10X_hg38_Human_bone_marrow_with_T_cell_depleted"]] = list()
# IReNA_result_all[["SNAREseq2_hg38_A549"]] = list()
# IReNA_result_all[["SNAREseq2_hg38_GM12878"]] = list()

### Run IReNA
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
    data.dir = paste0("I:/scDERI_Paper/3-benchmark/2-IReNA/", gsub(".qsave", "", file.name))
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

  ########## Step 0: Select motif database based on species ##########
  message("Do Tranfac201803_Hs_MotifTFsF...")
  motif1 = Tranfac201803_Hs_MotifTFsF
  candidate_tfs = unlist(strsplit(motif1$TFs,';'))
  gc()

  ########## Step 1: Construct modular pseudobulk expression profile ##########
  message("Convert count matrix to sparse format, create CellDataSet object")
  interest_scRNAseq_for_GRN = scRNAseq
  interest_RNAsparseMatrix_for_GRN = as(as.matrix(JoinLayers(interest_scRNAseq_for_GRN)@assays$RNA$counts), 'sparseMatrix')
  interest_cds = newCellDataSet(
    cellData = interest_RNAsparseMatrix_for_GRN,
    phenoData = new('AnnotatedDataFrame', data=JoinLayers(interest_scRNAseq_for_GRN)@meta.data),
    featureData = new('AnnotatedDataFrame', data=data.frame(gene_short_name = row.names(interest_RNAsparseMatrix_for_GRN), row.names = row.names(interest_RNAsparseMatrix_for_GRN))),
    lowerDetectionLimit = 0, # 0.1
    expressionFamily = negbinomial.size()) # Note: Do not view this interest_cds

  message("Estimate size factors and dispersion")
  interest_cds = estimateSizeFactors(interest_cds)
  interest_cds = estimateDispersions(interest_cds)

  message("Select highly variable genes")
  disp_table = dispersionTable(interest_cds)
  disp.genes = subset(
    disp_table,
    mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit
  )$gene_id

  message("Specify genes for subsequent analysise")
  interest_cds = setOrderingFilter(interest_cds, disp.genes)

  message("Dimensionality reduction")
  interest_cds = reduceDimension(interest_cds, method = 'DDRTree', max_components = 2)

  message("Order cells by pseudotime")
  interest_cds = orderCells(interest_cds)

  scRNAseq@meta.data[, "Pseudotime"] = 0
  for (cell in rownames(scRNAseq@meta.data)) {
    scRNAseq@meta.data[cell, "Pseudotime"] = pData(interest_cds)[cell, "Pseudotime"]
  }; gc()

  scRNAseq = scCustomize::Convert_Assay(seurat_object = scRNAseq, assay = "RNA", convert_to = "V3")
  expr_matrix = scRNAseq@assays$RNA@data
  expr_matrix[is.na(expr_matrix) | is.nan(expr_matrix) | is.infinite(expr_matrix)] = 0
  scRNAseq = SetAssayData(scRNAseq, layer = "data", new.data = expr_matrix)
  scRNAseq = scRNAseq[which(apply(expr_matrix, 1, var, na.rm = TRUE) > 0), ]

  message("Do irena_step1...")
  a = Sys.time(); exp = irena_step1(scRNAseq, 50, 4); print(Sys.time() - a); rm(a); gc()

  ########## Step 2: Infer GRNs ##########
  message("Do irena_step2...")
  a = Sys.time(); regulatory_relationships = irena_step2(as.matrix(exp), candidate_tfs = candidate_tfs); print(Sys.time() - a); rm(a); gc()

  message("Do regulatory_relationships")
  a = Sys.time(); regulatory_relationships = regulatory_relationships[abs(regulatory_relationships$Correlation)>0.8, ]; print(Sys.time() - a); rm(a); gc()

  message("Select gene tss for data...")
  gene_tss = hg38_gene_tss

  message("Predicting potential binding TFs for TGs using the R package motifmatchr...")
  gene_tss = gene_tss[gene_tss$symbol %in% rownames(exp), ]
  gene_tss_GR = GRanges(paste0(gene_tss[, 1], ':', gene_tss[,2 ], '-', gene_tss[, 3]))

  message("Do identify_region_tfs..."); gc()
  a = Sys.time(); fimo_regulation = identify_region_tfs(
    gene_tss_GR, gene_tss[, 6], motif1, BSdb = BSgenome.Hsapiens.UCSC.hg38
  ); print(Sys.time() - a); rm(a); gc()

  message("Do filter_regulation_fimo...")
  a = Sys.time(); regulatory_relationships = filter_regulation_fimo(fimo_regulation, regulatory_relationships); print(Sys.time() - a); rm(a)

  ########## Step 3: Decode gene regulatory networks ##########
  message("Do network_analysis...")
  IReNA_result = network_analysis(
    regulatory_relationships,
    exp[,1:2],
    TFFDR1 = 10,
    TFFDR2 = 10,
    ModuleFDR = 0.01
  )

  # View TFs regulating other modules
  IReNA_result$TF_module_regulation[1:3,]

  # View enriched TF network
  IReNA_result$TF_network[1:3,]

  # View simplified network
  IReNA_result$intramodular_network[1:3,]

  IReNA_result_all[[n]] = IReNA_result
  qs::qsave(IReNA_result, paste0(data.dir, "/IReNA_result.qsave"))
  save.image("IReNA_result.RData")
}
save.image("../IReNA_result.RData")

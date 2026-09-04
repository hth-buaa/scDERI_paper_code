###########################################
#                                         #
#    Run Pando on benchmarking dataset    #
#                                         #
###########################################


# ========== 1. Load Required Packages ==========
library(Seurat)
library(qs)
library(stringr)
library(Pando)
library(tidyverse)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Mmusculus.UCSC.mm10)
library(JASPAR2020)


# ========== 2. Set Working Directory ==========
setwd("/mnt/i/scDERI_Paper/3-benchmark/6-Pando/")


# ========== 3. Define Input Datasets ==========
# Note: Dataset A549 and GM12878 produce errors; only dataset BM runs successfully
list = c(
  "/mnt/i/scDERI_Paper/3-benchmark/0-downsampling/Downsampled-data/10X_hg38_Human_bone_marrow_with_T_cell_depleted/10X_hg38_Human_bone_marrow_with_T_cell_depleted",
  "/mnt/i/scDERI_Paper/3-benchmark/0-downsampling/Downsampled-data/SNAREseq2_hg38_A549/SNAREseq2_hg38_A549",
  "/mnt/i/scDERI_Paper/3-benchmark/0-downsampling/Downsampled-data/SNAREseq2_hg38_GM12878/SNAREseq2_hg38_GM12878"
)


# ========== 4. Load Motif Data ==========
# Load built-in motif database from Pando package
data('motifs', package = "Pando")


# ========== 5. Process Dataset (Working Dataset) ==========
# Note: Dataset 1 and 3 produce errors; only processing dataset 2
# name = list[1]
name = list[2]
# name = list[3]

# 5.1 Load Seurat object
qs = qread(file = paste0(name,'.qsave'))

# 5.2 Convert Assay5 objects to Assay objects if needed
if (grepl("Assay5", class(qs[["RNA"]]))) {qs[["RNA"]] = as(qs[["RNA"]], "Assay")}
if (grepl("Assay5", class(qs[["ATAC"]]))) {qs[["ATAC"]] = as(qs[["ATAC"]], "Assay")}

# 5.3 Select variable features
seurat_object = Seurat::FindVariableFeatures(qs, assay = 'RNA')

# 5.4 Extract gene annotations from EnsDb
#annotations = suppressWarnings(GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86, verbose = FALSE))
#seqlevelsStyle(annotations) = 'UCSC'
#Annotation(seurat_object@assays[["ATAC"]]) = annotations


# ========== 6. Pando GRN Analysis Pipeline ==========
# 6.1 Initialize GRN object and select candidate regions
seurat_object = initiate_grn(seurat_object, rna_assay = 'RNA', peak_assay = 'ATAC')

# 6.2 Scan candidate regions for TF binding motifs
# Scan candidate regions for TF binding motifs
# pfm = TFBSTools::getMatrixSet(
#   x = JASPAR2020,
#   opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE, species = "Homo sapiens")
# )
seurat_object = find_motifs(
  seurat_object,
  pfm = motifs,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

# 6.3 Infer gene regulatory network
seurat_object = infer_grn(seurat_object)

# 6.4 Find gene and regulatory modules
seurat_object = find_modules(seurat_object)

# 6.5 Extract information
modules = NetworkModules(seurat_object)

# 6.7 Visualize GRN using correlation-based UMAP
seurat_object = get_network_graph(seurat_object, umap_method = "corr")


# ========== 7. Save Results ==========
qsave(seurat_object, file = paste0(basename(name), '.qs'))

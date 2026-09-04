###########################################
#                                         #
#    Run scMEGA on benchmarking dataset   #
#                                         #
###########################################


# ========== 1. Load Required Packages ==========
library(Seurat)
library(qs)
library(stringr)
library(tidyverse)
library(stats)
suppressMessages(library(ArchR))
suppressMessages(library(Signac))
suppressMessages(library(scMEGA))
suppressMessages(library(Nebulosa))
suppressMessages(library(GenomeInfoDb))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg19))
suppressMessages(library(BSgenome.Mmusculus.UCSC.mm10))
suppressMessages(library(EnsDb.Hsapiens.v86))
suppressMessages(library(EnsDb.Hsapiens.v75))
suppressMessages(library(EnsDb.Mmusculus.v79))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(JASPAR2020))
suppressMessages(library(TFBSTools))
suppressMessages(library(igraph))
suppressMessages(library(ggraph))
suppressMessages(library(MOJITOO))


# ========== 2. Set Working Directory and Define Data Paths ==========
setwd("/mnt/i/scDERI_Paper/3-benchmark/7-scMEGA/")
list = c(
  "/mnt/i/scDERI_Paper/3-benchmark/0-downsampling/Downsampled-data/10X_hg38_Human_bone_marrow_with_T_cell_depleted/10X_hg38_Human_bone_marrow_with_T_cell_depleted",
  "/mnt/i/scDERI_Paper/3-benchmark/0-downsampling/Downsampled-data/SNAREseq2_hg38_A549/SNAREseq2_hg38_A549",
  "/mnt/i/scDERI_Paper/3-benchmark/0-downsampling/Downsampled-data/SNAREseq2_hg38_GM12878/SNAREseq2_hg38_GM12878"
)


# ========== 3. Load Motif Data ==========
# Get motif PFMs from JASPAR CORE database
pfm = getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)


# ========== 4. Load and Preprocess each Dataset ==========
# Run each dataset individually (by commenting and uncommenting the corresponding sections).
# name = list[1]
# name = list[2]
name = list[3]

qs = qread(file = paste0(name,'.qsave'))
Genome = 'hg38'

# 4.1 Separate ATAC and RNA assays
qs@active.assay = 'ATAC'
obj.atac = qs
qs@active.assay = 'RNA'
obj.rna = qs

# 4.2 Add gene annotations
annotations = suppressWarnings(GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86, verbose = FALSE))
seqlevelsStyle(annotations) = 'UCSC'
Annotation(obj.atac) = annotations

# 4.3 Convert Assay to Assay5 if needed
if (grepl("Assay", class(qs[["RNA"]]))) {qs[["RNA"]] = as(qs[["RNA"]], "Assay5")}
if (grepl("Assay", class(qs[["ATAC"]]))) {qs[["ATAC"]] = as(qs[["ATAC"]], "Assay5")}

# 4.4 RNA processing: SCTransform, PCA, UMAP
obj.rna = obj.rna %>%
  SCTransform(verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30, verbose = FALSE)

# 4.5 ATAC processing: TF-IDF, top features, SVD, UMAP
obj.atac = obj.atac %>%
  RunTFIDF() %>%
  FindTopFeatures() %>%
  RunSVD() %>%
  RunUMAP(reduction = 'lsi', dims = 2:30, verbose = FALSE)


# ========== 5. Multi-omic Integration with MOJITOO ==========
# 5.1 Create combined Seurat object
meta.data = obj.rna@meta.data %>% as.data.frame()
pbmc = CreateSeuratObject(
  counts = LayerData(obj.rna, assay = "RNA", layer = "counts"),
  assay = "RNA",
  meta.data = meta.data
)
pbmc[["ATAC"]] = CreateChromatinAssay(
  counts = obj.atac@assays$ATAC@counts,
  sep = c(":", "-"),
  min.cells = 1,
  genome = Genome
)

# 5.2 RNA processing
DefaultAssay(pbmc) = "RNA"
pbmc = pbmc %>%
  NormalizeData(verbose=F) %>%
  FindVariableFeatures(nfeatures = 3000, verbose = FALSE) %>%
  ScaleData(verbose=F) %>%
  RunPCA(npcs = 50, reduction.name = "RNA_PCA", verbose = FALSE)

# 5.3 ATAC processing
DefaultAssay(pbmc) = "ATAC"
pbmc = pbmc %>%
  RunTFIDF(verbose=F) %>%
  FindTopFeatures(min.cutoff = 'q0', verbose = FALSE) %>%
  RunSVD(verbose=F)

# 5.4 MOJITOO integration
pbmc = mojitoo(
  object = pbmc,
  reduction.list = list("RNA_PCA", "lsi"),
  dims.list = list(1:50, 2:50), # exclude 1st dimension of LSI
  reduction.name = 'MOJITOO',
  assay = "RNA"
)

# 5.5 UMAP on MOJITOO embedding
DefaultAssay(pbmc) = "RNA"
embedd = Embeddings(pbmc[["MOJITOO"]])
pbmc = RunUMAP(pbmc,
               reduction = "MOJITOO",
               reduction.name = "MOJITOO_UMAP",
               dims = 1:ncol(embedd),
               verbose = FALSE)


# ========== 6. Clustering and Trajectory Analysis ==========
# 6.1 Dataset-specific clustering resolution
if (basename(name) == "10X_hg38_Human_bone_marrow_with_T_cell_depleted") {
  resolution = 0.4
} else if (basename(name) == "SNAREseq2_hg38_A549") {
  resolution = 0.8
} else if (basename(name) == "SNAREseq2_hg38_GM12878") {
  resolution = 0.7
}

# 6.2 Find neighbors and clusters
pbmc = FindNeighbors(pbmc, reduction = "MOJITOO", dims = 1:30)
pbmc = FindClusters(pbmc, resolution = resolution, verbose = TRUE)

# 6.3 Visualize clusters and add trajectory based on dimension reduction visualization image
# DimPlot(pbmc, group.by = "celltype", shuffle = TRUE, label = TRUE, reduction = "MOJITOO_UMAP") + NoLegend()
# DimPlot(pbmc, group.by = paste0("RNA_snn_res.", resolution), shuffle = TRUE, label = TRUE, reduction = "MOJITOO_UMAP") + NoLegend()
DimPlot(pbmc, group.by = "seurat_clusters", shuffle = TRUE, label = TRUE, reduction = "MOJITOO_UMAP") + NoLegend()
if (basename(name) == "10X_hg38_Human_bone_marrow_with_T_cell_depleted") { # ***set trajectory by yourself****
  trajectory = c(2, 0, 1)
} else if (basename(name) == "SNAREseq2_hg38_A549") {
  trajectory = c(0, 1)
} else if (basename(name) == "SNAREseq2_hg38_GM12878") {
  trajectory = c(0, 1)
}
pbmc = AddTrajectory(
  object = pbmc,
  trajectory = trajectory,
  group.by = paste0("RNA_snn_res.", resolution),
  reduction = "MOJITOO_UMAP",
  dims = 1:2,
  use.all = TRUE
)
pbmc.t.cells = pbmc[, !is.na(pbmc$Trajectory)]


# ========== 7. Transcription Factor and Gene Selection ==========
# 7.1 Add motifs and run chromVAR
pbmc.t.cells = AddMotifs(
  object = pbmc.t.cells,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm,
  assay = "ATAC"
)
pbmc.t.cells = RunChromVAR(
  object = pbmc.t.cells,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  assay = "ATAC"
)

# 7.2 Select TFs and genes
sel.tfs = SelectTFs(object = pbmc.t.cells)
df.cor = sel.tfs$tfs
sel.genes = SelectGenes(object = pbmc.t.cells, labelTop1 = 0, labelTop2 = 0)
df.p2g = sel.genes$p2g

# ========== 8. Gene Regulatory Network Inference ==========
# 8.1 Gene regulatory network inference and visualization
tf.gene.cor = GetTFGeneCorrelation(object = pbmc.t.cells,
                                   tf.use = df.cor$tfs,
                                   gene.use = unique(df.p2g$gene),
                                   tf.assay = "chromvar",
                                   gene.assay = "RNA",
                                   trajectory.name = "Trajectory")

# 8.2 Generate heatmap
# ht = GRNHeatmap(tf.gene.cor, tf.timepoint = df.cor$time_point)
# ht


# ========== 9. Save Results ==========
time_point = df.cor$time_point
save(list = c("tf.gene.cor", "ht", "time_point"), file = paste0(basename(name), '.RData'))
save.image(paste0("FINAL_", basename(name), '.RData'))

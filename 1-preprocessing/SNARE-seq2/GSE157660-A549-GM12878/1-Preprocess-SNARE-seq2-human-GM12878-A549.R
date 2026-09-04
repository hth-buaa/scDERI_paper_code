#############################################################################################
#                                                                                           #
#   Preprocess the single-cell data set from SNARE-seq2 on human GM12878 and A549 samples   #
#                                                                                           #
#############################################################################################


# 1. Libraries
library(easypackages)
libs = c("qs", "Seurat", "Signac", "EnsDb.Hsapiens.v86")
libraries(libs)


# 2. Parameters
data.dir = "/mnt/i/scDERI_Paper/1-preprocessing/SNARE-seq2/GSE157660-A549-GM12878/rda/"
dir.create(data.dir)
setwd(data.dir)


# 3.1 Download https://ftp.ncbi.nlm.nih.gov/geo/series/GSE157nnn/GSE157660/suppl/GSE157660%5Fhuman.snare2%2Drna.matrix.tar.gz and unzip
# 3.2 Load the RNA data
rna_counts = Seurat::Read10X(data.dir = "../human.snare2-rna.matrix")
dim(rna_counts)
# [1] 29954 10762
rna_counts[1:4, 1:4]
# 4 x 4 sparse Matrix of class "dgCMatrix"
#          A549_cells_fixedfrozen_P1_CATCAAGTAGATCGCAACAAGCTA A549_cells_fixedfrozen_P1_TGGAACAACGGATTGCAGTACAAG
# A1BG                                                      .                                                  .
# A1BG-AS1                                                  .                                                  .
# A1CF                                                      .                                                  .
# A2M                                                       .                                                  .
#          A549_cells_fixedfrozen_P1_ACACGACCTCCGTCTAAGTACAAG A549_cells_fixedfrozen_P1_CATCAAGTATCCTGTAACAAGCTA
# A1BG                                                      .                                                  .
# A1BG-AS1                                                  .                                                  .
# A1CF                                                      .                                                  .
# A2M                                                       1                                                  .


# 4.1 Download https://ftp.ncbi.nlm.nih.gov/geo/series/GSE157nnn/GSE157660/suppl/GSE157660%5Fhuman.snare2%2Datac.matrix.tar.gz and unzip
# 4.2 Load the ATAC data
atac_counts = Seurat::Read10X(data.dir = "../human.snare2-atac.matrix")
dim(atac_counts)
# [1] 143101  21657
atac_counts[1:4, 1:4]
# 4 x 4 sparse Matrix of class "dgCMatrix"
#                    GM12878_cells_1pctPFA_4C_exp1_P2_TGGAACAAAATGTTGCCAAGACTA GM12878_cells_1pctPFA_4C_exp1_P5_AACGCTTATATCAGCACAAGACTA
# chr1:804776-805277                                                         .                                                         .
# chr1:826804-827305                                                         .                                                         .
# chr1:827350-827851                                                         .                                                         .
# chr1:828129-828630                                                         .                                                         .
#                    GM12878_cells_1pctPFA_4C_exp1_P3_AAGGTACAGACTAGTAGAGCTGAA GM12878_cells_1pctPFA_4C_exp1_P4_CCGTGAGATGGTGGTAGAGCTGAA
# chr1:804776-805277                                                         .                                                         .
# chr1:826804-827305                                                         .                                                         .
# chr1:827350-827851                                                         .                                                         .
# chr1:828129-828630                                                         .                                                         .


# 5. Unify cells between GEX and ATAC
common.cells = intersect(colnames(rna_counts), colnames(atac_counts))
rna_counts = rna_counts[, common.cells]
atac_counts = atac_counts[, common.cells]
dim(rna_counts)
# [1] 29954 10093
dim(atac_counts)
# [1] 143101  10093
identical(colnames(rna_counts), colnames(atac_counts))
# [1] TRUE


# 6. Get gene annotations for hg38
annotation = Signac::GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)
Signac::seqlevels(annotation) = paste0('chr', Signac::seqlevels(annotation))


# 7. Create a Seurat object containing the RNA data
sc = SeuratObject::CreateSeuratObject(
  counts = rna_counts,
  assay = "RNA"
)
dim(sc)
# [1] 29954 10093


# 8. Create ATAC assay and add it to the object
sc[["ATAC"]] = Signac::CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  annotation = annotation
)
dim(sc[["ATAC"]])
# [1] 143101  10093


# 9. Split the dataset into a list of two seurat objects (GM12878 and A549)
ifnb.list = Seurat::SplitObject(sc, split.by = "orig.ident")
names(ifnb.list)
# [1] "A549"    "GM12878"
A549 = ifnb.list[[1]]
GM = ifnb.list[[2]]


# 10. Quality control
DefaultAssay(A549) = "RNA"
A549[["percent.mt"]] = Seurat::PercentageFeatureSet(A549, pattern = "^MT-")
# qs::qsave(A549, "./A549_Seu_before_QC.qsave")

DefaultAssay(GM) = "RNA"
GM[["percent.mt"]] = Seurat::PercentageFeatureSet(GM, pattern = "^MT-")
# qs::qsave(GM, "./GM_Seu_before_QC.qsave")


# 11. Read metadata
png("./A549_RNA_ATAC_MT_before_QC.png", res = 150, width = 2000, height = 1000)
VlnPlot(object = A549, features = c("nCount_RNA", "nCount_ATAC", "percent.mt"), ncol = 4, pt.size = 0)
dev.off()

png("GM_RNA_ATAC_MT_before_QC.png", res = 150, width = 2000, height = 1000)
VlnPlot(object = GM, features = c("nCount_RNA", "nCount_ATAC", "percent.mt"), ncol = 4, pt.size = 0)
dev.off()


# 12. Filter out low quality cells
A549 = S4Vectors::subset(
  x = A549,
  subset = nCount_RNA < 15000 &
    nCount_RNA > 200 &
    nCount_ATAC < 3500 &
    nCount_ATAC > 200 &
    # nucleosome_signal < 1 &
    # TSS.enrichment > 2 &
    percent.mt < 5
)
dim(A549)
# [1] 29954  1088

GM = S4Vectors::subset(
  x = GM,
  subset = nCount_RNA < 30000 &
    nCount_RNA > 200 &
    nCount_ATAC < 10000 &
    nCount_ATAC > 200 &
    # nucleosome_signal < 1 &
    # TSS.enrichment > 2 &
    percent.mt < 2
)
dim(GM)
# [1] 29954  8876


# 13. Simpily cell names
cells = colnames(A549)
head(cells)
# [1] "A549_cells_fixedfrozen_P1_AACAACCAGCCAAGACACAAGCTA" "A549_cells_fixedfrozen_P1_ACAGCAGAACGCTCGAACAAGCTA"
# [3] "A549_cells_fixedfrozen_P1_AAGAGATCGCCAAGACAGTACAAG" "A549_cells_fixedfrozen_P1_TGGTGGTAGCGAGTAAAGTACAAG"
# [5] "A549_cells_fixedfrozen_P1_AAGACGGAGCCACATAACAAGCTA" "A549_cells_fixedfrozen_P1_CACCTTACGAGTTAGCAGTACAAG"
length(cells)
# [1] 1088
cells1 = sub("^.*_", "", cells)
length(cells1)
# [1] 1088
cells2 = unique(cells1)
length(cells2)
# [1] 1088
colnames(A549) = cells2
head(colnames(A549))
# [1] "AACAACCAGCCAAGACACAAGCTA" "ACAGCAGAACGCTCGAACAAGCTA" "AAGAGATCGCCAAGACAGTACAAG" "TGGTGGTAGCGAGTAAAGTACAAG"
# [5] "AAGACGGAGCCACATAACAAGCTA" "CACCTTACGAGTTAGCAGTACAAG"

qs::qsave(A549, "./SNAREseq2_hg38_A549.qsave")

cells = colnames(GM)
head(cells)
# [1] "GM12878_cells_1pctPFA_4C_exp1_P6_ACAAGCTAAACGTGATGAGCTGAA"  "GM12878_cells_1pctPFA_4C_exp1_P9_ACACAGAAAGCAGGAAGAGCTGAA"
# [3] "GM12878_cells_1pctPFA_4C_exp1_P10_GTGTTCTACGAACTTAGAGCTGAA" "GM12878_cells_1pctPFA_4C_exp1_P5_GCTAACGAATCCTGTAGAGCTGAA"
# [5] "GM12878_cells_1pctPFA_4C_exp1_P5_AAGGACACAACTCACCGAGCTGAA"  "GM12878_cells_1pctPFA_4C_exp1_P5_TTCACGCAACGCTCGAGAGCTGAA"
length(cells)
# [1] 8876
cells1 = sub("^.*_", "", cells)
length(cells1)
# [1] 8876
cells2 = unique(cells1)
length(cells2)
# [1] 7426
# Cell name simplification is skipped to avoid substantial cell loss.

qs::qsave(GM, "./SNAREseq2_hg38_GM12878.qsave")


# 14. Save
rm(cells, cells1, cells2)
save.image("./1-Preprocess-SNARE-seq2-human-GM12878-A549.RData")

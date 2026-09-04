#############################################################################################
#                                                                                           #
#    Preprocess the single-cell data set from 10X Genomics Multiome on human bone marrow    #
#                                                                                           #
#############################################################################################


# 1. Install HDF5 bash (run code below in wsl)
mkdir -p /mnt/i/scDERI_Paper/1-preprocessing/10X_hg38_Human_bone_marrow_with_T_cell_depleted/hdf5
cd /mnt/i/scDERI_Paper/1-preprocessing/10X_hg38_Human_bone_marrow_with_T_cell_depleted/hdf5
wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.6/src/hdf5-1.10.6.tar.gz
tar -xzvf hdf5-1.10.6.tar.gz
cd hdf5-1.10.6
./configure --prefix=/mnt/i/scDERI_Paper/1-preprocessing/10X_hg38_Human_bone_marrow_with_T_cell_depleted/hdf5/hdf5_1.10.6
make -j 8
make install


# 2. Libraries
dyn.load(x = "/mnt/i/scDERI_Paper/1-preprocessing/10X_hg38_Human_bone_marrow_with_T_cell_depleted/hdf5/hdf5_1.10.6/lib/libhdf5.so.103")

library(easypackages)
libs = c("qs", "hdf5r", "Seurat", "SeuratDisk", "Signac", "EnsDb.Hsapiens.v86")
libraries(libs)


# 3. Parameters
data.dir = "/mnt/i/scDERI_Paper/1-preprocessing/10X_hg38_Human_bone_marrow_with_T_cell_depleted/10X-human-bone-marrow/rda/"
dir.create(data.dir)
setwd(data.dir)


# 4. Convert h5ad to h5Seurat for GEX assay
SeuratDisk::Convert("../bm_multiome_rna.h5ad", dest = "h5seurat", overwrite = TRUE)
sc_rna = SeuratDisk::LoadH5Seurat("../bm_multiome_rna.h5seurat")
sc_rna
head(rownames(sc_rna))
head(colnames(sc_rna))


# 5. Convert gene names from integers to symbols
rownames(sc_rna@assays[["RNA"]]@counts) = read.table(paste0(data.dir, "../bm_multiome_rna.txt"))[[1]]
rownames(sc_rna@assays[["RNA"]]@data) = rownames(sc_rna@assays[["RNA"]]@counts)
rownames(sc_rna@assays[['RNA']]@meta.features) = rownames(sc_rna@assays[["RNA"]]@counts)


# 6. Simplify cell name
# Based on preliminary data observation,
# it was found that only the portion after the "#" in the cell names of the sequencing data matches the cell names in the _fragments.tsv.gz file.
# Therefore, cell names need to be streamlined.
# A small number of duplicate streamlined cell names were identified,
# and these were handled by retaining only the first occurrence of each duplicate name.
original_names = colnames(sc_rna)
simplified_names = sub(".*#", "", original_names)
dup_indices = which(duplicated(simplified_names))
if (length(dup_indices) > 0) {
  cat("Found ", length(dup_indices), " duplicate simplified cell names\n")
  cat("Duplicate simplified cell names:\n")
  print(unique(simplified_names[dup_indices]))

  keep_indices = !duplicated(simplified_names)
  sc_rna_unique = sc_rna[, keep_indices]
  colnames(sc_rna_unique) = simplified_names[keep_indices]

  cat("\nProcessing results:\n")
  cat("Original number of cells: ", length(original_names), "\n")
  cat("Number of cells after deduplication: ", ncol(sc_rna_unique), "\n")
  cat("Removed ", length(original_names) - ncol(sc_rna_unique), " duplicate cells\n")
} else {
  colnames(sc_rna) = simplified_names
  cat("No duplicate cell names found, successfully renamed all cells\n")
}


# 7. Convert h5ad to h5Seurat for ATAC assay
SeuratDisk::Convert("../bm_multiome_atac.h5ad", dest = "h5seurat", overwrite = TRUE)
sc_atac = SeuratDisk::LoadH5Seurat("../bm_multiome_atac.h5seurat")
sc_atac
head(rownames(sc_atac))
head(colnames(sc_atac))


# 8. Simplify cell name
original_names = colnames(sc_atac)
simplified_names = sub(".*#", "", original_names)
dup_indices = which(duplicated(simplified_names))
if (length(dup_indices) > 0) {
  cat("Found ", length(dup_indices), " duplicate simplified cell names\n")
  cat("Duplicate simplified cell names:\n")
  print(unique(simplified_names[dup_indices]))

  keep_indices = !duplicated(simplified_names)
  sc_atac_unique = sc_atac[, keep_indices]
  colnames(sc_atac_unique) = simplified_names[keep_indices]

  cat("\nProcessing results:\n")
  cat("Original number of cells: ", length(original_names), "\n")
  cat("Number of cells after deduplication: ", ncol(sc_atac_unique), "\n")
  cat("Removed ", length(original_names) - ncol(sc_atac_unique), " duplicate cells\n")
} else {
  colnames(sc_atac) = simplified_names
  cat("No duplicate cell names found, successfully renamed all cells\n")
}


# 9. Remove cells that are not present in both omics datasets
# Fortunately, all cells are present in the dual omics data.
table(sort(colnames(sc_rna_unique)) == sort(colnames(sc_atac_unique)))
# TRUE
# 7418


# 10. Get gene annotations for hg38
annotation = Signac::GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)
Signac::seqlevels(annotation) = paste0('chr', Signac::seqlevels(annotation))


# 8. Construct Seurat object composed of both ATAC and GEX
sc = sc_rna_unique
sc[['ATAC']] = Signac::CreateChromatinAssay(
  counts = sc_atac_unique@assays$RNA@counts,
  sep = c(":", "-"),
  annotation = annotation
)
sc
head(sc@meta.data)
dim(sc@meta.data)
colnames(sc@meta.data)
sc@meta.data[1:4, 1:4]
sc@assays[['RNA']]@counts[1:4, 1:4]
sc@assays[['ATAC']]@counts[1:4, 1:4]


# 9. Save
qs::qsave(sc, "./10X_hg38_Human_bone_marrow_with_T_cell_depleted.qsave")
save.image("./3-Preprocess-10X-human-bone-marrow.Rdata")

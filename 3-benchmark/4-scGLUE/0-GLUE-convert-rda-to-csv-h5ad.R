#####################################################################
#                                                                   #
#   Convert Seurat object to annData that can be accepted by GLUE   #
#                                                                   #
#####################################################################


# General libraries
dyn.load(x = "/mnt/i/scDERI_Paper/1-preprocessing/10X_hg38_Human_bone_marrow_with_T_cell_depleted/hdf5/hdf5_1.10.6/lib/libhdf5.so.103")
library(easypackages)
libs = c(
  "qs",
  "hdf5r",
  "Seurat",
  "Signac",
  "pbmcapply",
  "pbapply",
  "parallel",
  "EnsDb.Mmusculus.v79",
  "EnsDb.Hsapiens.v86",
  "dplyr",
  "ggplot2",
  "ggpubr",
  "SeuratDisk"
)
libraries(libs)


# Parameters
data.dir = "/mnt/i/scDERI_Paper/3-benchmark/4-scGLUE/"
setwd(data.dir)


# Conversion
message("Began converting Seurat objects")
pbmclapply(
  X = c(
    "/mnt/i/scDERI_Paper/3-benchmark/0-downsampling/Downsampled-data/10X_hg38_Human_bone_marrow_with_T_cell_depleted/10X_hg38_Human_bone_marrow_with_T_cell_depleted.qsave",
    "/mnt/i/scDERI_Paper/3-benchmark/0-downsampling/Downsampled-data/SNAREseq2_hg38_A549/SNAREseq2_hg38_A549.qsave",
    "/mnt/i/scDERI_Paper/3-benchmark/0-downsampling/Downsampled-data/SNAREseq2_hg38_GM12878/SNAREseq2_hg38_GM12878.qsave"
  ),
  mc.cores = detectCores() - 1,
  FUN = function(file.path) {
    message("----> Converting the file ", file.path, " ...")
    # file.path = "/mnt/i/scDERI_Paper/3-benchmark/0-downsampling/Downsampled-data/10X_hg38_Human_bone_marrow_with_T_cell_depleted/10X_hg38_Human_bone_marrow_with_T_cell_depleted.qsave"
    name = strsplit(file.path, split = "/") %>% `[[` (1) %>% tail(., n = 1) %>% gsub(".qsave", "", .)

    # Load data and generate scATAC-seq bed file
    obj = qs::qread(file.path)
    dim(obj)
    if (grepl("Assay5", class(obj[["RNA"]]))) {obj[["RNA"]] = as(obj[["RNA"]], "Assay")}
    if (grepl("Assay5", class(obj[["ATAC"]]))) {obj[["ATAC"]] = as(obj[["ATAC"]], "Assay")}


    # Convert the RNA and ATAC assays of Seurat object to Scanpy annData
    SaveH5Seurat(obj, filename = paste0(data.dir, name, ".h5Seurat"), overwrite = TRUE)
    Convert(paste0(data.dir, name, ".h5Seurat"), dest = "h5ad", assay = 'RNA', overwrite = TRUE)
    system(paste0("mv -f ", data.dir, name, ".h5ad", " ", data.dir, name, "-GEX.h5ad"))
    Convert(paste0(data.dir, name, ".h5Seurat"), dest = "h5ad", assay = 'ATAC', overwrite = TRUE)
    system(paste0("mv -f ", data.dir, name, ".h5ad", " ", data.dir, name, "-ATAC.h5ad"))


    # Convert the RNA and ATAC assays of Seurat object to csv files
    write.csv(obj[['RNA']]@counts, file = paste0(data.dir, name, "-GEX.csv"), quote = FALSE)
    write.csv(obj[['ATAC']]@counts, file = paste0(data.dir, name, "-ATAC.csv"), quote = FALSE)
  }
)
message("Finished converting Seurat objects")
save.image(paste0(data.dir, "1-GLUE-convert-rda-to-h5ad.RData"))

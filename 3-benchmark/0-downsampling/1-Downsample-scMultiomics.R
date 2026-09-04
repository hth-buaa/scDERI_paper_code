#####################################################################################################
#                                                                                                   #
#   Downsample data to incorporate 1,000 cells, 10,000 genes, and the first 10,000 variable peaks   #
#                                                                                                   #
#####################################################################################################


# 0. Task list
# (1) Remove the non-gene GEX features, which begin with "AC", "AL", or "LINC" followed by numbers, as well as mitochondrion genes starting with "MT" or "mt"
# (2) Remove peaks on non-standard chromosomes (only keep chr1-22, X, Y)
# (3) Select the top-ranked 10,000 variable genes
# (4) Select the top-ranked 10,000 variable peaks
# (5) Randomly select 1,000 cells (or all if less than 1,000)


# 1. Load libraries
library(easypackages)
libs = c("qs", "Seurat", "Signac", "pbapply", "dplyr", "stringr")
libraries(libs)


# 2. Setup paths
parent.dir = "/mnt/i/scDERI_Paper/3-benchmark/0-downsampling/"
data.dir = "/mnt/i/scDERI_Paper/3-benchmark/0-downsampling/Downsampled-data/"
dir.create("/mnt/i/scDERI_Paper/3-benchmark/0-downsampling/Downsampled-data/")
setwd(data.dir)

# Store the single-cell data prior to downsampling (Seurat object, _fragments.tsv.gz, _fragments.tsv.gz.tbi) (copy the preprocessed data here).
dir.create("/mnt/i/scDERI_Paper/3-benchmark/0-downsampling/10X_hg38_Human_bone_marrow_with_T_cell_depleted/")
dir.create("/mnt/i/scDERI_Paper/3-benchmark/0-downsampling/SNAREseq2_hg38_A549/")
dir.create("/mnt/i/scDERI_Paper/3-benchmark/0-downsampling/SNAREseq2_hg38_GM12878/")

system("cp /mnt/i/scDERI_Paper/1-preprocessing/10X_hg38_Human_bone_marrow_with_T_cell_depleted/10X-human-bone-marrow/rda/10X_hg38_Human_bone_marrow_with_T_cell_depleted.qsave                /mnt/i/scDERI_Paper/3-benchmark/0-downsampling/10X_hg38_Human_bone_marrow_with_T_cell_depleted/10X_hg38_Human_bone_marrow_with_T_cell_depleted.qsave")
system("cp /mnt/i/scDERI_Paper/1-preprocessing/10X_hg38_Human_bone_marrow_with_T_cell_depleted/10X-human-bone-marrow/rda/10X_hg38_Human_bone_marrow_with_T_cell_depleted_fragments.tsv.gz     /mnt/i/scDERI_Paper/3-benchmark/0-downsampling/10X_hg38_Human_bone_marrow_with_T_cell_depleted/10X_hg38_Human_bone_marrow_with_T_cell_depleted_fragments.tsv.gz")
system("cp /mnt/i/scDERI_Paper/1-preprocessing/10X_hg38_Human_bone_marrow_with_T_cell_depleted/10X-human-bone-marrow/rda/10X_hg38_Human_bone_marrow_with_T_cell_depleted_fragments.tsv.gz.tbi /mnt/i/scDERI_Paper/3-benchmark/0-downsampling/10X_hg38_Human_bone_marrow_with_T_cell_depleted/10X_hg38_Human_bone_marrow_with_T_cell_depleted_fragments.tsv.gz.tbi")

system("cp /mnt/i/scDERI_Paper/1-preprocessing/SNARE-seq2/GSE157660-A549-GM12878/rda/SNAREseq2_hg38_A549.qsave                       /mnt/i/scDERI_Paper/3-benchmark/0-downsampling/SNAREseq2_hg38_A549/SNAREseq2_hg38_A549.qsave")
system("cp /mnt/i/scDERI_Paper/1-preprocessing/SNARE-seq2/GSE157660-A549-GM12878/rda/SNAREseq2_hg38_A549_fragments.sorted.tsv.gz     /mnt/i/scDERI_Paper/3-benchmark/0-downsampling/SNAREseq2_hg38_A549/SNAREseq2_hg38_A549_fragments.tsv.gz")
system("cp /mnt/i/scDERI_Paper/1-preprocessing/SNARE-seq2/GSE157660-A549-GM12878/rda/SNAREseq2_hg38_A549_fragments.sorted.tsv.gz.tbi /mnt/i/scDERI_Paper/3-benchmark/0-downsampling/SNAREseq2_hg38_A549/SNAREseq2_hg38_A549_fragments.tsv.gz.tbi")

system("cp /mnt/i/scDERI_Paper/1-preprocessing/SNARE-seq2/GSE157660-A549-GM12878/rda/SNAREseq2_hg38_GM12878.qsave                       /mnt/i/scDERI_Paper/3-benchmark/0-downsampling/SNAREseq2_hg38_GM12878/SNAREseq2_hg38_GM12878.qsave")
system("cp /mnt/i/scDERI_Paper/1-preprocessing/SNARE-seq2/GSE157660-A549-GM12878/rda/SNAREseq2_hg38_GM12878_fragments.sorted.tsv.gz     /mnt/i/scDERI_Paper/3-benchmark/0-downsampling/SNAREseq2_hg38_GM12878/SNAREseq2_hg38_GM12878_fragments.tsv.gz")
system("cp /mnt/i/scDERI_Paper/1-preprocessing/SNARE-seq2/GSE157660-A549-GM12878/rda/SNAREseq2_hg38_GM12878_fragments.sorted.tsv.gz.tbi /mnt/i/scDERI_Paper/3-benchmark/0-downsampling/SNAREseq2_hg38_GM12878/SNAREseq2_hg38_GM12878_fragments.tsv.gz.tbi")

# Store the downsampling data (Seurat object, _fragments.tsv.gz, _fragments.tsv.gz.tbi).
dir.create("/mnt/i/scDERI_Paper/3-benchmark/0-downsampling/Downsampled-data/10X_hg38_Human_bone_marrow_with_T_cell_depleted/")
dir.create("/mnt/i/scDERI_Paper/3-benchmark/0-downsampling/Downsampled-data/SNAREseq2_hg38_A549/")
dir.create("/mnt/i/scDERI_Paper/3-benchmark/0-downsampling/Downsampled-data/SNAREseq2_hg38_GM12878/")


# 3. List input directories
dir.list = setdiff(list.dirs(parent.dir, full.names = FALSE, recursive = FALSE), "Downsampled-data")
length(dir.list)
dir.list
# [1] "10X_hg38_Human_bone_marrow_with_T_cell_depleted" "SNAREseq2_hg38_A549"                             "SNAREseq2_hg38_GM12878"


# 4. Main processing loop
pbapply::pblapply(dir.list, function(x) {
  # x = dir.list[1] # for test

  message(x)
  # dir.create(paste0(data.dir, x))
  rda.list = list.files(paste0(parent.dir, x), pattern = ".qsave")
  frag.list = list.files(paste0(parent.dir, x), pattern = "_fragments.tsv.gz")

  lapply(rda.list, function(y) {
    # y = rda.list[1] # for test

    # 4.1 Load data
    message("----> ", y)
    org = stringr::str_extract(y, stringr::regex("hg19|hg38|mm10|mm9", ignore_case = TRUE))
    if (grepl("^mm", org)) {
      mt.prefix <- "mt-"
    } else {
      mt.prefix <- "MT-"
    }
    obj = qs::qread(paste0(parent.dir, x, "/", y))
    message("----> Cells: ", ncol(obj), "\n",
            "----> Genes: ", nrow(obj[['RNA']]), "\n",
            "----> Peaks: ", nrow(obj[['ATAC']]), "\n")


    # 4.2 Random cell sampling (Task 5)
    # obj.save = obj # for test
    n.cells = min(1000, ncol(obj))
    obj = S4Vectors::subset(obj, cells = colnames(obj)[sample(x = 1:ncol(obj), size = n.cells)])
    message("----> ", ncol(obj), " cells were selected after downsampling")


    # 4.3 Filter non-gene and MT-genes features (Task 1)
    genes = setdiff(rownames(obj[['RNA']]), grep("(^AC|^AL|^LINC)[0-9]+", rownames(obj[['RNA']]), value = TRUE))
    message("----> ", length(genes), " genes after filtering out primer names")
    genes = grep(mt.prefix, genes, invert = TRUE, value = TRUE)
    message("----> ", length(genes), " genes after filtering out MT-genes")


    # 4.4 Select top variable genes (Task 3)
    DefaultAssay(obj) = "RNA"
    obj = Seurat::NormalizeData(obj)
    obj = Seurat::FindVariableFeatures(obj, selection.method = "vst")
    if ("meta.features" %in% slotNames(obj[["RNA"]])) {
      rna.meta = obj[['RNA']]@meta.features
      rna.meta = rna.meta[order(rna.meta$vst.variance.standardized, decreasing = TRUE),]
    } else if ("meta.data" %in% slotNames(obj[["RNA"]])) {
      rna.meta = obj[['RNA']]@meta.data
      rownames(rna.meta) = rownames(obj)
      rna.meta = rna.meta[order(rna.meta$vf_vst_counts_variance.standardized, decreasing = TRUE),]
    } else {stop("There is no meta.")}
    genes = rna.meta[rownames(rna.meta) %in% genes,] %>% head(n = 10000) %>% rownames
    message("----> ", length(genes), " genes after selecting highly variable genes")


    # 4.5 Filter out peaks on non-standard chromosomes (Task 2)
    peaks = rownames(obj[['ATAC']])
    peaks = grep("chr([0-9]|X|Y)", peaks, value = TRUE)
    message("----> ", length(peaks), " peaks after filtering out peaks on non-standard chromosomes")


    # 4.6 Select top variable peaks (Task 4)
    DefaultAssay(obj) = "ATAC"
    obj = Signac::FindTopFeatures(obj)
    atac.meta = obj[['ATAC']]@meta.features
    atac.meta = atac.meta[order(atac.meta$percentile, decreasing = TRUE), ]
    peaks = atac.meta[rownames(atac.meta) %in% peaks, ] %>% head(n = 10000) %>% rownames
    message("----> ", length(peaks), " peaks after selecting highly variable peaks")


    # 4.7 Create the subsetted Seurat object
    obj.sub = S4Vectors::subset(obj, features = c(genes, peaks))
    message("----> ", nrow(obj.sub[['RNA']]), " genes, ",
            nrow(obj.sub[['ATAC']]), " peaks, and ", ncol(obj.sub),
            " cells were retained after downsampling")


    # 4.8 Save the Seurat object and copy fragments file
    qs::qsave(obj.sub, paste0(data.dir, x, "/", y))
    if (length(frag.list) > 0) {
      system(paste0("cp ", data.dir, "../", x, "/", gsub(".qsave", "_fragments.tsv.gz ", y), data.dir, x))
      system(paste0("cp ", data.dir, "../", x, "/", gsub(".qsave", "_fragments.tsv.gz.tbi ", y), data.dir, x))
    }
    message("----> Finished copying ", gsub(".qsave", "_fragments.tsv.gz ", y, "\n\n"))
  })
})

#################################################
#                                               #
#       Apply scDERI on analysed dataset        #
#                                               #
#################################################
# Use R v4.4
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list = ls())
gc()

# Install scDERI from GitHub (run once)
# devtools::install_github("hth-buaa/scDERI")
library(scDERI)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(DDRTree)
library(BSgenome.Hsapiens.UCSC.hg38)

# OS configuration for path switching
sys = "wsl2"
# sys = "win11"

# Set working directory based on operating system.
if (sys == "wsl2") {
  Working_directory = "/mnt/i/scDERI_Paper/4-apply/"
} else if (sys == "win11") {
  Working_directory = "I:/scDERI_Paper/4-apply/"
} else {stop("Error.")}
setwd(Working_directory)




########## Part 0: Loading Data ##########
##### 0.1 Download raw data to "./raw_data/" #####
### Download scRNA-seq data of KP4 cells (a pancreatic ductal adenocarcinoma cell line).
# Data website: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM9343186
# Download GSM9343186_barcodes.tsv.gz website: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM9343186&format=file&file=GSM9343186%5Fbarcodes%2Etsv%2Egz
# Download GSM9343186_features.tsv.gz website: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM9343186&format=file&file=GSM9343186%5Ffeatures%2Etsv%2Egz
# Download GSM9343186_matrix.mtx.gz website: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM9343186&format=file&file=GSM9343186%5Fmatrix%2Emtx%2Egz
# In reality, GSM9343186_matrix.mtx.gz contains the count matrices for both scRNA-seq and scATAC-seq.
# Upon inspection, these datasets have been preprocessed by the uploader, such as through quality control.

### Download scATAC-seq data of KP4 cells (a pancreatic ductal adenocarcinoma cell line).
# Data website: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM9343187
# Download GSM9343187_atac_fragments.tsv.gz website: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM9343187&format=file&file=GSM9343187%5Fatac%5Ffragments%2Etsv%2Egz

### Rename raw files to standard 10x format
setwd(Working_directory)
folder_path = "./raw_data/"
file_mappings = list(
  c("GSM9343186_barcodes.tsv.gz", "barcodes.tsv.gz"),
  c("GSM9343186_features.tsv.gz", "features.tsv.gz"),
  c("GSM9343186_matrix.mtx.gz", "matrix.mtx.gz"),
  c("GSM9343187_atac_fragments.tsv.gz", "atac_fragments.tsv.gz")
)

### Batch rename files with existence check.
for (mapping in file_mappings) {
  old_name = mapping[1]
  new_name = mapping[2]
  old_path = file.path(folder_path, old_name)
  new_path = file.path(folder_path, new_name)
  if (file.exists(old_path)) {
    file.rename(from = old_path, to = new_path)
    message("Renamed: ", old_name, " -> ", new_name)
  } else {
    message("File not found, skipping: ", old_name)
  }
}
message("File list in folder after renaming:")
print(list.files(folder_path)); rm(folder_path, file_mappings, mapping, old_name, new_name, old_path, new_path)
# [1] "atac_fragments.tsv.gz"     "atac_fragments.tsv.gz.tbi" "barcodes.tsv.gz"           "features.tsv.gz"           "matrix.mtx.gz"



##### 0.2 Generate tabix index for fragments file #####
### Run the following 3 commands in a Linux terminal:
sudo apt-get install tabix
cd /mnt/i/scDERI_Paper/4-apply/raw_data
tabix -p bed atac_fragments.tsv.gz



##### 0.3 Construct Seurat objects from raw data #####
setwd(Working_directory)

### Read 10x multiome count matrix.
counts_data = Seurat::Read10X(data.dir = "./raw_data/")
class(counts_data) # list with two assay modalities
names(counts_data) # "Gene Expression", "Peaks"
dim(counts_data[["Gene Expression"]]) # 36601 genes × 9511 cells
counts_data[["Gene Expression"]][1:5, 1:5]
# 5 x 5 sparse Matrix of class "dgCMatrix"
#             AAACAGCCAAAGCTAA-1 AAACAGCCAATCATGT-1 AAACAGCCACCTCAGG-1 AAACAGCCACGAATCC-1 AAACAGCCAGGCTACT-1
# MIR1302-2HG                  .                  .                  .                  .                  .
# FAM138A                      .                  .                  .                  .                  .
# OR4F5                        .                  .                  .                  .                  .
# AL627309.1                   .                  .                  .                  .                  .
# AL627309.3                   .                  .                  .                  .                  .
dim(counts_data[["Peaks"]]) # 111870 peaks × 9511 cells
counts_data[["Peaks"]][1:5, 1:5]
# 5 x 5 sparse Matrix of class "dgCMatrix"
#                    AAACAGCCAAAGCTAA-1 AAACAGCCAATCATGT-1 AAACAGCCACCTCAGG-1 AAACAGCCACGAATCC-1 AAACAGCCAGGCTACT-1
# chr1:9735-10683                     .                  4                  .                  .                  .
# chr1:15776-16632                    .                  .                  .                  .                  .
# chr1:17121-17951                    .                  .                  .                  .                  .
# chr1:180650-181082                  .                  .                  .                  .                  2
# chr1:181185-181785                  .                  .                  .                  .                  .

### Create initial Seurat object with RNA assay.
data = SeuratObject::CreateSeuratObject(
  counts = counts_data[["Gene Expression"]],
  assay = "RNA",
  project = "KP4",
  min.cells = 0,
  min.features = 0
)

### Add cell type annotation.
data@meta.data$seurat_annotations = "KP4"

### Build gene annotation from EnsDb.
annotation = Signac::GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)
Signac::seqlevels(annotation) = paste0('chr', Signac::seqlevels(annotation))

### Add ATAC ChromatinAssay to Seurat object.
data[["ATAC"]] = Signac::CreateChromatinAssay(
  counts = counts_data[["Peaks"]],
  sep = c(":", "-"),
  annotation = annotation,
  min.cells = 0,
  min.features = 0
)

### Extract standalone scRNA-seq Seurat object.
scRNAseq_raw = SeuratObject::CreateSeuratObject(
  counts = SeuratObject::GetAssayData(data, assay = "RNA", layer = "counts"),
  assay = "RNA",
  meta.data = data@meta.data
)

### Extract standalone scATAC-seq Seurat object.
data@active.assay = "ATAC"
data@assays[["RNA"]] = NULL
scATACseq_raw = data
rm(data, annotation, counts_data)

### Verify assay class types.
print(class(scRNAseq_raw[["RNA"]])) # Assay5 (SeuratObject)
print(class(scATACseq_raw[["ATAC"]])) # ChromatinAssay (Signac)

### Save raw Seurat objects.
qs::qsave(scRNAseq_raw, file = "./scRNAseq_raw.qsave")
qs::qsave(scATACseq_raw, file = "./scATACseq_raw.qsave")




##### 0.4 Removing Data from Abnormal Chromosomes #####
setwd(Working_directory)
scATACseq_raw = qs::qread("./scATACseq_raw.qsave")

### Summarize peak distribution across chromosomes.
print(as.data.frame(table(sapply(strsplit(rownames(scATACseq_raw), "-"), `[`, 1))))
#          Var1  Freq
# 1  GL000194.1     2
# 2  GL000195.1    14
# 3  GL000205.2    20
# 4  GL000218.1     7
# 5  GL000219.1     6
# 6  KI270711.1     6
# 7  KI270713.1     6
# 8  KI270726.1     1
# 9  KI270728.1     2
# 10 KI270731.1     1
# 11       chr1 10488
# 12      chr10  5594
# 13      chr11  5846
# 14      chr12  5843
# 15      chr13  1977
# 16      chr14  3467
# 17      chr15  3154
# 18      chr16  3424
# 19      chr17  5510
# 20      chr18  1901
# 21      chr19  4101
# 22       chr2  8287
# 23      chr20  2883
# 24      chr21  1401
# 25      chr22  2061
# 26       chr3  7277
# 27       chr4  5200
# 28       chr5  7011
# 29       chr6  6483
# 30       chr7  6852
# 31       chr8  4673
# 32       chr9  5017
# 33       chrX  3329
# 34       chrY    26

### Retain only standard chromosomes; remove non-canonical contigs.
peaks_keep = as.character(Signac::seqnames(Signac::granges(scATACseq_raw))) %in%
  GenomeInfoDb::standardChromosomes(Signac::granges(scATACseq_raw))
print(table(peaks_keep)) # 65 non-standard peaks removed, 111805 peaks retained.

### Subset all assay slots by retained peaks.
atac_assay = scATACseq_raw[["ATAC"]]
atac_assay@counts = atac_assay@counts[peaks_keep, ]
atac_assay@data = atac_assay@data[peaks_keep, ]
if (nrow(atac_assay@scale.data) > 0) {atac_assay@scale.data = atac_assay@scale.data[peaks_keep, ]}
atac_assay@ranges = atac_assay@ranges[peaks_keep]
atac_assay@meta.features = atac_assay@meta.features[peaks_keep, , drop = FALSE]
if (length(atac_assay@bias) > 0) {atac_assay@bias = atac_assay@bias[peaks_keep]}
if (length(atac_assay@motifs) > 0) {atac_assay@motifs@data = atac_assay@motifs@data[peaks_keep, ]}

### Update Seurat object with filtered assay.
scATACseq_clean_chro = scATACseq_raw
scATACseq_clean_chro[["ATAC"]] = atac_assay; rm(peaks_keep, atac_assay, scATACseq_raw)

### Verify chromosome filtering result.
print(as.data.frame(table(sapply(strsplit(rownames(scATACseq_clean_chro), "-"), `[`, 1))))
#     Var1  Freq
# 1   chr1 10488
# 2  chr10  5594
# 3  chr11  5846
# 4  chr12  5843
# 5  chr13  1977
# 6  chr14  3467
# 7  chr15  3154
# 8  chr16  3424
# 9  chr17  5510
# 10 chr18  1901
# 11 chr19  4101
# 12  chr2  8287
# 13 chr20  2883
# 14 chr21  1401
# 15 chr22  2061
# 16  chr3  7277
# 17  chr4  5200
# 18  chr5  7011
# 19  chr6  6483
# 20  chr7  6852
# 21  chr8  4673
# 22  chr9  5017
# 23  chrX  3329
# 24  chrY    26

# Save chromosome-filtered ATAC object.
qs::qsave(scATACseq_clean_chro, file = "./scATACseq_clean_chro.qsave")



##### 0.5 Downsampling #####
##### 0.5.1 Gene downsampling: remove non-coding RNAs and mitochondrial genes, keep top 10000 highly variable genes (HVGs)
setwd(Working_directory)
scRNAseq_raw = qs::qread("./scRNAseq_raw.qsave")

genes = setdiff(rownames(scRNAseq_raw), grep("(^AC|^AL|^LINC)[0-9]+", rownames(scRNAseq_raw), value = TRUE))
genes = grep("MT-", genes, invert = TRUE, value = TRUE)
scRNAseq = S4Vectors::subset(scRNAseq_raw, features = genes)
print(scRNAseq)
# An object of class Seurat
# 23060 features across 9511 samples within 1 assay
# Active assay: RNA (23060 features, 0 variable features)
#  1 layer present: counts

### Normalize counts and calculate gene variability.
scRNAseq = Seurat::NormalizeData(scRNAseq)
scRNAseq = Seurat::FindVariableFeatures(scRNAseq, selection.method = "vst", nfeatures = nrow(scRNAseq))

### Select top 10000 most variable genes.
RNA_var_meta = scRNAseq@assays[["RNA"]]@meta.data
rownames(RNA_var_meta) = rownames(scRNAseq)
RNA_var_meta = RNA_var_meta[order(RNA_var_meta$vf_vst_counts_variance.standardized, decreasing = TRUE),]
genes = rownames(RNA_var_meta)[1:10000]
scRNAseq = S4Vectors::subset(scRNAseq, features = genes)
rm(RNA_var_meta, genes, scRNAseq_raw)
print(scRNAseq)
# An object of class Seurat
# 10000 features across 9511 samples within 1 assay
# Active assay: RNA (10000 features, 10000 variable features)
#  2 layers present: counts, data

### Save downsampled RNA object.
qs::qsave(scRNAseq, file = "./scRNAseq.qsave")


##### 0.5.2 Peak downsampling: keep top 10,000 most variable chromatin peaks.
setwd(Working_directory)
scATACseq_clean_chro = qs::qread("./scATACseq_clean_chro.qsave")

scATACseq = scATACseq_clean_chro
scATACseq = Signac::RunTFIDF(scATACseq)
scATACseq = Signac::FindTopFeatures(scATACseq, min.cutoff = "q0")

### Rank peaks by accessibility percentile and select top 10000.
ATAC_var_meta = scATACseq@assays[["ATAC"]]@meta.features
ATAC_var_meta = ATAC_var_meta[order(ATAC_var_meta$percentile, decreasing = TRUE), ]
peaks = rownames(ATAC_var_meta)[1:10000]

peaks_keep = rownames(scATACseq) %in% peaks
print(table(peaks_keep)) # 101805 peaks removed, 10000 peaks retained.

### Subset all assay slots by retained peaks.
atac_assay = scATACseq[["ATAC"]]
atac_assay@counts = atac_assay@counts[peaks_keep, ]
atac_assay@data = atac_assay@data[peaks_keep, ]
if (nrow(atac_assay@scale.data) > 0) {atac_assay@scale.data = atac_assay@scale.data[peaks_keep, ]}
atac_assay@ranges = atac_assay@ranges[peaks_keep]
atac_assay@meta.features = atac_assay@meta.features[peaks_keep, , drop = FALSE]
if (length(atac_assay@bias) > 0) {atac_assay@bias = atac_assay@bias[peaks_keep]}
if (length(atac_assay@motifs) > 0) {atac_assay@motifs@data = atac_assay@motifs@data[peaks_keep, ]}
atac_assay@var.features = character(0)

# Update Seurat object with filtered peaks.
scATACseq[["ATAC"]] = atac_assay; rm(peaks, atac_assay, ATAC_var_meta, scATACseq_clean_chro)
print(scATACseq)
# An object of class Seurat
# 10000 features across 9511 samples within 1 assay
# Active assay: ATAC (10000 features, 0 variable features)
#  2 layers present: counts, data

# Save downsampled ATAC object.
qs::qsave(scATACseq, file = "./scATACseq.qsave")




########## Part 1: Data Preprocessing ##########
##### 1.1 Identifying TFs #####
setwd(Working_directory)
scRNAseq = qs::qread("./scRNAseq.qsave")
options(timeout = 60)

### Retrieve TF prior networks from ChIP-Atlas database.
monitor = peakRAM::peakRAM({
  TFs = get_TFs_from_priors(
    scRNAseq = scRNAseq,
    species = "Homo sapiens",
    collection = c("CORE", "CNE", "PHYLOFACTS", "SPLICE", "POLII", "FAM", "PBM", "PBM_HOMEO", "PBM_HLH", "UNVALIDATED"),
    genome_assembly = "hg38",
    min_nonzero_freq = 0, min_nonzero_mean = 0, min_overall_mean = 0,
    min_nonzero_freq_quantile = 0, min_nonzero_mean_quantile = 0, min_overall_mean_quantile = 0,
    distance_from_TSS = 10,
    excluded_cell_types_in_priors = "NotExclude",
    prior_min_TG_number = 0, prior_top_n_keep = 100000, prior_top_q_keep = 1
  )
})
message("Get ", length(TFs), " TFs, ", print_run_stats(monitor))
# Get 237 TFs, running: 4 hours 11 mins 28 secs, peak RAM: 0 GB 948 MB

qs::qsave(TFs, file = "./1.1 Data Preprocessing - Get TFs From Database/TFs.qsave")

### Check per-TF download status.
download_information = readRDS("./1.1 Data Preprocessing - Get TFs From Database/ChIP_Atlas_Priors/download_information.rds")
print(download_information[download_information$status == "URL not found or download error", ])
#                                                             url
# ESR1   https://chip-atlas.dbcls.jp/data/hg38/target/ESR1.10.tsv
# FOXA1 https://chip-atlas.dbcls.jp/data/hg38/target/FOXA1.10.tsv
# KLF7   https://chip-atlas.dbcls.jp/data/hg38/target/KLF7.10.tsv
#                                                                                                  destination                          status          time
# ESR1   /mnt/i/scDERI_Paper/4-apply/1.1 Data Preprocessing - Get TFs From Database/ChIP_Atlas_Priors/ESR1.tsv URL not found or download error 15.33907 mins
# FOXA1 /mnt/i/scDERI_Paper/4-apply/1.1 Data Preprocessing - Get TFs From Database/ChIP_Atlas_Priors/FOXA1.tsv URL not found or download error 14.13492 mins
# KLF7   /mnt/i/scDERI_Paper/4-apply/1.1 Data Preprocessing - Get TFs From Database/ChIP_Atlas_Priors/KLF7.tsv URL not found or download error 51.54475 secs

### Manually download missing files using URLs from download_information$url, save to paths in download_information$destination.

### Re-run TF retrieval after manual file supplementation.
monitor = peakRAM::peakRAM({
  TFs = get_TFs_from_priors(
    scRNAseq = scRNAseq,
    species = "Homo sapiens",
    collection = c("CORE", "CNE", "PHYLOFACTS", "SPLICE", "POLII", "FAM", "PBM", "PBM_HOMEO", "PBM_HLH", "UNVALIDATED"),
    genome_assembly = "hg38",
    min_nonzero_freq = 0, min_nonzero_mean = 0, min_overall_mean = 0,
    min_nonzero_freq_quantile = 0, min_nonzero_mean_quantile = 0, min_overall_mean_quantile = 0,
    distance_from_TSS = 10,
    excluded_cell_types_in_priors = "NotExclude",
    prior_min_TG_number = 0, prior_top_n_keep = 100000, prior_top_q_keep = 1
  )
})
message("Get ", length(TFs), " TFs, ", print_run_stats(monitor))
# Get 240 TFs, running: 0 hours 4 mins 15 secs, peak RAM: 0 GB 937 MB

qs::qsave(TFs, file = "./1.1 Data Preprocessing - Get TFs From Database/TFs.qsave")



##### 1.2 Annotating Chromatin Fragments and Identifying TGs #####
##### 1.2.1 Annotating Chromatin Fragments and Identifying TGs
setwd(Working_directory)
scRNAseq = qs::qread("./scRNAseq.qsave")
scATACseq = qs::qread("./scATACseq.qsave")
library(dplyr)
monitor = peakRAM::peakRAM({
  results_identify_TGs_all = identify_TGs(
    genome_for_anno = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene,
    fragment_annotation_range = 5e+05, scRNAseq = scRNAseq, scATACseq = scATACseq,
    annoDb = "org.Hs.eg.db", flankDistance = 10000
  )
})
message(print_run_stats(monitor)) # running: 0 hours 0 mins 55 secs, peak RAM: 1 GB 92 MB
qs::qsave(results_identify_TGs_all, file = "./1.2 Data Preprocessing - Identify TGs By Annotation/results_identify_TGs_all.qsave")


##### 1.2.2 Filtering Annotation Results
setwd(Working_directory)
results_identify_TGs_all = qs::qread("./1.2 Data Preprocessing - Identify TGs By Annotation/results_identify_TGs_all.qsave")
scRNAseq = qs::qread("./scRNAseq.qsave")
scATACseq = qs::qread("./scATACseq.qsave")
monitor = peakRAM::peakRAM({
  results_identify_TGs = filter_peak_annotations(
    results_identify_TGs = results_identify_TGs_all, scRNAseq = scRNAseq, scATACseq = scATACseq,
    forced_feature_calculation = TRUE,
    distance_points = c(1, 500, 500000), min_quantile_points = c(0.35, 0.55, 1.55),
    distance_quantile_fit_start = list(a = 0.006667, b = 1.659167, c = 0.6003),
    one_peak_at_least_per_gene = TRUE, gene_quality_factor_min = 0.9,
    weight_gene_nonzero_freq = 0.7, weight_gene_nonzero_mean = 0.2, weight_gene_overall_mean = 0.1,
    min_peak_nonzero_freq_quantile_scale = 0.4, min_peak_nonzero_mean_scale = 1.2, min_peak_overall_mean_quantile_scale = 0.1
  )
})
message(nrow(results_identify_TGs[["peak_anno_promoter_all"]]), "/", length(unique(results_identify_TGs[["peak_anno_promoter"]][["peak"]])), " peaks, ",
        length(unique(results_identify_TGs[["peak_anno_promoter_all"]][["SYMBOL"]])), "/", length(results_identify_TGs[["TGs"]]), " TGs, ",
        nrow(results_identify_TGs[["peak_anno_promoter_all"]]), "/", nrow(unique(results_identify_TGs[[1]][,c("peak","SYMBOL")])), " peak-TG, peaks/TGs = ",
        round(length(unique(results_identify_TGs[["peak_anno_promoter"]][["peak"]]))/length(results_identify_TGs[["TGs"]]), digits = 2), ", ",
        print_run_stats(monitor))
# 3374/2760 peaks, 2760/2894 TGs, 3374/2949 peak-TG, peaks/TGs = 0.95, running: 0 hours 0 mins 6 secs, peak RAM: 1 GB 76 MB

qs::qsave(results_identify_TGs, file = "./1.2 Data Preprocessing - Identify TGs By Annotation/results_identify_TGs.qsave")



##### 1.3 Generating Counts Matrices for TGs/TFs Expression and TGs Activity (Seurat Object) #####
### Run scDERI_apply-survival.R in folder ./survival.

### Load data.
setwd(Working_directory)
TGs = qs::qread("./1.2 Data Preprocessing - Identify TGs By Annotation/results_identify_TGs.qsave")$TGs
TFs = qs::qread("./1.1 Data Preprocessing - Get TFs From Database/TFs.qsave")
scRNAseq = qs::qread("./scRNAseq.qsave")
scATACseq = qs::qread("./scATACseq.qsave")

### Retain genes whose high expression also shortens survival.
survival_info_TF = qs::qread("./survival/survival_info_TF.qsave")
print(table(survival_info_TF$features)) # high: 28, low: 7
TFs = intersect(TFs, rownames(survival_info_TF)[survival_info_TF$features == "high"])
print(length(TFs)) # 28

survival_info_TG = qs::qread("./survival/survival_info_TG.qsave")
print(table(survival_info_TG$features)) # high: 402, low: 112
TGs = intersect(TGs, rownames(survival_info_TG)[survival_info_TG$features == "high"])
print(length(TGs)) # 402

# Run get_expression_and_activity_matrix
monitor = peakRAM::peakRAM({
  Basic_Info = get_expression_and_activity_matrix(
    TGs = TGs, TFs = TFs, scRNAseq = scRNAseq, scATACseq = scATACseq,
    Ensdb = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86,
    flanking_region = 10000, path = paste0(Working_directory, "raw_data/atac_fragments.tsv.gz")
  )
})
message(length(Basic_Info[["TFs"]]), " TFs, ", length(Basic_Info[["TGs"]]), " TGs, ",
        length(Basic_Info[["scRNAseq_for_GRN"]]@meta.data$seurat_annotations), " cells, ", print_run_stats(monitor))
# 28 TFs, 386 TGs, 9511 cells, running: 0 hours 20 mins 55 secs, peak RAM: 3 GB 752 MB

qs::qsave(Basic_Info, file = "./1.3 Data Preprocessing - Get Expression And Activity Matrix/Basic_Info.qsave")




########## Part 2: Inferring Pseudotime, Branching, and Grouping Cells ##########
##### 2.1 Pseudotime Analysis and Cell Branch Partitioning #####
##### 2.1.1 Obtaining Relevant Data of Cell Types of Interest
setwd(Working_directory)
monitor = peakRAM::peakRAM({
  interest_cell_type_data = get_interest_cell_type_data(
    interest_cell_type = "KP4", Basic_Info = Basic_Info, ncores = parallel::detectCores() - 1
  )
})
message(print_run_stats(monitor)) # running: 0 hours 0 mins 1 secs, peak RAM: 0 GB 854 MB
qs::qsave(interest_cell_type_data, file = "./2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment/interest_cell_type_data.qsave")


##### 2.1.2 Performing Pseudotemporal Ordering and Branching Analysis of Cell Types of Interest Based on Transcriptomic Data
setwd(Working_directory)
interest_cell_type_data = qs::qread("./2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment/interest_cell_type_data.qsave")
library(DDRTree)
monitor = peakRAM::peakRAM({
  interest_cell_type_Branches = order_pseudotime_and_divide_branches(
    interest_cell_type_data = interest_cell_type_data, auto_param_selection = TRUE
  )
  # 1→2→3→11→4→5→6→7→8→9, 1→2→3→11→4→5→6→7→8→10
})
message(print_run_stats(monitor)) # running: 0 hours 7 mins 13 secs, peak RAM: 6 GB 489 MB
qs::qsave(interest_cell_type_Branches, file = "./2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment/interest_cell_type_Branches.qsave")



##### 2.2 Cell Grouping and Group-Wise Information Specification #####
##### 2.2.1 Retrieving Pseudotemporal Expression/Activity Vectors for Each Gene
setwd(Working_directory)
interest_cell_type_Branches = qs::qread("./2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment/interest_cell_type_Branches.qsave")
interest_cell_type_data = qs::qread("./2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment/interest_cell_type_data.qsave")

sparsity_cut = list()
sparsity_cut[["KP4"]] = as.data.frame(matrix(
  0, nrow = 3, ncol = length(interest_cell_type_Branches[[1]]),
  dimnames = list(c("TFE", "TGA", "TGE"), paste0("branch", seq_along(interest_cell_type_Branches[[1]])))
))
sparsity_cut[[1]]["TFE", ] = 0.005
sparsity_cut[[1]]["TGA", ] = 0.005
sparsity_cut[[1]]["TGE", ] = 0.005
print(sparsity_cut)
# $KP4
#     branch1 branch2
# TFE   0.005   0.005
# TGA   0.005   0.005
# TGE   0.005   0.005

monitor = peakRAM::peakRAM({
  interest_cell_type_genes_pseudotime_info = get_genes_pseudotime_info(
    interest_cell_type_Branches = interest_cell_type_Branches,
    interest_cell_type_data = interest_cell_type_data,
    alpha_cell = 1, alpha_TF = 0, alpha_TG = 0,
    ncores = parallel::detectCores() - 1, sparsity_cut = sparsity_cut
  )
})
message(print_run_stats(monitor)) # running: 0 hours 0 mins 16 secs, peak RAM: 0 GB 472 MB

print(unlist(interest_cell_type_genes_pseudotime_info[[1]][["Branches_n_cell"]])) # 7242 8277
print(sapply(interest_cell_type_genes_pseudotime_info[[1]][["Branches_TFs"]], length)) # 23 23
print(sapply(interest_cell_type_genes_pseudotime_info[[1]][["Branches_TGs"]], length)) # 380 380

qs::qsave(interest_cell_type_genes_pseudotime_info, file = "./2.2 Data Processing - Cell Grouping/interest_cell_type_genes_pseudotime_info.qsave")


##### 2.2.2 Grouping Cells Based on Pseudotime Analysis and Calculating Valid (Non-Zero) Expression and Activity Values per Group
setwd(Working_directory)
interest_cell_type_genes_pseudotime_info = qs::qread("./2.2 Data Processing - Cell Grouping/interest_cell_type_genes_pseudotime_info.qsave")

n_cell = round(mean(unlist(interest_cell_type_genes_pseudotime_info[[1]][["Branches_n_cell"]])))
fit_scale_factor = 1.2
points_x_for_fitting_nc_and_nmin = c(
  round(n_cell * 0.25 * fit_scale_factor),
  round(n_cell * fit_scale_factor),
  round(n_cell * 1.75 * fit_scale_factor)
)
print(points_x_for_fitting_nc_and_nmin) # 2328 9312 16296

n_min_expected = 3
points_y_for_fitting_nc_and_nmin = c(n_min_expected - 1, n_min_expected, n_min_expected + 1)
print(points_y_for_fitting_nc_and_nmin) # 2 3 4

monitor = peakRAM::peakRAM({
  interest_cell_type_group = cell_grouping(
    interest_cell_type_genes_pseudotime_info = interest_cell_type_genes_pseudotime_info,
    points_x_for_fitting_nc_and_nmin = points_x_for_fitting_nc_and_nmin,
    points_y_for_fitting_nc_and_nmin = points_y_for_fitting_nc_and_nmin,
    ncores = parallel::detectCores() - 1
  )
})
message(print_run_stats(monitor)) # running: 0 hours 4 mins 42 secs, peak RAM: 0 GB 362 MB

print(unlist(interest_cell_type_group[[1]][["Branches_n_group"]])) # 8 8
print(as.numeric(interest_cell_type_group[[1]][["n_min"]])) # 3 3

qs::qsave(interest_cell_type_group, file = "./2.2 Data Processing - Cell Grouping/interest_cell_type_group.qsave")




########## Part 3: Constructing Initial GRNs from TFBS PWMs ##########
##### 3.1 Constructing an Initial GRN for Comprehensive TG-TF Pairs #####
##### 3.1.1  Constructing an Initial GRN for Comprehensive TG-TF Pairs Based on Base PWM of TF Binding Sites
setwd(Working_directory)
interest_cell_type_genes_pseudotime_info = qs::qread("./2.2 Data Processing - Cell Grouping/interest_cell_type_genes_pseudotime_info.qsave")
results_identify_TGs = qs::qread("./1.2 Data Preprocessing - Identify TGs By Annotation/results_identify_TGs.qsave")
library(dplyr)
monitor = peakRAM::peakRAM({
  interest_cell_type_iGRN_all = get_iGRN_by_TFBS_pwm_by_JASPAR2024(
    interest_cell_type_genes_pseudotime_info = interest_cell_type_genes_pseudotime_info,
    results_identify_TGs = results_identify_TGs, genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
    min_score_for_matchPWM = "0%", species = "Homo sapiens",
    collection = c("CORE", "CNE", "PHYLOFACTS", "SPLICE", "POLII", "FAM", "PBM", "PBM_HOMEO", "PBM_HLH", "UNVALIDATED"),
    output_predicted_TFBS = FALSE, ncores = parallel::detectCores() - 1
  )
})
message(ncol(interest_cell_type_iGRN_all[[1]]), " TFs, ", nrow(interest_cell_type_iGRN_all[[1]]), " TGs, ",
        sum(interest_cell_type_iGRN_all[[1]] > 0), " pairs, ", print_run_stats(monitor))
# 23 TFs, 380 TGs, 8740 pairs, running: 0 hours 2 mins 10 secs, peak RAM: 0 GB 27 MB

qs::qsave(interest_cell_type_iGRN_all, file = "./3 get iGRN/interest_cell_type_iGRN_all.qsave")


##### 3.1.2 Filtering  the Initial GRN by ChIP-Atlas Priors
setwd(Working_directory)
interest_cell_type_iGRN_all = qs::qread("./3 get iGRN/interest_cell_type_iGRN_all.qsave")

results_identify_TGs = qs::qread("./1.2 Data Preprocessing - Identify TGs By Annotation/results_identify_TGs.qsave")
peak_anno_promoter = results_identify_TGs[["peak_anno_promoter"]]; rm(results_identify_TGs)
peak_anno_promoter = peak_anno_promoter[peak_anno_promoter$SYMBOL %in% rownames(interest_cell_type_iGRN_all[[1]]), ]
peak_TG_ratio = length(unique(peak_anno_promoter[["peak"]])) / length(unique(peak_anno_promoter[["SYMBOL"]]))
message("peak_TG_ratio = ", round(peak_TG_ratio, digits = 4)) # peak_TG_ratio = 1.0237

theta_i_top_quantile = 1
theta_i_top_n = nrow(interest_cell_type_iGRN_all[[1]])

monitor = peakRAM::peakRAM({
  interest_cell_type_iGRN = filter_iGRN_by_ChIP_Atlas_priors(
    interest_cell_type_iGRN = interest_cell_type_iGRN_all,
    theta_i_top_quantile = theta_i_top_quantile, theta_i_top_n = theta_i_top_n,
    TG_rescue_number = 5
  )
})
message("top_q=", round(theta_i_top_quantile, digits = 4),
        ", top_n=", theta_i_top_n, ", ",
        ncol(interest_cell_type_iGRN[[1]]), " TFs, ", nrow(interest_cell_type_iGRN[[1]]), " TGs, ",
        sum(interest_cell_type_iGRN[[1]] > 0), " pairs, ",
        round(sum(interest_cell_type_iGRN[[1]] > 0) / ncol(interest_cell_type_iGRN[[1]]), digits = 2), " TGs per TF, ",
        print_run_stats(monitor))
# top_q=1, top_n=380, 23 TFs, 368 TGs, 7327 pairs, 318.57 TGs per TF, running: 0 hours 0 mins 3 secs, peak RAM: 0 GB 69 MB

qs::qsave(interest_cell_type_iGRN, file = "./3 get iGRN/interest_cell_type_iGRN.qsave")



##### 3.2 Setting the Regulatory Threshold tau #####
setwd(Working_directory)
interest_cell_type_iGRN = qs::qread("./3 get iGRN/interest_cell_type_iGRN.qsave")
monitor = peakRAM::peakRAM({
  interest_cell_type_tau = get_cell_type_tau(
    interest_cell_type_iGRN = interest_cell_type_iGRN
  )
})
message("tau is ", round(interest_cell_type_tau[1], digits = 4), ", ", print_run_stats(monitor))
# tau is 0.0953, running: 0 hours 0 mins 0 secs, peak RAM: 0 GB 0 MB

qs::qsave(interest_cell_type_tau, file = "./3 get iGRN/interest_cell_type_tau.qsave")



##### 3.3 Constructing the iGRN for Each Branch #####
setwd(Working_directory)
interest_cell_type_iGRN = qs::qread("./3 get iGRN/interest_cell_type_iGRN.qsave")
interest_cell_type_group = qs::qread("./2.2 Data Processing - Cell Grouping/interest_cell_type_group.qsave")
monitor = peakRAM::peakRAM({
  interest_cell_type_branch_iGRN = get_branch_iGRN(
    interest_cell_type_iGRN = interest_cell_type_iGRN,
    interest_cell_type_group = interest_cell_type_group,
    ncores = parallel::detectCores() - 1
  )
})
message(paste(unlist(lapply(interest_cell_type_branch_iGRN[[1]], ncol)), collapse = "/"), " TFs, ",
        paste(unlist(lapply(interest_cell_type_branch_iGRN[[1]], nrow)), collapse = "/"), " TGs, ",
        paste(sapply(interest_cell_type_branch_iGRN[[1]], function(df) sum(df > 0, na.rm = TRUE)), collapse = "/"), " pairs, ", print_run_stats(monitor))
# 23/23 TFs, 368/368 TGs, 7327/7327 pairs, running: 0 hours 0 mins 3 secs, peak RAM: 0 GB 168 MB

qs::qsave(interest_cell_type_branch_iGRN, file = "./3 get iGRN/interest_cell_type_branch_iGRN.qsave")




########## Part 4: Developing DERI-Based GRN Regression Predictive Models for Each Branch ##########
##### 4.1 Partitioning Training Sets #####
setwd(Working_directory)
interest_cell_type_branch_iGRN = qs::qread("./3 get iGRN/interest_cell_type_branch_iGRN.qsave")
interest_cell_type_tau = qs::qread("./3 get iGRN/interest_cell_type_tau.qsave")
interest_cell_type_group = qs::qread("./2.2 Data Processing - Cell Grouping/interest_cell_type_group.qsave")
monitor = peakRAM::peakRAM({
  interest_cell_type_branch_training_set = select_training_set(
    interest_cell_type_branch_iGRN = interest_cell_type_branch_iGRN,
    interest_cell_type_tau  = interest_cell_type_tau, interest_cell_type_group = interest_cell_type_group,
    top_n_for_max_theta_in_each_regulon = 2, max_theta_retention_rate_in_tops = 1,
    max_theta_retention_rate_in_all = 1, retention_number_of_tops = 999999,
    refine = TRUE, keep_all_TF = TRUE, top_ratio_in_refine = 1,
    ncores = parallel::detectCores() - 1
  )
})
message(paste(lapply(interest_cell_type_branch_training_set[[1]], nrow), collapse = ","), " training TF-TG pairs, ", print_run_stats(monitor))
# 85,85,85 training TF-TG pairs, running: 0 hours 0 mins 1 secs, peak RAM: 0 GB 134 MB

qs::qsave(interest_cell_type_branch_training_set, file = "./4.1 Build Prediction Model - Select Training Set/interest_cell_type_branch_training_set.qsave")



##### 4.2 Build Predictive Model for Each Branch #####
##### 4.2.1 Setting Initial Values, Lower Bounds, and Upper Bounds for Parameters
setwd(Working_directory)
interest_cell_type_branch_training_set = qs::qread("./4.1 Build Prediction Model - Select Training Set/interest_cell_type_branch_training_set.qsave")
monitor = peakRAM::peakRAM({
  interest_cell_type_branch_init_params = set_init_params(
    interest_cell_type_branch_training_set = interest_cell_type_branch_training_set
  )
})
message(print_run_stats(monitor)) # running: 0 hours 0 mins 3 secs, peak RAM: 0 GB 15 MB
qs::qsave(interest_cell_type_branch_init_params, file = "./4.2 BUild Prediction Model/interest_cell_type_branch_init_params.qsave")


##### 4.2.2 Training models
setwd(Working_directory)
interest_cell_type_branch_training_set = qs::qread("./4.1 Build Prediction Model - Select Training Set/interest_cell_type_branch_training_set.qsave")
interest_cell_type_branch_init_params = qs::qread("./4.2 BUild Prediction Model/interest_cell_type_branch_init_params.qsave")
interest_cell_type_group = qs::qread("./2.2 Data Processing - Cell Grouping/interest_cell_type_group.qsave")
interest_cell_type_tau = qs::qread("./3 get iGRN/interest_cell_type_tau.qsave")
interest_cell_type_genes_pseudotime_info = qs::qread("./2.2 Data Processing - Cell Grouping/interest_cell_type_genes_pseudotime_info.qsave")
monitor = peakRAM::peakRAM({
  interest_cell_type_branch_model_train_all = model_train(
    interest_cell_type_branch_training_set = interest_cell_type_branch_training_set,
    interest_cell_type_branch_init_params = interest_cell_type_branch_init_params,
    interest_cell_type_group = interest_cell_type_group, interest_cell_type_tau = interest_cell_type_tau,
    interest_cell_type_genes_pseudotime_info = interest_cell_type_genes_pseudotime_info,
    max_epochs = 100, early_stop = 5, eps_theta = 1e-2, eps_loss = 1e-5, parallel_ga = FALSE, seed_ga = 123,
    popSize_ga = 50, maxiter_ga = 150, pcrossover_ga = 0.5, pmutation_ga = 0.5,
    iterations_grad = 30, ncores_grad = 1,
    ncores = parallel::detectCores() - 1
  )
})
message(print_run_stats(monitor)) # running: 0 hours 46 mins 26 secs, peak RAM: 0 GB 12 MB

qs::qsave(interest_cell_type_branch_model_train_all, file = "./4.2 BUild Prediction Model/interest_cell_type_branch_model_train_all.qsave")

interest_cell_type_branch_model_train = interest_cell_type_branch_model_train_all
for (cell_type in names(interest_cell_type_branch_model_train)) {
  df = interest_cell_type_branch_model_train[[cell_type]][["best_pred_result"]]
  bad_flag = with(df, theta_i < 0.011)
  rownames_filtered = rownames(df)[bad_flag]
  interest_cell_type_branch_model_train[[cell_type]][["best_pred_result"]] = df[!bad_flag, ]
  for (n in 1:(length(interest_cell_type_branch_model_train[[cell_type]]) - 1)) {
    lst = interest_cell_type_branch_model_train[[cell_type]][[paste0("branch", n)]]
    interest_cell_type_branch_model_train[[cell_type]][[paste0("branch", n)]] = lst[!names(lst) %in% rownames_filtered]
  }
}
message(nrow(interest_cell_type_branch_model_train_all[[1]][["best_pred_result"]]), "/",
        nrow(interest_cell_type_branch_model_train[[1]][["best_pred_result"]]), " training results across branches, ",
        paste(unlist(lapply(interest_cell_type_branch_model_train_all[[1]][2:length(interest_cell_type_branch_model_train_all[[1]])], length)), collapse = ","), "/",
        paste(unlist(lapply(interest_cell_type_branch_model_train[[1]][2:length(interest_cell_type_branch_model_train[[1]])], length)), collapse = ","), " training results per branches")
# 85/85 training results across branches, 85,85/85,85 training results per branches

qs::qsave(interest_cell_type_branch_model_train, file = "./4.2 BUild Prediction Model/interest_cell_type_branch_model_train.qsave")




########## Part 5: Predicting GRN ##########
##### 5.1 Getting Data for Inferring GRN #####
setwd(Working_directory)
interest_cell_type_branch_iGRN = qs::qread("./3 get iGRN/interest_cell_type_branch_iGRN.qsave")
interest_cell_type_group = qs::qread("./2.2 Data Processing - Cell Grouping/interest_cell_type_group.qsave")
ncores = parallel::detectCores() - 1
n_part = ncores
monitor = peakRAM::peakRAM({
  interest_cell_type_branch_data_for_GRN_infer = get_data_for_inferring_GRN(
    interest_cell_type_branch_iGRN = interest_cell_type_branch_iGRN,
    interest_cell_type_group = interest_cell_type_group,
    n_part = n_part, ncores = ncores
  )
})
message(print_run_stats(monitor)) # running: 0 hours 11 mins 34 secs, peak RAM: 0 GB 514 MB
# running: 0 hours 11 mins 17 secs, peak RAM: 0 GB 686 MB
sapply(interest_cell_type_branch_data_for_GRN_infer[[1]], dim)
#      branch1 branch2  GRN
# [1,]    7327    7327 7327
# [2,]      27      27    3
qs::qsave(interest_cell_type_branch_data_for_GRN_infer, file = "./5 Infer GRN/interest_cell_type_branch_data_for_GRN_infer.qsave")



##### 5.2 Inferring GRN for Each Branch #####
setwd(Working_directory)
interest_cell_type_branch_data_for_GRN_infer = qs::qread("./5 Infer GRN/interest_cell_type_branch_data_for_GRN_infer.qsave")
interest_cell_type_tau = qs::qread("./3 get iGRN/interest_cell_type_tau.qsave")
interest_cell_type_branch_model_train = qs::qread("./4.2 BUild Prediction Model/interest_cell_type_branch_model_train.qsave")
interest_cell_type_group = qs::qread("./2.2 Data Processing - Cell Grouping/interest_cell_type_group.qsave")
ncores = parallel::detectCores() - 1
interest_cell_type_pGRN = list()
cell_type = names(interest_cell_type_branch_data_for_GRN_infer)[1] # KP4
interest_cell_type_pGRN[[cell_type]] = list()
for (n in 1:(length(interest_cell_type_branch_data_for_GRN_infer[[cell_type]]) - 1)) {
  message("Infer cell type ", cell_type, " branch ", n, ".")
  each = ceiling(nrow(interest_cell_type_branch_data_for_GRN_infer[[1]][[paste0("branch", n)]]) / ncores)
  print(each) # 94 94
  gc()

  # Run model_train with runtime and peak memory monitoring.
  monitor = peakRAM::peakRAM({
    interest_cell_type_pGRN[[cell_type]][[paste0("branch", n)]] = infer_GRN(
      interest_cell_type_branch_data_for_GRN_infer = interest_cell_type_branch_data_for_GRN_infer,
      interest_cell_type_tau = interest_cell_type_tau,
      interest_cell_type_branch_model_train = interest_cell_type_branch_model_train,
      interest_cell_type_group = interest_cell_type_group,
      select_cell_type = cell_type, select_branch = n,
      each = each, ncores = ncores
    )[[cell_type]][[paste0("branch", n)]]
  })
  message("Cell type ", cell_type, " - Branch ", n, ": ",print_run_stats(monitor))
  # Cell type KP4 - Branch 1: running: 0 hours 7 mins 37 secs, peak RAM: 0 GB 200 MB
  # Cell type KP4 - Branch 2: running: 0 hours 5 mins 29 secs, peak RAM: 0 GB 200 MB

  # Check point
  save.image(paste0("./5 Infer GRN/5_Infer_GRN_branch", n, "middle.RData"))
  # load(paste0("./5 Infer GRN/5_Infer_GRN_branch", n, "middle.RData"))

  gc()
}
qs::qsave(interest_cell_type_pGRN, file = "./5 Infer GRN/interest_cell_type_pGRN.qsave")



##### 5.3 Selecting Best Prediction for each TF-TG #####
setwd(Working_directory)
interest_cell_type_pGRN = qs::qread("./5 Infer GRN/interest_cell_type_pGRN.qsave")
interest_cell_type_tau = qs::qread("./3 get iGRN/interest_cell_type_tau.qsave")
interest_cell_type_genes_pseudotime_info = qs::qread("./2.2 Data Processing - Cell Grouping/interest_cell_type_genes_pseudotime_info.qsave")

ncores = parallel::detectCores() - 1
interest_cell_type_branch_data_for_GRN_infer = qs::qread("./5 Infer GRN/interest_cell_type_branch_data_for_GRN_infer.qsave")
each = ceiling((max(sapply(
  interest_cell_type_branch_data_for_GRN_infer[[1]][1:(length(interest_cell_type_branch_data_for_GRN_infer[[1]]) - 1)],
  nrow
))) / ncores)
print(each) # 94

monitor = peakRAM::peakRAM({
  interest_cell_type_pGRN = select_best_prediction(
    interest_cell_type_pGRN = interest_cell_type_pGRN, interest_cell_type_tau = interest_cell_type_tau,
    interest_cell_type_genes_pseudotime_info = interest_cell_type_genes_pseudotime_info,
    fit_loss_quantile_threshold = 1, theta_p_quantile_threshold = 0, each = each, ncores = ncores
  )
})
message(nrow(interest_cell_type_pGRN[[1]][["pred_summary"]]), " pairs remaining, ", print_run_stats(monitor))
# 7327 pairs remaining, running: 0 hours 0 mins 42 secs, peak RAM: 0 GB 80 MB

qs::qsave(interest_cell_type_pGRN, file = "./5 Infer GRN/interest_cell_type_pGRN.qsave")



##### 5.4 Filtering GRN #####
interest_cell_type_pGRN = qs::qread("./5 Infer GRN/interest_cell_type_pGRN.qsave")
fit_loss_quantile_threshold = 1
theta_p_quantile_threshold = 0
fit_loss_upper = 0.2
theta_p_lower = qs::qread("./3 get iGRN/interest_cell_type_tau.qsave")[1]
monitor = peakRAM::peakRAM({
  interest_cell_type_pGRN = filter_GRN(
    interest_cell_type_pGRN = interest_cell_type_pGRN,
    fit_loss_quantile_threshold = fit_loss_quantile_threshold,
    theta_p_quantile_threshold = theta_p_quantile_threshold,
    fit_loss_upper = fit_loss_upper,
    theta_p_lower = theta_p_lower
  )
})
message(nrow(interest_cell_type_pGRN[[1]][["pred_summary_loose"]]), " loose pairs, ", print_run_stats(monitor))
# 6715 loose pairs, running: 0 hours 0 mins 0 secs, peak RAM: 0 GB 4 MB

qs::qsave(interest_cell_type_pGRN, file = "./5 Infer GRN/interest_cell_type_pGRN.qsave")



##### 5.5 Inferring eGRN #####
setwd(Working_directory)
scRNAseq = qs::qread(file = "./scRNAseq.qsave")
scATACseq = qs::qread("./scATACseq.qsave")
results_identify_TGs = qs::qread("./1.2 Data Preprocessing - Identify TGs By Annotation/results_identify_TGs.qsave")
interest_cell_type_iGRN = qs::qread("./3 get iGRN/interest_cell_type_iGRN.qsave")
interest_cell_type_pGRN = qs::qread("./5 Infer GRN/interest_cell_type_pGRN.qsave")
monitor = peakRAM::peakRAM({
  interest_cell_type_peGRN = infer_eGRN(
    scRNAseq = scRNAseq, scATACseq = scATACseq, results_identify_TGs = results_identify_TGs,
    interest_cell_type_iGRN = interest_cell_type_iGRN, interest_cell_type_pGRN = interest_cell_type_pGRN,
  )
})
message(print_run_stats(monitor)) # running: 0 hours 0 mins 21 secs, peak RAM: 2 GB 324 MB
qs::qsave(interest_cell_type_peGRN, file = "./5 Infer GRN/interest_cell_type_peGRN.qsave")



##### 5.6 Filtering peGRN #####
setwd(Working_directory)
interest_cell_type_peGRN = qs::qread("./5 Infer GRN/interest_cell_type_peGRN.qsave")

# keep only the top 50% of TGs per TF
min_peakTG_per_TF = round(length(unique(interest_cell_type_peGRN[["KP4"]][["loose"]][["TG"]])) * 0.5)

monitor = peakRAM::peakRAM({
  interest_cell_type_peGRN = filter_eGRN(
    interest_cell_type_peGRN = interest_cell_type_peGRN,
    theta_p_points = c(0, 0.5, 1),
    initial_quantile_points = c(5.79, 4.66, 4.23),
    theta_p_quantile_fit_start = list(a = 0.646435, b = 0.397096, c = 3.271714),
    rescue_peakTG = FALSE, peakTG_rescue_quantile_threshold = 0.00, min_TF_per_peakTG = 1,
    rescue_TF = TRUE, TF_rescue_quantile_threshold = 0.00, min_peakTG_per_TF = min_peakTG_per_TF,
    factor_balance = list(`FG-NL` = 0.2, `FG-theta_i-R` = 0.2,
                          `F-NZER-R` = 0.2, `F-NZEA-R` = 0.2, `F-EA-R` = 0.2,
                          `P-NZAR-R` = 0.2, `P-NZAA-R` = 0.2, `P-AA-R` = 0.2,
                          `G-NZER-R` = 0.2, `G-NZEA-R` = 0.2, `G-EA-R` = 0.2),
    factor_range = list(`FG-NL` = 0.43, `FG-theta_i-R` = 0.64,
                        `F-NZER-R` = 0.28, `F-NZEA-R` = 0.33, `F-EA-R` = 0.22,
                        `P-NZAR-R` = 0.19, `P-NZAA-R` = 0.67, `P-AA-R` = 0.13,
                        `G-NZER-R` = 0.21, `G-NZEA-R` = 0.12, `G-EA-R` = 0.26)
  )
})
message(print_run_stats(monitor), "\n",
        "[loose] ",
        length(unique(interest_cell_type_peGRN[[1]][["loose"]][["TF"]])), "/", length(unique(interest_cell_type_peGRN[[1]][["loose_filtered_simplified"]][["TF"]])), " TFs, ",
        length(unique(interest_cell_type_peGRN[[1]][["loose"]][["peak"]])), "/", length(unique(interest_cell_type_peGRN[[1]][["loose_filtered_simplified"]][["peak"]])), " peaks, ",
        length(unique(interest_cell_type_peGRN[[1]][["loose"]][["TG"]])), "/", length(unique(interest_cell_type_peGRN[[1]][["loose_filtered_simplified"]][["TG"]])), " TGs, ",
        nrow(unique(interest_cell_type_peGRN[[1]][["loose"]][, c("TF", "peak")])), "/", nrow(unique(interest_cell_type_peGRN[[1]][["loose_filtered_simplified"]][, c("TF", "peak")])), " TF-peak, ",
        nrow(unique(interest_cell_type_peGRN[[1]][["loose"]][, c("peak", "TG")])), "/", nrow(unique(interest_cell_type_peGRN[[1]][["loose_filtered_simplified"]][, c("peak", "TG")])), " peak-TG, ",
        nrow(unique(interest_cell_type_peGRN[[1]][["loose"]][, c("TF", "TG")])), "/", nrow(unique(interest_cell_type_peGRN[[1]][["loose_filtered_simplified"]][, c("TF", "TG")])), " TF-TG, ",
        nrow(unique(interest_cell_type_peGRN[[1]][["loose"]][, c("TF", "peak", "TG")])), "/", nrow(unique(interest_cell_type_peGRN[[1]][["loose_filtered_simplified"]][, c("TF", "peak", "TG")])), " TF-peak-TG, rescue ratio ",
        round(1 - sum(interest_cell_type_peGRN[[1]][["loose_filtered"]][["filter_keep"]]) / nrow(interest_cell_type_peGRN[[1]][["loose_filtered_simplified"]]), 3), "\n",
        "[strict] ",
        length(unique(interest_cell_type_peGRN[[1]][["strict"]][["TF"]])), "/", length(unique(interest_cell_type_peGRN[[1]][["strict_filtered_simplified"]][["TF"]])), " TFs, ",
        length(unique(interest_cell_type_peGRN[[1]][["strict"]][["peak"]])), "/", length(unique(interest_cell_type_peGRN[[1]][["strict_filtered_simplified"]][["peak"]])), " peaks, ",
        length(unique(interest_cell_type_peGRN[[1]][["strict"]][["TG"]])), "/", length(unique(interest_cell_type_peGRN[[1]][["strict_filtered_simplified"]][["TG"]])), " TGs, ",
        nrow(unique(interest_cell_type_peGRN[[1]][["strict"]][, c("TF", "peak")])), "/", nrow(unique(interest_cell_type_peGRN[[1]][["strict_filtered_simplified"]][, c("TF", "peak")])), " TF-peak, ",
        nrow(unique(interest_cell_type_peGRN[[1]][["strict"]][, c("peak", "TG")])), "/", nrow(unique(interest_cell_type_peGRN[[1]][["strict_filtered_simplified"]][, c("peak", "TG")])), " peak-TG, ",
        nrow(unique(interest_cell_type_peGRN[[1]][["strict"]][, c("TF", "TG")])), "/", nrow(unique(interest_cell_type_peGRN[[1]][["strict_filtered_simplified"]][, c("TF", "TG")])), " TF-TG, ",
        nrow(unique(interest_cell_type_peGRN[[1]][["strict"]][, c("TF", "peak", "TG")])), "/", nrow(unique(interest_cell_type_peGRN[[1]][["strict_filtered_simplified"]][, c("TF", "peak", "TG")])), " TF-peak-TG, rescue ratio ",
        round(1 - sum(interest_cell_type_peGRN[[1]][["strict_filtered"]][["filter_keep"]]) / nrow(interest_cell_type_peGRN[[1]][["strict_filtered_simplified"]]), 3))
# running: 0 hours 0 mins 3 secs, peak RAM: 0 GB 35 MB
# [loose] 23/23 TFs, 375/375 peaks, 367/367 TGs, 6885/3732 TF-peak, 377/377 peak-TG, 6715/3690 TF-TG, 6912/3744 TF-peak-TG, rescue ratio 0.908
# [strict] 23/23 TFs, 375/375 peaks, 367/367 TGs, 6881/3733 TF-peak, 377/377 peak-TG, 6711/3691 TF-TG, 6908/3745 TF-peak-TG, rescue ratio 0.908

qs::qsave(interest_cell_type_peGRN, file = "./5 Infer GRN/interest_cell_type_peGRN.qsave")



##### 5.7 Generate disease-associated TF-peak-TG network #####
setwd(Working_directory)

# Initialize disease-associated eGRN
interest_cell_type_peGRN = qs::qread("./5 Infer GRN/interest_cell_type_peGRN.qsave")
disease_associated_eGRN = interest_cell_type_peGRN[["KP4"]][["loose_filtered"]]
rm(interest_cell_type_peGRN)

# Get TFs, peaks, and TGs in disease-associated eGRN
eGRN = disease_associated_eGRN
disease_associated_eGRN = list()
disease_associated_eGRN[["TFs"]] = sort(unique(eGRN$TF))
disease_associated_eGRN[["peaks"]] = sort(unique(eGRN$peak))
disease_associated_eGRN[["TGs"]] = sort(unique(eGRN$TG))
disease_associated_eGRN[["TF_peaks_TGs"]] = eGRN

qs::qsave(disease_associated_eGRN, file = "./5 Infer GRN/disease_associated_eGRN.qsave")
# writexl::write_xlsx(disease_associated_eGRN[["TF_peaks_TGs"]], "./5 Infer GRN/disease_associated_eGRN.xlsx")




########## Part 6: Plot for Experiment ##########
##### 6.1 Peak-TG Distance Distribution Histogram #####
setwd(Working_directory)
library(dplyr)
library(ggplot2)

# Load promoter annotation data
results_identify_TGs = qs::qread("./1.2 Data Preprocessing - Identify TGs By Annotation/results_identify_TGs.qsave")
peak_anno_promoter = results_identify_TGs[["peak_anno_promoter"]]
peak_anno_promoter = peak_anno_promoter[, c("peak", "SYMBOL", "distanceToTSS")]
peak_anno_promoter[, "distanceToTSS"] = abs(peak_anno_promoter[, "distanceToTSS"])
peak_anno_promoter[, "SYMBOL"] = toupper(peak_anno_promoter[, "SYMBOL"])
dim(peak_anno_promoter) # [1] 2949    3

# Retain only results where distanceToTSS falls between the lower and upper quartiles
# distanceToTSS_q1 = quantile(peak_anno_promoter[, "distanceToTSS"], probs = 0.25)
# distanceToTSS_q3 = quantile(peak_anno_promoter[, "distanceToTSS"], probs = 0.75)

# Restrict peak-TG distance to (0, 15 kb)
distanceToTSS_q1 = 0
distanceToTSS_q3 = 15000
peak_anno_promoter = peak_anno_promoter[
  peak_anno_promoter$distanceToTSS > distanceToTSS_q1 & peak_anno_promoter$distanceToTSS < distanceToTSS_q3, , drop = FALSE
]

# Histogram parameters
n = length(peak_anno_promoter$distanceToTSS) # number of data points
binwidth = 500 # Histogram bin width

# Create data frame for plotting
df = data.frame(value = peak_anno_promoter$distanceToTSS)

# Histogram overlaid with a scaled density curve
All_Peak_TG_Distance_Distribution_Histogram = ggplot(df, aes(x = value)) +
  geom_histogram( # count-based histogram
    aes(y = after_stat(count)),
    binwidth = binwidth,
    boundary = 0,
    fill = "lightblue",
    color = "black",
    alpha = 0.7
  ) +
  geom_density( # density scaled to histogram counts
    aes(y = after_stat(density) * n * binwidth, color = "Density Curve"),
    linewidth = 1,
    adjust = 1,
    bounds = c(0, 15000),
    show.legend = FALSE
  ) +
  scale_y_continuous(name = "Count", breaks = scales::extended_breaks(n = 5)) +
  scale_x_continuous(breaks = seq(0, 13500, by = 4500)) +
  labs(x = "TG–peak regulatory distance", y = "Count") +
  theme_classic() +
  theme(
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
print(All_Peak_TG_Distance_Distribution_Histogram)

ggsave("./6 Plot for Experiment/6.1_All_Peak_TG_Distance_Distribution_Histogram.png",
       plot = All_Peak_TG_Distance_Distribution_Histogram,
       width = 3, height = 2, dpi = 300)

### Export Histogram Bin Counts
# Define 500-bp bins over [0, 15 kb]
bin_breaks = seq(distanceToTSS_q1, distanceToTSS_q3, by = binwidth)

# Build interval labels: [0,500], (500,1000], ...
bin_labels = c(
  sprintf("[%d,%d]", bin_breaks[1], bin_breaks[2]),
  sprintf("(%d,%d]", bin_breaks[2:(length(bin_breaks) - 1)], bin_breaks[3:length(bin_breaks)])
)

# Assign distances to bins, retaining all bins including zero counts
bin_factor = cut(df$value, breaks = bin_breaks, labels = bin_labels,
                  include.lowest = TRUE, right = TRUE)
hist_table = data.frame(
  "Distance interval" = bin_labels,
  "TG-peak regulation count" = as.integer(table(factor(bin_factor, levels = bin_labels))),
  check.names = FALSE
)

# Save to the same folder as the figure
openxlsx::write.xlsx(
  hist_table,
  "./6 Plot for Experiment/6.1_All_Peak_TG_Distance_Distribution_Histogram.xlsx",
  rowNames = FALSE
)



##### 6.2 Fitted Function Plot of Minimum Cell Count per Group vs. Branch Cell Number #####
setwd(Working_directory)
library(minpack.lm)
library(ggplot2)
library(ggrepel)
interest_cell_type_genes_pseudotime_info = qs::qread("./2.2 Data Processing - Cell Grouping/interest_cell_type_genes_pseudotime_info.qsave")
interest_cell_type_group = qs::qread("./2.2 Data Processing - Cell Grouping/interest_cell_type_group.qsave")

# Anchor points for fitting
points_x_for_fitting_nc_and_nmin = c(2328, 9312, 16296)
points_y_for_fitting_nc_and_nmin = c(2, 3, 4)

# Organizing data point coordinates
data = data.frame(x = points_x_for_fitting_nc_and_nmin, y = points_y_for_fitting_nc_and_nmin)

# Initial parameter values
a_init = 1
b_init = 1
c_init = stats::median(data$y)

# Fit y = c * log(a*x + b) with the Levenberg-Marquardt algorithm
fit = minpack.lm::nlsLM(
  y ~ c * log(a * x + b),
  data = data,
  start = list(a = a_init, b = b_init, c = c_init),
  control = minpack.lm::nls.lm.control(maxiter = 1000)
)

# Extract fitted parameters
a = stats::coef(fit)["a"]
b = stats::coef(fit)["b"]
c = stats::coef(fit)["c"]
options(scipen = 999)
message("Fitting results: a = ", a, ", b = ", b,", c = ", c, ".")
# Fitting results: a = 0.0000305752128550459, b = 1.03206672515629, c = 94.019683601527.

# Branch-level annotation points
x_annotate = unlist(interest_cell_type_genes_pseudotime_info[["KP4"]][["Branches_n_cell"]])
y_annotate = interest_cell_type_group[["KP4"]][["n_min"]]
label_annotate = sprintf("(%.0f, %.0f)", ceiling(x_annotate), ceiling(y_annotate))
label_annotate = paste0("Branch ", seq_along(label_annotate), ": ", label_annotate)
annotate_df = data.frame(
  x = x_annotate,
  y = y_annotate,
  label = label_annotate
)

# Fitted curve
x = seq(100, max(x_annotate) * 1.2, by = 1)
y = ceiling(c * log(a * x + b))
df = data.frame(x = x, y = y)

# Build plot
Fitted_Function_Plot_of_Minimum_Cell_Count_per_Group_vs._Branch_Cell_Number = ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
  ggplot2::geom_line(color = "blue", linewidth = 1) +
  ggplot2::geom_hline(yintercept = 0, color = "gray50") +
  ggplot2::geom_vline(xintercept = 0, color = "gray50") +
  ggplot2::geom_point(data = annotate_df, ggplot2::aes(x = x, y = y), color = "red", fill = NA, shape = 21, size = 2, stroke = 1.5) +
  ggrepel::geom_text_repel(
    data = annotate_df,
    ggplot2::aes(label = label),
    nudge_x = 500,
    nudge_y = -0.2 * max(df$y),
    force = 60,
    force_pull = 0,
    size = 3,
    segment.color = "gray60",
    box.padding = 0.5,
    max.overlaps = Inf
  ) + ggplot2::labs(
    title = "y = [c*log(a*x+b)]",
    subtitle = paste0("a = ", format(round(a, 5), nsmall = 4), ", b = ", format(round(b, 5), nsmall = 4), ", c = ", format(round(c, 5), nsmall = 4)),
    x = "Branch Cell Number",
    y = "Minimum Cell Count per Group When Grouping"
  ) + ggplot2::scale_x_continuous(labels = function(x) format(x, big.mark = ",")) +
  ggplot2::scale_y_continuous(labels = function(y) format(round(y, 1), nsmall = 1)) +
  ggplot2::coord_cartesian(xlim = c(0, max(x)), ylim = c(0, max(y))) +
  ggplot2::theme_minimal(base_family = "sans") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 12),
    plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 8, color = "gray40"),
    axis.text = ggplot2::element_text(size = 8),
    panel.grid.minor = ggplot2::element_blank()
  )
print(Fitted_Function_Plot_of_Minimum_Cell_Count_per_Group_vs._Branch_Cell_Number)

# Save plot
ggplot2::ggsave("./6 Plot for Experiment/6.2_Fitted_Function_Plot_of_Minimum_Cell_Count_per_Group_vs._Branch_Cell_Number.png",
                plot = Fitted_Function_Plot_of_Minimum_Cell_Count_per_Group_vs._Branch_Cell_Number,
                width = 3.5, height = 3.9, dpi = 300)



##### 6.3 Bar Chart of Branches Information #####
setwd(Working_directory)
library(Cairo)
interest_cell_type_group = qs::qread("./2.2 Data Processing - Cell Grouping/interest_cell_type_group.qsave")
interest_cell_type_genes_pseudotime_info = qs::qread("./2.2 Data Processing - Cell Grouping/interest_cell_type_genes_pseudotime_info.qsave")

branches = c("1→2→3→11→4→5→6→7→8→9", "1→2→3→11→4→5→6→7→8→10") # branch state sequences
n_bar = length(interest_cell_type_group[["KP4"]][["n_min"]]); n_bar # number of branches
x_length = max(unlist(interest_cell_type_group[["KP4"]][["Branches_T_k"]])) + 0.5; x_length # x-axis upper limit
times = 1.5 # plot scaling factor
colors = rainbow(max(unlist(interest_cell_type_group[["KP4"]][["Branches_n_group"]])), s = 1, v = 1, alpha = 0.3); colors # group fill colors

# Open PNG device
CairoPNG(
  filename = paste0("./6 Plot for Experiment/6.3_Branches_information_BarChart.png"),
  width = 0.8 * (x_length + 1.1) * times,
  height = (1.8 * n_bar + 1) * 0.75 * times,
  units = "in",
  res = 300
)

# Initialize empty plot
plot(NA, type = "n", xlim = c(-0.7, x_length - 0.4), ylim = c(-0.5, 1.7 * n_bar), xlab = "", ylab = "", main = "", axes = FALSE, asp = NA)

# Add bottom x-axis
axis(side = 1, at = seq(0, x_length + 1, by = 1), pos = -0.3, cex.axis = times * 1, padj = -0.2, tck = -0.032 * 0.75)

text(0, -0.1, "Pseudo ", cex = times * 1.2, adj = 1)
text(0, -0.5, "time ", cex = times * 1.2, adj = 1)

# Draw bars branch by branch
for (n in 1:n_bar) {
  v1 = c(0, interest_cell_type_group[["KP4"]][["Branches_Tslot_k"]][[n]]); v1 # time-slot boundaries
  v2_0 = c(0, interest_cell_type_group[["KP4"]][["Branches_n_k"]][[n]]); v2_0 # cell counts per group
  v2 = v2_0 / sum(v2_0) * sum(v1); v2 # scale cell counts to the time axis
  cum1 = cumsum(v1); cum1
  cum2 = cumsum(v2); cum2

  # Draw a trapezoid for each group
  for (i in 1:length(v1)) {
    x_poly = c(cum1[i], cum1[i + 1], cum2[i + 1], cum2[i])
    y_poly = c(n - 1 + 0.2 * (n - 1) + 0.6 * n, n - 1 + 0.2 * (n - 1) + 0.6 * n, n + 0.2 * (n - 1) + 0.6 * n, n + 0.2 * (n - 1) + 0.6 * n)
    polygon(x_poly, y_poly, col = colors[i], border = "gray45", lwd = times * 1)
  }

  # Draw outer rectangle
  polygon(
    c(0, cum2[length(v1)], cum2[length(v1)], 0),
    c(n - 1 + 0.2 * (n - 1) + 0.6 * n, n - 1 + 0.2 * (n - 1) + 0.6 * n, n + 0.2 * (n - 1) + 0.6 * n, n + 0.2 * (n - 1) + 0.6 * n),
    col = NULL,
    border = "black",
    lwd = times * 1
  )

  # Add labels inside each trapezoid
  for (i in 1:length(v1)) {
    # Time-slot label (bottom), with font size adapted to slot width
    mid_x1 = (cum1[i] + cum1[i+1]) / 2
    if (!is.na(round(v1, digits = 2)[i + 1])) {
      if (round(v1, digits = 2)[i + 1] < 0.25 && round(v1, digits = 2)[i + 1] > 0.15) {
        text(mid_x1, n - 0.85 + 0.2 * (n - 1) + 0.6 * n, labels = round(v1, digits = 1)[i + 1], col = "black", cex = times * 0.5)
      } else if (round(v1, digits = 2)[i + 1] < 0.15) {
        text(mid_x1, n - 0.85 + 0.2 * (n - 1) + 0.6 * n, labels = round(v1, digits = 1)[i + 1], col = "black", cex = times * 0.3)
      } else {
        text(mid_x1, n - 0.85 + 0.2 * (n - 1) + 0.6 * n, labels = round(v1, digits = 2)[i + 1], col = "black", cex = times * 0.8)
      }
    } else {
      text(mid_x1, n - 0.85 + 0.2 * (n - 1) + 0.6 * n, labels = round(v1, digits = 2)[i + 1], col = "black", cex = times * 0.8)
    }

    mid_x2 = (cum2[i] + cum2[i+1]) / 2 # cell count label (top)
    text(mid_x2, n - 0.15 + 0.2 * (n - 1) + 0.6 * n, labels = v2_0[i + 1], col = "black", cex = times * 0.8)

    x_center = mean(c(cum1[i], cum1[i + 1], cum2[i + 1], cum2[i])) # group number label (center)
    y_center = n - 0.5 + 0.2 * (n - 1) + 0.6 * n
    text(x_center, y_center, labels = i, col = "black", cex = times * 1.2, font = 2)
  }

  # Add left-side labels
  text(0, n - 1.2 + 0.2 * (n - 1) + 0.6 * n, "State ", col = "black", cex = times * 0.8, adj = 1)
  text(0, n - 0.85 + 0.2 * (n - 1) + 0.6 * n, "Period ", col = "black", cex = times * 0.8, adj = 1)
  text(0, n - 0.15 + 0.2 * (n - 1) + 0.6 * n, "Cells ", col = "black", cex = times * 0.8, adj = 1)
  text(0, n - 0.5 + 0.2 * (n - 1) + 0.6 * n, "Group ", col = "black", cex = times * 0.8, adj = 1)
  text(-0.75, n + 0 + 0.2 * (n - 1) + 0.6 * n, paste0("Branch", n), col = "black", cex = times * 1.0, adj = 1, srt = 90)

  # Add right-side summary labels
  text(sum(v1), n - 0.85 + 0.2 * (n - 1) + 0.6 * n, paste0(" ", round(sum(v1), digits = 2)), col = "black", cex = times * 0.8, adj = 0)
  text(sum(v1), n - 0.15 + 0.2 * (n - 1) + 0.6 * n, paste0(" ", sum(v2_0)), col = "black", cex = times * 0.8, adj = 0)

  # Add state divider lines
  states = unlist(strsplit(branches[n], "→"))
  pseudotime_info = interest_cell_type_genes_pseudotime_info[["KP4"]][["Branches"]][[n]]
  pseudotime_info$State = factor(pseudotime_info$State, levels = states)
  x_states = c(0, na.omit(as.vector(tapply(pseudotime_info$Pseudotime, pseudotime_info$State, max))))
  for (s in 1:(length(states) + 1)) {
    lines(
      x = rep(x_states[s], 2),
      y = c(n - 1 + 0.2 * (n - 1) + 0.6 * n - 0.4, n - 1 + 0.2 * (n - 1) + 0.6 * n),
      lty = "dashed",
      lwd = times * 1
    )
  }; rm(pseudotime_info)

  # Add state labels
  x_states_pos = (x_states[-length(x_states)] + x_states[-1]) / 2
  for (s in 1:length(x_states_pos)) {
    text(x_states_pos[s], n - 1.2 + 0.2 * (n - 1) + 0.6 * n, states[s], col = "black", cex = times * 0.8, adj = 0.5)
  }
}

# Close device
dev.off()



##### 6.4 Line Chart of Cell Number and Pseudotime Period in Branches #####
setwd(Working_directory)
library(Cairo)
library(ggplot2)
library(patchwork)
interest_cell_type_group = qs::qread("./2.2 Data Processing - Cell Grouping/interest_cell_type_group.qsave")
cell_counts = interest_cell_type_group[["KP4"]][["Branches_n_k"]]
pseudotimes = interest_cell_type_group[["KP4"]][["Branches_Tslot_k"]]

# Per-branch line chart of average pseudotime period per cell group
Line_Chart_of_Cell_Average_Pseudotime_Period_list = list()
for (i in seq_along(cell_counts)) {
  df = data.frame(
    Group = seq_along(cell_counts[[i]]),
    Cell_Average_Pseudotime_Period = pseudotimes[[i]] / cell_counts[[i]]
  )

  Line_Chart_of_Cell_Average_Pseudotime_Period_list[[i]] = ggplot(df, aes(x = Group, y = Cell_Average_Pseudotime_Period)) +
    geom_line(linewidth = 1, color = "grey50") +
    geom_point(size = 2, color = "grey50") +
    scale_x_continuous(breaks = seq(min(df$Group), max(df$Group), by = 1)) +
    scale_y_continuous(name = "Cell Average Pseudotime Period") +
    labs(title = paste("Branch", i), x = "Cell Group") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 8),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank()
    )

  ggplot2::ggsave(paste0("./6 Plot for Experiment/6.4_Line_Chart_of_Cell_Average_Pseudotime_Period_in_Branch_", i, ".png"),
                  plot = Line_Chart_of_Cell_Average_Pseudotime_Period_list[[i]], width = 3.5, height = 2.5, dpi = 300)
}

# Combined panel of per-branch charts
Line_Chart_of_Cell_Average_Pseudotime_Period1 =
  (Line_Chart_of_Cell_Average_Pseudotime_Period_list[[1]] | Line_Chart_of_Cell_Average_Pseudotime_Period_list[[2]]) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")

Line_Chart_of_Cell_Average_Pseudotime_Period1 = Line_Chart_of_Cell_Average_Pseudotime_Period1 +
  plot_annotation(
    title = "Cell Average Pseudotime Period",
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  )

ggplot2::ggsave("./6 Plot for Experiment/6.4_Line_Chart_of_Cell_Average_Pseudotime_Period1.png",
                plot = Line_Chart_of_Cell_Average_Pseudotime_Period1,
                width = 6, height = 2.8, dpi = 300)

# Overlaid chart of all branches
library(tidyverse)
df_combined = imap_dfr(cell_counts, ~ data.frame(
  Branch = factor(paste("Branch", .y)),
  Group = seq_along(.x),
  Cell_Average_Pseudotime_Period = pseudotimes[[.y]] / .x
))

Line_Chart_of_Cell_Average_Pseudotime_Period2 = ggplot(df_combined, aes(x = Group, y = Cell_Average_Pseudotime_Period, color = Branch)) +
  geom_line(aes(linewidth = Branch, alpha = Branch)) +
  geom_point(aes(size = Branch, alpha = Branch)) +
  scale_color_manual(values = c("Branch 1" = "#E64B35", "Branch 2" = "#3182BD")) +
  scale_linewidth_manual(values = c("Branch 1" = 1, "Branch 2" = 1)) +
  scale_size_manual(values = c("Branch 1" = 2, "Branch 2" = 2)) +
  scale_alpha_manual(values = c("Branch 1" = 0.5, "Branch 2" = 0.5)) +
  scale_x_continuous(breaks = ~ seq(1, max(.x), by = 1)) +
  scale_y_continuous(name = "Cell Average Pseudotime Period") +
  labs(title = "Cell Average Pseudotime Period", x = "Cell Group") +
  theme_minimal() + theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 6),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  )
print(Line_Chart_of_Cell_Average_Pseudotime_Period2)
ggsave("./6 Plot for Experiment/6.4_Line_Chart_of_Cell_Average_Pseudotime_Period2.png",
       plot = Line_Chart_of_Cell_Average_Pseudotime_Period2,
       width = 3.5, height = 2.5, dpi = 300)




##### 6.5 UpSet Plot of Genes Information in Branches #####
# TFs and TGs are identical across the two branches, so the UpSet plot is skipped
#
# setwd(Working_directory)
# library(UpSetR)
# library(grid)
# library(Cairo)
# interest_cell_type_genes_pseudotime_info = qs::qread("./2.2 Data Processing - Cell Grouping/interest_cell_type_genes_pseudotime_info.qsave")
#
# # Create UpSet plots for TFs and TGs
# for (gene_type in c("TFs", "TGs")) {
#   # Run separately
#   # gene_type = "TFs"
#   gene_type = "TGs"
#
#   # Get gene lists for each branch
#   venn_list = interest_cell_type_genes_pseudotime_info[["KP4"]][[paste0("Branches_", ifelse(gene_type == "TFs", "TFs", "TGs"))]]
#   names(venn_list) = paste0("Branch", 1:length(venn_list))
#
#   # Create binary matrix for UpSet plot
#   all_genes = unique(unlist(venn_list))
#   branch_names = names(venn_list)
#   binary_matrix = matrix(0, nrow = length(all_genes), ncol = length(branch_names), dimnames = list(all_genes, branch_names))
#   for (col_name in branch_names) {
#     binary_matrix[all_genes %in% venn_list[[col_name]], col_name] = 1
#   }
#
#   # Set scaling factor
#   times = 0.6
#
#   # Create PNG file
#   CairoPNG(
#     filename = paste0("./6 Plot for Experiment/6.5_Branches_", gene_type, "_UpSetPlot.png"),
#     width = (2 * length(venn_list) + 4) * times,
#     height = (1.5 * length(venn_list) + 4) * times,
#     units = "in",
#     res = 300
#   )
#
#   # Generate UpSet plot
#   upset(
#     data = as.data.frame(binary_matrix),
#     nsets = length(venn_list),
#     nintersects = NA,
#     sets = NULL,
#     keep.order = FALSE,
#     set.metadata = NULL,
#     intersections = NULL,
#     matrix.color = "gray23",
#     main.bar.color = "gray23",
#     mainbar.y.label = "Intersection Size",
#     mainbar.y.max = NULL,
#     sets.bar.color = "gray23",
#     sets.x.label = "Set Size",
#     point.size = 2.2,
#     line.size = 0.7,
#     mb.ratio = c(0.6, 0.4),
#     expression = NULL,
#     att.pos = NULL,
#     order.by = c("freq", "degree"),
#     decreasing = c(TRUE, FALSE),
#     show.numbers = "yes",
#     number.angles = 0,
#     group.by = "degree",
#     cutoff = NULL,
#     queries = NULL,
#     query.legend = "none",
#     shade.color = "gray88",
#     shade.alpha = 0.25,
#     matrix.dot.alpha = 0.5,
#     empty.intersections = NULL,
#     color.pal = 1,
#     boxplot.summary = NULL,
#     attribute.plots = NULL,
#     scale.intersections = "identity",
#     scale.sets = "identity",
#     text.scale = 2.5 * times * 0.95,
#     set_size.angles = 0,
#     set_size.show = TRUE,
#     set_size.numbers_size = NULL,
#     # set_size.scale_max = nrow(binary_matrix) + ifelse(gene_type == "TFs", 20, 1000)
#   )
#
#   # Add text label
#   grid.text(
#     label = paste("KP4", " Total ", gene_type, ":", nrow(binary_matrix)),
#     x = unit(0.5, "npc"),
#     y = unit(0.05, "npc"),
#     gp = gpar(fontsize = 20 * times, fontface = "bold"),
#     just = c("left", "top")
#   )
#
#   # Save plot
#   dev.off()
# }




########## Part 7: Plot for Analysis ##########
##### 7.1 Enrichment Analysis and Plot #####
##### 7.1.1 Enrichment Analysis
setwd(Working_directory)
library(dplyr)
library(clusterProfiler)
library(ggplot2)
library(org.Hs.eg.db)
library(enrichplot)

# Load disease-associated eGRN
disease_associated_eGRN = qs::qread("./5 Infer GRN/disease_associated_eGRN.qsave")
eGRN = disease_associated_eGRN$TF_peaks_TGs[, c("TF", "TG")]

# Build one regulon (TF plus its TGs) per TF
regulons = list()
tfs = sort(unique(eGRN$TF))
n_tfs = length(tfs)
time1 = Sys.time(); for (tf in sort(unique(eGRN$TF))) {
  pos = match(tf, tfs)
  running_time = format(Sys.time() - time1)
  running_time_spilt = strsplit(running_time, " ")[[1]]
  message(tf, " - ", pos, "/", n_tfs, " - running time: ", format(Sys.time() - time1), " - expected remaining time: ",
          round(as.numeric(running_time_spilt[1]) / pos * (n_tfs - pos), digits = 3), " ", running_time_spilt[2], " - expected running time: ",
          round(as.numeric(running_time_spilt[1]) / pos * n_tfs, digits = 3), " ", running_time_spilt[2], ".")
  regulons[[tf]] = list()
  regulons[[tf]][["SYMBOL"]] = unique(c(tf, eGRN[eGRN$TF == tf, "TG"]))
  regulons[[tf]][["ENTREZID"]] = bitr(regulons[[tf]][["SYMBOL"]], fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
}; message("Running time: ", format(Sys.time() - time1), ".")
# Running time: 10.8485 secs.

# Regulon size (TF plus TGs) per TF
for (tf in names(regulons)) {message(tf, ": ", length(regulons[[tf]][["SYMBOL"]]))}
# ATF6: 109
# BHLHE40: 182
# CEBPB: 168
# E2F1: 175
# E2F2: 182
# E2F7: 150
# ELF4: 166
# ETV6: 170
# FOSL1: 173
# FOXC1: 184
# IRF2: 140
# KLF10: 171
# KLF11: 183
# KLF5: 161
# KLF6: 180
# MAFF: 182
# MXI1: 179
# MYBL2: 164
# MYC: 179
# STAT1: 160
# TFAP4: 184
# TP63: 161
# TRPS1: 4

# Convert a "count/total" fraction string to numeric
convert_fraction = function(fraction_str) {
  sapply(strsplit(fraction_str, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
}

# Filtering thresholds
count_threshold = 3
generatio_threshold = 0
pq_threshold = 0.05

# Run GO/KEGG/MKEGG enrichment for each regulon
options(timeout = 99999)

enrich = list()
time1 = Sys.time(); for (tf in tfs) {
  # tf = tfs[1] # For test
  # names(enrich) # For test
  # match(tf, tfs) # For test

  pos = match(tf, tfs)
  running_time = format(Sys.time() - time1)
  running_time_spilt = strsplit(running_time, " ")[[1]]
  message(tf, " - ", pos, "/", n_tfs, " - running time: ", format(Sys.time() - time1), " - expected remaining time: ",
          round(as.numeric(running_time_spilt[1]) / pos * (n_tfs - pos), digits = 3), " ", running_time_spilt[2], " - expected running time: ",
          round(as.numeric(running_time_spilt[1]) / pos * n_tfs, digits = 3), " ", running_time_spilt[2], ".")

  message(tf, " - TG count: ", length(regulons[[tf]][["SYMBOL"]]) - 1)

  enrich[[tf]] = list()

  ### GO enrichment ###
  enrich_GO = enrichGO(
    gene = regulons[[tf]][["SYMBOL"]], OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "ALL",
    minGSSize = 2, readable = TRUE, pvalueCutoff = 0.05, qvalueCutoff = 0.05, pAdjustMethod = "BH"
  )

  enrich[[tf]][["GO"]] = enrich_GO
  if (!is.null(enrich_GO) & nrow(enrich_GO@result) > 0) {
    enrich_GO_filtered = clusterProfiler::simplify(enrich_GO)
    result = enrich_GO_filtered@result
    result = result[order(-result$FoldEnrichment), ]
    result = result[result$Count >= count_threshold & convert_fraction(result$GeneRatio) >= generatio_threshold, ]
    result = result[result$pvalue <= pq_threshold & result$p.adjust <= pq_threshold, ]
    enrich_GO_filtered@result = result; rm(result)
  } else {enrich_GO_filtered = enrich_empty}
  enrich[[tf]][["GO_filtered"]] = enrich_GO_filtered
  if (!is.null(enrich_GO)) {
    message(tf, " - GO Pathway count: ", nrow(enrich[[tf]][["GO"]]@result), "/", nrow(enrich[[tf]][["GO_filtered"]]@result))
  }

  # Template for an empty enrichment result
  enrich_empty = enrich_GO
  enrich_empty@result = enrich_empty@result[FALSE, ]
  enrich_empty@ontology = character(0)
  enrich_empty@gene = character(0)
  enrich_empty@universe = character(0)
  enrich_empty@geneSets = list()
  enrich[["empty"]] = enrich_empty

  ### KEGG enrichment ###
  enrich_KEGG = enrichKEGG(
    gene = regulons[[tf]][["ENTREZID"]], keyType = "kegg", organism = "hsa",
    minGSSize = 2, use_internal_data = FALSE, pvalueCutoff = 0.05, qvalueCutoff = 0.05, pAdjustMethod = "BH"
  )
  enrich[[tf]][["KEGG"]] = enrich_KEGG

  if (!is.null(enrich_KEGG)) {
    enrich_KEGG = setReadable(enrich_KEGG, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    enrich[[tf]][["KEGG"]] = enrich_KEGG
    enrich_KEGG_filtered = enrich_KEGG
    result = enrich_KEGG_filtered@result
    result = result[order(-result$FoldEnrichment), ]
    result = result[result$Count >= count_threshold & convert_fraction(result$GeneRatio) >= generatio_threshold, ]
    result = result[result$pvalue <= pq_threshold & result$p.adjust <= pq_threshold, ]
    enrich_KEGG_filtered@result = result; rm(result)
  } else {enrich_KEGG_filtered = NULL}
  enrich[[tf]][["KEGG_filtered"]] = enrich_KEGG_filtered
  if (!is.null(enrich_KEGG)) {
    message(tf, " - KEGG Pathway count: ", nrow(enrich[[tf]][["KEGG"]]@result), "/", nrow(enrich[[tf]][["KEGG_filtered"]]@result))
  }

  ### KEGG Module enrichment ###
  enrich_MKEGG = enrichMKEGG(
    gene = regulons[[tf]][["ENTREZID"]], keyType = "kegg", organism = "hsa",
    minGSSize = 2, pvalueCutoff = 0.05, qvalueCutoff = 0.05, pAdjustMethod = "BH"
  )
  enrich[[tf]][["MKEGG"]] = enrich_MKEGG

  if (!is.null(enrich_MKEGG)) {
    enrich_MKEGG = setReadable(enrich_MKEGG, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    enrich[[tf]][["MKEGG"]] = enrich_MKEGG
    enrich_MKEGG_filtered = enrich_MKEGG
    result = enrich_MKEGG_filtered@result
    result = result[order(-result$FoldEnrichment), ]
    result = result[result$Count >= count_threshold & convert_fraction(result$GeneRatio) >= generatio_threshold, ]
    result = result[result$pvalue <= pq_threshold & result$p.adjust <= pq_threshold, ]
    enrich_MKEGG_filtered@result = result; rm(result)
  } else {enrich_MKEGG_filtered = enrich_empty}
  enrich[[tf]][["MKEGG_filtered"]] = enrich_MKEGG_filtered
  if (!is.null(enrich_MKEGG)) {
    message(tf, " - MKEGG Pathway count: ", nrow(enrich[[tf]][["MKEGG"]]@result), "/", nrow(enrich[[tf]][["MKEGG_filtered"]]@result))
  }

  ### Combine all enrichment results ###
  enrich_all_filtered = enrich_GO_filtered

  # Add GO results
  if (nrow(enrich_all_filtered@result) > 0) {
    enrich_all_filtered@result$ID = paste0(
      enrich_all_filtered@result$ONTOLOGY, ".", enrich_all_filtered@result$ID
    )
  }
  enrich_all_filtered@result$ONTOLOGY = NULL

  # Add KEGG results
  if (!is.null(enrich_KEGG_filtered)) {
    # Update result
    result = enrich_KEGG_filtered@result
    result$ID = paste(result$category, result$subcategory, result$ID, sep = ".")
    result$category = NULL; result$subcategory = NULL
    enrich_all_filtered@result = rbind(enrich_all_filtered@result, result); rm(result)

    # Update ontology
    enrich_all_filtered@ontology = c(enrich_all_filtered@ontology, enrich_KEGG_filtered@ontology)

    # Update universe
    enrich_all_filtered@universe = unique(c(
      enrich_all_filtered@universe,
      bitr(enrich_KEGG_filtered@universe, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)$SYMBOL
    ))

    # Update gene sets from ENTREZID to SYMBOL
    geneSets = enrich_KEGG_filtered@geneSets
    entrez2symbol_map = bitr(unique(unlist(geneSets, use.names = FALSE)), fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
    map_vector = setNames(entrez2symbol_map$SYMBOL, entrez2symbol_map$ENTREZID)
    geneSets = lapply(geneSets, function(single_geneset) {na.omit(map_vector[single_geneset])}); rm(entrez2symbol_map, map_vector)
    enrich_all_filtered@geneSets = c(enrich_all_filtered@geneSets, geneSets); rm(geneSets)
  }

  # Add MKEGG results
  if (!is.null(enrich_MKEGG_filtered)) {
    # Update result
    enrich_all_filtered@result = rbind(enrich_all_filtered@result, enrich_MKEGG_filtered@result)

    # Update ontology
    enrich_all_filtered@ontology = c(enrich_all_filtered@ontology, enrich_MKEGG_filtered@ontology)

    # Update universe
    enrich_all_filtered@universe = unique(c(
      enrich_all_filtered@universe,
      bitr(enrich_MKEGG_filtered@universe, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)$SYMBOL
    ))

    # Update gene sets from ENTREZID to SYMBOL
    geneSets = enrich_MKEGG_filtered@geneSets
    entrez2symbol_map = bitr(unique(unlist(geneSets, use.names = FALSE)), fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
    map_vector = setNames(entrez2symbol_map$SYMBOL, entrez2symbol_map$ENTREZID)
    geneSets = lapply(geneSets, function(single_geneset) {na.omit(map_vector[single_geneset])}); rm(entrez2symbol_map, map_vector)
    enrich_all_filtered@geneSets = c(enrich_all_filtered@geneSets, geneSets); rm(geneSets)
  }

  # Final filtering
  enrich_all_filtered@result = enrich_all_filtered@result %>%
    dplyr::filter(pvalue < 0.05, p.adjust < 0.05, Count > 0) %>%
    dplyr::arrange(dplyr::desc(FoldEnrichment))
  enrich[[tf]][["all"]] = enrich_all_filtered

  message("Check point.")
  # save.image("./enrich_check_point.RData")
  # load("./enrich_check_point.RData")

  qs::qsave(enrich, file = "./enrich_check_point.qsave")
  # enrich = qs::qread("./enrich_check_point.qsave")
}; message("Running time: ", format(Sys.time() - time1), ".")
qs::qsave(enrich, file = "./7 Plot for Analysis/enrich.qsave")
# enrich = qs::qread("./7 Plot for Analysis/enrich.qsave")
# 54.51232 mins

# Further filter enrichment results for each regulon
enrichment = enrich
enrichment[["empty"]] = NULL
count_thres = 5
for (tf in names(enrichment)) {
  enrichment[[tf]] = list()
  enrichment[[tf]][["all"]] = enrich[[tf]][["all"]]
  enrichment[[tf]][["all_filtered"]] = enrichment[[tf]][["all"]]
  # rownames(enrichment[[tf]][["all_filtered"]]@result) = 1:nrow(enrichment[[tf]][["all_filtered"]]@result)
  enrichment[[tf]][["count_filtered"]] = enrichment[[tf]][["all_filtered"]]
  x = enrichment[[tf]][["count_filtered"]]@result
  x = x[x$Count >= count_thres, ]
  x = x[grepl(tf, x$geneID), ]
  if(nrow(x) == 0) {
    enrichment[[tf]] = NULL
    next
  }
  x$TF = tf
  rownames(x) = paste0(tf, "-", rownames(x))
  enrichment[[tf]][["count_filtered"]]@result = x

  enrichment[[tf]][["all"]] = NULL # Release memory
  enrichment[[tf]][["all_filtered"]] = NULL # Release memory
  # enrichment[[tf]][["count_filtered"]]@universe = "" # Release memory

  # message(tf, ": ", nrow(enrichment[[tf]][["count_filtered"]]))

}; rm(enrich, x)
qs::qsave(enrichment, file = "./7 Plot for Analysis/enrichment.qsave")
# enrichment = qs::qread("./7 Plot for Analysis/enrichment.qsave")

# Retained pathway count per TF
for (tf in names(enrichment)) {message(tf, ": ", nrow(enrichment[[tf]][["count_filtered"]]))}
# ATF6: 1
# BHLHE40: 1
# CEBPB: 7
# E2F1: 24
# E2F2: 3
# E2F7: 12
# FOSL1: 5
# FOXC1: 4
# IRF2: 2
# KLF10: 2
# KLF11: 3
# MAFF: 3
# MXI1: 1
# MYBL2: 7
# MYC: 12
# STAT1: 9
# TFAP4: 6
# TP63: 7

# Pool filtered pathways across regulons
enrichment[["count_filtered_all"]] = enrichment[[1]][["count_filtered"]]
enrichment[["count_filtered_all"]]@result = do.call(rbind, lapply(names(enrichment)[1:(length(enrichment) - 1)], function(tf) {
  df = enrichment[[tf]][["count_filtered"]]@result
  df
}))
dim(enrichment[["count_filtered_all"]]@result)
enrichment[["count_filtered_all"]]@result[, c("ID", "Description", "GeneRatio", "BgRatio", "FoldEnrichment")]

# L1: the only pathway explicitly annotated as pancreatic cancer
level_1 = c(
  "Human Diseases.Cancer: specific types.hsa05212" # KEGG Pancreatic cancer pathway
)

# L2: strongly associated with PDAC biology but not PDAC-specific
level_2 = c(
  "BP.GO:0002526", # acute inflammatory response; key component of the PDAC TME, matching the pro-inflammatory role of CEBPB
  "BP.GO:0006979", # oxidative stress; hallmark of the PDAC TME, linked to progression and chemoresistance
  "BP.GO:0030856", # epithelial cell differentiation; PDAC arises from ductal epithelium, but this is a generic process (borderline L2)
  "BP.GO:0042542", # response to H2O2; KRAS-driven ROS in the PDAC TME, and STAT1 can be activated by oxidative stress
  "BP.GO:0045639", # myeloid cell differentiation; MDSCs/TAMs underpin the immunosuppressive PDAC TME (borderline L2)
  "BP.GO:0048144", # fibroblast proliferation; CAF expansion forms the dense desmoplastic stroma characteristic of PDAC
  "BP.GO:0048145", # regulation of fibroblast proliferation; CAF/desmoplasia is a defining feature of PDAC
  "BP.GO:0050729", # negative regulation of inflammatory response; chronic inflammation is a PDAC hallmark and CEBPB drives inflammation-linked progression
  "Environmental Information Processing.Signal transduction.hsa04310", # Wnt/beta-catenin is aberrantly active in PDAC and cross-talks with FOSL1 (AP-1)
  "MF.GO:0045296" # cadherin binding; linked to adhesion, EMT and cadherin switch in PDAC invasion (borderline L2)
)

# L3: canonical pan-cancer hallmarks
level_3 = c(
  "BP.GO:0000079", # regulation of CDK activity; CDK4/6 are candidate therapeutic targets in PDAC
  "BP.GO:0000280", # nuclear division, core proliferation process
  "BP.GO:0000910", # cytokinesis, proliferation hallmark
  "BP.GO:0006260", # DNA replication, core proliferation process
  "BP.GO:0006261", # DNA-templated DNA replication, core proliferation process
  "BP.GO:0006935", # chemotaxis linked to invasion and metastasis; FOSL1/AP-1 implicated in PDAC metastasis (borderline L3)
  "BP.GO:0007051", # spindle organization, core mitotic process; chromosomal instability hallmark
  "BP.GO:0007059", # chromosome segregation, core mitotic process; chromosomal instability hallmark
  "BP.GO:0009410", # xenobiotic/drug response, linked to PDAC chemoresistance (borderline L3)
  "BP.GO:0010212", # ionizing radiation response, DNA damage response; genome instability hallmark
  "BP.GO:0032465", # regulation of cytokinesis, proliferation hallmark
  "BP.GO:0033044", # regulation of chromosome organization, linked to genome stability
  "BP.GO:0042770", # DNA damage signal transduction; KRAS-driven replication stress in PDAC
  "BP.GO:0044772", # mitotic cell cycle phase transition, core proliferation process
  "BP.GO:0044843", # G1/S transition, core cell-cycle checkpoint
  "BP.GO:0045787", # positive regulation of cell cycle, proliferation hallmark
  "BP.GO:0051054", # positive regulation of DNA metabolic process, core proliferation
  "BP.GO:0051302", # regulation of cell division, proliferation hallmark
  "BP.GO:0072331", # p53-class mediator signaling; TP63 is a p53 family member and p53 is frequently inactivated in PDAC
  "BP.GO:0090068", # positive regulation of cell cycle process, proliferation hallmark, FE=11.48
  "BP.GO:0090329", # regulation of DNA-templated DNA replication; E2F7 is an atypical repressive E2F, FE=15.6
  "BP.GO:0097193", # intrinsic apoptotic signaling, core cancer hallmark (evasion); CEBPB mediates apoptosis resistance in PDAC
  "BP.GO:0098813", # nuclear chromosome segregation, core mitotic process
  "BP.GO:0140014", # mitotic nuclear division, M-phase core; MYBL2 is a master mitotic transcription factor, FE=11.43
  "BP.GO:1901987", # regulation of cell cycle phase transition, core proliferation
  "BP.GO:1901990", # regulation of mitotic cell cycle phase transition, core cell-cycle process
  "BP.GO:2001233", # regulation of apoptotic signaling, core cancer hallmark (evasion)
  "Cellular Processes.Cell growth and death.hsa04110", # cell cycle, core cancer hallmark
  "Cellular Processes.Cell growth and death.hsa04218", # cellular senescence, linked to tumor suppression/escape, pan-cancer core process
  "Human Diseases.Cancer: overview.hsa05202", # transcriptional misregulation in cancer, pan-cancer feature
  "Human Diseases.Cancer: overview.hsa05206" # microRNAs in cancer, pan-cancer feature
)

table(level_1 %in% enrichment[["count_filtered_all"]]@result$ID) # 1 match
table(level_2 %in% enrichment[["count_filtered_all"]]@result$ID) # 10 matches
table(level_3 %in% enrichment[["count_filtered_all"]]@result$ID) # 31 matches

# Assign PDAC relevance levels (1 = highest, 3 = lowest)
enrichment[["count_filtered_all"]]@result$importance_level = "No" # not relevant
enrichment[["count_filtered_all"]]@result$importance_level[enrichment[["count_filtered_all"]]@result$ID %in% level_1] = 1 # L1: the only pathway explicitly tied to pancreatic cancer
enrichment[["count_filtered_all"]]@result$importance_level[enrichment[["count_filtered_all"]]@result$ID %in% level_2] = 2 # L2: strongly related to PDAC biology but not PDAC-specific
enrichment[["count_filtered_all"]]@result$importance_level[enrichment[["count_filtered_all"]]@result$ID %in% level_3] = 3 # L3: canonical pan-cancer hallmarks

enrichment[["count_filtered_all"]]@result = enrichment[["count_filtered_all"]]@result[enrichment[["count_filtered_all"]]@result$importance_level %in% c("1", "2", "3"), ]

# Synchronize per-regulon results with the retained pathway set
for (tf in setdiff(names(enrichment), "count_filtered_all")) {
  enrichment[[tf]][["count_filtered"]]@result = enrichment[[tf]][["count_filtered"]]@result[
    enrichment[[tf]][["count_filtered"]]@result$ID %in% enrichment[["count_filtered_all"]]@result$ID,
    ,
    drop = FALSE
  ]
  if (nrow(enrichment[[tf]][["count_filtered"]]@result) <= 0) {
    enrichment[[tf]] = NULL
  } else {
    enrichment[[tf]] = enrichment[[tf]][["count_filtered"]]
  }
}
names(enrichment)[names(enrichment) == "count_filtered_all"] = "all"

writexl::write_xlsx(enrichment[["all"]]@result, "./7 Plot for Analysis/7.1.1_enrichment_result.xlsx")
dim(enrichment[["all"]]@result) # [1] 66 14
# enrichment[["all"]]@result[, c("ID", "Description", "GeneRatio", "BgRatio", "FoldEnrichment", "importance_level")]

# Subset of L1/L2 pathways
enrichment[["L12"]] = enrichment[["all"]]
enrichment[["L12"]]@result = enrichment[["L12"]]@result[enrichment[["L12"]]@result$importance_level %in% c("1", "2"), ]

qs::qsave(enrichment, file = "./7 Plot for Analysis/enrichment.qsave")
# enrichment = qs::qread("./7 Plot for Analysis/enrichment.qsave")


##### 7.1.2 Enrichment Dot Plots
setwd(Working_directory)
library(enrichplot)
library(ggplot2)
library(ggh4x)
library(ggpubr)

enrichment = qs::qread("./7 Plot for Analysis/enrichment.qsave")

# Dot plot 1: pooled L1/L2 pathways
enrichment[["L12"]]@result$Description = paste0(enrichment[["L12"]]@result$TF, " - ", enrichment[["L12"]]@result$Description)

L12_Enrichment_Dotplot = dotplot(
  enrichment[["L12"]],
  x = "FoldEnrichment",
  color = "p.adjust",
  size = "Count",
  font.size = 7.5,
  showCategory = 10000,
  label_format = 31
)
print(L12_Enrichment_Dotplot)
ggsave(filename = "./7 Plot for Analysis/7.1.2_L12_Enrichment_Dotplot.png",
       plot = L12_Enrichment_Dotplot,
       width = 4, height = 4, dpi = 300)

# Dot plot 2: L1-L3 pathways faceted by regulon
Regulon_Enrichment_Dotplot = dotplot(
  enrichment[["all"]],
  x = "FoldEnrichment",
  color = "p.adjust",
  size = "Count",
  font.size = 7.5,
  showCategory = 10000,
  label_format = 31
)+ facet_wrap(~ TF, nrow = 3, scales = "free_y", dir = "v")
print(Regulon_Enrichment_Dotplot)
ggsave(filename = paste0("./7 Plot for Analysis/7.1.2_Regulon_L123_Enrichment_Dotplot.png"),
       plot = Regulon_Enrichment_Dotplot,
       width = 16, height = 10, dpi = 300)


##### 7.1.3 Gene-Pathway Heatmap
setwd(Working_directory)
library(dplyr)
library(ggplot2)
library(reshape2)

enrichment = qs::qread("./7 Plot for Analysis/enrichment.qsave")

# Clean pathway IDs
enrich_result = enrichment[["L12"]]
enrich_result@result = enrich_result@result[, c("ID", "geneID")]
enrich_result@result$ID = gsub("-", "_", gsub(":", ".", enrich_result@result$ID))
enrich_result@result$ID = gsub(" ", "_", enrich_result@result$ID)
enrich_result@result$ID = gsub(",", "", enrich_result@result$ID)
enrich_result@result = enrich_result@result %>% distinct(ID, geneID, .keep_all = TRUE)
enrich_result@result = as.data.frame(enrich_result@result %>% group_by(ID) %>% summarise(geneID = paste(geneID, collapse = "/"), .groups = 'drop'))

# Build binary gene-pathway matrix
genes = sort(unique(unlist(strsplit(enrich_result@result[["geneID"]], "/"))))
pathways = enrich_result@result$ID
pathways = gsub("-", "_", gsub(":", ".", pathways))

gene_and_pathway_matrix = data.frame(matrix(
  "No", nrow = length(genes), ncol = length(pathways), dimnames = list(genes, pathways)
))
for (pathway in pathways) {
  genes_in_pathway = enrich_result@result[enrich_result@result$ID == pathway, "geneID"]
  genes_in_pathway = unlist(strsplit(genes_in_pathway, "/"))
  for (gene in genes_in_pathway) {
    if (gene %in% genes) {gene_and_pathway_matrix[gene, pathway] = "Yes"}
  }
}

dim(gene_and_pathway_matrix)
class(gene_and_pathway_matrix)
gene_and_pathway_matrix[1:4,1:4]

# Shorten verbose KEGG pathway column names
colnames(gene_and_pathway_matrix)[colnames(gene_and_pathway_matrix) == "Environmental_Information_Processing.Signal_transduction.hsa04310"] = "hsa04310"
colnames(gene_and_pathway_matrix)[colnames(gene_and_pathway_matrix) == "Human_Diseases.Cancer._specific_types.hsa05212"] = "hsa05212"

# Convert to long format for ggplot2
gene_and_pathway_matrix$Gene = rownames(gene_and_pathway_matrix)
openxlsx::write.xlsx(as.data.frame(t(gene_and_pathway_matrix)),
                     file = "./7 Plot for Analysis/7.1.3_gene_and_pathway_matrix.xlsx",
                     rowNames = TRUE,
                     colNames = TRUE)
saveRDS(gene_and_pathway_matrix, file = "./7 Plot for Analysis/gene_and_pathway_matrix.rds")
data_long = melt(gene_and_pathway_matrix, id.vars = "Gene", variable.name = "Pathway", value.name = "Value")

# Build heatmap
times = 0.7
Gene_and_Pathway_Heatmap = ggplot(data_long, aes(x = Pathway, y = Gene, fill = as.factor(Value))) +
  geom_tile(color = "white", size = 0.5) +
  scale_fill_manual(values = c("No" = "grey95", "Yes" = "grey50"), name = "Value") +
  labs(x = "Pathway", y = "Gene") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6 * times),
    axis.text.y = element_text(size = 6 * times),
    panel.grid = element_blank(),
    legend.text = element_text(size = 6 * times),
    legend.title = element_text(size = 8 * times),
    axis.title.x = element_text(size = 10 * times),
    axis.title.y = element_text(size = 10 * times)
  )
print(Gene_and_Pathway_Heatmap)

# Save heatmap
ggsave(filename = "./7 Plot for Analysis/7.1.3_Gene_and_Pathway_Heatmap.png",
       plot = Gene_and_Pathway_Heatmap,
       width = 3.5 * times, height = 6 * times, dpi = 300)



##### 7.2 Network Plot of Disease-Associated eGRN #####
setwd(Working_directory)
disease_associated_eGRN = qs::qread("./5 Infer GRN/disease_associated_eGRN.qsave")
eGRN = disease_associated_eGRN[["TF_peaks_TGs"]][, c("TF", "peak", "TG", "FG-theta_p")]
colnames(eGRN)[4] = "theta_p"
rm(disease_associated_eGRN)

enrichment = qs::qread("./7 Plot for Analysis/enrichment.qsave")

library(dplyr)
library(tidyr)
library(igraph)
library(ggraph)
library(tidygraph)

# Plot a directed TF-peak-TG network; only peaks linked to >= 2 TFs or >= 2 TGs are labeled
plot_network = function(eGRN,
                        type_TF = 4,
                        type_TG = 3,
                        type_peak = 2) {
  # Build edge list
  edges_TF_to_peak = eGRN %>% select(from = TF, to = peak, weight = theta_p)
  edges_peak_to_TG = eGRN %>% select(from = peak, to = TG, weight = theta_p)
  edges = bind_rows(edges_TF_to_peak, edges_peak_to_TG)

  # Build node table from edge endpoints
  unique_nodes = unique(c(edges$from, edges$to))
  nodes = data.frame(
    name = unique_nodes,
    type = case_when(
      unique_nodes %in% eGRN$TF ~ "TF",
      unique_nodes %in% eGRN$TG ~ "TG",
      TRUE ~ "peak"
    )
  )

  # Compute per-peak connectivity and mark hub peaks
  peak_tf_count = eGRN %>% group_by(peak) %>% summarise(tf_count = n_distinct(TF))
  peak_tg_count = eGRN %>% group_by(peak) %>% summarise(tg_count = n_distinct(TG))
  peak_importance = full_join(peak_tf_count, peak_tg_count, by = "peak") %>% replace_na(list(tf_count = 0, tg_count = 0))
  peak_importance$is_important = with(peak_importance, tf_count >= 2 | tg_count >= 2)
  nodes = nodes %>% left_join(peak_importance, by = c("name" = "peak"))

  # Set label font size by node type
  nodes = nodes %>%
    mutate(font_size = case_when(
      type == "TF" ~ type_TF,
      type == "TG" ~ type_TG,
      type == "peak" ~ type_peak
    ))

  # Build tbl_graph object
  network = tbl_graph(nodes = nodes, edges = edges, directed = TRUE)

  # Render network
  NetworkPlot = ggraph(network, layout = "kk") +
    geom_edge_link(aes(edge_width = weight), color = "grey70", alpha = 0.6) +
    geom_node_point(aes(color = type), size = 1) +
    geom_node_text(
      aes(label = ifelse(type %in% c("TF", "TG") | (type == "peak" & is_important), name, ""), size = font_size),
      repel = TRUE, max.overlaps = 20
    ) +
    scale_edge_width(range = c(0.3, 2)) +
    scale_color_manual(values = c("TF" = "red", "peak" = "blue", "TG" = "green")) +
    scale_size_identity() +
    theme(plot.margin = margin(0, 0, 0, 0, "cm")) +
    labs(edge_width = "Edge Weight", color = "Node Type") +
    theme_graph()

  return(NetworkPlot)
}

### eGRN associated with L1/L2 pathways
L12_edges = enrichment[["L12"]]@result[, c("TF", "geneID")]
L12_pairs = data.frame(
  TF = rep(L12_edges$TF, lengths(strsplit(L12_edges$geneID, "/"))),
  geneID = unlist(strsplit(L12_edges$geneID, "/"))
)
eGRN_L12 = eGRN[
  paste(eGRN$TF, eGRN$TG) %in% paste(L12_pairs$TF, L12_pairs$geneID),
]; rm(L12_edges,L12_pairs)

L12_NetworkPlot = plot_network(eGRN = eGRN_L12, type_TF = 4, type_TG = 3, type_peak = 2)
print(L12_NetworkPlot)
ggsave("./7 Plot for Analysis/7.2_L12_NetworkPlot.png",
       plot = L12_NetworkPlot,
       width = 8, height = 6, dpi = 300)

openxlsx::write.xlsx(
  eGRN_L12,
  "./7 Plot for Analysis/7.2_L12_NetworkPlot.xlsx",
  rowNames = FALSE
)

### eGRN associated with L1-L3 pathways, one network per regulon
for (tf in setdiff(names(enrichment), c("all", "L12"))) {
  regulon_edges = enrichment[[tf]]@result[, c("TF", "geneID")]
  regulon_pairs = data.frame(
    TF = rep(regulon_edges$TF, lengths(strsplit(regulon_edges$geneID, "/"))),
    geneID = unlist(strsplit(regulon_edges$geneID, "/"))
  )
  eGRN_regulon = eGRN[
    paste(eGRN$TF, eGRN$TG) %in% paste(regulon_pairs$TF, regulon_pairs$geneID),
  ]; rm(regulon_edges,regulon_pairs)

  regulon_NetworkPlot = plot_network(eGRN = eGRN_regulon, type_TF = 4, type_TG = 3, type_peak = 2)
  # print(regulon_NetworkPlot)
  ggsave(paste0("./7 Plot for Analysis/7.2_", tf, "_L123_NetworkPlot.png"),
         plot = regulon_NetworkPlot,
         width = 6, height = 4.5, dpi = 300)
  openxlsx::write.xlsx(
    eGRN_regulon,
    paste0("./7 Plot for Analysis/7.2_", tf, "_L123_NetworkPlot.xlsx"),
    rowNames = FALSE
  )
}



##### 7.3 Connection Number Histogram of Disease-Associated eGRN Nodes #####
setwd(Working_directory)
disease_associated_eGRN = qs::qread("./5 Infer GRN/disease_associated_eGRN.qsave")
eGRN = disease_associated_eGRN[["TF_peaks_TGs"]][, c("TF", "peak", "TG")]

enrichment = qs::qread("./7 Plot for Analysis/enrichment.qsave")

# Subset eGRN to edges covered by L1-L3 pathways
L123_edges = enrichment[["all"]]@result[, c("TF", "geneID")]
L123_pairs = data.frame(
  TF = rep(L123_edges$TF, lengths(strsplit(L123_edges$geneID, "/"))),
  geneID = unlist(strsplit(L123_edges$geneID, "/"))
)
eGRN_L123 = eGRN[
  paste(eGRN$TF, eGRN$TG) %in% paste(L123_pairs$TF, L123_pairs$geneID),
]; rm(L123_edges,L123_pairs)
qs::qsave(eGRN_L123, file = "./7 Plot for Analysis/eGRN_L123.qsave")
# eGRN_L123 = qs::qread("./7 Plot for Analysis/eGRN_L123.qsave")

library(dplyr)
library(ggplot2)


##### 7.3.1 Bar plot: number of peaks per TF
# Peaks per TF, reflecting the regulatory breadth of each TF
tf_peak_count = eGRN_L123 %>%
  group_by(TF) %>%
  summarise(peak_count = n_distinct(peak))
saveRDS(tf_peak_count, "./7 Plot for Analysis/tf_peak_count.rds")
openxlsx::write.xlsx(
  as.data.frame(tf_peak_count),
  "./7 Plot for Analysis/7.3.1_tf_peak_count.xlsx",
  rowNames = FALSE
)

peak_number_per_TF_BarPlot = ggplot(tf_peak_count, aes(x = reorder(TF, -peak_count), y = peak_count)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
  geom_text(aes(label = peak_count), vjust = -0.1, size = 3, color = "black") +
  labs(x = "TF", y = "Number of regulated peaks") +
  scale_y_continuous(breaks = seq(0, max(tf_peak_count$peak_count) + 1, by = 5)) +
  theme_classic() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = 0.5, lineend = "butt"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
print(peak_number_per_TF_BarPlot)

ggsave("./7 Plot for Analysis/7.3.1_peak_number_per_TF_BarPlot_in_L123_eGRN.png", plot = peak_number_per_TF_BarPlot,
       width = nrow(tf_peak_count) * 0.2, height = 2.8, dpi = 300)


##### 7.3.2 Bar plot: number of TGs per TF
# TGs per TF, another measure of each TF's regulatory breadth
tf_tg_count = eGRN_L123 %>%
  group_by(TF) %>%
  summarise(tg_count = n_distinct(TG))
saveRDS(tf_tg_count, "./7 Plot for Analysis/tf_tg_count.rds")
openxlsx::write.xlsx(
  as.data.frame(tf_tg_count),
  "./7 Plot for Analysis/7.3.2_tf_tg_count.xlsx",
  rowNames = FALSE
)

TG_number_per_TF_BarPlot = ggplot(tf_tg_count, aes(x = reorder(TF, -tg_count), y = tg_count)) +
  geom_bar(stat = "identity", fill = "darkorange", alpha = 0.8) +
  geom_text(aes(label = tg_count), vjust = -0.1, size = 3, color = "black") +
  labs(x = "TF", y = "Number of regulated TGs") +
  scale_y_continuous(breaks = seq(0, max(tf_tg_count$tg_count), by = 5)) +
  theme_classic() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = 0.5, lineend = "butt"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
print(TG_number_per_TF_BarPlot)

ggsave("./7 Plot for Analysis/7.3.2_TG_number_per_TF_BarPlot_in_L123_eGRN.png", plot = TG_number_per_TF_BarPlot,
       width = nrow(tf_tg_count) * 0.2, height = 2.8, dpi = 300)


##### 7.3.3 Bar plot: distribution of TGs per peak
# Distribution of TG counts per peak, reflecting regulatory specificity of peaks
peak_tg_count = eGRN_L123 %>%
  group_by(peak) %>%
  summarise(tg_count = n_distinct(TG)) %>%
  group_by(tg_count) %>%
  summarise(peak_count = n())
openxlsx::write.xlsx(
  as.data.frame(peak_tg_count),
  "./7 Plot for Analysis/7.3.3_peak_tg_count.xlsx",
  rowNames = FALSE
)

# Complete missing count levels with zeros
all_tg_counts = data.frame(tg_count = 1:max(3, max(peak_tg_count$tg_count)))
peak_tg_count_complete = all_tg_counts %>%
  left_join(peak_tg_count, by = "tg_count") %>%
  mutate(peak_count = ifelse(is.na(peak_count), 0, peak_count))
TG_number_per_peak_BarPlot = ggplot(peak_tg_count_complete, aes(x = tg_count, y = peak_count)) +
  geom_bar(stat = "identity", fill = "purple", alpha = 0.8) +
  geom_text(aes(label = peak_count), vjust = -0.1, size = 3, color = "black") +
  scale_x_continuous(breaks = 1:max(3, max(peak_tg_count$tg_count))) +
  labs(x = "Number of TGs per peak", y = "Number of peaks") +
  # scale_y_continuous(breaks = seq(0, max(peak_tg_count_complete$peak_count) + 10, by = 45)) +
  theme_classic() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = 0.5, lineend = "butt")
  )
print(TG_number_per_peak_BarPlot)

ggsave("./7 Plot for Analysis/7.3.3_TG_number_per_peak_BarPlot_in_L123_eGRN.png",
       plot = TG_number_per_peak_BarPlot,
       width = max(3, max(peak_tg_count$tg_count)) * 0.8, height = 2.55, dpi = 300)


##### 7.3.4 Bar plot: distribution of TFs per TG
# Distribution of TF counts per TG, indicating potential TF co-regulation
tg_tf_count = eGRN_L123 %>%
  group_by(TG) %>%
  summarise(tf_count = n_distinct(TF)) %>%
  group_by(tf_count) %>%
  summarise(tg_count = n())
openxlsx::write.xlsx(
  as.data.frame(tg_tf_count),
  "./7 Plot for Analysis/7.3.4_tg_tf_count.xlsx",
  rowNames = FALSE
)

# Complete missing count levels with zeros
all_tf_counts = data.frame(tf_count = 1:max(3, max(tg_tf_count$tf_count)))
tg_tf_count_complete = all_tf_counts %>%
  left_join(tg_tf_count, by = "tf_count") %>%
  mutate(tg_count = ifelse(is.na(tg_count), 0, tg_count))
TF_number_per_TG_BarPlot = ggplot(tg_tf_count_complete, aes(x = tf_count, y = tg_count)) +
  geom_bar(stat = "identity", fill = "darkgreen", alpha = 0.8) +
  geom_text(aes(label = tg_count), vjust = -0.1, size = 3, color = "black") +
  scale_x_continuous(breaks = 1:max(3, max(tg_tf_count$tf_count))) +
  labs(x = "Number of TFs per TG", y = "Number of TGs") +
  # scale_y_continuous(breaks = seq(0, max(tg_tf_count_complete$tg_count) + 10, by = 10)) +
  theme_classic() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = 0.5, lineend = "butt")
  )
print(TF_number_per_TG_BarPlot)

ggsave("./7 Plot for Analysis/7.3.4_TF_number_per_TG_BarPlot_in_L123_eGRN.png",
       plot = TF_number_per_TG_BarPlot,
       width = max(3, max(tg_tf_count$tf_count)) * 0.25, height = 2.55, dpi = 300)


##### 7.3.5 Bar plot: distribution of peaks per TG
# Distribution of peak counts per TG, indicating potential enhancer co-regulation
tg_peak_count = eGRN_L123 %>%
  group_by(TG) %>%
  summarise(peak_count = n_distinct(peak)) %>%
  group_by(peak_count) %>%
  summarise(tg_count = n())
openxlsx::write.xlsx(
  as.data.frame(tg_peak_count),
  "./7 Plot for Analysis/7.3.5_tg_peak_count.xlsx",
  rowNames = FALSE
)

all_peak_counts = data.frame(peak_count = 1:max(3, max(tg_peak_count$peak_count)))
tg_peak_count_complete = all_peak_counts %>%
  left_join(tg_peak_count, by = "peak_count") %>%
  mutate(tg_count = ifelse(is.na(tg_count), 0, tg_count))
Peak_number_per_TG_BarPlot = ggplot(tg_peak_count_complete, aes(x = peak_count, y = tg_count)) +
  geom_bar(stat = "identity", fill = "brown", alpha = 0.8) +
  geom_text(aes(label = tg_count), vjust = -0.1, size = 3, color = "black") +
  scale_x_continuous(breaks = 1:max(3, max(tg_peak_count$peak_count))) +
  labs(x = "Number of peaks per TG", y = "Number of TGs") +
  # scale_y_continuous(breaks = seq(0, max(tg_peak_count_complete$tg_count) + 10, by = 30)) +
  theme_classic() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = 0.5, lineend = "butt")
  )
print(Peak_number_per_TG_BarPlot)

ggsave("./7 Plot for Analysis/7.3.5_Peak_number_per_TG_BarPlot_in_L123_eGRN.png",
       plot = Peak_number_per_TG_BarPlot,
       width = max(3, max(tg_peak_count$peak_count)) * 0.6, height = 2.55, dpi = 300)



##### 7.4 Peak-TG Distance Distribution Histogram #####
setwd(Working_directory)
library(dplyr)
library(ggplot2)

# Load disease-associated eGRN (L1-L3 edges)
eGRN_L123 = qs::qread("./7 Plot for Analysis/eGRN_L123.qsave")
eGRN_L123[, "TG"] = toupper(eGRN_L123[, "TG"])

# Load promoter annotation data
results_identify_TGs = qs::qread("./1.2 Data Preprocessing - Identify TGs By Annotation/results_identify_TGs.qsave")
peak_anno_promoter = results_identify_TGs[["peak_anno_promoter"]]
peak_anno_promoter = peak_anno_promoter[, c("peak", "SYMBOL", "distanceToTSS")]
peak_anno_promoter[, "distanceToTSS"] = abs(peak_anno_promoter[, "distanceToTSS"])
peak_anno_promoter[, "SYMBOL"] = toupper(peak_anno_promoter[, "SYMBOL"])
dim(peak_anno_promoter)

# Retain only peak-TG pairs present in the eGRN
peak_anno_promoter = dplyr::semi_join(peak_anno_promoter, eGRN_L123, by = c("peak" = "peak", "SYMBOL" = "TG"))

# Retain only results where distanceToTSS falls between the lower and upper quartiles
# distanceToTSS_q1 = quantile(peak_anno_promoter[, "distanceToTSS"], probs = 0.25)
# distanceToTSS_q3 = quantile(peak_anno_promoter[, "distanceToTSS"], probs = 0.75)

# Restrict peak-TG distance to (0, 15 kb)
distanceToTSS_q1 = 0
distanceToTSS_q3 = 15000
peak_anno_promoter = peak_anno_promoter[
  peak_anno_promoter$distanceToTSS > distanceToTSS_q1 & peak_anno_promoter$distanceToTSS < distanceToTSS_q3, , drop = FALSE
]

# Histogram parameters
n = length(peak_anno_promoter$distanceToTSS) # number of data points
binwidth = 500 # Histogram bin width

# Create data frame for plotting
df = data.frame(value = peak_anno_promoter$distanceToTSS)

# Histogram overlaid with a scaled density curve
Peak_TG_Distance_Distribution_Histogram = ggplot(df, aes(x = value)) +
  geom_histogram( # count-based histogram
    aes(y = after_stat(count)),
    binwidth = binwidth,
    boundary = 0,
    fill = "lightblue",
    color = "black",
    alpha = 0.7
  ) +
  geom_density( # density scaled to histogram counts
    aes(y = after_stat(density) * n * binwidth, color = "Density Curve"),
    linewidth = 1,
    adjust = 1,
    bounds = c(0, 15000),
    show.legend = FALSE
  ) +
  scale_y_continuous(name = "Count", breaks = scales::extended_breaks(n = 5)) +
  scale_x_continuous(breaks = seq(0, 13500, by = 4500)) +
  labs(x = "TG–peak regulatory distance", y = "Count") +
  theme_classic() +
  theme(
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
print(Peak_TG_Distance_Distribution_Histogram)

ggsave("./7 Plot for Analysis/7.4_Peak_TG_Distance_Distribution_Histogram_in_L123_eGRN.png",
       plot = Peak_TG_Distance_Distribution_Histogram,
       width = 3, height = 2, dpi = 300)

### Export Histogram Bin Counts
# Define 500-bp bins over [0, 15 kb]
bin_breaks = seq(distanceToTSS_q1, distanceToTSS_q3, by = binwidth)

# Build interval labels: [0,500], (500,1000], ...
bin_labels = c(
  sprintf("[%d,%d]", bin_breaks[1], bin_breaks[2]),
  sprintf("(%d,%d]", bin_breaks[2:(length(bin_breaks) - 1)], bin_breaks[3:length(bin_breaks)])
)

# Assign distances to bins, retaining all bins including zero counts
bin_factor = cut(df$value, breaks = bin_breaks, labels = bin_labels,
                  include.lowest = TRUE, right = TRUE)
hist_table = data.frame(
  "Distance interval" = bin_labels,
  "TG-peak regulation count" = as.integer(table(factor(bin_factor, levels = bin_labels))),
  check.names = FALSE
)

# Save to the same folder as the figure
openxlsx::write.xlsx(
  hist_table,
  "./7 Plot for Analysis/7.4_Peak_TG_Distance_Distribution_Histogram_in_L123_eGRN.xlsx",
  rowNames = FALSE
)



##### 7.5 Regulon Feature Heatmap #####
setwd(Working_directory)
library(ggplot2)
library(reshape2)
disease_associated_eGRN = qs::qread("./5 Infer GRN/disease_associated_eGRN.qsave")

# Subset full feature table to L1-L3 edges
eGRN_L123 = qs::qread("./7 Plot for Analysis/eGRN_L123.qsave")

x = disease_associated_eGRN[["TF_peaks_TGs"]]
eGRN_L123_feat = x[
  paste(x$TF, x$peak, x$TG, sep = "_") %in% paste(eGRN_L123$TF, eGRN_L123$peak, eGRN_L123$TG, sep = "_"),
  ,
  drop = FALSE
]

# Select and rename feature columns
eGRN_L123_feat = eGRN_L123_feat[, c(
  "TF",
  "G-EA", "P-AA", # "Expression", "Accessibility",
  "G-NZER", "P-NZAR", # "Expression_non_zero_ratio", "Accessibility_non_zero_ratio",
  "FG-theta_i", "FG-theta_p" # "theta_i", "theta_p"
)]
colnames(eGRN_L123_feat) = c(
  "TF",
  "TG Expression Average", "Peak Accessibility Average",
  "TG Expression Non-zero Ratio", "Peak Accessibility Non-zero Ratio",
  "TF–peak Match Strength", "Predicted TF–TG Regulatory Strength"
)

# Max-normalize each feature to [0, 1]
eGRN_L123_feat_norm = eGRN_L123_feat
eGRN_L123_feat_norm[2:ncol(eGRN_L123_feat_norm)] = lapply(eGRN_L123_feat_norm[2:ncol(eGRN_L123_feat_norm)], function(x) {
  x / max(x, na.rm = TRUE)
})
eGRN_L123_feat_norm["TG Expression Average"] = lapply(eGRN_L123_feat_norm["TG Expression Average"], function(x) {
  x / quantile(as.numeric(unlist(x)), probs = 0.8)
})

# Average features per TF
eGRN_norm_mean = aggregate(. ~ TF, data = eGRN_L123_feat_norm, FUN = mean, na.rm = TRUE)
rownames(eGRN_norm_mean) = eGRN_norm_mean$TF
eGRN_mean = aggregate(. ~ TF, data = eGRN_L123_feat, FUN = mean, na.rm = TRUE)
rownames(eGRN_mean) = eGRN_mean$TF

# Append total peak/TG count features
tf_peak_count = readRDS("./7 Plot for Analysis/tf_peak_count.rds")
tf_peak_count_norm = tf_peak_count
tf_peak_count_norm$peak_count = log(tf_peak_count_norm$peak_count) / max(log(tf_peak_count$peak_count * 1.2))

tf_tg_count = readRDS("./7 Plot for Analysis/tf_tg_count.rds")
tf_tg_count_norm = tf_tg_count
tf_tg_count_norm$tg_count = log(tf_tg_count_norm$tg_count) / max(log(tf_tg_count$tg_count * 1.2))

eGRN_mean$`Total Peak Count` = tf_peak_count$peak_count
eGRN_mean$`Total TG Count` = tf_tg_count$tg_count
eGRN_mean = eGRN_mean[, c("TF",
                          "Total Peak Count", "Total TG Count",
                          "TG Expression Average", "Peak Accessibility Average",
                          "TG Expression Non-zero Ratio", "Peak Accessibility Non-zero Ratio",
                          "TF–peak Match Strength", "Predicted TF–TG Regulatory Strength")]

eGRN_norm_mean$`Total Peak Count` = tf_peak_count_norm$peak_count
eGRN_norm_mean$`Total TG Count` = tf_tg_count_norm$tg_count
eGRN_norm_mean = eGRN_norm_mean[, c("TF",
                                    "Total Peak Count", "Total TG Count",
                                    "TG Expression Average", "Peak Accessibility Average",
                                    "TG Expression Non-zero Ratio", "Peak Accessibility Non-zero Ratio",
                                    "TF–peak Match Strength", "Predicted TF–TG Regulatory Strength")]

# Convert to long format
eGRN_norm_long = melt(eGRN_norm_mean, id.vars = "TF", variable.name = "Feature", value.name = "Value_norm")
eGRN_long = melt(eGRN_mean, id.vars = "TF", variable.name = "Feature", value.name = "Value")

# Attach raw mean values as tile labels
eGRN_norm_long$Value = round(eGRN_long$Value, digits = 2)

# Build heatmap: fill shows max-normalized values (red = 1, blue = 0, light yellow = 0.5);
# tile labels show the mean of the raw values
Regulon_Feature_Heatmap = ggplot(eGRN_norm_long, aes(x = Feature, y = TF)) +
  geom_tile(aes(fill = Value_norm), color = "white", size = 0.5) +
  geom_text(
    aes(label = sprintf("%.3g", Value), color = ifelse(Value_norm > 0.76 | Value_norm < 0.24, "white", "black")),
    size = 2, fontface = "bold"
  ) + scale_color_identity() +
  scale_fill_gradient2(
    low = "blue", mid = "lightyellow", high = "red",
    midpoint = 0.5,
    limits = c(0, 1),
    name = NULL,
    n.breaks = 2,
    labels = c("min", "max")
  ) + labs(x = "Feature", y = "Regulon (TF)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1, vjust = 1, size = 7),
    axis.text.y = element_text(size = 9),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 9),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )
print(Regulon_Feature_Heatmap)

ggsave("./7 Plot for Analysis/7.5_Regulon_Feature_Heatmap_in_L123_eGRN.png",
       plot = Regulon_Feature_Heatmap,
       width = 4, height = 3.5, dpi = 300)

### Export Heatmap Data Matrix
# Save matrix: rows = Regulon (TF), columns = features
openxlsx::write.xlsx(
  round(eGRN_mean[, -1, drop = FALSE], digits = 2),
  "./7 Plot for Analysis/7.5_Regulon_Feature_Heatmap_in_L123_eGRN.xlsx",
  rowNames = TRUE
)

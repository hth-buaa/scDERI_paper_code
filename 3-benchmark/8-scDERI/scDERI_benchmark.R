###################################################
#                                                 #
#       Run scDERI on benchmarking dataset        #
#                                                 #
###################################################

# Run in Win11, WSL2, Ubuntu 22.04.5 LTS, R v4.4.1.
# The value of parallel::detectCores() is 78.
# Allocated 224 GB of physical memory to the WSL2 instance.

# Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
# rm(list = ls())

# devtools::install_github("hth-buaa/scDERI")
library(scDERI); library(dplyr); library(qs)
library(TxDb.Hsapiens.UCSC.hg38.knownGene); library(EnsDb.Hsapiens.v86); library(DDRTree); library(BSgenome.Hsapiens.UCSC.hg38)

# Run this set of code three times, each time selecting a dataset (which in this context manifests as annotated vs. non-annotated).
dataset = "A549"
# dataset = "GM12878"
# dataset = "BM"

sys = "wsl2"
# sys = "win11"

if (dataset == "BM" & sys == "wsl2") {
  Working_directory = "/mnt/i/scDERI_Paper/3-benchmark/8-scDERI/10X_hg38_Human_bone_marrow_with_T_cell_depleted/"
} else if (dataset == "BM" & sys == "win11") {
  Working_directory = "I:/scDERI_Paper/3-benchmark/8-scDERI/10X_hg38_Human_bone_marrow_with_T_cell_depleted/"
} else if (dataset == "GM12878" & sys == "wsl2") {
  Working_directory = "/mnt/i/scDERI_Paper/3-benchmark/8-scDERI/SNAREseq2_hg38_GM12878/"
} else if (dataset == "GM12878" & sys == "win11") {
  Working_directory = "I:/scDERI_Paper/3-benchmark/8-scDERI/SNAREseq2_hg38_GM12878/"
} else if (dataset == "A549" & sys == "wsl2") {
  Working_directory = "/mnt/i/scDERI_Paper/3-benchmark/8-scDERI/SNAREseq2_hg38_A549/"
} else if (dataset == "A549" & sys == "win11") {
  Working_directory = "I:/scDERI_Paper/3-benchmark/8-scDERI/SNAREseq2_hg38_A549/"
} else {stop("Error.")}
setwd(Working_directory)





########## Part 0: Loading Data ##########
### Map dataset code to full directory name.
if (dataset == "BM") {
  dataset_full_name = "10X_hg38_Human_bone_marrow_with_T_cell_depleted"
} else if (dataset == "GM12878") {
  dataset_full_name = "SNAREseq2_hg38_GM12878"
} else if (dataset == "A549") {
  dataset_full_name = "SNAREseq2_hg38_A549"
} else {
  stop("Error: unrecognized dataset name.")
}

# Original downsampling: 1,000 cells + top 10k HVGs + top 10k HVPs.
# scDERI requires cell aggregation by pseudotime, which is infeasible with only 1,000 cells.
# Revised strategy: skip cell sampling; keep feature-level filtering (top 10k genes / peaks).

### Copy downsampled data files to working directory.
copy_from = c(
  paste0("../../0-downsampling/", dataset_full_name, "/", dataset_full_name, ".qsave"),
  paste0("../../0-downsampling/", dataset_full_name, "/", dataset_full_name, "_fragments.tsv.gz"),
  paste0("../../0-downsampling/", dataset_full_name, "/", dataset_full_name, "_fragments.tsv.gz.tbi")
)
file.copy(from = copy_from, to = Working_directory, overwrite = TRUE)

### Load Seurat object and report basic dimensions.
data = qs::qread(paste0("./", dataset_full_name, ".qsave"))
ncell = ncol(data)
ngene = nrow(data[['RNA']])
npeak = nrow(data[['ATAC']])
message("----> Cells: ", ncell, "\n",
        "----> Genes: ", ngene, "\n",
        "----> Peaks: ", npeak, "\n")

### Gene filtering: remove non-gene symbols and mitochondrial genes.
genes = setdiff(rownames(data[['RNA']]),
                grep("(^AC|^AL|^LINC)[0-9]+", rownames(data[['RNA']]), value = TRUE))
message("----> ", length(genes), " genes after filtering primer/pseudogene names")

genes = grep("MT-", genes, invert = TRUE, value = TRUE)
message("----> ", length(genes), " genes after filtering mitochondrial genes")

### Select top 10,000 highly variable genes (vst method).
Seurat::DefaultAssay(data) = "RNA"
data = Seurat::NormalizeData(data)
data = Seurat::FindVariableFeatures(data, selection.method = "vst", nfeatures = 10000)

# Extract variable feature metadata (slot name varies by Seurat version)
if ("meta.features" %in% slotNames(data[["RNA"]])) {
  rna.meta = data[['RNA']]@meta.features
  rna.meta = rna.meta[order(rna.meta$vst.variance.standardized, decreasing = TRUE), ]
} else if ("meta.data" %in% slotNames(data[["RNA"]])) {
  rna.meta = data[['RNA']]@meta.data
  rownames(rna.meta) = rownames(data)
  rna.meta = rna.meta[order(rna.meta$vf_vst_counts_variance.standardized, decreasing = TRUE), ]
} else {
  stop("Error: no meta feature slot found in RNA assay.")
}

library(dplyr)
genes = rna.meta[rownames(rna.meta) %in% genes, ] %>% head(n = 10000) %>% rownames
message("----> ", length(genes), " genes after selecting highly variable features")

### Filter peaks on standard chromosomes (chr1-22, chrX, chrY).
peaks = rownames(data[['ATAC']])
peaks = grep("chr([0-9]|X|Y)", peaks, value = TRUE)
message("----> ", length(peaks), " peaks after filtering non-standard chromosomes")

### Select top 10,000 most accessible peaks.
Seurat::DefaultAssay(data) = "ATAC"
data = Signac::FindTopFeatures(data)

atac.meta = data[['ATAC']]@meta.features
atac.meta = atac.meta[order(atac.meta$percentile, decreasing = TRUE), ]
peaks = atac.meta[rownames(atac.meta) %in% peaks, ] %>% head(n = 10000) %>% rownames
message("----> ", length(peaks), " peaks after selecting top variable features")

### Subset Seurat object to filtered genes and peaks.
data = S4Vectors::subset(data, features = c(genes, peaks))
message("----> ", dataset, ": ", nrow(data[['RNA']]), " genes, ",
        nrow(data[['ATAC']]), " peaks, and ", ncol(data),
        " cells retained after downsampling")
# ----> A549: 10000 genes, 10000 peaks, and 1088 cells retained after downsampling
# ----> GM12878: 10000 genes, 10000 peaks, and 8876 cells retained after downsampling
# ----> BM: 10000 genes, 10000 peaks, and 7418 cells retained after downsampling

### Assign dataset label and split into separate RNA / ATAC objects.
if (dataset == "BM") {
  data@meta.data$seurat_annotations = "BM"
} else if (dataset == "GM12878") {
  data@meta.data$seurat_annotations = "GM12878"
} else if (dataset == "A549") {
  data@meta.data$seurat_annotations = "A549"
} else {
  stop("Error: unrecognized dataset name.")
}

# Extract RNA assay as a standalone Seurat object
scRNAseq = SeuratObject::CreateSeuratObject(
  counts = SeuratObject::GetAssayData(data, assay = "RNA", layer = "counts"),
  assay = "RNA",
  meta.data = data@meta.data
)

# Extract ATAC assay as a standalone Seurat object
data@active.assay = "ATAC"
data@assays[["RNA"]] = NULL
scATACseq = data
rm(data)

if (dataset == "GM12878") {
  dedup_rename_cells = function(obj) {
    new_names = sub(".*_", "", colnames(obj))
    keep = !duplicated(new_names)
    obj = obj[, keep]
    colnames(obj) = new_names[keep]
    return(obj)
  }
  scRNAseq = dedup_rename_cells(scRNAseq)
  scATACseq = dedup_rename_cells(scATACseq)
}

### Save split objects to disk.
# saveRDS(scRNAseq, file = "./scRNAseq.rds")
# saveRDS(scATACseq, file = "./scATACseq.rds")
qs::qsave(scRNAseq, file = "./scRNAseq.qsave")
qs::qsave(scATACseq, file = "./scATACseq.qsave")




########## Part 1: Data Preprocessing ##########
##### 1.1 Identifying TFs #####
setwd(Working_directory)
scRNAseq = qs::qread("./scRNAseq.qsave")
options(timeout = 999999999)
if (dataset == "BM") {
  excluded_cell_types_in_priors = c("bone marrow", "NCI-H929", "BLaER1", "OCI-LY1", "OCI-LY3", "erythroblast", "SH-SY5Y", "NB4")
} else if (dataset == "GM12878") {
  excluded_cell_types_in_priors = "GM12878"
} else if (dataset == "A549") {
  excluded_cell_types_in_priors = "A549"
} else {stop("Error.")}
monitor = peakRAM::peakRAM({
  TFs = get_TFs_from_priors(
    scRNAseq = scRNAseq, species = "Homo sapiens", genome_assembly = "hg38",
    min_nonzero_freq = 0.032, min_nonzero_mean = 1.104, min_overall_mean = 0.046,
    min_nonzero_freq_quantile = 0.387, min_nonzero_mean_quantile = 0.382, min_overall_mean_quantile = 0.421,
    delay_time = 0.05, excluded_cell_types_in_priors = excluded_cell_types_in_priors, ignore_case_when_filtering_priors = TRUE,
    prior_min_TG_number = 5900, prior_top_n_keep = 30000, prior_top_q_keep = 1
  )
})
message(dataset, ": Get ", length(TFs), " TFs, ", print_run_stats(monitor))
saveRDS(TFs, file = "./1.1 Data Preprocessing - Get TFs From Database/TFs.rds")
# A549: Get 87 TFs, running: 0 hours 1 mins 18 secs, peak RAM: 0 GB 385 MB
# GM12878: Get 79 TFs, running: 0 hours 1 mins 35 secs, peak RAM: 0 GB 666 MB
# BM: Get 94 TFs, running: 0 hours 1 mins 29 secs, peak RAM: 0 GB 685 MB



##### 1.2 Annotating Chromatin Fragments and Identifying TGs #####
##### 1.2.1 Annotating Chromatin Fragments and Identifying TGs
setwd(Working_directory)
scRNAseq = qs::qread("./scRNAseq.qsave")
scATACseq = qs::qread("./scATACseq.qsave")
library(dplyr)
monitor = peakRAM::peakRAM({
  results_identify_TGs_all = identify_TGs(
    genome_for_anno = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene,
    scRNAseq = scRNAseq, scATACseq = scATACseq, annoDb = "org.Hs.eg.db", flankDistance = 500000
  )
})
message(dataset, ": ", print_run_stats(monitor))
saveRDS(results_identify_TGs_all, file = "./1.2 Data Preprocessing - Identify TGs By Annotation/results_identify_TGs_all.rds")
# A549: running: 0 hours 1 mins 15 secs, peak RAM: 0 GB 763 MB
# GM12878: running: 0 hours 0 mins 52 secs, peak RAM: 0 GB 748 MB
# BM: running: 0 hours 0 mins 54 secs, peak RAM: 1 GB 249 MB


##### 1.2.2 Filtering Annotation Results
setwd(Working_directory)
results_identify_TGs_all = readRDS("./1.2 Data Preprocessing - Identify TGs By Annotation/results_identify_TGs_all.rds")
scRNAseq = qs::qread("./scRNAseq.qsave")
scATACseq = qs::qread("./scATACseq.qsave")
monitor = peakRAM::peakRAM({
  results_identify_TGs = filter_peak_annotations(
    results_identify_TGs = results_identify_TGs_all, scRNAseq = scRNAseq, scATACseq = scATACseq, forced_feature_calculation = TRUE,
    distance_points = base::c(1, 500, 500000), min_quantile_points = base::c(0.35, 0.55, 1.55),
    distance_quantile_fit_start = base::list(a = 0.006667, b = 1.659167, c = 0.6003),
    one_peak_at_least_per_gene = TRUE, gene_quality_factor_min = 0.9,
    weight_gene_nonzero_freq = 0.7, weight_gene_nonzero_mean = 0.2, weight_gene_overall_mean = 0.1,
    min_peak_nonzero_freq_quantile_scale = 0.4, min_peak_nonzero_mean_scale = 1.2, min_peak_overall_mean_quantile_scale = 0.1
  )
})
message(dataset, ": ", nrow(results_identify_TGs[["peak_anno_promoter_all"]]), "/", length(unique(results_identify_TGs[["peak_anno_promoter"]][["peak"]])), " peaks, ",
        length(unique(results_identify_TGs[["peak_anno_promoter_all"]][["SYMBOL"]])), "/", length(results_identify_TGs[["TGs"]]), " TGs, ",
        nrow(results_identify_TGs[["peak_anno_promoter_all"]]), "/", nrow(unique(results_identify_TGs[[1]][,c("peak","SYMBOL")])), " peak-TG, peaks/TGs = ",
        round(length(unique(results_identify_TGs[["peak_anno_promoter"]][["peak"]]))/length(results_identify_TGs[["TGs"]]), digits = 2), ", ",
        print_run_stats(monitor))
saveRDS(results_identify_TGs, file = "./1.2 Data Preprocessing - Identify TGs By Annotation/results_identify_TGs.rds")
# A549: 3992/2176 peaks, 2297/3413 TGs, 3992/5099 peak-TG, peaks/TGs = 0.64, running: 0 hours 0 mins 12 secs, peak RAM: 0 GB 246 MB
# GM12878: 3714/2180 peaks, 2342/3655 TGs, 3714/5135 peak-TG, peaks/TGs = 0.6, running: 0 hours 0 mins 10 secs, peak RAM: 1 GB 112 MB
# BM: 4440/3411 peaks, 3598/4478 TGs, 4440/7182 peak-TG, peaks/TGs = 0.76, running: 0 hours 0 mins 14 secs, peak RAM: 1 GB 189 MB



##### 1.3 Generating Counts Matrices for TGs/TFs Expression and TGs Activity (Seurat Object) #####
setwd(Working_directory)
results_identify_TGs = readRDS("./1.2 Data Preprocessing - Identify TGs By Annotation/results_identify_TGs.rds")
TFs = readRDS("./1.1 Data Preprocessing - Get TFs From Database/TFs.rds")
scRNAseq = qs::qread("./scRNAseq.qsave")
scATACseq = qs::qread("./scATACseq.qsave")
if (dataset == "BM") {
  dataset_full_name = "10X_hg38_Human_bone_marrow_with_T_cell_depleted"
} else if (dataset == "GM12878") {
  dataset_full_name = "SNAREseq2_hg38_GM12878"
} else if (dataset == "A549") {
  dataset_full_name = "SNAREseq2_hg38_A549"
} else {stop("Error.")}
monitor = peakRAM::peakRAM({
  Basic_Info = get_expression_and_activity_matrix(
    TGs = results_identify_TGs$TGs, TFs = TFs, scRNAseq = scRNAseq, scATACseq = scATACseq,
    Ensdb = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86, flanking_region = 500000,
    path = paste0(Working_directory, dataset_full_name, "_fragments.tsv.gz")
  )
})
message(dataset, ": ", length(Basic_Info[["TFs"]]), " TFs, ", length(Basic_Info[["TGs"]]), " TGs, ",
        length(Basic_Info[["scRNAseq_for_GRN"]]@meta.data$seurat_annotations), " cells, ", print_run_stats(monitor))
qs::qsave(Basic_Info, file = "./1.3 Data Preprocessing - Get Expression And Activity Matrix/Basic_Info.qsave")
# A549: 85 TFs, 3147 TGs, 1088 cells, running: 0 hours 2 mins 34 secs, peak RAM: 1 GB 573 MB
# GM12878: 77 TFs, 3363 TGs, 7426 cells, running: 0 hours 5 mins 55 secs, peak RAM: 2 GB 39 MB
# BM: 93 TFs, 4159 TGs, 7418 cells, running: 1 hours 17 mins 23 secs, peak RAM: 25 GB 981 MB




########## Part 2: Inferring Pseudotime, Branching, and Grouping Cells ##########
##### 2.1 Pseudotime Analysis and Cell Branch Partitioning #####
##### 2.1.1 Obtaining Relevant Data of Cell Types of Interest
setwd(Working_directory)
Basic_Info = qs::qread("./1.3 Data Preprocessing - Get Expression And Activity Matrix/Basic_Info.qsave")
monitor = peakRAM::peakRAM({
  interest_cell_type_data = get_interest_cell_type_data(
    interest_cell_type = dataset, Basic_Info = Basic_Info,
    alpha_gene = 1, alpha_cell = 1, ncores = parallel::detectCores() - 1
  )
})
message(dataset, ": ", print_run_stats(monitor))
qs::qsave(interest_cell_type_data, file = "./2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment/interest_cell_type_data.qsave")
# A549: running: 0 hours 0 mins 3 secs, peak RAM: 0 GB 241 MB
# GM12878: running: 0 hours 0 mins 10 secs, peak RAM: 0 GB 852 MB
# BM: running: 0 hours 0 mins 14 secs, peak RAM: 3 GB 539 MB


##### 2.1.2 Performing Pseudotemporal Ordering and Branching Analysis of Cell Types of Interest Based on Transcriptomic Data
setwd(Working_directory)
library(DDRTree)
if (packageVersion("igraph") != "2.0.3") {
  remove.packages("igraph")
  install_version("igraph", version = "2.0.3")
}
interest_cell_type_data = qs::qread("./2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment/interest_cell_type_data.qsave")
monitor = peakRAM::peakRAM({
  interest_cell_type_Branches = order_pseudotime_and_divide_branches(
    interest_cell_type_data = interest_cell_type_data,
    auto_param_selection = FALSE, ncenter = 4 # ensure a single branch
  )
})
1
1
message(dataset, ": ", print_run_stats(monitor))
saveRDS(interest_cell_type_Branches, file = "./2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment/interest_cell_type_Branches.rds")
# A549: running: 0 hours 0 mins 47 secs, peak RAM: 2 GB 244 MB
# GM12878: running: 0 hours 6 mins 37 secs, peak RAM: 4 GB 208 MB
# BM: running: 0 hours 8 mins 16 secs, peak RAM: 4 GB 600 MB



##### 2.2 Cell Grouping and Group-Wise Information Specification #####
# Principles for Parameter Settings in Gene Filtering and Cell Grouping:
# First priority: n_min must be at least 2.
# Second priority: Given n_min ≥ 2, aim to create exactly 8 groups whenever possible; if not feasible, get as close to 8 groups as possible.
# Third priority: While targeting ~8 groups, set n_min as close to 6 or 7 as possible.
# Fourth priority: Once the group count and n_min are fixed, minimize the exclusion of TG.


##### 2.2.1 Retrieving Pseudotemporal Expression/Activity Vectors for Each Gene
setwd(Working_directory)
interest_cell_type_Branches = readRDS("./2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment/interest_cell_type_Branches.rds")
interest_cell_type_data = qs::qread("./2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment/interest_cell_type_data.qsave")
sparsity_cut = base::list()
sparsity_cut[[dataset]] = base::as.data.frame(base::matrix(
  0, nrow = 3, ncol = base::length(interest_cell_type_Branches[[dataset]]),
  dimnames = base::list(c("TFE", "TGA", "TGE"), base::paste0("branch", base::seq_along(interest_cell_type_Branches[[dataset]])))
))
sparsity_cut[[dataset]]["TFE", "branch1"] = 0
if (dataset == "A549") {
  sparsity_cut[[dataset]]["TGE", "branch1"] = 0.08
} else if (dataset == "GM12878") {
  sparsity_cut[[dataset]]["TGE", "branch1"] = 0.02
} else if (dataset == "BM") {
  sparsity_cut[[dataset]]["TGE", "branch1"] = 0.09
} else {stop("Error.")}
sparsity_cut[[dataset]]["TGA", "branch1"] = sparsity_cut[[dataset]]["TGE", "branch1"]
monitor = peakRAM::peakRAM({
  interest_cell_type_genes_pseudotime_info = suppressMessages(get_genes_pseudotime_info(
    interest_cell_type_Branches = interest_cell_type_Branches, interest_cell_type_data = interest_cell_type_data,
    alpha_cell = 1, alpha_TF = 0, alpha_TG = 0, ncores = parallel::detectCores() - 1, sparsity_cut = sparsity_cut
  ))
})
message(dataset, ": ", print_run_stats(monitor))
saveRDS(interest_cell_type_genes_pseudotime_info, file = "./2.2 Data Processing - Cell Grouping/interest_cell_type_genes_pseudotime_info.rds")
# A549: running: 0 hours 0 mins 5 secs, peak RAM: 0 GB 260 MB
# GM12878: running: 0 hours 0 mins 28 secs, peak RAM: 1 GB 869 MB
# BM: running: 0 hours 0 mins 52 secs, peak RAM: 6 GB 907 MB

# Check minimum nonzero frequency
for (i in seq_along(interest_cell_type_genes_pseudotime_info)) {
  for (j in 1:3) {
    for (k in seq_along(interest_cell_type_genes_pseudotime_info[[i]][[j]])) {
      df = interest_cell_type_genes_pseudotime_info[[i]][[j]][[k]]
      nonzero_freq = rowSums(df != 0) / ncol(df)
      print(head(sort(nonzero_freq)))
      min_nonzero_freq = min(nonzero_freq)
      cat("i: ", i, ", j: ", j, ", k: ", k, ", min_nonzero_freq: ", min_nonzero_freq, "\n")
    }
  }
}; rm(df, nonzero_freq)
# A549
#      PRDM1       JUND       E2F8      HNF1A      NFKB2      HNF4G
# 0.04503676 0.04595588 0.06158088 0.06158088 0.06433824 0.07169118
# i:  1 , j:  1 , k:  1 , min_nonzero_freq:  0.04503676
#       ZW10       INIP     LPCAT2     ZNF708       DCTD   C6orf141
# 0.08272059 0.08363971 0.08455882 0.08547794 0.08915441 0.09099265
# i:  1 , j:  2 , k:  1 , min_nonzero_freq:  0.08272059
#     DHCR24       RFX5      CDH23       JAG2    TMEM160     KCNJ14
# 0.08088235 0.08088235 0.08088235 0.08088235 0.08088235 0.08088235
# i:  1 , j:  3 , k:  1 , min_nonzero_freq:  0.08088235

# GM12878
#       KLF3      NR4A1     TCF7L2      GRHL2       KLF9       E2F1
# 0.04134123 0.04349583 0.04443846 0.05022893 0.05965527 0.07029356
# i:  1 , j:  1 , k:  1 , min_nonzero_freq:  0.04134123
#     EPHA3     ITGA1     MEIS1  ARHGAP42    TFAP2B      IPMK
# 0.1014005 0.1415298 0.1624024 0.1675195 0.1747913 0.1813897
# i:  1 , j:  2 , k:  1 , min_nonzero_freq:  0.1014005
#       ICOS      PRSS2     SCUBE2     THEMIS      TESK1       TLE1
# 0.02006464 0.02006464 0.02006464 0.02019930 0.02019930 0.02019930
# i:  1 , j:  3 , k:  1 , min_nonzero_freq:  0.02006464

# BM
#      GATA2       BATF      PRDM1       E2F1      GATA1       THRB
# 0.03248854 0.04826099 0.05109194 0.05244001 0.05769749 0.05971960
# i:  1 , j:  1 , k:  1 , min_nonzero_freq:  0.03248854
#     PLRG1      GUF1   TBC1D32    SRFBP1    PACRGL     GCFC2
# 0.6583985 0.6842815 0.6881909 0.6915611 0.6942572 0.6947964
# i:  1 , j:  2 , k:  1 , min_nonzero_freq:  0.6583985
#     LZTFL1     FAM76A      NDRG3     SLC9B2    TIMMDC1       OLR1
# 0.09005123 0.09032084 0.09032084 0.09032084 0.09059046 0.09099488
# i:  1 , j:  3 , k:  1 , min_nonzero_freq:  0.09005123

results_identify_TGs = readRDS("./1.2 Data Preprocessing - Identify TGs By Annotation/results_identify_TGs.rds")
peak_anno_promoter = results_identify_TGs[["peak_anno_promoter"]]
peak_anno_promoter = peak_anno_promoter[peak_anno_promoter$SYMBOL %in% Reduce(union, interest_cell_type_genes_pseudotime_info[[1]][["Branches_TGs"]]), ]
message(dataset, ": ", sapply(interest_cell_type_genes_pseudotime_info[[1]][["Branches_TFs"]], length), " TFs, ",
        length(unique(peak_anno_promoter[["peak"]])), " peaks, ",
        length(unique(results_identify_TGs[["peak_anno_promoter_all"]][["SYMBOL"]])), "/",
        length(unique(peak_anno_promoter[["SYMBOL"]])), " TGs, peaks/TGs = ",
        round(length(unique(peak_anno_promoter[["peak"]]))/length(unique(peak_anno_promoter[["SYMBOL"]])), digits = 2))
# A549: 84 TFs, 1358 peaks, 2297/1341 TGs, peaks/TGs = 1.01
# GM12878: 76 TFs, 1662 peaks, 2342/1824 TGs, peaks/TGs = 0.91
# BM: 92 TFs, 2120 peaks, 3598/2022 TGs, peaks/TGs = 1.05


##### 2.2.2 Grouping Cells Based on Pseudotime Values and Counts of Valid (Non-Zero) Expression and Activity Values to Calculate TFE_T, TGA_T, and TGE_T per Group
setwd(Working_directory)
interest_cell_type_genes_pseudotime_info = readRDS("./2.2 Data Processing - Cell Grouping/interest_cell_type_genes_pseudotime_info.rds")
n_cell = interest_cell_type_genes_pseudotime_info[[1]][["Branches_n_cell"]][[1]]
points_x_for_fitting_nc_and_nmin = round(c(100, n_cell * 1.1, n_cell * 2.2))
if (dataset == "BM") {
  points_y_for_fitting_nc_and_nmin = c(6, 7, 8)
} else if (dataset == "GM12878") {
  points_y_for_fitting_nc_and_nmin = c(5, 6, 7)
} else if (dataset == "A549") {
  points_y_for_fitting_nc_and_nmin = c(1, 2, 3)
} else {stop("Error.")}
monitor = peakRAM::peakRAM({
  interest_cell_type_group = cell_grouping(
    interest_cell_type_genes_pseudotime_info = interest_cell_type_genes_pseudotime_info,
    points_x_for_fitting_nc_and_nmin = points_x_for_fitting_nc_and_nmin,
    points_y_for_fitting_nc_and_nmin = points_y_for_fitting_nc_and_nmin,
    ncores = parallel::detectCores() - 1
  )
})
message(dataset, ": ", print_run_stats(monitor))
saveRDS(interest_cell_type_group, file = "./2.2 Data Processing - Cell Grouping/interest_cell_type_group.rds")
# A549: running: 0 hours 0 mins 6 secs, peak RAM: 0 GB 274 MB
# GM12878: running: 0 hours 4 mins 37 secs, peak RAM: 2 GB 530 MB
# BM: running: 0 hours 5 mins 39 secs, peak RAM: 4 GB 709 MB

message(dataset, ": ", paste(unlist(interest_cell_type_group[[1]][["Branches_n_group"]]), collapse = ", "))
# A549: 8
# GM12878: 8
# BM: 8

message(dataset, ": ", paste(unlist(interest_cell_type_group[[1]][["n_min"]]), collapse = ", "))
# A549: 2
# GM12878: 6
# BM: 7

message(dataset, ": ", paste(unlist(interest_cell_type_genes_pseudotime_info[[1]][["Branches_n_cell"]]), collapse = ", "))
# A549: 1088
# GM12878: 7426
# BM: 7418




########## Part 3: Constructing Initial GRNs from TFBS PWMs ##########
##### 3.1 Constructing an Initial GRN for Comprehensive TG-TF Pairs #####
##### 3.1.1 Constructing an Initial GRN for Comprehensive TG-TF Pairs Based on Base PWM of TF Binding Sites
setwd(Working_directory)
interest_cell_type_genes_pseudotime_info = readRDS("./2.2 Data Processing - Cell Grouping/interest_cell_type_genes_pseudotime_info.rds")
results_identify_TGs = base::readRDS("./1.2 Data Preprocessing - Identify TGs By Annotation/results_identify_TGs.rds")
suppressMessages(library(dplyr))
monitor = peakRAM::peakRAM({
  interest_cell_type_iGRN_all = get_iGRN_by_TFBS_pwm_by_JASPAR2024(
    interest_cell_type_genes_pseudotime_info = interest_cell_type_genes_pseudotime_info,
    results_identify_TGs = results_identify_TGs, genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
    min_score_for_matchPWM = "0%", species = "Homo sapiens", ncores = parallel::detectCores() - 1
  )
})
message(dataset, ": ", ncol(interest_cell_type_iGRN_all[[1]]), " TFs, ", nrow(interest_cell_type_iGRN_all[[1]]), " TGs, ",
        sum(interest_cell_type_iGRN_all[[1]] > 0), " pairs, ", print_run_stats(monitor))
saveRDS(interest_cell_type_iGRN_all, file = "./3 get iGRN/interest_cell_type_iGRN_all.rds")
# A549: 84 TFs, 1341 TGs, 112644 pairs, running: 0 hours 7 mins 23 secs, peak RAM: 0 GB 109 MB
# GM12878: 76 TFs, 1824 TGs, 138624 pairs, running: 0 hours 7 mins 24 secs, peak RAM: 0 GB 117 MB
# BM: 92 TFs, 2022 TGs, 186024 pairs, running: 0 hours 9 mins 42 secs, peak RAM: 0 GB 149 MB


##### 3.1.2  Filtering  the Initial GRN by Chip-Atlas Priors
setwd(Working_directory)
interest_cell_type_iGRN_all = readRDS("./3 get iGRN/interest_cell_type_iGRN_all.rds")

results_identify_TGs = readRDS("./1.2 Data Preprocessing - Identify TGs By Annotation/results_identify_TGs.rds")
peak_anno_promoter = results_identify_TGs[["peak_anno_promoter"]]
peak_anno_promoter = peak_anno_promoter[peak_anno_promoter$SYMBOL %in% rownames(interest_cell_type_iGRN_all[[1]]), ]
peak_TG_ratio = length(unique(peak_anno_promoter[["peak"]]))/length(unique(peak_anno_promoter[["SYMBOL"]]))
message(dataset, ": peak_TG_ratio = ", round(peak_TG_ratio, digits = 4))
# A549: peak_TG_ratio = 1.0127
# GM12878: peak_TG_ratio = 0.9112
# BM: peak_TG_ratio = 1.0485

expected_peak_number_per_TF = 50
theta_i_top_quantile = expected_peak_number_per_TF/peak_TG_ratio/nrow(interest_cell_type_iGRN_all[[1]])
theta_i_top_n = round(expected_peak_number_per_TF/peak_TG_ratio)

monitor = peakRAM::peakRAM({
  interest_cell_type_iGRN = filter_iGRN_by_ChIP_Atlas_priors(
    interest_cell_type_iGRN = interest_cell_type_iGRN_all,
    theta_i_top_quantile = theta_i_top_quantile, theta_i_top_n = theta_i_top_n, TG_rescue_number = 5
  )
})
message(dataset, ": top_q=", round(theta_i_top_quantile, digits = 4),
        ", top_n=", theta_i_top_n, ", ",
        ncol(interest_cell_type_iGRN[[1]]), " TFs, ", nrow(interest_cell_type_iGRN[[1]]), " TGs, ",
        sum(interest_cell_type_iGRN[[1]] > 0), " pairs, ",
        round(sum(interest_cell_type_iGRN[[1]] > 0) / ncol(interest_cell_type_iGRN[[1]]), digits = 2), " TGs per TF, ", print_run_stats(monitor))
saveRDS(interest_cell_type_iGRN, file = "./3 get iGRN/interest_cell_type_iGRN.rds")
# A549: top_q=0.0368, top_n=49, 84 TFs, 1319 TGs, 7977 pairs, 94.96 TGs per TF, running: 0 hours 0 mins 6 secs, peak RAM: 0 GB 204 MB
# GM12878: top_q=0.0301, top_n=55, 76 TFs, 1797 TGs, 9898 pairs, 130.24 TGs per TF, running: 0 hours 0 mins 9 secs, peak RAM: 0 GB 268 MB
# BM: top_q=0.0236, top_n=48, 92 TFs, 1978 TGs, 11016 pairs, 119.74 TGs per TF, running: 0 hours 0 mins 9 secs, peak RAM: 0 GB 269 MB

# Examine the distribution of the number of TGs associated with each TF
res = t(apply(interest_cell_type_iGRN[[1]], 2, function(x) {
  nz = x[x != 0]
  c(nonzero_freq = length(nz) / length(x),
    nonzero_count = length(nz),
    nonzero_min = if(length(nz)) min(nz) else NA)
}))
summary(res[,2])
# A549
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 49.00   53.00   74.50   94.96  116.50  272.00

# GM12878
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 55.00   90.75  112.00  130.24  170.25  370.00

# BM
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 48.00   55.75   76.50  119.74  134.75  569.00



##### 3.2 Setting the Regulatory Threshold tau #####
setwd(Working_directory)
interest_cell_type_iGRN = readRDS("./3 get iGRN/interest_cell_type_iGRN.rds")
monitor = peakRAM::peakRAM({
  interest_cell_type_tau = get_cell_type_tau(interest_cell_type_iGRN = interest_cell_type_iGRN)
})
message(dataset, ": tau is ", round(interest_cell_type_tau[1], digits = 4), ", ", print_run_stats(monitor))
saveRDS(interest_cell_type_tau, file = "./3 get iGRN/interest_cell_type_tau.rds")
# A549: tau is 0.2459, running: 0 hours 0 mins 0 secs, peak RAM: 0 GB 3 MB
# GM12878: tau is 0.1831, running: 0 hours 0 mins 0 secs, peak RAM: 0 GB 4 MB
# BM: tau is 0.1851, running: 0 hours 0 mins 0 secs, peak RAM: 0 GB 4 MB



##### 3.3 Constructing the iGRN for Each Branch #####
setwd(Working_directory)
interest_cell_type_iGRN = readRDS("./3 get iGRN/interest_cell_type_iGRN.rds")
interest_cell_type_group = readRDS("./2.2 Data Processing - Cell Grouping/interest_cell_type_group.rds")
monitor = peakRAM::peakRAM({
  interest_cell_type_branch_iGRN = get_branch_iGRN(
    interest_cell_type_iGRN = interest_cell_type_iGRN,
    interest_cell_type_group = interest_cell_type_group,
    ncores = parallel::detectCores() - 1
  )
})
message(dataset, ": ", paste(unlist(lapply(interest_cell_type_branch_iGRN[[1]], ncol)), collapse = "/"), " TFs, ",
        paste(unlist(lapply(interest_cell_type_branch_iGRN[[1]], nrow)), collapse = "/"), " TGs, ",
        paste(sapply(interest_cell_type_branch_iGRN[[1]], function(df) sum(df > 0, na.rm = TRUE)), collapse = "/"), " pairs, ", print_run_stats(monitor))
saveRDS(interest_cell_type_branch_iGRN, file = "./3 get iGRN/interest_cell_type_branch_iGRN.rds")
# A549: 84 TFs, 1319 TGs, 7977 pairs, running: 0 hours 0 mins 10 secs, peak RAM: 0 GB 506 MB
# GM12878: 76 TFs, 1797 TGs, 9898 pairs, running: 0 hours 0 mins 16 secs, peak RAM: 0 GB 506 MB
# BM: 92 TFs, 1978 TGs, 11016 pairs, running: 0 hours 0 mins 20 secs, peak RAM: 0 GB 506 MB




########## Part 4: Developing DERI-Based GRN Regression Predictive Models for Each Branch ##########
##### 4.1 Partitioning Training Sets #####
setwd(Working_directory)
interest_cell_type_branch_iGRN = readRDS("./3 get iGRN/interest_cell_type_branch_iGRN.rds")
interest_cell_type_tau = base::readRDS("./3 get iGRN/interest_cell_type_tau.rds")
interest_cell_type_group = readRDS("./2.2 Data Processing - Cell Grouping/interest_cell_type_group.rds")
monitor = peakRAM::peakRAM({
  interest_cell_type_branch_training_set = suppressMessages(select_training_set(
    interest_cell_type_branch_iGRN = interest_cell_type_branch_iGRN,
    interest_cell_type_tau  = interest_cell_type_tau, interest_cell_type_group = interest_cell_type_group,
    top_n_for_max_theta_in_each_regulon = 1, max_theta_retention_rate_in_tops = 1, max_theta_retention_rate_in_all = 1,
    retention_number_of_tops = 200000, refine = TRUE, keep_all_TF = TRUE, top_ratio_in_refine = 1, ncores = 1
  ))
})
message(dataset, ": ", paste(lapply(interest_cell_type_branch_training_set[[1]], nrow), collapse = ","), " training TF-TG pairs, ", print_run_stats(monitor))
saveRDS(interest_cell_type_branch_training_set, file = "./4.1 Build Prediction Model - Select Training Set/interest_cell_type_branch_training_set.rds")
# A549: 85,85 training TF-TG pairs, running: 0 hours 0 mins 6 secs, peak RAM: 0 GB 466 MB
# GM12878: 76,76 training TF-TG pairs, running: 0 hours 0 mins 9 secs, peak RAM: 0 GB 507 MB
# BM: 92,92 training TF-TG pairs, running: 0 hours 0 mins 11 secs, peak RAM: 0 GB 506 MB



##### 4.2 Build Predictive Model for Each Branch #####
##### 4.2.1 Setting Initial Values, Lower Bounds, and Upper Bounds for Parameters
setwd(Working_directory)
interest_cell_type_branch_training_set = readRDS("./4.1 Build Prediction Model - Select Training Set/interest_cell_type_branch_training_set.rds")
monitor = peakRAM::peakRAM({
  interest_cell_type_branch_init_params = suppressMessages(set_init_params(interest_cell_type_branch_training_set = interest_cell_type_branch_training_set))
})
message(dataset, ": ", print_run_stats(monitor))
saveRDS(interest_cell_type_branch_init_params, file = "./4.2 BUild Prediction Model/interest_cell_type_branch_init_params.rds")
# A549: running: 0 hours 0 mins 1 secs, peak RAM: 0 GB 10 MB
# GM12878: running: 0 hours 0 mins 1 secs, peak RAM: 0 GB 10 MB
# BM: running: 0 hours 0 mins 0 secs, peak RAM: 0 GB 6 MB



##### 4.2.2 Training models
setwd(Working_directory)
interest_cell_type_branch_training_set = readRDS("./4.1 Build Prediction Model - Select Training Set/interest_cell_type_branch_training_set.rds")
interest_cell_type_branch_init_params = readRDS("./4.2 BUild Prediction Model/interest_cell_type_branch_init_params.rds")
interest_cell_type_group = readRDS("./2.2 Data Processing - Cell Grouping/interest_cell_type_group.rds")
interest_cell_type_tau = readRDS("./3 get iGRN/interest_cell_type_tau.rds")
interest_cell_type_genes_pseudotime_info = readRDS("./2.2 Data Processing - Cell Grouping/interest_cell_type_genes_pseudotime_info.rds")
monitor = peakRAM::peakRAM({
  interest_cell_type_branch_model_train_all = model_train(
    interest_cell_type_branch_training_set = interest_cell_type_branch_training_set,
    interest_cell_type_branch_init_params = interest_cell_type_branch_init_params,
    interest_cell_type_group = interest_cell_type_group, interest_cell_type_tau = interest_cell_type_tau,
    interest_cell_type_genes_pseudotime_info = interest_cell_type_genes_pseudotime_info,
    max_epochs = 100, early_stop = 5, eps_theta = 1e-2, eps_loss = 1e-5, ncores = parallel::detectCores() - 1
  )
})
message(dataset, ": ", print_run_stats(monitor))
saveRDS(interest_cell_type_branch_model_train_all, file = "./4.2 BUild Prediction Model/interest_cell_type_branch_model_train_all.rds")
# A549: running: 0 hours 10 mins 27 secs, peak RAM: 0 GB 6 MB
# GM12878: running: 0 hours 12 mins 35 secs, peak RAM: 0 GB 11 MB
# BM: running: 0 hours 18 mins 41 secs, peak RAM: 0 GB 6 MB

# Remove training sets with poor training performance
interest_cell_type_branch_model_train = interest_cell_type_branch_model_train_all
for (cell_type in names(interest_cell_type_branch_model_train)) {
  df = interest_cell_type_branch_model_train[[cell_type]][["best_pred_result"]]
  bad_flag = with(df, theta_p < 0.05)
  rownames_filtered = rownames(df)[bad_flag]
  interest_cell_type_branch_model_train[[cell_type]][["best_pred_result"]] = df[!bad_flag, ]
  for (n in 1:(length(interest_cell_type_branch_model_train[[cell_type]]) - 1)) {
    lst = interest_cell_type_branch_model_train[[cell_type]][[paste0("branch", n)]]
    interest_cell_type_branch_model_train[[cell_type]][[paste0("branch", n)]] = lst[!names(lst) %in% rownames_filtered]
  }
}; rm(cell_type, df, bad_flag, rownames_filtered, n, lst)
message(
  dataset, ": ", nrow(interest_cell_type_branch_model_train_all[[1]][["best_pred_result"]]), "/",
  nrow(interest_cell_type_branch_model_train[[1]][["best_pred_result"]]), " training results across branches, ",
  paste(unlist(lapply(interest_cell_type_branch_model_train_all[[1]][2:length(interest_cell_type_branch_model_train_all[[1]])], length)), collapse = ","), "/",
  paste(unlist(lapply(interest_cell_type_branch_model_train[[1]][2:length(interest_cell_type_branch_model_train[[1]])], length)), collapse = ","), " training results per branches"
)
saveRDS(interest_cell_type_branch_model_train, file = "./4.2 BUild Prediction Model/interest_cell_type_branch_model_train.rds")
# A549: 85/82 training results across branches, 85/82 training results per branches
# GM12878: 76/73 training results across branches, 76/73 training results per branches
# BM: 92/91 training results across branches, 92/91 training results per branches



########## Part 5: Predicting GRN ##########
##### 5.1 Getting Data for Inferring GRN #####
setwd(Working_directory)
interest_cell_type_branch_iGRN = readRDS("./3 get iGRN/interest_cell_type_branch_iGRN.rds")
interest_cell_type_group = readRDS("./2.2 Data Processing - Cell Grouping/interest_cell_type_group.rds")
monitor = peakRAM::peakRAM({
  interest_cell_type_branch_data_for_GRN_infer = get_data_for_inferring_GRN(
    interest_cell_type_branch_iGRN = interest_cell_type_branch_iGRN,
    interest_cell_type_group = interest_cell_type_group,
    n_part = parallel::detectCores() - 1, ncores = parallel::detectCores() - 1
  )
})
message(dataset, ": ", print_run_stats(monitor))
saveRDS(interest_cell_type_branch_data_for_GRN_infer, file = "./5 Infer GRN/interest_cell_type_branch_data_for_GRN_infer.rds")
# A549: running: 0 hours 0 mins 37 secs, peak RAM: 0 GB 160 MB
# GM12878: running: 0 hours 0 mins 51 secs, peak RAM: 0 GB 193 MB
# BM: running: 0 hours 1 mins 30 secs, peak RAM: 0 GB 196 MB



##### 5.2 Inferring GRN for Each Branch #####
setwd(Working_directory)
interest_cell_type_branch_data_for_GRN_infer = readRDS("./5 Infer GRN/interest_cell_type_branch_data_for_GRN_infer.rds")
interest_cell_type_tau = readRDS("./3 get iGRN/interest_cell_type_tau.rds")
interest_cell_type_branch_model_train = readRDS("./4.2 BUild Prediction Model/interest_cell_type_branch_model_train.rds")
interest_cell_type_group = readRDS("./2.2 Data Processing - Cell Grouping/interest_cell_type_group.rds")
interest_cell_type_pGRN = list()
message("Infer cell type ", dataset, ".")
interest_cell_type_pGRN[[dataset]] = list()
ncores = parallel::detectCores() - 1
for (n in 1:(length(interest_cell_type_branch_data_for_GRN_infer[[dataset]]) - 1)) {
  message("Infer cell type ", dataset, " branch ", n, ".")
  each = ceiling(nrow(interest_cell_type_branch_data_for_GRN_infer[[1]][[paste0("branch", n)]]) / ncores)
  message("each: ", each)
  monitor = peakRAM::peakRAM({
    interest_cell_type_pGRN[[dataset]][[paste0("branch", n)]] = infer_GRN(
      interest_cell_type_branch_data_for_GRN_infer = interest_cell_type_branch_data_for_GRN_infer, interest_cell_type_tau = interest_cell_type_tau,
      interest_cell_type_branch_model_train = interest_cell_type_branch_model_train, interest_cell_type_group = interest_cell_type_group,
      select_cell_type = dataset, select_branch = n, each = each, ncores = ncores
    )[[dataset]][[paste0("branch", n)]]
  })
  message(dataset, " branch ", n, ": each is ", each, ", ", print_run_stats(monitor))
  # save.image(paste0("/mnt/i/scDERI_Paper/3-benchmark/8-scDERI/5_Infer_GRN_branch", n, "middle.RData")) # Check point
  # load(paste0("/mnt/i/scDERI_Paper/3-benchmark/8-scDERI/5_Infer_GRN_branch", 1, "middle.RData"))

  # A549 branch 1: each is 103, running: 0 hours 7 mins 45 secs, peak RAM: 0 GB 211 MB
  # GM12878 branch 1: each is 127, running: 0 hours 7 mins 34 secs, peak RAM: 0 GB 232 MB
  # BM branch 1: each is 142, running: 0 hours 9 mins 6 secs, peak RAM: 0 GB 261 MB
}
saveRDS(interest_cell_type_pGRN, file = "./5 Infer GRN/interest_cell_type_pGRN.rds")



##### 5.3 Selecting Best Prediction for each TF-TG #####
setwd(Working_directory)
interest_cell_type_branch_data_for_GRN_infer = readRDS("./5 Infer GRN/interest_cell_type_branch_data_for_GRN_infer.rds")
interest_cell_type_pGRN = readRDS("./5 Infer GRN/interest_cell_type_pGRN.rds")
interest_cell_type_tau = readRDS("./3 get iGRN/interest_cell_type_tau.rds")
interest_cell_type_genes_pseudotime_info = readRDS("./2.2 Data Processing - Cell Grouping/interest_cell_type_genes_pseudotime_info.rds")
ncores = parallel::detectCores() - 1
monitor = peakRAM::peakRAM({
  interest_cell_type_pGRN = select_best_prediction(
    interest_cell_type_pGRN = interest_cell_type_pGRN, interest_cell_type_tau = interest_cell_type_tau,
    interest_cell_type_genes_pseudotime_info = interest_cell_type_genes_pseudotime_info,
    fit_loss_quantile_threshold = 1, theta_p_quantile_threshold = 0,
    each = ceiling((max(sapply(
      interest_cell_type_branch_data_for_GRN_infer[[1]][1:(length(interest_cell_type_branch_data_for_GRN_infer[[1]])-1)],
      nrow
    ))) / ncores),
    ncores = ncores
  )
})
message(dataset, ": ", nrow(interest_cell_type_pGRN[[1]][["pred_summary"]]), " pairs remaining, ", print_run_stats(monitor))
saveRDS(interest_cell_type_pGRN, file = "./5 Infer GRN/interest_cell_type_pGRN.rds")
# A549: 7941 pairs remaining, running: 0 hours 0 mins 19 secs, peak RAM: 0 GB 80 MB
# GM12878: 9793 pairs remaining, running: 0 hours 0 mins 23 secs, peak RAM: 0 GB 80 MB
# BM: 10851 pairs remaining, running: 0 hours 0 mins 23 secs, peak RAM: 0 GB 111 MB



##### 5.4 Filtering GRN #####
setwd(Working_directory)
interest_cell_type_pGRN = readRDS("./5 Infer GRN/interest_cell_type_pGRN.rds")
monitor = peakRAM::peakRAM({
  interest_cell_type_pGRN = filter_GRN(
    interest_cell_type_pGRN = interest_cell_type_pGRN,
    fit_loss_quantile_threshold = 1, theta_p_quantile_threshold = 0, fit_loss_upper = 20, theta_p_lower = -0.01
  )
})
message(dataset, ": ", nrow(interest_cell_type_pGRN[[1]][["pred_summary_loose"]]), " loose pairs, ", print_run_stats(monitor))
saveRDS(interest_cell_type_pGRN, file = "./5 Infer GRN/interest_cell_type_pGRN.rds")
# A549: 7939 loose pairs, running: 0 hours 0 mins 0 secs, peak RAM: 0 GB 3 MB
# GM12878: 9791 loose pairs, running: 0 hours 0 mins 0 secs, peak RAM: 0 GB 4 MB
# BM: 10849 loose pairs, running: 0 hours 0 mins 0 secs, peak RAM: 0 GB 4 MB



##### 5.5 Inferring peGRN #####
setwd(Working_directory)
scRNAseq = qs::qread(file = "./scRNAseq.qsave")
scATACseq = qs::qread("./scATACseq.qsave")
results_identify_TGs = readRDS("./1.2 Data Preprocessing - Identify TGs By Annotation/results_identify_TGs.rds")
interest_cell_type_iGRN = readRDS("./3 get iGRN/interest_cell_type_iGRN.rds")
interest_cell_type_pGRN = readRDS("./5 Infer GRN/interest_cell_type_pGRN.rds")
monitor = peakRAM::peakRAM({
  interest_cell_type_peGRN = infer_eGRN(
    scRNAseq = scRNAseq, scATACseq = scATACseq, results_identify_TGs = results_identify_TGs,
    interest_cell_type_iGRN = interest_cell_type_iGRN, interest_cell_type_pGRN = interest_cell_type_pGRN,
  )
})
message(dataset, ": ", print_run_stats(monitor))
saveRDS(interest_cell_type_peGRN, file = "./5 Infer GRN/interest_cell_type_peGRN.rds")
# A549: running: 0 hours 0 mins 30 secs, peak RAM: 0 GB 388 MB
# GM12878: running: 0 hours 0 mins 43 secs, peak RAM: 0 GB 862 MB
# BM: running: 0 hours 0 mins 49 secs, peak RAM: 1 GB 378 MB



##### 5.6 Filtering peGRN #####
setwd(Working_directory)
interest_cell_type_peGRN = readRDS("./5 Infer GRN/interest_cell_type_peGRN.rds")
monitor = peakRAM::peakRAM({
  interest_cell_type_peGRN = filter_eGRN(
    interest_cell_type_peGRN = interest_cell_type_peGRN,
    theta_p_points = c(0, 0.5, 1), initial_quantile_points = c(5.79, 4.66, 4.23),
    rescue_peakTG = TRUE, peakTG_rescue_quantile_threshold = 0.25, min_TF_per_peakTG = 1,
    rescue_TF = TRUE, TF_rescue_quantile_threshold = 0.0, min_peakTG_per_TF = 7,
    factor_range = list(`FG-NL` = 0.43, `FG-theta_i-R` = 0.64,
                        `F-NZER-R` = 0.28, `F-NZEA-R` = 0.33, `F-EA-R` = 0.22,
                        `P-NZAR-R` = 0.19, `P-NZAA-R` = 0.67, `P-AA-R` = 0.13,
                        `G-NZER-R` = 0.21, `G-NZEA-R` = 0.12, `G-EA-R` = 0.26)
  )
})
message(dataset, ": ",
        length(unique(interest_cell_type_peGRN[[1]][["loose"]][["TF"]])), "/", length(unique(interest_cell_type_peGRN[[1]][["loose_filtered_simplified"]][["TF"]])), " TFs, ",
        length(unique(interest_cell_type_peGRN[[1]][["loose"]][["peak"]])), "/", length(unique(interest_cell_type_peGRN[[1]][["loose_filtered_simplified"]][["peak"]])), " peaks, ",
        length(unique(interest_cell_type_peGRN[[1]][["loose"]][["TG"]])), "/", length(unique(interest_cell_type_peGRN[[1]][["loose_filtered_simplified"]][["TG"]])), " TGs, ",
        nrow(unique(interest_cell_type_peGRN[[1]][["loose"]][, c("TF", "peak")])), "/", nrow(unique(interest_cell_type_peGRN[[1]][["loose_filtered_simplified"]][, c("TF", "peak")])), " TF-peak, ",
        nrow(unique(interest_cell_type_peGRN[[1]][["loose"]][, c("peak", "TG")])), "/", nrow(unique(interest_cell_type_peGRN[[1]][["loose_filtered_simplified"]][, c("peak", "TG")])), " peak-TG, ",
        nrow(unique(interest_cell_type_peGRN[[1]][["loose"]][, c("TF", "TG")])), "/", nrow(unique(interest_cell_type_peGRN[[1]][["loose_filtered_simplified"]][, c("TF", "TG")])), " TF-TG, ",
        nrow(unique(interest_cell_type_peGRN[[1]][["loose"]][, c("TF", "peak", "TG")])), "/", nrow(unique(interest_cell_type_peGRN[[1]][["loose_filtered_simplified"]][, c("TF", "peak", "TG")])), " TF-peak-TG, rescue ratio ",
        round(1 - sum(interest_cell_type_peGRN[[1]][["loose_filtered"]][["filter_keep"]]) / nrow(interest_cell_type_peGRN[[1]][["loose_filtered_simplified"]]), 3), ", ", print_run_stats(monitor))
saveRDS(interest_cell_type_peGRN, file = "./5 Infer GRN/interest_cell_type_peGRN.rds")
qs::qsave(interest_cell_type_peGRN, "./interest_cell_type_peGRN.qsave")
# A549: 84/84 TFs, 1341/1127 peaks, 1319/1114 TGs, 9440/2081 TF-peak, 1640/1324 peak-TG, 7939/1935 TF-TG, 9778/2128 TF-peak-TG, rescue ratio 0.316, running: 0 hours 0 mins 0 secs, peak RAM: 0 GB 44 MB
# GM12878: 76/76 TFs, 1642/1395 peaks, 1795/1498 TGs, 11316/2672 TF-peak, 2188/1769 peak-TG, 9791/2483 TF-TG, 11891/2763 TF-peak-TG, rescue ratio 0.3, running: 0 hours 0 mins 0 secs, peak RAM: 0 GB 55 MB
# BM: 92/92 TFs, 2083/1787 peaks, 1976/1725 TGs, 13071/3219 TF-peak, 2518/2043 peak-TG, 10849/3011 TF-TG, 13655/3309 TF-peak-TG, rescue ratio 0.257, running: 0 hours 0 mins 0 secs, peak RAM: 0 GB 60 MB

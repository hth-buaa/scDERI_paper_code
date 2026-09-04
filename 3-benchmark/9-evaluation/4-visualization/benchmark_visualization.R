################################################################
#                                                              #
#            Benchmarking Visualization Script                 #
#                                                              #
################################################################



##### 1. Environment Setup and Package Loading #####
# General libraries
library(Pando)
library(dplyr)
library(qs)
library(ggplot2)
library(Cairo)
library(data.table)
library(plyr)
library(ggpubr)

# Set parameters
parent.dir = "/mnt/i/scDERI_Paper/3-benchmark/9-evaluation/4-visualization/"
base.dir = "/mnt/i/scDERI_Paper/3-benchmark"
# setwd(parent.dir)



##### 2. Load Benchmark Results #####
# Panels to generate:
# 1. Chart of benchmarking pipeline
# 2. Numbers of TFs, regions per regulon, and genes per regulon for 8 methods
# 3. Bar plots of f-scores for TF coverage
# 4. Box plots of precision for region-to-gene linkages
# 5. Bar plots of f-scores for region-to-gene linkages

# Load TF coverage benchmark results
TF_coverage_benchmark = qs::qread(paste0(base.dir, "/9-evaluation/1-TF_coverage/TF_coverage_benchmark.qsave"))
TF_coverage_benchmark = TF_coverage_benchmark[sapply(TF_coverage_benchmark, sum) > 0]
print(TF_coverage_benchmark)
# $`10X_hg38_Human_bone_marrow_with_T_cell_depleted`
#                   prec recall    f.score
# scDERI     0.032608696    0.3 0.05882353
# STREAM     0.020000000    0.3 0.03750000
# scGLUE     0.009505703    0.5 0.01865672
# Pando      0.006802721    0.1 0.01273885
# scMEGA     0.034482759    0.2 0.05882353
# FigR       0.008928571    0.4 0.01746725
# DIRECT_NET 0.022988506    0.2 0.04123711
# IReNA      0.008620690    0.3 0.01675978
#
# $SNAREseq2_hg38_GM12878
#                  prec     recall    f.score
# scDERI     0.34210526 0.15950920 0.21757322
# STREAM     0.31578947 0.18404908 0.23255814
# scGLUE     0.17630854 0.39263804 0.24334601
# Pando      0.00000000 0.00000000 0.00000000
# scMEGA     0.16000000 0.02453988 0.04255319
# FigR       0.07471264 0.07975460 0.07715134
# DIRECT_NET 0.00000000 0.00000000 0.00000000
# IReNA      0.03448276 0.02453988 0.02867384
#
# $SNAREseq2_hg38_A549
#                  prec     recall    f.score
# scDERI     0.21428571 0.21176471 0.21301775
# STREAM     0.16788321 0.27058824 0.20720721
# scGLUE     0.09051724 0.49411765 0.15300546
# Pando      0.00000000 0.00000000 0.00000000
# scMEGA     0.04000000 0.01176471 0.01818182
# FigR       0.03319502 0.09411765 0.04907975
# DIRECT_NET 0.00000000 0.00000000 0.00000000
# IReNA      0.01020408 0.02352941 0.01423488

# Load TF-region benchmark results
TF_region_benchmark = qs::qread(paste0(base.dir, "/9-evaluation/2-TF-peak_relationship/TF_region_benchmark.qsave"))
TF_region_benchmark = TF_region_benchmark[sapply(TF_region_benchmark, sum) > 0]
print(TF_region_benchmark)
# $`10X_hg38_Human_bone_marrow_with_T_cell_depleted`
#         precision       recall      f.score
# scDERI 0.49206349 0.0003379604 0.0006754542
# STREAM 0.39259690 0.0037021453 0.0073279626
# scGLUE 0.09633912 0.0008693266 0.0017177718
# Pando  0.33333333 0.0001001904 0.0002003205
#
# $SNAREseq2_hg38_GM12878
#         precision      recall     f.score
# scDERI 0.36306967 0.001018087 0.002025024
# STREAM 0.34646050 0.002085584 0.004135767
# scGLUE 0.07400771 0.002496815 0.004784046
# Pando  0.00000000 0.000000000 0.000000000
#
# $SNAREseq2_hg38_A549
#        precision       recall     f.score
# scDERI 0.3418885 0.0006484421 0.001289745
# STREAM 0.3150941 0.0015431271 0.003061434
# scGLUE 0.0740586 0.0014694668 0.002850218
# Pando  0.0000000 0.0000000000 0.000000000

TF_region_benchmark_all = qs::qread(paste0(base.dir, "/9-evaluation/2-TF-peak_relationship/TF_region_benchmark_all.qsave"))

# Load region-gene benchmark results
region_TG_benchmark = qs::qread(paste0(base.dir, "/9-evaluation/3-peak-TG_relationship/region_TG_benchmark.qsave"))
region_TG_benchmark = region_TG_benchmark[sapply(region_TG_benchmark, sum) > 0]
print(region_TG_benchmark)
# $`10X_hg38_Human_bone_marrow_with_T_cell_depleted`
#        precision     recall    f.score
# scDERI 0.9955947 0.35207101 0.52018839
# STREAM 0.9765835 0.28951817 0.44662870
# scGLUE 1.0000000 0.20245139 0.33673111
# Pando  1.0000000 0.05029586 0.09577465
#
# $SNAREseq2_hg38_GM12878
#        precision    recall   f.score
# scDERI 0.9954777 0.3339056 0.5000748
# STREAM 0.9713428 0.3351931 0.4983980
# scGLUE 0.9990010 0.2171674 0.3567770
# Pando  0.0000000 0.0000000 0.0000000
#
# $SNAREseq2_hg38_A549
#        precision    recall   f.score
# scDERI 0.9833837 0.2562368 0.4065423
# STREAM 0.9664607 0.2131078 0.3492130
# scGLUE 0.9939024 0.1987315 0.3312328
# Pando  0.0000000 0.0000000 0.0000000

# Identify datasets present in all three benchmark results
dataset = Reduce("intersect", list(names(TF_coverage_benchmark), names(TF_region_benchmark), names(region_TG_benchmark)))
print(dataset)
# [1] "10X_hg38_Human_bone_marrow_with_T_cell_depleted" "SNAREseq2_hg38_GM12878"                          "SNAREseq2_hg38_A549"
dataset_simplify = c("BM", "GM12878", "A549")

# Define factor levels for methods
method_factor_levels = c("scDERI", "STREAM", "GLUE", "Pando", "scMEGA", "FigR", "DIRECT-NET", "IReNA")

# Set colors
colors = setNames(c("#C70039", "#614BC3", "#EA1179", "#1D5D9B", "#1A5D1A", "#FD8D14", "#068FFF", "darkgrey"), method_factor_levels)



##### 3. Load Regulatory Networks (Regulons) from Each Method #####
### 3.1 scDERI regulons
message("Loading scDERI regulons...")
scDERI_regulons = pbapply::pblapply(seq_along(dataset), function(i) {
  message("----> Conversion for the ", i, "-th dataset: ", dataset[i])
  if (!file.exists(paste0(base.dir, "/8-scDERI/", dataset[i], "/interest_cell_type_peGRN.qsave"))) {
    message("----> No eRegulon file on dataset: ", dataset[i])
    return(NULL)
  }
  qs::qread(paste0(base.dir, "/8-scDERI/", dataset[i], "/interest_cell_type_peGRN.qsave"))[[1]][["loose_filtered_simplified"]]
})
names(scDERI_regulons) = dataset
qs::qsave(scDERI_regulons, paste0(parent.dir, "scDERI_regulons.qsave"))
# scDERI_regulons = qs::qread(paste0(parent.dir, "scDERI_regulons.qsave"))
message("Number of datasets loaded: ", length(scDERI_regulons))
# Number of datasets loaded: 3


### 3.2 STREAM regulons
message("Loading STREAM regulons...")
STREAM_regulons = pbapply::pblapply(seq_along(dataset), function(i) {
  message("----> Conversion for the ", i, "-th dataset: ", dataset[i])
  if (!file.exists(paste0(base.dir, "/1-STREAM/", dataset[i], "/eRegulons.qsave"))) {
    message("----> No eRegulon file on dataset: ", dataset[i])
    return(NULL)
  }
  qs::qread(paste0(base.dir, "/1-STREAM/", dataset[i], "/eRegulons.qsave"))
})
names(STREAM_regulons) = dataset
qs::qsave(STREAM_regulons, paste0(parent.dir, "STREAM_regulons.qsave"))
# STREAM_regulons = qs::qread(paste0(parent.dir, "STREAM_regulons.qsave"))
message("Number of datasets loaded: ", length(STREAM_regulons))
# Number of datasets loaded: 3


### 3.3 IReNA regulons
message("Loading IReNA regulons...")
IReNA_regulons = pbapply::pblapply(seq_along(dataset), function(i) {
  message("----> Conversion for the ", i, "-th dataset: ", dataset[i])
  if (!file.exists(paste0(base.dir, "/2-IReNA/", dataset[i], "/IReNA_result.qsave"))) {
    message("----> No eRegulon file on dataset: ", dataset[i])
    return(NULL)
  }
  qs::qread(paste0(base.dir, "/2-IReNA/", dataset[i], "/IReNA_result.qsave"))
})
names(IReNA_regulons) = dataset
qs::qsave(IReNA_regulons, paste0(parent.dir, "IReNA_regulons.qsave"))
# IReNA_regulons = qs::qread(paste0(parent.dir, "IReNA_regulons.qsave"))
message("Number of datasets loaded: ", length(IReNA_regulons))
# Number of datasets loaded: 3


### 3.4 FigR regulons
message("Loading FigR regulons...")
FigR_regulons = pbapply::pblapply(seq_along(dataset), function(i) {
  message("----> Conversion for the ", i, "-th dataset: ", dataset[i])
  if (!file.exists(paste0(base.dir, "/3-FigR/", dataset[i], "/FigR_result.qsave"))) {
    message("----> No eRegulon file on dataset: ", dataset[i])
    return(NULL)
  }
  qs::qread(paste0(base.dir, "/3-FigR/", dataset[i], "/FigR_result.qsave"))
})
names(FigR_regulons) = dataset
qs::qsave(FigR_regulons, paste0(parent.dir, "FigR_regulons.qsave"))
# FigR_regulons = qs::qread(paste0(parent.dir, "FigR_regulons.qsave"))
message("Number of datasets loaded: ", length(FigR_regulons))
# Number of datasets loaded: 3


### 3.5 GLUE regulons
message("Loading GLUE regulons...")
GLUE_regulons = pbapply::pblapply(seq_along(dataset), function(i) {
  message("----> Conversion for the ", i, "-th dataset: ", dataset[i])
  if (!file.exists(paste0(base.dir, "/4-scGLUE/", dataset[i], "/draft_grn.csv"))) {
    message("----> No eRegulon file on dataset: ", dataset[i])
    return(NULL)
  }
  grn = read.csv(paste0(base.dir, "/4-scGLUE/", dataset[i], "/draft_grn.csv"), sep = ",")
  colnames(grn) = c("TF", "Gene", "Strength")
  gene2peak = read.csv(paste0(base.dir, "/4-scGLUE/", dataset[i], "/gene2peak.csv"), sep = ",")
  colnames(gene2peak) = c("Gene", "Region", "Weight")
  common.genes = intersect(unique(grn$Gene), unique(gene2peak$Gene))
  merge(grn[grn$Gene %in% common.genes,, drop = FALSE],
        gene2peak[gene2peak$Gene %in% common.genes,, drop = FALSE], by = "Gene")
})
names(GLUE_regulons) = dataset
qs::qsave(GLUE_regulons, paste0(parent.dir, "GLUE_regulons.qsave"))
# GLUE_regulons = qs::qread(paste0(parent.dir, "GLUE_regulons.qsave"))
message("Number of datasets loaded: ", length(GLUE_regulons))
# Number of datasets loaded: 3


### 3.6 DIRECT-NET regulons
message("Loading DIRECT-NET regulons...")
DIRECT_NET_regulons = pbapply::pblapply(seq_along(dataset), function(i) {
  message("----> Conversion for the ", i, "-th dataset: ", dataset[i])
  if (!file.exists(paste0(base.dir, "/5-DIRECT-NET/", dataset[i], "/", dataset[i], "_directnet_out.qs"))) {
    message("----> No eRegulon file on dataset: ", dataset[i])
    return(NULL)
  }
  qs::qread(paste0(base.dir, "/5-DIRECT-NET/", dataset[i], "/", dataset[i], "_directnet_out.qs"))$link
})
names(DIRECT_NET_regulons) = dataset
qs::qsave(DIRECT_NET_regulons, paste0(parent.dir, "DIRECT_NET_regulons.qsave"))
# DIRECT_NET_regulons = qs::qread(paste0(parent.dir, "DIRECT_NET_regulons.qsave"))
message("Number of datasets loaded: ", length(DIRECT_NET_regulons))
# Number of datasets loaded: 3


### 3.7 Pando regulons
message("Loading Pando regulons...")
Pando_regulons = pbapply::pblapply(seq_along(dataset), function(i) {
  message("----> Conversion for the ", i, "-th dataset: ", dataset[i])
  if (!file.exists(paste0(base.dir, "/6-Pando/", dataset[i], ".qs"))) {
    message("----> No eRegulon file on dataset: ", dataset[i])
    return(NULL)
  }
  dt = qs::qread(paste0(base.dir, "/6-Pando/", dataset[i], ".qs")) %>% NetworkModules
  dt@meta
})
names(Pando_regulons) = dataset
qs::qsave(Pando_regulons, paste0(parent.dir, "Pando_regulons.qsave"))
# Pando_regulons = qs::qread(paste0(parent.dir, "Pando_regulons.qsave"))
message("Number of datasets loaded: ", length(Pando_regulons))
# Number of datasets loaded: 3


### 3.8 scMEGA regulons
message("Loading scMEGA regulons...")
scMEGA_regulons = pbapply::pblapply(seq_along(dataset), function(i) {
  message("----> Conversion for the ", i, "-th dataset: ", dataset[i])
  if (!file.exists(paste0(base.dir, "/7-scMEGA/", dataset[i], ".RData"))) {
    message("----> No eRegulon file on dataset: ", dataset[i])
    return(NULL)
  }
  load(paste0(base.dir, "/7-scMEGA/", dataset[i], ".RData"))
  tf.gene.cor
})
names(scMEGA_regulons) = dataset
qs::qsave(scMEGA_regulons, paste0(parent.dir, "scMEGA_regulons.qsave"))
# scMEGA_regulons = qs::qread(paste0(parent.dir, "scMEGA_regulons.qsave"))
message("Number of datasets loaded: ", length(scMEGA_regulons))
# Number of datasets loaded: 3



##### 4. Analyze Number of TFs per Method #####
### 4.1 Count TFs for each method
# scDERI
nTFs_scDERI = sapply(dataset, function(x) {scDERI_regulons[[x]][["TF"]] %>% unique %>% length})
df_nTFs_scDERI = data.frame(nTFs = nTFs_scDERI, Method = rep("scDERI", length(nTFs_scDERI))); rm(nTFs_scDERI)
print(df_nTFs_scDERI)
#                                                 nTFs Method
# 10X_hg38_Human_bone_marrow_with_T_cell_depleted   92 scDERI
# SNAREseq2_hg38_GM12878                            76 scDERI
# SNAREseq2_hg38_A549                               84 scDERI

# STREAM
nTFs_STREAM = sapply(dataset, function(x) {sapply(STREAM_regulons[[x]], "[[", "TF") %>% unique %>% length})
df_nTFs_STREAM = data.frame(nTFs = nTFs_STREAM, Method = rep("STREAM", length(nTFs_STREAM))); rm(nTFs_STREAM)
print(df_nTFs_STREAM)
#                                                 nTFs Method
# 10X_hg38_Human_bone_marrow_with_T_cell_depleted  150 STREAM
# SNAREseq2_hg38_GM12878                            95 STREAM
# SNAREseq2_hg38_A549                              137 STREAM

# IReNA
nTFs_IReNA = sapply(dataset, function(x) {IReNA_regulons[[x]][["FOSF_RegMTF_Cor_EnTFsTarg"]][["TF"]] %>% as.character %>% unique %>% length})
df_nTFs_IReNA = data.frame(nTFs = nTFs_IReNA, Method = rep("IReNA", length(nTFs_IReNA))); rm()
print(df_nTFs_IReNA)
#                                                 nTFs Method
# 10X_hg38_Human_bone_marrow_with_T_cell_depleted  348  IReNA
# SNAREseq2_hg38_GM12878                           116  IReNA
# SNAREseq2_hg38_A549                              196  IReNA

# FigR
nTFs_FigR = sapply(dataset, function(x) {FigR_regulons[[x]]$Motif %>% unique %>% length})
df_nTFs_FigR = data.frame(nTFs = nTFs_FigR, Method = rep("FigR", length(nTFs_FigR))); rm(nTFs_FigR)
print(df_nTFs_FigR)
#                                                 nTFs Method
# 10X_hg38_Human_bone_marrow_with_T_cell_depleted  448   FigR
# SNAREseq2_hg38_GM12878                           174   FigR
# SNAREseq2_hg38_A549                              241   FigR

# GLUE
nTFs_GLUE = sapply(dataset, function(x) {unique(GLUE_regulons[[x]][, "TF"]) %>% length})
df_nTFs_GLUE = data.frame(nTFs = nTFs_GLUE, Method = rep("GLUE", length(nTFs_GLUE))); rm(nTFs_GLUE)
print(df_nTFs_GLUE)
#                                                 nTFs Method
# 10X_hg38_Human_bone_marrow_with_T_cell_depleted  526   GLUE
# SNAREseq2_hg38_GM12878                           363   GLUE
# SNAREseq2_hg38_A549                              464   GLUE

# DIRECT-NET
nTFs_DIRECT_NET = sapply(dataset, function(x) {unique(unlist(DIRECT_NET_regulons[[x]][, "TF"])) %>% length})
df_nTFs_DIRECT_NET = data.frame(nTFs = nTFs_DIRECT_NET, Method = rep("DIRECT-NET", length(nTFs_DIRECT_NET))); rm(nTFs_DIRECT_NET)
print(df_nTFs_DIRECT_NET)
#                                                 nTFs     Method
# 10X_hg38_Human_bone_marrow_with_T_cell_depleted   87 DIRECT-NET
# SNAREseq2_hg38_GM12878                             0 DIRECT-NET
# SNAREseq2_hg38_A549                                0 DIRECT-NET

# Pando
nTFs_Pando = sapply(dataset, function(x) {unique(unlist(Pando_regulons[[x]][, "tf"])) %>% length})
df_nTFs_Pando = data.frame(nTFs = nTFs_Pando, Method = rep("Pando", length(nTFs_Pando))); rm(nTFs_Pando)
print(df_nTFs_Pando)
#                                                 nTFs Method
# 10X_hg38_Human_bone_marrow_with_T_cell_depleted  147  Pando
# SNAREseq2_hg38_GM12878                             0  Pando
# SNAREseq2_hg38_A549                                0  Pando

# scMEGA
nTFs_scMEGA = sapply(dataset, function(x) {unique(unlist(scMEGA_regulons[[x]][, "tf"])) %>% length})
df_nTFs_scMEGA = data.frame(nTFs = nTFs_scMEGA, Method = rep("scMEGA", length(nTFs_scMEGA))); rm(nTFs_scMEGA)
print(df_nTFs_scMEGA)
#                                                 nTFs Method
# 10X_hg38_Human_bone_marrow_with_T_cell_depleted   58 scMEGA
# SNAREseq2_hg38_GM12878                            25 scMEGA
# SNAREseq2_hg38_A549                               25 scMEGA


### 4.2 Create box plot for number of TFs
nTFs_df = rbind(
  df_nTFs_scDERI,
  df_nTFs_STREAM,
  df_nTFs_GLUE,
  df_nTFs_Pando,
  df_nTFs_scMEGA,
  df_nTFs_FigR,
  df_nTFs_DIRECT_NET,
  df_nTFs_IReNA
)
nTFs_df$Method = factor(nTFs_df$Method, levels = method_factor_levels)
print(nTFs_df)
#                                                  nTFs     Method
# 10X_hg38_Human_bone_marrow_with_T_cell_depleted    92     scDERI
# SNAREseq2_hg38_GM12878                             76     scDERI
# SNAREseq2_hg38_A549                                84     scDERI
# 10X_hg38_Human_bone_marrow_with_T_cell_depleted1  150     STREAM
# SNAREseq2_hg38_GM128781                            95     STREAM
# SNAREseq2_hg38_A5491                              137     STREAM
# 10X_hg38_Human_bone_marrow_with_T_cell_depleted2  526       GLUE
# SNAREseq2_hg38_GM128782                           363       GLUE
# SNAREseq2_hg38_A5492                              464       GLUE
# 10X_hg38_Human_bone_marrow_with_T_cell_depleted3  147      Pando
# SNAREseq2_hg38_GM128783                             0      Pando
# SNAREseq2_hg38_A5493                                0      Pando
# 10X_hg38_Human_bone_marrow_with_T_cell_depleted4   58     scMEGA
# SNAREseq2_hg38_GM128784                            25     scMEGA
# SNAREseq2_hg38_A5494                               25     scMEGA
# 10X_hg38_Human_bone_marrow_with_T_cell_depleted5  448       FigR
# SNAREseq2_hg38_GM128785                           174       FigR
# SNAREseq2_hg38_A5495                              241       FigR
# 10X_hg38_Human_bone_marrow_with_T_cell_depleted6   87 DIRECT-NET
# SNAREseq2_hg38_GM128786                             0 DIRECT-NET
# SNAREseq2_hg38_A5496                                0 DIRECT-NET
# 10X_hg38_Human_bone_marrow_with_T_cell_depleted7  348      IReNA
# SNAREseq2_hg38_GM128787                           116      IReNA
# SNAREseq2_hg38_A5497                              196      IReNA

boxplots_nTFs = ggplot(nTFs_df, aes(x = Method, y = nTFs, fill = Method)) +
  geom_boxplot(lwd = 0.1) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.75, linetype = "dashed", color = "black", size = 0.2) +
  scale_fill_manual(values = colors) +
  coord_cartesian(ylim = c(NA, 510)) +
  labs(title = "TF Numbers") +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = NA, color = "black", size = 0.5) +
  theme(
    axis.title = element_text(size = 12, color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 9, angle = 45, hjust = 1),
    plot.title = element_text(size = 10, color = "black", hjust = 0.5, face = "bold"),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    # axis.text.x = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
    panel.grid.minor = element_line(color = "grey95", linewidth = 0.3),
    panel.background = element_rect(fill = "white", linewidth = 0.3),
    legend.position = "none"
  )
print(boxplots_nTFs)
qs::qsave(boxplots_nTFs, paste0(parent.dir, "boxplots_nTFs.qsave"))

CairoPNG(filename = paste0(parent.dir, "boxplots_nTFs.png"), width = 2.4, height = 3.2, units = "in", dpi = 300)
print(boxplots_nTFs)
dev.off()



##### 5. Analyze Number of Genes per Regulon #####
# 5.1 Count genes per regulon for each method
nTGsPerReg_scDERI = sapply(dataset, function(x) {
  if (is.null(scDERI_regulons[[x]])) {return(NULL)}
  GRN = unique(scDERI_regulons[[x]][, c("TF", "TG")])
  tapply(GRN$TG, GRN$TF, function(x) length(unique(x)))
}) %>% unlist
if (is.null(nTGsPerReg_scDERI)) {nTGsPerReg_scDERI = 0}
df_nTGsPerReg_scDERI = data.frame(
  nTG_per_regulon = nTGsPerReg_scDERI,
  Method = rep("scDERI", length(nTGsPerReg_scDERI))
); rm(nTGsPerReg_scDERI)

df_nTGsPerReg_scDERI$dataset = sub("(SNAREseq2_hg38_GM12878|SNAREseq2_hg38_A549|10X_hg38_Human_bone_marrow_with_T_cell_depleted).*", "\\1", rownames(df_nTGsPerReg_scDERI))
print(aggregate(nTG_per_regulon ~ dataset, data = df_nTGsPerReg_scDERI, FUN = mean))
#                                           dataset nTG_per_regulon
# 1 10X_hg38_Human_bone_marrow_with_T_cell_depleted        32.72826
# 2                             SNAREseq2_hg38_A549        23.03571
# 3                          SNAREseq2_hg38_GM12878        32.67105
df_nTGsPerReg_scDERI$dataset = NULL

rownames(df_nTGsPerReg_scDERI) = NULL
print(head(df_nTGsPerReg_scDERI))
#   nTG_per_regulon Method
# 1              40 scDERI
# 2              52 scDERI
# 3              32 scDERI
# 4              13 scDERI
# 5               7 scDERI
# 6               7 scDERI
print(dim(df_nTGsPerReg_scDERI))
# [1] 252   2

# STREAM
nTGsPerReg_STREAM = sapply(dataset, function(x) {
  if (is.null(STREAM_regulons[[x]])) {return(NULL)}
  sapply(STREAM_regulons[[x]], "[[", "genes") %>% sapply(., length)
}) %>% unlist
if (is.null(nTGsPerReg_STREAM)) {nTGsPerReg_STREAM = 0}
df_nTGsPerReg_STREAM = data.frame(
  nTG_per_regulon = nTGsPerReg_STREAM,
  Method = rep("STREAM", length(nTGsPerReg_STREAM))
); rm(nTGsPerReg_STREAM)

df_nTGsPerReg_STREAM$dataset = sub("(SNAREseq2_hg38_GM12878|SNAREseq2_hg38_A549|10X_hg38_Human_bone_marrow_with_T_cell_depleted).*", "\\1", rownames(df_nTGsPerReg_STREAM))
print(aggregate(nTG_per_regulon ~ dataset, data = df_nTGsPerReg_STREAM, FUN = mean))
#                                           dataset nTG_per_regulon
# 1 10X_hg38_Human_bone_marrow_with_T_cell_depleted       122.02465
# 2                             SNAREseq2_hg38_A549        57.07299
# 3                          SNAREseq2_hg38_GM12878        81.00000
df_nTGsPerReg_STREAM$dataset = NULL

rownames(df_nTGsPerReg_STREAM) = NULL
print(head(df_nTGsPerReg_STREAM))
#   nTG_per_regulon Method
# 1             526 STREAM
# 2             470 STREAM
# 3             453 STREAM
# 4             470 STREAM
# 5             440 STREAM
# 6             380 STREAM
print(dim(df_nTGsPerReg_STREAM))
# [1] 575   2

# IReNA
nTGsPerReg_IReNA = sapply(dataset, function(x) {
  if (is.null(IReNA_regulons[[x]])) {return(NULL)}
  IReNA_regulons[[x]][["FOSF_RegMTF_Cor_EnTFsTarg"]][, c("TF", "Target")] %>%
    dplyr::distinct(TF, Target, .keep_all = FALSE) %>%
    dplyr::count(TF, name = "Target") %>%
    tibble::column_to_rownames("TF")
}) %>% unlist
if (is.null(nTGsPerReg_IReNA)) {nTGsPerReg_IReNA = 0}
df_nTGsPerReg_IReNA = data.frame(
  nTG_per_regulon = nTGsPerReg_IReNA,
  Method = rep("IReNA", length(nTGsPerReg_IReNA))
); rm(nTGsPerReg_IReNA)

df_nTGsPerReg_IReNA$dataset = sub("(SNAREseq2_hg38_GM12878|SNAREseq2_hg38_A549|10X_hg38_Human_bone_marrow_with_T_cell_depleted).*", "\\1", rownames(df_nTGsPerReg_IReNA))
print(aggregate(nTG_per_regulon ~ dataset, data = df_nTGsPerReg_IReNA, FUN = mean))
#                                           dataset nTG_per_regulon
# 1 10X_hg38_Human_bone_marrow_with_T_cell_depleted       11.870690
# 2                             SNAREseq2_hg38_A549        3.586735
# 3                          SNAREseq2_hg38_GM12878        3.017241
df_nTGsPerReg_IReNA$dataset = NULL

rownames(df_nTGsPerReg_IReNA) = NULL
print(head(df_nTGsPerReg_IReNA))
#   nTG_per_regulon Method
# 1               7  IReNA
# 2               3  IReNA
# 3              45  IReNA
# 4              15  IReNA
# 5               4  IReNA
# 6              22  IReNA
print(dim(df_nTGsPerReg_IReNA))
# [1] 660   2

# FigR
nTGsPerReg_FigR = sapply(dataset, function(x) {
  if (is.null(FigR_regulons[[x]])) {return(NULL)}
  FigR_regulons[[x]][, c("Motif", "DORC")] %>%
    dplyr::distinct(Motif, DORC, .keep_all = FALSE) %>%
    dplyr::count(Motif, name = "DORC") %>%
    tibble::column_to_rownames("Motif")
}) %>% unlist
if (is.null(nTGsPerReg_FigR)) {nTGsPerReg_FigR = 0}
df_nTGsPerReg_FigR = data.frame(
  nTG_per_regulon = nTGsPerReg_FigR,
  Method = rep("FigR", length(nTGsPerReg_FigR))
); rm(nTGsPerReg_FigR)

df_nTGsPerReg_FigR$dataset = sub("(SNAREseq2_hg38_GM12878|SNAREseq2_hg38_A549|10X_hg38_Human_bone_marrow_with_T_cell_depleted).*", "\\1", rownames(df_nTGsPerReg_FigR))
print(aggregate(nTG_per_regulon ~ dataset, data = df_nTGsPerReg_FigR, FUN = mean))
#                                           dataset nTG_per_regulon
# 1 10X_hg38_Human_bone_marrow_with_T_cell_depleted        12.12723
# 2                             SNAREseq2_hg38_A549         9.39834
# 3                          SNAREseq2_hg38_GM12878        10.01149
df_nTGsPerReg_FigR$dataset = NULL

rownames(df_nTGsPerReg_FigR) = NULL
print(head(df_nTGsPerReg_FigR))
#   nTG_per_regulon Method
# 1               8   FigR
# 2              10   FigR
# 3               2   FigR
# 4               4   FigR
# 5              22   FigR
# 6               4   FigR
print(dim(df_nTGsPerReg_FigR))
# [1] 863   2

# GLUE
nTGsPerReg_GLUE = sapply(dataset, function(x) {
  if (is.null(GLUE_regulons[[x]])) {return(NULL)}
  split(GLUE_regulons[[x]], f = GLUE_regulons[[x]]$TF) %>% sapply(., "[[", "Gene") %>% sapply(., unique) %>% sapply(., length)
}) %>% unlist
if (is.null(nTGsPerReg_GLUE)) {nTGsPerReg_GLUE = 0}
df_nTGsPerReg_GLUE = data.frame(
  nTG_per_regulon = nTGsPerReg_GLUE,
  Method = rep("GLUE", length(nTGsPerReg_GLUE))
); rm(nTGsPerReg_GLUE)

df_nTGsPerReg_GLUE$dataset = sub("(SNAREseq2_hg38_GM12878|SNAREseq2_hg38_A549|10X_hg38_Human_bone_marrow_with_T_cell_depleted).*", "\\1", rownames(df_nTGsPerReg_GLUE))
print(aggregate(nTG_per_regulon ~ dataset, data = df_nTGsPerReg_GLUE, FUN = mean))
#                                           dataset nTG_per_regulon
# 1 10X_hg38_Human_bone_marrow_with_T_cell_depleted        99.21483
# 2                             SNAREseq2_hg38_A549        95.73060
# 3                          SNAREseq2_hg38_GM12878       131.51515
df_nTGsPerReg_GLUE$dataset = NULL

rownames(df_nTGsPerReg_GLUE) = NULL
print(head(df_nTGsPerReg_GLUE))
#   nTG_per_regulon Method
# 1             194   GLUE
# 2             106   GLUE
# 3             245   GLUE
# 4             258   GLUE
# 5              69   GLUE
# 6             264   GLUE
print(dim(df_nTGsPerReg_GLUE))
# [1] 1353    2

# DIRECT-NET
nTGsPerReg_DIRECT_NET = sapply(dataset, function(x) {
  if (is.null(DIRECT_NET_regulons[[x]])) {return(NULL)}
  split(DIRECT_NET_regulons[[x]], f = DIRECT_NET_regulons[[x]]$TF) %>% sapply(., "[[", "Gene") %>% sapply(., unique) %>% sapply(., length)
}) %>% unlist
if (is.null(nTGsPerReg_DIRECT_NET)) {nTGsPerReg_DIRECT_NET = 0}
df_nTGsPerReg_DIRECT_NET = data.frame(
  nTG_per_regulon = nTGsPerReg_DIRECT_NET,
  Method = rep("DIRECT-NET", length(nTGsPerReg_DIRECT_NET))
); rm(nTGsPerReg_DIRECT_NET)

df_nTGsPerReg_DIRECT_NET$dataset = sub("(SNAREseq2_hg38_GM12878|SNAREseq2_hg38_A549|10X_hg38_Human_bone_marrow_with_T_cell_depleted).*", "\\1", rownames(df_nTGsPerReg_DIRECT_NET))
print(aggregate(nTG_per_regulon ~ dataset, data = df_nTGsPerReg_DIRECT_NET, FUN = mean))
#                                           dataset nTG_per_regulon
# 1 10X_hg38_Human_bone_marrow_with_T_cell_depleted        6.942529
df_nTGsPerReg_DIRECT_NET$dataset = NULL

rownames(df_nTGsPerReg_DIRECT_NET) = NULL
print(head(df_nTGsPerReg_DIRECT_NET))
#   nTG_per_regulon     Method
# 1               6 DIRECT-NET
# 2               2 DIRECT-NET
# 3               2 DIRECT-NET
# 4              14 DIRECT-NET
# 5               1 DIRECT-NET
# 6               1 DIRECT-NET
print(dim(df_nTGsPerReg_DIRECT_NET))
# [1] 87  2

# Pando
nTGsPerReg_Pando = sapply(dataset, function(x) {
  if (is.null(Pando_regulons[[x]])) {return(NULL)}
  split(Pando_regulons[[x]], f = Pando_regulons[[x]]$tf) %>% sapply(., "[[", "target") %>% sapply(., unique) %>% sapply(., length)
}) %>% unlist
if (is.null(nTGsPerReg_Pando)) {nTGsPerReg_Pando = 0} # 新增
df_nTGsPerReg_Pando = data.frame(
  nTG_per_regulon = nTGsPerReg_Pando,
  Method = rep("Pando", length(nTGsPerReg_Pando))
); rm(nTGsPerReg_Pando)

df_nTGsPerReg_Pando$dataset = sub("(SNAREseq2_hg38_GM12878|SNAREseq2_hg38_A549|10X_hg38_Human_bone_marrow_with_T_cell_depleted).*", "\\1", rownames(df_nTGsPerReg_Pando))
print(aggregate(nTG_per_regulon ~ dataset, data = df_nTGsPerReg_Pando, FUN = mean))
#                                           dataset nTG_per_regulon
# 1 10X_hg38_Human_bone_marrow_with_T_cell_depleted        2.836735
df_nTGsPerReg_Pando$dataset = NULL

rownames(df_nTGsPerReg_Pando) = NULL
print(head(df_nTGsPerReg_Pando))
#   nTG_per_regulon Method
# 1               1  Pando
# 2               1  Pando
# 3               1  Pando
# 4               1  Pando
# 5               2  Pando
# 6               1  Pando
print(dim(df_nTGsPerReg_Pando))
# [1] 147   2

# scMEGA
nTGsPerReg_scMEGA = sapply(dataset, function(x) {
  if (is.null(scMEGA_regulons[[x]])) {return(NULL)}
  split(scMEGA_regulons[[x]], f = scMEGA_regulons[[x]]$tf) %>% lapply(., "[[", "gene") %>% lapply(., unique) %>% sapply(., length)
}) %>% unlist
if (is.null(nTGsPerReg_scMEGA)) {nTGsPerReg_scMEGA = 0}
nTGsPerReg_scMEGA = reshape2::melt(nTGsPerReg_scMEGA)
colnames(nTGsPerReg_scMEGA) = "nTG_per_regulon"
df_nTGsPerReg_scMEGA = data.frame(
  nTG_per_regulon = nTGsPerReg_scMEGA,
  Method = rep("scMEGA", length(nTGsPerReg_scMEGA))
); rm(nTGsPerReg_scMEGA)

df_nTGsPerReg_scMEGA$dataset = sub("(SNAREseq2_hg38_GM12878|SNAREseq2_hg38_A549|10X_hg38_Human_bone_marrow_with_T_cell_depleted).*", "\\1", rownames(df_nTGsPerReg_scMEGA))
print(aggregate(nTG_per_regulon ~ dataset, data = df_nTGsPerReg_scMEGA, FUN = mean))
#                                           dataset nTG_per_regulon
# 1 10X_hg38_Human_bone_marrow_with_T_cell_depleted             489
# 2                             SNAREseq2_hg38_A549             140
# 3                          SNAREseq2_hg38_GM12878             150
df_nTGsPerReg_scMEGA$dataset = NULL

rownames(df_nTGsPerReg_scMEGA) = NULL
print(head(df_nTGsPerReg_scMEGA))
#   nTG_per_regulon Method
# 1             489 scMEGA
# 2             489 scMEGA
# 3             489 scMEGA
# 4             489 scMEGA
# 5             489 scMEGA
# 6             489 scMEGA
print(dim(df_nTGsPerReg_scMEGA))
# [1] 108   2


### 5.2 Create box plot for number of genes per regulon
df_nTGsPerReg = as.data.frame(data.table::rbindlist(list(
  df_nTGsPerReg_scDERI,
  df_nTGsPerReg_STREAM,
  df_nTGsPerReg_GLUE,
  df_nTGsPerReg_Pando,
  df_nTGsPerReg_scMEGA,
  df_nTGsPerReg_FigR,
  df_nTGsPerReg_DIRECT_NET,
  df_nTGsPerReg_IReNA
), use.names = FALSE))
colnames(df_nTGsPerReg) = c("nTG_per_regulon", "Method")
df_nTGsPerReg$Method = factor(df_nTGsPerReg$Method, levels = method_factor_levels)
print(head(df_nTGsPerReg))
#   nTG_per_regulon Method
# 1              40 scDERI
# 2              52 scDERI
# 3              32 scDERI
# 4              13 scDERI
# 5               7 scDERI
# 6               7 scDERI
print(dim(df_nTGsPerReg))
# [1] 4045    2
print(aggregate(nTG_per_regulon ~ Method, data = df_nTGsPerReg, FUN = mean))
#       Method nTG_per_regulon
# 1     scDERI       29.480159
# 2     STREAM       95.561739
# 3       GLUE      106.685883
# 4      Pando        2.836735
# 5     scMEGA      329.740741
# 6       FigR       10.938586
# 7 DIRECT-NET        6.942529
# 8      IReNA        7.854545

boxplots_nTG_per_regulon = ggplot(df_nTGsPerReg, aes(x = Method, y = nTG_per_regulon, fill = Method)) +
  geom_boxplot(lwd = 0.1, outlier.size = 0.01) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.75, linetype = "dashed", color = "black", size = 0.2) +
  scale_fill_manual(values = colors) +
  coord_cartesian(ylim = c(NA, 510)) +
  labs(title = "TG Numbers per eRegulon") +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = NA, color = "black", size = 0.5) +
  theme(
    axis.title = element_text(size = 12, color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 9, angle = 45, hjust = 1),
    plot.title = element_text(size = 10, color = "black", hjust = 0.5, face = "bold"),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    # axis.text.x = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
    panel.grid.minor = element_line(color = "grey95", linewidth = 0.3),
    panel.background = element_rect(fill = "white", linewidth = 0.3),
    legend.position = "none"
  )
print(boxplots_nTG_per_regulon)
qs::qsave(boxplots_nTG_per_regulon, paste0(parent.dir, "boxplots_nTG_per_regulon.qsave"))

CairoPNG(filename = paste0(parent.dir, "boxplots_nTG_per_regulon.png"), width = 2.4, height = 3.2, units = "in", dpi = 300)
print(boxplots_nTG_per_regulon)
dev.off()



##### 6. Analyze Number of Regions per Regulon #####
### 6.1 Count regions per regulon for each method
nRegionsPerReg_scDERI = sapply(dataset, function(x) {
  if (is.null(scDERI_regulons[[x]])) {return(NULL)}
  GRN = unique(scDERI_regulons[[x]][, c("TF", "peak")])
  tapply(GRN$peak, GRN$TF, function(x) length(unique(x)))
}) %>% unlist
df_nRegionsPerReg_scDERI = data.frame(
  nRegion_per_regulon = nRegionsPerReg_scDERI,
  Method = rep("scDERI", length(nRegionsPerReg_scDERI))
); rm(nRegionsPerReg_scDERI)
colnames(df_nRegionsPerReg_scDERI) = c("nRegion_per_regulon", "Method")

df_nRegionsPerReg_scDERI$dataset = sub("(SNAREseq2_hg38_GM12878|SNAREseq2_hg38_A549|10X_hg38_Human_bone_marrow_with_T_cell_depleted).*", "\\1", rownames(df_nRegionsPerReg_scDERI))
print(aggregate(nRegion_per_regulon ~ dataset, data = df_nRegionsPerReg_scDERI, FUN = mean))
#                                           dataset nRegion_per_regulon
# 1 10X_hg38_Human_bone_marrow_with_T_cell_depleted            34.98913
# 2                             SNAREseq2_hg38_A549            24.77381
# 3                          SNAREseq2_hg38_GM12878            35.15789
df_nRegionsPerReg_scDERI$dataset = NULL

rownames(df_nRegionsPerReg_scDERI) = NULL
print(head(df_nRegionsPerReg_scDERI))
#   nRegion_per_regulon Method
# 1                  41 scDERI
# 2                  52 scDERI
# 3                  32 scDERI
# 4                  13 scDERI
# 5                   7 scDERI
# 6                   7 scDERI
print(dim(df_nRegionsPerReg_scDERI))
# [1] 252   2

# STREAM
nRegionsPerReg_STREAM = sapply(dataset, function(x) {
  sapply(STREAM_regulons[[x]], "[[", "peaks") %>% sapply(., length)
}) %>% unlist
df_nRegionsPerReg_STREAM = data.frame(
  nRegion_per_regulon = nRegionsPerReg_STREAM,
  Method = rep("STREAM", length(nRegionsPerReg_STREAM))
); rm(nRegionsPerReg_STREAM)
colnames(df_nRegionsPerReg_STREAM) = c("nRegion_per_regulon", "Method")

df_nRegionsPerReg_STREAM$dataset = sub("(SNAREseq2_hg38_GM12878|SNAREseq2_hg38_A549|10X_hg38_Human_bone_marrow_with_T_cell_depleted).*", "\\1", rownames(df_nRegionsPerReg_STREAM))
print(aggregate(nRegion_per_regulon ~ dataset, data = df_nRegionsPerReg_STREAM, FUN = mean))
#                                           dataset nRegion_per_regulon
# 1 10X_hg38_Human_bone_marrow_with_T_cell_depleted            75.28521
# 2                             SNAREseq2_hg38_A549            33.27737
# 3                          SNAREseq2_hg38_GM12878            39.61039
df_nRegionsPerReg_STREAM$dataset = NULL

rownames(df_nRegionsPerReg_STREAM) = NULL
print(head(df_nRegionsPerReg_STREAM))
#   nRegion_per_regulon Method
# 1                 328 STREAM
# 2                 313 STREAM
# 3                 300 STREAM
# 4                 297 STREAM
# 5                 276 STREAM
# 6                 250 STREAM
print(dim(df_nRegionsPerReg_STREAM))
# [1] 575   2

# IReNA
df_nRegionsPerReg_IReNA = data.frame(nRegion_per_regulon = NA, Method = "IReNA")
print(df_nRegionsPerReg_IReNA)
#   nRegion_per_regulon Method
# 1                  NA  IReNA

# FigR
df_nRegionsPerReg_FigR = data.frame(nRegion_per_regulon = NA, Method = "FigR")
print(df_nRegionsPerReg_FigR)
#   nRegion_per_regulon Method
# 1                  NA   FigR

# GLUE
nRegionsPerReg_GLUE = sapply(dataset, function(x) {
  if (is.null(GLUE_regulons[[x]])) {return(NULL)}
  split(GLUE_regulons[[x]], f = GLUE_regulons[[x]]$TF) %>% sapply(., "[[", "Region") %>% sapply(., unique) %>% sapply(., length)
}) %>% unlist
df_nRegionsPerReg_GLUE = data.frame(
  nRegion_per_regulon = nRegionsPerReg_GLUE,
  Method = rep("GLUE", length(nRegionsPerReg_GLUE))
); rm(nRegionsPerReg_GLUE)
colnames(df_nRegionsPerReg_GLUE) = c("nRegion_per_regulon", "Method")

df_nRegionsPerReg_GLUE$dataset = sub("(SNAREseq2_hg38_GM12878|SNAREseq2_hg38_A549|10X_hg38_Human_bone_marrow_with_T_cell_depleted).*", "\\1", rownames(df_nRegionsPerReg_GLUE))
print(aggregate(nRegion_per_regulon ~ dataset, data = df_nRegionsPerReg_GLUE, FUN = mean))
#                                           dataset nRegion_per_regulon
# 1 10X_hg38_Human_bone_marrow_with_T_cell_depleted            103.8783
# 2                             SNAREseq2_hg38_A549            163.8858
# 3                          SNAREseq2_hg38_GM12878            215.6446
df_nRegionsPerReg_GLUE$dataset = NULL

rownames(df_nRegionsPerReg_GLUE) = NULL
print(head(df_nRegionsPerReg_GLUE))
#   nRegion_per_regulon Method
# 1                 203   GLUE
# 2                 110   GLUE
# 3                 257   GLUE
# 4                 274   GLUE
# 5                  74   GLUE
# 6                 275   GLUE
print(dim(df_nRegionsPerReg_GLUE))
# [1] 1353    2

# DIRECT-NET
df_nRegionsPerReg_DIRECT_NET = data.frame(nRegion_per_regulon = NA, Method = "DIRECT-NET")
print(df_nRegionsPerReg_DIRECT_NET)
#   nRegion_per_regulon     Method
# 1                  NA DIRECT-NET

# Pando
nRegionsPerReg_Pando = sapply(dataset, function(x) {
  if (is.null(Pando_regulons[[x]])) {return(NULL)}
  split(Pando_regulons[[x]], f = Pando_regulons[[x]]$tf) %>% sapply(., "[[", "regions") %>% sapply(., unique) %>% sapply(., length)
}) %>% unlist
if (is.null(nRegionsPerReg_Pando)) {nRegionsPerReg_Pando = NA}
df_nRegionsPerReg_Pando = data.frame(
  nRegion_per_regulon = nRegionsPerReg_Pando,
  Method = rep("Pando", length(nRegionsPerReg_Pando))
); rm(nRegionsPerReg_Pando)

df_nRegionsPerReg_Pando$dataset = sub("(SNAREseq2_hg38_GM12878|SNAREseq2_hg38_A549|10X_hg38_Human_bone_marrow_with_T_cell_depleted).*", "\\1", rownames(df_nRegionsPerReg_Pando))
print(aggregate(nRegion_per_regulon ~ dataset, data = df_nRegionsPerReg_Pando, FUN = mean))
#                                           dataset nRegion_per_regulon
# 1 10X_hg38_Human_bone_marrow_with_T_cell_depleted            2.836735
df_nRegionsPerReg_Pando$dataset = NULL

rownames(df_nRegionsPerReg_Pando) = NULL
print(head(df_nRegionsPerReg_Pando))
#   nRegion_per_regulon Method
# 1                   1  Pando
# 2                   1  Pando
# 3                   1  Pando
# 4                   1  Pando
# 5                   2  Pando
# 6                   1  Pando
print(dim(df_nRegionsPerReg_Pando))
# [1] 147   2

# scMEGA
df_nRegionsPerReg_scMEGA = data.frame(nRegion_per_regulon = NA, Method = "scMEGA")
print(df_nRegionsPerReg_scMEGA)
#   nRegion_per_regulon Method
# 1                  NA scMEGA


### 6.2 Create box plot for number of regions per regulon
df_nRegionsPerReg = rbind(
  df_nRegionsPerReg_scDERI,
  df_nRegionsPerReg_STREAM,
  df_nRegionsPerReg_GLUE,
  df_nRegionsPerReg_Pando,
  df_nRegionsPerReg_scMEGA,
  df_nRegionsPerReg_FigR,
  df_nRegionsPerReg_DIRECT_NET,
  df_nRegionsPerReg_IReNA
)
df_nRegionsPerReg$Method = factor(df_nRegionsPerReg$Method, levels = method_factor_levels)
print(head(df_nRegionsPerReg))
#   nRegion_per_regulon Method
# 1                  41 scDERI
# 2                  52 scDERI
# 3                  32 scDERI
# 4                  13 scDERI
# 5                   7 scDERI
# 6                   7 scDERI
print(dim(df_nRegionsPerReg))
# [1] 2331    2

boxplots_nRegion_per_regulon = ggplot(df_nRegionsPerReg, aes(x = Method, y = nRegion_per_regulon, fill = Method)) +
  geom_boxplot(lwd = 0.1, outlier.size = 0.01) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.75, linetype = "dashed", color = "black", size = 0.2) +
  scale_fill_manual(values = colors) +
  coord_cartesian(ylim = c(NA, 510)) +
  labs(title = "Peak Numbers per eRegulon") +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = NA, color = "black", size = 0.5) +
  theme(
    axis.title = element_text(size = 12, color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 9, angle = 45, hjust = 1),
    plot.title = element_text(size = 10, color = "black", hjust = 0.5, face = "bold"),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    # axis.text.x = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
    panel.grid.minor = element_line(color = "grey95", linewidth = 0.3),
    panel.background = element_rect(fill = "white", linewidth = 0.3),
    legend.position = "none"
  )
print(boxplots_nRegion_per_regulon)
qs::qsave(boxplots_nRegion_per_regulon, paste0(parent.dir, "boxplots_nRegion_per_regulon.qsave"))

# CairoPNG(filename = paste0(parent.dir, "boxplots_nRegion_per_regulon.png"), width = 300, height = 300, res = 150)
CairoPNG(filename = paste0(parent.dir, "boxplots_nRegion_per_regulon.png"), width = 2.4, height = 3.2, units = "in", dpi = 300)
print(boxplots_nRegion_per_regulon)
dev.off()



##### 7. Visualize TF Coverage Benchmark Results #####
barplots_benchmark_TF_coverage_list = list()
for (i in seq_along(dataset)) {
  message("Creating TF coverage plot for dataset ", dataset[i])
  df = data.frame(
    Method = method_factor_levels,
    f.score = TF_coverage_benchmark[[i]][plyr::revalue(method_factor_levels, c("GLUE" = "scGLUE", "DIRECT-NET" = "DIRECT_NET")), "f.score"]
  )
  rownames(df) = df$Method
  df$Method = factor(df$Method, levels = method_factor_levels)

  barplots_benchmark_TF_coverage_list[[i]] = ggplot(df, aes(x = Method, y = f.score, fill = Method)) +
    geom_bar(stat = "identity", fill = colors, lwd = 0.1) +
    labs(title = paste0(dataset_simplify[i], " TF coverage")) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = NA, color = "black", size = 0.5) +
    theme(
      axis.title = element_text(size = 12, color = "black"),
      axis.text = element_text(size = 12, color = "black"),
      axis.text.x = element_text(size = 9, angle = 45, hjust = 1),
      plot.title = element_text(size = 10, color = "black", hjust = 0.5, face = "bold"),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      # axis.text.x = element_blank(),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
      panel.grid.minor = element_line(color = "grey95", linewidth = 0.3),
      panel.background = element_rect(fill = "white", linewidth = 0.3),
      legend.position = "none"
    )
  CairoPNG(filename = paste0(parent.dir, "barplots_benchmark_TF_coverage_", dataset_simplify[i], ".png"), width = 2.4, height = 3.2, units = "in", dpi = 300)
  print(barplots_benchmark_TF_coverage_list[[i]])
  dev.off()
}
names(barplots_benchmark_TF_coverage_list) = dataset
qs::qsave(barplots_benchmark_TF_coverage_list, paste0(parent.dir, "barplots_benchmark_TF_coverage_list.qsave"))

# barplots_benchmark_TF_coverage = ggpubr::ggarrange(plotlist = barplots_benchmark_TF_coverage_list, nrow = 1)
# barplots_benchmark_TF_coverage
# qs::qsave(barplots_benchmark_TF_coverage, paste0(parent.dir, "barplots_benchmark_TF_coverage.qsave"))

# CairoPNG(filename = paste0(parent.dir, "barplots_benchmark_TF_coverage.png"), width = 1200, height = 400, res = 150)
# barplots_benchmark_TF_coverage
# dev.off()



##### 8. Visualize TF-to-Region Linkage Benchmark Results #####
barplots_benchmark_TF_region_list = list()
for (i in seq_along(dataset)) {
  message("Creating TF-region plot for dataset ", dataset[i])
  df = data.frame(
    Method = as.character(TF_region_benchmark_all[[i]]$method),
    precision = as.numeric(TF_region_benchmark_all[[i]][, "precision"])
  )
  df = dplyr::bind_rows(df, data.frame(Method = c("scMEGA", "FigR", "DIRECT-NET", "IReNA"), precision = c(0, 0, 0, 0)))
  df$Method[df$Method == "scGLUE"] = "GLUE"
  df$Method = factor(df$Method, levels = method_factor_levels)
  df$precision = as.numeric(df$precision)

  barplots_benchmark_TF_region_list[[i]] = ggplot(df, aes(x = Method, y = precision, fill = Method)) +
    geom_boxplot(lwd = 0.1, outlier.size = 0.01) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.75, linetype = "dashed", color = "black", size = 0.2) +
    scale_fill_manual(values = colors) +
    labs(title = paste0(dataset_simplify[i], " TF-peak")) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = NA, color = "black", size = 0.5) +
    theme(
      axis.title = element_text(size = 12, color = "black"),
      axis.text = element_text(size = 12, color = "black"),
      axis.text.x = element_text(size = 9, angle = 45, hjust = 1),
      plot.title = element_text(size = 10, color = "black", hjust = 0.5, face = "bold"),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      # axis.text.x = element_blank(),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
      panel.grid.minor = element_line(color = "grey95", linewidth = 0.3),
      panel.background = element_rect(fill = "white", linewidth = 0.3),
      legend.position = "none"
    )
  CairoPNG(filename = paste0(parent.dir, "barplots_benchmark_TF_region_", dataset_simplify[i], ".png"), width = 2.4, height = 3.2, units = "in", dpi = 300)
  print(barplots_benchmark_TF_region_list[[i]])
  dev.off()
}
names(barplots_benchmark_TF_region_list) = dataset
qs::qsave(barplots_benchmark_TF_region_list, paste0(parent.dir, "barplots_benchmark_TF_region_list.qsave"))

# barplots_benchmark_TF_region = ggpubr::ggarrange(plotlist = barplots_benchmark_TF_region_list, nrow = 1)
# barplots_benchmark_TF_region
# qs::qsave(barplots_benchmark_TF_region, paste0(parent.dir, "barplots_benchmark_TF_region.qsave"))

# CairoPNG(filename = paste0(parent.dir, "barplots_benchmark_TF_region.png"), width = 1200, height = 400, res = 150)
# barplots_benchmark_TF_region
# dev.off()



##### 9. Visualize Region-to-Gene Linkage Benchmark Results #####
barplots_benchmark_region_TG_list = list()
for (i in seq_along(dataset)) {
  message("Creating region-gene plot for dataset ", dataset[i])
  region_TG_benchmark_plus = rbind(
    region_TG_benchmark[[i]],
    matrix(0, nrow = 4, ncol = ncol(region_TG_benchmark[[i]]),
           dimnames = list(c("IReNA", "FigR", "DIRECT_NET", "scMEGA"), NULL))
  )

  df = data.frame(
    Method = method_factor_levels,
    f.score = as.numeric(region_TG_benchmark_plus[plyr::revalue(method_factor_levels, c("GLUE" = "scGLUE", "DIRECT-NET" = "DIRECT_NET")), "f.score"])
  )
  df$Method = factor(df$Method, method_factor_levels)

  barplots_benchmark_region_TG_list[[i]] = ggplot(df, aes(x = Method, y = f.score)) +
    geom_bar(stat = "identity", fill = colors, lwd = 0.1) +
    labs(title = paste0(dataset_simplify[i], " peak-TG")) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = NA, color = "black", size = 0.5) +
    theme(
      axis.title = element_text(size = 12, color = "black"),
      axis.text = element_text(size = 12, color = "black"),
      axis.text.x = element_text(size = 9, angle = 45, hjust = 1),
      plot.title = element_text(size = 10, color = "black", hjust = 0.5, face = "bold"),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      # axis.text.x = element_blank(),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
      panel.grid.minor = element_line(color = "grey95", linewidth = 0.3),
      panel.background = element_rect(fill = "white", linewidth = 0.3),
      legend.position = "none"
    )
  CairoPNG(filename = paste0(parent.dir, "barplots_benchmark_region_TG_", dataset_simplify[i], ".png"), width = 2.4, height = 3.2, units = "in", dpi = 300)
  print(barplots_benchmark_region_TG_list[[i]])
  dev.off()
}
names(barplots_benchmark_region_TG_list) = dataset
qs::qsave(barplots_benchmark_region_TG_list, paste0(parent.dir, "barplots_benchmark_region_TG_list.qsave"))

# barplots_benchmark_region_TG = ggpubr::ggarrange(plotlist = barplots_benchmark_region_TG_list, nrow = 1, ncol = 3)
# barplots_benchmark_region_TG
# qs::qsave(barplots_benchmark_region_TG, paste0(parent.dir, "barplots_benchmark_region_TG.qsave"))

# CairoPNG(filename = paste0(parent.dir, "barplots_benchmark_region_TG.png"), width = 1200, height = 400, res = 150)
# barplots_benchmark_region_TG
# dev.off()



##### 10. Combine All Panels into Final Figure #####
# CairoPNG(filename = paste0(parent.dir, "benchmark_plots.png"), width = 2400, height = 1800, res = 150)
# ggpubr::ggarrange(
#   ggpubr::ggarrange(
#     NA,
#     ggpubr::ggarrange(boxplots_nTFs, boxplots_nTG_per_regulon, boxplots_nRegion_per_regulon, ncol = 1, labels = LETTERS[2:4]),
#     ggpubr::ggarrange(plotlist = barplots_benchmark_TF_coverage_list, nrow = 1, labels = LETTERS[5:7]),
#     nrow = 1, widths = c(1, 1, 2)
#   ),
#   ggpubr::ggarrange(
#     ggpubr::ggarrange(plotlist = barplots_benchmark_TF_region_list, nrow = 1, labels = LETTERS[8:10]),
#     ggpubr::ggarrange(plotlist = barplots_benchmark_region_TG_list, nrow = 1, labels = LETTERS[11:13]),
#     nrow = 1, widths = c(1, 1)
#   ),
#   ncol = 1
# )
# dev.off()

save.image(paste0(parent.dir, "benchmark-visualization.RData"))
# load("/mnt/i/scDERI_Paper/3-benchmark/9-evaluation/4-visualization/benchmark-visualization.RData")

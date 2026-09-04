############################################################
#                                                          #
#   Apply scDERI on analysed dataset - survival analysis   #
#                                                          #
############################################################


# ================================================
# 1. Download TCGA PAAD Dataset from UCSC Xena
# ================================================
# Access the GDC TCGA Pancreatic Cancer (PAAD) dataset (18 datasets) in UCSC Xena
# (https://xenabrowser.net/datapages/?cohort=GDC%20TCGA%20Pancreatic%20Cancer%20(PAAD)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443).
# 1. Navigate to the "gene expression RNAseq - STAR TPM (n=183)" section in GDC Hub (https://xenabrowser.net/datapages/?dataset=TCGA-PAAD.star_tpm.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443)
# and download the TPM expression matrix (https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-PAAD.star_tpm.tsv.gz).
# 2. Navigate to the "phenotype - Phenotype (n=196)" section in GDC Hub (https://xenabrowser.net/datapages/?dataset=TCGA-PAAD.clinical.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443)
# and download the clinical information for each sample (https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-PAAD.clinical.tsv.gz).
# 3. Navigate to the "phenotype - survival data (n=195)" section in GDC Hub (https://xenabrowser.net/datapages/?dataset=TCGA-PAAD.survival.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443)
# and download the survival information for each sample (https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-PAAD.survival.tsv.gz).
# Download all files to I:\scDERI_Paper\4-apply\survivaland extract them to this folder for viewing in Excel.


# ================================================
# 2. Environment Setup
# ================================================
# Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
library(data.table)
library(survival)
library(survminer)
library(org.Hs.eg.db)
library(stringr)
library(survRM2)
library(parallel)

sys = "wsl2"
# sys = "win11"
if (sys == "wsl2") {
  Working_directory = "/mnt/i/scDERI_Paper/4-apply/"
} else if (sys == "win11") {
  Working_directory = "I:/scDERI_Paper/4-apply/"
} else {stop("Error.")}
setwd(Working_directory)


# ================================================
# 3. Load TCGA PAAD Data
# ================================================
# Load TPM expression matrix
tpm = data.table::fread("./survival/TCGA-PAAD.star_tpm.tsv.gz")
tpm[1:5, 1:5]
#            Ensembl_ID TCGA-IB-A6UG-01A TCGA-3A-A9IU-01A TCGA-FB-A545-01A TCGA-IB-AAUQ-01A
#                <char>            <num>            <num>            <num>            <num>
# 1: ENSG00000000003.15         5.148837       5.72554213         4.667376        4.4337407
# 2:  ENSG00000000005.6         0.000000       0.06170708         0.000000        0.3015876
# 3: ENSG00000000419.13         5.806968       6.45111218         6.630772        6.1337464
# 4: ENSG00000000457.14         2.621876       3.25945328         2.963567        2.6911545
# 5: ENSG00000000460.17         1.457909       2.13405661         2.222588        1.7894789

# Load survival data (OS: overall survival, 1=death, 0=censored)
survival = data.table::fread("./survival/TCGA-PAAD.survival.tsv.gz")
head(survival)
#              sample OS.time    OS     _PATIENT
#              <char>   <num> <int>       <char>
# 1: TCGA-HZ-7922-01A       4     0 TCGA-HZ-7922
# 2: TCGA-HZ-A8P1-01A       7     0 TCGA-HZ-A8P1
# 3: TCGA-IB-AAUM-01A       8     0 TCGA-IB-AAUM
# 4: TCGA-RL-AAAS-01A       9     0 TCGA-RL-AAAS
# 5: TCGA-US-A77G-01A      12     1 TCGA-US-A77G
# 6: TCGA-FB-AAQ3-01A      31     1 TCGA-FB-AAQ3

# Load clinical data
clinical = data.table::fread("./survival/TCGA-PAAD.clinical.tsv.gz")
dim(clinical) # [1] 196  89


# ================================================
# 4. Gene ID Conversion (Ensembl to Symbol)
# ================================================
# Extract Ensembl IDs (remove version numbers)
ensembl_ids = str_split(tpm$Ensembl_ID, "\\.", simplify = TRUE)[,1]
tpm$Ensembl_ID = ensembl_ids

# Convert Ensembl IDs to gene symbols using clusterProfiler
gene_symbol = clusterProfiler::bitr(
  ensembl_ids,
  fromType = "ENSEMBL",
  toType = c("SYMBOL", "ENTREZID"),
  OrgDb = "org.Hs.eg.db"
)

# Merge symbol information into TPM matrix
tpm = merge(tpm, gene_symbol, by.x = "Ensembl_ID", by.y = "ENSEMBL", all.x = TRUE)
tpm = tpm %>% dplyr::relocate(tail(names(.), 2), .after = 1)
tpm = na.omit(tpm, cols = "SYMBOL")
tpm = tpm[rowSums(tpm[, -(1:3)]) > 0, ]
tpm = aggregate(tpm[, 4:ncol(tpm)], by = list(SYMBOL = tpm$SYMBOL), FUN = mean)


# ================================================
# 5. Load scDERI-Predicted genes
# ================================================
# For test
# interest_cell_type_peGRN = qs::qread("./5 Infer GRN/interest_cell_type_peGRN.qsave")
# TFs = sort(unique(interest_cell_type_peGRN[["KP4"]][["loose_filtered_simplified"]][["TF"]]))
# TGs = sort(unique(interest_cell_type_peGRN[["KP4"]][["loose_filtered_simplified"]][["TG"]]))

TGs = qs::qread("./1.2 Data Preprocessing - Identify TGs By Annotation/results_identify_TGs.qsave")$TGs
TFs = qs::qread("./1.1 Data Preprocessing - Get TFs From Database/TFs.qsave")

TFs = TFs[TFs %in% tpm$SYMBOL]
TGs = TGs[TGs %in% tpm$SYMBOL]


# ================================================
# 6. Survival Analysis for Each gene
# ================================================
survival_info = list()

for (d in c("./survival/TF", "./survival/TG")) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}; rm(d)

# set p value threshold
p_value_threshold = 0.05
# p_value_threshold = 0.01 # For test

gene_type_set = c("TF", "TG")
# gene_type_set = "TF" # For test
# gene_type_set = "TG" # For test

draw = c(TRUE, FALSE)
names(draw) = gene_type_set

t1 = Sys.time(); for (gene_type in gene_type_set) {
  survival_info[[gene_type]] = list()
  if (gene_type == "TF") {genes = TFs} else {genes = TGs}
  survival_info[[gene_type]][["result"]] = data.frame(matrix(
    NA, nrow = length(genes), ncol = 3,
    dimnames = list(genes, c("rmst_area_ratio_of_group_high_to_low", "p_of_rmst_area_ratio_of_group_high_to_low", "features"))
  ))
  n_gene = length(genes)

  t2 = Sys.time()
  for (target_gene in genes) {
    survival_info[[gene_type]][[target_gene]] = list()
    survival_info[[gene_type]][[target_gene]][[gene_type]] = target_gene

    # Extract expression data for target gene
    gene_expression = as.numeric(tpm[tpm$SYMBOL == target_gene, -(1:3)])
    survival_info[[gene_type]][[target_gene]][["gene_expression"]] = gene_expression

    # Extract sample IDs
    sample_ids = colnames(tpm)[-(1:3)]
    survival_info[[gene_type]][[target_gene]][["sample_ids"]] = sample_ids

    # Create expression data frame
    expr_df = data.frame(sample = sample_ids, expression = gene_expression)

    # Merge with survival data
    surv_data = merge(survival, expr_df, by = "sample")

    # Group samples by median expression
    surv_data$group = ifelse(surv_data$expression > median(surv_data$expression, na.rm = TRUE), "High", "Low")
    survival_info[[gene_type]][[target_gene]][["survival_data"]] = surv_data

    # Calculate restricted mean survival time (RMST)
    # Set truncation time tau to 85th percentile of survival times
    tau = quantile(surv_data$OS.time, 0.85, na.rm = TRUE)
    survival_info[[gene_type]][[target_gene]][["stage_time_tau_for_RMST_calculation"]] = tau

    # Perform RMST analysis
    # rmst_result = survRM2::rmst2(
    #   time = surv_data$OS.time,
    #   status = surv_data$OS,
    #   arm = ifelse(surv_data$group == "High", 1, 0),
    #   tau = tau,
    #   alpha = p_value_threshold
    # )
    rmst_result = tryCatch({
      survRM2::rmst2(
        time = surv_data$OS.time,
        status = surv_data$OS,
        arm = ifelse(surv_data$group == "High", 1, 0),
        tau = tau,
        alpha = p_value_threshold
      )
    }, error = function(e) {
      return(NULL)
    })
    if (is.null(rmst_result)) {
      survival_info[[gene_type]][[target_gene]] = NULL
      next
    }

    survival_info[[gene_type]][[target_gene]][["RMST_result"]] = rmst_result

    # Extract maximum p-value from RMST results
    p_max = max(rmst_result$unadjusted.result[, "p"])

    # Determine if subsequent analysis is needed (p < p_value_threshold)
    if (p_max < p_value_threshold) {
      subsequent_analysis = TRUE
    } else {
      subsequent_analysis = FALSE
    }

    time_period = format(Sys.time() - t2)
    time_period_spilt = strsplit(time_period, " ")[[1]]
    position =  match(target_gene, genes)
    message(gene_type, " - ", target_gene, " - ", position, "/", n_gene, " - ", time_period, " - ",
            "Expected ", round(as.numeric(time_period_spilt[1]) / position * (n_gene - position), digits = 3), " ", time_period_spilt[2], " remaining",
            " - draw: ", subsequent_analysis)
    survival_info[[gene_type]][[target_gene]][["subsequent_analysis"]] = subsequent_analysis

    if (subsequent_analysis) {
      # Calculate RMST ratio (High/Low)
      rmst_high = rmst_result$RMST.arm1$rmst["Est."]
      rmst_low = rmst_result$RMST.arm0$rmst["Est."]
      rmst_ratio = rmst_high / rmst_low

      # Store results
      survival_info[[gene_type]][[target_gene]][["rmst_area_ratio_of_group_high_to_low"]] = rmst_ratio
      survival_info[[gene_type]][["result"]][target_gene, "rmst_area_ratio_of_group_high_to_low"] = rmst_ratio
      survival_info[[gene_type]][[target_gene]][["p_of_rmst_area_ratio_of_group_high_to_low"]] = p_max
      survival_info[[gene_type]][["result"]][target_gene, "p_of_rmst_area_ratio_of_group_high_to_low"] = p_max

      # Determine survival association direction
      if (rmst_ratio < 1) {
        survival_info[[gene_type]][["result"]][target_gene, "features"] = "high" # High expression of this TF is detrimental to survival.
      } else {
        survival_info[[gene_type]][["result"]][target_gene, "features"] = "low" # Low expression of this TF is detrimental to survival.
      }

      # Perform Kaplan-Meier survival analysis
      fit = survfit(Surv(OS.time, OS) ~ group, data = surv_data)
      survival_info[[gene_type]][[target_gene]][["survival_fit"]] = fit

      if (draw[gene_type]) {
        # Generate survival plot
        survival_curve_plot = ggsurvplot(
          fit, data = surv_data,
          pval = TRUE,
          risk.table = TRUE,
          palette = c("red", "blue"),
          legend.title = target_gene,
          xlab = "Time (Days)"
        )
        survival_info[[gene_type]][[target_gene]][["survival_curve_plot"]] = survival_curve_plot

        # Save survival plot
        ggsave(paste0("./survival/", gene_type, "/", target_gene, "_survival_curve_plot.png"),
               plot = survminer::arrange_ggsurvplots(list(survival_curve_plot), print = TRUE, ncol = 1, nrow = 1),
               width = 5, height = 5, dpi = 300)
      }
    }
  }
};t2 = Sys.time(); print(t2 - t1)
# 7.687771 mins

qs::qsave(survival_info, "./survival/survival_info.qsave")
# survival_info = qs::qread("./survival/survival_info.qsave")

survival_info_TF = survival_info[["TF"]][["result"]]
survival_info_TF = na.omit(survival_info_TF)
qs::qsave(survival_info_TF, "./survival/survival_info_TF.qsave")
# survival_info_TF = qs::qread("./survival/survival_info_TF.qsave")

survival_info_TG = survival_info[["TG"]][["result"]]
survival_info_TG = na.omit(survival_info_TG)
qs::qsave(survival_info_TG, "./survival/survival_info_TG.qsave")
# survival_info_TG = qs::qread("./survival/survival_info_TG.qsave")

table(survival_info_TF$features)
# high  low
#   28    7

table(survival_info_TG$features)
# high  low
#  402  112

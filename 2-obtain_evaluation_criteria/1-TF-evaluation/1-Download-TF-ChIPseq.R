#####################################################################################
#                                                                                   #
#   Download ChIP-seq data to validate the identified TFs and TF-region relations   #
#                                                                                   #
#####################################################################################


# ======================== 1. Setup and Parameters ========================
library(dplyr)
library(readxl)
library(openxlsx)
work_dir <- "/mnt/i/scDERI_Paper/2-obtain_evaluation_criteria/1-TF-evaluation/"
setwd(work_dir)


# ======================== 2. ENCODE Metadata Acquisition & TF Peak Selection ========================
# Retrieve metadata of TF ChIP-seq peak files from the ENCODE database (https://www.encodeproject.org/search/?type=Files);
# for each TF, retain the record with the largest file size for benchmarking
# All metadata TSV files must be manually downloaded to the current working directory and named by dataset

# ---------- 2.1 A549 Cell Line ----------
# Filtering criteria:
#   - Format: output category = annotation;
#             Output type = IDR ranked peaks, peaks and background as input for IDR
#             File format = bed;
#             File format type = idr_ranked_peak, narrowPeak
#   - Assay: Assay title = TF ChIP-seq
#   - Biosample: Biosample = A549
#   - Analysis: Genome assembly = GRCh38
#   - Provenance: Project = ENCODE
#   - Quality: Status = released, archived, revoked
# Filter results page: https://www.encodeproject.org/report/?type=File&output_category=annotation&output_type=peaks+and+background+as+input+for+IDR&output_type=IDR+ranked+peaks&file_format=bed&file_format_type=idr_ranked_peak&file_format_type=narrowPeak&file_type=bed+idr_ranked_peak&file_type=bed+narrowPeak&assay_title=TF+ChIP-seq&biosample_ontology.term_name=A549&biosample_ontology.organ_slims=lung&assembly=GRCh38&award.project=ENCODE&status=released
# TSV metadata download link: https://www.encodeproject.org/report.tsv?type=File&output_category=annotation&file_format=bed&file_format_type=idr_ranked_peak&file_format_type=narrowPeak&biosample_ontology.term_name=A549&assembly=GRCh38&biosample_ontology.organ_slims=lung&award.project=ENCODE&status=released&status=revoked&status=archived&output_type=IDR+ranked+peaks&output_type=peaks+and+background+as+input+for+IDR&file_type=bed+narrowPeak&file_type=bed+idr_ranked_peak&assay_title=TF+ChIP-seq
# Local saved filename: A549_TF_ChIP_seq_ENCODE_all.tsv

# ---------- 2.2 GM12878 Cell Line ----------
# Filtering criteria:
#   - Format: output category = annotation;
#             Output type = IDR ranked peaks, peaks and background as input for IDR
#             File format = bed;
#             File format type = idr_ranked_peak, narrowPeak
#   - Assay: Assay title = TF ChIP-seq
#   - Biosample: Biosample = GM12878
#   - Analysis: Genome assembly = GRCh38
#   - Provenance: Project = ENCODE
#   - Quality: Status = released, archived, revoked
# Filter results page: https://www.encodeproject.org/report/?type=File&output_category=annotation&file_format=bed&output_type=IDR+ranked+peaks&output_type=peaks+and+background+as+input+for+IDR&file_format_type=idr_ranked_peak&assay_title=TF+ChIP-seq&biosample_ontology.term_name=GM12878&file_format_type=narrowPeak&assembly=GRCh38&award.project=ENCODE&status=released&status=archived&status=revoked
# TSV metadata download link: https://www.encodeproject.org/report.tsv?type=File&output_category=annotation&file_format=bed&output_type=IDR+ranked+peaks&output_type=peaks+and+background+as+input+for+IDR&file_format_type=idr_ranked_peak&assay_title=TF+ChIP-seq&biosample_ontology.term_name=GM12878&file_format_type=narrowPeak&assembly=GRCh38&award.project=ENCODE&status=released&status=archived&status=revoked
# Local saved filename: GM12878_TF_ChIP_seq_ENCODE_all.tsv

# ---------- 2.3 Bone Marrow Mixed Cell Lines ----------
# Filtering criteria:
#   - Format: output category = annotation;
#             Output type = peaks, peaks and background as input for IDR, IDR ranked peaks
#             File format = bed;
#             File format type = idr_ranked_peak, narrowPeak
#   - Assay: Assay title = TF ChIP-seq
#   - Biosample: Biosample = BLaER1, NCI-H929, OCI-LY1, OCI-LY3, erythroblast, SH-SY5Y, NB4
#   - Analysis: Genome assembly = GRCh38
#   - Provenance: Project = ENCODE
#   - Quality: Status = released, archived, revoked
# Filter results page: https://www.encodeproject.org/report/?type=File&file_format=bed&output_category=annotation&output_type=peaks&output_type=IDR+ranked+peaks&output_type=peaks+and+background+as+input+for+IDR&file_format_type=narrowPeak&file_format_type=idr_ranked_peak&file_type=bed+narrowPeak&file_type=bed+idr_ranked_peak&assay_title=TF+ChIP-seq&biosample_ontology.term_name=BLaER1&biosample_ontology.term_name=NCI-H929&biosample_ontology.term_name=OCI-LY1&biosample_ontology.term_name=OCI-LY3&biosample_ontology.term_name=erythroblast&biosample_ontology.term_name=SH-SY5Y&biosample_ontology.term_name=NB4&assembly=GRCh38&award.project=ENCODE&status=archived&status=released&status=revoked
# TSV metadata download link: https://www.encodeproject.org/report.tsv?type=File&file_format=bed&output_category=annotation&output_type=peaks&output_type=IDR+ranked+peaks&output_type=peaks+and+background+as+input+for+IDR&file_format_type=narrowPeak&file_format_type=idr_ranked_peak&file_type=bed+narrowPeak&file_type=bed+idr_ranked_peak&assay_title=TF+ChIP-seq&biosample_ontology.term_name=BLaER1&biosample_ontology.term_name=NCI-H929&biosample_ontology.term_name=OCI-LY1&biosample_ontology.term_name=OCI-LY3&biosample_ontology.term_name=erythroblast&biosample_ontology.term_name=SH-SY5Y&biosample_ontology.term_name=NB4&assembly=GRCh38&award.project=ENCODE&status=archived&status=released&status=revoked
# Local saved filename: Bone_Marrow_TF_ChIP_seq_ENCODE_all.tsv

# Read metadata by dataset; retain the record with the largest file size for each TF
TF_ChIP_seq_ENCODE <- list()
for (dataset in c("A549", "GM12878", "Bone-marrow")) {
  TF_ChIP_seq_ENCODE[[paste0(dataset, "-hg38-ChIP-seq")]] <-
    data.table::fread(paste0("./", dataset, "_TF_ChIP_seq_ENCODE_all.tsv")) %>%
    as.data.frame() %>%
    dplyr::group_by(`Target label`) %>%
    dplyr::slice_max(`File size`, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    as.data.frame()
}

# View names of the output datasets
names(TF_ChIP_seq_ENCODE)
# [1] "A549-hg38-ChIP-seq"        "GM12878-hg38-ChIP-seq"     "Bone-marrow-hg38-ChIP-seq"

# Export the filtered metadata as an Excel file
writexl::write_xlsx(
  x = TF_ChIP_seq_ENCODE,
  path = "./encode-chipseq-tf-enhancer-benchmarking.xlsx"
)


# ======================== 3. Load Excel File Information ========================
message("\nBegan reading the Excel sheets one by one ...")
sheet_names <- readxl::excel_sheets("encode-chipseq-tf-enhancer-benchmarking.xlsx")
sheet_names
# [1] "A549-hg38-ChIP-seq"        "GM12878-hg38-ChIP-seq"     "Bone-marrow-hg38-ChIP-seq"
message(length(sheet_names), " sheets exist in this Excel file")
# 3 sheets exist in this Excel file


# ======================== 4. Initialize Data Storage ========================
data_info <- list()
multithreaded_download = TRUE
options(timeout = 10000)


# ======================== 5. Process Each Sheet ========================
for(i in sheet_names) {
  # c

  # 5.1 Read sheet
  message("----> Reading the sheet ", i, " ...")
  data_info[[i]] = readxl::read_excel("encode-chipseq-tf-enhancer-benchmarking.xlsx", sheet = i)
  message("----> Found ", nrow(data_info[[i]]), " BED files in sheet ", i)
  # ----> Found 85 BED files in sheet A549-hg38-ChIP-seq
  # ----> Found 163 BED files in sheet GM12878-hg38-ChIP-seq
  # ----> Found 10 BED files in sheet Bone-marrow-hg38-ChIP-seq

  # 5.2 Create output directory
  dir.create(paste0(work_dir, i))


  # ================= 6. Download Files =================
  message("\n----> Downloading datasets from sheet ", i, " ...")
  if (multithreaded_download) {
    library(parallel)
    n_cores <- detectCores() - 1
    message("Using ", n_cores, " cores for parallel download")
    cl <- makeCluster(n_cores, type = "FORK")
    result <- parLapply(cl, 1:nrow(data_info[[i]]), function(j) {
      tryCatch({
        message("--------> Downloading ChIP-seq data of ID ", j, " ...")

        url <- data_info[[i]]$`Download URL`[j]
        target <- data_info[[i]]$`Target label`[j]
        bed.file <- strsplit(url, split = "/")
        suffix <- bed.file[[1]][length(bed.file[[1]])]

        download.result <- try(download.file(
          url = paste0("https://www.encodeproject.org", url),
          destfile = paste0(work_dir, i, "/", target, "_", suffix),
          quiet = FALSE
        ))

        if (inherits(download.result, "try-error")) {
          message("--------> Failed to download file for ", target)
          return(paste("Failed:", target))
        } else {
          message("--------> Finished downloading ChIP-seq data of target ", target)
          return(paste("Success:", target))
        }
      }, error = function(e) {
        message("--------> Error in processing ID ", j, ": ", e$message)
        return(paste("Error:", j))
      })
    })
    stopCluster(cl)
    success_count <- sum(grepl("Success:", result))
    fail_count <- sum(grepl("Failed:|Error:", result))
    message("\nDownload summary: ", success_count, " succeeded, ", fail_count, " failed")
  } else {
    lapply(1:nrow(data_info[[i]]), function(j) {
      # j=1
      message("--------> Downloading ChIP-seq data of ID ", j, "/", nrow(data_info[[i]])," ...")
      url <- data_info[[i]]$`Download URL`[j]
      target <- data_info[[i]]$`Target label`[j]
      bed.file <- strsplit(url, split = "/")
      suffix <- bed.file[[1]][length(bed.file[[1]])]
      messag <- try(download.file(
        url = paste0("https://www.encodeproject.org", url),
        destfile = paste0(work_dir, i, "/", target, "_", suffix)
      ))
      if (class(messag) == "try-error") {message("--------> Failed to download file for ", target)}
      message("--------> Finished downloading ChIP-seq data of target ", target, "\n")
    })
  }

  for (j in 1:nrow(data_info[[i]])) {
    url <- data_info[[i]]$`Download URL`[j]
    target <- data_info[[i]]$`Target label`[j]
    bed.file <- strsplit(url, split = "/")
    suffix <- bed.file[[1]][length(bed.file[[1]])]
    if (file.exists(paste0(work_dir, i, "/", suffix))) {
      file.rename(paste0(work_dir, i, "/", suffix), paste0(work_dir, i, "/", toupper(target), "_", suffix))
    }
  }


  # ================= 7. Unzip Files =================
  gz_files = list.files(path = paste0(work_dir, i), pattern = "\\.bed\\.gz$", full.names = TRUE)
  for (f in gz_files) {R.utils::gunzip(f, overwrite = TRUE)}
  message("--------> unziped .bed.gz files in dir: ", i)
}


# ======================== 8. Save Results ========================
# write.xlsx(data_info, paste0(work_dir, "ATAC-GEX-TF-occupancy-validation.xlsx"))
message("Finished reading the Excel sheets one by one.")
# save.image(paste0(work_dir, "1-Download-TF-ChiPseq.RData"))

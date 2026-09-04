################################################################################
#                                                                              #
#   Select Hi-C data to validate the identified TFs and TF-region relations    #
#                                                                              #
################################################################################


# ======================== 1. Setup and Parameters ========================
library(dplyr)
library(readxl)
library(openxlsx)
work_dir <- "/mnt/i/scDERI_Paper/2-obtain_evaluation_criteria/2-HiC-evaluation/"
setwd(work_dir)


# ======================== 2. ENCODE Metadata Acquisition & TF Peak Selection ========================
# Retrieve metadata of TF ChIP-seq peak files from the ENCODE database (https://www.encodeproject.org/search/?type=Files);
# for each TF, retain the record with the largest file size for benchmarking
# All metadata TSV files must be manually downloaded to the current working directory and named by dataset

# ---------- 2.1 A549 Cell Line ----------
# Filtering criteria:
#   - Format: output category = annotation;
#             Output type = contact matrix;
#             File format = hic
#   - Assay: Assay title = in situ Hi-C
#   - Biosample: Biosample = A549
#   - Analysis: Genome assembly = GRCh38
#   - Provenance: Project = ENCODE
#   - Quality: Status = released
# Filter results page: https://www.encodeproject.org/search/?type=File&output_category=annotation&output_type=contact+matrix&file_format=hic&assay_title=in+situ+Hi-C&biosample_ontology.term_name=A549&assembly=GRCh38&award.project=ENCODE&status=released
# TSV metadata download link: https://www.encodeproject.org/report.tsv?type=File&output_category=annotation&output_type=contact+matrix&file_format=hic&assay_title=in+situ+Hi-C&biosample_ontology.term_name=A549&assembly=GRCh38&award.project=ENCODE&status=released
# Local saved filename: A549_hic_ENCODE_all.tsv

# ---------- 2.2 GM12878 Cell Line ----------
# Filtering criteria:
#   - Format: output category = annotation;
#             Output type = contact matrix;
#             File format = hic
#   - Assay: Assay title = in situ Hi-C
#   - Biosample: Biosample = GM2878
#   - Analysis: Genome assembly = GRCh38
#   - Provenance: Project = ENCODE
#   - Quality: Status = released
# Filter results page: https://www.encodeproject.org/report/?type=File&output_category=annotation&file_format=hic&file_type=hic&output_type=contact+matrix&biosample_ontology.term_name=GM12878&assembly=GRCh38&award.project=ENCODE&assay_title=in+situ+Hi-C&sort=-file_size
# TSV metadata download link: https://www.encodeproject.org/report.tsv?type=File&output_category=annotation&file_format=hic&file_type=hic&output_type=contact+matrix&biosample_ontology.term_name=GM12878&assembly=GRCh38&award.project=ENCODE&assay_title=in+situ+Hi-C&sort=-file_size
# Local saved filename: GM12878_hic_ENCODE_all.tsv

# ---------- 2.3 Bone Marrow Mixed Cell Lines ----------
# Filtering criteria:
#   - Format: output category = annotation;
#             Output type = contact matrix;
#             File format = hic
#   - Assay: Assay title = in situ Hi-C
#   - Analysis: Genome assembly = GRCh38
#   - Provenance: Project = ENCODE
#   - Quality: Status = released
# Filter results page: https://www.encodeproject.org/report/?type=File&output_category=annotation&file_format=hic&output_type=contact+matrix&assembly=GRCh38&award.project=ENCODE&assay_title=in+situ+Hi-C&sort=-file_size&biosample_ontology.organ_slims=bone+marrow&status=released
# TSV metadata download link: https://www.encodeproject.org/report.tsv?type=File&output_category=annotation&file_format=hic&output_type=contact+matrix&assembly=GRCh38&award.project=ENCODE&assay_title=in+situ+Hi-C&sort=-file_size&biosample_ontology.organ_slims=bone+marrow&status=released
# Local saved filename: Bone_Marrow_hic_ENCODE_all.tsv

# Read metadata by dataset; retain the record with the largest file size for each dataset
hic_ENCODE <- list()
for (dataset in c("A549", "GM12878", "Bone-marrow")) {
  hic_ENCODE[[paste0(dataset, "-hg38-HiC")]] <-
    data.table::fread(paste0("./", dataset, "_hic_ENCODE_all.tsv")) %>%
    as.data.frame() %>%
    dplyr::slice_max(`File size`, n = 1, with_ties = FALSE) %>%
    as.data.frame()
}

# View names of the output datasets
names(hic_ENCODE)
# [1] "A549-hg38-HiC"        "GM12878-hg38-HiC"     "Bone-marrow-hg38-HiC"

# Export the filtered metadata as an Excel file
writexl::write_xlsx(
  x = hic_ENCODE,
  path = "./encode-hic-enhancer-gene-benchmarking.xlsx"
)


# ======================== 3. Load Excel File Information ========================
message("\nBegan reading the Excel sheets one by one ...")
sheet_names <- readxl::excel_sheets("encode-hic-enhancer-gene-benchmarking.xlsx")
message(length(sheet_names), " sheets exist in this Excel file")


# ======================== 4. Initialize Data Storage ========================
data_info <- list()
options(timeout = 60*60*24)


# ======================== 5. Process Each Sheet ========================
for(i in sheet_names) {
  # i = sheet_names[2]

  # 5.1 Read sheet
  message("----> Reading the sheet ", i, " ...")
  data_info[[i]] = readxl::read_excel("encode-hic-enhancer-gene-benchmarking.xlsx", sheet = i)
  message("----> Found ", nrow(data_info[[i]]), " Hi-C files in sheet ", i)

  # 5.2 Create output directory
  dir.create(paste0(work_dir, i))


  # ================= 6. Download Files =================
  message("\n----> Downloading datasets from sheet ", i, " ...")
  hic.file <- strsplit(data_info[[i]]$`Download URL`, split = "/")
  suffix <- hic.file[[1]][length(hic.file[[1]])]
  messag <- try(download.file(
    url = paste0("https://www.encodeproject.org", data_info[[i]]$`Download URL`),
    destfile = paste0(work_dir, i, "/", suffix)
  ))
  if (class(messag) == "try-error") {
    message("--------> Failed to download file for ", suffix)
  }
  message("--------> Finished downloading Hi-C data ", suffix, "\n")
}


# ======================== 7. Save Results ========================
# write.xlsx(data_info, paste0(work_dir, "ATAC-GEX-chrom-interaction-validation.xlsx"))
message("Finished reading the Excel sheets one by one\n")
# save.image(paste0(work_dir, "1-Download-chrom-interact.RData"))

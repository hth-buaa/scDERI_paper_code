##############################################################
#                                                            #
#   Generate All Subdirectories for Whole-Process Analysis   #
#                                                            #
##############################################################


# --- 1. Define Root Path ---
# root_path = "I:/scDERI_Paper/"
root_path = "/mnt/i/scDERI_Paper/"

# --- 2. Define Relative Paths for All Subdirectories ---
sub_dirs = c(
  # 1-preprocessing
  "1-preprocessing",
  "1-preprocessing/10X_hg38_Human_bone_marrow_with_T_cell_depleted",
  "1-preprocessing/10X_hg38_Human_bone_marrow_with_T_cell_depleted/10X-human-bone-marrow",
  "1-preprocessing/10X_hg38_Human_bone_marrow_with_T_cell_depleted/10X-human-bone-marrow/rda",
  "1-preprocessing/10X_hg38_Human_bone_marrow_with_T_cell_depleted/10X-human-bone-marrow/sc_analysis",
  "1-preprocessing/10X_hg38_Human_bone_marrow_with_T_cell_depleted/hdf5",
  "1-preprocessing/SNARE-seq2",
  "1-preprocessing/SNARE-seq2/GSE157660-A549-GM12878",
  "1-preprocessing/SNARE-seq2/GSE157660-A549-GM12878/rda",
  "1-preprocessing/SNARE-seq2/GSE157660-A549-GM12878/human.snare2-rna.matrix",
  "1-preprocessing/SNARE-seq2/GSE157660-A549-GM12878/human.snare2-atac.matrix",

  # 2-obtain_evaluation_criteria
  "2-obtain_evaluation_criteria",
  "2-obtain_evaluation_criteria/1-TF-evaluation",
  "2-obtain_evaluation_criteria/1-TF-evaluation/A549-hg38-ChIP-seq",
  "2-obtain_evaluation_criteria/1-TF-evaluation/Bone-marrow-hg38-ChIP-seq",
  "2-obtain_evaluation_criteria/1-TF-evaluation/GM12878-hg38-ChIP-seq",
  "2-obtain_evaluation_criteria/2-HiC-evaluation",
  "2-obtain_evaluation_criteria/2-HiC-evaluation/A549-hg38-ChIP-seq",
  "2-obtain_evaluation_criteria/2-HiC-evaluation/Bone-marrow-hg38-ChIP-seq",
  "2-obtain_evaluation_criteria/2-HiC-evaluation/GM12878-hg38-ChIP-seq",

  # 3-benchmark
  "3-benchmark",
  "3-benchmark/0-downsampling",
  "3-benchmark/0-downsampling/10X_hg38_Human_bone_marrow_with_T_cell_depleted",
  "3-benchmark/0-downsampling/SNAREseq2_hg38_A549",
  "3-benchmark/0-downsampling/SNAREseq2_hg38_GM12878",
  "3-benchmark/0-downsampling/Downsampled-data",
  "3-benchmark/0-downsampling/Downsampled-data/10X_hg38_Human_bone_marrow_with_T_cell_depleted",
  "3-benchmark/0-downsampling/Downsampled-data/SNAREseq2_hg38_A549",
  "3-benchmark/0-downsampling/Downsampled-data/SNAREseq2_hg38_GM12878",
  "3-benchmark/1-STREAM",
  "3-benchmark/1-STREAM/10X_hg38_Human_bone_marrow_with_T_cell_depleted",
  "3-benchmark/1-STREAM/SNAREseq2_hg38_A549",
  "3-benchmark/1-STREAM/SNAREseq2_hg38_GM12878",
  "3-benchmark/2-IReNA",
  "3-benchmark/2-IReNA/10X_hg38_Human_bone_marrow_with_T_cell_depleted",
  "3-benchmark/2-IReNA/SNAREseq2_hg38_A549",
  "3-benchmark/2-IReNA/SNAREseq2_hg38_GM12878",
  "3-benchmark/3-FigR",
  "3-benchmark/3-FigR/10X_hg38_Human_bone_marrow_with_T_cell_depleted",
  "3-benchmark/3-FigR/SNAREseq2_hg38_A549",
  "3-benchmark/3-FigR/SNAREseq2_hg38_GM12878",
  "3-benchmark/4-scGLUE",
  "3-benchmark/4-scGLUE/10X_hg38_Human_bone_marrow_with_T_cell_depleted",
  "3-benchmark/4-scGLUE/SNAREseq2_hg38_A549",
  "3-benchmark/4-scGLUE/SNAREseq2_hg38_GM12878",
  "3-benchmark/4-scGLUE/ChIPseq-data",
  "3-benchmark/4-scGLUE/gtf-files",
  "3-benchmark/5-DIRECT-NET",
  "3-benchmark/5-DIRECT-NET/10X_hg38_Human_bone_marrow_with_T_cell_depleted",
  "3-benchmark/5-DIRECT-NET/SNAREseq2_hg38_A549",
  "3-benchmark/5-DIRECT-NET/SNAREseq2_hg38_GM12878",
  "3-benchmark/6-Pando",
  "3-benchmark/6-Pando/10X_hg38_Human_bone_marrow_with_T_cell_depleted",
  "3-benchmark/6-Pando/SNAREseq2_hg38_A549",
  "3-benchmark/6-Pando/SNAREseq2_hg38_GM12878",
  "3-benchmark/7-scMEGA",
  "3-benchmark/7-scMEGA/10X_hg38_Human_bone_marrow_with_T_cell_depleted",
  "3-benchmark/7-scMEGA/SNAREseq2_hg38_A549",
  "3-benchmark/7-scMEGA/SNAREseq2_hg38_GM12878",
  "3-benchmark/8-scDERI",
  "3-benchmark/8-scDERI/10X_hg38_Human_bone_marrow_with_T_cell_depleted",
  "3-benchmark/8-scDERI/SNAREseq2_hg38_A549",
  "3-benchmark/8-scDERI/SNAREseq2_hg38_GM12878",
  "3-benchmark/9-evaluation",
  "3-benchmark/9-evaluation/1-TF_coverage",
  "3-benchmark/9-evaluation/2-TF-peak_relationship",
  "3-benchmark/9-evaluation/3-peak-TG_relationship",
  "3-benchmark/9-evaluation/4-visualization",

  # 4-apply
  "4-apply",
  "4-apply/raw_data",
  "4-apply/survival",
  "4-apply/survival/TF",
  "4-apply/survival/TG",
)

# --- Concatenate root path with subdirectories to generate full path vector ---
full_paths = file.path(root_path, sub_dirs)

# --- Recursively Create Directories ---
lapply(full_paths, dir.create, recursive = TRUE, showWarnings = FALSE)

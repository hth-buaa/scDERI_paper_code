##########################################################################
#                                                                        #
#   Generate fragments file of human GM12878 and A549 single-cell data   #
#                                                                        #
##########################################################################

data.dir = "/mnt/i/scDERI_Paper/1-preprocessing/SNARE-seq2/GSE157660-A549-GM12878/rda/"
setwd(data.dir)


# 1. Load required libraries
library(qs)
library(tidyr)
library(stringr)
library(dplyr)


# 2. Define dataset list
dataset = c(
  "SNAREseq2_hg38_GM12878",
  "SNAREseq2_hg38_A549"
)


# 3. Process each dataset
for (i in dataset) {
  # i = dataset[1] # for test

  # 3.1 Load preprocessed object
  obj = qread(paste0("./", i, ".qsave"))

  max(obj@assays[["ATAC"]]@counts)
  # [1] 4

  dim(obj@assays[["ATAC"]]@counts)
  # [1] 143101   8876

  # [1] 143101   1088

  class(obj@assays[["ATAC"]]@counts)
  # [1] "dgCMatrix"
  # attr(,"package")
  # [1] "Matrix"

  Matrix::nnzero(obj@assays[["ATAC"]]@counts)
  # [1] 14283243

  # [1] 1103693

  rownames(obj@assays[["ATAC"]]@counts)[1:3]
  # [1] "chr1-804776-805277" "chr1-826804-827305" "chr1-827350-827851"

  colnames(obj@assays[["ATAC"]]@counts)[1:3]
  # [1] "GM12878_cells_1pctPFA_4C_exp1_P6_ACAAGCTAAACGTGATGAGCTGAA"  "GM12878_cells_1pctPFA_4C_exp1_P9_ACACAGAAAGCAGGAAGAGCTGAA"
  # [3] "GM12878_cells_1pctPFA_4C_exp1_P10_GTGTTCTACGAACTTAGAGCTGAA"

  # [1] "AACAACCAGCCAAGACACAAGCTA" "ACAGCAGAACGCTCGAACAAGCTA" "AAGAGATCGCCAAGACAGTACAAG"

  # 3.2 Extract ATAC counts matrix
  ATACcounts = obj@assays[["ATAC"]]@counts
  ATACcounts[1:3,1:3]
  # 3 x 3 sparse Matrix of class "dgCMatrix"
  #                    GM12878_cells_1pctPFA_4C_exp1_P6_ACAAGCTAAACGTGATGAGCTGAA GM12878_cells_1pctPFA_4C_exp1_P9_ACACAGAAAGCAGGAAGAGCTGAA
  # chr1-804776-805277                                                         .                                                         .
  # chr1-826804-827305                                                         .                                                         .
  # chr1-827350-827851                                                         .                                                         .
  #                    GM12878_cells_1pctPFA_4C_exp1_P10_GTGTTCTACGAACTTAGAGCTGAA
  # chr1-804776-805277                                                          .
  # chr1-826804-827305                                                          .
  # chr1-827350-827851                                                          .

  # 3 x 3 sparse Matrix of class "dgCMatrix"
  #                    AACAACCAGCCAAGACACAAGCTA ACAGCAGAACGCTCGAACAAGCTA AAGAGATCGCCAAGACAGTACAAG
  # chr1-804776-805277                        .                        .                        .
  # chr1-826804-827305                        .                        .                        .
  # chr1-827350-827851                        .                        1                        .

  # 3.3 Convert sparse matrix to data frame
  sparse_triplets = summary(ATACcounts)
  result_df = data.frame(
    peak = rownames(ATACcounts)[sparse_triplets$i],
    cell = colnames(ATACcounts)[sparse_triplets$j],
    value = sparse_triplets$x
  )
  head(result_df)
  #                   peak                                                      cell value
  # 1 chr1-1064063-1064564 GM12878_cells_1pctPFA_4C_exp1_P6_ACAAGCTAAACGTGATGAGCTGAA     1
  # 2 chr1-1247646-1248147 GM12878_cells_1pctPFA_4C_exp1_P6_ACAAGCTAAACGTGATGAGCTGAA     1
  # 3 chr1-1308841-1309342 GM12878_cells_1pctPFA_4C_exp1_P6_ACAAGCTAAACGTGATGAGCTGAA     2
  # 4 chr1-1691909-1692410 GM12878_cells_1pctPFA_4C_exp1_P6_ACAAGCTAAACGTGATGAGCTGAA     1
  # 5 chr1-1724210-1724711 GM12878_cells_1pctPFA_4C_exp1_P6_ACAAGCTAAACGTGATGAGCTGAA     2
  # 6 chr1-2412425-2412926 GM12878_cells_1pctPFA_4C_exp1_P6_ACAAGCTAAACGTGATGAGCTGAA     2

  #                   peak                     cell value
  # 1 chr1-1693286-1693787 AACAACCAGCCAAGACACAAGCTA     1
  # 2 chr1-2133094-2133595 AACAACCAGCCAAGACACAAGCTA     2
  # 3 chr1-3543470-3543971 AACAACCAGCCAAGACACAAGCTA     1
  # 4 chr1-6785705-6786206 AACAACCAGCCAAGACACAAGCTA     1
  # 5 chr1-7893821-7894322 AACAACCAGCCAAGACACAAGCTA     2
  # 6 chr1-8876735-8877236 AACAACCAGCCAAGACACAAGCTA     1

  # 3.4 Process and clean data
  result_df1 = separate(result_df, peak, into = c("chrom", "start", "end"), sep = "-")
  head(result_df1)
  #   chrom   start     end                                                      cell value
  # 1  chr1 1064063 1064564 GM12878_cells_1pctPFA_4C_exp1_P6_ACAAGCTAAACGTGATGAGCTGAA     1
  # 2  chr1 1247646 1248147 GM12878_cells_1pctPFA_4C_exp1_P6_ACAAGCTAAACGTGATGAGCTGAA     1
  # 3  chr1 1308841 1309342 GM12878_cells_1pctPFA_4C_exp1_P6_ACAAGCTAAACGTGATGAGCTGAA     2
  # 4  chr1 1691909 1692410 GM12878_cells_1pctPFA_4C_exp1_P6_ACAAGCTAAACGTGATGAGCTGAA     1
  # 5  chr1 1724210 1724711 GM12878_cells_1pctPFA_4C_exp1_P6_ACAAGCTAAACGTGATGAGCTGAA     2
  # 6  chr1 2412425 2412926 GM12878_cells_1pctPFA_4C_exp1_P6_ACAAGCTAAACGTGATGAGCTGAA     2

  #   chrom   start     end                     cell value
  # 1  chr1 1693286 1693787 AACAACCAGCCAAGACACAAGCTA     1
  # 2  chr1 2133094 2133595 AACAACCAGCCAAGACACAAGCTA     2
  # 3  chr1 3543470 3543971 AACAACCAGCCAAGACACAAGCTA     1
  # 4  chr1 6785705 6786206 AACAACCAGCCAAGACACAAGCTA     1
  # 5  chr1 7893821 7894322 AACAACCAGCCAAGACACAAGCTA     2
  # 6  chr1 8876735 8877236 AACAACCAGCCAAGACACAAGCTA     1

  result_df2 = arrange(result_df1, chrom, start, end)
  head(result_df2)
  #   chrom   start     end                                                       cell value
  # 1  chr1 1000056 1000557 GM12878_cells_1pctPFA_4C_exp1_P12_GATAGACAAACAACCACAAGACTA     2
  # 2  chr1 1000056 1000557 GM12878_cells_1pctPFA_4C_exp1_P12_AGCCATGCGAGCTGAACAAGACTA     2
  # 3  chr1 1000056 1000557  GM12878_cells_1pctPFA_4C_exp1_P5_AAGGTACATCTTCACAGAGCTGAA     2
  # 4  chr1 1000056 1000557  GM12878_cells_1pctPFA_4C_exp1_P5_AGCAGGAACGACACACGAGCTGAA     2
  # 5  chr1 1000056 1000557  GM12878_cells_1pctPFA_4C_exp1_P5_CGACTGGATATCAGCAGAGCTGAA     1
  # 6  chr1 1000056 1000557  GM12878_cells_1pctPFA_4C_exp1_P5_CGACTGGAATTGAGGACAAGACTA     1

  #   chrom   start     end                     cell value
  # 1  chr1 1000056 1000557 GCTCGGTACTAAGGTCAGTACAAG     1
  # 2  chr1 1000056 1000557 ACACAGAACATACCAAAGTACAAG     2
  # 3  chr1 1000056 1000557 CGGATTGCACATTGGCACAAGCTA     2
  # 4  chr1 1000056 1000557 CATACCAAAACAACCAAGTACAAG     1
  # 5  chr1 1000056 1000557 ATAGCGACCTAAGGTCCTGTAGCC     2
  # 6  chr1 1000056 1000557 ACACAGAACCATCCTCCTGTAGCC     4

  # 3.5 Save fragments file
  write.table(
    x = result_df2,
    file = paste0("./", i, "_fragments.tsv"),
    sep = "\t",
    col.names = FALSE,
    row.names = FALSE,
    quote = FALSE
  )
}


# 4 Process and index each dataset
# 4.1 Process and index GM12878 dataset
# system("bgzip /mnt/i/scDERI_Paper/1-preprocessing/SNARE-seq2/GSE157660-A549-GM12878/rda/SNAREseq2_hg38_GM12878_fragments.tsv")
# system("gunzip -c /mnt/i/scDERI_Paper/1-preprocessing/SNARE-seq2/GSE157660-A549-GM12878/rda/SNAREseq2_hg38_GM12878_fragments.tsv.gz | sort -k1,1 -k2,2n | bgzip > /mnt/i/scDERI_Paper/1-preprocessing/SNARE-seq2/GSE157660-A549-GM12878/rda/SNAREseq2_hg38_GM12878_fragments.sorted.tsv.gz")
# system("tabix -p bed /mnt/i/scDERI_Paper/1-preprocessing/SNARE-seq2/GSE157660-A549-GM12878/rda/SNAREseq2_hg38_GM12878_fragments.sorted.tsv.gz")

system("sort -k1,1 -k2,2n /mnt/i/scDERI_Paper/1-preprocessing/SNARE-seq2/GSE157660-A549-GM12878/rda/SNAREseq2_hg38_GM12878_fragments.tsv | bgzip > /mnt/i/scDERI_Paper/1-preprocessing/SNARE-seq2/GSE157660-A549-GM12878/rda/SNAREseq2_hg38_GM12878_fragments.tsv.gz")
system("tabix -p bed /mnt/i/scDERI_Paper/1-preprocessing/SNARE-seq2/GSE157660-A549-GM12878/rda/SNAREseq2_hg38_GM12878_fragments.tsv.gz")

# 4.2. Process and index A549 dataset
# system("bgzip /mnt/i/scDERI_Paper/1-preprocessing/SNARE-seq2/GSE157660-A549-GM12878/rda/SNAREseq2_hg38_A549_fragments.tsv")
# system("gunzip -c /mnt/i/scDERI_Paper/1-preprocessing/SNARE-seq2/GSE157660-A549-GM12878/rda/SNAREseq2_hg38_A549_fragments.tsv.gz | sort -k1,1 -k2,2n | bgzip > /mnt/i/scDERI_Paper/1-preprocessing/SNARE-seq2/GSE157660-A549-GM12878/rda/SNAREseq2_hg38_A549_fragments.sorted.tsv.gz")
# system("tabix -p bed /mnt/i/scDERI_Paper/1-preprocessing/SNARE-seq2/GSE157660-A549-GM12878/rda/SNAREseq2_hg38_A549_fragments.sorted.tsv.gz")

system("sort -k1,1 -k2,2n /mnt/i/scDERI_Paper/1-preprocessing/SNARE-seq2/GSE157660-A549-GM12878/rda/SNAREseq2_hg38_A549_fragments.tsv | bgzip > /mnt/i/scDERI_Paper/1-preprocessing/SNARE-seq2/GSE157660-A549-GM12878/rda/SNAREseq2_hg38_A549_fragments.tsv.gz")
system("tabix -p bed /mnt/i/scDERI_Paper/1-preprocessing/SNARE-seq2/GSE157660-A549-GM12878/rda/SNAREseq2_hg38_A549_fragments.tsv.gz")

# save.image("./2-Preprocess-SNARE-seq2-human-GM12878-A549.RData")

# Up to this point, we have obtained the files _fragments.sorted.tsv.gz and _fragments.sorted.tsv.gz.tbi, which can be used with Signac.

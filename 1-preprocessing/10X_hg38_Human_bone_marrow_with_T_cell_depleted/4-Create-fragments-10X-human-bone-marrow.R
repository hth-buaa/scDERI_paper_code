#####################################################################
#                                                                   #
#    Create fragments file of human bone marrow single cell data    #
#                                                                   #
#####################################################################


# 1. Load libraries
library(data.table)
library(qs)
library(dplyr)
data.dir = "/mnt/i/scDERI_Paper/1-preprocessing/10X_hg38_Human_bone_marrow_with_T_cell_depleted/10X-human-bone-marrow/rda/"
setwd(data.dir)


# 2. Read fragment files
frag1 = fread("../BM_Tcelldep_Rep1_atac_fragments.tsv.gz")
frag2 = fread("../BM_Tcelldep_Rep2_atac_fragments.tsv.gz")


# 3. Standardize column names
colnames(frag1) = c("chrom", "start", "end", "barcode", "count")
colnames(frag2) = c("chrom", "start", "end", "barcode", "count")
head(frag1)
#     chrom start   end            barcode count
#    <char> <int> <int>             <char> <int>
# 1:   chr1 10071 10204 ACTAACTCATTGTGTG-1     1
# 2:   chr1 10072 10339 GCACCTAAGGCATGAA-1     1
# 3:   chr1 10073 10309 CATTCATTCTAAGTGC-1     1
# 4:   chr1 10078 10279 TTTCAGTTCCGCAAAT-1     1
# 5:   chr1 10079 10321 AAGGTATAGGTCGAGG-1     1
# 6:   chr1 10079 10326 ACCAAACTCAAGACTC-1     1
head(frag2)
#     chrom start   end            barcode count
#    <char> <int> <int>             <char> <int>
# 1:   chr1 10001 10346 CTAACCCTCGCAACAT-1     1
# 2:   chr1 10065 10419 TGGTAAACATAGTCAT-1     1
# 3:   chr1 10073 10203 ACCTTCCAGCGTGCGT-1     1
# 4:   chr1 10073 10210 GAGTCAAAGGTCGATT-1     1
# 5:   chr1 10073 10321 AGACCCGGTGACATGC-1     1
# 6:   chr1 10077 10273 GCCTACTTCCAAACAC-1     2


# 4. Load preprocessed object
obj = qs::qread("./10X_hg38_Human_bone_marrow_with_T_cell_depleted.qsave")


# 5. Extract and intersect cell barcodes
cells1 = unique(frag1$barcode)
head(cells1)
# [1] "ACTAACTCATTGTGTG-1" "GCACCTAAGGCATGAA-1" "CATTCATTCTAAGTGC-1" "TTTCAGTTCCGCAAAT-1" "AAGGTATAGGTCGAGG-1" "ACCAAACTCAAGACTC-1"

cells2 = unique(frag2$barcode)
head(cells2)
# [1] "CTAACCCTCGCAACAT-1" "TGGTAAACATAGTCAT-1" "ACCTTCCAGCGTGCGT-1" "GAGTCAAAGGTCGATT-1" "AGACCCGGTGACATGC-1" "GCCTACTTCCAAACAC-1"

cells = intersect(cells1, cells2)


# 6. Get cell IDs from RNA and ATAC assays
cs1 = colnames(obj@assays$RNA@counts)
head(cs1)
# [1] "AAACATGCACATAGCC-1" "AAACATGCACCTATAG-1" "AAACATGCATTGCGGT-1" "AAACCAACACAGCCTG-1" "AAACCGAAGACAAAGT-1" "AAACCGAAGATGCCTG-1"

cs2 = colnames(obj@assays$ATAC@counts)
head(cs2)
# [1] "AAACATGCACATAGCC-1" "AAACATGCACCTATAG-1" "AAACATGCATTGCGGT-1" "AAACCAACACAGCCTG-1" "AAACCGAAGACAAAGT-1" "AAACCGAAGATGCCTG-1"
# cs0 = intersect(cs1, cs2)

cs = unique(c(cs1, cs2))
head(cs); rm(cs1, cs2)
# [1] "AAACATGCACATAGCC-1" "AAACATGCACCTATAG-1" "AAACATGCATTGCGGT-1" "AAACCAACACAGCCTG-1" "AAACCGAAGACAAAGT-1" "AAACCGAAGATGCCTG-1"

table(cs %in% cells1)
# FALSE  TRUE
#    48  7370
table(cs %in% cells2)
# TRUE
# 7418
table(cs %in% cells)
# FALSE  TRUE
#    48  7370


# 7. Filter fragments by shared cells
frag11 = frag1[frag1$barcode %in% cs, ]
head(frag11)
#     chrom start   end            barcode count
#    <char> <int> <int>             <char> <int>
# 1:   chr1 10072 10339 GCACCTAAGGCATGAA-1     1
# 2:   chr1 10078 10279 TTTCAGTTCCGCAAAT-1     1
# 3:   chr1 10079 10326 ACCAAACTCAAGACTC-1     1
# 4:   chr1 10084 10334 GTTCTCATCGTAATCA-1     1
# 5:   chr1 10084 10339 GATTACGGTGGATTGC-1     3
# 6:   chr1 10085 10185 CTTGTAAAGGTCCGTA-1     1

frag21 = frag2[frag2$barcode %in% cs, ]
head(frag21)
#     chrom start   end            barcode count
#    <char> <int> <int>             <char> <int>
# 1:   chr1 10001 10346 CTAACCCTCGCAACAT-1     1
# 2:   chr1 10065 10419 TGGTAAACATAGTCAT-1     1
# 3:   chr1 10073 10203 ACCTTCCAGCGTGCGT-1     1
# 4:   chr1 10073 10210 GAGTCAAAGGTCGATT-1     1
# 5:   chr1 10073 10321 AGACCCGGTGACATGC-1     1
# 6:   chr1 10077 10273 GCCTACTTCCAAACAC-1     2


# 8. Merge and sum counts
frag11_renamed = frag11 %>% rename(count_1 = count)
frag21_renamed = frag21 %>% rename(count_2 = count)
merged_df = full_join(frag11_renamed, frag21_renamed, by = c("chrom", "start", "end", "barcode"))
head(merged_df)
#     chrom start   end            barcode count_1 count_2
#    <char> <int> <int>             <char>   <int>   <int>
# 1:   chr1 10072 10339 GCACCTAAGGCATGAA-1       1      NA
# 2:   chr1 10078 10279 TTTCAGTTCCGCAAAT-1       1      NA
# 3:   chr1 10079 10326 ACCAAACTCAAGACTC-1       1      NA
# 4:   chr1 10084 10334 GTTCTCATCGTAATCA-1       1      NA
# 5:   chr1 10084 10339 GATTACGGTGGATTGC-1       3      NA
# 6:   chr1 10085 10185 CTTGTAAAGGTCCGTA-1       1      NA

merged_df1 = merged_df %>%
  mutate(count = case_when(
    !is.na(count_1) & !is.na(count_2) ~ count_1 + count_2,
    !is.na(count_1) & is.na(count_2) ~ count_1,
    is.na(count_1) & !is.na(count_2) ~ count_2
  )) %>% select(-count_1, -count_2)
head(merged_df1); rm(merged_df)
#     chrom start   end            barcode count
#    <char> <int> <int>             <char> <int>
# 1:   chr1 10072 10339 GCACCTAAGGCATGAA-1     1
# 2:   chr1 10078 10279 TTTCAGTTCCGCAAAT-1     1
# 3:   chr1 10079 10326 ACCAAACTCAAGACTC-1     1
# 4:   chr1 10084 10334 GTTCTCATCGTAATCA-1     1
# 5:   chr1 10084 10339 GATTACGGTGGATTGC-1     3
# 6:   chr1 10085 10185 CTTGTAAAGGTCCGTA-1     1


# 9. Sort and filter
merged_df2 = merged_df1[order(chrom, start, end)]
merged_df2 = merged_df2[grepl("chr", merged_df2$chrom), ]
head(merged_df2); rm(merged_df1)
#     chrom start   end            barcode count
#    <char> <int> <int>             <char> <int>
# 1:   chr1 10001 10346 CTAACCCTCGCAACAT-1     1
# 2:   chr1 10065 10419 TGGTAAACATAGTCAT-1     1
# 3:   chr1 10072 10339 GCACCTAAGGCATGAA-1     1
# 4:   chr1 10073 10203 ACCTTCCAGCGTGCGT-1     1
# 5:   chr1 10073 10210 GAGTCAAAGGTCGATT-1     1
# 6:   chr1 10073 10321 AGACCCGGTGACATGC-1     1


# 10. Save and compress
data.table::fwrite(
  merged_df2,
  file = "./10X_hg38_Human_bone_marrow_with_T_cell_depleted_fragments.tsv",
  sep = "\t",
  col.names = FALSE
)

system("bgzip        /mnt/i/scDERI_Paper/1-preprocessing/10X_hg38_Human_bone_marrow_with_T_cell_depleted/10X-human-bone-marrow/rda/10X_hg38_Human_bone_marrow_with_T_cell_depleted_fragments.tsv")
system("tabix -p bed /mnt/i/scDERI_Paper/1-preprocessing/10X_hg38_Human_bone_marrow_with_T_cell_depleted/10X-human-bone-marrow/rda/10X_hg38_Human_bone_marrow_with_T_cell_depleted_fragments.tsv.gz")

# save.image("./4-Preprocess-10X-human-bone-marrow.Rdata")

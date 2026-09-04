###############################################################################
#                                                                             #
#   Evaluate the peak-TG relationship prediction performance of each method   #
#                                                                             #
###############################################################################

# Evaluation steps:
# 1. For each method, extract predicted peak-TG relationship
# 2. Compare with ground truth peak-peak relationship from Hi-C (for each dataset, calculate the proportion of region-to-gene linkages against Hi-C contact pairs)
# 3. Calculate precision, recall, and F1-score for each method


##### 1. Environment Setup and Parameter Configuration #####
# Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
# rm(list = ls())


##### 2. Path and Package Configuration #####
base.dir = "/mnt/i/scDERI_Paper/3-benchmark"
parent.dir = "/mnt/i/scDERI_Paper/3-benchmark/9-evaluation/3-peak-TG_relationship/"
eval.dir = "/mnt/i/scDERI_Paper/2-obtain_evaluation_criteria/2-HiC-evaluation/"

# Load required packages
library(data.table); library(dplyr); library(Signac); library(Pando); library(strawr)

# Define dataset names
data.ls = c(
  "10X_hg38_Human_bone_marrow_with_T_cell_depleted",
  "SNAREseq2_hg38_GM12878",
  "SNAREseq2_hg38_A549"
)
message("There are ", length(data.ls), " datasets")


##### 3. Define Helper Functions #####
# 3.1 Collapse transcripts to longest transcript per gene
# This function collapses multiple transcripts of a gene to a single "longest transcript" representation, taking the full genomic span from the earliest start to the latest end of all transcripts
CollapseToLongestTranscript = function(
    ranges # GRanges object containing transcript annotations for genes
) {
  range.df = as.data.table(x = ranges)
  range.df$strand = as.character(x = range.df$strand)
  range.df$strand = ifelse(test = range.df$strand == "*", yes = "+", no = range.df$strand)
  collapsed = range.df[, .(unique(seqnames), min(start), max(end), strand[[1]], gene_biotype[[1]], gene_name[[1]]), "gene_id"]
  colnames(x = collapsed) = c("gene_id", "seqnames", "start", "end", "strand", "gene_biotype", "gene_name")
  collapsed$gene_name = make.unique(names = collapsed$gene_name)
  gene.ranges = makeGRangesFromDataFrame(df = collapsed, keep.extra.columns = TRUE)
  return(gene.ranges)
}

# 3.2 Preload gene annotations for different genomes
Annot = list()
# library(EnsDb.Mmusculus.v79)
# Annot[["mm10"]] = GetGRangesFromEnsDb(EnsDb.Mmusculus.v79) %>% CollapseToLongestTranscript
# library(EnsDb.Hsapiens.v75)
# Annot[["hg19"]] = GetGRangesFromEnsDb(EnsDb.Hsapiens.v75) %>% CollapseToLongestTranscript
library(EnsDb.Hsapiens.v86)
Annot[["hg38"]] = GetGRangesFromEnsDb(EnsDb.Hsapiens.v86) %>% CollapseToLongestTranscript

# 3.3 Function to get genomic locations for gene symbols
get_gene_loc = function(
    symbols, # Vector of gene symbols to get genomic locations for
    org = "hg38" # Genome assembly: "hg38", "hg19", or "mm10"
) {
  # Get annotation for the specified genome assembly
  annot = Annot[[org]]
  message("Loaded genomic loci for genes in ", org, " assembly")

  # Create named vector: gene_name -> genomic_location
  named.loci = setNames(paste0("chr", GRangesToString(annot)), annot$gene_name)

  # Extract locations for input symbols
  subset.loci = named.loci[symbols]
  names(subset.loci) = NULL
  subset.loci
}


##### 4. Extract Predicted Region-Gene Relationships from Each Method #####
# Goal: Extract region-gene (peak-target gene) relationships predicted by each method
# Each method's results are stored in different formats and locations

# 4.1 scDERI region-gene extraction
message("Extracting region-gene relationships from scDERI...")
work.dir = paste0(base.dir, "/8-scDERI/"); setwd(work.dir)
scDERI_region_TG_links = pbapply::pblapply(seq_along(data.ls), function(i) {
  message("----> Conversion for the ", i, "-th dataset: ", data.ls[i])
  if (!file.exists(paste0(work.dir, data.ls[i], "/interest_cell_type_peGRN.qsave"))) {
    message("----> No eRegulon file on dataset: ", data.ls[i])
    return(NULL)
  }
  res = qs::qread(paste0(work.dir, data.ls[i], "/interest_cell_type_peGRN.qsave"))[[1]][["loose_filtered_simplified"]][, c("peak", "TG")] %>% unique
  res1 = do.call(rbind, strsplit(res$peak, "-"))
  res2 = data.frame(
    seqnames = factor(res1[, 1]),
    start = as.numeric(res1[, 2]),
    end = as.numeric(res1[, 3]),
    width = as.numeric(res1[, 3]) - as.numeric(res1[, 2]) + 1,
    strand = factor("*", levels = c("*", "+", "-")),
    # gene = res[, "SYMBOL"],
    gene = res[, "TG"],
    row.names = NULL
  )
  dt = data.frame(
    Region = paste0(res2$seqnames, "-", res2$start, "-", res2$end),
    Gene = res2$gene,
    Gene.loc = get_gene_loc(symbols = res2$gene, org = strsplit(data.ls[i], split = "_")[[1]][2])
  )
  return(dt[!duplicated(dt) & !is.na(dt$Gene.loc), , drop = FALSE])
})
names(scDERI_region_TG_links) = data.ls
qs::qsave(scDERI_region_TG_links, paste0(parent.dir, "scDERI_region_TG_links.qsave"))
# scDERI_region_TG_links = qs::qread(paste0(parent.dir, "scDERI_region_TG_links.qsave"))

# 4.2 STREAM region-gene extraction
message("Extracting region-gene relationships from STREAM...")
work.dir = paste0(base.dir, "/1-STREAM/"); setwd(work.dir)
STREAM_region_TG_links = pbapply::pblapply(seq_along(data.ls), function(i) {
  message("----> Conversion for the ", i, "-th dataset: ", data.ls[i])
  if (!file.exists(paste0(work.dir, data.ls[i], "/eRegulons.qsave"))) {
    message("----> No eRegulon file on dataset: ", data.ls[i])
    return(NULL)
  }
  res = qs::qread(paste0(work.dir, data.ls[i], "/eRegulons.qsave")) %>%
    lapply(., "[[", "links") %>% lapply(., as.data.frame) %>% Reduce("rbind", .)
  dt = data.frame(
    Region = paste0(res$seqnames, "-", res$start, "-", res$end),
    Gene = res$gene,
    Gene.loc = get_gene_loc(symbols = res$gene, org = strsplit(data.ls[i], split = "_")[[1]][2])
  )
  return(dt[!duplicated(dt) & !is.na(dt$Gene.loc), , drop = FALSE])
})
names(STREAM_region_TG_links) = data.ls
qs::qsave(STREAM_region_TG_links, paste0(parent.dir, "STREAM_region_TG_links.qsave"))
# STREAM_region_TG_links = qs::qread(paste0(parent.dir, "STREAM_region_TG_links.qsave"))

# 4.3 Pando region-gene extraction
message("Extracting region-gene relationships from Pando...")
work.dir = paste0(base.dir, "/6-Pando/"); setwd(work.dir)
Pando_region_TG_links = pbapply::pblapply(seq_along(data.ls), function(i) {
  message("----> Conversion for the ", i, "-th dataset: ", data.ls[i])
  if (!file.exists(paste0(work.dir, data.ls[i], ".qs"))) {
    message("----> No eRegulon file on dataset: ", data.ls[i])
    return(NULL)
  }
  grn = qs::qread(paste0(work.dir, data.ls[i], ".qs")) %>% NetworkModules
  dt = data.frame(
    Region = grn@meta$regions,
    Gene = grn@meta$target,
    Gene.loc = get_gene_loc(symbols = grn@meta$target, org = strsplit(data.ls[i], split = "_")[[1]][2])
  )
  return(dt[!duplicated(dt) & !is.na(dt$Gene.loc), , drop = FALSE])
  message("----> Wrote eRegulons to file: ", "eRegulons_", data.ls[i], ".qsave")
})
names(Pando_region_TG_links) = data.ls
Pando_region_TG_links = Pando_region_TG_links[!sapply(Pando_region_TG_links, is.null)]
qs::qsave(Pando_region_TG_links, paste0(parent.dir, "Pando_region_TG_links.qsave"))
# Pando_region_TG_links = qs::qread(paste0(parent.dir, "Pando_region_TG_links.qsave"))

# 4.4 GLUE region-gene extraction
message("Extracting region-gene relationships from GLUE...")
work.dir = paste0(base.dir, "/4-scGLUE/"); setwd(work.dir)
scGLUE_region_TG_links = pbapply::pblapply(seq_along(data.ls), function(i) {
  message("----> Conversion for the ", i, "-th dataset: ", data.ls[i])
  if (!file.exists(paste0(work.dir, data.ls[i], "/draft_grn.csv"))) {
    message("----> No eRegulon file on dataset: ", data.ls[i])
    return(NULL)
  }
  grn = read.csv(paste0(work.dir, data.ls[i], "/draft_grn.csv"), sep = ",")
  gene2peak = read.csv(paste0(work.dir, data.ls[i], "/gene2peak.csv"), sep = ",")
  colnames(gene2peak) = c("Gene", "Region", "Weight")
  colnames(grn) = c("TF", "Gene", "importance")
  res = gene2peak[gene2peak$Gene %in% grn$Gene,, drop = FALSE]
  dt = data.frame(
    Region = res$Region,
    Gene = res$Gene,
    Gene.loc = get_gene_loc(symbols = res$Gene, org = strsplit(data.ls[i], split = "_")[[1]][2])
  )
  return(dt[!duplicated(dt) | is.na(dt$Gene.loc), , drop = FALSE])
})
names(scGLUE_region_TG_links) = data.ls
scGLUE_region_TG_links = scGLUE_region_TG_links[!sapply(scGLUE_region_TG_links, is.null)]
qs::qsave(scGLUE_region_TG_links, paste0(parent.dir, "scGLUE_region_TG_links.qsave"))
# scGLUE_region_TG_links = qs::qread(paste0(parent.dir, "scGLUE_region_TG_links.qsave"))

# 4.5 Compile all method results
method.list = list(
  scDERI = scDERI_region_TG_links,
  STREAM = STREAM_region_TG_links,
  scGLUE = scGLUE_region_TG_links,
  Pando = Pando_region_TG_links
)


##### 5. Load Ground Truth from Hi-C Data #####
# 5.1 Dataset correspondence table
chipseq.hic = tibble::tibble(
  Dataset = c("10X_hg38_Human_bone_marrow_with_T_cell_depleted", "SNAREseq2_hg38_A549", "SNAREseq2_hg38_GM12878"),
  `ChIP.seq` = c("Bone-marrow-hg38-ChIP-seq", "A549-hg38-ChIP-seq", "GM12878-hg38-ChIP-seq"),
  `Hi.C` = c("Bone-marrow-hg38-HiC", "A549-hg38-HiC", "GM12878-hg38-HiC")
)

# 5.2 Function to read Hi-C contact matrices
read_hic = function(
    fname, # Path to .hic file
    org = "hg38", # Genome assembly
    unit = "BP", # Unit for binning
    binsize = 2500000, # Bin size in base pairs
    top.peaks = 100 # Number of top interactions to extract per chromosome
) {
  # Determine chromosome list based on organism
  if (grepl("^mm", org)) {chr.ll = c(1:40, "X", "Y")} else {chr.ll = c(1:23, "X", "Y")}

  # Test chromosome name format
  messag = try(straw(norm = "NONE", fname = fname, chr1loc = "1", chr2loc = "1", unit = unit, binsize = binsize))
  if (class(messag) == "try-error") {
    message("--------> Chromosome names lack a prefix 'chr'")
    message("--------> Adding 'chr' prefix to chromosome names")
    chr.ll = paste0("chr", chr.ll)
  }

  # Extract top interactions for each chromosome
  df = pbapply::pblapply(chr.ll, function(x) {
    message("Processing chromosome: ", x)
    chr.hic = try(straw(norm = "NONE", fname = fname, x, x, unit = unit, binsize = binsize) %>% arrange(desc(counts)) %>% head(n = top.peaks))
    if (class(chr.hic) == "try-error") {return(data.frame(peak1 = NA, peak2 = NA, counts = NA))}
    lapply(1:nrow(chr.hic), function(i) {
      y = chr.hic[i,]
      data.frame(
        peak1 = paste0(x, "-", format(y$x + 1, scientific = FALSE), "-", format(y$x + binsize, scientific = FALSE)),
        peak2 = paste0(x, "-", format(y$y + 1, scientific = FALSE), "-", format(y$y + binsize, scientific = FALSE)),
        counts = y$counts
      )
    }) %>% Reduce("rbind", .)
  }) %>% na.omit %>% Reduce("rbind", .)
  df = df[!is.na(df$counts),, drop = FALSE]

  # Add 'chr' prefix if needed
  if (!class(messag) == "try-error") {
    df$peak1 = paste0("chr", df$peak1)
    df$peak2 = paste0("chr", df$peak2)
  }
  return(df)
}

# 5.3 Load Hi-C contact matrices
message("Loading Hi-C contact matrices...")
top.peaks = 100
HiC_contact_matrix_list = pbmcapply::pbmclapply(1:nrow(chipseq.hic), mc.cores = parallel::detectCores() - 1, function(i) {
  message("----> Loading Hi-C data for dataset: ", chipseq.hic$Dataset[i])
  if (chipseq.hic$Hi.C[i] == "N/A") { return(NULL) }
  dir.ll = strsplit(chipseq.hic$Hi.C[i], split = ", ") %>% `[[` (1)
  pbapply::pblapply(dir.ll, function(x) {
    message("Processing directory: ", x)
    org = strsplit(x, split = "-")[[1]][length(strsplit(x, split = "-")[[1]])-1]
    hic.files = list.files(paste0(eval.dir, x), pattern = ".hic")
    dir_region_TG_links = lapply(hic.files, function(y) {
      read_hic(fname = paste0(eval.dir, x, "/", y), org = org, top.peaks = top.peaks)
    }) %>% Reduce("rbind", .)
  }) %>% Reduce("rbind", .)
})
names(HiC_contact_matrix_list) = chipseq.hic$Dataset
qs::qsave(HiC_contact_matrix_list, paste0(parent.dir, "HiC_contact_matrix_list.qsave"))
# HiC_contact_matrix_list = qs::qread(paste0(parent.dir, "HiC_contact_matrix_list.qsave"))


##### 6. Define Evaluation Functions #####
# Calculate overlap metrics between predicted and Hi-C interactions
cal_overlap_regions = function(
    r1, # Predicted region-gene relationships (data frame with Region, Gene, Gene.loc)
    r2 # Hi-C bin-to-bin contacts (data frame with peak1, peak2, counts)
) {
  if (is.null(r1) | is.null(r2)) {return(c(precision = 0, recall = 0, f.score = 0))}

  # Extract primary region from semicolon-separated list
  r1$Region = strsplit(r1$Region, split = ";") %>% sapply(., "[", 1)
  r1 = r1[!is.na(r1$Gene.loc),, drop = FALSE]

  # Find overlaps in both orientations
  forward1 = findOverlaps(StringToGRanges(r1$Region, sep = c(":", "-")), StringToGRanges(r2$peak1, sep = c(":", "-")))
  forward2 = findOverlaps(StringToGRanges(r1$Gene.loc, sep = c(":", "-")), StringToGRanges(r2$peak2, sep = c(":", "-")))
  reverse1 = findOverlaps(StringToGRanges(r1$Region, sep = c(":", "-")), StringToGRanges(r2$peak2, sep = c(":", "-")))
  reverse2 = findOverlaps(StringToGRanges(r1$Gene.loc, sep = c(":", "-")), StringToGRanges(r2$peak1, sep = c(":", "-")))

  # Combine forward and reverse matches
  hits1 = inner_join(as.data.frame(forward1), as.data.frame(forward2))
  hits2 = inner_join(as.data.frame(reverse1), as.data.frame(reverse2))
  hits = rbind(hits1, hits2)
  hits = hits[!duplicated(hits), , drop = FALSE]

  # Calculate precision, recall, and F1-score
  prec = length(unique(hits$queryHits)) / nrow(r1)
  recall = length(unique(hits$subjectHits)) / nrow(r2)
  if (prec + recall == 0) {f.score = 0} else {f.score = 2 * prec * recall / (prec + recall)}
  return(c(precision = prec, recall = recall, f.score = f.score))
}


##### 7. Evaluate All Methods #####
# 7.1 Calculate evaluation metrics for all methods and datasets
region_TG_benchmark = pbmcapply::pbmclapply(data.ls, mc.cores = parallel::detectCores() - 1, function(chip.data) {
  message("----> Evaluating dataset: ", chip.data)
  data.scores = lapply(seq_along(method.list), function(i) {
    message("Method: ", names(method.list)[i])
    method = method.list[[i]]
    if (!chip.data %in% names(method)) {return(c(prec = 0, recall = 0, f.score = 0))}
    method.scores = cal_overlap_regions(r1 = method[[chip.data]], r2 = HiC_contact_matrix_list[[chip.data]])
    method.scores[is.na(method.scores)] = 0.0
    method.scores
  }) %>% Reduce("rbind", .)
  rownames(data.scores) = names(method.list)
  data.scores
})
names(region_TG_benchmark) = data.ls

# 7.2 Save and display results
print(region_TG_benchmark)
qs::qsave(region_TG_benchmark, paste0(parent.dir, "region_TG_benchmark.qsave"))
# save.image("3_peak-TG_relationship_benchmark.RData")
# region_TG_benchmark = qs::qread(paste0(parent.dir, "region_TG_benchmark.qsave"))

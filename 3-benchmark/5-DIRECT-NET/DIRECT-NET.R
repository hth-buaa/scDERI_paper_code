################################################
#                                              #
#    Run DIRECT-NET on benchmarking dataset    #
#                                              #
################################################


# ========== 1. Setup and Library Loading ==========
options(stringsAsFactors = FALSE)
library(tibble)
library(rlang)
library(dplyr)
library(presto)
library(JASPAR2016)
library(spatstat.sparse)
library(DIRECTNET)
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(patchwork)
library(ggplot2)
library(qs)
library(BSgenome.Hsapiens.UCSC.hg38)


# ========== 2. Create Promoter Regions File ==========
# Import GTF annotation
# download Basic gene annotation (CHR) https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.basic.annotation.gtf.gz from https://www.gencodegenes.org/human/
gtf_data = BiocIO::import("/mnt/i/scDERI_Paper/3-benchmark/5-DIRECT-NET/gencode.v49.chr_patch_hapl_scaff.basic.annotation.gtf.gz")

# Extract gene features
genes_gtf = subset(gtf_data, type == "gene")

# Define promoter regions (-2000/+200 for + strand, -200/+2000 for - strand)
promoter_regions = data.frame(
  Chrom = seqnames(genes_gtf),
  Starts = ifelse(strand(genes_gtf) == "+", start(genes_gtf) - 2000, start(genes_gtf) - 200),
  Ends = ifelse(strand(genes_gtf) == "+", end(genes_gtf) + 200, end(genes_gtf) + 2000),
  gene_id = genes_gtf$gene_id,
  gene_name = genes_gtf$gene_name,
  stringsAsFactors = FALSE
)

# Ensure start positions are positive
promoter_regions$Starts = pmax(1, promoter_regions$Starts)

# Format for DIRECT-NET input
genome.info = promoter_regions %>%
  dplyr::select(Chrom, Starts, Ends, genes = gene_name) %>%
  distinct()

# Save promoter regions file
write.table(genome.info,
            file = "/mnt/i/scDERI_Paper/3-benchmark/5-DIRECT-NET/hg38.promoter.regions.txt",
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
rm(list = ls())


# ========== 3. Define Input Files ==========
# Note: Dataset A549 and GM12878 produce errors; only dataset BM runs successfully
args = c(
  "/mnt/i/scDERI_Paper/3-benchmark/0-downsampling/Downsampled-data/10X_hg38_Human_bone_marrow_with_T_cell_depleted/10X_hg38_Human_bone_marrow_with_T_cell_depleted.qsave",
  "/mnt/i/scDERI_Paper/3-benchmark/0-downsampling/Downsampled-data/SNAREseq2_hg38_A549/SNAREseq2_hg38_A549.qsave",
  "/mnt/i/scDERI_Paper/3-benchmark/0-downsampling/Downsampled-data/SNAREseq2_hg38_GM12878/SNAREseq2_hg38_GM12878.qsave"
)


# ========== 4. Main Processing Loop ==========
# file.path = args[1]
# file.path = args[2]
# file.path = args[3]
for (n in 1:length(args)) {
  file.path = args[n]
  file.name = strsplit(file.path, split = "\\/") %>% `[[` (1) %>% tail(n = 1)
  data.dir = paste0("/mnt/i/scDERI_Paper/3-benchmark/5-DIRECT-NET/", gsub(".qsave", "", file.name))
  dir.create(data.dir)

  # ========== 5. Run DIRECTNET ==========
  setwd(data.dir)
  output_name = paste0(tools::file_path_sans_ext(file.name), "_directnet_out.qs")
  directnet_run = function(file.path, output_path, output_name){
    start_time = Sys.time()


    # 5.1 Load and process genome annotation
    genome.info = read.table(file = "/mnt/i/scDERI_Paper/3-benchmark/5-DIRECT-NET/hg38.promoter.regions.txt")
    names(genome.info) = c("Chrom","Starts","Ends","genes")
    genes = lapply(genome.info$genes, function(x) strsplit(x,"[|]")[[1]][1])
    genes = lapply(genes, function(x) strsplit(x,"[.]")[[1]][1])
    genes = unlist(genes)
    genome.info$genes = genes
    unik = !duplicated(genes) # filter out different transcript
    genome.info = genome.info[unik,]


    # 5.2 Load and preprocess Seurat object
    examp_3_3_3 = qread(file.path)
    # Convert Assay5 to Assay objects if needed
    if (grepl("Assay5", class(examp_3_3_3[["RNA"]]))) {examp_3_3_3[["RNA"]] = as(examp_3_3_3[["RNA"]], "Assay")}
    if (grepl("Assay5", class(examp_3_3_3[["ATAC"]]))) {examp_3_3_3[["ATAC"]] = as(examp_3_3_3[["ATAC"]], "Assay")}


    # 5.3 RNA processing
    DefaultAssay(examp_3_3_3) = "RNA"
    examp_3_3_3 = SCTransform(examp_3_3_3, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')


    # 5.4 ATAC processing
    DefaultAssay(examp_3_3_3) = "ATAC"
    examp_3_3_3 = RunTFIDF(examp_3_3_3)
    examp_3_3_3 = FindTopFeatures(examp_3_3_3, min.cutoff = 'q0')
    examp_3_3_3 = RunSVD(examp_3_3_3)
    examp_3_3_3 = RunUMAP(examp_3_3_3, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")


    # 5.5 Multi-modal integration
    examp_3_3_3 = FindMultiModalNeighbors(examp_3_3_3, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
    examp_3_3_3 = RunUMAP(examp_3_3_3, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
    examp_3_3_3 = FindClusters(examp_3_3_3, graph.name = "wsnn", algorithm = 3, verbose = FALSE)


    # 5.6 Identify cluster markers
    DefaultAssay(examp_3_3_3)="SCT"
    markers_all = presto:::wilcoxauc.Seurat(X = examp_3_3_3, group_by = 'wsnn_res.0.8', assay = 'data', seurat_assay = 'SCT')
    auc_threshold = 0.5
    markers_all = markers_all[which(markers_all$auc > auc_threshold), , drop = FALSE ]
    markers_all = markers_all[which(markers_all$padj < 0.05), , drop = FALSE]
    markers = data.frame(gene = markers_all$feature, group = markers_all$group)


    # 5.7 Run DIRECTNET (network inference)
    examp_3_3_3 = Run_DIRECT_NET(
      examp_3_3_3, peakcalling = FALSE, k_neigh = 50, atacbinary = TRUE, max_overlap = 0.5,
      size_factor_normalize = FALSE, genome.info = genome.info, focus_markers = markers_all$feature,
      nthread = parallel::detectCores() - 1
    )
    saveRDS(examp_3_3_3, file = "seuratobj.rds")

    # Extract DIRECT-NET results
    direct.net_result = Misc(examp_3_3_3, slot = 'direct.net')
    direct.net_result = as.data.frame(do.call(cbind,direct.net_result)) # links for markers

    # Standardize function type labels
    direct.net_result$function_type = gsub("HF", "HC", direct.net_result$function_type)
    direct.net_result$function_type = gsub("Rest", "MC", direct.net_result$function_type)
    direct.net_result$function_type = gsub("LF", "LC", direct.net_result$function_type)


    # 5.8 Identify differentially accessible peaks
    DefaultAssay(examp_3_3_3) = 'ATAC'
    focused_markers = markers
    groups = unique(focused_markers$group)
    da_peaks_list = list()
    Idents(examp_3_3_3) = examp_3_3_3$wsnn_res.0.8
    for (i in 1:length(groups)) {
      print(i)
      da_peaks = FindMarkers(
        object = examp_3_3_3,
        min.pct = 0.1,
        logfc.threshold = 0.25,
        ident.1 = groups[i],
        group.by = "wsnn_res.0.8",
        test.use = 'LR',
        only.pos = TRUE
      )
      da_peaks_list[[i]] = da_peaks
    }


    # 5.9 Generate regulatory networks
    # CRE-gene connections
    CREs_Gene = generate_CRE_Gene_links(direct.net_result, markers = focused_markers)

    # Find focused CREs which is overlapped with DA
    Focused_CREs = generate_CRE(L_G_record = CREs_Gene$distal, P_L_G_record = CREs_Gene$promoter, da_peaks_list)

    # Detect TFs for distal CREs
    L_TF_record = generate_peak_TF_links(peaks_bed_list = Focused_CREs$distal, species="Homo sapiens", genome = BSgenome.Hsapiens.UCSC.hg38, markers = focused_markers)

    # detect TFs for Promoters
    P_L_TF_record = generate_peak_TF_links(peaks_bed_list = Focused_CREs$promoter, species="Homo sapiens", genome = BSgenome.Hsapiens.UCSC.hg38, markers = focused_markers)


    # 5.10 Build network for visualization
    # Detect CRE-TF connections
    network_links = generate_links_for_Cytoscape(L_G_record = Focused_CREs$L_G_record, L_TF_record, P_L_G_record = Focused_CREs$P_L_G_record, P_L_TF_record,groups)

    # Ouput Node attribute
    Node_attribute = generate_node_for_Cytoscape(network_links, markers = focused_markers)


    # 5.11 Compile and save results
    list_output=list()
    list_output[[1]] = network_links
    list_output[[2]] = Node_attribute
    names(list_output) = c("link","node")
    qsave(list_output, paste(output_path, output_name, sep = "/"))

    end_time = Sys.time()
    message("Running time: ")
    print(end_time - start_time)
  }

  # Execute DIRECTNET analysis
  directnet_run(file.path = file.path, output_path = data.dir, output_name = output_name)
  save.image("directnet.RData")
}

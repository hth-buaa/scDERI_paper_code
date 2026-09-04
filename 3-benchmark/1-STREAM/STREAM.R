##################################################
#                                                #
#       Run STREAM on benchmarking dataset       #
#                                                #
##################################################


# General libraries
library(easypackages)
libs = c(
  "qs",
  "Seurat",
  "Signac",
  "pbmcapply",
  "pbapply",
  "parallel",
  # "EnsDb.Mmusculus.v79",
  "EnsDb.Hsapiens.v86",
  "dplyr",
  "ggplot2",
  "ggpubr",
  "igraph"
)
libraries(libs)

# devtools::install_github("OSU-BMBL/STREAM")
library(stream)


# Parameters
args = c(
  "/mnt/i/scDERI_Paper/3-benchmark/0-downsampling/Downsampled-data/10X_hg38_Human_bone_marrow_with_T_cell_depleted/",
  "/mnt/i/scDERI_Paper/3-benchmark/0-downsampling/Downsampled-data/SNAREseq2_hg38_A549/",
  "/mnt/i/scDERI_Paper/3-benchmark/0-downsampling/Downsampled-data/SNAREseq2_hg38_GM12878/"
)
# file.path = args[1]
# file.path = args[2]
# file.path = args[3]
for (n in 1:length(args)) {
  # Get org
  # n=2
  file.path = args[n]
  file.name = strsplit(file.path, split = "\\/") %>% `[[` (1) %>% tail(n = 1)
  data.dir = paste0("/mnt/i/scDERI_Paper/3-benchmark/1-STREAM/", gsub(".qsave", "", file.name))
  dir.create(data.dir)
  data.dir
  org = strsplit(file.name, split = "_") %>% `[[` (1) %>% `[` (2)
  print(org)


  # Get obj
  setwd(data.dir)
  obj = qs::qread(paste0(file.path, "/", file.name, ".qsave"))
  print(dim(obj))
  class(obj)
  if (grepl("Assay5", class(obj[["RNA"]]))) {obj[["RNA"]] = as(obj[["RNA"]], "Assay")}
  if (grepl("Assay5", class(obj[["ATAC"]]))) {obj[["ATAC"]] = as(obj[["ATAC"]], "Assay")}


  # Fix bugs, probably because the underlying code has changed after one year of development of my replica of STREAM.
  # Download https://figshare.com/ndownloader/files/38794185 to /mnt/i/scDERI_Paper/3-benchmark/1-STREAM/TFBS.list.rda
  run_stream_fixed = function(obj = NULL, qubic.path = "/users/PAS1475/liyang/software/QUBIC2/qubic",
                              candidate.TFs = NULL, peak.assay = "ATAC", var.genes = 3000,
                              top.peaks = 3000, ifPutativeTFs = FALSE, min.cells = 10,
                              out.dir = "./", org = "hg38", top.ngenes = 5, c.cutoff = 1,
                              n.blocks = 100, distance = 5e+05, filter_peaks_for_cicero = FALSE,
                              ifWeighted = TRUE, cicero.covar = -Inf, cicero.method = "covariance",
                              signac.score = -Inf, signac.method = "pearson", signac.pval = Inf,
                              intra.cutoff = 1, inter.cutoff = 0.8, peak.cutoff = 0.8,
                              score.cutoff = 1, KL = "min.exp", quantile.cutoff = 4, BlockOverlap = 0.5,
                              Extension = 0.9, submod.step = 30, min.eRegs = 100, url.link = "https://figshare.com/ndownloader/files/38794185") {
    # Quick assignment for testing
    # obj = obj; qubic.path = "/mnt/i/scDERI_Paper/3-benchmark/1-STREAM/qubic2-master/qubic";
    # candidate.TFs = NULL; peak.assay = "ATAC"; var.genes = 3000;
    # top.peaks = 3000; ifPutativeTFs = FALSE; min.cells = 10;
    # out.dir = "./"; org = org; top.ngenes = 5; c.cutoff = 1;
    # n.blocks = 100; distance = 5e+05; filter_peaks_for_cicero = FALSE;
    # ifWeighted = TRUE; cicero.covar = -Inf; cicero.method = "covariance";
    # signac.score = -Inf; signac.method = "pearson"; signac.pval = Inf;
    # intra.cutoff = 1; inter.cutoff = 0.8; peak.cutoff = 0.8;
    # score.cutoff = 1; KL = "min.exp"; quantile.cutoff = 4; BlockOverlap = 0.5;
    # Extension = 0.9; submod.step = 30; min.eRegs = 100; url.link = "https://figshare.com/ndownloader/files/38794185"

    start.time = Sys.time()
    message("Began running STREAM2 at time: \n", start.time)
    if (!dir.exists(out.dir)) {
      message("Creating the directory: ", out.dir, " to save the intermediate or final results ...")
      dir.create(out.dir)
    }
    stream:::quiet(peaks.use = rownames(obj[[peak.assay]])[as.vector(
      as.character(BSgenome::seqnames(obj[[peak.assay]]@ranges)) %in% GenomeInfoDb::standardChromosomes(obj[[peak.assay]]@ranges)
    )])
    obj = subset(obj, features = c(rownames(obj[["RNA"]])[!grepl("\\.", rownames(obj[["RNA"]]))], peaks.use))
    message("Filtered out features not a gene symbol or belonging to non-standard chromosomes,\n", "leading to ", nrow(obj[["RNA"]]), " genes and ", nrow(obj[[peak.assay]]), " enhancers.")
    qs::qsave(obj, paste0(out.dir, "Obj_filtered.qsave"))
    message("Saved the Seurat object after filtering to file: ", paste0(out.dir, "Obj_filtered.qsave."))
    set.seed(1234)
    stream:::quiet(if (org == "mm10") {
      require(BSgenome.Mmusculus.UCSC.mm10)
      org.gs = BSgenome.Mmusculus.UCSC.mm10
    } else if (org == "mm9") {
      require(BSgenome.Mmusculus.UCSC.mm9)
      org.gs = BSgenome.Mmusculus.UCSC.mm9
    } else if (org == "hg19") {
      require(BSgenome.Hsapiens.UCSC.hg19)
      org.gs = BSgenome.Hsapiens.UCSC.hg19
    } else if (org == "hg38") {
      require(BSgenome.Hsapiens.UCSC.hg38)
      org.gs = BSgenome.Hsapiens.UCSC.hg38
    } else {
      stop("Wrong organism assembly is input\n", "Please specify the correct organism assembly, e.g., hg38, hg19, mm10, or mm9")
    })
    message("Loaded full genome sequences for ", org)
    if (is.null(Signac::Annotation(object = obj[[peak.assay]]))) {
      message("Annotation slot does not exist in the ", peak.assay, " assay")
      annotations = Signac::GetGRangesFromEnsDb(org_to_DB(org = org))
      GenomeInfoDb::seqlevels(annotations) = paste0("chr", GenomeInfoDb::seqlevels(annotations))
      GenomeInfoDb::genome(annotations) = org
      GenomeInfoDb::seqlevelsStyle(annotations) = "UCSC"
      Signac::Annotation(object = obj[[peak.assay]]) = annotations
      message("Attach an annotation slot to the ", peak.assay, " assay")
    }
    Seurat::DefaultAssay(obj) = peak.assay
    obj = Signac::RegionStats(object = obj, assay = peak.assay, genome = org.gs)
    message("Computed the GC content, region lengths, and dinucleotide base frequencies", " for regions in the ", peak.assay, " assay and add them to the feature metadata.")
    qs::qsave(obj, paste0(out.dir, "Obj_quality_ctr.qsave"))
    links.df = stream:::filter_nearby_genes(obj = obj, peak.assay = peak.assay)
    qs::qsave(links.df, paste0(out.dir, "Potential_enhancer_gene_pairs_within_1mb.qsave"))
    message("Wrote potential enhancer-gene pairs to file ", out.dir, "Potential_enhancer_gene_pairs_within_1mb.qsave")
    obj = subset(x = obj, features = c(unique(links.df$peak), unique(links.df$gene)))
    # obj = obj[, Matrix::colSums(obj[["RNA"]]@counts) > 0 & Matrix::colSums(obj[[peak.assay]]@counts) > 0]
    obj = obj[, Matrix::colSums(GetAssayData(obj, assay = "RNA", layer = "counts")) > 0 &
                 Matrix::colSums(GetAssayData(obj, assay = peak.assay, layer = "counts")) > 0]
    message("Removed the cells with zero count of RNA or ATAC, \n", "leading to ", ncol(obj), " cells.")
    message("Filtered the enhancers that are beyond 500 kb from any genes in the Seurat object.\n",
            length(unique(links.df$peak)), " enhancers and ", length(unique(links.df$gene)), " genes were retained.")
    Seurat::DefaultAssay(obj) = "RNA"
    suppressWarnings(obj = Seurat::SCTransform(obj, verbose = FALSE, variable.features.n = var.genes))
    Seurat::DefaultAssay(obj) = peak.assay
    obj = Signac::RunTFIDF(obj)
    obj = Signac::FindTopFeatures(obj, min.cutoff = "q5")
    obj = subset(x = obj, features = c(obj[[peak.assay]]@var.features, obj[["SCT"]]@var.features))
    message("Saved the Seurat object composed of ", nrow(obj[["RNA"]]), " genes (", length(obj[["SCT"]]@var.features), " variable genes) and ",
            nrow(obj[[peak.assay]]), " enhancers (", length(obj[[peak.assay]]@var.features), " top-ranked enhancers) to file: ", out.dir, "Obj_var_genes_top_enhs.qsave.")
    qs::qsave(obj, paste0(out.dir, "Obj_var_genes_top_enhs.qsave"))
    options(timeout = 2000)
    # load(url(url.link))
    load("/mnt/i/scDERI_Paper/3-benchmark/1-STREAM/TFBS.list.rda")
    message("Loaded TF binding sites from JASPAR 2022.")
    TF.CRE.pairs = stream:::find_TFBS(peaks = rownames(obj[[peak.assay]]), TFBS.list = TFBS.list, org = org, candidate.TFs = candidate.TFs)
    qs::qsave(TF.CRE.pairs, paste0(out.dir, "TF_binding_sites_on_enhs.qsave"))
    bound.TFs = TF.CRE.pairs$CRE
    binding.CREs = TF.CRE.pairs$TF
    rm(TF.CRE.pairs)
    message("Identified ", length(binding.CREs), " TFs binding ", length(bound.TFs), " enhancers,\n", "which were saved to file: ", out.dir, "TF_binding_sites_on_enhs.qsave.")
    if (ifPutativeTFs) {
      stream:::quiet(require(JASPAR2022))
      pfm = TFBSTools::getMatrixSet(x = JASPAR2022, opts = list(collection = "CORE", tax_group = "vertebrates", all_versions = FALSE))
      puta.TF.peak.pairs = run_motifmatchr(pfm = pfm, peaks = rownames(obj[[peak.assay]]), org.gs = org.gs, candidate.TFs = candidate.TFs)
      qs::qsave(puta.TF.peak.pairs, paste0(out.dir, "Putative_TF_binding_sites_on_enhs.qsave"))
      puta.bound.TFs = puta.TF.peak.pairs$puta.bound.TFs
      puta.binding.peaks = puta.TF.peak.pairs$puta.binding.peaks
      rm(puta.TF.peak.pairs)
      message("Identified ", length(puta.binding.peaks), " TFs putatively binding ", length(puta.bound.TFs), " enhancers,\n", "which were saved to file: ", out.dir, "Putative_TF_binding_sites_on_enhs.qsave.")
    }
    message("\nBegan LTMG modeling ===================")
    LTMG.obj = stream:::call_LTMG(obj = obj)

    LTMG_discrete = LTMG.obj@LTMG@LTMG_discrete
    message("End LTMG modeling ===================")
    LTMG.obj = IRISFGM::CalBinaryMultiSignal(LTMG.obj)
    qs::qsave(LTMG.obj, paste0(out.dir, "LTMG_obj.qsave"))
    write.table(
      cbind(o = rownames(LTMG.obj@LTMG@LTMG_BinaryMultisignal), as.data.frame(LTMG.obj@LTMG@LTMG_BinaryMultisignal)),
      file = paste0(out.dir, "rna_counts.txt.chars"), quote = FALSE, sep = "\t", row.names = FALSE
    )
    message("Finished LTMG modeling for ", nrow(LTMG.obj@LTMG@LTMG_discrete), " genes across ", ncol(LTMG.obj@LTMG@LTMG_discrete),
            " cells.\n", "There are ", length(unique(range(LTMG.obj@LTMG@LTMG_discrete))), " transcriptional regulatory states (TRSs) identified.")
    message("\nBegan QUBIC biclustering ===================")
    system(paste0(qubic.path, " -i rna_counts.txt.chars -d -f ", BlockOverlap, " -o ", n.blocks, " -c ", Extension, "-k ", min.cells))
    message("End QUBIC biclustering ===================\n")
    qubic.cells = readLines(paste0(out.dir, "rna_counts.txt.chars.blocks")) %>% grep(pattern = "^ Conds ", value = TRUE) %>% strsplit(., split = " ")
    CoCond_cell = pbapply::pblapply(seq_along(qubic.cells), function(i) {
      x = qubic.cells[[i]]
      cell_name = x[-(1:3)]
      data.frame(cell_name = cell_name, Condition = rep(i, length(cell_name)))
    }) %>% Reduce("rbind", .)
    message("Read co-conditioned cells (n = ", length(unique(CoCond_cell$cell_name)), ") in ", length(unique(CoCond_cell$Condition)), " biclusters")
    qubic.genes = readLines(paste0(out.dir, "rna_counts.txt.chars.blocks")) %>% grep(pattern = "^ Genes ", value = TRUE) %>% strsplit(., split = " ")
    CoReg_gene = pbapply::pblapply(seq_along(qubic.genes), function(i) {
      x = qubic.genes[[i]]
      Gene = x[-(1:3)] %>% strsplit(., split = "_") %>% sapply(., "[", 1) %>% unique
      data.frame(Gene = Gene, Condition = rep(i, length(Gene)))
    }) %>% Reduce("rbind", .)
    CoReg_gene = CoReg_gene[CoReg_gene$Condition %in% unique(CoCond_cell$Condition), ]
    message("Read co-regulated genes (n = ", length(unique(CoReg_gene$Gene)), ") in ", length(unique(CoReg_gene$Condition)), " biclusters")
    rna.dis = subset(x = LTMG_discrete, rownames(LTMG_discrete) %in% rownames(Seurat::GetAssayData(obj, assay = "RNA", slot = "data")))
    atac.dis = stream:::binarize(Seurat::GetAssayData(object = obj, slot = "data", assay = peak.assay))
    link.ratio = sapply(apply(links.df, 2, unique), length)
    if (length(unique(CoReg_gene$Condition)) == 0) {
      stop("No bicluster was discovered by IRIS-FGM\n", "Please perform quality control or denoising before analysis")
    } else if (length(unique(CoReg_gene$Condition)) < 2) {
      block.list = list(unique(CoReg_gene$Gene))
    } else {
      block.list = split(CoReg_gene, f = CoReg_gene$Condition) %>% sapply(., "[[", "Gene") %>% head(n = n.blocks)
    }
    obj.list = stream:::subset_object(CoCond_cell = CoCond_cell, object = obj,
                                      links.df = links.df, atac.dis = atac.dis, n.blocks = n.blocks,
                                      max.peaks = round(nrow(rna.dis) * link.ratio[1]/link.ratio[2]),
                                      peak.assay = peak.assay)
    flags = sapply(obj.list, stream:::is.not.null)
    obj.list = obj.list[flags]
    block.list = block.list[flags]
    if (length(block.list) < 1) {stop("No IRIS-FGM bicluster passed quality control\n", "Please Please perform quality control or denoising before analysis")}

    G.list = stream:::build_graph(obj.list = obj.list, obj = obj, rna.dis = rna.dis,
                                  atac.dis = atac.dis, ifWeighted = ifWeighted, distance = distance,
                                  cicero.covar = cicero.covar, cicero.method = cicero.method,
                                  org.gs = org.gs, peak.assay = peak.assay, signac.score = signac.score,
                                  signac.pval = signac.pval, signac.method = signac.method,
                                  filter_peaks_for_cicero = filter_peaks_for_cicero, min.cells = min.cells,
                                  out.dir = out.dir)
    flags = sapply(G.list, stream:::is.not.null)
    G.list = G.list[flags]
    if (length(G.list) < 1) {stop("No heterogeneous graph was constructed\n", "Please perform quality control or denoising before analysis")}
    obj.list = obj.list[flags]
    block.list = block.list[flags]
    qs::qsave(obj.list, paste0(out.dir, "Obj_list.qsave"))
    qs::qsave(G.list, paste0(out.dir, "Graph_list.qsave"))
    if (length(G.list) < 1) {stop("No heterogeneous graph was constructed")}
    message("List containing ", length(obj.list), " Seurat objects were saved to file: ", out.dir, "Obj_list.qsave.")
    message("Constructed ", length(G.list), " heterogeneous graphs,\n", "which were saved to file: ", out.dir, "Graph_list.qsave.")
    TFGene.pairs = stream:::find_TF_gene(G.list = G.list, bound.TFs = bound.TFs, binding.CREs = binding.CREs)
    qs::qsave(TFGene.pairs, paste0(out.dir, "TF_target_pairs.qsave"))
    message("Saved the TF-target pairs to file: ", out.dir, "TF_target_pairs.qsave.")
    if (ifPutativeTFs) {
      puta.TF.gene.pairs = stream:::find_TF_gene(G.list = G.list, bound.TFs = puta.bound.TFs, binding.CREs = puta.binding.peaks)
      qs::qsave(puta.TF.gene.pairs, paste0(out.dir, "Putative_TF_target_pairs.qsave"))
      message("Saved the putative TF-target pairs to file: ", out.dir, "Putative_TF_target_pairs.qsave.")
    }
    if (ifPutativeTFs) {
      if (grepl("^mm", org)) {tax.id = 10090} else {tax.id = 9606}
      stream:::quiet(require(JASPAR2022))
      pfm = TFBSTools::getMatrixSet(x = JASPAR2022, opts = list(collection = "CORE", tax_group = "vertebrates", species = tax.id, all_versions = FALSE))
      stream:::quiet(require(motifmatchr))
      stream:::quiet(obj = Signac::AddMotifs(object = obj, genome = org.gs, pfm = pfm))
      motif.obj = subset(obj, features = intersect(rownames(obj[[peak.assay]]@motifs), rownames(obj[[peak.assay]])))
    } else {motif.obj = NULL}
    seeds = stream:::SFP_seeding(motif.obj = motif.obj, G.list = G.list,
                                 obj.list = obj.list, block.list = block.list, bound.TFs = bound.TFs,
                                 binding.CREs = binding.CREs, TFGene.pairs = TFGene.pairs,
                                 score.cutoff = score.cutoff, quantile.cutoff = quantile.cutoff,
                                 ifPutativeTFs = ifPutativeTFs, rna.dis = rna.dis, atac.dis = atac.dis,
                                 KL = KL, peak.assay = peak.assay)
    if (length(seeds) < 1) {stop("No seeds was identified\n", "Please perform quality control or denoising before analysis")}
    message(length(seeds), " seeds are identified for hybrid biclustering,\n", "which were saved to file: ", out.dir, "Seeds.qsave.")
    qs::qsave(seeds, paste0(out.dir, "Seeds.qsave"))
    rna.list = stream:::get_matrix_list(m = rna.dis, obj.list = obj.list, assay = "RNA")
    atac.list = stream:::get_matrix_list(m = atac.dis, obj.list = obj.list, assay = peak.assay)
    rm(obj.list)
    message("Biclustering based on Tier-1 TF-enhancer binding annotations ...")
    HBCs = stream:::hybrid_biclust(seeds = seeds[(sapply(seeds, "[[", "Tier") == 1)], rna.list = rna.list, atac.list = atac.list,
                                   top.ngenes = top.ngenes, bound.TFs = bound.TFs, binding.CREs = binding.CREs,
                                   G.list = G.list, TFGene.pairs = TFGene.pairs, peak.cutoff = peak.cutoff,
                                   c.cutoff = c.cutoff, KL = KL, org.gs = org.gs, intra.cutoff = intra.cutoff,
                                   inter.cutoff = inter.cutoff, quantile.cutoff = quantile.cutoff,
                                   rna.dis = rna.dis, atac.dis = atac.dis, min.cells = min.cells)
    if (ifPutativeTFs) {
      message("Biclustering based on Tier-2 motif enrichment ...")
      puta.HBCs = stream:::hybrid_biclust(seeds = seeds[(sapply(seeds, "[[", "Tier") == 2)], rna.list = rna.list, atac.list = atac.list,
                                          top.ngenes = top.ngenes, bound.TFs = puta.bound.TFs,
                                          binding.CREs = puta.binding.peaks, G.list = G.list,
                                          TFGene.pairs = puta.TF.gene.pairs, peak.cutoff = peak.cutoff,
                                          c.cutoff = c.cutoff, KL = KL, org.gs = org.gs, intra.cutoff = intra.cutoff,
                                          inter.cutoff = inter.cutoff, quantile.cutoff = quantile.cutoff,
                                          rna.dis = rna.dis, atac.dis = atac.dis, min.cells = min.cells)
      HBCs = c(HBCs, puta.HBCs)
    }
    gene.range = sapply(HBCs, "[[", "genes") %>% sapply(., length) %>% range
    peak.range = sapply(HBCs, "[[", "peaks") %>% sapply(., length) %>% range
    cell.range = sapply(HBCs, "[[", "cells") %>% sapply(., length) %>% range
    message("Identified ", length(HBCs), " hybrid biclusters (HBCs),\n",
            "which were saved to file: ", out.dir, "HBCs.qsave.\n",
            "HBCs contain ", gene.range[1], "-", gene.range[2], " genes, ",
            peak.range[1], "-", peak.range[2], " enhancers, ", cell.range[1],
            "-", cell.range[2], " cells, and ", length(unique(sapply(HBCs, "[[", "TF"))), " TFs.")
    qs::qsave(HBCs, paste0(out.dir, "HBCs.qsave"))
    rm(TFGene.pairs)
    rm(atac.list)
    rm(atac.dis)
    rm(bound.TFs)
    if (length(HBCs) <= min.eRegs) {submod.HBCs = HBCs} else {
      message("Performing submodular optimization ...")
      sim.m = stream:::compute_sim(HBCs = HBCs)
      submod.obj = stream:::sub_mod(HBCs = HBCs, sim.m = sim.m, rna.list = rna.list,
                                    G.list = G.list, links.df = links.df, block.list = block.list,
                                    n.cells = ncol(rna.dis), obj = obj, min.eRegs = min.eRegs,
                                    step = submod.step, peak.assay = peak.assay)
      rm(sim.m)
      submod.HBCs = submod.obj$eRegs
      message("HBC scores after submodular optimization were written to file: ", out.dir, "Submodular_scores.qsave.")
      qs::qsave(submod.obj$obj, paste0(out.dir, "Submodular_scores.qsave"))
      end.time = Sys.time()
      message("Finished running STREAM2 at time: \n", end.time, "\n\n", "Total running time: ", end.time - start.time)
      return(submod.obj$eRegs)
    }
    message(length(submod.HBCs), " enhancer regulons (eRegulons) were discovered")
    return(submod.HBCs)
  }


  # Run STREAM
  en.regs = run_stream_fixed(
    obj = obj,
    qubic.path = "/mnt/i/scDERI_Paper/3-benchmark/1-STREAM/qubic2-master/qubic", # Search qubic with your browser and install it in the specified directory for reading.
    org = org
  )
  # Output Structure: List of eRegulons
  # Outer list: All identified eRegulons (e.g., MEF2A, MYB, MEIS2)
  # Inner list: Each eRegulon contains 11 core components:
  # 1. terminal (int): ID of the IRIS-FGM algorithm block for prediction tracking
  # 2. Tier (num): Confidence tier (1=high, 2=lower) for TF-enhancer relationships
  # 3. TF (chr): Core transcription factor name
  # 4. genes (chr): Vector of putative target genes
  # 5. peaks (chr): Vector of target enhancer/accessible region coordinates
  # 6. cells (chr): Vector of cell identifiers where the eRegulon is active
  # 7. atac.ratio (num): Proportion of cells with accessible peaks in active cells
  # 8. score (int): eRegulon score, often correlates with target gene count
  # 9. weight (num): Algorithm-computed weight for network analysis
  # 10. links (GRanges): Enhancer-gene connections as GRanges object
  # 11. seed (num): Random seed for reproducibility
  save.image(paste0(data.dir, "/", gsub(".qsave", "", file.name), ".RData"))
  qs::qsave(en.regs, paste0(data.dir, "/eRegulons.qsave"))
}

# https://pypi.org/project/scglue/#files
# https://anaconda.org/channels/bioconda/packages/scglue/overview
# https://github.com/gao-lab/GLUE
# https://scglue.readthedocs.io/zh-cn/latest/install.html

conda create -n scglue python=3.9 -y
conda activate scglue
# conda install -c conda-forge -c bioconda scglue -y
conda install -c conda-forge -c bioconda scglue pytorch-gpu -y



####################################################
#                                                  #
#               Stage 1: Preprocessing             #
#                                                  #
####################################################
import os
os.chdir("/mnt/i/scDERI_Paper/3-benchmark/4-scGLUE/")
os.getcwd()


# Notes
# 1. GLUE depends on numpy 1.24 or less
# 2. Bugs exist in "Check integration diagnostics"
# 3. Bugs may occur in "scanpy.pp.highly_variable_genes" when setting flavor as "seurat_v3".
#    Seting flavor as "seurat" can solve this issue


# Get parameters
import sys
name_arr = prefix.split('_')
org = "hg38"
print ("Parsed parameters for preprocessing")


# Modules
import sys
import anndata as ad
import networkx as nx
import scanpy as sc
import scglue
import os
import dill
import pickle
import pandas as pd
import numpy as np
from matplotlib import rcParams
from pybedtools import BedTool
print ("Imported modules for preprocessing")


# Parameters
prefix = "/mnt/i/scDERI_Paper/3-benchmark/4-scGLUE/10X_hg38_Human_bone_marrow_with_T_cell_depleted" # Process each dataset separately
# prefix = "/mnt/i/scDERI_Paper/3-benchmark/4-scGLUE/SNAREseq2_hg38_A549"                           # Process each dataset separately
# prefix = "/mnt/i/scDERI_Paper/3-benchmark/4-scGLUE/SNAREseq2_hg38_GM12878"                        # Process each dataset separately
rna_name = prefix + "-GEX.h5ad"
atac_name = prefix + "-ATAC.h5ad"
data_dir = prefix + "/"
if not os.path.exists(data_dir):
    os.makedirs(data_dir)
gtf_dir = "/mnt/i/scDERI_Paper/3-benchmark/4-scGLUE/gtf-files/"
print ("Preprared input files")


scglue.plot.set_publication_params()
rcParams["figure.figsize"] = (4, 4)
print ("Parsed parameters")


# Prepare GTF file
# Download to gtf_dir (Only the GTF of hg38 was used in my experiment.)
# https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.chr_patch_hapl_scaff.annotation.gtf.gz
# https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.chr_patch_hapl_scaff.annotation.gtf.gz
# https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.chr_patch_hapl_scaff.annotation.gtf.gz
# https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M1/gencode.vM1.annotation.gtf.gz
org2gtf = {
  "hg38": "gencode.v44.chr_patch_hapl_scaff.annotation.gtf.gz",
  "hg19": "gencode.v19.chr_patch_hapl_scaff.annotation.gtf.gz",
  "mm10": "gencode.vM25.chr_patch_hapl_scaff.annotation.gtf.gz",
  "mm9": "gencode.vM1.annotation.gtf.gz"
}


gtf_file = org2gtf[org]
print ("Determined the GTF file")


# Read data
rna = sc.read_h5ad(rna_name)
rna
atac = sc.read_h5ad(atac_name)
atac
print ("Loaded GEX and ATAC modalities")


# Value format conversion
rna_raw_copy = rna.raw.to_adata()
rna_raw_copy.X = rna_raw_copy.X.astype(int)
rna.raw = rna_raw_copy
rna = rna.raw.to_adata()
rna.layers['counts'] = rna
print ("Assign raw counts to annData of GEX")


atac_raw_copy = atac.raw.to_adata()
atac_raw_copy.X = atac_raw_copy.X.astype(int)
atac.raw = atac_raw_copy
atac = atac.raw.to_adata()
atac.layers['counts'] = atac
print ("Assign raw counts to annData of ATAC")


# Preprocess scRNA-seq data
rna.X, rna.X.data
rna.layers["counts"] = rna.X.copy()
sc.pp.log1p(rna)
sc.pp.highly_variable_genes(rna, n_top_genes = 2000, flavor = "seurat") # the authors of scglue forgot this
sc.pp.normalize_total(rna)
sc.pp.scale(rna)
sc.tl.pca(rna, n_comps = 20, svd_solver = "auto")
sc.pp.neighbors(rna, metric = "cosine")
sc.tl.umap(rna)
print ("Preprocessed GEX data")


# Preprocess scATAC-seq data
atac.X, atac.X.data
atac.var_names = atac.var_names.str.replace('.', '-')
scglue.data.lsi(atac, n_components = 20, n_iter = 15)
atac.obsm['X_lsi'][np.isnan(atac.obsm['X_lsi'])] = 0
sc.pp.neighbors(atac, use_rep = "X_lsi", metric = "cosine")
sc.tl.umap(atac)
print ("Preprocessed ATAC data")


# Obtain genomic coordinates
rna.var.head()
scglue.data.get_gene_annotation(
    rna, gtf = gtf_dir + gtf_file,
    gtf_by = "gene_name"
)
print ("Got GEX annotations")


atac.var_names[:5]
split = atac.var_names.str.split(r"[:-]")
atac.var["chrom"] = split.map(lambda x: x[0])
atac.var["chromStart"] = split.map(lambda x: x[1]).astype(int)
atac.var["chromEnd"] = split.map(lambda x: x[2]).astype(int)
print ("Got ATAC annotations")


# Graph construction
rna.var['hgnc_id'] = rna.var.index
keep_rna = pd.notna(rna.var['hgnc_id'])
rna = rna[:, keep_rna]
rna = rna[:, ~rna.var['chrom'].isnull()] # The authors of GLUE did not think of this step
guidance = scglue.genomics.rna_anchored_guidance_graph(rna, atac)
guidance
scglue.graph.check_graph(guidance, [rna, atac])
atac.var.head()
print ("Constructed guidance graph")


# Save preprocessed data files
with open(data_dir + "rna-pp.pkl", 'wb') as f:
    pickle.dump(rna, f)
with open(data_dir + "atac-pp.pkl", 'wb') as f:
    pickle.dump(atac, f)
nx.write_graphml(guidance, data_dir + "guidance.graphml.gz")
print ("Constructed guidance graph")


####################################################
#                                                  #
#           Stage 2: Model training                #
#                                                  #
####################################################


from itertools import chain

import anndata as ad
import itertools
import networkx as nx
import pandas as pd
import scanpy as sc
import scglue
import seaborn as sns
from matplotlib import rcParams
scglue.plot.set_publication_params()
rcParams["figure.figsize"] = (4, 4)
print ("Imported modules for training")


import os
work_dir = data_dir
os.chdir(work_dir)
os.getcwd()
print ("Working directory: " + work_dir)


# Read preprocessed data
import pickle
rna = pickle.load(open('rna-pp.pkl', 'rb'))
atac = pickle.load(open('atac-pp.pkl', 'rb'))
guidance = nx.read_graphml("guidance.graphml.gz")
rna
atac
guidance
print ("Read preprocessed data")


# Configure data
scglue.models.configure_dataset(
    rna, "NB", use_highly_variable=True,
    use_layer="counts", use_rep="X_pca"
)
print ("Configured GEX data")


scglue.models.configure_dataset(
    atac, "NB", use_highly_variable=True,
    use_rep="X_lsi"
)
print ("Configured ATAC data")


guidance_hvf = guidance.subgraph(chain(
    rna.var.query("highly_variable").index,
    atac.var.query("highly_variable").index
)).copy()
guidance_hvf
print ("Subsetted guidance graph")


# Train GLUE model
glue = scglue.models.fit_SCGLUE( # Run for two days.
    {"rna": rna, "atac": atac}, guidance_hvf,
    fit_kws={"directory": "glue"}
)
with open(os.path.join(os.getcwd(), "glue_model.pkl"), 'wb') as file: 
  pickle.dump(glue, file)
glue.save("glue.dill")
glue
print ("Trained GLUE model")


# # Check integration diagnostics
# dx = scglue.models.integration_consistency(
#     glue, {"rna": rna, "atac": atac}, guidance_hvf
# )
# dx
# print ("Checked integration diagnostics")
# 
# 
# _ = sns.lineplot(x="n_meta", y="consistency", data=dx).axhline(y=0.05, c="darkred", ls="--")
# 
# 
# # Apply model for cell and feature embedding
# rna.obsm["X_glue"] = glue.encode_data("rna", rna)
# atac.obsm["X_glue"] = glue.encode_data("atac", atac)
# print ("Applied model for cell and feature embedding")
# 
# 
# combined = ad.concat([rna, atac])
# combined
# print ("Combined GEX and ATAC")
# 
# 
# sc.pp.neighbors(combined, use_rep="X_glue", metric="cosine")
# sc.tl.leiden(combined, resolution=1.0)
# sc.tl.umap(combined)
# sc.pl.umap(combined, color=[ "leiden"], wspace=0.65)
# print ("Finished clustering on GEX+ATAC")


feature_embeddings = glue.encode_graph(guidance_hvf)
feature_embeddings = pd.DataFrame(feature_embeddings, index=glue.vertices)
feature_embeddings.iloc[:5, :5]


rna.varm["X_glue"] = feature_embeddings.reindex(rna.var_names).to_numpy()
atac.varm["X_glue"] = feature_embeddings.reindex(atac.var_names).to_numpy()
rna.varm["X_glue"]
print("Calculated GLUE embedding")


import pickle
pickle.dump(rna, open('rna-emb.pkl', 'wb'))
pickle.dump(atac, open('atac-emb.pkl', 'wb'))
nx.write_graphml(guidance_hvf, "guidance-hvf.graphml.gz")
print ("Saved GLUE embeddings")


####################################################
#                                                  #
#          Stage 3: Regulatory inference           #
#                                                  #
####################################################


import os
import anndata as ad
import networkx as nx
import numpy as np
import pandas as pd
import scglue
import seaborn as sns
# !pip install ipython
from IPython import display
from matplotlib import rcParams
from networkx.algorithms.bipartite import biadjacency_matrix
from networkx.drawing.nx_agraph import graphviz_layout


# Read intermediate results
import os
import pickle
rna = pickle.load(open('rna-emb.pkl', 'rb'))
atac = pickle.load(open('atac-emb.pkl', 'rb'))
guidance_hvf = nx.read_graphml("guidance-hvf.graphml.gz")
print ("Read intermediate results")


rna.var["name"] = rna.var_names
atac.var["name"] = atac.var_names
rna.var["name"]


genes = rna.var.query("highly_variable").index
peaks = atac.var.query("highly_variable").index
genes


# Cis-regulatory inference with GLUE feature embeddings
import numpy as np
features = pd.Index(np.concatenate([rna.var_names, atac.var_names]))
feature_embeddings = np.concatenate([rna.varm["X_glue"], atac.varm["X_glue"]])
features
print ("Cis-regulatory inference with GLUE feature embeddings")


skeleton = guidance_hvf.edge_subgraph(
    e for e, attr in dict(guidance_hvf.edges).items()
    if attr["type"] == "fwd"
).copy()
skeleton


reginf = scglue.genomics.regulatory_inference(
    features, feature_embeddings,
    skeleton=skeleton, random_state=0
)
reginf


gene2peak = reginf.edge_subgraph(
    e for e, attr in dict(reginf.edges).items()
    if attr["qval"] < 0.05
)
gene2peak
print ("Calcualted peak-to-gene linkages")


# Write the edges to the CSV file
import csv
with open("gene2peak.csv", mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['Source', 'Target', 'Weight'])  # Header row

    for edge in gene2peak.edges(data=True):
        source, target, attributes = edge
        weight = attributes['weight']
        writer.writerow([source, target, weight])
print ("Wrote peak-to-gene linakges to csv file")


# Construct TF-gene regulatory network from inferred cis-regulatory regions
# Draft a coexpression-based network
# Download the bed files to "/mnt/i/scDERI_Paper/3-benchmark/4-scGLUE/ChIPseq-data/"
# http://download.gao-lab.org/GLUE/cisreg/ENCODE-TF-ChIP-hg38.bed.gz
# http://download.gao-lab.org/GLUE/cisreg/ENCODE-TF-ChIP-hg19.bed.gz
# http://download.gao-lab.org/GLUE/cisreg/ENCODE-TF-ChIP-mm10.bed.gz
# http://download.gao-lab.org/GLUE/cisreg/ENCODE-TF-ChIP-mm9.bed.gz
org2motif = {
    # "hg38": "JASPAR2022-hg38.bed.gz",
    "hg38": "ENCODE-TF-ChIP-hg38.bed.gz",
    # "hg19": "JASPAR2022-hg19.bed.gz",
    "hg19": "ENCODE-TF-ChIP-hg19.bed.gz",
    # "mm10": "JASPAR2022-mm10.bed.gz",
    "mm10": "ENCODE-TF-ChIP-mm10.bed.gz",
    # "mm9": "JASPAR2022-mm9.bed.gz"
    "mm9": "ENCODE-TF-ChIP-mm9.bed.gz"
}
motif_bed = scglue.genomics.read_bed("/mnt/i/scDERI_Paper/3-benchmark/4-scGLUE/ChIPseq-data/" + org2motif[org])
motif_bed.head()
print ("Drafted coexpression-based network")


tfs = pd.Index(motif_bed["name"]).intersection(rna.var_names)
tfs.size


# !conda install loompy
import loompy
rna[:, np.union1d(genes, tfs)].write_loom("rna.loom")
np.savetxt("tfs.txt", tfs, fmt="%s")
rna[:, np.union1d(genes, tfs)]


# !pip install pyscenic
import pyscenic
import os
pyscenic.__version__
os.system("which pyscenic")


with loompy.connect("rna.loom") as ds:
    print("Row attributes:", ds.ra.keys())
    print("Column attributes:", ds.ca.keys())


# conda create -n pyscenic python=3.8.20 -y
# conda activate pyscenic
# conda install pyscenic -y
# which pyscenic
os.system("/home/hth/miniconda3/envs/pyscenic/bin/pyscenic grn rna.loom tfs.txt -o draft_grn.csv --seed 0 --num_workers 8 --cell_id_attribute obs_names --gene_attribute var_names")
print ("Finished running pySCENIC")

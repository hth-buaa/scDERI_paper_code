#####################################################
#                                                   #
#   Extract the gene symbols and save them in txt   #
#                                                   #
#####################################################


# 1. Parameters
data_dir = "/mnt/i/scDERI_Paper/1-preprocessing/10X_hg38_Human_bone_marrow_with_T_cell_depleted/10X-human-bone-marrow/"


# 2. Modules
import scanpy as sc


# 3. Load annData
print("Loading h5ad file...")
adata = sc.read_h5ad(data_dir + "bm_multiome_rna.h5ad")
print(f"Data loaded successfully! Contains {adata.n_vars} genes and {adata.n_obs} cells.")


# 4. Gene symbol conversion
gene_names = adata.var_names.tolist()
print(f"Successfully extracted {len(gene_names)} gene symbols.")
print("Top 10 gene examples:")
for i, gene in enumerate(gene_names[:10]):
  print(f"{i+1}. {gene}")


# 5. Save to file
output_file = data_dir + 'bm_multiome_rna.txt'
with open(output_file, 'w') as f:
  for item in gene_names:
    f.write("%s\n" % item)
print(f"Gene symbols have been saved to: {output_file}.")

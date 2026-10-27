import scanpy as sc
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sc.settings.verbosity = 1
sc.set_figure_params(dpi=100, facecolor="white")

adata = sc.read_h5ad("data/adata_raw.h5ad")
adata.var_names_make_unique()
print("Loaded:", adata)

# Note: this data is already log2(TPM/10+1) normalized per GEO description, not raw counts
print("X min/max:", adata.X.min(), adata.X.max())

# Basic QC metrics
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None)
print(adata.obs[["n_genes_by_counts", "total_counts", "pct_counts_mt"]].describe())

# Filter very low-quality cells/genes (data pre-filtered already but double check)
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
print("After filtering:", adata.shape)

# Keep raw layer
adata.layers["lognorm"] = adata.X.copy()

# HVG selection (data already log-transformed, use seurat flavor on log data)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
adata.raw = adata
adata_hvg = adata[:, adata.var.highly_variable].copy()

sc.pp.scale(adata_hvg, max_value=10)
sc.tl.pca(adata_hvg, svd_solver="arpack", n_comps=40)
sc.pp.neighbors(adata_hvg, n_neighbors=15, n_pcs=30)
sc.tl.leiden(adata_hvg, resolution=0.5, key_added="leiden", flavor="igraph", n_iterations=2, directed=False)
sc.tl.umap(adata_hvg)

print("Number of Leiden clusters:", adata_hvg.obs["leiden"].nunique())
print(adata_hvg.obs["leiden"].value_counts())

# Transfer results back to full adata for record-keeping
adata.obsm["X_umap"] = adata_hvg.obsm["X_umap"]
adata.obs["leiden"] = adata_hvg.obs["leiden"]
adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"]

adata.write("data/adata_processed.h5ad")
adata_hvg.write("data/adata_hvg.h5ad")
print("Saved processed data.")

# Plot UMAP colored by Leiden cluster and by known cell label
fig, axes = plt.subplots(1, 2, figsize=(16, 7))
sc.pl.umap(adata, color="leiden", ax=axes[0], show=False, title="Leiden Clusters (unsupervised)")
sc.pl.umap(adata, color="cell_label", ax=axes[1], show=False, title="Known Cell Labels (ground truth)")
plt.tight_layout()
plt.savefig("data/umap_comparison.png", dpi=150)
print("Saved umap_comparison.png")

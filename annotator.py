import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

adata = sc.read_h5ad("data/adata_processed.h5ad")
adata.var_names_make_unique()
adata.X = adata.layers["lognorm"]

# ---- Marker gene reference panel (curated, standard TME + melanoma markers) ----
marker_ref = {
    "T cell": ["CD3D", "CD3E", "IL7R", "TCF7", "CCL5"],
    "B cell": ["MS4A1", "CD79A", "CD79B", "BANK1", "CD19"],
    "NK cell": ["NKG7", "GNLY", "KLRD1", "GZMB", "PRF1"],
    "Macrophage": ["CD68", "LYZ", "AIF1", "FCER1G", "ITGAX"],
    "CAF": ["COL1A1", "COL1A2", "COL3A1", "PDGFRB", "FAP"],
    "Endothelial": ["PECAM1", "VWF", "CLDN5", "PLVAP", "ENG"],
    "Malignant (melanocytic)": ["MLANA", "DCT", "TYR", "PMEL", "MITF", "PRAME", "MCAM"],
}

# Known BRAFi/targeted-therapy drug-resistance program markers (from melanoma literature,
# e.g. Tirosh et al. 2016 "resistance program", Hugo et al., Sun et al. AXL-program)
resistance_markers = ["AXL", "NGFR", "EGFR", "WNT5A", "SERPINE1", "ROR2"]

# Score each cluster's average expression per marker set (simple, automated annotator)
def score_clusters(adata, marker_dict, groupby="leiden"):
    scores = pd.DataFrame(index=sorted(adata.obs[groupby].unique(), key=int))
    for label, genes in marker_dict.items():
        genes_present = [g for g in genes if g in adata.var_names]
        if not genes_present:
            scores[label] = np.nan
            continue
        sub = adata[:, genes_present]
        vals = np.asarray(sub.X.mean(axis=1)).flatten() if hasattr(sub.X, "mean") else sub.X.mean(axis=1)
        adata.obs["_score_tmp"] = vals
        scores[label] = adata.obs.groupby(groupby)["_score_tmp"].mean()
    return scores

cluster_scores = score_clusters(adata, marker_ref)
print("=== Cluster marker scores (mean log-expression) ===")
print(cluster_scores.round(2))

# Automated annotation: assign each cluster its top-scoring cell type
auto_annotation = cluster_scores.idxmax(axis=1)
print("\n=== Automated Cluster Annotations ===")
print(auto_annotation)

adata.obs["auto_celltype"] = adata.obs["leiden"].map(auto_annotation)

# Resistance program score (malignant clusters only)
res_score = score_clusters(adata, {"Resistance_program": resistance_markers})
print("\n=== Resistance-program score per cluster ===")
print(res_score.round(3))

# Combine into a summary table
summary = cluster_scores.round(2).copy()
summary["Resistance_score"] = res_score["Resistance_program"].round(3)
summary["Auto_annotation"] = auto_annotation
summary["N_cells"] = adata.obs["leiden"].value_counts()
summary.to_csv("data/cluster_annotation_summary.csv")
print("\nSaved cluster_annotation_summary.csv")

adata.write("data/adata_annotated.h5ad")

# Plot: UMAP colored by auto annotation + resistance score overlay on malignant cells
fig, axes = plt.subplots(1, 2, figsize=(16, 7))
sc.pl.umap(adata, color="auto_celltype", ax=axes[0], show=False, title="Automated Cell Type Annotation")

# highlight resistance score only for malignant-labeled cells
malig_mask = adata.obs["auto_celltype"] == "Malignant (melanocytic)"
adata.obs["resistance_score_plot"] = np.nan
genes_present = [g for g in resistance_markers if g in adata.var_names]
sub = adata[:, genes_present]
vals = np.asarray(sub.X.mean(axis=1)).flatten()
adata.obs.loc[malig_mask, "resistance_score_plot"] = vals[malig_mask.values]

sc.pl.umap(adata, color="resistance_score_plot", ax=axes[1], show=False,
           title="Resistance-program score\n(malignant cells only)", cmap="Reds")
plt.tight_layout()
plt.savefig("data/annotation_and_resistance.png", dpi=150)
print("Saved annotation_and_resistance.png")

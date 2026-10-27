import scanpy as sc
import pandas as pd
import numpy as np

adata = sc.read_h5ad("data/adata_annotated.h5ad")

# Map ground truth labels to same categories as our annotator for comparison
gt_map = {
    "Malignant": "Malignant (melanocytic)",
    "T cell": "T cell",
    "B cell": "B cell",
    "NK": "NK cell",
    "Macrophage": "Macrophage",
    "CAF": "CAF",
    "Endothelial": "Endothelial",
    "N/A": "N/A",
    "Unresolved": "Unresolved"
}
adata.obs["gt_mapped"] = adata.obs["cell_label"].map(gt_map)

valid = adata.obs[adata.obs["gt_mapped"].isin(
    ["Malignant (melanocytic)","T cell","B cell","NK cell","Macrophage","CAF","Endothelial"])].copy()
valid["auto_celltype"] = valid["auto_celltype"].astype(str)
valid["gt_mapped"] = valid["gt_mapped"].astype(str)

acc = (valid["auto_celltype"] == valid["gt_mapped"]).mean()
print(f"Per-cell annotation accuracy (auto vs ground truth): {acc*100:.1f}%  (n={len(valid)} cells)")

# Confusion-style crosstab
conf = pd.crosstab(valid["gt_mapped"], valid["auto_celltype"])
print("\nConfusion matrix (rows=ground truth, cols=predicted):")
print(conf)
conf.to_csv("data/confusion_matrix.csv")

# Malignant subcluster resistance summary
malig = adata.obs[adata.obs["auto_celltype"] == "Malignant (melanocytic)"]
print(f"\nTotal malignant cells identified: {len(malig)}")
print("\nMalignant cells per Leiden subcluster (patient-driven substructure):")
print(malig["leiden"].value_counts())

# recompute resistance score for malignant cells
adata.X = adata.layers["lognorm"]
resistance_markers = ["AXL", "NGFR", "EGFR", "WNT5A", "SERPINE1", "ROR2"]
genes_present = [g for g in resistance_markers if g in adata.var_names]
sub = adata[malig.index, genes_present]
vals = np.asarray(sub.X.mean(axis=1)).flatten()
thresh = pd.Series(vals, index=malig.index)
high_res_thresh = thresh.quantile(0.9)
n_high_res = (thresh >= high_res_thresh).sum()
print(f"\nHigh-resistance-score malignant cells (top 10%, score >= {high_res_thresh:.2f}): {n_high_res} cells")
print(f"This represents {n_high_res/len(malig)*100:.1f}% of malignant cells and {n_high_res/len(adata)*100:.2f}% of all cells")

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

adata = sc.read_h5ad("data/adata_annotated.h5ad")

# Figure 3: Confusion matrix heatmap
conf = pd.read_csv("data/confusion_matrix.csv", index_col=0)
plt.figure(figsize=(7,6))
sns.heatmap(conf, annot=True, fmt="d", cmap="Blues")
plt.title("Annotator Accuracy: Predicted vs Ground-Truth Cell Type")
plt.ylabel("Ground truth (published labels)")
plt.xlabel("Automated annotation (this study)")
plt.tight_layout()
plt.savefig("data/fig3_confusion_matrix.png", dpi=150)
plt.close()

# Figure 4: Resistance score distribution among malignant cells
adata.X = adata.layers["lognorm"]
resistance_markers = ["AXL", "NGFR", "EGFR", "WNT5A", "SERPINE1", "ROR2"]
genes_present = [g for g in resistance_markers if g in adata.var_names]
malig = adata[adata.obs["auto_celltype"] == "Malignant (melanocytic)"]
sub = malig[:, genes_present]
vals = np.asarray(sub.X.mean(axis=1)).flatten()

plt.figure(figsize=(7,5))
plt.hist(vals, bins=40, color="firebrick", alpha=0.8)
thresh = np.quantile(vals, 0.9)
plt.axvline(thresh, color="black", linestyle="--", label=f"Top 10% threshold ({thresh:.2f})")
plt.xlabel("Resistance-program score (mean log-expr of AXL/NGFR/EGFR/WNT5A/SERPINE1/ROR2)")
plt.ylabel("Number of malignant cells")
plt.title("Distribution of Drug-Resistance Signature Score\nAcross Malignant Melanoma Cells")
plt.legend()
plt.tight_layout()
plt.savefig("data/fig4_resistance_distribution.png", dpi=150)
plt.close()

# Figure 5: Cluster sizes / composition bar chart
comp = adata.obs["auto_celltype"].value_counts()
plt.figure(figsize=(8,5))
comp.plot(kind="bar", color=sns.color_palette("tab10", len(comp)))
plt.ylabel("Number of cells")
plt.title("Cell Type Composition (Automated Annotation)")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.savefig("data/fig5_composition.png", dpi=150)
plt.close()

print("Saved fig3, fig4, fig5")

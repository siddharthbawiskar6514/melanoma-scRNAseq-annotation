import scanpy as sc, pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, seaborn as sns
from scipy.stats import kruskal
from sklearn.metrics import classification_report

adata = sc.read_h5ad("data/adata_annotated.h5ad")

# 1. Per-class precision / recall
gt_map = {"Malignant": "Malignant (melanocytic)", "T cell": "T cell",
          "B cell": "B cell", "NK": "NK cell", "Macrophage": "Macrophage",
          "CAF": "CAF", "Endothelial": "Endothelial"}
gt = adata.obs["cell_label"].map(gt_map)
m = gt.notna()
print("=== Per-class performance ===")
print(classification_report(gt[m].astype(str),
      adata.obs.loc[m, "auto_celltype"].astype(str), zero_division=0))

# 2. Resistance score by patient (malignant cells only)
genes = [g for g in ["AXL", "NGFR", "EGFR", "WNT5A", "SERPINE1", "ROR2"]
         if g in adata.var_names]
malig = adata[adata.obs["auto_celltype"] == "Malignant (melanocytic)"].copy()
score = np.asarray(malig[:, genes].X.mean(axis=1)).flatten()
df = pd.DataFrame({"tumor": malig.obs["tumor"].astype(str).values,
                   "res_score": score})
df = df[df.groupby("tumor")["res_score"].transform("size") >= 20]

print("\n=== Resistance score by patient (malignant cells) ===")
print(df.groupby("tumor")["res_score"].describe().round(2))

stat, p = kruskal(*[g.values for _, g in df.groupby("tumor")["res_score"]])
print(f"\nKruskal-Wallis across patients: H = {stat:.1f}, p = {p:.2e}")

grand = df["res_score"].mean()
ss_total = ((df["res_score"] - grand) ** 2).sum()
ss_between = df.groupby("tumor")["res_score"].apply(
    lambda s: len(s) * (s.mean() - grand) ** 2).sum()
print(f"Share of score variance explained by patient identity: "
      f"{ss_between / ss_total:.0%}")

thr = df["res_score"].quantile(0.9)
top = df[df["res_score"] >= thr]
print(f"\nTop-10% cells (score >= {thr:.2f}) per patient:")
print(top["tumor"].value_counts())

plt.figure(figsize=(9, 4))
sns.boxplot(data=df, x="tumor", y="res_score")
plt.title("Resistance-signature score by patient (malignant cells)")
plt.xlabel("Patient (tumour ID)")
plt.tight_layout()
plt.savefig("data/fig6_resistance_by_patient.png", dpi=150)
print("\nSaved fig6_resistance_by_patient.png")

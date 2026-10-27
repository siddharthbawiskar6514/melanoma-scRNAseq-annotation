import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad

print("Loading raw text file...")
df = pd.read_csv("data/GSE72056.txt", sep="\t", index_col=0, low_memory=False)
print("Shape:", df.shape)
print(df.index[:6].tolist())

# First 3 rows are metadata
malignant_key = "malignant(1=no,2=yes,0=unresolved)"
celltype_key = "non-malignant cell type (1=T,2=B,3=Macro.4=Endo.,5=CAF;6=NK)"

tumor = df.loc["tumor"]
malignant = df.loc[malignant_key]
celltype = df.loc[celltype_key]

expr = df.drop(index=["tumor", malignant_key, celltype_key])
print("Expression matrix shape (genes x cells):", expr.shape)

expr = expr.apply(pd.to_numeric, errors="coerce")
print("Any NaN after coercion:", expr.isna().sum().sum())

adata = ad.AnnData(X=expr.T.values.astype(np.float32))
adata.obs_names = expr.columns.tolist()
adata.var_names = expr.index.tolist()
adata.obs["tumor"] = tumor.values
adata.obs["malignant_flag"] = malignant.values.astype(int)
adata.obs["celltype_flag"] = celltype.values.astype(int)

malig_map = {0: "unresolved", 1: "non-malignant", 2: "malignant"}
adata.obs["malignant"] = adata.obs["malignant_flag"].map(malig_map)

ct_map = {0: "N/A", 1: "T cell", 2: "B cell", 3: "Macrophage", 4: "Endothelial", 5: "CAF", 6: "NK"}
adata.obs["celltype_known"] = adata.obs["celltype_flag"].map(ct_map)

def combine(row):
    if row["malignant"] == "malignant":
        return "Malignant"
    elif row["malignant"] == "non-malignant":
        return row["celltype_known"]
    else:
        return "Unresolved"
adata.obs["cell_label"] = adata.obs.apply(combine, axis=1)

print(adata)
print(adata.obs["cell_label"].value_counts())

adata.write("data/adata_raw.h5ad")
print("Saved adata_raw.h5ad")

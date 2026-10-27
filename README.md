# Marker-based Cell Type Annotation and Resistance-Signature Scoring in Melanoma scRNA-seq (GSE72056)

A Python (Scanpy) pipeline that clusters single cells from melanoma tumours, labels the clusters automatically using curated marker genes, checks those labels against the published annotations, and scores malignant cells for a literature-derived drug-resistance gene signature.

## Question

Bulk RNA-seq averages expression over all cells in a tumour. Can single-cell data be used to (1) identify the cell types present with a simple, transparent marker-based annotator, and (2) look at how a known drug-resistance gene program varies across malignant cells?

## Data

- **GSE72056**, Tirosh et al., *Science* 2016: 4,645 single cells from 19 melanoma patients.
- Source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72056
- Values are log2(TPM/10 + 1), not raw counts. The dataset ships with the authors' labels (malignant / non-malignant and cell type), which are used here only for validation.
- The data file is **not** included in this repository. Download it from GEO and save it as `data/GSE72056.txt`.

## Methods

| Step | Script | What it does |
|------|--------|--------------|
| 1 | `load_data.py` | Separates metadata rows from the expression matrix, builds an AnnData object with the published labels attached |
| 2 | `pipeline.py` | QC filtering, 2,000 highly variable genes, scaling, PCA (40 components), neighbour graph, Leiden clustering (resolution 0.5), UMAP |
| 3 | `annotator.py` | Scores each cluster against marker panels for 7 cell types (T, B, NK, macrophage, CAF, endothelial, malignant) and assigns the top-scoring label; scores malignant cells for a 6-gene resistance signature (AXL, NGFR, EGFR, WNT5A, SERPINE1, ROR2) |
| 4 | `validate.py` | Compares automated labels with the published labels (accuracy, confusion matrix) |
| 5 | `final_figures.py` | Confusion-matrix heatmap, resistance-score histogram, cell-type composition bar chart |

## Results

- Automated annotation agreed with the published labels for **93.2%** of confidently labelled cells (see `validate.py` output and the confusion matrix).
- 1,189 cells were annotated as malignant.
- Resistance-signature scores vary across malignant cells (see `fig4_resistance_distribution.png`).

*(Add your figures here, for example `![UMAP](figures/umap_comparison.png)`.)*

## Limitations

- **Resistance "high" group is defined by a percentile cut-off.** Cells in the top 10% of the signature score were flagged as high. That is a fixed fraction by construction, so it does not by itself show that a distinct rare subpopulation exists. Evidence for that would need, for example, a bimodal score distribution, a per-patient breakdown, or an independent threshold.
- **Patient effects.** Malignant cells in this dataset tend to group by patient. No batch correction was applied, so some clustering structure and signature variation may reflect patient identity rather than resistance biology.
- **Signature genes come from the literature** and were not validated on this dataset.
- **Accuracy is dominated by the largest classes.** Per-class precision and recall should be reported alongside overall accuracy.
- Leiden clustering is a classical method; no deep-learning approach was tested.

## How to run

```bash
pip install -r requirements.txt
mkdir -p data
# download GSE72056 from GEO and save as data/GSE72056.txt
python load_data.py
python pipeline.py
python annotator.py
python validate.py
python final_figures.py
```

## Note on tooling

This pipeline was built with AI assistance (Claude) for code generation and drafting. The dataset, research question, and interpretation were chosen and reviewed by the author, who can explain each step.

## Author

Siddharth Bawiskar, B.Tech Biotechnology, Parul University

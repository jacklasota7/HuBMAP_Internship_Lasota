# Integrating and Interpreting Single-Cell Datasets
**Pittsburgh Supercomputing Center (PSC) — HuBMAP Project**  
PSC / Carnegie Mellon University

Analysis and integration of single-cell RNA sequencing (scRNAseq) data from the HuBMAP Consortium using standardized, reproducible pipelines. The goal is to transform raw biological data into clean, structured datasets ready for exploration, visualization, and cell-type annotation. As well as, building accessible tools for interpreting the results.

---

## Project Goals

| Goal | Description |
|---|---|
| **Objective** | Process large-scale biological data through pipelines and provide accessible tools for its interpretation |
| **Standardization** | Uniformly process datasets so results are comparable across tissue types |
| **Products** | Clean, structured `.h5ad` datasets ready for visualization, exploration, or annotation |
| **Environment** | All computation runs in high-performance UNIX environments (PSC Bridges-2) |

---

## Tools & Dependencies

| Tool | Role |
|---|---|
| `scanpy` | Core single-cell analysis toolkit: QC, normalization, clustering, visualization |
| `anndata` | Primary data structure (`.h5ad`) for storing cell × gene matrices with metadata |
| `pandas` | DataFrame operations for cell and gene-level metadata |
| `matplotlib` | Plotting backend for Scanpy visualizations |
| Jupyter Notebook | Interactive analysis environment |
| Salmon-Alevin | Upstream sequencing pipeline (FASTQ → `.h5ad`) |

---

## Data Source

Data is pulled from the **[HuBMAP Consortium Data Portal](https://portal.hubmapconsortium.org/)**, which provides centralized access to human organ datasets. Datasets are searchable and browsable by organ, cell type, or molecule, and are available in both raw and processed formats. This project currently works with heart and kidney tissue datasets.

---

## Pipeline

### Stage 1: Raw Sequencing (Salmon-Alevin)

Raw sequencing data arrives in FASTQ format. The **Salmon-Alevin** pipeline handles:

- Quality control of raw reads
- Quantification of gene expression per cell
- Initial filtering of low-quality cells and genes
- Output as `.h5ad` format (AnnData-compatible)

This stage runs upstream on PSC infrastructure and produces the raw `.h5ad` files that serve as input to the notebook pipeline below.

---

### Stage 2: Data Loading & Exploration (`heart.ipynb`)

The notebook reads raw `.h5ad` files using AnnData and performs an initial structural inspection before any processing.

```python
import anndata as ad
heart = ad.read_h5ad("HT_raw.h5ad")
```

**Key AnnData components examined:**

| Attribute | Contents |
|---|---|
| `heart.X` | Sparse cell × gene count matrix (mostly zeros in raw data) |
| `heart.obs` | Cell-level metadata (e.g., `predicted_label` for cell type annotations) |
| `heart.var` | Gene-level metadata (e.g., `hugo_symbol`) |

Pandas is used to interrogate `.obs`. For example, `value_counts()` on `predicted_label` gives the distribution of cell types present in the dataset.

---

### Stage 3: Preprocessing & Quality Control

Initial filtering removes cells and genes that don't meet minimum thresholds:

```python
sc.pp.filter_cells(heart, min_genes=200)   # drop cells with too few detected genes
sc.pp.filter_genes(heart, min_cells=3)     # drop genes expressed in too few cells
```

**QC metrics are then calculated for three gene categories:**

| Category | Gene Prefix | Biological Significance |
|---|---|---|
| Mitochondrial | `MT-` | High % suggests dying or low-quality cells |
| Ribosomal | `RPS`, `RPL` | Flags cells with unusual translational activity |
| Hemoglobin | `^HB[^(P)]` | Flags potential red blood cell contamination |

```python
heart.var["mt"] = heart.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(heart, qc_vars=["mt"], inplace=True)
```

QC distributions are visualized via violin plots (`n_genes_by_counts`, `total_counts`) and scatter plots colored by mitochondrial percentage. Cells failing QC thresholds are filtered out before downstream analysis.

---

### Stage 4: Normalization & Feature Selection

Raw counts are preserved before normalization by copying to a separate layer:

```python
heart.layers["counts"] = heart.X.copy()   # preserve raw counts
sc.pp.normalize_total(heart)               # scale each cell to equal total counts
sc.pp.log1p(heart)                         # log-transform for variance stabilization
```

Highly variable genes (HVGs) are then identified. HVGs are the genes that carry the most biological signal and will be used in downstream dimensionality reduction:

```python
sc.pp.highly_variable_genes(heart, min_mean=0.0125, max_mean=3, min_disp=0.5)
```

---

### Stage 5: Dimensionality Reduction

**PCA** captures the major axes of variation across cells and reduces the feature space before clustering:

```python
sc.tl.pca(heart)
sc.pl.pca_variance_ratio(heart, n_pcs=50, log=True)   # determine how many PCs to use
```

PCA biplots are generated to check whether technical factors (e.g., mitochondrial %, total counts) are driving variation. This variation would be a signal that further normalization may be needed.

**UMAP** then projects cells into 2D for visual inspection of cluster structure:

```python
sc.pp.neighbors(heart)
sc.tl.umap(heart)
sc.pl.umap(heart, color=["total_counts"], size=2)
```

---

### Stage 6: Clustering

The processed pipeline runs **Leiden clustering**, which groups cells based on their neighborhood graph. Clusters are visualized on the UMAP:

```python
sc.pl.umap(heart_processed, color="leiden")
```

After the full pipeline runs, filtering reduces the dataset to roughly **1% of original cells** and **25% of original genes** which is a dramatic reduction, reflecting how much of the raw data is low-quality/uninformative.

---

### Stage 7: Data Integration (Multi-Dataset)

When working across tissue types (heart + kidney) or datasets (e.g., adding PBMCs to a reference), two integration strategies are used:

**Ingest Mapping** —> projects new cells into an existing reference UMAP and transfers cell-type annotations. Works well for clearly distinct cell types (T cells, B cells) but leaves some batch separation in subtler populations (monocytes).

**BBKNN** —> rebuilds the nearest-neighbor graph across batches, correcting for batch effects more aggressively. Improves mixing of monocytes and dendritic cells, though some rare clusters (Megakaryocytes) can be absorbed.

Both methods are used in combination, and results are compared visually on shared UMAPs.

---

## Data Outputs

| File | Format | Description |
|---|---|---|
| `HT_raw.h5ad` | AnnData / HDF5 | Raw heart scRNAseq data from HuBMAP portal |
| `HT_processed.h5ad` | AnnData / HDF5 | Fully processed data with Leiden clusters and UMAP embeddings |
| `*.json` (metadata) | JSON | Dataset-level metadata loaded alongside `.h5ad` for context |

---

## User Interface

Processed datasets and metadata feed into a **Shiny app** that provides tissue-type overviews and enables gene and cell-level comparisons without requiring direct notebook access. Thus, making results accessible to non-technical collaborators.

---

## Acknowledgements

Thanks to the HuBMAP HIVE teams, the CMU Tools Component team, and the PSC IEC team.

Special thanks to:
- **Matt Ruffalo** (CMU) —> Principal Investigator, Systems Scientist
- **Penny Cuda** (CMU) —> Research Programmer, HIVE Tools Component
- **Xiang Li** (PSC) —> Senior Bioinformatics Support Specialist

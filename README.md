# BinarySPA

**Binary-SPA** is a marker-based cell type annotation framework for spatial transcriptomics analysis.

---

## Installation

```r
# install.packages("devtools")
devtools::install_github("HonghaoNU/BinarySPA")
```

---

## Input Requirements

Binary-SPA accepts one of the following expression inputs:

1. Expression Data

10x Genomics HDF5 file (.h5)

Example: Xenium cell_feature_matrix.h5

Must contain a gene-by-cell count matrix

If multiple assays are present, the default assay "Gene Expression" is used (configurable via assay_name)


2. Expression matrix

A gene-by-cell matrix (matrix or dgCMatrix)

Row names: gene symbols

Column: cell IDs

---

## Quick Start Example

```r
library(BinarySPA)

res <- run_binary_spa(
  h5_file = "cell_feature_matrix.h5",
  marker_file = "markers.csv",
  delta_thresh = 0.15
)
```

---

## Output

`run_binary_spa()` returns:

- a CSV file containing cell-level annotations
- a CSV file summarizing cell-type proportions
- an .rds file bundling the two outputs above


## Method Overview

1. Load gene-by-cell matrix from 10x HDF5
2. Convert to cell-by-gene format
3. Align marker matrix to expression matrix
4. Compute binary marker scores
5. Assign cell types based on delta threshold

---

## Reproducibility

Tested with:

- R >= 4.2
- Seurat
- Matrix
- hdf5r

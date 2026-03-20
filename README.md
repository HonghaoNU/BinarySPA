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

Binary-SPA requires an expression matrix and a marker file.

---

## 1. Expression data (10x Genomics HDF5)

- File format: `.h5`  
- Example: `cell_feature_matrix.h5` (Xenium / Visium / 10x output)  
- Must contain a **gene-by-cell count matrix**  
- If multiple assays are present, the default assay **"Gene Expression"** is used  
  (configurable via `assay_name`)

---

## 2. Expression matrix

- A gene-by-cell matrix (`matrix` or `dgCMatrix`)  
- Row names: gene symbols  
- Column names: cell IDs  

---

## 3. Marker file (CSV required)

The marker file must be a **CSV file** containing a binary marker matrix.

- Rows = cell types  
- Columns = marker genes  
- Values must be `0` or `1`

### Format

- Use `1` to indicate that a gene is a marker for a cell type  
- Use `0` when the gene is not a marker  
- Each cell type can have multiple marker genes  
- A gene can be a marker for multiple cell types  
- All non-marker entries should be `0`  
- Column names should be gene symbols  
- Row names should be cell type labels
- Make sure the list of cell types is comprehensive, since Binary-SPA will assign every cell to one of the provided cell types  
- The marker file should be saved as a `.csv` file
  
### Example (markers.csv)

```
cell_type,CD3D,CD3E,MS4A1,LYZ,GAD1
T_cell,1,1,0,0,0
B_cell,0,0,1,0,0
Myeloid,0,0,0,1,0
Neuron,0,0,0,0,1
Mixed,1,0,1,0,0
```
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
2. Compute binary marker scores
3. Assign cell types

---

## Environment and Dependencies

Tested with:

- R >= 4.2
- Seurat
- Matrix
- hdf5r

---

## Reference

If you use Binary-SPA, please cite:

Binary-SPA: A marker-based framework for cell type annotation in spatial transcriptomics
https://www.biorxiv.org/content/10.64898/2026.03.17.712369v1

---

## License

Binary-SPA is free for academic and non-profit research use.

Commercial use requires a commercial license.


For licensing inquiries: honghao.bi@northwestern.edu, ji.peng@northwestern.edu

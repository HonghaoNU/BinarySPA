## Data source

The COAD Xenium data used in this example was downloaded from SPATCH:

https://spatch.pku-genomics.org/#/dataset/xenium

The dataset was published in:

Ren P, Zhang R, Wang Y, Zhang P, Luo C, Wang S, Li X, Zhang Z, Zhao Y, He Y, Li Y, Gao Z, Zhang X, Zhao Y, Liu Z, Meng Y, Zhang Z, Zeng Z. Systematic benchmarking of high-throughput subcellular spatial transcriptomics platforms across human tumors. Nature Communications. 2025;16:9232. https://pmc.ncbi.nlm.nih.gov/articles/PMC12534522/

## Required input

BinarySPA requires a 10x-style HDF5 count matrix.

If the downloaded SPATCH file is an AnnData `.h5ad` file, place it here:

```txt
example/data/adata.h5ad


## Installation

Create and activate a conda environment with R 4.2 and the required dependencies:

```bash
conda create -p /path/to/binaryspa -c conda-forge -c bioconda \
  r-base=4.2 \
  r-devtools \
  r-seurat \
  r-matrix \
  r-hdf5r \
  r-remotes \
  r-irkernel \
  git

conda activate /path/to/binaryspa

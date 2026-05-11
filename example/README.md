## Data source

The COAD Xenium data used in this example was downloaded from SPATCH:

https://spatch.pku-genomics.org/#/dataset/xenium

The dataset was published in:

Ren P, Zhang R, Wang Y, Zhang P, Luo C, Wang S, Li X, Zhang Z, Zhao Y, He Y, Li Y, Gao Z, Zhang X, Zhao Y, Liu Z, Meng Y, Zhang Z, Zeng Z. Systematic benchmarking of high-throughput subcellular spatial transcriptomics platforms across human tumors. Nature Communications. 2025;16:9232. https://pmc.ncbi.nlm.nih.gov/articles/PMC12534522/

## Required input

The SPATCH download was converted to a 10x-style HDF5 file named:

## Convert SPATCH h5ad to 10x HDF5

If the downloaded SPATCH file is an AnnData `.h5ad` file, place it here:

```txt
example/data/adata.h5ad

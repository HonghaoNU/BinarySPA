import anndata as ad
import h5py
import numpy as np
import scipy.sparse as sp

input_file = "example/data/adata.h5ad"
output_file = "example/data/example.h5"

adata = ad.read_h5ad(input_file)

# Use raw counts if available; otherwise use X
if adata.raw is not None:
    X = adata.raw.X
    genes = adata.raw.var_names.astype(str).to_numpy()
else:
    X = adata.X
    genes = adata.var_names.astype(str).to_numpy()

cells = adata.obs_names.astype(str).to_numpy()

# Seurat::Read10X_h5 expects genes x cells.
# AnnData is usually cells x genes, so transpose.
X = X.T

if not sp.issparse(X):
    X = sp.csc_matrix(X)
else:
    X = X.tocsc()

data = X.data
indices = X.indices
indptr = X.indptr
shape = np.array(X.shape, dtype=np.uint64)

gene_ids = genes
gene_names = genes
feature_types = np.array(["Gene Expression"] * len(genes), dtype="S")

with h5py.File(output_file, "w") as f:
    grp = f.create_group("matrix")

    grp.create_dataset("data", data=data)
    grp.create_dataset("indices", data=indices)
    grp.create_dataset("indptr", data=indptr)
    grp.create_dataset("shape", data=shape)

    grp.create_dataset("barcodes", data=cells.astype("S"))

    features = grp.create_group("features")
    features.create_dataset("id", data=gene_ids.astype("S"))
    features.create_dataset("name", data=gene_names.astype("S"))
    features.create_dataset("feature_type", data=feature_types)
    features.create_dataset("genome", data=np.array([""] * len(genes), dtype="S"))

print(f"Saved: {output_file}")
print(f"Genes x cells: {X.shape[0]} x {X.shape[1]}")

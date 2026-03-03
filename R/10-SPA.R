#' Perform Binary-SPA label transfer using Seurat
#'
#' Uses high-confidence labels (non-NA and not "Others") as reference and
#' transfers labels to NA-labeled cells using Seurat anchors.
#' Final metadata column added: "Binary_SPA".
#'
#' @param seu Seurat object containing labels.
#' @param label_col Metadata column with clean labels (default "Binary_CellType").
#' @param npcs Number of PCs (default 30).
#' @param dims Dimensions for anchors/UMAP (default 1:30).
#' @param n_neighbors UMAP neighbors for reference (default 5).
#' @param verbose Print progress (default TRUE).
#'
#' @return Seurat object with Binary_SPA added.
#' @export
binary_spa_transfer_labels <- function(
    seu,
    label_col = "Binary_CellType",
    npcs = 30,
    dims = 1:30,
    n_neighbors = 5,
    verbose = TRUE
) {

  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat package is required.", call. = FALSE)
  }
  if (!inherits(seu, "Seurat")) {
    stop("`seu` must be a Seurat object.", call. = FALSE)
  }
  if (!(label_col %in% colnames(seu@meta.data))) {
    stop("Label column not found in Seurat metadata: ", label_col, call. = FALSE)
  }

  lab <- seu@meta.data[[label_col]]
  if (is.null(lab)) stop("Label column is NULL: ", label_col, call. = FALSE)

  ref_cells <- rownames(seu@meta.data)[!is.na(lab) & lab != "Others"]
  query_cells <- rownames(seu@meta.data)[is.na(lab)]

  if (verbose) {
    message(sprintf("Reference cells: %d | Query cells: %d",
                    length(ref_cells), length(query_cells)))
  }

  # If no query, just copy labels over
  if (length(query_cells) == 0L) {
    seu$Binary_SPA <- as.character(lab)
    if (verbose) message("No query cells. Binary_SPA copied from ", label_col, ".")
    return(seu)
  }

  # Build reference/query objects safely
  ref <- seu[, ref_cells, drop = FALSE]
  query <- seu[, query_cells, drop = FALSE]


  # ---- Reference processing ----
  ref <- Seurat::NormalizeData(ref, verbose = FALSE)
  ref <- Seurat::FindVariableFeatures(ref, verbose = FALSE)
  ref <- Seurat::ScaleData(ref, verbose = FALSE)
  ref <- Seurat::RunPCA(ref, npcs = npcs, verbose = FALSE)
  ref <- Seurat::RunUMAP(ref, dims = dims, n.neighbors = n_neighbors,
                         return.model = TRUE, verbose = FALSE)

  # ---- Query processing ----
  query <- Seurat::NormalizeData(query, verbose = FALSE)
  query <- Seurat::FindVariableFeatures(query, verbose = FALSE)
  query <- Seurat::ScaleData(query, verbose = FALSE)
  query <- Seurat::RunPCA(query, npcs = npcs, verbose = FALSE)



  # ---- Anchors ----
  anchors <- Seurat::FindTransferAnchors(
    reference = ref,
    query = query,
    normalization.method = "LogNormalize",
    reference.reduction = "pca",
    dims = dims,
    verbose = FALSE
  )


  # ---- Map query ----
  query <- Seurat::MapQuery(
    anchorset = anchors,
    reference = ref,
    query = query,
    refdata = list(Binary_SPA = ref[[label_col, drop = TRUE]]),
    reference.reduction = "pca",
    reduction.model = "umap",
    verbose = FALSE
  )

  # ---- Merge back into full object ----
  seu$Binary_SPA <- NA_character_
  seu$Binary_SPA[ref_cells] <- as.character(ref[[label_col, drop = TRUE]])
  seu$Binary_SPA[query_cells] <- as.character(query$predicted.Binary_SPA)

  if (verbose) message("Binary_SPA label transfer complete.")
  seu
}

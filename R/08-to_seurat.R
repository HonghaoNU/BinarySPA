#' Create a Seurat object from Binary-SPA counts
#'
#' Creates a Seurat object from a genes x cells count matrix returned by
#' binary_spa_read_expr().
#'
#' @param counts A genes x cells count matrix (dgCMatrix / Matrix / matrix).
#' @param project Seurat project name (default "BinarySPA").
#' @param assay Assay name for Seurat (default "RNA").
#' @param min.cells Minimum cells expressing a gene (passed to CreateSeuratObject).
#' @param min.features Minimum features expressed in a cell (passed to CreateSeuratObject).
#'
#' @return A Seurat object.
#' @export
binary_spa_counts_to_seurat <- function(
  counts,
  project = "BinarySPA",
  assay = "RNA",
  min.cells = 0,
  min.features = 0
) {

  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat package is required to create a Seurat object.", call. = FALSE)
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Matrix package is required.", call. = FALSE)
  }

  if (!(inherits(counts, "Matrix") || is.matrix(counts))) {
    stop("`counts` must be a matrix or Matrix object (genes x cells).", call. = FALSE)
  }
  if (is.null(rownames(counts)) || is.null(colnames(counts))) {
    stop("`counts` must have rownames (genes) and colnames (cells).", call. = FALSE)
  }

  # Ensure sparse dgCMatrix for Seurat
  if (!inherits(counts, "dgCMatrix")) {
    counts <- Matrix::Matrix(counts, sparse = TRUE)
    counts <- methods::as(counts, "dgCMatrix")
  }

  Seurat::CreateSeuratObject(
    counts = counts,
    project = project,
    assay = assay,
    min.cells = min.cells,
    min.features = min.features
  )
}

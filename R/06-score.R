#' Compute Binary-SPA score matrix by matrix multiplication
#'
#' Aligns marker genes with a binary expression matrix and computes scores:
#' (cells x genes) %*% t(cell_types x genes) -> (cells x cell_types)
#'
#' @param marker_mat Marker matrix (cell_types x genes). Row names = cell types.
#' @param binary_mat Binary expression matrix (cells x genes), typically dgCMatrix of 0/1.
#' @param verbose Print alignment summary (default TRUE).
#'
#' @return Numeric matrix of scores (cells x cell_types).
#' @export
binary_spa_score_mm <- function(marker_mat, binary_mat, verbose = TRUE) {

  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Matrix package is required.", call. = FALSE)
  }

  if (is.null(rownames(marker_mat)) || is.null(colnames(marker_mat))) {
    stop("`marker_mat` must have rownames (cell types) and colnames (genes).", call. = FALSE)
  }
  if (is.null(rownames(binary_mat)) || is.null(colnames(binary_mat))) {
    stop("`binary_mat` must have rownames (cells) and colnames (genes).", call. = FALSE)
  }

  # 1) Align on shared genes (use marker order by default)
  common_genes <- intersect(colnames(marker_mat), colnames(binary_mat))
  if (length(common_genes) == 0L) {
    stop("No overlapping genes between marker_mat and binary_mat.", call. = FALSE)
  }

  # keep marker gene order (stable, reproducible)
  common_genes <- colnames(marker_mat)[colnames(marker_mat) %in% common_genes]

  if (verbose && (length(common_genes) < ncol(marker_mat) || length(common_genes) < ncol(binary_mat))) {
    message(sprintf(
      "Aligning on %d shared genes (dropping %d from marker_mat, %d from binary_mat).",
      length(common_genes),
      ncol(marker_mat) - length(common_genes),
      ncol(binary_mat) - length(common_genes)
    ))
  }

  M <- marker_mat[, common_genes, drop = FALSE]   # cell_types x genes
  B <- binary_mat[, common_genes, drop = FALSE]   # cells x genes

  # 2) Ensure numeric marker storage; keep B sparse
  M <- as.matrix(M)
  storage.mode(M) <- "double"

  # 3) Multiply: (cells x genes) %*% (genes x cell_types) = (cells x cell_types)
  #    Use Matrix methods; result will usually be a dense matrix (OK: cell_types is small)
  S <- B %*% Matrix::t(M)

  # 4) Reattach dimnames (Matrix multiplication should keep these, but enforce)
  rownames(S) <- rownames(B)  # cells
  colnames(S) <- rownames(M)  # cell types

  S
}

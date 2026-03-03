#' Transpose expression matrix to cells x genes
#'
#' Takes a genes x cells count matrix (e.g., output of binary_spa_read_expr())
#' and returns a cells x genes matrix for Binary-SPA downstream steps.
#'
#' @param counts A genes x cells expression matrix (dgCMatrix or matrix)
#'
#' @return dgCMatrix cells x genes
#' @export
binary_spa_counts_to_cells_x_genes <- function(counts) {

  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Matrix package is required.", call. = FALSE)
  }

  if (!(inherits(counts, "dgCMatrix") || is.matrix(counts) || inherits(counts, "Matrix"))) {
    stop(
      "`counts` must be a matrix or Matrix object (genes x cells). Got: ",
      paste(class(counts), collapse = ", "),
      call. = FALSE
    )
  }

  if (is.null(rownames(counts)) || is.null(colnames(counts))) {
    stop("`counts` must have rownames (genes) and colnames (cells).", call. = FALSE)
  }

  # Ensure sparse, then transpose
  if (!inherits(counts, "dgCMatrix")) {
    counts <- Matrix::Matrix(counts, sparse = TRUE)
  }

  out <- Matrix::t(counts)

  # Now: rows=cells, cols=genes
  out
}

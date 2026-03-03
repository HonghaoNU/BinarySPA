#' Binarize expression matrix
#'
#' Converts a cells x genes expression matrix into a binary 0/1 matrix.
#'
#' @param expr_sub A cells x genes expression matrix.
#' @param threshold Numeric cutoff (default 0).
#'
#' @return dgCMatrix (cells x genes), values 0 or 1.
#' @export
binary_spa_binarize <- function(expr_sub, threshold = 0) {

  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Matrix package is required.", call. = FALSE)
  }

  if (!(inherits(expr_sub, "Matrix") || is.matrix(expr_sub))) {
    stop("`expr_sub` must be matrix-like.", call. = FALSE)
  }

  if (is.null(rownames(expr_sub)) || is.null(colnames(expr_sub))) {
    stop("Matrix must have row and column names.", call. = FALSE)
  }

  # Replace NA with 0 (safety)
  expr_sub[is.na(expr_sub)] <- 0

  # Create numeric 0/1 sparse matrix
  binary_mat <- Matrix::Matrix((expr_sub > threshold) * 1, sparse = TRUE)
  binary_mat <- methods::as(binary_mat, "dgCMatrix")

  rownames(binary_mat) <- rownames(expr_sub)
  colnames(binary_mat) <- colnames(expr_sub)

  binary_mat
}

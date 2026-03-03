#' Read expression input (H5 or matrix)
#'
#' Automatically recognize input type:
#' - .h5 file → read using Seurat::Read10X_h5
#' - matrix / dgCMatrix → use directly
#'
#' @param input Either:
#'   - Path to .h5 file
#'   - A matrix (genes x cells)
#' @param assay_name Assay name if .h5 contains multiple assays
#'
#' @return dgCMatrix (genes x cells)
#' @export
binary_spa_read_expr <- function(input, assay_name = "Gene Expression") {

  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Matrix package is required.", call. = FALSE)
  }

  # -------------------------
  # Case 1: input is file path
  # -------------------------
  if (is.character(input) && length(input) == 1L) {

    if (!file.exists(input)) {
      stop("File not found: ", input, call. = FALSE)
    }

    # Detect .h5
    if (grepl("\\.h5$", input, ignore.case = TRUE)) {

      if (!requireNamespace("Seurat", quietly = TRUE)) {
        stop("Seurat package is required to read .h5 files.", call. = FALSE)
      }

      counts <- Seurat::Read10X_h5(filename = input)

      if (is.list(counts)) {

        nm <- names(counts)

        if (!is.null(nm) && assay_name %in% nm) {
          counts <- counts[[assay_name]]
        } else if (!is.null(nm) && "Gene Expression" %in% nm) {
          counts <- counts[["Gene Expression"]]
        } else {
          counts <- counts[[1]]
        }
      }

    } else {
      stop("Unsupported file format. Only .h5 currently supported.", call. = FALSE)
    }

  }

  # -------------------------
  # Case 2: input is matrix
  # -------------------------
  else if (is.matrix(input) || inherits(input, "Matrix")) {
    counts <- input
  }

  else {
    stop("Unsupported input type. Provide .h5 file or matrix.", call. = FALSE)
  }

  # -------------------------
  # Ensure sparse dgCMatrix
  # -------------------------
  if (!inherits(counts, "dgCMatrix")) {
    counts <- Matrix::Matrix(counts, sparse = TRUE)
  }

  # -------------------------
  # Validate dimensions
  # -------------------------
  if (is.null(rownames(counts)) || is.null(colnames(counts))) {
    stop("Expression matrix must have gene (rownames) and cell (colnames) names.",
         call. = FALSE)
  }

  counts
}

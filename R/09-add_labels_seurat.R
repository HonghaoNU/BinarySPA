#' Add Binary-SPA clean labels (named vector) to a Seurat object
#'
#' Adds `res$cell_annotations` (named character vector: names=cell IDs,
#' values=cell types for Clean cells, NA for Mixed) into Seurat metadata.
#'
#' @param sce A Seurat object.
#' @param cell_annotations Named character vector. Names must be cell IDs.
#' @param col_name Metadata column name to create (default "Binary_CellType").
#' @param verbose Print overlap summary (default TRUE).
#'
#' @return Seurat object with added metadata column.
#' @export
binary_spa_add_labels_to_seurat <- function(
    sce,
    cell_annotations,
    col_name = "Binary_CellType",
    verbose = TRUE
) {

  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat package is required.", call. = FALSE)
  }
  if (!inherits(sce, "Seurat")) {
    stop("`sce` must be a Seurat object.", call. = FALSE)
  }

  if (is.null(cell_annotations) || !is.atomic(cell_annotations)) {
    stop("`cell_annotations` must be a named atomic vector (character recommended).", call. = FALSE)
  }
  if (is.null(names(cell_annotations))) {
    stop("`cell_annotations` must be a NAMED vector: names = cell IDs.", call. = FALSE)
  }

  seurat_cells <- colnames(sce)
  overlap <- intersect(seurat_cells, names(cell_annotations))

  if (verbose) {
    message(sprintf(
      "Seurat cells: %d | annotation cells: %d | overlap: %d",
      length(seurat_cells), length(cell_annotations), length(overlap)
    ))
  }
  if (length(overlap) == 0L) {
    stop("No overlapping cell IDs between Seurat object and cell_annotations.", call. = FALSE)
  }

  # Create vector aligned to Seurat cell order
  labels <- rep(NA_character_, length(seurat_cells))
  names(labels) <- seurat_cells
  labels[overlap] <- as.character(cell_annotations[overlap])

  # Add to metadata
  sce[[col_name]] <- labels

  sce
}

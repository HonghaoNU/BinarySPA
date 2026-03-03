#' Read marker file and align markers to expression genes
#'
#' Reads a marker CSV with a `cell_type` column, sets rownames to cell types,
#' normalizes marker gene column names (replaces "_" with "-"), intersects genes
#' with the expression matrix, and returns a cleaned marker matrix.
#'
#' Expected marker file format:
#' - One column named `cell_type`
#' - Other columns are gene names, with marker weights / 0-1 / binary values
#'
#' @param markers Path to marker CSV file.
#' @param expr_mat_t Expression matrix in cells x genes format (e.g. output of
#'   binary_spa_counts_to_cells_x_genes()).
#' @param out_csv Optional path to write cleaned marker table (default NULL = no write).
#' @param verbose Print common gene summary (default TRUE).
#'
#' @return A data.frame marker_df_use: cell_types x common_genes
#' @export
binary_spa_read_markers_align <- function(
  markers,
  expr_mat_t,
  out_csv = NULL,
  verbose = TRUE
) {

  if (length(markers) != 1L || is.na(markers) || !nzchar(markers)) {
    stop("`markers` must be a single, non-empty file path.", call. = FALSE)
  }
  if (!file.exists(markers)) {
    stop("Marker file not found: ", markers, call. = FALSE)
  }

  if (is.null(colnames(expr_mat_t))) {
    stop("`expr_mat_t` must have gene names in colnames (cells x genes).", call. = FALSE)
  }

  marker_df <- utils::read.csv(markers, check.names = FALSE)

  if (!("cell_type" %in% colnames(marker_df))) {
    stop("Marker file must contain a `cell_type` column.", call. = FALSE)
  }

  if (anyDuplicated(marker_df$cell_type)) {
    stop("Duplicate values found in `cell_type` column; cell types must be unique.", call. = FALSE)
  }

  rownames(marker_df) <- marker_df$cell_type
  marker_df$cell_type <- NULL

  # Normalize gene naming: "_" -> "-" (common Xenium/10x mismatch)
  colnames(marker_df) <- gsub("_", "-", colnames(marker_df))

  common_genes <- intersect(colnames(expr_mat_t), colnames(marker_df))

  if (verbose) {
    cat("\nNumber of common genes:", length(common_genes), "\n")
    print(head(common_genes))
  }

  if (length(common_genes) == 0L) {
    stop(
      "No common genes found between expr_mat_t colnames and marker_df colnames.\n",
      "Example expr genes: ", paste(utils::head(colnames(expr_mat_t), 10), collapse = ", "), "\n",
      "Example marker genes: ", paste(utils::head(colnames(marker_df), 10), collapse = ", "),
      call. = FALSE
    )
  }

  marker_df_use <- marker_df[, common_genes, drop = FALSE]
  marker_df_use[is.na(marker_df_use)] <- 0

  # Optionally write cleaned marker table
  if (!is.null(out_csv)) {
    marker_df_use_out <- marker_df_use
    marker_df_use_out$cell_type <- rownames(marker_df_use_out)
    marker_df_use_out <- marker_df_use_out[, c("cell_type", setdiff(colnames(marker_df_use_out), "cell_type"))]
    utils::write.csv(marker_df_use_out, file = out_csv, row.names = FALSE)
  }

  marker_df_use
}

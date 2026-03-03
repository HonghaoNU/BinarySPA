#' Subset expression matrix to marker genes
#'
#' Creates expr_sub by keeping only genes shared between the expression matrix
#' (cells x genes) and the marker matrix (cell_types x genes). Gene order is
#' aligned to the marker matrix columns (or optionally expression columns).
#'
#' @param expr_mat_t Expression matrix in cells x genes format.
#' @param marker_df_use Marker matrix (cell_types x genes), typically from
#'   binary_spa_read_markers_align().
#' @param align_to Which gene order to use: "marker" (default) or "expr".
#' @param verbose Print overlap summary (default TRUE).
#'
#' @return A list with:
#'   - expr_sub: cells x common_genes
#'   - common_genes: character vector of genes used
#' @export
binary_spa_make_expr_sub <- function(
  expr_mat_t,
  marker_df_use,
  align_to = c("marker", "expr"),
  verbose = TRUE
) {
  align_to <- match.arg(align_to)

  if (is.null(colnames(expr_mat_t))) {
    stop("`expr_mat_t` must have gene names in colnames (cells x genes).", call. = FALSE)
  }
  if (is.null(colnames(marker_df_use))) {
    stop("`marker_df_use` must have gene names in colnames.", call. = FALSE)
  }

  expr_genes <- colnames(expr_mat_t)
  marker_genes <- colnames(marker_df_use)

  common_genes <- intersect(expr_genes, marker_genes)

  if (verbose) {
    cat("\nNumber of overlapping genes:", length(common_genes), "\n")
    print(head(common_genes))
  }

  if (length(common_genes) == 0L) {
    stop("No overlapping genes found between expr_mat_t and marker_df_use.", call. = FALSE)
  }

  if (align_to == "marker") {
    common_genes <- marker_genes[marker_genes %in% common_genes]
  } else {
    common_genes <- expr_genes[expr_genes %in% common_genes]
  }

  expr_sub <- expr_mat_t[, common_genes, drop = FALSE]

  list(expr_sub = expr_sub, common_genes = common_genes)
}

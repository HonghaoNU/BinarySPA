#' Run Binary-SPA End-to-End Annotation
#'
#' @description
#' Executes the complete Binary-SPA workflow starting from:
#' - Xenium cell_feature_matrix.h5
#' - Marker matrix CSV
#'
#' Returns annotated Seurat object and optionally saves outputs.
#'
#' @param h5_file Path to Xenium cell_feature_matrix.h5
#' @param marker_file Path to marker matrix CSV
#' @param delta_thresh Delta score threshold for label assignment (default = 0.15)
#' @param subset_n Optional number of cells to subset for testing (default = NULL)
#' @param output_prefix Prefix for output files (default = "Binary_SPA")
#' @param save_outputs Logical; whether to write RDS and CSV files (default = TRUE)
#'
#' @return A list containing:
#'   - seu
#'   - score_matrix
#'   - annotation_table
#'   - proportion_table
#'
#' @export
run_binary_spa <- function(
  h5_file,
  marker_file,
  delta_thresh = 0.15,
  subset_n = NULL,
  output_prefix = "Binary_SPA",
  save_outputs = TRUE
) {

  message("Step 1: Reading expression matrix...")
  counts_gxc <- binary_spa_read_expr(h5_file)

  if (!is.null(subset_n)) {
    cells_keep <- colnames(counts_gxc)[seq_len(min(subset_n, ncol(counts_gxc)))]
    counts_gxc <- counts_gxc[, cells_keep, drop = FALSE]

    message("Subset to ", ncol(counts_gxc), " cells.")
  }

  message("Step 2: Converting to cells x genes...")
  expr_mat_t <- binary_spa_counts_to_cells_x_genes(counts_gxc)

  message("Step 3: Aligning marker matrix...")
  marker_df_use <- binary_spa_read_markers_align(
    markers = marker_file,
    expr_mat_t = expr_mat_t
  )

  message("Step 4: Subsetting expression to markers...")
  sub <- binary_spa_make_expr_sub(expr_mat_t, marker_df_use)
  expr_sub <- sub$expr_sub

  message("Step 5: Binarizing expression...")
  binary_mat <- binary_spa_binarize(expr_sub)

  message("Step 6: Scoring marker matrix...")
  score_mat <- binary_spa_score_mm(marker_df_use, binary_mat)

  message("Step 7: Assigning labels...")
  res <- binary_spa_assign_labels_delta(
    score_mat,
    delta_thresh = delta_thresh
  )

  message("Step 8: Creating Seurat object...")
  seu <- binary_spa_counts_to_seurat(counts_gxc)

  seu <- binary_spa_add_labels_to_seurat(
    seu,
    res$cell_annotations,
    col_name = "Binary_CellType"
  )

  seu <- binary_spa_transfer_labels(
    seu,
    label_col = "Binary_CellType"
  )

  # Export cell annotations
  cell_annot <- data.frame(
    cell_id = colnames(seu),
    Binary_SPA = as.character(seu$Binary_SPA),
    Binary_SPA_confidence = as.character(seu$Binary_SPA_confidence),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  # Export proportions
  tab <- table(seu$Binary_SPA, useNA = "ifany")
  prop <- prop.table(tab)

  prop_df <- data.frame(
    Binary_SPA = names(tab),
    n_cells = as.integer(tab),
    proportion = as.numeric(prop),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  prop_df <- prop_df[order(prop_df$n_cells, decreasing = TRUE), ]

  if (save_outputs) {

    saveRDS(seu, file = paste0(output_prefix, "_seurat.rds"))

    write.csv(
      cell_annot,
      file = paste0(output_prefix, "_cell_annotations.csv"),
      row.names = FALSE
    )

    write.csv(
      prop_df,
      file = paste0(output_prefix, "_celltype_proportions.csv"),
      row.names = FALSE
    )

    message("Outputs written to disk.")
  }

  message("Binary-SPA completed successfully.")

  return(list(
    seu = seu,
    score_matrix = score_mat,
    annotation_table = cell_annot,
    proportion_table = prop_df
  ))
}

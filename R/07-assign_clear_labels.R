#' Assign Binary-SPA labels using normalized top1-top2 margin (delta)
#'
#' Takes a score matrix (cells x cell_types) and assigns each cell the top
#' cell type if delta >= delta_thresh; otherwise marks "Mixed".
#'
#' Cell types with all-zero scores are automatically removed.
#'
#' @param score_mat Score matrix (cells x cell_types).
#' @param delta_thresh Threshold for delta (default 0.15).
#' @param write_prefix If not NULL, write CSV outputs.
#' @param verbose Print messages.
#'
#' @return list(annot_final, cell_annotations, norm_scores)
#' @export
binary_spa_assign_labels_delta <- function(
  score_mat,
  delta_thresh = 0.15,
  write_prefix = NULL,
  verbose = TRUE
) {

  score_mat <- as.matrix(score_mat)

  if (is.null(rownames(score_mat)) || is.null(colnames(score_mat))) {
    stop("`score_mat` must have rownames (cells) and colnames (cell types).",
         call. = FALSE)
  }

  # ---- ALWAYS drop all-zero cell types (columns) ----
  keep_types <- colSums(score_mat != 0, na.rm = TRUE) > 0
  removed <- sum(!keep_types)

  score_mat <- score_mat[, keep_types, drop = FALSE]

  if (verbose && removed > 0) {
    message(removed, " cell types removed (all-zero columns).")
  }

  if (ncol(score_mat) == 0L) {
    stop("All cell types were removed (no non-zero scores).",
         call. = FALSE)
  }

  # ---- Normalize per cell type (column-wise max) ----
  type_max <- apply(score_mat, 2, max, na.rm = TRUE)
  norm_scores <- sweep(score_mat, 2, type_max, "/")
  norm_scores[!is.finite(norm_scores)] <- NA_real_

  # ---- Top1 / Top2 per cell ----
  get_top2 <- function(v) {
    ord <- order(v, decreasing = TRUE, na.last = NA)
    t1 <- if (length(ord) >= 1) ord[1] else NA_integer_
    t2 <- if (length(ord) >= 2) ord[2] else NA_integer_
    s1 <- if (!is.na(t1)) v[t1] else NA_real_
    s2 <- if (!is.na(t2)) v[t2] else NA_real_
    c(t1 = as.numeric(t1), t2 = as.numeric(t2), s1 = s1, s2 = s2, delta = s1 - s2)
  }

  top2_mat <- vapply(
    seq_len(nrow(norm_scores)),
    function(i) get_top2(norm_scores[i, ]),
    FUN.VALUE = c(t1 = 0, t2 = 0, s1 = 0, s2 = 0, delta = 0)
  )
  colnames(top2_mat) <- rownames(norm_scores)

  celltypes <- colnames(norm_scores)

  top1_idx  <- as.integer(top2_mat["t1", ])
  top2_idx  <- as.integer(top2_mat["t2", ])

  top1_type <- ifelse(is.na(top1_idx), NA_character_, celltypes[top1_idx])
  top2_type <- ifelse(is.na(top2_idx), NA_character_, celltypes[top2_idx])

  top1_norm <- as.numeric(top2_mat["s1", ])
  top2_norm <- as.numeric(top2_mat["s2", ])
  delta     <- as.numeric(top2_mat["delta", ])

  status <- ifelse(!is.na(delta) & delta >= delta_thresh, "Clean", "Mixed")

  annot_final <- data.frame(
    cell          = rownames(norm_scores),
    celltype      = top1_type,
    top1_norm     = top1_norm,
    top2_type     = top2_type,
    top2_norm     = top2_norm,
    purity_margin = delta,
    status        = status,
    stringsAsFactors = FALSE,
    row.names = NULL,
    check.names = FALSE
  )

  cell_annotations <- stats::setNames(
    ifelse(status == "Clean", top1_type, NA_character_),
    annot_final$cell
  )

  if (!is.null(write_prefix)) {

    delta_tag <- sprintf("delta%03d", round(delta_thresh * 100))

    out1 <- paste0(write_prefix, "_cell_annotations_Binary_", delta_tag, ".csv")
    out2 <- paste0(write_prefix, "_cell_annotations_simple_cleanOnly_", delta_tag, ".csv")

    utils::write.csv(annot_final, out1, row.names = FALSE)

    utils::write.csv(
      data.frame(cell = names(cell_annotations),
                 celltype = as.character(cell_annotations),
                 stringsAsFactors = FALSE),
      out2,
      row.names = FALSE
    )

    if (verbose) message("Wrote:\n- ", out1, "\n- ", out2)
  }

  list(
    annot_final = annot_final,
    cell_annotations = cell_annotations,
    norm_scores = norm_scores
  )
}

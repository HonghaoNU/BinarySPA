#' Compute Binary-SPA per-cell prediction confidence scores
#'
#' Standalone companion to \code{binary_spa_transfer_labels}. Performs the same
#' Seurat label-transfer workflow and returns a tidy data.frame of predicted
#' labels with their confidence scores (max prediction score from MapQuery),
#' and optionally per-class scores. Does NOT modify the input Seurat object
#' and does NOT depend on any prior call to \code{binary_spa_transfer_labels}.
#'
#' Reference cells (those with non-NA, non-"Others" labels) are reported with
#' \code{Binary_SPA_Score = 1} by convention; only query cells receive model-
#' derived scores.
#'
#' @param seu Seurat object containing labels.
#' @param label_col Metadata column with clean labels (default "Binary_CellType").
#' @param npcs Number of PCs (default 30).
#' @param dims Dimensions for anchors/UMAP (default 1:30).
#' @param n_neighbors UMAP neighbors for reference (default 5).
#' @param include_per_class Logical; if TRUE, append per-class prediction scores
#'   for query cells (NA for reference cells). Default FALSE.
#' @param query_only Logical; if TRUE, return only predicted cells, excluding
#'   reference cells. Default FALSE.
#' @param file Optional path to write a CSV. If NULL (default) nothing is written.
#' @param verbose Print progress (default TRUE).
#'
#' @return A data.frame with columns: \code{cell}, \code{Binary_SPA},
#'   \code{Binary_SPA_Score}, \code{Binary_SPA_Source}, and (optionally) one
#'   \code{score_<class>} column per cell type.
#' @export
binary_spa_confidence_scores <- function(
    seu,
    label_col = "Binary_CellType",
    npcs = 30,
    dims = 1:30,
    n_neighbors = 5,
    include_per_class = FALSE,
    query_only = FALSE,
    file = NULL,
    verbose = TRUE
) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat package is required.", call. = FALSE)
  }
  if (!inherits(seu, "Seurat")) {
    stop("`seu` must be a Seurat object.", call. = FALSE)
  }
  if (!(label_col %in% colnames(seu@meta.data))) {
    stop("Label column not found in Seurat metadata: ", label_col, call. = FALSE)
  }

  lab <- seu@meta.data[[label_col]]
  if (is.null(lab)) stop("Label column is NULL: ", label_col, call. = FALSE)

  ref_cells   <- rownames(seu@meta.data)[!is.na(lab) & lab != "Others"]
  query_cells <- rownames(seu@meta.data)[is.na(lab)]

  if (verbose) {
    message(sprintf("Reference cells: %d | Query cells: %d",
                    length(ref_cells), length(query_cells)))
  }

  # No query cells: return reference-only table early.
  if (length(query_cells) == 0L) {
    out <- data.frame(
      cell              = ref_cells,
      Binary_SPA        = as.character(lab[ref_cells]),
      Binary_SPA_Score  = 1,
      Binary_SPA_Source = "reference",
      stringsAsFactors  = FALSE,
      row.names         = NULL
    )
    if (isTRUE(query_only)) out <- out[0, , drop = FALSE]
    if (!is.null(file)) {
      utils::write.csv(out, file = file, row.names = FALSE)
      if (isTRUE(verbose)) message("Wrote ", nrow(out), " rows to: ", file)
    }
    return(out)
  }

  ref   <- seu[, ref_cells,   drop = FALSE]
  query <- seu[, query_cells, drop = FALSE]

  # ---- Reference processing ----
  ref <- Seurat::NormalizeData(ref, verbose = FALSE)
  ref <- Seurat::FindVariableFeatures(ref, verbose = FALSE)
  ref <- Seurat::ScaleData(ref, verbose = FALSE)
  ref <- Seurat::RunPCA(ref, npcs = npcs, verbose = FALSE)
  ref <- Seurat::RunUMAP(ref, dims = dims, n.neighbors = n_neighbors,
                         return.model = TRUE, verbose = FALSE)

  # ---- Query processing ----
  query <- Seurat::NormalizeData(query, verbose = FALSE)
  query <- Seurat::FindVariableFeatures(query, verbose = FALSE)
  query <- Seurat::ScaleData(query, verbose = FALSE)
  query <- Seurat::RunPCA(query, npcs = npcs, verbose = FALSE)

  # ---- Anchors ----
  anchors <- Seurat::FindTransferAnchors(
    reference = ref,
    query = query,
    normalization.method = "LogNormalize",
    reference.reduction = "pca",
    dims = dims,
    verbose = FALSE
  )

  # ---- Map query ----
  query <- Seurat::MapQuery(
    anchorset = anchors,
    reference = ref,
    query = query,
    refdata = list(Binary_SPA = ref[[label_col, drop = TRUE]]),
    reference.reduction = "pca",
    reduction.model = "umap",
    verbose = FALSE
  )

  # ---- Assemble scores ----
  ref_df <- data.frame(
    cell              = ref_cells,
    Binary_SPA        = as.character(ref[[label_col, drop = TRUE]]),
    Binary_SPA_Score  = 1,
    Binary_SPA_Source = "reference",
    stringsAsFactors  = FALSE,
    row.names         = NULL
  )

  query_df <- data.frame(
    cell              = query_cells,
    Binary_SPA        = as.character(query$predicted.Binary_SPA),
    Binary_SPA_Score  = as.numeric(query$predicted.Binary_SPA.score),
    Binary_SPA_Source = "predicted",
    stringsAsFactors  = FALSE,
    row.names         = NULL
  )

  out <- rbind(ref_df, query_df)

  # ---- Per-class scores (query cells only) ----
  if (isTRUE(include_per_class)) {
    score_assay <- "prediction.score.Binary_SPA"
    if (score_assay %in% Seurat::Assays(query)) {
      score_mat <- as.matrix(Seurat::GetAssayData(query, assay = score_assay))
      # rows = classes, cols = query cells -> transpose to cells x classes
      pc_df <- as.data.frame(t(score_mat), stringsAsFactors = FALSE,
                             check.names = FALSE)
      pc_df$cell <- rownames(pc_df)
      class_cols <- setdiff(colnames(pc_df), "cell")
      colnames(pc_df)[match(class_cols, colnames(pc_df))] <-
        paste0("score_", class_cols)
      out <- merge(out, pc_df, by = "cell", all.x = TRUE, sort = FALSE)
    } else if (isTRUE(verbose)) {
      warning("Per-class score assay 'prediction.score.Binary_SPA' not found; ",
              "skipping per-class export.", call. = FALSE)
    }
  }

  # ---- Optional filter ----
  if (isTRUE(query_only)) {
    out <- out[out$Binary_SPA_Source == "predicted", , drop = FALSE]
  }

  if (!is.null(file)) {
    utils::write.csv(out, file = file, row.names = FALSE)
    if (isTRUE(verbose)) message("Wrote ", nrow(out), " rows to: ", file)
  }

  out
}

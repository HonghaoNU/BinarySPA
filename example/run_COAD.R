
library(BinarySPA)

h5_file <- "example/data/example.h5"
marker_file <- "example/markers_COAD.csv"

if (!file.exists(h5_file)) stop("Missing file: ", h5_file)
if (!file.exists(marker_file)) stop("Missing file: ", marker_file)

dir.create("example/results", recursive = TRUE, showWarnings = FALSE)

res <- run_binary_spa(
  h5_file = h5_file,
  marker_file = marker_file,
  delta_thresh = 0.15,
  subset_n = 1000,
  output_prefix = "example/results/COAD_BinarySPA",
  save_outputs = TRUE
)

head(res$annotation_table)
res$proportion_table

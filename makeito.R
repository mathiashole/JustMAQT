#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(ape)
  library(RColorBrewer)
  library(readr)
  library(viridisLite)
})

# ---- Initialize variables ----
tree_file <- NULL
ids_file <- NULL
header_file <- NULL
# keywords <- c()
output_file <- NULL

# Plot
keywords <- c()
heatmap_file <- NULL

# Defoult colors
continuous_palette <- "viridis"
discrete_palette <- "Dark2"

# ---- Parsing arguments ----
args <- commandArgs(trailingOnly = TRUE)
i <- 1
while (i <= length(args)) {
  if (args[i] %in% c("--tree", "-t")) {
    tree_file <- args[i + 1]
    i <- i + 2
  } else if (args[i] %in% c("--ids", "-id")) {
    ids_file <- args[i + 1]
    i <- i + 2
  } else if (args[i] %in% c("--header", "-h")) {
    header_file <- args[i + 1]
    i <- i + 2
  } else if (args[i] %in% c("--keywords", "-k")) {
    # collect all keywords up to the next flag
    j <- i + 1
    while (j <= length(args) && !startsWith(args[j], "--") && !startsWith(args[j], "-")) {
      keywords <- c(keywords, args[j])
      j <- j + 1
    }
    i <- j
  } else if (args[i] %in% c("--out", "-o")) {
    output_file <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--heatmap") {
    heatmap_files <- c(heatmap_files, args[i + 1]); i <- i + 2
  } else if (args[i] == "--continuous-palette") {
    continuous_palette <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--discrete-palette") {
    discrete_palette <- args[i + 1]; i <- i + 2
  } else {
    stop(paste("Unknown argument:", args[i]))
  }
}

# ---- Validation ----
if (is.null(tree_file) && is.null(ids_file)) stop("Error: You must provide either --tree or --ids")
if (!is.null(tree_file) && !is.null(ids_file)) stop("Error: Provide only one of --tree or --ids, not both")
if (length(keywords) == 0 && length(heatmap_files) == 0) stop("Error: You must provide at least one keyword with --keywords or one heatmap with --heatmap")
if (is.null(header_file)) stop("Error: You must provide --header")
if (is.null(output_file)) stop("Error: You must provide --out")

# ---- Gets id ----
ids <- c()
if (!is.null(tree_file)) {
  tree <- read.tree(tree_file)
  ids <- tree$tip.label
} else {
  ids <- readLines(ids_file)
}

# ---- Read header ----
header <- readLines(header_file)

# ---- Colors palette ----
pal <- RColorBrewer::brewer.pal(max(3, length(keywords)), "Dark2")
names(pal) <- keywords

# ---- Find matches ----
matches <- lapply(ids, function(id) {
  found <- keywords[sapply(keywords, function(k) str_detect(id, k))]
  if (length(found) == 0) return(NULL)
  
  # Generate a row for each keyword that matches
  do.call(rbind, lapply(found, function(keyword) {
    data.frame(
      ID = id,
      symbol = 2,      # circle
      size = 10,
      color = pal[keyword],
      fill = 1,
      position = -1,
      label = keyword,
      stringsAsFactors = FALSE
    )
  }))
})

data_block <- do.call(rbind, matches)

# ---- Heatmap mode ----
if (!is.null(heatmap_file)) {
  df <- read_tsv(heatmap_file, show_col_types = FALSE)
  df[is.na(df)] <- "X"
  
  out_lines <- c(
    header,
    "DATA",
    apply(df, 1, function(x) paste(x, collapse = "\t"))
  )
}

# ---- Create final file ----
writeLines(c(
  header,
  "DATA",
  apply(data_block, 1, paste, collapse = ",")
), con = output_file)

cat("File save in:", output_file, "\n")


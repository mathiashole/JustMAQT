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
barplot_file <- NULL

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
    heatmap_file <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--barplot") {
    barplot_file <- args[i + 1]
    i <- i + 2
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

modes_selected <- sum(length(keywords) > 0, !is.null(heatmap_file), !is.null(barplot_file))
if(modes_selected == 0) stop("Error: You must provide either --keywords or --heatmap or --barplot")
if(modes_selected > 1) stop("Error: Choose only one mode: keywords OR heatmap OR barplot")

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


if(length(keywords) > 0) {

# ---- Colors palette ----
if (!is.null(discrete_palette) && str_detect(discrete_palette, "#")) {
  pal <- unlist(strsplit(discrete_palette, "\\s+"))
  if (length(pal) < length(keywords)) {
      stop("Error: Not enough colors provided for all keywords.")
    }
    pal <- pal[1:length(keywords)]
  } else {
    # Usar RColorBrewer con la paleta especificada
    pal <- RColorBrewer::brewer.pal(
      max(3, length(keywords)), 
      discrete_palette
    )
  }
  names(pal) <- keywords

# pal <- RColorBrewer::brewer.pal(max(3, length(keywords)), "Dark2")
# names(pal) <- keywords

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

  out_lines <- c(
    header,
    "DATA",
    apply(data_block, 1, paste, collapse = ",")
  )

}

# ---- Heatmap mode ----
if (!is.null(heatmap_file)) {
  # Detect file extension
  ext <- tools::file_ext(heatmap_file)
  
  if (ext %in% c("csv", "CSV")) {
    df <- read.csv(heatmap_file, stringsAsFactors = FALSE)
  } else {
    # Default to TSV
    df <- readr::read_tsv(heatmap_file, show_col_types = FALSE)
  }

  # df <- read_tsv(heatmap_file, show_col_types = FALSE)
  df[is.na(df)] <- "X"
  
  # Extract column names (except first one, whiche is ID)
  col_labels <- colnames(df)[-1]

  # Create the dynamic FIELD line
  field_labels_line <- paste("FIELD_LABELS", paste(col_labels, collapse = "\t"))

  # Replace in the header if FIELD_LABELS exists, if not add it
  header_mod <- gsub("^FIELD_LABELS.*", field_labels_line, header)
  if (identical(header, header_mod)) {
    #If you didn't find FIELD_LABELS in the header, we add it at the end
    header_mod <- c(header, field_labels_line)
  }

  out_lines <- c(
    header_mod,
    # "DATA",
    apply(df, 1, function(x) paste(x, collapse = "\t"))
  )
}

# ---- Procesar barplot ----
if (!is.null(barplot_file)) {
  # Detect file extension
  ext <- tools::file_ext(barplot_file)
  
  if (ext %in% c("csv", "CSV")) {
    df <- read.csv(barplot_file, stringsAsFactors = FALSE)
  } else {
    df <- readr::read_tsv(barplot_file, show_col_types = FALSE)
  }
  # Minimum validation: must have 2 or 3 columns
  if (ncol(df) < 2 || ncol(df) > 3) {
    stop("Barplot file must have 2 or 3 columns: ID,value[,label]")
  }
  # Craft output
  out_lines <- c(
    header,
    "DATA",
    apply(df, 1, function(x) paste(x, collapse = ","))
  )
  
}

# ---- Create final file ----

writeLines(out_lines, con = output_file)
cat("File saved in:", output_file, "\n")

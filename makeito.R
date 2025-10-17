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
output_file <- NULL

# Plot
keywords <- c()
keywords_shape <- c()
heatmap_file <- NULL
barplot_file <- NULL
binary_file <- NULL
multibarplot_file <- NULL
multibar_type <- "default"  # options: default, aligned, stacked

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
  } else if (args[i] == "--multibarplot") {
    multibarplot_file <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--multibar-type") {
    multibar_type <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--binary") {
    binary_file <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--continuous-palette") {
    continuous_palette <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--discrete-palette") {
    discrete_palette <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--keywords-shape") {
    # collect all shape codes up to the next flag
    j <- i + 1
    while (j <= length(args) && !startsWith(args[j], "--") && !startsWith(args[j], "-")) {
      keywords_shape <- c(keywords_shape, args[j])
      j <- j + 1
    }
    i <- j
  } else {
    stop(paste("Unknown argument:", args[i]))
  }
}

# ---- Validation ----
if (is.null(tree_file) && is.null(ids_file)) stop("Error: You must provide either --tree or --ids")
if (!is.null(tree_file) && !is.null(ids_file)) stop("Error: Provide only one of --tree or --ids, not both")
if (is.null(header_file)) stop("Error: You must provide --header")
if (is.null(output_file)) stop("Error: You must provide --out")

modes_selected <- sum(length(keywords) > 0, length(keywords_shape) > 0, !is.null(heatmap_file), !is.null(barplot_file), !is.null(binary_file))
if(modes_selected == 0) stop("Error: You must provide either --keywords or --heatmap or --barplot or --binary")
if (modes_selected > 1 && modes_selected != (length(keywords) > 0) + (length(keywords_shape) > 0)) {
  stop("Error: Only one mode allowed, except combining --keywords and --keywords-shape")
}

if (!tolower(multibar_type) %in% c("default", "aligned", "stacked")) {
  stop("Error: --multibar-type must be one of 'default', 'aligned', or 'stacked'.")
}

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
if (length(keywords) > 0) {
  if (!is.null(discrete_palette) && str_detect(discrete_palette, "#")) {
    pal <- unlist(strsplit(discrete_palette, "\\s+"))
    if (length(pal) < length(keywords)) stop("Error: Not enough colors provided for keywords.")
    pal <- pal[1:length(keywords)]
  } else {
    pal <- brewer.pal(max(3, length(keywords)), discrete_palette)
  }
  names(pal) <- keywords
} else {
  pal <- c()
}

# Parse --keywords-shape (format keyword:shape)
shape_map <- list()
if (length(keywords_shape) > 0) {
  for (entry in keywords_shape) {
    if (!str_detect(entry, ":")) stop("Each --keywords-shape must be in format keyword:shape (e.g., Retro:3)")
    parts <- str_split(entry, ":", simplify = TRUE)
    kw <- parts[1]; sh <- as.numeric(parts[2])
    shape_map[[kw]] <- sh
  }
}

# Combine all relevant keywords (for matching IDs)
all_keywords <- unique(c(keywords, names(shape_map)))

# ---- Find matches ----
matches <- lapply(ids, function(id) {
  found <- keywords[sapply(keywords, function(k) str_detect(id, k))]
  if (length(found) == 0) return(NULL)
  
  # Generate a row for each keyword that matches
  do.call(rbind, lapply(found, function(keyword) {
    shape <- if (!is.null(shape_map[[keyword]])) shape_map[[keyword]] else 2
    color <- if (!is.null(pal[keyword])) pal[keyword] else "#000000"
    data.frame(
      ID = id,
      # symbol = 2,      # circle
      symbol = shape,
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

# ---- Procesar multi-barplot ----
if (!is.null(multibarplot_file)) {
  # Detect file extension
  ext <- tools::file_ext(multibarplot_file)

  if (ext %in% c("csv", "CSV")) {
    df <- read.csv(multibarplot_file, stringsAsFactors = FALSE)
  } else {
    # Default to TSV
    df <- readr::read_tsv(multibarplot_file, show_col_types = FALSE)
  }

  # Validation: must have at least 3 columns (ID + 2+ values)
  if (ncol(df) < 3) {
    stop("Multi-barplot file must have at least 3 columns: ID and at least 2 numeric fields.")
  }

  # Replace NAs with 0
  df[is.na(df)] <- 0

  # Extract column names except ID of header
  col_labels <- colnames(df)[-1]
  n_fields <- length(col_labels)

  # # Assign automatic colors (Dark2 or Set2 palette)
  # pal <- brewer.pal(min(max(3, n_fields), 8), "Set2")[1:n_fields]
  # Generate color palette dynamically
  if (!is.null(discrete_palette) && str_detect(discrete_palette, "#")) {
    pal <- unlist(strsplit(discrete_palette, "\\s+"))
    if (length(pal) < n_fields) stop("Error: Not enough colors provided for fields.")
    pal <- pal[1:n_fields]
  } else {
    pal <- brewer.pal(min(max(3, n_fields), 8), discrete_palette)
  }

  # Generate dynamic header lines
  field_labels_line <- paste("FIELD_LABELS", paste(col_labels, collapse = ","))
  field_colors_line <- paste("FIELD_COLORS", paste(pal, collapse = ","))

  # Replace existing lines or append if missing
  if (any(grepl("^FIELD_LABELS", header))) {
    header <- gsub("^FIELD_LABELS.*", field_labels_line, header)
  } else {
    header <- c(header, field_labels_line)
  }

  if (any(grepl("^FIELD_COLORS", header))) {
    header <- gsub("^FIELD_COLORS.*", field_colors_line, header)
  } else {
    header <- c(header, field_colors_line)
  }

  # ---- Add multibar layout configuration ----
  if (tolower(multibar_type) == "aligned") {
    align_line <- "ALIGN_FIELDS,1"
    side_line <- "SIDE_STACKED,0"
  }

  # Build final dataset lines
  out_lines <- c(
    header,
    "DATA",
    apply(df, 1, function(x) paste(x, collapse = ","))
  )
}

# ---- Create final file ----
writeLines(out_lines, con = output_file)
cat("File saved in:", output_file, "\n")

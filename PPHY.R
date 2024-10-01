#!/usr/bin/env Rscript

# Load functions from functions.R
source("quick_off_init.R")

required_packages <- c("ggtree", "treeio", "ape", "RColorBrewer", "optparse", "cluster", "factoextra", "ggplot2") # , "DECIPHER"
manage_packages(required_packages)

# Define command line options
option_list <- list(
  make_option(c("--config"), type = "character", default = NULL, help = "Archivo de configuración", metavar = "FILE"),
  make_option(c("--phy"), type = "character", default = NULL, help = "Archivo de filogenia", metavar = "FILE"),
  make_option(c("--keyword"), type = "character", default = NULL, help = "Palabras para colorear separadas por espacios", metavar = "WORDS"),
  make_option(c("--alignment"), type = "character", default = NULL, help = "Archivo de alineamiento", metavar = "FILE"),
  make_option(c("--root"), type = "character", default = NULL, help = "Raíz del árbol", metavar = "NODE"),
  make_option(c("--layout"), type = "character", default = "rectangular", help = "Tipo de layout del árbol (rectangular, equal_angle, daylight, circular, roundrect)", metavar = "LAYOUT"),
   # make_option(c("--cluster"), type = "logical", default = FALSE, help = "Activar clustering automático", action = "store_true")
  make_option(c("--genotype"), type = "character", default = NULL, help = "Archivo de genotipo", metavar = "FILE"),
  make_option(c("--countineous"), type = "character", default = NULL, help = "Archivo de valor continuos", metavar = "FILE"),
  make_option(c("--cluster"), type = "character", default = "no", help = "Clustering automático: AUTO, número para especificar k, o NO para desactivar", metavar = "CLUSTER")
)

############
# FUNCTION #
############

# Parse arguments
opt <- parse_args(OptionParser(option_list = option_list))

if (!is.null(opt$config)) {
  config <- read.csv(opt$config, stringsAsFactors = FALSE, header = TRUE)
  
  for (i in seq_len(nrow(config))) {
    arg_name <- config$argument[i]
    arg_value <- config$value[i]
    
    if (is.null(opt[[arg_name]])) {
      opt[[arg_name]] <- arg_value
    }
  }
}

# Check required arguments
if (is.null(opt$phy)) {
  stop("You must specify a phylogeny file with --phy")
}

# Read phylogenetic tree
tree <- read.tree(opt$phy)

# Check if the tree contains tags
if (length(tree$tip.label) == 0) {
  stop("The tree does not contain tip tags.")
}

# Create a dataframe with the names and colors information
data <- data.frame(label = tree$tip.label)
data$color <- "grey"  # Default color

# Generate a color palette for the words
if (!is.null(opt$keyword)) {
  palabras <- unlist(strsplit(opt$keyword, " "))
  n_palabras <- length(palabras)
  colores <- brewer.pal(min(n_palabras, 12), "Dark2")  # Limit maximum palette size
  
  for (i in 1:n_palabras) {
    palabra <- palabras[i]
    data$color[grepl(palabra, data$label)] <- colores[i]
  }
}

# Root the tree si se especifica
if (!is.null(opt$root)) {
  tree <- root(tree, outgroup = opt$root)
}

# Function to read genotype data from different formats without using the tools library
read_data_table <- function(file) {
  file_ext <- tolower(sub(".*\\.", "", basename(file)))
  
  if (file_ext == "csv") {
    genotype_data <- read.csv(file, row.names = 1)
  } else if (file_ext %in% c("tsv", "txt")) {
    genotype_data <- read.delim(file, row.names = 1, sep = "\t")
  # } else if (file_ext %in% c("xlsx", "xls")) {
  #   library(readxl)
  #   genotype_data <- read_excel(file, col_names = TRUE)
  #   row.names(genotype_data) <- genotype_data[[1]]
  #   genotype_data <- genotype_data[-1]
  } else {
    stop("Unsupported file format. Please provide a csv, tsv, txt, or Excel file.")
  }
  
  return(genotype_data)
}

# Function to apply alignment to the tree plot if alignment file is specified
apply_alignment <- function(p, alignment_file, layout_option) {

  # Check if layout is "daylight" or "equal_angle"
  if (layout_option %in% c("daylight", "equal_angle")) {
    warning("Alignment cannot be applied with 'daylight' or 'equal_angle' layout. The function will not be executed.")
    return(p)  # Return the original plot without applying alignment
  }

  if (!is.null(alignment_file)) {
    tryCatch({
      p <- msaplot(p, fasta = alignment_file)
      # al_read <- readDNAStringSet(alignment_file, format = "fasta")
      # BrowseSeqs(al_read, highlight=0)
      
    }, error = function(e) {
      stop("Error applying alignment: Check if the IDs in the alignment file match the tree tip labels.")
    })
  }
  return(p)
}

# Function to create a gheatmap
# need debug if there are empty data
plot_genotype_heatmap <- function(tree_plot, genotype_file, alignment_file, offset = 4, width = 0.17, 
                          colnames_angle = -45, hjust = 0) {

  # Check if both alignment_file and genotype_file are provided or if neither is provided
  if ((!is.null(alignment_file) && !is.null(genotype_file))) {
    warning("The function cannot be executed if both an alignment file and a genotype file are provided, or if the genotype file is missing.")
    return(tree_plot)
  }

  if (!is.null(genotype_file)) {
    tryCatch({
      # Read the genotype file into phenotype_data (assuming it's a data frame)
      phenotype_data <- read_data_table(genotype_file)

      # Get the unique values ​​of the genotype
      unique_phenotypes <- unique(unlist(phenotype_data))
      # Get number of columns
      number_cols <- length(phenotype_data)
      # Create number width
      width <- number_cols * width
      
      # Create the color palette using RColorBrewer
      if (length(unique_phenotypes) <= 12) { 
          # Check if the number of unique phenotypes is less than or equal to 12
          colors <- brewer.pal(length(unique_phenotypes), "Set3")
      } else if (length(unique_phenotypes) <= 20) { 
          # If the number of unique phenotypes is greater than 12
          palette_tandem <- c(brewer.pal(12, "Set3"), brewer.pal(8, "Dark2"))
          colors <- palette_tandem[1:length(unique_phenotypes)]
      } else {
          # If the number of unique phenotypes exceeds 20
          stop("The number of unique genotypes exceeds the number of colors in the selected palette.")
      }
      
      # Generate the heat map
      p <- gheatmap(tree_plot, phenotype_data, offset = offset, width = width, 
                    colnames_angle = colnames_angle, hjust = hjust) +
        scale_fill_manual(breaks = unique_phenotypes, values = colors, name = "Phenotype")
    }, error = function(e) {
      stop("Error applying genotype plot on phylogenetics tree: Check if the IDs in the genotype file match the tree tip labels.")
    })
  }

  return(p)
}

# NEED ADVICED TO SCALE
# Function to create a continuous value heatmap
plot_continuous_heatmap <- function(tree_plot, contineous_file, alignment_file, offset = 6.5, width = 0.17, 
                          colnames_angle = -45, hjust = 0){
  
  # Check if both alignment_file and genotype_file are provided or if neither is provided
  if ((!is.null(alignment_file) && !is.null(contineous_file))) {
    warning("The function cannot be executed if both an alignment file and a contineous value file are provided, or if the contineous file is missing.")
    return(tree_plot)
  }

    if (!is.null(contineous_file)) {
    tryCatch({
      # Read the contineous data file (assuming it's a data frame)
      contineous_data <- read_data_table(contineous_file)
      # Get number of columns
      number_cols <- length(contineous_data)
      # Create number width
      width <- number_cols * width

      p <- gheatmap(tree_plot, contineous_data, offset = offset, width = width, colnames_angle = colnames_angle,
                        hjust = hjust, colnames_offset_y = 0.25) +  
        scale_fill_viridis_c(option="B", name="continuous\nvalue")

    }, error = function(e) {
      stop("Error applying contineous plot on phylogenetic: Check if the IDs in the contineous file match the tree tip labels.")
    })
  }

  return(p)

}

## NEED work on this function
check_both_heatmaps <- function(p, genotype_file, contineous_file, alignment_file){
  # Check if both contineous_file and genotype_file are provided
  if ((!is.null(genotype_file) && !is.null(contineous_file))) {
    newscale_package <- "ggnewscale"
    manage_packages(newscale_package)

    # Get the names of the passed arguments
    opt_names <- names(opt)

    # Identify the order of arguments
    genotype_index <- which(opt_names == "genotype")
    contineous_index <- which(opt_names == "countineous")
    # Make offset variable on order to put args
    # Check which file was specified first
    if (genotype_index < contineous_index) {
      # If --genotype was specified first
      p1 <- plot_genotype_heatmap(p, genotype_file, alignment_file)
      p2 <- p1 + new_scale_fill()
      p2 <- plot_continuous_heatmap(p2, contineous_file, alignment_file)
    } else if (genotype_index > contineous_index) {
      # If --contineous was specified first
      p1 <- plot_continuous_heatmap(p, contineous_file, alignment_file)
      p2 <- p1 + new_scale_fill()
      p2 <- plot_genotype_heatmap(p2, genotype_file, alignment_file)
    }
      # Return the last generated graph
      return(p2)


  } else if (!is.null(genotype_file)) {
    # If there is only genotype
    return(plot_genotype_heatmap(p, genotype_file, alignment_file))
  } else if (!is.null(contineous_file)) {
    # If there is only continuous variable
    return(plot_continuous_heatmap(p, contineous_file, alignment_file))
  } else {
    stop("Debe proporcionar al menos uno de los archivos: genotype o contineous.")
  }
}

save_phylogenetic_plot <- function(p, genotype_file = NULL, continuous_file = NULL, phy_file, alignment, pdf_width, pdf_height) {
  
  # Case 1: Both genotype and continuous files are provided
  if (!is.null(genotype_file) && !is.null(continuous_file)) {
  
    # 1. Phylogeny with genotype heatmap
    p_genotype <- plot_genotype_heatmap(p, genotype_file, alignment)
    output_genotype <- sub("\\..+$", "_genotype.pdf", phy_file)  # Create output file for genotype heatmap
    ggsave(output_genotype, plot = p_genotype, device = "pdf", width = pdf_width, height = pdf_height)  # Save plot
    
    # 2. Phylogeny with continuous heatmap
    p_continuous <- plot_continuous_heatmap(p, continuous_file, alignment)
    output_continuous <- sub("\\..+$", "_continuous.pdf", phy_file)  # Create output file for continuous heatmap
    ggsave(output_continuous, plot = p_continuous, device = "pdf", width = pdf_width, height = pdf_height)  # Save plot
    
    # 3. Phylogeny with both heatmaps
    p_both <- check_both_heatmaps(p, genotype_file, continuous_file, alignment)  # Generate phylogeny with both heatmaps
    output_both <- sub("\\..+$", "_both_heatmaps.pdf", phy_file)  # Create output file for both heatmaps
    ggsave(output_both, plot = p_both, device = "pdf", width = pdf_width, height = pdf_height)  # Save plot

    # Print the names of the files that were saved
    cat("The graphs have been saved to", output_genotype, output_continuous, "and", output_both, "\n")
    
  } else if (!is.null(genotype_file) || !is.null(continuous_file)) {
    
    # Case 2: Only one of genotype or continuous files is provided
    p <- check_both_heatmaps(p, genotype_file, continuous_file, alignment)  # Generate plot with the provided heatmap
    output_pdf <- sub("\\..+$", ".pdf", phy_file)  # Create output file for the provided heatmap
    ggsave(output_pdf, plot = p, device = "pdf", width = pdf_width, height = pdf_height)  # Save plot
    cat("The graph with the provided heatmap has been saved to", output_pdf, "\n")
    
  } else {
    
    # Case 3: Neither genotype nor continuous files are provided
    output_pdf <- sub("\\..+$", ".pdf", phy_file)  # Create output file for the basic phylogeny plot
    ggsave(output_pdf, plot = p, device = "pdf", width = pdf_width, height = pdf_height)  # Save the basic phylogeny plot
    cat("The base phylogeny graph has been saved to", output_pdf, "\n")
    
  }
}

# Set constant dimensions for the PDF
pdf_width <- 20  
pdf_height <- 20 

#############
# EXECUTION #
#############

# Convert cluster option to lowercase to handle case insensitivity
cluster_option <- tolower(opt$cluster)

# Handle clustering based on the option provided
if (cluster_option == "auto") {
  cat("Automatic clustering activated.\n")
  
  # Calculate distance matrix
  dist_matrix <- cophenetic(tree)
  data_matrix <- as.matrix(dist_matrix)
  
  # Perform hierarchical clustering
  hc <- hclust(as.dist(dist_matrix), method = "average")
  
  optclus <- sapply(2:9, function(x) summary(silhouette(cutree(hc, k = x), data_matrix))$avg.width)
  optnclust <- which(optclus == max(optclus))
  
  cat("Optimal number of clusters k=", optnclust, "\n")

  
  # Perform cut based on k
  clusters <- cutree(hc, k = optnclust)
  
  # Obtain leaf labels
  tip_labels <- tree$tip.label
  
  # Create a list of groups
  groups <- list()
  for (i in 1:optnclust) {
    group_name <- paste("Cluster", i)
    groups[[group_name]] <- tip_labels[clusters == i]
  }
  
  # Identify nodes for each clade
  nodes <- lapply(groups, function(tips) {
    MRCA(tree, tips)
  })
  
  num_colors <- length(groups)
  colors <- brewer.pal(num_colors, "Set3")
  
  p <- ggtree(tree, layout = opt$layout)
  
  # Add geom_hilight and geom_cladelabel based on nodes
  for (i in seq_along(nodes)) {
    p <- p + geom_hilight(node = nodes[[i]], fill = colors[i], alpha = 0.2) +
      geom_cladelabel(node = nodes[[i]], label = names(nodes)[i], 
                      color = colors[i], offset = .1, barsize = 2,
                      fontsize = 8, align = TRUE, alpha = 1)
  }
  
  p <- p %<+% data + 
    geom_tiplab(aes(color = I(color))) +  # Paint species labels
    geom_text2(aes(subset = !isTip, label = label), hjust = -.3) +  # Show bootstrap values
    geom_treescale() +
    theme(legend.position = "none")  # Hide color legend

  # Apply alignment if specified
  p <- apply_alignment(p, opt$alignment, opt$layout)

  # # Apply genotype heatmap plot on phylogenetic
  # p <- plot_genotype_heatmap(p, opt$genotype, opt$alignment)

  # # Apply contineous heatmap plot on phylogenetic
  # p <- plot_continuous_heatmap(p, opt$countineous, opt$alignment)

  # Save the graph to a PDF file with specified dimensions
  # output_pdf <- sub("\\..+$", ".pdf", opt$phy)
  # ggsave(output_pdf, plot = p, device = "pdf", width = pdf_width, height = pdf_height)

  save_phylogenetic_plot(
  p, 
  genotype_file = opt$genotype, 
  continuous_file = opt$countineous, 
  phy_file = opt$phy, 
  alignment = opt$alignment, 
  pdf_width = pdf_width, 
  pdf_height = pdf_height
  )
  
  # # Show success message
  # cat("The graph has been saved to", output_pdf, "\n")

  
} else if (grepl("^[0-9]+$", cluster_option)) {
  k <- as.numeric(cluster_option)
  cat("Automatic clustering activated with k =", k, "\n")
  
  # Calculate distance matrix
  dist_matrix <- cophenetic(tree)
  
  # Perform hierarchical clustering
  hc <- hclust(as.dist(dist_matrix), method = "average")
  
  # Cut the tree into k clusters
  clusters <- cutree(hc, k = k)
  
  # Obtain leaf labels
  tip_labels <- tree$tip.label
  
  # Create a list of groups
  groups <- list()
  for (i in 1:k) {
    group_name <- paste("Cluster", i)
    groups[[group_name]] <- tip_labels[clusters == i]
  }
  
  # Identify nodes for each clade
  nodes <- lapply(groups, function(tips) {
    MRCA(tree, tips)
  })
  
  num_colors <- length(groups)
  colors <- brewer.pal(num_colors, "Set3")
  
  p <- ggtree(tree, layout = opt$layout)
  
  # Add geom_hilight and geom_cladelabel based on nodes
  for (i in seq_along(nodes)) {
    p <- p + geom_hilight(node = nodes[[i]], fill = colors[i], alpha = 0.2) +
      geom_cladelabel(node = nodes[[i]], label = names(nodes)[i], 
                      color = colors[i], offset = .1, barsize = 2,
                      fontsize = 8, align = TRUE, alpha = 1)
  }
  
  p <- p %<+% data + 
    geom_tiplab(aes(color = I(color))) +  # Color the species labels
    geom_text2(aes(subset = !isTip, label = label), hjust = -.3) +  # Show bootstrap values
    geom_treescale() +
    theme(legend.position = "none")  # Hide color legend

  # Apply alignment if specified
  p <- apply_alignment(p, opt$alignment, opt$layout)

  # # Apply genotype heatmap plot on phylogenetic
  # p <- plot_genotype_heatmap(p, opt$genotype, opt$alignment)

  # # Apply contineous heatmap plot on phylogenetic
  # p <- plot_continuous_heatmap(p, opt$countineous, opt$alignment)

# Save the graph to a PDF file with specified dimensions
  # output_pdf <- sub("\\..+$", ".pdf", opt$phy)
  # ggsave(output_pdf, plot = p, device = "pdf", width = pdf_width, height = pdf_height)

  save_phylogenetic_plot(
  p, 
  genotype_file = opt$genotype, 
  continuous_file = opt$countineous, 
  phy_file = opt$phy, 
  alignment = opt$alignment, 
  pdf_width = pdf_width, 
  pdf_height = pdf_height
  )
  
  # # Show success message
  # cat("The graph has been saved to", output_pdf, "\n")

} else {
  cat("Automatic clustering disabled.\n")
  
  # Create ggtree object with specified layout
  p <- ggtree(tree, layout = opt$layout)
  
  # Add colors to tree branches
  p <- p %<+% data + 
    geom_tiplab(aes(color = I(color))) +  # Paint species labels
    geom_text2(aes(subset = !isTip, label = label), hjust = -.3) + # Show bootstrap values
    geom_treescale() + 
    theme(legend.position = "none")  # Hide color legend

  # Apply alignment if specified
  p <- apply_alignment(p, opt$alignment, opt$layout)

  save_phylogenetic_plot(
  p, 
  genotype_file = opt$genotype, 
  continuous_file = opt$countineous, 
  phy_file = opt$phy, 
  alignment = opt$alignment, 
  pdf_width = pdf_width, 
  pdf_height = pdf_height
  )
  
  # # Show success message
  # cat("The graph has been saved to", output_pdf, "\n")
}

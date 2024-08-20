#!/usr/bin/env Rscript

# Load functions from functions.R
source("quick_off_init.R")

required_packages <- c("ggtree", "treeio", "ape", "RColorBrewer", "optparse", "cluster", "factoextra") # , "DECIPHER"
manage_packages(required_packages)

# Define command line options
option_list <- list(
  make_option(c("--phy"), type = "character", default = NULL, help = "Archivo de filogenia", metavar = "FILE"),
  make_option(c("--keyword"), type = "character", default = NULL, help = "Palabras para colorear separadas por espacios", metavar = "WORDS"),
  make_option(c("--alignment"), type = "character", default = NULL, help = "Archivo de alineamiento", metavar = "FILE"),
  make_option(c("--root"), type = "character", default = NULL, help = "Raíz del árbol", metavar = "NODE"),
  make_option(c("--layout"), type = "character", default = "rectangular", help = "Tipo de layout del árbol (rectangular, equal_angle, daylight, circular, roundrect)", metavar = "LAYOUT"),
   # make_option(c("--cluster"), type = "logical", default = FALSE, help = "Activar clustering automático", action = "store_true")
  make_option(c("--genotype"), type = "character", default = NULL, help = "Archivo de genotipo", metavar = "FILE"),
  make_option(c("--cluster"), type = "character", default = "no", help = "Clustering automático: AUTO, número para especificar k, o NO para desactivar", metavar = "CLUSTER")
)

############
# FUNCTION #
############

# Parse arguments
opt <- parse_args(OptionParser(option_list = option_list))

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
plot_gheatmap <- function(tree_plot, genotype_file, alignment_file, offset = 5, width = 0.5, font_size = 3, 
                          colnames_angle = -45, hjust = 0, color_palette = "Set3") {
  # Check if alignment file
  if (alignment_file) {
    warning("Genotype cannot be applied with alignment file. The function will not be executed.")
    return(p)
  }

  if (!is.null(genotype_file)) {
    tryCatch({                          
      # Get the unique values ​​of the genotype
      unique_phenotypes <- unique(unlist(phenotype_data))
      
      # Create the color palette using RColorBrewer
      colors <- brewer.pal(length(unique_phenotypes), color_palette)
      
      # Generate the heat map
      p <- gheatmap(tree_plot, phenotype_data, offset = offset, width = width, font.size = font_size, 
                    colnames_angle = colnames_angle, hjust = hjust) +
        scale_fill_manual(breaks = unique_phenotypes, values = colors, name = "Phenotype")
    }, error = function(e) {
      stop("Error applying alignment: Check if the IDs in the genotype file match the tree tip labels.")
    })
  }  
  return(p)
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
    theme(legend.position = "none")  # Hide color legend

# Apply alignment if specified
p <- apply_alignment(p, opt$alignment, opt$layout)

# Apply genotype heatmap plot on phylogenetic
p <- plot_gheatmap(p, opt$genotype, opt$alignment)

  # Save the graph to a PDF file with specified dimensions
  output_pdf <- sub("\\..+$", ".pdf", opt$phy)
  ggsave(output_pdf, plot = p, device = "pdf", width = pdf_width, height = pdf_height)
  
  # Show success message
  cat("The graph has been saved to", output_pdf, "\n")

  
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
    theme(legend.position = "none")  # Hide color legend

  # Apply alignment if specified
  p <- apply_alignment(p, opt$alignment, opt$layout)

  # Apply genotype heatmap plot on phylogenetic
p <- plot_gheatmap(p, opt$genotype, opt$alignment)

# Save the graph to a PDF file with specified dimensions
  output_pdf <- sub("\\..+$", ".pdf", opt$phy)
  ggsave(output_pdf, plot = p, device = "pdf", width = pdf_width, height = pdf_height)
  
  # Show success message
  cat("The graph has been saved to", output_pdf, "\n")

} else {
  cat("Automatic clustering disabled.\n")
  
  # Create ggtree object with specified layout
  p <- ggtree(tree, layout = opt$layout)
  
  # Add colors to tree branches
  p <- p %<+% data + 
    geom_tiplab(aes(color = I(color))) +  # Paint species labels
    geom_text2(aes(subset = !isTip, label = label), hjust = -.3) + # Show bootstrap values
    theme(legend.position = "none")  # Hide color legend

  # Apply alignment if specified
  p <- apply_alignment(p, opt$alignment, opt$layout)

  # Apply genotype heatmap plot on phylogenetic
  p <- plot_gheatmap(p, opt$genotype, opt$alignment)
  
  # Save the graph to a PDF file with specified dimensions
  output_pdf <- sub("\\..+$", ".pdf", opt$phy)
  ggsave(output_pdf, plot = p, device = "pdf", width = pdf_width, height = pdf_height)
  
  # Show success message
  cat("The graph has been saved to", output_pdf, "\n")
}


# # Example of using the --cluster option
# if (opt$cluster) {
#   cat("Automatic clustering activated.\n")
  
#   # Assuming you already have the distance matrix calculated
#   dist_matrix <- cophenetic(tree)
#   data_matrix <- as.matrix(dist_matrix)
  
#   # Perform hierarchical clustering
#   hc <- hclust(as.dist(dist_matrix), method = "average")
  
#   optclus <- sapply(2:9, function(x) summary(silhouette(cutree(hc, k = x), data_matrix))$avg.width)
#   optnclust <- which(optclus == max(optclus))
  
#   print(optnclust)
  
#   # Perform cut based on k
#   clusters <- cutree(hc, k = optnclust)
  
#   # Obtain leaf labels
#   tip_labels <- tree$tip.label
  
#   # Create a list of groups
#   groups <- list()
#   for (i in 1:optnclust) {
#     group_name <- paste("Cluster", i)
#     groups[[group_name]] <- tip_labels[clusters == i]
#   }
  
#   # Identify nodes for each clade
#   nodes <- lapply(groups, function(tips) {
#     MRCA(tree, tips)
#   })
  
#   num_colors <- length(groups)
#   colors <- brewer.pal(num_colors, "Set3")
  
#   p <- ggtree(tree, layout = opt$layout)
  
#   # Add geom_hilight and geom_cladelabel based on nodes
#   for (i in seq_along(nodes)) {
#     p <- p + geom_hilight(node = nodes[[i]], fill = colors[i], alpha = 0.2) +
#       geom_cladelabel(node = nodes[[i]], label = names(nodes)[i], 
#                       color = colors[i], offset = .1, barsize = 2,
#                       fontsize = 8, align = TRUE, alpha = 0.2)
#   }
  
#   p <- p %<+% data + 
#     geom_tiplab(aes(color = I(color))) +  # Paint species labels
#     geom_text2(aes(subset = !isTip, label = label), hjust = -.3) +  # Show bootstrap values
#     theme(legend.position = "none")  # Hide color legend

#   # Save the graph to a PDF file with specified dimensions
#   output_pdf <- sub("\\..+$", ".pdf", opt$phy)
#   ggsave(output_pdf, plot = p, device = "pdf", width = pdf_width, height = pdf_height)
  
#   # Show success message
#   cat("The graph has been saved to", output_pdf, "\n")
  
# } else {
#   cat("Automatic clustering disabled.\n")
  
#   # Create ggtree object with specified layout
#   p <- ggtree(tree, layout = opt$layout)
  
#   # Add colors to tree branches
#   p <- p %<+% data + 
#     geom_tiplab(aes(color = I(color))) +  # Paint species labels
#     geom_text2(aes(subset = !isTip, label = label), hjust = -.3) + # Show bootstrap values
#     theme(legend.position = "none")  # Hide color legend
  
#   # Save the graph to a PDF file with specified dimensions
#   output_pdf <- sub("\\..+$", ".pdf", opt$phy)
#   ggsave(output_pdf, plot = p, device = "pdf", width = pdf_width, height = pdf_height)
  
#   # Show success message
#   cat("The graph has been saved to", output_pdf, "\n")
# }




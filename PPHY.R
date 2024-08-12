#!/usr/bin/env Rscript

# Cargar paquetes necesarios
library(ggtree)
library(treeio)
library(ape)
library(RColorBrewer)
library(optparse)

## culster

library(cluster)
library(factoextra)


# Definir opciones de línea de comandos
option_list <- list(
  make_option(c("--phy"), type = "character", default = NULL, help = "Archivo de filogenia", metavar = "FILE"),
  make_option(c("--keyword"), type = "character", default = NULL, help = "Palabras para colorear separadas por espacios", metavar = "WORDS"),
  make_option(c("--alignment"), type = "character", default = NULL, help = "Archivo de alineamiento", metavar = "FILE"),
  make_option(c("--root"), type = "character", default = NULL, help = "Raíz del árbol", metavar = "NODE"),
  make_option(c("--layout"), type = "character", default = "rectangular", help = "Tipo de layout del árbol (rectangular, equal_angle, daylight)", metavar = "LAYOUT"),
  make_option(c("--width"), type = "numeric", default = 10, help = "Ancho del PDF en pulgadas", metavar = "WIDTH"),
  make_option(c("--height"), type = "numeric", default = 10, help = "Alto del PDF en pulgadas", metavar = "HEIGHT"),
  make_option(c("--cluster"), type = "logical", default = FALSE, help = "Activar clustering automático", action = "store_true")
)

# Parsear argumentos
opt <- parse_args(OptionParser(option_list = option_list))

# Verificar argumentos requeridos
if (is.null(opt$phy)) {
  stop("Debe especificar un archivo de filogenia con --phy")
}

# Leer árbol filogenético
tree <- read.tree(opt$phy)

# Verificar si el árbol contiene etiquetas
if (length(tree$tip.label) == 0) {
  stop("El árbol no contiene etiquetas de tips.")
}

# Crear un dataframe con la información de los nombres y los colores
data <- data.frame(label = tree$tip.label)
data$color <- "grey"  # Color por defecto

# Generar una paleta de colores para las palabras
if (!is.null(opt$keyword)) {
  palabras <- unlist(strsplit(opt$keyword, " "))
  n_palabras <- length(palabras)
  colores <- brewer.pal(min(n_palabras, 12), "Dark2")  # Limitar al tamaño máximo de la paleta
  
  for (i in 1:n_palabras) {
    palabra <- palabras[i]
    data$color[grepl(palabra, data$label)] <- colores[i]
  }
}

# Root the tree si se especifica
if (!is.null(opt$root)) {
  tree <- root(tree, outgroup = opt$root)
}


# Ejemplo de uso de la opción --cluster
if (opt$cluster) {
  cat("Clustering automático activado.\n")
  
  # Suponiendo que ya tienes la matriz de distancias calculada
  dist_matrix <- cophenetic(tree)
  data_matrix <- as.matrix(dist_matrix)
  
  # Realizar agrupamiento jerárquico
  hc <- hclust(as.dist(dist_matrix), method = "average")
  
  optclus <- sapply(2:9, function(x) summary(silhouette(cutree(hc, k = x), data_matrix))$avg.width)
  optnclust <- which(optclus == max(optclus))
  
  print(optnclust)
  
  # Realizar corte basado en k
  clusters <- cutree(hc, k = optnclust)
  
  # Obtener etiquetas de las hojas
  tip_labels <- tree$tip.label
  
  # Crear lista de grupos
  groups <- list()
  for (i in 1:optnclust) {
    group_name <- paste("Cluster", i)
    groups[[group_name]] <- tip_labels[clusters == i]
  }
  
  # Identificar nodos para cada clado
  nodes <- lapply(groups, function(tips) {
    MRCA(tree, tips)
  })
  
  num_colors <- length(groups)
  colors <- brewer.pal(num_colors, "Set3")
  
  p <- ggtree(tree, layout = opt$layout)
  
  # Añadir geom_hilight y geom_cladelabel basado en los nodos
  for (i in seq_along(nodes)) {
    p <- p + geom_hilight(node = nodes[[i]], fill = colors[i], alpha = 0.2) +
      geom_cladelabel(node = nodes[[i]], label = names(nodes)[i], 
                      color = colors[i], offset = .1, barsize = 2,
                      fontsize = 5, align = TRUE, alpha = 0.5)
  }
  
  p <- p %<+% data + 
    geom_tiplab(aes(color = I(color))) +  # Pintar las etiquetas de las especies
    geom_text2(aes(subset = !isTip, label = label), hjust = -.3) +  # Mostrar los valores de bootstrap
    theme(legend.position = "none")  # Ocultar leyenda de colores
  # Guardar el gráfico en un archivo PDF con dimensiones especificadas
  output_pdf <- sub("\\..+$", ".pdf", opt$phy)
  ggsave(output_pdf, plot = p, device = "pdf", width = opt$width*3, height = opt$height*3)
  
  # Mostrar mensaje de éxito
  cat("El gráfico se ha guardado en", output_pdf, "\n")
  
} else {
  cat("Clustering automático desactivado.\n")
  
  # Crear objeto ggtree con layout especificado
  p <- ggtree(tree, layout = opt$layout)
  
  # Añadir colores a las ramas del árbol
  p <- p %<+% data + 
    geom_tiplab(aes(color = I(color))) +  # Pintar las etiquetas de las especies
    geom_text2(aes(subset = !isTip, label = label), hjust = -.3) +  # Mostrar los valores de bootstrap
    theme(legend.position = "none")  # Ocultar leyenda de colores
  
  # Guardar el gráfico en un archivo PDF con dimensiones especificadas
  output_pdf <- sub("\\..+$", ".pdf", opt$phy)
  ggsave(output_pdf, plot = p, device = "pdf", width = opt$width*3, height = opt$height*3)
  
  # Mostrar mensaje de éxito
  cat("El gráfico se ha guardado en", output_pdf, "\n")
}


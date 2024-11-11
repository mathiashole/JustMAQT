
# Function to check if a given input is a list of numeric values
is_numeric_list <- function(input) {
  # Split the input by spaces
  elements <- strsplit(input, " ")[[1]]
  
  # Check if there is more than one element and all elements are numeric
  length(elements) > 1 && all(sapply(elements, function(x) !is.na(as.numeric(x))))
}

# Function to manage packages, installing if not available
manage_packages <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}
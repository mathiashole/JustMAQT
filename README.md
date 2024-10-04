# :tanabata_tree: JustMAQT

[![R](https://img.shields.io/badge/R-276DC3?style=for-the-badge&logo=R&logoColor=white&labelColor=101010)](https://www.r-project.org/about.html)
![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/mathiashole/JustMAQT?style=for-the-badge&labelColor=101010&color=white)
![GitHub](https://img.shields.io/github/license/mathiashole/GScissors?color=%23179287&style=for-the-badge&logoColor=white&labelColor=101010)

`JustMAQT` is an R script for customizing and visualizing phylogenetic trees.

## :book: Features

-   Plotting phylogenetic trees using input from a variety of file formats.
-   Customizable tree layouts (rectangular, circular, equal_angle, etc.).
-   Color-coding tree branches based on specified keywords.
-   Option to root the tree at a specified ID.
-   Clustering options 
-   Heatmaps option for continuous values or genotypes.
-   Command-line interface.

## :wrench: Usage

To use the phylogenetic tree plotting script, follow these steps:

1. Download or clone the script from the repository.
2. Ensure the required R packages are installed: `ggtree`, `treeio`, `ape`, `RColorBrewer`, `optparse`, `cluster`, `factoextra`, and `ggplot2`.
3. **Option 1: Using the terminal**  
   Open a terminal, navigate to the folder containing the script, and run the script with the appropriate arguments.
4. **Option 2: Using RStudio**  
   Open the script directly in RStudio and run it from there, passing the arguments manually or setting them as variables within the script.

## :hammer: in progress ...

## :bulb: R Quick Examples

```{r, eval = FALSE}
args <- c("PPHY.R",
          "--phy", "/path/phylogeny_file",
          "--keyword", shQuote("species1 species2 species4"),
          "--layout", "roundrect",
          "--genotype", "/path/genotype_data_file",
          "--countineous", "/path/contineous_data_file")

# Execution command
system2("Rscript", args = args)
```

## :sparkling_heart: Contributing

- :octocat: [Pull requests](https://github.com/mathiashole/JustMAQT/pulls) and :star2: stars are always welcome.
- For major changes, please open an [issue](https://github.com/mathiashole/JustMAQT/issues) first to discuss what you would like to change.
- Please make sure to update tests as appropriate.

## :mega: Contact

:bird: [Mathias](https://twitter.com/joaquinmangino)

## License
MIT &copy; [Mathias Mangino](https://github.com/mathiashole)
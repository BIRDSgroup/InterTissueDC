# Set working directory ---------------------------------------------------

working_directory <- "/data/kailash/DDP"
setwd(working_directory)
.libPaths(c("./libraries", .libPaths()))
.libPaths()

# Loading libraries -------------------------------------------------------


if (!require("pacman")) install.packages("pacman")

# Package names
packages <- c("remotes", "ggplot2", "dplyr", "BiocManager", "readxl", "limma", "edgeR",
              "BiocGenerics", "Biobase", "S4Vectors", "IRanges", "cli", "AnnotationDbi",
              "annotate", "genefilter", "BiocParallel", "sva", "CellCODE", "BRETIGEA",
              "Matrix", "ggdendro", "pheatmap", "viridis", "RColorBrewer",
              "psych", "ggpubr", "GGally", "ggpmisc", "ggfortify", "corrplot", "reshape2",
              "ellipsis", "lubridate", "pillar", "S4Vectors", "splus2R", "vctrs", "VennDiagram",
              "xts", "purrr", "stringr")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages], INSTALL_opts = "--no-multiarch",
                   dependencies = TRUE)
}

remotes::install_github("DesiQuintans/librarian", INSTALL_opts = "--no-multiarch")
remotes::install_github("mchikina/CellCODE", INSTALL_opts = "--no-multiarch")

BiocPackages <- c("limma", "edgeR", "BiocGenerics", "Biobase", "S4Vectors",
                  "IRanges", "AnnotationDbi", "annotate", "genefilter",
                  "BiocParallel", "sva", "GSALightning")
BioC_installed_packages <- BiocPackages %in% rownames(installed.packages())
if (any(BioC_installed_packages == FALSE)) {
  BiocManager::install(BiocPackages[!BioC_installed_packages], INSTALL_opts = "--no-multiarch",
                       dependencies = TRUE, ask = FALSE, update = FALSE)
}

# ipak function: install and load multiple R packages.
# check to see if packages are installed. Install them if they are not,
# then load them into the R session.

# Source: https://stackoverflow.com/a/24183297/8281802

ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE, INSTALL_opts = "--no-multiarch")
  sapply(pkg, require, character.only = TRUE)
}

# usage
packages <- c(packages, BiocPackages)
ipak(packages)

# Set working directory ---------------------------------------------------

.libPaths(c("libraries/", .libPaths()))
.libPaths()

libDirectory <- "libraries/"

# Loading libraries -------------------------------------------------------


if (!require("pacman")){install.packages("pacman", lib = libDirectory, repos = "http://cran.us.r-project.org")} 

# Package names
packages <- c("tidyverse", "remotes", "ggplot2", "dplyr", "BiocManager", "readxl", "limma", "edgeR",
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
                   dependencies = TRUE, lib = libDirectory, repos = "http://cran.us.r-project.org")
}

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager", repos = "http://cran.us.r-project.org")

# BiocManager::install("sva", lib = libDirectory)

if (!require("librarian")){remotes::install_github("DesiQuintans/librarian", INSTALL_opts = "--no-multiarch", lib = libDirectory)}
if (!require("CellCODE")){remotes::install_github("mchikina/CellCODE", INSTALL_opts = "--no-multiarch", lib = libDirectory)}

BiocPackages <- c("limma", "edgeR", "BiocGenerics", "Biobase", "S4Vectors",
                  "IRanges", "AnnotationDbi", "annotate", "genefilter",
                  "BiocParallel", "sva", "GSALightning")
BioC_installed_packages <- BiocPackages %in% rownames(installed.packages())
if (any(BioC_installed_packages == FALSE)) {
  BiocManager::install(BiocPackages[!BioC_installed_packages], INSTALL_opts = "--no-multiarch",
                       dependencies = TRUE, ask = FALSE, update = FALSE, lib = libDirectory)
}

# ipak function: install and load multiple R packages.
# check to see if packages are installed. Install them if they are not,
# then load them into the R session.

# Source: https://stackoverflow.com/a/24183297/8281802

ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE, INSTALL_opts = "--no-multiarch", lib = libDirectory, repos = "http://cran.us.r-project.org")
  sapply(pkg, require, character.only = TRUE)
}

# usage
packages <- c(packages, BiocPackages)
ipak(packages)


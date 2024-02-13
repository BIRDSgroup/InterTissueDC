# Making the scatter plot of pVals before and after cell type correction

rm(list = ls())

# Loading libraries -------------------------------------------------------

.libPaths(c("libraries/", .libPaths()))
source("scripts/libraries.R")
source("/utils.R")

library(ggplot2)
library(tidyverse)
library(ggpubr)

# Looking into DC plots ---------------------------------------------------

make_looking_into_dc_plots <- function(BR_before_ctc, BR_after_ctc){
  before_ctc <- read.csv(paste0("results/before_ctc_between_region/", BR_before_ctc[1], BR_before_ctc[2], "Mapped_DCEdgesWithoutNA_fdr_0.01.csv"),
                    stringsAsFactors = FALSE)
  after_ctc <- read.csv(paste0("results/after_ctc_between_region/", BR_after_ctc[1], BR_after_ctc[2], "_ct_correctionMapped_DCEdgesWithoutNA_fdr_0.01.csv"),
                      stringsAsFactors = FALSE)
  
  # Filtering only common DC-edges between the comparisons
  
  before_ctc$dc_edge <- paste(before_ctc$a, before_ctc$b, sep = "-")
  after_ctc$dc_edge <- paste(after_ctc$a, after_ctc$b, sep = "-")
  common_dc_edges <- intersect(before_ctc$dc_edge, after_ctc$dc_edge)
  
  message("Number of DC edges before CTC ", dim(before_ctc)[1])
  message("Number of DC edges after CTC ", dim(after_ctc)[1])
  message("Number of common DC edges ", length(common_dc_edges))
  
  before_ctc_mm <- before_ctc %>% filter(m1 == "M", m2 == "M") %>% filter(dc_edge %in% common_dc_edges)
  before_ctc_nmnm <- before_ctc %>% filter(m1 == "NM", m2 == "NM") %>% filter(dc_edge %in% common_dc_edges)
  after_ctc_mm <- after_ctc %>% filter(m1 == "M", m2 == "M") %>% filter(dc_edge %in% common_dc_edges)
  after_ctc_nmnm <- after_ctc %>% filter(m1 == "NM", m2 == "NM") %>% filter(dc_edge %in% common_dc_edges)
  
  rm(before_ctc, after_ctc)
  
  p1 <- ggplot(before_ctc_mm, aes(x = z)) +
    ggtitle(paste0("M-M edges before cell-type correction ", "(n=", dim(before_ctc_mm)[1], ")")) +
    geom_histogram( colour = "black", fill = "blue") +
    theme_bw(base_size = 11) +
    theme(axis.text = element_text(size = 10, face = "bold"),
          axis.title = element_text(size = 10,face = "bold"),
          plot.title = element_text(hjust = 0.5))
  
  
  p2 <- ggplot(before_ctc_nmnm, aes(x = z)) + 
    geom_histogram(colour = "black", fill = "blue") +
    ggtitle(paste0("NM-NM edges before cell-type correction ", "(n=", dim(before_ctc_nmnm)[1], ")")) +
    theme_bw(base_size = 11) +
    theme(axis.text = element_text(size = 10, face = "bold"),
          axis.title = element_text(size = 10,face = "bold"),
          plot.title = element_text(hjust = 0.5))
  
  p3 <- ggplot(after_ctc_mm, aes(x = z)) +
    ggtitle(paste0("M-M edges after cell-type correction ", "(n=", dim(after_ctc_mm)[1], ")")) + 
    geom_histogram(colour = "black", fill = "blue") +
    theme_bw(base_size = 11) +
    theme(axis.text = element_text(size = 10, face = "bold"),
          axis.title = element_text(size = 10,face = "bold"),
          plot.title = element_text(hjust = 0.5))
  
  p4 <- ggplot(after_ctc_nmnm, aes(x = z)) +
    ggtitle(paste0("NM-NM edges after cell-type correction", "(n=", dim(after_ctc_nmnm)[1], ")")) + 
    geom_histogram(colour = "black", fill = "blue") +
    theme_bw(base_size = 11) +
    theme(axis.text = element_text(size = 10, face = "bold"),
          axis.title = element_text(size = 10,face = "bold"),
          plot.title = element_text(hjust = 0.5))
  
  figure <- ggarrange(p1, p3, p2, p4, ncol = 2, nrow = 2)
  ggsave(path = "plots/",
         device = "pdf",
         dpi = 300,
         width = 18,
         height = 10,
         filename = paste0("looking_into_dc_", BR_after_ctc[1], BR_after_ctc[2], "common_edges", ".pdf"),
         plot = annotate_figure(figure,
                                top = text_grob(paste("Looking into ", BR_after_ctc[1], "-", BR_after_ctc[2], " DC edges", sep = ""),
                                                face = "bold", size = 14)))
}

make_looking_into_dc_plots(BR_before_ctc <- c("BM_10", "BM_22"),
                           BR_after_ctc <- c("BM10", "BM22"))

make_looking_into_dc_plots(BR_before_ctc <- c("BM_10", "BM_36"),
                           BR_after_ctc <- c("BM10", "BM36"))

make_looking_into_dc_plots(BR_before_ctc <- c("BM_10", "BM_44"),
                           BR_after_ctc <- c("BM10", "BM44"))

make_looking_into_dc_plots(BR_before_ctc <- c("BM_22", "BM_36"),
                           BR_after_ctc <- c("BM22", "BM36"))

make_looking_into_dc_plots(BR_before_ctc <- c("BM_22", "BM_44"),
                           BR_after_ctc <- c("BM22", "BM44"))

make_looking_into_dc_plots(BR_before_ctc <- c("BM_36", "BM_44"),
                           BR_after_ctc <- c("BM36", "BM44"))

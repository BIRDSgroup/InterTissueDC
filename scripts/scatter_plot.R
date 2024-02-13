# Kailash, can you please add the example for a gene pair based on Fig 1b of Cell Systems paper.
# Can you do it for 2 gene pairs- ZNF81-ZKSCAN1 & SH2D5-ZKSCAN1, both from PHG_IFG (BM36-44)
# Use cell-type corrected data

# Braak based colouring

rm(list = ls())

# Loading necessary packages ----------------------------------------------

.libPaths(c("libraries/", .libPaths()))
source("scripts/libraries.R")
source("scripts/utils.R")

library(tidyverse)
library(ggplot2)
library(ggpubr)

# Code --------------------------------------------------------------------

metadata <- read.csv("data/rna_metadata_onlyBAM_noExclude_noDup_Nov2018.csv", stringsAsFactors = FALSE)
mapping_file <- read.csv("data/Mapped_h37_ens_to_gene.csv", stringsAsFactors = FALSE)

# Input to be specified ---------------------------------------------------

BR <- c("BM36", "BM44")
gene_pair <- c("ZNF81", "ZKSCAN1")
ct_corrected <- TRUE

# Prepping the data -------------------------------------------------------

if(ct_corrected == TRUE){
  ct_file_name <- "_cellTypeCorrected"
} else {
  ct_file_name = ""
}

ad1 <- read.csv(paste("results/AD", BR[1], "_", BR[1], BR[2], ct_file_name, ".csv", sep = ""), stringsAsFactors = FALSE, row.names = 1)
ad2 <- read.csv(paste("results/AD", BR[2], "_", BR[1], BR[2], ct_file_name, ".csv", sep = ""), stringsAsFactors = FALSE, row.names = 1)
ctl1 <- read.csv(paste("results/CTL", BR[1], "_", BR[1], BR[2], ct_file_name, ".csv", sep = ""), stringsAsFactors = FALSE, row.names = 1)
ctl2 <- read.csv(paste("results/CTL", BR[2], "_", BR[1], BR[2], ct_file_name, ".csv", sep = ""), stringsAsFactors = FALSE, row.names = 1)

ad1 <- data.frame(ad1)
ad2 <- data.frame(ad2)
ctl1 <- data.frame(ctl1)
ctl2 <- data.frame(ctl2)

ens_ids <- ensid_sampleid[,1]
colnames(ad1) <- colnames(ad2) <- colnames(ctl1) <- colnames(ctl2) <- ens_ids

# ad_ids <- read.csv(paste("results/AD_ids", BR[1], BR[2], ".csv", sep = ""), stringsAsFactors = FALSE)
# ctl_ids <- read.csv(paste("results/CTL_ids", BR[1], BR[2], ".csv", sep = ""), stringsAsFactors = FALSE)

ens <- unlist(lapply(gene_pair, function(s) mapping_file$gene_id[match(s, mapping_file$gene_name)]))

gene_ctl <- cbind(ctl1[, ens[1], drop = FALSE], ctl2[, ens[2], drop = FALSE])
# gene_ctl <- cbind(gene_ctl, ctl_ids$braak_score)
gene_ad <- cbind(ad1[, ens[1], drop = FALSE], ad2[, ens[2], drop = FALSE])
# gene_ad <- cbind(gene_ad, ad_ids$braak_score)

colnames(gene_ad) = colnames(gene_ctl) =  c("g1", "g2")

head(gene_ad)
head(gene_ctl)

# Making the plots --------------------------------------------------------

eps <- 0.025
ylim_min <- min(min(gene_ad$g2), min(gene_ctl$g2)) - eps
ylim_max <- max(max(gene_ad$g2), max(gene_ctl$g2)) + eps

p1 <- ggplot(gene_ctl, aes(g1, g2)) +
  geom_point(size = 2) + theme_bw(base_size = 20) +
  theme(axis.text = element_text(size = 11, face = "bold"),
        axis.title = element_text(size = 11,face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5)) +
  ggtitle("Control") +
  geom_smooth(method = 'lm', colour = "red") +
  stat_cor(method = "spearman") +
  labs(x = gene_pair[1], y = gene_pair[2]) +
  geom_point(data = data.frame(x = mean(gene_ctl %>% pull(g1)),
                               y = mean(gene_ctl %>% pull(g2)))
             , aes(x = x, y = y), shape=23, fill="blue", color="darkred", size=5) +
  ylim(ylim_min, ylim_max)

p2 <- ggplot(gene_ad, aes(g1, g2)) +
  geom_point(size = 2) + theme_bw(base_size = 20) +
  theme(axis.text = element_text(size = 11, face = "bold"),
        axis.title = element_text(size = 11,face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5)) +
  ggtitle("AD") +
  geom_smooth(method = 'lm', colour = "red") +
  stat_cor(method = "spearman") +
  labs(x = gene_pair[1], y = gene_pair[2]) +
  geom_point(data = data.frame(x = mean(gene_ad %>% pull(g1)),
                               y = mean(gene_ad %>% pull(g2)))
             , aes(x = x, y = y), shape=23, fill="blue", color="darkred", size=5) +
  ylim(ylim_min, ylim_max)

figure <- ggarrange(p1, p2, ncol = 2)

ggsave(path = "plots/",
       device = "pdf",
       dpi = 300,
       width = 9,
       height = 5,
       filename = paste("scatter_plot", BR[1], BR[2], gene_pair[1], gene_pair[2], ct_file_name, ".pdf", sep = "_"),
       plot = annotate_figure(figure,
                              top = text_grob(paste(BR[1], "-", BR[2], ct_file_name, sep = ""),
                                              face = "bold", size = 16)))

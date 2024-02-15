# Making the scatter plot of pVals before and after cell type correction

rm(list = ls())

# Loading libraries -------------------------------------------------------

.libPaths(c("libraries/", .libPaths()))
source("scripts/libraries.R")
source("scripts/utils.R")

library(ggplot2)
library(tidyverse)
library(ggpubr)

# Looking into DC plots ---------------------------------------------------

make_looking_into_dc_plots <- function(BR_before_ctc, BR_after_ctc){
  before_ctc <- read.csv(paste0("results/before_ctc_between_region/", BR_before_ctc[1], BR_before_ctc[2], "Mapped_DCEdgesWithoutNA_fdr_0.01.csv"),
                         stringsAsFactors = FALSE)
  
  analysis_type <- "between_region/"
  
  # Variables loaded - rna_metadata_processed, ensid_sampleid, clinical, genes_mapped
  
  read_cell_type_corrected_gExp <- function(br, index, condition){
    if(BR_after_ctc[1] != BR_after_ctc[2]){
      df <- read.csv(paste("results/", condition, BR_after_ctc[index], "_", BR_after_ctc[1], BR_after_ctc[2], "_cellTypeCorrected.csv",
                           sep = ""), stringsAsFactors = FALSE, row.names = 1)
    } else if(BR_after_ctc[1] == BR_after_ctc[2]){
      df <- read.csv(paste("results/", BR_after_ctc[1], "_AD_CTL_cellTypeCorrected.csv", sep = ""),
                     stringsAsFactors = FALSE, row.names = 1)
      barcodes <- na.omit(sapply(strsplit(row.names(df),"_"), function(x){x[3]}))
      ind_id <- rna_metadata_processed$individualIdentifier[match(barcodes, rna_metadata_processed$barcode)]
      # Noticed that some barcodes do not have indId mapping
      ind_id[which(is.na(ind_id))] <- 0
      
      sample_group <- clinical$DxCondition[match(ind_id, clinical$individualIdentifier)]
      
      if(condition == "AD"){
        df <- df[row.names(df)[which(sample_group == "Definite AD")],]
      }else if(condition == "CTL"){
        df <- df[row.names(df)[which(sample_group == "Control")],]
      }
    }
    return(as.matrix(df))
  }
  
  ens_genes <- ensid_sampleid[,1]
  AD1 <- read_cell_type_corrected_gExp(br = BR_after_ctc, index = 1, condition = "AD")
  AD2 <- read_cell_type_corrected_gExp(br = BR_after_ctc, index = 2, condition = "AD")
  CTL1 <- read_cell_type_corrected_gExp(br = BR_after_ctc, index = 1, condition = "CTL")
  CTL2 <- read_cell_type_corrected_gExp(br = BR_after_ctc, index = 2, condition = "CTL")
  
  colnames(AD1) <- colnames(AD2) <- colnames(CTL1) <- colnames(CTL2) <- ens_genes
  
  # BR is the normalized gene expression values which is corrected for batch effects
  # It is a matrix with dimension (#genes X #samples)
  
  # Finding number of AD and control samples in BR1 and BR2
  n.ad <- dim(AD1)[1]
  n.ctl <- dim(CTL1)[1]
  
  union12 <- paste(before_ctc$a, before_ctc$b)
  edges <- strsplit(union12, " ")
  mat  <- matrix(unlist(edges), ncol = 2, byrow = TRUE)
  df   <- as.data.frame(mat, stringsAsFactors = FALSE)
  colnames(df) <- c("a", "b")
  
  edges <- df
  edges$combined <- paste(edges$a, edges$b)
  # Initializing list of pvalues vector
  p <- rep(NA, dim(edges)[1])
  z <- rep(NA, length(p))
  cor_CTL <- rep(NA, length(p))
  cor_AD <- rep(NA, length(p))
  
  # Potential improvement: Use purrr instead of this while loop.
  i <- 1
  while(i <= dim(edges)[1]){
    
    if(i%%100000 == 0){
      cat("\n",i,"out of", dim(edges)[1],"times \n") # Printing out iteration number
    }
    
    cor1 <- cor( CTL1[ , edges$a[i] ], CTL2[ , edges$b[i] ], method = "spearman" )
    cor2 <- cor( AD1[ , edges$a[i] ], AD2[ , edges$b[i] ], method = "spearman" )
    
    cor_CTL[i] <- cor1
    cor_AD[i] <- cor2
    
    rtest_vals <- r.test(n = n.ctl, n2 = n.ad, r12 = cor1, r34 = cor2)
    
    p[i] <- rtest_vals$p
    # To keep info regarding LOC and GOC
    # LOC - Loss of correlation wrt Control, GOC - Gain of correlation wrt Control
    z[i] <- rtest_vals$z * sign(cor2 - cor1)
    
    i <- i + 1
  }
  
  
  # Print DC edges
  dcres <- cbind(edges, p, z, cor_CTL, cor_AD)
  rm(edges)
  
  mapping_file <- read.csv(paste("data/", 
                                 "Mapped_h37_ens_to_gene.csv", 
                                 sep = ""), 
                           stringsAsFactors = FALSE)
  
  cell_types <- c("ast", "mic", "neu", "oli")
  number_mg <- 1000
  number_mg_per_ct <- rep(number_mg, length(cell_types))
  
  ens1 <- dcres$a
  ens2 <- dcres$b
  map1 <- match(ens1, mapping_file$gene_id)
  map2 <- match(ens2, mapping_file$gene_id)
  g1 <- mapping_file$gene_name[map1]
  g2 <- mapping_file$gene_name[map2]
  
  dcres <- cbind(dcres, g1, g2)
  
  mg <- c()
  ct <- c()
  
  for(i in 1:length(cell_types)){
    ct_markers <- head(markers_df_human_brain$markers[which(markers_df_brain$cell == cell_types[i])],
                       n = number_mg_per_ct[i])
    ct <- c(ct, rep(cell_types[i], number_mg_per_ct[i]))
    mg <- c(mg, ct_markers)
  }
  
  
  marker_label_g1 <- rep(0, dim(dcres)[1])
  marker_label_g2 <- rep(0, dim(dcres)[1])
  
  match1 <- match(toupper(dcres$g1), mg)
  match2 <- match(toupper(dcres$g2), mg)
  
  marker_label_g1[which(is.na(match1))] = "NM"
  marker_label_g1[which(marker_label_g1 == 0)] = "M"
  marker_label_g2[which(is.na(match2))] = "NM"
  marker_label_g2[which(marker_label_g2 == 0)] = "M"
  
  dcres$m1 <- marker_label_g1
  dcres$m2 <- marker_label_g2
  
  dcres$c1 <- ct[match1]
  dcres$c2 <- ct[match2]
  after_ctc <- dcres
  
  rm(dcres)
  # Filtering only common DC-edges between the comparisons
  
  before_ctc$dc_edge <- paste(before_ctc$a, before_ctc$b, sep = "-")
  after_ctc$dc_edge <- paste(after_ctc$a, after_ctc$b, sep = "-")
  
  common_edges <- intersect(before_ctc$dc_edge, after_ctc$dc_edge)
  
  before_ctc_mm <- before_ctc %>% filter(m1 == "M", m2 == "M") %>% filter(dc_edge %in% common_edges)
  before_ctc_nmnm <- before_ctc %>% filter(m1 == "NM", m2 == "NM") %>% filter(dc_edge %in% common_edges)
  after_ctc_mm <- after_ctc %>% filter(m1 == "M", m2 == "M") %>% filter(dc_edge %in% common_edges)
  after_ctc_nmnm <- after_ctc %>% filter(m1 == "NM", m2 == "NM") %>% filter(dc_edge %in% common_edges)
  
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
         filename = paste0("looking_into_dc_", BR_after_ctc[1], BR_after_ctc[2], "_dc_before_ctc_edges", ".pdf"),
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

# Author: B P Kailash
# Purpose: Find differentially correlated edges between two brain regions

rm(list = ls())

.libPaths(c("libraries/", .libPaths()))

library(tidyverse)
library(optparse)

# INPUT to be specified by the user ---------------------------------------

# Brain regions (BR) under study BM10, BM22, BM44, BM36
# BR <- c(BR1, BR2) # Specify two brain regions at a time
# False discovery rate value (fdr_val)

option_list <- list(
  make_option(c('--br1'), help = 'Brain region 1', default = ""),
  make_option(c('--br2'), help = 'Brain region 2', default = ""),
  make_option(c('--fdr'), help = "fdr threshold", default = 0.01)
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

BR <- c(opt$br1, opt$br2)
message("**Brain region combination is ", print(BR), "**")
fdr_val <- opt$fdr
message("**FDR threshold is ", fdr_val, "**")

file_name <- paste("Edges", toString(fdr_val), sep = "")

# Loading necessary packages ----------------------------------------------

message("**Loading packages**")
source("scripts/libraries.R")
message("**Finished loading packages**")

# The following are the file names of the original gene expression files
# available on synapse https://www.synapse.org/#!Synapse:syn7391833 
# AMP-AD_MSBB_MSSM_BM_10.normalized.sex_race_age_RIN_PMI_exonicRate_rRnaRate_batch_adj.tsv
# AMP-AD_MSBB_MSSM_BM_22.normalized.sex_race_age_RIN_PMI_exonicRate_rRnaRate_batch_adj.tsv
# AMP-AD_MSBB_MSSM_BM_36.normalized.sex_race_age_RIN_PMI_exonicRate_rRnaRate_batch_adj.tsv
# AMP-AD_MSBB_MSSM_BM_44.normalized.sex_race_age_RIN_PMI_exonicRate_rRnaRate_batch_adj.tsv
# max(theAboveMatrix - GeneExpMatrixInFolder) is very small - hence they are same matrices

# OUTPUT files generated --------------------------------------------------

# 1) AD and Control matrices for the two brain regions under study with genes in ENS_ID
# Only common AD, CTL samples across the two brain regions are considered
# 2) Corrrelated edges - UnionCorrelatedEdgesBM36BM44.csv - genes in ENS_ID format
# All the gene pairs that are identified as significant are recorded in this file
# 3) DC files
#   a) p_values generated in r.test (note this is not adjusted)
#   b) DC edges in ENS_ID format, along with adjusted p values and z values

# DC Code --------------------------------------------------------------------

message("**Loading metadata and functions**")
source("scripts/utils.R")
message("**Finished loading metadata and functions**")

# Variables loaded - rna_metadata_processed, ensid_sampleid, clinical, genes_mapped

read_cell_type_corrected_gExp <- function(br, index, condition){
  if(BR[1] != BR[2]){
    df <- read.csv(paste("results/", condition, BR[index], "_", BR[1], BR[2], "_cellTypeCorrected.csv",
                         sep = ""), stringsAsFactors = FALSE, row.names = 1)
  } else if(BR[1] == BR[2]){
    df <- read.csv(paste("results/", BR[1], "_AD_CTL_cellTypeCorrected.csv", sep = ""),
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
AD1 <- read_cell_type_corrected_gExp(br = BR, index = 1, condition = "AD")
AD2 <- read_cell_type_corrected_gExp(br = BR, index = 2, condition = "AD")
CTL1 <- read_cell_type_corrected_gExp(br = BR, index = 1, condition = "CTL")
CTL2 <- read_cell_type_corrected_gExp(br = BR, index = 2, condition = "CTL")

colnames(AD1) <- colnames(AD2) <- colnames(CTL1) <- colnames(CTL2) <- ens_genes
# BR is the normalized gene expression values which is corrected for batch effects
# It is a matrix with dimension (#genes X #samples)

# Finding number of AD and control samples in BR1 and BR2
n.ad <- dim(AD1)[1]
n.ctl <- dim(CTL1)[1]

# How many values should cor_AD12 contain?
# It should contain {#genes}C{1} * {#genes}C{1} number of values, here in this case,
# It should be 23201*23201 = 538286401 correlation values
# Spearman is used here instead of Pearson, because Spearman works better for non-linear data.

message("**Beginning calculating spearman correlations within AD, and within CTL**")
if(!file.exists(paste("results/",BR[1], BR[2], "CTCorrectedADCorrelatedEdges.csv", sep = "")) &&
   !file.exists(paste("results/",BR[1], BR[2], "CTCorrectedCTLCorrelatedEdges.csv", sep = ""))){
  ########### Gene-Gene Correlation analysis
  # Rows are from BR1, columns are from BR2
  cor_AD12 <- cor(AD1, AD2, method = "spearman")
  message("**Finished calculating spearman correlations within AD**")
  cor_CTL12 <- cor(CTL1, CTL2, method = "spearman")
  message("**Finished calculating spearman correlations within CTL**")
  
  stopifnot(colnames(cor_AD12) == rownames(cor_AD12), colnames(cor_AD12) == ens_genes)
  stopifnot(colnames(cor_CTL12) == rownames(cor_CTL12), rownames(cor_CTL12) == ens_genes)
  
  stopifnot(nrow(AD2) == n.ad)
  p_AD12 <- pvalfromtcor(cor_AD12, n.ad)
  stopifnot(nrow(CTL2) == n.ctl)
  p_CTL12 <- pvalfromtcor(cor_CTL12, n.ctl)
  
  message("**Benjamini Hochberg adjusting p values**")
  # Benjamini-Hochberg pval adjustment
  A_AD12 <- matrix(p.adjust(p_AD12, method = "BH"), 
                   nrow = dim(p_AD12)[1], ncol = dim(p_AD12)[2])
  A_CTL12 <- matrix(p.adjust(p_CTL12, method = "BH"), 
                    nrow = dim(p_CTL12)[1], ncol = dim(p_CTL12)[2])
  
  # cbind(x,y) - Adds y after columns of x
  
  # The following property of the indices is used to obtain the correlated edges
  stopifnot(A_AD12[23202] == A_AD12[1,2])
  stopifnot(row(A_AD12)[23202] == 1)
  stopifnot(col(A_AD12)[23202] == 2)
  
  indices <- which(A_AD12 <= fdr_val)
  cor_AD12_edges <- cbind(ens_genes[row(A_AD12)[indices]], 
                          ens_genes[col(A_AD12)[indices]],
                          as.double(cor_AD12[indices]),
                          as.double(p_AD12[indices]),
                          as.double(A_AD12[indices]))
  colnames(cor_AD12_edges) <- c("g1-BR1", "g2-BR2", "r_ad", "p-val", "p-adj")
  
  indices <- which(A_CTL12 <= fdr_val)
  cor_CTL12_edges <- cbind(ens_genes[row(A_CTL12)[indices]],
                           ens_genes[col(A_CTL12)[indices]],
                           as.double(cor_CTL12[indices]),
                           as.double(p_CTL12[indices]),
                           as.double(A_CTL12[indices]))
  colnames(cor_CTL12_edges) <- c("g1-BR1", "g2-BR2", "r_ctl", "p-val", "p-adj")
  rm(indices)
  
  # hist(as.double(cor_AD12_edges[,2]))
  # hist(as.double(cor_CTL12_edges[,2]))
  cor_AD12_edges <- data.frame(cor_AD12_edges, stringsAsFactors = FALSE)
  cor_CTL12_edges <- data.frame(cor_CTL12_edges, stringsAsFactors = FALSE)
  
  cor_AD12_edges$combined <- paste(cor_AD12_edges$g1.BR1,cor_AD12_edges$g2.BR2)
  cor_CTL12_edges$combined <- paste(cor_CTL12_edges$g1.BR1,cor_CTL12_edges$g2.BR2)
  
  write.csv(cor_AD12_edges, paste("results/", BR[1], BR[2], "CTCorrectedADCorrelatedEdges.csv", sep = ""))
  write.csv(cor_CTL12_edges, paste("results/", BR[1], BR[2], "CTCorrectedCTLCorrelatedEdges.csv", sep = ""))
  message("**Number of significantly correlated edges in AD ", dim(cor_AD12_edges)[1], "**")
  message("**Number of significantly correlated edges in CTL ", dim(cor_CTL12_edges)[1], "**")
  
  # Find union of cor_ad_34_edges, cor_ctl_34_edges
  # NOTE: A particular edge in union12 need not have a signficant edge in both AD and CTl
  # but they definitely will have significance in either AD or CTL.
  union12 <- union(cor_AD12_edges$combined, cor_CTL12_edges$combined)
  stopifnot(union12 != '')
  message("**Number of significantly correlated edges in the union ", length(union12), "**")
} else {
  
  message("**Pre-computed spearman correlations within AD, and within CTL are being used**")
  
}

message("** Finding union edges **")
if(!file.exists(paste("results/","UnionCorrelatedEdges_ct_correction",BR[1],BR[2], ".csv", sep = ""))){
  # For a given edge a-b
  edges <- strsplit(union12, " ")
  mat  <- matrix(unlist(edges), ncol = 2, byrow = TRUE)
  df   <- as.data.frame(mat, stringsAsFactors = FALSE)
  colnames(df) <- c("a", "b")
  
  edges <- df
  
  edges$combined <- paste(edges$a, edges$b)
  union_ad <- cor_AD12_edges[match(edges$combined, cor_AD12_edges$combined),c("r_ad", "p.val", "p.adj")]
  colnames(union_ad) <- c("r_ad", "p.val.ad", "p.adj.ad")
  union_ctl <- cor_CTL12_edges[match(edges$combined, cor_CTL12_edges$combined), c("r_ctl", "p.val", "p.adj")]
  colnames(union_ctl) <- c("r_ctl", "p.val.ctl", "p.adj.ctl")
  
  ### Removing correlation files to reduce memory
  
  rm(cor_AD12_edges, cor_CTL12_edges)
  
  edges <- cbind(dplyr::select(edges, c(a,b)), union_ad, union_ctl)
  
  write.csv(edges, paste("results/","UnionCorrelatedEdges_ct_correction",BR[1],BR[2], ".csv", sep = ""))
} else {
  
  message("**Pre-computed UnionCorrelatedEdges are being used**")
  edges <- read.csv(paste("results/", "UnionCorrelatedEdges_ct_correction",BR[1],BR[2], ".csv", sep = ""))
}
message("** Finished finding/loading Union correlated edges **")

# DC Analysis -------------------------------------------------------------

message("** Beginning DC analysis **")

# r.test DC edges ---------------------------------------------------------

if(!file.exists(paste("results/", "p_vals_without_adjustment_ct_correction",BR[1],BR[2],toString(file_name),".rda",sep="")) &&
   !file.exists(paste("results/", "DCEdges_ct_correction",BR[1],BR[2],".csv",sep=""))){
  
  # Initializing list of pvalues vector
  p <- rep(NA, dim(edges)[1])
  z <- rep(NA, length(p))
  cor_CTL <- rep(NA, length(p))
  cor_AD <- rep(NA, length(p))
  
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
 
  # NOTE the following file contains p values before adjustment
  save(p, file = paste("results/", "p_vals_without_adjustment_ct_correction",BR[1],BR[2],toString(file_name),".rda",sep=""))
  
  p12 <- p
  p_adjusted <- p.adjust(p12, method = 'BH')
  
  # Print DC edges
  dcres <- cbind(edges, p_adjusted, z, cor_CTL, cor_AD)
  rm(edges)
  r <- which(p_adjusted <= fdr_val)
  dc_rtest <- dcres[r,]
  write.csv(dc_rtest, paste("results/", "DCEdges_ct_correction",BR[1],BR[2],".csv",sep=""))
} else {
  dc_rtest <- read.csv(paste("results/", "DCEdges_ct_correction",BR[1],BR[2],".csv",sep=""))
}

message("**Number of DC edges ", dim(dc_rtest)[1], "**")

rm(list = ls()) 

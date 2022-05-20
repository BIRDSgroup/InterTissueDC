
# Author: B P Kailash
# Purpose: Find differentially correlated edges between two brain regions

rm(list = ls())

# Loading necessary packages ----------------------------------------------

.libPaths(c("/data/kailash/DDP/libraries/", .libPaths()))
source("/data/kailash/DDP/src/main/scripts/libraries.R")

# INPUT to be specified by the user ---------------------------------------

server = "IBSE"
BR <- c("BM36", "BM44")
fdr_val <- 0.01

# Brain regions (BR) under study BM10, BM22, BM44, BM36
# BR <- c(BR1, BR2) # Specify two brain regions at a time
# False discovery rate value (fdr_val)
# Server can be either "IBSE" or "BIRDS"

# Setting file directories ------------------------------------------------

# Path to the original gene expression files (path_geneExp)
# Path to the data folder
# Path to the results folder

# The following are the file names of the original gene expression files
# available on synapse https://www.synapse.org/#!Synapse:syn7391833 
# AMP-AD_MSBB_MSSM_BM_10.normalized.sex_race_age_RIN_PMI_exonicRate_rRnaRate_batch_adj.tsv
# AMP-AD_MSBB_MSSM_BM_22.normalized.sex_race_age_RIN_PMI_exonicRate_rRnaRate_batch_adj.tsv
# AMP-AD_MSBB_MSSM_BM_36.normalized.sex_race_age_RIN_PMI_exonicRate_rRnaRate_batch_adj.tsv
# AMP-AD_MSBB_MSSM_BM_44.normalized.sex_race_age_RIN_PMI_exonicRate_rRnaRate_batch_adj.tsv

# http://10.24.6.51:8787/ - IBSE server

if(server == "IBSE"){
  # The following is for runnning on IBSE server Actinium
  
  setwd("/data/kailash/DDP/src/main/scripts/")
  
  path_geneExp <- "/data/srivatsan/Original_Expression_Files/Original_"
  path_data <- "/data/kailash/DDP/data/"
  path_results <- "/data/kailash/DDP/results/"
  
  file_name <- paste("Edges",toString(fdr_val),sep = "")
}

# http://10.21.226.107:8787/ - BIRDS server

if(server == "BIRDS"){
  # The following is for running on BIRDS server
  
  setwd("/data/users/be15b007/DDP/src/main/scripts/")
  
  path_geneExp <- "/data/users/be15b007/DDP/data/Original_"
  path_data <- "/data/users/be15b007/DDP/data/"
  path_results <- "/data/users/be15b007/DDP/results/"
  
  file_name <- paste("Edges",toString(fdr_val),sep="")
}

# OUTPUT files generated --------------------------------------------------

# 1) AD and Control matrices for the two brain regions under study, genes in ENS_ID
# Only common AD, CTL samples across the two brain regions are considered
# 2) Correlated edges - UnionCorrelatedEdgesBM36BM44.csv - genes in ENS_ID format
# All the gene pairs that are identified as significant are recorded in this file
# 3) DC files
#   a) p_values generated in r.test (note this is not adjusted)
#   b) DC edges in ENS_ID format, along with adjusted p values and z values

# DC Code --------------------------------------------------------------------

source("/data/kailash/DDP/src/main/scripts/utils.R")
# Variables loaded - rna_metadata_processed, ensid_sampleid, clinical, genes_mapped

# BR is the normalized gene expression values which is corrected for batch effects
# It is a matrix with dimension (#genes X #samples)
BR1 <- readBMfile(BR[1])
BR2 <- readBMfile(BR[2])
ens_genes <- colnames(BR1)
stopifnot(ens_genes == colnames(BR2))

BR1$sample_gp <- find_sample_group(BR[1])
BR2$sample_gp <- find_sample_group(BR[2])

prev_rownames_BR1 <- row.names(BR1)
prev_rownames_BR2 <- row.names(BR2)

row.names(BR1) <- find_individual_id(BR[1])
row.names(BR2) <- find_individual_id(BR[2])

ctl_ids_BR1 <- row.names(BR1)[which(BR1$sample_gp == "Control")]
ad_ids_BR1 <- row.names(BR1)[which(BR1$sample_gp == "Definite AD")]

ctl_ids_BR2 <- row.names(BR2)[which(BR2$sample_gp == "Control")]
ad_ids_BR2 <- row.names(BR2)[which(BR2$sample_gp == "Definite AD")]

BR1 <- subset(BR1,select = -sample_gp)
BR2 <- subset(BR2, select = -sample_gp)

ad_ids_BR1_2 <- intersect(ad_ids_BR1, ad_ids_BR2)
write.csv(ad_ids_BR1_2, paste(path_results, "AD_ids", BR[1], BR[2], ".csv", sep = ""), 
          row.names = FALSE)
ctl_ids_BR1_2 <- intersect(ctl_ids_BR1, ctl_ids_BR2)
write.csv(ctl_ids_BR1_2, paste(path_results, "CTL_ids", BR[1], BR[2], ".csv", sep = ""), 
          row.names = FALSE)

# Finding number of AD and control samples in BR1 and BR2
n.ad <- length(ad_ids_BR1_2)
n.ctl <- length(ctl_ids_BR1_2)

AD1 <- BR1[ad_ids_BR1_2,]
AD2 <- BR2[ad_ids_BR1_2,]
CTL1 <- BR1[ctl_ids_BR1_2,]
CTL2 <- BR2[ctl_ids_BR1_2,]

write.csv(AD1, paste(path_results, "AD", BR[1], "_", BR[1], BR[2], ".csv", sep = ""))
write.csv(CTL1, paste(path_results, "CTL", BR[1], "_", BR[1], BR[2], ".csv", sep = ""))
write.csv(AD2, paste(path_results, "AD", BR[2], "_", BR[1], BR[2], ".csv", sep = ""))
write.csv(CTL2, paste(path_results, "CTL", BR[2], "_", BR[1], BR[2], ".csv", sep = ""))

# How many values should cor_AD12 contain?
# It should contain {#genes}C{1} * {#genes}C{1} number of values, here in this case,
# It should be 23201*23201 = 538286401 correlation values
# Spearman is used here instead of Pearson, because Spearman works better for non-linear data.

########### Gene-Gene Correlation analysis
# Rows are from BR1, columns are from BR2
cor_AD12 <- cor(AD1, AD2, method = "spearman")
cor_CTL12 <- cor(CTL1, CTL2, method = "spearman")

stopifnot(colnames(cor_AD12) == rownames(cor_AD12), colnames(cor_AD12) == ens_genes)
stopifnot(colnames(cor_CTL12) == rownames(cor_CTL12), rownames(cor_CTL12) == ens_genes)

stopifnot(nrow(AD2) == n.ad)
p_AD12 <- pvalfromtcor(cor_AD12, n.ad)
stopifnot(nrow(CTL2) == n.ctl)
p_CTL12 <- pvalfromtcor(cor_CTL12, n.ctl)

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

write.csv(cor_AD12_edges, paste(path_results,BR[1], BR[2], "ADCorrelatedEdges.csv", sep = ""))
write.csv(cor_CTL12_edges, paste(path_results,BR[1], BR[2], "CTLCorrelatedEdges.csv", sep = ""))

# Find union of cor_ad_34_edges, cor_ctl_34_edges
# NOTE: A particular edge in union12 need not have a signficant edge in both AD and CTl
# but they definitely will have significance in either AD or CTL.
union12 <- union(cor_AD12_edges$combined, cor_CTL12_edges$combined)
stopifnot(union12 != '')

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

edges <- cbind(dplyr::select(edges, c(a,b)), union_ad, union_ctl)

write.csv(edges, paste(path_results,"UnionCorrelatedEdges",BR[1],BR[2], ".csv", sep = ""))

# DC Analysis -------------------------------------------------------------

# From Brain region 1
t1 <- AD1[,edges[,1]]
t2 <- CTL1[,edges[,1]]

# From Brain region 2
t3 <- AD2[,edges[,2]]
t4 <- CTL2[,edges[,2]]

# Zip fastener(x,y)
# col1 is col1 of x
# col2 is col1 of y
# and so on

res_AD12 <- zipFastener(t1,t3)
res_CTL12 <- zipFastener(t2,t4)

res_AD12 <- as.data.frame(res_AD12)
res_CTL12 <- as.data.frame(res_CTL12)

# r.test DC edges ---------------------------------------------------------

# Initializing list of pvalues vector
p <- rep(NA,length(colnames(res_AD12))/2 + 1)
z <- rep(NA, length(p))
cor_CTL <- rep(NA, length(p))
cor_AD <- rep(NA, length(p))

i <- 1
j <- 1

while(i <= length(colnames(res_AD12))){
  cat("\n",i,"out of", length(colnames(res_AD12)),"times \n") # Printing out iteration number
  temp1 <- res_CTL12[,c(i,i+1)]
  temp2 <- res_AD12[,c(i,i+1)]
  colnames(temp1) <- colnames(temp2) <- c('a','b')
  
  cor1 <- cor(temp1[,'a'], temp1[,'b'], method = "spearman")
  cor2 <- cor(temp2[,'a'], temp2[,'b'], method = "spearman")
  
  cor_CTL[j] <- cor1
  cor_AD[j] <- cor2
  
  rtest_vals <- r.test(n = n.ctl, n2 = n.ad, r12 = cor1, r34 = cor2)
  
  p[j] <- rtest_vals$p
  z[j] <- rtest_vals$z * sign(cor2 - cor1)
  
  # To keep info regarding LOC and GOC
  # LOC - Loss of correlation wrt Control, GOC - Gain of correlation wrt Control
  
  # i + 2 because we are taking the bm36 and bm44 adjacent columns gene exp values together
  j <- j + 1
  i <- i + 2
}

# Removing NAs in the end of the vector
indices <- which(!is.na(p))
stopifnot(length(is.na(p)) > 2)
p <- p[indices]
z <- z[indices]
cor_CTL <- cor_CTL[indices]
cor_AD <- cor_AD[indices]
# NOTE the following file contains p values before adjustment
save(p, file = paste(path_results,"p_vals_without_adjustment",BR[1],BR[2],toString(file_name),".rda",sep=""))

p12 <- p
p_adjusted <- p.adjust(p12, method = 'BH')

# Print DC edges
dcres <- cbind(edges, p_adjusted, z, cor_CTL, cor_AD)
rm(edges)
r <- which(p_adjusted <= fdr_val)
dc_rtest <- dcres[r,]
write.csv(dc_rtest, paste(path_results,"DCEdges",BR[1],BR[2],".csv",sep=""))

cat("Number of significantly correlated edges in AD", dim(cor_AD12_edges)[1], "\n")
cat("Number of significantly correlated edges in CTL", dim(cor_CTL12_edges)[1], "\n")
cat("Number of significantly correlated edges in the union", length(union12), "\n")
cat("Number of DC edges ", dim(dc_rtest)[1], "\n")

rm(list = ls()) 
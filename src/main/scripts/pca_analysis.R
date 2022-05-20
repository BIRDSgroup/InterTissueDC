rm(list = ls())

# Loading necessary packages ----------------------------------------------

.libPaths(c("/data/kailash/DDP/libraries/", .libPaths()))
source("/data/kailash/DDP/src/main/scripts/libraries.R")

# Global variables --------------------------------------------------------

path_geneExp <- "/data/srivatsan/Original_Expression_Files/Original_"
path_data <- "/data/kailash/DDP/data/"
path_results <- "/data/kailash/DDP/results/"
cell_type_correction <- TRUE

# Loading data ------------------------------------------------------------

setwd(path_results)
source("/data/kailash/DDP/src/main/scripts/utils.R")

if(cell_type_correction == FALSE){
  ct_name <- ""
  bm10 <- readBMfile("BM10")
  bm22 <- readBMfile("BM22")
  bm36 <- readBMfile("BM36")
  bm44 <- readBMfile("BM44")
  bm1 <- rbind(bm10, bm22, bm36, bm44)
  bm <- bm1
  
  barcodes_bm10 <- sapply(strsplit(rownames(bm10),"_"), function(x){x[3]})
  barcodes_bm22 <- sapply(strsplit(rownames(bm22),"_"), function(x){x[3]})
  barcodes_bm36 <- sapply(strsplit(rownames(bm36),"_"), function(x){x[3]})
  barcodes_bm44 <- sapply(strsplit(rownames(bm44),"_"), function(x){x[3]})
  barcodes <- as.numeric(c(barcodes_bm10, barcodes_bm22, barcodes_bm36, barcodes_bm44))
  stopifnot(length(barcodes) == dim(bm10)[1] + dim(bm22)[1] + dim(bm36)[1] + dim(bm44)[1])
  stopifnot(!duplicated(barcodes))
  
  barcodes[which(!barcodes %in% rna_metadata_processed$barcode)]
  barcodes_bm10[which(!barcodes_bm10 %in% rna_metadata_processed$barcode)]
  barcodes_bm22[which(!barcodes_bm22 %in% rna_metadata_processed$barcode)]
  barcodes_bm36[which(!barcodes_bm36 %in% rna_metadata_processed$barcode)]
  barcodes_bm44[which(!barcodes_bm44 %in% rna_metadata_processed$barcode)]
  
} else{
  ct_name <- "cell_type_corrected"
  bm10 <- read.csv(paste(path_results, "BM10_AD_CTL_cellTypeCorrected.csv", sep = ""), stringsAsFactors = FALSE, row.names = 1)
  bm22 <- read.csv(paste(path_results, "BM22_AD_CTL_cellTypeCorrected.csv", sep = ""), stringsAsFactors = FALSE, row.names = 1)
  bm36 <- read.csv(paste(path_results, "BM36_AD_CTL_cellTypeCorrected.csv", sep = ""), stringsAsFactors = FALSE, row.names = 1)
  bm44 <- read.csv(paste(path_results, "BM44_AD_CTL_cellTypeCorrected.csv", sep = ""), stringsAsFactors = FALSE, row.names = 1)
  bm1 <- rbind(bm10, bm22, bm36, bm44)
  colnames(bm1) <- ensid_sampleid[,1]
  bm <- data.frame(bm1)
  
  barcodes_bm10 <- sapply(strsplit(rownames(bm10),"_"), function(x){x[3]})
  barcodes_bm22 <- sapply(strsplit(rownames(bm22),"_"), function(x){x[3]})
  barcodes_bm36 <- sapply(strsplit(rownames(bm36),"_"), function(x){x[3]})
  barcodes_bm44 <- sapply(strsplit(rownames(bm44),"_"), function(x){x[3]})
  barcodes <- as.numeric(c(barcodes_bm10, barcodes_bm22, barcodes_bm36, barcodes_bm44))
  stopifnot(length(barcodes) == dim(bm10)[1] + dim(bm22)[1] + dim(bm36)[1] + dim(bm44)[1])
  stopifnot(!duplicated(barcodes))
  
  barcodes[which(!barcodes %in% rna_metadata_processed$barcode)]
  barcodes_bm10[which(!barcodes_bm10 %in% rna_metadata_processed$barcode)]
  barcodes_bm22[which(!barcodes_bm22 %in% rna_metadata_processed$barcode)]
  barcodes_bm36[which(!barcodes_bm36 %in% rna_metadata_processed$barcode)]
}

bm_prcomp <- prcomp(bm)

# rna_metadata_processed is loaded in utils.R
metadata_mapped <- rna_metadata_processed[match(barcodes, rna_metadata_processed$barcode),]
metadata_pca <- metadata_mapped[,c("RACE.inferred", "SEX.inferred", "RIN", "batch", "DxCondition", "SEX", "PlaqueMean", "Action", "Apo1", "Apo2", "bbscore", "CDR", "AOD", "PMI", "BrodmannArea")]
dat <- data.frame(cbind(metadata_pca, bm))

# Expectation: No clustering

autoplot(bm_prcomp, data = dat)
autoplot(bm_prcomp, data = dat, colour = "SEX.inferred")
autoplot(bm_prcomp, data = dat, colour = "RACE.inferred")
autoplot(bm_prcomp, data = dat, colour = "RIN")
autoplot(bm_prcomp, data = dat, colour = "batch")
autoplot(bm_prcomp, data = dat, colour = "SEX")
autoplot(bm_prcomp, data = dat, colour = "AOD")
autoplot(bm_prcomp, data = dat, colour = "PMI")

# Expected clustering
autoplot(bm_prcomp, data = dat, colour = "BrodmannArea", frame.type = "norm")
autoplot(bm_prcomp, data = dat, colour = "PlaqueMean")
autoplot(bm_prcomp, data = dat, colour = "Apo1")
autoplot(bm_prcomp, data = dat, colour = "CDR")
autoplot(bm_prcomp, data = dat, colour = "DxCondition")

# PCA on top 50% of genes showing variance --------------------------------

dat2 <- bm
dat2 <- dat2[,order(apply(bm,2,var))]
dat2 <- dat2[,floor(length(colnames)/2):length(colnames(bm))]
dat2_prcomp <- prcomp(dat2)
barcodes2 <- sapply(strsplit(rownames(dat2),"_"), function(x){x[3]})
RNA_metadata_mapped2 <- data.frame(rna_metadata)[match(barcodes2,rna_metadata$barcode),]
metadata_pca2 <- metadata_mapped2[,c("RACE.inferred", "SEX.inferred", "RIN", "batch", "DxCondition", "SEX", "PlaqueMean", "Action", "Apo1", "Apo2", "bbscore", "CDR", "AOD", "PMI", "BrodmannArea")]
dat2 <- data.frame(cbind(metadata_pca2, dat2))
autoplot(dat2_prcomp, data=dat2, colour = "BrodmannArea", frame.type = "norm", x = 1, y = 2)
autoplot(dat2_prcomp, data=dat2, colour = "BrodmannArea", frame.type = "norm", x = 1, y = 3)
autoplot(dat2_prcomp, data=dat2, colour = "BrodmannArea", frame.type = "norm", x = 2, y = 3)


# PCA on scaled data ------------------------------------------------------

dat3 <- scale(bm)
dat3_prcomp <- prcomp(dat3)
p1 <- autoplot(dat3_prcomp, data=dat, colour = "BrodmannArea", frame.type = "norm")
p2 <- autoplot(dat3_prcomp, data=dat, colour = "DxCondition", frame.type = "norm")
p3 <- autoplot(dat3_prcomp, data=dat, colour = "RIN", frame.type = "norm")
p4 <- autoplot(dat3_prcomp, data=dat, colour = "SEX", frame.type = "norm")

ggsave(path = path_results,
       device = "png",
       dpi = 300,
       width = 18,
       height = 10,
       filename = paste("PCA_scaled_combined_BM_DxCondition_SEX_RIN", ct_name, ".png"),
       plot = ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2))

# PCA scaled, top 50 ------------------------------------------------------

dat4 <- bm[,order(apply(bm,2,var))]
dat4 <- dat4[,floor(length(colnames(bm))/2):length(colnames(bm))]
dat4_prcomp <- prcomp(dat4)
autoplot(dat4_prcomp, data=dat2, colour = "BrodmannArea", frame.type = "norm", x = 1, y = 2)
autoplot(dat4_prcomp, data=dat2, colour = "BrodmannArea", frame.type = "norm", x = 2, y = 3)
+ dark_theme_gray(base_family = "Fira Sans Condensed Light", base_size = 14)

# Boxplots to check how data was normalized -------------------------------

library(ggplot2)

dat6 <- data.frame(cbind(bm10[,1], bm10[,2], bm10[,3], bm10[,4], bm10[,5]))
colnames(dat6) <- colnames(bm10)[1:5]
rownames(dat6) <- rownames(bm10)
boxplot(dat6, names = colnames(bm10)[1:5], xlab = "Gene", ylab = "GeneExp", main = "BM10")
ggplot(data.frame(bm10[,1]), aes(y = colnames(bm10)[1])) + geom_boxplot()

bm10_n <- bm10[,order(apply(bm10,2,var))]
bm22_n <- bm22[,order(apply(bm22,2,var))]
bm36_n <- bm36[,order(apply(bm36,2,var))]
bm44_n <- bm44[,order(apply(bm44,2,var))]

# Boxplot of unordered top 5 genes
box1 <- boxplot(bm10_n[,1], bm10_n[,2], bm10_n[,3], bm10_n[,4], bm10_n[,5], names = colnames(bm10_n)[1:5], xlab = "Gene", ylab = "GeneExp", main = "BM10")
box2 <- boxplot(bm22_n[,1], bm22_n[,2], bm22_n[,3], bm22_n[,4], bm22_n[,5], names = colnames(bm10_n)[1:5], xlab = "Gene", ylab = "GeneExp", main = "BM22")
box3 <- boxplot(bm36_n[,1], bm36_n[,2], bm36_n[,3], bm36_n[,4], bm36_n[,5], names = colnames(bm10_n)[1:5], xlab = "Gene", ylab = "GeneExp", main = "BM36")
box4 <- boxplot(bm44_n[,1], bm44_n[,2], bm44_n[,3], bm44_n[,4], bm44_n[,5], names = colnames(bm10_n)[1:5], xlab = "Gene", ylab = "GeneExp", main = "BM44")

# Boxplot of top 5 genes with highest variance
box5 <- boxplot(bm10[,6], bm10[,7], bm10[,8], bm10[,9], bm10[,10], names = colnames(bm10)[1:5], xlab = "Gene", ylab = "GeneExp", main = "BM10")
box6 <- boxplot(bm22[,6], bm22[,7], bm22[,8], bm22[,9], bm22[,10], names = colnames(bm10)[1:5], xlab = "Gene", ylab = "GeneExp", main = "BM22")
box7 <- boxplot(bm36[,6], bm36[,7], bm36[,8], bm36[,9], bm36[,10], names = colnames(bm10)[1:5], xlab = "Gene", ylab = "GeneExp", main = "BM36")
box8 <- boxplot(bm44[,6], bm44[,7], bm44[,8], bm44[,9], bm44[,10], names = colnames(bm10)[1:5], xlab = "Gene", ylab = "GeneExp", main = "BM44")
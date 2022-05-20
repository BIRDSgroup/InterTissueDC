# Cell-type gene expression correction using a linear model

rm(list = ls())

# Loading necessary packages ----------------------------------------------

.libPaths(c("/data/kailash/DDP/libraries/", .libPaths()))
source("/data/kailash/DDP/src/main/scripts/libraries.R")

# Global variables --------------------------------------------------------

path_geneExp <- "/data/srivatsan/Original_Expression_Files/Original_"
path_data <- "/data/kailash/DDP/data/"
path_results <- "/data/kailash/DDP/results/"

source("/data/kailash/DDP/src/main/scripts/utils.R")

# Functions ---------------------------------------------------------------

make_w_matrix <- function(markerGenes, cellTypes, numberMarkersPerCellType){
  mg <- c()
  for (i in 1:length(cellTypes)){
    mg_ct <- head(markerGenes$markers[which(markerGenes$cell == cellTypes[i])], n = numberMarkersPerCellType[i])
    mg <- c(mg, mg_ct)
  }
  
  w <- matrix(0L, nrow = sum(numberMarkersPerCellType), ncol = length(cellTypes))
  colnames(w) <- cellTypes
  rownames(w) <- mg
  
  sum_number_markers <- rep(0, length(numberMarkersPerCellType))
  for (i in 1:length(sum_number_markers)){
    if (i == 1){
      sum_number_markers[i] = numberMarkersPerCellType[i]
    } else {
      sum_number_markers[i] = numberMarkersPerCellType[i] + sum_number_markers[i-1]
    }
  }
  
  for (i in 1:length(cellTypes)){
    if (i == 1){
      w[,cellTypes[i]][1:numberMarkersPerCellType[i]] = rep(1, numberMarkersPerCellType[i])
    } 
    else{
      w[,cellTypes[i]][(sum_number_markers[i-1]+1) : sum_number_markers[i]] = rep(1, numberMarkersPerCellType[i])
    }
  }
  
  return(list(w, mg))
}

read_BR_AD_CTL <- function(BrainRegion){
  # Getting the expression matrix
  BR <- readBMfile(BrainRegion)
  colnames(BR) <- genes_mapped
  BR$sample_gp <- find_sample_group(BrainRegion)
  BR_AD_Ctl <- BR[c(which(BR$sample_gp == "Control"), which(BR$sample_gp == "Definite AD")),]
  return(BR_AD_Ctl)
}

get_SPV_BrainRegion <- function(BrainRegion, MarkerGenes){
  
  BR_AD_Ctl <- read_BR_AD_CTL(BrainRegion)
  sample_gp <- BR_AD_Ctl$sample_gp
  BR_AD_Ctl <- subset(BR_AD_Ctl, select = -c(sample_gp))
  
  mg_new <- intersect(mg, colnames(BR_AD_Ctl))
  w_new <- w[mg_new,]
  
  # Finding E_m
  # NOTE: length(intersect(mg, colnames(bm36))) = 4181 != 5500
  # NOTE: length(which(is.na(colnames(bm36) == ""))) = 1029
  
  E_m <- t(BR_AD_Ctl[,mg_new])
  colnames(E_m) <- rownames(BR_AD_Ctl)
  spvs <- as.data.frame(my_getAllSPVs(brainRegion = BrainRegion, data = E_m,
                                      grp = sample_gp, dataTag = w_new,
                                      method = "mixed", plot = T, mix.par = 0.3))
  return(spvs)
}

readCellTypeFreq <- function(BR, whichBR, condition){
  file <- read.csv(paste(path_results, BR[whichBR], "_CellCODE_80mg_4cellTypes_", BR[1], BR[2], "_", 
                         condition, ".csv", sep = ""), stringsAsFactors = FALSE)
  row.names(file) <- file[,1]
  file <- file[, -1]
  return(file)
}

# X is the cell type frequencies - training data - #samples X # cellTypes
# y is the gene expression data - training data - #samples X 1 gene
# fit(X,y) - linear model
# We are predicting the gene's expression based on the model
# Then we correct for the gene's expression based on the prediction

predict_gExp <- function(gExp, cf){
  
  # Rows are samples, Columns are genes and cell-types
  colnames(gExp) <- gene_list <- ensid_sampleid[,1]
  stopifnot(row.names(gExp) == row.names(cf))
  df <- cbind(gExp, cf)
  gExp[!is.na(gExp)] <- NA # Making the original gene expression matrix into all NA
  df2 <- gExp # Cell-type corrected gene expression data will be added in place
  
  count <- 0
  for(i in gene_list){
    model <- lm(paste0(i, "~", "ast + mic + neu + oli"), data = df)
    pred <- predict(model, df)
    ct_correction <- pred - rep(coef(model)["(Intercept)"], length(pred))
    df2[,i] <- ct_correction
    count = count + 1
    cat(paste(count, "of", length(gene_list), "\n"))
  }
  return(t(df2))
}

# Getting the spvs ------------------------------------------------------

# CellCODE 80 mg, 4 cell types seems to be the best performing model
ip <- make_w_matrix(markers_df_human_brain, c("ast", "mic", "neu", "oli"), c(20,20,20,20))
w <- ip[[1]]
mg <- ip[[2]]

# Produces a heatmap plot of all pairwise correlations. The colors are green to white to red 
# with white being 0 correlation.

BM10_spv <- get_SPV_BrainRegion("BM10", MarkerGenes = mg)
write.csv(BM10_spv,
          paste(path_results, "BM10_CellCODE_80mg_4cellTypes.csv", sep = ""))

BM22_spv <- get_SPV_BrainRegion("BM22", MarkerGenes = mg)
write.csv(BM22_spv,
          paste(path_results, "BM22_CellCODE_80mg_4cellTypes.csv", sep = ""))

BM36_spv <- get_SPV_BrainRegion("BM36", MarkerGenes = mg)
write.csv(BM36_spv,
          paste(path_results, "BM36_CellCODE_80mg_4cellTypes.csv", sep = ""))

BM44_spv <- get_SPV_BrainRegion("BM44", MarkerGenes = mg)
write.csv(BM44_spv,
          paste(path_results, "BM44_CellCODE_80mg_4cellTypes.csv", sep = ""))

# Loading AD, CTL gene files to be cell-type corrected -------------------

bm10_geneExp <- read_BR_AD_CTL("BM10")
bm10_sample_gp <- bm10_geneExp$sample_gp
bm10_geneExp <- subset(bm10_geneExp, select = -c(sample_gp))

bm22_geneExp <- read_BR_AD_CTL("BM22")
bm22_sample_gp <- bm22_geneExp$sample_gp
bm22_geneExp <- subset(bm22_geneExp, select = -c(sample_gp))

bm36_geneExp <- read_BR_AD_CTL("BM36")
bm36_sample_gp <- bm36_geneExp$sample_gp
bm36_geneExp <- subset(bm36_geneExp, select = -c(sample_gp))

bm44_geneExp <- read_BR_AD_CTL("BM44")
bm44_sample_gp <- bm44_geneExp$sample_gp
bm44_geneExp <- subset(bm44_geneExp, select = -c(sample_gp))

# Cell-type correction on whole matrix ------------------------------------

bm10_pred <- predict_gExp(bm10_geneExp, BM10_spv)
bm22_pred <- predict_gExp(bm22_geneExp, BM22_spv)
bm36_pred <- predict_gExp(bm36_geneExp, BM36_spv)
bm44_pred <- predict_gExp(bm44_geneExp, BM44_spv) # 8 26
write.csv(matrix(rexp(10), 2), "random.csv") # 6:05 startime 6 46 end time for BM36

# Linear regression model - cell type correction --------------------------

bm10_corrected <- bm10_geneExp - t(bm10_pred)
bm22_corrected <- bm22_geneExp - t(bm22_pred)
bm36_corrected <- bm36_geneExp - t(bm36_pred)
bm44_corrected <- bm44_geneExp - t(bm44_pred)

# Writing the files -------------------------------------------------------

write.csv(bm10_corrected, paste(path_results, "BM10_AD_CTL_cellTypeCorrected.csv", sep = ""))
write.csv(bm22_corrected, paste(path_results, "BM22_AD_CTL_cellTypeCorrected.csv", sep = ""))
write.csv(bm36_corrected, paste(path_results, "BM36_AD_CTL_cellTypeCorrected.csv", sep = ""))
write.csv(bm44_corrected, paste(path_results, "BM44_AD_CTL_cellTypeCorrected.csv", sep = ""))

# Finding common samples --------------------------------------------------

barcodes_bm10 <- as.vector(na.omit(sapply(strsplit(row.names(bm10_corrected),"_"), function(x){x[3]})))
barcodes_bm22 <- as.vector(na.omit(sapply(strsplit(row.names(bm22_corrected),"_"), function(x){x[3]})))
barcodes_bm36 <- as.vector(na.omit(sapply(strsplit(row.names(bm36_corrected),"_"), function(x){x[3]})))
barcodes_bm44 <- as.vector(na.omit(sapply(strsplit(row.names(bm44_corrected),"_"), function(x){x[3]})))

mdata <- read.csv(paste(path_data, "rna_metadata_onlyBAM_noExclude_noDup_Nov2018.csv", sep = ""), stringsAsFactors = FALSE)

return_individual_identifier <- function(barcodes){
  individual_id <- mdata$individualIdentifier[match(barcodes, mdata$barcode)]
  individual_id[which(is.na(individual_id))] <- 0
  return(individual_id)
}

# Convert to individual ids and then find intersection
row.names(bm10_corrected) <- return_individual_identifier(barcodes_bm10)
row.names(bm22_corrected) <- return_individual_identifier(barcodes_bm22)
row.names(bm36_corrected) <- return_individual_identifier(barcodes_bm36)
row.names(bm44_corrected) <- return_individual_identifier(barcodes_bm44)

bm10_corrected_ad <- bm10_corrected[which(mdata$DxCondition[match(row.names(bm10_corrected), mdata$individualIdentifier)] == "Definite AD"),]
bm10_corrected_ctl <- bm10_corrected[which(mdata$DxCondition[match(row.names(bm10_corrected), mdata$individualIdentifier)] == "Control"),]
stopifnot(dim(bm10_corrected_ad)[1] + dim(bm10_corrected_ctl)[1] == dim(bm10_corrected)[1])

bm22_corrected_ad <- bm22_corrected[which(mdata$DxCondition[match(row.names(bm22_corrected), mdata$individualIdentifier)] == "Definite AD"),]
bm22_corrected_ctl <- bm22_corrected[which(mdata$DxCondition[match(row.names(bm22_corrected), mdata$individualIdentifier)] == "Control"),]
stopifnot(dim(bm22_corrected_ad)[1] + dim(bm22_corrected_ctl)[1] == dim(bm22_corrected)[1])

bm36_corrected_ad <- bm36_corrected[which(mdata$DxCondition[match(row.names(bm36_corrected), mdata$individualIdentifier)] == "Definite AD"),]
bm36_corrected_ctl <- bm36_corrected[which(mdata$DxCondition[match(row.names(bm36_corrected), mdata$individualIdentifier)] == "Control"),]
stopifnot(dim(bm36_corrected_ad)[1] + dim(bm36_corrected_ctl)[1] == dim(bm36_corrected)[1])

bm44_corrected_ad <- bm44_corrected[which(mdata$DxCondition[match(row.names(bm44_corrected), mdata$individualIdentifier)] == "Definite AD"),]
bm44_corrected_ctl <- bm44_corrected[which(mdata$DxCondition[match(row.names(bm44_corrected), mdata$individualIdentifier)] == "Control"),]
stopifnot(dim(bm44_corrected_ad)[1] + dim(bm44_corrected_ctl)[1] == dim(bm44_corrected)[1])

bm10_1022_ad_geneExp_corrected <- bm10_corrected_ad[intersect(row.names(bm10_corrected_ad), row.names(bm22_corrected_ad)),]
bm22_1022_ad_geneExp_corrected <- bm22_corrected_ad[intersect(row.names(bm10_corrected_ad), row.names(bm22_corrected_ad)),]
bm10_1022_ctl_geneExp_corrected <- bm10_corrected_ctl[intersect(row.names(bm10_corrected_ctl), row.names(bm22_corrected_ctl)),]
bm22_1022_ctl_geneExp_corrected <- bm22_corrected_ctl[intersect(row.names(bm10_corrected_ctl), row.names(bm22_corrected_ctl)),]

bm10_1036_ad_geneExp_corrected <- bm10_corrected_ad[intersect(row.names(bm10_corrected_ad), row.names(bm36_corrected_ad)),]
bm36_1036_ad_geneExp_corrected <- bm36_corrected_ad[intersect(row.names(bm10_corrected_ad), row.names(bm36_corrected_ad)),]
bm10_1036_ctl_geneExp_corrected <- bm10_corrected_ctl[intersect(row.names(bm10_corrected_ctl), row.names(bm36_corrected_ctl)),]
bm36_1036_ctl_geneExp_corrected <- bm36_corrected_ctl[intersect(row.names(bm10_corrected_ctl), row.names(bm36_corrected_ctl)),]

bm10_1044_ad_geneExp_corrected <- bm10_corrected_ad[intersect(row.names(bm10_corrected_ad), row.names(bm44_corrected_ad)),]
bm44_1044_ad_geneExp_corrected <- bm44_corrected_ad[intersect(row.names(bm10_corrected_ad), row.names(bm44_corrected_ad)),]
bm10_1044_ctl_geneExp_corrected <- bm10_corrected_ctl[intersect(row.names(bm10_corrected_ctl), row.names(bm44_corrected_ctl)),]
bm44_1044_ctl_geneExp_corrected <- bm44_corrected_ctl[intersect(row.names(bm10_corrected_ctl), row.names(bm44_corrected_ctl)),]

bm22_2236_ad_geneExp_corrected <- bm22_corrected_ad[intersect(row.names(bm22_corrected_ad), row.names(bm36_corrected_ad)),]
bm36_2236_ad_geneExp_corrected <- bm36_corrected_ad[intersect(row.names(bm22_corrected_ad), row.names(bm36_corrected_ad)),]
bm22_2236_ctl_geneExp_corrected <- bm22_corrected_ctl[intersect(row.names(bm22_corrected_ctl), row.names(bm36_corrected_ctl)),]
bm36_2236_ctl_geneExp_corrected <- bm36_corrected_ctl[intersect(row.names(bm22_corrected_ctl), row.names(bm36_corrected_ctl)),]

bm22_2244_ad_geneExp_corrected <- bm22_corrected_ad[intersect(row.names(bm22_corrected_ad), row.names(bm44_corrected_ad)),]
bm44_2244_ad_geneExp_corrected <- bm44_corrected_ad[intersect(row.names(bm22_corrected_ad), row.names(bm44_corrected_ad)),]
bm22_2244_ctl_geneExp_corrected <- bm22_corrected_ctl[intersect(row.names(bm22_corrected_ctl), row.names(bm44_corrected_ctl)),]
bm44_2244_ctl_geneExp_corrected <- bm44_corrected_ctl[intersect(row.names(bm22_corrected_ctl), row.names(bm44_corrected_ctl)),]

bm36_3644_ad_geneExp_corrected <- bm36_corrected_ad[intersect(row.names(bm36_corrected_ad), row.names(bm44_corrected_ad)),]
bm44_3644_ad_geneExp_corrected <- bm44_corrected_ad[intersect(row.names(bm36_corrected_ad), row.names(bm44_corrected_ad)),]
bm36_3644_ctl_geneExp_corrected <- bm36_corrected_ctl[intersect(row.names(bm36_corrected_ctl), row.names(bm44_corrected_ctl)),]
bm44_3644_ctl_geneExp_corrected <- bm44_corrected_ctl[intersect(row.names(bm36_corrected_ctl), row.names(bm44_corrected_ctl)),]

# Writing the files -------------------------------------------------------

write.csv(bm10_1022_ad_geneExp_corrected, paste(path_results, "ADBM10_BM10BM22_cellTypeCorrected.csv", sep = ""))
write.csv(bm22_1022_ad_geneExp_corrected, paste(path_results, "ADBM22_BM10BM22_cellTypeCorrected.csv", sep = ""))
write.csv(bm10_1022_ctl_geneExp_corrected, paste(path_results, "CTLBM10_BM10BM22_cellTypeCorrected.csv", sep = ""))
write.csv(bm22_1022_ctl_geneExp_corrected, paste(path_results, "CTLBM22_BM10BM22_cellTypeCorrected.csv", sep = ""))

write.csv(bm10_1036_ad_geneExp_corrected, paste(path_results, "ADBM10_BM10BM36_cellTypeCorrected.csv", sep = ""))
write.csv(bm36_1036_ad_geneExp_corrected, paste(path_results, "ADBM36_BM10BM36_cellTypeCorrected.csv", sep = ""))
write.csv(bm10_1036_ctl_geneExp_corrected, paste(path_results, "CTLBM10_BM10BM36_cellTypeCorrected.csv", sep = ""))
write.csv(bm36_1036_ctl_geneExp_corrected, paste(path_results, "CTLBM36_BM10BM36_cellTypeCorrected.csv", sep = ""))

write.csv(bm10_1044_ad_geneExp_corrected, paste(path_results, "ADBM10_BM10BM44_cellTypeCorrected.csv", sep = ""))
write.csv(bm44_1044_ad_geneExp_corrected, paste(path_results, "ADBM44_BM10BM44_cellTypeCorrected.csv", sep = ""))
write.csv(bm10_1044_ctl_geneExp_corrected, paste(path_results, "CTLBM10_BM10BM44_cellTypeCorrected.csv", sep = ""))
write.csv(bm44_1044_ctl_geneExp_corrected, paste(path_results, "CTLBM44_BM10BM44_cellTypeCorrected.csv", sep = ""))

write.csv(bm22_2236_ad_geneExp_corrected, paste(path_results, "ADBM22_BM22BM36_cellTypeCorrected.csv", sep = ""))
write.csv(bm36_2236_ad_geneExp_corrected, paste(path_results, "ADBM36_BM22BM36_cellTypeCorrected.csv", sep = ""))
write.csv(bm22_2236_ctl_geneExp_corrected, paste(path_results, "CTLBM22_BM22BM36_cellTypeCorrected.csv", sep = ""))
write.csv(bm36_2236_ctl_geneExp_corrected, paste(path_results, "CTLBM36_BM22BM36_cellTypeCorrected.csv", sep = ""))

write.csv(bm22_2244_ad_geneExp_corrected, paste(path_results, "ADBM22_BM22BM44_cellTypeCorrected.csv", sep = ""))
write.csv(bm44_2244_ad_geneExp_corrected, paste(path_results, "ADBM44_BM22BM44_cellTypeCorrected.csv", sep = ""))
write.csv(bm22_2244_ctl_geneExp_corrected, paste(path_results, "CTLBM22_BM22BM44_cellTypeCorrected.csv", sep = ""))
write.csv(bm44_2244_ctl_geneExp_corrected, paste(path_results, "CTLBM44_BM22BM44_cellTypeCorrected.csv", sep = ""))

write.csv(bm36_3644_ad_geneExp_corrected, paste(path_results, "ADBM36_BM36BM44_cellTypeCorrected.csv", sep = ""))
write.csv(bm44_3644_ad_geneExp_corrected, paste(path_results, "ADBM44_BM36BM44_cellTypeCorrected.csv", sep = ""))
write.csv(bm36_3644_ctl_geneExp_corrected, paste(path_results, "CTLBM36_BM36BM44_cellTypeCorrected.csv", sep = ""))
write.csv(bm44_3644_ctl_geneExp_corrected, paste(path_results, "CTLBM44_BM36BM44_cellTypeCorrected.csv", sep = ""))

rm(list = ls())

# Using BM44 gene-coefficients, making scatter plots
# with other CellCODE models

# Global variables --------------------------------------------------------

path_results <- "/data/kailash/DDP/results/"
path_data <- "/data/kailash/DDP/data/"
path_geneExp <- "/data/srivatsan/Original_Expression_Files/Original_"
mapping_file <- read.csv(paste(path_data, 
                               "Mapped_h37_ens_to_gene.csv", 
                               sep = ""), 
                         stringsAsFactors = FALSE)
source("/data/kailash/DDP/src/main/scripts/utils.R")

# Loading libraries -------------------------------------------------------

.libPaths(c("/data/kailash/DDP/libraries/", .libPaths()))
source("/data/kailash/DDP/src/main/scripts/libraries.R")

# Functions ---------------------------------------------------------------

read_scale_HGNC_BM <- function(BR){
  file <- read.csv(paste(path_geneExp,BR,'.csv',sep=""),stringsAsFactors = FALSE)
  row.names(file) <- file[,1]
  file <- as.data.frame(t(file[,-1]))
  # svds function scales the gene expression data before deconvolving
  mm = apply(file,2,mean)
  tmp = sweep(file,2,mm, "-")
  file <- tmp
  # Changing genes from ENS_ID format to HGNC format
  genes <- toupper(mapping_file$gene_name[match(colnames(file), mapping_file$gene_id)])
  colnames(file) <- genes
  return(file)
}

read_cellCODE_est <- function(BR){
  ctfs <- read.csv(paste(path_results, BR, "_CellCODE_80mg_4cellTypes.csv", sep = ""),
                   stringsAsFactors = FALSE)
  row.names(ctfs) <- ctfs[,1]
  ctfs <- as.data.frame(ctfs[,-1])
  return(ctfs)
}

predict_ctf <- function(cellType, geneExp, geneCoef){
  mg <- geneCoef$genes[which(geneCoef$X.1 == cellType)]
  Em <- geneExp[, which(genes %in% mg)]
  umg <- geneCoef$su[match(mg, u$genes)]
  ugenes <- geneCoef$genes[match(mg, u$genes)]
  Em <- Em[, match(mg, colnames(Em))]
  stopifnot(mg == colnames(Em))
  stopifnot(ugenes == colnames(Em))
  Em <- data.matrix(Em)
  umg <- as.vector(umg)
  ctfWtEst <- umg %*% t(Em)
  ctfWtEst <- -ctfWtEst
  return(ctfWtEst)
}

make_corr_plot <- function(geneWtPrediction, cellCODEpred, cellType){
  geneWtPrediction <- geneWtPrediction[match(row.names(cellCODEpred), colnames(geneWtPrediction))]
  cellCODEpredCellType <- cellCODEpred[, cellType]
  scatter_data <- cbind(geneWtPrediction, cellCODEpredCellType)
  p <- ggplot(scatter_data, aes(geneWtPrediction, cellCODEpredCellType)) + geom_point() +
    xlab(paste("Estimating", cellType, "proportions using gene weight coeff")) +
    ylab(paste("CellCODE", cellType, "estimate")) +
    theme_bw(base_size = 20) +
    geom_smooth(method = 'lm', colour = "red") +
    stat_cor(method = "spearman") +
    theme(axis.text = element_text(size = 11),
          axis.title = element_text(size = 11),
          plot.title = element_text(size = 14, hjust = 0.5))
  return(p)
}

# Loading data ------------------------------------------------------------

u <- read.csv(paste(path_results, "BM44", "U_SVD_CellCODE.csv", sep = ""),
              stringsAsFactors = FALSE)
cellTypes <- c("ast", "mic", "neu", "oli")

# Checking if BM44 gene coefficients can predict BM44 CellCODE est --------

gExp <- read_scale_HGNC_BM("BM44")
genes <- colnames(gExp)
cellcodeEst <- read_cellCODE_est("BM44")

prediction <- predict_ctf(cellTypes[1], gExp, u)
p1 <- make_corr_plot(prediction, cellcodeEst, cellTypes[1])

prediction <- predict_ctf(cellTypes[2], gExp, u)
p2 <- make_corr_plot(prediction, cellcodeEst, cellTypes[2])

prediction <- predict_ctf(cellTypes[3], gExp, u)
p3 <- make_corr_plot(prediction, cellcodeEst, cellTypes[3])

prediction <- predict_ctf(cellTypes[4], gExp, u)
p4 <- make_corr_plot(prediction, cellcodeEst, cellTypes[4])

figure <- ggarrange(p1, p2, p3, p4, nrow = 2, ncol = 2)
title <- expression(atop(bold("Predicting BM44 ctfs using BM44 CellCODE gene coefficients")))
annotate_figure(figure, top=text_grob(title))
ggsave(paste(path_results, "pred_bm44ctfs_using_bm44genCoef.png", sep = ""))

# Using BM44 model to predict BM36 ctfs -----------------------------------

gExp <- read_scale_HGNC_BM("BM36")
genes <- colnames(gExp)
cellcodeEst <- read_cellCODE_est("BM36")

prediction1 <- predict_ctf(cellTypes[1], gExp, u)
p1 <- make_corr_plot(prediction1, cellcodeEst, cellTypes[1])

prediction2 <- predict_ctf(cellTypes[2], gExp, u)
p2 <- make_corr_plot(prediction2, cellcodeEst, cellTypes[2])

prediction3 <- predict_ctf(cellTypes[3], gExp, u)
p3 <- make_corr_plot(prediction3, cellcodeEst, cellTypes[3])

prediction4 <- predict_ctf(cellTypes[4], gExp, u)
p4 <- make_corr_plot(prediction4, cellcodeEst, cellTypes[4])

stopifnot(colnames(prediction1) == colnames(prediction2),
          colnames(prediction3) == colnames(prediction4))
prediction <- cbind(matrix(prediction1),
                    matrix(prediction2),
                    matrix(prediction3),
                    matrix(prediction4))
colnames(prediction) <- c("ast", "mic", "neu", "oli")
row.names(prediction) <- colnames(prediction1)
head(prediction)

write.csv(data.frame(prediction), paste(path_results,
                                        "BM36_CellCODE_80mg_4cellTypes_using_BM44geneCoef.csv",
                                        sep = ""))
  
figure <- ggarrange(p1, p2, p3, p4, nrow = 2, ncol = 2)
title <- expression(atop(bold("Predicting BM36 ctfs using BM44 CellCODE gene coefficients")))
annotate_figure(figure, top=text_grob(title))
ggsave(paste(path_results, "pred_bm36ctfs_using_bm44genCoef.png", sep = ""))

# Using BM44 model to predict BM22 ctfs -----------------------------------

gExp <- read_scale_HGNC_BM("BM22")
genes <- colnames(gExp)
cellcodeEst <- read_cellCODE_est("BM22")

prediction <- predict_ctf(cellTypes[1], gExp, u)
p1 <- make_corr_plot(prediction, cellcodeEst, cellTypes[1])

prediction <- predict_ctf(cellTypes[2], gExp, u)
p2 <- make_corr_plot(prediction, cellcodeEst, cellTypes[2])

prediction <- predict_ctf(cellTypes[3], gExp, u)
p3 <- make_corr_plot(prediction, cellcodeEst, cellTypes[3])

prediction <- predict_ctf(cellTypes[4], gExp, u)
p4 <- make_corr_plot(prediction, cellcodeEst, cellTypes[4])

figure <- ggarrange(p1, p2, p3, p4, nrow = 2, ncol = 2)
title <- expression(atop(bold("Predicting BM22 ctfs using BM44 CellCODE gene coefficients")))
annotate_figure(figure, top=text_grob(title))
ggsave(paste(path_results, "pred_bm22ctfs_using_bm44genCoef.png", sep = ""))

# Using BM44 model to predict BM10 ctfs -----------------------------------

gExp <- read_scale_HGNC_BM("BM10")
genes <- colnames(gExp)
cellcodeEst <- read_cellCODE_est("BM10")

prediction <- predict_ctf(cellTypes[1], gExp, u)
p1 <- make_corr_plot(prediction, cellcodeEst, cellTypes[1])

prediction <- predict_ctf(cellTypes[2], gExp, u)
p2 <- make_corr_plot(prediction, cellcodeEst, cellTypes[2])

prediction <- predict_ctf(cellTypes[3], gExp, u)
p3 <- make_corr_plot(prediction, cellcodeEst, cellTypes[3])

prediction <- predict_ctf(cellTypes[4], gExp, u)
p4 <- make_corr_plot(prediction, cellcodeEst, cellTypes[4])

figure <- ggarrange(p1, p2, p3, p4, nrow = 2, ncol = 2)
title <- expression(atop(bold("Predicting BM10 ctfs using BM44 CellCODE gene coefficients")))
annotate_figure(figure, top=text_grob(title))
ggsave(paste(path_results, "pred_bm10ctfs_using_bm44genCoef.png", sep = ""))

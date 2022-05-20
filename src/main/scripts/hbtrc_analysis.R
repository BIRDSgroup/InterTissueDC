# Cell-type gene expression correction using a linear model

rm(list = ls())

# Loading necessary packages ----------------------------------------------

setwd("F:/Brain-Disease-networks-version-2/")
.libPaths(c("libraries/", .libPaths()))

# Global variables --------------------------------------------------------

path_geneExp <- "data/hbtrc/"
path_data <- "data/"
path_results <- "results/"
path_geneExp <- "data/hbtrc/"

BR <- c("BA9", "BA17")
fdr_val <- 0.01

file_name <- paste("Edges",toString(fdr_val),"ct_corrected", sep = "")

# Installing libraries ----------------------------------------------------

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

# Utils -------------------------------------------------------------------

# Given list of corr vals r, return p-vals ----------------------------------------

pvalfromtcor <- function(r, n)
{
  tstat = r * sqrt( (n-2) / (1 - r^2) )
  p = 2 * pt(-abs(tstat), n-2) # 2* because two tailed test
  rm(tstat) #reduce memory footprint
  stopifnot(length(p) == length(r))
  p[abs(r)==1] = 0
  return(p)
}

# Zip fastener ------------------------------------------------------------

# Zip fastener(x,y)
# col1 is col1 of x
# col2 is col1 of y
# and so on

zipFastener <- function(df1, df2, along=2)
{
  if(!is.element(along, c(1,2))){
    stop("along must be 1 or 2 for rows and columns
         respectively")
  }
  if(along==1 & (ncol(df1)!= ncol(df2))) {
    stop ("the no. of columns has to be equal to merge
          them by zip feeding")
  }
  if(along==2 & (nrow(df1)!= nrow(df2))) {
    stop ("the no. of rows has to be equal to merge them by
          zip feeding")
  }
  
  d1 <- dim(df1)[along]
  d2 <- dim(df2)[along]
  i1 <- 1:d1           # index vector 1
  i2 <- 1:d2 + d1      # index vector 2
  
  if(d1==d2) {
    dMax <- d1
  } else if (d1 > d2) {
    length(i2) <- length(i1)    # make vectors same length,
    dMax <- d1                  # fill blanks with NAs
  } else  if(d1 < d2){
    length(i1) <- length(i2)    # make vectors same length,
    dMax <- d2                  # fill blanks with NAs
  }
  
  index <- as.vector(matrix(c(i1, i2), ncol=dMax, byrow=T))
  index <- index[!is.na(index)]         # remove NAs
  
  if(along==1){
    colnames(df2) <- colnames(df1)   # keep 1st colnames
    res <- rbind(df1,df2)[ index, ]  # reorder data frame
  }
  if(along==2) res <- cbind(df1,df2)[ , index]
  
  return(res)
}

# Read gene expression data -----------------------------------------------

readBMfile <- function(BR){
  file <- read.csv(paste(path_geneExp, BR,'.csv',sep=""),stringsAsFactors = FALSE)
  row.names(file) <- file[,1]
  return(as.data.frame(t(file[,-1])))
}

# Getting the CELLCODE cell-type estimates --------------------------------

svdIntGrp_U <- function (dat, grp, cutoff = 0.3) 
{
  pp = f.pvalue(dat, grp)
  datr = resid(dat, grp)
  iismall = order(pp)[1:floor(nrow(dat) * cutoff)]
  for (i in 1:2) {
    svdres = svds(dat[-iismall, ])
    sv = svdres$v[, 1, drop = F]
    mod = cbind(grp, sv, sv * grp[, 2])
    mod0 = cbind(grp[, 1], sv)
    pp = f.pvalue(dat[, ], mod, mod0)
    iismall = order(pp)[1:floor(nrow(dat) * cutoff)]
  }
  svdres = svds(dat[-iismall, ])
  geneWeights <- rep(0, nrow(dat))
  geneWeights[-iismall] = svdres$u[,1]
  return(geneWeights)
}

mean_center <- function(data){
  mm=apply(data,1,mean)
  tmp=sweep(data,1,mm, "-")
  return(tmp)
}

# Difference between getALLSPVs and my_getALLSPVs - row.names(SPVs)=colnames(data)
my_getAllSPVs <-
  function(brainRegion = "", data, grp, dataTag, method=c("mixed", "raw", "residual", "SVA"), plot=F, mix.par=0.3){
    cm=intersect(rownames(data), rownames(dataTag)) # data - genes X samples, dataTag - genes X cellTypes
    dataTag=dataTag[cm,, drop=F]
    if (is.null(dim(grp))) {
      grp = model.matrix(~1 + grp)
      # Intercept 1, grpAD - reference, so second column name would be grpCTL an indicator variable
    }
    if(is.null(method)){
      method="mixed"
    }
    method=match.arg(method)
    
    
    #build all SPVs
    SPVs=matrix(nrow=ncol(data), ncol=0) # individuals  X #CellTypes
    
    # Following removes 
    vv=apply(data,1,var) # row-wise i.e per gene variance
    iiuse=which(vv[cm]>0)
    # Removing genes which has zero variance from being considered as marker gene
    dataTag=dataTag[iiuse,,drop=F]
    
    # ncol(dataTag) = #CellTypes
    for (i in 1:ncol(dataTag)) {
      
      # Choosing the genes that have been supplied as 
      # marker for one of the cell type
      genes=rownames(dataTag)[(which(dataTag[,i]!=0))]
      
      if (method=="residual"){
        datar=resid(data,grp)	
        svdres=svds(datar[genes,])
        sv=svdres$v[,1]
      }
      else if (method=="mixed"){     
        
        sv=svdIntGrp(data[genes,], grp, cutoff=mix.par)
        su=svdIntGrp_U(data[genes,], grp, cutoff=mix.par)
        su=cbind(su, genes, colnames(dataTag)[i])
        ppcell=f.pvalue(t(sv), grp)  
        message(paste(colnames(dataTag)[i],"f p.value=",ppcell))
        if(min (ppcell)<0.05){
          message(paste("cell proportions may not be constant in ", colnames(dataTag)[i] )) 
        }
        
      }
      else if (method=="SVA"){
        
        svdres=irwsva.build(data[genes,], grp, cbind(rep(1,ncol(data))),n.sv=1, B=5)
        sv=svdres$sv
      }
      
      else if (method=="raw"){
        svdres=svds(data[genes,])
        sv=svdres$v[,1]
        
        ppcell=f.pvalue(t(sv), grp)  
        message(paste(colnames(dataTag)[i],"f p.value=",ppcell))
      }
      
      cc=cor(sv, t(data[genes,]))
      
      if (mean(cc)<0){
        sv=-sv
      }
      
      SPVs=cbind(SPVs, sv)
      if(i == 1){
        Us=su
      }
      else{
        Us=rbind(Us,su)
      }
    }
    
    if(brainRegion != ""){write.csv(Us, paste(path_results, brainRegion, "U_SVD_CellCODE.csv", sep = ""))}
    
    if (plot){
      parlast=par(no.readonly=T)
      datatmp=rbind(data[rownames(dataTag),], t(SPVs))
      datatmp=resid(datatmp, grp)
      labGenes=c(apply(dataTag, 1, function(x){which(x>0)[1]}), rep(ncol(dataTag)+1, ncol(SPVs)))
      mycol=c(rainbow(ncol(dataTag))[sample(ncol(dataTag))], "black")
      
      corplot(datatmp, usecordist=T, ColSideColors=mycol[labGenes], labRow=c(rep("", nrow(dataTag)), colnames(dataTag)), density.info="none", cexRow=1, dendrogram="column", key=F)
      par(mfg=c(1,par()$mfcol[2]))
      legend("topright", fill=mycol, legend=colnames(dataTag), ncol=1, xpd=NA, cex=1, inset=-0.1)
      par(parlast)
    }
    colnames(SPVs)=colnames(dataTag)
    row.names(SPVs)=colnames(data)
    SPVs
  }

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
  file <- read.csv("data/hbtrc/Sample_information.csv", stringsAsFactors = FALSE)
  stopifnot(row.names(BR) == file$ï..Samples)
  BR$sample_gp <- file$case_control
  BR_AD_Ctl <- BR[c(which(BR$sample_gp == "control"), which(BR$sample_gp == "affected")),]
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
  gene_list <- colnames(gExp) <- gene_temp_names
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

png(filename="images/BA9_CellCODE_80mg_4ct_heatmap.png", width = 800, height = 800)
BA9_SPV <- get_SPV_BrainRegion("BA9", MarkerGenes = mg)
dev.off()
write.csv(BA9_SPV,
          paste(path_results, "BA9_CellCODE_80mg_4cellTypes.csv", sep = ""))

png(filename="images/BA17_CellCODE_80mg_4ct_heatmap.png", width = 800, height = 800)
BA17_SPV <- get_SPV_BrainRegion("BA17", MarkerGenes = mg)
dev.off()
write.csv(BA17_SPV,
          paste(path_results, "BA17_CellCODE_80mg_4cellTypes.csv", sep = ""))

# Loading AD, CTL gene files to be cell-type corrected -------------------

BA9_geneExp <- read_BR_AD_CTL("BA9")
BA9_sample_gp <- BA9_geneExp$sample_gp
BA9_geneExp <- subset(BA9_geneExp, select = -c(sample_gp))

BA17_geneExp <- read_BR_AD_CTL("BA17")
BA17_sample_gp <- BA17_geneExp$sample_gp
BA17_geneExp <- subset(BA17_geneExp, select = -c(sample_gp))

# Changing gene names -----------------------------------------------------

stopifnot(length(unique(colnames(BA9_geneExp))) == length(colnames(BA9_geneExp)))
stopifnot(length(unique(colnames(BA17_geneExp))) == length(colnames(BA17_geneExp)))
stopifnot(colnames(BA17_geneExp) == colnames(BA9_geneExp))

gene_temp_names <- unlist(lapply(1:length(colnames(BA9_geneExp)), 
                                 function(s){c(paste(c("gene", s), collapse = ""))}))
gene_list <- cbind(colnames(BA9_geneExp), gene_temp_names)

# Cell-type correction on whole matrix ------------------------------------

BA9_pred <- predict_gExp(BA9_geneExp, BA9_SPV)
BA17_pred <- predict_gExp(BA17_geneExp, BA17_SPV)

# Linear regression model - cell type correction --------------------------

BA9_corrected <- BA9_geneExp - t(BA9_pred)
BA17_corrected <- BA17_geneExp - t(BA17_pred)

# Writing the files -------------------------------------------------------

write.csv(BA9_corrected, paste(path_results, "BA9_AD_CTL_cellTypeCorrected.csv", sep = ""))
write.csv(BA17_corrected, paste(path_results, "BA17_AD_CTL_cellTypeCorrected.csv", sep = ""))

mdata <- read.csv("data/hbtrc/Sample_information.csv", stringsAsFactors = FALSE)

BA9_corrected_ad <- BA9_corrected[which(mdata$case_control[match(row.names(BA9_corrected), mdata[,1])] == "affected"),]
BA9_corrected_ctl <- BA9_corrected[which(mdata$case_control[match(row.names(BA9_corrected), mdata[,1])] == "control"),]

BA17_corrected_ad <- BA17_corrected[which(mdata$case_control[match(row.names(BA17_corrected), mdata[,1])] == "affected"),]
BA17_corrected_ctl <- BA17_corrected[which(mdata$case_control[match(row.names(BA17_corrected), mdata[,1])] == "control"),]

stopifnot(dim(BA17_corrected_ad)[1] + dim(BA17_corrected_ctl)[1] == dim(BA17_corrected)[1])

BA9_917_ad_geneExp_corrected <- BA9_corrected_ad[intersect(row.names(BA9_corrected_ad), row.names(BA17_corrected_ad)),]
BA17_917_ad_geneExp_corrected <- BA17_corrected_ad[intersect(row.names(BA9_corrected_ad), row.names(BA17_corrected_ad)),]
BA9_917_ctl_geneExp_corrected <- BA9_corrected_ctl[intersect(row.names(BA9_corrected_ctl), row.names(BA17_corrected_ctl)),]
BA17_917_ctl_geneExp_corrected <- BA17_corrected_ctl[intersect(row.names(BA9_corrected_ctl), row.names(BA17_corrected_ctl)),]

# Writing the files -------------------------------------------------------

write.csv(BA9_917_ad_geneExp_corrected, paste(path_results, "ADBA9_BA9BA17_cellTypeCorrected.csv", sep = ""))
write.csv(BA17_917_ad_geneExp_corrected, paste(path_results, "ADBA17_BA9BA17_cellTypeCorrected.csv", sep = ""))
write.csv(BA9_917_ctl_geneExp_corrected, paste(path_results, "CTLBA9_BA9BA17_cellTypeCorrected.csv", sep = ""))
write.csv(BA17_917_ctl_geneExp_corrected, paste(path_results, "CTLBA17_BA9BA17_cellTypeCorrected.csv", sep = ""))

# DC code -----------------------------------------------------------------

read_cell_type_corrected_gExp <- function(br, index, condition){
  if(BR[1] != BR[2]){
    df <- read.csv(paste(path_results, condition, BR[index], "_", BR[1], BR[2], "_cellTypeCorrected.csv",
                         sep = ""), stringsAsFactors = FALSE, row.names = 1)
  } else if(BR[1] == BR[2]){
    df <- read.csv(paste(path_results, BR[1], "_AD_CTL_cellTypeCorrected.csv", sep = ""),
                   stringsAsFactors = FALSE, row.names = 1)
    if(condition == "AD"){
      df <- df[row.names(df)[which(find_sample_group(BR[1]) == "Definite AD")],]
    }else if(condition == "CTL"){
      df <- df[row.names(df)[which(find_sample_group(BR[1]) == "Control")],]
    }
  }
  return(as.data.frame(df))
}

ens_genes <- gene_temp_names

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

########### Gene-Gene Correlation analysis
# Rows are from BR1, columns are from BR2
cor_AD12 <- cor(AD1, AD2, method = "spearman")
cor_CTL12 <- cor(CTL1, CTL2, method = "spearman")

stopifnot(colnames(cor_AD12) == rownames(cor_AD12), colnames(cor_AD12) == ens_genes)
stopifnot(colnames(cor_CTL12) == rownames(cor_CTL12), rownames(cor_CTL12) == ens_genes)

# Garbage collection ------------------------------------------------------

gc()
memory.limit(9999999999)
gc() 

# Code continues ----------------------------------------------------------

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

cor_AD12_edges <- data.frame(cor_AD12_edges, stringsAsFactors = FALSE)
cor_CTL12_edges <- data.frame(cor_CTL12_edges, stringsAsFactors = FALSE)

cor_AD12_edges$combined <- paste(cor_AD12_edges$g1.BR1,cor_AD12_edges$g2.BR2)
cor_CTL12_edges$combined <- paste(cor_CTL12_edges$g1.BR1,cor_CTL12_edges$g2.BR2)

write.csv(cor_AD12_edges, paste(path_results,BR[1], BR[2], "CTCorrectedADCorrelatedEdges.csv", sep = ""))
write.csv(cor_CTL12_edges, paste(path_results,BR[1], BR[2], "CTCorrectedCTLCorrelatedEdges.csv", sep = ""))

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

write.csv(edges, paste(path_results,"UnionCorrelatedEdges_ct_correction",BR[1],BR[2], ".csv", sep = ""))

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

rm(t1,t2,t3,t4)

# r.test DC edges ---------------------------------------------------------

# Initializing list of pvalues vector
p <- rep(NA,length(colnames(res_AD12))/2 + 1)
z <- rep(NA, length(p))
cor_CTL <- rep(NA, length(p))
cor_AD <- rep(NA, length(p))

i <- 1
j <- 1
while(i <= length(colnames(res_AD12)))
{
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
save(p, file = paste(path_results,"p_vals_without_adjustment_ct_correction",BR[1],BR[2],toString(file_name),".rda",sep=""))

p12 <- p
p_adjusted <- p.adjust(p12, method = 'BH')

# Print DC edges
dcres <- cbind(edges, p_adjusted, z, cor_CTL, cor_AD)
r <- which(p_adjusted <= fdr_val)
dc_rtest <- dcres[r,]
write.csv(dc_rtest, paste(path_results,"DCEdges_ct_correction",BR[1],BR[2],".csv",sep=""))

cat("Number of significantly correlated edges in AD", dim(cor_AD12_edges)[1], "\n")
cat("Number of significantly correlated edges in CTL", dim(cor_CTL12_edges)[1], "\n")
cat("Number of significantly correlated edges in the union", length(union12), "\n")
cat("Number of DC edges ", dim(dc_rtest)[1], "\n")

rm(list = ls())

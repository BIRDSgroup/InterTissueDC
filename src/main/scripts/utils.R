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
  file <- read.csv(paste(path_geneExp,BR,'.csv',sep=""),stringsAsFactors = FALSE)
  row.names(file) <- file[,1]
  return(as.data.frame(t(file[,-1])))
}

# ENSEMBL to HGNC mapping -------------------------------------------------
ensid_sampleid <- read.csv(paste(path_data, "msbb_metadata_mapping_and_group_assignment.csv",
                                 sep = ""), stringsAsFactors = FALSE)

find_genes_mapped <- function(){
  ens_ids <- ensid_sampleid[,1]
  
  gencode_mapping <- read.csv(paste(path_data, "Mapped_h37_ens_to_gene.csv", sep = ""), 
                              stringsAsFactors = FALSE)
  genes_mapped <- gencode_mapping$gene_name[match(ens_ids, gencode_mapping$gene_id)]
  return(genes_mapped)
}

genes_mapped <- find_genes_mapped()

# Finding sample group ----------------------------------------------------

rna_metadata_processed <- read.csv(paste(path_data, 
                                         "rna_metadata_onlyBAM_noExclude_noDup_Nov2018.csv",
                                         sep = ""), stringsAsFactors = FALSE)

find_individual_id <- function(BR){
  barcodes <- na.omit(sapply(strsplit(ensid_sampleid[,BR],"_"), function(x){x[3]}))
  ind_id <- rna_metadata_processed$individualIdentifier[match(barcodes, rna_metadata_processed$barcode)]
  # Noticed that some barcodes do not have indId mapping
  ind_id[which(is.na(ind_id))] <- 0
  return(ind_id)
}

clinical <- read.csv(paste(path_data, "MSBB_clinical.csv", sep = ""),
                     stringsAsFactors = FALSE)

find_sample_group <- function(BrainRegion){
  barcodes <- na.omit(sapply(strsplit(ensid_sampleid[,BrainRegion],"_"), function(x){x[3]}))
  ind_id <- rna_metadata_processed$individualIdentifier[match(barcodes, rna_metadata_processed$barcode)]
  return(clinical$DxCondition[match(ind_id, clinical$individualIdentifier)])
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
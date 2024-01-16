# To identity DEG in cell-type corrected data
# using Mann Whitney Wilcoxon test (a.k.a Mann Whitney U test and Wilcoxon Rank Sum test)

rm(list = ls())

# Loading global variables ------------------------------------------------

option_list <- list(
  make_option(c('--workingDir'), help = "working directory", default = ""),
  make_option(c('--br'), help = 'Brain region', default = ""),
  make_option(c('--ctc'), help = "either TRUE or FALSE - TRUE if cell-type corrected", default = FALSE)
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

wdir <- opt$workingDir
setwd(wdir)
BR <- opt$br
ct_corrected <- opt$ctc

if(ct_corrected == TRUE){
  ct_file_name = "cellTypeCorrected"
} else{
  ct_file_name = ""
}

source("scripts/utils.R")

# Loading libraries -------------------------------------------------------

.libPaths(c("libraries/", .libPaths()))
source("scripts/libraries.R")

# Loading data files ------------------------------------------------------

if(ct_corrected == TRUE){
  ad_ctl_gExp <- read.csv(paste("results/", BR, "_AD_CTL_cellTypeCorrected.csv", sep = ""),
                          stringsAsFactors = FALSE, row.names = 1)
  colnames(ad_ctl_gExp) <- ensid_sampleid$ENS_ID
} else{
  bm <- read.csv(paste("data/",BR,'.csv',sep=""),stringsAsFactors = FALSE, row.names = 1)
  sample_gp <- find_sample_group(BR)
  ad_ctl_gExp <- bm[,which(sample_gp %in% c("Control", "Definite AD"))]
  ad_ctl_gExp <- t(ad_ctl_gExp)
}

barcodes <- sapply(strsplit(row.names(ad_ctl_gExp),"_"), function(x){x[3]})
mdata <- read.csv(paste(path_data, "rna_metadata_onlyBAM_noExclude_noDup_Nov2018.csv", sep = ""),
                     stringsAsFactors = FALSE)
sample_condition <- mdata$DxCondition[match(barcodes, mdata$barcode)]

# Differential gene expression --------------------------------------------

# wilcoxTest(eset, fac, tests = c("unpaired", "paired"))
# eset - The expression matrix. Each row is a gene, and each column is a subject/sample.
# The gene names must be presented as the row names.
# fac - Subject labels, for unpaired T-test, either a factor or something that can be coerced into a factor 
# (e.g. 0 and 1, Experiment and Control). For paired T-test, fac must be an integer vector of 1,-1,2,-2,...,
# where each number represents a pair, and the sign represents the conditions.
# tests - The tests to performed. Can be either the default "unpaired" for unpaired T-tests or "paired" for paired T-tests.

deg <- wilcoxTest(t(ad_ctl_gExp), sample_condition, tests = "unpaired")
deg <- as.data.frame(deg)

# Adding average expression -----------------------------------------------

deg$AvgExpn <- apply(ad_ctl_gExp, 2, mean)

# Calculating logFC -------------------------------------------------------

# logFC = log((geneAD - geneCTL)/(geneCTL))

ad_gExp <- ad_ctl_gExp[which(sample_condition == "Definite AD"),]
ctl_gExp <- ad_ctl_gExp[which(sample_condition == "Control"),]

ad_avg <- apply(ad_gExp, 2, mean)
ctl_avg <- apply(ctl_gExp, 2, mean)

FC <- ad_avg - ctl_avg
FC["ENSG00000080709"]

deg$logFC <- FC

# Adding mapped gene symbol -----------------------------------------------

mapping_file <- read.csv(paste("data/Mapped_h37_ens_to_gene.csv", sep = ""),
                         stringsAsFactors = FALSE)
deg$gene <- mapping_file$gene_name[match(row.names(deg), mapping_file$gene_id)]

write.csv(deg, paste("results/", BR, ct_file_name, "_DEG_WilcoxonRankSumTest.csv", sep = ""))

# Author: B P Kailash
# Purpose: Analyze the DC relations observed - Part 1

# 1) Find the sample size of the analysis
# 2) Convert all DC edges from ensids to HGNC symbols
# 3) Add M, NM, cell type labels to each DC edge

rm(list = ls())

# Global variables --------------------------------------------------------

source("scripts/libraries.R")
source("scripts/utils.R")

library(purrr)
library(BRETIGEA)

mapping_file <- read.csv(paste("data/", 
                               "Mapped_h37_ens_to_gene.csv", 
                               sep = ""), 
                         stringsAsFactors = FALSE)

cell_types <- c("ast", "mic", "neu", "oli")
number_mg <- 1000
number_mg_per_ct <- rep(number_mg, length(cell_types))

# Finding total sample size -----------------------------------------------

colnames(ensid_sampleid)

bm10_samples <- ensid_sampleid$BM10[which(ensid_sampleid$BM10!= "")]
bm22_samples <- ensid_sampleid$BM22[which(ensid_sampleid$BM22!= "")]
bm36_samples <- ensid_sampleid$BM36[which(ensid_sampleid$BM36!= "")]
bm44_samples <- ensid_sampleid$BM44[which(ensid_sampleid$BM44!= "")]

get_individual_id <- function(samples){
  return(rna_metadata_processed$individualIdentifier[
    match(unlist(sapply(samples, function(x){unlist(strsplit(x, "_"))[3]})),
          rna_metadata_processed$barcode)]
  )
}

all_samples <- unique(c(get_individual_id(bm10_samples),
                        get_individual_id(bm22_samples),
                        get_individual_id(bm36_samples),
                        get_individual_id(bm44_samples)))
Dx <- rna_metadata_processed$DxCondition[match(all_samples, rna_metadata_processed$individualIdentifier)]

paste("Number of CTL individuals in this study: ", length(which(Dx == "Control")))
paste("Number of AD individuals in this study: ", length(which(Dx == "Definite AD")))

# Finding individual brain region sample size -----------------------------

get_sample_dx_class <- function(BrainRegion, brain_region_samples){
  samples <- unlist(map(str_split(brain_region_samples, "\\."), 2))
  sample_dx_class <- rna_metadata_processed$DxCondition[match(samples, rna_metadata_processed$sampleIdentifier)]
  cat(paste("Number of AD samples for ", BrainRegion, " is: ", length(which(sample_dx_class == "Definite AD")), "\n"))
  cat(paste("Number of Control samples for ", BrainRegion, " is: ", length(which(sample_dx_class == "Control")), "\n"))
}

get_sample_dx_class("BM10", bm10_samples)
get_sample_dx_class("BM22", bm22_samples)
get_sample_dx_class("BM36", bm36_samples)
get_sample_dx_class("BM44", bm44_samples)

# Finding paired brain-region sample sizes --------------------------------

get_paired_sample_size <- function(BR){
  ad_ids_BR1_2 <- read.csv(paste("results/", "AD_ids", BR[1], BR[2], ".csv", sep = ""),
                           stringsAsFactors = FALSE)
  ctl_ids_BR1_2 <- read.csv(paste("results/", "CTL_ids", BR[1], BR[2], ".csv", sep = ""),
                            stringsAsFactors = FALSE)
  cat(paste("Number of AD samples for ", BR[1], "-", BR[2], " is : ", length(ad_ids_BR1_2$x), "\n"))
  cat(paste("Number of CTL samples for ", BR[1], "-", BR[2], " is : ", length(ctl_ids_BR1_2$x)))
}

get_paired_sample_size(c("BM_10", "BM_22"))
get_paired_sample_size(c("BM_10", "BM_36"))
get_paired_sample_size(c("BM_10", "BM_44"))
get_paired_sample_size(c("BM_22", "BM_36"))
get_paired_sample_size(c("BM_22", "BM_44"))
get_paired_sample_size(c("BM_36", "BM_44"))

get_paired_sample_size(c("BM_10", "BM_10"))
get_paired_sample_size(c("BM_22", "BM_22"))
get_paired_sample_size(c("BM_36", "BM_36"))
get_paired_sample_size(c("BM_44", "BM_44"))

# Convert significant correlation edges to HGNC mapping ------------------

ensembl_to_approved_gene_symbol_dc_edges <- function(MappingFile, PathResults = "results/", BR, ct_corrected){
  
  if(BR[1] == BR[2]){
    analysis_type <- "within_region/"
  } else if(BR[1] != BR[2]){
    analysis_type <- "between_region/"
  }
  
  if(ct_corrected == TRUE){
    ct_file_name <- "_ct_correction"
    DCEdges <- read.csv(paste("results/", "after_ctc_", analysis_type, "DCEdges", ct_file_name, BR[1],BR[2], ".csv", sep = ""),
                                        stringsAsFactors = FALSE)
    outputFileName <- paste(PathResults, "after_ctc_", analysis_type, BR[1],BR[2], ct_file_name, "Mapped_DCEdgesWithoutNA", ".csv",sep="")
  } else {
    ct_file_name = ""
    DCEdges <- read.csv(paste("results/", "before_ctc_", analysis_type, "DCEdges", ct_file_name, BR[1],BR[2], ".csv", sep = ""),
                                        stringsAsFactors = FALSE)
    outputFileName <- paste(PathResults, "before_ctc_", analysis_type, BR[1],BR[2], ct_file_name, "Mapped_DCEdgesWithoutNA", ".csv",sep="")
  }
  
  ens1 <- DCEdges$a
  ens2 <- DCEdges$b
  map1 <- match(ens1, MappingFile$gene_id)
  map2 <- match(ens2, MappingFile$gene_id)
  g1 <- MappingFile$gene_name[map1]
  g2 <- MappingFile$gene_name[map2]
  
  geneMapped <- cbind(DCEdges, g1, g2)
  count_before_mapping <- dim(geneMapped)[[1]]
  
  geneMapped <- geneMapped[-union(which(is.na(geneMapped$g1)), which(is.na(geneMapped$g2))),]
  count_after_mapping <- dim(geneMapped)[[1]]
  
  cat(paste("Number of significant DC edges before mapping for: ", BR[1], "-", BR[2], count_before_mapping), sep = "")
  cat(paste("\nNumber of DC edges after mapping for: ", BR[1], "-", BR[2], count_after_mapping), sep = "")
  cat("\nNumber of edges lost to mapping are: ", count_before_mapping - count_after_mapping, "\n")
  write.csv(geneMapped , file = outputFileName, row.names = FALSE)
}

# ********* 1.1

# ct_corrected = FALSE

## Inter-region
ensembl_to_approved_gene_symbol_dc_edges(MappingFile = mapping_file,
                                     PathResults = "results/",
                                     BR = c("BM_10", "BM_22"), ct_corrected = FALSE)

ensembl_to_approved_gene_symbol_dc_edges(MappingFile = mapping_file,
                                     PathResults = "results/",
                                     BR = c("BM_10", "BM_36"), ct_corrected = FALSE)

ensembl_to_approved_gene_symbol_dc_edges(MappingFile = mapping_file,
                                     PathResults = "results/",
                                     BR = c("BM_10", "BM_44"), ct_corrected = FALSE)

ensembl_to_approved_gene_symbol_dc_edges(MappingFile = mapping_file,
                                     PathResults = "results/",
                                     BR = c("BM_22", "BM_36"), ct_corrected = FALSE)

ensembl_to_approved_gene_symbol_dc_edges(MappingFile = mapping_file,
                                     PathResults = "results/",
                                     BR = c("BM_22", "BM_44"), ct_corrected = FALSE)

ensembl_to_approved_gene_symbol_dc_edges(MappingFile = mapping_file,
                                     PathResults = "results/",
                                     BR = c("BM_36", "BM_44"), ct_corrected = FALSE)

# ********* 1.2
# ct_corrected = TRUE

## Inter-region
ensembl_to_approved_gene_symbol_dc_edges(MappingFile = mapping_file,
                                     PathResults = "results/",
                                     BR = c("BM10", "BM22"), ct_corrected = TRUE)

ensembl_to_approved_gene_symbol_dc_edges(MappingFile = mapping_file,
                                     PathResults = "results/",
                                     BR = c("BM10", "BM36"), ct_corrected = TRUE)

ensembl_to_approved_gene_symbol_dc_edges(MappingFile = mapping_file,
                                     PathResults = "results/",
                                     BR = c("BM10", "BM44"), ct_corrected = TRUE)

ensembl_to_approved_gene_symbol_dc_edges(MappingFile = mapping_file,
                                     PathResults = "results/",
                                     BR = c("BM22", "BM36"), ct_corrected = TRUE)

ensembl_to_approved_gene_symbol_dc_edges(MappingFile = mapping_file,
                                     PathResults = "results/",
                                     BR = c("BM22", "BM44"), ct_corrected = TRUE)

ensembl_to_approved_gene_symbol_dc_edges(MappingFile = mapping_file,
                                     PathResults = "results/",
                                     BR = c("BM36", "BM44"), ct_corrected = TRUE)

# Labeling significantly correlated edges ---------------------------------

cell_types <- c("ast", "mic", "neu", "oli")
number_mg <- 1000
number_mg_per_ct <- rep(number_mg, length(cell_types))

label_dc_edges <- function(CellTypes, NumberOfMarkersPerCellType, BR, CTCorrection){
  
  if(BR[1] == BR[2]){
    analysis_type <- "within_region/"
  } else if(BR[1] != BR[2]){
    analysis_type <- "between_region/"
  }
  
  if(CTCorrection == TRUE){
    ct_file_name <- "_ct_correction"
    inputFileName <- paste("results/", "after_ctc_", analysis_type, BR[1], BR[2], ct_file_name, "Mapped_DCEdgesWithoutNA", ".csv", sep = "")
    outPutFileName <- paste("results/", "after_ctc_", analysis_type, BR[1],BR[2], ct_file_name, "Mapped_DCEdgesWithoutNA", ".csv",sep="")
  } else {
    ct_file_name = ""
    inputFileName <- paste("results/", "before_ctc_", analysis_type, BR[1], BR[2], ct_file_name, "Mapped_DCEdgesWithoutNA", ".csv", sep = "")
    outPutFileName <- paste("results/", "before_ctc_", analysis_type, BR[1],BR[2], ct_file_name, "Mapped_DCEdgesWithoutNA", ".csv",sep="")
  }
  
  if(length(CellTypes) != length(NumberOfMarkersPerCellType)){
    print("Length of the cell types and the number of markers list is not the same")
    break
  }
  
  mg <- c()
  ct <- c()
  
  for(i in 1:length(CellTypes)){
    ct_markers <- head(markers_df_human_brain$markers[which(markers_df_brain$cell == CellTypes[i])],
                       n = NumberOfMarkersPerCellType[i])
    ct <- c(ct, rep(CellTypes[i], NumberOfMarkersPerCellType[i]))
    mg <- c(mg, ct_markers)
  }
  
  correlated_edges <- read.csv(inputFileName, stringsAsFactors = FALSE)
  
  marker_label_g1 <- rep(0, dim(correlated_edges)[1])
  marker_label_g2 <- rep(0, dim(correlated_edges)[1])
  
  match1 <- match(toupper(correlated_edges$g1), mg)
  match2 <- match(toupper(correlated_edges$g2), mg)
  
  marker_label_g1[which(is.na(match1))] = "NM"
  marker_label_g1[which(marker_label_g1 == 0)] = "M"
  marker_label_g2[which(is.na(match2))] = "NM"
  marker_label_g2[which(marker_label_g2 == 0)] = "M"
  
  correlated_edges$m1 <- marker_label_g1
  correlated_edges$m2 <- marker_label_g2
  
  correlated_edges$c1 <- ct[match1]
  correlated_edges$c2 <- ct[match2]
  
  write.csv(correlated_edges, outPutFileName)
}

# ********* 3.1
# ct_corrected = FALSE

## Inter-region

label_dc_edges(BR = c("BM_10",  "BM_22"), CellTypes = cell_types, NumberOfMarkersPerCellType = number_mg_per_ct, CTCorrection = FALSE)
label_dc_edges(BR = c("BM_10",  "BM_36"), CellTypes = cell_types, NumberOfMarkersPerCellType = number_mg_per_ct, CTCorrection = FALSE)
label_dc_edges(BR = c("BM_10",  "BM_44"), CellTypes = cell_types, NumberOfMarkersPerCellType = number_mg_per_ct, CTCorrection = FALSE)
label_dc_edges(BR = c("BM_22",  "BM_36"), CellTypes = cell_types, NumberOfMarkersPerCellType = number_mg_per_ct, CTCorrection = FALSE)
label_dc_edges(BR = c("BM_22",  "BM_44"), CellTypes = cell_types, NumberOfMarkersPerCellType = number_mg_per_ct, CTCorrection = FALSE)
label_dc_edges(BR = c("BM_36",  "BM_44"), CellTypes = cell_types, NumberOfMarkersPerCellType = number_mg_per_ct, CTCorrection = FALSE)

# ********* 3.2
# ct_corrected = TRUE

## Inter-region

label_dc_edges(BR = c("BM10",  "BM22"), CellTypes = cell_types, NumberOfMarkersPerCellType = number_mg_per_ct, CTCorrection = TRUE)
label_dc_edges(BR = c("BM10",  "BM36"), CellTypes = cell_types, NumberOfMarkersPerCellType = number_mg_per_ct, CTCorrection = TRUE)
label_dc_edges(BR = c("BM10",  "BM44"), CellTypes = cell_types, NumberOfMarkersPerCellType = number_mg_per_ct, CTCorrection = TRUE)
label_dc_edges(BR = c("BM22",  "BM36"), CellTypes = cell_types, NumberOfMarkersPerCellType = number_mg_per_ct, CTCorrection = TRUE)
label_dc_edges(BR = c("BM22",  "BM44"), CellTypes = cell_types, NumberOfMarkersPerCellType = number_mg_per_ct, CTCorrection = TRUE)
label_dc_edges(BR = c("BM36",  "BM44"), CellTypes = cell_types, NumberOfMarkersPerCellType = number_mg_per_ct, CTCorrection = TRUE)

# Adding Chr and Gene_type labels to correlated edges ---------------------

add_chr_geneType_dc <- function(BR, CTCorrection){
  
  if(BR[1] == BR[2]){
    analysis_type <- "within_region/"
  } else if(BR[1] != BR[2]){
    analysis_type <- "between_region/"
  }
  
  if(CTCorrection == TRUE){
    ct_file_name <- "_ct_correction"
    inputFileName <- paste("results/", "after_ctc_", analysis_type, BR[1], BR[2], ct_file_name, "Mapped_DCEdgesWithoutNA", ".csv", sep = "")
    outPutFileName <- paste("results/", "after_ctc_", analysis_type, BR[1],BR[2], ct_file_name, "Mapped_DCEdgesWithoutNA", ".csv",sep="")
  } else {
    ct_file_name = ""
    inputFileName <- paste("results/", "before_ctc_", analysis_type, BR[1], BR[2], ct_file_name, "Mapped_DCEdgesWithoutNA", ".csv", sep = "")
    outPutFileName <- paste("results/", "before_ctc_", analysis_type, BR[1],BR[2], ct_file_name, "Mapped_DCEdgesWithoutNA", ".csv",sep="")
  }
  
  correlated_edges <- read.csv(inputFileName, stringsAsFactors = FALSE)
  correlated_edges$gene_type_1 <- mapping_file$gene_type[match(correlated_edges$g1, mapping_file$gene_name)]
  correlated_edges$gene_type_2 <- mapping_file$gene_type[match(correlated_edges$g2, mapping_file$gene_name)]
  correlated_edges$chr1 <- mapping_file$chr[match(correlated_edges$g1, mapping_file$gene_name)]
  correlated_edges$chr2 <- mapping_file$chr[match(correlated_edges$g2, mapping_file$gene_name)]
  write.csv(correlated_edges, outPutFileName)
}

# ********* 5.1
# ct_corrected = FALSE

## Inter-region

add_chr_geneType_dc(BR = c("BM_10",  "BM_22"), CTCorrection = FALSE)
add_chr_geneType_dc(BR = c("BM_10",  "BM_36"), CTCorrection = FALSE)
add_chr_geneType_dc(BR = c("BM_10",  "BM_44"), CTCorrection = FALSE)
add_chr_geneType_dc(BR = c("BM_22",  "BM_36"), CTCorrection = FALSE)
add_chr_geneType_dc(BR = c("BM_22",  "BM_44"), CTCorrection = FALSE)
add_chr_geneType_dc(BR = c("BM_36",  "BM_44"), CTCorrection = FALSE)

# ********* 5.2
# ct_corrected = TRUE

## Inter-region

add_chr_geneType_dc(BR = c("BM10",  "BM22"), CTCorrection = TRUE)
add_chr_geneType_dc(BR = c("BM10",  "BM36"), CTCorrection = TRUE)
add_chr_geneType_dc(BR = c("BM10",  "BM44"), CTCorrection = TRUE)
add_chr_geneType_dc(BR = c("BM22",  "BM36"), CTCorrection = TRUE)
add_chr_geneType_dc(BR = c("BM22",  "BM44"), CTCorrection = TRUE)
add_chr_geneType_dc(BR = c("BM36",  "BM44"), CTCorrection = TRUE)

######################################################################################################

## Within-region

ensembl_to_approved_gene_symbol_dc_edges(MappingFile = mapping_file,
                                         PathResults = "results/",
                                         BR = c("BM_10", "BM_10"), ct_corrected = FALSE)

ensembl_to_approved_gene_symbol_dc_edges(MappingFile = mapping_file,
                                         PathResults = "results/",
                                         BR = c("BM_22", "BM_22"), ct_corrected = FALSE)

ensembl_to_approved_gene_symbol_dc_edges(MappingFile = mapping_file,
                                         PathResults = "results/",
                                         BR = c("BM_36", "BM_36"), ct_corrected = FALSE)

ensembl_to_approved_gene_symbol_dc_edges(MappingFile = mapping_file,
                                         PathResults = "results/",
                                         BR = c("BM_44", "BM_44"), ct_corrected = FALSE)

## Within-region

ensembl_to_approved_gene_symbol_dc_edges(MappingFile = mapping_file,
                                         PathResults = "results/",
                                         BR = c("BM10", "BM10"), ct_corrected = TRUE)

ensembl_to_approved_gene_symbol_dc_edges(MappingFile = mapping_file,
                                         PathResults = "results/",
                                         BR = c("BM22", "BM22"), ct_corrected = TRUE)

ensembl_to_approved_gene_symbol_dc_edges(MappingFile = mapping_file,
                                         PathResults = "results/",
                                         BR = c("BM36", "BM36"), ct_corrected = TRUE)
ensembl_to_approved_gene_symbol_dc_edges(MappingFile = mapping_file,
                                         PathResults = "results/",
                                         BR = c("BM44", "BM44"), ct_corrected = TRUE)

## Within-region

label_dc_edges(BR = c("BM_10",  "BM_10"), CellTypes = cell_types, NumberOfMarkersPerCellType = number_mg_per_ct, CTCorrection = FALSE)
label_dc_edges(BR = c("BM_22",  "BM_22"), CellTypes = cell_types, NumberOfMarkersPerCellType = number_mg_per_ct, CTCorrection = FALSE)
label_dc_edges(BR = c("BM_36",  "BM_36"), CellTypes = cell_types, NumberOfMarkersPerCellType = number_mg_per_ct, CTCorrection = FALSE)
label_dc_edges(BR = c("BM_44",  "BM_44"), CellTypes = cell_types, NumberOfMarkersPerCellType = number_mg_per_ct, CTCorrection = FALSE)

## Within-region

label_dc_edges(BR = c("BM10",  "BM10"), CellTypes = cell_types, NumberOfMarkersPerCellType = number_mg_per_ct, CTCorrection = TRUE)
label_dc_edges(BR = c("BM22",  "BM22"), CellTypes = cell_types, NumberOfMarkersPerCellType = number_mg_per_ct, CTCorrection = TRUE)
label_dc_edges(BR = c("BM36",  "BM36"), CellTypes = cell_types, NumberOfMarkersPerCellType = number_mg_per_ct, CTCorrection = TRUE)
label_dc_edges(BR = c("BM44",  "BM44"), CellTypes = cell_types, NumberOfMarkersPerCellType = number_mg_per_ct, CTCorrection = TRUE)

## Within-region

add_chr_geneType_dc(BR = c("BM_10",  "BM_10"), CTCorrection = FALSE)
add_chr_geneType_dc(BR = c("BM_22",  "BM_22"), CTCorrection = FALSE)
add_chr_geneType_dc(BR = c("BM_36",  "BM_36"), CTCorrection = FALSE)
add_chr_geneType_dc(BR = c("BM_44",  "BM_44"), CTCorrection = FALSE)

## Within-region

add_chr_geneType_dc(BR = c("BM10",  "BM10"), CTCorrection = TRUE)
add_chr_geneType_dc(BR = c("BM22",  "BM22"), CTCorrection = TRUE)
add_chr_geneType_dc(BR = c("BM36",  "BM36"), CTCorrection = TRUE)
add_chr_geneType_dc(BR = c("BM44",  "BM44"), CTCorrection = TRUE)

######################################################################################################

# Writing significant DC results

write_significant_dc <- function(BR, CTCorrection, fdr_val = 0.01){
  
  if(BR[1] == BR[2]){
    analysis_type <- "within_region/"
  } else if(BR[1] != BR[2]){
    analysis_type <- "between_region/"
  }
  
  if(CTCorrection == TRUE){
    ct_file_name <- "_ct_correction"
    inputFileName <- paste("results/", "after_ctc_", analysis_type, BR[1], BR[2], ct_file_name, "Mapped_DCEdgesWithoutNA", ".csv", sep = "")
    outPutFileName <- paste("results/", "after_ctc_", analysis_type, BR[1],BR[2], ct_file_name, "Mapped_DCEdgesWithoutNA_fdr_", fdr_val, ".csv",sep="")
  } else {
    ct_file_name = ""
    inputFileName <- paste("results/", "before_ctc_", analysis_type, BR[1], BR[2], ct_file_name, "Mapped_DCEdgesWithoutNA", ".csv", sep = "")
    outPutFileName <- paste("results/", "before_ctc_", analysis_type, BR[1],BR[2], ct_file_name, "Mapped_DCEdgesWithoutNA_fdr_", fdr_val, ".csv",sep="")
  }
  
  dcedges <- read.csv(inputFileName, stringsAsFactors = FALSE)
  r <- which(dcedges$p_adjusted <= fdr_val)
  write.csv(dcedges[r,], outPutFileName)
  message("**Number of DC edges ", dim(dcedges[r,])[1], "**")
}

# ********* 7.1
# ct_corrected = FALSE

## Between-region

write_significant_dc(BR = c("BM_10",  "BM_22"), CTCorrection = FALSE)
write_significant_dc(BR = c("BM_10",  "BM_36"), CTCorrection = FALSE)
write_significant_dc(BR = c("BM_10",  "BM_44"), CTCorrection = FALSE)
write_significant_dc(BR = c("BM_22",  "BM_36"), CTCorrection = FALSE)
write_significant_dc(BR = c("BM_22",  "BM_44"), CTCorrection = FALSE)
write_significant_dc(BR = c("BM_36",  "BM_44"), CTCorrection = FALSE)

## Within-region

write_significant_dc(BR = c("BM_10",  "BM_10"), CTCorrection = FALSE)
write_significant_dc(BR = c("BM_22",  "BM_22"), CTCorrection = FALSE)
write_significant_dc(BR = c("BM_36",  "BM_36"), CTCorrection = FALSE)
write_significant_dc(BR = c("BM_44",  "BM_44"), CTCorrection = FALSE)

# ********* 7.2
# ct_corrected = TRUE

## Inter-region

write_significant_dc(BR = c("BM10",  "BM22"), CTCorrection = TRUE)
write_significant_dc(BR = c("BM10",  "BM36"), CTCorrection = TRUE)
write_significant_dc(BR = c("BM10",  "BM44"), CTCorrection = TRUE)
write_significant_dc(BR = c("BM22",  "BM36"), CTCorrection = TRUE)
write_significant_dc(BR = c("BM22",  "BM44"), CTCorrection = TRUE)
write_significant_dc(BR = c("BM36",  "BM44"), CTCorrection = TRUE)

## Within-region

write_significant_dc(BR = c("BM10",  "BM10"), CTCorrection = TRUE)
write_significant_dc(BR = c("BM22",  "BM22"), CTCorrection = TRUE)
write_significant_dc(BR = c("BM36",  "BM36"), CTCorrection = TRUE)
write_significant_dc(BR = c("BM44",  "BM44"), CTCorrection = TRUE)

# Count M-M, NM-NM, M-NM, NM-M interactions -------------------------------

count_interactions <- function(BR, ct_correction, fdr_val){
  
  if(BR[1] == BR[2]){
    analysis_type <- "within_region/"
  } else if(BR[1] != BR[2]){
    analysis_type <- "between_region/"
  }
  
  if(ct_correction == TRUE){
    ct_file_name <- "_ct_correction"
    inputFileName <- paste("results/", "after_ctc_", analysis_type, BR[1],BR[2], ct_file_name, "Mapped_DCEdgesWithoutNA_fdr_", fdr_val, ".csv",sep="")
  } else {
    ct_file_name = ""
    inputFileName <- paste("results/", "before_ctc_", analysis_type, BR[1],BR[2], ct_file_name, "Mapped_DCEdgesWithoutNA_fdr_", fdr_val, ".csv",sep="")
  }
  
  dc <- read.csv(inputFileName, stringsAsFactors = FALSE)
  interaction = paste(dc$m1, dc$m2, sep = "")
  cat(paste("Number of M-M interactions for"), BR[1], BR[2], "is", length(which(interaction == "MM")))
  cat(paste("\nNumber of M-NM interactions is"), BR[1], BR[2], "is", length(which(interaction == "MNM")))
  cat(paste("\nNumber of NM-M interactions is"), BR[1], BR[2], "is",length(which(interaction == "NMM")))
  cat(paste("\nNumber of NM-NM interactions is"), BR[1], BR[2], "is",length(which(interaction == "NMNM")))
}

# NOTE that there are no duplicates among the 4000 markers for the 4 cell types
length(unique(markers_df_human_brain$markers[which(markers_df_human_brain$cell %in% c("ast", "mic", "neu", "oli"))]))

# ********* 7.1
# ct_corrected = FALSE

## Inter-region

count_interactions(c("BM_10", "BM_22"), ct_correction =  FALSE, fdr_val = 0.01)
count_interactions(c("BM_10", "BM_36"), ct_correction =  FALSE, fdr_val = 0.01)
count_interactions(c("BM_10", "BM_44"), ct_correction =  FALSE, fdr_val = 0.01)
count_interactions(c("BM_22", "BM_36"), ct_correction =  FALSE, fdr_val = 0.01)
count_interactions(c("BM_22", "BM_44"), ct_correction =  FALSE, fdr_val = 0.01)
count_interactions(c("BM_36", "BM_44"), ct_correction =  FALSE, fdr_val = 0.01)
         
# ********* 7.2          
# ct_corrected = TRUE

## Inter-region

count_interactions(c("BM10", "BM22"), ct_correction =  TRUE, fdr_val = 0.01)
count_interactions(c("BM10", "BM36"), ct_correction =  TRUE, fdr_val = 0.01)
count_interactions(c("BM10", "BM44"), ct_correction =  TRUE, fdr_val = 0.01)
count_interactions(c("BM22", "BM36"), ct_correction =  TRUE, fdr_val = 0.01)
count_interactions(c("BM22", "BM44"), ct_correction =  TRUE, fdr_val = 0.01)
count_interactions(c("BM36", "BM44"), ct_correction =  TRUE, fdr_val = 0.01)


## Within-region

count_interactions(c("BM_10", "BM_10"), ct_correction =  FALSE, fdr_val = 0.01)
count_interactions(c("BM_22", "BM_22"), ct_correction =  FALSE, fdr_val = 0.01)
count_interactions(c("BM_36", "BM_36"), ct_correction =  FALSE, fdr_val = 0.01)
count_interactions(c("BM_44", "BM_44"), ct_correction =  FALSE, fdr_val = 0.01)

## Within-region

count_interactions(c("BM10", "BM10"), ct_correction =  TRUE, fdr_val = 0.01)
count_interactions(c("BM22", "BM22"), ct_correction =  TRUE, fdr_val = 0.01)
count_interactions(c("BM36", "BM36"), ct_correction =  TRUE, fdr_val = 0.01)
count_interactions(c("BM44", "BM44"), ct_correction =  TRUE, fdr_val = 0.01)


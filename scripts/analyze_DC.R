# Author: B P Kailash
# Purpose: Analyze the DC relations observed

# What does Part 1 of the code do?
# 1) Find the sample size of the analysis
# 2) Convert all DC edges from ensids to HGNC symbols
# 3) Add M, NM, cell type labels to each DC edge
# 4) Add braak scores to the sample ids
# 5) Find the dlpfc snRNAseq DEG within a cell type 
# 6) Find the intersection of dlpfc in the BM10 combinations

rm(list = ls())

# Global variables --------------------------------------------------------

source("scripts/libraries.R")
source("scripts/utils.R")

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

# Convert significant correlation edges to HGNC mapping ------------------

ensembl_to_approved_gene_symbol_corr <- function(MappingFile, PathResults, BR, ct_corrected, fdr_val){
  
  if(ct_corrected == TRUE){
    ct_file_name <- "_ct_correction"
    UnionCorrelatedEdgeFile <- read.csv(paste("results/", "UnionCorrelatedEdges", ct_file_name, BR[1],BR[2], "_fdr_", fdr_val, ".csv", sep = ""),
                                        stringsAsFactors = FALSE)
  } else {
    ct_file_name = ""
    UnionCorrelatedEdgeFile <- read.csv(paste("results/", "UnionCorrelatedEdges", ct_file_name, BR[1],BR[2], "_fdr", fdr_val, ".csv", sep = ""),
                                        stringsAsFactors = FALSE)
  }
  
  
  ens1 <- UnionCorrelatedEdgeFile$a
  ens2 <- UnionCorrelatedEdgeFile$b
  map1 <- match(ens1, MappingFile$gene_id)
  map2 <- match(ens2, MappingFile$gene_id)
  g1 <- MappingFile$gene_name[map1]
  g2 <- MappingFile$gene_name[map2]
  
  geneMapped <- cbind(UnionCorrelatedEdgeFile, g1, g2)
  count_before_mapping <- dim(geneMapped)[[1]]
  
  geneMapped <- geneMapped[-union(which(is.na(geneMapped$g1)), which(is.na(geneMapped$g2))),]
  count_after_mapping <- dim(geneMapped)[[1]]
  
  cat(paste("Number of significant correlated edges before mapping for: ", BR[1], "-", BR[2], count_before_mapping), sep = "")
  cat(paste("\nNumber of DC edges after mapping for: ", BR[1], "-", BR[2], count_after_mapping), sep = "")
  cat("\nNumber of edges lost to mapping are: ", count_before_mapping - count_after_mapping, "\n")
  write.csv(geneMapped , file = paste(PathResults, BR[1],BR[2], ct_file_name, "Mapped_CorrelatedEdgesWithoutNA_fdr_", fdr_val, ".csv",sep=""),
            row.names = FALSE)
}

# ct_corrected = FALSE

## Inter-region
ensembl_to_approved_gene_symbol_corr(MappingFile = mapping_file, 
                                     PathResults = "results/",
                                     BR = c("BM_10", "BM_22"), ct_corrected = FALSE, fdr_val = 0.01)

ensembl_to_approved_gene_symbol_corr(MappingFile = mapping_file, 
                                     PathResults = "results/", 
                                     BR = c("BM_10", "BM_36"), ct_corrected = FALSE, fdr_val = 0.01)

ensembl_to_approved_gene_symbol_corr(MappingFile = mapping_file, 
                                     PathResults = "results/", 
                                     BR = c("BM_10", "BM_44"), ct_corrected = FALSE, fdr_val = 0.01)

ensembl_to_approved_gene_symbol_corr(MappingFile = mapping_file, 
                                     PathResults = "results/", 
                                     BR = c("BM_22", "BM_36"), ct_corrected = FALSE, fdr_val = 0.01)

ensembl_to_approved_gene_symbol_corr(MappingFile = mapping_file, 
                                     PathResults = "results/", 
                                     BR = c("BM_22", "BM_44"), ct_corrected = FALSE, fdr_val = 0.01)

ensembl_to_approved_gene_symbol_corr(MappingFile = mapping_file, 
                                     PathResults = "results/", 
                                     BR = c("BM_36", "BM_44"), ct_corrected = FALSE, fdr_val = 0.01)

## Within-region

ensembl_to_approved_gene_symbol_corr(MappingFile = mapping_file, 
                                     PathResults = "results/", 
                                     BR = c("BM_10", "BM_10"), ct_corrected = FALSE, fdr_val = 0.01)

ensembl_to_approved_gene_symbol_corr(MappingFile = mapping_file, 
                                     PathResults = "results/", 
                                     BR = c("BM_22", "BM_22"), ct_corrected = FALSE, fdr_val = 0.01)

ensembl_to_approved_gene_symbol_corr(MappingFile = mapping_file, 
                                     PathResults = "results/", 
                                     BR = c("BM_36", "BM_36"), ct_corrected = FALSE, fdr_val = 0.01)

ensembl_to_approved_gene_symbol_corr(MappingFile = mapping_file, 
                                     PathResults = "results/", 
                                     BR = c("BM_44", "BM_44"), ct_corrected = FALSE, fdr_val = 0.01)

# ct_corrected = TRUE

## Inter-region
ensembl_to_approved_gene_symbol_corr(MappingFile = mapping_file, 
                                     PathResults = "results/",
                                     BR = c("BM10", "BM22"), ct_corrected = TRUE, fdr_val = 0.01)

ensembl_to_approved_gene_symbol_corr(MappingFile = mapping_file, 
                                     PathResults = "results/", 
                                     BR = c("BM10", "BM36"), ct_corrected = TRUE, fdr_val = 0.01)

ensembl_to_approved_gene_symbol_corr(MappingFile = mapping_file, 
                                     PathResults = "results/", 
                                     BR = c("BM10", "BM44"), ct_corrected = TRUE, fdr_val = 0.01)

ensembl_to_approved_gene_symbol_corr(MappingFile = mapping_file, 
                                     PathResults = "results/", 
                                     BR = c("BM22", "BM36"), ct_corrected = TRUE, fdr_val = 0.01)

ensembl_to_approved_gene_symbol_corr(MappingFile = mapping_file, 
                                     PathResults = "results/", 
                                     BR = c("BM22", "BM44"), ct_corrected = TRUE, fdr_val = 0.01)

ensembl_to_approved_gene_symbol_corr(MappingFile = mapping_file, 
                                     PathResults = "results/", 
                                     BR = c("BM36", "BM44"), ct_corrected = TRUE, fdr_val = 0.01)

## Within-region

ensembl_to_approved_gene_symbol_corr(MappingFile = mapping_file, 
                                     PathResults = "results/", 
                                     BR = c("BM10", "BM10"), ct_corrected = TRUE, fdr_val = 0.01)

ensembl_to_approved_gene_symbol_corr(MappingFile = mapping_file, 
                                     PathResults = "results/", 
                                     BR = c("BM22", "BM22"), ct_corrected = TRUE, fdr_val = 0.01)

ensembl_to_approved_gene_symbol_corr(MappingFile = mapping_file, 
                                     PathResults = "results/", 
                                     BR = c("BM36", "BM36"), ct_corrected = TRUE, fdr_val = 0.01)

ensembl_to_approved_gene_symbol_corr(MappingFile = mapping_file, 
                                     PathResults = "results/", 
                                     BR = c("BM44", "BM_44"), ct_corrected = TRUE, fdr_val = 0.01)

# Convert DC edges to HGNC mapping ----------------------------------------

ensembl_to_approved_gene_symbol_dc <- function(DCEdgeFile, MappingFile, PathResults, BR, ct_corrected, fdr_val){
  
  if(ct_corrected == TRUE){
    ct_file_name <- "_ct_correction"
  } else {
    ct_file_name = ""
  }
  
  DCEdgeFile <- read.csv(paste("results/", "DCEdges", ct_file_name, BR[1], BR[2], "fdr_", fdr_val, ".csv", sep = ""),
                         stringsAsFactors = FALSE)
  ens1 <- DCEdgeFile$a
  ens2 <- DCEdgeFile$b
  map1 <- match(ens1, MappingFile$gene_id)
  map2 <- match(ens2, MappingFile$gene_id)
  g1 <- MappingFile$gene_name[map1]
  g2 <- MappingFile$gene_name[map2]
  geneMapped <- cbind(DCEdgeFile, g1, g2)
  colnames(geneMapped) <- c(colnames(DCEdgeFile), "g1", "g2")
  dc_count_before_mapping <- dim(geneMapped)[[1]]
  cat(paste("Number of DC edges before mapping for ", BR[1], "-", BR[2], dc_count_before_mapping), sep = "")
  geneMapped <- geneMapped[-union(which(is.na(geneMapped$g1)),which(is.na(geneMapped$g2))),]
  dc_count_after_mapping <- dim(geneMapped)[[1]]
  cat(paste("\nNumber of DC edges after mapping for ", BR[1], "-", BR[2], dc_count_after_mapping), sep = "")
  # Losing a lot of edges because of mapping not being available.
  cat("\nNumber of edges lost to mapping are", dc_count_before_mapping - dc_count_after_mapping, "\n")
  write.csv(geneMapped , file = paste(PathResults,BR[1],BR[2], ct_file_name, "Mapped_DCEdgesWithoutNA_fdr_", fdr_val, ".csv",sep=""),
            row.names = FALSE)
  
}

# ct_corrected = FALSE

## Inter-region
ensembl_to_approved_gene_symbol_dc(
  MappingFile = mapping_file, 
  PathResults = "results/", 
  BR = c("BM_10", "BM_22"), ct_corrected = FALSE, fdr_val = 0.01)

ensembl_to_approved_gene_symbol_dc(
  MappingFile = mapping_file, 
  PathResults = "results/", 
  BR = c("BM_10", "BM_36"), ct_corrected = FALSE, fdr_val = 0.01)

ensembl_to_approved_gene_symbol_dc(
  MappingFile = mapping_file, 
  PathResults = "results/", 
  BR = c("BM_10", "BM_44"), ct_corrected = FALSE, fdr_val = 0.01)

ensembl_to_approved_gene_symbol_dc(
  MappingFile = mapping_file, 
  PathResults = "results/", 
  BR = c("BM_22", "BM_36"), ct_corrected = FALSE, fdr_val = 0.01)

ensembl_to_approved_gene_symbol_dc(
  MappingFile = mapping_file, 
  PathResults = "results/", 
  BR = c("BM_22", "BM_44"), ct_corrected = FALSE, fdr_val = 0.01)

ensembl_to_approved_gene_symbol_dc(
  MappingFile = mapping_file, 
  PathResults = "results/", 
  BR = c("BM_36", "BM_44"), ct_corrected = FALSE, fdr_val = 0.01)

## Within-region

ensembl_to_approved_gene_symbol_dc(
  MappingFile = mapping_file, 
  PathResults = "results/", 
  BR = c("BM_10", "BM_10"), ct_corrected = FALSE, fdr_val = 0.01)

ensembl_to_approved_gene_symbol_dc(
  MappingFile = mapping_file, 
  PathResults = "results/", 
  BR = c("BM_22", "BM_22"), ct_corrected = FALSE, fdr_val = 0.01)

ensembl_to_approved_gene_symbol_dc(
  MappingFile = mapping_file, 
  PathResults = "results/", 
  BR = c("BM_36", "BM_36"), ct_corrected = FALSE, fdr_val = 0.01)

ensembl_to_approved_gene_symbol_dc(
  MappingFile = mapping_file, 
  PathResults = "results/", 
  BR = c("BM_44", "BM_44"), ct_corrected = FALSE, fdr_val = 0.01)

# ct_corrected = FALSE

## Inter-region
ensembl_to_approved_gene_symbol_dc(
  MappingFile = mapping_file, 
  PathResults = "results/", 
  BR = c("BM10", "BM22"), ct_corrected = TRUE, fdr_val = 0.01)

ensembl_to_approved_gene_symbol_dc(
  MappingFile = mapping_file, 
  PathResults = "results/", 
  BR = c("BM10", "BM36"), ct_corrected = TRUE, fdr_val = 0.01)

ensembl_to_approved_gene_symbol_dc(
  MappingFile = mapping_file, 
  PathResults = "results/", 
  BR = c("BM10", "BM44"), ct_corrected = TRUE, fdr_val = 0.01)

ensembl_to_approved_gene_symbol_dc(
  MappingFile = mapping_file, 
  PathResults = "results/", 
  BR = c("BM22", "BM36"), ct_corrected = TRUE, fdr_val = 0.01)

ensembl_to_approved_gene_symbol_dc(
  MappingFile = mapping_file, 
  PathResults = "results/", 
  BR = c("BM22", "BM44"), ct_corrected = TRUE, fdr_val = 0.01)

ensembl_to_approved_gene_symbol_dc(
  MappingFile = mapping_file, 
  PathResults = "results/", 
  BR = c("BM36", "BM44"), ct_corrected = TRUE, fdr_val = 0.01)

## Within-region

ensembl_to_approved_gene_symbol_dc(
  MappingFile = mapping_file, 
  PathResults = "results/", 
  BR = c("BM10", "BM10"), ct_corrected = TRUE, fdr_val = 0.01)

ensembl_to_approved_gene_symbol_dc(
  MappingFile = mapping_file, 
  PathResults = "results/", 
  BR = c("BM22", "BM22"), ct_corrected = TRUE, fdr_val = 0.01)

ensembl_to_approved_gene_symbol_dc(
  MappingFile = mapping_file, 
  PathResults = "results/", 
  BR = c("BM36", "BM36"), ct_corrected = TRUE, fdr_val = 0.01)

ensembl_to_approved_gene_symbol_dc(
  MappingFile = mapping_file, 
  PathResults = "results/", 
  BR = c("BM44", "BM44"), ct_corrected = TRUE, fdr_val = 0.01)

# Labeling significantly correlated edges ---------------------------------

label_correlated_edges <- function(CellTypes, NumberOfMarkersPerCellType, BR, CTCorrection, fdr_val){
  
  if(CTCorrection == TRUE){
    ct_file_name <- "_ct_correction"
  } else {
    ct_file_name = ""
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
  
  correlated_edges <- read.csv(paste("results/", BR[1], BR[2], ct_file_name, "Mapped_CorrelatedEdgesWithoutNA_fdr_", fdr_val, ".csv", sep = ""),
                               stringsAsFactors = FALSE)
  
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
  
  write.csv(correlated_edges, paste("results/", BR[1],BR[2], ct_file_name, "Mapped_CorrelatedEdgesWithoutNA_fdr_", fdr_val, ".csv",sep=""))
}

# ct_corrected = FALSE

## Inter-region

label_correlated_edges(BR = c("BM_10",  "BM_22"), CellTypes = cell_types, NumberOfMarkersPerCellType = number_mg_per_ct, CTCorrection = FALSE, fdr_val = 0.01)
label_correlated_edges(BR = c("BM_10",  "BM_36"), CellTypes = cell_types, NumberOfMarkersPerCellType = number_mg_per_ct, CTCorrection = FALSE, fdr_val = 0.01)
label_correlated_edges(BR = c("BM_10",  "BM_44"), CellTypes = cell_types, NumberOfMarkersPerCellType = number_mg_per_ct, CTCorrection = FALSE, fdr_val = 0.01)
label_correlated_edges(BR = c("BM_22",  "BM_36"), CellTypes = cell_types, NumberOfMarkersPerCellType = number_mg_per_ct, CTCorrection = FALSE, fdr_val = 0.01)
label_correlated_edges(BR = c("BM_22",  "BM_44"), CellTypes = cell_types, NumberOfMarkersPerCellType = number_mg_per_ct, CTCorrection = FALSE, fdr_val = 0.01)
label_correlated_edges(BR = c("BM_36",  "BM_44"), CellTypes = cell_types, NumberOfMarkersPerCellType = number_mg_per_ct, CTCorrection = FALSE, fdr_val = 0.01)

## Within-region

label_correlated_edges(BR = c("BM_10",  "BM_10"), CellTypes = cell_types, NumberOfMarkersPerCellType = number_mg_per_ct, CTCorrection = FALSE, fdr_val = 0.01)
label_correlated_edges(BR = c("BM_22",  "BM_22"), CellTypes = cell_types, NumberOfMarkersPerCellType = number_mg_per_ct, CTCorrection = FALSE, fdr_val = 0.01)
label_correlated_edges(BR = c("BM_36",  "BM_36"), CellTypes = cell_types, NumberOfMarkersPerCellType = number_mg_per_ct, CTCorrection = FALSE, fdr_val = 0.01)
label_correlated_edges(BR = c("BM_44",  "BM_44"), CellTypes = cell_types, NumberOfMarkersPerCellType = number_mg_per_ct, CTCorrection = FALSE, fdr_val = 0.01)

# ct_corrected = TRUE

## Inter-region

label_correlated_edges(BR = c("BM10",  "BM22"), CellTypes = cell_types, NumberOfMarkersPerCellType = number_mg_per_ct, CTCorrection = TRUE, fdr_val = 0.01)
label_correlated_edges(BR = c("BM10",  "BM36"), CellTypes = cell_types, NumberOfMarkersPerCellType = number_mg_per_ct, CTCorrection = TRUE, fdr_val = 0.01)
label_correlated_edges(BR = c("BM10",  "BM44"), CellTypes = cell_types, NumberOfMarkersPerCellType = number_mg_per_ct, CTCorrection = TRUE, fdr_val = 0.01)
label_correlated_edges(BR = c("BM22",  "BM36"), CellTypes = cell_types, NumberOfMarkersPerCellType = number_mg_per_ct, CTCorrection = TRUE, fdr_val = 0.01)
label_correlated_edges(BR = c("BM22",  "BM44"), CellTypes = cell_types, NumberOfMarkersPerCellType = number_mg_per_ct, CTCorrection = TRUE, fdr_val = 0.01)
label_correlated_edges(BR = c("BM36",  "BM44"), CellTypes = cell_types, NumberOfMarkersPerCellType = number_mg_per_ct, CTCorrection = TRUE, fdr_val = 0.01)

## Within-region

label_correlated_edges(BR = c("BM10",  "BM10"), CellTypes = cell_types, NumberOfMarkersPerCellType = number_mg_per_ct, CTCorrection = TRUE, fdr_val = 0.01)
label_correlated_edges(BR = c("BM22",  "BM22"), CellTypes = cell_types, NumberOfMarkersPerCellType = number_mg_per_ct, CTCorrection = TRUE, fdr_val = 0.01)
label_correlated_edges(BR = c("BM36",  "BM36"), CellTypes = cell_types, NumberOfMarkersPerCellType = number_mg_per_ct, CTCorrection = TRUE, fdr_val = 0.01)
label_correlated_edges(BR = c("BM44",  "BM44"), CellTypes = cell_types, NumberOfMarkersPerCellType = number_mg_per_ct, CTCorrection = TRUE, fdr_val = 0.01)

# Label significant DC edges --------------------------------

label_DC_edges <- function(CellTypes, NumberOfMarkersPerCellType, BR, CTCorrection, fdr_val){
  
  if(CTCorrection == TRUE){
    ct_file_name <- "_ct_correction"
  } else {
    ct_file_name = ""
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
  
  DC_edges <- read.csv(paste("results/", BR[1], BR[2], ct_file_name, "Mapped_DCEdgesWithoutNA_fdr_", fdr_val, ".csv", sep = ""),
                       stringsAsFactors = FALSE)
  
  marker_label_g1 <- rep(0, dim(DC_edges)[1])
  marker_label_g2 <- rep(0, dim(DC_edges)[1])
  
  match1 <- match(toupper(DC_edges$g1), mg)
  match2 <- match(toupper(DC_edges$g2), mg)
  
  marker_label_g1[which(is.na(match1))] = "NM"
  marker_label_g1[which(marker_label_g1 == 0)] = "M"
  marker_label_g2[which(is.na(match2))] = "NM"
  marker_label_g2[which(marker_label_g2 == 0)] = "M"
  
  DC_edges$m1 <- marker_label_g1
  DC_edges$m2 <- marker_label_g2
  
  DC_edges$c1 <- ct[match1]
  DC_edges$c2 <- ct[match2]
  
  write.csv(DC_edges, paste("results/", BR[1], BR[2], ct_file_name, "Mapped_DCEdgesWithoutNA_fdr_", fdr_val, ".csv", sep = ""))
}

cell_types <- c("ast", "mic", "neu", "oli")
number_mg <- 1000
number_mg_per_ct <- rep(number_mg, length(cell_types))

# ct_corrected = FALSE

## Inter-region

label_DC_edges(BR = c("BM_10",  "BM_22"), CellTypes = cell_types, NumberOfMarkersPerCellType = number_mg_per_ct, CTCorrection = FALSE, fdr_val = 0.01)
label_DC_edges(BR = c("BM_10",  "BM_36"), CellTypes = cell_types, NumberOfMarkersPerCellType = number_mg_per_ct, CTCorrection = FALSE, fdr_val = 0.01)
label_DC_edges(BR = c("BM_10",  "BM_44"), CellTypes = cell_types, NumberOfMarkersPerCellType = number_mg_per_ct, CTCorrection = FALSE, fdr_val = 0.01)
label_DC_edges(BR = c("BM_22",  "BM_36"), CellTypes = cell_types, NumberOfMarkersPerCellType = number_mg_per_ct, CTCorrection = FALSE, fdr_val = 0.01)
label_DC_edges(BR = c("BM_22",  "BM_44"), CellTypes = cell_types, NumberOfMarkersPerCellType = number_mg_per_ct, CTCorrection = FALSE, fdr_val = 0.01)
label_DC_edges(BR = c("BM_36",  "BM_44"), CellTypes = cell_types, NumberOfMarkersPerCellType = number_mg_per_ct, CTCorrection = FALSE, fdr_val = 0.01)

## Within-region

label_DC_edges(BR = c("BM_10",  "BM_10"), CellTypes = cell_types, NumberOfMarkersPerCellType = number_mg_per_ct, CTCorrection = FALSE, fdr_val = 0.01)
label_DC_edges(BR = c("BM_22",  "BM_22"), CellTypes = cell_types, NumberOfMarkersPerCellType = number_mg_per_ct, CTCorrection = FALSE, fdr_val = 0.01)
label_DC_edges(BR = c("BM_36",  "BM_36"), CellTypes = cell_types, NumberOfMarkersPerCellType = number_mg_per_ct, CTCorrection = FALSE, fdr_val = 0.01)
label_DC_edges(BR = c("BM_44",  "BM_44"), CellTypes = cell_types, NumberOfMarkersPerCellType = number_mg_per_ct, CTCorrection = FALSE, fdr_val = 0.01)

# ct_corrected = TRUE

## Inter-region

label_DC_edges(BR = c("BM10",  "BM22"), CellTypes = cell_types, NumberOfMarkersPerCellType = number_mg_per_ct, CTCorrection = TRUE, fdr_val = 0.01)
label_DC_edges(BR = c("BM10",  "BM36"), CellTypes = cell_types, NumberOfMarkersPerCellType = number_mg_per_ct, CTCorrection = TRUE, fdr_val = 0.01)
label_DC_edges(BR = c("BM10",  "BM44"), CellTypes = cell_types, NumberOfMarkersPerCellType = number_mg_per_ct, CTCorrection = TRUE, fdr_val = 0.01)
label_DC_edges(BR = c("BM22",  "BM36"), CellTypes = cell_types, NumberOfMarkersPerCellType = number_mg_per_ct, CTCorrection = TRUE, fdr_val = 0.01)
label_DC_edges(BR = c("BM22",  "BM44"), CellTypes = cell_types, NumberOfMarkersPerCellType = number_mg_per_ct, CTCorrection = TRUE, fdr_val = 0.01)
label_DC_edges(BR = c("BM36",  "BM44"), CellTypes = cell_types, NumberOfMarkersPerCellType = number_mg_per_ct, CTCorrection = TRUE, fdr_val = 0.01)

## Within-region

label_DC_edges(BR = c("BM10",  "BM10"), CellTypes = cell_types, NumberOfMarkersPerCellType = number_mg_per_ct, CTCorrection = TRUE, fdr_val = 0.01)
label_DC_edges(BR = c("BM22",  "BM22"), CellTypes = cell_types, NumberOfMarkersPerCellType = number_mg_per_ct, CTCorrection = TRUE, fdr_val = 0.01)
label_DC_edges(BR = c("BM36",  "BM36"), CellTypes = cell_types, NumberOfMarkersPerCellType = number_mg_per_ct, CTCorrection = TRUE, fdr_val = 0.01)
label_DC_edges(BR = c("BM44",  "BM44"), CellTypes = cell_types, NumberOfMarkersPerCellType = number_mg_per_ct, CTCorrection = TRUE, fdr_val = 0.01)

# Adding Chr and Gene_type labels to correlated edges ---------------------

add_chr_geneType_corr <- function(BR, CTCorrection, fdr_val){
  
  if(CTCorrection == TRUE){
    ct_file_name <- "_ct_correction"
  } else {
    ct_file_name = ""
  }
  
  correlated_edges <- read.csv(paste("results/", BR[1], BR[2], ct_file_name, "Mapped_CorrelatedEdgesWithoutNA_fdr_", fdr_val, ".csv", sep = ""),
                               stringsAsFactors = FALSE)
  correlated_edges$gene_type_1 <- mapping_file$gene_type[match(correlated_edges$g1, mapping_file$gene_name)]
  correlated_edges$gene_type_2 <- mapping_file$gene_type[match(correlated_edges$g2, mapping_file$gene_name)]
  correlated_edges$chr1 <- mapping_file$chr[match(correlated_edges$g1, mapping_file$gene_name)]
  correlated_edges$chr2 <- mapping_file$chr[match(correlated_edges$g2, mapping_file$gene_name)]
  write.csv(correlated_edges, paste("results/", BR[1],BR[2], ct_file_name, "Mapped_CorrelatedEdgesWithoutNA_fdr_", fdr_val, ".csv",sep=""))
}

# ct_corrected = FALSE

## Inter-region

add_chr_geneType_corr(BR = c("BM_10",  "BM_22"), CTCorrection = FALSE, fdr_val = 0.01)
add_chr_geneType_corr(BR = c("BM_10",  "BM_36"), CTCorrection = FALSE, fdr_val = 0.01)
add_chr_geneType_corr(BR = c("BM_10",  "BM_44"), CTCorrection = FALSE, fdr_val = 0.01)
add_chr_geneType_corr(BR = c("BM_22",  "BM_36"), CTCorrection = FALSE, fdr_val = 0.01)
add_chr_geneType_corr(BR = c("BM_22",  "BM_44"), CTCorrection = FALSE, fdr_val = 0.01)
add_chr_geneType_corr(BR = c("BM_36",  "BM_44"), CTCorrection = FALSE, fdr_val = 0.01)

## Within-region

add_chr_geneType_corr(BR = c("BM_10",  "BM_10"), CTCorrection = FALSE, fdr_val = 0.01)
add_chr_geneType_corr(BR = c("BM_22",  "BM_22"), CTCorrection = FALSE, fdr_val = 0.01)
add_chr_geneType_corr(BR = c("BM_36",  "BM_36"), CTCorrection = FALSE, fdr_val = 0.01)
add_chr_geneType_corr(BR = c("BM_44",  "BM_44"), CTCorrection = FALSE, fdr_val = 0.01)

# ct_corrected = TRUE

## Inter-region

add_chr_geneType_corr(BR = c("BM10",  "BM22"), CTCorrection = TRUE, fdr_val = 0.01)
add_chr_geneType_corr(BR = c("BM10",  "BM36"), CTCorrection = TRUE, fdr_val = 0.01)
add_chr_geneType_corr(BR = c("BM10",  "BM44"), CTCorrection = TRUE, fdr_val = 0.01)
add_chr_geneType_corr(BR = c("BM22",  "BM36"), CTCorrection = TRUE, fdr_val = 0.01)
add_chr_geneType_corr(BR = c("BM22",  "BM44"), CTCorrection = TRUE, fdr_val = 0.01)
add_chr_geneType_corr(BR = c("BM36",  "BM44"), CTCorrection = TRUE, fdr_val = 0.01)

## Within-region

add_chr_geneType_corr(BR = c("BM10",  "BM10"), CTCorrection = TRUE, fdr_val = 0.01)
add_chr_geneType_corr(BR = c("BM22",  "BM22"), CTCorrection = TRUE, fdr_val = 0.01)
add_chr_geneType_corr(BR = c("BM36",  "BM36"), CTCorrection = TRUE, fdr_val = 0.01)
add_chr_geneType_corr(BR = c("BM44",  "BM44"), CTCorrection = TRUE, fdr_val = 0.01)

# Adding chr, gene type to dc edges ---------------------------------------

add_chr_geneType_dc <- function(BR, CTCorrection, fdr_val){
  
  if(CTCorrection == TRUE){
    ct_file_name <- "_ct_correction"
  } else {
    ct_file_name = ""
  }
  
  DC_edges <- read.csv(paste("results/", BR[1], BR[2], ct_file_name, "Mapped_DCEdgesWithoutNA_fdr_", fdr_val, ".csv", sep = ""),
                       stringsAsFactors = FALSE)
  DC_edges$gene_type_1 <- mapping_file$gene_type[match(DC_edges$g1, mapping_file$gene_name)]
  DC_edges$gene_type_2 <- mapping_file$gene_type[match(DC_edges$g2, mapping_file$gene_name)]
  DC_edges$chr1 <- mapping_file$chr[match(DC_edges$g1, mapping_file$gene_name)]
  DC_edges$chr2 <- mapping_file$chr[match(DC_edges$g2, mapping_file$gene_name)]
  write.csv(DC_edges, paste("results/", BR[1], BR[2], ct_file_name, "Mapped_DCEdgesWithoutNA_fdr_", fdr_val, ".csv", sep = ""))
}

# ct_corrected = FALSE

## Inter-region

add_chr_geneType_dc(BR = c("BM_10",  "BM_22"), CTCorrection = FALSE, fdr_val = 0.01)
add_chr_geneType_dc(BR = c("BM_10",  "BM_36"), CTCorrection = FALSE, fdr_val = 0.01)
add_chr_geneType_dc(BR = c("BM_10",  "BM_44"), CTCorrection = FALSE, fdr_val = 0.01)
add_chr_geneType_dc(BR = c("BM_22",  "BM_36"), CTCorrection = FALSE, fdr_val = 0.01)
add_chr_geneType_dc(BR = c("BM_22",  "BM_44"), CTCorrection = FALSE, fdr_val = 0.01)
add_chr_geneType_dc(BR = c("BM_36",  "BM_44"), CTCorrection = FALSE, fdr_val = 0.01)

## Within-region

add_chr_geneType_dc(BR = c("BM_10",  "BM_10"), CTCorrection = FALSE, fdr_val = 0.01)
add_chr_geneType_dc(BR = c("BM_22",  "BM_22"), CTCorrection = FALSE, fdr_val = 0.01)
add_chr_geneType_dc(BR = c("BM_36",  "BM_36"), CTCorrection = FALSE, fdr_val = 0.01)
add_chr_geneType_dc(BR = c("BM_44",  "BM_44"), CTCorrection = FALSE, fdr_val = 0.01)

# ct_corrected = TRUE

## Inter-region

add_chr_geneType_dc(BR = c("BM10",  "BM22"), CTCorrection = TRUE, fdr_val = 0.01)
add_chr_geneType_dc(BR = c("BM10",  "BM36"), CTCorrection = TRUE, fdr_val = 0.01)
add_chr_geneType_dc(BR = c("BM10",  "BM44"), CTCorrection = TRUE, fdr_val = 0.01)
add_chr_geneType_dc(BR = c("BM22",  "BM36"), CTCorrection = TRUE, fdr_val = 0.01)
add_chr_geneType_dc(BR = c("BM22",  "BM44"), CTCorrection = TRUE, fdr_val = 0.01)
add_chr_geneType_dc(BR = c("BM36",  "BM44"), CTCorrection = TRUE, fdr_val = 0.01)

## Within-region

add_chr_geneType_dc(BR = c("BM10",  "BM10"), CTCorrection = TRUE, fdr_val = 0.01)
add_chr_geneType_dc(BR = c("BM22",  "BM22"), CTCorrection = TRUE, fdr_val = 0.01)
add_chr_geneType_dc(BR = c("BM36",  "BM36"), CTCorrection = TRUE, fdr_val = 0.01)
add_chr_geneType_dc(BR = c("BM44",  "BM44"), CTCorrection = TRUE, fdr_val = 0.01)

# Count M-M, NM-NM, M-NM, NM-M interactions -------------------------------

count_interactions <- function(BR, ct_correction, fdr_val){
  if(ct_correction == TRUE){
    ct_file_name <- "_ct_correction"
  } else {
    ct_file_name = ""
  }
  dc <- read.csv(paste("results/", BR[1], BR[2], ct_file_name,
                       "Mapped_DCEdgesWithoutNA_fdr_", fdr_val, ".csv", sep = ""), 
                 stringsAsFactors = FALSE)
  interaction = paste(dc$m1, dc$m2, sep = "")
  cat(paste("Number of M-M interactions for"), BR[1], BR[2], "is", length(which(interaction == "MM")))
  cat(paste("\nNumber of M-NM interactions is"), BR[1], BR[2], "is", length(which(interaction == "MNM")))
  cat(paste("\nNumber of NM-M interactions is"), BR[1], BR[2], "is",length(which(interaction == "NMM")))
  cat(paste("\nNumber of NM-NM interactions is"), BR[1], BR[2], "is",length(which(interaction == "NMNM\n")))
}

# NOTE that there are no duplicates among the 4000 markers for the 4 cell types
length(unique(markers_df_human_brain$markers[which(markers_df_human_brain$cell %in% c("ast", "mic", "neu", "oli"))]))

# ct_corrected = FALSE

## Inter-region

count_interactions(c("BM_10", "BM_22"), ct_correction =  FALSE, fdr_val = 0.01)
count_interactions(c("BM_10", "BM_36"), ct_correction =  FALSE, fdr_val = 0.01)
count_interactions(c("BM_10", "BM_44"), ct_correction =  FALSE, fdr_val = 0.01)
count_interactions(c("BM_22", "BM_36"), ct_correction =  FALSE, fdr_val = 0.01)
count_interactions(c("BM_22", "BM_44"), ct_correction =  FALSE, fdr_val = 0.01)
count_interactions(c("BM_36", "BM_44"), ct_correction =  FALSE, fdr_val = 0.01)

## Within-region

count_interactions(c("BM_10", "BM_10"), ct_correction =  FALSE, fdr_val = 0.01)
count_interactions(c("BM_22", "BM_22"), ct_correction =  FALSE, fdr_val = 0.01)
count_interactions(c("BM_36", "BM_36"), ct_correction =  FALSE, fdr_val = 0.01)
count_interactions(c("BM_44", "BM_44"), ct_correction =  FALSE, fdr_val = 0.01)
                   
# ct_corrected = TRUE

## Inter-region

count_interactions(c("BM10", "BM22"), ct_correction =  TRUE, fdr_val = 0.01)
count_interactions(c("BM10", "BM36"), ct_correction =  TRUE, fdr_val = 0.01)
count_interactions(c("BM10", "BM44"), ct_correction =  TRUE, fdr_val = 0.01)
count_interactions(c("BM22", "BM36"), ct_correction =  TRUE, fdr_val = 0.01)
count_interactions(c("BM22", "BM44"), ct_correction =  TRUE, fdr_val = 0.01)
count_interactions(c("BM36", "BM44"), ct_correction =  TRUE, fdr_val = 0.01)

## Within-region

count_interactions(c("BM10", "BM10"), ct_correction =  TRUE, fdr_val = 0.01)
count_interactions(c("BM22", "BM22"), ct_correction =  TRUE, fdr_val = 0.01)
count_interactions(c("BM36", "BM36"), ct_correction =  TRUE, fdr_val = 0.01)
count_interactions(c("BM44", "BM44"), ct_correction =  TRUE, fdr_val = 0.01)

# Adding braak score to the sample ids ------------------------------------

metadata <- read.csv(paste("data/", "rna_metadata_onlyBAM_noExclude_noDup_Nov2018.csv", sep = ""), stringsAsFactors = FALSE)

label_with_braak <- function(BR){
  ad_ids <- read.csv(paste("results/", "AD_ids", BR[1], BR[2], ".csv", sep = ""), stringsAsFactors = FALSE)
  ctl_ids <- read.csv(paste("results/", "CTL_ids", BR[1], BR[2], ".csv", sep = ""), stringsAsFactors = FALSE)
  
  ad_ids$braak_score <- metadata$bbscore[match(ad_ids$x, metadata$individualIdentifier)]
  ctl_ids$braak_score <- metadata$bbscore[match(ctl_ids$x, metadata$individualIdentifier)]
  
  write.csv(ad_ids, paste("results/", "AD_ids", BR[1], BR[2], ".csv", sep = ""))
  write.csv(ctl_ids, paste("results/", "CTL_ids", BR[1], BR[2], ".csv", sep = ""))
}

label_with_braak(c("BM10", "BM22"))
label_with_braak(c("BM10", "BM36"))
label_with_braak(c("BM10", "BM44"))
label_with_braak(c("BM22", "BM36"))
label_with_braak(c("BM22", "BM44"))
label_with_braak(c("BM36", "BM44"))

# Find DLPFC scRNAseq DEG -------------------------------------------------

# Purpose: Find cell intrinsic DEG (ciDEG) from DE results
# Data used is Supplementary table 2 from the following paper
# https://www.nature.com/articles/s41586-019-1195-2
# Single-cell transcriptomic analysis of Alzheimer's disease

# Loading necessary libraries

library(readxl)

# Loading data

# NOTE: The following data is from the following paper
# https://www.biorxiv.org/content/10.1101/2020.04.15.041376v1.full
# Deciphering cellular transcriptional alterations in Alzheimer's disease brains

# Finding the independent mixed model cell-intrinsic DEGs
find_deg <- function(df, cellType){
  df <- as.data.frame(df[,1:9])
  colnames(df) <- df[1,]
  df <- df[-1,]
  row.names(df) <- df[,1]
  df <- df[,-1]
  df$cellType <- rep(cellType, length(row.names(df)))
  return(df[which(df$DEGs.Ind.Mix.models == TRUE),])
} 

# DLPFC snRNAseq DEG in DC edges ------------------------------------------

ex_df <- read_excel(paste("data/", "41586_2019_1195_MOESM4_ESM.xlsx", sep = ""), sheet = 2)
in_df <- read_excel(paste("data/", "41586_2019_1195_MOESM4_ESM.xlsx", sep = ""), sheet = 3)
ast_df <- read_excel(paste("data/", "41586_2019_1195_MOESM4_ESM.xlsx", sep = ""), sheet = 4)
oli_df <- read_excel(paste("data/", "41586_2019_1195_MOESM4_ESM.xlsx", sep = ""), sheet = 5)
mic_df <- read_excel(paste("data/", "41586_2019_1195_MOESM4_ESM.xlsx", sep = ""), sheet = 7)

ast_deg <- row.names(find_deg(ast_df, "Ast"))
mic_deg <- row.names(find_deg(mic_df, "Mic"))
ex_deg <- row.names(find_deg(ex_df, "Ex"))
in_deg <- row.names(find_deg(in_df, "In"))
oli_deg <- row.names(find_deg(oli_df, "Oli"))

stopifnot(length(unique(ast_deg)) == length(ast_deg))
stopifnot(length(unique(mic_deg)) == length(mic_deg))
stopifnot(length(unique(ex_deg)) == length(ex_deg))
stopifnot(length(unique(in_deg)) == length(in_deg))
stopifnot(length(unique(oli_deg)) == length(oli_deg))
dlpfc_scRNAseq_deg <- unique(c(ex_deg, in_deg, ast_deg, oli_deg, mic_deg))
# stopifnot(length(dlpfc_scRNAseq_deg) == 1031)

snDEG <- dlpfc_scRNAseq_deg
write.csv(dlpfc_scRNAseq_deg, paste("results/", "dlpfc_snRNAseq_deg.csv", sep = ""), 
          row.names = FALSE)

find_ciDEG_in_DC <- function(BR, ct_correction){
  
  if(ct_correction == TRUE){
    ct_file_name <- "_ct_correction"
  } else {
    ct_file_name = ""
  }
  
  dc_edges <- read.csv(paste("results/", BR[1], BR[2], ct_file_name, "Mapped_DCEdgesWithoutNA_fdr_", fdr_val, ".csv", sep = ""),
                       stringsAsFactors = FALSE)
  # Since BM10 is the region scRNA-seq data is available for
  dc_genes <- unique(dc_edges$g1)
  
  # Cell-type specific overlap for BM10 genes
  ast_overlap_g1 <- intersect(ast_deg, dc_genes)
  mic_overlap_g1 <- intersect(mic_deg, dc_genes)
  neu_overlap_g1 <- intersect(unique(c(ex_deg, in_deg)), dc_genes)
  oli_overlap_g1 <- intersect(oli_deg, dc_genes)
  total_overlap <- intersect(snDEG, dc_genes)
  
  library(BRETIGEA)
  marker_genes <- markers_df_human_brain$markers[which(markers_df_human_brain$cell %in%
                                                         c("ast", "mic", "neu", "oli"))]
  ast_markers <- markers_df_human_brain$markers[which(markers_df_human_brain$cell == "ast")]
  mic_markers <- markers_df_human_brain$markers[which(markers_df_human_brain$cell == "mic")]
  neu_markers <- markers_df_human_brain$markers[which(markers_df_human_brain$cell == "neu")]
  oli_markers <- markers_df_human_brain$markers[which(markers_df_human_brain$cell == "oli")]
  
  ast_overlap_g1_m <- intersect(ast_overlap_g1, ast_markers)
  mic_overlap_g1_m <- intersect(mic_overlap_g1, mic_markers)
  neu_overlap_g1_m <- intersect(neu_overlap_g1, neu_markers)
  oli_overlap_g1_m <- intersect(oli_overlap_g1, oli_markers)
  total_overlap_m <- intersect(total_overlap, marker_genes)
  
  cat("For g1", "\n",
      "Ast", length(ast_overlap_g1), "(", length(ast_overlap_g1_m), ")", "\n",
      "Mic", length(mic_overlap_g1), "(", length(mic_overlap_g1_m), ")", "\n",
      "Ex and In Neu", length(neu_overlap_g1), "(", length(neu_overlap_g1_m), ")", "\n",
      "Oli", length(oli_overlap_g1), "(", length(oli_overlap_g1_m), ")","\n",
      "Total", length(total_overlap), "(", length(total_overlap_m), ")", "\n"
  )
  
  # Finding the overlap of the marker genes with unique BM10 DC genes
  # a + b + c + d = Number of BM10 DC genes
  # b + d = Number of Marker genes among BM10 DC genes
  
  function_round_to_2 <- function(x){
    return(format(round(x, 2), nsmall = 2))
  }
  
  ast_b_d <- intersect(ast_markers, dc_genes)
  mic_b_d <- intersect(mic_markers, dc_genes)
  neu_b_d <- intersect(neu_markers, dc_genes)
  oli_b_d <- intersect(oli_markers, dc_genes)
  total_b_d <- intersect(marker_genes, dc_genes)
  
  cat("d/(b+d)%", "\n",
      "Ast", function_round_to_2(length(ast_overlap_g1_m)*100/length(ast_b_d)), "\n",
      "Mic", function_round_to_2(length(mic_overlap_g1_m)*100/length(mic_b_d)),"\n",
      "Neu", function_round_to_2(length(neu_overlap_g1_m)*100/length(neu_b_d)),"\n",
      "Oli", function_round_to_2(length(oli_overlap_g1_m)*100/length(oli_b_d)),"\n",
      "Total", function_round_to_2(length(total_overlap_m)*100/length(total_b_d)), "\n")
  
  universe <- intersect(dc_genes, unique(c(marker_genes, total_overlap)))
  
  ast_b_d <- intersect(ast_markers, universe)
  mic_b_d <- intersect(mic_markers, universe)
  neu_b_d <- intersect(neu_markers, universe)
  oli_b_d <- intersect(oli_markers, universe)
  total_b_d <- intersect(marker_genes, universe)
  
  cat("d/(b+d)%", "\n",
      "Ast", function_round_to_2(length(ast_overlap_g1_m)*100/length(ast_b_d)), "\n",
      "Mic", function_round_to_2(length(mic_overlap_g1_m)*100/length(mic_b_d)),"\n",
      "Neu", function_round_to_2(length(neu_overlap_g1_m)*100/length(neu_b_d)),"\n",
      "Oli", function_round_to_2(length(oli_overlap_g1_m)*100/length(oli_b_d)),"\n",
      "Total", function_round_to_2(length(total_overlap_m)*100/length(total_b_d)), "\n")
  
}

# NOTE: We are only considering the BM10 combinations because the
# dlpfc DEG are BM10 region's DEG
find_ciDEG_in_DC(c("BM10", "BM22"), ct_correction = FALSE)
find_ciDEG_in_DC(c("BM10", "BM36"), ct_correction = FALSE)
find_ciDEG_in_DC(c("BM10", "BM44"), ct_correction = FALSE)

find_ciDEG_in_DC(c("BM10", "BM22"), ct_correction = TRUE)
find_ciDEG_in_DC(c("BM10", "BM36"), ct_correction = TRUE)
find_ciDEG_in_DC(c("BM10", "BM44"), ct_correction = TRUE)

markers_in_ci_deg <- markers_df_human_brain[which(markers_df_human_brain$markers %in% dlpfc_scRNAseq_deg),]
markers_in_ci_deg_4_ct <- markers_in_ci_deg[which(markers_in_ci_deg$cell %in% c("ast", "mic", "neu", "oli")),]
dim(markers_in_ci_deg_4_ct)

# Getting DEG results of the community of interest ------------------------

community <- read.csv(paste("results/", "Community -715 genes-DE values.csv", sep = ""),
                      stringsAsFactors = FALSE)
bm10_community_genes <- community$BM10[which(community$BM10 != "")]
bm36_community_genes <- community$BM36[which(community$BM36 != "")]

bm10_deg <- read.csv(paste("results/", "BM10_DEG_results_limma-voom.csv", sep = ""),
                     stringsAsFactors = FALSE)
bm36_deg <- read.csv(paste("results/", "BM36_DEG_results_limma-voom.csv", sep = ""),
                     stringsAsFactors = FALSE)

bm10_deg_community <- bm10_deg[which(bm10_deg$GeneSymbol %in% bm10_community_genes),]
bm36_deg_community <- bm36_deg[which(bm36_deg$GeneSymbol %in% bm36_community_genes),]

write.csv(bm10_deg_community, paste("results/", "bm10_deg_community_genes.csv",sep = ""))
write.csv(bm36_deg_community, paste("results/", "bm36_deg_community_genes.csv",sep = ""))

which(bm10_deg_community$adj.P.Val < 0.05)
length(which(bm36_deg_community$adj.P.Val < 0.05))

bm10_community_ciDEG <- rbind(ex_deg[which(row.names(ex_deg) %in% bm10_community_genes),],
                              in_deg[which(row.names(in_deg) %in% bm10_community_genes),],
                              mic_deg[which(row.names(mic_deg) %in% bm10_community_genes),],
                              ast_deg[which(row.names(ast_deg) %in% bm10_community_genes),],
                              oli_deg[which(row.names(oli_deg) %in% bm10_community_genes),],
                              opc_deg[which(row.names(opc_deg) %in% bm10_community_genes),])

bm36_community_ciDEG <- rbind(ex_deg[which(row.names(ex_deg) %in% bm36_community_genes),],
                              in_deg[which(row.names(in_deg) %in% bm36_community_genes),],
                              mic_deg[which(row.names(mic_deg) %in% bm36_community_genes),],
                              ast_deg[which(row.names(ast_deg) %in% bm36_community_genes),],
                              oli_deg[which(row.names(oli_deg) %in% bm36_community_genes),],
                              opc_deg[which(row.names(opc_deg) %in% bm36_community_genes),])

write.csv(bm10_community_ciDEG, paste("results/", "bm10_community_genes_ciDEG.csv", sep = ""))
write.csv(bm36_community_ciDEG, paste("results/", "bm36_community_genes_ciDEG.csv", sep = ""))
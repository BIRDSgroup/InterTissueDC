# Author: B P Kailash
# Purpose: Analyze the DC relations observed - Part 2

# 4) Add braak scores to the sample ids
# 5) Find the dlpfc snRNAseq DEG within a cell type 
# 6) Find the intersection of dlpfc in the BM10 combinations

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

label_with_braak(c("BM_10", "BM_22"))
label_with_braak(c("BM_10", "BM_36"))
label_with_braak(c("BM_10", "BM_44"))
label_with_braak(c("BM_22", "BM_36"))
label_with_braak(c("BM_22", "BM_44"))
label_with_braak(c("BM_36", "BM_44"))

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
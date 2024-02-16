library(tidyverse)

find_overlap_with_within_dc_per_br_all_gene_pairs <- function(BR){
  within_region <- read_csv(paste0("results/after_ctc_within_region/", BR, BR, "_ct_correctionMapped_DCEdgesWithoutNA.csv"))
  dim(within_region)
  all_dc_edges <-  read_csv(paste0("results/", BR, "_BR2_allDC.csv"))
  dim(all_dc_edges)
  
  all_dc_edges$edge <- paste(all_dc_edges$BR1, all_dc_edges$BR2)
  
  within_region$edge <- paste(within_region$g1, within_region$g2)
  all_dc_edges_ab <- left_join(all_dc_edges, within_region, by = "edge")
  
  within_region$edge <- paste(within_region$g2, within_region$g1)
  all_dc_edges_ba <- left_join(all_dc_edges, within_region, by = "edge")
  
  all_dc_edges <- rbind(all_dc_edges_ab, all_dc_edges_ba)
  
  
  write_csv(na.omit(all_dc_edges), paste0("results/", BR, "_BR2_overlap_with_within_DC.csv"))
}

find_overlap_with_within_dc_per_br_all_gene_pairs(BR = "BM10")
find_overlap_with_within_dc_per_br_all_gene_pairs(BR = "BM22")
find_overlap_with_within_dc_per_br_all_gene_pairs(BR = "BM36")
find_overlap_with_within_dc_per_br_all_gene_pairs(BR = "BM44")


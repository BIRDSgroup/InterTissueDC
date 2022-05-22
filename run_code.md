# How to obtain Differentially Correlated (DC) edges after Cell-type correction (CTC)
## Workflow

This document details the analysis workflow. A walkthrough of the code used in this study is presented as a tutorial to ensure reproducibility.

This markdown file gives a step by step guide to obtain all the results obtained in this study for a particular brain region combination BM36 - BM44. The process is identical for studying any of the five other brain region combinations. 

## Server details

All the code used in this study was run on the IBSE Actinium server (http://10.24.6.51:8787) with the following specifications. 

* Processor : 4x Intel Xeon E7-4850V4 family (16 Cores, 32 Threads, Processor
base frequency : 2.1GHz)
* Cache : 40MB
* Memory : 32x 32GB DDR4 2400MHz ECC LRDIMM (1TB total)
* Disk Space: 10x 4TB 3.5” 7.2K RPM SATA Hard Disk Drive
* Disk Space: 2x 1.2TB 2.5” MLC SATA Solid-State Drive
* Architecture: x86_64

Details can also be found at session_info.txt.

## run find_DC.R

The find_dc code identifies the significantly correlated edges, and significantly differentially correlated edges across brain tissues/regions.

### INPUT
The code requires the following input from the user. The name of the server (IBSE), brain region combination under study - BM10-BM22, BM10-BM36, BM10-BM44, BM22-BM36, BM22-BM44, or BM36-BM44. NOTE: The order cannot change, e.g. BM44-BM36 != BM36-BM44, and the FDR cutoff for identfying statistically significant correlations.

```{r}
server = "IBSE"
BR <- c("BM36", "BM44")
fdr_val <- 0.01
```

### OUTPUT
The code generates the following output files, e.g. for BM36-BM44.
1. CTL_idsBM36BM44.csv
2. AD_idsBM36BM44.csv
3. ADBM36_BM36BM44.csv
4. CTLBM36_BM36BM44.csv
5. ADBM44_BM36BM44.csv
6. CTLBM44_BM36BM44.csv
7. BM36BM44ADCorrelatedEdges.csv
8. BM36BM44CTLCorrelatedEdges.csv
9. UnionCorrelatedEdgesBM36BM44.csv
10. p_vals_without_adjustmentBM36BM44Edges0.01.rda
11. DCEdgesBM36BM44.csv

**NOTE:** Please find below the estimate running time for the find_dc code on RStudio in a server with previously mentioned specifications.

### Time taken

| Brain regions | Time taken (hrs) |
| ----------- | ----------- |
| BM10-BM22 | 9.5 |
| BM10-BM36 | 1 |
| BM10-BM44 | 0.6 |
| BM22-BM36 | 15 |
| BM222-BM44 | 11 |
| BM36-BM44 | 0.33 |

## run cell_type_correction.R

This part of the code predicts the cell-type frequencies using the CellCODE 80 mg 4 ct model (see BRETIGEAvsCellCODEanalysis.R to learn how we chose this model). Cellular deconvolution models are generated for each brain region. Based on the CellCODE estimated cell-type frequencies, each gene's expression is corrected for to reflect only the cell-intrinsic changes.

Load the read_BR_AD_CTL, readCellTypeFreq and predict_gExp functions.

### CODE OUTLINE
```{r}
# CellCODE 80 mg, 4 cell types seems to be the best performing model
ip <- make_w_matrix(markers_df_human_brain, c("ast", "mic", "neu", "oli"), c(20,20,20,20))
w <- ip[[1]]
mg <- ip[[2]]

BM36_spv <- get_SPV_BrainRegion("BM36", MarkerGenes = mg)
BM44_spv <- get_SPV_BrainRegion("BM44", MarkerGenes = mg)

bm36_geneExp <- read_BR_AD_CTL("BM36")
bm36_sample_gp <- bm36_geneExp$sample_gp
bm36_geneExp <- subset(bm36_geneExp, select = -c(sample_gp))

bm44_geneExp <- read_BR_AD_CTL("BM44")
bm44_sample_gp <- bm44_geneExp$sample_gp
bm44_geneExp <- subset(bm44_geneExp, select = -c(sample_gp))

bm36_pred <- predict_gExp(bm36_geneExp, BM36_spv)
bm44_pred <- predict_gExp(bm44_geneExp, BM44_spv)

bm36_corrected <- bm36_geneExp - t(bm36_pred)
bm44_corrected <- bm44_geneExp - t(bm44_pred)
```

### OUTPUT

1. BM36_CellCODE_80mg_4cellTypes.csv
2. BM44_CellCODE_80mg_4cellTypes.csv
3. BM36_AD_CTL_cellTypeCorrected.csv
4. BM44_AD_CTL_cellTypeCorrected.csv
5. ADBM36_BM36BM44_cellTypeCorrected.csv
6. ADBM44_BM36BM44_cellTypeCorrected.csv
7. CTLBM36_BM36BM44_cellTypeCorrected.csv
8. CTLBM44_BM36BM44_cellTypeCorrected.csv

### Time taken (estimate)
1.5 hours.

## run find_DC_after_ct_correction.R

### INPUT

```{r}
server = "IBSE"
BR <- c("BM36", "BM44")
fdr_val <- 0.01
```

### OUTPUT

1. BM36BM44CTCorrectedADCorrelatedEdges.csv
2. BM36BM44CTCorrectedCTLCorrelatedEdges.csv
3. UnionCorrelatedEdges_ct_correctionBM36BM44.csv
4. p_vals_without_adjustment_ct_correctionBM36BM44Edges0.01.rda
5. DCEdges_ct_correctionBM36BM44.csv

### Time taken

| Brain regions | Time taken (mins) |
| ----------- | ----------- |
| BM10-BM22 | 20 |
| BM10-BM36 | 19 |
| BM10-BM44 | 18 |
| BM22-BM36 | 21 |
| BM222-BM44 | 18 |
| BM36-BM44 | 19 |

## run analyze_dc.R
The user needs to specify
* path_results
* path_data
* cell_types
* number_mg  

Then the following functions also present in the code needs to be loaded
* ensembl_to_approved_gene_symbol_corr
* ensembl_to_approved_gene_symbol_dc
* label_DC_edges
* count_interactions
* label_with_braak

then run as specified below.

```{r}
path_results <- "/data/kailash/DDP/results/"
path_data <- "/data/kailash/DDP/data/"
mapping_file <- read.csv(paste(path_data, 
                               "Mapped_h37_ens_to_gene.csv", 
                               sep = ""), 
                         stringsAsFactors = FALSE)


ensembl_to_approved_gene_symbol_corr(MappingFile = mapping_file, 
                                     PathResults = path_results, 
                                     BR = c("BM36", "BM44"), ct_corrected = FALSE)
ensembl_to_approved_gene_symbol_corr(MappingFile = mapping_file, 
                                     PathResults = path_results, 
                                     BR = c("BM36", "BM44"), ct_corrected = TRUE)

ensembl_to_approved_gene_symbol_dc(
  MappingFile = mapping_file, 
  PathResults = path_results, 
  BR = c("BM36", "BM44"), ct_corrected = FALSE)
ensembl_to_approved_gene_symbol_dc(
  MappingFile = mapping_file, 
  PathResults = path_results, 
  BR = c("BM36", "BM44"), ct_corrected = TRUE)

cell_types <- c("ast", "mic", "neu", "oli")
number_mg <- 1000
number_mg_per_ct <- rep(number_mg, length(cell_types))

label_DC_edges(BR = c("BM36",  "BM44"), CellTypes = cell_types, NumberOfMarkersPerCellType = number_mg_per_ct, CTCorrection = FALSE)
label_DC_edges(BR = c("BM36",  "BM44"), CellTypes = cell_types, NumberOfMarkersPerCellType = number_mg_per_ct, CTCorrection = TRUE)
count_interactions(c("BM36", "BM44"), ct_correction =  FALSE)

metadata <- read.csv(paste(path_data, "rna_metadata_onlyBAM_noExclude_noDup_Nov2018.csv", sep = ""), stringsAsFactors = FALSE)
label_with_braak(c("BM36", "BM44"))
```

We are using all the 1000 marker genes for each cell type (ast, mic, neu, oli) identified in the BRETIGEA study to label the DC edges.

### OUTPUT

1. BM36BM44Mapped_CorrelatedEdgesWithoutNA.csv
2. BM36BM44_ct_correctionMapped_CorrelatedEdgesWithoutNA.csv
3. BM36BM44Mapped_DCEdgesWithoutNA.csv
4. BM36BM44_ct_correctionMapped_DCEdgesWithoutNA.csv
5. CTL_idsBM36BM44.csv with braak labelling
6. AD_idsBM36BM44.csv with braak labelling

## run_find_deg.R

### INPUT

```r
path_data <- "/data/kailash/DDP/data/"
path_geneExp <- "/data/srivatsan/Original_Expression_Files/Original_"
path_results <- "/data/kailash/DDP/results/"
BR <- "BM36"
ct_corrected <- TRUE
````

### OUTPUT

Change input parameters to obtain the following output files
1. BM36_DEG_WilcoxonRankSumTest.csv
2. BM44_DEG_WilcoxonRankSumTest.csv
3. BM36cellTypeCorrected_DEG_WilcoxonRankSumTest.csv
4. BM44cellTypeCorrected_DEG_WilcoxonRankSumTest.csv

###### Please feel free to raise an issue on GitHub for clarifications.

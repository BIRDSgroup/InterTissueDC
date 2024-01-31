# How to run the code? - Description of the workflow

This document details the analysis workflow. A walkthrough of the code used in this study is presented as a tutorial to ensure reproducibility.

## run find_DC.R

The find_dc code identifies the significantly correlated edges, and significantly differentially correlated edges across brain tissues/regions.

### INPUT
The code requires the following input from the user. The brain region combination under study - e.g. BM_36, BM_44. NOTE: The order cannot change, e.g. BM44-BM36 != BM36-BM44, and the FDR cutoff for identfying statistically significant correlations.

```
Rscript scripts/find_DC.R --br1 BM_10 --br2 BM_22 --fdr 0.01
# bsub -n 8 -R "rusage[mem=16000]" -R span[hosts=1] -W 12:00 -eo %J.err -oo %J.out

Rscript scripts/find_DC.R --br1 BM_10 --br2 BM_36 --fdr 0.01
# bsub -n 8 -R "rusage[mem=16000]" -R span[hosts=1] -W 12:00 -eo %J.err -oo %J.out

Rscript scripts/find_DC.R --br1 BM_10 --br2 BM_44 --fdr 0.01
# bsub -n 8 -R "rusage[mem=16000]" -R span[hosts=1] -W 12:00 -eo %J.err -oo %J.out

Rscript scripts/find_DC.R --br1 BM_22 --br2 BM_36 --fdr 0.01
# bsub -n 8 -R "rusage[mem=16000]" -R span[hosts=1] -W 12:00 -eo %J.err -oo %J.out

Rscript scripts/find_DC.R --br1 BM_22 --br2 BM_44 --fdr 0.01
# bsub -n 8 -R "rusage[mem=16000]" -R span[hosts=1] -W 12:00 -eo %J.err -oo %J.out

Rscript scripts/find_DC.R --br1 BM_36 --br2 BM_44 --fdr 0.01
# bsub -n 8 -R "rusage[mem=16000]" -R span[hosts=1] -W 12:00 -eo %J.err -oo %J.out

Rscript scripts/find_DC.R --br1 BM_10 --br2 BM_10 --fdr 0.01
# bsub -n 8 -R "rusage[mem=16000]" -R span[hosts=1] -W 50:00 -eo %J.err -oo %J.out

Rscript scripts/find_DC.R --br1 BM_22 --br2 BM_22 --fdr 0.01
# bsub -n 8 -R "rusage[mem=16000]" -R span[hosts=1] -W 50:00 -eo %J.err -oo %J.out

Rscript scripts/find_DC.R --br1 BM_36 --br2 BM_36 --fdr 0.01
# bsub -n 10 -R "rusage[mem=16000]" -R span[hosts=1] -W 50:00 -eo %J.err -oo %J.out

Rscript scripts/find_DC.R --br1 BM_44 --br2 BM_44 --fdr 0.01
# bsub -n 8 -R "rusage[mem=16000]" -R span[hosts=1] -W 60:00 -eo %J.err -oo %J.out
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

## run cell_type_correction.R

This part of the code predicts the cell-type frequencies using the CellCODE 80 mg 4 ct model (see BRETIGEAvsCellCODEanalysis.R to learn how we chose this model). Cellular deconvolution models are generated for each brain region. Based on the CellCODE estimated cell-type frequencies, each gene's expression is corrected for to reflect only the cell-intrinsic changes.

### RUN
```{r}
Rscript scripts/cell_type_correction.R
```

### OUTPUT example for BM36, BM44 combination.

1. BM36_CellCODE_80mg_4cellTypes.csv
2. BM44_CellCODE_80mg_4cellTypes.csv
3. BM36_AD_CTL_cellTypeCorrected.csv
4. BM44_AD_CTL_cellTypeCorrected.csv
5. ADBM36_BM36BM44_cellTypeCorrected.csv
6. ADBM44_BM36BM44_cellTypeCorrected.csv
7. CTLBM36_BM36BM44_cellTypeCorrected.csv
8. CTLBM44_BM36BM44_cellTypeCorrected.csv

## run find_DC_after_ct_correction.R

The code requires the following input from the user. The brain region combination under study - e.g. BM_36, BM_44. NOTE: The order cannot change, e.g. BM44-BM36 != BM36-BM44, and the FDR cutoff for identfying statistically significant correlations.

```
Rscript scripts/find_dc_after_ct_correction.R --br1 BM10 --br2 BM22 --fdr 0.01
# bsub -n 8 -R "rusage[mem=16000]" -R span[hosts=1] -W 12:00 -eo %J.err -oo %J.out

Rscript scripts/find_dc_after_ct_correction.R --br1 BM10 --br2 BM36 --fdr 0.01
# bsub -n 8 -R "rusage[mem=16000]" -R span[hosts=1] -W 12:00 -eo %J.err -oo %J.out

Rscript scripts/find_dc_after_ct_correction.R --br1 BM10 --br2 BM44 --fdr 0.01
# bsub -n 8 -R "rusage[mem=16000]" -R span[hosts=1] -W 12:00 -eo %J.err -oo %J.out

Rscript scripts/find_dc_after_ct_correction.R --br1 BM22 --br2 BM36 --fdr 0.01
# bsub -n 8 -R "rusage[mem=16000]" -R span[hosts=1] -W 12:00 -eo %J.err -oo %J.out

Rscript scripts/find_dc_after_ct_correction.R --br1 BM22 --br2 BM44 --fdr 0.01
# bsub -n 8 -R "rusage[mem=16000]" -R span[hosts=1] -W 12:00 -eo %J.err -oo %J.out

Rscript scripts/find_dc_after_ct_correction.R --br1 BM36 --br2 BM44 --fdr 0.01
# bsub -n 8 -R "rusage[mem=16000]" -R span[hosts=1] -W 12:00 -eo %J.err -oo %J.out

Rscript scripts/find_dc_after_ct_correction.R --br1 BM10 --br2 BM10 --fdr 0.01
# bsub -n 8 -R "rusage[mem=16000]" -R span[hosts=1] -W 36:00 -eo %J.err -oo %J.out

Rscript scripts/find_dc_after_ct_correction.R --br1 BM22 --br2 BM22 --fdr 0.01
# bsub -n 8 -R "rusage[mem=16000]" -R span[hosts=1] -W 36:00 -eo %J.err -oo %J.out

Rscript scripts/find_dc_after_ct_correction.R --br1 BM36 --br2 BM36 --fdr 0.01
# bsub -n 10 -R "rusage[mem=16000]" -R span[hosts=1] -W 36:00 -eo %J.err -oo %J.out

Rscript scripts/find_dc_after_ct_correction.R --br1 BM44 --br2 BM44 --fdr 0.01
# bsub -n 8 -R "rusage[mem=16000]" -R span[hosts=1] -W 36:00 -eo %J.err -oo %J.out
```

### OUTPUT example for BM36, BM44 combination.

1. BM36BM44CTCorrectedADCorrelatedEdges.csv
2. BM36BM44CTCorrectedCTLCorrelatedEdges.csv
3. UnionCorrelatedEdges_ct_correctionBM36BM44.csv
4. p_vals_without_adjustment_ct_correctionBM36BM44Edges0.01.rda
5. DCEdges_ct_correctionBM36BM44.csv

## run_find_deg.R
Find differentially expressed genes between "Definite AD" and "Control" for each brain region.

```
Rscript scripts/find_deg.R --br BM_10 --ctc FALSE
Rscript scripts/find_deg.R --br BM_22 --ctc FALSE
Rscript scripts/find_deg.R --br BM_36 --ctc FALSE
Rscript scripts/find_deg.R --br BM_44 --ctc FALSE

Rscript scripts/find_deg.R --br BM10 --ctc TRUE
Rscript scripts/find_deg.R --br BM22 --ctc TRUE
Rscript scripts/find_deg.R --br BM36 --ctc TRUE
Rscript scripts/find_deg.R --br BM44 --ctc TRUE
```
### OUTPUT example for BM36

Change input parameters to obtain the following output files
1. BM36_DEG_WilcoxonRankSumTest.csv
2. BM36cellTypeCorrected_DEG_WilcoxonRankSumTest.csv

## run analyze_dc.R

The user needs to specify. We are using all the 1000 marker genes for each cell type (ast, mic, neu, oli) identified in the BRETIGEA study to label the DC edges.

* cell_types (default: cell_types <- c("ast", "mic", "neu", "oli"))
* number_mg  (default: number_mg <- 1000)

then run as specified below.

```
Rscript scripts/analyze_DC.R
```

### OUTPUT example for BM36, BM44 combination.

1. BM36BM44Mapped_CorrelatedEdgesWithoutNA.csv
2. BM36BM44_ct_correctionMapped_CorrelatedEdgesWithoutNA.csv
3. BM36BM44Mapped_DCEdgesWithoutNA.csv
4. BM36BM44_ct_correctionMapped_DCEdgesWithoutNA.csv
5. CTL_idsBM36BM44.csv with braak labelling
6. AD_idsBM36BM44.csv with braak labelling

###### Please feel free to raise an issue on GitHub for clarifications.

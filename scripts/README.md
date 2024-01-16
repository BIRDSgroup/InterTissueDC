## run BRETIGEAvsCellCODEanalysis.R

* Analyzing different models of BRETIGEA and CellCODE, CellCODE 80 marker gene (mg), four cell types (ct) (ast, mic, neu, oli) was chosen to be the best model - see CellCODE analysis, BRETIGEA analysis and Making the ggpairs plot sections of the code.
* Go to section **Getting the models** in the code
   + load functions make_w_matrix, get_SPV_BrainRegion present in the code
   + run the following lines of code

```{r}
ip <- make_w_matrix(markers_df_human_brain, c("ast", "mic", "neu", "oli"), c(20,20,20,20))
w <- ip[[1]]
MG <- ip[[2]]

BM36_spv <- get_SPV_BrainRegion("BM36", MarkerGenes = MG)
write.csv(BM36_spv,
          paste(path_results, "BM36_CellCODE_80mg_4cellTypes.csv", sep = ""))
BM44_spv <- get_SPV_BrainRegion("BM44", MarkerGenes = mg)
write.csv(BM44_spv,
          paste(path_results, "BM44_CellCODE_80mg_4cellTypes.csv", sep = ""))
```

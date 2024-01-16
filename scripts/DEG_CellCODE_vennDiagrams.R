rm(list = ls())

# Loading libraries -------------------------------------------------------

path_data <- "/data/kailash/DDP/data/"
  
.libPaths(c("/data/kailash/DDP/libraries/", .libPaths()))

source("/data/kailash/DDP/src/main/scripts/libraries.R")
source("/data/kailash/DDP/src/main/scripts/utils.R")

# CellCODE Venn Diagram ---------------------------------------------------

path_results <- "/data/kailash/DDP/results/"
setwd(path_results)

f1 <- read.csv("BM10U_SVD_CellCODE.csv", row.names = 1, stringsAsFactors = FALSE)
f1$BrainRegion <- rep("BM10", dim(f1)[1])
f2 <- read.csv("BM22U_SVD_CellCODE.csv", row.names = 1, stringsAsFactors = FALSE)
f2$BrainRegion <- rep("BM22", dim(f2)[1])
f3 <- read.csv("BM36U_SVD_CellCODE.csv", row.names = 1, stringsAsFactors = FALSE)
f3$BrainRegion <- rep("BM36", dim(f3)[1])
f4 <- read.csv("BM44U_SVD_CellCODE.csv", row.names = 1, stringsAsFactors = FALSE)
f4$BrainRegion <- rep("BM44", dim(f4)[1])

# Violin plot -------------------------------------------------------------

f <- rbind(f1, f2, f3, f4)
head(f)

ggplot(f, aes(x = BrainRegion, y = su, color = BrainRegion)) + geom_violin() +
  labs(x = "Brain Region", y = "CellCODE gene coefficient") +
  geom_jitter(position=position_jitter(0.2))


# Scatter plots -----------------------------------------------------------

plot_scatter <- function(f1, f2){
  scatter1 <- data.frame(x1 = f1, x2 = f2)
  plot.new()
  ggplot(scatter1, aes(x = x1.su, y = x2.su)) + geom_point() +
    labs(x = scatter1$x1.BrainRegion[1] , y = scatter1$x2.BrainRegion[1]) +
    title("Scatter plot of the CellCODE gene coefficients") +
    scale_x_continuous(limits = c(-1.0, 0.8)) +
    scale_y_continuous(limits = c(-1.0, 0.8))
}

plot_scatter(f1, f2)
plot_scatter(f1, f3)
plot_scatter(f1, f4)
plot_scatter(f2, f3)
plot_scatter(f2, f4)
plot_scatter(f3, f4)

dim(f1)
dim(f2)
dim(f3)
dim(f4)

length(which(f1$su != 0))
length(which(f2$su != 0))
length(which(f3$su != 0))
length(which(f4$su != 0))

f1[which(f1$su != 0),]
f2[which(f2$su != 0),]
f3[which(f3$su != 0),]
f4[which(f4$su != 0),]

x1 <- f1$genes[which(f1$su != 0)]
x2 <- f2$genes[which(f2$su != 0)]
x3 <- f3$genes[which(f3$su != 0)]
x4 <- f4$genes[which(f4$su != 0)]

length(which(x1 %in% x2))

library(VennDiagram)
library(RColorBrewer)

colours <- brewer.pal(n = 8, name = 'Dark2')
  
venn.diagram(
  x = list(x1, x2, x3, x4),
  category.names = c("BM10" , "BM22 " , "BM36", "BM44"),
  filename = "CellCODE_coeffs_genes.png",
  output = TRUE,
  
  fill = c(alpha("#1B9E77",0.3), alpha("#D95F02",0.3), alpha("#7570B3",0.3), alpha("#E7298A",0.3)),
  fontfamily = "sans",
  cat.default.pos = "outer",
  cat.fontfamily = "sans",
  cat.col = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A"),
)

count_cell_types <- function(f){
  f <- f[which(f$su != 0),]
  cat("Number of ast markers used ", length(which(f$cellType == "ast")), "\n")
  cat("Number of mic markers used ", length(which(f$cellType == "mic")), "\n")
  cat("Number of neu markers used ", length(which(f$cellType == "neu")), "\n")
  cat("Number of oli markers used ", length(which(f$cellType == "oli")), "\n")
}

count_cell_types(f1)
count_cell_types(f2)
count_cell_types(f3)
count_cell_types(f4)

# DGE analysis Venn Diagram -----------------------------------------------

path_results <- "/data/kailash/DDP/results/"


# Reading all the DEG files -----------------------------------------------

ct_file_name <- "cellTypeCorrected"
bm10_deg_limma <- read.csv(paste(path_results, "BM10", ct_file_name, "_DGE_results_limma_moderated_t_test.csv", sep = ""),
                           stringsAsFactors = FALSE)
bm22_deg_limma <- read.csv(paste(path_results, "BM22", ct_file_name, "_DGE_results_limma_moderated_t_test.csv", sep = ""),
                           stringsAsFactors = FALSE)
bm36_deg_limma <- read.csv(paste(path_results, "BM36", ct_file_name, "_DGE_results_limma_moderated_t_test.csv", sep = ""),
                           stringsAsFactors = FALSE)
bm44_deg_limma <- read.csv(paste(path_results, "BM44", ct_file_name, "_DGE_results_limma_moderated_t_test.csv", sep = ""),
                           stringsAsFactors = FALSE)

ct_file_name <- "_ct_corrected"
bm10_deg_limma_voom <- read.csv(paste(path_results, "BM10", ct_file_name, "_DGE_results_limma-voom.csv", sep = ""),
                                stringsAsFactors = FALSE)
bm22_deg_limma_voom <- read.csv(paste(path_results, "BM22", ct_file_name, "_DGE_results_limma-voom.csv", sep = ""),
                                stringsAsFactors = FALSE)
bm36_deg_limma_voom <- read.csv(paste(path_results, "BM36", ct_file_name, "_DGE_results_limma-voom.csv", sep = ""),
                                stringsAsFactors = FALSE)
bm44_deg_limma_voom <- read.csv(paste(path_results, "BM44", ct_file_name, "_DGE_results_limma-voom.csv", sep = ""),
                                stringsAsFactors = FALSE)

ct_file_name = "cellTypeCorrected"
bm10_ct_deg_wilcoxon <- read.csv(paste(path_results, "BM10", ct_file_name, "_DEG_WilcoxonRankSumTest.csv", sep = ""),
                                 stringsAsFactors = FALSE)
bm22_ct_deg_wilcoxon <- read.csv(paste(path_results, "BM22", ct_file_name, "_DEG_WilcoxonRankSumTest.csv", sep = ""),
                                 stringsAsFactors = FALSE)
bm36_ct_deg_wilcoxon <- read.csv(paste(path_results, "BM36", ct_file_name, "_DEG_WilcoxonRankSumTest.csv", sep = ""),
                                 stringsAsFactors = FALSE)
bm44_ct_deg_wilcoxon <- read.csv(paste(path_results, "BM44", ct_file_name, "_DEG_WilcoxonRankSumTest.csv", sep = ""),
                                 stringsAsFactors = FALSE)


# Checking number of DEG --------------------------------------------------

print_number_of_deg_limma <- function(ct_file_name, fdr){
  bm10_deg_limma <- read.csv(paste(path_results, "BM10", ct_file_name, "_DGE_results_limma_moderated_t_test.csv", sep = ""),
                             stringsAsFactors = FALSE)
  bm22_deg_limma <- read.csv(paste(path_results, "BM22", ct_file_name, "_DGE_results_limma_moderated_t_test.csv", sep = ""),
                             stringsAsFactors = FALSE)
  bm36_deg_limma <- read.csv(paste(path_results, "BM36", ct_file_name, "_DGE_results_limma_moderated_t_test.csv", sep = ""),
                             stringsAsFactors = FALSE)
  bm44_deg_limma <- read.csv(paste(path_results, "BM44", ct_file_name, "_DGE_results_limma_moderated_t_test.csv", sep = ""),
                             stringsAsFactors = FALSE)
  
  cat("BM10 ", length(which(bm10_deg_limma$adj.P.Val <= fdr)), "\n")
  cat("BM22 ", length(which(bm22_deg_limma$adj.P.Val <= fdr)), "\n")
  cat("BM36 ", length(which(bm36_deg_limma$adj.P.Val <= fdr)), "\n")
  cat("BM44 ", length(which(bm44_deg_limma$adj.P.Val <= fdr)), "\n")
  
}

print_number_of_deg_limma(fdr = 0.1, ct_file_name = "cellTypeCorrected")
print_number_of_deg_limma(fdr = 0.2, ct_file_name = "cellTypeCorrected")
print_number_of_deg_limma(fdr = 0.1, ct_file_name = "")

# Checking number of DEG limma-voom ---------------------------------------

print_number_of_deg_limma_voom <- function(ct_file_name, fdr){
  bm10_deg_limma_voom <- read.csv(paste(path_results, "BM10", ct_file_name, "_DGE_results_limma-voom.csv", sep = ""),
                             stringsAsFactors = FALSE)
  bm22_deg_limma_voom <- read.csv(paste(path_results, "BM22", ct_file_name, "_DGE_results_limma-voom.csv", sep = ""),
                             stringsAsFactors = FALSE)
  bm36_deg_limma_voom <- read.csv(paste(path_results, "BM36", ct_file_name, "_DGE_results_limma-voom.csv", sep = ""),
                             stringsAsFactors = FALSE)
  bm44_deg_limma_voom <- read.csv(paste(path_results, "BM44", ct_file_name, "_DGE_results_limma-voom.csv", sep = ""),
                             stringsAsFactors = FALSE)
  
  cat("BM10 ", length(which(bm10_deg_limma_voom$adj.P.Val <= fdr)), "\n")
  cat("BM22 ", length(which(bm22_deg_limma_voom$adj.P.Val <= fdr)), "\n")
  cat("BM36 ", length(which(bm36_deg_limma_voom$adj.P.Val <= fdr)), "\n")
  cat("BM44 ", length(which(bm44_deg_limma_voom$adj.P.Val <= fdr)), "\n")
  
}

print_number_of_deg_limma_voom(fdr = 0.1, ct_file_name = "_ct_corrected")
print_number_of_deg_limma_voom(fdr = 0.2, ct_file_name = "_ct_corrected")
print_number_of_deg_limma_voom(fdr = 0.1, ct_file_name = "")

# Checking number of DEG --------------------------------------------------

# NOTE: The q-values are the BH adjusted p-values

print_number_of_deg_wilcoxon <- function(ct_file_name, fdr){
  bm10_ct_deg_wilcoxon <- read.csv(paste(path_results, "BM10", ct_file_name, "_DEG_WilcoxonRankSumTest.csv", sep = ""),
                                   stringsAsFactors = FALSE)
  bm22_ct_deg_wilcoxon <- read.csv(paste(path_results, "BM22", ct_file_name, "_DEG_WilcoxonRankSumTest.csv", sep = ""),
                                   stringsAsFactors = FALSE)
  bm36_ct_deg_wilcoxon <- read.csv(paste(path_results, "BM36", ct_file_name, "_DEG_WilcoxonRankSumTest.csv", sep = ""),
                                   stringsAsFactors = FALSE)
  bm44_ct_deg_wilcoxon <- read.csv(paste(path_results, "BM44", ct_file_name, "_DEG_WilcoxonRankSumTest.csv", sep = ""),
                                   stringsAsFactors = FALSE)
  
  cat("BM10 ", length(which(bm10_ct_deg_wilcoxon$q.value.up.regulated.in.Control <= fdr)) +
      length(which(bm10_ct_deg_wilcoxon$q.value.up.regulated.in.Definite.AD <= fdr)), "\n")
  cat("BM22 ", length(which(bm22_ct_deg_wilcoxon$q.value.up.regulated.in.Control <= fdr)) +
      length(which(bm22_ct_deg_wilcoxon$q.value.up.regulated.in.Definite.AD <= fdr)), "\n")
  cat("BM36 ", length(which(bm36_ct_deg_wilcoxon$q.value.up.regulated.in.Control <= fdr)) +
      length(which(bm36_ct_deg_wilcoxon$q.value.up.regulated.in.Definite.AD <= fdr)), "\n")
  cat("BM44 ", length(which(bm44_ct_deg_wilcoxon$q.value.up.regulated.in.Control <= fdr)) +
      length(which(bm44_ct_deg_wilcoxon$q.value.up.regulated.in.Definite.AD <= fdr)), "\n")
}

print_number_of_deg_wilcoxon(fdr = 0.1, ct_file_name = "cellTypeCorrected")
print_number_of_deg_wilcoxon(fdr = 0.2, ct_file_name = "cellTypeCorrected")
print_number_of_deg_wilcoxon(fdr = 0.05, ct_file_name = "cellTypeCorrected")
print_number_of_deg_wilcoxon(fdr = 0.05, ct_file_name = "")

# Loading wang4 -----------------------------------------------------------

fdr <- 0.05

wang4 <- read.csv("/data/kailash/DDP/results/1-s2.0-S0896627320308618-mmc3.csv",
                  stringsAsFactors = FALSE, skip = 1, header = TRUE)
colnames(wang4) <- wang4[1,]
wang4 <- wang4[-1,]
head(wang4)
wang4_cerad <- wang4[which(wang4$Trait == "CERAD"),]

wang4_bm10 <- wang4_cerad[which(wang4_cerad$adj.P.Val <= fdr & wang4_cerad$Region == "BM10"),]
wang4_bm22 <- wang4_cerad[which(wang4_cerad$adj.P.Val <= fdr & wang4_cerad$Region == "BM22"),]
wang4_bm36 <- wang4_cerad[which(wang4_cerad$adj.P.Val <= fdr & wang4_cerad$Region == "BM36"),]
wang4_bm44 <- wang4_cerad[which(wang4_cerad$adj.P.Val <= fdr & wang4_cerad$Region == "BM44"),]

# Venn Diagram ------------------------------------------------------------

draw_venn_diagram <- function(x1, x2, x3, x4, brainRegion, ct_correction){
  venn.diagram(
    x = list(c(x1$X[which(x1$q.value.up.regulated.in.Control <= fdr)],
               x1$X[which(x1$q.value.up.regulated.in.Definite.AD <= fdr)]),
             x2$Gene[which(x2$adj.P.Val <= fdr)],
             x3,
             x4
             ),
    category.names = c("Wilcoxon Rank Sum Test" , "Limma on normalized data"),
    filename = paste(brainRegion, ct_correction, "_DEG_Overlap.png", sep = ""),
    output = TRUE,
  
    fill = c(alpha("#1B9E77",0.3), alpha("#D95F02",0.3)),
    fontfamily = "sans",
    cat.default.pos = "text",
    cat.fontfamily = "sans",
    cat.dist = c(-0.05, -0.05),
    cat.col = c("#1B9E77", "#D95F02"),
  
    cex = 0.9,
    main = paste(brainRegion, ct_correction),
    main.fontfamily = "sans",
    fontface = "bold"
  )
  }

draw_venn_diagram(bm10_ct_deg_wilcoxon, bm10_deg_limma, bm10_deg_limma_voom, wang4_bm10,
                  "BM10", "_ct_corrected_")

draw_venn_diagram(bm22_ct_deg_wilcoxon, bm22_deg_limma, bm22_deg_limma_voom, wang4_bm22,
                  "BM22", "_ct_corrected_")

draw_venn_diagram(bm36_ct_deg_wilcoxon, bm36_deg_limma, bm36_deg_limma_voom, wang4_bm36,
                  "BM36", "_ct_corrected_")

draw_venn_diagram(bm44_ct_deg_wilcoxon, bm44_deg_limma, bm44_deg_limma_voom, wang4_bm44, 
                  "BM44", "_ct_corrected_")

# Finding number of marker genes within DEGs ------------------------------

get_overlap_deg_mg_wilcoxon <- function(x1){
  degs <- c(x1$gene[which(x1$q.value.up.regulated.in.Control <= fdr)],
            x1$gene[which(x1$q.value.up.regulated.in.Definite.AD <= fdr)])
  
  l1 <- length(degs)
  l2 <- length(which(degs %in% markers_df_human_brain$markers))
  cat(l1, l2)
}

get_overlap_deg_mg_wilcoxon(bm10_deg_wilcoxon)
get_overlap_deg_mg_wilcoxon(bm22_deg_wilcoxon)
get_overlap_deg_mg_wilcoxon(bm36_deg_wilcoxon)
get_overlap_deg_mg_wilcoxon(bm44_deg_wilcoxon)

get_overlap_deg_mg_wilcoxon(bm10_ct_deg_wilcoxon)
get_overlap_deg_mg_wilcoxon(bm22_ct_deg_wilcoxon)
get_overlap_deg_mg_wilcoxon(bm36_ct_deg_wilcoxon)
get_overlap_deg_mg_wilcoxon(bm44_ct_deg_wilcoxon)

get_overlap_deg_mg_limma <- function(x2){
  degs <- c(x2$GeneSymbol[which(x2$adj.P.Val <= fdr)])
  
  l1 <- length(degs)
  l2 <- length(which(degs %in% markers_df_human_brain$markers))
  cat(l1, l2)
}

get_overlap_deg_mg_limma(bm10_deg_limma)
get_overlap_deg_mg_limma(bm22_deg_limma)
get_overlap_deg_mg_limma(bm36_deg_limma)
get_overlap_deg_mg_limma(bm44_deg_limma)

get_overlap_deg_mg_limma(bm10_ct_deg_limma)
get_overlap_deg_mg_limma(bm22_ct_deg_limma)
get_overlap_deg_mg_limma(bm36_ct_deg_limma)
get_overlap_deg_mg_limma(bm44_ct_deg_limma)

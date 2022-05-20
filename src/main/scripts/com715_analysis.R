# Analyzing Community 715

rm(list = ls())

# Loading libraries -------------------------------------------------------

.libPaths(c("/data/kailash/DDP/libraries/", .libPaths()))
source("/data/kailash/DDP/src/main/scripts/libraries.R")
source("/data/kailash/DDP/src/main/scripts/utils.R")

# Loading datasets --------------------------------------------------------

path_results <- "/data/kailash/DDP/results/"
path_data <- "/data/kailash/DDP/data/"
path_script <- "/data/kailash/DDP/src/main/scripts/"
path_geneExp <- "/data/srivatsan/Original_Expression_Files/Original_"

source(file = paste(path_script,  "utils.R", sep = ""))

bm10 <- readBMfile("BM10")
bm36 <- readBMfile("BM36")
mapping_file <- read.csv(paste(path_data, 
                               "Mapped_h37_ens_to_gene.csv", 
                               sep = ""), 
                         stringsAsFactors = FALSE)

row.names(bm10) <- find_individual_id("BM10")
row.names(bm36) <- find_individual_id("BM36")
colnames(bm10) <- colnames(bm36) <- ensid_sampleid[,1]

com715_genes <- data.frame(read_excel(paste(path_data, "com715.xlsx", sep = ""), sheet = 1))
com715_go <- data.frame(read_excel(paste(path_data, "com715.xlsx", sep = ""), sheet = 2))
com715_annot <- data.frame(read_excel(paste(path_data, "com715.xlsx", sep = ""), sheet = 3))
com715_genes$BM10[which(com715_genes$BM10 == "SEPTIN3")] <- "SEPT3"
com715_annot$BM10[which(com715_annot$BM10 == "SEPTIN3")] <- "SEPT3"

com715_bm10_genes <- ensid_sampleid$ENS_ID[match(com715_genes$BM10[which(!is.na(com715_genes$BM10))],
                                                 genes_mapped)]
com715_bm36_genes <- ensid_sampleid$ENS_ID[match(com715_genes$BM36[which(!is.na(com715_genes$BM36))],
                                                 genes_mapped)]
com715_bm10_annot <- com715_annot$BM10.annotation[which(!is.na(com715_annot$BM10.annotation))]
com715_bm36_annot <- com715_annot$BM36.annotation[which(!is.na(com715_annot$BM36.annotation))]

com715_bm10_annotgenes <- na.omit(mapping_file$gene_name[match(com715_bm10_genes[which(com715_bm10_annot == "Annotated")], 
                                                       mapping_file$gene_id)])
com715_bm10_NAnnotgenes <- na.omit(mapping_file$gene_name[match(com715_bm10_genes[which(com715_bm10_annot == "Not_Annotated")], 
                                                        mapping_file$gene_id)])
com715_bm36_annotgenes <- na.omit(mapping_file$gene_name[match(com715_bm36_genes[which(com715_bm36_annot == "Annotated")], 
                                                       mapping_file$gene_id)])
com715_bm36_NAnnotgenes <- na.omit(mapping_file$gene_name[match(com715_bm36_genes[which(com715_bm36_annot == "Not_Annotated")], 
                                                        mapping_file$gene_id)])

bm10_com715 <- bm10[,which(colnames(bm10) %in% com715_bm10_genes)]
bm36_com715 <- bm36[,which(colnames(bm36) %in% com715_bm36_genes)]

stopifnot(dim(bm10_com715)[2] <= length(com715_bm10_genes))
stopifnot(dim(bm36_com715)[2] <= length(com715_bm36_genes))

dim(bm10_com715)
dim(bm36_com715)

# Finding AD, CTL common samples --------------------------------------------------

bm10_sample_gp <- find_sample_group("BM10")
bm36_sample_gp <- find_sample_group("BM36")

ctl_ids_bm10 <- row.names(bm10_com715)[which(bm10_sample_gp == "Control")]
ad_ids_bm10 <- row.names(bm10_com715)[which(bm10_sample_gp == "Definite AD")]

ctl_ids_bm36 <- row.names(bm36_com715)[which(bm36_sample_gp == "Control")]
ad_ids_bm36 <- row.names(bm36_com715)[which(bm36_sample_gp == "Definite AD")]

ad_ids_BR1_2 <- intersect(ad_ids_bm10, ad_ids_bm36)
ctl_ids_BR1_2 <- intersect(ctl_ids_bm10, ctl_ids_bm36)

ad_bm10 <- bm10_com715[ad_ids_BR1_2,]
ad_bm36 <- bm36_com715[ad_ids_BR1_2,]
ctl_bm10 <- bm10_com715[ctl_ids_BR1_2,]
ctl_bm36 <- bm36_com715[ctl_ids_BR1_2,]

dim(ctl_bm10)
dim(ad_bm10)
dim(ctl_bm36)
dim(ad_bm36)

AD1 <- ad_bm10
AD2 <- ad_bm36
CTL1 <- ctl_bm10
CTL2 <- ctl_bm36

n.ad <- dim(AD1)[1]
n.ctl <- dim(CTL1)[1]

ad_bm10_genes <- colnames(ad_bm10)
ad_bm36_genes <- colnames(ad_bm36)
ctl_bm10_genes <- colnames(ctl_bm10)
ctl_bm36_genes <- colnames(ctl_bm36)

# Creating all possible edges file ----------------------------------------

fdr_val <- 1

cor_AD12 <- cor(AD1, AD2, method = "spearman")
cor_CTL12 <- cor(CTL1, CTL2, method = "spearman")

p_AD12 <- pvalfromtcor(cor_AD12, n.ad)
p_CTL12 <- pvalfromtcor(cor_CTL12, n.ctl)

# Benjamini-Hochberg pval adjustment
A_AD12 <- matrix(p.adjust(p_AD12, method = "BH"), 
                 nrow = dim(p_AD12)[1], ncol = dim(p_AD12)[2])
A_CTL12 <- matrix(p.adjust(p_CTL12, method = "BH"), 
                  nrow = dim(p_CTL12)[1], ncol = dim(p_CTL12)[2])

# cbind(x,y) - Adds y after columns of x

indices <- which(A_AD12 <= fdr_val)
cor_AD12_edges <- cbind(ad_bm10_genes[row(A_AD12)[indices]], 
                        ad_bm36_genes[col(A_AD12)[indices]],
                        as.double(cor_AD12[indices]),
                        as.double(p_AD12[indices]),
                        as.double(A_AD12[indices]))
colnames(cor_AD12_edges) <- c("g1-BR1", "g2-BR2", "r_ad", "p-val", "p-adj")

rm(indices)

indices <- which(A_CTL12 <= fdr_val)
cor_CTL12_edges <- cbind(ctl_bm10_genes[row(A_CTL12)[indices]],
                         ctl_bm36_genes[col(A_CTL12)[indices]],
                         as.double(cor_CTL12[indices]),
                         as.double(p_CTL12[indices]),
                         as.double(A_CTL12[indices]))
colnames(cor_CTL12_edges) <- c("g1-BR1", "g2-BR2", "r_ctl", "p-val", "p-adj")
rm(indices)

# hist(as.double(cor_AD12_edges[,2]))
# hist(as.double(cor_CTL12_edges[,2]))
cor_AD12_edges <- data.frame(cor_AD12_edges, stringsAsFactors = FALSE)
cor_CTL12_edges <- data.frame(cor_CTL12_edges, stringsAsFactors = FALSE)

cor_AD12_edges$combined <- paste(cor_AD12_edges$g1.BR1,cor_AD12_edges$g2.BR2)
cor_CTL12_edges$combined <- paste(cor_CTL12_edges$g1.BR1,cor_CTL12_edges$g2.BR2)

write.csv(cor_AD12_edges, paste(path_results, "BM10-BM36_com715_ADCorrelatedEdges.csv", sep = ""))
write.csv(cor_CTL12_edges, paste(path_results, "BM10-BM36_com715_CTLCorrelatedEdges.csv", sep = ""))

# Find union of cor_ad_34_edges, cor_ctl_34_edges
# NOTE: A particular edge in union12 need not have a signficant edge in both AD and CTl
# but they definitely will have significance in either AD or CTL.
union12 <- union(cor_AD12_edges$combined, cor_CTL12_edges$combined)
stopifnot(union12 != '')

# For a given edge a-b
edges <- strsplit(union12, " ")
mat  <- matrix(unlist(edges), ncol = 2, byrow = TRUE)
df   <- as.data.frame(mat, stringsAsFactors = FALSE)
colnames(df) <- c("a", "b")

edges <- df

edges$combined <- paste(edges$a, edges$b)
union_ad <- cor_AD12_edges[match(edges$combined, cor_AD12_edges$combined),c("r_ad", "p.val", "p.adj")]
colnames(union_ad) <- c("r_ad", "p.val.ad", "p.adj.ad")
union_ctl <- cor_CTL12_edges[match(edges$combined, cor_CTL12_edges$combined), c("r_ctl", "p.val", "p.adj")]
colnames(union_ctl) <- c("r_ctl", "p.val.ctl", "p.adj.ctl")

edges <- cbind(select(edges, c(a,b)), union_ad, union_ctl)

write.csv(edges, paste(path_results,"UnionCorrelatedEdges_BM10_BM36_com715", ".csv", sep = ""))

rm(df)

# Finding z scores --------------------------------------------------------

# From Brain region 1
t1 <- ad_bm10[,edges[,1]]
t2 <- ctl_bm10[,edges[,1]]

# From Brain region 2
t3 <- ad_bm36[,edges[,2]]
t4 <- ctl_bm36[,edges[,2]]

# Zip fastener(x,y)
# col1 is col1 of x
# col2 is col1 of y
# and so on

res_AD12 <- zipFastener(t1,t3)
res_CTL12 <- zipFastener(t2,t4)

res_AD12 <- as.data.frame(res_AD12)
res_CTL12 <- as.data.frame(res_CTL12)

rm(t1,t2,t3,t4)

# r.test DC edges ---------------------------------------------------------

# Initializing list of p-values vector
p <- rep(NA,length(colnames(res_AD12))/2 + 1)
z <- rep(NA, length(p))
cor_CTL <- rep(NA, length(p))
cor_AD <- rep(NA, length(p))

i <- 1
j <- 1

while(i <= length(colnames(res_AD12))){
  cat("\n",i,"out of", length(colnames(res_AD12)),"times \n") # Printing out iteration number
  temp1 <- res_CTL12[,c(i,i+1)]
  temp2 <- res_AD12[,c(i,i+1)]
  colnames(temp1) <- colnames(temp2) <- c('a','b')
  
  cor1 <- cor(temp1[,'a'], temp1[,'b'], method = "spearman")
  cor2 <- cor(temp2[,'a'], temp2[,'b'], method = "spearman")
  
  cor_CTL[j] <- cor1
  cor_AD[j] <- cor2
  
  rtest_vals <- r.test(n = n.ctl, n2 = n.ad, r12 = cor1, r34 = cor2)
  
  p[j] <- rtest_vals$p
  z[j] <- rtest_vals$z * sign(cor2 - cor1)
  
  # To keep info regarding LOC and GOC
  # LOC - Loss of correlation wrt Control, GOC - Gain of correlation wrt Control
  
  # i + 2 because we are taking the bm36 and bm44 adjacent columns gene exp values together
  j <- j + 1
  i <- i + 2
}

rm(i,j,temp1,temp2, rtest_vals)

# Removing NAs in the end of the vector
indices <- which(!is.na(p))
stopifnot(length(is.na(p)) > 2)
p <- p[indices]
z <- z[indices]
cor_CTL <- cor_CTL[indices]
cor_AD <- cor_AD[indices]
# NOTE the following file contains p values before adjustment
save(p, file = paste(path_results,"p_vals_without_adjustment_com715_BM10_BM36",".rda",sep=""))

p12 <- p
p_adjusted <- p.adjust(p12, method = 'BH')

rm(p, indices)

# Print DC edges
dcres <- cbind(edges, p_adjusted, z, cor_CTL, cor_AD)
rm(edges)
r <- which(p_adjusted <= fdr_val)
dc_rtest <- dcres[r,]

dc_rtest$BM10_annot <- com715_bm10_annot[match(dc_rtest$a, com715_bm10_genes)]
dc_rtest$BM36_annot <- com715_bm36_annot[match(dc_rtest$b, com715_bm36_genes)]

dc_rtest$g1 <- mapping_file$gene_name[match(dc_rtest$a, mapping_file$gene_id)]
dc_rtest$g2 <- mapping_file$gene_name[match(dc_rtest$b, mapping_file$gene_id)]

dc_rtest$BM10_go <- com715_go$BM10.annotation[match(dc_rtest$g1, com715_go$BM10)]
dc_rtest$BM36_go <- com715_go$BM36.annotation[match(dc_rtest$g2, com715_go$BM36)]

write.csv(dc_rtest, paste(path_results,"DCEdges_com715_BM10_BM36",".csv",sep=""))

# Making the Z-matrix and pheatmap ----------------------------------------

z_df <- dc_rtest[, c("a", "b", "z")]

# Reference for the following line of code: https://stackoverflow.com/a/54471567/8281802
z_matrix <- with(z_df, tapply(z, list(a, b), sum))
row.names(z_matrix) <- mapping_file$gene_name[match(row.names(z_matrix), mapping_file$gene_id)]
colnames(z_matrix) <- mapping_file$gene_name[match(colnames(z_matrix), mapping_file$gene_id)]
dim(z_matrix)
stopifnot(dim(dc_rtest)[1] == dim(z_matrix)[1] * dim(z_matrix)[2])

annotation_y <- data.frame(com715_annot$BM10.annotation[match(row.names(z_matrix), com715_annot$BM10)])
row.names(annotation_y) <- row.names(z_matrix)
colnames(annotation_y) <- "Annotation"
dim(annotation_y)

annotation_x <- data.frame(com715_annot$BM36.annotation[match(colnames(z_matrix), com715_annot$BM36)])
row.names(annotation_x) <- colnames(z_matrix)
colnames(annotation_x) <- "Annotation"
dim(annotation_x)

# creat colours for each group
# newCols <- colorRampPalette(grDevices::rainbow(length(unique(annotation_x[,1]))))
# annoCol <- newCols(length(unique(annotation_x[,1])))
annoCol <- c("orchid", "coral")
names(annoCol) <- unique(annotation_x[,1])
annoCol <- list(Annotation = annoCol)

# Var1        <- c("navy", "darkgreen")
# names(Var1) <- c("Annotated", "Not_Annotated")
# anno_colors <- list(Var1 = Var1)

p1 <- pheatmap(z_matrix,
              fontsize = 5,
              method = "spearman",
              cluster_rows = T,
              clustering_method = "complete",
              annotation_col =  annotation_x,
              annotation_row = annotation_y,
              color = colorRampPalette(c("orchid", "white", "coral"))(50),
              # color = viridis(50),
              treeheight_row = 40,
              treeheight_col = 40,
              annotation_colors =  annoCol
)

ggsave(plot = p1,
       units = "in",
       height = 10, # 6 - 8, try different sizes
       width = 15,
       filename = "/data/kailash/DDP/results/com715_pheatmap_all.pdf"
)

annotation_x_p2 <- data.frame(annotation_x[com715_bm36_annotgenes,],
                              row.names = com715_bm36_annotgenes,
                              stringsAsFactors = FALSE)
colnames(annotation_x_p2) <- "Annotation"

annotation_y_p2 <- data.frame(annotation_y[com715_bm10_annotgenes,],
                              row.names = com715_bm10_annotgenes,
                              stringsAsFactors = FALSE)
colnames(annotation_y_p2) <-  "Annotation"

# A - A
p2 <- pheatmap(z_matrix[com715_bm10_annotgenes, com715_bm36_annotgenes],
               fontsize = 8,
               method = "spearman",
               cluster_rows = T,
               clustering_method = "complete",
               color = colorRampPalette(c("orchid", "white", "coral"))(50),
               # color = viridis(50),
               treeheight_row = 40,
               treeheight_col = 40,
               annotation_col = annotation_x_p2,
               annotation_row = annotation_y_p2,
               annotation_colors =  annoCol
)

ggsave(plot = p2,
       units = "in",
       height = 10, # 6 - 8, try different sizes
       width = 15,
       filename = "/data/kailash/DDP/results/com715_pheatmap_a_a.pdf"
)

annotation_x_p3 <- data.frame(annotation_x[com715_bm36_NAnnotgenes,],
                              row.names = com715_bm36_NAnnotgenes,
                              stringsAsFactors = FALSE)
colnames(annotation_x_p3) <- "Annotation"

annotation_y_p3 <- data.frame(annotation_y[com715_bm10_annotgenes,],
                              row.names = com715_bm10_annotgenes,
                              stringsAsFactors = FALSE)
colnames(annotation_y_p3) <-  "Annotation"

# A-NA
p3 <- pheatmap(z_matrix[com715_bm10_annotgenes, com715_bm36_NAnnotgenes],
               fontsize = 8,
               method = "spearman",
               cluster_rows = T,
               clustering_method = "complete",
               color = colorRampPalette(c("orchid", "white", "coral"))(50),
               # color = viridis(50),
               treeheight_row = 40,
               treeheight_col = 40,
               annotation_col =  annotation_x_p3,
               annotation_row = annotation_y_p3,
               annotation_colors =  annoCol
)

ggsave(plot = p3,
       units = "in",
       height = 10, # 6 - 8, try different sizes
       width = 15,
       filename = "/data/kailash/DDP/results/com715_pheatmap_a_na.pdf"
)

annotation_x_p4 <- data.frame(annotation_x[com715_bm36_annotgenes,],
                              row.names = com715_bm36_annotgenes,
                              stringsAsFactors = FALSE)
colnames(annotation_x_p4) <- "Annotation"

annotation_y_p4 <- data.frame(annotation_y[com715_bm10_NAnnotgenes,],
                              row.names = com715_bm10_NAnnotgenes,
                              stringsAsFactors = FALSE)
colnames(annotation_y_p4) <-  "Annotation"

# NA-A
p4 <- pheatmap(z_matrix[com715_bm10_NAnnotgenes, com715_bm36_annotgenes],
               fontsize = 8,
               method = "spearman",
               cluster_rows = T,
               clustering_method = "complete",
               color = colorRampPalette(c("orchid", "white", "coral"))(50),
               # color = viridis(50),
               treeheight_row = 40,
               treeheight_col = 40,
               annotation_col =  annotation_x_p4,
               annotation_row = annotation_y_p4,
               annotation_colors =  annoCol
)

ggsave(plot = p4,
       units = "in",
       height = 10, # 6 - 8, try different sizes
       width = 15,
       filename = "/data/kailash/DDP/results/com715_pheatmap_na_a.pdf"
)

annotation_x_p5 <- data.frame(annotation_x[com715_bm36_NAnnotgenes,],
                              row.names = com715_bm36_NAnnotgenes,
                              stringsAsFactors = FALSE)
colnames(annotation_x_p5) <- "Annotation"

annotation_y_p5 <- data.frame(annotation_y[com715_bm10_NAnnotgenes,],
                              row.names = com715_bm10_NAnnotgenes,
                              stringsAsFactors = FALSE)
colnames(annotation_y_p5) <-  "Annotation"

# NA-NA
p5 <- pheatmap(z_matrix[com715_bm10_NAnnotgenes, com715_bm36_NAnnotgenes],
               fontsize = 8,
               method = "spearman",
               cluster_rows = T,
               clustering_method = "complete",
               color = colorRampPalette(c("orchid", "white", "coral"))(50),
               # color = viridis(50),
               treeheight_row = 40,
               treeheight_col = 40,
               annotation_col = annotation_x_p5,
               annotation_row = annotation_y_p5,
               annotation_colors =  annoCol
)

ggsave(plot = p5,
       units = "in",
       height = 10, # 6 - 8, try different sizes
       width = 15,
       filename = "/data/kailash/DDP/results/com715_pheatmap_na_na.pdf"
)

# Pheatmap of annotated-annotated -----------------------------------------

z_df <- z_matrix[unique(com715_go$BM10), unique(com715_go$BM36)]
row.names(z_df) <- unique(com715_go$BM10)
colnames(z_df) <- unique(com715_go$BM36)
dim(z_df)

annotation_x <- data.frame(dc_rtest$BM10_go[match(row.names(z_df), dc_rtest$g1)],
                           stringsAsFactors = FALSE)
row.names(annotation_x) <- row.names(z_df)
colnames(annotation_x) <- "GO_BP"

annotation_y <- data.frame(dc_rtest$BM36_go[match(colnames(z_df), dc_rtest$g2)],
                           stringsAsFactors = FALSE)
row.names(annotation_y) <- colnames(z_df)
colnames(annotation_y) <- "GO_BP"

annot_df_color <- data.frame(rbind(annotation_x, annotation_y))
all_GO_BP <- unique(annot_df_color[,1])

# creat colours for each group
# newCols <- colorRampPalette(colorRamps::primary.colors(length(unique(all_GO_BP))))
# annoCol <- newCols(length(unique(all_GO_BP)))
# annoCol <- magma(length(unique(all_GO_BP)))
annoCol <- brewer.pal(n = length(unique(all_GO_BP)), name = "Paired")
names(annoCol) <- unique(all_GO_BP)
annoCol <- list(GO_BP = annoCol)

p6 <- pheatmap(z_df,
              fontsize = 10,
              method = "spearman",
              cluster_rows = T,
              clustering_method = "complete",
              annotation_row = annotation_x,
              annotation_col = annotation_y,
              color = colorRampPalette(c("orchid", "white", "coral"))(50),
              # color = viridis(50),
              treeheight_row = 40,
              treeheight_col = 40,
              annotation_colors = annoCol
)

ggsave(plot = p6,
       units = "in",
       height = 10, # 6 - 8, try different sizes
       width = 15,
       filename = "/data/kailash/DDP/results/com715_pheatmap_annotated_annotated_edges.pdf"
       )

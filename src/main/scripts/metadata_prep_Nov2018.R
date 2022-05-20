rm(list = ls())

# Prepping the metadata file ----------------------------------------------

path_data <- "/data/kailash/DDP/data/"

rna_metadata <- read.csv(paste(path_data, "MSBB_RNAseq_covariates_November2018Update.csv", sep = ""),stringsAsFactors = FALSE)
rna_metadata <- rna_metadata[which(rna_metadata$Action != "Exclude") , ]
rna_metadata <- rna_metadata[which(rna_metadata$fileType == "bam") , ]
duplicated_barcodes <- rna_metadata$barcode[duplicated(rna_metadata$barcode)]
unique_barcodes <- unique(rna_metadata$barcode)

fileName_remove <- c()
i <- 1
while(i <= length(duplicated_barcodes)){
  dat <- rna_metadata[rna_metadata$barcode == duplicated_barcodes[i],]
  minFile <- dat$fileName[dat$Mapped == min(dat$Mapped)]
  fileName_remove <- c(fileName_remove, minFile)
  i = i + 1
}

rna_metadata <- rna_metadata[!rna_metadata$fileName %in% fileName_remove,]
duplicated_barcodes <- rna_metadata$barcode[duplicated(rna_metadata$barcode)]

cat(duplicated_barcodes)
i <- 1
fileName_remove <- c()
while(i <= length(duplicated_barcodes)){
  dat <- rna_metadata[rna_metadata$barcode == duplicated_barcodes[i],]
  minFile <- dat$fileName[dat$Mapped == min(dat$Mapped)]
  fileName_remove <- c(fileName_remove, minFile)
  i = i + 1
}
rna_metadata <- rna_metadata[!rna_metadata$fileName %in% fileName_remove,]

duplicated_barcodes <- rna_metadata$barcode[duplicated(rna_metadata$barcode)]

# Adding clinical variables to the metadata file --------------------------

c_data <- read.csv(paste(path_data, "MSBB_clinical.csv", sep = ""), stringsAsFactors = FALSE)
m <- read.csv(paste(path_data, "rna_metadata_onlyBAM_noExclude_noDup_Nov2018.csv", sep = ""),
              stringsAsFactors = FALSE)
combined_metadata <- cbind(m, subset(c_data[match(m$individualIdentifier, c_data$individualIdentifier),], select = -individualIdentifier))

dim(c_data)
dim(m)
dim(combined_metadata)

library(ggplot2)
library(ggpubr)
p1 <- ggplot(combined_metadata, aes(bbscore)) + geom_histogram()
p2 <- ggplot(c_data, aes(bbscore)) + geom_histogram()
ggarrange(p1, p2, nrow = 2)

# Writing the prepped metadata file ---------------------------------------

write.csv(combined_metadata, paste(path_data, "rna_metadata_onlyBAM_noExclude_noDup_Nov2018.csv", sep = ""))
rna_metadata <- read.csv("rna_metadata_onlyBAM_noExclude_noDup_Nov2018.csv", stringsAsFactors = FALSE)
dim(rna_metadata)
colnames(rna_metadata)
rna_metadata$barcode
rna_metadata$bbscore

# Prepping the GenCode mapping --------------------------------------------

setwd("/data/kailash/DDP/data/")

# hg37 --------------------------------------------------------------------

h37 <- read.table("./gencode.v19.annotation.gff3", sep = "\t",
                  stringsAsFactors = FALSE)
dim(h37) # 26,15,566 X 9
unique(h37$V3)
h37 <- h37[which(h37$V3 == "gene"),]
dim(h37) # 57,820 X 9
unique(h37$V1) # chromsome 1 - 22, chr X, chr Y, chr M
unique(h37$V2) # HAVANA, ENSEMBL

v9 <- strsplit(as.character(h37$V9),split = ";")
v9 <- do.call(rbind, v9)
colnames(v9) <- c("ID1", "gene_id", "transcript_id", "gene_type", "gene_status", "gene_name", "transcript_type",
                  "transcript_status", "transcript_name", "level", "havana_gene", "ID2")
v9 <- data.frame(v9, stringsAsFactors = FALSE)
dim(v9) # 57,820 X 12
v9 <- v9[, c(2, 4, 5, 6)]

df1 <- strsplit(as.character(v9[,1]), split = "=")
df1 <- do.call(rbind, df1)[,2]
df1 <- strsplit(as.character(df1), split = "[.]")
df1 <- do.call(rbind, df1)[,1]

split_by_equal_to <- function(df, i){
  x <- strsplit(as.character(df[,i]), split = "=")
  return(do.call(rbind, x)[,2])
}

df2 <- split_by_equal_to(v9, 2)
df3 <- split_by_equal_to(v9, 3)
df4 <- split_by_equal_to(v9, 4)

rm(v9)
h37_map <- data.frame(cbind(h37$V1, df1, df2, df3, df4), stringsAsFactors = FALSE)
colnames(h37_map) <- c("chr", "gene_id", "gene_type", "gene_status", "gene_name")

head(h37_map)
dim(h37_map) # 57,820 X 5

# Checking with our ENS_IDS -----------------------------------------------

ens_ids <- read.csv("./ens_ids.csv", stringsAsFactors = FALSE)
ens_ids <- ens_ids$x
head(ens_ids)
# Checking if all ENS_IDS are unique
stopifnot(length(unique(ens_ids)) == length(ens_ids)) # 23,201

length(ens_ids) - length(which(ens_ids %in% h37_map$gene_id)) # 304 ENS_IDS are not present in h37 mapping
ens_ids_no_map <- ens_ids[-which(ens_ids %in% h37_map$gene_id)]
write.csv(ens_ids_no_map, "ens_ids_no_map.csv")


# Checking if there exists mapping for ens_ids_no_map in hg38 --------------------------------

h38 <- read.table("./gencode.v32.annotation.gff3", stringsAsFactors = FALSE, sep = "\t")
h38 <- h38[which(h38$V3 == "gene"),]
dim(h38) # 60,609 X 9

v9 <- strsplit(as.character(h38$V9), split = ";")
v9 <- do.call(rbind, v9)
v9 <- v9[,c(1,2,3,4)]

df1 <- strsplit(as.character(v9[,1]), split = "=")
df1 <- do.call(rbind, df1)[,2]
df1 <- strsplit(as.character(df1), split = "[.]")
df1 <- do.call(rbind, df1)[,1]

split_by_equal_to <- function(df, i){
  x <- strsplit(as.character(df[,i]), split = "=")
  return(do.call(rbind, x)[,2])
}

df2 <- split_by_equal_to(v9, 2)
df3 <- split_by_equal_to(v9, 3)
df4 <- split_by_equal_to(v9, 4)

rm(v9)
h38_map <- data.frame(cbind(h38$V1, df1, df3, df4), stringsAsFactors = FALSE)
colnames(h38_map) <- c("chr", "gene_id", "gene_type", "gene_name")

head(h38_map)
dim(h38_map) # 60,609 X 4

which(h38_map$gene_id %in% ens_ids_no_map) # 0
# SO these 304 ENS_IDS are not present in hg38 mapping as well.

# Checking if those 304 ENS_IDS are present in my mapping -----------------

my_mapping <- read.csv("biomart_biotools_combined_mapping_file.csv",
                       stringsAsFactors = FALSE)
which(my_mapping$Gene.stable.ID %in% ens_ids_no_map) # 304

stopifnot(unique(my_mapping[match(ens_ids_no_map, my_mapping$Gene.stable.ID),]$HGNC.symbol) == "")

# Working with h37 mapping ------------------------------------------------

h37_ens <- h37_map[which(h37_map$gene_id %in% ens_ids),] # Hay in needle analogy
dim(h37_ens) # 22,909 X 5

duplicates1 <- h37_ens[which(duplicated(h37_ens$gene_id, fromLast = TRUE)),]
duplicates2 <- h37_ens[which(duplicated(h37_ens$gene_id, fromLast = FALSE)),]

# They are exact duplicates, so I need not worry, but are from different chromosomes
stopifnot(duplicates1[,2:length(colnames(duplicates1))] == duplicates2[,2:length(colnames(duplicates2))])

# NOTE: match returns the first match, so there is a possibility of 
# an ENSEMBL ID mapping to multiple HGNC, which needs to be identified and
# corrected for
# But it doesn't matter here if I just get the first match, because the duplicates
# have the same values
# Multiple ENS_ID pointing to the same HGNC is fine, but multiple HGNC pointing 
# to the same ENS_ID is not - we fixed that, now we are looking at duplicate HGNC

# Have to remove one set of the 12 duplicates, because
# this will increase unnecessarily the number of HGNC duplicates as well
h37_ens <- h37_ens[-which(duplicated(h37_ens$gene_id, fromLast = TRUE)),]
dim(h37_ens) # 22,897 X 5
stopifnot(length(which(ens_ids %in% h37_map$gene_id)) == dim(h37_ens)[1])
stopifnot(length(unique(h37_ens$gene_id)) == dim(h37_ens)[1])
write.csv(h37_ens, "Mapped_h37_ens_to_gene.csv")

# Finding number of duplicates HGNC symbols

length(unique(h37_ens$gene_name))

duplicates3 <- h37_ens[which(duplicated(h37_ens$gene_name, fromLast = TRUE)),]
duplicates3 <- duplicates3[order(row.names(duplicates3)), ]
duplicates4 <- h37_ens[which(duplicated(h37_ens$gene_name, fromLast = FALSE)),]
duplicates4 <- duplicates4[order(row.names(duplicates4)), ]

length(unique(duplicates3$gene_id))
length(unique(duplicates4$gene_id))
length(which(duplicates3$gene_status == "NOVEL"))
length(which(duplicates4$gene_status == "NOVEL"))
intersect(duplicates3, duplicates4)
dim(duplicates3) # 81 X 5
dim(duplicates4) # 81 X 5

write.csv(duplicates3, "duplicates3.csv")
write.csv(duplicates4, "duplicates4.csv")

intersect(row.names(duplicates3),row.names(duplicates4))
length(union(duplicates3$gene_name, duplicates4$gene_name)) # 27 genes are duplicated many times

# Check overlap with markers_df_human_brain -------------------------------

library(BRETIGEA)
overlap_genes_with_BRETIGEA <- as.character(intersect(markers_df_human_brain[,1], union(duplicates3$gene_name, duplicates4$gene_name)))
markers_df_human_brain[which(markers_df_human_brain$markers %in% overlap_genes_with_BRETIGEA),]
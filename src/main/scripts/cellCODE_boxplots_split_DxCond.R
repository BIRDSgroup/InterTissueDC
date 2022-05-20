rm(list = ls())

# Global variables --------------------------------------------------------

path_results <- "/data/kailash/DDP/results/"
path_data <- "/data/kailash/DDP/data/"
BR <- "BM44"

# Loading libraries ------------------------------------------------------

.libPaths(c("/data/kailash/DDP/libraries/", .libPaths()))
source("/data/kailash/DDP/src/main/scripts/libraries.R")
source("/data/kailash/DDP/src/main/scripts/utils.R")

# Loading data ------------------------------------------------------------

ctfs <- read.csv(paste(path_results, BR, "_CellCODE_80mg_4cellTypes.csv", sep = ""),
                 stringsAsFactors = FALSE)
row.names(ctfs) <- ctfs$X
ctfs <- subset(ctfs, select = -X)
mdata <- read.csv(paste(path_data, "rna_metadata_onlyBAM_noExclude_noDup_Nov2018.csv", sep = ""),
                  stringsAsFactors = FALSE)

barcodes <- sapply(strsplit(row.names(ctfs), "_"), function(x){x[3]})
sample_gp <- mdata$DxCondition[match(barcodes, mdata$barcode)]

nSamples <- length(sample_gp)

ctfs_boxplot <- do.call("rbind", list(data.frame(cbind(ctfs$ast, rep("Ast", nSamples))),
                                      data.frame(cbind(ctfs$mic, rep("Mic", nSamples))),
                                      data.frame(cbind(ctfs$neu, rep("Neu", nSamples))),
                                      data.frame(cbind(ctfs$oli, rep("Oli", nSamples)))))
ctfs_boxplot$sample_gp <- sample_gp
ctfs_boxplot <- data.frame(ctfs_boxplot, stringsAsFactors = FALSE)
ctfs_boxplot$X1 <- as.numeric(as.character(ctfs_boxplot$X1))

ggplot(ctfs_boxplot, aes(x = X2, fill = sample_gp, y = X1)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position = position_dodge(0.8), cex = 0.5) +
  ylab("CellCODE cell-type frequency estimates") +
  xlab("Cell-types") + theme_bw() +
  labs(title = BR) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        axis.text.x = element_text(face = "bold", colour = "black"),
        legend.title = element_text(face = "bold")) +
  stat_compare_means(method = "wilcox.test") +
  ylim(-0.5, 0.5)

compare_means(X1 ~ sample_gp, ctfs_boxplot, group.by = "X2",
              p.adjust.method = "BH")

ggsave(paste(BR, "Cell-type-frequencies boxplot split by sample group.png", sep = ""))


### Import packages
library(dplyr)
library(reshape2)
library(ggplot2)
library(readxl)
library(edgeR)
library(clusterProfiler)
library(org.Hs.eg.db)

### Import datasets
breast_raw_data <- read.table("../GSE131512_cancerTPM.txt", 
                              header = T, 
                              row.names = 1, 
                              stringsAsFactors = F)
normal_raw_data <- read.table("../GSE131512_normalTPM.txt", 
                              header = T, 
                              row.names = 1, 
                              stringsAsFactors = F)

total_raw_data <- bind_cols(breast_raw_data, normal_raw_data)
rownames(total_raw_data) <- rownames(breast_raw_data)

metadate <- read_xlsx("../GSE131512_metaData.xlsx")
metadate$Type <- c(rep("Cancer", 96), rep("Normal", 32))


### Initial Filter
##    Filter those have small means
filter_data <- total_raw_data[which(rowSums(total_raw_data) > 128 * 20), ]
##    Filter those have too many 0
filter_data <- filter_data[which(rowSums(filter_data == 0) < 128 * 0.2), ]

{
  edgeR_list <- DGEList(counts = filter_data, group = metadate$Type)
  edgeR_model <- estimateGLMCommonDisp(edgeR_list, verbose = TRUE)
  edgeR_fit <- glmFit(edgeR_model)
  edgeR_lrt <- glmLRT(edgeR_fit)
  topTags(edgeR_lrt)
  edgeR_selected_feature <- edgeR_lrt$table[which((edgeR_lrt$table$PValue < 1e-05 / nrow(edgeR_lrt$table) 
                                                  & ((edgeR_lrt$table$logFC > 2) | (edgeR_lrt$table$logFC < 0.5)))
                                                  & (edgeR_lrt$table$logCPM > 3)
  ),]
  diff_data <- total_raw_data[which(rownames(total_raw_data) %in% rownames(edgeR_selected_feature)),]
  #rownames(diff_data) <- rownames(filter_data)[as.numeric(rownames(diff_data))]
}

##    Convert ENSG into gene names
gene_ID_matched <- mapIds(org.Hs.eg.db, 
                          keys = rownames(diff_data),
                          column = "SYMBOL", 
                          keytype = "ENSEMBL")
##    Some not matched automaticall7, find them manually on Ensembl
gene_ID_matched[7] <- "MT-RNR2"
gene_ID_matched[8] <- "MTND4P12"
gene_ID_matched[10] <- "AL353644.10"

### Merge into a final dataset
rownames(diff_data) <- gene_ID_matched
processed_data <- as.data.frame(t(diff_data))
processed_data <- bind_cols(processed_data, metadate)

write.csv(processed_data, "Processed data.csv", 
          quote = F,
          row.names = F)

##  Plotting
###   Volcano plt
volvano_df <- as.data.frame(edgeR_lrt$table)
volvano_df$Significance <- "No"
volvano_df$Significance[which(rownames(volvano_df) %in% rownames(edgeR_selected_feature))] <- "Yes"

p1 <- ggplot(data = volvano_df, 
             aes(x = logFC, y = -1*log10(PValue), col = Significance)) +  
  geom_point(size = 1) + 
  scale_color_manual(values=c("black", "red")) +
  xlim(c(-4, 4)) +
  geom_vline(xintercept = c(-2, 2), lty = 4,col = "black", lwd = 0.5) +
  geom_hline(yintercept = -log10(1e-05 / nrow(edgeR_lrt$table)), lty = 4, col = "black", lwd = 0.5) +
  labs(x = "log2(fold change)", y = "-log10 (p-value)", title = "Volcano Plot of Differential Genes") + 
  theme_bw()
p1
pdf("../Plot/Data process/Volcano plot.pdf", 6, 5)
p1
dev.off()





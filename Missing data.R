##  Import libraries
library(ggplot2)
library(dplyr)

## Import dataset
all_data <- read.csv("Processed data.csv", header = T, stringsAsFactors = F)

## Missing data summary
gene_data <- all_data[,1:10]
missing_number <- colSums(gene_data == 0)
missing_number_df <- data.frame(gene = names(missing_number), 
                                count = missing_number, 
                                percentage = missing_number / nrow(all_data))
p1 <- ggplot(data = missing_number_df, aes(x = gene, y = count, fill = gene)) + 
  geom_col() +
  xlab("Gene ID") +
  ylab("Count of Missing Data") + 
  scale_fill_discrete(name = "Gene ID") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 7, angle = 45, vjust = 0.7))

pdf("Plot/Missing data/Count of missing data.pdf", 5, 4)
p1
dev.off()


##  Relation to sample types
missing_type <- data.frame()
for(i in 1:10){
  loc <- which(gene_data[,i] == 0)
  temp <- data.frame(Gene = colnames(gene_data)[i], 
                     Type = all_data$Type[loc],
                     Gender = all_data$Gender[loc],
                     Age = all_data$Age[loc]
                     )
  missing_type <- rbind(missing_type, temp)
}

p2 <- ggplot(data = missing_type, aes(x = Age, col = Type)) + 
  geom_density() + 
  ylab("Density of missing data") + 
  theme_bw()

pdf("Plot/Missing data/Missing on age.pdf", 5, 4)
p2
dev.off()
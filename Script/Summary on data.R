##  Import libraries
library(ggplot2)
library(dplyr)
library(pheatmap)
library(psych)

## Import dataset
all_data <- read.csv("../Processed data.csv", header = T, stringsAsFactors = F)

all_data %>%
  filter(Type == 'Normal') %>%
  select(-Type, -Sample.name, -Age, -Gender) %>%
  colMeans()


pheatmap(t(log(all_data[,1:10] + 1)),
         scale = "row",
         cluster_cols = F)

corr.test(all_data[,1:10], method = "pearson")

##  Summary on genes
p1 <- ggplot(data = all_data, aes(x = RAC2)) +
  geom_density() +
  scale_x_continuous(trans = "log") +
  theme_bw()
p1

p2 <- ggplot(data = all_data, aes(x = PRDM6)) +
  geom_density() +
  scale_x_continuous(trans = "log") +
  theme_bw()
p2



##  Import libraries
library(tidyverse)
library(gam)
library(splines)
library(splines2)
library(pheatmap)
library(pROC)
library(RColorBrewer)

##  Import dataset
all_data <- read.csv("../Processed data.csv", header = T, stringsAsFactors = F)
all_data$Type <- as.factor(all_data$Type)

regression_data <- all_data[, c(1:10, 12)]

regression_data <- regression_data[order(regression_data$DCAKD),]

##  Example: Linear analysis on PRDM6 (the best) / DCAKD (the most related)
model_rough <- lm(PRDM6 ~ ., data = regression_data)
summary(model_rough)

model_linear <- lm(PRDM6 ~ DCAKD, data = regression_data)
#summary(model_linear)
model_spline <- lm(PRDM6 ~ bSpline(DCAKD, df = 4), data = regression_data)

plot(PRDM6 ~ DCAKD, data = regression_data)
lines(regression_data$DCAKD, fitted(model_linear), col="red")
lines(regression_data$DCAKD, fitted(model_spline), col="green")
legend(x="topleft", legend = c("Linear","Spline"),fill=c("red", "green"))

all_data %>% 
  ggplot(aes(x = DCAKD, y = PRDM6, col = Type)) + 
  geom_smooth() + 
  theme_bw() + 
  facet_grid(Type ~ .)
  
  
  
  
  

  
  
  


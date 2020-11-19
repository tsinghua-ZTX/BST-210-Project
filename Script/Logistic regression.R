##  Import libraries
library(tidyverse)
library(pheatmap)
library(pROC)
library(RColorBrewer)
library(ResourceSelection)

##  Import dataset
all_data <- read.csv("../Processed data.csv", header = T, stringsAsFactors = F)
all_data$Type <- as.factor(all_data$Type)

regression_data <- all_data[, c(1:10, 14)]

##  Rough Regression: single gene model
##    Logistic models
{
  PRDM6_mod <- glm(Type ~ PRDM6, family = binomial(), data = regression_data) 
  CAMK2A_mod <- glm(Type ~ CAMK2A, family = binomial(), data = regression_data) 
  RAC2_mod <- glm(Type ~ RAC2, family = binomial(), data = regression_data) 
  KRAS_mod <- glm(Type ~ KRAS, family = binomial(), data = regression_data) 
  DCAKD_mod <- glm(Type ~ DCAKD, family = binomial(), data = regression_data) 
  TMSB4X_mod <- glm(Type ~ TMSB4X, family = binomial(), data = regression_data) 
  MT.RNR2_mod <- glm(Type ~ MT.RNR2, family = binomial(), data = regression_data) 
  MTND4P12_mod <- glm(Type ~ MTND4P12, family = binomial(), data = regression_data) 
  RN7SL2_mod <- glm(Type ~ RN7SL2, family = binomial(), data = regression_data) 
  AL353644.10_mod <- glm(Type ~ AL353644.10, family = binomial(), data = regression_data) 
}

##    ROC curve
{
  PRDM6_roc <- roc(regression_data$Type, predict(PRDM6_mod),
                   levels = c("Normal", "Cancer"))
  CAMK2A_roc <- roc(regression_data$Type, predict(CAMK2A_mod),
                    levels = c("Normal", "Cancer"))
  RAC2_roc <- roc(regression_data$Type, predict(RAC2_mod),
                  levels = c("Normal", "Cancer"))
  KRAS_roc <- roc(regression_data$Type, predict(KRAS_mod),
                  levels = c("Normal", "Cancer"))
  DCAKD_roc <- roc(regression_data$Type, predict(DCAKD_mod),
                   levels = c("Normal", "Cancer"))
  TMSB4X_roc <- roc(regression_data$Type, predict(TMSB4X_mod),
                    levels = c("Normal", "Cancer"))
  MT.RNR2_roc <- roc(regression_data$Type, predict(MT.RNR2_mod),
                     levels = c("Normal", "Cancer"))
  MTND4P12_roc <- roc(regression_data$Type, predict(MTND4P12_mod),
                      levels = c("Normal", "Cancer"))
  RN7SL2_roc <- roc(regression_data$Type, predict(RN7SL2_mod),
                    levels = c("Normal", "Cancer"))
  AL353644.10_roc <- roc(regression_data$Type, predict(AL353644.10_mod),
                         levels = c("Normal", "Cancer"))
}

##    ROC plot
{
  roc_color <- brewer.pal(9,"Set1")
  
  pdf("../Plot/Logistic regression/ROC plot 1.pdf", 6 ,5)
  plot(PRDM6_roc, col = roc_color[1])
  plot.roc(CAMK2A_roc, add = TRUE, col = roc_color[2]) 
  plot.roc(RAC2_roc, add = TRUE, col = roc_color[3]) 
  plot.roc(KRAS_roc, add = TRUE, col = roc_color[4]) 
  plot.roc(DCAKD_roc, add = TRUE, col = roc_color[5]) 
  legend("bottomright", 
         legend = names(regression_data)[1:5],
         col = roc_color[1:5], 
         lwd = 2)
  dev.off()
  
  pdf("../Plot/Logistic regression/ROC plot 2.pdf", 6 ,5)
  plot(TMSB4X_roc, col = roc_color[6])
  plot.roc(MT.RNR2_roc, add = TRUE, col = roc_color[7]) 
  plot.roc(MTND4P12_roc, add = TRUE, col = roc_color[8]) 
  plot.roc(RN7SL2_roc, add = TRUE, col = roc_color[9]) 
  plot.roc(DCAKD_roc, add = TRUE, col = roc_color[1]) 
  legend("bottomright",
         legend = names(regression_data)[c(6:9, 1)],
         col = roc_color[c(6:9, 1)], lwd = 2)
  dev.off()
}

##    H-L GOF test
{
  hoslem.test(all_data$Type, fitted(PRDM6_mod))
  hoslem.test(all_data$Type, fitted(CAMK2A_mod))
  hoslem.test(all_data$Type, fitted(RAC2_mod))
  hoslem.test(all_data$Type, fitted(KRAS_mod))
  hoslem.test(all_data$Type, fitted(DCAKD_mod))
  hoslem.test(all_data$Type, fitted(TMSB4X_mod))
  hoslem.test(all_data$Type, fitted(MT.RNR2_mod))
  hoslem.test(all_data$Type, fitted(MTND4P12_mod))
  hoslem.test(all_data$Type, fitted(RN7SL2_mod))
  hoslem.test(all_data$Type, fitted(AL353644.10_mod))
}

##    Age related
{
  age_mod <- glm(Type ~ Age, family = binomial(), data = all_data) 
  summary(age_mod)
  
  ##  Linear model (attempt)
  test_model <- lm(Age ~ RAC2, data = all_data)
  summary(test_model)
}


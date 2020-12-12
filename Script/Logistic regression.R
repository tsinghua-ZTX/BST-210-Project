##  Import libraries
library(tidyverse)
library(caret)
library(glmnet)
library(lars)
library(pheatmap)
library(pROC)
library(RColorBrewer)
library(ResourceSelection)

##  Import dataset
{
  all_data <- read.csv("../Processed data.csv", header = T, stringsAsFactors = F)
  all_data$Type <- factor(all_data$Type, levels = c("Normal", "Cancer"))
  
  # Data including only genes and types
  regression_data <- all_data[, c(1:10, 14)]
}

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
         legend = names(regression_data)[c(6:10)],
         col = roc_color[c(6:9, 1)], lwd = 2)
  dev.off()
}

##    H-L GOF test (unnecessary)
{
  # hoslem.test(all_data$Type, fitted(PRDM6_mod))
  # hoslem.test(all_data$Type, fitted(CAMK2A_mod))
  # hoslem.test(all_data$Type, fitted(RAC2_mod))
  # hoslem.test(all_data$Type, fitted(KRAS_mod))
  # hoslem.test(all_data$Type, fitted(DCAKD_mod))
  # hoslem.test(all_data$Type, fitted(TMSB4X_mod))
  # hoslem.test(all_data$Type, fitted(MT.RNR2_mod))
  # hoslem.test(all_data$Type, fitted(MTND4P12_mod))
  # hoslem.test(all_data$Type, fitted(RN7SL2_mod))
  # hoslem.test(all_data$Type, fitted(AL353644.10_mod))
}

##    Age related
{
  age_mod <- glm(Type ~ Age, family = binomial(), data = all_data) 
  summary(age_mod)
  
  ##  Linear model (attempt)
  test_model <- lm(Age ~ RAC2, data = all_data)
  summary(test_model)
}


##  Model Selection: stepwise and LASSO
{
  ##    Both stepwise regression (default)
  {
    full_model <- glm(Type ~ ., family = binomial(), data = regression_data) 
    both_step <- step(full_model)
    summary(both_step)
    # KRAS and AL353644.10 
  }
  
  ##    Backward stepwise regression
  {
    back_step <- step(full_model, direction = "backward")
    summary(back_step)
    # KRAS and AL353644.10 
    # Results similar to "both" method
  }
  
  ##    Special model: without KRAS
  {
    KRAS_lack_model <- glm(Type ~ ., family = binomial(), data = regression_data[, -4]) 
    both_step <- step(KRAS_lack_model)
    summary(both_step)
  }
  
  ##    LASSO
  {
    set.seed(1234)
    x <- as.matrix(regression_data[,1:10])
    y <- regression_data[,11]
    
    lar1 <- lars(x, as.matrix(as.numeric(y)))
    lar1
    plot(lar1) 
    summary(lar1)
    lar1$beta[10, ]
    coef <- coef.lars(lar1, mode="step", s=10) 
    predict(lar1, x, s = 10)$fit
    
    lasso_pred <- data.frame(LASSO = predict(lar1, x, s = 10)$fit, Type = all_data$Type)
    ggplot(data = lasso_pred, aes(y = LASSO, x = Type, col = Type)) + 
      geom_boxplot() + 
      theme_bw()
    
    # fit <- glmnet(x, y, alpha = 1, family = 'binomial')
    # fit_cv <- cv.glmnet(x, y, alpha = 1, family = 'binomial', type.measure = 'class')
    # plot(fit_cv)
    # pdf("../Plot/Logistic regression/LASSO_Lambda_CV.pdf", 6, 5)
    # plot(fit_cv)
    # dev.off()
    # coef(fit_cv)
    ##  PRDM6, DCAKD, TMSB4X are chosen
  }
  
  
}


##  Model 1: KRAS and AL353644.10
{
  model_1 <- glm(Type ~ KRAS + AL353644.10, family = binomial(), data = regression_data) 
  summary(model_1)
  model_1_roc <- roc(regression_data$Type, 
                     predict(model_1),
                     levels = c("Normal", "Cancer"))
  # AUROC is 1: too good to be true; also, there seems to be a HUGE collinearity in model 1
  # Why?
  # Shown in point plot
  KA_plot <- regression_data %>%
    ggplot(aes(x = KRAS, y = AL353644.10, col = Type)) +
    geom_point() + 
    labs(x = "KRAS (TPM)", y = "AL353644.10 (TPM)") + 
    scale_color_brewer(type = "div", palette = "Set1") + 
    theme_bw()
  pdf("../Plot/Logistic regression/KA_plot.pdf", 6, 5)
  KA_plot
  dev.off()
  
  ##    Strangely, in cancer groups, there's a linear relation between KRAS and AL353644.10, but not in 
  ##  the normal group; 
  
  # Analysis only on cancer group
  KA_cancer <- lm(KRAS ~ AL353644.10, data = regression_data[which(regression_data$Type == "Cancer"), ])
  summary(KA_cancer)
  # \beta = 0.031548, with p-value <2e-16; a strong linear relationship.
  
  KA_normal <- lm(KRAS ~ AL353644.10, data = regression_data[which(regression_data$Type == "Normal"), ])
  summary(KA_normal)
  # \beta = -5.773e-05, with p-value = 0.75; no linear significance.
}

# Since the too-good performance, maybe we should consider other genes
##  Model 2: PRDM6 and RAC2, similar as model 1
{
  model_2 <- glm(Type ~ PRDM6 + RAC2, family = binomial(), data = regression_data) 
  summary(model_2)
  model_2_roc <- roc(regression_data$Type, 
                     predict(model_2),
                     levels = c("Normal", "Cancer"))
  PR_plot <- regression_data %>%
    ggplot(aes(x = PRDM6, y = RAC2, col = Type)) +
    geom_point() + 
    labs(x = "PRDM6 (TPM)", y = "RAC2 (TPM)") + 
    scale_color_brewer(type = "div", palette = "Set1") + 
    theme_bw()
  pdf("../Plot/Logistic regression/PR_plot.pdf", 6, 5)
  PR_plot
  dev.off()  
}

##  Model 3 (LASSO): PRDM6, DCAKD and TMSB4X
{
  model_3_roc <- roc(regression_data$Type, 
                     as.vector(predict(lar1, x, s = 10)$fit),
                     levels = c("Normal", "Cancer"))
  PDT_ROC <- plot(model_3_roc, 
                  col = roc_color[1],
                  print.auc = T,
                  print.auc.x = 0.2,
                  print.auc.y = 0.05)
  pdf("../Plot/Logistic regression/PDT_ROC.pdf", 5, 5)
  plot(model_3_roc, 
       col = roc_color[1],
       print.auc = T,
       print.auc.x = 0.25,
       print.auc.y = 0.05)
  dev.off()    
}


#   In general, all models could have a very good performance in estimating the binary outcome 
# by logistic regression; 

##  Model 4: 1 + age
{
  model_4 <- glm(Type ~ KRAS + AL353644.10 + Age, family = binomial(), data = all_data) 
  summary(model_4)
}


##  test

ggplot(data = all_data) + 
  geom_point(aes(x = DCAKD, y = CAMK2A, col = Type)) + 
  theme_bw()







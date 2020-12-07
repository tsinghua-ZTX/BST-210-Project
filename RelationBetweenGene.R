##  Import libraries
library(tidyverse)
library(gam)
library(splines)
library(splines2)
library(pheatmap)
library(pROC)
library(RColorBrewer)

##  Import dataset
all_data <- read.csv("Processed data.csv", header = T, stringsAsFactors = F)
all_data$Type <- as.factor(all_data$Type)

regression_data <- all_data[, c(1:10, 12)]

#Linear analysis on PRDM6/DCAKD 
model_rough_PRDM6 <- lm(PRDM6 ~ ., data = regression_data)
summary(model_rough_PRDM6)

model_linear1_PRDM6 <- lm(PRDM6 ~ DCAKD, data = regression_data)
summary(model_linear1_PRDM6)

plot(PRDM6 ~ DCAKD, data = regression_data)
lines(regression_data$DCAKD, fitted(model_linear1_PRDM6), col="red")

all_data %>% 
  ggplot(aes(x = DCAKD, y = PRDM6, col = Type)) + 
  geom_smooth() + 
  theme_bw() + 
  facet_grid(Type ~ .)

#Linear analysis on PRDM6/CAMK2A
model_linear2_PRDM6 <- lm(PRDM6 ~ CAMK2A, data = regression_data)
summary(model_linear2_PRDM6)

plot(PRDM6 ~ CAMK2A, data = regression_data)
lines(regression_data$CAMK2A, fitted(model_linear2_PRDM6), col="red")

all_data %>% 
  ggplot(aes(x = CAMK2A, y = PRDM6, col = Type)) + 
  geom_smooth() + 
  theme_bw() + 
  facet_grid(Type ~ .)

#Linear analysis on PRDM6/RAC2
model_linear3_PRDM6 <- lm(PRDM6 ~ RAC2, data = regression_data)
summary(model_linear3_PRDM6)

plot(PRDM6 ~ RAC2, data = regression_data)
lines(regression_data$RAC2, fitted(model_linear3_PRDM6), col="red")

all_data %>% 
  ggplot(aes(x = RAC2, y = PRDM6, col = Type)) + 
  geom_smooth() + 
  theme_bw() + 
  facet_grid(Type ~ .)

#Linear analysis on PRDM6/KRAS
model_linear4_PRDM6 <- lm(PRDM6 ~ KRAS, data = regression_data)
summary(model_linear4_PRDM6)

plot(PRDM6 ~ KRAS, data = regression_data)
lines(regression_data$KRAS, fitted(model_linear4_PRDM6), col="red")

all_data %>% 
  ggplot(aes(x = KRAS, y = PRDM6, col = Type)) + 
  geom_smooth() + 
  theme_bw() + 
  facet_grid(Type ~ .)

#Linear analysis on PRDM6/TMSB4X
model_linear5_PRDM6 <- lm(PRDM6 ~ TMSB4X, data = regression_data)
summary(model_linear5_PRDM6)

plot(PRDM6 ~ TMSB4X, data = regression_data)
lines(regression_data$TMSB4X, fitted(model_linear5_PRDM6), col="red")

all_data %>% 
  ggplot(aes(x = TMSB4X, y = PRDM6, col = Type)) + 
  geom_smooth() + 
  theme_bw() + 
  facet_grid(Type ~ .)

#Linear analysis on PRDM6/MT.RNR2
model_linear6_PRDM6 <- lm(PRDM6 ~ MT.RNR2, data = regression_data)
summary(model_linear6_PRDM6)

plot(PRDM6 ~ MT.RNR2, data = regression_data)
lines(regression_data$MT.RNR2, fitted(model_linear6_PRDM6), col="red")

all_data %>% 
  ggplot(aes(x = MT.RNR2, y = PRDM6, col = Type)) + 
  geom_smooth() + 
  theme_bw() + 
  facet_grid(Type ~ .)

#Linear analysis on PRDM6/MTND4P12
model_linear7_PRDM6 <- lm(PRDM6 ~ MTND4P12, data = regression_data)
summary(model_linear7_PRDM6)

plot(PRDM6 ~ MTND4P12, data = regression_data)
lines(regression_data$MTND4P12, fitted(model_linear7_PRDM6), col="red")

all_data %>% 
  ggplot(aes(x = MTND4P12, y = PRDM6, col = Type)) + 
  geom_smooth() + 
  theme_bw() + 
  facet_grid(Type ~ .)

#Linear analysis on PRDM6/RN7SL2
model_linear8_PRDM6 <- lm(PRDM6 ~ RN7SL2, data = regression_data)
summary(model_linear8_PRDM6)

plot(PRDM6 ~ RN7SL2, data = regression_data)
lines(regression_data$RN7SL2, fitted(model_linear8_PRDM6), col="red")

all_data %>% 
  ggplot(aes(x = RN7SL2, y = PRDM6, col = Type)) + 
  geom_smooth() + 
  theme_bw() + 
  facet_grid(Type ~ .)

#Linear analysis on PRDM6/AL353644.10
model_linear9_PRDM6 <- lm(PRDM6 ~ AL353644.10, data = regression_data)
summary(model_linear9_PRDM6)

plot(PRDM6 ~ AL353644.10, data = regression_data)
lines(regression_data$AL353644.10, fitted(model_linear9_PRDM6), col="red")

all_data %>% 
  ggplot(aes(x = AL353644.10, y = PRDM6, col = Type)) + 
  geom_smooth() + 
  theme_bw() + 
  facet_grid(Type ~ .)

#Linear analysis on CAMK2A
model_rough_CAMK2A <- lm(CAMK2A ~ ., data = regression_data)
summary(model_rough_CAMK2A)
 

#Linear analysis on CAMK2A/RAC2
model_linear1_CAMK2A <- lm(CAMK2A ~ RAC2, data = regression_data)
summary(model_linear1_CAMK2A)

plot(CAMK2A ~ RAC2, data = regression_data)
lines(regression_data$RAC2, fitted(model_linear1_CAMK2A), col="red")

all_data %>% 
  ggplot(aes(x = RAC2, y = CAMK2A, col = Type)) + 
  geom_smooth() + 
  theme_bw() + 
  facet_grid(Type ~ .)

#Linear analysis on CAMK2A/KRAS
model_linear2_CAMK2A <- lm(CAMK2A ~ KRAS, data = regression_data)
summary(model_linear2_CAMK2A)

plot(CAMK2A ~ KRAS, data = regression_data)
lines(regression_data$KRAS, fitted(model_linear2_CAMK2A), col="red")

all_data %>% 
  ggplot(aes(x = KRAS, y = CAMK2A, col = Type)) + 
  geom_smooth() + 
  theme_bw() + 
  facet_grid(Type ~ .)

#Linear analysis on CAMK2A/DCAKD
model_linear3_CAMK2A <- lm(CAMK2A ~ DCAKD, data = regression_data)
summary(model_linear3_CAMK2A)

plot(CAMK2A ~ DCAKD, data = regression_data)
lines(regression_data$DCAKD, fitted(model_linear3_CAMK2A), col="red")

all_data %>% 
  ggplot(aes(x = DCAKD, y = CAMK2A, col = Type)) + 
  geom_smooth() + 
  theme_bw() + 
  facet_grid(Type ~ .)

#Linear analysis on CAMK2A/TMSB4X
model_linear4_CAMK2A <- lm(CAMK2A ~ TMSB4X, data = regression_data)
summary(model_linear4_CAMK2A)

plot(CAMK2A ~ TMSB4X, data = regression_data)
lines(regression_data$TMSB4X, fitted(model_linear4_CAMK2A), col="red")

all_data %>% 
  ggplot(aes(x = TMSB4X, y = CAMK2A, col = Type)) + 
  geom_smooth() + 
  theme_bw() + 
  facet_grid(Type ~ .)

#Linear analysis on CAMK2A/MT.RNR2
model_linear5_CAMK2A <- lm(CAMK2A ~ MT.RNR2, data = regression_data)
summary(model_linear5_CAMK2A)

plot(CAMK2A ~ MT.RNR2, data = regression_data)
lines(regression_data$MT.RNR2, fitted(model_linear5_CAMK2A), col="red")

all_data %>% 
  ggplot(aes(x = MT.RNR2, y = CAMK2A, col = Type)) + 
  geom_smooth() + 
  theme_bw() + 
  facet_grid(Type ~ .)

#Linear analysis on CAMK2A/MTND4P12
model_linear6_CAMK2A <- lm(CAMK2A ~ MTND4P12, data = regression_data)
summary(model_linear6_CAMK2A)

plot(CAMK2A ~ MTND4P12, data = regression_data)
lines(regression_data$MTND4P12, fitted(model_linear6_CAMK2A), col="red")

all_data %>% 
  ggplot(aes(x = MTND4P12, y = CAMK2A, col = Type)) + 
  geom_smooth() + 
  theme_bw() + 
  facet_grid(Type ~ .)

#Linear analysis on CAMK2A/RN7SL2
model_linear7_CAMK2A <- lm(CAMK2A ~ RN7SL2, data = regression_data)
summary(model_linear7_CAMK2A)

plot(CAMK2A ~ RN7SL2, data = regression_data)
lines(regression_data$RN7SL2, fitted(model_linear7_CAMK2A), col="red")

all_data %>% 
  ggplot(aes(x = RN7SL2, y = CAMK2A, col = Type)) + 
  geom_smooth() + 
  theme_bw() + 
  facet_grid(Type ~ .)

#Linear analysis on CAMK2A/AL353644.10
model_linear8_CAMK2A <- lm(CAMK2A ~ AL353644.10, data = regression_data)
summary(model_linear8_CAMK2A)

plot(CAMK2A ~ AL353644.10, data = regression_data)
lines(regression_data$AL353644.10, fitted(model_linear8_CAMK2A), col="red")

all_data %>% 
  ggplot(aes(x = AL353644.10, y = CAMK2A, col = Type)) + 
  geom_smooth() + 
  theme_bw() + 
  facet_grid(Type ~ .)

#Linear analysis on RAC2
model_rough_RAC2 <- lm(RAC2 ~ ., data = regression_data)
summary(model_rough_RAC2)

#Linear analysis on RAC2/KRAS
model_linear1_RAC2 <- lm(RAC2 ~ KRAS, data = regression_data)
summary(model_linear1_RAC2)

plot(RAC2 ~ KRAS, data = regression_data)
lines(regression_data$KRAS, fitted(model_linear1_RAC2), col="red")

all_data %>% 
  ggplot(aes(x = KRAS, y = RAC2, col = Type)) + 
  geom_smooth() + 
  theme_bw() + 
  facet_grid(Type ~ .)


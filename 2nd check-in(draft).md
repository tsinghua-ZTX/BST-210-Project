## Second project check-in

###1. Group Identification

**Group 15**

Yaoqi Zhang 

Yu Zhang 

Tianxiao Zhao

### 2. Background knowledge

&emsp;a. 

- There are many researches studying on similar problems. Those researches mostly used the RNAseq results from tissue as the source of expression data, but in our study we are focusing on extracellular RNAs, which are present in serum and thus become a much better biomarker than those in previous studies.
  - Until now, there is no widely accepted exRNA biomarkers on breast cancers, also letting us be interest in this project to find out if any exRNA could have such effect.
- Breast cancer is also a hot topic in cancer-related field. Previous studies gave a much similar analysis method as we've proposed in 1st check-in: after differential analysis, chosen genes were used as variables in logistic regression models in order to achieve a classification between BC and normal samples. Meanwhile, some further approaches, like machine learning and SVM, were also used in studies. 
  - If possible and necessary, those further methods may be applied in this project, but only as supplementary.

&emsp;b.

- We've reviewed the review(s) of our 1st Project Check-In, and we are glad to have such helpful reviews on our project.
  - Most reviews concern with the number of variables (genes), which is also what we've concerned. Some great ideas like volcano plot and clustering by RNA types are very helpful for us in choosing highly variable features for the analysis. 
  - We are now planning to dealing with variables in two ways:
    - By differential analysis, we shall select the highly variable features considering both fold-change and FDR. After getting this new dataset with much less but highly variable genes, logistic regression and downstream analysis would be done.
    - Calculate the percentage of RNA types, use these secondary variables as another dataset, apply logistic regression and to see if that can reveal more information. Further clustering on variables would also be possible if needed.
  - Multicollinearity is another concern. We believe that would be a key part in both two approaches above. Advice like correlation test and clustering would be helpful, and we are sure that those would be in our project.

### 3. Analysis Plan

- Firstly, use edgeR, a biostatistical packages, to select the highly variable features between two groups. 
- Visualize the overall structure of the data, determine if outliers exist. Also, check how differential genes distribute over samples.
- With the new data, apply logistic regression to build a classfication model, thus try to figure out which genes contributes the most and the most significant to the occurence of breast cancer.
  - Logistic regression on each gene to determine the independent effect of each gene; also, to figure out the order of importance;
  - Forward selection to build a multi-gene regression model; 
  - Use published BC biomarkers as controls, estimate their effect by logistic regression as comparisons to our selected genes
- To avoid collinearity, check the relations between remain variables by correlation test, VIF and other methods and to decide if further filtration is needed.
- Use statistical methods, test whether gender or age is a confounder or effect modifiers; if so, adjust our model by stratification on those variables to build a covariable-included model.
  - Add age into each one-gene model to see if age plays as a confounder or effect modifiers
  - For gender, check if genes are gender-related (difference in both cancer-female/normal-female and normal-female/normal-male), or not (difference in cancer-female/normal-female, but not in normal-female/normal-male)
- After determing the biomarkers from the logistic model, apply correlation test, clustering analysis or even biological research to find out the relationship among those genes. 
- Calculate the percentage of each exRNA types as a RNA-type variables, use them to see if serum of BC samples has a different exRNA type pattern. If so, identify how it works and why (try our best). 

### 4. Missing Data

_There is only missing in X, but not in Y. So, the following analysis is just dealing with X, but not with MCAR, MAR or MNAR_

&emsp;a.

- Generally, since we are using a well-published dataset, there is no typical missing data in our selected data:
  - For every samples, we know their type (BC or normal), as well as their age and gender;
- Still, if we exam all 60675 genes, there are lots of 0 value in the huge matrix. Most of them are because those genes' expressions are too low and the sequencing method did not catch them. It is nearly impossible to exam each gene carefully to decide whether it's MCAR, MAR or MNAR. So, we've decided to do a first round of filtration, as said in the modeling part, to avoid those low-quality varaibles. 
  - For gene TPM (X variables), in many cases 0 does not mean missing. For some specific genes, they did have no transcription for BC samples;

&emsp;b. 

- For most genes, only few samples have 0 value while other most samples are normal. Some genes like DCAKD have more 0 values, but these genes are high related to BC and all these 0 value are in BC samples. So, right now without further analysis, we have no confidence to do imputation on these variable values. 

  <img src="/Users/zhaotianxiao/Desktop/Project/2nd check-in/Plot/Missing data/Count of missing data.jpg" alt="Count of missing data" style="zoom:50%;" />

- Due to this special occasion, our first plan is to consider all 0 in X as not missing and continue applying logistic regression. If any problem or abnormality occurs, we shall return on these 0 values and to see if imputation is needed.



### 5. Modeling _**all plots and tables are just drafts_

&emsp;a. <u>Linear, flexible/additive or other methods (LASSO, ridge) from this topic</u>

- Our analysis plan does not involve these methods, since the outcome Y is a binary term.
- One of our thought is to explore if there is any genes' expression in exRNA that can reflect samples' age, and whether those genes actually involve in BC. It may seem not-so-related, though, while it maybe still worth a try, by applying LASSO.



&emsp;b. <u>Logistic</u>

- Before constructing any model, we've first done a general filtration to remove those genes that have less expression level (sum <  128 * 20 among all samples), and those genes that have too much 0 values (#0 > 128 * 20%). After these filtrations, the number of genes dropped from 60675 to __1873__. After differential analysis, a first-verion of 10 genes were selected into the rough model part. Those include 4 down-regulated genes and 6 up-regulated genes.

  - The volcano plot below gives the criteria on finding the top high-variance genes: 

    <img src="/Users/zhaotianxiao/Desktop/Project/2nd check-in/Plot/Data process/Volcano plot.jpg" alt="avatar" style="zoom:50%;" />

  - Here are the gene ID of 10 selected genes:

  | Gene symbol ID | Gene Ensembl ID | Regulation |
  | -------------- | --------------- | ---------- |
  | PRDM6          | ENSG00000061455 | Up         |
  | CAMK2A         | ENSG00000070808 | Down       |
  | RAC2           | ENSG00000128340 | Down       |
  | KRAS           | ENSG00000133703 | Down       |
  | DCAKD          | ENSG00000172992 | Up         |
  | TMSB4X         | ENSG00000205542 | Up         |
  | MT-RNR2        | ENSG00000210082 | UP         |
  | MTND4P12       | ENSG00000247627 | Down       |
  | RN7SL2         | ENSG00000274012 | Up         |
  | AL353644.10    | ENSG00000277048 | Up         |

- The first step is to build up a rough model: for 10 genes, build up a simple logistic regression model that contains only the chosen gene alone.

  - For each model, the coefficient and p-value are listed below:

  | Gene ID       | Mean TPM   | Estimate of Coefficient | Odds Ratio | Std. Error | P-value  |
  | ------------- | ---------- | ----------------------- | ---------- | ---------- | -------- |
  | PRDM6         | 36.46756   | 0.16101                 | 1.174697   | 0.03915    | 3.92e-05 |
  | CAMK2A        | 86.35928   | -0.17646                | 0.8382323  | 0.04215    | 2.84e-05 |
  | RAC2          | 300.68722  | -0.029634               | 0.9708008  | 0.007329   | 5.27e-05 |
  | KRAS          | 102.15116  | -0.18815                | 0.8284904  | 0.04249    | 9.52e-06 |
  | DCAKD         | 24.03341   | 0.13591                 | 1.145579   | 0.02865    | 2.09e-06 |
  | TMSB4X        | 107.00828  | 0.010045                | 1.010096   | 0.002136   | 2.57e-06 |
  | MT-RNR2*      | 6496.68522 | 1.836e-04               | 1.000184   | 3.826e-05  | 1.6e-06  |
  | MTND4P12      | 76.93664   | -0.025312               | 0.9750057  | 0.008608   | 0.00328  |
  | RN7SL2        | 683.43801  | 0.0015311               | 1.001532   | 0.0003503  | 1.24e-05 |
  | AL353644.10** | 5918.46639 | 1.975e-04               | 1.000198   | 3.999e-05  | 7.90e-07 |

  - It can be seen that the coefficient is greatly affected by the mean TPM level. Mitochondrial gene* and pseudogene** have a much higher expression than normal genes, which tells us a further separation of gene types would be essential.
  - The interpretation of coefficients is not as usual models we've learnt. 
    - Since the varaibles are TPM, which had already normalized within samples, a one-unit increase in variables indicates a relatively higher expression among all genome, but not an absolute value. So, it can be said that for each gene above, a 0.0001% relative more expression among your serum exRNA would bring the odds ratio calculated in the table.

- To evaluate models, we've calculated and plotted the ROC curve for each genes:

  <img src="/Users/zhaotianxiao/Desktop/Project/2nd check-in/Plot/Logistic regression/ROC plot 1.jpg" alt="ROC plot 1" style="zoom: 50%;" />

  <img src="/Users/zhaotianxiao/Desktop/Project/2nd check-in/Plot/Logistic regression/ROC plot 2.jpg" alt="ROC plot 2" style="zoom: 50%;" />

  - We here can see that most genes give a relatively good classification performance in the ROC cruve. Though, there are still place for improvement.

- Also, Hosmer-Lemeshow Goodness of Fit (GOF) tests for each model have been applied. Results all showed that those models provide a good fit to the data. 
- Our next step is to consider apply multiple genes into the regression model using forward selection method. Furthermore, we would find out some reported BC-related genes as comparison.



&emsp;c. <u>Poisson</u>

- No; right now, there is no plan to apply Poisson model in our project: 
  - Our outcome Y is a binomial variable, where Poisson regression is not suitable here;
  - Our variables X are gene expression matrix normalized by TPM. No Poisson regression method could be done on them, as well.



&emsp;d. <u>Survival Analysis</u>

- No, as our data don’t warrant it



### 6. Writing

- __Abstract__

  Breast cancer is one of the three most common cancers worldwide<sup>1</sup>. Extracellular RNA (exRNAs) has been widely studied through high-throughput sequencing as a noninvasive biomarker for liquid biopsy. Here, we applied a series of statistical studies, including logistic regression, confounding and effect modification analysis, on a published exRNA expression data<sup>2</sup> to find out how exRNA expressions are associated with breast cancer, given the information of age and gender. As a result, _TBC_

- __Introduction and Motivation__ 

  With high incidence and mortality, breast cancer is considered as one of the most common cancers among women. One of our major goals is to detect cancer earlier with less cost in order to prevent deaths from breast cancer. 

  In the past years, developed methods like mammograms and ultrasound had been successfully used in early screening of breast cancers. But still, an effective and inexpensive appoarch of biopsy is needed. Liquid biopsy, detecting including circulating tumor cells, exosomes, extracellular DNAs and extracellular RNAs (exRNAs), has nowadays become a rapid-developing in vitro diagnostic means because of its accessibility and commercial advantages. 

  In this study, we assessed the assocaition between the exRNA expression in serum samples and breast cancer, using multiple statistical methods based on a published exRNA expression dataset<sup>2</sup>. In general, high variance genes were selected by differential analysis and logistic regression were then used as the classification model. Patients' age and gender were collected as covariables and also considered in the adjusted model. _TBC_

- __Data description__

  We used a established exRNA expression dataset from a published study<sup>2</sup>. The dataset included the expression of totally 60675 genes among 96 serum samples with breast cancer and 32 normal serum samples as controls. All expression data are normalized into TPM (Transcripts Per Million reads) so that they can be compared between samples. The gene names are all present in the format of Ensembl gene ID. Metadata including age and gender are also present in the dataset. _TBC_



[1] Harbeck N, Gnant M. Breast cancer. Lancet. 2017 Mar 18;389(10074):1134-1150. doi: 10.1016/S0140-6736(16)31891-8. Epub 2016 Nov 17. PMID: 27865536.

[2] Zixu Z, Sheng Z. Extracellular RNA in a single droplet of human serum reflects physiologic and disease states. PNAS. 2019 Sep 116 (38) 19200-19208; DOI: 10.1073/pnas.1908252116



- __Outline of Materials and Methods__
  - Data Collection and Description
    - As in data description
  - Selection of High Variance Features
  - Model for logistic regression
  - Adjustment on Covariables
  - Gene Categories and exRNA Types
  - Classification Analysis
  - Data viability



### 7. Appendix

- Data process.R

```R
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
pdf("Plot/Data process/Volcano plot.pdf", 6, 5)
p1
dev.off()
```

- Missing data.R

```R
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
```



- Logistic regression.R

```R
##  Import libraries
library(tidyverse)
library(pheatmap)
library(pROC)
library(RColorBrewer)
library(ResourceSelection)

##  Import dataset
all_data <- read.csv("Processed data.csv", header = T, stringsAsFactors = F)
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
  
  pdf("Plot/Logistic regression/ROC plot 1.pdf", 6 ,5)
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
  
  pdf("Plot/Logistic regression/ROC plot 2.pdf", 6 ,5)
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
```

### 8. Intentions

Right now, there is no intention to try to reach publication involving work/results from this project. 


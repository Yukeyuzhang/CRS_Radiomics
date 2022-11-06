####################################################################################################################
### title: "Validation-Model test"
### author: "Kzzhu"
### time: "22.07.04"
#################################################################################################################### 
### 1.Packages------------------------------------------------------------------
{
  suppressMessages({
    library(readxl)
    library(stats)
    library(glmnet)
    library(ROCR)
    library(ggpubr) 
    library(cowplot)
    library(data.table)
    library(class)
    library(kknn)
    library(kernlab)
    library(caret)
    library(MASS)
    library(reshape2)
    library(ggplot2)
    library(rpart) 
    library(partykit)
    library(vcd)
    library(rmda)
    library(patchwork)
    library(rms)
    library(rmda)
    library(PredictABEL)
    library(nricens)
    library(pheatmap)
    library(pROC)
    library(gmodels)
    library(ResourceSelection)
  })
  rm(list = ls())
  options(stringsAsFactors = F)
  work_dir <- "E:/Radiomics"
  setwd(work_dir)
}

### 2.Load model/data----------------------------------------------------
{
  load("E:/fitEos.RData")
  load("E:/fitEM.RData")
  load("E:/fitCV.RData")
  load("E:/fitRad.RData")
  load("E:/fitCom.RData")
}
All <- read.table("E:/All.csv", header = T, sep = ",")
data <- read.table("E:/Training.csv", header = T, sep = ",")
data <- read.table("E:/Validation1.csv", header = T, sep = ",")
data <- read.table("E:/Validation2.csv", header = T, sep = ",")
data <- read.table("E:/Validation3.csv", header = T, sep = ",")
dim(All)

### 3.Radscore-----------------------------------------------------------------------
# radiomics-ROC
{
  X_Rad <- All[, 32:dim(All)[2]]
  X_Rad <- data.matrix(X_Rad)
  Y_Rad <- All["Diagnosis"]
  colnames(Y_Rad) <- "Diagnosis"
  Y_Rad <- data.matrix(Y_Rad)
  Rad <- cbind(Y_Rad, X_Rad)
  Rad <- as.data.frame(Rad)
}
{
  Y_Radpredict <- predict(fitRad, X_Rad, s = fitCV$lambda.1se, type = "link")
  KZRad <- cbind(Rad$Diagnosis, Y_Radpredict)
  pred_Rad <- ROCR::prediction(KZRad[, 2], KZRad[, 1])
  auc_Rad <- performance(pred_Rad, "auc")@y.values[[1]]
  perf_Rad <- performance(pred_Rad, "tpr", "fpr")
  plot(perf_Rad, lwd = 2, colorize = FALSE, col = "red", add = T)
  lines(c(0, 1), c(0, 1), lwd = 2, col = "gray", lty = 4 )
  text(0.8, 0.2, cex = c(2, 2), labels = paste0("AUC: ", round(auc_Rad, 3)), col = "red")
}

All <- read.table("E:/All.csv", header = T, sep = ",")
dim(All)

### 4.Radscore-----------------------------------------------------------------------
# radiomics-ROC
{
  X_Rad <- data[, 32:dim(data)[2]]
  X_Rad <- data.matrix(X_Rad)
  Y_Rad <- data["Diagnosis"]
  colnames(Y_Rad) <- "Diagnosis"
  Y_Rad <- data.matrix(Y_Rad)
  Rad <- cbind(Y_Rad, X_Rad)
  Rad <- as.data.frame(Rad)
}
{
  Y_Radpredict <- predict(fitRad, X_Rad, s = fitCV$lambda.1se, type = "link")
  KZRad <- cbind(Rad$Diagnosis, Y_Radpredict)
  pred_Rad <- ROCR::prediction(KZRad[, 2], KZRad[, 1])
  auc_Rad <- performance(pred_Rad, "auc")@y.values[[1]]
  perf_Rad <- performance(pred_Rad, "tpr", "fpr")
  plot(perf_Rad, lwd = 2, colorize = FALSE, col = "red", add = T)
  lines(c(0, 1), c(0, 1), lwd = 2, col = "gray", lty = 4 )
  text(0.8, 0.2, cex = c(2, 2), labels = paste0("AUC: ", round(auc_Rad, 3)), col = "red")
}

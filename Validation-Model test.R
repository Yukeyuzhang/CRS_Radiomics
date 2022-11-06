####################################################################################################################
### title: "Validation"
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

### 2.Load models------------------------------------------------------------------
{
  load("E:/fitEos.RData")
  load("E:/fitEM.RData")
  load("E:/fitCV.RData")
  load("E:/fitRad.RData")
  load("E:/fitCom.RData")
}
Validate <- read.table("E:/Validation1.csv", header = T, sep = ",")
Validate <- read.table("E:/Validation2.csv", header = T, sep = ",")
Validate <- read.table("E:/Validation3.csv", header = T, sep = ",")
dim(Validate)

### 3.Models--------------------------------------------------------------
# Eos-ROC
{
  Y_Eospredict <- predict(fitEos, Validate, type = "link")
  KZEos <- cbind(Validate$Diagnosis, Y_Eospredict) 
  pred_Rad <- ROCR::prediction(KZEos[, 2], KZEos[, 1])
  auc_Rad <- performance(pred_Rad, "auc")@y.values[[1]]
  perf_Rad <- performance(pred_Rad, "tpr", "fpr")
  plot(perf_Rad, lwd = 2, colorize = FALSE, col = "#2258AA")
  lines(c(0, 1), c(0, 1), lwd = 2, col = "gray", lty = 4 )
  text(0.8, 0.3, cex = c(2, 2), labels = paste0("AUC: ", round(auc_Rad, 3)), col = "#2258AA")
}
{
  rocobj <- plot.roc(KZEos[, 1], KZEos[, 2],
                     main = "ROC",
                     print.thres = "best",
                     percent = T, ci = T, print.auc = T)
  Y_Eosclass <- as.numeric(Y_Eospredict > 0.3291612)
}
table(Y_Eosclass, Validate$Diagnosis)
{
  zhu <- table(Y_Eosclass, Validate$Diagnosis)
  ACC <- (zhu[1,1] + zhu[2,2])/(zhu[1,1] + zhu[1,2] + zhu[2,1] + zhu[2,2])
  Sensitivity <- zhu[2,2]/(zhu[2,2] + zhu[1,2]) 
  Specificity <- zhu[1,1]/(zhu[1,1] + zhu[2,1])
  PPV <- zhu[2,2]/(zhu[2,2] + zhu[2,1])
  NPV <- zhu[1,1]/(zhu[1,1] + zhu[1,2])
  Zhu <- rbind(ACC, Sensitivity, Specificity, PPV, NPV)
}

# EM-ROC
{
  Y_EMpredict <- predict(fitEM, Validate, type = "link")
  KZEM <- cbind(Validate$Diagnosis, Y_EMpredict)
  pred_Rad <- ROCR::prediction(KZEM[, 2], KZEM[, 1])
  auc_Rad <- performance(pred_Rad, "auc")@y.values[[1]]
  perf_Rad <- performance(pred_Rad, "tpr", "fpr")
  plot(perf_Rad, lwd = 2, colorize = FALSE, col = "#2258AA")
  lines(c(0, 1), c(0, 1), lwd = 2, col = "gray", lty = 4 )
  text(0.8, 0.3, cex = c(2, 2), labels = paste0("AUC: ", round(auc_Rad, 3)), col = "#2258AA")
}
{
  rocobj <- plot.roc(KZEM[, 1], KZEM[, 2],
                     main = "ROC",
                     print.thres = "best",
                     percent = T, ci = T, print.auc = T)
  Y_EMclass <- as.numeric(Y_EMpredict > 0.3603362)
}
table(Y_EMclass, Validate$Diagnosis)
{
  zhu <- table(Y_EMclass, Validate$Diagnosis)
  ACC <- (zhu[1,1] + zhu[2,2])/(zhu[1,1] + zhu[1,2] + zhu[2,1] + zhu[2,2])
  Sensitivity <- zhu[2,2]/(zhu[2,2] + zhu[1,2]) 
  Specificity <- zhu[1,1]/(zhu[1,1] + zhu[2,1])
  PPV <- zhu[2,2]/(zhu[2,2] + zhu[2,1])
  NPV <- zhu[1,1]/(zhu[1,1] + zhu[1,2])
  Zhu <- rbind(ACC, Sensitivity, Specificity, PPV, NPV)
}

# radiomics-ROC
{
  X_Rad <- Validate[, 32:dim(Validate)[2]]
  X_Rad <- data.matrix(X_Rad)
  Y_Rad <- Validate["Diagnosis"]
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
{
  rocobj <- plot.roc(KZRad[, 1], KZRad[, 2],
                     main = "ROC",
                     print.thres = "best",
                     percent = T, ci = T, print.auc = T)
  Y_Radclass <- as.numeric(Y_Radpredict > -0.2877397)
}
table(Y_Radclass, Validate$Diagnosis)
{
  zhu <- table(Y_Radclass, Validate$Diagnosis)
  ACC <- (zhu[1,1] + zhu[2,2])/(zhu[1,1] + zhu[1,2] + zhu[2,1] + zhu[2,2])
  Sensitivity <- zhu[2,2]/(zhu[2,2] + zhu[1,2]) 
  Specificity <- zhu[1,1]/(zhu[1,1] + zhu[2,1])
  PPV <- zhu[2,2]/(zhu[2,2] + zhu[2,1])
  NPV <- zhu[1,1]/(zhu[1,1] + zhu[1,2])
  Zhu <- rbind(ACC, Sensitivity, Specificity, PPV, NPV)
}

# Com-ROC
Radscore <- as.numeric(Y_Radpredict)
Validate2 <- cbind(Validate, Radscore)
ddist <- datadist(Validate2)
options(datadist = 'ddist')
{
  Y_Compredict <- predict(fitCom, Validate2, type = "link")
  KZCom <- cbind(Validate$Diagnosis, Y_Compredict)
  pred_Rad <- ROCR::prediction(KZCom[, 2], KZCom[, 1])
  auc_Rad <- performance(pred_Rad, "auc")@y.values[[1]]
  perf_Rad <- performance(pred_Rad, "tpr", "fpr")
  plot(perf_Rad, lwd = 2, colorize = FALSE, col = "black",add = T)
  lines(c(0, 1), c(0, 1), lwd = 2, col = "gray", lty = 4 )
  text(0.8, 0.1, cex = c(2, 2), labels = paste0("AUC: ", round(auc_Rad, 3)), col = "black")
}
{
  rocobj <- plot.roc(KZCom[, 1], KZCom[, 2],
                     main = "ROC",
                     print.thres = "best",
                     percent = T, ci = T, print.auc = T)
  Y_Comclass <- as.numeric(Y_Compredict > 0.4401986)
}
table(Y_Comclass, Validate$Diagnosis)
{
  zhu <- table(Y_Comclass, Validate$Diagnosis)
  ACC <- (zhu[1,1] + zhu[2,2])/(zhu[1,1] + zhu[1,2] + zhu[2,1] + zhu[2,2])
  Sensitivity <- zhu[2,2]/(zhu[2,2] + zhu[1,2]) 
  Specificity <- zhu[1,1]/(zhu[1,1] + zhu[2,1])
  PPV <- zhu[2,2]/(zhu[2,2] + zhu[2,1])
  NPV <- zhu[1,1]/(zhu[1,1] + zhu[1,2])
  Zhu <- rbind(ACC, Sensitivity, Specificity, PPV, NPV)
}

### 4.Boxplot----------------------------------------------------------
{
  KZCom <- as.data.frame(KZCom)
  colnames(KZCom) <- c('event','prob')
  KZCom$event <- as.factor(KZCom$event)
  p.Com <- ggplot(KZCom, aes(x = event, y = prob, fill = event)) +
    stat_boxplot(geom = "errorbar", size = 1, width = 0.4, col = c("#2258AA","#FF0A0A")) +
    geom_boxplot(size = 1, outlier.fill = "white", outlier.size = 0.001,
                 outlier.color = "white", col = c("#2258AA","#FF0A0A"),
                 fill= "white") +
    geom_jitter(aes(fill = event), width = 0.15, shape = 20, size = 1.9) +
    aes(color = event)+
    labs(x = "Diagnosis", y = "Predict") +
    theme_classic() +
    scale_color_manual(values = c("#2258AA","#FF0A0A")) +
    ylim(0, 1) +
    theme(axis.title.x = element_text(face = "bold", size = 28),
          axis.title.y = element_text(face = "bold", size = 28))+
    stat_compare_means()
  p.Com
}

### 5.Confusion----------------------------------------------------------
{ 
  table(Y_Comclass, Validate$Diagnosis)
  a1 <- table(Y_Comclass, Validate$Diagnosis)
  pheatmap(a1,
           cluster_rows = F,
           cluster_cols = F,
           border_color = "#7C858D",
           color = colorRampPalette(c("#E4EBF7","#2258AA"))(100))
}

### 6.Delong test----------------------------------------------------------
roc1 <- roc(Validate$Diagnosis, Y_Eospredict)
roc1 <- roc(Validate$Diagnosis, Y_Radpredict)
roc1 <- roc(Validate$Diagnosis, Y_EMpredict)
roc2 <- roc(Validate$Diagnosis, Y_Compredict)
roc.test(roc1, roc2, method = 'delong')

### 7.Calibration curve-----------------------------------------------------
{
  Y_Comfitted <- predict(fitCom, Validate2, type = "link")
  Com.fitted <- lrm(Diagnosis ~ Y_Comfitted, data = Validate2, x = T, y = T)
  cal1 <- calibrate(Com.fitted, 
                    cmethod = 'hare', method = 'boot', B = 1000, 
                    data = data)
  plot(cal1, xlim = c(0, 1.0), ylim = c(0, 1.0))
  hoslem.test(Validate$Diagnosis, Y_Compredict, g = 10)
}

### 8.DCA curve-------------------------------------------------------------
{
  DCAX <- as.matrix(Validate2[c("Eosinophil","Radscore")])
  DCAY <- Validate2$Diagnosis
  DCAdata <- as.data.frame(cbind(DCAY, DCAX))
  colnames(DCAdata) <- c("Diagnosis","Eosinophil","Radscore")
  Rad <- decision_curve(Diagnosis ~ Radscore,
                        data = DCAdata,
                        family = binomial(link ='logit'),
                        thresholds = seq(0, 0.75, by = 0.01),
                        confidence.intervals = 0.95,
                        study.design =  'cohort')
  Cli <- decision_curve(Diagnosis ~ Eosinophil,
                        data = DCAdata,
                        family = binomial(link ='logit'),
                        thresholds = seq(0, 0.75, by = 0.01),
                        confidence.intervals = 0.95,
                        study.design =  'cohort')
  Com <- decision_curve(Diagnosis ~ Eosinophil + Radscore,
                        data = DCAdata,
                        family = binomial(link ='logit'),
                        thresholds = seq(0, 0.75, by = 0.01),
                        confidence.intervals = 0.95,
                        study.design = 'cohort')
  List <- list(Cli, Rad, Com)
  plot_decision_curve(List, curve.names = c('Eos','Rad', "Com"),
                      lwd = 2.5,
                      cost.benefit.axis = FALSE,
                      col = c('#2258aa','#ff0a0a','#28393E'),
                      confidence.intervals = FALSE,
                      standardize = FALSE)
}

### 9.NRI and IDI-----------------------------------------------------------
{
  z.simple = as.matrix(subset(DCAdata, select = c(Radscore)))
  z.complex = as.matrix(subset(DCAdata, select = c(Eosinophil, Radscore)))
  
  msimple = glm(DCAY~., binomial(logit), data.frame(z.simple), x = T)
  mcomplex = glm(DCAY~., binomial(logit), data.frame(z.complex), x = T)
  
  pRad <- msimple$fitted.values
  pCom <- mcomplex$fitted.values
  reclassification(data = DCAdata,
                   cOutcome = 1,
                   predrisk1 = pRad, 
                   predrisk2 = pCom,
                   cutoff = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))
}
{
  z.simple = as.matrix(subset(DCAdata, select = c(Eosinophil)))
  z.complex = as.matrix(subset(DCAdata, select = c(Eosinophil, Radscore)))
  
  msimple = glm(DCAY~., binomial(logit), data.frame(z.simple), x = T)
  mcomplex = glm(DCAY~., binomial(logit), data.frame(z.complex), x = T)
  
  pRad <- msimple$fitted.values
  pCom <- mcomplex$fitted.values
  reclassification(data = DCAdata,
                   cOutcome = 1,
                   predrisk1 = pRad, 
                   predrisk2 = pCom,
                   cutoff = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))
}

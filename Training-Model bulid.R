####################################################################################################################
### title: "Model Build"
### author: "Kzzhu"
### time: "22.07.04"
####################################################################################################################
### 1.Package----------------------------------------------------------
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
    library(e1071)
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
    library(mRMRe)
  })
  rm(list = ls())
  options(stringsAsFactors = F)
  work_dir <- "E:/Radiomics"
  setwd(work_dir)
}

### 2.Train-val--------------------------------------------------------
Train <- read.table("E:/Training.csv", header = T, sep = ",")
dim(Train)

### 3.Clinical model---------------------------------------------------
# Coefficient
fitEos <- glm(Diagnosis ~ Eosinophil, data = Train, family = binomial())
fitEos$coefficients
# ROC
{
  Y_Eospredict <- predict(fitEos, Train, type = "link")
  KZEos <- cbind(Train$Diagnosis, Y_Eospredict) 
  pred_Rad <- ROCR::prediction(KZEos[, 2], KZEos[, 1])
  auc_Rad <- performance(pred_Rad, "auc")@y.values[[1]]
  perf_Rad <- performance(pred_Rad, "tpr", "fpr")
  plot(perf_Rad, lwd = 2, colorize = FALSE, col = "#2258AA")
  lines(c(0, 1), c(0, 1), lwd = 2, col = "gray", lty = 4 )
  text(0.8, 0.15, cex = c(2, 2), labels = paste0("Eos: ", round(auc_Rad, 3)), col = "#2258AA")
}
{
  rocobj <- plot.roc(KZEos[, 1], KZEos[, 2],
                     main = "ROC",
                     print.thres = "best",
                     percent = T, ci = T, print.auc = T)
  Eosparameters <- coords(rocobj, "best")
  Eoscutoff <-  Eosparameters[1, 1]
  Eoscutoff
  Y_Eosclass <- as.numeric(Y_Eospredict > Eoscutoff)
}
table(Y_Eosclass, Train$Diagnosis)
{
  zhu <- table(Y_Eosclass, Train$Diagnosis)
  ACC <- (zhu[1,1] + zhu[2,2])/(zhu[1,1] + zhu[1,2] + zhu[2,1] + zhu[2,2])
  Sensitivity <- zhu[2,2]/(zhu[2,2] + zhu[1,2]) 
  Specificity <- zhu[1,1]/(zhu[1,1] + zhu[2,1])
  PPV <- zhu[2,2]/(zhu[2,2] + zhu[2,1])
  NPV <- zhu[1,1]/(zhu[1,1] + zhu[1,2])
  Zhu <- rbind(ACC, Sensitivity, Specificity, PPV, NPV)
}

### 4.CT score model---------------------------------------------------
# Coefficient
fitEM <- glm(Diagnosis ~ EM, data = Train, family = binomial())
fitEM$coefficients
# ROC
{
  Y_EMpredict <- predict(fitEM, Train, type = "link")
  KZEM <- cbind(Train$Diagnosis, Y_EMpredict)
  pred_Rad <- ROCR::prediction(KZEM[, 2], KZEM[, 1])
  auc_Rad <- performance(pred_Rad, "auc")@y.values[[1]]
  perf_Rad <- performance(pred_Rad, "tpr", "fpr")
  plot(perf_Rad, lty = 2, lwd = 2, colorize = FALSE, col = "#FF0A0A",add=T)
  lines(c(0, 1), c(0, 1), lwd = 2, col = "gray", lty = 4 )
  text(0.8, 0.35, cex = c(2, 2), labels = paste0("EM: ", round(auc_Rad, 3)), col = "#FF0A0A")
}
{
  rocobj <- plot.roc(KZEM[, 1], KZEM[, 2],
                     main = "ROC",
                     print.thres = "best",
                     percent = T, ci = T, print.auc = T)
  EMparameters <- coords(rocobj, "best")
  EMcutoff <-  EMparameters[1, 1]
  EMcutoff
  Y_EMclass <- as.numeric(Y_EMpredict > Eoscutoff)
}
table(Y_EMclass, Train$Diagnosis)
{
  zhu <- table(Y_EMclass, Train$Diagnosis)
  ACC <- (zhu[1,1] + zhu[2,2])/(zhu[1,1] + zhu[1,2] + zhu[2,1] + zhu[2,2])
  Sensitivity <- zhu[2,2]/(zhu[2,2] + zhu[1,2]) 
  Specificity <- zhu[1,1]/(zhu[1,1] + zhu[2,1])
  PPV <- zhu[2,2]/(zhu[2,2] + zhu[2,1])
  NPV <- zhu[1,1]/(zhu[1,1] + zhu[1,2])
  Zhu <- rbind(ACC, Sensitivity, Specificity, PPV, NPV)
}

### 5.mRMR-------------------------------------------------------------
{
  X_mRMR <- Train[, 32:dim(Train)[2]]
  Y_mRMR <- Train[, "Diagnosis"]
  colnames(Y_mRMR) <- "Diagnosis"
  mRMR <- as.data.frame(cbind(Y_mRMR, X_mRMR))
  
  target_indices <- which(names(Y_mRMR)=='Diagnosis')
  for (m in which(sapply(mRMR, class) != "numeric")){
    mRMR[, m] = as.numeric(mRMR[, m])
  }
  mRMR <- scale(mRMR, center = F, scale = T)
  Data <- mRMR.data(data = data.frame(mRMR))
  mrmr <- mRMR.ensemble(data = Data
                        , target_indices = target_indices
                        , feature_count = 30
                        , solution_count = 1)
  index <- mrmr@filters[[as.character(mrmr@target_indices)]]
  mRMR_feature <- mRMR[,index]
}
mRMR_feature <- scale(mRMR_feature)

### 6.Radiomics model--------------------------------------------------
{
  X_Rad <- mRMR_feature
  X_Rad <- data.matrix(X_Rad)
  Y_Rad <- Train["Diagnosis"]
  colnames(Y_Rad) <- "Diagnosis"
  Y_Rad <- data.matrix(Y_Rad)
  Rad <- cbind(Y_Rad, X_Rad)
  Rad <- as.data.frame(Rad)
}

# LASSO
# set.seed
{
  fitCV <- cv.glmnet(X_Rad, Y_Rad, 
                     family = "binomial", 
                     type.measure = "auc",
                     nfolds = 10)
  plot(fitCV)
}
{
  fitRad <- glmnet(X_Rad, Y_Rad, family = "binomial", alpha = 1)
  print(fitRad)
  plot(fitRad, xvar = "lambda", label = TRUE)
  abline(v = log(fitCV$lambda.1se), lty = 3, lwd = 2, col = "black")

  Coefficients <- coef(fitRad, s = fitCV$lambda.1se)
  Active.Index <- which(Coefficients != 0)
  Active.Coefficients <- Coefficients[Active.Index]
  
  get_coefficients <- function(fitRad, lambda.1se, Coefficients, Active.Index, Active.Coefficients){
    kz <- data.frame(rownames(Coefficients)[Active.Index], Active.Coefficients)
    kz <- data.table('feature' = rownames(Coefficients)[Active.Index],
                     'coef' = Active.Coefficients)
    kz$expcoef <- exp(kz$coef)
    return(kz[order(expcoef)])}
  Feature_Weight <- get_coefficients(fitRad, fitCV$lambda.1se, Coefficients, Active.Index, Active.Coefficients)
  get_coefficients(fitRad, fitCV$lambda.1se, Coefficients, Active.Index, Active.Coefficients)
}
{
  Feature_Weight <- Feature_Weight[2:dim(Feature_Weight)[1], ]
  FW_plot <- ggplot(data = Feature_Weight,
                    mapping = aes(x = feature, y = coef)) +
    geom_bar(stat = "identity", fill = "#2258AA") +
    coord_flip() +
    labs(x = " ", y = "Coefficients") +
    theme_classic() #+
  FW_plot
}
{
  Y_Radpredict <- predict(fitRad, X_Rad, s = fitCV$lambda.1se, type = "link")
  KZRad <- cbind(Rad$Diagnosis, Y_Radpredict)
  pred_Rad <- ROCR::prediction(KZRad[, 2], KZRad[, 1])
  auc_Rad <- performance(pred_Rad, "auc")@y.values[[1]]
  perf_Rad <- performance(pred_Rad, "tpr", "fpr")
  plot(perf_Rad, lwd = 2, colorize = FALSE, col = "#FF0A0A", add = T)
  lines(c(0, 1), c(0, 1), lwd = 2, col = "gray", lty = 4)
  text(0.8, 0.25, cex = c(2, 2), labels = paste0("Rad: ", round(auc_Rad, 3)), col = "#FF0A0A")
}
{
  rocobj <- plot.roc(KZRad[, 1], KZRad[, 2],
                     main = "ROC",
                     print.thres = "best",
                     percent = T, ci = T, print.auc = T)
  Radparameters <- coords(rocobj, "best")
  Radcutoff <- Radparameters[1, 1]
  Radcutoff
  Y_Radclass <- as.numeric(Y_Radpredict > Radcutoff)
}
table(Y_Radclass, Train$Diagnosis)
{
  zhu <- table(Y_Radclass, Train$Diagnosis)
  ACC <- (zhu[1,1] + zhu[2,2])/(zhu[1,1] + zhu[1,2] + zhu[2,1] + zhu[2,2])
  Sensitivity <- zhu[2,2]/(zhu[2,2] + zhu[1,2]) 
  Specificity <- zhu[1,1]/(zhu[1,1] + zhu[2,1])
  PPV <- zhu[2,2]/(zhu[2,2] + zhu[2,1])
  NPV <- zhu[1,1]/(zhu[1,1] + zhu[1,2])
  Zhu <- rbind(ACC, Sensitivity, Specificity, PPV, NPV)
}

### 7.Combined model----------------------------------------------------
Radscore <- as.numeric(Y_Radpredict)
Train2 <- cbind(Train, Radscore)
fitCom <- glm(Diagnosis ~ Eosinophil + Radscore, data = Train2, family = binomial())
fitCom$coefficients
# ROC
{
  Y_Compredict <- predict(fitCom, Train, type = "link")
  KZCom <- cbind(Train$Diagnosis, Y_Compredict)
  pred_Rad <- ROCR::prediction(KZCom[, 2], KZCom[, 1])
  auc_Rad <- performance(pred_Rad, "auc")@y.values[[1]]
  perf_Rad <- performance(pred_Rad, "tpr", "fpr")
  plot(perf_Rad, lwd = 2, colorize = FALSE, col = "black", add = T)
  lines(c(0, 1), c(0, 1), lwd = 2, col = "gray", lty = 4 )
  text(0.8, 0.05, cex = c(2, 2), labels = paste0("Com: ", round(auc_Rad, 3)), col = "black")
}
{
  rocobj <- plot.roc(KZCom[, 1], KZCom[, 2],
                     main = "ROC",
                     print.thres = "best",
                     percent = T, ci = T, print.auc = T)
  Comparameters <- coords(rocobj, "best")
  Comcutoff <-  Comparameters[1, 1]
  Comcutoff
  Y_Comclass <- as.numeric(Y_Compredict > Comcutoff)
}
table(Y_Comclass, Train$Diagnosis)
{
  zhu <- table(Y_Comclass, Train$Diagnosis)
  ACC <- (zhu[1,1] + zhu[2,2])/(zhu[1,1] + zhu[1,2] + zhu[2,1] + zhu[2,2])
  Sensitivity <- zhu[2,2]/(zhu[2,2] + zhu[1,2]) 
  Specificity <- zhu[1,1]/(zhu[1,1] + zhu[2,1])
  PPV <- zhu[2,2]/(zhu[2,2] + zhu[2,1])
  NPV <- zhu[1,1]/(zhu[1,1] + zhu[1,2])
  Zhu <- rbind(ACC, Sensitivity, Specificity, PPV, NPV)
}

# Boxplot---------------------------------------------------------------
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

# Confusion-------------------------------------------------------------
{ 
  table(Y_Comclass, Train$Diagnosis)
  a1 <- table(Y_Comclass, Train$Diagnosis)
  pheatmap(a1,
           cluster_rows = F,
           cluster_cols = F,
           border_color = "#7C858D",
           color = colorRampPalette(c("#E4EBF7","#2258AA"))(100))
}

### 8.Delong test-------------------------------------------------------
roc1 <- roc(Train$Diagnosis, Y_Radpredict)
roc1 <- roc(Train$Diagnosis, Y_Eospredict)
roc1 <- roc(Train$Diagnosis, Y_EMpredict)
roc2 <- roc(Train$Diagnosis, Y_Compredict)
roc.test(roc1, roc2, method = 'delong')

### 9.Nomogram----------------------------------------------------------
{
  ddist <- datadist(Train2)
  options(datadist = 'ddist')
  fitCom <- lrm(Diagnosis ~ Eosinophil + Radscore, 
                 data = Train2, x = T, y = T)
  nom <- nomogram(fitCom
                  , fun = plogis
                  , fun.at = c(.001, .01, .05, seq(.1, .9, by = .1), .95, .99, .999)
                  , lp = F
                  , funlabel = "Diagnosis")
  plot(nom)
}

### 10.Calibration curve------------------------------------------------
{
  Y_Comfitted <- predict(fitCom, Train, type = "link")
  Com.fitted <- lrm(Diagnosis ~ Y_Comfitted, data = Train, x = T, y = T)
  cal1 <- calibrate(Com.fitted, 
                    cmethod = 'hare', method = 'boot', B = 1000, 
                    data = data)
  plot(cal1, xlim = c(0, 1.0), ylim = c(0, 1.0))
  hoslem.test(Train$Diagnosis, Y_Compredict, g = 10)
}

### 11.DCA curve--------------------------------------------------------
{
  DCAX <- as.matrix(Train2[c("Eosinophil","Radscore")])
  DCAY <- Train2$Diagnosis
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
                      cost.benefit.axis = F,
                      col = c('#2258aa','#ff0a0a','#28393E'),
                      confidence.intervals = F,
                      standardize = F)
}

### 12.NRI and IDI------------------------------------------------------
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

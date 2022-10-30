####################################################################################################################
### title: "Figures"
### author: "Kzzhu"
### time: "2022.07.04"
####################################################################################################################
### 1.Packages----------------------------------------------------------------------------------------------------
{
  suppressMessages({
    library(readxl)
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
    library(riskRegression)
  })
  rm(list = ls())
  options(stringsAsFactors = F)
  work_dir <- "E:/Radiomics"
  setwd(work_dir)
}

### 2.Data-----------------------------------------------------------------------
{
  load("E:/fitEos.RData")
  load("E:/fitEM.RData")
  load("E:/fitCV.RData")
  load("E:/fitRad.RData")
  load("E:/fitCom.RData")
  load("E:/fitCom2.RData")
}
train <- read.table("E:/Training.csv", header = T, sep = ",")
validate1 <- read.table("E:/Validation1.csv", header = T, sep = ",")
validate2 <- read.table("E:/Validation2.csv", header = T, sep = ",")
validate3 <- read.table("E:/Validation3.csv", header = T, sep = ",")

### 3.Calibration curve-----------------------------------------------------------
{
  Y_Comfitted <- predict(fitCom, train, type = "response")
  Com.fitted1 <- lrm(Diagnosis ~ Y_Comfitted, data = train, x = T, y = T)
  cal1 <- calibrate(Com.fitted1, B = 1000, 
                    data = train)
  plot(cal1, xlim = c(0, 1.0), ylim = c(0, 1.0))
  
  Y_Comfitted <- predict(fitCom, validate1, type = "response")
  Com.fitted2 <- lrm(Diagnosis ~ Y_Comfitted, data = validate1, x = T, y = T)
  cal2 <- calibrate(Com.fitted2, B = 1000, 
                    data = validate1)
  plot(cal2, xlim = c(0, 1.0), ylim = c(0, 1.0))
  
  Y_Comfitted <- predict(fitCom, validate2, type = "response")
  Com.fitted3 <- lrm(Diagnosis ~ Y_Comfitted, data = validate2, x = T, y = T)
  cal3 <- calibrate(Com.fitted3, B = 1000, 
                    data = validate2)
  plot(cal3, xlim = c(0, 1.0), ylim = c(0, 1.0))
  
  Y_Comfitted <- predict(fitCom, validate4, type = "response")
  Com.fitted4 <- lrm(Diagnosis ~ Y_Comfitted, data = validate4, x = T, y = T)
  cal4 <- calibrate(Com.fitted4, B = 1000, 
                    data = validate4)
  plot(cal4, xlim = c(0, 1.0), ylim = c(0, 1.0))
  
  plot(1, type = "n",
       xlim = c(0,1),
       ylim = c(0,1),
       xlab = "Predicted",
       ylab = "Dignosis",
       legend = F,
       subtitles = F)
  abline(0, 1, col = "black", lty = 2, lwd = 1)
}

### 4.DCA curve-------------------------------------------------------------------
{
  DCAX <- as.matrix(data[c("Eosinophil","Radscore")])
  DCAY <- Y_data
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
  plot_decision_curve(List, curve.names = c('Cli','Rad', "Com"),
                      lwd = 2.5,
                      cost.benefit.axis = FALSE,
                      col = c('#2258aa','#ff0a0a','#28393E'),
                      confidence.intervals = FALSE, 
                      standardize = FALSE)
}

### 5.waterfall-------------------------------------------------------------------
{
# Radcutoff
  data$id = 1:dim(data)[1]
  data$Radscore = data$Radscore
  data$Diagnosis = as.factor(data$Diagnosis)
  ggplot(data, aes(x = reorder(id, Radscore), y = Radscore-Radcutoff, fill = Diagnosis))+
    geom_bar(stat = "identity", width = 0.6) +
    xlab('Patients') + ylab("Rad-Score") +
    theme(axis.text.x = element_blank()) +
    theme(panel.grid = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black")) +
    theme(axis.title = element_text(size = 28))+
    ylim(-2.5, 2.5)+
    scale_fill_manual(values = c("#2258AA","#FF0A0A"),
                      labels = c("NCRSwNP", "ECRSwNP"))
}


### 6.Boxplot---------------------------------------------------------------------
{
  KZ_Com <- as.data.frame(KZ_Com)
  colnames(KZ_Com) <- c('event','prob')
  KZ_Com$event <- as.factor(KZ_Com$event)
  p.Com <- ggplot(KZ_Com, aes(x = event, y = prob, fill = event)) +
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

### 7.Confusion-------------------------------------------------------------------
{ 
  table(Y_Comclass, data$Diagnosis)
  a1 <- table(Y_Comclass, data$Diagnosis)
  pheatmap(a1,
           cluster_rows = F,
           cluster_cols = F,
           border_color = "#7C858D",
           color = colorRampPalette(c("#E4EBF7","#2258AA"))(100))
}

### 8.ROC curve-------------------------------------------------------------------
{
  Y_Compredict_train <- predict(fitCom, train, type = "response")
  KZ_Com_train <- cbind(train$Diagnosis, Y_Compredict_train) 
  pred_Com_train <- ROCR::prediction(KZ_Com_train[, 2], KZ_Com_train[, 1])
  auc_Com_train <- performance(pred_Com_train, "auc")@y.values[[1]]
  perf_Com_train <- performance(pred_Com_train, "tpr", "fpr")
  plot(perf_Com_train, lwd = 2.5, colorize = FALSE, col = "#2258aa")
  lines(c(0, 1), c(0, 1), lwd = 2, col = "gray", lty = 4 )
  text(0.7, 0.3, cex = c(2, 2), labels = paste0("Training cohort = ", round(auc_Com_train, 3)), col = "#2258aa")
  
  Y_Compredict_validate <- predict(fitCom, validate1, type = "response")
  KZ_Com_validate <- cbind(validate$Diagnosis, Y_Compredict_validate) 
  pred_Com_validate <- ROCR::prediction(KZ_Com_validate[, 2], KZ_Com_validate[, 1])
  auc_Com_validate <- performance(pred_Com_validate, "auc")@y.values[[1]]
  perf_Com_validate <- performance(pred_Com_validate, "tpr", "fpr")
  plot(perf_Com_validate, lwd = 2.5, colorize = FALSE, col = "#ff0a0a", add = T)
  # lines(c(0, 1), c(0, 1), lwd = 2, col = "gray", lty = 4 )
  text(0.7, 0.22, cex = c(2, 2), labels = paste0("Internal validation cohort = ", round(auc_Com_validate, 3)), col = "#ff0a0a")
  
  Y_Compredict_test1 <- predict(fitCom, validate2, type = "response")
  KZ_Com_test1 <- cbind(validate2$Diagnosis, Y_Compredict_test1) 
  pred_Com_test1 <- ROCR::prediction(KZ_Com_test1[, 2], KZ_Com_test1[, 1])
  auc_Com_test1 <- performance(pred_Com_test1, "auc")@y.values[[1]]
  perf_Com_test1 <- performance(pred_Com_test1, "tpr", "fpr")
  plot(perf_Com_test1, lwd = 2.5, colorize = FALSE, col = "#eea00c", add = T)
  # lines(c(0, 1), c(0, 1), lwd = 2, col = "gray", lty = 4 )
  text(0.7, 0.14, cex = c(2, 2), labels = paste0("External validation cohort-1 = ", round(auc_Com_test1, 3)), col = "#eea00c")
  
  Y_Compredict_test2 <- predict(fitCom, validate3, type = "response")
  KZ_Com_test2 <- cbind(validate3$Diagnosis, Y_Compredict_test2) 
  pred_Com_test2 <- ROCR::prediction(KZ_Com_test2[, 2], KZ_Com_test2[, 1])
  auc_Com_test2 <- performance(pred_Com_test2, "auc")@y.values[[1]]
  perf_Com_test2 <- performance(pred_Com_test2, "tpr", "fpr")
  plot(perf_Com_test2, lwd = 2.5, colorize = FALSE, col = "#059e45", add = T)
  # lines(c(0, 1), c(0, 1), lwd = 2, col = "gray", lty = 4 )
  text(0.7, 0.06, cex = c(2, 2), labels = paste0("External validation cohort-2 = ", round(auc_Com_test2, 3)), col = "#059e45")
}

### 9.NRI,IDI------------------------------------------------------------------
{
  z.simple = as.matrix(subset(DCAdata, select = c(Eosinophil)))
  z.complex = as.matrix(subset(DCAdata, select = c(Eosinophil, Radscore)))
  
  msimple = glm(DCAY~., binomial(logit), data.frame(z.simple), x = T)
  mcomplex = glm(DCAY~., binomial(logit), data.frame(z.complex), x = T)
  
  pCli <- msimple$fitted.values
  pCom <- mcomplex$fitted.values
  reclassification(data = DCAdata,
                   cOutcome = 1,
                   predrisk1 = pCli, 
                   predrisk2 = pCom,
                   cutoff = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))
}
{
# Clicutoff
# Comcutoff
  NRI <- data.frame(Cli = pCli, Com = pCom, 
                    event = as.factor(data$Diagnosis))
  p.NRI <- ggplot(data = NRI, aes(x = Cli, y = Com, color = event)) +
    geom_point(size = 2) +
    scale_discrete_manual(values=c("#2258AA","#FF0A0A"),
                          aesthetics = 'colour', labels = c("NON", 'EOS')) +
    labs(x = "Clinical Model Predict", y = "Combined Model Predict") +
    xlim(0, 1) +
    ylim(0, 1) +
    theme_bw() +
    geom_hline(aes(yintercept = Comcutoff2), linetype = "dashed") +
    geom_vline(aes(xintercept = Clicutoff2), linetype = "dashed")
  p.NRI 
}
{
  IDI <- data.frame(event = as.factor(data$Diagnosis), 
                    Cli = pCli, 
                    Com = pCom,
                    id = 1:dim(data)[1])
  IDI_Eos <- subset(IDI, subset = (event == 1)) 
  Cli_Eos <- mean(IDI_Eos$Cli)
  Com_Eos <- mean(IDI_Eos$Com)
  IDI_Non <- subset(IDI, subset = (event == 0)) 
  Cli_Non <- mean(IDI_Non$Cli)
  Com_Non <- mean(IDI_Non$Com)
  p.IDI <- ggplot(data = IDI) +
    geom_point(aes(x = reorder(id, Cli), y = Cli, color = event),
               size = 2, shape = 17) +
    # geom_point(aes(x = reorder(id, Com), y = Com, color = event),
    #                size = 2) +
    scale_discrete_manual(values = c("#2258AA","#FF0A0A"),
                          aesthetics = 'colour',
                          labels = c("NON", "EOS")) +
    labs(x = "Patient", y = "Model Predict") +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    geom_hline(aes(yintercept = Cli_Eos), linetype = "dashed") +
    geom_hline(aes(yintercept = Cli_Non), linetype = "dashed")
  # geom_hline(aes(yintercept = Com_Eos), linetype = "dashed") +
  # geom_hline(aes(yintercept = Com_Non), linetype = "dashed")
  p.IDI
}

### 13.Heatmap------------------------------------------------------------------
data <- read_excel("E:/Train_Heatmap.xls", sheet = "Sheet1")
data1 <- data[, 6:dim(data)[2]]
mydata_cor <- cor(data1, method = "spearman") 
data2 <- mydata_cor[19:48,1:18]
mydata_cor2 <- rcorr(as.matrix(data1))
pvalue <- mydata_cor2$P[19:48,1:18]

data_mark = as.matrix(pvalue) 
for(i in 1:dim(pvalue)[1]){    
  for(j in 1:dim(pvalue)[2]){        
    if(pvalue[i,j] <= 0.001 && pvalue[i,j] > 0){data_mark[i,j]="***"}            
    else if(pvalue[i,j] <= 0.01 && pvalue[i,j] > 0.001){data_mark[i,j]="**"} 
    else if(pvalue[i,j] <= 0.05 && pvalue[i,j] > 0.01){data_mark[i,j]="*"}            
    else {data_mark[i,j] = " "}}}
# * 0.05>=p>0.01; ** 0.01>=p>0.001; *** 0.001>=p

corrplot(data2) 
pheatmap(data2,
         cluster_rows = T,
         cluster_cols = T,
         border_color = "#7C858D",
         display_numbers = data_mark, 
         color = colorRampPalette(c("#0F3BA9", "white", "#FF0A0A"))(100)
)


p <- data2 %>% data.frame()
phr <- hclust(dist(p)) %>% ggtree(layout = "rectangular", branch.length="none")
phc <- hclust(dist(t(p))) %>% ggtree() + layout_dendrogram()
p$Rad <- rownames(p)
p1 <- gather(p, 1:18, key = "Cli", value = 'correlation')
pp <- ggplot(p1, aes(x = Cli, y = Rad)) +
  geom_tile(aes(fill = correlation)) +
  theme_minimal() +
  scale_y_discrete(position = "right") +
  xlab(NULL) +
  ylab(NULL) +
  scale_fill_gradient2(low = "#0532C7", high = "red", mid = "white", midpoint = 0.01) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
pp %>%
  insert_left(phr, width = .2) %>%
  insert_top(phc, height = .2)

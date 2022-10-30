####################################################################################################################
### title: "Center patients"
### author: "Kzzhu"
### time: "2022.06.04"
####################################################################################################################
### 1.Packages---------------------------------------------------------------
{
  suppressMessages({
    library(readxl)
    library(stats)
    })
  rm(list = ls())
  options(stringsAsFactors = F)
  work_dir <- "E:/Radiomics"
  setwd(work_dir)
}

### 2.Readcsv----------------------------------------------------------------
{
  NCRSwNP <- read.table("E:/Non.csv", header = T, sep = ",")
  ECRSwNP <- read.table("E:/Eos.csv", header = T, sep = ",")
  EOS <- ECRSwNP[, 1:dim(ECRSwNP)[2]]
  NON <- NCRSwNP[, 1:dim(NCRSwNP)[2]]
}

### 3.Group------------------------------------------------------------------
set.seed(n)
{
  ind1 <- sample(2, nrow(EOS), replace = T, prob = c(0.8, 0.2))
  train_eos <- EOS[ind1 == 1, ]
  validate_eos <- EOS[ind1 == 2, ]
  
  ind2 <- sample(2, nrow(NON), replace = T, prob = c(0.8, 0.2))
  train_non <- NON2[ind2 == 1, ]
  validate_non <- NON2[ind2 == 2, ]
  
  train <- rbind(train_eos, train_non)
  validate <- rbind(validate_eos, validate_non)
}
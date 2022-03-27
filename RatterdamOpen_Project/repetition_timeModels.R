# Repetition Models
# WH Oct 24 2021
# 2022-03-23 revamping this script for repetition manuscript
# Running models to test whether there's any effect of time
# without commenting here on what that time signal looks like 
# Time interaction models 


library(lme4)
library(ggplot2)
library(lmtest)
library(stringr)
library(emmeans)
library(tidyr)
library(ggpubr)
library(splines)

# Load data and recode factors
alleypath <- "E:\\Ratterdam\\R_data_repetition\\2022-03-23_AlleySuperpopDirVisitFiltered.csv"
alleydf <- read.csv(alleypath,header=TRUE)

alleydf$PrevDir <- as.factor(alleydf$PrevDir)
alleydf$CurrDir <- as.factor(alleydf$CurrDir)
alleydf$NextDir <- as.factor(alleydf$NextDir)
alleydf$RetroEgo <- as.factor(alleydf$RetroEgo)
alleydf$ProspEgo <- as.factor(alleydf$ProspEgo)
alleydf$Repeating <- as.factor(alleydf$Repeating)
alleydf$Traversal <- as.factor(alleydf$Traversal)
alleydf$CellID <- as.factor(alleydf$CellID)
alleydf$FieldNum <- as.factor(alleydf$FieldNum)
alleydf$FieldID <- as.factor(alleydf$FieldID)
alleydf$Alleys <- as.factor(alleydf$Alleys)


alleydf = alleydf[alleydf$Traversal=='True',]


###
### Time
###
savepath <- 'E:\\Ratterdam\\temp\\TimeCD_Model\\'

sigs <- c()
sigreps <- c()
signonreps <- c()
run <- 0
nreps <- 0
nnonreps <- 0

repOrNot <- c()
sigs <- c()
orientation <- c()

alpha  = 0.05/3

rmse <- c()

rmse_base <- c()
rmse_alt <- c()

set.seed(123)
for(o in c("V","H")){
  oriendf <- subset(alleydf, Orientation==o)
  
  for(fid in unique(oriendf$FieldID)){
    u<-try({
      run <- run + 1
      
      field <- subset(oriendf, FieldID==fid)
      basemod <- glm(Rate + 1 ~ CurrDir,family='Gamma',data=field)
      mod <- glm(Rate + 1 ~ CurrDir + ns(StartTimes,3),family='Gamma',data=field)
      s<-summary(mod)
      anysig <- FALSE
      
      if(unique(field$Repeating)=="True"){
        nreps <- nreps + 1
      }
      else if(unique(field$Repeating)=="False"){
        nnonreps <- nnonreps + 1
      }
      
      s <- summary(mod)
      
      p1 <- s$coefficients["ns(StartTimes, 3)1","Pr(>|t|)"]
      p2 <- s$coefficients["ns(StartTimes, 3)2","Pr(>|t|)"]
      p3 <- s$coefficients["ns(StartTimes, 3)3","Pr(>|t|)"]
      
      # save labels
      
      orientation <- c(orientation, o)
      
      if((p1<alpha)|(p2<alpha)|(p3<alpha)){
        
        sigs <- c(sigs, TRUE)
        
      }
      else{
        sigs <- c(sigs, FALSE)
      }
      
      if(unique(field$Repeating)=="True"){
        repOrNot <- c(repOrNot, TRUE)
      }
      else if(unique(field$Repeating)=="False"){
        repOrNot <- c(repOrNot, FALSE)
      }
      
      rmse_base <- c(rmse_base, sqrt(mean((field$Rate-basemod$fitted.values)^2)))
      rmse_alt <- c(rmse_alt, sqrt(mean((field$Rate-mod$fitted.values)^2)))
      
      
      
      
    },silent=FALSE)
  }
  
}



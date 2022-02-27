# Repetition Project WH 2022-02-25
# Independence of Temporal Dynamics Among Repeating Fields
# Test model with interaction between field ID and time 

library(lme4)
library(ggplot2)
library(lmtest)
library(stringr)
library(emmeans)
library(tidyr)
library(ggpubr)
library(splines)
library(car)

# 211210 has 30% field overlap threshold and slightly looser traversal thresholds 
#alleypath <- "E:\\Ratterdam\\R_data_repetition\\20220215-140515_superPopAlleyBehaviorResponse_1.5vfilt_PosInFieldNormedFR.csv"
allleypath <- "E:\\Ratterdam\\R_data_repetition\\220222_AlleySuperpopDirVisitFiltered.csv"

alleydf <- read.csv(alleypath,header=TRUE)


alleydf <- alleydf[is.finite(alleydf$Rate),]
alleydf <- alleydf[alleydf$Traversal=="True",]
alleydf <- alleydf[alleydf$Reward=="False",]


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
alleydf$NumFields <- as.numeric(alleydf$NumFields)

startTimeKnots <- 3

alpha <- 0.05

total_models_run <- 0

sigInts <- c()

for(o in c('V','H')){
  oriendf <- subset(alleydf, Orientation==o)
  
  for(cid in unique(oriendf$CellID)){
    cell <- subset(oriendf, CellID==cid)
    if(unique(cell$Repeating)=='True'){
      
      
      try({
      
      m <- glm(Rate + 1 ~ FieldID:ns(StartTimes,3),data=cell,family='Gamma')
      s <- summary(m)
      if(sum(tail(s$coefficients[,'Pr(>|t|)']<0.05,-1))>=1){
        sigInts <- c(sigInts, cid)
      }
      total_models_run <- total_models_run + 1
      
      },silent=FALSE)
    }
    
    
    
  }
  
}
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
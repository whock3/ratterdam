# GLMs to test egocentric vs allocentric tuning within single fields

library(lme4)
library(ggplot2)
library(lmtest)
library(stringr)
library(emmeans)
library(tidyr)
library(ggpubr)
library(splines)

# Load data and recode factors
alleypath <- "E:\\Ratterdam\\R_data_repetition\\211005_AlleySuperpopDirVisitFiltered.csv"

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

for(o in c('V','H')){
  oriendf <- subset(alleydf, Orientation==o)
  
  for(fid in unique(oriendf$FieldID)){
    
    field_ <- subset(oriendf, FieldID == fid)
    
    n <- sample(nrow(field_))
    shuff.field <- data.frame(field_)
    shuff.field$CurrDir <- shuff.field$CurrDir[n]
    shuff.field$PrevDir <- shuff.field$PrevDir[n]
    shuff.field$NextDir <- shuff.field$NextDir[n]
    
    if(shuffle==TRUE){
      field <- shuff.field
      
    }
    else if(shuffle==FALSE){
      field <- field_
    }
    
    
    mrun <- try({

      m1 <- glm(Rate+1 ~ ns(StartTimes,startTimeKnots)+CurrDir+PrevDir+NextDir+ProspEgo+RetroEgo,
               family='Gamma', data=field)
      
      s <- summary(m1)
      print(s$coefficients)
      
      
      
    },silent=TRUE)
    
  }
  
}

      
# Repetition Project - Detecting and Analyzing Multicollinearity
# WH 2022-01-13

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
alleypath <- "E:\\Ratterdam\\R_data_repetition\\211220_AlleySuperpopDirVisitFiltered.csv"

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
alleydf$NumFields <- as.numeric(alleydf$NumFields)


vifs <- data.frame(matrix(nrow=0,ncol=3))
colnames(vifs) <- cols



m <- glm(Rate+1 ~ CurrDir+PrevDir+NextDir,
family='Gamma',
data=field)
al <- alias(m)
if("Complete" %in% names(al)){
  
  field2 <- data.frame(field) # copy
  
  alidx <- which(al$Complete!=0, arr.ind=T)
  
  #this is not R-thonic but it works
  for(i in dim(alidx)[1]){
  
    #there are, by definition, two factor levels that are collinear for every
    # instance of multicollinearity 
    mc_A <- row.names(al$Complete)[alidx[i,1]]
    mc_B <- colnames(al$Complete)[alidx[i,2]]
    
    for(mc in c(mc_A, mc_B)){
      
      factor_level <- substr(dl,nchar(dl),nchar(dl))
      factor_type <- substr(dl,1,nchar(dl)-1)
      
      if(factor_type == 'NextDir'){
        field2 <- droplevels(subset(field, NextDir!=factor_level))
      }
      else if(factor_type == 'PrevDir'){
        field2 <- droplevels(subset(field, PrevDir!=factor_level))
        
      }
      else if(factor_type == 'CurrDir'){
        field2 <- droplevels(subset(field, CurrDir!=factor_level))
      }
    }
    
  }

}


total_unaliased <- 0
total_aliased <- 0

vif_curr <- c()
vif_prev <- c()
vif_next <- c()

for(o in c('V','H')){
  oriendf <- subset(alleydf, Orientation==o)
  
  for(fid in unique(oriendf$FieldID)){
    print(fid)
    
    try({
    field <- subset(oriendf, FieldID == fid)
    
    m <- glm(Rate + 1 ~ CurrDir + PrevDir + NextDir,family='Gamma',data=field)
    al <- alias(m)
    print(al)
    if(!("Complete" %in% names(al))){
      
      total_unaliased <- total_unaliased + 1
      
      v <- vif(m) # will error out if there are any cases of complete MC
      
      #sometimes vif() returns a single value for each var,
      # sometimes 2 values (gvif and gvif^(1/(2*df)))
      if(length(v)==9){
        cvif <- v['CurrDir',3]
        pvif <- v['PrevDir',3]
        nvif <- v['NextDir',3]
      }
      else if(length(v)==3){
        cvif <- v['CurrDir']
        pvif <- v['PrevDir']
        nvif <- v['NextDir']
      }
      
      vif_curr <- c(vif_curr, cvif)
      vif_prev <- c(vif_prev, pvif)
      vif_next <- c(vif_next, nvif)
      
    }
    else if("Complete" %in% names(al)){
      
      total_aliased <- total_aliased + 1
      
    }
    
    },silent=TRUE)
    
  }
  
}





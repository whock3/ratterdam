# 22-02-22 Current Direction Model Comparison
# WH Repetition project
# GLM model comparison testing effect of current direction
# Compare base model (time + intercept) to model with current direction as well
# This analysis was contained in repetition_GLM_ModelComparison, but there
# were many other comparions being made and the script got confusing. Just CD here. 



library(lme4)
library(ggplot2)
library(lmtest)
library(stringr)
library(emmeans)
library(tidyr)
library(ggpubr)
library(splines)


alleypath <- "E:\\Ratterdam\\R_data_repetition\\220222_AlleySuperpopDirVisitFiltered.csv"

alleydf <- read.csv(alleypath,header=TRUE)

alleydf$PrevDir <- as.factor(alleydf$PrevDir)
alleydf$CurrDir <- as.factor(alleydf$CurrDir)
alleydf$NextDir <- as.factor(alleydf$NextDir)
alleydf$RetroEgo <- as.factor(alleydf$RetroEgo)
alleydf$ProspEgo <- as.factor(alleydf$ProspEgo)
alleydf$Repeating <- as.factor(alleydf$Repeating)
alleydf$Traversal <- as.factor(alleydf$Traversal)
alleydf$Reward <- as.factor(alleydf$Reward)
alleydf$CellID <- as.factor(alleydf$CellID)
alleydf$FieldNum <- as.factor(alleydf$FieldNum)
alleydf$FieldID <- as.factor(alleydf$FieldID)
alleydf$Alleys <- as.factor(alleydf$Alleys)
alleydf$NumFields <- as.numeric(alleydf$NumFields)

m1_rmse <- c()
m2_rmse <- c()

m1_aic <- c()
m2_aic <- c()

lrCurr_pvals <- c() # keep all pvalues from lrtest(base, base+cd) regardless
# if thats the best model because want to plot rmse colored
# by whether cd was helpful for fig 3 sfn2021. fig 5 then gets
# into whats the best model. 

startTimeKnots = 3

sigP <- c() # save whether CD is sig, for figures

repOrNot <- c()

for(o in c('V','H')){
  oriendf <- subset(alleydf, Orientation==o)
  
  for(fid in unique(oriendf$FieldID)){
    
    field <- subset(oriendf, FieldID == fid)
    
    mrun <- try({
       
      m1 <- glm(Rate+1 ~ ns(StartTimes,startTimeKnots), family='Gamma', data=field)
      m2 <- glm(Rate+1 ~ CurrDir + ns(StartTimes,startTimeKnots), family = 'Gamma', data=field)
     
      m1_rmse <- c(m1_rmse, sqrt(mean((field$Rate-m1$fitted.values)^2)))
      m2_rmse <- c(m2_rmse, sqrt(mean((field$Rate-m2$fitted.values)^2)))
      
      m1_aic <- c(m1_aic, m1$aic)
      m2_aic <- c(m2_aic, m2$aic)
     
      lrCurr <- lrtest(m1,m2)
      
      lrCurr_pvals <- c(lrCurr_pvals, lrCurr[2,"Pr(>Chisq)"])
      
      if((lrCurr[2,"Pr(>Chisq)"] < 0.05)&(!is.nan(lrCurr[2,"Pr(>Chisq)"]))){
        
        sigP <- c(sigP, 1)
      }
      else{
        sigP <- c(sigP, 0)
      }
      
      if(unique(field$Repeating)=='True'){
        
        repOrNot <- c(repOrNot,TRUE)
      }
      else if(unique(field$Repeating)=='False'){
        repOrNot <- c(repOrNot, FALSE)
      }

      
      
    },silent=FALSE) 
  }
}




















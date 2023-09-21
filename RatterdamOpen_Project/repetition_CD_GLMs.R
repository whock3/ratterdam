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


alleypath <- "E:\\Ratterdam\\R_data_repetition\\2022-04-05_AlleySuperpopDirVisitFiltered.csv"

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

failed_units <- c()

fids <- c() # save fieldIDs to cross reference with other methods to see overlap in which are significant

lrCurr_pvals <- c() # keep all pvalues from lrtest(base, base+cd) regardless
# if thats the best model because want to plot rmse colored
# by whether cd was helpful for fig 3 sfn2021. fig 5 then gets
# into whats the best model. 

startTimeKnots = 3

shuffle <- TRUE
nshuffles <- 1000


sigP <- c() # save whether CD is sig, for figures

repOrNot <- c()

total_current_responsive <- c()


for(s in 1:nshuffles){
  
  print(s)

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
      
      mrun <- tryCatch({
         
        m1 <- glm(Rate+1 ~ ns(StartTimes,startTimeKnots), family='Gamma', data=field)
        m2 <- glm(Rate+1 ~ CurrDir + ns(StartTimes,startTimeKnots), family = 'Gamma', data=field)
       
        m1_rmse <- c(m1_rmse, sqrt(mean((field$Rate-m1$fitted.values)^2)))
        m2_rmse <- c(m2_rmse, sqrt(mean((field$Rate-m2$fitted.values)^2)))
        
        fids <- c(fids, fid)
        
        m1_aic <- c(m1_aic, m1$aic)
        m2_aic <- c(m2_aic, m2$aic)
       
        lrCurr <- lrtest(m1,m2)
        
        lrCurr_pvals <- c(lrCurr_pvals, lrCurr[2,"Pr(>Chisq)"])
        
        #dropping pvalue to 0.05/2 2022-04-14
        if((lrCurr[2,"Pr(>Chisq)"] < 0.05/2)&(!is.nan(lrCurr[2,"Pr(>Chisq)"]))){
          
          sigP <- c(sigP, 1)
        }
        else{
          sigP <- c(sigP, 0)
        }
        
        if(unique(field$Repeating)=='TRUE'){
          
          repOrNot <- c(repOrNot,TRUE)
        }
        else if(unique(field$Repeating)=='FALSE'){
          repOrNot <- c(repOrNot, FALSE)
        }
  
        
        
      },
      error=function(e){
        failed_units <<- c(failed_units, c(fid, o))
        
      }
      ,silent=FALSE) 
      
    
    }
  }
  
  total_current_responsive <- c(total_current_responsive, sum(sigP)/length(sigP))
  
}



# debugging 4/11/22 why glm analysis and MW test have different sample sizes
# seems to be model complexity in GLM causes some fields to fail out of test

# field = subset(alleydf, FieldID==47 & Orientation=='H')
# 
# plot(field$StartTimes, field$Rate)
# m1 <- glm(Rate+1 ~ ns(StartTimes,3), family='Gamma', data=field)
# m2 <- glm(Rate+1 ~ CurrDir, family = 'Gamma', data=field)
# m3 <- glm(Rate+1 ~ CurrDir + ns(StartTimes,3), family = 'Gamma', data=field)
# 
# 


# > failed_units
# "14"  "V"   
# "125" "V"   
# "166" "V"   
# "278" "V"   
# "301" "V"   
# "310" "V"   
# "47"  "H"   
# "122" "H"   
# "140" "H"   
# "143" "H"   
# "154" "H"   
# "169" "H"   
# "286" "H"   
# "321" "H"   
# "323" "H"   
# "347" "H"  
 











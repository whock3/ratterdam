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
alleypath <- "E:\\Ratterdam\\R_data_repetition\\2022-04-05_AlleySuperpopDirVisitFiltered.csv"
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


###
### Time
###
savepath <- 'E:\\Ratterdam\\repetition_manuscript\\Figure6_TemporalDynamics\\2022-04-20_timeGLMResults.csv'

sigs <- c()
sigreps <- c()
signonreps <- c()
run <- 0
nreps <- 0
nnonreps <- 0

repOrNot <- c()
sigs <- c()
orientation <- c()

alpha  = 0.05/6 # 6 = 3 * 2. 3 for time knots, 2 for orientation filtering. 

rmse <- c()

rmse_base <- c()
rmse_alt <- c()

shuffle <- TRUE
nshuffles <- 1000

set.seed(123)

total_time_responsive <- c()

for(s in 1:nshuffles){
  
  print(s)
  for(o in c("V","H")){
    oriendf <- subset(alleydf, Orientation==o)
    
    for(fid in unique(oriendf$FieldID)){
      u<-try({
        run <- run + 1
        
        # shuffle or not 
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
          
          sigs <- c(sigs, 1)
          
        }
        else{
          sigs <- c(sigs, 0)
        }
        
        if(unique(field$Repeating)==TRUE){
          repOrNot <- c(repOrNot, TRUE)
        }
        else if(unique(field$Repeating)==FALSE){
          repOrNot <- c(repOrNot, FALSE)
        }
        
        rmse_base <- c(rmse_base, sqrt(mean((field$Rate-basemod$fitted.values)^2)))
        rmse_alt <- c(rmse_alt, sqrt(mean((field$Rate-mod$fitted.values)^2)))
        
        
        
        
      },silent=FALSE)
    }
    
  }
  
  total_time_responsive <- c(total_time_responsive, sum(sigs)/length(sigs))
  
}



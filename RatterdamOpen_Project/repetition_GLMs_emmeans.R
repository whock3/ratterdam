# Repetition Project
# 22-01-08 WH 
# Script to create GLM per field and test effects of factors 
# via post-hoc comparisons between levels within factor, via emmeans
# Time will be modeled as a FE or RE (using glmer) depending on further thought
#

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
alleypath <- "E:\\Ratterdam\\R_data_repetition\\20220215-140515_superPopAlleyBehaviorResponse_1.5vfilt_PosInFieldNormedFR.csv"


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
total_current <- 0 
total_next <- 0
total_previous <- 0

current_responsive <- c()
next_responsive <- c()
previous_responsive <- c()

vif_thresh = 5


for(o in c('V','H')){
  oriendf <- subset(alleydf, Orientation==o)
  
  for(fid in unique(oriendf$FieldID)){
    print(fid)
    
    field <- subset(oriendf, FieldID == fid)
    try({
    m <- glm(Rate + 1 ~ CurrDir + PrevDir + NextDir + ns(StartTimes,3),
             family='Gamma',
             data=field)
    
    al <- alias(m)
    if(!("Complete" %in% names(al))){
    
    v <- vif(m) # will error out if there are any cases of complete MC
    
    #sometimes vif() returns a single value for each var,
    # sometimes 2 values (gvif and gvif^(1/(2*df)))
    if(length(v)==12){
      cvif <- v['CurrDir',3]
      pvif <- v['PrevDir',3]
      nvif <- v['NextDir',3]
    }
    else if(length(v)==4){
      cvif <- v['CurrDir']
      pvif <- v['PrevDir']
      nvif <- v['NextDir']
    }
    
    total_models_run <- total_models_run + 1
    
    #Previous Direction
    if(pvif < vif_thresh){
      total_previous <- total_previous + 1
      em_m <- emmeans(m,"PrevDir")
      pwc <- summary(pairs(em_m))
      sig <- pwc[1][pwc[6]<alpha]
      if(length(sig)>=1){
        previous_responsive <- c(previous_responsive, fid)
      }
    }
    
    #Current Direction
    if(cvif < vif_thresh){
    total_current <- total_current + 1
    em_m <- emmeans(m,"CurrDir")
    pwc <- summary(pairs(em_m))
    sig <- pwc[1][pwc[6]<alpha]
    if(length(sig)>=1){
      current_responsive <- c(current_responsive, fid)
    }
    }

    

    #Next Direction
    if(nvif < vif_thresh){
    total_next <- total_next + 1
    em_m <- emmeans(m,"NextDir")
    pwc = summary(pairs(em_m))
    sig <- pwc[1][pwc[6]<alpha]
    if(length(sig)>=1){
      next_responsive <- c(next_responsive, fid)
    }
    }
    
    
    
    }
    
    },silent=TRUE)
    
  }
  
}

print(length(current_responsive)/total_current)
print(length(previous_responsive)/total_previous)
print(length(next_responsive)/total_next)



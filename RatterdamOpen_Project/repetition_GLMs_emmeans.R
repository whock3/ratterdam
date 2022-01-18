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

startTimeKnots <- 3

alpha <- 0.05

total_models_run <- 0
total_current <- 0 
total_next <- 0
total_previous <- 0

current_responsive <- c()
next_responsive <- c()
previous_responsive <- c()

for(o in c('V','H')){
  oriendf <- subset(alleydf, Orientation==o)
  
  for(fid in unique(oriendf$FieldID)){
    print(fid)
    
    field <- subset(oriendf, FieldID == fid)
    
    try({
      rm(em_m, pwc, sig)
    })
    
    try({
    total_models_run <- total_models_run + 1
    m <- glm(Rate + 1 ~ CurrDir + PrevDir + NextDir,family='Gamma',data=field)
    
    try({
    #Current Direction
    total_current <- total_current + 1
    em_m <- emmeans(m,"CurrDir")
    pwc <- summary(pairs(em_m))
    sig <- pwc[1][pwc[6]<alpha]
    if(length(sig)>=1){
      current_responsive <- c(current_responsive, fid)
    }
    },silent=FALSE)
    
    try({
    #Next Direction
    total_next <- total_next + 1
    em_m <- emmeans(m,"NextDir")
    pwc = summary(pairs(em_m))
    sig <- pwc[1][pwc[6]<alpha]
    if(length(sig)>=1){
      next_responsive <- c(next_responsive, fid)
    }
    },silent=FALSE)
    
    try({
    #Previous Direction
    total_previous <- total_previous + 1
    em_m <- emmeans(m,"PrevDir")
    pwc <- summary(pairs(em_m))
    sig <- pwc[1][pwc[6]<alpha]
    if(length(sig)>=1){
      previous_responsive <- c(previous_responsive, fid)
    }
    },silent=FALSE)
    
    },silent=FALSE)
  
    
    
  }
  
}


# Shuffling LMER Analysis
# WH Mid March 2021
# Run LMER analysis on data with shuffled texture lables within an alley
# Compare distribution of effects (e.g. coefficients) to empirical

# imports
library(lme4)
library(splines)
library(ggplot2)
library(lmtest)
library(stringr)


setwd("E:\\UserData\\Documents\\GitHub\\ratterdam\\")
source("E:\\UserData\\Documents\\GitHub\\ratterdam\\Beltway_Project\\cicheck.R")
source("E:\\UserData\\Documents\\GitHub\\ratterdam\\Beltway_Project\\glmer_fx.R")
source("E:\\UserData\\Documents\\GitHub\\ratterdam\\Beltway_Project\\ratterdam_RunLMER.R")


# load cell and define key variables

cellname <- "TT1\\cl-maze1.6"
alley=11

celldf <- df[which(df$name==cellname & df$alley==alley),]
celldf_copy <- data.frame(celldf)

ntrials <- length(unique(celldf$trial))
nshuffles <- 1000

# Create trial list

trialList <- vector(mode="numeric", length=ntrials)

for (i in 0:ntrials-1){
  
  trialList[i+1] = unique(df[df$trial==i,]$texture)
  
}


# shuffle loop

for (i in 1:nshuffles){
  
  
  shuffledTrials <- trialList[sample(1:ntrials)]
  
  
  for (j in 1:ntrials){
    
  celldf_copy[celldf$trial==j-1,]$texture <- shuffledTrials[j]
  
      
  }
  
  celldf_copy <- lmer_routine(celldf_copy)
  
  any_nonoverlap <- FALSE
  
  for (a in unique(celldf$alley)){
    
    nonoverlap <- checkAlleyNonoverlap(celldf, a)
    if (nonoverlap == TRUE){
      
      any_nonoverlap <- TRUE
    }
      
    
  }
  
  
  
  
  
}










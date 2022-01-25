# Repetition Project
# GLMs with posthoc emmeans contrasts using dummified trajectory codes 
# where each three-label route (previous, current, next dirs) is represented
# as a unique categorical label. E.g. "NNE"



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


codes <- c()
for(i in seq(1:length(1:nrow(alleydf))))
  {codes <- c(codes, sprintf("%s%s%s",alleydf[i,"PrevDir"],alleydf[i,"CurrDir"],alleydf[i,"NextDir"]))}

alleydf$Code <- as.factor(codes)

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

alpha <- 0.05  # emmeans pairs does a multiple comparisons correction internally


traj_responsive <- c()
total_models_run <- 0
rep_coding <- c()
sf_coding <- c()

for(o in c('V','H')){
  oriendf <- subset(alleydf, Orientation==o)
  
  for(fid in unique(oriendf$FieldID)){
    print(fid)
    
    field <- subset(oriendf, FieldID == fid)
    try({
      total_models_run <- total_models_run + 1
      m <- glm(Rate + 1 ~ Code + ns(StartTimes, 3),family='Gamma',data=field)
      em_m <- emmeans(m,"Code")
      pwc <- summary(pairs(em_m))
      sig <- pwc[1][pwc[6]<alpha]
      if(length(sig)>=1){
        traj_responsive <- c(traj_responsive, c(c(fid,o)))
        if(unique(field$Repeating)=="True"){
          rep_coding <- c(rep_coding, TRUE)
        }
        else if(unique(field$NumFields)==1){
          sf_coding <- c(sf_coding, TRUE)
        }
      }
      else{
        if(unique(field$Repeating)=="True"){
          rep_coding <- c(rep_coding, FALSE)
        }
        else if(unique(field$NumFields)==1){
          sf_coding <- c(sf_coding, FALSE)
        }
      }
   
      
    },silent=FALSE)
    
  }
  
}
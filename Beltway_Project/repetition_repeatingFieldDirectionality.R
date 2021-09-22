#Repetition Project Will Hockeimer 9/22/21
# Looking at directionality across repeating fields, specifically if they
# share the same directionality 

library(lme4)
library(ggplot2)
library(lmtest)
library(stringr)
library(emmeans)
library(tidyr)
library(ggpubr)
library(splines)

# Load data and recode factors
alleypath <- "E:\\Ratterdam\\R_data_repetition\\20210921-183738_superPopAlleyBehaviorResponse_1.5vfilt.csv"
interpath <- "E:\\Ratterdam\\R_data_repetition\\20210921-183739_superPopInterBehaviorResponse_1.5vfilt.csv"

alleydf <- read.csv(alleypath,header=TRUE)
interdf <- read.csv(interpath, header=TRUE)

alleydf$PrevDir <- as.factor(alleydf$PrevDir)
alleydf$CurrDir <- as.factor(alleydf$CurrDir)
alleydf$NextDir <- as.factor(alleydf$NextDir)
alleydf$RetroEgo <- as.factor(alleydf$RetroEgo)
alleydf$ProspEgo <- as.factor(alleydf$ProspEgo)
alleydf$Repeating <- as.factor(alleydf$Repeating)

interdf$PrevDir <- as.factor(interdf$PrevDir)
interdf$CurrentEgo <- as.factor(interdf$CurrentEgo)
interdf$NextDir <- as.factor(interdf$NextDir)
interdf$RetroEgo <- as.factor(interdf$RetroEgo)
interdf$ProspEgo <- as.factor(interdf$ProspEgo)
interdf$Repeating <- as.factor(interdf$Repeating)


# First and simplest way: two way anova, current direction interacting field 
horizsig <- c()
horiztotal <- 0
vertsig <- c()
verttotal <- 0

for(cid in unique(alleydf$CellID)){
  
  fdf <- subset(alleydf, CellID==cid)
  
  #check if its repeating and if multiple fields have data. this isnt always the
  # case if, e.g. all the passes were turnarounds or invalid due to bad sampling of field
  if((unique(fdf$Repeating)=='True') & length(unique(fdf$FieldNum))>1){
   
     #need a try block because sometimes you dont get a valid interaction, e.g. if
    # all visits to a given field are 0 it wont estimate interaction (?)
    
    # Do vertical alleys
    try({
      vertfdf <- subset(fdf, Orientation=='V')
      mod <- aov(Rate~CurrDir*FieldNum, data=vertfdf)
      
      # In R you apparently fit with aov (or lm) first, then call anova() to get results. summary() seems
      # to give same answer up to rounding. 
      # see: https://stackoverflow.com/questions/40823310/when-should-i-use-aov-and-when-anova
      p <- anova(mod)['CurrDir:FieldNum','Pr(>F)']
      print(p)
      if(p < 0.05/2){
        vertsig <- c(vertsig, cid)
      }
      
      verttotal <- verttotal + 1 #increment counter for number of cells actually analyzed for mult comp purposes 
      
    },silent=FALSE)
    
    #Do horizontal alleys
    try({
      horizfdf <- subset(fdf, Orientation=='H')
      mod <- aov(Rate~CurrDir*FieldNum, data=horizfdf)
      
      # In R you apparently fit with aov (or lm) first, then call anova() to get results. summary() seems
      # to give same answer up to rounding. 
      # see: https://stackoverflow.com/questions/40823310/when-should-i-use-aov-and-when-anova
      p <- anova(mod)['CurrDir:FieldNum','Pr(>F)']
      print(p)
      if(p < 0.05/2){
        horizsig <- c(horizsig, cid)
      }
      
      horiztotal <- horiztotal + 1 #increment counter for number of cells actually analyzed for mult comp purposes 
      
    },silent=FALSE)
  }
}

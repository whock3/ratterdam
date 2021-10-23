## (G)LMM Models and Model Comparisons For Testing Effect of Behavioral Variables
## On Individual Field Firing.
## WH Repetition Project 18-10-21

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



repFull <- 0
repCD <- 0
repPD <- 0
repND <- 0

nonrepFull <- 0
nonrepCD <- 0
nonrepPD <- 0
nonrepND <- 0

msigCD <- 0 # how many models have better fit w current dir, use to see if
            # num models with better fit of P+C+N is sig  by binomial 
actualRunRep <- 0
actualRunNonrep <- 0


o <- 'H'
oriendf <- subset(alleydf, Orientation==o)

for(fid in unique(oriendf$FieldID)){
  
  field <- subset(oriendf, FieldID == fid)
  mrun <- try({
    m1 <- glm(Rate+1 ~ 1, family='Gamma', data=field)
    m2 <- glm(Rate+1 ~ CurrDir, family = 'Gamma', data=field)
    m3 <- glm(Rate+1 ~ PrevDir + CurrDir, family='Gamma',data=field)
    m4 <- glm(Rate+1 ~ CurrDir + NextDir, family='Gamma',data=field)
    m5 <- glm(Rate+1 ~ PrevDir + CurrDir + NextDir, family='Gamma',data=field)
    
  
    #keep track of how many fields were tested by group
    if(unique(field$Repeating)=='True'){
      actualRunRep <- actualRunRep + 1
    }
    else if(unique(field$Repeating)=='False'){
      actualRunNonrep <- actualRunNonrep + 1
    }
    

      
    lrCurr <- lrtest(m1,m2)
    lrPrev <- lrtest(m2,m3)
    lrNext <- lrtest(m2,m4)
    lrFullA <- lrtest(m3,m5)
    lrFullB <- lrtest(m4,m5)
    bfAdj <- 5
    
    # Case checks for repeating fields
    if(unique(field$Repeating)=='True'){
      
      #check full model against partial C+N and partial P+C, see if its better than both
      if((lrFullA[2,"Pr(>Chisq)"]<(0.05/bfAdj))&(lrFullB[2,"Pr(>Chisq)"]<(0.05/bfAdj))){
        repFull <- repFull + 1
      }
      #otherwise drop down to test C vs C+P
      else if(lrPrev[2,"Pr(>Chisq)"]<(0.05/bfAdj)){
        repPD <- repPD + 1
      }
      #otherwise drop down to test C vs C+N
      else if(lrNext[2,"Pr(>Chisq)"]<(0.05/bfAdj)){
        repND <- repND + 1
      }
      #otherwise drop down to test null vs C
      else if(lrCurr[2,"Pr(>Chisq)"]<(0.05/bfAdj)){
        repCD <- repCD + 1
      }
    }
    
    
    # Case checks for nonrepeating fields
    if(unique(field$Repeating)=='False'){
      
      #check full model against partial C+N and partial P+C, see if its better than both
      if((lrFullA[2,"Pr(>Chisq)"]<(0.05/bfAdj))&(lrFullB[2,"Pr(>Chisq)"]<(0.05/bfAdj))){
        nonrepFull <- nonrepFull + 1
      }
      #otherwise drop down to test C vs C+P
      else if(lrPrev[2,"Pr(>Chisq)"]<(0.05/bfAdj)){
        nonrepPD <- nonrepPD + 1
      }
      #otherwise drop down to test C vs C+N
      else if(lrNext[2,"Pr(>Chisq)"]<(0.05/bfAdj)){
        nonrepND <- nonrepND + 1
      }
      #otherwise drop down to test null vs C
      else if(lrCurr[2,"Pr(>Chisq)"]<(0.05/bfAdj)){
        nonrepCD <- nonrepCD + 1
      }
    }
      
   
    
  },silent=TRUE) 
}
sprintf("%s Alleys",o)

print(sprintf("Repeating fields best model is P+C+N: %s",repFull/actualRunRep))
print(sprintf("Repeating fields best mode is Current Dir: %s",repCD/actualRunRep))
print(sprintf("Repeating fields best model is Previous Dir: %s",repPD/actualRunRep))
print(sprintf("Repeating fields best model is Next Dir: %s",repND/actualRunRep))

print(sprintf("Non-repeating fields best model is P+C+N: %s",nonrepFull/actualRunNonrep))
print(sprintf("Non-repeating fields best mode is Current Dir: %s",nonrepCD/actualRunNonrep))
print(sprintf("Non-repeating fields best model is Previous Dir: %s",nonrepPD/actualRunNonrep))
print(sprintf("Non-repeating fields best model is Next Dir: %s",nonrepND/actualRunNonrep))















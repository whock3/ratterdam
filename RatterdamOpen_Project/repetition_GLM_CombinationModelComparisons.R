# Repetition Project
# 22-01-03 WH 
# Script to create all GLM models using combinations of CD,ND,PD
# and test their relative improvements over one another using LRT.
# Comparisons made between nested versions of models. 
# eg. Base > CD > CD + PD > CD + PD + ND etc


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


alleydf <- subset(alleydf, Repeating=='True')

startTimeKnots = 3
shuffle <- FALSE
nruns <- 1 # 1 if not shuffling as the GLM LRT is (practically) deterministic. otherwise num shuffles you want

bf <- 0.05/15 # 15 LRT comparisons we're doing

## lists for pvalues of each field and model comparison for each run 
# (real data gets 1 run, shuffle gets n runs, e.g, 1000 runs)
# single variables
p_base_versus_CD <- c()
p_base_versus_PD <- c()
p_base_versus_ND <- c()

# variable pairs
p_CD_versus_CD_PD <- c() 
p_CD_versus_CD_ND <- c()

p_PD_versus_PD_CD <- c()
p_PD_versus_PD_ND <- c()

p_ND_versus_ND_CD <- c()
p_ND_versus_ND_PD <- c()


# variable triplets
p_CD_PD_versus_CD_PD_ND <- c() 
p_CD_ND_versus_CD_ND_PD <- c()

p_PD_CD_versus_PD_CD_ND <- c()
p_PD_ND_versus_PD_ND_CD <- c()

p_ND_CD_versus_ND_CD_PD <- c()
p_ND_PD_versus_ND_PD_CD <- c()

## lists of proportion of fields significant by each LRT for a given run 
nsig_base_versus_CD <- c()
nsig_base_versus_PD <- c()
nsig_base_versus_ND <- c()

# variable pairs
nsig_CD_versus_CD_PD <- c() 
nsig_CD_versus_CD_ND <- c()

nsig_PD_versus_PD_CD <- c()
nsig_PD_versus_PD_ND <- c()

nsig_ND_versus_ND_CD <- c()
nsig_ND_versus_ND_PD <- c()


# variable triplets
nsig_CD_PD_versus_CD_PD_ND <- c() 
nsig_CD_ND_versus_CD_ND_PD <- c()

nsig_PD_CD_versus_PD_CD_ND <- c()
nsig_PD_ND_versus_PD_ND_CD <- c()

nsig_ND_CD_versus_ND_CD_PD <- c()
nsig_ND_PD_versus_ND_PD_CD <- c()


for(n in 1:nruns){
  print(n)
  for(o in c('V','H')){
    oriendf <- subset(alleydf, Orientation==o)
    
    for(fid in unique(oriendf$FieldID)){
      
      field_ <- subset(oriendf, FieldID == fid)
      
      n <- sample(nrow(field_))
      shuff.field <- data.frame(field_)
      shuff.field$Rate <- shuff.field$Rate[n]
      
      if(shuffle==TRUE){
        field <- shuff.field
        
      }
      else if(shuffle==FALSE){
        field <- field_
      }
  
      ### Create models
      # wrap everything in try blocks because some fields won't have enough vars for
      ## (order of terms shouldn't matter so some models may be redundant, except if 
      ## you're running a post-hoc ANOVA with an unbalanced design. But since I don't 
      ## know all the inner mechanics of glm() I will proceed as if order matters)
      try({
      basem <- glm(Rate+1 ~ ns(StartTimes,startTimeKnots), 
                   family='Gamma', data=field)
      
      cd <- glm(Rate+1 ~ ns(StartTimes,startTimeKnots) 
                + CurrDir, 
                family='Gamma', data=field)
      
      nd <- glm(Rate+1 ~ ns(StartTimes,startTimeKnots) 
                + NextDir, 
                family='Gamma', data=field)
      
      pd <- glm(Rate+1 ~ ns(StartTimes,startTimeKnots) 
                + PrevDir, 
                family='Gamma', data=field)
      
      cd_pd <- glm(Rate+1 ~ ns(StartTimes,startTimeKnots) 
                   + CurrDir 
                   + PrevDir, 
                   family='Gamma', data=field)
      
      cd_nd <- glm(Rate+1 ~ ns(StartTimes,startTimeKnots) 
                   + CurrDir 
                   + NextDir, 
                   family='Gamma', data=field)
      
      pd_cd <- glm(Rate+1 ~ ns(StartTimes,startTimeKnots) 
                   + PrevDir 
                   + CurrDir, 
                   family='Gamma', data=field)
      
      pd_nd <- glm(Rate+1 ~ ns(StartTimes,startTimeKnots) 
                   + PrevDir 
                   + NextDir, 
                   family='Gamma', data=field)
      
      nd_cd <- glm(Rate+1 ~ ns(StartTimes,startTimeKnots) 
                   + NextDir 
                   + CurrDir, 
                   family='Gamma', data=field)
      
      nd_pd <- glm(Rate+1 ~ ns(StartTimes,startTimeKnots) 
                   + NextDir 
                   + PrevDir, 
                   family='Gamma', data=field)
      
      cd_pd_nd <- glm(Rate+1 ~ ns(StartTimes,startTimeKnots) 
                      + CurrDir 
                      + PrevDir
                      + NextDir, 
                      family='Gamma', data=field)
      
      cd_nd_pd <- glm(Rate+1 ~ ns(StartTimes,startTimeKnots) 
                   + CurrDir 
                   + NextDir
                   + PrevDir, 
                   family='Gamma', data=field)
      
      pd_cd_nd <- glm(Rate+1 ~ ns(StartTimes,startTimeKnots) 
                     + PrevDir 
                     + CurrDir
                     + NextDir, 
                     family='Gamma', data=field)
      
      pd_nd_cd <- glm(Rate+1 ~ ns(StartTimes,startTimeKnots) 
                     + PrevDir 
                     + NextDir
                     + CurrDir, 
                     family='Gamma', data=field)
      
      nd_cd_pd <- glm(Rate+1 ~ ns(StartTimes,startTimeKnots) 
                      + NextDir 
                      + CurrDir
                      + PrevDir, 
                      family='Gamma', data=field)
      
      nd_pd_cd <- glm(Rate+1 ~ ns(StartTimes,startTimeKnots) 
                      + NextDir 
                      + PrevDir
                      + CurrDir, 
                      family='Gamma', data=field)
      
      
      ## Run LRTs on different models. Start with each individual direction and 
      # add other terms in all possible combinations in a tree-like manner.
      # Eg base > CD > CD + PD > CD + PD + ND
      #              > CD + ND > CD + ND + PD
      
      # single variables
      lrt_base_versus_CD <- lrtest(basem, cd)
      lrt_base_versus_PD <- lrtest(basem, pd)
      lrt_base_versus_ND <- lrtest(basem, nd)
      
      # variable pairs
      lrt_CD_versus_CD_PD <- lrtest(cd, cd_pd)
      lrt_CD_versus_CD_ND <- lrtest(cd, cd_nd)
      
      lrt_PD_versus_PD_CD <- lrtest(pd, pd_cd)
      lrt_PD_versus_PD_ND <- lrtest(pd, pd_nd)
      
      lrt_ND_versus_ND_CD <- lrtest(nd,nd_cd)
      lrt_ND_versus_ND_PD <- lrtest(nd,nd_pd)
      
      
      # variable triplets
      lrt_CD_PD_versus_CD_PD_ND <- lrtest(cd_pd, cd_pd_nd)
      lrt_CD_ND_versus_CD_ND_PD <- lrtest(cd_nd, cd_nd_pd)
      
      lrt_PD_CD_versus_PD_CD_ND <- lrtest(pd_cd, pd_cd_nd)
      lrt_PD_ND_versus_PD_ND_CD <- lrtest(pd_nd, pd_nd_cd)
      
      lrt_ND_CD_versus_ND_CD_PD <- lrtest(nd_cd, nd_cd_pd)
      lrt_ND_PD_versus_ND_PD_CD <- lrtest(nd_pd, nd_pd_cd)
      
      
      ## Save p-values (can apply multiple comparison correction post-hoc) 
      
      # single variables
      p_base_versus_CD <- c(p_base_versus_CD, lrt_base_versus_CD[2,"Pr(>Chisq)"])
      p_base_versus_PD <- c(p_base_versus_PD, lrt_base_versus_PD[2,"Pr(>Chisq)"])
      p_base_versus_ND <- c(p_base_versus_ND, lrt_base_versus_ND[2,"Pr(>Chisq)"])
      
      # variable pairs
      p_CD_versus_CD_PD <- c(p_CD_versus_CD_PD, lrt_CD_versus_CD_PD[2,"Pr(>Chisq)"])
      p_CD_versus_CD_ND <- c(p_CD_versus_CD_ND, lrt_CD_versus_CD_ND[2,"Pr(>Chisq)"])
      
      p_PD_versus_PD_CD <- c(p_PD_versus_PD_CD, lrt_PD_versus_PD_CD[2,"Pr(>Chisq)"])
      p_PD_versus_PD_ND <- c(p_PD_versus_PD_ND, lrt_PD_versus_PD_ND[2,"Pr(>Chisq)"])
      
      p_ND_versus_ND_CD <- c(p_ND_versus_ND_CD, lrt_ND_versus_ND_CD[2,"Pr(>Chisq)"])
      p_ND_versus_ND_PD <- c(p_ND_versus_ND_PD, lrt_ND_versus_ND_PD[2,"Pr(>Chisq)"])
      
      
      # variable triplets
      p_CD_PD_versus_CD_PD_ND <- c(p_CD_PD_versus_CD_PD_ND, lrt_CD_PD_versus_CD_PD_ND[2,"Pr(>Chisq)"])
      p_CD_ND_versus_CD_ND_PD <- c(p_CD_ND_versus_CD_ND_PD, lrt_CD_ND_versus_CD_ND_PD[2,"Pr(>Chisq)"])
      
      p_PD_CD_versus_PD_CD_ND <- c(p_PD_CD_versus_PD_CD_ND, lrt_PD_CD_versus_PD_CD_ND[2,"Pr(>Chisq)"])
      p_PD_ND_versus_PD_ND_CD <- c(p_PD_ND_versus_PD_ND_CD, lrt_PD_ND_versus_PD_ND_CD[2,"Pr(>Chisq)"])
      
      p_ND_CD_versus_ND_CD_PD <- c(p_ND_CD_versus_ND_CD_PD, lrt_ND_CD_versus_ND_CD_PD[2,"Pr(>Chisq)"])
      p_ND_PD_versus_ND_PD_CD <- c(p_ND_PD_versus_ND_PD_CD, lrt_ND_PD_versus_ND_PD_CD[2,"Pr(>Chisq)"])
      
      },silent=FALSE) # try per field
      
    }
    
  }
      nsig_base_versus_CD <- c(nsig_base_versus_CD, sum(na.omit(p_base_versus_CD) < bf)/length(na.omit(p_base_versus_CD)))
      nsig_base_versus_PD <- c(nsig_base_versus_PD, sum(na.omit(p_base_versus_PD) < bf)/length(na.omit(p_base_versus_PD)))
      nsig_base_versus_ND <- c(nsig_base_versus_ND, sum(na.omit(p_base_versus_ND) < bf)/length(na.omit(p_base_versus_ND)))
      
      # variable pairs
      nsig_CD_versus_CD_PD <- c(nsig_CD_versus_CD_PD, sum(na.omit(p_CD_versus_CD_PD) < bf)/length(na.omit(p_CD_versus_CD_PD))) 
      nsig_CD_versus_CD_ND <- c(nsig_CD_versus_CD_ND, sum(na.omit(p_CD_versus_CD_ND) < bf)/length(na.omit(p_CD_versus_CD_ND)))
      
      nsig_PD_versus_PD_CD <- c(nsig_PD_versus_PD_CD, sum(na.omit(p_PD_versus_PD_CD) < bf)/length(na.omit(p_PD_versus_PD_CD)))
      nsig_PD_versus_PD_ND <- c(nsig_PD_versus_PD_ND, sum(na.omit(p_PD_versus_PD_ND) < bf)/length(na.omit(p_PD_versus_PD_ND)))
      
      nsig_ND_versus_ND_CD <- c(nsig_ND_versus_ND_CD, sum(na.omit(p_ND_versus_ND_CD) < bf)/length(na.omit(p_ND_versus_ND_CD)))
      nsig_ND_versus_ND_PD <- c(nsig_ND_versus_ND_PD, sum(na.omit(p_ND_versus_ND_PD) < bf)/length(na.omit(p_ND_versus_ND_PD)))
      
      
      # variable triplets
      nsig_CD_PD_versus_CD_PD_ND <- c(nsig_CD_PD_versus_CD_PD_ND, sum(na.omit(p_CD_PD_versus_CD_PD_ND) < bf)/length(na.omit(p_CD_PD_versus_CD_PD_ND))) 
      nsig_CD_ND_versus_CD_ND_PD <- c(nsig_CD_ND_versus_CD_ND_PD, sum(na.omit(p_CD_ND_versus_CD_ND_PD) < bf)/length(na.omit(p_CD_ND_versus_CD_ND_PD)))
      
      nsig_PD_CD_versus_PD_CD_ND <- c(nsig_PD_CD_versus_PD_CD_ND, sum(na.omit(p_PD_CD_versus_PD_CD_ND) < bf)/length(na.omit(p_PD_CD_versus_PD_CD_ND)))
      nsig_PD_ND_versus_PD_ND_CD <- c(nsig_PD_ND_versus_PD_ND_CD, sum(na.omit(p_PD_ND_versus_PD_ND_CD) < bf)/length(na.omit(p_PD_ND_versus_PD_ND_CD))) 
      
      nsig_ND_CD_versus_ND_CD_PD <- c(nsig_ND_CD_versus_ND_CD_PD, sum(na.omit(p_ND_CD_versus_ND_CD_PD) < bf)/length(na.omit(p_ND_CD_versus_ND_CD_PD)))
      nsig_ND_PD_versus_ND_PD_CD <- c(nsig_ND_PD_versus_ND_PD_CD, sum(na.omit(p_ND_PD_versus_ND_PD_CD) < bf)/length(na.omit(p_ND_PD_versus_ND_PD_CD)))
  
  
  
  
}

# Print results 
print("Base versus CD")
print(sum(na.omit(p_base_versus_CD) < bf)/length(na.omit(p_base_versus_CD)))

print("Base versus PD")
print(sum(na.omit(p_base_versus_PD) < bf)/length(na.omit(p_base_versus_PD)))

print("Base versus ND")
print(sum(na.omit(p_base_versus_ND) < bf)/length(na.omit(p_base_versus_ND)))

print("CD versus CD+PD")
print(sum(na.omit(p_CD_versus_CD_PD) < bf)/length(na.omit(p_CD_versus_CD_PD)))

print("CD versus CD+ND")
print(sum(na.omit(p_CD_versus_CD_ND) < bf)/length(na.omit(p_CD_versus_CD_ND)))

print("PD versus PD+CD")
print(sum(na.omit(p_PD_versus_PD_CD) < bf)/length(na.omit(p_PD_versus_PD_CD)))

print("PD versus PD+ND")
print(sum(na.omit(p_PD_versus_PD_ND) < bf)/length(na.omit(p_PD_versus_PD_ND)))

print("ND versus ND+CD")
print(sum(na.omit(p_ND_versus_ND_CD) < bf)/length(na.omit(p_ND_versus_ND_CD)))

print("ND versus ND+PD")
print(sum(na.omit(p_ND_versus_ND_PD) < bf)/length(na.omit(p_ND_versus_ND_PD)))

print("CD+PD versus CD+PD+ND")
print(sum(na.omit(p_CD_PD_versus_CD_PD_ND) < bf)/length(na.omit(p_CD_PD_versus_CD_PD_ND)))

print("CD+ND versus CD+ND+PD")
print(sum(na.omit(p_CD_ND_versus_CD_ND_PD) < bf)/length(na.omit(p_CD_ND_versus_CD_ND_PD)))

print("PD+CD versus PD+CD+ND")
print(sum(na.omit(p_PD_CD_versus_PD_CD_ND) < bf)/length(na.omit(p_PD_CD_versus_PD_CD_ND)))

print("PD+ND versus PD+ND+CD")
print(sum(na.omit(p_PD_ND_versus_PD_ND_CD) < bf)/length(na.omit(p_PD_ND_versus_PD_ND_CD)))

print("ND+CD versus ND+CD+PD")
print(sum(na.omit(p_ND_CD_versus_ND_CD_PD) < bf)/length(na.omit(p_ND_CD_versus_ND_CD_PD)))

print("ND+PD versus ND+PD+CD")
print(sum(na.omit(p_ND_PD_versus_ND_PD_CD) < bf)/length(na.omit(p_ND_PD_versus_ND_PD_CD)))



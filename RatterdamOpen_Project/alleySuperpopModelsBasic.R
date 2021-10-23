# 10-3-21 LMER Script for Alley Superpopulation
# Doing single model for whole dataset with cell/fields as RE

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

# mod <- lmer(Rate ~ CurrDir + 
#               (CurrDir|CellID/FieldID)
#             ,data=valleydf)

# 
# mod <- lmer(Rate ~ CurrDir + 
#                     PrevDir + 
#                     NextDir + 
#                     (CurrDir|FieldID) +
#                     (PrevDir|FieldID) +
#                     (NextDir|FieldID)
#                     , data=valleydf)

# ALLOCENTRIC
oalleydf <- subset(alleydf, Orientation=='H')

repsigCD <- 0
repsigPD <- 0
repsigND <- 0

nonrepsigCD <- 0
nonrepsigPD <- 0
nonrepsigND <- 0

for(fid in unique(oalleydf$FieldID)){
  field <- subset(oalleydf, FieldID==fid)
  mrun <- try({
  mod <- lm(Rate ~ CurrDir + NextDir + PrevDir, data=field)
  a <- anova(mod)
  },silent=TRUE)
  
  if(class(mrun)!='try-error'){
  x <- try({
  #check current direction 
  if(a['CurrDir','Pr(>F)'] < 0.05){
    if(unique(field$Repeating)=='True'){
      repsigCD <- repsigCD + 1
    }
    else if(unique(field$Repeating)=='False'){
      nonrepsigCD <- nonrepsigCD + 1
    }
  }
  },silent=TRUE)
  
  x <- try({
  #check previous direction 
  if(a['PrevDir','Pr(>F)'] < 0.05){
    if(unique(field$Repeating)=='True'){
      repsigPD <- repsigPD + 1
    }
    else if(unique(field$Repeating)=='False'){
      nonrepsigPD <- nonrepsigPD + 1
    }
  }
  },silent=TRUE)
  
  
  x<-try({
  #check next direction 
  if(a['NextDir','Pr(>F)'] < 0.05){
    if(unique(field$Repeating)=='True'){
      repsigND <- repsigND + 1
    }
    else if(unique(field$Repeating)=='False'){
      nonrepsigND <- nonrepsigND + 1
    }
  }
  },silent=TRUE)
  
  }
}

print(repsigCD)
print(repsigPD)
print(repsigND)

print(nonrepsigCD)
print(nonrepsigPD)
print(nonrepsigND)



# EGOCENTRIC


repsigCD <- 0
repsigPD <- 0
repsigND <- 0

nonrepsigCD <- 0
nonrepsigPD <- 0
nonrepsigND <- 0

for(fid in unique(oalleydf$FieldID)){
  field <- subset(oalleydf, FieldID==fid)
  mrun <- try({
  mod <- lm(Rate ~ CurrDir + RetroEgo + ProspEgo, data=field)
  a <- anova(mod)
  },silent=TRUE)
  
  if(class(mrun)!='try-error'){
  x <- try({
    #check current direction 
    if(a['CurrDir','Pr(>F)'] < 0.05){
      if(unique(field$Repeating)=='True'){
        repsigCD <- repsigCD + 1
      }
      else if(unique(field$Repeating)=='False'){
        nonrepsigCD <- nonrepsigCD + 1
      }
    }
  },silent=TRUE)
  
  x <- try({
    #check previous direction 
    if(a['RetroEgo','Pr(>F)'] < 0.05){
      if(unique(field$Repeating)=='True'){
        repsigPD <- repsigPD + 1
      }
      else if(unique(field$Repeating)=='False'){
        nonrepsigPD <- nonrepsigPD + 1
      }
    }
  },silent=TRUE)
  
  
  x<-try({
    #check next direction 
    if(a['ProspEgo','Pr(>F)'] < 0.05){
      if(unique(field$Repeating)=='True'){
        repsigND <- repsigND + 1
      }
      else if(unique(field$Repeating)=='False'){
        nonrepsigND <- nonrepsigND + 1
      }
    }
  },silent=TRUE)
  
  }
}

print(repsigCD)
print(repsigPD)
print(repsigND)

print(nonrepsigCD)
print(nonrepsigPD)
print(nonrepsigND)


#### Repeating vs Nonrepeating tested via MSE across fields in each pop

replr <- c()
repmse <- c()
repaic <- c()

nonreplr <- c()
nonrepmse <- c()
nonrepaic <- c()



oalleydf <- subset(alleydf, Orientation=='H')

for(fid in unique(oalleydf$FieldID)){
  field <- subset(oalleydf, FieldID==fid)
  
  x <- try({
    gmod <- glm(Rate+1 ~ PrevDir + CurrDir + NextDir, data=field, family='Gamma')
    gnull <- glm(Rate+1 ~ CurrDir, family='Gamma', data=field)
    },silent=TRUE)
  
  if(class(x)!='try-error'){
    
    if(unique(field$Repeating)=='True'){
      replr <- c(replr, lrtest(gmod,gnull)[2,'Pr(>Chisq)'])
      repmse <- c(repmse, sqrt(mean(gmod$residuals^2)))
      repaic <- c(repaic, gmod$aic-gnull$aic)
    }
    else if(unique(field$Repeating)=='False'){
      nonreplr <- c(nonreplr, lrtest(gmod,gnull)[2,'Pr(>Chisq)'])
      nonrepmse <- c(nonrepmse, sqrt(mean(gmod$residuals^2)))
      nonrepaic <- c(nonrepaic, gmod$aic - gnull$aic)
      
    }
  }
}


# m2r <- length(replr[replr<0.05])/length(replr)
# m2n <- length(nonreplr[nonreplr<0.05])/length(nonreplr)

m2r <- repaic
m2n <- nonrepaic



length(repaic[repaic<0])/length(repaic)
length(nonrepaic[nonrepaic<0])/length(nonrepaic)

boxplot(repaic,nonrepaic,ylim=c(-50,10))








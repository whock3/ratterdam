# Repetition Project William Hockeimer 9/20/21
# Field LMER Models without temporal regressor 
# Looking at prospective and retrospective coding in allocentric framework 

library(lme4)
library(ggplot2)
library(lmtest)
library(stringr)
library(emmeans)
library(tidyr)
library(ggpubr)
library(splines)
# define data paths and load in data
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


# Allocentric frame modeling 

#subset dfs
repalleydf <- subset(alleydf, Repeating == 'True')
repinterdf <- subset(interdf, Repeating == 'True')

nonrepalleydf <- subset(alleydf, Repeating == 'False')
nonrepinterdf <- subset(interdf, Repeating = 'False')

#initalize data vectors. Making comparisons between prosp and retro encoding
# at alleys and intersections between repeating and nonrepeating cells 
repAlleyProspAllo <- c()
repAlleyProspAllo_count <- 0
repAlleyRetroAllo <- c()
repAlleyRetroAllo_count <- 0

repInterProspAllo <- c()
repInterProspAllo_count <- 0
repInterRetroAllo <- c()
repInterRetroAllo_count <- 0

nonrepAlleyProspAllo <- c()
nonrepAlleyProspAllo_count <- 0
nonrepAlleyRetroAllo <- c()
nonrepAlleyRetroAllo_count <- 0

nonrepInterProspAllo <- c()
nonrepInterProspAllo_count <- 0
nonrepInterRetroAllo <- c()
nonrepInterRetroAllo_count <- 0

# Repeating fields in alley, allo frame prosp/retro coding
for(fid in unique(repalleydf$FieldID)){
  fdf <- subset(repalleydf, FieldID == fid)
  
  try({
    p <- aov(Rate ~ PrevDir, data=fdf)
    repAlleyRetroAllo <- c(repAlleyRetroAllo, eta_sq(p)['etasq']$etasq)
    if(summary(p)[[1]][["Pr(>F)"]][[1]] < 0.05/4){
      repAlleyRetroAllo_count <- repAlleyRetroAllo_count + 1
    }
  },silent=TRUE)
  
  try({
    n <- aov(Rate ~ NextDir, data=fdf)
    repAlleyProspAllo <- c(repAlleyProspAllo, eta_sq(n)['etasq']$etasq)
    if(summary(n)[[1]][["Pr(>F)"]][[1]] < 0.05/4){
      repAlleyProspAllo_count <- repAlleyProspAllo_count + 1
    }
  },silent=TRUE)
  
}

# Repeating fields in intersection, allo frame prosp/retro coding
for(fid in unique(repinterdf$FieldID)){
  fdf <- subset(repinterdf, FieldID == fid)
  
  try({
    p <- aov(Rate ~ PrevDir, data=fdf)
    repInterRetroAllo <- c(repInterRetroAllo, eta_sq(p)['etasq']$etasq)
    if(summary(p)[[1]][["Pr(>F)"]][[1]] < 0.05/4){
      repInterRetroAllo_count <- repInterRetroAllo_count + 1
    }
  },silent=TRUE)
  
  
  try({
    n <- aov(Rate ~ NextDir, data=fdf)
    repInterProspAllo <- c(repInterProspAllo, eta_sq(n)['etasq']$etasq)
    if(summary(n)[[1]][["Pr(>F)"]][[1]] < 0.05/4){
      repInterProspAllo_count <- repInterProspAllo_count + 1
    }
  },silent=TRUE)
  
}

# Non-repeating fields in alley, allo frame prosp/retro coding
for(fid in unique(nonrepalleydf$FieldID)){
  fdf <- subset(nonrepalleydf, FieldID == fid)
  
  try({
    p <- aov(Rate ~ PrevDir, data=fdf)
    nonrepAlleyRetroAllo <- c(nonrepAlleyRetroAllo, eta_sq(p)['etasq']$etasq)
    if(summary(p)[[1]][["Pr(>F)"]][[1]] < 0.05/4){
      nonrepAlleyRetroAllo_count <- nonrepAlleyRetroAllo_count + 1
    }
  },silent=TRUE)
  
  try({
    n <- aov(Rate ~ NextDir, data=fdf)
    nonrepAlleyProspAllo <- c(nonrepAlleyProspAllo, eta_sq(n)['etasq']$etasq)
    if(summary(n)[[1]][["Pr(>F)"]][[1]] < 0.05/4){
      nonrepAlleyProspAllo_count <- nonrepAlleyProspAllo_count + 1
    }
  },silent=TRUE)
  
}

# Non-repeating fields in intersection, allo frame prosp/retro coding
for(fid in unique(nonrepinterdf$FieldID)){
  fdf <- subset(nonrepinterdf, FieldID == fid)
  
  try({
    p <- aov(Rate ~ PrevDir, data=fdf)
    nonrepInterRetroAllo <- c(nonrepInterRetroAllo, eta_sq(p)['etasq']$etasq)
    if(summary(p)[[1]][["Pr(>F)"]][[1]] < 0.05/4){
      nonrepInterRetroAllo_count <- nonrepInterRetroAllo_count + 1
    }
  },silent=TRUE)
  
  try({
    n <- aov(Rate ~ NextDir, data=fdf)
    nonrepInterProspAllo <- c(nonrepInterProspAllo, eta_sq(n)['etasq']$etasq)
    if(summary(n)[[1]][["Pr(>F)"]][[1]] < 0.05/4){
      nonrepInterProspAllo_count <- nonrepInterProspAllo_count + 1
    }
  },silent=TRUE)
  
}


# Run binomial test on counts, see num fields w sig effects (broken  down by groups as above)
count <- repAlleyRetroAllo_count
numfields <- length(unique(repalleydf$FieldID))
bt <- binom.test(count,numfields, p=0.05, alternative='greater')
sprintf("Repeating fields in alleys, retrospective coding: %s, p=%s",count/numfields,  bt$p.value)

count <- repAlleyProspAllo_count
numfields <- length(unique(repalleydf$FieldID))
bt <- binom.test(count,numfields, p=0.05, alternative='greater')
sprintf("Repeating fields in alleys, prospective coding: %s, p=%s",count/numfields,  bt$p.value)

count <- repInterRetroAllo_count
numfields <- length(unique(repinterdf$FieldID))
bt <- binom.test(count,numfields, p=0.05, alternative='greater')
sprintf("Repeating fields in intersections, retrospective coding: %s, p=%s",count/numfields,  bt$p.value)

count <- repInterProspAllo_count
numfields <- length(unique(repinterdf$FieldID))
bt <- binom.test(count,numfields, p=0.05, alternative='greater')
sprintf("Repeating fields in intersections, prospective coding: %s, p=%s",count/numfields,  bt$p.value)

count <- nonrepAlleyRetroAllo_count
numfields <- length(unique(nonrepalleydf$FieldID))
bt <- binom.test(count,numfields, p=0.05, alternative='greater')
sprintf("Non-repeating fields in alleys, retrospective coding: %s, p=%s",count/numfields,  bt$p.value)

count <- nonrepAlleyProspAllo_count
numfields <- length(unique(nonrepalleydf$FieldID))
bt <- binom.test(count,numfields, p=0.05, alternative='greater')
sprintf("Non-repeating fields in alleys, prospective coding: %s, p=%s",count/numfields,  bt$p.value)

count <- nonrepInterRetroAllo_count
numfields <- length(unique(nonrepinterdf$FieldID))
bt <- binom.test(count,numfields, p=0.05, alternative='greater')
sprintf("Non-repeating fields in intersections, retrospective coding: %s, p=%s",count/numfields,  bt$p.value)

count <- nonrepInterProspAllo_count
numfields <- length(unique(nonrepinterdf$FieldID))
bt <- binom.test(count,numfields, p=0.05, alternative='greater')
sprintf("Non-repeating fields in intersections, prospective coding: %s, p=%s",count/numfields,  bt$p.value)

      
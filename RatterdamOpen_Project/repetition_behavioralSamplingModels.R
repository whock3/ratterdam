# GLMs to test the effect of behavioral sampling on firing rates.


library(lme4)
library(ggplot2)
library(lmtest)
library(stringr)
library(emmeans)
library(tidyr)
library(ggpubr)
library(splines)


alleypath <- "E:\\Ratterdam\\R_data_repetition\\2022-04-05_AlleySuperpopDirVisitFiltered.csv"

alleydf <- read.csv(alleypath,header=TRUE)

alleydf$PrevDir <- as.factor(alleydf$PrevDir)
alleydf$CurrDir <- as.factor(alleydf$CurrDir)
alleydf$NextDir <- as.factor(alleydf$NextDir)
alleydf$RetroEgo <- as.factor(alleydf$RetroEgo)
alleydf$ProspEgo <- as.factor(alleydf$ProspEgo)
alleydf$Repeating <- as.factor(alleydf$Repeating)
alleydf$Traversal <- as.factor(alleydf$Traversal)
alleydf$Reward <- as.factor(alleydf$Reward)
alleydf$CellID <- as.factor(alleydf$CellID)
alleydf$FieldNum <- as.factor(alleydf$FieldNum)
alleydf$FieldID <- as.factor(alleydf$FieldID)
alleydf$Alleys <- as.factor(alleydf$Alleys)
alleydf$NumFields <- as.numeric(alleydf$NumFields)

rdf <- subset(alleydf, Rat=='R781' & Day == 'D3')

rdf$StartTimes <- (rdf$StartTimes - min(rdf$StartTimes))/1e6

wins <- c()


window <- 15*60
offset <- 1*60
wins <- c()
begin <- 0
stop <- FALSE

while(stop==FALSE){
  a <- begin
  b <- begin + window
  if( b < ceiling(max(rdf$StartTimes))){
    wins <- append(wins, c(a,b))
    begin <- begin + offset
  }
  else {
    stop = TRUE
  }

}


# Repetition Project - Detecting and Analyzing Multicollinearity
# WH 2022-01-13

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


vifs <- data.frame(matrix(nrow=0,ncol=3))
colnames(vifs) <- cols

total <- 0

m <- glm(Rate+1 ~ CurrDir+PrevDir+NextDir,
family='Gamma',
data=field)
al <- alias(m)
"Complete" %in% names(al)
alidx = which(al$Complete!=0, arr.ind=T)
row.names(al$Complete)[alidx[1,1]]
colnames(al$Complete)[alidx[1,2]]

#do this for each entry in alidx
factor_level <- substr(dl,nchar(dl),nchar(dl))
factor_type <- substr(dl,1,nchar(dl)-1)

field2 <- droplevels(subset(field, CurrDir!="W" & NextDir!="S"))





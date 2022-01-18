# Repetition Project - Using bootstrapped lasso (Bolasso)
# Purpose is to use bolasso to select best regressors among
# those related to direction, while accounting for multicollinearity,
# and use bootstrapm estimates of variable importance to test which, if any,
# variables each field cares about
# WH 2022-01-18

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



for(o in c("V", "H")){
  
  oriendf <- subset(alleydf, Orientation == o)
  
  for(fid in unique(oriendf$FieldID)){
    
    field <- subset(oriendf, FieldID == fid)
    
    x <-model.matrix(Rate ~ -1+CurrDir+NextDir+PrevDir,data=droplevels(field))[,-1]
    gm <- gglasso(x=x,y=field$Rate,lambda=1)
   
    
  }

  
}



create_factor <- function(nb_lvl, n= 100 ){
  factor(sample(letters[1:nb_lvl],n, replace = TRUE))}

df <- data.frame(var1 = create_factor(5), 
                 var2 = create_factor(5), 
                 var3 = create_factor(5), 
                 var4 = create_factor(5),
                 var5 = rnorm(100),
                 y = rnorm(100))

y <- df$y
x <- model.matrix( ~ ., dplyr::select(df, -y))[, -1]
groups <- c(rep(1:4, each = 4), 5)
fit <- gglasso(x = x, y = y, group = groups, lambda = 1)
fit$beta
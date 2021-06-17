

# LMER Routine - Beltway Rate Remapping Project
# WH March 20 2021

# Functions to take a dataframe, run a (G)LMER model on it, and return statistics
# such as confidence intervals, coefficients, 

# imports
library(lme4)
library(splines)
library(ggplot2)
library(lmtest)
library(stringr)


setwd("E:\\UserData\\Documents\\GitHub\\ratterdam\\")
source("E:\\UserData\\Documents\\GitHub\\ratterdam\\Beltway_Project\\cicheck.R")
source("E:\\UserData\\Documents\\GitHub\\ratterdam\\Beltway_Project\\glmer_fx.R")

lmer_routine <- function(celldf, nalleys){
  
  nsplineknots <- 6 # 6 if the alley includes intersections and is ~24bins. 12bins if its the expanded alley that only includes a portion of intersection

  bonf <- 0.05/(3*nsplineknots*nalleys) #*8 for reassignment analysis
  cipct <- 1-bonf
  z <- qnorm(cipct)
  
  # run model
  modi <- lmer_feTrial(celldf)

  # Reorder stimulus factors to make other comparisons
  celldf$texture <- factor(celldf$texture, levels = c("B", "A", "C"))
  modreord <- lmer_feTrial(celldf)
  celldf$texture <- factor(celldf$texture, levels = c("A","B","C"))
  
  
  # fit for model with int
  celldf$fit <- predict(modi, newdata=celldf, re.form=NA)
  
  # for interaction model, Wald CI. txt A as default
  # ci <- confint(modi, method='Wald', level=cipct)
  # ci <- ci[-c(1,2),]
  # fe <- fixef(modi)
  # cidf <- data.frame("low"=ci[,1], "up"=ci[,2],"fe"=fe)
  # cidf <- cidf[-c(1),]
  
  # do CI on model with B as default and concat 
  # cireord <- confint(modreord, method='Wald', level=cipct)
  # cireord <- cireord[-c(1,2),]
  # fereord <- fixef(modreord)
  # cireorddf <- data.frame("low"=cireord[,1], "up"=cireord[,2],"fe"=fereord)
  # ciall <- rbind(cidf, cireorddf[c("textureC", "ns(spatialBin, 6)1:textureC","ns(spatialBin, 6)2:textureC","ns(spatialBin, 6)3:textureC","ns(spatialBin, 6)4:textureC","ns(spatialBin, 6)5:textureC","ns(spatialBin, 6)6:textureC"),])
  # ndf <- data.frame(ciall,cellname,alley)
  # wdf <- rbind(wdf, ndf)
  
  Designmat <- model.matrix(rate ~ ns(spatialBin, 6)*texture+reward, celldf)
  predvar <- diag(Designmat %*% vcov(modi) %*% t(Designmat))
  celldf$fitCI <- sqrt(predvar)*z
  
  
  return(celldf)
}

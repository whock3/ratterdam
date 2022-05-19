# Repetition Manuscript 
# Same/different corridor regression analysis
# Comparing field sampling bias with directionality product
# WH 2022-05-17

library(lme4)
library(ggplot2)
library(lmtest)
library(stringr)
library(emmeans)
library(tidyr)
library(ggpubr)
library(splines)

cpath <- "E:\\Ratterdam\\repetition_manuscript\\Figure5_RepeatingFieldNonconservation\\corridordf.csv"
cdf <- read.csv(cpath,header=TRUE)

cdf_share <- subset(cdf, shareness==1)
mod <- glm(fieldB ~ fieldA + prods,data = cdf_share)
print(summary(mod))

cdf_nonshare <- subset(cdf, shareness==0)
mod <- glm(fieldB ~ fieldA + prods, data = cdf_nonshare)
print(summary(mod))
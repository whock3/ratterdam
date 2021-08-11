# Repetition Project
# Will Hockeimer July 6 2021
# Examining role of directionality on fields

library(lme4)
library(ggplot2)
library(lmtest)
library(stringr)
library(emmeans)
library(tidyr)
library(ggpubr)
library(splines)
# define data paths and load in data
path <- "E:\\Ratterdam\\R_data_repetition\\20210810-164908_R765DFD8_timedirLMER_1.5vfilt.csv"
code <- "R859D2"
df <- read.csv(path,header=TRUE)
df <- df[!(df$dirC==0),]
df <- df[!(df$dirP1==0),]
df <- df[!(df$dirM1==0),]

ts <- str_replace(Sys.time()," ","_")
ts <- str_replace_all(ts, ":", "_")

df$field <- as.factor(df$field)
df$dirC <- as.factor(df$dirC)
df$dirM1 <- as.factor(df$dirM1)
df$dirP1 <- as.factor(df$dirP1)
df$epoch <- as.factor(df$epoch)


df$dirC <- plyr::revalue(df$dirC, c("1"="N","2"="E","3"="S","4"="W"))
df$dirM1 <- plyr::revalue(df$dirM1, c("1"="N","2"="E","3"="S","4"="W"))
df$dirP1 <- plyr::revalue(df$dirP1, c("1"="N","2"="E","3"="S","4"="W"))


bonf <- 0.05
cipct <- 1-bonf
z <- qnorm(cipct)

unitname <- 'TT15\\cl-maze1.4'
celldf <- subset(df, unit == unitname)

f = '0'
fdf <- subset(celldf, field == f)
mod <- lm(rate ~ ns(timestamp,4)*dirC, data=fdf)

fdf$fit <- predict(mod, newdata=fdf, re.form=NA)

# Designmat <- model.matrix(rate ~ ns(visitIdx, 3)*dirC, fdf)
# predvar <- diag(Designmat %*% vcov(mod) %*% t(Designmat))
# fdf$fitCI <- sqrt(predvar)*z

p <- ggplot(data=fdf, aes(x=timestamp, y=rate, color=dirC))+
  geom_smooth(method='loess')+
  geom_point()+
  ggtitle(sprintf("Cell %s Field %s", unitname, f))
print(p)
 







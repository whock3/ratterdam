# Repetition Project
# Will Hockeimer July 6 2021
# Examining role of directionality on fields
# 
# Sep 2021 - Directional LMER using segmented field traversal\
# csv created by taking all turns through a field on a given traversal and
# grouping by heading (with multiple headings in the same direction concatenated
# if theyre consecutive or separated if theyre separated by a different heading)
# 

library(lme4)
library(ggplot2)
library(lmtest)
library(stringr)
library(emmeans)
library(tidyr)
library(ggpubr)
library(splines)
# define data paths and load in data
path <- "E:\\Ratterdam\\R_data_repetition\\20210903-102901_R781D3_segmentedDirections_1.5vfilt.csv"
code <- "R781D3"
df <- read.csv(path,header=TRUE)

ts <- str_replace(Sys.time()," ","_")
ts <- str_replace_all(ts, ":", "_")

df <- df[df$Direction != 'X',]

df$Field <- as.factor(df$Field)
df$Direction <- as.factor(df$Direction)

bonf <- 0.05
cipct <- 1-bonf
z <- qnorm(cipct)

unitname <- 'TT2\\cl-maze1.2'
celldf <- subset(df, Unit == unitname)

f = '1'
fdf <- subset(celldf, Field == f)
mod <- lm(Rate ~ ns(StartTime,3)*Direction, data=fdf)

fdf$fit <- predict(mod, newdata=fdf, re.form=NA)

# Designmat <- model.matrix(rate ~ ns(visitIdx, 3)*dirC, fdf)
# predvar <- diag(Designmat %*% vcov(mod) %*% t(Designmat))
# fdf$fitCI <- sqrt(predvar)*z

p <- ggplot(data=fdf, aes(x=StartTime, y=Rate, color=Direction))+
  geom_smooth(method='loess')+
  geom_point()+
  ggtitle(sprintf("Cell %s Field %s", unitname, f))
print(p)
 
print(anova(mod))







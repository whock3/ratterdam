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

# define data paths and load in data
path <- "E:\\Ratterdam\\R_data_repetition\\20210713-182039_R808D6_dirLMER_1.5vfilt.csv"
savepath <- "E:\\Ratterdam\\R_data_repetition\\210624-onward_turnLM\\"
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

results <- data.frame()
for(unitname in unique(df$unit)){
  print(unitname)
  celldf <- subset(df, unit==unitname)
  
  u <- try({
    res <- aov(rate ~ dirC*field, data=celldf)
    resdf <- summary(res)[[1]]
    p <- resdf["dirC:field","Pr(>F)"]
    results  <- rbind(results, c(unitname, p))
    print("completed")
  },silent=FALSE)
  
}



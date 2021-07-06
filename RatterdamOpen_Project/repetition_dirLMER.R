# Repetition Project
# Will Hockeimer July 6 2021
# Examining role of directionality on fields

library(lme4)
library(splines)
library(ggplot2)
library(lmtest)
library(stringr)

# define data paths and load in data
path <- "E:\\Ratterdam\\R_data_repetition\\20210706-142623_R781D3_dirLMER_1.5vfilt.csv"
savepath <- "E:\\Ratterdam\\R_data_repetition\\210624-onward_turnLM\\"
code <- "R781D3"
df <- read.csv(path,header=TRUE)
df <- df[!(df$direction==0),]
df$field <- as.factor(df$field)
df$direction <- as.factor(df$direction)
df$epoch <- as.factor(df$epoch)
ts <- str_replace(Sys.time()," ","_")
ts <- str_replace_all(ts, ":", "_")


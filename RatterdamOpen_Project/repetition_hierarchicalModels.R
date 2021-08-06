# Repetition Project
# WH 2021-07-29
# Hierarchical Models
# Goal is to create series of (g)lmers with increasing model complexity
# and perform model comparison to see what level of complexity (w.r.t route
# coding) best capture variance in data

# Load Libraries and data ---- 

library(lme4)
library(lmtest)
library(stringr)
library(emmeans)
library(tidyverse)
library(ggpubr)

path <- "E:\\Ratterdam\\R_data_repetition\\20210729-161818_R859D2_dirLMER_1.5vfilt.csv"
savepath <- "E:\\Ratterdam\\R_data_repetition\\210624-onward_turnLM\\"
code <- "R859D1"
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


runHierarchicalModel_ <- function(fdf){
  
  mod_dir <- lmer(rate ~ dirC + (1|epoch), data=fdf)
  
  cf <- confint(mod_dir, method='Wald', level=0.95)
  
  # 4th row should be direction Avs B comparison. If they nonoverlap zero, product is positive
  if(sign(cf[4,"2.5 %"]*cf[4,"97.5 %"])==1){
    
    mod_prev <- lmer(rate ~ dirC*dirM1 + (1|epoch), data=fdf)
    mod_next <- lmer(rate ~ dirC*dirP1 + (1|epoch), data=fdf)
    mod_traj <- lmer(rate ~ dirC*dirP1*dirM1 + (1|epoch), data=fdf)
    
    lr_dir_prev <- lrtest(mod_dir, mod_prev)
    p_dir_prev <- lr_dir_prev$`Pr(>Chisq)`[2] # 1st value is NA
    
    lr_dir_next <- lrtest(mod_dir, mod_next)
    p_dir_next <- lr_dir_next$`Pr(>Chisq)`[2]
    
    lr_prev_traj <- lrtest(mod_prev, mod_traj)
    p_prev_traj <- lr_prev_traj$`Pr(>Chisq)`[2]
    
    
    lr_next_traj <- lrtest(mod_next, mod_traj)
    p_next_traj <- lr_next_traj$`Pr(>Chisq)`[2]
    
    
    
    
  }
  
  
  
  
 
  
}



for(cellname in unique(df$unit)){
  celldf <- subset(df, unit == cellname)
  for(field in unique(celldf$field)){
    
    fdf <- subset(celldf, field == field)
    hierarchical_model(fdf)
    
  }
  
}
  
  
  

celldf <- subset(df, unit == 'TT1\\cl-maze1.1')
fdf <- subset(celldf, field =='0')

# Create Models ----



# Repetition Project
# 22-01-08 WH 
# Script to create GLM per field and test effects of factors 
# via post-hoc comparisons between levels within factor, via emmeans
# Time will be modeled as a FE or RE (using glmer) depending on further thought
#

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
#alleypath <- "E:\\Ratterdam\\R_data_repetition\\20220215-140515_superPopAlleyBehaviorResponse_1.5vfilt_PosInFieldNormedFR.csv"
alleypath <- "E:\\Ratterdam\\R_data_repetition\\2022-03-23_AlleySuperpopDirVisitFiltered.csv"

alleydf <- read.csv(alleypath,header=TRUE)


alleydf <- alleydf[is.finite(alleydf$Rate),]
alleydf <- alleydf[alleydf$Traversal=="True",]
alleydf <- alleydf[alleydf$Reward=="False",]


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

startTimeKnots <- 3

alpha <- 0.05

total_models_run <- 0
total_current <- 0 
total_next <- 0
total_previous <- 0

current_responsive <- c()
next_responsive <- c()
previous_responsive <- c()

# heres a page that lists refs for different choices of threshold
# https://quantifyinghealth.com/vif-threshold
vif_thresh = 5

repOrNot <- c()

for(o in c('V','H')){
  oriendf <- subset(alleydf, Orientation==o)
  
  for(fid in unique(oriendf$FieldID)){
    print(fid)
    
    field <- subset(oriendf, FieldID == fid)
    
    dirs <- unique(field$CurrDir)
    
    # for(dir in dirs){
    #   dfield <- field[field$CurrDir==dir,]
    #   try({
      # m <- glm(Rate + 1 ~ PrevDir + NextDir + ns(StartTimes,3),
      #          family='Gamma',
      #          data=dfield)
      
    try({
    
    m <- glm(Rate + 1 ~ CurrDir + PrevDir + NextDir + ns(StartTimes,3),
             family='Gamma',
             data=field)
    
    
    al <- alias(m)
    if(!("Complete" %in% names(al))){
    
    v <- vif(m) # will error out if there are any cases of complete MC
    
    #sometimes vif() returns a single value for each var,
    # sometimes 2 values (gvif and gvif^(1/(2*df)))
    if(length(v)==12){
      cvif <- v['CurrDir',3]
      pvif <- v['PrevDir',3]
      nvif <- v['NextDir',3]
    }
    else if(length(v)==4){
      cvif <- v['CurrDir']
      pvif <- v['PrevDir']
      nvif <- v['NextDir']
    }
    
    total_models_run <- total_models_run + 1

    #Previous Direction
    if(pvif < vif_thresh){
      total_previous <- total_previous + 1
      em_m <- emmeans(m,"PrevDir")
      pwc <- summary(pairs(em_m))
      sig <- pwc[1][pwc[6]<alpha]
      if(length(sig)>=1){
        previous_responsive <- c(previous_responsive, fid)
      }
    }
    
    #Current Direction
    if(cvif < vif_thresh){
    total_current <- total_current + 1
    em_m <- emmeans(m,"CurrDir")
    pwc <- summary(pairs(em_m))
    sig <- pwc[1][pwc[6]<alpha]
    if(length(sig)>=1){
      current_responsive <- c(current_responsive, fid)
    }
    }

    

    #Next Direction
    if(nvif < vif_thresh){
    total_next <- total_next + 1
    em_m <- emmeans(m,"NextDir")
    pwc = summary(pairs(em_m))
    sig <- pwc[1][pwc[6]<alpha]
    if(length(sig)>=1){
      next_responsive <- c(next_responsive, fid)
    }
    }
    
    
    
    }
    
    },silent=TRUE)
      
    }
    
  
  
}
print("Current")
print(length(current_responsive))
print(total_current)
print(length(current_responsive)/total_current)

print("Previous")
print(length(previous_responsive))
print(total_previous)
print(length(previous_responsive)/total_previous)

print("Next")
print(length(next_responsive))
print(total_next)
print(length(next_responsive)/total_next)




##### Break down by repeating or non
# 
# total_rep <- c()
# total_nonrep <- c()
# 
# for(cid in unique(alleydf$CellID)){
#   
#   if(unique(alleydf[alleydf$CellID==cid,'Repeating'])=="True"){
#     total_rep <- c(total_rep, cid)
#   }
#   else if(unique(alleydf[alleydf$CellID==cid,'Repeating'])=="False"){
#     total_nonrep <- c(total_nonrep, cid)
#   }
# }
# 
# repResponsive <- c()
# nonrepResponsive <- c()
# 
# for(fid in current_responsive){
#   repStatus <- unique(alleydf[alleydf$FieldID==fid,"Repeating"])
#   if(repStatus=="True"){
#     repResponsive <- c(repResponsive, fid)
#   }
#   else if(repStatus == "False"){
#     nonrepResponsive <- c(nonrepResponsive, fid)
#     
#   }
#   
# }
# 
# print(length(repResponsive)/length(total_rep))
# print(length(nonrepResponsive)/length(total_nonrep))
    



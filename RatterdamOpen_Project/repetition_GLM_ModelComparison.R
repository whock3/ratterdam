## (G)LMM Models and Model Comparisons For Testing Effect of Behavioral Variables
## On Individual Field Firing.
## WH Repetition Project 21-10-18

library(lme4)
library(ggplot2)
library(lmtest)
library(stringr)
library(emmeans)
library(tidyr)
library(ggpubr)
library(splines)

# Load data and recode factors

# 211207 has 30% field overlap threshold and revised traversal thresholds as of 12-6-21
#alleypath <- "E:\\Ratterdam\\R_data_repetition\\211207_AlleySuperpopDirVisitFiltered.csv"

# 211206 has 45% field overlap threshold and revised traversal thresholds as of 12-6-21
#alleypath <- "E:\\Ratterdam\\R_data_repetition\\211206_AlleySuperpopDirVisitFiltered.csv"


# 211208 has 10-45% overlap. I.e. we want just the penumbra. Traversal thresholds same
# except length thresh is 0.1 not 0.25 (because field pieces so small you dont want
# to miss everything)
#alleypath <- "E:\\Ratterdam\\R_data_repetition\\211208_AlleySuperpopDirVisitFiltered.csv"


# 211210 has 30% field overlap threshold and slightly looser traversal thresholds 
# alleypath <- "E:\\Ratterdam\\R_data_repetition\\211220_AlleySuperpopDirVisitFiltered.csv"

alleypath <- "E:\\Ratterdam\\R_data_repetition\\220218_AlleySuperpopDirVisitFiltered.csv"



alleydf <- read.csv(alleypath,header=TRUE)

alleydf$PrevDir <- as.factor(alleydf$PrevDir)
alleydf$CurrDir <- as.factor(alleydf$CurrDir)
alleydf$NextDir <- as.factor(alleydf$NextDir)
alleydf$RetroEgo <- as.factor(alleydf$RetroEgo)
alleydf$ProspEgo <- as.factor(alleydf$ProspEgo)
alleydf$Repeating <- as.factor(alleydf$Repeating)
alleydf$Traversal <- as.factor(alleydf$Traversal)
alleydf$Reward <- as.factor(alleydf$Reward)
alleydf$CellID <- as.factor(alleydf$CellID)
alleydf$FieldNum <- as.factor(alleydf$FieldNum)
alleydf$FieldID <- as.factor(alleydf$FieldID)
alleydf$Alleys <- as.factor(alleydf$Alleys)
alleydf$NumFields <- as.numeric(alleydf$NumFields)

m1_rmse <- c()
m2_rmse <- c()
m3_rmse <- c()
m4_rmse <- c()
m5_rmse <- c()

m1_aic <- c()
m2_aic <- c()
m3_aic <- c()
m4_aic <- c()
m5_aic <- c()


lrCurr_pvals <- c() # keep all pvalues from lrtest(base, base+cd) regardless
                    # if thats the best model because want to plot rmse colored
                    # by whether cd was helpful for fig 3 sfn2021. fig 5 then gets
                    # into whats the best model. 

msigCD <- 0 # how many models have better fit w current dir, use to see if
            # num models with better fit of P+C+N is sig  by binomial 

startTimeKnots = 3

nshuffle <- 1
shuffle <- FALSE

shuff_rep_CD <- c()
shuff_rep_ND <- c()
shuff_rep_PD <- c()
shuff_rep_Full <- c()

shuff_nonrep_CD <- c()
shuff_nonrep_ND <- c()
shuff_nonrep_PD <- c()
shuff_nonrep_Full <- c()

sigP <- c() # save whether CD is sig, for figures

if(shuffle==FALSE){
  for(s in 1:nshuffle){
    
    
    repFull <- 0
    repCD <- 0
    repPD <- 0
    repND <- 0
    
    nonrepFull <- 0
    nonrepCD <- 0
    nonrepPD <- 0
    nonrepND <- 0
    
    actualRunRep <- 0
    actualRunNonrep <- 0
    
    repOrNot <- c() # keep track of which fields are repeating for slicing purposes 
    
    
    print(s)
    for(o in c('V','H')){
      oriendf <- subset(alleydf, Orientation==o)
      
      for(fid in unique(oriendf$FieldID)){
        
        field_ <- subset(oriendf, FieldID == fid)
        
        n <- sample(nrow(field_))
        shuff.field <- data.frame(field_)
        shuff.field$CurrDir <- shuff.field$CurrDir[n]
        shuff.field$PrevDir <- shuff.field$PrevDir[n]
        shuff.field$NextDir <- shuff.field$NextDir[n]
        
        if(shuffle==TRUE){
          field <- shuff.field
          
        }
        else if(shuffle==FALSE){
          field <- field_
        }
        
        
        mrun <- try({
          # 21-10-24 changed models to include main effect of time expressed as 
          # ns(StartTimes,nknots). Before the m1 ws y ~ 1 and other models were same minus that time term
          # m1 <- glm(Rate+1 ~ ns(StartTimes,startTimeKnots), family='Gamma', data=field)
          # m2 <- glm(Rate+1 ~ CurrDir + ns(StartTimes,startTimeKnots), family = 'Gamma', data=field)
          # m3 <- glm(Rate+1 ~ PrevDir + CurrDir + ns(StartTimes,startTimeKnots), family='Gamma',data=field)
          # m4 <- glm(Rate+1 ~ CurrDir + NextDir + ns(StartTimes,startTimeKnots), family='Gamma',data=field)
          # m5 <- glm(Rate+1 ~ PrevDir + CurrDir + NextDir + ns(StartTimes,startTimeKnots), family='Gamma',data=field)
          
          m1 <- glm(Rate+1 ~ ns(StartTimes,startTimeKnots), family='Gamma', data=field)
          m2 <- glm(Rate+1 ~ CurrDir + ns(StartTimes,startTimeKnots), family = 'Gamma', data=field)
          m3 <- glm(Rate+1 ~ NextDir + CurrDir + ns(StartTimes,startTimeKnots), family='Gamma',data=field)
          m4 <- glm(Rate+1 ~ CurrDir + PrevDir + ns(StartTimes,startTimeKnots), family='Gamma',data=field)
          m5 <- glm(Rate+1 ~ PrevDir + CurrDir + NextDir + ns(StartTimes,startTimeKnots), family='Gamma',data=field)
          
          m1_rmse <- c(m1_rmse, sqrt(mean((field$Rate-m1$fitted.values)^2)))
          m2_rmse <- c(m2_rmse, sqrt(mean((field$Rate-m2$fitted.values)^2)))
          m3_rmse <- c(m3_rmse, sqrt(mean((field$Rate-m3$fitted.values)^2)))
          m4_rmse <- c(m4_rmse, sqrt(mean((field$Rate-m4$fitted.values)^2)))
          m5_rmse <- c(m5_rmse, sqrt(mean((field$Rate-m5$fitted.values)^2)))
          
          m1_aic <- c(m1_aic, m1$aic)
          m2_aic <- c(m2_aic, m2$aic)
          m3_aic <- c(m3_aic, m3$aic)
          m4_aic <- c(m4_aic, m4$aic)
          m5_aic <- c(m5_aic, m5$aic)
 
          lrCurr <- lrtest(m1,m2)
          lrPrev <- lrtest(m2,m3)
          lrNext <- lrtest(m2,m4)
          lrFullA <- lrtest(m3,m5)
          lrFullB <- lrtest(m4,m5)
          
          lrCurr_pvals <- c(lrCurr_pvals, lrCurr[2,"Pr(>Chisq)"])
          
          if((lrCurr[2,"Pr(>Chisq)"] < 0.05)&(!is.nan(lrCurr[2,"Pr(>Chisq)"]))){
            
            sigP <- c(sigP, 1)
          }
          else{
            sigP <- c(sigP, 0)
          }
          
          bfAdj <- 5 # adjust for multiple comparisons via Bonferroni correction
          nfields <- as.numeric(unique(as.character(field$NumFields)))
          # Case checks for repeating fields
          if(unique(field$Repeating)=='True'){
            actualRunRep <- actualRunRep + 1
            repOrNot <- c(repOrNot, TRUE)
              if((lrFullA[2,"Pr(>Chisq)"]<(0.05/bfAdj))&(lrFullB[2,"Pr(>Chisq)"]<(0.05/bfAdj))){
                repFull <- repFull + 1
              }
              else if((lrPrev[2,"Pr(>Chisq)"]<(0.05/bfAdj))&(lrNext[2,"Pr(>Chisq)"]>(0.05/bfAdj))){
                repPD <- repPD + 1
                
              }
              else if((lrPrev[2,"Pr(>Chisq)"]>(0.05/bfAdj))&(lrNext[2,"Pr(>Chisq)"]<(0.05/bfAdj))){
                repND <- repND + 1 
              }
              else if(lrCurr[2,"Pr(>Chisq)"]<(0.05/bfAdj)){
                repCD <- repCD + 1 
                
              }
              
          }
          
          
          # Case checks for nonrepeating fields
          else if((unique(field$Repeating)=='False')){
            actualRunNonrep <- actualRunNonrep + 1
            repOrNot <- c(repOrNot, FALSE)
              if((lrFullA[2,"Pr(>Chisq)"]<(0.05/bfAdj))&(lrFullB[2,"Pr(>Chisq)"]<(0.05/bfAdj))){
                nonrepFull <- nonrepFull + 1
              }
              else if((lrPrev[2,"Pr(>Chisq)"]<(0.05/bfAdj))&(lrNext[2,"Pr(>Chisq)"]>(0.05/bfAdj))){
                nonrepPD <- nonrepPD + 1

              }
              else if((lrPrev[2,"Pr(>Chisq)"]>(0.05/bfAdj))&(lrNext[2,"Pr(>Chisq)"]<(0.05/bfAdj))){
                nonrepND <- nonrepND + 1 
              }
              else if(lrCurr[2,"Pr(>Chisq)"]<(0.05/bfAdj)){
                nonrepCD <- nonrepCD + 1 
              }
              
            
            
          }
          
          
        },silent=FALSE) 
      }
    }
    
    shuff_rep_CD <- c(shuff_rep_CD, repCD/actualRunRep)
    shuff_rep_ND <- c(shuff_rep_ND, repND/actualRunRep)
    shuff_rep_PD <- c(shuff_rep_PD, repPD/actualRunRep)
    shuff_rep_Full <- c(shuff_rep_Full, repFull/actualRunRep)
    
    shuff_nonrep_CD <- c(shuff_nonrep_CD, nonrepCD/actualRunNonrep)
    shuff_nonrep_ND <- c(shuff_nonrep_ND, nonrepND/actualRunNonrep)
    shuff_nonrep_PD <- c(shuff_nonrep_PD, nonrepPD/actualRunNonrep)
    shuff_nonrep_Full <- c(shuff_nonrep_Full, nonrepFull/actualRunNonrep)
    
  }
}


total <- actualRunRep + actualRunNonrep
print(sprintf("Repeating fields best model is P+C+N: %s (%s/%s)",repFull/actualRunRep,repFull,actualRunRep))
print(sprintf("Repeating fields best mode is Current Dir: %s (%s/%s)",repCD/actualRunRep,repCD,actualRunRep))
print(sprintf("Repeating fields best model is Previous Dir: %s (%s/%s)",repPD/actualRunRep,repPD,actualRunRep))
print(sprintf("Repeating fields best model is Next Dir: %s (%s/%s)",repND/actualRunRep,repND,actualRunRep))

print(sprintf("Non-repeating fields best model is P+C+N: %s (%s/%s)",nonrepFull/actualRunNonrep,nonrepFull,actualRunNonrep))
print(sprintf("Non-repeating fields best mode is Current Dir: %s (%s/%s)",nonrepCD/actualRunNonrep,nonrepCD,actualRunNonrep))
print(sprintf("Non-repeating fields best model is Previous Dir: %s (%s/%s)",nonrepPD/actualRunNonrep,nonrepPD,actualRunNonrep))
print(sprintf("Non-repeating fields best model is Next Dir: %s (%s/%s)",nonrepND/actualRunNonrep,nonrepND,actualRunNonrep))

counts <- c(repFull/actualRunRep,
            repCD/actualRunRep,
            repPD/actualRunRep,
            repND/actualRunRep,
            nonrepFull/actualRunNonrep,
            nonrepCD/actualRunNonrep,
            nonrepPD/actualRunNonrep,
            nonrepND/actualRunNonrep
            )
barnames <- c("Rep Full",
              "Rep CD",
              "Rep PD",
              "Rep ND",
              "Nonrep Full",
              "Nonrep CD",
              "Nonrep PD",
              "Nonrep ND")

barplot(counts, names.arg=barnames,
        col=c('red','red','red','red','grey','grey','grey','grey'),
        main='Repeating versus Single-Fielded Non-repeating')
abline(h = 0.05, col="black", lwd=3, lty=2)












# Repetition Models
# WH Oct 24 2021
# Time interaction models 

library(lme4)
library(ggplot2)
library(lmtest)
library(stringr)
library(emmeans)
library(tidyr)
library(ggpubr)
library(splines)

# Load data and recode factors
alleypath <- "E:\\Ratterdam\\R_data_repetition\\211005_AlleySuperpopDirVisitFiltered.csv"
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



###
### Time * Current Direction
###
savepath <- 'E:\\Ratterdam\\temp\\TimeCD_Model\\'

sigInters <- c()
set.seed(123)
for(o in c("V","H")){
  oriendf <- subset(alleydf, Orientation==o)
  
  for(fid in unique(oriendf$FieldID)){
    u<-try({
    field <- subset(oriendf, FieldID==fid)
    mod <- glm(Rate + 1 ~ CurrDir*ns(StartTimes,3),family='Gamma',data=field)
    s<-summary(mod)
    anysig <- FALSE
    
    #gen labels for saving, relative field number and repeating or not
    relfid <- unique(field$FieldNum)
    repeating <- unique(field$Repeating)
    
    # this checks to see if a given interaction term is present in the coefficients.
    # a given coef may not be there due to, e.g. singularities
    # obviously hardcoded for CurrDir*ns(StartTimes,3)
    if(sprintf("%s:ns(StartTimes, 3)1",row.names(s$coefficients)[2]) %in% row.names(s$coefficients)){
      pval <- s$coefficients[sprintf("%s:ns(StartTimes, 3)1",row.names(s$coefficients)[2]),4]
      if((!is.nan(pval)&(pval<(0.05/3)))){
        sigInters <- c(sigInters, fid)
        anysig <- TRUE
      }
      
    }
    if(sprintf("%s:ns(StartTimes, 3)2",row.names(s$coefficients)[2]) %in% row.names(s$coefficients)){
      pval <- s$coefficients[sprintf("%s:ns(StartTimes, 3)2",row.names(s$coefficients)[2]),4]
      if((!is.nan(pval)&(pval<(0.05/3)))){
        sigInters <- c(sigInters, fid)
        anysig <- TRUE
      }
      
    }
    if(sprintf("%s:ns(StartTimes, 3)3",row.names(s$coefficients)[2]) %in% row.names(s$coefficients)){
      pval <- s$coefficients[sprintf("%s:ns(StartTimes, 3)3",row.names(s$coefficients)[2]),4]
      if((!is.nan(pval)&(pval<(0.05/3)))){
        sigInters <- c(sigInters, fid)
        anysig <-TRUE
      }
    }
    
    
    field$fit <- predict(mod,newdata=field,re.form=NA,type='response')
    title_rat <- unique(oriendf[oriendf$FieldID==fid,]$Rat)
    title_day <- unique(oriendf[oriendf$FieldID==fid,]$Day)
    title_cell <- unique(oriendf[oriendf$FieldID==fid,]$CellName)
    p <- ggplot(data=field,aes(x=((StartTimes/1e6)/60)-((StartTimes[1]/1e6)/60),col=CurrDir))+
      geom_point(aes(y=Rate))+
      geom_line(aes(y=fit, col=CurrDir),size=1.5)+
      ggtitle(sprintf("%s %s %s Field %s  %s  Section. repeating = %s, sig interaction = %s", title_rat, 
                                                                                              title_day, 
                                                                                              title_cell, 
                                                                                              relfid, 
                                                                                              o,
                                                                                              repeating, 
                                                                                              anysig
              ))+
      theme(plot.title = element_text(size = 10, face = "bold"))+
      labs(x="Time (min)")
    
    ggsave(sprintf("%s%s_%s_%s_%s_%s.png",savepath,title_rat,title_day,str_replace(title_cell,"\\\\","-"),relfid,o))
    
    },silent=FALSE)
    }
    
}



##
## Time alone
##


sigTime <- c()
set.seed(123)
for(o in c("V","H")){
  oriendf <- subset(alleydf, Orientation==o)
  
  for(fid in unique(oriendf$FieldID)){
    u<-try({
      field <- subset(oriendf, FieldID==fid)
      mod <- glm(Rate + 1 ~ ns(StartTimes,3),family='Gamma',data=field)
      s<-summary(mod)
      anysig <- FALSE
      
      print(fid)
      print(s)
      
      
    },silent=TRUE)
  }
}
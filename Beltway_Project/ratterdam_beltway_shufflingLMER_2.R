# Shuffling LMER Analysis
# WH Mid March 2021
# Run LMER analysis on data with shuffled texture lables within an alley
# Compare distribution of effects (e.g. coefficients) to empirical

# imports
library(lme4)
library(splines)
library(ggplot2)
library(lmtest)
library(stringr)


setwd("E:\\UserData\\Documents\\GitHub\\ratterdam\\")
source("E:\\UserData\\Documents\\GitHub\\ratterdam\\Beltway_Project\\cicheck.R")
source("E:\\UserData\\Documents\\GitHub\\ratterdam\\Beltway_Project\\glmer_fx.R")
source("E:\\UserData\\Documents\\GitHub\\ratterdam\\Beltway_Project\\ratterdam_RunLMER.R")

save <- FALSE # toggle to save multipage pdfs and wald csvs
database <- "E:\\Ratterdam\\R_data_beltway\\"
datapath <- sprintf("%s%s",database,"20210602-131553_lmerPassingPop_1.5vfilt_0.5stepsmooth_24bins_2R_3qual.csv")
df <- read.csv(datapath,header=TRUE)
df$texture <- as.factor(df$texture)
df$reward <- as.factor(df$reward)

nshuffles <- 1000

shuffResults_cells <- integer(nshuffles) # initialize 0s of size=nshuffles. These will be incremented according to
                                          # how many cells per shuffle (i.e index in vector) pass the null lmer test
shuffResults_cellsName <- vector("list",100)

for (cellname in unique(df$name)){
  
  print(cellname)

  # load cell and define key variables
  celldf <- df[df$name==cellname,]
  
  nalleys <- length(unique(celldf$alley))
  

  #loop over alleys 
  for (a in unique(celldf$alley)){
    u <- try({
    print(a)
    celldf_a <- celldf[celldf$alley==a,]
    celldf_a_copy <- data.frame(celldf_a) # used to shuffle, keep original trial structure in celldf_a 
    
    
    celldf_a <- lmer_routine(celldf_a, nalleys) # nalleys is number of alleys w activity and used for bonferroni corr
    #empirical <- calcAlleyMaxSD(celldf_a, a)
    
    
    # Create trial list
    ntrials <- length(unique(celldf_a$trial))
    trialList <- data.frame(texture=factor())
    for (i in 0:ntrials-1){
      stim <- unique(celldf_a[celldf_a$trial==i,]$texture)
      trialList <- rbind(trialList, paste(stim))
    }
    
    shuffresults <- 0

    # shuffle loop
    for (i in 1:nshuffles){
      
      shuffledTrials <- data.frame(trialList[sample(1:ntrials),])
      
      for (j in 1:ntrials){
        celldf_a_copy[celldf_a$trial==j-1,]$texture <- shuffledTrials[j,]

      }
  
        
        celldf_a_copy <- lmer_routine(celldf_a_copy, nalleys)
        celldf_a_copy_no_r <- data.frame(celldf_a_copy[celldf_a_copy$reward=="0",])
        
        outcome <- checkAlleyNonoverlap(celldf_a_copy_no_r,a)
        if(outcome==TRUE){
          shuffResults_cells[i] <- shuffResults_cells[i] + 1
          if(i%%10==0){
            shuffResults_cellsName[[i/10]] <- c(shuffResults_cellsName[[i/10]],sprintf("%s_%s", cellname,a))
          }
          
        }

    }
    #print(shuffresults/nshuffles)
    # pct95 <- quantile(sort(shuffledData),0.95)
    # title <- sprintf("Cell %s Alley %s",cellname,a)
    # hist(shuffledData, main=title, xlab="Max SD Between Conditions", ylab="Frequency")
    # abline(v=c(empirical, pct95), col=c("red", "black"))
    
    },silent=FALSE)
    
  }


}










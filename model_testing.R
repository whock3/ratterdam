## Thru 10/28/20 beofre SZ meeting. 

#Params for whole analysis. Like number spatial bins,etc
#nbin <-12 

# cellID = 1
# alleyID = 5
# name <- unique(df[df$cell==cellID,]$name)
# print(name)
# celldf <- subset(df, cell == cellID | alley == alleyID)
# mod <- glmer_alley(cellID, alleyID, df)
# celldf$fit <- exp(predict(mod, newdata=celldf,re.form=NA))

# dfA <- celldf[celldf$texture=='A',]
# dfB <- celldf[celldf$texture=='B',]
# dfC <- celldf[celldf$texture=='C',]
# 
# seA <- 1:nbin # this is not how to preallocate properly. but 'right' way didnt let me do vector math so idk
# seB <- 1:nbin
# seC <- 1:nbin
# 
# 
# for (bin in seq(1,nbin)){
#   subA <- dfA[dfA$spatialBin==bin-1,]
#   subB <- dfB[dfB$spatialBin==bin-1,]
#   subC <- dfC[dfC$spatialBin==bin-1,]
#   
#   seA[bin] <- sum((subA$rate-subA$fit)^2)/(nrow(subA)-2)
#   seB[bin] <- sd(subB)/sqrt(length(subB))
#   seC[bin] <- sd(subC)/sqrt(length(subC))
# }
# 
# 
# (g <- ggplot(NULL, aes(spatialBin, fit))+
#       geom_line(data=dfA, col='red')+
#       geom_line(data=dfB, col='blue')+
#       geom_line(data=dfC, col='green')+
#       geom_ribbon(data=dfA,aes(y = fit, ymin=fit-1.96*seA,
#                       ymax=fit+1.96*seA),alpha=0.1,fill='red')+
#       geom_ribbon(data=dfB,aes(y = fit, ymin=fit-1.96*seB,
#                                ymax=fit+1.96*seB),alpha=0.1,fill='blue')+
#       geom_ribbon(data=dfC,aes(y = fit, ymin=fit-1.96*seC,
#                                ymax=fit+1.96*seC),alpha=0.1,fill='green')+
#       ggtitle(sprintf("%s alley %s", name, alleyID ))
# 
# )
# print(g)

## 10/29/20

# imports
library(lme4)
library(splines)
library(ggplot2)
library(lmtest)
setwd("E:\\UserData\\Documents\\GitHub\\ratterdam\\")
source("glmer_fx.R")

# Choose where figs will go
df <- read.csv("E:\\Ratterdam\\R_data\\R781BRD3.csv",header=TRUE)
figbasepath <- "E:\\Ratterdam\\R_data\\results\\"

for(cellID in unique(df$cell)){

  
  # Select and set up data
  celldf <- subset(df, cell == cellID)
  celldf$rate <- log(celldf$rate+1)
  cellname <- unique(df[df$cell==cellID,]$name)
  alleys <- unique(celldf$alley)
  
  pdf(paste(figbasepath,"R781BRD3_",cellID,".pdf",sep=""),onefile=TRUE)
  
  for(alleyID in alleys){
    
    print(alleyID)
    
    celldf <- subset(df, cell == cellID & alley == alleyID)
    
    # run models
    modm <- lmer_alley_main(celldf)
    modi <- lmer_alley_int(celldf)
    modn <- lmer_alley_none(celldf)
    
    # Reorder stimulus factors to make other comparisons
    celldf$texture <- factor(celldf$texture, levels = c("B", "A", "C"))
    modreord <- lmer_alley_int(celldf)
    celldf$texture <- factor(celldf$texture, levels = c("A","B","C"))
    
    
    # likelihood test modi vs modm and modn vs modm
    likeli_main <- lrtest(modn, modm)
    likeli_int <- lrtest(modm, modi)
    plrmain <- format(round(likeli_main$`Pr(>Chisq)`, 2), nsmall=3)[2]
    plrint <- format(round(likeli_int$`Pr(>Chisq)`, 2), nsmall=3)[2]
    
    # fit for model with int
    celldf$fit <- predict(modi, newdata=celldf, re.form=NA)
    
    # for interaction model, Wald CI. txt A as default
    ci <- confint(modi, method='Wald', level=0.995)
    ci <- ci[-c(1,2),]
    fe <- fixef(modi)
    cidf <- data.frame("low"=ci[,1], "up"=ci[,2],"fe"=fe)
    cidf <- cidf[-c(1),]
    
    # do CI on model with B as default and concat 
    cireord <- confint(modreord, method='Wald', level=0.995)
    cireord <- cireord[-c(1,2),]
    fereord <- fixef(modreord)
    cireorddf <- data.frame("low"=cireord[,1], "up"=cireord[,2],"fe"=fereord)
    ciall <- rbind(cidf, cireorddf[c("textureA", "ns(spatialBin, 3)1:textureC","ns(spatialBin, 3)2:textureC","ns(spatialBin, 3)3:textureC"),])
    
    
    # ggplot of wald CIs
    x = seq(1,nrow(ciall))
    waldplot <- ggplot(data=ciall, aes(x=x,y=fe))+
      geom_point()+
      geom_errorbar(aes(ymin=low,ymax=up))+
      scale_x_discrete(limits=row.names(ciall))+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
      geom_hline(yintercept=0)+
      ggtitle(sprintf("Cell %s Alley %s 99.5pct Wald CI",cellname, alleyID))
    print(waldplot)
    
    # calc CI of fits. CI = 95%
    # set up design matrix then multiply mat*var-cov mat * mat-1 to get var
    # then sqrt and mult by crit value to get CI of given pct 
    Designmat <- model.matrix(rate ~ ns(spatialBin, 3)*texture + reward, celldf)
    predvar <- diag(Designmat %*% vcov(modi) %*% t(Designmat))
    celldf$fitCI <- sqrt(predvar)*2.58 #99.5
    
    # ggplot of splines with 95% CI of fits
    celldf <- celldf[celldf$reward=="0",]
    p <- ggplot(data=celldf, aes(x=spatialBin))+
      geom_line(aes(y=fit, color=texture,linetype=factor(reward)))+
      geom_ribbon(aes(y=fit, ymin=fit-fitCI, ymax=fit+fitCI, fill=texture), alpha=0.2)+
      ggtitle(sprintf("Cell %s, Alley %s", cellname, alleyID))
    print(p)
  }
  
  dev.off()

}


  
  
  
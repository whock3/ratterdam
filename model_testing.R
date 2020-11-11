# imports
library(lme4)
library(splines)
library(ggplot2)
library(lmtest)
library(stringr)

setwd("E:\\UserData\\Documents\\GitHub\\ratterdam\\")
source("glmer_fx.R")

# Read in data
code <- "R886BRD2"
database <- "E:\\Ratterdam\\R_data\\csv_vfilt2_p5stepsmooth_12bins\\"
datapath <- sprintf("%s%s%s",database,code,".csv")
df <- read.csv(datapath,header=TRUE)

# Select output, create timestamp 
figbasepath <- "E:\\Ratterdam\\R_data\\graphs\\"
ts <- str_replace(Sys.time()," ","_")
ts <- str_replace_all(ts, ":", "_")


if(code=="R859BRD5"){
  
  df<-df[!df$alley==3,]
  df<-df[!df$alley==5,]
  df<-df[!df$alley==7,]
}



wdf <- data.frame(matrix(ncol=5))
colnames(wdf) <- c("cellname", "alleyID", "low", "up", "fe")

for(cellID in unique(df$cell)){

  
  # Select and set up data
  celldf <- subset(df, cell == cellID)
  celldf$rate <- log(celldf$rate+1)
  cellname <- unique(df[df$cell==cellID,]$name)
  alleys <- unique(celldf$alley)
  nalleys <-length(alleys)
  
  bonf <- 0.05/(9*nalleys) # 9 comparisons in alley (3 stim * 3 texture) * # alleys
  cipct <- 1-bonf
  z <- qnorm(cipct)
  
  pdf(paste(figbasepath,ts,"_",code,"_",cellID,".pdf",sep=""),onefile=TRUE)
  
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
    ci <- confint(modi, method='Wald', level=cipct)
    ci <- ci[-c(1,2),]
    fe <- fixef(modi)
    cidf <- data.frame("low"=ci[,1], "up"=ci[,2],"fe"=fe)
    cidf <- cidf[-c(1),]
    
    # do CI on model with B as default and concat 
    cireord <- confint(modreord, method='Wald', level=cipct)
    cireord <- cireord[-c(1,2),]
    fereord <- fixef(modreord)
    cireorddf <- data.frame("low"=cireord[,1], "up"=cireord[,2],"fe"=fereord)
    ciall <- rbind(cidf, cireorddf[c("textureC", "ns(spatialBin, 3)1:textureC","ns(spatialBin, 3)2:textureC","ns(spatialBin, 3)3:textureC"),])
    ndf <- data.frame(ciall,cellname,alleyID)
    wdf <- rbind(wdf, ndf)
    
    # ggplot of wald CIs
    x = seq(1,nrow(ciall))
    waldplot <- ggplot(data=ciall, aes(x=x,y=fe))+
      geom_point()+
      geom_errorbar(aes(ymin=low,ymax=up))+
      scale_x_discrete(limits=row.names(ciall))+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
      geom_hline(yintercept=0)+
      ggtitle(sprintf("Cell %s Alley %s %s Wald CI",cellname, alleyID, cipct))
    print(waldplot)
    
    # calc CI of fits. CI = 95%
    # set up design matrix then multiply mat*var-cov mat * mat-1 to get var
    # then sqrt and mult by crit value to get CI of given pct 
    Designmat <- model.matrix(rate ~ ns(spatialBin, 3)*texture + reward, celldf)
    predvar <- diag(Designmat %*% vcov(modi) %*% t(Designmat))
    celldf$fitCI <- sqrt(predvar)*z
    
    # ggplot of splines with CI of fits
    celldf <- celldf[celldf$reward=="0",]
    p <- ggplot(data=celldf, aes(x=spatialBin))+
      geom_line(aes(y=fit, color=texture,linetype=factor(reward)))+
      geom_ribbon(aes(y=fit, ymin=fit-fitCI, ymax=fit+fitCI, fill=texture), alpha=0.2)+
      ggtitle(sprintf("Cell %s, Alley %s", cellname, alleyID))
    print(p)
  }
  
  dev.off()

}

# Write Wald CI data to csv
dpath <- "E:\\Ratterdam\\R_data\\"
write.csv(wdf,paste(dpath,ts,"_",code,"_wald.csv",sep=""),row.names=TRUE)
  
  
  
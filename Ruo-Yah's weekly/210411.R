#  Ratterdam Repetition Project
#  Assessing relationship between rate and turns
#  LMER Analysis

library(lme4)
library(splines)
library(ggplot2)
library(lmtest)
library(stringr)

# define data paths and load in data
path <- "C:\\Users\\Ruo-Yah Lai\\Desktop\\My folder\\College\\Junior\\K lab research\\Graphs\\LMER\\20210321-163252_R781D2_1vfilt_.csv"
df <- read.csv(path,header=TRUE)

# Select output, create timestamp 
save <- TRUE
figbasepath <- "C:\\Users\\Ruo-Yah Lai\\Desktop\\My folder\\College\\Junior\\K lab research\\Graphs\\LMER\\"
ts <- str_replace_all(Sys.time(),"-","")
ts <- str_replace(ts," ","-")
ts <- str_replace_all(ts, ":", "")

exp <- "R781D2"

if(save==TRUE){
  pdf(paste(figbasepath,ts,"_","nextTurnEgo","_",exp,".pdf",sep=""),onefile=TRUE)
}

for(cellID in unique(df$cell)){
  print(cellID)
  if (cellID == "TT3\\cl-maze1.4"){
    next
  }
  cellname <- cellID
  celldf <- subset(df, cell == cellname)
  
  # get rid of epoch 0, it's a py artefact from bisecting approach to finding 
  # epoch that needs to be fixed
  celldf <- celldf[!celldf$epoch==0,]
  celldf <- celldf[!celldf$nextTurnEgo==0,]
  
  #mod <- lmer(rate ~ ns(nextTurnEgo,nsplineknots):epoch + (1+ns(direction,nsplineknots)|field) + (1+epoch|field), data=celldf)
  #mod <- lmer(rate ~ as.factor(prevTurnAllo) + as.factor(nextTurnAllo) + (1+as.factor(prevTurnAllo)|field) + (1+as.factor(nextTurnAllo)|field), data=celldf)
  #celldf$fit <- predict(mod, newdata=celldf, re.form=NA)
  
  #Designmat <- model.matrix(rate ~ ns(direction,nsplineknots):epoch + (ns(direction,nsplineknots)|field) + (epoch|field), celldf)
  #Designmat <- model.matrix(rate ~ as.factor(prevTurnAllo) + as.factor(nextTurnAllo), data=celldf)
  #predvar <- diag(Designmat %*% vcov(mod) %*% t(Designmat))
  #celldf$fitCI <- sqrt(predvar)*z
  
  # interaction plot w spline fits
  p <- ggplot(celldf, aes(x=nextTurnEgo, y=rate))+
    geom_violin(aes(x=as.factor(nextTurnEgo), y=rate))+
    geom_boxplot(aes(x=as.factor(nextTurnEgo), y=rate), width=0.05)+
    ggtitle(sprintf("%s %s", exp, cellname))
    #geom_point(size=2,alpha=0.5)
  #ggplot(celldf, aes(x=as.factor(prevTurnAllo), y=rate, fill=as.factor(prevTurnEgo)))
  
  print(p)
  
  
  #individual plots by group 
  (mm_plot <- ggplot(celldf, aes(x = nextTurnEgo, y = rate)) +
      facet_wrap(~field, nrow=2) +   # a panel for each mountain range
      geom_violin(aes(x=as.factor(nextTurnEgo), y=rate))+
      geom_boxplot(aes(x=as.factor(nextTurnEgo), y=rate), width=0.1)+
      #geom_point(alpha = 0.5)
      ggtitle(sprintf("%s %s", exp, cellname))
  )
  print(mm_plot)
}

if(save==TRUE){
  dev.off()
}
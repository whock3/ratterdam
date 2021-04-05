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
  pdf(paste(figbasepath,ts,"_","alloTurn","_",exp,".pdf",sep=""),onefile=TRUE)
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
  celldf <- celldf[!celldf$prevTurnAllo==0,]
  celldf <- celldf[!celldf$nextTurnAllo==0,]
  
  #define params
  nsplineknots <- 6
  nepochs <- length(unique(celldf$epoch))
  
  
  bonf <- 0.17/(nepochs*nsplineknots) 
  cipct <- 1-bonf
  z <- qnorm(cipct)
  
  #mod <- lmer(rate ~ ns(nextTurnEgo,nsplineknots):epoch + (1+ns(direction,nsplineknots)|field) + (1+epoch|field), data=celldf)
  mod <- lmer(rate ~ as.factor(prevTurnAllo) + as.factor(nextTurnAllo) + (1+as.factor(prevTurnAllo)|field) + (1+as.factor(nextTurnAllo)|field), data=celldf)
  celldf$fit <- predict(mod, newdata=celldf, re.form=NA)
  
  #Designmat <- model.matrix(rate ~ ns(direction,nsplineknots):epoch + (ns(direction,nsplineknots)|field) + (epoch|field), celldf)
  Designmat <- model.matrix(rate ~ as.factor(prevTurnAllo) + as.factor(nextTurnAllo), data=celldf)
  predvar <- diag(Designmat %*% vcov(mod) %*% t(Designmat))
  celldf$fitCI <- sqrt(predvar)*z
  
  # interaction plot w spline fits
  p <- ggplot(celldf, aes(x=prevTurnAllo, y=rate, color=as.factor(nextTurnAllo)))+
    geom_line(aes(prevTurnAllo, fit, group=as.factor(nextTurnAllo)))+
    geom_ribbon(aes(y=fit, ymin=fit-fitCI, ymax=fit+fitCI, fill=as.factor(nextTurnAllo)), alpha=0.2)+
    ggtitle(sprintf("%s %s", exp, cellname))+
    geom_point(size=2,alpha=0.5)
  
  print(p)
  
  
  #individual plots by group 
  (mm_plot <- ggplot(celldf, aes(x = prevTurnAllo, y = rate, color=as.factor(nextTurnAllo))) +
      facet_wrap(~field, nrow=2) +   # a panel for each mountain range
      geom_point(alpha = 0.5) +
      geom_ribbon(aes(y=fit, ymin=fit-fitCI, ymax=fit+fitCI, fill=as.factor(nextTurnAllo)), alpha=0.2) + 
      geom_line(data = cbind(celldf, pred = predict(mod)), aes(y = pred), size = 1)  # adding predicted line from mixed model 
  )
  print(mm_plot)
}

if(save==TRUE){
  dev.off()
}
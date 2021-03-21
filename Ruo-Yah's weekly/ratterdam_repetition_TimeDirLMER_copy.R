#  Ratterdam Repetition Project
#  Assessing confound/relationship between directionality and time
#  LMER Analysis, Early March 2021
# Model Fr ~ direction (continuous) : epoch (factor) with each x|field 

library(lme4)
library(splines)
library(ggplot2)
library(lmtest)
library(stringr)

# define data paths and load in data
path <- "C:\\Users\\Ruo-Yah Lai\\Desktop\\My folder\\College\\Junior\\K lab research\\Graphs\\LMER\\20210321-163252_R781D2_1vfilt_.csv"
df <- read.csv(path,header=TRUE)

exp <- "R781D2"
cellname <- "TT9\\cl-maze1.4"
celldf <- subset(df, cell == cellname)

# get rid of epoch 0, it's a py artefact from bisecting approach to finding 
# epoch that needs to be fixed
celldf <- celldf[!celldf$epoch==0,]

#define params
nsplineknots <- 6
nepochs <- length(unique(celldf$epoch))

histdir <- vector(mode="numeric", length=nrow(celldf))

for (i in 1:length(histdir)){
  
  dir = celldf$direction[i]
  
  if (dir > 0 && dir <= 90){
    histdir[i] <- 1
  }
  else if (dir > 90 && dir <= 180){
    histdir[i] <- 2
  }
  else if (dir > 180 && dir <= 270){
    histdir[i] <- 3
  }
  else if (dir > 270 && dir <= 360){
    histdir[i] <- 4
  }
  else {}
  
}

celldf$histdir <- histdir


bonf <- 0.17/(nepochs*nsplineknots) 
cipct <- 1-bonf
z <- qnorm(cipct)

#mod <- lmer(rate ~ ns(prevTurnEgo,nsplineknots):epoch + (1+ns(direction,nsplineknots)|field) + (1+epoch|field), data=celldf)
mod <- lmer(rate ~ ns(direction,nsplineknots):epoch + (1+ns(direction,nsplineknots)|field) + (1+epoch|field), data=celldf)
celldf$fit <- predict(mod, newdata=celldf, re.form=NA)

#Designmat <- model.matrix(rate ~ ns(direction,nsplineknots):epoch + (ns(direction,nsplineknots)|field) + (epoch|field), celldf)
Designmat <- model.matrix(rate ~ ns(direction,nsplineknots):epoch, celldf)
predvar <- diag(Designmat %*% vcov(mod) %*% t(Designmat))
celldf$fitCI <- sqrt(predvar)*z

# interaction plot w spline fits
p <- ggplot(celldf, aes(x=direction, y=rate, color=as.factor(epoch)))+
  geom_line(aes(direction, fit, group=epoch))+
  geom_ribbon(aes(y=fit, ymin=fit-fitCI, ymax=fit+fitCI, fill=as.factor(epoch)), alpha=0.2)+
  ggtitle(sprintf("%s %s", exp, cellname))+
  geom_point(size=2,alpha=0.5)

print(p)


#individual plots by group 
#(mm_plot <- ggplot(celldf, aes(x = direction, y = rate, color = as.factor(epoch))) +
#    facet_wrap(~field, nrow=2) +   # a panel for each mountain range
#    geom_point(alpha = 0.5) +
#    geom_line(data = cbind(celldf, pred = predict(mod)), aes(y = pred), size = 1)  # adding predicted line from mixed model 
#)
#print(mm_plot)
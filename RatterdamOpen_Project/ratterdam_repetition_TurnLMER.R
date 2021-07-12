#  Ratterdam Repetition Project
#  Assessing role of turns on firing rate
#  LMER Analysis, Late June 2021

library(lme4)
library(splines)
library(ggplot2)
library(lmtest)
library(stringr)

# define data paths and load in data
path <- "E:\\Ratterdam\\R_data_repetition\\20210701-170138_R781D3_TurnLMER_1.5vfilt_.csv"
savepath <- "E:\\Ratterdam\\R_data_repetition\\210624-onward_turnLM\\"
code <- "R781D3"
df <- read.csv(path,header=TRUE)
ts <- str_replace(Sys.time()," ","_")
ts <- str_replace_all(ts, ":", "_")

df$field <- as.factor(df$field)

df$alloPrevPre <- as.factor(df$alloPrevPre)
df$egoPre <- as.factor(df$egoPre)
df$alloPrevPost <- as.factor(df$alloPrevPost)
df$alloCurrPre <- as.factor(df$alloCurrPre)
df$egoCurr <- as.factor(df$egoCurr)
df$alloCurrPost <- as.factor(df$alloCurrPost)
df$alloNextPre <- as.factor(df$alloNextPre)
df$egoNext <- as.factor(df$egoNext)
df$alloNextPost <-as.factor(df$alloNextPost)

df$alloPrevPre <- plyr::revalue(df$alloPrevPre, c("1"="N","2"="E","3"="S","4"="W"))
df$alloPrevPost <- plyr::revalue(df$alloPrevPost, c("1"="N","2"="E","3"="S","4"="W"))
df$alloCurrPre <- plyr::revalue(df$alloCurrPre, c("1"="N","2"="E","3"="S","4"="W"))
df$alloCurrPost <- plyr::revalue(df$alloCurrPost, c("1"="N","2"="E","3"="S","4"="W"))
df$alloNextPre <- plyr::revalue(df$alloNextPre, c("1"="N","2"="E","3"="S","4"="W"))
df$alloNextPost <- plyr::revalue(df$alloNextPost, c("1"="N","2"="E","3"="S","4"="W"))
df$egoPre <- plyr::revalue(df$egoPre, c("1"="F","2"="R","3"="B","4"="L"))
df$egoCurr <- plyr::revalue(df$egoCurr, c("1"="F","2"="R","3"="B","4"="L"))
df$egoNext <- plyr::revalue(df$egoNext, c("1"="F","2"="R","3"="B","4"="L"))



for (cellname in unique(df$cell)){
  print(cellname)
  pdf(paste(savepath,ts,"_",code,"_",str_replace(cellname,"\\\\","-"),".pdf",sep=""),onefile=TRUE)

  celldf <- subset(df, cell == cellname)
  
  u<-try({
  # 1- prev turn, all variables model
  mod<-lm(rate ~ alloPrevPre*field + alloPrevPost*field + egoPre*field, data=celldf)
  celldf$fit <- predict(mod, newdata=celldf, re.form=NA)
  celldf$err <- NA
  z <- qnorm(0.99)
  
  for(a in c("N","E","S","W")){
    for(b in c("N","E","S","W")){
      for(c in c("F","R","B","L")){
        for(f in seq(0,length(unique(celldf$field))-1)){
          group <- subset(celldf,celldf$alloPrevPre == a & celldf$alloPrevPost == b & celldf$egoPre==c & celldf$field==f)
          if(nrow(group)>1){
            celldf[celldf$alloPrevPre == a & celldf$alloPrevPost == b & celldf$egoPre==c & celldf$field==f,14] <- sd(group$rate)/sqrt(nrow(group))
          }
          
        }
      }
    }
  }
  
  p1<-ggplot(data=celldf, aes(x=alloPrevPre,y=fit,group=egoPre,fill=alloPrevPost))+
    facet_wrap(~field)+
    geom_bar(position='dodge',stat='identity')+
    geom_errorbar(aes(ymin=fit-err,ymax=fit+err),position='dodge')+
    geom_text(position = position_dodge(width=1), aes(group=egoPre, label=egoPre, hjust=0), angle=90)+
    ggtitle(sprintf("%s Prev Turn, All Vars, Interacting Field", unique(celldf$cell)))+
    scale_color_brewer(palette="Dark2")
  print(p1)
  
  },silent=FALSE)
  
  u<-try({
  # 2- curr turn, all variables model
  
  mod<-lm(rate ~ alloCurrPre*field + alloCurrPost*field + egoCurr*field, data=celldf)
  celldf$fit <- predict(mod, newdata=celldf, re.form=NA)
  celldf$err <- NA
  z <- qnorm(0.99)
  
  for(a in c("N","E","S","W")){
    for(b in c("N","E","S","W")){
      for(c in c("F","R","B","L")){
        for(f in seq(0,length(unique(celldf$field))-1)){
          group <- subset(celldf,celldf$alloCurrPre == a & celldf$alloCurrPost == b & celldf$egoCurr==c & celldf$field==f)
          if(nrow(group)>1){
            celldf[celldf$alloCurrPre == a & celldf$alloCurrPost == b & celldf$egoCurr==c & celldf$field==f,14] <- sd(group$rate)/sqrt(nrow(group))
          }
          
        }
      }
    }
  }
  
  p1<-ggplot(data=celldf, aes(x=alloCurrPre,y=fit,group=egoCurr,fill=alloCurrPost))+
    facet_wrap(~field)+
    geom_bar(position='dodge',stat='identity')+
    geom_errorbar(aes(ymin=fit-err,ymax=fit+err),position='dodge')+
    geom_text(position = position_dodge(width=1), aes(group=egoCurr, label=egoCurr, hjust=0), angle=90)+
    ggtitle(sprintf("%s Closest Turn, All Vars, Interacting Field", unique(celldf$cell)))+
    scale_color_brewer(palette="Dark2")
  print(p1)
  
  },silent=FALSE)
  
  u<-try({
  
  # 3 - next turn all vars
  
  mod<-lm(rate ~ alloNextPre*field + alloNextPost*field + egoNext*field, data=celldf)
  celldf$fit <- predict(mod, newdata=celldf, re.form=NA)
  celldf$err <- NA
  z <- qnorm(0.99)
  
  for(a in c("N","E","S","W")){
    for(b in c("N","E","S","W")){
      for(c in c("F","R","B","L")){
        for(f in seq(0,length(unique(celldf$field))-1)){
          group <- subset(celldf,celldf$alloNextPre == a & celldf$alloNextPost == b & celldf$egoNext==c & celldf$field==f)
          if(nrow(group)>1){
            celldf[celldf$alloNextPre == a & celldf$alloNextPost == b & celldf$egoNext==c & celldf$field==f,14] <- sd(group$rate)/sqrt(nrow(group))
          }
          
        }
      }
    }
  }
  
  p1<-ggplot(data=celldf, aes(x=alloNextPre,y=fit,group=egoNext,fill=alloNextPost))+
    facet_wrap(~field)+
    geom_bar(position='dodge',stat='identity')+
    geom_errorbar(aes(ymin=fit-err,ymax=fit+err),position='dodge')+
    geom_text(position = position_dodge(width=1), aes(group=egoNext,label=egoNext, hjust=0), angle=90)+
    ggtitle(sprintf("%s Next Turn, All Vars, Interacting Field", unique(celldf$cell)))+
    scale_color_brewer(palette="Dark2")
  print(p1)
  
  },silent=FALSE)

  dev.off()
  
  
}


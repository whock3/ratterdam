# Repetition Direction and Turn analysis
# WH 7/8/21 - 

library(lme4)
library(ggplot2)
library(lmtest)
library(stringr)
library(emmeans)
library(tidyr)
library(ggpubr)

# define data paths and load in data
path <- "E:\\Ratterdam\\R_data_repetition\\20210708-151540_R859D2_dirLMER_1.5vfilt.csv"
savepath <- "E:\\Ratterdam\\R_data_repetition\\210624-onward_turnLM\\"
code <- "R859D2"
df <- read.csv(path,header=TRUE)
df <- df[!(df$dirC==0),]
df <- df[!(df$dirP1==0),]
df <- df[!(df$dirM1==0),]

ts <- str_replace(Sys.time()," ","_")
ts <- str_replace_all(ts, ":", "_")

df$field <- as.factor(df$field)
df$dirC <- as.factor(df$dirC)
df$dirM1 <- as.factor(df$dirM1)
df$dirP1 <- as.factor(df$dirP1)
df$epoch <- as.factor(df$epoch)


df$dirC <- plyr::revalue(df$dirC, c("1"="N","2"="E","3"="S","4"="W"))
df$dirM1 <- plyr::revalue(df$dirM1, c("1"="N","2"="E","3"="S","4"="W"))
df$dirP1 <- plyr::revalue(df$dirP1, c("1"="N","2"="E","3"="S","4"="W"))


 

for(cellname in unique(df$unit)){
  print(cellname)
  celldf <- subset(df, unit == cellname)
  adjP <- 0.05/length(unique(celldf$field)) # bonferroni correction by field #
  
  pdf(paste(savepath,ts,"_",code,"_",str_replace(cellname,"\\\\","-"),".pdf",sep=""),onefile=TRUE)
  
  for(f in unique(celldf$field)){
    print(f)
    u<-try({
    fdf = subset(celldf, field==f)
    
    # ANOVA - three way interaction plus temporal epoch
    amod <- aov(rate ~ dirM1*dirC*dirP1+epoch, data = fdf)
    ramod <- drop_na(summary(amod)[[1]]) # need to access inside the s-class to get the table
    
    # GLM - three way interaction plus temporal epoch
    gmod <- glm(rate+1 ~ dirM1*dirC*dirP1+epoch, family="Gamma", data = fdf)
    rgmod <- drop_na(anova(gmod, test="F")) # run anova to get main effects and interactions (i.e. not all level contrasts)
    
    asig <- row.names(ramod[ramod$`Pr(>F)`<adjP,])
    gsig <- row.names(rgmod[rgmod$`Pr(>F)`<adjP,])
    asig <- str_c(asig,collapse=',')
    gsig <- str_c(gsig,collapse=',')
    
    asig <- c(asig,".") # need something else geom_label complains if empty
    gsig <- c(gsig,".") # ibid
    
    p1 <- ggplot(fdf, aes(x=dirM1, y=rate))+
      geom_violin(color='black',fill="slategray1")+
      geom_boxplot(width=0.1,color='navy')+
      theme(text=element_text(size=10))+
      ggtitle(sprintf("Field %s, Prev Dir, All ANOVA effects: %s", f, asig))
      
    p2 <- ggplot(fdf, aes(x=dirC, y=rate))+
      geom_violin(color='black',fill="slategray1")+
      geom_boxplot(width=0.1,color='navy')+
      theme(text=element_text(size=10))+
      ggtitle(sprintf("Field %s, Curr Dir, All GLM effects: %s", f, gsig))
  
    p3 <- ggplot(fdf, aes(x=dirP1, y=rate))+
      geom_violin(color='black',fill="slategray1")+
      geom_boxplot(width=0.1,color='navy')+
      theme(text=element_text(size=10))+
      ggtitle(sprintf("Field %s, Next Dir", f))
    
    print(ggarrange(p1,p2,p3,nrow=3,ncol=1))
    },silent=FALSE)
  
  }
  dev.off()
}


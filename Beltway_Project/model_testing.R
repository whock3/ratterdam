# imports
library(lme4)
library(splines)
library(ggplot2)
library(lmtest)
library(stringr)


setwd("E:\\UserData\\Documents\\GitHub\\ratterdam\\")
source("E:\\UserData\\Documents\\GitHub\\ratterdam\\Beltway_Project\\cicheck.R")
source("E:\\UserData\\Documents\\GitHub\\ratterdam\\Beltway_Project\\glmer_fx.R")

# Read in data
code <- "R859BRD5"
save <- FALSE # toggle to save multipage pdfs and wald csvs
database <- "E:\\Ratterdam\\R_data_beltway\\"
datapath <- sprintf("%s%s",database,"20210520-181551_R859BRD5_1.5vfilt_0.5stepsmooth_24bins_2R_3qual.csv")
df <- read.csv(datapath,header=TRUE)

# Select output, create timestamp 
figbasepath <- "E:\\Ratterdam\\R_data_beltway\\graphs\\210317_reassignmentAnalysis\\"
ts <- str_replace(Sys.time()," ","_")
ts <- str_replace_all(ts, ":", "_")


# if(code=="R859BRD5"){
# 
#   df<-df[!df$alley==3,]
#   df<-df[!df$alley==5,]
#   df<-df[!df$alley==7,]
# }


# set up dfs that will store data. first is wald info, second is ci overlap check outcome
# for all cells/alleys (for wald and fitted curves)
wdf <- data.frame(matrix(ncol=5))
colnames(wdf) <- c("cellname", "alleyID", "low", "up", "fe")

cicheckdf <- data.frame(matrix(ncol=4))
colnames(cicheckdf) <- c("cell", "alley", "wald", "fit")

nsplineknots <- 6 # set number spline knots 

for(cellID in unique(df$cell)){

  
  # Select and set up data
  celldf <- subset(df, cell == cellID)
  celldf$rate <- log(celldf$rate+1)
  cellname <- unique(df[df$cell==cellID,]$name)
  alleys <- unique(celldf$alley)
  nalleys <-length(alleys)
  
  bonf <- 0.17/(3*nsplineknots*nalleys*5) #*8 for reassignment analysis
  cipct <- 1-bonf
  z <- qnorm(cipct)
  
  if(save==TRUE){
  pdf(paste(figbasepath,ts,"_",code,"_",str_replace(cellname,"\\\\","-"),".pdf",sep=""),onefile=TRUE)
  }
  
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
    ciall <- rbind(cidf, cireorddf[c("textureC", "ns(spatialBin, 6)1:textureC","ns(spatialBin, 6)2:textureC","ns(spatialBin, 6)3:textureC","ns(spatialBin, 6)4:textureC","ns(spatialBin, 6)5:textureC","ns(spatialBin, 6)6:textureC"),])
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
    
    if(save==TRUE){
    print(waldplot)
    }
    
    u<- try({
    celldf <- lmer_routine(celldf, nalleys)
    celldf_no_r <- data.frame(celldf[celldf$reward=="0",])
    
    outcome <- checkAlleyNonoverlap(celldf_no_r,alleyID)
    if(outcome==TRUE){
      print("Pass")
      
    }
    },silent=FALSE)
    
    
    # Check non-overlap between walds and 0, and fits and each other
    # fitoverlap <- checkAlleyNonoverlap(celldf, alleyID)
    # 
    # waldoverlap <- checkWaldNonoverlap(ndf)
    # r <- data.frame("cell"=cellname, "alley"=alleyID,"wald"=waldoverlap,"fit"=fitoverlap)
    # cicheckdf <- rbind(cicheckdf, r)
    # 
    
    # ggplot of splines with CI of fits
    celldf <- celldf[celldf$reward=="0",]
    p <- ggplot(data=celldf, aes(x=spatialBin))+
      geom_line(aes(y=fit, color=texture,linetype=factor(reward)))+
      geom_ribbon(aes(y=fit, ymin=fit-fitCI, ymax=fit+fitCI, fill=texture), alpha=0.2)+
      ggtitle(sprintf("Cell %s, Alley %s", cellname, alleyID))+
      scale_color_manual(values=c("red", "blue", "darkgreen"))+
      scale_fill_manual(values=c("red", "blue", "darkgreen"))
    
    if(save==TRUE){
    print(p)
    }
  }
  
  if(save==TRUE){
  dev.off()
  }
  
    
}


  
 
  


  
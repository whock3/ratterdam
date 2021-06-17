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
path <- "E:\\Ratterdam\\R_data_repetition\\20210301-205647_R859D2_1vfilt_.csv"
savepath <- "E:\\Ratterdam\\R_data_repetition\\210325_timeDirLMER\21-06-11\\"
df <- read.csv(path,header=TRUE)

# get rid of epoch 0, it's a py artefact from bisecting approach to finding 
# epoch that needs to be fixed
df <- df[!df$epoch==0,]
df$epoch <- as.factor(df$epoch)


exp <- "R859D2"

for (cellname in unique(df$cell)){
  print(cellname)
    celldf <- subset(df, cell == cellname)
    
    #define params
    nsplineknots <- 9
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
    
    mod <- NULL
    
    u <- try(mod <- lmer(rate ~ histdir*field + (1+histdir|epoch) + (1+field|epoch), data=celldf))
    if (!is.null(mod)){
      celldf$fit <- predict(mod, newdata=celldf, re.form=NA)
      
      #Designmat <- model.matrix(rate ~ ns(direction,nsplineknots):epoch + (ns(direction,nsplineknots)|field) + (epoch|field), celldf)
      # Designmat <- model.matrix(rate ~ histdir*field, celldf)
      # predvar <- diag(Designmat %*% vcov(mod) %*% t(Designmat))
      # celldf$fitCI <- sqrt(predvar)*z
      
      # interaction plot w spline fits
      # p <- ggplot(celldf, aes(x=histdir, y=rate, color=as.factor(epoch)))+
      #   geom_line(aes(histdir, fit, group=epoch))+
      #   geom_ribbon(aes(y=fit, ymin=fit-fitCI, ymax=fit+fitCI, fill=as.factor(epoch)), alpha=0.2)+
      #   ggtitle(sprintf("%s %s", exp, cellname))+
      #   geom_point(size=2,alpha=0.5)
      # 
      # print(p)
      # 
      # 
      # #individual plots by group 
      # (mm_plot <- ggplot(celldf, aes(x = histdir, y = rate, color = as.factor(epoch))) +
      #     facet_wrap(~field, nrow=2) +   
      #     geom_point(alpha = 0.5) +
      #     geom_bar(aes(histdir, fit, group=epoch))+ 
      #     #geom_errorbar(aes(y=fit, ymin=fit-fitCI, ymax=fit+fitCI, fill=as.factor(epoch)), alpha=0.2)+
      #     ggtitle(sprintf("%s %s", exp, cellname))# adding predicted CIs from mixed model
      # )
      # print(mm_plot)
      
     
      title <- sprintf("%s_%s", exp, cellname)
      title <- gsub("\\\\","-",title)
      
      #png(file=sprintf("%s%s_overall.png",savepath,title))
      # p <- ggplot(celldf, aes(x=histdir, y=fit, fill=epoch))+
      #       geom_bar(position='dodge',stat='identity')+
      #       geom_errorbar(aes(ymin=fit-(sd(rate)/sqrt(length(rate))), ymax=fit+(sd(rate)/sqrt(length(rate)))),position='dodge')+
      #       ggtitle(title)
      # print(p)
      # #dev.off()
      
      
      #png(file=sprintf("%s%s_fields.png",savepath,title))
      fp <- ggplot(celldf, aes(x=field, y=fit, fill=as.factor(histdir)))+
        facet_wrap(~epoch, nrow=3)+
        geom_bar(position='dodge',stat='identity')+
        geom_errorbar(aes(ymin=fit-(sd(rate)/sqrt(length(rate))), ymax=fit+(sd(rate)/sqrt(length(rate)))),position='dodge')+
        ggtitle(title)
      
      print(fp)
      #dev.off()
    
    }
    


}

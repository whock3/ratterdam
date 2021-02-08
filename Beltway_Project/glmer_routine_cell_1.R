
cellID <- 12
celldf <- subset(df, cell == cellID)
alleys <- unique(celldf$alley)
cellname <- unique(subset(df, cell==cellID, select=name))

for(alleyID in alleys){
  
  print(sprintf("============ Cell: %s, Cell ID: %s, alley ID: %s =============== ", cellname, cellID, alleyID))
  mod <- glmer_alley(cellID, alleyID, df)
  # title <- sprintf("Cell %s alley %s", cellname, alleyID)
  # print(sjPlot::plot_model(mod, type='int', title=title))
  print(summary(mod))
  ci = confint(mod, parm='beta_', level = 0.97, method='Wald')
  print("Confidence Intervals (Wald)")
  print(ci)
  # 
  # assemble into df with wald ci w intercept removed
  fe <- fixef(mod)
  newdf <- data.frame("low"=ci[,1], "up"=ci[,2],"fe"=fe)
  newdf <- newdf[-c(1),]
  
  x = seq(1,nrow(newdf))
  waldplot <- ggplot(data=newdf, aes(x=x,y=fe))+
              geom_point()+
              geom_errorbar(aes(ymin=low,ymax=up))+
              scale_x_discrete(limits=row.names(newdf))+
              theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
              geom_hline(yintercept=0)+
              ggtitle(sprintf("Cell %s Alley %s Wald CI 1.5/98.5",cellname, alleyID))
  print(waldplot)
}
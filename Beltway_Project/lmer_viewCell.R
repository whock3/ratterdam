# script for looking at a single cell from, e.g. a lmer loop




cellname <- "TT13\\cl-maze1.1"
alley <- 7

celldf <- df[df$name==cellname,]
celldf_a <- celldf[celldf$alley==alley,]
nalleys <- length(unique(celldf$alley))
celldf_a <- lmer_routine(celldf_a, nalleys) # nalleys is number of alleys w activity and used for bonferroni corr


celldf_a <- celldf_a[celldf_a$reward=="0",]
p <- ggplot(data=celldf_a, aes(x=spatialBin))+
       geom_line(aes(y=fit, color=texture,linetype=factor(reward)))+
       geom_ribbon(aes(y=fit, ymin=fit-fitCI, ymax=fit+fitCI, fill=texture), alpha=0.2)+
       ggtitle(sprintf("Cell %s, Alley %s", cellname, alley))+
       scale_color_manual(values=c("red", "blue", "darkgreen"))+
       scale_fill_manual(values=c("red", "blue", "darkgreen"))
print(p)
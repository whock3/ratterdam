# script for looking at a single cell from, e.g. a lmer loop
#




cellname <- "TT10\\cl-maze1.6"
alley <- 16

celldf <- df[df$name==cellname,]
celldf_a <- celldf[celldf$alley==alley,]
nalleys <- length(unique(celldf$alley))
celldf_a <- lmer_routine(celldf_a, nalleys) # nalleys is number of alleys w activity and used for bonferroni corr


# plot model results with CI broken down by reward
p <- ggplot(data=celldf_a, aes(x=spatialBin))+
       geom_line(aes(y=fit, color=texture))+
       geom_ribbon(data=celldf_a[celldf_a$reward=="1",], aes(y=fit, x=spatialBin, ymin=fit-fitCI, ymax=fit+fitCI, fill=texture), alpha=0.2, inherit.aes=FALSE)+
       geom_ribbon(data=celldf_a[celldf_a$reward=="0",], aes(y=fit, x=spatialBin, ymin=fit-fitCI, ymax=fit+fitCI, fill=texture), alpha=0.2, inherit.aes=FALSE)+
       ggtitle(sprintf("Cell %s, Alley %s", cellname, alley))+
       scale_color_manual(values=c("red", "blue", "darkgreen"))+
       scale_fill_manual(values=c("red", "blue", "darkgreen"))
print(p)


# plot raw data
p <- ggplot(data=celldf_a, aes(x=spatialBin))+
        geom_smooth(aes(y=rate, color=texture,linetype=factor(reward)))+
        ggtitle(sprintf("Cell %s, Alley %s", cellname, alley))+
        scale_color_manual(values=c("red", "blue", "darkgreen"))+
        scale_fill_manual(values=c("red", "blue", "darkgreen"))
print(p)

#plot fits trying to see ns(trial,2) effect on data 

p <- ggplot(data=celldf_a, aes(x=spatialBin))+
        geom_line(aes(y=fit, color=texture,linetype=factor(reward)))+
        ggtitle(sprintf("Cell %s, Alley %s", cellname, alley))+
        scale_color_manual(values=c("red", "blue", "darkgreen"))+
        scale_fill_manual(values=c("red", "blue", "darkgreen"))
print(p)





p <- qplot(nullcells$`X[[i]]`)+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(p)


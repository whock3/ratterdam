# Script to run 1x shuffle on a given cell and visualize the the spline fit
# of texture response across alleys with CIs. Used to see how cells are passing
# null test so often.

cellname <- "TT10\\cl-maze1.6"
alley <- 16

nshuffles <- 1000

celldf <- df[df$name==cellname,]
nalleys <- length(unique(celldf$alley)) # for alpha adjustment, even though we're only looking at 1 alley here
celldf_a <- celldf[celldf$alley==alley,]
celldf_a_copy <- data.frame(celldf_a) # used to shuffle, keep original trial structure in celldf_a 


celldf_a <- lmer_routine(celldf_a, nalleys) # nalleys is number of alleys w activity and used for bonferroni corr

# create trial list
ntrials <- length(unique(celldf_a$trial))
trialList <- data.frame(texture=factor())
for (i in 0:ntrials-1){
  stim <- unique(celldf_a[celldf_a$trial==i,]$texture)
  trialList <- rbind(trialList, paste(stim))
}

# Shuffle trials
shuffledTrials <- data.frame(trialList[sample(1:ntrials),])
for (j in 1:ntrials){
  celldf_a_copy[celldf_a$trial==j-1,]$texture <- shuffledTrials[j,]
}

# Run shuffle through lmer and check CI overlap (non R trials)
celldf_a_copy <- lmer_routine(celldf_a_copy, nalleys)
celldf_a_copy_no_r <- data.frame(celldf_a_copy[celldf_a_copy$reward=="0",])
outcome <- checkAlleyNonoverlap(celldf_a_copy_no_r,alley)
print(outcome)
#Visualize real data and single shuffle result
p1 <- ggplot(data=celldf_a[celldf_a$reward=="0",], aes(x=spatialBin))+
  geom_line(aes(y=fit, color=texture,linetype=factor(reward)))+
  geom_ribbon(aes(y=fit, ymin=fit-fitCI, ymax=fit+fitCI, fill=texture), alpha=0.2)+
  ggtitle(sprintf("Cell %s, Alley %s, Real Data", cellname, alley))+
  scale_color_manual(values=c("red", "blue", "darkgreen"))+
  scale_fill_manual(values=c("red", "blue", "darkgreen"))
print(p1)

p2 <- ggplot(data=celldf_a_copy_no_r, aes(x=spatialBin))+
  geom_line(aes(y=fit, color=texture,linetype=factor(reward)))+
  geom_ribbon(aes(y=fit, ymin=fit-fitCI, ymax=fit+fitCI, fill=texture), alpha=0.2)+
  ggtitle(sprintf("Cell %s, Alley %s, single Shuffle", cellname, alley))+
  scale_color_manual(values=c("red", "blue", "darkgreen"))+
  scale_fill_manual(values=c("red", "blue", "darkgreen"))
print(p2)

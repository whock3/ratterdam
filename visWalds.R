
passcounts <- vector()
nonpasscounts <- vector()
cellcounts <- 0
names <- unique(wdf$cellname)
for(targ in c(7,8,9,10,11,12,14,15,16)){
  passvals <- vector()
  nonpassvals <- vector()
for(name in names){
  cellpass <- FALSE
  cdf <- wdf[wdf$cellname==name,]
  alleys <- unique(cdf$alleyID)
    for(a in alleys){
      cadf <- cdf[cdf$alleyID==a,]
      if(cadf[targ,3] > 0 | cadf[targ,4] < 0){
        passvals <- append(passvals, cadf[targ,5])
        cellpass <- TRUE
      }
      else{
        nonpassvals <- append(nonpassvals, cadf[targ,5])
      }
    }
  
}

passcounts <- append(passcounts, length(passvals))
nonpasscounts <- append(nonpasscounts, length(nonpassvals))

nv <- data.frame('value'=abs(nonpassvals))
pv <- data.frame('value'=abs(passvals))
p <- ggplot(nv, aes(x=value))+
      geom_histogram(aes(y=1*(..density..)/sum(..density..)),color='grey', binwidth=0.1, alpha=0.8)+
      geom_histogram(data=pv, aes(y=1*(..density..)/sum(density)), color='firebrick',fill='firebrick', binwidth=0.1, alpha=0.8)+
      theme(axis.text.x = element_text(size=14), 
            axis.text.y = element_text(size=14),
            axis.title.y = element_text(size=16),
            axis.title.x = element_text(size=16))+
      xlab("Absolute Wald Coefficient")+ylab("Frequency Density")+
      geom_vline(aes(xintercept=0), color="black", linetype="dashed", size=1.5)+
      #geom_density(data=nv,alpha=0.4, fill='grey')+
    #geom_density(data=pv,alpha=0.4, fill='firebrick')+
    ggtitle(targ)
    

print(p)

}



p <- ggplot(cdf, aes(fill=type, y=counts, x=cond))+
      geom_bar(position='stack', stat='identity')+
      scale_fill_manual(values=c("grey","firebrick"))+
    scale_x_continuous(breaks = c(seq(1,9)), labels= labels, vjust=0.1)+
      theme(axis.text.x = element_text(size=20, angle=45), 
            axis.text.y = element_text(size=20))
      
print(p)

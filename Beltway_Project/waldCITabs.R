
passcounts <- vector()
nonpasscounts <- vector()
cellcounts <- 0
names <- unique(wdf$cellname)
passvals <- vector()
nonpassvals <- vector()

for(name in names){
  cellpass <- FALSE
  for(targ in c(7,8,9,10,11,12,14,15,16)){
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
  
  if(cellpass==TRUE){
    cellcounts <- cellcounts + 1
  }
}
print(datapath)
print(cellcounts)
print(length(names))
print(cellcounts/length(names))
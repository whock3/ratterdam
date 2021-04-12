# For alley

singleCIbounds <- function(celldf, alley, b, stim){
  fit <- unique(celldf[(celldf$alley==alley & celldf$spatialBin==b & celldf$texture==stim),]$fit)
  ci <- unique(celldf[(celldf$alley==alley & celldf$spatialBin==b & celldf$texture==stim),]$fitCI)
  upper <- fit+ci
  lower <- fit-ci
  return(c(upper,lower))
  
}

checkSingleNonoverlap <- function(ciX, ciY){
  # ciX,Y are lists where its [upper, lower] ci bound
  
  #length will be 0 if there is no sample, i.e. no sampling at that bin for that txt
  if (length(ciX) > 0 & length(ciY) > 0){
    if(ciX[1] < ciY[2]){
      nonoverlap = TRUE
      
    }
    else if(ciX[2] > ciY[1]){
      nonoverlap = TRUE
    }
    else{
      nonoverlap = FALSE
    }
  }
  else{
    nonoverlap = FALSE
  }
  
  return(nonoverlap)
}

checkAlleyNonoverlap <- function(celldf, alley){
  bins <- seq(0,max(celldf$spatialBin))
  
  alleynonoverlap <- FALSE
  
  for(b in bins){
    cia <- singleCIbounds(celldf, alley, b, 'A')
    cib <- singleCIbounds(celldf, alley, b, 'B')
    cic <- singleCIbounds(celldf, alley, b, 'C')
    abo <- checkSingleNonoverlap(cia, cib)
    bco <- checkSingleNonoverlap(cib, cic)
    cao <- checkSingleNonoverlap(cic,cia)
    if(abo == TRUE || bco == TRUE || cao == TRUE){
      alleynonoverlap = TRUE
    }
    
  }
  return(alleynonoverlap)
  
}

checkWaldNonoverlap <- function(walddf){
  conditions = c("ns(spatialBin, 5)1:textureB",
                 "ns(spatialBin, 5)2:textureB",
                 "ns(spatialBin, 5)3:textureB",
                 "ns(spatialBin, 5)4:textureB",
                 "ns(spatialBin, 5)5:textureB",
                 "ns(spatialBin, 5)1:textureC",
                 "ns(spatialBin, 5)2:textureC",
                 "ns(spatialBin, 5)3:textureC",
                 "ns(spatialBin, 5)4:textureC",
                 "ns(spatialBin, 5)5:textureC",
                 "ns(spatialBin, 5)1:textureC1",
                 "ns(spatialBin, 5)2:textureC1",
                 "ns(spatialBin, 5)3:textureC1",
                 "ns(spatialBin, 5)4:textureC1",
                 "ns(spatialBin, 5)5:textureC1"
                 )
  nonoverlap <- FALSE
  for(cond in conditions){
    low <- walddf[cond,'low']
    up <- walddf[cond, 'up']
    if(low > 0 | up < 0){
      nonoverlap <- TRUE
    }
    
  }
  return(nonoverlap)
  
  
}


calcAlleyMaxSD <- function(celldf, alley) {

  celldf <- celldf[celldf$alley==alley,]
  nbins <- length(unique(celldf$spatialBin))
  
  overallMaxDist <- 0 
  
  for (b in 0:(nbins-1)){
  
  aci=unique(celldf[celldf$spatialBin==b & celldf$texture=='A',]$fitCI)
  afit=unique(celldf[celldf$spatialBin==b & celldf$texture=='A',]$fit)
  
  bci=unique(celldf[celldf$spatialBin==b & celldf$texture=='B',]$fitCI)
  bfit=unique(celldf[celldf$spatialBin==b & celldf$texture=='B',]$fit)
  
  cci=unique(celldf[celldf$spatialBin==b & celldf$texture=='C',]$fitCI)
  cfit=unique(celldf[celldf$spatialBin==b & celldf$texture=='C',]$fit)
  
  ab_dist = mean(c((abs(afit-bfit)/aci),(abs(afit-bfit)/bci)))
  bc_dist = mean(c((abs(bfit-cfit)/bci),(abs(bfit-cfit)/cci)))
  ca_dist = mean(c((abs(cfit-afit)/cci),(abs(cfit-afit)/aci)))
  
  dists <- c(ab_dist, bc_dist, ca_dist)
  dists <- dists[!is.na(dists)]
  
  
  if (length(dists > 0)){
  
    max_dist <- max(dists)
    
    if (max_dist > overallMaxDist){
      
      overallMaxDist <- max_dist
      }
  
    }
  
  }
  
  
  return(overallMaxDist)
}
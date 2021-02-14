# For alley

singleCIbounds <- function(celldf, alley, b, stim){
  fit <- unique(celldf[(celldf$alley==alley & celldf$spatialBin==b & celldf$texture==stim),]$fit)
  ci <- unique(celldf[(celldf$alley==alley & celldf$spatialBin==b & celldf$texture==stim),]$fitCI)
  upper <- fit+ci
  lower <- fit-ci
  return(c(upper,lower))
  
}

checkSingleNonoverlap <- function(ciA, ciB){
  # ciA,B are lists where its [upper, lower] ci bound
  if(ciA[1] < ciB[2]){
    nonoverlap = TRUE
    
  }
  else if(ciA[2] > ciB[1]){
    nonoverlap = TRUE
  }
  else{
    nonoverlap = FALSE
  }
  
  return(nonoverlap)
}

checkAlleyNonoverlap <- function(celldf, alley){
  #celldf may or may not be alreayd filtered by alley
  # (possibly redundant) filtering done in fx
  celldf <- celldf[celldf$alley==alley,]
  bins <- seq(0,max(celldf$spatialBin))
  
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
    else{
      alleynonoverlap = FALSE
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
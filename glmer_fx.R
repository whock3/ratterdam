glmer_cell <- function(cellID, df){
  
  celldf <- subset(df, cell == cellID)
  nalleys <- length(unique(celldf$alley))
  
  
  # If only 1 alley with sufficient activity for analysis, do not include alley as RE
  if(nalleys > 20){
    mod <- glmer(log(rate+1)+1 ~ spatialBin*texture
                 + (1|alley)
                 + (1|trial)
                 + (1|reward), 
                 family=Gamma, 
                 data = celldf, 
                 control = glmerControl(optimizer='Nelder_Mead'))
  }
  else{
    mod <- glmer(log(rate+1)+1 ~ spatialBin*texture
                 + (1|trial)
                 + (1|reward), 
                 family=Gamma, 
                 data = celldf, 
                 control = glmerControl(optimizer='Nelder_Mead'))
  }
  
  return(mod)
  
}

glmer_alley <- function(cellID, alleyID, df){
  celldf <- subset(df, cell == cellID | alley == alleyID)
  mod <- glmer(log(rate+1)+1 ~ ns(spatialBin,3)*texture
               + (1|trial)
               + (1|reward), 
               family=Gamma, 
               data = celldf, 
               control = glmerControl(optimizer='Nelder_Mead'),
               nAGQ=0
               )
  
  return(mod)
}
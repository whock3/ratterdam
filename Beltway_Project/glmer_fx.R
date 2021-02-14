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
  mod <- glmer(rate+1 ~ ns(spatialBin,3)*texture+reward
               + (1|trial),
               #+ (1|reward), 
               family=Gamma, 
               data = celldf, 
               control = glmerControl(optimizer='Nelder_Mead'),
               nAGQ=0
               )
  
  return(mod)
}

glmer_alley_int <- function(df){
  mod <- glmer(rate+1 ~ ns(spatialBin,3)*texture+reward
               + (1|trial),
               family=Gamma, 
               data = df, 
               control = glmerControl(optimizer='Nelder_Mead'),
               nAGQ=0
  )
  
  return(mod)
}

lmer_alley_int <- function(mydf){
  mod <- lmer(rate ~ ns(spatialBin,6)*texture+reward
               + (1|trial),
               data = mydf
  )
  
  return(mod)
}

lmer_alley_main <- function(mydf){
  mod <- lmer(rate ~ ns(spatialBin,3)+texture+reward
              + (1|trial),
              #+ (1|reward),
              data = mydf
  )
  
  return(mod)
}

lmer_alley_none <- function(mydf){
  mod <- lmer(rate ~ ns(spatialBin,3)+reward
              + (1|trial),
              #+ (1|reward),
              data = mydf
  )
  
  return(mod)
}

predFun <- function(fittedmodel, df){
  fit <- exp(predict(fittedmodel, newdata=df, re.form = NA))
  return(fit)
  
}
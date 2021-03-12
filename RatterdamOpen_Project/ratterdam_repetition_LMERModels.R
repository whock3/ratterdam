# G/LMER Models for Repetition

repLMER_timeDir <- function(mydf){
  mod <- lmer(rate ~ ns(spatialBin,6)*texture
              + (1|trial),
              data = mydf
  )
  
  return(mod)
}

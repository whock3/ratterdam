# Case checks for repeating fields
if(unique(field$Repeating)=='True'){
  
  #check full model against partial C+N and partial P+C, see if its better than both
  if((lrFullA[2,"Pr(>Chisq)"]<(0.05/bfAdj))&(lrFullB[2,"Pr(>Chisq)"]<(0.05/bfAdj))){
    repFull <- repFull + 1
  }
  #otherwise drop down to test C vs C+P
  else if(lrPrev[2,"Pr(>Chisq)"]<(0.05/bfAdj)){
    repPD <- repPD + 1
  }
  #otherwise drop down to test C vs C+N
  else if(lrNext[2,"Pr(>Chisq)"]<(0.05/bfAdj)){
    repND <- repND + 1
  }
  #otherwise drop down to test null vs C
  else if(lrCurr[2,"Pr(>Chisq)"]<(0.05/bfAdj)){
    repCD <- repCD + 1
  }
}


# Case checks for nonrepeating fields
if(unique(field$Repeating)=='False'){
  
  #check full model against partial C+N and partial P+C, see if its better than both
  if((lrFullA[2,"Pr(>Chisq)"]<(0.05/bfAdj))&(lrFullB[2,"Pr(>Chisq)"]<(0.05/bfAdj))){
    nonrepFull <- nonrepFull + 1
  } 
  #otherwise drop down to test C vs C+P
  else if(lrPrev[2,"Pr(>Chisq)"]<(0.05/bfAdj)){
    nonrepPD <- nonrepPD + 1
  }
  #otherwise drop down to test C vs C+N
  else if(lrNext[2,"Pr(>Chisq)"]<(0.05/bfAdj)){
    nonrepND <- nonrepND + 1
  }
  #otherwise drop down to test null vs C
  else if(lrCurr[2,"Pr(>Chisq)"]<(0.05/bfAdj)){
    nonrepCD <- nonrepCD + 1
  }
}
 

#### Generate OU process under H0 ##########
simDiscretizedOrnsteinH0<-function(dmark, nbMarkers, sizeChro){ 
  
  # generate a discretized OU process on marker locations
  X = rep(0,nbMarkers)
  # X contains my test statistic computed along the genome
  X[1]=rnorm(1,0,1);
  
 for (i in 1:nbMarkers){
  epsilon= sqrt(1-exp(-4*dmark))*rnorm(1,0,1);
  X[i+1]=exp(-2*dmark)*X[i]+ epsilon;
 }
  
 return(X)  
  
}


#### Generate OU process under H1 ##########
simDiscretizedOrnsteinH1<-function(nbTRN, dmark, nbMarkers, sizeChro, posMarkers, posQTL, effectQTL, Psi0, errorVariance){ 
  
  # generate a discretized OU process on marker locations
  X = rep(0,nbMarkers)
  Z = rep(0,nbMarkers)
  # X contains my test statistic computed along the genome
  X[1]=rnorm(1,0,1);
  Z[1]=X[1] + sum((sqrt(nbTRN)*effectQTL/sqrt(errorVariance)) * exp(-2*abs(posMarkers[1]-posQTL ))) 
  
  for (i in 1:(nbMarkers-1)){
    epsilon= sqrt(1-exp(-4*dmark))*rnorm(1,0,1);
    X[i+1]=exp(-2*dmark)*X[i]+ epsilon;
    Z[i+1]=X[i+1] + sum((sqrt(nbTRN)*effectQTL/sqrt(errorVariance)) * sqrt(Psi0) * exp(-2*abs(posMarkers[i+1]-posQTL ))) 
  }
  # Z is my asymptotic test statistic at markers
   
  return(Z)
 
}
 

#######################  Smart Formula #################### 

computeAsymptProcessSmartFormula<-function(testStatOnMarkers, dmark, nbMarkers){
   
  currentMax=testStatOnMarkers[1]^2
  
  for (i in 1:(nbMarkers-1)){ 
    # deal with the ith marker interval
    
    currentMax=max(currentMax,testStatOnMarkers[(i+1)]^2)
    
    ratio=testStatOnMarkers[(i+1)]/testStatOnMarkers[(i)]
    
    if ( (ratio > exp(-2*dmark)) && (ratio < exp(2*dmark)) ) { 
      
      numSuppStat = testStatOnMarkers[(i+1)]^2 + testStatOnMarkers[i]^2 - 2*exp(-2*dmark)*testStatOnMarkers[(i+1)]*testStatOnMarkers[i]
      denomSuppStat = 1- exp(-4*dmark)
      suppStat = numSuppStat/denomSuppStat
      currentMax= max(currentMax,suppStat)
      
    }#end if
    
  }# end for
  
  return(currentMax)
  
   
}
  
  
 
 

#### Generate process under H0 according to Theorem 3.9 ##########
simDiscretizedProcessH0<-function(dmark, nbMarkers, sizeChro, Psi1, Psi0){ 
  
  # generate a discretized OU process on marker locations
  X = rep(0,nbMarkers)
  # X contains my test statistic computed along the genome
  X[1]=rnorm(1,0,sqrt(Psi1/Psi0));
  
 for (i in 1:nbMarkers){
  epsilon= sqrt((Psi1/Psi0)*(1-exp(-4*dmark)))*rnorm(1,0,1);
  X[i+1]=exp(-2*dmark)*X[i]+ epsilon;
 }
  
 return(X)  
  
}


#### Generate process under H1    ##########
simDiscretizedProcessH1<-function(nbTRN, dmark, nbMarkers, sizeChro, posMarkers, posQTL, effectQTL, Psi2, Psi1, Psi0, errorVariance){ 
  
  # generate a discretized OU process on marker locations
  X = rep(0,nbMarkers)
  Z = rep(0,nbMarkers)
  #Term =(Psi2-Psi1)/sqrt(Psi0)
  Term=(Psi2-Psi1)/sqrt(Psi1-0.57721^2)
  
  # X contains my test statistic computed along the genome
  X[1]=rnorm(1,0,sqrt(Psi1/Psi0));
  Z[1]=X[1] + Term * sum((sqrt(nbTRN)*effectQTL/sqrt(errorVariance)) * exp(-2*abs(posMarkers[1]-posQTL ))) 
  
  for (i in 1:(nbMarkers-1)){
     
    epsilon= sqrt((Psi1/Psi0)*(1-exp(-4*dmark)))*rnorm(1,0,1);
    X[i+1]=exp(-2*dmark)*X[i]+ epsilon;  
    Z[i+1]=X[i+1] + Term * sum((sqrt(nbTRN)*effectQTL/sqrt(errorVariance)) * exp(-2*abs(posMarkers[i+1]-posQTL ))) 
    
  }
  # Z is my asymptotic test statistic at markers
   
  return(Z)
 
}

 

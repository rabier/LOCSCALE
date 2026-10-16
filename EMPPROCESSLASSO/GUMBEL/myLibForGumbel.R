# My functions for Gumbel family


#### Generate Genomes ##########
generateGenomes<-function(nbInd, nbMarkersAndQTL,sizeChro){ 
  
  myGenome = matrix(0, nrow = nbInd, ncol = nbMarkersAndQTL)
  nbRecomb<- 0

  for (i in 1:nbInd){
  
    myRandBinomial=rbinom(1,1,0.5)
  
    if (myRandBinomial==0){
      myRandBinomial=-1	
    }
  
    myGenome[i,1]=myRandBinomial
    indice=1
    myGenome[i,2:nbMarkersAndQTL]=myGenome[i,1]
    posRecomb = rexp(1, 1)
    # generate exponential distribution with rate 1
  
    while (posRecomb<sizeChro){
    
      nbRecomb = nbRecomb+1
      # a recombination event happened
      # find the first marker on the left side of the recombination location
      indleftMarkerRecomb<-max(which(posRecomb-posMarkersAndQTL >0 ))
   
      if (myGenome[i,indleftMarkerRecomb]==1){
        # we switch the sign 
        myGenome[i,(indleftMarkerRecomb+1):nbMarkersAndQTL]=-1
      } else {
        myGenome[i,(indleftMarkerRecomb+1):nbMarkersAndQTL]=1	
      }
    
      posRecomb = posRecomb+rexp(1, 1)
    # draw the next recombination event
    
   }#end while
  
  }#end for
  return (myGenome)
}

#### Generate Phenotypes ##########
  
generateGumbelPhenotypes<-function(myGenome, nbTRN, mu, effectQTL, indQTL, errorVariance){ 
  
  #requires library(ordinal).
   
  myY = mu + myGenome[ , indQTL] %*% effectQTL  +  sqrt(errorVariance)*rgumbel(nbTRN, location = 0, scale = 1, max = TRUE)
  
  # So, Variance(Y)= errorVariance*pi^2/6
  return(myY)
  
}

##### Compute the statistical test along the genome #####
 
# Compute the test statistic along the genome using the interpolation ... 
 
###################### Gumbel distribution ###########################
 
## First, only on markers

computeTestStatIMOnlyOnMarkersGumbel<-function(myGenome,myY,dmark,nbMarkers,nbTRN){
  
  testStatOnMarkers = rep(0,nbMarkers)
  # statTest contains my test statistic computed along the genome
  
  # I start by computing the test statistic at markers , I fill my statTest vector
  for (k in 1:nbMarkers) { 
    
    Term =  (myY-mean(myY))*pi / sqrt(6*var(myY))[1]  + 0.5772156
    #testStatOnMarkers[k] = myGenome[ , indMarker[k]]  %*%  (1-exp(-(myY-mean(myY))/sqrt(var(myY)[1])))  /sqrt(nbTRN)
    
    testStatOnMarkers[k] = myGenome[ , indMarker[k]]  %*%  (1-exp(-Term) ) /sqrt(nbTRN)
    
    
  } 
  
  return(testStatOnMarkers)
} 

 

##### Compute the statistical test along the genome For the Gumbel Distribution
##### but using a maximum of one test within  marker intervals #####

computeTestStatIMSmartFormulaGumbel<-function(myGenome,nbTRN,dmark,nbMarkers,indMarker,posMarkers,myY){
   
  testStatOnMarkers = rep(0,nbMarkers)
  # statTest contains my test statistic computed along the genome
  
  # I start by computing the test statistic at markers , I fill my statTest vector
  for (k in 1:nbMarkers) { 
    
   Term =  (myY-mean(myY))*pi / sqrt(6*var(myY))[1]  + 0.5772156
   testStatOnMarkers[k] = myGenome[ , indMarker[k]]  %*%  (1-exp(-Term) ) /sqrt(nbTRN)
    
  } 
  # 0.5772156 is the euler mascheroni constant
  
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
 
#  No Plugin below
 

computeTestStatIMSmartFormulaGumbelNoPlugin<-function(myGenome,nbTRN,dmark,nbMarkers,indMarker,posMarkers,myY,mu,errorVariance){
  
  testStatOnMarkers = rep(0,nbMarkers)
  # statTest contains my test statistic computed along the genome
  
  # I start by computing the test statistic at markers , I fill my statTest vector
  for (k in 1:nbMarkers) { 
    
    testStatOnMarkers[k] = myGenome[ , indMarker[k]]  %*%  (1-exp(-(myY-mu)/sqrt(errorVariance)))  /sqrt(nbTRN)
    
  } 
  
  
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


 
############# Logistic distribution now ############################
 

##### Compute the statistical test along the genome For the Logistic Distribution
##### but using a maximum of one test within  marker intervals #####
 
computeTestStatIMSmartFormulaLogistic<-function(myGenome,nbTRN,dmark,nbMarkers,indMarker,posMarkers,myY){
  
  testStatOnMarkers = rep(0,nbMarkers)
  # statTest contains my test statistic computed along the genome
  
  #var(myY)= (pi^2 /3)* errorVariance
  #errorVariance =  var(myY) *3/pi^2
  
  
  # I start by computing the test statistic at markers , I fill my statTest vector
  for (k in 1:nbMarkers) { 
    
    myStat=0
    
    for (i in 1:nbTRN) { 
      
      num= myGenome[ i, indMarker[k]] %*% (exp(-( (myY[i]-mean(myY))*pi / sqrt(3*var(myY) )))-1 )
      denom= (1 + exp(-( (myY[i]-mean(myY))*pi / sqrt(3*var(myY) ))) ) * sqrt(nbTRN/3)
      myStat= myStat + num/denom
    }
    
    testStatOnMarkers[k] = myStat
    
  } 
    
  
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


#### No plugin now ##########
 
computeTestStatIMSmartFormulaLogisticNoPlugin<-function(myGenome,nbTRN,dmark,nbMarkers,indMarker,posMarkers,myY,mu,errorVariance){
  
  testStatOnMarkers = rep(0,nbMarkers)
  # statTest contains my test statistic computed along the genome
   
  
  # I start by computing the test statistic at markers , I fill my statTest vector
  for (k in 1:nbMarkers) { 
    
    myStat=0
    
    for (i in 1:nbTRN) { 
      
      num= myGenome[ i, indMarker[k]] %*% (exp(-( (myY[i]-mu) /sqrt(errorVariance)))-1 )
      denom= (1 + exp(-( (myY[i]-mu)  / sqrt(errorVariance))) ) * sqrt(nbTRN/3)
      myStat= myStat + num/denom
    }
    
    testStatOnMarkers[k] = myStat
    
  } 
  
  
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


  


################################ END STAND BY #################################
 

##########################################################


r <- function(posLeft,posRight){
  
  myResult =  0.5*(1-exp(-2*abs(posRight-posLeft))) 
  
  return(myResult)
  
  
}



Q <- function(alleleLeft,alleleRight,pos,posLeft,posRight){
  
  
  if ( (alleleLeft==1) && (alleleRight==1)) {
    
    myQ = (1-r(posLeft,pos))*(1-r(pos,posRight))/(1-r(posLeft,posRight))
    
  } else if ( (alleleLeft==1) && (alleleRight==-1)) {
    
    
    myQ = (1-r(posLeft,pos))*r(pos,posRight)/r(posLeft,posRight)
    
    
  } else if ( (alleleLeft==-1) && (alleleRight==1)) {
    
    myQ = r(posLeft,pos)*(1-r(pos,posRight))/r(posLeft,posRight)
    
  } else  {
    
    
    # (alleleLeft==-1) && (alleleRight==-1))
    
    myQ = r(posLeft,pos)*r(pos,posRight)/(1-r(posLeft,posRight))
    
  }
  
  return(myQ)
  
}







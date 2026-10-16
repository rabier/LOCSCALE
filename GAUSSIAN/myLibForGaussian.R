# My functions for Gaussian family

#### Generate Genomes ##########
generateGenomes<-function(nbInd, nbMarkersAndQTL){ 
  
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

generateGaussianPhenotypes<-function(mu, effectQTL, indQTL, errorVariance, myGenome){ 
   
  myY = mu + genomeTRN[ , indQTL] %*% effectQTL  + rnorm(nbTRN,0,sqrt(errorVariance))
     
  return(myY)

}

##### Compute the maximum of the process  #######
##### but using a maximum of one test within  marker intervals #####
#### according to Azais formula #####
 
computeTestStatIMSmartFormula<-function(dmark,nbMarkers,myGenome,myY){
  
  testStatOnMarkers = rep(0,nbMarkers)
  # statTest contains my test statistic computed along the genome
   
  # I start by computing the test statistic at markers
  for (k in 1:nbMarkers) { 
    testStatOnMarkers[k] = ( myGenome[ , indMarker[k]] %*% (myY-mean(myY)) ) / sqrt(sum((myY-mean(myY))^2))
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







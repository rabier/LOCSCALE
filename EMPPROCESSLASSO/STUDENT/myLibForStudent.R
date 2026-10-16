# My functions for the Student family


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
 
generateStudentPhenotypes<-function(myGenome, nbTRN, mu, effectQTL, indQTL, errorVariance, df){ 
  
  # df is the numer of degrees of freedom of the student distribution
  myY = mu + myGenome[ , indQTL] %*% effectQTL  +  sqrt(errorVariance)*rt(nbTRN, df, 0)
  
  # So, Variance(Y)= errorVariance*df/(df-2)
  return(myY)
  
}
 

### Student distribution

# First function, only on markers 
computeTestStatIMStudentOnlyOnMarkers<-function(myGenome,myY,dmark,nbMarkers,df){
  
  testStatOnMarkers = rep(0,nbMarkers)
  # statTest contains my test statistic computed along the genome
  
  # I start by computing the test statistic at markers, I fill my statTest vector
  for (k in 1:nbMarkers) { 
    
    myStat=0
    
    for (i in 1:nbTRN) { 
      num=( myGenome[ i, indMarker[k]] %*% (myY[i]-mean(myY))/ (sqrt(var(myY) )) ) * sqrt(df*(df+1)*(df+3)/(df-2))/sqrt(nbTRN) 
      denom=df +  ( (myY[i]-mean(myY))^2 / var(myY) )* df/(df-2)
      myStat= myStat + num/denom
    }
    
    testStatOnMarkers[k] = myStat
    
  }# TO DO Check that it works with this tets statistic ...
  
  return(testStatOnMarkers)
  
}    

 

##### Compute the statistical test along the genome For the Student Distribution
##### but using a maximum of one test within  marker intervals #####
 

computeTestStatIMSmartFormulaStudent<-function(myGenome,nbTRN,sizeChro,dmark,nbMarkers,indMarker,posMarkers,myY,df){
   
  testStatOnMarkers = rep(0,nbMarkers)
  # statTest contains my test statistic computed along the genome
  
  # I start by computing the test statistic at markers, I fill my statTest vector
  for (k in 1:nbMarkers) { 
    
    myStat=0
    
    for (i in 1:nbTRN) { 
      num=( myGenome[ i, indMarker[k]] %*% (myY[i]-mean(myY))/ (sqrt(var(myY) )) ) * sqrt(df*(df+1)*(df+3)/(df-2))/sqrt(nbTRN) 
      denom=df +  ( (myY[i]-mean(myY))^2 / var(myY) )* df/(df-2)
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


# IDEM but without the plugin, so still for Student distribution
 
computeTestStatIMSmartFormulaStudentNoPlugin<-function(myGenome,nbTRN,dmark,nbMarkers,myY,df,mu,errorVariance){
  
  # I am using mu and errorVariance since I am not using a plug in 
  
  testStatOnMarkers = rep(0,nbMarkers)
  # statTest contains my test statistic computed along the genome
  
  # I start by computing the test statistic at markers , I fill my statTest vector
  for (k in 1:nbMarkers) { 
     
    myStat=0
    
    for (i in 1:nbTRN) { 
    
      num=(myGenome[i, indMarker[k]] %*% (myY[i]-mu)/sqrt(errorVariance)) *  sqrt((df+1)*(df+3)) /sqrt(nbTRN)
      denom=df +  ( (myY[i]-mu)^2 / errorVariance )
      myStat= myStat + num/denom
      
    }
    
    testStatOnMarkers[k] =  myStat
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







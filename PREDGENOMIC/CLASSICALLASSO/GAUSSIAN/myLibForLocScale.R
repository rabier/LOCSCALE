# My functions for my location scale family code


#### Generate Genomes under Haldane ##########
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

generateGaussianPhenotypes<-function(myGenome, nbTRN, mu, effectQTL, indQTL, errorVariance){ 
   
  myY = mu + myGenome[ , indQTL] %*% effectQTL  + rnorm(nbTRN,0,sqrt(errorVariance))
     
  return(myY)

}
 


generateStudentPhenotypes<-function(myGenome, nbTRN, mu, effectQTL, indQTL, errorVariance, df){ 
  
  # df is the numer of degrees of freedom of the student distribution
  myY = mu + myGenome[ , indQTL] %*% effectQTL  +  sqrt(errorVariance)*rt(nbTRN, df, 0)
  
  # So, Variance(Y)= errorVariance*df/(df-2)
  return(myY)
  
}

 

generateGumbelPhenotypes<-function(myGenome, nbTRN, mu, effectQTL, indQTL, errorVariance){ 
  
  #requires library(ordinal).
   
  myY = mu + myGenome[ , indQTL] %*% effectQTL  +  sqrt(errorVariance)*rgumbel(nbTRN, location = 0, scale = 1, max = TRUE)
  
  # So, Variance(Y)= errorVariance*pi^2/6
  # ok it works also for E(Y), I have checked 
  return(myY)
  
}



generateLogisticPhenotypes<-function(myGenome, nbTRN, mu, effectQTL, indQTL, errorVariance){ 
  
  #requires library(stats)
   
  myY = mu + myGenome[ , indQTL] %*% effectQTL  +  sqrt(errorVariance)*rlogis(nbTRN, location = 0, scale = 1)
  
  # So, Variance(Y)= errorVariance*pi^2/3
  # ok it works also for E(Y), I have checked 
  return(myY)
  
}

 
 

##### Compute the statistical test along the genome #####

# First function, only on markers 
computeTestStatIMOnlyOnMarkers<-function(myGenome,myY,dmark,nbMarkers){
   
  testStatOnMarkers = rep(0,nbMarkers)
  # statTest contains my test statistic at the different markers along the genome
  for (k in 1:nbMarkers) { 
    testStatOnMarkers[k] = ( myGenome[ , indMarker[k]] %*% (myY-mean(myY)) ) / sqrt(sum((myY-mean(myY))^2))
  }
   
  return(testStatOnMarkers)
  
}
## end function, only on markers

# Compute the test statitic along the genome using the interpolation ... 

computeTestStatIM<-function(myGenome,myY,sizeChro,dtest,dmark,nbMarkers,indMarker,posMarkers){

  testedPositions=seq(0,sizeChro,by=dtest)

  nbTests=length(testedPositions)

  testStat = rep(0,nbTests)
  # statTest contains my test statistic computed along the genome

  indMarkerTestStat<-seq(1,nbTests,by=(dmark/dtest))

  indMarkerTestStat<-as.integer(round(indMarkerTestStat))

  # I start by computing the test statistic at markers , I fill my statTest vector
  for (k in 1:nbMarkers) { 
    # a remettre
    testStat[indMarkerTestStat[k]] = ( myGenome[ , indMarker[k]] %*% (myY-mean(myY)) ) / sqrt(sum((myY-mean(myY))^2))
    #testStat[indMarkerTestStat[k]] = ( myGenome[ , indMarker[k]] %*% (myY-0.57721) ) / sqrt(sum((myY-mean(myY))^2))
    # end a remettre 
    
    # j essaie juste un truc c est ici que ca joue pour la robustesse!
    #testStat[indMarkerTestStat[k]] = ( myGenome[ , indMarker[k]] %*% (myY) ) / sqrt(length(myY))
    ### en d j essaie 
    
    
  }

  # Now, I am handling the test statistics between markers 

  for (i in 1:(nbMarkers-1)){ 
   # deal with the ith marker interval
  
   posFirstTestWithin<-posMarkers[i]+dtest
   posLastTestWithin<-posMarkers[(i+1)]-dtest
   seqTestPosMarkerInterval <- seq(posFirstTestWithin,posLastTestWithin,by=dtest)
   nbTestMarkerInterval <- length(seqTestPosMarkerInterval) 
   
   for (j in 1:nbTestMarkerInterval){ 
   #jth test after marker i
    
    currentTestPos <- seqTestPosMarkerInterval[j]
    leftMarkerPos <- posMarkers[i]
    rightMarkerPos <- posMarkers[(i+1)]
    
    currentTestInd <- as.integer(round(indMarkerTestStat[i] + j))
    indLeftMarker<-indMarkerTestStat[i]
    indRightMarker<-indMarkerTestStat[(i+1)]
    
    alpha=0
    beta=0
    alpha = Q(1,1,currentTestPos,leftMarkerPos,rightMarkerPos) -  Q(-1,1,currentTestPos,leftMarkerPos,rightMarkerPos)
    beta = Q(1,1,currentTestPos,leftMarkerPos,rightMarkerPos) -  Q(1,-1,currentTestPos,leftMarkerPos,rightMarkerPos)
    num = alpha * testStat[indLeftMarker] +  beta * testStat[indRightMarker] 
    denom = sqrt( alpha^2 + beta^2 + 2*alpha*beta*exp(-2*dmark) )
    testStat[currentTestInd] = num /denom
    
   }
  
  
  }

  return(testStat)
  
}



##### Compute the statistical test along the genome
##### but using a maximum of one test within  marker intervals #####
# Fix it: place indmarker in argument 

computeTestStatIMSmartFormula<-function(myGenome,myY,dmark,nbMarkers){
  
  testStatOnMarkers = rep(0,nbMarkers)
  # statTest contains my test statistic computed along the genome
   
  # I start by computing the test statistic at markers , I fill my statTest vector
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

################################ STAND BY #################################

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

computeTestStatIMSmartFormulaStudent<-function(myGenome,nbTRN,sizeChro,dtest,dmark,nbMarkers,indMarker,posMarkers,myY,df){
   
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
# Fix it: place indmarker in argument 


computeTestStatIMSmartFormulaStudentNoPlugin<-function(myGenome,nbTRN,dmark,nbMarkers,myY,mu,errorVariance){
  
  # I am using mu since I am not using a plug in 
   
  
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
  }# TO DO Check that it works with this tets statistic ...
  
  
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
   #testStatOnMarkers[k] = myGenome[ , indMarker[k]]  %*%  (1-exp(-(myY-mean(myY))/sqrt(var(myY)[1])))  /sqrt(nbTRN)
 
   testStatOnMarkers[k] = myGenome[ , indMarker[k]]  %*%  (1-exp(-Term) ) /sqrt(nbTRN)
    
   
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
 
# WAIT No Plugin now 
 

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
 
 
  
computeTestStatIMLogisticOnlyOnMarkers<-function(myGenome,myY,dmark,nbMarkers,nbTRN){
    
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
    return(testStatOnMarkers)
    
}    



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







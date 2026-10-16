# Here are my functions useful for analyzing real data

 
##### new function for our real data based on Gaussian data ##################


##### Compute the statistical test along the genome #####

############## Gaussian assumption ###############

# Compute the test statistic along the genome using the interpolation ... 

computeTestStatIM<-function(myGenome,myY,sizeChro,dtest,posMarkers,nbMarkers){
  
  testedPositions=seq(posMarkers[1],posMarkers[1]+sizeChro,by=dtest)
  nbTests=length(testedPositions)
  testStat = rep(0,nbTests)
  print(nbTests)
  print(length(testStat))
  # statTest contains my test statistic computed along the genome
  
  indMarkerTestStat<- posMarkers/dtest - posMarkers[1]/dtest +1
  print(indMarkerTestStat)
  indMarkerTestStat<-as.integer(round(indMarkerTestStat))
  print(indMarkerTestStat)
  
  # I start by computing the test statistic at markers , I fill my statTest vector
  for (k in 1:nbMarkers) { 
    
    testStat[indMarkerTestStat[k]] = ( myGenome[ ,  k] %*% (myY-mean(myY)) ) / sqrt(sum((myY-mean(myY))^2))
    
  }
  print(length(testStat))
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
      denom = sqrt( alpha^2 + beta^2 + 2*alpha*beta*exp(-2*(rightMarkerPos-leftMarkerPos)))
      testStat[currentTestInd] = num /denom
      
    }
    
    
  }
  print(length(testStat))
  return(testStat)
  
}


 
 
computeTestStatIMSmartFormula<-function(myGenome,myY,posMarkers,nbMarkers){
  
  testStatOnMarkers = rep(0,nbMarkers)
  # statTest contains my test statistic computed along the genome
  
  # I start by computing the test statistic at markers , I fill my statTest vector
  for (k in 1:nbMarkers) { 
    testStatOnMarkers[k] = ( myGenome[ , k] %*% (myY-mean(myY)) ) / sqrt(sum((myY-mean(myY))^2))
  }
  
  currentMax=testStatOnMarkers[1]^2
  # Now a maximum of one test between markers
  
  for (i in 1:(nbMarkers-1)){ 
    # deal with the ith marker interval
    
    currentMax=max(currentMax,testStatOnMarkers[(i+1)]^2)
    
    ratio=testStatOnMarkers[(i+1)]/testStatOnMarkers[(i)]
    # dmark is the distance between markers of interest 
    # (i.e. the marker interval under study)
    
    print(i)
    dmark=posMarkers[(i+1)]- posMarkers[(i)]
    
    if ( (ratio > exp(-2*dmark)) && (ratio < exp(2*dmark)) ) { 
      
      numSuppStat = testStatOnMarkers[(i+1)]^2 + testStatOnMarkers[i]^2 - 2*exp(-2*dmark)*testStatOnMarkers[(i+1)]*testStatOnMarkers[i]
      denomSuppStat = 1- exp(-4*dmark)
      suppStat = numSuppStat/denomSuppStat
      currentMax= max(currentMax,suppStat)
      
    }#end if
    
  }# end for
  
  return(c(currentMax,testStatOnMarkers))
  
}

###################
  
################## Student assumption ##############

### Student with interpolation ######


computeTestStatStudentIM<-function(myGenome,myY,sizeChro,dtest,posMarkers,nbMarkers,df){
  
  testedPositions=seq(posMarkers[1],posMarkers[1]+sizeChro,by=dtest)
  nbTests=length(testedPositions)
  testStat = rep(0,nbTests)
  print(nbTests)
  print(length(testStat))
  # statTest contains my test statistic computed along the genome
  
  indMarkerTestStat<- posMarkers/dtest - posMarkers[1]/dtest +1
  print(indMarkerTestStat)
  indMarkerTestStat<-as.integer(round(indMarkerTestStat))
  print(indMarkerTestStat)
  
  
  # I start by computing the test statistic at markers, I fill my statTest vector
  for (k in 1:nbMarkers) { 
    
    myStat=0
    
    for (i in 1:nbTRN) { 
      num=( myGenome[ i, k] %*% (myY[i]-mean(myY))/ (sqrt(var(myY) )) ) * sqrt(df*(df+1)*(df+3)/(df-2))/sqrt(nbTRN) 
      denom=df +  ( (myY[i]-mean(myY))^2 / var(myY) )* df/(df-2)
      myStat= myStat + num/denom
    }
    
    testStat[indMarkerTestStat[k]] = myStat
    
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
      denom = sqrt( alpha^2 + beta^2 + 2*alpha*beta*exp(-2*(rightMarkerPos-leftMarkerPos)))
      testStat[currentTestInd] = num /denom
      
    }
    
    
  }
  print(length(testStat))
  return(testStat)
  
}

#### End Student with interpolation ########

### Smart Formula for Student #######


computeTestStatIMSmartFormulaStudent<-function(myGenome,nbTRN,nbMarkers,posMarkers,myY,df){
  
  testStatOnMarkers = rep(0,nbMarkers)
  # statTest contains my test statistic computed along the genome
  
  # I start by computing the test statistic at markers, I fill my statTest vector
  for (k in 1:nbMarkers) { 
    
    myStat=0
    
    for (i in 1:nbTRN) { 
      num=( myGenome[ i, k] %*% (myY[i]-mean(myY))/ (sqrt(var(myY) )) ) * sqrt(df*(df+1)*(df+3)/(df-2))/sqrt(nbTRN) 
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
    
    dmark=posMarkers[(i+1)]- posMarkers[(i)]
    
    if ( (ratio > exp(-2*dmark)) && (ratio < exp(2*dmark)) ) { 
      
      numSuppStat = testStatOnMarkers[(i+1)]^2 + testStatOnMarkers[i]^2 - 2*exp(-2*dmark)*testStatOnMarkers[(i+1)]*testStatOnMarkers[i]
      denomSuppStat = 1- exp(-4*dmark)
      suppStat = numSuppStat/denomSuppStat
      currentMax= max(currentMax,suppStat)
      
    }#end if
    
  }# end for
  
  return(currentMax)
  
}
  

########## End Student #####################

  
##################### Useful functions #####################


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
 
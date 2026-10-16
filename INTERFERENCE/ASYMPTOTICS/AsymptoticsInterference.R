
# Deal with asymptotics under the interference model 

source("myLibForInterference.R")
 
#### Generate OU process under H1 at marker locations ##########
simDiscretizedOrnsteinInterferenceH1<-function(nbTRN, dmark, nbMarkers, posMarkers, posMarkersAndQTL, indQTL, effectQTL, Psi0, errorVariance){ 
   
  # generate a discretized OU process on marker locations
  X = rep(0,nbMarkers)
  Z = rep(0,nbMarkers)
  # Z contains my test statistic computed along the genome
  
for (i in 1:nbMarkers){  
  # deal with ith marker
  Shift<-0
  
  for (k in 1:nbQTLs){
    # deal with kth QTL
  
   myInd<-indQTL[k]
   myPosQTL<-posMarkersAndQTL[myInd]
   myIndFlankingMarkerLeftSideQTL<-myInd-1
   myIndFlankingMarkerRightSideQTL<-myInd+1
   myPosMarkerLeftSideQTL<-posMarkersAndQTL[myIndFlankingMarkerLeftSideQTL]
   myPosMarkerRightSideQTL<-posMarkersAndQTL[myIndFlankingMarkerRightSideQTL]
  
   coeff<- sqrt(Psi0) *sqrt(nbTRN)*effectQTL[k]/sqrt(errorVariance)  
   myAlphaTilde<-alphaTilde(posQTL[k],myPosMarkerLeftSideQTL,myPosMarkerRightSideQTL)
   myBetaTilde<-betaTilde(posQTL[k],myPosMarkerLeftSideQTL,myPosMarkerRightSideQTL)
     
   if (myPosQTL>posMarkers[i]){
   
     Term<-exp(-2*abs(posMarkers[i]-myPosMarkerLeftSideQTL)) * 
      (myAlphaTilde + myBetaTilde* exp(-2*abs(myPosMarkerRightSideQTL-myPosMarkerLeftSideQTL)))
   } else {
   
     Term<-exp(-2*abs(posMarkers[i]-myPosMarkerRightSideQTL)) *  
       (myAlphaTilde * exp(-2*abs(myPosMarkerRightSideQTL-myPosMarkerLeftSideQTL)) + myBetaTilde)     
   
   }
   Shift<-coeff*Term + Shift
     
  }# end k
     
     if (i==1){
       #first marker
       X[1]=rnorm(1,0,1);
       Z[1]=X[1] + Shift
     }else{
       epsilon= sqrt(1-exp(-4*dmark))*rnorm(1,0,1);
       X[i]=exp(-2*dmark)*X[(i-1)]+ epsilon;
       Z[i]=X[i] + Shift
     } # Z is my asymptotic test statistic at markers
        
  } # end i
    
  return(Z)
 
} # end function

 
  
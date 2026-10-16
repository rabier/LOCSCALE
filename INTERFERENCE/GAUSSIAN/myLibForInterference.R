# Functions for interference

#### Generate Genomes under interference ##########

# We assume that the QTLs belong to different marker interval 
generateGenomesInterference<-function(genomeTRN, nbInd, posMarkersAndQTL,indQTL, nbQTLs){ 
  
# input: genomeTRN contains a genome generated under Haldane
# output: myGenomeTRNInterf will contain a genome under interference  
# Recall that at marker locations we have Haldane model, so no changes needed at markers  
# so we will update only QTL alleles   
  
myGenomeTRNInterf<-genomeTRN
  
for (i in 1:nbInd){   
  # individual i
  
  
  for (k in 1:nbQTLs){
    # we deal with the k th QTL
  
     myInd<-indQTL[k]
     myPosQTL<-posMarkersAndQTL[myInd]
     myIndMarkerLeft<-myInd-1
     myIndMarkerRight<-myInd+1
     myPosMarkerLeft<-posMarkersAndQTL[myIndMarkerLeft]
     myPosMarkerRight<-posMarkersAndQTL[myIndMarkerRight]
     
     if ( (genomeTRN[i,myIndMarkerLeft]==1) & (genomeTRN[i,myIndMarkerRight]==1) ) {
       
       #print("11")
       myGenomeTRNInterf[i,myInd]=1
       
     } else if ( (genomeTRN[i,myIndMarkerLeft]==-1) & (genomeTRN[i,myIndMarkerRight]==-1) ) {
  
       #print("-1-1")
       myGenomeTRNInterf[i,myInd]=-1
    
     } else if ( (genomeTRN[i,myIndMarkerLeft]==1) & (genomeTRN[i,myIndMarkerRight]==-1) ) {
       
       #print("1-1")
       # generate random number from uniform distribution
       if (runif(1, 0, 1)< (myPosMarkerRight-myPosQTL)/(myPosMarkerRight-myPosMarkerLeft)){
         myGenomeTRNInterf[i,myInd]=1
       } else {
         myGenomeTRNInterf[i,myInd]=-1
       }
       
     } else {
       #print("-11")
       
       #we have ( (genomeTRN[i,myIndMarkerLeft]==-1) & (genomeTRN[i,myIndMarkerRight]==1) )
       # generate random number from uniform distribution
       if (runif(1, 0, 1)< (myPosQTL-myPosMarkerLeft)/(myPosMarkerRight-myPosMarkerLeft)){
         myGenomeTRNInterf[i,myInd]=1
       } else {
         myGenomeTRNInterf[i,myInd]=-1
       }
     
     } # end if 
     
    
  } # end for k
} # end for i
 
  return(myGenomeTRNInterf)
}


alphaTilde<-function(t, tleft, tright){ 
  #compute alphaTilde for interference
  result<-(t-tleft)/(tright-tleft)
  return(result)
}  
  

betaTilde<-function(t, tleft, tright){ 
  #compute betaTilde for interference
  result<-(tright-t)/(tright-tleft)
  return(result)
}  


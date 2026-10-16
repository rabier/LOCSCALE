#  Empirical Process Lasso 
# I consider here Logistic case and focus on QTL mapping
# You can choose below the setting for the QTL effects

library(glmnet)

source("myLibForLogistic.R")
set.seed(12)

 
calcLASSO = function (Y,X) {
  
  fit1<-cv.glmnet(X[,],Y,intercept=FALSE)
  
  MyCoeff<-coef(fit1,s=fit1$lambda.1se)
  #lambda.1se
  #lambda.min 
   
  return (MyCoeff)
}
 
##################### Input parameters ####################

sizeChro=4
# sizeChro: chromosome length in Morgan

dmark=0.001
# dmark: distance between markers, we suppose equally spaced 
 
nbTRN=2000
# number of training individuals

nbMarkers= sizeChro/dmark +1
# nbMarkers: number of markers

# QTL locations (same locations as in our AdaptSgenoLasso paper)
posQTL=c(0.65001, 1.50001, 2.35001, 2.75001, 3.10001, 3.75001)


# QTL effects : same QTL effect as in our AdaptSgenoLasso paper 

# First setting: 6 large QTLs
effectQTL=c(0.3794,0.3794,0.3794,0.3794,0.3794,0.3794)

# Second setting: 3 large QTLs and 3 small QTLs
#effectQTL=c(0.1897,0.3794,0.1897,0.3794,0.1897,0.3794)
# effectQTL: vector of QTL effects

# Third setting: 1 large QTL and 5 small QTLs
#effectQTL=c(0.1897,0.3794,0.1897,0.18974,0.1897,0.1897)

nbQTLs=length(posQTL)

mu=0
# mu is the global mean
errorVariance=1
# errorVariance is the environmental variance

nbSample=100

############# End of input parameters ###############

# INITIALIZATION
nbMarkersAndQTL = nbMarkers + nbQTLs
posMarkers = seq(0,sizeChro,dmark)
posMarkersAndQTL = sort(c(posMarkers,posQTL))
indQTL<- which(posMarkersAndQTL %in% posQTL)
# posMarkersAndQTL[indQTL] gives the QTL locations 
indMarker<- which(posMarkersAndQTL %in% posMarkers)
# END INITIALIZATION

########### Calculate the regression matrix for the Empirical Process Lasso
V<-matrix(0,nrow=nbMarkers,ncol=nbMarkers)

for (j in 1:nbMarkers){
  for (i in 1:j){
    V[i,j]=exp(-2*(posMarkersAndQTL[indMarker[j]] -posMarkersAndQTL[indMarker[i]]))
  }    
}

matSigma<-(V+t(V)) - diag(1,nrow=nbMarkers,ncol=nbMarkers)
Abis<-chol(matSigma)
A<-t(Abis)

########### End Regression matrix 


##### Setting the scene for counting QTLs in QTL mapping ######

# Set up for counting QTL correctly identified 
#### if I allow of having a distance of either 1cM or 2 cM 
## on the left side and right side
  
nbInd1cM<-0.01/dmark

# 1cM
indQTLSmart<-indQTL
for (j in 1:nbQTLs){
  Left<-seq(indQTL[j]-nbInd1cM,indQTL[j]-1,1)
  Right<-seq(indQTL[j]+1,indQTL[j]+nbInd1cM,1)
  indQTLSmart<-c(indQTLSmart,Left,Right)
}
indQTLSmart<-sort(indQTLSmart)
 
# 2cM
nbInd2cM<-0.02/dmark
indQTLSmart2cM<-indQTL
for (j in 1:nbQTLs){
  Left<-seq(indQTL[j]-nbInd2cM,indQTL[j]-1,1)
  Right<-seq(indQTL[j]+1,indQTL[j]+nbInd2cM,1)
  indQTLSmart2cM<-c(indQTLSmart2cM,Left,Right)
}
indQTLSmart2cM<-sort(indQTLSmart2cM)
   
L1Ratio<-rep(0,nbSample)
L1RatioSmart<-rep(0,nbSample)
L1RatioSmart2cM<-rep(0,nbSample)
 
################# End setting the scene ###############

# Now let us compute the EmpProcessLasso ...

for (indsimul in 1:nbSample){

  genomeTRN = generateGenomes(nbTRN,nbMarkersAndQTL,sizeChro)
  # genomeTRN : genome of the TRN individuals of size nbTRN times nbMarkersAndQTL
  yTRN = generateLogisticPhenotypes(genomeTRN, nbTRN, mu, effectQTL, indQTL, errorVariance)
  myObservedTestStatAtMarkers<-computeTestStatIMLogisticOnlyOnMarkers(genomeTRN,yTRN,dmark,nbMarkers,nbTRN)
  YEmpLasso<-solve(A)%*%myObservedTestStatAtMarkers 
  CoeffLasso<-calcLASSO(YEmpLasso,t(A))
  indNonNulCoeffLasso <-which(CoeffLasso!=0)
  TruePositives<-which(indNonNulCoeffLasso %in% indQTL )
  TruePositivesSmart<-which(indNonNulCoeffLasso %in% indQTLSmart )
  TruePositivesSmart2cM<-which( indNonNulCoeffLasso  %in% indQTLSmart2cM )  
  
  L1Ratio[indsimul]<- sum(abs(CoeffLasso[indNonNulCoeffLasso[TruePositives]]))  / sum(abs(CoeffLasso)) 
  L1RatioSmart[indsimul]<- sum(abs(CoeffLasso[indNonNulCoeffLasso[TruePositivesSmart]]))  / sum(abs(CoeffLasso)) 
  L1RatioSmart2cM[indsimul]<- sum(abs(CoeffLasso[indNonNulCoeffLasso[TruePositivesSmart2cM]]))  / sum(abs(CoeffLasso)) 
  
}

 
mean(L1Ratio[!is.na(L1Ratio)])

#L1Ratio(0.01)
mean(L1RatioSmart[!is.na(L1RatioSmart)])

#L1Ratio(0.02)
mean(L1RatioSmart2cM[!is.na(L1RatioSmart2cM)])





 





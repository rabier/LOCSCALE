#  EmpProcessLasso for Student case,  in the context of genomic prediction  

library(glmnet)

source("myLibForStudent.R")
set.seed(12)

calcLASSO = function (Y,X) {
  
  fit1<-cv.glmnet(X[,],Y,intercept=FALSE)
  
  #MyCoeff<-coef(fit1,s=fit1$lambda.min )
  #lambda.1se
  #lambda.min 
   
  return (fit1$lambda.1se)
  
}
  
sizeChro=4
# sizeChro: chromosome length in Morgan

dmark=0.001
# dmark: distance between markers, we suppose equally spaced 
 
nbTRN=2000
# number of training individuals

nbTST=100
# number of test individuals

nbMarkers= sizeChro/dmark +1
# nbMarkers: number of markers

# Same QTL location as in our AdaptSgenoLasso paper
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

# Student distribution
df=4

### End of input parameters ###

# INITIALIZATION
nbMarkersAndQTL = nbMarkers + nbQTLs
posMarkers = seq(0,sizeChro,dmark)
posMarkersAndQTL = sort(c(posMarkers,posQTL))
indQTL<- which(posMarkersAndQTL %in% posQTL)
# posMarkersAndQTL[indQTL] gives our QTL locations 
indMarker<- which(posMarkersAndQTL %in% posMarkers)
# END INITIALIZATION

########### Calculate the regression matrix for our Empirical Lasso
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
   
############# Genomic prediction ###########

accuracyPredGenomic<-rep(0,nbSample)
 
# In the Student case 
Psi0= (df+1)/(df+3)

# Now let us compute the Lasso ...

for (indsimul in 1:nbSample){
   
  genomeTRN = generateGenomes(nbTRN,nbMarkersAndQTL,sizeChro)
  # genomeTRN : genome of the TRN individuals of size nbTRN times nbMarkersAndQTL
  yTRN =generateStudentPhenotypes(genomeTRN, nbTRN, mu, effectQTL, indQTL, errorVariance, df)
   
  # for TEST now
  genomeTST = generateGenomes(nbTST,nbMarkersAndQTL,sizeChro)
  yTST=generateStudentPhenotypes(genomeTST, nbTST, mu, effectQTL, indQTL, errorVariance, df)
  # end for TEST
  
  myObservedTestStatAtMarkers<-computeTestStatIMStudentOnlyOnMarkers(genomeTRN,yTRN,dmark,nbMarkers,df)
  YEmpLasso<-solve(A)%*%myObservedTestStatAtMarkers 
  optimalLambda<-calcLASSO(YEmpLasso,t(A))
    
  # Genomic prediction now
  myFit<-glmnet(t(A),YEmpLasso,intercept=FALSE, alpha=1, lambda=optimalLambda)
  MyCoeff<-coef(myFit)
  indNotNull<-which(MyCoeff!=0)
  YHat<- genomeTST[,indNotNull] %*% (MyCoeff[indNotNull,]/sqrt(Psi0*nbTRN))
 
  
  # Precision de prediction pour le LASSO , en cours 
  accuracyPredGenomic[indsimul]=cor(YHat,yTST)
  
  # End Genomic prediction now
  
}

 
mean(accuracyPredGenomic[!is.na(accuracyPredGenomic)])
 



 





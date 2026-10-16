#  Classical Lasso ... 
# I consider here logistic case

library(glmnet)

library(stats)


source("myLibForLocScale.R")
set.seed(20)

 

calcLASSO = function (Y,X) {
  
  fit1<-cv.glmnet(X[,],Y,intercept=FALSE)
  
  #MyCoeff<-coef(fit1,s=fit1$lambda.min )
  #lambda.1se
  #lambda.min 
  
  return (fit1$lambda.1se)
  
}

 
### end fonctions interessantes ####

sizeChro=4
# sizeChro: chromosome length in Morgan

dmark=0.001
# dmark: distance between markers, we suppose equally spaced 
 

nbTRN=1000
# number of training individuals

nbTST=100

nbMarkers= sizeChro/dmark +1
# nbMarkers: number of markers

# Same QTL location as in our AdaptSgenoLasso paper
posQTL=c(0.65001, 1.50001, 2.35001, 2.75001, 3.10001, 3.75001)

nbQTLs=length(posQTL)

# Same QTL effect as in our AdaptSgenoLasso paper 
effectQTL=c(0.1897,0.3794,0.1897,0.3794,0.1897,0.3794)
# effectQTL: vector of QTL effects

mu=0
# mu is the global mean
errorVariance=1
# errorVariance is the environmental variance

nbSample=100
### End of input parameters ###

# INITIALIZATION
nbMarkersAndQTL = nbMarkers + nbQTLs
posMarkers = seq(0,sizeChro,dmark)
posMarkersAndQTL = sort(c(posMarkers,posQTL))
indQTL<- which(posMarkersAndQTL %in% posQTL)
# posMarkersAndQTL[indQTL] gives our QTL locations 
indMarker<- which(posMarkersAndQTL %in% posMarkers)
# END INITIALIZATION
 
accuracyPredGenomic<-rep(0,nbSample)

# In the logistic case 


# End set up for counting QTL correctly identified 

 
# Now let us compute the Lasso ...

for (indsimul in 1:nbSample){
  print(indsimul)
  genomeTRN = generateGenomes(nbTRN,nbMarkersAndQTL,sizeChro)
  # genomeTRN : genome of the TRN individuals of size nbTRN times nbMarkersAndQTL
 
  yTRN = generateLogisticPhenotypes(genomeTRN, nbTRN, mu, effectQTL, indQTL, errorVariance)
  # for TEST now
  genomeTST = generateGenomes(nbTST,nbMarkersAndQTL,sizeChro)
  yTST=generateLogisticPhenotypes(genomeTST, nbTST, mu, effectQTL, indQTL, errorVariance)
  # end for TEST  
   
  optimalLambda<-calcLASSO(yTRN,genomeTRN)
  
  # Genomic prediction now
  myFit<-glmnet(genomeTRN,yTRN,intercept=TRUE, alpha=1, lambda=optimalLambda)
  MyCoeff<-coef(myFit)
  indNotNull<-which(MyCoeff!=0)
  YHat<- genomeTST[,indNotNull] %*% MyCoeff[indNotNull,]
  
  # Precision de prediction pour le LASSO , en cours 
  accuracyPredGenomic[indsimul]=cor(YHat,yTST)
  
   
}
 
mean(accuracyPredGenomic[!is.na(accuracyPredGenomic)])

 


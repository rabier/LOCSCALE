# Source code for Git in the Gaussian setting, Table 3 and Table 7 of the manuscript
 
source("myLibForLogistic.R")
set.seed(12)

library(stats)

######################### Input parameters #############################

sizeChro=1
# sizeChro: chromosome length in Morgan

nbQTLs=2
# nbQTLs: number of QTLs

#dmark: distance between markers
# for 6 markers
dmark=0.20 
# for 21 markers
#dmark=0.05
# for 51 markers
#dmark=0.02
# for 101 markers
#dmark=0.01
# for 1001 markers
#dmark=0.1
#dmark: distance between markers
#we suppose markers equally spaced
 

nbTRN=100
# nbTRN: number of training individuals

nbMarkers= sizeChro/dmark +1
# nbMarkers: number of markers

posQTL=c(0.75001, 0.90011)
# posQTL: vector of QTL locations 
# QTL locations must not match exactly marker locations

# effectQTL: vector of QTL effects
# H0
#effectQTL=c(0,0)
# I set the QTL effects to 0

#H1
effectQTL=c(4/sqrt(nbTRN),0) 
# Be careful, in the case  you have only one QTL, 
# you need to specify in this code
# two QTLs (nbQTLs=2) and set the second effect (or first effect) to 0
# As soon as you have at least 2 QTLs, it is always fine 

mu=0
# mu is the global mean
errorVariance=1
# errorVariance is the environmental variance
 

thresholdAzais=6.76
# Threshold for QTL mapping, to be replaced by the value associated to the marker density
# It is the threshold according to Azais et al formula

# Recall that in this tutorial sizeChro=1 and we consider 6 markers

thresholdDavies=6.92
thresholdFeingold=6.78
 
 
#################### End of input parameters ######################

nbMarkersAndQTL = nbMarkers + nbQTLs
posMarkers = seq(0,sizeChro,dmark)
posMarkersAndQTL = sort(c(posMarkers,posQTL))

indQTL<- which(posMarkersAndQTL %in% posQTL)
# posMarkersAndQTL[indQTL] gives our QTL locations 
indMarker<- which(posMarkersAndQTL %in% posMarkers)
 
countRejections<-0
nbSample=10000 
mySmartMax<-rep(0,nbSample)
mySmartMaxNoPlugin<-rep(0,nbSample)

for (i in 1:nbSample){
#  print(i)
  genomeTRN = generateGenomes(nbTRN,nbMarkersAndQTL)
  # genomeTRN : genome of the TRN individuals of size nbTRN times nbMarkersAndQTL
  yTRN = generateLogisticPhenotypes(genomeTRN, nbTRN, mu, effectQTL, indQTL, errorVariance)

  mySmartMax[i]<-computeTestStatIMSmartFormulaLogistic(genomeTRN,nbTRN ,dmark,nbMarkers,indMarker,posMarkers,yTRN)
  # compute the smart formula for the max of the process, according to Azais et al (2014)
  # i.e. tests only on markers and a signle test within marker intervals
      
  mySmartMaxNoPlugin[i]<-computeTestStatIMSmartFormulaLogisticNoPlugin(genomeTRN,nbTRN,dmark,nbMarkers,indMarker,posMarkers,yTRN,mu,errorVariance)
  # same as above but we do not use a plug in (ie we use true value of parameters)

  }

print("Number of rejections")

print("Azais formula")
print("Plugin")
sum(mySmartMax>thresholdAzais)
print("No Plugin")
sum(mySmartMaxNoPlugin>thresholdAzais)
 
print("Davies formula")
print("Plugin")
sum(mySmartMax>thresholdDavies)
print("No Plugin")
sum(mySmartMaxNoPlugin>thresholdDavies)

print("Feingold formula")
print("Plugin")
sum(mySmartMax>thresholdFeingold)
print("No Plugin")
sum(mySmartMaxNoPlugin>thresholdFeingold)

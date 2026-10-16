# Source code for Git 
# Robustness: we believe that it is Gaussian but it is Logistic

# We use a test statistic without any plugin 
source("myLibForGaussianNoPlugin.R")
source("myLibForLogistic.R")
set.seed(12)

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
 

for (i in 1:nbSample){
#  print(i)
  genomeTRN = generateGenomes(nbTRN,nbMarkersAndQTL,sizeChro)
  # genomeTRN : genome of the TRN individuals of size nbTRN times nbMarkersAndQTL
    
  yTRN =generateLogisticPhenotypes(genomeTRN, nbTRN, mu, effectQTL, indQTL, errorVariance)
  # yTRN : vector of size nbTRN containing the phenotypes Y
   
  mySmartMax[i]<-computeTestStatIMSmartFormulaNoPlugin(dmark,nbMarkers,genomeTRN,yTRN)
  # compute the smart formula for the max of the process, according to Azais et al (2014)
  # i.e. tests only on markers and a signle test within marker intervals
  
}

print("Number of rejections")

print("Azais formula")
sum(mySmartMax>thresholdAzais)
 
print("Davies formula")
sum(mySmartMax>thresholdDavies)

print("Feingold formula")
sum(mySmartMax>thresholdFeingold)

 

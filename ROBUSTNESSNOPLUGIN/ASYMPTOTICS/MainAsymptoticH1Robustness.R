# My main code for the asymptotic process for robustness
# We believe that it is Gaussian but it is either Student, Logistic or Gumbel
 
# Be careful, you have to set the value of Psi1 and Psi2 to the value
# corresponding to the family under study
# Psi0 is always equal to 1 because of Gaussian

# This tutorial is ready for Sudent
 
# Robustness
 
source("myLibForAsymptotics.R")
source("myLibForAsymptoticsRobustness.R")

set.seed(15)

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

# For Student distribution, specify the number df of degrees of freedom
df=4

#################### End of input parameters ######################

nbMarkersAndQTL = nbMarkers + nbQTLs
posMarkers = seq(0,sizeChro,dmark)
posMarkersAndQTL = sort(c(posMarkers,posQTL))
indQTL<- which(posMarkersAndQTL %in% posQTL)
# posMarkersAndQTL[indQTL] gives our QTL locations 
indMarker<- which(posMarkersAndQTL %in% posMarkers)
  
# number of paths of the asymptotic process 
nbPaths=100000
mySmartMaxAsympt<-rep(0,nbPaths)
myProcessAsymptAtMarkers<-rep(0,nbMarkers)
 
 
# For Student  
Psi1=df/(df-2)
Psi2=2/(df-2)
Psi0=1


# For Gumbel , cela marche pas pour linstant
#Psi1=(pi^2)/6 + 0.57721^2
#Psi2=(pi^2)/6 + 0.57721^2 -1
#Psi0=1

# For Logistic
#Psi1=(pi^2)/3
#Psi2=(pi^2)/3-1
#Psi0=1
 
 
for (i in 1:nbPaths){
 
  myProcessAsymptAtMarkers<-simDiscretizedProcessH1(nbTRN, dmark, nbMarkers, sizeChro, posMarkers, posQTL, effectQTL, Psi2, Psi1, Psi0, errorVariance) 
   
  mySmartMaxAsympt[i]<-computeAsymptProcessSmartFormula(myProcessAsymptAtMarkers, dmark, nbMarkers)
  
}  
 
print("AsymptoticProcessAzais")
sum(mySmartMaxAsympt>thresholdAzais)

print("AsymptoticProcessDavis")
sum(mySmartMaxAsympt>thresholdDavies)

print("AsymptoticProcessFeingold")
sum(mySmartMaxAsympt>thresholdFeingold)

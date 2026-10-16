# My main code for the interference phenomenon under Gaussian setting
 
source("myLibForGaussian.R")
source("myLibForInterference.R")
  
######################### Input parameters #############################

set.seed(12)
 
#sizeChro=1
# sizeChro: chromosome length in Morgan
#sizeChro=4
sizeChro=7
 
nbQTLs=2
# nbQTLs: number of QTLs

#dmark=0.2
dmark=0.25
# dmark: distance between markers, we suppose equally spaced
 
nbTRN=1000
#nbTRN=100
# nbTRN: number of training individuals

nbMarkers= sizeChro/dmark +1
# nbMarkers: number of markers

# posQTL: vector of QTL locations 
# QTL locations must not match exactly marker locations
# Be careful under interference QTL must not belong to same marker intervals
#posQTL=c(0.10001, 0.750011)
posQTL=c(0.10001, 2.800011)

#effectQTL=c(0,0) 
effectQTL=c(0,4/sqrt(nbTRN))
# effectQTL: vector of QTL effects

mu=0
# mu is the global mean

errorVariance=1
# errorVariance is the environmental variance


#thresholdAzais=6.76
#thresholdAzais=9.01
thresholdAzais=10.12
# Threshold for QTL mapping, to be replaced by the value associated to the marker density
# It is the threshold according to Azais et al formula
 

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
  #print(i)
  
  # we begin by generating genomes under Haldane
  genomeTRN = generateGenomes(nbTRN,nbMarkersAndQTL,sizeChro)
  # now for QTL Location, we generate under interference
  newGenomeTRN=generateGenomesInterference(genomeTRN,nbTRN, posMarkersAndQTL,indQTL, nbQTLs)
  # genomeTRN : genome of the TRN individuals of size nbTRN times nbMarkersAndQTL

  yTRN =generateGaussianPhenotypes(mu, effectQTL, indQTL, errorVariance,newGenomeTRN,nbTRN)
  # yTRN : vector of size nbTRN containing the phenotypes Y
  mySmartMax[i]<-computeTestStatIMSmartFormula(dmark,nbMarkers,newGenomeTRN,yTRN)
  # compute the smart formula for the max of the process
} 


print("Number of rejections")

print("Azais formula")
sum(mySmartMax>thresholdAzais)
#sum(myIMStat>6.76)


# map 4
#sum(mySmartMax>8.25)
#sum(myIMStat>8.25)
 
## Check empirical level

# NBTRN=200 , 5.1% donc niveau respecté pour T=1 d=0.2 Gaussienne avec seuil 6.76

# NBTRN=1000, 4.8% de rejet donc niveau respecté pour T=4 d=0.5 Gaussienne
# map 4 du papier interference. avec un seuil de 8.25
 
### End Check the empirical level for Gaussian #######
 
### Check the power ######### 
# pour la puissance empirique en interference j ai pas exactement les resultats de TEST

# map 4 du papier interference avec t*=0.25
# NBTRN=1000, 68.2% de rejet 

# 70 % pour autre simul
# 69.4 % avec 10000 traj
 

# NBTRN=1000,   68.31% de rejet

# nBTRN = 10 000, 69.84% de rejet

# new 75.09 % avec 10000 traj  en asymptotique , en accord avec le papier de Test 74.20%

# Si je fais Haldane, j ai du 52.61% en empirique avec nBTRN = 10 000


# map 4 du papier interference avec t*=0.60
# 79.2 % avec 1000  traj en empirique , nbTRN=1000


# 76.91% de power avec nBTRBN=1000
# 77.62% en empirique  avec nBTRBN=1000 et 10000 samples

# New 

# 77.93% en empirique  avec nBTRBN=1000 et 10000 samples
# 80.3% en empirique  avec nBTRBN=2000 et 1000 samples
# 77.45% en empirique  avec nBTRBN=5000 et 10000 samples

# 82.30% en processus asymptotique avec 10 000 paths , en accord avec le papier ds Test 82.01%

### End Check the power #########
  
 
  ################################# Asymptotic process ###########
  # L asymptotique marche ! je trouve comme dans le papier TEST
 
  nbPaths=100
  mySmartMaxAsympt<-rep(0,nbPaths)
  myProcessAsymptAtMarkers<-rep(0,nbMarkers)
  myMaxProcessAsympt<-rep(0,nbSample)
  mySmartMaxAsympt<-rep(0,nbSample)
  
  Psi0= 1
  
  
  for (i in 1:nbPaths){
    
    myProcessAsymptAtMarkers<-simDiscretizedOrnsteinInterferenceH1(nbTRN, dmark, nbMarkers, posMarkers, posMarkersAndQTL, indQTL, effectQTL, Psi0, errorVariance)  
    
   # myObservedProcessAsympt<-computeAsymptProcess(myProcessAsymptAtMarkers, dmark, dtest, nbMarkers, posMarkers)
   # myMaxProcessAsympt[i]<-max(myObservedProcessAsympt^2)
    
    mySmartMaxAsympt[i]<-computeAsymptProcessSmartFormula(myProcessAsymptAtMarkers, dmark, nbMarkers)
    
  }  
  
  
  # To plot the last path
  #plot(myObservedProcessAsympt)
 
  #sum(myMaxProcessAsympt>8.25)

  print("Asymptotic")  
  sum(mySmartMaxAsympt>6.76)
    
# Check asymptotics   
  

 

# My main code for the interference phenomenon under Student setting
 
source("myLibForStudent.R")
source("myLibForInterference.R")
  
######################### Input parameters #############################
 
set.seed(15)
 
sizeChro=1
# sizeChro: chromosome length in Morgan
#sizeChro=4
#sizeChro=7
 
nbQTLs=2
# nbQTLs: number of QTLs

dmark=0.2
#dmark=0.25
# dmark: distance between markers, we suppose equally spaced
 
nbTRN=1000
#nbTRN=100
# nbTRN: number of training individuals

nbMarkers= sizeChro/dmark +1
# nbMarkers: number of markers

# posQTL: vector of QTL locations 
# QTL locations must not match exactly marker locations
# Be careful under interference QTL must not belong to same marker intervals
posQTL=c(0.10001, 0.750011)
#posQTL=c(0.10001, 2.800011)

#effectQTL=c(0,0) 
effectQTL=c(0,4/sqrt(nbTRN))
# effectQTL: vector of QTL effects

mu=0
# mu is the global mean

errorVariance=1
# errorVariance is the environmental variance

thresholdAzais=6.76
#thresholdAzais=9.01
#thresholdAzais=10.12
# Threshold for QTL mapping, to be replaced by the value associated to the marker density
# It is the threshold according to Azais et al formula
 
# Number of degrees of freedom
df=4

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
  #print(i)
  
  # we begin by generating genomes under Haldane
  genomeTRN = generateGenomes(nbTRN,nbMarkersAndQTL,sizeChro)
  # now for QTL Location, we generate under interference
  newGenomeTRN=generateGenomesInterference(genomeTRN,nbTRN, posMarkersAndQTL,indQTL, nbQTLs)
  # genomeTRN : genome of the TRN individuals of size nbTRN times nbMarkersAndQTL

  yTRN =generateStudentPhenotypes(newGenomeTRN, nbTRN, mu, effectQTL, indQTL, errorVariance, df)
  # yTRN : vector of size nbTRN containing the phenotypes Y
   
  mySmartMax[i]<-computeTestStatIMSmartFormulaStudent(newGenomeTRN,nbTRN,sizeChro,dmark,nbMarkers,indMarker,posMarkers,yTRN,df)
  # compute the smart formula for the max of the process, according to Azais et al (2014)
  # i.e. tests only on markers and a signle test within marker intervals
  mySmartMaxNoPlugin[i]<-computeTestStatIMSmartFormulaStudentNoPlugin(newGenomeTRN,nbTRN,dmark,nbMarkers,yTRN,df,mu,errorVariance)
  # same as above but we do not use a plug in (ie we use true value of parameters)
  
} 


print("Number of rejections")

print("Azais formula")
print("Plugin")
sum(mySmartMax>thresholdAzais)
print("No Plugin")
sum(mySmartMaxNoPlugin>thresholdAzais)
 
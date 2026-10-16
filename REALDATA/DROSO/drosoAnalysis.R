# Analysis of Drosophila melanogaster data

source("myLibsForRealData.R")

# Analysis of real data of Drosophila melanogaster bristle

library(qtl)
library(ordinal)
library(stats)


data(bristleX)
summary(bristleX)
# cf. p34 and 35 of qtl manual

# Genetic map
bristleX$geno$X$map

length(bristleX$geno$X$map)
# 17 intervals

posMarkers =c(bristleX$geno$X$map)/100
# marker positions in Morgan

sizeChro=posMarkers[length(posMarkers)]-posMarkers[1]
 
dtest=0.00001


yTRN<-bristleX$pheno$Abfemale
nbTRN<-length(yTRN)
# 92


genomeTRN<-bristleX$geno$X$data
nbMarkers<-dim(genomeTRN)[2]
# 17

length(posMarkers)
# 17

# Removing marker 14 because too many missing data
genomeTRNFiltered<-genomeTRN[,-14]
posMarkersFiltered<-posMarkers[-14]

nbMarkersFiltered<-dim(genomeTRNFiltered)[2]
# 16
length(posMarkersFiltered)
# 16

################
genomeTRNNewCode<-genomeTRNFiltered
# we consider the code with +1 and -1

for (i in 1:dim(genomeTRNFiltered)[1]){
  for (j in 1:dim(genomeTRNFiltered)[2]){
    
    if (genomeTRNFiltered[i,j]==2){
      genomeTRNNewCode[i,j]<--1
    }
    
  } 
}


# To sum up,  we have only removed individuals with missing data 
# and we have kept all markers

#Last, we have used a new encoding to have alleles codes +1 and -1

########### Gaussian analysis #####################

# Compute the max thanks to the interpolation
myScoreGaussian<-computeTestStatIM(genomeTRNNewCode,yTRN,sizeChro,dtest,posMarkersFiltered,nbMarkersFiltered)

max(myScoreGaussian^2)
### 54.64473

# Compute the max thanks to Azais formula
mySmartMax<-computeTestStatIMSmartFormula(genomeTRNNewCode,yTRN,posMarkersFiltered,nbMarkersFiltered)
max(mySmartMax)

# As expected, the two maximum are the same !
  
#### Argmax now

myArgMax=which.max(myScoreGaussian^2)
seq(posMarkers[1],posMarkers[1]+sizeChro,by=dtest)[myArgMax]
# 0.06
# QTL estimated at  6 cM

############# Student below #############

# df is the number of degrees of freedom

###### Student with 4 degrees of freedom  #######
df=4

# using smart formula 
mySmartMaxStudentDf4<-computeTestStatIMSmartFormulaStudent(genomeTRNNewCode,nbTRN,nbMarkersFiltered,posMarkersFiltered,yTRN,df)
#[1] 77.9203 

# using interpolations 
myScoreStudentDf4<-computeTestStatStudentIM(genomeTRNNewCode,yTRN,sizeChro,dtest,posMarkersFiltered,nbMarkersFiltered,df)

max(myScoreStudentDf4^2)
#[1] 77.9203 

myArgMax=which.max(myScoreStudentDf4^2)
seq(posMarkers[1],posMarkers[1]+sizeChro,by=dtest)[myArgMax]
# [1] 0.06
 

###### Student with 20 degrees of freedom  #######
df=20

# using smart formula 
mySmartMaxStudentDf20<-computeTestStatIMSmartFormulaStudent(genomeTRNNewCode,nbTRN,nbMarkersFiltered,posMarkersFiltered,yTRN,df)
#[1] 60.38718

# using interpolations 
myScoreStudentDf20<-computeTestStatStudentIM(genomeTRNNewCode,yTRN,sizeChro,dtest,posMarkersFiltered,nbMarkersFiltered,df)

max(myScoreStudentDf20^2)
#[1] 60.38718

myArgMax=which.max(myScoreStudentDf20^2)
seq(posMarkers[1],posMarkers[1]+sizeChro,by=dtest)[myArgMax]
# [1] 0.06
 
###### Student with 60 degrees of freedom  #######
df=60

# using smart formula 
mySmartMaxStudentDf60<-computeTestStatIMSmartFormulaStudent(genomeTRNNewCode,nbTRN,nbMarkersFiltered,posMarkersFiltered,yTRN,df)
#[1] 56.21366

# using interpolations 
myScoreStudentDf60<-computeTestStatStudentIM(genomeTRNNewCode,yTRN,sizeChro,dtest,posMarkersFiltered,nbMarkersFiltered,df)

max(myScoreStudentDf60^2)
#[1] 56.21366

myArgMax=which.max(myScoreStudentDf60^2)
seq(posMarkers[1],posMarkers[1]+sizeChro,by=dtest)[myArgMax]
# [1] 0.06
 

################# Graphs #######################

pdf("DrosoProfile.pdf")
plot(seq(posMarkers[1],posMarkers[1]+sizeChro,by=dtest),myScoreGaussian^2,col="green",ylim=c(0,80), xlab="Location in Morgan", ylab="Square of Score statistic")

lines(seq(posMarkers[1],posMarkers[1]+sizeChro,by=dtest),myScoreStudentDf4^2,col="red")

lines(seq(posMarkers[1],posMarkers[1]+sizeChro,by=dtest),myScoreStudentDf20^2,col="orange")

lines(seq(posMarkers[1],posMarkers[1]+sizeChro,by=dtest),myScoreStudentDf60^2,col="black")

abline(h=8.64,col="pink",lty=2)

dev.off()

################# End Graphs #######################



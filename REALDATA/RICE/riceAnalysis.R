# Analysis of Huang data

source("myLibsForRealData.R")

# Genetic map
posMarkers=c(0.0,19.2,35.3,40.1,44.8,60.1,75.6,90.6,94.4,97.7,132.0,134.5,158.0,166.2,179.4,212.5,215.1,224.3)/100
# Marker locations in Morgan

nbMarkers=length(posMarkers)
sizeChro=posMarkers[length(posMarkers)]-posMarkers[1]
dtest=0.00001

myHuangTable<-read.table("HuangData", head=FALSE, sep="")
dim(myHuangTable)
# [1] 106  19

# Filtering for missing data
myHuangTableFiltered<-myHuangTable[!(apply(myHuangTable, 1, function(y) any(y == 0))),]
dim(myHuangTableFiltered)
# [1] 55 19

nbTRN<-dim(myHuangTableFiltered)[1]
## 55 TRN individuals
yTRN<-myHuangTableFiltered$V19 
# phenotype of interest

# Genomes 
genomeTRN<-myHuangTableFiltered[,1:18]
genomeTRNNewCode<-genomeTRN 
# we consider the code with +1 and -1

for (i in 1:dim(genomeTRN)[1]){
  for (j in 1:dim(genomeTRN)[2]){
    
    if (genomeTRN[i,j]==2){
      genomeTRNNewCode[i,j]<--1
    }
    
  } 
}


# To sum up,  we have only removed individuals with missing data 
# and we have kept all markers

#Last, we have used a new encoding to have alleles codes +1 and -1

########### Gaussian analysis #####################

# Compute the max thanks to the interpolation
myScoreGaussian<-computeTestStatIM(genomeTRNNewCode,yTRN,sizeChro,dtest,posMarkers,nbMarkers)
# myScoreGaussian[1:224301] 

max(myScoreGaussian^2)
### 23.58633 

# Compute the max thanks to Azais formula
mySmartMax<-computeTestStatIMSmartFormula(genomeTRNNewCode,yTRN,posMarkers,nbMarkers)
max(mySmartMax)
### 23.58633

# As expected, the two maximum are the same !

# Threshold 8.36 according Azais formula
#  So 23.58>8.36, we reject H0 at the 5% significance level
 
#### Argmax now

myArgMax=which.max(myScoreGaussian^2)
#[1] 197279
seq(posMarkers[1],posMarkers[1]+sizeChro,by=dtest)[myArgMax]
# 1.97278
# QTL estimated at  1.973 Morgan

# If I want to plot the interpolation
plot(seq(0,sizeChro,by=dtest),myScoreGaussian^2)
 
############# Student below #############

# df is the number of degrees of freedom

###### Student with 4 degrees of freedom  #######
df=4

# using smart formula 
mySmartMaxStudentDf4<-computeTestStatIMSmartFormulaStudent(genomeTRNNewCode,nbTRN,nbMarkers,posMarkers,yTRN,df)

# using interpolations 
myScoreStudentDf4<-computeTestStatStudentIM(genomeTRNNewCode,yTRN,sizeChro,dtest,posMarkers,nbMarkers,df)
max(myScoreStudentDf4^2)
# 30.57138

myArgMax=which.max(myScoreStudentDf4^2)
seq(posMarkers[1],posMarkers[1]+sizeChro,by=dtest)[myArgMax]
# [1] 1.9998

###### Student with 20 degrees of freedom  #######
df=20

# using smart formula 
mySmartMaxStudentDf20<-computeTestStatIMSmartFormulaStudent(genomeTRNNewCode,nbTRN,nbMarkers,posMarkers,yTRN,df)
# [1] 25.42155

# using interpolations 
myScoreStudentDf20<-computeTestStatStudentIM(genomeTRNNewCode,yTRN,sizeChro,dtest,posMarkers,nbMarkers,df)
max(myScoreStudentDf20^2)
# [1] 25.42155

myArgMax=which.max(myScoreStudentDf20^2)
seq(posMarkers[1],posMarkers[1]+sizeChro,by=dtest)[myArgMax]
# [1] 1.97728

df=60
# using smart formula 
mySmartMaxStudentDf60<-computeTestStatIMSmartFormulaStudent(genomeTRNNewCode,nbTRN,nbMarkers,posMarkers,yTRN,df)

# using interpolations 
myScoreStudentDf60<-computeTestStatStudentIM(genomeTRNNewCode,yTRN,sizeChro,dtest,posMarkers,nbMarkers,df)
max(myScoreStudentDf60^2)
# [1] 23.94233
 
myArgMax=which.max(myScoreStudentDf60^2)
seq(posMarkers[1],posMarkers[1]+sizeChro,by=dtest)[myArgMax]
# [1] 1.97426
 
################# Graphs #######################
 
pdf("HuangProfile.pdf")
plot(seq(posMarkers[1],posMarkers[1]+sizeChro,by=dtest),myScoreGaussian^2,col="green",ylim=c(0,35), xlab="Location in Morgan", ylab="Square of Score statistic")

lines(seq(posMarkers[1],posMarkers[1]+sizeChro,by=dtest),myScoreStudentDf4^2,col="red")

lines(seq(posMarkers[1],posMarkers[1]+sizeChro,by=dtest),myScoreStudentDf20^2,col="orange")

lines(seq(posMarkers[1],posMarkers[1]+sizeChro,by=dtest),myScoreStudentDf60^2,col="black")

abline(h=8.64,col="pink",lty=2)

dev.off()


################# End Graphs #######################

  

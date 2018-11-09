#Code for comparing across brain regions: Vikram's data
###Megan Hagenauer, Feb 11 2014


#Loading a package that specializes in multilevel modeling:
install.packages("lme4")
install.packages("Matrix")
install.packages("Rcpp")
library(Matrix)
library(Rcpp)
library(lme4)
library(nlme)


#Loading the standard packages
library(gdata)
library(fields)
library(stats)
library(car)
library(affy)
library(preprocessCore)
library(multtest)
library(All)
library(ConsensusClusterPlus)
library(plyr)
library(lattice)
library(datasets)
library(RColorBrewer)
library(grDevices)



#1: Quick and dirty: Looking for overlapping top genes by reading in the p-values and betas associated with each of the main effects in different regions.  Note: For Vikram's data, this may mean going back and running an identical model for each brain region. 


#1.A.1: Change the Working Directory to that for LM2 for each brain region. 
           #E.g, Pritzker960 Data Release ->RegionName->I am interested in diagnosis

#1.A.2: Read in the data and assign it a brain region identifier.


ReadInData<-function(filename, Region, Variable){
	
LM2testOutput<-as.matrix(read.table(filename, header=T, row.names=2, sep=","))
print("Dimensions:")
print(dim(LM2testOutput))
print("Colnames:")
print(colnames(LM2testOutput))
print("Example:")
print(LM2testOutput[c(1:3), c(1:5)])

#Double check that the probe name is the first column, unadjusted p-values is the second column, and the beta is the 5th column

TempPvalues<-as.matrix(as.numeric(LM2testOutput[,2]))
row.names(TempPvalues)<-row.names(LM2testOutput)
print("Pvalues are Numeric?")
print(is.numeric(TempPvalues))

TempBetas<-as.matrix(as.numeric(LM2testOutput[,5]))
row.names(TempBetas)<-row.names(LM2testOutput)
print("Betas are Numeric?")
print(is.numeric(TempBetas))

assign(paste(Region, Variable, "Pvalues", sep="_"), TempPvalues, envir=as.environment(1))
rm(TempPvalues)

assign(paste(Region, Variable, "Betas", sep="_"), TempBetas, envir=as.environment(1))
rm(TempBetas)
}



#Changed working directory to AAA

ReadInData("LM2testOutputAgev2.csv", "AAA", "Age")
ReadInData("LM2testOutputRNAIntegrityv2.csv", "AAA", "RNAIntegrity")
ReadInData("LM2testOutputGenderv2.csv", "AAA", "Gender")
ReadInData("LM2testOutputHoursFinalv2.csv", "AAA", "HoursFinal")
ReadInData("LM2testOutputMDDv2.csv", "AAA", "MDD")
ReadInData("LM2testOutputBrainPHv2.csv", "AAA", "PH")

AAA_PvalueOutputSummary<-read.csv("PvalueOutputSummary.csv", header=TRUE)


AAA_NormNoOutliers<-as.matrix(read.table("Quantile Normalized filtered data No Outliers.txt", header=T, row.names=1, sep="\t"))
dim(AAA_NormNoOutliers)

AAA_SampleInfo<-as.matrix(read.csv("SampleInfoNoOutliers.csv",  header=T))

cbind(colnames(AAA_NormNoOutliers), AAA_SampleInfo[,1])

NewColnamesForNormNoOutliers<- matrix(0, length(AAA_NormNoOutliers[1,]), 1)

for(i in 1:length(AAA_NormNoOutliers[1,])){
NewColnamesForNormNoOutliers[i]<-paste("AAA", AAA_SampleInfo[i,5], AAA_SampleInfo[i,1], sep="_")
}

colnames(AAA_NormNoOutliers)<-NewColnamesForNormNoOutliers




#Changed working directory to AB

ReadInData("LM2testOutputAgev2.csv", "AB", "Age")
ReadInData("LM2testOutputRNAIntegrityv2.csv", "AB", "RNAIntegrity")
ReadInData("LM2testOutputGenderv2.csv", "AB", "Gender")
ReadInData("LM2testOutputHoursFinalv2.csv", "AB", "HoursFinal")
ReadInData("LM2testOutputMDDv2.csv", "AB", "MDD")
ReadInData("LM2testOutputBrainPHv2.csv", "AB", "PH")

AB_PvalueOutputSummary<-read.csv("PvalueOutputSummary.csv", header=TRUE)


AB_NormNoOutliers<-as.matrix(read.table("Quantile Normalized filtered data No Outliers.txt", header=T, row.names=1, sep="\t"))
dim(AB_NormNoOutliers)

AB_SampleInfo<-as.matrix(read.csv("SampleInfoNoOutliers.csv",  header=T))

cbind(colnames(AB_NormNoOutliers), AB_SampleInfo[,1])

NewColnamesForNormNoOutliers<- matrix(0, length(AB_NormNoOutliers[1,]), 1)

for(i in 1:length(AB_NormNoOutliers[1,])){
NewColnamesForNormNoOutliers[i]<-paste("AB", AB_SampleInfo[i,5], AB_SampleInfo[i,1], sep="_")
}

colnames(AB_NormNoOutliers)<-NewColnamesForNormNoOutliers




#Changed working directory to AHA

ReadInData("LM2testOutputAgev2.csv", "AHA", "Age")
ReadInData("LM2testOutputRNAIntegrityv2.csv", "AHA", "RNAIntegrity")
ReadInData("LM2testOutputGenderv2.csv", "AHA", "Gender")
ReadInData("LM2testOutputHoursFinalv2.csv", "AHA", "HoursFinal")
ReadInData("LM2testOutputMDDv2.csv", "AHA", "MDD")
ReadInData("LM2testOutputBrainPHv2.csv", "AHA", "PH")

AHA_PvalueOutputSummary<-read.csv("PvalueOutputSummary.csv", header=TRUE)


AHA_NormNoOutliers<-as.matrix(read.table("Quantile Normalized filtered data No Outliers.txt", header=T, row.names=1, sep="\t"))
dim(AHA_NormNoOutliers)

AHA_SampleInfo<-as.matrix(read.csv("SampleInfoNoOutliers.csv",  header=T))

cbind(colnames(AHA_NormNoOutliers), AHA_SampleInfo[,1])

NewColnamesForNormNoOutliers<- matrix(0, length(AHA_NormNoOutliers[1,]), 1)

for(i in 1:length(AHA_NormNoOutliers[1,])){
NewColnamesForNormNoOutliers[i]<-paste("AHA", AHA_SampleInfo[i,5], AHA_SampleInfo[i,1], sep="_")
}

colnames(AHA_NormNoOutliers)<-NewColnamesForNormNoOutliers




#Changed working directory to Basal

ReadInData("LM2testOutputAgev2.csv", "Basal", "Age")
ReadInData("LM2testOutputRNAIntegrityv2.csv", "Basal", "RNAIntegrity")
ReadInData("LM2testOutputGenderv2.csv", "Basal", "Gender")
ReadInData("LM2testOutputHoursFinalv2.csv", "Basal", "HoursFinal")
ReadInData("LM2testOutputMDDv2.csv", "Basal", "MDD")
ReadInData("LM2testOutputBrainPHv2.csv", "Basal", "PH")

Basal_PvalueOutputSummary<-read.csv("PvalueOutputSummary.csv", header=TRUE)



Basal_NormNoOutliers<-as.matrix(read.table("Quantile Normalized filtered data No Outliers.txt", header=T, row.names=1, sep="\t"))
dim(Basal_NormNoOutliers)

Basal_SampleInfo<-as.matrix(read.csv("SampleInfoNoOutliers.csv",  header=T))

cbind(colnames(Basal_NormNoOutliers), Basal_SampleInfo[,1])

NewColnamesForNormNoOutliers<- matrix(0, length(Basal_NormNoOutliers[1,]), 1)

for(i in 1:length(Basal_NormNoOutliers[1,])){
NewColnamesForNormNoOutliers[i]<-paste("Basal", Basal_SampleInfo[i,5], Basal_SampleInfo[i,1], sep="_")
}

colnames(Basal_NormNoOutliers)<-NewColnamesForNormNoOutliers




#Changed working directory to Central

ReadInData("LM2testOutputAgev2.csv", "Central", "Age")
ReadInData("LM2testOutputRNAIntegrityv2.csv", "Central", "RNAIntegrity")
ReadInData("LM2testOutputGenderv2.csv", "Central", "Gender")
ReadInData("LM2testOutputHoursFinalv2.csv", "Central", "HoursFinal")
ReadInData("LM2testOutputMDDv2.csv", "Central", "MDD")
ReadInData("LM2testOutputBrainPHv2.csv", "Central", "PH")

Central_PvalueOutputSummary<-read.csv("PvalueOutputSummary.csv", header=TRUE)


Central_NormNoOutliers<-as.matrix(read.table("Quantile Normalized filtered data No Outliers.txt", header=T, row.names=1, sep="\t"))
dim(Central_NormNoOutliers)

Central_SampleInfo<-as.matrix(read.csv("SampleInfoNoOutliers.csv",  header=T))

cbind(colnames(Central_NormNoOutliers), Central_SampleInfo[,1])

NewColnamesForNormNoOutliers<- matrix(0, length(Central_NormNoOutliers[1,]), 1)

for(i in 1:length(Central_NormNoOutliers[1,])){
NewColnamesForNormNoOutliers[i]<-paste("Central", Central_SampleInfo[i,5], Central_SampleInfo[i,1], sep="_")
}

colnames(Central_NormNoOutliers)<-NewColnamesForNormNoOutliers





#Changed working directory to CO

ReadInData("LM2testOutputAgev2.csv", "CO", "Age")
ReadInData("LM2testOutputRNAIntegrityv2.csv", "CO", "RNAIntegrity")
ReadInData("LM2testOutputGenderv2.csv", "CO", "Gender")
ReadInData("LM2testOutputHoursFinalv2.csv", "CO", "HoursFinal")
ReadInData("LM2testOutputMDDv2.csv", "CO", "MDD")
ReadInData("LM2testOutputBrainPHv2.csv", "CO", "PH")

CO_PvalueOutputSummary<-read.csv("PvalueOutputSummary.csv", header=TRUE)


CO_NormNoOutliers<-as.matrix(read.table("Quantile Normalized filtered data No Outliers.txt", header=T, row.names=1, sep="\t"))
dim(CO_NormNoOutliers)

CO_SampleInfo<-as.matrix(read.csv("SampleInfoNoOutliers.csv",  header=T))

cbind(colnames(CO_NormNoOutliers), CO_SampleInfo[,1])

NewColnamesForNormNoOutliers<- matrix(0, length(CO_NormNoOutliers[1,]), 1)

for(i in 1:length(CO_NormNoOutliers[1,])){
NewColnamesForNormNoOutliers[i]<-paste("CO", CO_SampleInfo[i,5], CO_SampleInfo[i,1], sep="_")
}

colnames(CO_NormNoOutliers)<-NewColnamesForNormNoOutliers




#Changed working directory to Lateral

ReadInData("LM2testOutputAgev2.csv", "Lateral", "Age")
ReadInData("LM2testOutputRNAIntegrityv2.csv", "Lateral", "RNAIntegrity")
ReadInData("LM2testOutputGenderv2.csv", "Lateral", "Gender")
ReadInData("LM2testOutputHoursFinalv2.csv", "Lateral", "HoursFinal")
ReadInData("LM2testOutputMDDv2.csv", "Lateral", "MDD")
ReadInData("LM2testOutputBrainPHv2.csv", "Lateral", "PH")

Lateral_PvalueOutputSummary<-read.csv("PvalueOutputSummary.csv", header=TRUE)


Lateral_NormNoOutliers<-as.matrix(read.table("Quantile Normalized filtered data No Outliers.txt", header=T, row.names=1, sep="\t"))
dim(Lateral_NormNoOutliers)

Lateral_SampleInfo<-as.matrix(read.csv("SampleInfoNoOutliers.csv",  header=T))

cbind(colnames(Lateral_NormNoOutliers), Lateral_SampleInfo[,1])

NewColnamesForNormNoOutliers<- matrix(0, length(Lateral_NormNoOutliers[1,]), 1)

for(i in 1:length(Lateral_NormNoOutliers[1,])){
NewColnamesForNormNoOutliers[i]<-paste("Lateral", Lateral_SampleInfo[i,5], Lateral_SampleInfo[i,1], sep="_")
}

colnames(Lateral_NormNoOutliers)<-NewColnamesForNormNoOutliers






#Changed working directory to Medial

ReadInData("LM2testOutputAgev2.csv", "Medial", "Age")
ReadInData("LM2testOutputRNAIntegrityv2.csv", "Medial", "RNAIntegrity")
ReadInData("LM2testOutputGenderv2.csv", "Medial", "Gender")
ReadInData("LM2testOutputHoursFinalv2.csv", "Medial", "HoursFinal")
ReadInData("LM2testOutputMDDv2.csv", "Medial", "MDD")
ReadInData("LM2testOutputBrainPHv2.csv", "Medial", "PH")

Medial_PvalueOutputSummary<-read.csv("PvalueOutputSummary.csv", header=TRUE)

Medial_NormNoOutliers<-as.matrix(read.table("Quantile Normalized filtered data No Outliers.txt", header=T, row.names=1, sep="\t"))
dim(Medial_NormNoOutliers)

Medial_SampleInfo<-as.matrix(read.csv("SampleInfoNoOutliers.csv",  header=T))

cbind(colnames(Medial_NormNoOutliers), Medial_SampleInfo[,1])

NewColnamesForNormNoOutliers<- matrix(0, length(Medial_NormNoOutliers[1,]), 1)

for(i in 1:length(Medial_NormNoOutliers[1,])){
NewColnamesForNormNoOutliers[i]<-paste("Medial", Medial_SampleInfo[i,5], Medial_SampleInfo[i,1], sep="_")
}

colnames(Medial_NormNoOutliers)<-NewColnamesForNormNoOutliers







#Changed working directory to PAC

ReadInData("LM2testOutputAgev2.csv", "PAC", "Age")
ReadInData("LM2testOutputRNAIntegrityv2.csv", "PAC", "RNAIntegrity")
ReadInData("LM2testOutputGenderv2.csv", "PAC", "Gender")
ReadInData("LM2testOutputHoursFinalv2.csv", "PAC", "HoursFinal")
ReadInData("LM2testOutputMDDv2.csv", "PAC", "MDD")
ReadInData("LM2testOutputBrainPHv2.csv", "PAC", "PH")

PAC_PvalueOutputSummary<-read.csv("PvalueOutputSummary.csv", header=TRUE)


PAC_NormNoOutliers<-as.matrix(read.table("Quantile Normalized filtered data No Outliers.txt", header=T, row.names=1, sep="\t"))
dim(PAC_NormNoOutliers)

PAC_SampleInfo<-as.matrix(read.csv("SampleInfoNoOutliers.csv",  header=T))

cbind(colnames(PAC_NormNoOutliers), PAC_SampleInfo[,1])

NewColnamesForNormNoOutliers<- matrix(0, length(PAC_NormNoOutliers[1,]), 1)

for(i in 1:length(PAC_NormNoOutliers[1,])){
NewColnamesForNormNoOutliers[i]<-paste("PAC", PAC_SampleInfo[i,5], PAC_SampleInfo[i,1], sep="_")
}

colnames(PAC_NormNoOutliers)<-NewColnamesForNormNoOutliers






#Changed working directory to PL

ReadInData("LM2testOutputAgev2.csv", "PL", "Age")
ReadInData("LM2testOutputRNAIntegrityv2.csv", "PL", "RNAIntegrity")
ReadInData("LM2testOutputGenderv2.csv", "PL", "Gender")
ReadInData("LM2testOutputHoursFinalv2.csv", "PL", "HoursFinal")
ReadInData("LM2testOutputMDDv2.csv", "PL", "MDD")
ReadInData("LM2testOutputBrainPHv2.csv", "PL", "PH")

PL_PvalueOutputSummary<-read.csv("PvalueOutputSummary.csv", header=TRUE)

PL_NormNoOutliers<-as.matrix(read.table("Quantile Normalized filtered data No Outliers.txt", header=T, row.names=1, sep="\t"))
dim(PL_NormNoOutliers)

PL_SampleInfo<-as.matrix(read.csv("SampleInfoNoOutliers.csv",  header=T))

cbind(colnames(PL_NormNoOutliers), PL_SampleInfo[,1])

NewColnamesForNormNoOutliers<- matrix(0, length(PL_NormNoOutliers[1,]), 1)

for(i in 1:length(PL_NormNoOutliers[1,])){
NewColnamesForNormNoOutliers[i]<-paste("PL", PL_SampleInfo[i,5], PL_SampleInfo[i,1], sep="_")
}

colnames(PL_NormNoOutliers)<-NewColnamesForNormNoOutliers










###B: Determine which probes are expressed in *all 10 brain regions* (this is true for 960 and Vikram's data) - I believe this would be a combination of the find statement and a really complex AND statement

AAA_IndicesForOverlappingGenes<- matrix(c(

row.names(AAA_MDD_Pvalues)%in%row.names(AB_MDD_Pvalues),
row.names(AAA_MDD_Pvalues)%in%row.names(AHA_MDD_Pvalues),
row.names(AAA_MDD_Pvalues)%in%row.names(Basal_MDD_Pvalues),
row.names(AAA_MDD_Pvalues)%in%row.names(Central_MDD_Pvalues),
row.names(AAA_MDD_Pvalues)%in%row.names(CO_MDD_Pvalues),
row.names(AAA_MDD_Pvalues)%in%row.names(Lateral_MDD_Pvalues),
row.names(AAA_MDD_Pvalues)%in%row.names(Medial_MDD_Pvalues),
row.names(AAA_MDD_Pvalues)%in%row.names(PAC_MDD_Pvalues),
row.names(AAA_MDD_Pvalues)%in%row.names(PL_MDD_Pvalues)
), 
nrow=length(AAA_MDD_Pvalues), ncol=9)

AAA_IndicesForOverlappingGenes[c(1:3),]

AAA_IndicesForOverlappingGenes_All<-matrix(T, length(AAA_MDD_Pvalues[,1]), 1)

for(i in 1:length(AAA_MDD_Pvalues[,1])){	
AAA_IndicesForOverlappingGenes_All[i,1]<-sum(AAA_IndicesForOverlappingGenes[i,]==T)==9
}

sum(AAA_IndicesForOverlappingGenes_All[,1]==T)
#[1]  15788
#Out of 18863 original probes?  That's crazy!  




AAA_Age_Pvalues_OverlappingGenes<-AAA_Age_Pvalues[(AAA_IndicesForOverlappingGenes_All[,1]==T),]
AAA_RNAIntegrity_Pvalues_OverlappingGenes<-AAA_RNAIntegrity_Pvalues[(AAA_IndicesForOverlappingGenes_All[,1]==T),]
AAA_Gender_Pvalues_OverlappingGenes<-AAA_Gender_Pvalues[(AAA_IndicesForOverlappingGenes_All[,1]==T),]
AAA_HoursFinal_Pvalues_OverlappingGenes<-AAA_HoursFinal_Pvalues[(AAA_IndicesForOverlappingGenes_All[,1]==T),]
AAA_MDD_Pvalues_OverlappingGenes<-AAA_MDD_Pvalues[(AAA_IndicesForOverlappingGenes_All[,1]==T),]
AAA_PH_Pvalues_OverlappingGenes<-AAA_PH_Pvalues[(AAA_IndicesForOverlappingGenes_All[,1]==T),]

AAA_NormNoOutliers_OverlappingGenes<-AAA_NormNoOutliers[(AAA_IndicesForOverlappingGenes_All[,1]==T),]


length(AAA_MDD_Pvalues_OverlappingGenes)
AAA_MDD_Pvalues_OverlappingGenes[c(1:5)]



#Ideally I should really make this code neater so that AAA is included in the Overlapping genes matrix, so that I can line the code up between regions and compare the total number of overlapping genes across regions.

AAA_IndicesForOverlappingGenes2<- matrix(c(
row.names(AAA_MDD_Pvalues)%in%row.names(AAA_MDD_Pvalues),
row.names(AAA_MDD_Pvalues)%in%row.names(AB_MDD_Pvalues),
row.names(AAA_MDD_Pvalues)%in%row.names(AHA_MDD_Pvalues),
row.names(AAA_MDD_Pvalues)%in%row.names(Basal_MDD_Pvalues),
row.names(AAA_MDD_Pvalues)%in%row.names(Central_MDD_Pvalues),
row.names(AAA_MDD_Pvalues)%in%row.names(CO_MDD_Pvalues),
row.names(AAA_MDD_Pvalues)%in%row.names(Lateral_MDD_Pvalues),
row.names(AAA_MDD_Pvalues)%in%row.names(Medial_MDD_Pvalues),
row.names(AAA_MDD_Pvalues)%in%row.names(PAC_MDD_Pvalues),
row.names(AAA_MDD_Pvalues)%in%row.names(PL_MDD_Pvalues)
), 
nrow=length(AAA_MDD_Pvalues), ncol=10)


AAA_IndicesForOverlappingGenes2[c(1:3),]

AAA_IndicesForOverlappingGenes_Any<-matrix(T, length(AAA_MDD_Pvalues[,1]), 1)

for(i in 1:length(AAA_MDD_Pvalues[,1])){	
AAA_IndicesForOverlappingGenes_Any[i,1]<-sum(AAA_IndicesForOverlappingGenes2[i,]==T)>1
}

sum(AAA_IndicesForOverlappingGenes_Any[,1]==T)
[1] 18774


TotalOverlappingAAAwOther<-apply(AAA_IndicesForOverlappingGenes2, 2, sum)
names(TotalOverlappingAAAwOther)<-c("AAA", "AB", "AHA", "Basal", "Central", "CO", "Lateral", "Medial", "PAC", "PL")

sum((AAA_IndicesForOverlappingGenes_Any[,1]==T)&(AAA_IndicesForOverlappingGenes_All[,1]==F))
[1] 2986

18774-15788
[1] 2986

AAA_NormNoOutliers_OverlappingGenesNotAll<-AAA_NormNoOutliers[(AAA_IndicesForOverlappingGenes_Any[,1]==T)&(AAA_IndicesForOverlappingGenes_All[,1]==F),]


boxplot(data.frame(AAA_NormNoOutliers_OverlappingGenesNotAll), cex=0.25, las=3, par(cex.axis=0.75))
#Very few of these genes look like they are not noise.
#The entire box for their signal distribution tends to fall below 7.5

#What if we made things slightly more stringent:

AAA_IndicesForOverlappingGenes_Atleast3<-matrix(T, length(AAA_MDD_Pvalues[,1]), 1)

for(i in 1:length(AAA_MDD_Pvalues[,1])){	
AAA_IndicesForOverlappingGenes_Atleast3[i,1]<-sum(AAA_IndicesForOverlappingGenes2[i,]==T)>2
}

sum(AAA_IndicesForOverlappingGenes_Atleast3[,1]==T)
[1] 18660
18660-15788
[1] 2872

AAA_NormNoOutliers_OverlappingGenes_Atleast3<-AAA_NormNoOutliers[(AAA_IndicesForOverlappingGenes_Atleast3[,1]==T)&(AAA_IndicesForOverlappingGenes_All[,1]==F),]

boxplot(data.frame(AAA_NormNoOutliers_OverlappingGenes_Atleast3), cex=0.25, las=3, par(cex.axis=0.75))
#Looks about the same. So maybe these probes are just noise, perhaps they are just really low-level expression genes.





AB_IndicesForOverlappingGenes<- matrix(c(

row.names(AB_MDD_Pvalues)%in%row.names(AAA_MDD_Pvalues),
row.names(AB_MDD_Pvalues)%in%row.names(AHA_MDD_Pvalues),
row.names(AB_MDD_Pvalues)%in%row.names(Basal_MDD_Pvalues),
row.names(AB_MDD_Pvalues)%in%row.names(Central_MDD_Pvalues),
row.names(AB_MDD_Pvalues)%in%row.names(CO_MDD_Pvalues),
row.names(AB_MDD_Pvalues)%in%row.names(Lateral_MDD_Pvalues),
row.names(AB_MDD_Pvalues)%in%row.names(Medial_MDD_Pvalues),
row.names(AB_MDD_Pvalues)%in%row.names(PAC_MDD_Pvalues),
row.names(AB_MDD_Pvalues)%in%row.names(PL_MDD_Pvalues)
), 
nrow=length(AB_MDD_Pvalues), ncol=9)

AB_IndicesForOverlappingGenes[c(1:3),]

AB_IndicesForOverlappingGenes_All<-matrix(T, length(AB_MDD_Pvalues[,1]), 1)

for(i in 1:length(AB_MDD_Pvalues[,1])){	
AB_IndicesForOverlappingGenes_All[i,1]<-sum(AB_IndicesForOverlappingGenes[i,]==T)==9
}

sum(AB_IndicesForOverlappingGenes_All[,1]==T)
#[1] 15788
#Out of 19250 original probes.  Interesting.


AB_Age_Pvalues_OverlappingGenes<-AB_Age_Pvalues[(AB_IndicesForOverlappingGenes_All[,1]==T),]
AB_RNAIntegrity_Pvalues_OverlappingGenes<-AB_RNAIntegrity_Pvalues[(AB_IndicesForOverlappingGenes_All[,1]==T),]
AB_Gender_Pvalues_OverlappingGenes<-AB_Gender_Pvalues[(AB_IndicesForOverlappingGenes_All[,1]==T),]
AB_HoursFinal_Pvalues_OverlappingGenes<-AB_HoursFinal_Pvalues[(AB_IndicesForOverlappingGenes_All[,1]==T),]
AB_MDD_Pvalues_OverlappingGenes<-AB_MDD_Pvalues[(AB_IndicesForOverlappingGenes_All[,1]==T),]
AB_PH_Pvalues_OverlappingGenes<-AB_PH_Pvalues[(AB_IndicesForOverlappingGenes_All[,1]==T),]

AB_NormNoOutliers_OverlappingGenes<-AB_NormNoOutliers[(AB_IndicesForOverlappingGenes_All[,1]==T),]


length(AB_MDD_Pvalues_OverlappingGenes)
AB_MDD_Pvalues_OverlappingGenes[c(1:5)]



#Ideally I should really make this code neater so that AAA is included in the Overlapping genes matrix, so that I can line the code up between regions and compare the total number of overlapping genes across regions.

AB_IndicesForOverlappingGenes2<- matrix(c(
row.names(AB_MDD_Pvalues)%in%row.names(AAA_MDD_Pvalues),
row.names(AB_MDD_Pvalues)%in%row.names(AB_MDD_Pvalues),
row.names(AB_MDD_Pvalues)%in%row.names(AHA_MDD_Pvalues),
row.names(AB_MDD_Pvalues)%in%row.names(Basal_MDD_Pvalues),
row.names(AB_MDD_Pvalues)%in%row.names(Central_MDD_Pvalues),
row.names(AB_MDD_Pvalues)%in%row.names(CO_MDD_Pvalues),
row.names(AB_MDD_Pvalues)%in%row.names(Lateral_MDD_Pvalues),
row.names(AB_MDD_Pvalues)%in%row.names(Medial_MDD_Pvalues),
row.names(AB_MDD_Pvalues)%in%row.names(PAC_MDD_Pvalues),
row.names(AB_MDD_Pvalues)%in%row.names(PL_MDD_Pvalues)
), 
nrow=length(AB_MDD_Pvalues), ncol=10)


AB_IndicesForOverlappingGenes2[c(1:3),]

AB_IndicesForOverlappingGenes_Any<-matrix(T, length(AB_MDD_Pvalues[,1]), 1)

for(i in 1:length(AB_MDD_Pvalues[,1])){	
AB_IndicesForOverlappingGenes_Any[i,1]<-sum(AB_IndicesForOverlappingGenes2[i,]==T)>1
}

sum(AB_IndicesForOverlappingGenes_Any[,1]==T)
[1] 19201
19201-15788
[1] 3413

TotalOverlappingABwOther<-apply(AB_IndicesForOverlappingGenes2, 2, sum)
names(TotalOverlappingABwOther)<-c("AAA", "AB", "AHA", "Basal", "Central", "CO", "Lateral", "Medial", "PAC", "PL")

sum((AB_IndicesForOverlappingGenes_Any[,1]==T)&(AB_IndicesForOverlappingGenes_All[,1]==F))
[1] 3413


AB_NormNoOutliers_OverlappingGenesNotAll<-AB_NormNoOutliers[(AB_IndicesForOverlappingGenes_Any[,1]==T)&(AB_IndicesForOverlappingGenes_All[,1]==F),]


boxplot(data.frame(AB_NormNoOutliers_OverlappingGenesNotAll), cex=0.25, las=3, par(cex.axis=0.75))







AHA_IndicesForOverlappingGenes<- matrix(c(

row.names(AHA_MDD_Pvalues)%in%row.names(AAA_MDD_Pvalues),
row.names(AHA_MDD_Pvalues)%in%row.names(AB_MDD_Pvalues),
row.names(AHA_MDD_Pvalues)%in%row.names(Basal_MDD_Pvalues),
row.names(AHA_MDD_Pvalues)%in%row.names(Central_MDD_Pvalues),
row.names(AHA_MDD_Pvalues)%in%row.names(CO_MDD_Pvalues),
row.names(AHA_MDD_Pvalues)%in%row.names(Lateral_MDD_Pvalues),
row.names(AHA_MDD_Pvalues)%in%row.names(Medial_MDD_Pvalues),
row.names(AHA_MDD_Pvalues)%in%row.names(PAC_MDD_Pvalues),
row.names(AHA_MDD_Pvalues)%in%row.names(PL_MDD_Pvalues)
), 
nrow=length(AHA_MDD_Pvalues), ncol=9)

AHA_IndicesForOverlappingGenes[c(1:3),]

AHA_IndicesForOverlappingGenes_All<-matrix(T, length(AHA_MDD_Pvalues[,1]), 1)

for(i in 1:length(AHA_MDD_Pvalues[,1])){	
AHA_IndicesForOverlappingGenes_All[i,1]<-sum(AHA_IndicesForOverlappingGenes[i,]==T)==9
}

sum(AHA_IndicesForOverlappingGenes_All[,1]==T)
#[1] 15788
#Out of 20593 original probes


AHA_Age_Pvalues_OverlappingGenes<-AHA_Age_Pvalues[(AHA_IndicesForOverlappingGenes_All[,1]==T),]
AHA_RNAIntegrity_Pvalues_OverlappingGenes<-AHA_RNAIntegrity_Pvalues[(AHA_IndicesForOverlappingGenes_All[,1]==T),]
AHA_Gender_Pvalues_OverlappingGenes<-AHA_Gender_Pvalues[(AHA_IndicesForOverlappingGenes_All[,1]==T),]
AHA_HoursFinal_Pvalues_OverlappingGenes<-AHA_HoursFinal_Pvalues[(AHA_IndicesForOverlappingGenes_All[,1]==T),]
AHA_MDD_Pvalues_OverlappingGenes<-AHA_MDD_Pvalues[(AHA_IndicesForOverlappingGenes_All[,1]==T),]
AHA_PH_Pvalues_OverlappingGenes<-AHA_PH_Pvalues[(AHA_IndicesForOverlappingGenes_All[,1]==T),]


AHA_NormNoOutliers_OverlappingGenes<-AHA_NormNoOutliers[(AHA_IndicesForOverlappingGenes_All[,1]==T),]



length(AHA_MDD_Pvalues_OverlappingGenes)
AHA_MDD_Pvalues_OverlappingGenes[c(1:5)]




AHA_IndicesForOverlappingGenes2<- matrix(c(
row.names(AHA_MDD_Pvalues)%in%row.names(AAA_MDD_Pvalues),
row.names(AHA_MDD_Pvalues)%in%row.names(AB_MDD_Pvalues),
row.names(AHA_MDD_Pvalues)%in%row.names(AHA_MDD_Pvalues),
row.names(AHA_MDD_Pvalues)%in%row.names(Basal_MDD_Pvalues),
row.names(AHA_MDD_Pvalues)%in%row.names(Central_MDD_Pvalues),
row.names(AHA_MDD_Pvalues)%in%row.names(CO_MDD_Pvalues),
row.names(AHA_MDD_Pvalues)%in%row.names(Lateral_MDD_Pvalues),
row.names(AHA_MDD_Pvalues)%in%row.names(Medial_MDD_Pvalues),
row.names(AHA_MDD_Pvalues)%in%row.names(PAC_MDD_Pvalues),
row.names(AHA_MDD_Pvalues)%in%row.names(PL_MDD_Pvalues)
), 
nrow=length(AHA_MDD_Pvalues), ncol=10)


AHA_IndicesForOverlappingGenes2[c(1:3),]

AHA_IndicesForOverlappingGenes_Any<-matrix(T, length(AHA_MDD_Pvalues[,1]), 1)

for(i in 1:length(AHA_MDD_Pvalues[,1])){	
AHA_IndicesForOverlappingGenes_Any[i,1]<-sum(AHA_IndicesForOverlappingGenes2[i,]==T)>1
}

sum(AHA_IndicesForOverlappingGenes_Any[,1]==T)
[1] 19788
19788-15788
[1] 4000

TotalOverlappingAHAwOther<-apply(AHA_IndicesForOverlappingGenes2, 2, sum)
names(TotalOverlappingAHAwOther)<-c("AAA", "AB", "AHA", "Basal", "Central", "CO", "Lateral", "Medial", "PAC", "PL")

sum((AHA_IndicesForOverlappingGenes_Any[,1]==T)&(AHA_IndicesForOverlappingGenes_All[,1]==F))
[1] 4000


AHA_NormNoOutliers_OverlappingGenesNotAll<-AHA_NormNoOutliers[(AHA_IndicesForOverlappingGenes_Any[,1]==T)&(AHA_IndicesForOverlappingGenes_All[,1]==F),]


boxplot(data.frame(AHA_NormNoOutliers_OverlappingGenesNotAll), cex=0.25, las=3, par(cex.axis=0.75))
#Wow - the signal for those is really low - the boxes don't even reach 7.2 in general









Basal_IndicesForOverlappingGenes<- matrix(c(

row.names(Basal_MDD_Pvalues)%in%row.names(AAA_MDD_Pvalues),
row.names(Basal_MDD_Pvalues)%in%row.names(AB_MDD_Pvalues),
row.names(Basal_MDD_Pvalues)%in%row.names(AHA_MDD_Pvalues),
row.names(Basal_MDD_Pvalues)%in%row.names(Central_MDD_Pvalues),
row.names(Basal_MDD_Pvalues)%in%row.names(CO_MDD_Pvalues),
row.names(Basal_MDD_Pvalues)%in%row.names(Lateral_MDD_Pvalues),
row.names(Basal_MDD_Pvalues)%in%row.names(Medial_MDD_Pvalues),
row.names(Basal_MDD_Pvalues)%in%row.names(PAC_MDD_Pvalues),
row.names(Basal_MDD_Pvalues)%in%row.names(PL_MDD_Pvalues)
), 
nrow=length(Basal_MDD_Pvalues), ncol=9)

Basal_IndicesForOverlappingGenes[c(1:3),]

Basal_IndicesForOverlappingGenes_All<-matrix(T, length(Basal_MDD_Pvalues[,1]), 1)

for(i in 1:length(Basal_MDD_Pvalues[,1])){	
Basal_IndicesForOverlappingGenes_All[i,1]<-sum(Basal_IndicesForOverlappingGenes[i,]==T)==9
}

sum(Basal_IndicesForOverlappingGenes_All[,1]==T)
#[1] 15788
#Out of 18731 original probes


Basal_Age_Pvalues_OverlappingGenes<-Basal_Age_Pvalues[(Basal_IndicesForOverlappingGenes_All[,1]==T),]
Basal_RNAIntegrity_Pvalues_OverlappingGenes<-Basal_RNAIntegrity_Pvalues[(Basal_IndicesForOverlappingGenes_All[,1]==T),]
Basal_Gender_Pvalues_OverlappingGenes<-Basal_Gender_Pvalues[(Basal_IndicesForOverlappingGenes_All[,1]==T),]
Basal_HoursFinal_Pvalues_OverlappingGenes<-Basal_HoursFinal_Pvalues[(Basal_IndicesForOverlappingGenes_All[,1]==T),]
Basal_MDD_Pvalues_OverlappingGenes<-Basal_MDD_Pvalues[(Basal_IndicesForOverlappingGenes_All[,1]==T),]
Basal_PH_Pvalues_OverlappingGenes<-Basal_PH_Pvalues[(Basal_IndicesForOverlappingGenes_All[,1]==T),]


Basal_NormNoOutliers_OverlappingGenes<-Basal_NormNoOutliers[(Basal_IndicesForOverlappingGenes_All[,1]==T),]


length(Basal_MDD_Pvalues_OverlappingGenes)
Basal_MDD_Pvalues_OverlappingGenes[c(1:5)]






Basal_IndicesForOverlappingGenes2<- matrix(c(
row.names(Basal_MDD_Pvalues)%in%row.names(AAA_MDD_Pvalues),
row.names(Basal_MDD_Pvalues)%in%row.names(AB_MDD_Pvalues),
row.names(Basal_MDD_Pvalues)%in%row.names(AHA_MDD_Pvalues),
row.names(Basal_MDD_Pvalues)%in%row.names(Basal_MDD_Pvalues),
row.names(Basal_MDD_Pvalues)%in%row.names(Central_MDD_Pvalues),
row.names(Basal_MDD_Pvalues)%in%row.names(CO_MDD_Pvalues),
row.names(Basal_MDD_Pvalues)%in%row.names(Lateral_MDD_Pvalues),
row.names(Basal_MDD_Pvalues)%in%row.names(Medial_MDD_Pvalues),
row.names(Basal_MDD_Pvalues)%in%row.names(PAC_MDD_Pvalues),
row.names(Basal_MDD_Pvalues)%in%row.names(PL_MDD_Pvalues)
), 
nrow=length(Basal_MDD_Pvalues), ncol=10)


Basal_IndicesForOverlappingGenes2[c(1:3),]

Basal_IndicesForOverlappingGenes_Any<-matrix(T, length(Basal_MDD_Pvalues[,1]), 1)

for(i in 1:length(Basal_MDD_Pvalues[,1])){	
Basal_IndicesForOverlappingGenes_Any[i,1]<-sum(Basal_IndicesForOverlappingGenes2[i,]==T)>1
}

sum(Basal_IndicesForOverlappingGenes_Any[,1]==T)
[1] 18714
18714-15788
[1] 2926

TotalOverlappingBasalwOther<-apply(Basal_IndicesForOverlappingGenes2, 2, sum)
names(TotalOverlappingBasalwOther)<-c("AAA", "AB", "AHA", "Basal", "Central", "CO", "Lateral", "Medial", "PAC", "PL")

sum((Basal_IndicesForOverlappingGenes_Any[,1]==T)&(Basal_IndicesForOverlappingGenes_All[,1]==F))
[1] 2926


Basal_NormNoOutliers_OverlappingGenesNotAll<-Basal_NormNoOutliers[(Basal_IndicesForOverlappingGenes_Any[,1]==T)&(Basal_IndicesForOverlappingGenes_All[,1]==F),]


boxplot(data.frame(Basal_NormNoOutliers_OverlappingGenesNotAll), cex=0.25, las=3, par(cex.axis=0.75))
#The signal on some of these is actually pretty decent - boxes are more centered around 7.5






Central_IndicesForOverlappingGenes<- matrix(c(

row.names(Central_MDD_Pvalues)%in%row.names(AAA_MDD_Pvalues),
row.names(Central_MDD_Pvalues)%in%row.names(AB_MDD_Pvalues),
row.names(Central_MDD_Pvalues)%in%row.names(AHA_MDD_Pvalues),
row.names(Central_MDD_Pvalues)%in%row.names(Basal_MDD_Pvalues),
row.names(Central_MDD_Pvalues)%in%row.names(CO_MDD_Pvalues),
row.names(Central_MDD_Pvalues)%in%row.names(Lateral_MDD_Pvalues),
row.names(Central_MDD_Pvalues)%in%row.names(Medial_MDD_Pvalues),
row.names(Central_MDD_Pvalues)%in%row.names(PAC_MDD_Pvalues),
row.names(Central_MDD_Pvalues)%in%row.names(PL_MDD_Pvalues)
), 
nrow=length(Central_MDD_Pvalues), ncol=9)

Central_IndicesForOverlappingGenes[c(1:3),]

Central_IndicesForOverlappingGenes_All<-matrix(T, length(Central_MDD_Pvalues[,1]), 1)

for(i in 1:length(Central_MDD_Pvalues[,1])){	
Central_IndicesForOverlappingGenes_All[i,1]<-sum(Central_IndicesForOverlappingGenes[i,]==T)==9
}

sum(Central_IndicesForOverlappingGenes_All[,1]==T)
#[1] 15788
#Out of 18704 original probes


Central_Age_Pvalues_OverlappingGenes<-Central_Age_Pvalues[(Central_IndicesForOverlappingGenes_All[,1]==T),]
Central_RNAIntegrity_Pvalues_OverlappingGenes<-Central_RNAIntegrity_Pvalues[(Central_IndicesForOverlappingGenes_All[,1]==T),]
Central_Gender_Pvalues_OverlappingGenes<-Central_Gender_Pvalues[(Central_IndicesForOverlappingGenes_All[,1]==T),]
Central_HoursFinal_Pvalues_OverlappingGenes<-Central_HoursFinal_Pvalues[(Central_IndicesForOverlappingGenes_All[,1]==T),]
Central_MDD_Pvalues_OverlappingGenes<-Central_MDD_Pvalues[(Central_IndicesForOverlappingGenes_All[,1]==T),]
Central_PH_Pvalues_OverlappingGenes<-Central_PH_Pvalues[(Central_IndicesForOverlappingGenes_All[,1]==T),]

Central_NormNoOutliers_OverlappingGenes<-Central_NormNoOutliers[(Central_IndicesForOverlappingGenes_All[,1]==T),]



length(Central_MDD_Pvalues_OverlappingGenes)
Central_MDD_Pvalues_OverlappingGenes[c(1:5)]




Central_IndicesForOverlappingGenes2<- matrix(c(
row.names(Central_MDD_Pvalues)%in%row.names(AAA_MDD_Pvalues),
row.names(Central_MDD_Pvalues)%in%row.names(AB_MDD_Pvalues),
row.names(Central_MDD_Pvalues)%in%row.names(AHA_MDD_Pvalues),
row.names(Central_MDD_Pvalues)%in%row.names(Basal_MDD_Pvalues),
row.names(Central_MDD_Pvalues)%in%row.names(Central_MDD_Pvalues),
row.names(Central_MDD_Pvalues)%in%row.names(CO_MDD_Pvalues),
row.names(Central_MDD_Pvalues)%in%row.names(Lateral_MDD_Pvalues),
row.names(Central_MDD_Pvalues)%in%row.names(Medial_MDD_Pvalues),
row.names(Central_MDD_Pvalues)%in%row.names(PAC_MDD_Pvalues),
row.names(Central_MDD_Pvalues)%in%row.names(PL_MDD_Pvalues)
), 
nrow=length(Central_MDD_Pvalues), ncol=10)


Central_IndicesForOverlappingGenes2[c(1:3),]

Central_IndicesForOverlappingGenes_Any<-matrix(T, length(Central_MDD_Pvalues[,1]), 1)

for(i in 1:length(Central_MDD_Pvalues[,1])){	
Central_IndicesForOverlappingGenes_Any[i,1]<-sum(Central_IndicesForOverlappingGenes2[i,]==T)>1
}

sum(Central_IndicesForOverlappingGenes_Any[,1]==T)
[1] 18632
18632-15788
[1] 2844

TotalOverlappingCentralwOther<-apply(Central_IndicesForOverlappingGenes2, 2, sum)
names(TotalOverlappingCentralwOther)<-c("AAA", "AB", "AHA", "Basal", "Central", "CO", "Lateral", "Medial", "PAC", "PL")

sum((Central_IndicesForOverlappingGenes_Any[,1]==T)&(Central_IndicesForOverlappingGenes_All[,1]==F))
[1] 2844


Central_NormNoOutliers_OverlappingGenesNotAll<-Central_NormNoOutliers[(Central_IndicesForOverlappingGenes_Any[,1]==T)&(Central_IndicesForOverlappingGenes_All[,1]==F),]


boxplot(data.frame(Central_NormNoOutliers_OverlappingGenesNotAll), cex=0.25, las=3, par(cex.axis=0.75))
#The signal on some of these is actually pretty decent - boxes are more centered around 7.5





CO_IndicesForOverlappingGenes<- matrix(c(

row.names(CO_MDD_Pvalues)%in%row.names(AAA_MDD_Pvalues),
row.names(CO_MDD_Pvalues)%in%row.names(AB_MDD_Pvalues),
row.names(CO_MDD_Pvalues)%in%row.names(AHA_MDD_Pvalues),
row.names(CO_MDD_Pvalues)%in%row.names(Basal_MDD_Pvalues),
row.names(CO_MDD_Pvalues)%in%row.names(Central_MDD_Pvalues),
row.names(CO_MDD_Pvalues)%in%row.names(Lateral_MDD_Pvalues),
row.names(CO_MDD_Pvalues)%in%row.names(Medial_MDD_Pvalues),
row.names(CO_MDD_Pvalues)%in%row.names(PAC_MDD_Pvalues),
row.names(CO_MDD_Pvalues)%in%row.names(PL_MDD_Pvalues)
), 
nrow=length(CO_MDD_Pvalues), ncol=9)

CO_IndicesForOverlappingGenes[c(1:3),]

CO_IndicesForOverlappingGenes_All<-matrix(T, length(CO_MDD_Pvalues[,1]), 1)

for(i in 1:length(CO_MDD_Pvalues[,1])){	
CO_IndicesForOverlappingGenes_All[i,1]<-sum(CO_IndicesForOverlappingGenes[i,]==T)==9
}

sum(CO_IndicesForOverlappingGenes_All[,1]==T)
#[1] 15788
#Out of 18826 original probes


CO_Age_Pvalues_OverlappingGenes<-CO_Age_Pvalues[(CO_IndicesForOverlappingGenes_All[,1]==T),]
CO_RNAIntegrity_Pvalues_OverlappingGenes<-CO_RNAIntegrity_Pvalues[(CO_IndicesForOverlappingGenes_All[,1]==T),]
CO_Gender_Pvalues_OverlappingGenes<-CO_Gender_Pvalues[(CO_IndicesForOverlappingGenes_All[,1]==T),]
CO_HoursFinal_Pvalues_OverlappingGenes<-CO_HoursFinal_Pvalues[(CO_IndicesForOverlappingGenes_All[,1]==T),]
CO_MDD_Pvalues_OverlappingGenes<-CO_MDD_Pvalues[(CO_IndicesForOverlappingGenes_All[,1]==T),]
CO_PH_Pvalues_OverlappingGenes<-CO_PH_Pvalues[(CO_IndicesForOverlappingGenes_All[,1]==T),]

CO_NormNoOutliers_OverlappingGenes<-CO_NormNoOutliers[(CO_IndicesForOverlappingGenes_All[,1]==T),]



length(CO_MDD_Pvalues_OverlappingGenes)
CO_MDD_Pvalues_OverlappingGenes[c(1:5)]



CO_IndicesForOverlappingGenes2<- matrix(c(
row.names(CO_MDD_Pvalues)%in%row.names(AAA_MDD_Pvalues),
row.names(CO_MDD_Pvalues)%in%row.names(AB_MDD_Pvalues),
row.names(CO_MDD_Pvalues)%in%row.names(AHA_MDD_Pvalues),
row.names(CO_MDD_Pvalues)%in%row.names(Basal_MDD_Pvalues),
row.names(CO_MDD_Pvalues)%in%row.names(Central_MDD_Pvalues),
row.names(CO_MDD_Pvalues)%in%row.names(CO_MDD_Pvalues),
row.names(CO_MDD_Pvalues)%in%row.names(Lateral_MDD_Pvalues),
row.names(CO_MDD_Pvalues)%in%row.names(Medial_MDD_Pvalues),
row.names(CO_MDD_Pvalues)%in%row.names(PAC_MDD_Pvalues),
row.names(CO_MDD_Pvalues)%in%row.names(PL_MDD_Pvalues)
), 
nrow=length(CO_MDD_Pvalues), ncol=10)


CO_IndicesForOverlappingGenes2[c(1:3),]

CO_IndicesForOverlappingGenes_Any<-matrix(T, length(CO_MDD_Pvalues[,1]), 1)

for(i in 1:length(CO_MDD_Pvalues[,1])){	
CO_IndicesForOverlappingGenes_Any[i,1]<-sum(CO_IndicesForOverlappingGenes2[i,]==T)>1
}

sum(CO_IndicesForOverlappingGenes_Any[,1]==T)
[1] 18753
18753-15788
[1] 2965

TotalOverlappingCOwOther<-apply(CO_IndicesForOverlappingGenes2, 2, sum)
names(TotalOverlappingCOwOther)<-c("AAA", "AB", "AHA", "Basal", "Central", "CO", "Lateral", "Medial", "PAC", "PL")

sum((CO_IndicesForOverlappingGenes_Any[,1]==T)&(CO_IndicesForOverlappingGenes_All[,1]==F))
[1] 2965


CO_NormNoOutliers_OverlappingGenesNotAll<-CO_NormNoOutliers[(CO_IndicesForOverlappingGenes_Any[,1]==T)&(CO_IndicesForOverlappingGenes_All[,1]==F),]


boxplot(data.frame(CO_NormNoOutliers_OverlappingGenesNotAll), cex=0.25, las=3, par(cex.axis=0.75))
#Very low - box centers around 7.





Lateral_IndicesForOverlappingGenes<- matrix(c(

row.names(Lateral_MDD_Pvalues)%in%row.names(AAA_MDD_Pvalues),
row.names(Lateral_MDD_Pvalues)%in%row.names(AB_MDD_Pvalues),
row.names(Lateral_MDD_Pvalues)%in%row.names(AHA_MDD_Pvalues),
row.names(Lateral_MDD_Pvalues)%in%row.names(Basal_MDD_Pvalues),
row.names(Lateral_MDD_Pvalues)%in%row.names(Central_MDD_Pvalues),
row.names(Lateral_MDD_Pvalues)%in%row.names(CO_MDD_Pvalues),
row.names(Lateral_MDD_Pvalues)%in%row.names(Medial_MDD_Pvalues),
row.names(Lateral_MDD_Pvalues)%in%row.names(PAC_MDD_Pvalues),
row.names(Lateral_MDD_Pvalues)%in%row.names(PL_MDD_Pvalues)
), 
nrow=length(Lateral_MDD_Pvalues), ncol=9)

Lateral_IndicesForOverlappingGenes[c(1:3),]

Lateral_IndicesForOverlappingGenes_All<-matrix(T, length(Lateral_MDD_Pvalues[,1]), 1)

for(i in 1:length(Lateral_MDD_Pvalues[,1])){	
Lateral_IndicesForOverlappingGenes_All[i,1]<-sum(Lateral_IndicesForOverlappingGenes[i,]==T)==9
}

sum(Lateral_IndicesForOverlappingGenes_All[,1]==T)
#[1] 15788
#Out of 18603 original probes


Lateral_Age_Pvalues_OverlappingGenes<-Lateral_Age_Pvalues[(Lateral_IndicesForOverlappingGenes_All[,1]==T),]
Lateral_RNAIntegrity_Pvalues_OverlappingGenes<-Lateral_RNAIntegrity_Pvalues[(Lateral_IndicesForOverlappingGenes_All[,1]==T),]
Lateral_Gender_Pvalues_OverlappingGenes<-Lateral_Gender_Pvalues[(Lateral_IndicesForOverlappingGenes_All[,1]==T),]
Lateral_HoursFinal_Pvalues_OverlappingGenes<-Lateral_HoursFinal_Pvalues[(Lateral_IndicesForOverlappingGenes_All[,1]==T),]
Lateral_MDD_Pvalues_OverlappingGenes<-Lateral_MDD_Pvalues[(Lateral_IndicesForOverlappingGenes_All[,1]==T),]
Lateral_PH_Pvalues_OverlappingGenes<-Lateral_PH_Pvalues[(Lateral_IndicesForOverlappingGenes_All[,1]==T),]

Lateral_NormNoOutliers_OverlappingGenes<-Lateral_NormNoOutliers[(Lateral_IndicesForOverlappingGenes_All[,1]==T),]



length(Lateral_MDD_Pvalues_OverlappingGenes)
Lateral_MDD_Pvalues_OverlappingGenes[c(1:5)]




Lateral_IndicesForOverlappingGenes2<- matrix(c(
row.names(Lateral_MDD_Pvalues)%in%row.names(AAA_MDD_Pvalues),
row.names(Lateral_MDD_Pvalues)%in%row.names(AB_MDD_Pvalues),
row.names(Lateral_MDD_Pvalues)%in%row.names(AHA_MDD_Pvalues),
row.names(Lateral_MDD_Pvalues)%in%row.names(Basal_MDD_Pvalues),
row.names(Lateral_MDD_Pvalues)%in%row.names(Central_MDD_Pvalues),
row.names(Lateral_MDD_Pvalues)%in%row.names(CO_MDD_Pvalues),
row.names(Lateral_MDD_Pvalues)%in%row.names(Lateral_MDD_Pvalues),
row.names(Lateral_MDD_Pvalues)%in%row.names(Medial_MDD_Pvalues),
row.names(Lateral_MDD_Pvalues)%in%row.names(PAC_MDD_Pvalues),
row.names(Lateral_MDD_Pvalues)%in%row.names(PL_MDD_Pvalues)
), 
nrow=length(Lateral_MDD_Pvalues), ncol=10)


Lateral_IndicesForOverlappingGenes2[c(1:3),]

Lateral_IndicesForOverlappingGenes_Any<-matrix(T, length(Lateral_MDD_Pvalues[,1]), 1)

for(i in 1:length(Lateral_MDD_Pvalues[,1])){	
Lateral_IndicesForOverlappingGenes_Any[i,1]<-sum(Lateral_IndicesForOverlappingGenes2[i,]==T)>1
}

sum(Lateral_IndicesForOverlappingGenes_Any[,1]==T)
[1] 18343
18343-15788
[1] 2555

TotalOverlappingLateralwOther<-apply(Lateral_IndicesForOverlappingGenes2, 2, sum)
names(TotalOverlappingLateralwOther)<-c("AAA", "AB", "AHA", "Basal", "Central", "CO", "Lateral", "Medial", "PAC", "PL")

sum((Lateral_IndicesForOverlappingGenes_Any[,1]==T)&(Lateral_IndicesForOverlappingGenes_All[,1]==F))
[1] 2555


Lateral_NormNoOutliers_OverlappingGenesNotAll<-Lateral_NormNoOutliers[(Lateral_IndicesForOverlappingGenes_Any[,1]==T)&(Lateral_IndicesForOverlappingGenes_All[,1]==F),]


boxplot(data.frame(Lateral_NormNoOutliers_OverlappingGenesNotAll), cex=0.25, las=3, par(cex.axis=0.75))
#Very low - boxes *tightly* center around 7. 






Medial_IndicesForOverlappingGenes<- matrix(c(

row.names(Medial_MDD_Pvalues)%in%row.names(AAA_MDD_Pvalues),
row.names(Medial_MDD_Pvalues)%in%row.names(AB_MDD_Pvalues),
row.names(Medial_MDD_Pvalues)%in%row.names(AHA_MDD_Pvalues),
row.names(Medial_MDD_Pvalues)%in%row.names(Basal_MDD_Pvalues),
row.names(Medial_MDD_Pvalues)%in%row.names(Central_MDD_Pvalues),
row.names(Medial_MDD_Pvalues)%in%row.names(CO_MDD_Pvalues),
row.names(Medial_MDD_Pvalues)%in%row.names(Lateral_MDD_Pvalues),
row.names(Medial_MDD_Pvalues)%in%row.names(PAC_MDD_Pvalues),
row.names(Medial_MDD_Pvalues)%in%row.names(PL_MDD_Pvalues)
), 
nrow=length(Medial_MDD_Pvalues), ncol=9)

Medial_IndicesForOverlappingGenes[c(1:3),]

Medial_IndicesForOverlappingGenes_All<-matrix(T, length(Medial_MDD_Pvalues[,1]), 1)

for(i in 1:length(Medial_MDD_Pvalues[,1])){	
Medial_IndicesForOverlappingGenes_All[i,1]<-sum(Medial_IndicesForOverlappingGenes[i,]==T)==9
}

sum(Medial_IndicesForOverlappingGenes_All[,1]==T)
#[1] 15788
#Out of 19820 original probes


Medial_Age_Pvalues_OverlappingGenes<-Medial_Age_Pvalues[(Medial_IndicesForOverlappingGenes_All[,1]==T),]
Medial_RNAIntegrity_Pvalues_OverlappingGenes<-Medial_RNAIntegrity_Pvalues[(Medial_IndicesForOverlappingGenes_All[,1]==T),]
Medial_Gender_Pvalues_OverlappingGenes<-Medial_Gender_Pvalues[(Medial_IndicesForOverlappingGenes_All[,1]==T),]
Medial_HoursFinal_Pvalues_OverlappingGenes<-Medial_HoursFinal_Pvalues[(Medial_IndicesForOverlappingGenes_All[,1]==T),]
Medial_MDD_Pvalues_OverlappingGenes<-Medial_MDD_Pvalues[(Medial_IndicesForOverlappingGenes_All[,1]==T),]
Medial_PH_Pvalues_OverlappingGenes<-Medial_PH_Pvalues[(Medial_IndicesForOverlappingGenes_All[,1]==T),]

Medial_NormNoOutliers_OverlappingGenes<-Medial_NormNoOutliers[(Medial_IndicesForOverlappingGenes_All[,1]==T),]



length(Medial_MDD_Pvalues_OverlappingGenes)
Medial_MDD_Pvalues_OverlappingGenes[c(1:5)]


Medial_IndicesForOverlappingGenes2<- matrix(c(
row.names(Medial_MDD_Pvalues)%in%row.names(AAA_MDD_Pvalues),
row.names(Medial_MDD_Pvalues)%in%row.names(AB_MDD_Pvalues),
row.names(Medial_MDD_Pvalues)%in%row.names(AHA_MDD_Pvalues),
row.names(Medial_MDD_Pvalues)%in%row.names(Basal_MDD_Pvalues),
row.names(Medial_MDD_Pvalues)%in%row.names(Central_MDD_Pvalues),
row.names(Medial_MDD_Pvalues)%in%row.names(CO_MDD_Pvalues),
row.names(Medial_MDD_Pvalues)%in%row.names(Lateral_MDD_Pvalues),
row.names(Medial_MDD_Pvalues)%in%row.names(Medial_MDD_Pvalues),
row.names(Medial_MDD_Pvalues)%in%row.names(PAC_MDD_Pvalues),
row.names(Medial_MDD_Pvalues)%in%row.names(PL_MDD_Pvalues)
), 
nrow=length(Medial_MDD_Pvalues), ncol=10)


Medial_IndicesForOverlappingGenes2[c(1:3),]

Medial_IndicesForOverlappingGenes_Any<-matrix(T, length(Medial_MDD_Pvalues[,1]), 1)

for(i in 1:length(Medial_MDD_Pvalues[,1])){	
Medial_IndicesForOverlappingGenes_Any[i,1]<-sum(Medial_IndicesForOverlappingGenes2[i,]==T)>1
}

sum(Medial_IndicesForOverlappingGenes_Any[,1]==T)
[1] 19665
19665-15788
[1] 3877

TotalOverlappingMedialwOther<-apply(Medial_IndicesForOverlappingGenes2, 2, sum)
names(TotalOverlappingMedialwOther)<-c("AAA", "AB", "AHA", "Basal", "Central", "CO", "Lateral", "Medial", "PAC", "PL")

sum((Medial_IndicesForOverlappingGenes_Any[,1]==T)&(Medial_IndicesForOverlappingGenes_All[,1]==F))
[1] 3877


Medial_NormNoOutliers_OverlappingGenesNotAll<-Medial_NormNoOutliers[(Medial_IndicesForOverlappingGenes_Any[,1]==T)&(Medial_IndicesForOverlappingGenes_All[,1]==F),]


boxplot(data.frame(Medial_NormNoOutliers_OverlappingGenesNotAll), cex=0.25, las=3, par(cex.axis=0.75))
#Very low - boxes center around 7. 





PAC_IndicesForOverlappingGenes<- matrix(c(

row.names(PAC_MDD_Pvalues)%in%row.names(AAA_MDD_Pvalues),
row.names(PAC_MDD_Pvalues)%in%row.names(AB_MDD_Pvalues),
row.names(PAC_MDD_Pvalues)%in%row.names(AHA_MDD_Pvalues),
row.names(PAC_MDD_Pvalues)%in%row.names(Basal_MDD_Pvalues),
row.names(PAC_MDD_Pvalues)%in%row.names(Central_MDD_Pvalues),
row.names(PAC_MDD_Pvalues)%in%row.names(CO_MDD_Pvalues),
row.names(PAC_MDD_Pvalues)%in%row.names(Lateral_MDD_Pvalues),
row.names(PAC_MDD_Pvalues)%in%row.names(Medial_MDD_Pvalues),
row.names(PAC_MDD_Pvalues)%in%row.names(PL_MDD_Pvalues)
), 
nrow=length(PAC_MDD_Pvalues), ncol=9)

PAC_IndicesForOverlappingGenes[c(1:3),]

PAC_IndicesForOverlappingGenes_All<-matrix(T, length(PAC_MDD_Pvalues[,1]), 1)

for(i in 1:length(PAC_MDD_Pvalues[,1])){	
PAC_IndicesForOverlappingGenes_All[i,1]<-sum(PAC_IndicesForOverlappingGenes[i,]==T)==9
}

sum(PAC_IndicesForOverlappingGenes_All[,1]==T)
#[1] 15788
#Out of 15788 original probes.


PAC_Age_Pvalues_OverlappingGenes<-PAC_Age_Pvalues[(PAC_IndicesForOverlappingGenes_All[,1]==T),]
PAC_RNAIntegrity_Pvalues_OverlappingGenes<-PAC_RNAIntegrity_Pvalues[(PAC_IndicesForOverlappingGenes_All[,1]==T),]
PAC_Gender_Pvalues_OverlappingGenes<-PAC_Gender_Pvalues[(PAC_IndicesForOverlappingGenes_All[,1]==T),]
PAC_HoursFinal_Pvalues_OverlappingGenes<-PAC_HoursFinal_Pvalues[(PAC_IndicesForOverlappingGenes_All[,1]==T),]
PAC_MDD_Pvalues_OverlappingGenes<-PAC_MDD_Pvalues[(PAC_IndicesForOverlappingGenes_All[,1]==T),]
PAC_PH_Pvalues_OverlappingGenes<-PAC_PH_Pvalues[(PAC_IndicesForOverlappingGenes_All[,1]==T),]


PAC_NormNoOutliers_OverlappingGenes<-PAC_NormNoOutliers[(PAC_IndicesForOverlappingGenes_All[,1]==T),]



length(PAC_MDD_Pvalues_OverlappingGenes)
PAC_MDD_Pvalues_OverlappingGenes[c(1:5)]



PAC_IndicesForOverlappingGenes2<- matrix(c(
row.names(PAC_MDD_Pvalues)%in%row.names(AAA_MDD_Pvalues),
row.names(PAC_MDD_Pvalues)%in%row.names(AB_MDD_Pvalues),
row.names(PAC_MDD_Pvalues)%in%row.names(AHA_MDD_Pvalues),
row.names(PAC_MDD_Pvalues)%in%row.names(Basal_MDD_Pvalues),
row.names(PAC_MDD_Pvalues)%in%row.names(Central_MDD_Pvalues),
row.names(PAC_MDD_Pvalues)%in%row.names(CO_MDD_Pvalues),
row.names(PAC_MDD_Pvalues)%in%row.names(Lateral_MDD_Pvalues),
row.names(PAC_MDD_Pvalues)%in%row.names(Medial_MDD_Pvalues),
row.names(PAC_MDD_Pvalues)%in%row.names(PAC_MDD_Pvalues),
row.names(PAC_MDD_Pvalues)%in%row.names(PL_MDD_Pvalues)
), 
nrow=length(PAC_MDD_Pvalues), ncol=10)


PAC_IndicesForOverlappingGenes2[c(1:3),]

PAC_IndicesForOverlappingGenes_Any<-matrix(T, length(PAC_MDD_Pvalues[,1]), 1)

for(i in 1:length(PAC_MDD_Pvalues[,1])){	
PAC_IndicesForOverlappingGenes_Any[i,1]<-sum(PAC_IndicesForOverlappingGenes2[i,]==T)>1
}

sum(PAC_IndicesForOverlappingGenes_Any[,1]==T)
[1] 18352
18352-15788
[1] 2564

TotalOverlappingPACwOther<-apply(PAC_IndicesForOverlappingGenes2, 2, sum)
names(TotalOverlappingPACwOther)<-c("AAA", "AB", "AHA", "Basal", "Central", "CO", "Lateral", "Medial", "PAC", "PL")

sum((PAC_IndicesForOverlappingGenes_Any[,1]==T)&(PAC_IndicesForOverlappingGenes_All[,1]==F))
[1] 2564


PAC_NormNoOutliers_OverlappingGenesNotAll<-PAC_NormNoOutliers[(PAC_IndicesForOverlappingGenes_Any[,1]==T)&(PAC_IndicesForOverlappingGenes_All[,1]==F),]


boxplot(data.frame(PAC_NormNoOutliers_OverlappingGenesNotAll), cex=0.25, las=3, par(cex.axis=0.75))
#Like Basal, many of these are quite high, centered on 7.5 or higher.  I'm going to guess that this is because many of the genes in Basal are found in PAC (since they are so similar)





PL_IndicesForOverlappingGenes<- matrix(c(

row.names(PL_MDD_Pvalues)%in%row.names(AAA_MDD_Pvalues),
row.names(PL_MDD_Pvalues)%in%row.names(AB_MDD_Pvalues),
row.names(PL_MDD_Pvalues)%in%row.names(AHA_MDD_Pvalues),
row.names(PL_MDD_Pvalues)%in%row.names(Basal_MDD_Pvalues),
row.names(PL_MDD_Pvalues)%in%row.names(Central_MDD_Pvalues),
row.names(PL_MDD_Pvalues)%in%row.names(CO_MDD_Pvalues),
row.names(PL_MDD_Pvalues)%in%row.names(Lateral_MDD_Pvalues),
row.names(PL_MDD_Pvalues)%in%row.names(Medial_MDD_Pvalues),
row.names(PL_MDD_Pvalues)%in%row.names(PAC_MDD_Pvalues)
), 
nrow=length(PL_MDD_Pvalues), ncol=9)

PL_IndicesForOverlappingGenes[c(1:3),]

PL_IndicesForOverlappingGenes_All<-matrix(T, length(PL_MDD_Pvalues[,1]), 1)

for(i in 1:length(PL_MDD_Pvalues[,1])){	
PL_IndicesForOverlappingGenes_All[i,1]<-sum(PL_IndicesForOverlappingGenes[i,]==T)==9
}

sum(PL_IndicesForOverlappingGenes_All[,1]==T)
#[1] 15788
#Out of 19026 original probes


PL_Age_Pvalues_OverlappingGenes<-PL_Age_Pvalues[(PL_IndicesForOverlappingGenes_All[,1]==T),]
PL_RNAIntegrity_Pvalues_OverlappingGenes<-PL_RNAIntegrity_Pvalues[(PL_IndicesForOverlappingGenes_All[,1]==T),]
PL_Gender_Pvalues_OverlappingGenes<-PL_Gender_Pvalues[(PL_IndicesForOverlappingGenes_All[,1]==T),]
PL_HoursFinal_Pvalues_OverlappingGenes<-PL_HoursFinal_Pvalues[(PL_IndicesForOverlappingGenes_All[,1]==T),]
PL_MDD_Pvalues_OverlappingGenes<-PL_MDD_Pvalues[(PL_IndicesForOverlappingGenes_All[,1]==T),]
PL_PH_Pvalues_OverlappingGenes<-PL_PH_Pvalues[(PL_IndicesForOverlappingGenes_All[,1]==T),]

PL_NormNoOutliers_OverlappingGenes<-PL_NormNoOutliers[(PL_IndicesForOverlappingGenes_All[,1]==T),]




length(PL_MDD_Pvalues_OverlappingGenes)
PL_MDD_Pvalues_OverlappingGenes[c(1:5)]



PL_IndicesForOverlappingGenes2<- matrix(c(
row.names(PL_MDD_Pvalues)%in%row.names(AAA_MDD_Pvalues),
row.names(PL_MDD_Pvalues)%in%row.names(AB_MDD_Pvalues),
row.names(PL_MDD_Pvalues)%in%row.names(AHA_MDD_Pvalues),
row.names(PL_MDD_Pvalues)%in%row.names(Basal_MDD_Pvalues),
row.names(PL_MDD_Pvalues)%in%row.names(Central_MDD_Pvalues),
row.names(PL_MDD_Pvalues)%in%row.names(CO_MDD_Pvalues),
row.names(PL_MDD_Pvalues)%in%row.names(Lateral_MDD_Pvalues),
row.names(PL_MDD_Pvalues)%in%row.names(Medial_MDD_Pvalues),
row.names(PL_MDD_Pvalues)%in%row.names(PAC_MDD_Pvalues),
row.names(PL_MDD_Pvalues)%in%row.names(PL_MDD_Pvalues)
), 
nrow=length(PL_MDD_Pvalues), ncol=10)


PL_IndicesForOverlappingGenes2[c(1:3),]

PL_IndicesForOverlappingGenes_Any<-matrix(T, length(PL_MDD_Pvalues[,1]), 1)

for(i in 1:length(PL_MDD_Pvalues[,1])){	
PL_IndicesForOverlappingGenes_Any[i,1]<-sum(PL_IndicesForOverlappingGenes2[i,]==T)>1
}

sum(PL_IndicesForOverlappingGenes_Any[,1]==T)
[1] 18950
18950-15788
[1] 3162

TotalOverlappingPLwOther<-apply(PL_IndicesForOverlappingGenes2, 2, sum)
names(TotalOverlappingPLwOther)<-c("AAA", "AB", "AHA", "Basal", "Central", "CO", "Lateral", "Medial", "PAC", "PL")

sum((PL_IndicesForOverlappingGenes_Any[,1]==T)&(PL_IndicesForOverlappingGenes_All[,1]==F))
[1] 3162


PL_NormNoOutliers_OverlappingGenesNotAll<-PL_NormNoOutliers[(PL_IndicesForOverlappingGenes_Any[,1]==T)&(PL_IndicesForOverlappingGenes_All[,1]==F),]


boxplot(data.frame(PL_NormNoOutliers_OverlappingGenesNotAll), cex=0.25, las=3, par(cex.axis=0.75))
#Wow, these are really really low - box centers on 6.5!!



TotalOverlappingProbes10Nuclei<-matrix(
c(
TotalOverlappingAAAwOther,
TotalOverlappingABwOther,
TotalOverlappingAHAwOther,
TotalOverlappingBasalwOther,
TotalOverlappingCentralwOther,
TotalOverlappingCOwOther,
TotalOverlappingLateralwOther,
TotalOverlappingMedialwOther,
TotalOverlappingPACwOther,
TotalOverlappingPLwOther
), nrow=10, ncol=10
)


row.names(TotalOverlappingProbes10Nuclei)<-c("AAA", "AB", "AHA", "Basal", "Central", "CO", "Lateral", "Medial", "PAC", "PL")
colnames(TotalOverlappingProbes10Nuclei)<-c("AAA", "AB", "AHA", "Basal", "Central", "CO", "Lateral", "Medial", "PAC", "PL")


write.csv(TotalOverlappingProbes10Nuclei, "TotalOverlappingProbes10Nuclei.csv")







AAA_Age_Betas_OverlappingGenes<-AAA_Age_Betas[(AAA_IndicesForOverlappingGenes_All[,1]==T),]
AAA_RNAIntegrity_Betas_OverlappingGenes<-AAA_RNAIntegrity_Betas[(AAA_IndicesForOverlappingGenes_All[,1]==T),]
AAA_Gender_Betas_OverlappingGenes<-AAA_Gender_Betas[(AAA_IndicesForOverlappingGenes_All[,1]==T),]
AAA_HoursFinal_Betas_OverlappingGenes<-AAA_HoursFinal_Betas[(AAA_IndicesForOverlappingGenes_All[,1]==T),]
AAA_MDD_Betas_OverlappingGenes<-AAA_MDD_Betas[(AAA_IndicesForOverlappingGenes_All[,1]==T),]
AAA_PH_Betas_OverlappingGenes<-AAA_PH_Betas[(AAA_IndicesForOverlappingGenes_All[,1]==T),]


AB_Age_Betas_OverlappingGenes<-AB_Age_Betas[(AB_IndicesForOverlappingGenes_All[,1]==T),]
AB_RNAIntegrity_Betas_OverlappingGenes<-AB_RNAIntegrity_Betas[(AB_IndicesForOverlappingGenes_All[,1]==T),]
AB_Gender_Betas_OverlappingGenes<-AB_Gender_Betas[(AB_IndicesForOverlappingGenes_All[,1]==T),]
AB_HoursFinal_Betas_OverlappingGenes<-AB_HoursFinal_Betas[(AB_IndicesForOverlappingGenes_All[,1]==T),]
AB_MDD_Betas_OverlappingGenes<-AB_MDD_Betas[(AB_IndicesForOverlappingGenes_All[,1]==T),]
AB_PH_Betas_OverlappingGenes<-AB_PH_Betas[(AB_IndicesForOverlappingGenes_All[,1]==T),]


AHA_Age_Betas_OverlappingGenes<-AHA_Age_Betas[(AHA_IndicesForOverlappingGenes_All[,1]==T),]
AHA_RNAIntegrity_Betas_OverlappingGenes<-AHA_RNAIntegrity_Betas[(AHA_IndicesForOverlappingGenes_All[,1]==T),]
AHA_Gender_Betas_OverlappingGenes<-AHA_Gender_Betas[(AHA_IndicesForOverlappingGenes_All[,1]==T),]
AHA_HoursFinal_Betas_OverlappingGenes<-AHA_HoursFinal_Betas[(AHA_IndicesForOverlappingGenes_All[,1]==T),]
AHA_MDD_Betas_OverlappingGenes<-AHA_MDD_Betas[(AHA_IndicesForOverlappingGenes_All[,1]==T),]
AHA_PH_Betas_OverlappingGenes<-AHA_PH_Betas[(AHA_IndicesForOverlappingGenes_All[,1]==T),]

Basal_Age_Betas_OverlappingGenes<-Basal_Age_Betas[(Basal_IndicesForOverlappingGenes_All[,1]==T),]
Basal_RNAIntegrity_Betas_OverlappingGenes<-Basal_RNAIntegrity_Betas[(Basal_IndicesForOverlappingGenes_All[,1]==T),]
Basal_Gender_Betas_OverlappingGenes<-Basal_Gender_Betas[(Basal_IndicesForOverlappingGenes_All[,1]==T),]
Basal_HoursFinal_Betas_OverlappingGenes<-Basal_HoursFinal_Betas[(Basal_IndicesForOverlappingGenes_All[,1]==T),]
Basal_MDD_Betas_OverlappingGenes<-Basal_MDD_Betas[(Basal_IndicesForOverlappingGenes_All[,1]==T),]
Basal_PH_Betas_OverlappingGenes<-Basal_PH_Betas[(Basal_IndicesForOverlappingGenes_All[,1]==T),]


Central_Age_Betas_OverlappingGenes<-Central_Age_Betas[(Central_IndicesForOverlappingGenes_All[,1]==T),]
Central_RNAIntegrity_Betas_OverlappingGenes<-Central_RNAIntegrity_Betas[(Central_IndicesForOverlappingGenes_All[,1]==T),]
Central_Gender_Betas_OverlappingGenes<-Central_Gender_Betas[(Central_IndicesForOverlappingGenes_All[,1]==T),]
Central_HoursFinal_Betas_OverlappingGenes<-Central_HoursFinal_Betas[(Central_IndicesForOverlappingGenes_All[,1]==T),]
Central_MDD_Betas_OverlappingGenes<-Central_MDD_Betas[(Central_IndicesForOverlappingGenes_All[,1]==T),]
Central_PH_Betas_OverlappingGenes<-Central_PH_Betas[(Central_IndicesForOverlappingGenes_All[,1]==T),]

CO_Age_Betas_OverlappingGenes<-CO_Age_Betas[(CO_IndicesForOverlappingGenes_All[,1]==T),]
CO_RNAIntegrity_Betas_OverlappingGenes<-CO_RNAIntegrity_Betas[(CO_IndicesForOverlappingGenes_All[,1]==T),]
CO_Gender_Betas_OverlappingGenes<-CO_Gender_Betas[(CO_IndicesForOverlappingGenes_All[,1]==T),]
CO_HoursFinal_Betas_OverlappingGenes<-CO_HoursFinal_Betas[(CO_IndicesForOverlappingGenes_All[,1]==T),]
CO_MDD_Betas_OverlappingGenes<-CO_MDD_Betas[(CO_IndicesForOverlappingGenes_All[,1]==T),]
CO_PH_Betas_OverlappingGenes<-CO_PH_Betas[(CO_IndicesForOverlappingGenes_All[,1]==T),]


Lateral_Age_Betas_OverlappingGenes<-Lateral_Age_Betas[(Lateral_IndicesForOverlappingGenes_All[,1]==T),]
Lateral_RNAIntegrity_Betas_OverlappingGenes<-Lateral_RNAIntegrity_Betas[(Lateral_IndicesForOverlappingGenes_All[,1]==T),]
Lateral_Gender_Betas_OverlappingGenes<-Lateral_Gender_Betas[(Lateral_IndicesForOverlappingGenes_All[,1]==T),]
Lateral_HoursFinal_Betas_OverlappingGenes<-Lateral_HoursFinal_Betas[(Lateral_IndicesForOverlappingGenes_All[,1]==T),]
Lateral_MDD_Betas_OverlappingGenes<-Lateral_MDD_Betas[(Lateral_IndicesForOverlappingGenes_All[,1]==T),]
Lateral_PH_Betas_OverlappingGenes<-Lateral_PH_Betas[(Lateral_IndicesForOverlappingGenes_All[,1]==T),]


Medial_Age_Betas_OverlappingGenes<-Medial_Age_Betas[(Medial_IndicesForOverlappingGenes_All[,1]==T),]
Medial_RNAIntegrity_Betas_OverlappingGenes<-Medial_RNAIntegrity_Betas[(Medial_IndicesForOverlappingGenes_All[,1]==T),]
Medial_Gender_Betas_OverlappingGenes<-Medial_Gender_Betas[(Medial_IndicesForOverlappingGenes_All[,1]==T),]
Medial_HoursFinal_Betas_OverlappingGenes<-Medial_HoursFinal_Betas[(Medial_IndicesForOverlappingGenes_All[,1]==T),]
Medial_MDD_Betas_OverlappingGenes<-Medial_MDD_Betas[(Medial_IndicesForOverlappingGenes_All[,1]==T),]
Medial_PH_Betas_OverlappingGenes<-Medial_PH_Betas[(Medial_IndicesForOverlappingGenes_All[,1]==T),]

PAC_Age_Betas_OverlappingGenes<-PAC_Age_Betas[(PAC_IndicesForOverlappingGenes_All[,1]==T),]
PAC_RNAIntegrity_Betas_OverlappingGenes<-PAC_RNAIntegrity_Betas[(PAC_IndicesForOverlappingGenes_All[,1]==T),]
PAC_Gender_Betas_OverlappingGenes<-PAC_Gender_Betas[(PAC_IndicesForOverlappingGenes_All[,1]==T),]
PAC_HoursFinal_Betas_OverlappingGenes<-PAC_HoursFinal_Betas[(PAC_IndicesForOverlappingGenes_All[,1]==T),]
PAC_MDD_Betas_OverlappingGenes<-PAC_MDD_Betas[(PAC_IndicesForOverlappingGenes_All[,1]==T),]
PAC_PH_Betas_OverlappingGenes<-PAC_PH_Betas[(PAC_IndicesForOverlappingGenes_All[,1]==T),]


PL_Age_Betas_OverlappingGenes<-PL_Age_Betas[(PL_IndicesForOverlappingGenes_All[,1]==T),]
PL_RNAIntegrity_Betas_OverlappingGenes<-PL_RNAIntegrity_Betas[(PL_IndicesForOverlappingGenes_All[,1]==T),]
PL_Gender_Betas_OverlappingGenes<-PL_Gender_Betas[(PL_IndicesForOverlappingGenes_All[,1]==T),]
PL_HoursFinal_Betas_OverlappingGenes<-PL_HoursFinal_Betas[(PL_IndicesForOverlappingGenes_All[,1]==T),]
PL_MDD_Betas_OverlappingGenes<-PL_MDD_Betas[(PL_IndicesForOverlappingGenes_All[,1]==T),]
PL_PH_Betas_OverlappingGenes<-PL_PH_Betas[(PL_IndicesForOverlappingGenes_All[,1]==T),]


#UPDATE: Change nrow to match # of overlapping probes:******************************************



AllRegions_NormNoOutliers_OverlappingGenes<-cbind(
AB_NormNoOutliers_OverlappingGenes,
Central_NormNoOutliers_OverlappingGenes,
CO_NormNoOutliers_OverlappingGenes,
Lateral_NormNoOutliers_OverlappingGenes,
PAC_NormNoOutliers_OverlappingGenes,
AAA_NormNoOutliers_OverlappingGenes,
AHA_NormNoOutliers_OverlappingGenes,
Basal_NormNoOutliers_OverlappingGenes,
Medial_NormNoOutliers_OverlappingGenes,
PL_NormNoOutliers_OverlappingGenes)

write.csv(AllRegions_NormNoOutliers_OverlappingGenes, "AllRegions_NormNoOutliers_OverlappingGenes.csv")

heatmap(cor(AllRegions_NormNoOutliers_OverlappingGenes))


AllRegions_NormNoOutliers_OverlappingGenesv2<-cbind(row.names(AllRegions_NormNoOutliers_OverlappingGenes), AllRegions_NormNoOutliers_OverlappingGenes)
colnames(AllRegions_NormNoOutliers_OverlappingGenesv2)[1]<-"X"

AllRegions_NormNoOutliers_OverlappingGenesv3<- join(as.data.frame(AllRegions_NormNoOutliers_OverlappingGenesv2), ProbeInfoOutput, by = "X")
write.csv(AllRegions_NormNoOutliers_OverlappingGenesv3, "AllRegions_NormNoOutliers_OverlappingGenesv3.csv")





Age_Pvalues_OverlappingGenes<-matrix(c(

AAA_Age_Pvalues_OverlappingGenes,
AB_Age_Pvalues_OverlappingGenes,
AHA_Age_Pvalues_OverlappingGenes,
Basal_Age_Pvalues_OverlappingGenes,
Central_Age_Pvalues_OverlappingGenes,
CO_Age_Pvalues_OverlappingGenes,
Lateral_Age_Pvalues_OverlappingGenes,
Medial_Age_Pvalues_OverlappingGenes,
PAC_Age_Pvalues_OverlappingGenes,
PL_Age_Pvalues_OverlappingGenes),
nrow=15788, ncol=10)

row.names(Age_Pvalues_OverlappingGenes)<-names(AAA_Age_Pvalues_OverlappingGenes)
colnames(Age_Pvalues_OverlappingGenes)<-c("AAA", "AB", "AHA", "Basal", "Central", "CO", "Lateral", "Medial", "PAC", "PL")

Age_Pvalues_OverlappingGenes[c(1:3), c(1:10)]


RNAIntegrity_Pvalues_OverlappingGenes<-matrix(c(

AAA_RNAIntegrity_Pvalues_OverlappingGenes,
AB_RNAIntegrity_Pvalues_OverlappingGenes,
AHA_RNAIntegrity_Pvalues_OverlappingGenes,
Basal_RNAIntegrity_Pvalues_OverlappingGenes,
Central_RNAIntegrity_Pvalues_OverlappingGenes,
CO_RNAIntegrity_Pvalues_OverlappingGenes,
Lateral_RNAIntegrity_Pvalues_OverlappingGenes,
Medial_RNAIntegrity_Pvalues_OverlappingGenes,
PAC_RNAIntegrity_Pvalues_OverlappingGenes,
PL_RNAIntegrity_Pvalues_OverlappingGenes),
nrow=15788, ncol=10)

row.names(RNAIntegrity_Pvalues_OverlappingGenes)<-names(AAA_RNAIntegrity_Pvalues_OverlappingGenes)
colnames(RNAIntegrity_Pvalues_OverlappingGenes)<-c("AAA", "AB", "AHA", "Basal", "Central", "CO", "Lateral", "Medial", "PAC", "PL")

RNAIntegrity_Pvalues_OverlappingGenes[c(1:3), c(1:10)]



Gender_Pvalues_OverlappingGenes<-matrix(c(

AAA_Gender_Pvalues_OverlappingGenes,
AB_Gender_Pvalues_OverlappingGenes,
AHA_Gender_Pvalues_OverlappingGenes,
Basal_Gender_Pvalues_OverlappingGenes,
Central_Gender_Pvalues_OverlappingGenes,
CO_Gender_Pvalues_OverlappingGenes,
Lateral_Gender_Pvalues_OverlappingGenes,
Medial_Gender_Pvalues_OverlappingGenes,
PAC_Gender_Pvalues_OverlappingGenes,
PL_Gender_Pvalues_OverlappingGenes),
nrow=15788, ncol=10)

row.names(Gender_Pvalues_OverlappingGenes)<-names(AAA_Gender_Pvalues_OverlappingGenes)
colnames(Gender_Pvalues_OverlappingGenes)<-c("AAA", "AB", "AHA", "Basal", "Central", "CO", "Lateral", "Medial", "PAC", "PL")

Gender_Pvalues_OverlappingGenes[c(1:3), c(1:10)]



HoursFinal_Pvalues_OverlappingGenes<-matrix(c(

AAA_HoursFinal_Pvalues_OverlappingGenes,
AB_HoursFinal_Pvalues_OverlappingGenes,
AHA_HoursFinal_Pvalues_OverlappingGenes,
Basal_HoursFinal_Pvalues_OverlappingGenes,
Central_HoursFinal_Pvalues_OverlappingGenes,
CO_HoursFinal_Pvalues_OverlappingGenes,
Lateral_HoursFinal_Pvalues_OverlappingGenes,
Medial_HoursFinal_Pvalues_OverlappingGenes,
PAC_HoursFinal_Pvalues_OverlappingGenes,
PL_HoursFinal_Pvalues_OverlappingGenes),
nrow=15788, ncol=10)

row.names(HoursFinal_Pvalues_OverlappingGenes)<-names(AAA_HoursFinal_Pvalues_OverlappingGenes)
colnames(HoursFinal_Pvalues_OverlappingGenes)<-c("AAA", "AB", "AHA", "Basal", "Central", "CO", "Lateral", "Medial", "PAC", "PL")

HoursFinal_Pvalues_OverlappingGenes[c(1:3), c(1:10)]



MDD_Pvalues_OverlappingGenes<-matrix(c(

AAA_MDD_Pvalues_OverlappingGenes,
AB_MDD_Pvalues_OverlappingGenes,
AHA_MDD_Pvalues_OverlappingGenes,
Basal_MDD_Pvalues_OverlappingGenes,
Central_MDD_Pvalues_OverlappingGenes,
CO_MDD_Pvalues_OverlappingGenes,
Lateral_MDD_Pvalues_OverlappingGenes,
Medial_MDD_Pvalues_OverlappingGenes,
PAC_MDD_Pvalues_OverlappingGenes,
PL_MDD_Pvalues_OverlappingGenes),
nrow=15788, ncol=10)

row.names(MDD_Pvalues_OverlappingGenes)<-names(AAA_MDD_Pvalues_OverlappingGenes)
colnames(MDD_Pvalues_OverlappingGenes)<-c("AAA", "AB", "AHA", "Basal", "Central", "CO", "Lateral", "Medial", "PAC", "PL")

MDD_Pvalues_OverlappingGenes[c(1:3), c(1:10)]



PH_Pvalues_OverlappingGenes<-matrix(c(

AAA_PH_Pvalues_OverlappingGenes,
AB_PH_Pvalues_OverlappingGenes,
AHA_PH_Pvalues_OverlappingGenes,
Basal_PH_Pvalues_OverlappingGenes,
Central_PH_Pvalues_OverlappingGenes,
CO_PH_Pvalues_OverlappingGenes,
Lateral_PH_Pvalues_OverlappingGenes,
Medial_PH_Pvalues_OverlappingGenes,
PAC_PH_Pvalues_OverlappingGenes,
PL_PH_Pvalues_OverlappingGenes),
nrow=15788, ncol=10)

row.names(PH_Pvalues_OverlappingGenes)<-names(AAA_PH_Pvalues_OverlappingGenes)
colnames(PH_Pvalues_OverlappingGenes)<-c("AAA", "AB", "AHA", "Basal", "Central", "CO", "Lateral", "Medial", "PAC", "PL")

PH_Pvalues_OverlappingGenes[c(1:3), c(1:10)]



Age_Betas_OverlappingGenes<-matrix(c(

AAA_Age_Betas_OverlappingGenes,
AB_Age_Betas_OverlappingGenes,
AHA_Age_Betas_OverlappingGenes,
Basal_Age_Betas_OverlappingGenes,
Central_Age_Betas_OverlappingGenes,
CO_Age_Betas_OverlappingGenes,
Lateral_Age_Betas_OverlappingGenes,
Medial_Age_Betas_OverlappingGenes,
PAC_Age_Betas_OverlappingGenes,
PL_Age_Betas_OverlappingGenes),
nrow=15788, ncol=10)

row.names(Age_Betas_OverlappingGenes)<-names(AAA_Age_Betas_OverlappingGenes)
colnames(Age_Betas_OverlappingGenes)<-c("AAA", "AB", "AHA", "Basal", "Central", "CO", "Lateral", "Medial", "PAC", "PL")

Age_Betas_OverlappingGenes[c(1:3), c(1:10)]


RNAIntegrity_Betas_OverlappingGenes<-matrix(c(

AAA_RNAIntegrity_Betas_OverlappingGenes,
AB_RNAIntegrity_Betas_OverlappingGenes,
AHA_RNAIntegrity_Betas_OverlappingGenes,
Basal_RNAIntegrity_Betas_OverlappingGenes,
Central_RNAIntegrity_Betas_OverlappingGenes,
CO_RNAIntegrity_Betas_OverlappingGenes,
Lateral_RNAIntegrity_Betas_OverlappingGenes,
Medial_RNAIntegrity_Betas_OverlappingGenes,
PAC_RNAIntegrity_Betas_OverlappingGenes,
PL_RNAIntegrity_Betas_OverlappingGenes),
nrow=15788, ncol=10)

row.names(RNAIntegrity_Betas_OverlappingGenes)<-names(AAA_RNAIntegrity_Betas_OverlappingGenes)
colnames(RNAIntegrity_Betas_OverlappingGenes)<-c("AAA", "AB", "AHA", "Basal", "Central", "CO", "Lateral", "Medial", "PAC", "PL")

RNAIntegrity_Betas_OverlappingGenes[c(1:3), c(1:10)]



Gender_Betas_OverlappingGenes<-matrix(c(

AAA_Gender_Betas_OverlappingGenes,
AB_Gender_Betas_OverlappingGenes,
AHA_Gender_Betas_OverlappingGenes,
Basal_Gender_Betas_OverlappingGenes,
Central_Gender_Betas_OverlappingGenes,
CO_Gender_Betas_OverlappingGenes,
Lateral_Gender_Betas_OverlappingGenes,
Medial_Gender_Betas_OverlappingGenes,
PAC_Gender_Betas_OverlappingGenes,
PL_Gender_Betas_OverlappingGenes),
nrow=15788, ncol=10)


row.names(Gender_Betas_OverlappingGenes)<-names(AAA_Gender_Betas_OverlappingGenes)
colnames(Gender_Betas_OverlappingGenes)<-c("AAA", "AB", "AHA", "Basal", "Central", "CO", "Lateral", "Medial", "PAC", "PL")

Gender_Betas_OverlappingGenes[c(1:3), c(1:10)]



HoursFinal_Betas_OverlappingGenes<-matrix(c(

AAA_HoursFinal_Betas_OverlappingGenes,
AB_HoursFinal_Betas_OverlappingGenes,
AHA_HoursFinal_Betas_OverlappingGenes,
Basal_HoursFinal_Betas_OverlappingGenes,
Central_HoursFinal_Betas_OverlappingGenes,
CO_HoursFinal_Betas_OverlappingGenes,
Lateral_HoursFinal_Betas_OverlappingGenes,
Medial_HoursFinal_Betas_OverlappingGenes,
PAC_HoursFinal_Betas_OverlappingGenes,
PL_HoursFinal_Betas_OverlappingGenes),
nrow=15788, ncol=10)

row.names(HoursFinal_Betas_OverlappingGenes)<-names(AAA_HoursFinal_Betas_OverlappingGenes)
colnames(HoursFinal_Betas_OverlappingGenes)<-c("AAA", "AB", "AHA", "Basal", "Central", "CO", "Lateral", "Medial", "PAC", "PL")

HoursFinal_Betas_OverlappingGenes[c(1:3), c(1:10)]



MDD_Betas_OverlappingGenes<-matrix(c(

AAA_MDD_Betas_OverlappingGenes,
AB_MDD_Betas_OverlappingGenes,
AHA_MDD_Betas_OverlappingGenes,
Basal_MDD_Betas_OverlappingGenes,
Central_MDD_Betas_OverlappingGenes,
CO_MDD_Betas_OverlappingGenes,
Lateral_MDD_Betas_OverlappingGenes,
Medial_MDD_Betas_OverlappingGenes,
PAC_MDD_Betas_OverlappingGenes,
PL_MDD_Betas_OverlappingGenes),
nrow=15788, ncol=10)

row.names(MDD_Betas_OverlappingGenes)<-names(AAA_MDD_Betas_OverlappingGenes)
colnames(MDD_Betas_OverlappingGenes)<-c("AAA", "AB", "AHA", "Basal", "Central", "CO", "Lateral", "Medial", "PAC", "PL")

MDD_Betas_OverlappingGenes[c(1:3), c(1:10)]



PH_Betas_OverlappingGenes<-matrix(c(

AAA_PH_Betas_OverlappingGenes,
AB_PH_Betas_OverlappingGenes,
AHA_PH_Betas_OverlappingGenes,
Basal_PH_Betas_OverlappingGenes,
Central_PH_Betas_OverlappingGenes,
CO_PH_Betas_OverlappingGenes,
Lateral_PH_Betas_OverlappingGenes,
Medial_PH_Betas_OverlappingGenes,
PAC_PH_Betas_OverlappingGenes,
PL_PH_Betas_OverlappingGenes),
nrow=15788, ncol=10)

row.names(PH_Betas_OverlappingGenes)<-names(AAA_PH_Betas_OverlappingGenes)
colnames(PH_Betas_OverlappingGenes)<-c("AAA", "AB", "AHA", "Basal", "Central", "CO", "Lateral", "Medial", "PAC", "PL")

PH_Betas_OverlappingGenes[c(1:3), c(1:10)]




#Change the working directory:

png("Age_Betas_OverlappingGenes.png")
heatmap(cor(Age_Betas_OverlappingGenes))
dev.off()

png("RNAIntegrity_Betas_OverlappingGenes.png")
heatmap(cor(RNAIntegrity_Betas_OverlappingGenes))
dev.off()

png("Gender_Betas_OverlappingGenes.png")
heatmap(cor(Gender_Betas_OverlappingGenes))
dev.off()

png("HoursFinal_Betas_OverlappingGenes.png")
heatmap(cor(HoursFinal_Betas_OverlappingGenes))
dev.off()

png("MDD_Betas_OverlappingGenes.png")
heatmap(cor(MDD_Betas_OverlappingGenes))
dev.off()

png("PH_Betas_OverlappingGenes.png")
heatmap(cor(PH_Betas_OverlappingGenes))
dev.off()



temp<-cor(MDD_Betas_OverlappingGenes)
write.csv(temp, "MDDBetasRegionCorrelations.csv")

temp<-cor(Age_Betas_OverlappingGenes)
write.csv(temp, "AgeBetasRegionCorrelations.csv")

temp<-cor(Gender_Betas_OverlappingGenes)
write.csv(temp, "GenderBetasRegionCorrelations.csv")

temp<-cor(PH_Betas_OverlappingGenes)
write.csv(temp, "PHBetasRegionCorrelations.csv")

temp<-cor(HoursFinal_Betas_OverlappingGenes)
write.csv(temp, "HoursFinalBetasRegionCorrelations.csv")

temp<-cor(RNAIntegrity_Betas_OverlappingGenes)
write.csv(temp, "RNAIntegrityBetasRegionCorrelations.csv")


#Code for running PCA to do Regional Comparisons of the Betas for Various Effects

pcaMDD<-prcomp(t(MDD_Betas_OverlappingGenes))
tmp<-pcaMDD$x[,1:10]
rownames(tmp)<-colnames(MDD_Betas_OverlappingGenes)
write.csv(tmp, "MDD_PCA.csv")

tmp<-pcaMDD$rotation[,1:10]
write.csv(tmp, "MDD_PCA_Eigenvectors.csv")

png("MDD_PC1vsPC2.png")
plot(pcaMDD$x[,1]~pcaMDD$x[,2], xlab="PC2", ylab="PC1", main="PC1 vs. PC2")
dev.off()

png("MDD_PCA_ScreePlot.png")
plot(summary(pcaMDD)$importance[2,]~(c(1:length(summary(pcaMDD)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC#", ylab="Proportion of Variance Explained", col=2)
dev.off()

png("MDD_PCA_ScreePlot2.png")
plot(summary(pcaMDD)$importance[3,]~(c(1:length(summary(pcaMDD)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC#", ylab="Cumulative Proportion of Variance Explained", col=2)
dev.off()


pcaGender<-prcomp(t(Gender_Betas_OverlappingGenes))
tmp<-pcaGender$x[,1:10]
rownames(tmp)<-colnames(Gender_Betas_OverlappingGenes)
write.csv(tmp, "Gender_PCA.csv")

tmp<-pcaGender$rotation[,1:10]
write.csv(tmp, "Gender_PCA_Eigenvectors.csv")

png("Gender_PC1vsPC2.png")
plot(pcaGender$x[,1]~pcaGender$x[,2], xlab="PC2", ylab="PC1", main="PC1 vs. PC2")
dev.off()

png("Gender_PCA_ScreePlot.png")
plot(summary(pcaGender)$importance[2,]~(c(1:length(summary(pcaGender)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC#", ylab="Proportion of Variance Explained", col=2)
dev.off()

png("Gender_PCA_ScreePlot2.png")
plot(summary(pcaGender)$importance[3,]~(c(1:length(summary(pcaGender)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC#", ylab="Cumulative Proportion of Variance Explained", col=2)
dev.off()


pcaAge<-prcomp(t(Age_Betas_OverlappingGenes))
tmp<-pcaAge$x[,1:10]
rownames(tmp)<-colnames(Age_Betas_OverlappingGenes)
write.csv(tmp, "Age_PCA.csv")

tmp<-pcaAge$rotation[,1:10]
write.csv(tmp, "Age_PCA_Eigenvectors.csv")

png("Age_PC1vsPC2.png")
plot(pcaAge$x[,1]~pcaAge$x[,2], xlab="PC2", ylab="PC1", main="PC1 vs. PC2")
dev.off()

png("Age_PCA_ScreePlot.png")
plot(summary(pcaAge)$importance[2,]~(c(1:length(summary(pcaAge)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC#", ylab="Proportion of Variance Explained", col=2)
dev.off()

png("Age_PCA_ScreePlot2.png")
plot(summary(pcaAge)$importance[3,]~(c(1:length(summary(pcaAge)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC#", ylab="Cumulative Proportion of Variance Explained", col=2)
dev.off()


pcaPH<-prcomp(t(PH_Betas_OverlappingGenes))
tmp<-pcaPH$x[,1:10]
rownames(tmp)<-colnames(PH_Betas_OverlappingGenes)
write.csv(tmp, "PH_PCA.csv")

tmp<-pcaPH$rotation[,1:10]
write.csv(tmp, "PH_PCA_Eigenvectors.csv")

png("PH_PC1vsPC2.png")
plot(pcaPH$x[,1]~pcaPH$x[,2], xlab="PC2", ylab="PC1", main="PC1 vs. PC2")
dev.off()

png("PH_PCA_ScreePlot.png")
plot(summary(pcaPH)$importance[2,]~(c(1:length(summary(pcaPH)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC#", ylab="Proportion of Variance Explained", col=2)
dev.off()

png("PH_PCA_ScreePlot2.png")
plot(summary(pcaPH)$importance[3,]~(c(1:length(summary(pcaPH)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC#", ylab="Cumulative Proportion of Variance Explained", col=2)
dev.off()


pcaHoursFinal<-prcomp(t(HoursFinal_Betas_OverlappingGenes))
tmp<-pcaHoursFinal$x[,1:10]
rownames(tmp)<-colnames(HoursFinal_Betas_OverlappingGenes)
write.csv(tmp, "HoursFinal_PCA.csv")

tmp<-pcaHoursFinal$rotation[,1:10]
write.csv(tmp, "HoursFinal_PCA_Eigenvectors.csv")

png("HoursFinal_PC1vsPC2.png")
plot(pcaHoursFinal$x[,1]~pcaHoursFinal$x[,2], xlab="PC2", ylab="PC1", main="PC1 vs. PC2")
dev.off()

png("HoursFinal_PCA_ScreePlot.png")
plot(summary(pcaHoursFinal)$importance[2,]~(c(1:length(summary(pcaHoursFinal)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC#", ylab="Proportion of Variance Explained", col=2)
dev.off()

png("HoursFinal_PCA_ScreePlot2.png")
plot(summary(pcaHoursFinal)$importance[3,]~(c(1:length(summary(pcaHoursFinal)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC#", ylab="Cumulative Proportion of Variance Explained", col=2)
dev.off()



pcaRNAIntegrity<-prcomp(t(RNAIntegrity_Betas_OverlappingGenes))
tmp<-pcaRNAIntegrity$x[,1:10]
rownames(tmp)<-colnames(RNAIntegrity_Betas_OverlappingGenes)
write.csv(tmp, "RNAIntegrity_PCA.csv")

tmp<-pcaRNAIntegrity$rotation[,1:10]
write.csv(tmp, "RNAIntegrity_PCA_Eigenvectors.csv")

png("RNAIntegrity_PC1vsPC2.png")
plot(pcaRNAIntegrity$x[,1]~pcaRNAIntegrity$x[,2], xlab="PC2", ylab="PC1", main="PC1 vs. PC2")
dev.off()

png("RNAIntegrity_PCA_ScreePlot.png")
plot(summary(pcaRNAIntegrity)$importance[2,]~(c(1:length(summary(pcaRNAIntegrity)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC#", ylab="Proportion of Variance Explained", col=2)
dev.off()

png("RNAIntegrity_PCA_ScreePlot2.png")
plot(summary(pcaRNAIntegrity)$importance[3,]~(c(1:length(summary(pcaRNAIntegrity)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC#", ylab="Cumulative Proportion of Variance Explained", col=2)
dev.off()


#Consensus Clustering

#Normally we median center by gene before running the clustering. Does that actually make sense in this case?  (since the betas are already centered around 0 for the most part?)
MDD_Betas_OverlappingGenesMedianCentered<-sweep(MDD_Betas_OverlappingGenes, 1, apply(MDD_Betas_OverlappingGenes, 1, median, na.rm=T), FUN="-")

ClusterMDDAttempt1<-ConsensusClusterPlus(MDD_Betas_OverlappingGenesMedianCentered, maxK=6, reps=20, clusterAlg="km", title="ConsensusCluster_MDDAttempt1", distance="euclidean", plot="png", writeTable=TRUE)
ClusterMDDAttempt2<-ConsensusClusterPlus(MDD_Betas_OverlappingGenes, maxK=6, reps=20, clusterAlg="km", title="ConsensusCluster_MDDAttempt2", distance="euclidean", plot="png", writeTable=TRUE)

#Both Jun and Wilkerson suggest choosing the top most variable genes (variable betas?)

MDD_Betas_OverlappingGenesStDev<-apply(MDD_Betas_OverlappingGenesMedianCentered, 1, sd)

png("Histogram_MDD_Betas_StDev.png")
hist(MDD_Betas_OverlappingGenesStDev)
dev.off()

quantile(MDD_Betas_OverlappingGenesStDev, probs=0.75)
# 75%  0.2374042 

MDD_Betas_OverlappingGenesMedianCenterMostVar<-MDD_Betas_OverlappingGenesMedianCentered[(MDD_Betas_OverlappingGenesStDev>quantile(MDD_Betas_OverlappingGenesStDev, probs=0.75)),]
dim(MDD_Betas_OverlappingGenesMedianCenterMostVar)
#[1] 3947   10


ClusterMDDAttempt3<-ConsensusClusterPlus(MDD_Betas_OverlappingGenesMedianCenterMostVar, maxK=6, reps=20, clusterAlg="km", title="ConsensusCluster_MDDAttempt2", distance="euclidean", plot="png", writeTable=TRUE)
#It seems to me like narrowing the analysis to just the top most variable genes simply produces a loss of information


#This is just me trying to figure out how I could make a plot of all of the regional beta correlations that is colored to match the magnitude of the correlation.
#So far I have totally failed - I think I will just come back to this problem later. :(

image(cor(MDD_Betas_OverlappingGenesMedianCenterMostVar))
cor(MDD_Betas_OverlappingGenesMedianCenterMostVar)[1,]

png("MDD_Betas_1_CorrelationSummaryMostVar.png", width = 1200, height = 1200)
pairs(MDD_Betas_OverlappingGenesMedianCenterMostVar,
	gap=0,
	par(col=rainbow(100)[round(((0.3+1)/2)*100, digits=0)]),
	diag.panel=function(x,...){
		par(new=TRUE)
		hist(x,
			col="light blue",
			probability = TRUE,
			axes=FALSE,
			main="")
		lines(density(x), 
			col="red",
			lwd=3)
		rug(x)
		})
dev.off()

(rainbow(100)[c(1:100)])


val2col<-function(z, zlim, col = heat.colors(12), breaks){
 if(!missing(breaks)){
  if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
 }
 if(missing(breaks) & !missing(zlim)){
  zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
  zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
  breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1)) 
 }
 if(missing(breaks) & missing(zlim)){
  zlim <- range(z, na.rm=TRUE)
  zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
  zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
  breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
 }
 CUT <- cut(z, breaks=breaks)
 colorlevels <- col[match(CUT, levels(CUT))] # assign colors to heights for each point
 return(colorlevels)
}

val2col(cor(MDD_Betas_OverlappingGenesMedianCenterMostVar)[1,], zlim=1, col=heat.colors(10), breaks=11)

CUT <- cut(cor(MDD_Betas_OverlappingGenesMedianCenterMostVar)[1,], breaks=9)

rainbow(100)[round(((0.3+1)/2)*100, digits=0)]





# Average correlation between betas:

Average_Correlation_Betas_AllRegions<-(cor(Age_Betas_OverlappingGenes)+cor(RNAIntegrity_Betas_OverlappingGenes)+cor(Gender_Betas_OverlappingGenes)+cor(HoursFinal_Betas_OverlappingGenes)+cor(MDD_Betas_OverlappingGenes)+cor(PH_Betas_OverlappingGenes))/6

png("Average_Correlation_Betas_AllRegions.png")
heatmap(Average_Correlation_Betas_AllRegions)
dev.off()


png("Age_Pvalues_OverlappingGenes.png")
heatmap(cor(log(Age_Pvalues_OverlappingGenes)))
dev.off()

png("RNAIntegrity_Pvalues_OverlappingGenes.png")
heatmap(cor(log(RNAIntegrity_Pvalues_OverlappingGenes)))
dev.off()

png("Gender_Pvalues_OverlappingGenes.png")
heatmap(cor(log(Gender_Pvalues_OverlappingGenes)))
dev.off()

png("HoursFinal_Pvalues_OverlappingGenes.png")
heatmap(cor(log(HoursFinal_Pvalues_OverlappingGenes)))
dev.off()

png("MDD_Pvalues_OverlappingGenes.png")
heatmap(cor(log(MDD_Pvalues_OverlappingGenes)))
dev.off()

png("PH_Pvalues_OverlappingGenes.png")
heatmap(cor(log(PH_Pvalues_OverlappingGenes)))
dev.off()


#Examples of correlations between individual brain regions:

colnames(MDD_Betas_OverlappingGenes)

#Two Regions that supposedly correlate well:

png("MDD_Betas_AHAvsPL.png")
plot(MDD_Betas_OverlappingGenes[,10]~MDD_Betas_OverlappingGenes[,3])
temp<-lm(MDD_Betas_OverlappingGenes[,10]~MDD_Betas_OverlappingGenes[,3])
abline(temp, col=2)
mtext(paste("R-Squared=", round(summary.lm(temp)$r.squared, digits=3), sep=" "))
dev.off()

png("MDD_Betas_PACvsBasal.png")
plot(MDD_Betas_OverlappingGenes[,9]~MDD_Betas_OverlappingGenes[,4])
temp<-lm(MDD_Betas_OverlappingGenes[,9]~MDD_Betas_OverlappingGenes[,4])
abline(temp, col=2)
mtext(paste("R-Squared=", round(summary.lm(temp)$r.squared, digits=3), sep=" "))
dev.off()


png("MDD_Betas_CentralvsCO.png")
plot(MDD_Betas_OverlappingGenes[,5]~MDD_Betas_OverlappingGenes[,6])
temp<-lm(MDD_Betas_OverlappingGenes[,5]~MDD_Betas_OverlappingGenes[,6])
abline(temp, col=2)
mtext(paste("R-Squared=", round(summary.lm(temp)$r.squared, digits=3), sep=" "))
dev.off()

png("MDD_Betas_AAAvsLateral.png")
plot(MDD_Betas_OverlappingGenes[,1]~MDD_Betas_OverlappingGenes[,7])
temp<-lm(MDD_Betas_OverlappingGenes[,1]~MDD_Betas_OverlappingGenes[,7])
abline(temp, col=2)
mtext(paste("R-Squared=", round(summary.lm(temp)$r.squared, digits=3), sep=" "))
dev.off()


#Making  nicer, high res figures for the paper:

pdf("MDD_Betas_BasalVsPAC.pdf", height=5.5, width=5)
plot(MDD_Betas_OverlappingGenes[,4]~MDD_Betas_OverlappingGenes[,9], cex.lab=1.4, cex.axis=1, font.lab=2, ylab="Basal", xlab="PAC")
temp<-lm(MDD_Betas_OverlappingGenes[,4]~MDD_Betas_OverlappingGenes[,9])
abline(temp, col=2, lwd=4)
dev.off()

pdf("MDD_Betas_BasalVsCentral.pdf", height=5.5, width=5)
plot(MDD_Betas_OverlappingGenes[,4]~MDD_Betas_OverlappingGenes[,5], cex.lab=1.4, cex.axis=1, font.lab=2, ylab="Basal", xlab="Central")
temp<-lm(MDD_Betas_OverlappingGenes[,4]~MDD_Betas_OverlappingGenes[,5])
abline(temp, col=2, lwd=4)
dev.off()

pdf("MDD_Betas_COVsCentral.pdf", height=5.5, width=5)
plot(MDD_Betas_OverlappingGenes[,6]~MDD_Betas_OverlappingGenes[,5], cex.lab=1.4, cex.axis=1, font.lab=2, ylab="CO", xlab="Central")
temp<-lm(MDD_Betas_OverlappingGenes[,6]~MDD_Betas_OverlappingGenes[,5])
abline(temp, col=2, lwd=4)
dev.off()


#Two Regions that supposedly correlate badly:

png("MDD_Betas_AHAvsMedial.png")
plot(MDD_Betas_OverlappingGenes[,8]~MDD_Betas_OverlappingGenes[,3])
temp<-lm(MDD_Betas_OverlappingGenes[,8]~MDD_Betas_OverlappingGenes[,3])
abline(temp, col=2)
mtext(paste("R-Squared=", round(summary.lm(temp)$r.squared, digits=3), sep=" "))
dev.off()

png("MDD_Betas_PACvsLateral.png")
plot(MDD_Betas_OverlappingGenes[,9]~MDD_Betas_OverlappingGenes[,7])
temp<-lm(MDD_Betas_OverlappingGenes[,9]~MDD_Betas_OverlappingGenes[,7])
abline(temp, col=2)
mtext(paste("R-Squared=", round(summary.lm(temp)$r.squared, digits=3), sep=" "))
dev.off()

png("MDD_Betas_PACvsCO.png")
plot(MDD_Betas_OverlappingGenes[,9]~MDD_Betas_OverlappingGenes[,6])
temp<-lm(MDD_Betas_OverlappingGenes[,9]~MDD_Betas_OverlappingGenes[,6])
abline(temp, col=2)
mtext(paste("R-Squared=", round(summary.lm(temp)$r.squared, digits=3), sep=" "))
dev.off()


png("MDD_Betas_AAAvsMedial.png")
plot(MDD_Betas_OverlappingGenes[,1]~MDD_Betas_OverlappingGenes[,8])
temp<-lm(MDD_Betas_OverlappingGenes[,1]~MDD_Betas_OverlappingGenes[,8])
abline(temp, col=2)
mtext(paste("R-Squared=", round(summary.lm(temp)$r.squared, digits=3), sep=" "))
dev.off()


png("MDD_Betas_ABvsMedial.png")
plot(MDD_Betas_OverlappingGenes[,2]~MDD_Betas_OverlappingGenes[,8])
temp<-lm(MDD_Betas_OverlappingGenes[,2]~MDD_Betas_OverlappingGenes[,8])
abline(temp, col=2)
mtext(paste("R-Squared=", round(summary.lm(temp)$r.squared, digits=3), sep=" "))
dev.off()




AmygdalaRegions<-c("AAA", "AB","AHA", "Basal", "Central", "CO", "Lateral", "Medial", "PAC", "PL") 

for(i in 1:10){
for(j in 1:10){
png(paste("MDD_Betas_", AmygdalaRegions[i], "vs", AmygdalaRegions[j], ".png", sep=""))
plot(MDD_Betas_OverlappingGenes[,i]~MDD_Betas_OverlappingGenes[,j])
temp<-lm(MDD_Betas_OverlappingGenes[,i]~MDD_Betas_OverlappingGenes[,j])
abline(temp, col=2)
mtext(paste("R-Squared=", round(summary.lm(temp)$r.squared, digits=3), sep=" "))
dev.off()
}
}

png("MDD_Betas_1_CorrelationSummary.png", width = 1200, height = 1200)
pairs(MDD_Betas_OverlappingGenes,
	gap=0,
	diag.panel=function(x,...){
		par(new=TRUE)
		hist(x,
			col="light blue",
			probability = TRUE,
			axes=FALSE,
			main="")
		lines(density(x), 
			col="red",
			lwd=3)
		rug(x)
		})
dev.off()




for(i in 1:10){
for(j in 1:10){
png(paste("Gender_Betas_", AmygdalaRegions[i], "vs", AmygdalaRegions[j], ".png", sep=""))
plot(Gender_Betas_OverlappingGenes[,i]~Gender_Betas_OverlappingGenes[,j])
temp<-lm(Gender_Betas_OverlappingGenes[,i]~Gender_Betas_OverlappingGenes[,j])
abline(temp, col=2)
mtext(paste("R-Squared=", round(summary.lm(temp)$r.squared, digits=3), sep=" "))
dev.off()
}
}


png("Gender_Betas_1_CorrelationSummary.png", width = 1200, height = 1200)
pairs(Gender_Betas_OverlappingGenes,
	gap=0,
	diag.panel=function(x,...){
		par(new=TRUE)
		hist(x,
			col="light blue",
			probability = TRUE,
			axes=FALSE,
			main="")
		lines(density(x), 
			col="red",
			lwd=3)
		rug(x)
		})
dev.off()






for(i in 1:10){
for(j in 1:10){
png(paste("Age_Betas_", AmygdalaRegions[i], "vs", AmygdalaRegions[j], ".png", sep=""))
plot(Age_Betas_OverlappingGenes[,i]~Age_Betas_OverlappingGenes[,j])
temp<-lm(Age_Betas_OverlappingGenes[,i]~Age_Betas_OverlappingGenes[,j])
abline(temp, col=2)
mtext(paste("R-Squared=", round(summary.lm(temp)$r.squared, digits=3), sep=" "))
dev.off()
}
}



png("Age_Betas_1_CorrelationSummary.png", width = 1200, height = 1200)
pairs(Age_Betas_OverlappingGenes,
	gap=0,
	diag.panel=function(x,...){
		par(new=TRUE)
		hist(x,
			col="light blue",
			probability = TRUE,
			axes=FALSE,
			main="")
		lines(density(x), 
			col="red",
			lwd=3)
		rug(x)
		})
dev.off()





for(i in 1:10){
for(j in 1:10){
png(paste("PH_Betas_", AmygdalaRegions[i], "vs", AmygdalaRegions[j], ".png", sep=""))
plot(PH_Betas_OverlappingGenes[,i]~PH_Betas_OverlappingGenes[,j])
temp<-lm(PH_Betas_OverlappingGenes[,i]~PH_Betas_OverlappingGenes[,j])
abline(temp, col=2)
mtext(paste("R-Squared=", round(summary.lm(temp)$r.squared, digits=3), sep=" "))
dev.off()
}
}


png("PH_Betas_1_CorrelationSummary.png", width = 1200, height = 1200)
pairs(PH_Betas_OverlappingGenes,
	gap=0,
	diag.panel=function(x,...){
		par(new=TRUE)
		hist(x,
			col="light blue",
			probability = TRUE,
			axes=FALSE,
			main="")
		lines(density(x), 
			col="red",
			lwd=3)
		rug(x)
		})
dev.off()





for(i in 1:10){
for(j in 1:10){
png(paste("HoursFinal_Betas_", AmygdalaRegions[i], "vs", AmygdalaRegions[j], ".png", sep=""))
plot(HoursFinal_Betas_OverlappingGenes[,i]~HoursFinal_Betas_OverlappingGenes[,j])
temp<-lm(HoursFinal_Betas_OverlappingGenes[,i]~HoursFinal_Betas_OverlappingGenes[,j])
abline(temp, col=2)
mtext(paste("R-Squared=", round(summary.lm(temp)$r.squared, digits=3), sep=" "))
dev.off()
}
}


png("HoursFinal_Betas_1_CorrelationSummary.png", width = 1200, height = 1200)
pairs(HoursFinal_Betas_OverlappingGenes,
	gap=0,
	diag.panel=function(x,...){
		par(new=TRUE)
		hist(x,
			col="light blue",
			probability = TRUE,
			axes=FALSE,
			main="")
		lines(density(x), 
			col="red",
			lwd=3)
		rug(x)
		})
dev.off()




for(i in 1:10){
for(j in 1:10){
png(paste("RNAIntegrity_Betas_", AmygdalaRegions[i], "vs", AmygdalaRegions[j], ".png", sep=""))
plot(RNAIntegrity_Betas_OverlappingGenes[,i]~RNAIntegrity_Betas_OverlappingGenes[,j])
temp<-lm(RNAIntegrity_Betas_OverlappingGenes[,i]~RNAIntegrity_Betas_OverlappingGenes[,j])
abline(temp, col=2)
mtext(paste("R-Squared=", round(summary.lm(temp)$r.squared, digits=3), sep=" "))
dev.off()
}
}

png("RNAIntegrity_Betas_1_CorrelationSummary.png", width = 1200, height = 1200)
pairs(RNAIntegrity_Betas_OverlappingGenes,
	gap=0,
	diag.panel=function(x,...){
		par(new=TRUE)
		hist(x,
			col="light blue",
			probability = TRUE,
			axes=FALSE,
			main="")
		lines(density(x), 
			col="red",
			lwd=3)
		rug(x)
		})
dev.off()











#Calculating Fisher's p-value across brain regions for each variable:

#Change working directory to Annotation

ProbeInfoOutput<-read.csv("ProbeInfoOutput.csv")
colnames(ProbeInfoOutput)


#Change working directory back to regional comparison:

Age_Summary_CrossRegionalResults<-cbind(
row.names(Age_Pvalues_OverlappingGenes),
Age_Pvalues_OverlappingGenes[, c(1:10)], 
apply(Age_Pvalues_OverlappingGenes, 1, median),
apply(log(Age_Pvalues_OverlappingGenes), 1, mean),
apply(log(Age_Pvalues_OverlappingGenes), 1, function(x) -2*sum(x)))

Age_Summary_CrossRegionalResults<-cbind(
Age_Summary_CrossRegionalResults,
pchisq(as.numeric(Age_Summary_CrossRegionalResults[,14]), 20, lower.tail=FALSE),
Age_Betas_OverlappingGenes[, c(1:10)],
apply(Age_Betas_OverlappingGenes, 1, median),
apply(Age_Betas_OverlappingGenes, 1, mean)
)

Age_Summary_CrossRegionalResults<-cbind(
Age_Summary_CrossRegionalResults,
apply(Age_Betas_OverlappingGenes, 1, function(x) sd(x)/sqrt(10))
)

colnames(Age_Summary_CrossRegionalResults)[1]<-"X"
colnames(Age_Summary_CrossRegionalResults)[c(12:15)]<-c("Median Pvalue", "Mean Log Pvalue", "Chi-Square Stat", "Fishers Pvalue")
colnames(Age_Summary_CrossRegionalResults)[c(26:28)]<-c("Median Beta", "Mean Beta", "SE Beta")

Age_Summary_CrossRegionalResults[c(1:3), c(1:27)]

Age_Summary_CrossRegionalResultsv2<- join(as.data.frame(Age_Summary_CrossRegionalResults), ProbeInfoOutput, by = "X")
write.csv(Age_Summary_CrossRegionalResultsv2, "Age_Summary_CrossRegionalResults.csv")



RNAIntegrity_Summary_CrossRegionalResults<-cbind(
row.names(RNAIntegrity_Pvalues_OverlappingGenes),
RNAIntegrity_Pvalues_OverlappingGenes[, c(1:10)], 
apply(RNAIntegrity_Pvalues_OverlappingGenes, 1, median),
apply(log(RNAIntegrity_Pvalues_OverlappingGenes), 1, mean),
apply(log(RNAIntegrity_Pvalues_OverlappingGenes), 1, function(x) -2*sum(x)))

RNAIntegrity_Summary_CrossRegionalResults<-cbind(
RNAIntegrity_Summary_CrossRegionalResults,
pchisq(as.numeric(RNAIntegrity_Summary_CrossRegionalResults[,14]), 20, lower.tail=FALSE),
RNAIntegrity_Betas_OverlappingGenes[, c(1:10)],
apply(RNAIntegrity_Betas_OverlappingGenes, 1, median),
apply(RNAIntegrity_Betas_OverlappingGenes, 1, mean)
)

RNAIntegrity_Summary_CrossRegionalResults<-cbind(
RNAIntegrity_Summary_CrossRegionalResults,
apply(RNAIntegrity_Betas_OverlappingGenes, 1, function(x) sd(x)/sqrt(10))
)

colnames(RNAIntegrity_Summary_CrossRegionalResults)[1]<-"X"
colnames(RNAIntegrity_Summary_CrossRegionalResults)[c(12:15)]<-c("Median Pvalue", "Mean Log Pvalue", "Chi-Square Stat", "Fishers Pvalue")
colnames(RNAIntegrity_Summary_CrossRegionalResults)[c(26:28)]<-c("Median Beta", "Mean Beta", "SE Beta")

RNAIntegrity_Summary_CrossRegionalResults[c(1:3), c(1:27)]

RNAIntegrity_Summary_CrossRegionalResultsv2<- join(as.data.frame(RNAIntegrity_Summary_CrossRegionalResults), ProbeInfoOutput, by = "X")
write.csv(RNAIntegrity_Summary_CrossRegionalResultsv2, "RNAIntegrity_Summary_CrossRegionalResults.csv")





Gender_Summary_CrossRegionalResults<-cbind(
row.names(Gender_Pvalues_OverlappingGenes),
Gender_Pvalues_OverlappingGenes[, c(1:10)], 
apply(Gender_Pvalues_OverlappingGenes, 1, median),
apply(log(Gender_Pvalues_OverlappingGenes), 1, mean),
apply(log(Gender_Pvalues_OverlappingGenes), 1, function(x) -2*sum(x)))

Gender_Summary_CrossRegionalResults<-cbind(
Gender_Summary_CrossRegionalResults,
pchisq(as.numeric(Gender_Summary_CrossRegionalResults[,14]), 20, lower.tail=FALSE),
Gender_Betas_OverlappingGenes[, c(1:10)],
apply(Gender_Betas_OverlappingGenes, 1, median),
apply(Gender_Betas_OverlappingGenes, 1, mean)
)

Gender_Summary_CrossRegionalResults<-cbind(
Gender_Summary_CrossRegionalResults,
apply(Gender_Betas_OverlappingGenes, 1, function(x) sd(x)/sqrt(10))
)

colnames(Gender_Summary_CrossRegionalResults)[1]<-"X"
colnames(Gender_Summary_CrossRegionalResults)[c(12:15)]<-c("Median Pvalue", "Mean Log Pvalue", "Chi-Square Stat", "Fishers Pvalue")
colnames(Gender_Summary_CrossRegionalResults)[c(26:28)]<-c("Median Beta", "Mean Beta", "SE Beta")

Gender_Summary_CrossRegionalResults[c(1:3), c(1:27)]

Gender_Summary_CrossRegionalResultsv2<- join(as.data.frame(Gender_Summary_CrossRegionalResults), ProbeInfoOutput, by = "X")
write.csv(Gender_Summary_CrossRegionalResultsv2, "Gender_Summary_CrossRegionalResults.csv")


HoursFinal_Summary_CrossRegionalResults<-cbind(
row.names(HoursFinal_Pvalues_OverlappingGenes),
HoursFinal_Pvalues_OverlappingGenes[, c(1:10)], 
apply(HoursFinal_Pvalues_OverlappingGenes, 1, median),
apply(log(HoursFinal_Pvalues_OverlappingGenes), 1, mean),
apply(log(HoursFinal_Pvalues_OverlappingGenes), 1, function(x) -2*sum(x)))

HoursFinal_Summary_CrossRegionalResults<-cbind(
HoursFinal_Summary_CrossRegionalResults,
pchisq(as.numeric(HoursFinal_Summary_CrossRegionalResults[,14]), 20, lower.tail=FALSE),
HoursFinal_Betas_OverlappingGenes[, c(1:10)],
apply(HoursFinal_Betas_OverlappingGenes, 1, median),
apply(HoursFinal_Betas_OverlappingGenes, 1, mean)
)

HoursFinal_Summary_CrossRegionalResults<-cbind(
HoursFinal_Summary_CrossRegionalResults,
apply(HoursFinal_Betas_OverlappingGenes, 1, function(x) sd(x)/sqrt(10))
)

colnames(HoursFinal_Summary_CrossRegionalResults)[1]<-"X"
colnames(HoursFinal_Summary_CrossRegionalResults)[c(12:15)]<-c("Median Pvalue", "Mean Log Pvalue", "Chi-Square Stat", "Fishers Pvalue")
colnames(HoursFinal_Summary_CrossRegionalResults)[c(26:28)]<-c("Median Beta", "Mean Beta", "SE Beta")

HoursFinal_Summary_CrossRegionalResults[c(1:3), c(1:27)]

HoursFinal_Summary_CrossRegionalResultsv2<- join(as.data.frame(HoursFinal_Summary_CrossRegionalResults), ProbeInfoOutput, by = "X")
write.csv(HoursFinal_Summary_CrossRegionalResultsv2, "HoursFinal_Summary_CrossRegionalResults.csv")


MDD_Summary_CrossRegionalResults<-cbind(
row.names(MDD_Pvalues_OverlappingGenes),
MDD_Pvalues_OverlappingGenes[, c(1:10)], 
apply(MDD_Pvalues_OverlappingGenes, 1, median),
apply(log(MDD_Pvalues_OverlappingGenes), 1, mean),
apply(log(MDD_Pvalues_OverlappingGenes), 1, function(x) -2*sum(x)))

MDD_Summary_CrossRegionalResults<-cbind(
MDD_Summary_CrossRegionalResults,
pchisq(as.numeric(MDD_Summary_CrossRegionalResults[,14]), 20, lower.tail=FALSE),
MDD_Betas_OverlappingGenes[, c(1:10)],
apply(MDD_Betas_OverlappingGenes, 1, median),
apply(MDD_Betas_OverlappingGenes, 1, mean)
)

MDD_Summary_CrossRegionalResults<-cbind(
MDD_Summary_CrossRegionalResults,
apply(MDD_Betas_OverlappingGenes, 1, function(x) sd(x)/sqrt(10))
)

colnames(MDD_Summary_CrossRegionalResults)[1]<-"X"
colnames(MDD_Summary_CrossRegionalResults)[c(12:15)]<-c("Median Pvalue", "Mean Log Pvalue", "Chi-Square Stat", "Fishers Pvalue")
colnames(MDD_Summary_CrossRegionalResults)[c(26:28)]<-c("Median Beta", "Mean Beta", "SE Beta")

MDD_Summary_CrossRegionalResults[c(1:3), c(1:27)]

MDD_Summary_CrossRegionalResultsv2<- join(as.data.frame(MDD_Summary_CrossRegionalResults), ProbeInfoOutput, by = "X")
write.csv(MDD_Summary_CrossRegionalResultsv2, "MDD_Summary_CrossRegionalResults.csv")


#Cluster1:  the AAA, AHA, Lateral, and PL together (medium PC1, medium PC2)

MDD_Summary_Cluster1<-cbind(
row.names(MDD_Pvalues_OverlappingGenes),
MDD_Pvalues_OverlappingGenes[, c(1, 3, 7, 10)], 
apply(MDD_Pvalues_OverlappingGenes[, c(1, 3, 7, 10)], 1, median),
apply(log(MDD_Pvalues_OverlappingGenes[, c(1, 3, 7, 10)]), 1, mean),
apply(log(MDD_Pvalues_OverlappingGenes[, c(1, 3, 7, 10)]), 1, function(x) -2*sum(x)))

MDD_Summary_Cluster1<-cbind(
MDD_Summary_Cluster1,
pchisq(as.numeric(MDD_Summary_Cluster1[,8]), 8, lower.tail=FALSE),
MDD_Betas_OverlappingGenes[, c(1, 3, 7, 10)],
apply(MDD_Betas_OverlappingGenes[, c(1, 3, 7, 10)], 1, median),
apply(MDD_Betas_OverlappingGenes[, c(1, 3, 7, 10)], 1, mean),
apply(MDD_Betas_OverlappingGenes[, c(1, 3, 7, 10)], 1, function(x) sd(x)/sqrt(10))
)

colnames(MDD_Summary_Cluster1)[1]<-"X"
colnames(MDD_Summary_Cluster1)[c(6:9)]<-c("Median Pvalue", "Mean Log Pvalue", "Chi-Square Stat", "Fishers Pvalue")
colnames(MDD_Summary_Cluster1)[c(14:16)]<-c("Median Beta", "Mean Beta", "SE Beta")

MDD_Summary_Cluster1v2<- join(as.data.frame(MDD_Summary_Cluster1), ProbeInfoOutput, by = "X")
write.csv(MDD_Summary_Cluster1v2, "MDD_Summary_Cluster1v2.csv")


png("Histogram Cluster1 Fisher Pvals.png")
hist(as.numeric(MDD_Summary_Cluster1[,9]), breaks=100, col=2)
abline(a=158, b=0)
dev.off()


#Cluster2: just the AB (high PC1, high PC2) - correlates with both the Basal/PAC cluster and the medial amygdala

#Cluster3: the Basal and PAC (high PC1, low PC2)

MDD_Summary_Cluster3<-cbind(
row.names(MDD_Pvalues_OverlappingGenes),
MDD_Pvalues_OverlappingGenes[, c(4,9)], 
apply(MDD_Pvalues_OverlappingGenes[, c(4,9)], 1, median),
apply(log(MDD_Pvalues_OverlappingGenes[, c(4,9)]), 1, mean),
apply(log(MDD_Pvalues_OverlappingGenes[, c(4,9)]), 1, function(x) -2*sum(x)))

MDD_Summary_Cluster3<-cbind(
MDD_Summary_Cluster3,
pchisq(as.numeric(MDD_Summary_Cluster3[,6]), 4, lower.tail=FALSE),
MDD_Betas_OverlappingGenes[, c(4,9)],
apply(MDD_Betas_OverlappingGenes[, c(4,9)], 1, median),
apply(MDD_Betas_OverlappingGenes[, c(4,9)], 1, mean),
apply(MDD_Betas_OverlappingGenes[, c(4,9)], 1, function(x) sd(x)/sqrt(10))
)

colnames(MDD_Summary_Cluster3)[1]<-"X"
colnames(MDD_Summary_Cluster3)[c(4:7)]<-c("Median Pvalue", "Mean Log Pvalue", "Chi-Square Stat", "Fishers Pvalue")
colnames(MDD_Summary_Cluster3)[c(10:12)]<-c("Median Beta", "Mean Beta", "SE Beta")

MDD_Summary_Cluster3v2<- join(as.data.frame(MDD_Summary_Cluster3), ProbeInfoOutput, by = "X")
write.csv(MDD_Summary_Cluster3v2, "MDD_Summary_Cluster3v2.csv")

png("Histogram Cluster3 Fisher Pvals.png")
hist(as.numeric(MDD_Summary_Cluster3[,7]), breaks=100, col=2)
abline(a=158, b=0)
dev.off()



Cluster4: the Central and CO amygdala (low PC1, medium-high PC2)

MDD_Summary_Cluster4<-cbind(
row.names(MDD_Pvalues_OverlappingGenes),
MDD_Pvalues_OverlappingGenes[, c(5,6)], 
apply(MDD_Pvalues_OverlappingGenes[, c(5,6)], 1, median),
apply(log(MDD_Pvalues_OverlappingGenes[, c(5,6)]), 1, mean),
apply(log(MDD_Pvalues_OverlappingGenes[, c(5,6)]), 1, function(x) -2*sum(x)))

MDD_Summary_Cluster4<-cbind(
MDD_Summary_Cluster4,
pchisq(as.numeric(MDD_Summary_Cluster4[,6]), 4, lower.tail=FALSE),
MDD_Betas_OverlappingGenes[, c(5,6)],
apply(MDD_Betas_OverlappingGenes[, c(5,6)], 1, median),
apply(MDD_Betas_OverlappingGenes[, c(5,6)], 1, mean),
apply(MDD_Betas_OverlappingGenes[, c(5,6)], 1, function(x) sd(x)/sqrt(10))
)

colnames(MDD_Summary_Cluster4)[1]<-"X"
colnames(MDD_Summary_Cluster4)[c(4:7)]<-c("Median Pvalue", "Mean Log Pvalue", "Chi-Square Stat", "Fishers Pvalue")
colnames(MDD_Summary_Cluster4)[c(10:12)]<-c("Median Beta", "Mean Beta", "SE Beta")

MDD_Summary_Cluster4v2<- join(as.data.frame(MDD_Summary_Cluster4), ProbeInfoOutput, by = "X")
write.csv(MDD_Summary_Cluster4v2, "MDD_Summary_Cluster4v2.csv")


png("Histogram Cluster4 Fisher Pvals.png")
hist(as.numeric(MDD_Summary_Cluster4[,7]), breaks=100, col=2)
abline(a=158, b=0)
dev.off()








PH_Summary_CrossRegionalResults<-cbind(
row.names(PH_Pvalues_OverlappingGenes),
PH_Pvalues_OverlappingGenes[, c(1:10)], 
apply(PH_Pvalues_OverlappingGenes, 1, median),
apply(log(PH_Pvalues_OverlappingGenes), 1, mean),
apply(log(PH_Pvalues_OverlappingGenes), 1, function(x) -2*sum(x)))

PH_Summary_CrossRegionalResults<-cbind(
PH_Summary_CrossRegionalResults,
pchisq(as.numeric(PH_Summary_CrossRegionalResults[,14]), 20, lower.tail=FALSE),
PH_Betas_OverlappingGenes[, c(1:10)],
apply(PH_Betas_OverlappingGenes, 1, median),
apply(PH_Betas_OverlappingGenes, 1, mean)
)

PH_Summary_CrossRegionalResults<-cbind(
PH_Summary_CrossRegionalResults,
apply(PH_Betas_OverlappingGenes, 1, function(x) sd(x)/sqrt(10))
)

colnames(PH_Summary_CrossRegionalResults)[1]<-"X"
colnames(PH_Summary_CrossRegionalResults)[c(12:15)]<-c("Median Pvalue", "Mean Log Pvalue", "Chi-Square Stat", "Fishers Pvalue")
colnames(PH_Summary_CrossRegionalResults)[c(26:28)]<-c("Median Beta", "Mean Beta", "SE Beta")

PH_Summary_CrossRegionalResults[c(1:3), c(1:27)]

PH_Summary_CrossRegionalResultsv2<- join(as.data.frame(PH_Summary_CrossRegionalResults), ProbeInfoOutput, by = "X")
write.csv(PH_Summary_CrossRegionalResultsv2, "PH_Summary_CrossRegionalResults.csv")



MDDpvaluesAllRegions<-matrix(0,10,5)
row.names(MDDpvaluesAllRegions)<-c("AAA", "AB", "AHA", "Basal", "Central", "CO", "Lateral", "Medial", "PAC", "PL")
colnames(MDDpvaluesAllRegions)<-c("pval<0.001", "Random: pval<0.001", "pval<0.01", "Random: pval<0.01", "BH pval<0.10")

MDDpvaluesAllRegions[1,]<-as.numeric(AAA_PvalueOutputSummary[1, c(2:6)])
MDDpvaluesAllRegions[2,]<-as.numeric(AB_PvalueOutputSummary[1, c(2:6)])
MDDpvaluesAllRegions[3,]<-as.numeric(AHA_PvalueOutputSummary[1, c(2:6)])
MDDpvaluesAllRegions[4,]<-as.numeric(Basal_PvalueOutputSummary[1, c(2:6)])
MDDpvaluesAllRegions[5,]<-as.numeric(Central_PvalueOutputSummary[1, c(2:6)])
MDDpvaluesAllRegions[6,]<-as.numeric(CO_PvalueOutputSummary[1, c(2:6)])
MDDpvaluesAllRegions[7,]<-as.numeric(Lateral_PvalueOutputSummary[1, c(2:6)])
MDDpvaluesAllRegions[8,]<-as.numeric(Medial_PvalueOutputSummary[1, c(2:6)])
MDDpvaluesAllRegions[9,]<-as.numeric(PAC_PvalueOutputSummary[1, c(2:6)])
MDDpvaluesAllRegions[10,]<-as.numeric(PL_PvalueOutputSummary[1, c(2:6)])

write.csv(MDDpvaluesAllRegions, "MDDpvaluesAllRegions.csv")

AgepvaluesAllRegions<-matrix(0,10,5)
row.names(AgepvaluesAllRegions)<-c("AAA", "AB", "AHA", "Basal", "Central", "CO", "Lateral", "Medial", "PAC", "PL")
colnames(AgepvaluesAllRegions)<-c("pval<0.001", "Random: pval<0.001", "pval<0.01", "Random: pval<0.01", "BH pval<0.10")

AgepvaluesAllRegions[1,]<-as.numeric(AAA_PvalueOutputSummary[2, c(2:6)])
AgepvaluesAllRegions[2,]<-as.numeric(AB_PvalueOutputSummary[2, c(2:6)])
AgepvaluesAllRegions[3,]<-as.numeric(AHA_PvalueOutputSummary[2, c(2:6)])
AgepvaluesAllRegions[4,]<-as.numeric(Basal_PvalueOutputSummary[2, c(2:6)])
AgepvaluesAllRegions[5,]<-as.numeric(Central_PvalueOutputSummary[2, c(2:6)])
AgepvaluesAllRegions[6,]<-as.numeric(CO_PvalueOutputSummary[2, c(2:6)])
AgepvaluesAllRegions[7,]<-as.numeric(Lateral_PvalueOutputSummary[2, c(2:6)])
AgepvaluesAllRegions[8,]<-as.numeric(Medial_PvalueOutputSummary[2, c(2:6)])
AgepvaluesAllRegions[9,]<-as.numeric(PAC_PvalueOutputSummary[2, c(2:6)])
AgepvaluesAllRegions[10,]<-as.numeric(PL_PvalueOutputSummary[2, c(2:6)])

write.csv(AgepvaluesAllRegions, "AgepvaluesAllRegions.csv")

RNAIntegritypvaluesAllRegions<-matrix(0,10,5)
row.names(RNAIntegritypvaluesAllRegions)<-c("AAA", "AB", "AHA", "Basal", "Central", "CO", "Lateral", "Medial", "PAC", "PL")
colnames(RNAIntegritypvaluesAllRegions)<-c("pval<0.001", "Random: pval<0.001", "pval<0.01", "Random: pval<0.01", "BH pval<0.10")

RNAIntegritypvaluesAllRegions[1,]<-as.numeric(AAA_PvalueOutputSummary[3, c(2:6)])
RNAIntegritypvaluesAllRegions[2,]<-as.numeric(AB_PvalueOutputSummary[3, c(2:6)])
RNAIntegritypvaluesAllRegions[3,]<-as.numeric(AHA_PvalueOutputSummary[3, c(2:6)])
RNAIntegritypvaluesAllRegions[4,]<-as.numeric(Basal_PvalueOutputSummary[3, c(2:6)])
RNAIntegritypvaluesAllRegions[5,]<-as.numeric(Central_PvalueOutputSummary[3, c(2:6)])
RNAIntegritypvaluesAllRegions[6,]<-as.numeric(CO_PvalueOutputSummary[3, c(2:6)])
RNAIntegritypvaluesAllRegions[7,]<-as.numeric(Lateral_PvalueOutputSummary[3, c(2:6)])
RNAIntegritypvaluesAllRegions[8,]<-as.numeric(Medial_PvalueOutputSummary[3, c(2:6)])
RNAIntegritypvaluesAllRegions[9,]<-as.numeric(PAC_PvalueOutputSummary[3, c(2:6)])
RNAIntegritypvaluesAllRegions[10,]<-as.numeric(PL_PvalueOutputSummary[3, c(2:6)])

write.csv(RNAIntegritypvaluesAllRegions, "RNAIntegritypvaluesAllRegions.csv")


BrainPHpvaluesAllRegions<-matrix(0,10,5)
row.names(BrainPHpvaluesAllRegions)<-c("AAA", "AB", "AHA", "Basal", "Central", "CO", "Lateral", "Medial", "PAC", "PL")
colnames(BrainPHpvaluesAllRegions)<-c("pval<0.001", "Random: pval<0.001", "pval<0.01", "Random: pval<0.01", "BH pval<0.10")

BrainPHpvaluesAllRegions[1,]<-as.numeric(AAA_PvalueOutputSummary[4, c(2:6)])
BrainPHpvaluesAllRegions[2,]<-as.numeric(AB_PvalueOutputSummary[4, c(2:6)])
BrainPHpvaluesAllRegions[3,]<-as.numeric(AHA_PvalueOutputSummary[4, c(2:6)])
BrainPHpvaluesAllRegions[4,]<-as.numeric(Basal_PvalueOutputSummary[4, c(2:6)])
BrainPHpvaluesAllRegions[5,]<-as.numeric(Central_PvalueOutputSummary[4, c(2:6)])
BrainPHpvaluesAllRegions[6,]<-as.numeric(CO_PvalueOutputSummary[4, c(2:6)])
BrainPHpvaluesAllRegions[7,]<-as.numeric(Lateral_PvalueOutputSummary[4, c(2:6)])
BrainPHpvaluesAllRegions[8,]<-as.numeric(Medial_PvalueOutputSummary[4, c(2:6)])
BrainPHpvaluesAllRegions[9,]<-as.numeric(PAC_PvalueOutputSummary[4, c(2:6)])
BrainPHpvaluesAllRegions[10,]<-as.numeric(PL_PvalueOutputSummary[4, c(2:6)])

write.csv(BrainPHpvaluesAllRegions, "BrainPHpvaluesAllRegions.csv")



GenderpvaluesAllRegions<-matrix(0,10,5)
row.names(GenderpvaluesAllRegions)<-c("AAA", "AB", "AHA", "Basal", "Central", "CO", "Lateral", "Medial", "PAC", "PL")
colnames(GenderpvaluesAllRegions)<-c("pval<0.001", "Random: pval<0.001", "pval<0.01", "Random: pval<0.01", "BH pval<0.10")

GenderpvaluesAllRegions[1,]<-as.numeric(AAA_PvalueOutputSummary[5, c(2:6)])
GenderpvaluesAllRegions[2,]<-as.numeric(AB_PvalueOutputSummary[5, c(2:6)])
GenderpvaluesAllRegions[3,]<-as.numeric(AHA_PvalueOutputSummary[5, c(2:6)])
GenderpvaluesAllRegions[4,]<-as.numeric(Basal_PvalueOutputSummary[5, c(2:6)])
GenderpvaluesAllRegions[5,]<-as.numeric(Central_PvalueOutputSummary[5, c(2:6)])
GenderpvaluesAllRegions[6,]<-as.numeric(CO_PvalueOutputSummary[5, c(2:6)])
GenderpvaluesAllRegions[7,]<-as.numeric(Lateral_PvalueOutputSummary[5, c(2:6)])
GenderpvaluesAllRegions[8,]<-as.numeric(Medial_PvalueOutputSummary[5, c(2:6)])
GenderpvaluesAllRegions[9,]<-as.numeric(PAC_PvalueOutputSummary[5, c(2:6)])
GenderpvaluesAllRegions[10,]<-as.numeric(PL_PvalueOutputSummary[5, c(2:6)])

write.csv(GenderpvaluesAllRegions, "GenderpvaluesAllRegions.csv")



HoursFinalpvaluesAllRegions<-matrix(0,10,5)
row.names(HoursFinalpvaluesAllRegions)<-c("AAA", "AB", "AHA", "Basal", "Central", "CO", "Lateral", "Medial", "PAC", "PL")
colnames(HoursFinalpvaluesAllRegions)<-c("pval<0.001", "Random: pval<0.001", "pval<0.01", "Random: pval<0.01", "BH pval<0.10")

HoursFinalpvaluesAllRegions[1,]<-as.numeric(AAA_PvalueOutputSummary[6, c(2:6)])
HoursFinalpvaluesAllRegions[2,]<-as.numeric(AB_PvalueOutputSummary[6, c(2:6)])
HoursFinalpvaluesAllRegions[3,]<-as.numeric(AHA_PvalueOutputSummary[6, c(2:6)])
HoursFinalpvaluesAllRegions[4,]<-as.numeric(Basal_PvalueOutputSummary[6, c(2:6)])
HoursFinalpvaluesAllRegions[5,]<-as.numeric(Central_PvalueOutputSummary[6, c(2:6)])
HoursFinalpvaluesAllRegions[6,]<-as.numeric(CO_PvalueOutputSummary[6, c(2:6)])
HoursFinalpvaluesAllRegions[7,]<-as.numeric(Lateral_PvalueOutputSummary[6, c(2:6)])
HoursFinalpvaluesAllRegions[8,]<-as.numeric(Medial_PvalueOutputSummary[6, c(2:6)])
HoursFinalpvaluesAllRegions[9,]<-as.numeric(PAC_PvalueOutputSummary[6, c(2:6)])
HoursFinalpvaluesAllRegions[10,]<-as.numeric(PL_PvalueOutputSummary[6, c(2:6)])

write.csv(HoursFinalpvaluesAllRegions, "HoursFinalpvaluesAllRegions.csv")


#######################################

#Code for a failed project to compare relative expression levels for genes across subnuclei - despite the fact that the subnuclei were in different microarray processing batches!!!

#1.A.2: Read in the data and assign it a brain region identifier.


ReadInData2<-function(filename, Region, LM_ID){
	
LM_Intercept<-as.matrix(read.table(filename, header=T, row.names=2, sep=","))
print("Dimensions:")
print(dim(LM_Intercept))
print("Colnames:")
print(colnames(LM_Intercept))
print("Example:")
print(LM_Intercept[c(1:3), c(1:5)])

#Double check that the probe name is the first column, unadjusted p-values is the second column, and the beta is the 5th column

TempIntercept<-as.matrix(as.numeric(LM_Intercept[,2]))
row.names(TempIntercept)<-row.names(LM_Intercept)
print("Intercepts are Numeric?")
print(is.numeric(TempIntercept))
print(TempIntercept[c(1:3),1])

TempInterceptRank<-as.matrix(as.numeric(LM_Intercept[,4]))
row.names(TempInterceptRank)<-row.names(LM_Intercept)
print("Intercept Ranks are Numeric?")
print(is.numeric(TempInterceptRank))
print(TempInterceptRank[c(1:3),1])

assign(paste(Region, LM_ID, "Intercepts", sep="_"), TempIntercept, envir=as.environment(1))
TempIntercept<<-TempIntercept

assign(paste(Region, LM_ID, "InterceptRanks", sep="_"), TempInterceptRank, envir=as.environment(1))
rm(TempInterceptRank)
}




#Apparently the columns and rows are labeled incorrectly in the LM5 intercept output, so I had to make a special read in data function for it:

ReadInData4<-function(filename, Region, LM_ID){
	
LM_Intercept<-as.matrix(read.table(filename, header=T, row.names=1, sep=","))
print("Dimensions:")
print(dim(LM_Intercept))
print("Colnames:")
print(colnames(LM_Intercept))
print("Example:")
print(LM_Intercept[c(1:3), c(1:5)])

#Double check that the probe name is the first column, unadjusted p-values is the second column, and the beta is the 5th column

TempIntercept2<-as.matrix(as.numeric(LM_Intercept[,1]))
row.names(TempIntercept2)<-row.names(TempIntercept)
print("Intercepts are Numeric?")
print(is.numeric(TempIntercept2))
print(TempIntercept2[c(1:3),1])

TempInterceptRank<-as.matrix(as.numeric(LM_Intercept[,3]))
row.names(TempInterceptRank)<-row.names(TempIntercept)
print("Intercept Ranks are Numeric?")
print(is.numeric(TempInterceptRank))
print(TempInterceptRank[c(1:3),1])

assign(paste(Region, LM_ID, "Intercepts", sep="_"), TempIntercept2, envir=as.environment(1))
rm(TempIntercept2)

assign(paste(Region, LM_ID, "InterceptRanks", sep="_"), TempInterceptRank, envir=as.environment(1))
rm(TempInterceptRank)
}





ReadInData3<-function(filename, Region, LM_ID, RankOrSimple){
	
LM_Residuals<-as.matrix(read.table(filename, header=T, row.names=1, sep=","))
print("Dimensions:")
print(dim(LM_Residuals))
print("Colnames:")
print(colnames(LM_Residuals))
print("Example:")
print(LM_Residuals[c(1:3), c(1:5)])

print("Intercepts are Numeric?")
print(is.numeric(LM_Residuals))
print(LM_Residuals[c(1:3),c(1:3)])
LM_Residuals<<-LM_Residuals

assign(paste(Region, LM_ID, "Residuals", RankOrSimple, sep="_"), LM_Residuals, envir=as.environment(1))

}


ReadInData5<-function(filename, Region, LM_ID, RankOrSimple){
	
LM_Residuals2<-as.matrix(read.table(filename, header=T, row.names=1, sep=","))
print("Dimensions:")
print(dim(LM_Residuals2))
colnames(LM_Residuals2)<-colnames(LM_Residuals)
row.names(LM_Residuals2)<-row.names(LM_Residuals)
print("Colnames:")
print(colnames(LM_Residuals2))
print("Example:")
print(LM_Residuals2[c(1:3), c(1:5)])

print("Intercepts are Numeric?")
print(is.numeric(LM_Residuals2))
print(LM_Residuals2[c(1:3),c(1:3)])

assign(paste(Region, LM_ID, "Residuals", RankOrSimple, sep="_"), LM_Residuals2, envir=as.environment(1))
rm(LM_Residuals2)

}



#Changed working directory to AAA

ReadInData2("LM4InterceptOutputv2.csv", "AAA", "LM4")
ReadInData4("LM5InterceptOutputv2.csv", "AAA", "LM5")

ReadInData3("LM4ResidualsCleanedSemiOptimal.csv", "AAA", "LM4", "Simple")
ReadInData3("LM4ResidualsRank.csv", "AAA", "LM4", "Ranked")

ReadInData5("LM5ResidualsCleanedMedian.csv", "AAA", "LM5", "Simple")
ReadInData5("LM5ResidualsRank.csv", "AAA", "LM5", "Ranked")



#Changed working directory to AB

ReadInData2("LM4InterceptOutputv2.csv", "AB", "LM4")
ReadInData4("LM5InterceptOutputv2.csv", "AB", "LM5")

ReadInData3("LM4ResidualsCleanedSemiOptimal.csv", "AB", "LM4", "Simple")
ReadInData3("LM4ResidualsRank.csv", "AB", "LM4", "Ranked")

ReadInData5("LM5ResidualsCleanedMedian.csv", "AB", "LM5", "Simple")
ReadInData5("LM5ResidualsRank.csv", "AB", "LM5", "Ranked")


#Changed working directory to AHA

ReadInData2("LM4InterceptOutputv2.csv", "AHA", "LM4")
ReadInData4("LM5InterceptOutputv2.csv", "AHA", "LM5")

ReadInData3("LM4ResidualsCleanedSemiOptimal.csv", "AHA", "LM4", "Simple")
ReadInData3("LM4ResidualsRank.csv", "AHA", "LM4", "Ranked")

ReadInData5("LM5ResidualsCleanedMedian.csv", "AHA", "LM5", "Simple")
ReadInData5("LM5ResidualsRank.csv", "AHA", "LM5", "Ranked")


#Changed working directory to Basal

ReadInData2("LM4InterceptOutputv2.csv", "Basal", "LM4")
ReadInData4("LM5InterceptOutputv2.csv", "Basal", "LM5")

ReadInData3("LM4ResidualsCleanedSemiOptimal.csv", "Basal", "LM4", "Simple")
ReadInData3("LM4ResidualsRank.csv", "Basal", "LM4", "Ranked")

ReadInData5("LM5ResidualsCleanedMedian.csv", "Basal", "LM5", "Simple")
ReadInData5("LM5ResidualsRank.csv", "Basal", "LM5", "Ranked")


#Changed working directory to Central

ReadInData2("LM4InterceptOutputv2.csv", "Central", "LM4")
ReadInData4("LM5InterceptOutputv2.csv", "Central", "LM5")

ReadInData3("LM4ResidualsCleanedSemiOptimal.csv", "Central", "LM4", "Simple")
ReadInData3("LM4ResidualsRank.csv", "Central", "LM4", "Ranked")

ReadInData5("LM5ResidualsCleanedMedian.csv", "Central", "LM5", "Simple")
ReadInData5("LM5ResidualsRank.csv", "Central", "LM5", "Ranked")


#Changed working directory to CO

ReadInData2("LM4InterceptOutputv2.csv", "CO", "LM4")
ReadInData4("LM5InterceptOutputv2.csv", "CO", "LM5")

ReadInData3("LM4ResidualsCleanedSemiOptimal.csv", "CO", "LM4", "Simple")
ReadInData3("LM4ResidualsRank.csv", "CO", "LM4", "Ranked")

ReadInData5("LM5ResidualsCleanedMedian.csv", "CO", "LM5", "Simple")
ReadInData5("LM5ResidualsRank.csv", "CO", "LM5", "Ranked")


#Changed working directory to Lateral

ReadInData2("LM4InterceptOutputv2.csv", "Lateral", "LM4")
ReadInData4("LM5InterceptOutputv2.csv", "Lateral", "LM5")

ReadInData3("LM4ResidualsCleanedSemiOptimal.csv", "Lateral", "LM4", "Simple")
ReadInData3("LM4ResidualsRank.csv", "Lateral", "LM4", "Ranked")

ReadInData5("LM5ResidualsCleanedMedian.csv", "Lateral", "LM5", "Simple")
ReadInData5("LM5ResidualsRank.csv", "Lateral", "LM5", "Ranked")


#Changed working directory to Medial

ReadInData2("LM4InterceptOutputv2.csv", "Medial", "LM4")
ReadInData4("LM5InterceptOutputv2.csv", "Medial", "LM5")

ReadInData3("LM4ResidualsCleanedSemiOptimal.csv", "Medial", "LM4", "Simple")
ReadInData3("LM4ResidualsRank.csv", "Medial", "LM4", "Ranked")

ReadInData5("LM5ResidualsCleanedMedian.csv", "Medial", "LM5", "Simple")
ReadInData5("LM5ResidualsRank.csv", "Medial", "LM5", "Ranked")


#Changed working directory to PAC

ReadInData2("LM4InterceptOutputv2.csv", "PAC", "LM4")
ReadInData4("LM5InterceptOutputv2.csv", "PAC", "LM5")

ReadInData3("LM4ResidualsCleanedSemiOptimal.csv", "PAC", "LM4", "Simple")
ReadInData3("LM4ResidualsRank.csv", "PAC", "LM4", "Ranked")

ReadInData5("LM5ResidualsCleanedMedian.csv", "PAC", "LM5", "Simple")
ReadInData5("LM5ResidualsRank.csv", "PAC", "LM5", "Ranked")


#Changed working directory to PL

ReadInData2("LM4InterceptOutputv2.csv", "PL", "LM4")
ReadInData4("LM5InterceptOutputv2.csv", "PL", "LM5")

ReadInData3("LM4ResidualsCleanedSemiOptimal.csv", "PL", "LM4", "Simple")
ReadInData3("LM4ResidualsRank.csv", "PL", "LM4", "Ranked")

ReadInData5("LM5ResidualsCleanedMedian.csv", "PL", "LM5", "Simple")
ReadInData5("LM5ResidualsRank.csv", "PL", "LM5", "Ranked")


#Narrowing things down to just the genes found in all 10 nuclei:

AAA_LM4_InterceptRanks_OverlappingGenes<- AAA_LM4_InterceptRanks[(AAA_IndicesForOverlappingGenes_All[,1]==T),]                     
AAA_LM4_Intercepts_OverlappingGenes<-AAA_LM4_Intercepts[(AAA_IndicesForOverlappingGenes_All[,1]==T),]                        
AAA_LM4_Residuals_Ranked_OverlappingGenes<-AAA_LM4_Residuals_Ranked[(AAA_IndicesForOverlappingGenes_All[,1]==T),]                      
AAA_LM4_Residuals_Simple_OverlappingGenes<-AAA_LM4_Residuals_Simple[(AAA_IndicesForOverlappingGenes_All[,1]==T),]                      
AAA_LM5_InterceptRanks_OverlappingGenes<-AAA_LM5_InterceptRanks[(AAA_IndicesForOverlappingGenes_All[,1]==T),]                        
AAA_LM5_Intercepts_OverlappingGenes<-AAA_LM5_Intercepts[(AAA_IndicesForOverlappingGenes_All[,1]==T),]                           
AAA_LM5_Residuals_Ranked_OverlappingGenes<-AAA_LM5_Residuals_Ranked[(AAA_IndicesForOverlappingGenes_All[,1]==T),]                     
AAA_LM5_Residuals_Simple_OverlappingGenes<-AAA_LM5_Residuals_Simple[(AAA_IndicesForOverlappingGenes_All[,1]==T),]    



AB_LM4_InterceptRanks_OverlappingGenes<- AB_LM4_InterceptRanks[(AB_IndicesForOverlappingGenes_All[,1]==T),]                     
AB_LM4_Intercepts_OverlappingGenes<-AB_LM4_Intercepts[(AB_IndicesForOverlappingGenes_All[,1]==T),]                        
AB_LM4_Residuals_Ranked_OverlappingGenes<-AB_LM4_Residuals_Ranked[(AB_IndicesForOverlappingGenes_All[,1]==T),]                      
AB_LM4_Residuals_Simple_OverlappingGenes<-AB_LM4_Residuals_Simple[(AB_IndicesForOverlappingGenes_All[,1]==T),]                      
AB_LM5_InterceptRanks_OverlappingGenes<-AB_LM5_InterceptRanks[(AB_IndicesForOverlappingGenes_All[,1]==T),]                        
AB_LM5_Intercepts_OverlappingGenes<-AB_LM5_Intercepts[(AB_IndicesForOverlappingGenes_All[,1]==T),]                           
AB_LM5_Residuals_Ranked_OverlappingGenes<-AB_LM5_Residuals_Ranked[(AB_IndicesForOverlappingGenes_All[,1]==T),]                     
AB_LM5_Residuals_Simple_OverlappingGenes<-AB_LM5_Residuals_Simple[(AB_IndicesForOverlappingGenes_All[,1]==T),]    


AHA_LM4_InterceptRanks_OverlappingGenes<- AHA_LM4_InterceptRanks[(AHA_IndicesForOverlappingGenes_All[,1]==T),]                     
AHA_LM4_Intercepts_OverlappingGenes<-AHA_LM4_Intercepts[(AHA_IndicesForOverlappingGenes_All[,1]==T),]                        
AHA_LM4_Residuals_Ranked_OverlappingGenes<-AHA_LM4_Residuals_Ranked[(AHA_IndicesForOverlappingGenes_All[,1]==T),]                      
AHA_LM4_Residuals_Simple_OverlappingGenes<-AHA_LM4_Residuals_Simple[(AHA_IndicesForOverlappingGenes_All[,1]==T),]                      
AHA_LM5_InterceptRanks_OverlappingGenes<-AHA_LM5_InterceptRanks[(AHA_IndicesForOverlappingGenes_All[,1]==T),]                        
AHA_LM5_Intercepts_OverlappingGenes<-AHA_LM5_Intercepts[(AHA_IndicesForOverlappingGenes_All[,1]==T),]                           
AHA_LM5_Residuals_Ranked_OverlappingGenes<-AHA_LM5_Residuals_Ranked[(AHA_IndicesForOverlappingGenes_All[,1]==T),]                     
AHA_LM5_Residuals_Simple_OverlappingGenes<-AHA_LM5_Residuals_Simple[(AHA_IndicesForOverlappingGenes_All[,1]==T),]    


Basal_LM4_InterceptRanks_OverlappingGenes<- Basal_LM4_InterceptRanks[(Basal_IndicesForOverlappingGenes_All[,1]==T),]                     
Basal_LM4_Intercepts_OverlappingGenes<-Basal_LM4_Intercepts[(Basal_IndicesForOverlappingGenes_All[,1]==T),]                        
Basal_LM4_Residuals_Ranked_OverlappingGenes<-Basal_LM4_Residuals_Ranked[(Basal_IndicesForOverlappingGenes_All[,1]==T),]                      
Basal_LM4_Residuals_Simple_OverlappingGenes<-Basal_LM4_Residuals_Simple[(Basal_IndicesForOverlappingGenes_All[,1]==T),]                      
Basal_LM5_InterceptRanks_OverlappingGenes<-Basal_LM5_InterceptRanks[(Basal_IndicesForOverlappingGenes_All[,1]==T),]                        
Basal_LM5_Intercepts_OverlappingGenes<-Basal_LM5_Intercepts[(Basal_IndicesForOverlappingGenes_All[,1]==T),]                           
Basal_LM5_Residuals_Ranked_OverlappingGenes<-Basal_LM5_Residuals_Ranked[(Basal_IndicesForOverlappingGenes_All[,1]==T),]                     
Basal_LM5_Residuals_Simple_OverlappingGenes<-Basal_LM5_Residuals_Simple[(Basal_IndicesForOverlappingGenes_All[,1]==T),]    

Central_LM4_InterceptRanks_OverlappingGenes<- Central_LM4_InterceptRanks[(Central_IndicesForOverlappingGenes_All[,1]==T),]                     
Central_LM4_Intercepts_OverlappingGenes<-Central_LM4_Intercepts[(Central_IndicesForOverlappingGenes_All[,1]==T),]                        
Central_LM4_Residuals_Ranked_OverlappingGenes<-Central_LM4_Residuals_Ranked[(Central_IndicesForOverlappingGenes_All[,1]==T),]                      
Central_LM4_Residuals_Simple_OverlappingGenes<-Central_LM4_Residuals_Simple[(Central_IndicesForOverlappingGenes_All[,1]==T),]                      
Central_LM5_InterceptRanks_OverlappingGenes<-Central_LM5_InterceptRanks[(Central_IndicesForOverlappingGenes_All[,1]==T),]                        
Central_LM5_Intercepts_OverlappingGenes<-Central_LM5_Intercepts[(Central_IndicesForOverlappingGenes_All[,1]==T),]                           
Central_LM5_Residuals_Ranked_OverlappingGenes<-Central_LM5_Residuals_Ranked[(Central_IndicesForOverlappingGenes_All[,1]==T),]                     
Central_LM5_Residuals_Simple_OverlappingGenes<-Central_LM5_Residuals_Simple[(Central_IndicesForOverlappingGenes_All[,1]==T),]    

CO_LM4_InterceptRanks_OverlappingGenes<- CO_LM4_InterceptRanks[(CO_IndicesForOverlappingGenes_All[,1]==T),]                     
CO_LM4_Intercepts_OverlappingGenes<-CO_LM4_Intercepts[(CO_IndicesForOverlappingGenes_All[,1]==T),]                        
CO_LM4_Residuals_Ranked_OverlappingGenes<-CO_LM4_Residuals_Ranked[(CO_IndicesForOverlappingGenes_All[,1]==T),]                      
CO_LM4_Residuals_Simple_OverlappingGenes<-CO_LM4_Residuals_Simple[(CO_IndicesForOverlappingGenes_All[,1]==T),]                      
CO_LM5_InterceptRanks_OverlappingGenes<-CO_LM5_InterceptRanks[(CO_IndicesForOverlappingGenes_All[,1]==T),]                        
CO_LM5_Intercepts_OverlappingGenes<-CO_LM5_Intercepts[(CO_IndicesForOverlappingGenes_All[,1]==T),]                           
CO_LM5_Residuals_Ranked_OverlappingGenes<-CO_LM5_Residuals_Ranked[(CO_IndicesForOverlappingGenes_All[,1]==T),]                     
CO_LM5_Residuals_Simple_OverlappingGenes<-CO_LM5_Residuals_Simple[(CO_IndicesForOverlappingGenes_All[,1]==T),]    


Lateral_LM4_InterceptRanks_OverlappingGenes<- Lateral_LM4_InterceptRanks[(Lateral_IndicesForOverlappingGenes_All[,1]==T),]                     
Lateral_LM4_Intercepts_OverlappingGenes<-Lateral_LM4_Intercepts[(Lateral_IndicesForOverlappingGenes_All[,1]==T),]                        
Lateral_LM4_Residuals_Ranked_OverlappingGenes<-Lateral_LM4_Residuals_Ranked[(Lateral_IndicesForOverlappingGenes_All[,1]==T),]                      
Lateral_LM4_Residuals_Simple_OverlappingGenes<-Lateral_LM4_Residuals_Simple[(Lateral_IndicesForOverlappingGenes_All[,1]==T),]                      
Lateral_LM5_InterceptRanks_OverlappingGenes<-Lateral_LM5_InterceptRanks[(Lateral_IndicesForOverlappingGenes_All[,1]==T),]                        
Lateral_LM5_Intercepts_OverlappingGenes<-Lateral_LM5_Intercepts[(Lateral_IndicesForOverlappingGenes_All[,1]==T),]                           
Lateral_LM5_Residuals_Ranked_OverlappingGenes<-Lateral_LM5_Residuals_Ranked[(Lateral_IndicesForOverlappingGenes_All[,1]==T),]                     
Lateral_LM5_Residuals_Simple_OverlappingGenes<-Lateral_LM5_Residuals_Simple[(Lateral_IndicesForOverlappingGenes_All[,1]==T),]    


Medial_LM4_InterceptRanks_OverlappingGenes<- Medial_LM4_InterceptRanks[(Medial_IndicesForOverlappingGenes_All[,1]==T),]                     
Medial_LM4_Intercepts_OverlappingGenes<-Medial_LM4_Intercepts[(Medial_IndicesForOverlappingGenes_All[,1]==T),]                        
Medial_LM4_Residuals_Ranked_OverlappingGenes<-Medial_LM4_Residuals_Ranked[(Medial_IndicesForOverlappingGenes_All[,1]==T),]                      
Medial_LM4_Residuals_Simple_OverlappingGenes<-Medial_LM4_Residuals_Simple[(Medial_IndicesForOverlappingGenes_All[,1]==T),]                      
Medial_LM5_InterceptRanks_OverlappingGenes<-Medial_LM5_InterceptRanks[(Medial_IndicesForOverlappingGenes_All[,1]==T),]                        
Medial_LM5_Intercepts_OverlappingGenes<-Medial_LM5_Intercepts[(Medial_IndicesForOverlappingGenes_All[,1]==T),]                           
Medial_LM5_Residuals_Ranked_OverlappingGenes<-Medial_LM5_Residuals_Ranked[(Medial_IndicesForOverlappingGenes_All[,1]==T),]                     
Medial_LM5_Residuals_Simple_OverlappingGenes<-Medial_LM5_Residuals_Simple[(Medial_IndicesForOverlappingGenes_All[,1]==T),]    


PAC_LM4_InterceptRanks_OverlappingGenes<- PAC_LM4_InterceptRanks[(PAC_IndicesForOverlappingGenes_All[,1]==T),]                     
PAC_LM4_Intercepts_OverlappingGenes<-PAC_LM4_Intercepts[(PAC_IndicesForOverlappingGenes_All[,1]==T),]                        
PAC_LM4_Residuals_Ranked_OverlappingGenes<-PAC_LM4_Residuals_Ranked[(PAC_IndicesForOverlappingGenes_All[,1]==T),]                      
PAC_LM4_Residuals_Simple_OverlappingGenes<-PAC_LM4_Residuals_Simple[(PAC_IndicesForOverlappingGenes_All[,1]==T),]                      
PAC_LM5_InterceptRanks_OverlappingGenes<-PAC_LM5_InterceptRanks[(PAC_IndicesForOverlappingGenes_All[,1]==T),]                        
PAC_LM5_Intercepts_OverlappingGenes<-PAC_LM5_Intercepts[(PAC_IndicesForOverlappingGenes_All[,1]==T),]                           
PAC_LM5_Residuals_Ranked_OverlappingGenes<-PAC_LM5_Residuals_Ranked[(PAC_IndicesForOverlappingGenes_All[,1]==T),]                     
PAC_LM5_Residuals_Simple_OverlappingGenes<-PAC_LM5_Residuals_Simple[(PAC_IndicesForOverlappingGenes_All[,1]==T),]    

PL_LM4_InterceptRanks_OverlappingGenes<- PL_LM4_InterceptRanks[(PL_IndicesForOverlappingGenes_All[,1]==T),]                     
PL_LM4_Intercepts_OverlappingGenes<-PL_LM4_Intercepts[(PL_IndicesForOverlappingGenes_All[,1]==T),]                        
PL_LM4_Residuals_Ranked_OverlappingGenes<-PL_LM4_Residuals_Ranked[(PL_IndicesForOverlappingGenes_All[,1]==T),]                      
PL_LM4_Residuals_Simple_OverlappingGenes<-PL_LM4_Residuals_Simple[(PL_IndicesForOverlappingGenes_All[,1]==T),]                      
PL_LM5_InterceptRanks_OverlappingGenes<-PL_LM5_InterceptRanks[(PL_IndicesForOverlappingGenes_All[,1]==T),]                        
PL_LM5_Intercepts_OverlappingGenes<-PL_LM5_Intercepts[(PL_IndicesForOverlappingGenes_All[,1]==T),]                           
PL_LM5_Residuals_Ranked_OverlappingGenes<-PL_LM5_Residuals_Ranked[(PL_IndicesForOverlappingGenes_All[,1]==T),]                     
PL_LM5_Residuals_Simple_OverlappingGenes<-PL_LM5_Residuals_Simple[(PL_IndicesForOverlappingGenes_All[,1]==T),]    

 



#Combining data together from different regions into a single matrix

LM4_InterceptRanks_OverlappingGenes<-matrix(c(

AAA_LM4_InterceptRanks_OverlappingGenes,
AB_LM4_InterceptRanks_OverlappingGenes,
AHA_LM4_InterceptRanks_OverlappingGenes,
Basal_LM4_InterceptRanks_OverlappingGenes,
Central_LM4_InterceptRanks_OverlappingGenes,
CO_LM4_InterceptRanks_OverlappingGenes,
Lateral_LM4_InterceptRanks_OverlappingGenes,
Medial_LM4_InterceptRanks_OverlappingGenes,
PAC_LM4_InterceptRanks_OverlappingGenes,
PL_LM4_InterceptRanks_OverlappingGenes),
nrow=15788, ncol=10)

row.names(LM4_InterceptRanks_OverlappingGenes)<-names(AAA_LM4_InterceptRanks_OverlappingGenes)
colnames(LM4_InterceptRanks_OverlappingGenes)<-c("AAA", "AB", "AHA", "Basal", "Central", "CO", "Lateral", "Medial", "PAC", "PL")


LM5_InterceptRanks_OverlappingGenes<-matrix(c(

AAA_LM5_InterceptRanks_OverlappingGenes,
AB_LM5_InterceptRanks_OverlappingGenes,
AHA_LM5_InterceptRanks_OverlappingGenes,
Basal_LM5_InterceptRanks_OverlappingGenes,
Central_LM5_InterceptRanks_OverlappingGenes,
CO_LM5_InterceptRanks_OverlappingGenes,
Lateral_LM5_InterceptRanks_OverlappingGenes,
Medial_LM5_InterceptRanks_OverlappingGenes,
PAC_LM5_InterceptRanks_OverlappingGenes,
PL_LM5_InterceptRanks_OverlappingGenes),
nrow=15788, ncol=10)

row.names(LM5_InterceptRanks_OverlappingGenes)<-names(AAA_LM5_InterceptRanks_OverlappingGenes)
colnames(LM5_InterceptRanks_OverlappingGenes)<-c("AAA", "AB", "AHA", "Basal", "Central", "CO", "Lateral", "Medial", "PAC", "PL")


LM4_Intercepts_OverlappingGenes<-matrix(c(

AAA_LM4_Intercepts_OverlappingGenes,
AB_LM4_Intercepts_OverlappingGenes,
AHA_LM4_Intercepts_OverlappingGenes,
Basal_LM4_Intercepts_OverlappingGenes,
Central_LM4_Intercepts_OverlappingGenes,
CO_LM4_Intercepts_OverlappingGenes,
Lateral_LM4_Intercepts_OverlappingGenes,
Medial_LM4_Intercepts_OverlappingGenes,
PAC_LM4_Intercepts_OverlappingGenes,
PL_LM4_Intercepts_OverlappingGenes),
nrow=15788, ncol=10)

row.names(LM4_Intercepts_OverlappingGenes)<-names(AAA_LM4_Intercepts_OverlappingGenes)
colnames(LM4_Intercepts_OverlappingGenes)<-c("AAA", "AB", "AHA", "Basal", "Central", "CO", "Lateral", "Medial", "PAC", "PL")



LM5_Intercepts_OverlappingGenes<-matrix(c(

AAA_LM5_Intercepts_OverlappingGenes,
AB_LM5_Intercepts_OverlappingGenes,
AHA_LM5_Intercepts_OverlappingGenes,
Basal_LM5_Intercepts_OverlappingGenes,
Central_LM5_Intercepts_OverlappingGenes,
CO_LM5_Intercepts_OverlappingGenes,
Lateral_LM5_Intercepts_OverlappingGenes,
Medial_LM5_Intercepts_OverlappingGenes,
PAC_LM5_Intercepts_OverlappingGenes,
PL_LM5_Intercepts_OverlappingGenes),
nrow=15788, ncol=10)

row.names(LM5_Intercepts_OverlappingGenes)<-names(AAA_LM5_Intercepts_OverlappingGenes)
colnames(LM5_Intercepts_OverlappingGenes)<-c("AAA", "AB", "AHA", "Basal", "Central", "CO", "Lateral", "Medial", "PAC", "PL")


plot(LM5_InterceptRanks_OverlappingGenes[,1]~LM4_InterceptRanks_OverlappingGenes[,1])



#Note - I'm looking over this code again in an effort to determine how best to progress (July 18 2014) and it looks like the ranking was performed after the linear models were run - i.e., the intercepts for each gene for each region were ranked in the output, and then I also ranked the residuals by subject... which thinking about it is a rather silly thing to do for residuals. I think what would make more sense would be to go back to NormNoOutliers for each region, convert the gene data to ranks for each subject, re-run the linear models, and then output the intercepts with the SE for the regression as a manner of determining differences between regions. We could also try running a large multilevel model that includes subject and region as a term.

#Code to remove all Residuals_Ranked files in order to clear up some working memory!
rm(list=ls(pattern="*Residuals_Ranked*"))
# *nice* so much more efficient than doing that by hand...





png("LM5_InterceptRanks_OverlappingGenes_CorrSummary.png", width = 1200, height = 1200)
pairs(LM5_InterceptRanks_OverlappingGenes,
	gap=0,
	diag.panel=function(x,...){
		par(new=TRUE)
		hist(x,
			col="light blue",
			probability = TRUE,
			axes=FALSE,
			main="")
		lines(density(x), 
			col="red",
			lwd=3)
		rug(x)
		})
dev.off()



png("LM4_InterceptRanks_OverlappingGenes_CorrSummary.png", width = 1200, height = 1200)
pairs(LM4_InterceptRanks_OverlappingGenes,
	gap=0,
	diag.panel=function(x,...){
		par(new=TRUE)
		hist(x,
			col="light blue",
			probability = TRUE,
			axes=FALSE,
			main="")
		lines(density(x), 
			col="red",
			lwd=3)
		rug(x)
		})
dev.off()



png("LM4_Intercepts_OverlappingGenes_CorrSummary.png", width = 1200, height = 1200)
pairs(LM4_Intercepts_OverlappingGenes,
	gap=0,
	diag.panel=function(x,...){
		par(new=TRUE)
		hist(x,
			col="light blue",
			probability = TRUE,
			axes=FALSE,
			main="")
		lines(density(x), 
			col="red",
			lwd=3)
		rug(x)
		})
dev.off()



png("LM5_Intercepts_OverlappingGenes_CorrSummary.png", width = 1200, height = 1200)
pairs(LM5_Intercepts_OverlappingGenes,
	gap=0,
	diag.panel=function(x,...){
		par(new=TRUE)
		hist(x,
			col="light blue",
			probability = TRUE,
			axes=FALSE,
			main="")
		lines(density(x), 
			col="red",
			lwd=3)
		rug(x)
		})
dev.off()


png("Heatmap_LM4_Intercepts_OverlappingGenes.png")
heatmap(cor(LM4_Intercepts_OverlappingGenes))
dev.off()

write.csv(cor(LM4_Intercepts_OverlappingGenes), "LM4Intercepts_CorrSummary.csv")


png("Heatmap_LM5_Intercepts_OverlappingGenes.png")
heatmap(cor(LM5_Intercepts_OverlappingGenes))
dev.off()

write.csv(cor(LM5_Intercepts_OverlappingGenes), "LM5Intercepts_CorrSummary.csv")


png("Heatmap_LM4_InterceptRanks_OverlappingGenes.png")
heatmap(cor(LM4_InterceptRanks_OverlappingGenes))
dev.off()

write.csv(cor(LM4_InterceptRanks_OverlappingGenes), "LM4InterceptRanks_CorrSummary.csv")


png("Heatmap_LM5_InterceptRanks_OverlappingGenes.png")
heatmap(cor(LM5_InterceptRanks_OverlappingGenes))
dev.off()

write.csv(cor(LM5_InterceptRanks_OverlappingGenes), "LM5InterceptRanks_CorrSummary.csv")



#O.k., it seems to me like the non-ranked data is unusuable becuase there are clearly differences in the general distribution of the signal values
#Also, it seems like median centering the data so as to better detect correlations between regions might be the way to go


LM5_InterceptRanks_OverlappingGenes_MedianCentered<-sweep(LM5_InterceptRanks_OverlappingGenes, 1, apply(LM5_InterceptRanks_OverlappingGenes, 1, median))
plot(sort(LM5_InterceptRanks_OverlappingGenes_MedianCentered))

LM4_InterceptRanks_OverlappingGenes_MedianCentered<-sweep(LM4_InterceptRanks_OverlappingGenes, 1, apply(LM4_InterceptRanks_OverlappingGenes, 1, median))
plot(sort(LM4_InterceptRanks_OverlappingGenes_MedianCentered))

LM5_Intercepts_OverlappingGenes_MedianCentered<-sweep(LM5_Intercepts_OverlappingGenes, 1, apply(LM5_Intercepts_OverlappingGenes, 1, median))
plot(sort(LM5_Intercepts_OverlappingGenes_MedianCentered))

LM4_Intercepts_OverlappingGenes_MedianCentered<-sweep(LM4_Intercepts_OverlappingGenes, 1, apply(LM4_Intercepts_OverlappingGenes, 1, median))
plot(sort(LM4_Intercepts_OverlappingGenes_MedianCentered))



write.csv(LM4_Intercepts_OverlappingGenes, "LM4Intercepts.csv")
write.csv(LM5_Intercepts_OverlappingGenes, "LM5Intercepts.csv")
write.csv(LM4_InterceptRanks_OverlappingGenes, "LM4InterceptRanks.csv")
write.csv(LM5_InterceptRanks_OverlappingGenes, "LM5InterceptRanks.csv")

LargeInterceptMatrix<-cbind(row.names(LM4_Intercepts_OverlappingGenes), LM4_Intercepts_OverlappingGenes, LM5_Intercepts_OverlappingGenes, LM4_InterceptRanks_OverlappingGenes, LM5_InterceptRanks_OverlappingGenes )

colnames(LargeInterceptMatrix)[1]<-"X"

LargeInterceptMatrix[c(1:3), c(1:3)]

LargeInterceptMatrixOutput<-join(as.data.frame(LargeInterceptMatrix), ProbeInfoOutput, by = "X")

write.csv(LargeInterceptMatrixOutput, "LargeInterceptMatrixOutput.csv")






png("Heatmap_LM4_Intercepts_OverlappingGenesCntrd.png")
heatmap(cor(LM4_Intercepts_OverlappingGenes_MedianCentered))
dev.off()

write.csv(cor(LM4_Intercepts_OverlappingGenes_MedianCentered), "LM4InterceptsCtrd_CorrSummary.csv")

png("Heatmap_LM5_Intercepts_OverlappingGenesCntrd.png")
heatmap(cor(LM5_Intercepts_OverlappingGenes_MedianCentered))
dev.off()

write.csv(cor(LM5_Intercepts_OverlappingGenes_MedianCentered), "LM5InterceptsCtrd_CorrSummary.csv")


png("Heatmap_LM4_InterceptRanks_OverlappingGenesCntrd.png")
heatmap(cor(LM4_InterceptRanks_OverlappingGenes_MedianCentered))
dev.off()

write.csv(cor(LM4_InterceptRanks_OverlappingGenes_MedianCentered), "LM4InterceptRanksCtrd_CorrSummary.csv")



png("Heatmap_LM5_InterceptRanks_OverlappingGenesCntrd.png")
heatmap(cor(LM5_InterceptRanks_OverlappingGenes_MedianCentered))
dev.off()

write.csv(cor(LM5_InterceptRanks_OverlappingGenes_MedianCentered), "LM5InterceptRanksCtrd_CorrSummary.csv")



png("LM4_Intercepts_OverlappingGenesCntrd_CorrSummary.png", width = 1200, height = 1200)
pairs(LM4_Intercepts_OverlappingGenes_MedianCentered,
	gap=0,
	diag.panel=function(x,...){
		par(new=TRUE)
		hist(x,
			col="light blue",
			probability = TRUE,
			axes=FALSE,
			main="")
		lines(density(x), 
			col="red",
			lwd=3)
		rug(x)
		})
dev.off()



png("LM5_Intercepts_OverlappingGenesCntrd_CorrSummary.png", width = 1200, height = 1200)
pairs(LM5_Intercepts_OverlappingGenes_MedianCentered,
	gap=0,
	diag.panel=function(x,...){
		par(new=TRUE)
		hist(x,
			col="light blue",
			probability = TRUE,
			axes=FALSE,
			main="")
		lines(density(x), 
			col="red",
			lwd=3)
		rug(x)
		})
dev.off()



png("LM4_InterceptRanks_OverlappingGenesCntrd_CorrSummary.png", width = 1200, height = 1200)
pairs(LM4_InterceptRanks_OverlappingGenes_MedianCentered,
	gap=0,
	diag.panel=function(x,...){
		par(new=TRUE)
		hist(x,
			col="light blue",
			probability = TRUE,
			axes=FALSE,
			main="")
		lines(density(x), 
			col="red",
			lwd=3)
		rug(x)
		})
dev.off()



png("LM5_InterceptRanks_OverlappingGenesCntrd_CorrSummary.png", width = 1200, height = 1200)
pairs(LM5_InterceptRanks_OverlappingGenes_MedianCentered,
	gap=0,
	diag.panel=function(x,...){
		par(new=TRUE)
		hist(x,
			col="light blue",
			probability = TRUE,
			axes=FALSE,
			main="")
		lines(density(x), 
			col="red",
			lwd=3)
		rug(x)
		})
dev.off()


pcaLM5_InterceptRanks<-prcomp(t(LM5_InterceptRanks_OverlappingGenes_MedianCentered))
tmp<-pcaLM5_InterceptRanks$x[,1:10]
rownames(tmp)<-colnames(LM5_InterceptRanks_OverlappingGenes_MedianCentered)
write.csv(tmp, "LM5_InterceptRanksCntrd_PCA.csv")

tmp<-pcaLM5_InterceptRanks$rotation[,1:10]
write.csv(tmp, "pcaLM5_InterceptRanksCntrd_Eigenvectors.csv")

png("LM5_InterceptRanksCntrd_PC1vsPC2.png")
plot(pcaLM5_InterceptRanks$x[,1]~pcaLM5_InterceptRanks$x[,2], xlab="PC2", ylab="PC1", main="PC1 vs. PC2")
dev.off()

png("LM5_InterceptRanksCntrd_PCA_ScreePlot.png")
plot(summary(pcaLM5_InterceptRanks)$importance[2,]~(c(1:length(summary(pcaLM5_InterceptRanks)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC#", ylab="Proportion of Variance Explained", col=2)
dev.off()

png("LM5_InterceptRanksCntrd_PCA_ScreePlot2.png")
plot(summary(pcaLM5_InterceptRanks)$importance[3,]~(c(1:length(summary(pcaLM5_InterceptRanks)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC#", ylab="Cumulative Proportion of Variance Explained", col=2)
dev.off()

png("PCALM5_InterceptRanksCntrd_vs_PCAMDD.png")
plot(pcaLM5_InterceptRanks$x[,1]~pcaMDD$x[,1])
dev.off()

png("EigenVector1_LM5_InterceptRanksCntrd_vs_MDD.png")
plot(pcaLM5_InterceptRanks$rotation[,1]~pcaMDD$rotation[,1])
dev.off()


RNAConcbyRegion<-as.matrix(read.csv("RNAConcbyRegion.csv", header=T, row.names=1))
RNAIntegbyRegion<-as.matrix(read.csv("RNAIntegbyRegion.csv", header=T, row.names=1))

png("Heatmap of RNAConcByRegion.png")
heatmap(cor(RNAConcbyRegion))
dev.off()

png("Heatmap of RNAConcBySubject.png")
heatmap(cor(t(RNAConcbyRegion)))
dev.off()

temp<-cor(RNAConcbyRegion)
write.csv(temp, "RNAConcbyRegion_CorrSummary.csv")

temp<-cor(t(RNAConcbyRegion))
write.csv(temp, "RNAConcbySubject_CorrSummary.csv")


temp<-t(RNAConcbyRegion)
png("RNAConcbySubject_CorrSummary.png", width = 2400, height = 2400)
pairs(temp,
	gap=0,
	diag.panel=function(x,...){
		par(new=TRUE)
		hist(x,
			col="light blue",
			probability = TRUE,
			axes=FALSE,
			main="")
		lines(density(x), 
			col="red",
			lwd=3)
		rug(x)
		})
dev.off()


colnames(RNAConcbyRegion)
#[1] "Conc.AAA" "Conc.AB"  "Conc.AHA" "Conc.B"   "Conc.CE"  "Conc.CO"  "Conc.L"   "Conc.M"   "Conc.PAC" "Conc.PL" 

RNAConcbyRegionDissect1<-RNAConcbyRegion[, c(1, 2, 4, 6, 7, 9, 10)]

RNAConcbyRegionDissect2<-RNAConcbyRegion[, c(3, 5, 8)]


temp<-t(RNAConcbyRegionDissect1)


png("RNAConcbySubjectDissect1_CorrSummary.png", width = 2400, height = 2400)
pairs(temp,
	gap=0,
	diag.panel=function(x,...){
		par(new=TRUE)
		hist(x,
			col="light blue",
			probability = TRUE,
			axes=FALSE,
			main="")
		lines(density(x), 
			col="red",
			lwd=3)
		rug(x)
		})
dev.off()

temp<-cor(t(RNAConcbyRegionDissect1))
write.csv(temp, "RNAConcbySubjectDissect1_CorrSummary.csv")


png("Heatmap of RNAIntegbyRegion.png")
heatmap(cor(RNAIntegbyRegion))
dev.off()

temp<-cor(RNAIntegbyRegion)
write.csv(temp, "RNAIntegbyRegion_CorrSummary.csv")


DissectionOrder<-as.matrix(read.csv("Subject Order for LCM for R.csv", header=T, row.names=1))


for(i in 1:10){
png(paste(colnames(RNAConcbyRegion)[i], "vs", "DissectionOrder1", ".png", sep=""))
plot(RNAConcbyRegion[,i]~DissectionOrder[,1])
temp<-lm(RNAConcbyRegion[,i]~DissectionOrder[,1])
abline(temp, col=2)
mtext(paste("R-Squared=", round(summary.lm(temp)$r.squared, digits=3), sep=" "))
dev.off()
}

for(i in 1:10){
png(paste(colnames(RNAConcbyRegion)[i], "vs", "DissectionOrder2", ".png", sep=""))
plot(RNAConcbyRegion[,i]~DissectionOrder[,2])
temp<-lm(RNAConcbyRegion[,i]~DissectionOrder[,2])
abline(temp, col=2)
mtext(paste("R-Squared=", round(summary.lm(temp)$r.squared, digits=3), sep=" "))
dev.off()
}


#************Re-quantile normalizing data across brain regions*********************

#Alright, since it looks like only the rank data is standing up to the litmus test of correctly identifying appropriate regional variations in our litmus test genes, I'm going to try going back and quantile normalizing across all brain regions before running linear models. 

dim(AllRegions_NormNoOutliers_OverlappingGenes)
[1] 15788   404

AllRegions_NormNoOutliers_OverlappingGenes[c(1:3), c(1:3)]
#The data structure:
             AB_CTRL_4235 AB_CTRL_4623 AB_CTRL_4638
ILMN_2055271     7.721411     7.187390     7.854929
ILMN_1814316     8.046176     8.848606     8.849257
ILMN_2359168     8.068351     8.867992     8.724161

#"If you make use of quantile normalization please cite Bolstad et al, Bioinformatics (2003)."

AllRegions_NormNoOutliers_OverlappingGenes_ReNorm<-matrix(0, nrow=15788, ncol=404)

AllRegions_NormNoOutliers_OverlappingGenes_ReNorm<-normalize.quantiles(AllRegions_NormNoOutliers_OverlappingGenes)
row.names(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm)<-row.names(AllRegions_NormNoOutliers_OverlappingGenes)
colnames(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm)<-colnames(AllRegions_NormNoOutliers_OverlappingGenes)


dim(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm)
[1] 15788   404

AllRegions_NormNoOutliers_OverlappingGenes_ReNorm[c(1:3), c(1:3)]

#What the data looks like now:
             AB_CTRL_4235 AB_CTRL_4623 AB_CTRL_4638
ILMN_2055271     7.614650     7.127730     7.726372
ILMN_1814316     7.905285     8.616169     8.616726
ILMN_2359168     7.925292     8.633456     8.503379


#Output the quantile normalized data:

AllRegions_NormNoOutliers_OverlappingGenes_ReNorm2<-cbind(row.names(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm), AllRegions_NormNoOutliers_OverlappingGenes_ReNorm)

colnames(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm2)[1]<-"X"

#What the data looks like now:
AllRegions_NormNoOutliers_OverlappingGenes_ReNorm2[c(1:3), c(1:3)]
             X              AB_CTRL_4235       AB_CTRL_4623      
ILMN_2055271 "ILMN_2055271" "7.61465003703849" "7.12772956433538"
ILMN_1814316 "ILMN_1814316" "7.90528511396674" "8.61616931892008"
ILMN_2359168 "ILMN_2359168" "7.92529218233563" "8.63345634763011"

AllRegions_NormNoOutliers_OverlappingGenes_ReNorm2<- join(as.data.frame(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm2), ProbeInfoOutput, by = "X")

write.table(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm2, "AllRegions_NormNoOutliers_OverlappingGenes_ReNorm2.txt", sep="\t")

rm( AllRegions_NormNoOutliers_OverlappingGenes_ReNorm2)

#Double-check: a boxplot illustrating the (now identical) signal distributions for each sample:

boxplot(data.frame(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm), main="Boxplot of normalized signal values per sample (1 box= all filtered probes)", xlab="Sample ID", ylab="Quantile Normalized Signal")
#Wow -that is ugly as heck, but it gets the point across.

#Visualize the sample-sample correlations using a heatmap:
png("Sample Sample Correlations Heatmap_All10RegionsReNorm.png")
image(cor(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm), main="Visualizing the correlations between entire samples (by index#)", xlab="Red=Less correlated, Light yellow=Highly correlated")
dev.off()
#That's interesting - samples from particular brain regions don't necessarily really seem to be more correlated than samples between brain regions
#There also doesn't necessarily appear to be a correlation between samples across brain regions.

#I should probably test that. 

colnames(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm)[]

AllRegions_SampleInfo<-rbind(
	cbind(matrix("AB", nrow=length(AB_SampleInfo[,1]), ncol=1), AB_SampleInfo),
	cbind(matrix("Central", nrow=length(Central_SampleInfo[,1]), ncol=1), Central_SampleInfo),	
	cbind(matrix("CO", nrow=length(CO_SampleInfo[,1]), ncol=1), CO_SampleInfo),	
	cbind(matrix("Lateral", nrow=length(Lateral_SampleInfo[,1]), ncol=1), Lateral_SampleInfo),	
	cbind(matrix("PAC", nrow=length(PAC_SampleInfo[,1]), ncol=1), PAC_SampleInfo),			
	cbind(matrix("AAA", nrow=length(AAA_SampleInfo[,1]), ncol=1), AAA_SampleInfo),
	cbind(matrix("AHA", nrow=length(AHA_SampleInfo[,1]), ncol=1), AHA_SampleInfo),	
	cbind(matrix("Basal", nrow=length(Basal_SampleInfo[,1]), ncol=1), Basal_SampleInfo),	
	cbind(matrix("Medial", nrow=length(Medial_SampleInfo[,1]), ncol=1), Medial_SampleInfo),	
	cbind(matrix("PL", nrow=length(PL_SampleInfo[,1]), ncol=1), PL_SampleInfo)		
)

colnames(AllRegions_SampleInfo)[1]<-"Region"

dim(AllRegions_SampleInfo)
[1] 404  22

AllRegions_SampleInfo[c(1:3), c(1:3)]
     Region X      Sample.Processing.Order
[1,] "AB"   "4235" " 1"                   
[2,] "AB"   "4623" " 2"                   
[3,] "AB"   "4638" " 3" 

#Double Checking that the sample order matches that of AllRegions_NormNoOutliers_OverlappingGenes_ReNorm
cbind(colnames(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm), AllRegions_SampleInfo[, c(1:2)] )
#First attempt didn't work - fixed code and came back to this and it looked fine.

#Moving some of the sample info into named variables:

colnames(AllRegions_SampleInfo)
 [1] "Region"                  "X"                       "Sample.Processing.Order"
 [4] "GeneralChip"             "LocationOnChip"          "Sample.Group"           
 [7] "Subject.Number"          "Cohort"                  "Gender"                 
[10] "Age"                     "Age.of.RNA..years."      "TOD..hrs."              
[13] "Suicide"                 "Overdose"                "pH"                     
[16] "RNAConc"                 "RNAIntegrity"            "Hours.Cold"             
[19] "Hours.Ice"               "PMI"                     "HoursColdPlusIce"       
[22] "PMICorrected"  

AllRegions_Region<-as.factor(AllRegions_SampleInfo[,1])
AllRegions_SubjectID<-as.factor(AllRegions_SampleInfo[,2])
AllRegions_Diagnosis<-as.factor(AllRegions_SampleInfo[,6])
AllRegions_Gender<-as.factor(AllRegions_SampleInfo[,9])

AllRegions_Age<-as.numeric(AllRegions_SampleInfo[,10])
AllRegions_TOD<-as.numeric(AllRegions_SampleInfo[,12])

AllRegions_Suicide<-as.factor(AllRegions_SampleInfo[,13])
AllRegions_BrainPH<-as.numeric(AllRegions_SampleInfo[,15])
AllRegions_PMICorrected<-as.numeric(AllRegions_SampleInfo[,22])



#A quick double-check that everything was coded properly:

#The previously identified top pH-related gene: HSD17B14
plot(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm[row.names(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm)=="ILMN_1809483", ]~AllRegions_BrainPH)
#Wow - without taking into account region, that data is a mess. Subject is pretty obvious in the data too.

#Another attempt to double check things: gender gene XIST
plot(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm[row.names(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm)=="ILMN_1764573", ]~AllRegions_Gender)
#Alright, that worked, so I think the data is functional, it just is messy due to regional variation.


#RNA Integrity is going to be weird because it differs by a little bit for each brain region. It might make sense to replace it with average RNA Integrity:

AverageRNAIntegrity_AllSubjects<-read.csv("AverageRNAIntegrity_AllSubjects.csv", header=TRUE, sep=",")

colnames(AverageRNAIntegrity_AllSubjects)[1]<-"Subject.Number"

AllRegions_SampleInfo2<-join(as.data.frame(AllRegions_SampleInfo), as.data.frame(AverageRNAIntegrity_AllSubjects), by="Subject.Number")

AllRegions_SampleInfo2[c(1:3), ]
colnames(AllRegions_SampleInfo2)

AllRegions_RNAIntegrity<-as.numeric(AllRegions_SampleInfo2[,23])

#I should really median-center the variables before running models with them too:

AllRegions_AgeCentered<-AllRegions_Age-50
hist(AllRegions_AgeCentered)
AllRegions_RNAIntegrityCentered<-AllRegions_RNAIntegrity-5
hist(AllRegions_RNAIntegrityCentered)
AllRegions_BrainPHCentered<-AllRegions_BrainPH-6.8
hist(AllRegions_BrainPHCentered)
AllRegions_PMICorrectedCentered<-AllRegions_PMICorrected-22.5
hist(AllRegions_PMICorrectedCentered)

#If I want to do this analysis properly, I will probably need to run a multilevel model that takes into account the fact that the same subjects are included in each brain region or the effect of various subject variables will be inflated. I will also need to include interaction terms with region for all variables (ouch!). However, in the meantime I just want to test whether things are generally working properly. 

LinearModel3<-function(i){lm(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm[i,]~AllRegions_Diagnosis+AllRegions_AgeCentered+AllRegions_RNAIntegrityCentered+AllRegions_BrainPHCentered+AllRegions_Gender+AllRegions_PMICorrectedCentered+AllRegions_Region)}

##3. Then run this code (Ctrl+R) to get the p-values:


#***CODE TO RUN FOR STEP 17****

#Determining how many variables are in the linear model (since LM breaks factors down into dummy variables)
NumberofXsLM3<-length(summary.lm(LinearModel3(1))$coefficients[,1])
NameofXsLM3<-dimnames(summary.lm(LinearModel3(1))$coefficients)[1][[1]]

#Running the linear model:
LM3Betas<-matrix(0, length(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm[,1]), NumberofXsLM3)
colnames(LM3Betas)<-NameofXsLM3
row.names(LM3Betas)<-row.names(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm)

LM3pvalues<-matrix(0, length(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm[,1]), NumberofXsLM3)
colnames(LM3pvalues)<-NameofXsLM3
row.names(LM3pvalues)<-row.names(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm)

LM3SE<-matrix(0, length(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm[,1]), NumberofXsLM3)
colnames(LM3SE)<-NameofXsLM3
row.names(LM3SE)<-row.names(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm)

LM3Tstat<-matrix(0, length(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm[,1]), NumberofXsLM3)
colnames(LM3Tstat)<-NameofXsLM3
row.names(LM3Tstat)<-row.names(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm)



for(i in 1:length(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm[,1])){
	TempLM3<-LinearModel3(i)
	LM3Betas[i,]<-summary.lm(TempLM3)$coefficients[,1]
	LM3SE[i,]<-summary.lm(TempLM3)$coefficients[,2]
	LM3Tstat[i,]<-summary.lm(TempLM3)$coefficients[,3]
	LM3pvalues[i,]<-summary.lm(TempLM3)$coefficients[,4]
}


#Outputting a histogram illustrating the raw p-values in comparison to what would be expected by chance:
for (i in 1:NumberofXsLM3){
png(paste(paste("17 Histogram of Raw Pvalues Using LM3 for", NameofXsLM3[i], sep="  "), "png", sep="."))	
hist(LM3pvalues[,i], breaks=100, col=i, main=paste("Raw P-values using LM3 for", NameofXsLM3[i], sep="  "), xlab="Unadjusted p-value", ylab="Count")
abline(a=(length(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm[,1])/100), b=0)
dev.off()		
}		



#Outputting the raw pvalues and LM related statistics:

LM3BetasOutput<-cbind(row.names(LM3Betas), LM3Betas) 
colnames(LM3BetasOutput)[1]<-"X"

LM3BetasOutput2<-join(as.data.frame(LM3BetasOutput), ProbeInfoOutput, by="X")

LM3SEOutput<-cbind(LM3Betas, LM3BetasOutput2[ , c(18:25)])
LM3TstatOutput<-cbind(LM3Tstat, LM3BetasOutput2[ , c(18:25)])
LM3pvaluesOutput<-cbind(LM3pvalues, LM3BetasOutput2[ , c(18:25)])

write.csv(LM3BetasOutput2, "LM3Betas.csv")
write.csv(LM3SEOutput, "LM3SE.csv")
write.csv(LM3TstatOutput, "LM3Tstat.csv")
write.csv(LM3pvaluesOutput, "LM3pvaluesRAW.csv")


#It looks like the results are not quite the same as what I got by comparing effects (pH, gender, diagnosis, etc) across different brain regions using Fisher's p-value. Let's see:

colnames(LM3Betas)
[1] "(Intercept)"                     "AllRegions_DiagnosisMDD"        
 [3] "AllRegions_AgeCentered"          "AllRegions_RNAIntegrityCentered"
 [5] "AllRegions_BrainPHCentered"      "AllRegions_GenderM"             
 [7] "AllRegions_PMICorrectedCentered" "AllRegions_RegionAB"            
 [9] "AllRegions_RegionAHA"            "AllRegions_RegionBasal"         
[11] "AllRegions_RegionCentral"        "AllRegions_RegionCO"            
[13] "AllRegions_RegionLateral"        "AllRegions_RegionMedial"        
[15] "AllRegions_RegionPAC"            "AllRegions_RegionPL" 

colnames(PH_Summary_CrossRegionalResults)
 [1] "X"               "AAA"             "AB"              "AHA"             "Basal"          
 [6] "Central"         "CO"              "Lateral"         "Medial"          "PAC"            
[11] "PL"              "Median Pvalue"   "Mean Log Pvalue" "Chi-Square Stat" "Fishers Pvalue" 
[16] "AAA"             "AB"              "AHA"             "Basal"           "Central"        
[21] "CO"              "Lateral"         "Medial"          "PAC"             "PL"             
[26] "Median Beta"     "Mean Beta"       "SE Beta"   


png("PH LM3 Beta vs MedianBeta.png")
plot(LM3Betas[,5]~PH_Summary_CrossRegionalResults[,26])
dev.off()

png("MDD LM3 Beta vs MedianBeta.png")
plot(LM3Betas[,2]~MDD_Summary_CrossRegionalResults[,26])
dev.off()

#There seems to be a reasonably good correlation between the two sets of results
#When examining the litmus test genes, it seems like the betas associated with the central lateral and basal nuclei have similar relative relationships to those found in my intercept rank calculations and to previously published results. Next I need to see what happens when I add regional interaction terms. 


LinearModel6<-function(i){lm(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm[i,]~AllRegions_Diagnosis+AllRegions_AgeCentered+AllRegions_RNAIntegrityCentered+AllRegions_BrainPHCentered+AllRegions_Gender+AllRegions_PMICorrectedCentered+AllRegions_Region+AllRegions_Region*AllRegions_Diagnosis+AllRegions_Region*AllRegions_AgeCentered+AllRegions_Region*AllRegions_RNAIntegrityCentered+AllRegions_Region*AllRegions_BrainPHCentered+AllRegions_Region*AllRegions_Gender+AllRegions_Region*AllRegions_PMICorrectedCentered
	)}

##3. Then run this code (Ctrl+R) to get the p-values:


#***CODE TO RUN FOR STEP 17****

#Determining how many variables are in the linear model (since LM breaks factors down into dummy variables)
NumberofXsLM6<-length(summary.lm(LinearModel6(1))$coefficients[,1])
NameofXsLM6<-dimnames(summary.lm(LinearModel6(1))$coefficients)[1][[1]]

#Running the linear model:
LM6Betas<-matrix(0, length(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm[,1]), NumberofXsLM6)
colnames(LM6Betas)<-NameofXsLM6
row.names(LM6Betas)<-row.names(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm)

LM6pvalues<-matrix(0, length(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm[,1]), NumberofXsLM6)
colnames(LM6pvalues)<-NameofXsLM6
row.names(LM6pvalues)<-row.names(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm)

LM6SE<-matrix(0, length(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm[,1]), NumberofXsLM6)
colnames(LM6SE)<-NameofXsLM6
row.names(LM6SE)<-row.names(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm)

LM6Tstat<-matrix(0, length(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm[,1]), NumberofXsLM6)
colnames(LM6Tstat)<-NameofXsLM6
row.names(LM6Tstat)<-row.names(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm)



for(i in 1:length(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm[,1])){
	TempLM6<-LinearModel6(i)
	LM6Betas[i,]<-summary.lm(TempLM6)$coefficients[,1]
	LM6SE[i,]<-summary.lm(TempLM6)$coefficients[,2]
	LM6Tstat[i,]<-summary.lm(TempLM6)$coefficients[,3]
	LM6pvalues[i,]<-summary.lm(TempLM6)$coefficients[,4]
}


#Outputting a histogram illustrating the raw p-values in comparison to what would be expected by chance:
for (i in 1:NumberofXsLM6){
png(paste(paste("17 Histogram of Raw Pvalues Using LM6 for", NameofXsLM6[i], sep="  "), "png", sep="."))	
hist(LM6pvalues[,i], breaks=100, col=i, main=paste("Raw P-values using LM6 for", NameofXsLM6[i], sep="  "), xlab="Unadjusted p-value", ylab="Count")
abline(a=(length(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm[,1])/100), b=0)
dev.off()		
}		



#Outputting the raw pvalues and LM related statistics:

LM6BetasOutput<-cbind(row.names(LM6Betas), LM6Betas) 
colnames(LM6BetasOutput)[1]<-"X"

LM6BetasOutput2<-join(as.data.frame(LM6BetasOutput), ProbeInfoOutput, by="X")

LM6SEOutput<-cbind(LM6Betas, LM6BetasOutput2[ , c(18:25)])
LM6TstatOutput<-cbind(LM6Tstat, LM6BetasOutput2[ , c(18:25)])
LM6pvaluesOutput<-cbind(LM6pvalues, LM6BetasOutput2[ , c(18:25)])

write.csv(LM6BetasOutput2, "LM6Betas.csv")
write.csv(LM6SEOutput, "LM6SE.csv")
write.csv(LM6TstatOutput, "LM6Tstat.csv")
write.csv(LM6pvaluesOutput, "LM6pvaluesRAW.csv")

#Note: When running this in the future, we need to decide which brain region will be the intercept because the interaction terms mean that the overall main effects will be determined *at the interecept brain region* - i.e. currently "Diagnosis" is actually "Diagnosis effects in the AAA for males with average pre- and post-mortem parameters.  Similarly, p-values for effect of region actually indicate how much a region differs from the AAA.

#Hmm... Interestingly enough, the model with all of the interaction terms does a *worse* job of predicting the litmus test genes than the model without the intercept terms... except for in the Central, for which it does better. Um. Hmm.

#I'm going to try removing all of the interaction terms for RNAIntegrity because I think they are just adding noise. I'm also going to relevel region so that the Lateral nucleus is the intercept, because it will make it easier to compare the p-values for the litmus test genes.
is.factor(AllRegions_Region)
#eek - I just got false for that, so I went back and made it a factor and then re-ran LM3 and LM6 (although it seemed like my original results reflected this varible being treated as a factor with AAA as the intercept...I'm a little confused as to why)
AllRegions_Region<-relevel(AllRegions_Region, ref="Lateral")






#testing just basic ANOVA and repeated measures ANOVA first:
summary.aov(aov(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm[1,]~AllRegions_Region))

#vs.

summary.aov(aov(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm[1,]~AllRegions_Region+ Error(AllRegions_SubjectID/AllRegions_Region)))
#I got an error message for this, probably because it isn't balanced. I could either plug in median values for the missing subjects that were removed as outliers or just move on to mixed models. 

#instead:

lmer(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm[1,]~AllRegions_Region +  (1 | AllRegions_SubjectID)+(AllRegions_Region | AllRegions_SubjectID)))
#I got an error message for that too related to the fact that some of the subjects are missing data for particular regions. O.k. back to troubleshooting. Also, I'm not sure that the interaction term is really appropriate for this dataset! 
lmer(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm[1,]~AllRegions_Region +  (1 | AllRegions_SubjectID))
#that actually worked. Interesting.  But it takes a really long time to  process. And there doesn't seem to be a p-value available...?
summary(lmer(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm[1,]~AllRegions_Region +  (1 | AllRegions_SubjectID)))
#O.k. I found a long write-up online regarding how to use this package lme4... and it looks like the question of how to calculate a p-value is actually fairly complex because this is not just a least squares fit but a RLME which involves fitting multiple iterations of the model.  Since the data is not overly complex, it seems like this might be a good reason to just fill in the NA subject values with median or mean data and then run a repeated measures ANOVA.

#/Another version of mixed effects modeling. This version didn't work. It says that the the variable lengths differ, but that doesn't seem to be true. Confused.
summary(lme(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm[1,]~AllRegions_Region, random = ~ 1 | AllRegions_SubjectID))
length(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm[1,])
length(AllRegions_Region)
length(AllRegions_SubjectID)

#I think for the current regional analysis what might make sense instead is to first:
## 1) Just use repeated measures anova to remove subject-variable related effects because for the current analysis we don't really care about effects of pH, gender, etc, we just want them out of the analysis.  That said, if any of these variables are having different sorts of effects on each of the 10 nuclei they may artificially create "regional differences" that are actually regional differences in response to subject variables.  That would be an argument for using a different strategy with co-variates.  This would be a good reason (besides IP) for only running the analysis on control data. Perhaps we can see how the litmus test results turn out? 
## 2) To do that, I need to create a new dataset that has some non-data included for the outlier subjects so that the design is balanced. (either average or median data???)

#Figuring out who is missing:
table(AllRegions_SubjectID, AllRegions_Region)

#AAA is missing 4939
#AHA is missing 4905 & 4968
#Basal is missing 3204, 4590
#CO is missing 3204, 4235, 4839, 4968
#Medial is missing 3204, 4232, 4424, 
#PAC is missing 4424, 4638, 4652
#PL is missing 4310

#hmm... I should probably doublecheck what this looks like just for controls since they are such a small sample size:

AllRegions_SubjectID_Controls<-AllRegions_SubjectID[AllRegions_Diagnosis=="CTRL"]
AllRegions_Region_Controls<-AllRegions_Region[AllRegions_Diagnosis=="CTRL"]
AllRegions_NormNoOutliers_OverlappingGenes_ReNorm_Controls<-AllRegions_NormNoOutliers_OverlappingGenes_ReNorm[,AllRegions_Diagnosis=="CTRL"]

length(AllRegions_SubjectID_Controls)
AllRegions_SubjectID_Controls[c(1:3)]

length(AllRegions_Region_Controls)
AllRegions_Region_Controls[c(1:3)]

dim(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm_Controls)


table(AllRegions_SubjectID_Controls, AllRegions_Region_Controls)
#I don't understand why this code isn't working correctly but basically it looks like the missing control values are:

#AHA=4905, 4968
#CO=4235, 4839, 4968
#PAC=4638, 4652

#... which means that CO is missing ~1/5 of its control sample size (because the original n was 16)

#Note to self: R crashed at this point. **The workspace that was previously saved is from *before* I re-ran LM3 and LM6 with region defined as a factor releveled with lateral as the intercept. I have redefined region again, but I *have not* re-run LM3 and LM6 because I don't know if I need the results in the workspace at this point so it seems like a waste of time. 


AllRegions_SubjectID_Controls_wMissingNA<-as.factor(c(as.character(AllRegions_SubjectID_Controls), "4905", "4968", "4235", "4839", "4968", "4638", "4652"))
AllRegions_Region_Controls_wMissingNA<-as.factor(c(as.character(AllRegions_Region_Controls), "AHA", "AHA", "CO", "CO", "CO", "PAC", "PAC"))
AllRegions_Region_Controls_wMissingNA<-relevel(AllRegions_Region_Controls_wMissingNA, ref="Lateral")

AllRegions_NormNoOutliers_OverlappingGenes_ReNorm_Controls_wMissingNA<-cbind(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm_Controls, matrix(NA, length(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm_Controls[,1]),7))

length(AllRegions_SubjectID_Controls_wMissingNA)
length(AllRegions_Region_Controls_wMissingNA)

is.matrix(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm_Controls_wMissingNA)
dim(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm_Controls_wMissingNA)

#test case again:

summary.aov(aov(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm_Controls_wMissingNA[1,]~AllRegions_Region_Controls_wMissingNA))

#vs.

summary.aov(aov(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm_Controls_wMissingNA[1,]~AllRegions_Region_Controls_wMissingNA+ Error(AllRegions_SubjectID_Controls_wMissingNA/AllRegions_Region_Controls_wMissingNA)))

#got an error - I haven't added median values in here

temp<-matrix(NA, length(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm_Controls[,1]),7)

for(i in 1:length(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm_Controls[,1])){
temp[i, ]<-rep(mean(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm_Controls[i, ]), 7)
}

AllRegions_NormNoOutliers_OverlappingGenes_ReNorm_Controls_wMissingAve<-cbind(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm_Controls, temp)

#test case again:

summary.aov(aov(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm_Controls_wMissingAve[1,]~AllRegions_Region_Controls_wMissingNA))

#vs.

summary.aov(aov(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm_Controls_wMissingAve[1,]~AllRegions_Region_Controls_wMissingNA+ Error(AllRegions_Region_Controls_wMissingNA/AllRegions_SubjectID_Controls_wMissingNA)))
#Still an error term. Something is wrong with how I am running repeated measures ANOVA. Sigh.

table(AllRegions_Region_Controls_wMissingNA, AllRegions_SubjectID_Controls_wMissingNA)
#Ah - well that explains it.  Apparently when I coded in the new data I didn't realize that Region was being coded as numbers (factor). The same seems to be true for subject. 
                                     AllRegions_SubjectID_Controls_wMissingNA
AllRegions_Region_Controls_wMissingNA 5 7 14 15 16 17 21 24 25 27 28 31 34 36 41 42 4235 4638
                                  1   1 1  1  1  1  1  1  1  1  1  1  1  1  1  1  1    0    0
                                  10  1 1  1  1  1  1  1  1  1  1  1  1  1  1  1  1    0    0
                                  2   1 1  1  1  1  1  1  1  1  1  1  1  1  1  1  1    0    0
                                  3   1 1  1  1  1  1  1  1  1  1  1  1  1  1  1  1    0    0
                                  4   1 1  1  1  1  1  1  1  1  1  1  0  1  0  1  1    0    0
                                  5   1 1  1  1  1  1  1  1  1  1  1  1  1  1  1  1    0    0
                                  6   1 1  1  1  1  1  1  1  1  1  1  1  1  1  1  1    0    0
                                  7   0 1  1  1  1  1  1  1  1  1  0  1  1  0  1  1    0    0
                                  8   1 1  1  1  1  1  1  1  1  1  1  1  1  1  1  1    0    0
                                  9   1 1  1  0  0  1  1  1  1  1  1  1  1  1  1  1    0    0
                                  AHA 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0    0    0
                                  CO  0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0    1    0
                                  PAC 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0    0    1
                                     AllRegions_SubjectID_Controls_wMissingNA
AllRegions_Region_Controls_wMissingNA 4652 4839 4905 4968
                                  1      0    0    0    0
                                  10     0    0    0    0
                                  2      0    0    0    0
                                  3      0    0    0    0
                                  4      0    0    0    0
                                  5      0    0    0    0
                                  6      0    0    0    0
                                  7      0    0    0    0
                                  8      0    0    0    0
                                  9      0    0    0    0
                                  AHA    0    0    1    1
                                  CO     0    1    0    1
                                  PAC    1    0    0    0

#Alright, I fixed the problem by going back and forcing the factors ID and Region to be characters when combining them with the new datat, then forcing the data back into factor format. *phew*

levels(AllRegions_Region_Controls_wMissingNA)

 [1] "Lateral" "AAA"     "AB"      "AHA"     "Basal"   "Central" "CO"      "Medial"  "PAC"    
[10] "PL"  


table(AllRegions_Region_Controls_wMissingNA, AllRegions_SubjectID_Controls_wMissingNA)

                                     AllRegions_SubjectID_Controls_wMissingNA
AllRegions_Region_Controls_wMissingNA 4235 4375 4623 4638 4652 4664 4706 4744 4754 4830 4839
                              Lateral    1    1    1    1    1    1    1    1    1    1    1
                              AAA        1    1    1    1    1    1    1    1    1    1    1
                              AB         1    1    1    1    1    1    1    1    1    1    1
                              AHA        1    1    1    1    1    1    1    1    1    1    1
                              Basal      1    1    1    1    1    1    1    1    1    1    1
                              Central    1    1    1    1    1    1    1    1    1    1    1
                              CO         1    1    1    1    1    1    1    1    1    1    1
                              Medial     1    1    1    1    1    1    1    1    1    1    1
                              PAC        1    1    1    1    1    1    1    1    1    1    1
                              PL         1    1    1    1    1    1    1    1    1    1    1
                                     AllRegions_SubjectID_Controls_wMissingNA
AllRegions_Region_Controls_wMissingNA 4905 4945 4968 5000 5008
                              Lateral    1    1    1    1    1
                              AAA        1    1    1    1    1
                              AB         1    1    1    1    1
                              AHA        1    1    1    1    1
                              Basal      1    1    1    1    1
                              Central    1    1    1    1    1
                              CO         1    1    1    1    1
                              Medial     1    1    1    1    1
                              PAC        1    1    1    1    1
                              PL         1    1    1    1    1

#O.k. back to attempting to run repeated measures ANOVA:
summary.aov(aov(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm_Controls_wMissingAve[1,]~AllRegions_Region_Controls_wMissingNA+ Error(AllRegions_Region_Controls_wMissingNA/AllRegions_SubjectID_Controls_wMissingNA)))
#Error in 1L:object$rank : argument of length 0

#Oops - I think I reversed the syntax. Let's try again:
aov(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm_Controls_wMissingAve[1,]~AllRegions_Region_Controls_wMissingNA+ Error(AllRegions_SubjectID_Controls_wMissingNA/AllRegions_Region_Controls_wMissingNA))
#O.k., that seems to be working maybe.  
#From: http://ww2.coastal.edu/kingw/statistics/R-tutorials/repeated.html


temp<-summary(aov(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm_Controls_wMissingAve[1,]~AllRegions_Region_Controls_wMissingNA+ Error(AllRegions_SubjectID_Controls_wMissingNA/AllRegions_Region_Controls_wMissingNA)))
#Heck yeah, that worked. It is a little less significant than before (when I didn't include a repeated measures component). I wonder if that is because of the added fake interpolated data (Averages)


temp:

Error: AllRegions_SubjectID_Controls_wMissingNA
          			Df Sum Sq Mean Sq F value Pr(>F)
Residuals 15 0.7506 0.05004               

Error: AllRegions_SubjectID_Controls_wMissingNA:AllRegions_Region_Controls_wMissingNA
                                      								 Df Sum Sq Mean Sq F value   Pr(>F)    
AllRegions_Region_Controls_wMissingNA   9  4.044  0.4494   6.361 1.72e-07 ***
						Residuals                             135  9.536  0.0706                     
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#These are the df for the Region and Residuals:
unlist(temp)[6][[1]]
Error: AllRegions_SubjectID_Controls_wMissingNA:AllRegions_Region_Controls_wMissingNA.Df1 
                                                                                        9
unlist(temp)[7][[1]]
Error: AllRegions_SubjectID_Controls_wMissingNA:AllRegions_Region_Controls_wMissingNA.Df2 
                                                                                      135 
 #Extracting the F-value:                                                              
unlist(temp)[12][[1]]
Error: AllRegions_SubjectID_Controls_wMissingNA:AllRegions_Region_Controls_wMissingNA.F value1 
                                                                                      6.361142 
#Extracting the p-value:                                                                                   
unlist(temp)[14][[1]]                                                                              
Error: AllRegions_SubjectID_Controls_wMissingNA:AllRegions_Region_Controls_wMissingNA.Pr(>F)1 
                                                                                 1.719329e-07 
                                                                                 

#Running repeated measures ANOVA:

Control_RepeatedMeasuresANOVAresults<-matrix(NA, nrow=length(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm_Controls_wMissingAve[,1]), ncol=4)
row.names(Control_RepeatedMeasuresANOVAresults)<-row.names(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm_Controls_wMissingAve)
colnames(Control_RepeatedMeasuresANOVAresults)<-c("Df1", "Df2", "F", "p-value")
head(Control_RepeatedMeasuresANOVAresults)

 for(i in 1:length(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm_Controls_wMissingAve[,1])){
	temp<-summary(aov(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm_Controls_wMissingAve[i,]~AllRegions_Region_Controls_wMissingNA+ Error(AllRegions_SubjectID_Controls_wMissingNA/AllRegions_Region_Controls_wMissingNA)))
	
	Control_RepeatedMeasuresANOVAresults[i,1]<-unlist(temp)[6][[1]]
	Control_RepeatedMeasuresANOVAresults[i,2]<-unlist(temp)[7][[1]]
	Control_RepeatedMeasuresANOVAresults[i,3]<-unlist(temp)[12][[1]]
	Control_RepeatedMeasuresANOVAresults[i,4]<-unlist(temp)[14][[1]]    
}

head(Control_RepeatedMeasuresANOVAresults)



#Outputting a histogram illustrating the raw p-values in comparison to what would be expected by chance:
png("Histogram of Raw Pvalues for Control_Repeated Measures ANOVA.png")
hist(Control_RepeatedMeasuresANOVAresults[,4], breaks=100, col=2,  xlab="Unadjusted p-value", ylab="Count")
abline(a=(length(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm_Controls_wMissingAve[,1])/100), b=0)
dev.off()		
		

#Adding a calculation of the average values and Z-scores for each region to the output:

Control_AverageValuesByRegion<-matrix(NA, nrow=length(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm_Controls_wMissingAve[,1]), ncol=10)
row.names(Control_AverageValuesByRegion)<-row.names(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm_Controls_wMissingAve)
colnames(Control_AverageValuesByRegion)<-levels(AllRegions_Region_Controls_wMissingNA)
head(Control_AverageValuesByRegion)


temp<-tapply(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm_Controls_wMissingAve[1,], AllRegions_Region_Controls_wMissingNA, mean)
str(temp)

 for(i in 1:length(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm_Controls_wMissingAve[,1])){

Control_AverageValuesByRegion[i,]<-tapply(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm_Controls_wMissingAve[i,], AllRegions_Region_Controls_wMissingNA, mean)
}

head(Control_AverageValuesByRegion)

Control_AverageValuesFor10Regions<-apply(Control_AverageValuesByRegion, 1, mean)


Control_StDevByRegion<-matrix(NA, nrow=length(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm_Controls_wMissingAve[,1]), ncol=10)
row.names(Control_StDevByRegion)<-row.names(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm_Controls_wMissingAve)
colnames(Control_StDevByRegion)<-levels(AllRegions_Region_Controls_wMissingNA)
head(Control_StDevByRegion)


temp<-tapply(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm_Controls_wMissingAve[1,], AllRegions_Region_Controls_wMissingNA, sd)
str(temp)

 for(i in 1:length(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm_Controls_wMissingAve[,1])){

Control_StDevByRegion[i,]<-tapply(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm_Controls_wMissingAve[i,], AllRegions_Region_Controls_wMissingNA, sd)
}

head(Control_StDevByRegion)


Control_ZscoreByRegion<-matrix(NA, nrow=length(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm_Controls_wMissingAve[,1]), ncol=10)
row.names(Control_ZscoreByRegion)<-row.names(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm_Controls_wMissingAve)
colnames(Control_ZscoreByRegion)<-levels(AllRegions_Region_Controls_wMissingNA)
head(Control_ZscoreByRegion)




for(i in 1:length(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm_Controls_wMissingAve[,1])){

Control_ZscoreByRegion[i,]<-(Control_AverageValuesByRegion[i,]-Control_AverageValuesFor10Regions[i])/Control_StDevByRegion[i,]

}

Control_RepeatedMeasuresANOVAresults_Output<-cbind(row.names(Control_RepeatedMeasuresANOVAresults), Control_RepeatedMeasuresANOVAresults, Control_AverageValuesByRegion, Control_AverageValuesFor10Regions, Control_StDevByRegion, Control_ZscoreByRegion)
colnames(Control_RepeatedMeasuresANOVAresults_Output)[1]<-"X"

for(i in 1:10){
	j<-i+5
	colnames(Control_RepeatedMeasuresANOVAresults_Output)[j]<-paste("AVE", levels(AllRegions_Region_Controls_wMissingNA)[i], sep="_")
}

for(i in 1:10){
	j<-i+16
	colnames(Control_RepeatedMeasuresANOVAresults_Output)[j]<-paste("STDEV", levels(AllRegions_Region_Controls_wMissingNA)[i], sep="_")
}

for(i in 1:10){
	j<-i+26
	colnames(Control_RepeatedMeasuresANOVAresults_Output)[j]<-paste("ZScore", levels(AllRegions_Region_Controls_wMissingNA)[i], sep="_")
}



head(Control_RepeatedMeasuresANOVAresults_Output)


#Outputting the raw pvalues and LM related statistics:

Control_RepeatedMeasuresANOVAresults_Output2<-join(as.data.frame(Control_RepeatedMeasuresANOVAresults_Output), ProbeInfoOutput, by="X")

write.csv(Control_RepeatedMeasuresANOVAresults_Output2, "Control_RepeatedMeasuresANOVAresults_Output2.csv")
 
png("CorrelationSummary_ZscoresByRegion.png", width = 1200, height = 1200)
pairs(Control_ZscoreByRegion,
	gap=0,
	diag.panel=function(x,...){
		par(new=TRUE)
		hist(x,
			col="light blue",
			probability = TRUE,
			axes=FALSE,
			main="")
		lines(density(x), 
			col="red",
			lwd=3)
		rug(x)
		})
dev.off()
     
       
 CorrelationCoefficientZscoresbyRegion_Controls<-cor(Control_ZscoreByRegion) 
 write.csv( CorrelationCoefficientZscoresbyRegion_Controls, " CorrelationCoefficientZscoresbyRegion_Controls.csv")
 
 
# For SFN poster:
 
 colnames(Control_ZscoreByRegion)
 [1] "Lateral" "AAA"     "AB"      "AHA"     "Basal"   "Central" "CO"     
 [8] "Medial"  "PAC"     "PL"
 
Control_ZscoreByRegionReordered<-Control_ZscoreByRegion[, c(1,2,6,7, 8, 4, 3, 10,9,5)]

 colnames(Control_ZscoreByRegionReordered)  
 [1] "Lateral" "AAA"     "Central" "CO"      "Medial"  "AHA"     "AB"     
 [8] "PL"      "PAC"     "Basal" 
 
 png("CorrelationSummary_ZscoresByRegionSFN.png", width = 1200, height = 1200)
pairs(Control_ZscoreByRegionReordered, 
	gap=0, upper.panel = NULL,
	diag.panel=function(x,...){
		par(new=TRUE)
		hist(x,
			col="light blue",
			probability = TRUE,
			axes=FALSE,
			main="")
		lines(density(x), 
			col="red",
			lwd=3)
		rug(x)
		})
dev.off()
 
  CorrelationCoefficientZscoresbyRegion_ControlsSFN<-cor(Control_ZscoreByRegionReordered) 
 write.csv( CorrelationCoefficientZscoresbyRegion_ControlsSFN, " CorrelationCoefficientZscoresbyRegion_ControlsSFN.csv")
 
 
 
 
 
 
 
 png("HeatmapCorrelation_ZscoresByRegion_Controls.png")      
 heatmap(cor(Control_ZscoreByRegion)) 
 dev.off()
 
  
 dim(Control_ZscoreByRegion)  
             
 summary.lm(lm(Control_ZscoreByRegion[,1]~Control_ZscoreByRegion[,2]))  
 
 

 
pcaControl_ZscoreByRegion<-prcomp(t(Control_ZscoreByRegion))
tmp<-pcaControl_ZscoreByRegion$x[,1:10]
rownames(tmp)<-colnames(Control_ZscoreByRegion)
write.csv(tmp, "Control_ZscoreByRegion_PCA.csv")

tmp<-pcaControl_ZscoreByRegion$rotation[,1:10]
write.csv(tmp, "Control_ZscoreByRegion_PCA_Eigenvectors.csv")

png("Control_ZscoreByRegion_PC1vsPC2.png")
plot(pcaControl_ZscoreByRegion$x[,1]~pcaControl_ZscoreByRegion$x[,2], xlab="PC2", ylab="PC1", main="PC1 vs. PC2")
dev.off()

png("Control_ZscoreByRegion_PCA_ScreePlot.png")
plot(summary(pcaControl_ZscoreByRegion)$importance[2,]~(c(1:length(summary(pcaControl_ZscoreByRegion)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC#", ylab="Proportion of Variance Explained", col=2)
dev.off()

png("Control_ZscoreByRegion_PCA_ScreePlot2.png")
plot(summary(pcaControl_ZscoreByRegion)$importance[3,]~(c(1:length(summary(pcaControl_ZscoreByRegion)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC#", ylab="Cumulative Proportion of Variance Explained", col=2)
dev.off()




#Consensus Clustering

#Normally we median center by gene before running the clustering. Does that actually make sense in this case?  (since the z-scores are already centered around 0 for the most part?)
#I skipped this step: MDD_Betas_OverlappingGenesMedianCentered<-sweep(MDD_Betas_OverlappingGenes, 1, apply(MDD_Betas_OverlappingGenes, 1, median, na.rm=T), FUN="-")

ClusterMDDAttempt1<-ConsensusClusterPlus(Control_ZscoreByRegion, maxK=6, reps=20, clusterAlg="km", title="ConsensusCluster_Control_ZscoreByRegion", distance="euclidean", plot="png", writeTable=TRUE)


#Both Jun and Wilkerson suggest choosing the top most variable genes (variable betas?)

MDD_Betas_OverlappingGenesStDev<-apply(MDD_Betas_OverlappingGenesMedianCentered, 1, sd)

png("Histogram_MDD_Betas_StDev.png")
hist(MDD_Betas_OverlappingGenesStDev)
dev.off()

quantile(MDD_Betas_OverlappingGenesStDev, probs=0.75)
# 75%  0.2374042 

MDD_Betas_OverlappingGenesMedianCenterMostVar<-MDD_Betas_OverlappingGenesMedianCentered[(MDD_Betas_OverlappingGenesStDev>quantile(MDD_Betas_OverlappingGenesStDev, probs=0.75)),]
dim(MDD_Betas_OverlappingGenesMedianCenterMostVar)
#[1] 3947   10


ClusterMDDAttempt3<-ConsensusClusterPlus(MDD_Betas_OverlappingGenesMedianCenterMostVar, maxK=6, reps=20, clusterAlg="km", title="ConsensusCluster_MDDAttempt2", distance="euclidean", plot="png", writeTable=TRUE)
#It seems to me like narrowing the analysis to just the top most variable genes simply produces a loss of information


 
 
 
 
 

         
         
         
                                
#Alright, it seems to me like the control-only data is still pretty noisy. I'm going to try using the full dataset and see how the z-scores relate, although diagnosis is going to influence the results then too (since the average looks at the center of the data instead of a defined intercept)... :(


#Adding a calculation of the average values and Z-scores for each region to the output:

All_AverageValuesByRegion<-matrix(NA, nrow=length(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm[,1]), ncol=10)
row.names(All_AverageValuesByRegion)<-row.names(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm)
colnames(All_AverageValuesByRegion)<-levels(AllRegions_Region)
head(All_AverageValuesByRegion)


temp<-tapply(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm[1,], AllRegions_Region, mean)
str(temp)

 for(i in 1:length(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm[,1])){

All_AverageValuesByRegion[i,]<-tapply(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm[i,], AllRegions_Region, mean)
}

head(All_AverageValuesByRegion)

All_AverageValuesFor10Regions<-apply(All_AverageValuesByRegion, 1, mean)


All_StDevByRegion<-matrix(NA, nrow=length(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm[,1]), ncol=10)
row.names(All_StDevByRegion)<-row.names(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm)
colnames(All_StDevByRegion)<-levels(AllRegions_Region)
head(All_StDevByRegion)


temp<-tapply(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm[1,], AllRegions_Region, sd)
str(temp)

 for(i in 1:length(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm[,1])){

All_StDevByRegion[i,]<-tapply(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm[i,], AllRegions_Region, sd)
}

head(All_StDevByRegion)


All_ZscoreByRegion<-matrix(NA, nrow=length(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm[,1]), ncol=10)
row.names(All_ZscoreByRegion)<-row.names(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm)
colnames(All_ZscoreByRegion)<-levels(AllRegions_Region)
head(All_ZscoreByRegion)




for(i in 1:length(AllRegions_NormNoOutliers_OverlappingGenes_ReNorm[,1])){

All_ZscoreByRegion[i,]<-(All_AverageValuesByRegion[i,]-All_AverageValuesFor10Regions[i])/All_StDevByRegion[i,]

}

#Outputting it:
Control_All_RepeatedMeasuresANOVAresults_Output<-cbind(Control_RepeatedMeasuresANOVAresults_Output, All_AverageValuesByRegion, All_AverageValuesFor10Regions, All_StDevByRegion, All_ZscoreByRegion)                                                                            

colnames(Control_All_RepeatedMeasuresANOVAresults_Output)


for(i in 1:10){
	j<-i+36
	colnames(Control_All_RepeatedMeasuresANOVAresults_Output)[j]<-paste("All AVE", levels(AllRegions_Region)[i], sep="_")
}

for(i in 1:10){
	j<-i+47
	colnames(Control_All_RepeatedMeasuresANOVAresults_Output)[j]<-paste("All STDEV", levels(AllRegions_Region)[i], sep="_")
}

for(i in 1:10){
	j<-i+57
	colnames(Control_All_RepeatedMeasuresANOVAresults_Output)[j]<-paste("All ZScore", levels(AllRegions_Region)[i], sep="_")
}

colnames(Control_All_RepeatedMeasuresANOVAresults_Output)


Control_All_RepeatedMeasuresANOVAresults_Output2<-join(as.data.frame(Control_All_RepeatedMeasuresANOVAresults_Output), ProbeInfoOutput, by="X")

write.csv(Control_All_RepeatedMeasuresANOVAresults_Output2, "Control_All_RepeatedMeasuresANOVAresults_Output2.csv")
 


png("CorrelationSummary_ZscoresByRegion_All.png", width = 1200, height = 1200)
pairs(All_ZscoreByRegion,
	gap=0,
	diag.panel=function(x,...){
		par(new=TRUE)
		hist(x,
			col="light blue",
			probability = TRUE,
			axes=FALSE,
			main="")
		lines(density(x), 
			col="red",
			lwd=3)
		rug(x)
		})
dev.off()


for(i in 1:10){
	
png(paste("Control vs All Z scores_", levels(AllRegions_Region)[i], ".png", sep=""))
j<-i+26
k<-i+57
plot(Control_All_RepeatedMeasuresANOVAresults_Output[,j]~Control_All_RepeatedMeasuresANOVAresults_Output[,k], xlab="All Z score", ylab="Control Z score", main=levels(AllRegions_Region)[i])
dev.off()
}

#Alright, looking at that it seems like using all of the data is definitely not providing a stronger indication of regional differences (when not accounting for diagnosis), it seems to be muddling the picture. 




#Note: http://jaredknowles.com/journal/2013/11/25/getting-started-with-mixed-effect-models-in-r.   Apparently instead of comparing the multilevel models to a simple lm it is better to use glm (generalized linear model) which produces a model fit through maximum likelihood estimation because glm has AIC in its output (a measure of model fit). The model formula specification is the same as for lm (although it looks like : may be used as an operator to specify interaction terms instead of *)  

#The lmer function with the familiar formula interface, but now group level variables are specified using a special syntax: e.g. (1|school) tells lmer to fit a linear model with a varying-intercept group effect using the variable school.
#we can fit nested group effect terms through the following syntax: lmer(extro ~ open + agree + social + (1 | school/class), data = lmm.data) Here the (1|school/class) says that we want to fit a mixed effect term for varying intercepts 1| by schools, and for classes that are nested within schools.
#But, what if we want to explore the effect of different student level indicators as they vary across classrooms. Instead of fitting unique models by school (or school/class) we can fit a varying slope model. Here we modify our random effect term to include variables before the grouping terms: (1 + open|school/class) tells R to fit a varying slope and varying intercept model for schools and classes nested within schools, and to allow the slope of the open variable to vary by school.


######1) At the same time, it might be nice to determine which of the genes are specific to each brain region (basically the opposite, but using an OR statement)

###C: Use Fisher's method to examine the consistency of p-values across regions (similar to the circadian paper) for diagnosis, as well as for various confounds potentially (age, gender, PMI, pH)

###D: Examine the consistency of the direction of the betas across the 10 brain regions for genes with consistently low p-values. (Average beta +/- SE?)

###D: Run brain region correlations for p-values and/or betas (maybe just for top genes???) and determine which regions are most similar (matrix of Rsquared or p-values for correlations?). Note: this analysis is pretty sensitive to outliers.

###E: Run PCA on the matrix of p-values and/or betas and determine if any of the brain regions cluster together.


#Upside: Avoids issues related to missing subjects in different brain regions as well as differing quantile normalization and confounds
#Downside: There will be a large number of genes excluded due to lack of detection in all 10 brain regions. We could potentially run this using individual pairs of regions of brain regions as well, similar to what is discussed in 2B.



#2: Looking for regional similarities across control subjects:

 
###A:For the Pritzker 960 dataset, there are different dissection effects across brain regions, so there is no way to really compare across brain regions without removing the influence of confounds first. Due to issues of limited range/sample size, I am inclined to apply this correction over the entire data set instead of limiting the correction to controls. The code for this can be found in the AAA analysis I did for Kevin, but basically it includes running a linear model with the major confounds and then subtracting out the effects of the confounds but not diagnosis. This may be a good place to start with Vikram's data too.


###B: For the first analysis, in order to maximize the number of genes included in regional comparisons, it would be nice to simply compare sample-sample correlations across regions.  Ideally, this will show 1) A nice difference in sample-sample correlations between regions, 2) A subject-dependent signature (although perhaps not, since we would have removed the influence of age, gender, and post-mortem factors). To do this, we would need to:

#####1)  Determine which probes are shared across both brain regions in the cleaned data (using code similar to 1.A.1 above), and perhaps at the same time determine which probes are exclusive to each brain region.

#####2)  Quantile normalize across the brain regions. This could be done over the entire data set or just over the controls.

#####3)  Calculate the sample-sample correlations and visualize them in a matrix.

#####4)  Determine whether there is a significant difference in sample-sample correlations for same brain region vs. different brain region (and the magnitude of this difference - Z-score?)

#####5)  Determine whether the sample-sample correlation across brain regions is on average higher for the same subject vs. different subjects.

#####6)  Run PCA and see if the brain regions segregate in PC1 or PC2.

#Upside: This may help identify brain region circuits (not just individual markers for regions)
#Downside:  Since the regions were run separately, it is possible that "more similar" regions may simply be regions that were run at similar times or with similar techniques.



###C: We could also try running a large-scale linear model. To do this, would need to:

#####1)  Determine which probes are shared across *all* brain regions in the cleaned data (using code identical to 1.A.1 above).

#####2)  Quantile normalize across all brain regions. This could be done over the entire data set or just over the controls.

#####3)  Replace any missing samples (outlier subjects) with median values for all genes.

#####4)  Run PCA and see if the brain regions segregate in PC1 or PC2.

#####5)  Run a linear model for each gene: y~Region1+Region2+Region3+Region4+Region5+Region6+Region7+Region8+Region9+Region10. Note: It may be worth it to run this for all subjects as well (not just controls) and add diagnosis as a term to identify genes that are consistently affected by diagnosis (and perhaps some region*diagnosis interaction terms?).  The p-values will be inflated for this analysis (since the data has been cleaned) but probably not by too much, since the df used in the model (10, maybe 14 with confounds) will be small relative to the number actually included in the analysis (960!)

#Upside: This is a powerful analysis.
#Downside: There will be a large number of genes excluded due to lack of detection in all 10 brain regions. Also, this analysis is better at identifying strong markers for individual brain regions than at identifying regional circuits.  ... and since the regions were run separately, this is likely to include technical artifacts. 




#Looking at cell type distribution across the different Amygdala nuclei


cbind(row.names(Control_RepeatedMeasuresANOVAresults), Control_RepeatedMeasuresANOVAresults, Control_AverageValuesByRegion, Control_AverageValuesFor10Regions, Control_StDevByRegion, Control_ZscoreByRegion

colnames(Control_All_RepeatedMeasuresANOVAresults_Output2)
#Col 68 is reannotated symbol
#Col 5-15 are control averages for each region
#Col 27-37 are control z-scores for each region

[1] "X"                                 "Df1"                               "Df2"                              
 [4] "F"                                 "p-value"                           "AVE_Lateral"                      
 [7] "AVE_AAA"                           "AVE_AB"                            "AVE_AHA"                          
[10] "AVE_Basal"                         "AVE_Central"                       "AVE_CO"                           
[13] "AVE_Medial"                        "AVE_PAC"                           "AVE_PL"                           
[16] "Control_AverageValuesFor10Regions" "STDEV_Lateral"                     "STDEV_AAA"                        
[19] "STDEV_AB"                          "STDEV_AHA"                         "STDEV_Basal"                      
[22] "STDEV_Central"                     "STDEV_CO"                          "STDEV_Medial"                     
[25] "STDEV_PAC"                         "STDEV_PL"                          "ZScore_Lateral"                   
[28] "ZScore_AAA"                        "ZScore_AB"                         "ZScore_AHA"                       
[31] "ZScore_Basal"                      "ZScore_Central"                    "ZScore_CO"                        
[34] "ZScore_Medial"                     "ZScore_PAC"                        "ZScore_PL"                        
[37] "All AVE_Lateral"                   "All AVE_AAA"                       "All AVE_AB"                       
[40] "All AVE_AHA"                       "All AVE_Basal"                     "All AVE_Central"                  
[43] "All AVE_CO"                        "All AVE_Medial"                    "All AVE_PAC"                      
[46] "All AVE_PL"                        "All_AverageValuesFor10Regions"     "All STDEV_Lateral"                
[49] "All STDEV_AAA"                     "All STDEV_AB"                      "All STDEV_AHA"                    
[52] "All STDEV_Basal"                   "All STDEV_Central"                 "All STDEV_CO"                     
[55] "All STDEV_Medial"                  "All STDEV_PAC"                     "All STDEV_PL"                     
[58] "All ZScore_Lateral"                "All ZScore_AAA"                    "All ZScore_AB"                    
[61] "All ZScore_AHA"                    "All ZScore_Basal"                  "All ZScore_Central"               
[64] "All ZScore_CO"                     "All ZScore_Medial"                 "All ZScore_PAC"                   
[67] "All ZScore_PL"                     "SYMBOLREANNOTATED"                 "ENTREZREANNOTATED"                
[70] "ProbeQuality"                      "CodingZone"                        "REPEATMASK"                       
[73] "OVERLAPPINGSNP"                    "CONTROLREPORTERGROUPNAME"          "X.1"   




OligodendrocyteGenes<-read.csv("OligodendrocyteGenes_CahoyZhang_forR.csv")
dim(OligodendrocyteGenes)
[1] 80  2

EndothelialGenes<-read.csv("EndothelialGenes_CahoyZhang_forR.csv")
dim(EndothelialGenes)
[1] 84  2

AstrocyteGenes<-read.csv("AstrocyteGenes_CahoyZhang_forR.csv")
dim(AstrocyteGenes)
[1] 96  2

NeuronGenes<-read.csv("NeuronGenes_CahoyZhang_forR.csv")
dim(NeuronGenes)
[1] 113   2

PericyteGenes<-read.csv("Zhang_PericyteGenes.csv")
dim(PericyteGenes)
[1] 40  1

MicrogliaGenes<-read.csv("Zhang_MicrogliaGenes.csv")
dim(MicrogliaGenes)
[1] 40  1

GlutamateNeuronGenes<-read.csv("Sugino_GlutamateNeuronGenes.csv")
dim(GlutamateNeuronGenes)
[1] 67  1

GABANeuronGenes<-read.csv("Sugino_GABANeuronGenes.csv")
dim(GABANeuronGenes)
[1] 32  1

RBCGenes<-read.csv("RBCGenes_GeneCard_forR.csv")
dim(RBCGenes)
[1] 17  1

OligodendrocyteControlSummary<-Control_All_RepeatedMeasuresANOVAresults_Output2[Control_All_RepeatedMeasuresANOVAresults_Output2[,68]%in%OligodendrocyteGenes[,1],]
dim(OligodendrocyteControlSummary)
[1] 66 75

EndothelialControlSummary<-Control_All_RepeatedMeasuresANOVAresults_Output2[Control_All_RepeatedMeasuresANOVAresults_Output2[,68]%in%EndothelialGenes[,1],]
dim(EndothelialControlSummary)
[1] 43 75

AstrocyteControlSummary<-Control_All_RepeatedMeasuresANOVAresults_Output2[Control_All_RepeatedMeasuresANOVAresults_Output2[,68]%in%AstrocyteGenes[,1],]
dim(AstrocyteControlSummary)
[1] 70 75

NeuronControlSummary<-Control_All_RepeatedMeasuresANOVAresults_Output2[Control_All_RepeatedMeasuresANOVAresults_Output2[,68]%in%NeuronGenes[,1],]
dim(NeuronControlSummary)
[1] 80 75

PericyteControlSummary<-Control_All_RepeatedMeasuresANOVAresults_Output2[Control_All_RepeatedMeasuresANOVAresults_Output2[,68]%in%PericyteGenes[,1],]
dim(PericyteControlSummary)
[1] 41 75

MicrogliaControlSummary<-Control_All_RepeatedMeasuresANOVAresults_Output2[Control_All_RepeatedMeasuresANOVAresults_Output2[,68]%in%MicrogliaGenes[,1],]
dim(MicrogliaControlSummary)
[1] 14 75

GlutamateNeuronControlSummary<-Control_All_RepeatedMeasuresANOVAresults_Output2[Control_All_RepeatedMeasuresANOVAresults_Output2[,68]%in%GlutamateNeuronGenes[,1],]
dim(GlutamateNeuronControlSummary)
[1] 58 75

GABANeuronControlSummary<-Control_All_RepeatedMeasuresANOVAresults_Output2[Control_All_RepeatedMeasuresANOVAresults_Output2[,68]%in%GABANeuronGenes[,1],]
dim(GABANeuronControlSummary)
[1] 31 75

RBCControlSummary<-Control_All_RepeatedMeasuresANOVAresults_Output2[Control_All_RepeatedMeasuresANOVAresults_Output2[,68]%in%RBCGenes[,1],]
dim(RBCControlSummary)
[1]  7 75


OligodendrocyteControlAverages<-matrix(data=as.numeric(as.matrix(OligodendrocyteControlSummary[,c(6:15)])), ncol=10)
colnames(OligodendrocyteControlAverages)<-colnames(OligodendrocyteControlSummary[,c(6:15)])
row.names(OligodendrocyteControlAverages)<-OligodendrocyteControlSummary[,68]
png("Heatmap_OligodendrocyteByRegion.csv")
heatmap(cor(t(OligodendrocyteControlAverages)))
dev.off()

OligodendrocyteIndex1<-apply(OligodendrocyteControlAverages, 2, mean)

OligodendrocyteControlAverages2<-matrix(data=as.numeric(as.matrix(OligodendrocyteControlSummary[,c(58:67)])), ncol=10)
colnames(OligodendrocyteControlAverages2)<-colnames(OligodendrocyteControlSummary[,c(58:67)])
row.names(OligodendrocyteControlAverages2)<-OligodendrocyteControlSummary[,68]
png("Heatmap_OligodendrocyteByRegion2.csv")
heatmap(cor(t(OligodendrocyteControlAverages2)))
dev.off()

OligodendrocyteIndex2<-apply(OligodendrocyteControlAverages2, 2, mean)



MicrogliaControlAverages<-matrix(data=as.numeric(as.matrix(MicrogliaControlSummary[,c(6:15)])), ncol=10)
colnames(MicrogliaControlAverages)<-colnames(MicrogliaControlSummary[,c(6:15)])
row.names(MicrogliaControlAverages)<-MicrogliaControlSummary[,68]
png("Heatmap_MicrogliaByRegion.csv")
heatmap(cor(t(MicrogliaControlAverages)))
dev.off()

MicrogliaIndex1<-apply(MicrogliaControlAverages, 2, mean)

MicrogliaControlAverages2<-matrix(data=as.numeric(as.matrix(MicrogliaControlSummary[,c(58:67)])), ncol=10)
colnames(MicrogliaControlAverages2)<-colnames(MicrogliaControlSummary[,c(58:67)])
row.names(MicrogliaControlAverages2)<-MicrogliaControlSummary[,68]
png("Heatmap_MicrogliaByRegion2.csv")
heatmap(cor(t(MicrogliaControlAverages2)))
dev.off()

MicrogliaIndex2<-apply(MicrogliaControlAverages2, 2, mean)



plot(OligodendrocyteIndex1~MicrogliaIndex1)
#Strong relationship - regions that have more oligodendrocytes have more microglia, like I found in the preliminary analysis.
plot(OligodendrocyteIndex2~MicrogliaIndex2)
#Not quite as pretty



EndothelialControlAverages<-matrix(data=as.numeric(as.matrix(EndothelialControlSummary[,c(6:15)])), ncol=10)
colnames(EndothelialControlAverages)<-colnames(EndothelialControlSummary[,c(6:15)])
row.names(EndothelialControlAverages)<-EndothelialControlSummary[,68]
png("Heatmap_EndothelialByRegion.csv")
heatmap(cor(t(EndothelialControlAverages)))
dev.off()

EndothelialIndex1<-apply(EndothelialControlAverages, 2, mean)

EndothelialControlAverages2<-matrix(data=as.numeric(as.matrix(EndothelialControlSummary[,c(58:67)])), ncol=10)
colnames(EndothelialControlAverages2)<-colnames(EndothelialControlSummary[,c(58:67)])
row.names(EndothelialControlAverages2)<-EndothelialControlSummary[,68]
png("Heatmap_EndothelialByRegion2.csv")
heatmap(cor(t(EndothelialControlAverages2)))
dev.off()

EndothelialIndex2<-apply(EndothelialControlAverages2, 2, mean)




AstrocyteControlAverages<-matrix(data=as.numeric(as.matrix(AstrocyteControlSummary[,c(6:15)])), ncol=10)
colnames(AstrocyteControlAverages)<-colnames(AstrocyteControlSummary[,c(6:15)])
row.names(AstrocyteControlAverages)<-AstrocyteControlSummary[,68]
png("Heatmap_AstrocyteByRegion.csv")
heatmap(cor(t(AstrocyteControlAverages)))
dev.off()

AstrocyteIndex1<-apply(AstrocyteControlAverages, 2, mean)

AstrocyteControlAverages2<-matrix(data=as.numeric(as.matrix(AstrocyteControlSummary[,c(58:67)])), ncol=10)
colnames(AstrocyteControlAverages2)<-colnames(AstrocyteControlSummary[,c(58:67)])
row.names(AstrocyteControlAverages2)<-AstrocyteControlSummary[,68]
png("Heatmap_AstrocyteByRegion2.csv")
heatmap(cor(t(AstrocyteControlAverages2)))
dev.off()

AstrocyteIndex2<-apply(AstrocyteControlAverages2, 2, mean)



NeuronControlAverages<-matrix(data=as.numeric(as.matrix(NeuronControlSummary[,c(6:15)])), ncol=10)
colnames(NeuronControlAverages)<-colnames(NeuronControlSummary[,c(6:15)])
row.names(NeuronControlAverages)<-NeuronControlSummary[,68]
png("Heatmap_NeuronByRegion.csv")
heatmap(cor(t(NeuronControlAverages)))
dev.off()

NeuronIndex1<-apply(NeuronControlAverages, 2, mean)

NeuronControlAverages2<-matrix(data=as.numeric(as.matrix(NeuronControlSummary[,c(58:67)])), ncol=10)
colnames(NeuronControlAverages2)<-colnames(NeuronControlSummary[,c(58:67)])
row.names(NeuronControlAverages2)<-NeuronControlSummary[,68]
png("Heatmap_NeuronByRegion2.csv")
heatmap(cor(t(NeuronControlAverages2)))
dev.off()

NeuronIndex2<-apply(NeuronControlAverages2, 2, mean)



PericyteControlAverages<-matrix(data=as.numeric(as.matrix(PericyteControlSummary[,c(6:15)])), ncol=10)
colnames(PericyteControlAverages)<-colnames(PericyteControlSummary[,c(6:15)])
row.names(PericyteControlAverages)<-PericyteControlSummary[,68]
png("Heatmap_PericyteByRegion.csv")
heatmap(cor(t(PericyteControlAverages)))
dev.off()

PericyteIndex1<-apply(PericyteControlAverages, 2, mean)

PericyteControlAverages2<-matrix(data=as.numeric(as.matrix(PericyteControlSummary[,c(58:67)])), ncol=10)
colnames(PericyteControlAverages2)<-colnames(PericyteControlSummary[,c(58:67)])
row.names(PericyteControlAverages2)<-PericyteControlSummary[,68]
png("Heatmap_PericyteByRegion2.csv")
heatmap(cor(t(PericyteControlAverages2)))
dev.off()

PericyteIndex2<-apply(PericyteControlAverages2, 2, mean)



RBCControlAverages<-matrix(data=as.numeric(as.matrix(RBCControlSummary[,c(6:15)])), ncol=10)
colnames(RBCControlAverages)<-colnames(RBCControlSummary[,c(6:15)])
row.names(RBCControlAverages)<-RBCControlSummary[,68]
png("Heatmap_RBCByRegion.csv")
heatmap(cor(t(RBCControlAverages)))
dev.off()

RBCIndex1<-apply(RBCControlAverages, 2, mean)

RBCControlAverages2<-matrix(data=as.numeric(as.matrix(RBCControlSummary[,c(58:67)])), ncol=10)
colnames(RBCControlAverages2)<-colnames(RBCControlSummary[,c(58:67)])
row.names(RBCControlAverages2)<-RBCControlSummary[,68]
png("Heatmap_RBCByRegion2.csv")
heatmap(cor(t(RBCControlAverages2)))
dev.off()

RBCIndex2<-apply(RBCControlAverages2, 2, mean)



GlutamateNeuronControlAverages<-matrix(data=as.numeric(as.matrix(GlutamateNeuronControlSummary[,c(6:15)])), ncol=10)
colnames(GlutamateNeuronControlAverages)<-colnames(GlutamateNeuronControlSummary[,c(6:15)])
row.names(GlutamateNeuronControlAverages)<-GlutamateNeuronControlSummary[,68]
png("Heatmap_GlutamateNeuronByRegion.csv")
heatmap(cor(t(GlutamateNeuronControlAverages)))
dev.off()

GlutamateNeuronIndex1<-apply(GlutamateNeuronControlAverages, 2, mean)

GlutamateNeuronControlAverages2<-matrix(data=as.numeric(as.matrix(GlutamateNeuronControlSummary[,c(58:67)])), ncol=10)
colnames(GlutamateNeuronControlAverages2)<-colnames(GlutamateNeuronControlSummary[,c(58:67)])
row.names(GlutamateNeuronControlAverages2)<-GlutamateNeuronControlSummary[,68]
png("Heatmap_GlutamateNeuronByRegion2.csv")
heatmap(cor(t(GlutamateNeuronControlAverages2)))
dev.off()

GlutamateNeuronIndex2<-apply(GlutamateNeuronControlAverages2, 2, mean)



GABANeuronControlAverages<-matrix(data=as.numeric(as.matrix(GABANeuronControlSummary[,c(6:15)])), ncol=10)
colnames(GABANeuronControlAverages)<-colnames(GABANeuronControlSummary[,c(6:15)])
row.names(GABANeuronControlAverages)<-GABANeuronControlSummary[,68]
png("Heatmap_GABANeuronByRegion.csv")
heatmap(cor(t(GABANeuronControlAverages)))
dev.off()

GABANeuronIndex1<-apply(GABANeuronControlAverages, 2, mean)

GAD1Signal<-Control_All_RepeatedMeasuresANOVAresults_Output2[Control_All_RepeatedMeasuresANOVAresults_Output2[,68]%in%"GAD1", c(6:15)]
GAD2Signal<-Control_All_RepeatedMeasuresANOVAresults_Output2[Control_All_RepeatedMeasuresANOVAresults_Output2[,68]%in%"GAD2", c(6:15)]
GADSignal<-rbind(GAD1Signal, GAD2Signal)
GADIndex1<-apply(GADSignal, 2, mean)

plot(GADIndex1~GABANeuronIndex1)


GABANeuronControlAverages2<-matrix(data=as.numeric(as.matrix(GABANeuronControlSummary[,c(58:67)])), ncol=10)
colnames(GABANeuronControlAverages2)<-colnames(GABANeuronControlSummary[,c(58:67)])
row.names(GABANeuronControlAverages2)<-GABANeuronControlSummary[,68]
png("Heatmap_GABANeuronByRegion2.csv")
heatmap(cor(t(GABANeuronControlAverages2)))
dev.off()

GABANeuronIndex2<-apply(GABANeuronControlAverages2, 2, mean)

Control_All_RepeatedMeasuresANOVAresults_Output2[Control_All_RepeatedMeasuresANOVAresults_Output2[,68]=="GAD1",]



CellTypeIndices1<-rbind(NeuronIndex1, GABANeuronIndex1, GlutamateNeuronIndex1, AstrocyteIndex1, RBCIndex1, PericyteIndex1, EndothelialIndex1, MicrogliaIndex1, OligodendrocyteIndex1)

write.csv(CellTypeIndices1, "CellTypeIndices1.csv")

CellTypeIndices2<-rbind(NeuronIndex2, GABANeuronIndex2, GlutamateNeuronIndex2, AstrocyteIndex2, RBCIndex2, PericyteIndex2, EndothelialIndex2, MicrogliaIndex2, OligodendrocyteIndex2)


write.csv(CellTypeIndices2, "CellTypeIndices2.csv")

CellTypeIndices1_CorrelationMatrixbyCellType<-cor(t(CellTypeIndices1))

write.csv(CellTypeIndices1_CorrelationMatrixbyCellType, "CellTypeIndices1_CorrelationMatrixbyCellType.csv")


png("Heatmap_CellTypeIndices1_CorrMatrixbyCellType.png")
heatmap(cor(t(CellTypeIndices1)))
dev.off()

png("Scatterplots_CellTypeIndices1_CorrMatrixbyCellType.png", width = 1200, height = 1200)
pairs(t(CellTypeIndices1),
	gap=0,
	diag.panel=function(x,...){
		par(new=TRUE)
		hist(x,
			col="light blue",
			probability = TRUE,
			axes=FALSE,
			main="")
		lines(density(x), 
			col="red",
			lwd=3)
		rug(x)
		})
dev.off()


pca_CellTypeIndices1byCellType<-prcomp(t(CellTypeIndices1))
tmp<-pca_CellTypeIndices1byCellType$x[,1:9]
rownames(tmp)<-colnames(CellTypeIndices1)
write.csv(tmp, "pca_CellTypeIndices1byCellType.csv")

tmp<-pca_CellTypeIndices1byCellType$rotation[,1:9]
write.csv(tmp, "pca_CellTypeIndices1byCellType_Eigenvectors.csv")

png("pca_CellTypeIndices1byCellType_PC1vsPC2.png")
plot(pca_CellTypeIndices1byCellType$x[,1]~pca_CellTypeIndices1byCellType$x[,2], xlab="PC2", ylab="PC1", main="PC1 vs. PC2")
dev.off()

png("pca_CellTypeIndices1byCellType_ScreePlot.png")
plot(summary(pca_CellTypeIndices1byCellType)$importance[2,]~(c(1:length(summary(pca_CellTypeIndices1byCellType)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC#", ylab="Proportion of Variance Explained", col=2)
dev.off()



CellTypeIndices1_CorrelationMatrixbyRegion<-cor(CellTypeIndices1)

write.csv(CellTypeIndices1_CorrelationMatrixbyRegion, "CellTypeIndices1_CorrelationMatrixbyRegion.csv")


png("Heatmap_CellTypeIndices1_CorrMatrixbyRegion.png")
heatmap(cor(CellTypeIndices1))
dev.off()


png("Scatterplots_CellTypeIndices1_CorrMatrixbyRegion.png", width = 1200, height = 1200)
pairs(CellTypeIndices1,
	gap=0,
	diag.panel=function(x,...){
		par(new=TRUE)
		hist(x,
			col="light blue",
			probability = TRUE,
			axes=FALSE,
			main="")
		lines(density(x), 
			col="red",
			lwd=3)
		rug(x)
		})
dev.off()


pca_CellTypeIndices1byRegion<-prcomp(CellTypeIndices1)
tmp<-pca_CellTypeIndices1byRegion$x[,1:10]
rownames(tmp)<-rownames(CellTypeIndices1)
write.csv(tmp, "pca_CellTypeIndices1byRegion.csv")

tmp<-pca_CellTypeIndices1byRegion$rotation[,1:9]
write.csv(tmp, "pca_CellTypeIndices1byRegion_Eigenvectors.csv")

png("pca_CellTypeIndices1byRegion_PC1vsPC2.png")
plot(pca_CellTypeIndices1byRegion$x[,1]~pca_CellTypeIndices1byRegion$x[,2], xlab="PC2", ylab="PC1", main="PC1 vs. PC2")
dev.off()

png("pca_CellTypeIndices1byRegion_ScreePlot.png")
plot(summary(pca_CellTypeIndices1byRegion)$importance[2,]~(c(1:length(summary(pca_CellTypeIndices1byRegion)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC#", ylab="Proportion of Variance Explained", col=2)
dev.off()
#something is wrong - I should probably wait to debug until tomorrow when not tired.







CellTypeIndices2_CorrelationMatrixbyCellType<-cor(t(CellTypeIndices2))

write.csv(CellTypeIndices2_CorrelationMatrixbyCellType, "CellTypeIndices1_CorrelationMatrixbyCellType.csv")



png("Heatmap_CellTypeIndices2_CorrMatrixbyCellType.png")
heatmap(cor(t(CellTypeIndices2)))
dev.off()

CellTypeIndices2_CorrelationMatrixbyRegion<-cor(CellTypeIndices2)

write.csv(CellTypeIndices2_CorrelationMatrixbyRegion, "CellTypeIndices1_CorrelationMatrixbyRegion.csv")


png("Heatmap_CellTypeIndices2_CorrMatrixbyRegion.png")
heatmap(cor(CellTypeIndices2))
dev.off()



#Pulling out cell type specific genes from the MDD results:

dim(MDD_Summary_CrossRegionalResultsv2)
[1] 15788    36

colnames(MDD_Summary_CrossRegionalResultsv2)
[29] "SYMBOLREANNOTATED"


MDD_MeanStatsByNucleus<-apply(matrix(as.numeric(MDD_Summary_CrossRegionalResults[, c(2:27)]), nrow=15788, ncol=26), 2, mean)
names(MDD_MeanStatsByNucleus)<-colnames(MDD_Summary_CrossRegionalResults[, c(2:27)])

MDD_MedianStatsByNucleus<-apply(matrix(as.numeric(MDD_Summary_CrossRegionalResults[, c(2:27)]), nrow=15788, ncol=26), 2, median)
names(MDD_MedianStatsByNucleus)<-colnames(MDD_Summary_CrossRegionalResults[, c(2:27)])

MDD_1stQuantileStatsByNucleus<-apply(matrix(as.numeric(MDD_Summary_CrossRegionalResults[, c(2:27)]), nrow=15788, ncol=26), 2, function(x) quantile(x, 0.25))
names(MDD_1stQuantileStatsByNucleus)<-colnames(MDD_Summary_CrossRegionalResults[, c(2:27)])

MDD_10QuantileStatsByNucleus<-apply(matrix(as.numeric(MDD_Summary_CrossRegionalResults[, c(2:27)]), nrow=15788, ncol=26), 2, function(x) quantile(x, 0.10))
names(MDD_10QuantileStatsByNucleus)<-colnames(MDD_Summary_CrossRegionalResults[, c(2:27)])

temp<-rbind(MDD_MeanStatsByNucleus, MDD_MedianStatsByNucleus, MDD_1stQuantileStatsByNucleus, MDD_10QuantileStatsByNucleus)

write.csv(temp, "MDD_StatsDistributionByNucleus.csv")

 [1] "X"                        "AAA"                      "AB"                       "AHA"                     
 [5] "Basal"                    "Central"                  "CO"                       "Lateral"                 
 [9] "Medial"                   "PAC"                      "PL"                       "Median Pvalue"           
[13] "Mean Log Pvalue"          "Chi-Square Stat"          "Fishers Pvalue"           "AAA"                     
[17] "AB"                       "AHA"                      "Basal"                    "Central"                 
[21] "CO"                       "Lateral"                  "Medial"                   "PAC"                     
[25] "PL"                       "Median Beta"              "Mean Beta"                "SE Beta" 


MDD_CellTypeMeanStatsByNucleus<-matrix(0, nrow=9, ncol=26)
colnames(MDD_CellTypeMeanStatsByNucleus)<-colnames(MDD_Summary_CrossRegionalResults[, c(2:27)])
row.names(MDD_CellTypeMeanStatsByNucleus)<-c("Oligodendrocyte", "Endothelial", "Astrocyte", "CorticalNeuron", "Pericyte", "Microglia", "CorticalGlutamate", "CorticalGABA", "RBC")

MDD_CellTypeMedianStatsByNucleus<-matrix(0, nrow=9, ncol=26)
colnames(MDD_CellTypeMedianStatsByNucleus)<-colnames(MDD_Summary_CrossRegionalResults[, c(2:27)])
row.names(MDD_CellTypeMedianStatsByNucleus)<-c("Oligodendrocyte", "Endothelial", "Astrocyte", "CorticalNeuron", "Pericyte", "Microglia", "CorticalGlutamate", "CorticalGABA", "RBC")


OligodendrocyteMDDSummary<-MDD_Summary_CrossRegionalResultsv2[MDD_Summary_CrossRegionalResultsv2[,29]%in%OligodendrocyteGenes[,1],]
dim(OligodendrocyteMDDSummary)
[1] 66 36

write.csv(OligodendrocyteMDDSummary, "OligodendrocyteMDDSummary.csv")


EndothelialMDDSummary<-MDD_Summary_CrossRegionalResultsv2[MDD_Summary_CrossRegionalResultsv2[,29]%in%EndothelialGenes[,1],]
dim(EndothelialMDDSummary)
[1] 43 36

write.csv(EndothelialMDDSummary, "EndothelialMDDSummary.csv")

AstrocyteMDDSummary<-MDD_Summary_CrossRegionalResultsv2[MDD_Summary_CrossRegionalResultsv2[,29]%in%AstrocyteGenes[,1],]
dim(AstrocyteMDDSummary)
[1] 70 36

write.csv(AstrocyteMDDSummary, "AstrocyteMDDSummary.csv")

NeuronMDDSummary<-MDD_Summary_CrossRegionalResultsv2[MDD_Summary_CrossRegionalResultsv2[,29]%in%NeuronGenes[,1],]
dim(NeuronMDDSummary)
[1] 80 36

write.csv(NeuronMDDSummary, "NeuronMDDSummary.csv")

PericyteMDDSummary<-MDD_Summary_CrossRegionalResultsv2[MDD_Summary_CrossRegionalResultsv2[,29]%in%PericyteGenes[,1],]
dim(PericyteMDDSummary)
[1] 41 36

write.csv(PericyteMDDSummary, "PericyteMDDSummary.csv")

MicrogliaMDDSummary<-MDD_Summary_CrossRegionalResultsv2[MDD_Summary_CrossRegionalResultsv2[,29]%in%MicrogliaGenes[,1],]
dim(MicrogliaMDDSummary)
[1] 14 36

write.csv(MicrogliaMDDSummary, "MicrogliaMDDSummary.csv")

GlutamateNeuronMDDSummary<-MDD_Summary_CrossRegionalResultsv2[MDD_Summary_CrossRegionalResultsv2[,29]%in%GlutamateNeuronGenes[,1],]
dim(GlutamateNeuronMDDSummary)
[1] 58 36

write.csv(GlutamateNeuronMDDSummary, "GlutamateNeuronMDDSummary.csv")

GABANeuronMDDSummary<-MDD_Summary_CrossRegionalResultsv2[MDD_Summary_CrossRegionalResultsv2[,29]%in%GABANeuronGenes[,1],]
dim(GABANeuronMDDSummary)
[1] 31 36

write.csv(GABANeuronMDDSummary, "GABANeuronMDDSummary.csv")

RBCMDDSummary<-MDD_Summary_CrossRegionalResultsv2[MDD_Summary_CrossRegionalResultsv2[,29]%in%RBCGenes[,1],]
dim(RBCMDDSummary)
[1]  7 36

write.csv(RBCMDDSummary, "RBCMDDSummary.csv")



     
      

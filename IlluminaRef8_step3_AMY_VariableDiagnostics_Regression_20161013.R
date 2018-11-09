library(gdata)
library(fields)
library(stats)
library(car)
library(affy)
library(preprocessCore)
library(multtest)
library(pylr)
library(lmPerm)


#read in the demographic data and the expression signal data (which has been median centered, checked for outliers, checked for gender incongruities, and had missing subjects removed):

setwd("/Users/eaurbach/Dropbox/University of Michigan/Akil Lab/Postmortem Human Microarray/Freeze3/Freeze3 Data Reworked/All Regions - Variable Diagnostics-Regression-CleanSignalofConfounds-InteractionsDiagnosis/Amy")

AmySubjDemographics<-as.data.frame(read.csv("Amy_SubjectDemographics_noNAs_OrderedBySubject.csv",header=T))
row.names(AmySubjDemographics)<-AmySubjDemographics$subject_id
AmySubjDemographics<-AmySubjDemographics[,-1] #weird read-in thing where it gave us an extra column of subject IDs
AmyDetectedSignalData<-as.matrix(read.csv("Amy_PvalDetected_MedianCenteredData_noNAs_OrderedBySubject.csv",header=T,row.names=1))
colnames(AmyDetectedSignalData) <- as.numeric(gsub("[^0-9]", "", colnames(AmyDetectedSignalData)))
AmyUnfilteredSignalData<-as.matrix(read.csv("Amy_Unfiltered_MedianCenteredData_noNAs_OrderedBySubject.csv",header=T,row.names=1))
colnames(AmyUnfilteredSignalData) <- as.numeric(gsub("[^0-9]", "", colnames(AmyUnfilteredSignalData)))
colnames(AmyUnfilteredSignalData) <- as.numeric(gsub("[^0-9]", "", colnames(AmyUnfilteredSignalData))) #remove the X from the column names and make them numeric for ease of matching
colnames(AmyDetectedSignalData) == row.names(AmySubjDemographics) # retruns a matrix of TRUEs
colnames(AmyUnfilteredSignalData) == row.names(AmySubjDemographics) # retruns a matrix of TRUEs

#make sure to have data structures for the annotation linking probe ID and gene symbol
GeneNamesAmydetected<-as.data.frame(read.csv("AmydetectedIlluminaDecodeSorted22177Abbreviated.csv",header=T,row.names=16))
row.names(GeneNamesAmydetected)==row.names(AmyDetectedSignalData) #all true
GeneNamesAmyunfiltered<-as.data.frame(read.csv("AmyUnfilteredIlluminaDecodeSorted22177Abbreviated.csv",header=T,row.names=16))
row.names(GeneNamesAmyunfiltered)==row.names(AmyUnfilteredSignalData) #all true


#extract the relevant subject variables to examine during diagnostics and use in the "basic" linear model:
Diagnosis <- as.factor(AmySubjDemographics$Disease)
Diagnosis <- relevel(Diagnosis, ref = "C")
str(Diagnosis)

AgonalFactor <- as.numeric(AmySubjDemographics$agonal_factor_d002_b003)

BrainPH <- as.numeric(AmySubjDemographics$ph_a002_b002) 

Gender <- as.factor(AmySubjDemographics$gender_f002)
Gender <- relevel(Gender, ref = "M")
str(Gender)

Age <- as.numeric(AmySubjDemographics$age_f003)

PMI <- as.numeric(AmySubjDemographics$f016_hoursfinal)


SuicideIndices<-which(AmySubjDemographics$f021_mannerofdeath=="Suicide")
NotSuicideIndices<-which(AmySubjDemographics$f021_mannerofdeath!="Suicide")
Suicide<-matrix(nrow=length(row.names(AmySubjDemographics)))
Suicide[SuicideIndices]<-1
Suicide[NotSuicideIndices]<-0



#*********** Characterizing subject demographics and looking for relationships between subject variables:******

FactorVariables <- cbind(Gender, Diagnosis) #note: using cbind to make a matrix out of factor variables eliminates the labels associated with the levels to replace them with the numerical values--so this is now a matrix of 0s, 1s, 2s, 3s, and 4s. In order to re-establish the labels with gender and diagnosis, I have to make variables to identify them.
ContinuousVariables <- cbind(BrainPH, AgonalFactor, PMI, Age)
row.names(ContinuousVariables) <- row.names(AmySubjDemographics) #since these variables were numeric, they didn't have row names automatically and were just numbered.
row.names(FactorVariables)<- row.names(AmySubjDemographics)

#make histograms of the continuous variables and store them in a sub-folder called "Continuous Variable Histograms"
Amy_mainWD <- getwd()
Amy_hist <- paste(Amy_mainWD, "/Continuous Variable Histograms", sep = "")
dir.create(Amy_hist)
setwd(Amy_hist)
for (i in 1:length(colnames(ContinuousVariables))) {
	png(paste(paste("Histogram of", colnames(ContinuousVariables)[i], sep = " "), ".png", sep = ""))
	hist(ContinuousVariables[, i], col = i + 1, main = paste("Histogram of", colnames(ContinuousVariables)[i], sep = " "), xlab = paste(colnames(ContinuousVariables)[i]))
	dev.off()
}
setwd(Amy_mainWD)


#make scatter plots of the continuous variables and store them in a sub-folder called "Continuous Variable Scatterplots"
Amy_scatter <- paste(Amy_mainWD, "/Continuous Variable Scatterplots", sep = "")
dir.create(Amy_scatter)
setwd(Amy_scatter)
for (i in 1:length(colnames(ContinuousVariables))) {
	for (j in 1:length(colnames(ContinuousVariables))) {
		png(paste(paste("Scatterplot of", colnames(ContinuousVariables)[i], "vs.", colnames(ContinuousVariables)[j], sep = " "), 
			".png", sep = ""))

		plot.default(ContinuousVariables[, i], ContinuousVariables[, j], main = paste("Scatterplot of", colnames(ContinuousVariables)[i], 
			"vs.", colnames(ContinuousVariables)[j], sep = " "), xlab = paste(colnames(ContinuousVariables)[i]), ylab = paste(colnames(ContinuousVariables)[j]), 
			col = "maroon")
		RegressionLine <- lm(ContinuousVariables[, j] ~ ContinuousVariables[, i])
		abline(RegressionLine, col = "black")
		mtext(paste("p-value=", round(summary.lm(RegressionLine)$coefficients[8], digits = 4)))

		dev.off()
	}
}
setwd(Amy_mainWD)


#Using boxplots to visually examine the relationships between the continuous subject variables and categorical subject variables:
Amy_box <- paste(Amy_mainWD, "/Continuous vs Categorical Variable Boxplots", sep = "")
dir.create(Amy_box)
setwd(Amy_box)
for (i in 1:length(colnames(ContinuousVariables))) {
	for (j in 1:length(colnames(FactorVariables))) {
		png(paste(paste("Boxplot of", colnames(ContinuousVariables)[i], "vs.", colnames(FactorVariables)[j], sep = " "), ".png", 
			sep = ""))

		if (colnames(FactorVariables)[j] == "Gender") {
			boxplot(ContinuousVariables[, i] ~ FactorVariables[, j], main = paste("Boxplot of", colnames(ContinuousVariables)[i], "vs.", 
				colnames(FactorVariables)[j], sep = " "), xlab = paste(colnames(FactorVariables)[j]), ylab = paste(colnames(ContinuousVariables)[i]), 
				names = levels(Gender))
				
			RegressionLine <- lm(ContinuousVariables[, i] ~ Gender)
			mtext(paste("Female p-value=", round(summary.lm(RegressionLine)$coefficients[2, 4], digits = 4)))

			
		} else if (colnames(FactorVariables)[j] == "Diagnosis") {
			boxplot(ContinuousVariables[, i] ~ FactorVariables[, j], main = paste("Boxplot of", colnames(ContinuousVariables)[i], "vs.", 
				colnames(FactorVariables)[j], sep = " "), xlab = paste(colnames(FactorVariables)[j]), ylab = paste(colnames(ContinuousVariables)[i]), 
				names = levels(Diagnosis))

			RegressionLine <- lm(ContinuousVariables[, i] ~ Diagnosis)
			mtext(paste(paste(paste("BP p-value=", round(summary.lm(RegressionLine)$coefficients[2, 4], digits = 4)),"MD p-value=", round(summary.lm(RegressionLine)$coefficients[3, 4], digits = 4),sep = "   "),"Schiz p-value=", round(summary.lm(RegressionLine)$coefficients[4, 4], digits = 4),sep = "   "))
			
		} else if (colnames(FactorVariables)[j] == "Suicide") {
			boxplot(ContinuousVariables[, i] ~ FactorVariables[, j], main = paste("Boxplot of", colnames(ContinuousVariables)[i], "vs.", 
				colnames(FactorVariables)[j], sep = " "), xlab = paste(colnames(FactorVariables)[j]), ylab = paste(colnames(ContinuousVariables)[i]), 
				names = cbind("No Suicide", "Suicide"))

			RegressionLine <- lm(ContinuousVariables[, i] ~ Suicide)
			mtext(paste("Suicide p-value=", round(summary.lm(RegressionLine)$coefficients[2, 4], digits = 4)))

		} else {
			print("Something's wrong with the contents of FactorVariables")
			stop()
		}

		dev.off()
	}
}
setwd(Amy_mainWD)


###Outputting statistical relationships between Continuous and Factor Variables
#note: this code is stolen straight from Megan

Amy_subjectstats <- paste(Amy_mainWD, "/Stats between Subject Variables", sep = "")
dir.create(Amy_subjectstats)
setwd(Amy_subjectstats)

#Creating a text file of contingency tables to visually examine the relationships between categorical subject variables:
CrossTabsIV<-file("Cross Tabs Between Subject Factors.txt")
out<-c(
capture.output(

summary(Diagnosis),
summary(Gender),

for (i in 1:length(FactorVariables[1,])){
for(j in 1:length(FactorVariables[1,])){
ContingencyTable<-table(FactorVariables[,i],FactorVariables[,j])
print(paste(colnames(FactorVariables)[i], "vs", colnames(FactorVariables)[j], sep="  "))
print(ContingencyTable)
print(paste("p-value=", chisq.test(ContingencyTable)$p.value))	
}		
}
)
)
cat(out, file="Cross Tabs Between Subject Factors.txt", sep="\n", append=TRUE)
close(CrossTabsIV)
rm(out)





StatisticalRelationshipsIV<-file("Statistical Relationships between Subject Variables.txt")
out<-c(

capture.output(
#Calculating the variance inflation factor (vif) to determine which subject variables are highly related to other subject variables in the data set. Most important, of course, is whether any of the subject variables strongly correlate with Diagnosis. 
vif(lm(AmyDetectedSignalData[1,]~BrainPH + AgonalFactor +PMI+Diagnosis+Gender+Age))



#It might be nice to run a hierarchical cluster amongst these variables also...
),

#Using linear regression to examine the statistical relationships between the continuous subject variables:

capture.output(
for (i in 1:length(ContinuousVariables[1,])){
for(j in 1:length(ContinuousVariables[1,])){
print(paste(colnames(ContinuousVariables)[i], "vs", colnames(ContinuousVariables)[j], sep="  "))
print(summary.lm(lm(ContinuousVariables[,i]~ContinuousVariables[,j])))
}		
}
),

#Using anova to examine the statistical relationships between the continuous subject variables and categorical subject variables:

capture.output(
for (i in 1:length(ContinuousVariables[1,])){
for(j in 1:length(FactorVariables[1,])){
print(paste(colnames(ContinuousVariables)[i], "vs", colnames(FactorVariables)[j], sep="  "))
print(summary(aov(ContinuousVariables[,i]~FactorVariables[,j])))		
}		
}
),


#Using chi-square to examine the statistical relationships between the categorical subject variables:

capture.output(
for (i in 1:length(FactorVariables[1,])){
for(j in 1:length(FactorVariables[1,])){
ContingencyTable<-table(FactorVariables[,i],FactorVariables[,j])
print(paste(colnames(FactorVariables)[i], "vs", colnames(FactorVariables)[j], sep="  "))
print(chisq.test(ContingencyTable))		
}		
}
)

)
cat(out, file="Statistical Relationships between Subject Variables.txt", sep="\n", append=TRUE)
close(StatisticalRelationshipsIV)
rm(out)


#Flagging variables that are collinear with other subject variables (p<0.05 for simple bivariate relationships):
FlaggedRelationshipsBetweenIV<-file("Flagged Relationships Between Subject Variables.txt")
out<-c(

#Using linear regression to examine the statistical relationships between the continuous subject variables:
capture.output(
for (i in 1:length(ContinuousVariables[1,])){
for(j in 1:length(ContinuousVariables[1,])){
if(summary.lm(lm(ContinuousVariables[,i]~ContinuousVariables[,j]))$coefficient[8]<0.05){
print(paste(colnames(ContinuousVariables)[i], "vs", colnames(ContinuousVariables)[j], "p-value=", summary.lm(lm(ContinuousVariables[,i]~ContinuousVariables[,j]))$coefficient[8], sep="  "))}else{}
}		
}
),

#Using anova to examine the statistical relationships between the continuous subject variables and categorical subject variables:
capture.output(
for (i in 1:length(ContinuousVariables[1,])){
for(j in 1:length(FactorVariables[1,])){
if(summary(aov(ContinuousVariables[,i]~FactorVariables[,j]))[[1]][["Pr(>F)"]][1]<0.05){
print(paste(colnames(ContinuousVariables)[i], "vs", colnames(FactorVariables)[j], "p-value=", summary(aov(ContinuousVariables[,i]~FactorVariables[,j]))[[1]][["Pr(>F)"]][1], sep="  "))	
}else{}		
}		
}
),

#Using chi-square to examine the statistical relationships between the categorical subject variables:
capture.output(
for (i in 1:length(FactorVariables[1,])){
for(j in 1:length(FactorVariables[1,])){
ContingencyTable<-table(FactorVariables[,i], FactorVariables[,j])
if(chisq.test(ContingencyTable)$p.value<0.05){
print(paste(colnames(FactorVariables)[i], "vs", colnames(FactorVariables)[j], "p-value=", chisq.test(ContingencyTable)$p.value, sep="  "))	
}else{}		
}		
}
)
)
cat(out, file="Flagged Relationships Between Subject Variables.txt", sep="\n", append=TRUE)
close(FlaggedRelationshipsBetweenIV)
rm(out)


setwd(Amy_mainWD)

Amy_descstats <- paste(Amy_mainWD, "/Descriptive Statistics for Important Covariates", sep = "")
dir.create(Amy_descstats)
setwd(Amy_descstats)


DescriptiveStats<-file("Descriptive Stats for Important Covariates.txt")
out<-c(print("sample sizes by diagnosis:"),
capture.output(
table(AmySubjDemographics$Disease)
),print("number of males/females by diagnosis:"),
capture.output(print("sums (C, BP, MD, SC)"),
table(AmySubjDemographics$gender_f002[AmySubjDemographics$Disease=="C"]),
table(AmySubjDemographics$gender_f002[AmySubjDemographics$Disease=="BP"]),
table(AmySubjDemographics$gender_f002[AmySubjDemographics$Disease=="MD"]),
table(AmySubjDemographics$gender_f002[AmySubjDemographics$Disease=="SC"])
),print("mean, median, and standard deviation of age by diagnosis:"),
capture.output(
print("mean (C, BP, MD, SC)"),
mean(AmySubjDemographics$age_f003[AmySubjDemographics$Disease=="C"]),
mean(AmySubjDemographics$age_f003[AmySubjDemographics$Disease=="BP"]),
mean(AmySubjDemographics$age_f003[AmySubjDemographics$Disease=="MD"]),
mean(AmySubjDemographics$age_f003[AmySubjDemographics$Disease=="SC"]),
print("median (C, BP, MD, SC)"),
median(AmySubjDemographics$age_f003[AmySubjDemographics$Disease=="C"]),
median(AmySubjDemographics$age_f003[AmySubjDemographics$Disease=="BP"]),
median(AmySubjDemographics$age_f003[AmySubjDemographics$Disease=="MD"]),
median(AmySubjDemographics$age_f003[AmySubjDemographics$Disease=="SC"]),
print("sd (C, BP, MD, SC)"),
sd(AmySubjDemographics$age_f003[AmySubjDemographics$Disease=="C"]),
sd(AmySubjDemographics$age_f003[AmySubjDemographics$Disease=="BP"]),
sd(AmySubjDemographics$age_f003[AmySubjDemographics$Disease=="MD"]),
sd(AmySubjDemographics$age_f003[AmySubjDemographics$Disease=="SC"])
),print("mean, median, and standard deviation of PMI by diagnosis:"),
capture.output(
print("mean (C, BP, MD, SC)"),
mean(AmySubjDemographics$f016_hoursfinal[AmySubjDemographics$Disease=="C"]),
mean(AmySubjDemographics$f016_hoursfinal[AmySubjDemographics$Disease=="BP"]),
mean(AmySubjDemographics$f016_hoursfinal[AmySubjDemographics$Disease=="MD"]),
mean(AmySubjDemographics$f016_hoursfinal[AmySubjDemographics$Disease=="SC"]),
print("median (C, BP, MD, SC)"),
median(AmySubjDemographics$f016_hoursfinal[AmySubjDemographics$Disease=="C"]),
median(AmySubjDemographics$f016_hoursfinal[AmySubjDemographics$Disease=="BP"]),
median(AmySubjDemographics$f016_hoursfinal[AmySubjDemographics$Disease=="MD"]),
median(AmySubjDemographics$f016_hoursfinal[AmySubjDemographics$Disease=="SC"]),
print("sd (C, BP, MD, SC)"),
sd(AmySubjDemographics$f016_hoursfinal[AmySubjDemographics$Disease=="C"]),
sd(AmySubjDemographics$f016_hoursfinal[AmySubjDemographics$Disease=="BP"]),
sd(AmySubjDemographics$f016_hoursfinal[AmySubjDemographics$Disease=="MD"]),
sd(AmySubjDemographics$f016_hoursfinal[AmySubjDemographics$Disease=="SC"])
),print("mean, median, and standard deviation of agonal factor by diagnosis:"),
capture.output(
print("mean (C, BP, MD, SC)"),
mean(AmySubjDemographics$agonal_factor_d002_b003[AmySubjDemographics$Disease=="C"]),
mean(AmySubjDemographics$agonal_factor_d002_b003[AmySubjDemographics$Disease=="BP"]),
mean(AmySubjDemographics$agonal_factor_d002_b003[AmySubjDemographics$Disease=="MD"]),
mean(AmySubjDemographics$agonal_factor_d002_b003[AmySubjDemographics$Disease=="SC"]),
print("median (C, BP, MD, SC)"),
median(AmySubjDemographics$agonal_factor_d002_b003[AmySubjDemographics$Disease=="C"]),
median(AmySubjDemographics$agonal_factor_d002_b003[AmySubjDemographics$Disease=="BP"]),
median(AmySubjDemographics$agonal_factor_d002_b003[AmySubjDemographics$Disease=="MD"]),
median(AmySubjDemographics$agonal_factor_d002_b003[AmySubjDemographics$Disease=="SC"]),
print("sd (C, BP, MD, SC)"),
sd(AmySubjDemographics$agonal_factor_d002_b003[AmySubjDemographics$Disease=="C"]),
sd(AmySubjDemographics$agonal_factor_d002_b003[AmySubjDemographics$Disease=="BP"]),
sd(AmySubjDemographics$agonal_factor_d002_b003[AmySubjDemographics$Disease=="MD"]),
sd(AmySubjDemographics$agonal_factor_d002_b003[AmySubjDemographics$Disease=="SC"])
),print("mean, median, and standard deviation of pH by diagnosis:"),
capture.output(
print("mean (C, BP, MD, SC)"),
mean(AmySubjDemographics$ph_a002_b002[AmySubjDemographics$Disease=="C"]),
mean(AmySubjDemographics$ph_a002_b002[AmySubjDemographics$Disease=="BP"]),
mean(AmySubjDemographics$ph_a002_b002[AmySubjDemographics$Disease=="MD"]),
mean(AmySubjDemographics$ph_a002_b002[AmySubjDemographics$Disease=="SC"]),
print("median (C, BP, MD, SC)"),
median(AmySubjDemographics$ph_a002_b002[AmySubjDemographics$Disease=="C"]),
median(AmySubjDemographics$ph_a002_b002[AmySubjDemographics$Disease=="BP"]),
median(AmySubjDemographics$ph_a002_b002[AmySubjDemographics$Disease=="MD"]),
median(AmySubjDemographics$ph_a002_b002[AmySubjDemographics$Disease=="SC"]),
print("sd (C, BP, MD, SC)"),
sd(AmySubjDemographics$ph_a002_b002[AmySubjDemographics$Disease=="C"]),
sd(AmySubjDemographics$ph_a002_b002[AmySubjDemographics$Disease=="BP"]),
sd(AmySubjDemographics$ph_a002_b002[AmySubjDemographics$Disease=="MD"]),
sd(AmySubjDemographics$ph_a002_b002[AmySubjDemographics$Disease=="SC"])
)
)
cat(out, file="Descriptive Stats for Important Covariates.txt", sep="\n", append=TRUE)
close(DescriptiveStats)
rm(out)

mean(AmySubjDemographics$age_f003[AmySubjDemographics$Disease=="MD"])





#********************UNFILTERED Characterizing the principal components of variation (PCA) in the microarray signal data********
#note: again, code largely stolen from Megan

Amy_pca <- paste(Amy_mainWD, "/Principal Components Analysis", sep = "")
dir.create(Amy_pca)
setwd(Amy_pca)

Amy_pca_unf <- paste(Amy_pca, "/Unfiltered", sep = "")
dir.create(Amy_pca_unf)
setwd(Amy_pca_unf)

# # #Visualize the sample-sample correlations using a heatmap:
png("Sample Sample Correlations Heatmap.png")
image(cor(AmyUnfilteredSignalData), main="Visualizing the correlations between entire samples (by index #)", xlab="Red=Less correlated, Light yellow=Highly correlated")
dev.off()
# #Note that the heatmap can be tailored to focus on a certain level of correlation by using the command zlim=c(lower limit, upper limit)

#note that this generates a basically unreadable boxplot bc everything is on top of itself...
png("Boxplot Sample Sample Correlations.png")
boxplot(data.frame(cor(AmyUnfilteredSignalData)), cex=0.25, las=3, par(cex.axis=0.75))
Median10thQuantile<-median(apply((cor(AmyUnfilteredSignalData)), 1, quantile, 0.1))
abline(a=Median10thQuantile, b=0, col=2)
dev.off()

plot(AmyUnfilteredSignalData[,1]~ AmyUnfilteredSignalData[,2])
#Note that the normalized signal data (which is all median centered around 0 for each gene) has made it so that the sample-sample correlations are sometimes even negative!

#run the PCA itself
PCA_AmyUnfilteredSignalData<-prcomp(t(AmyUnfilteredSignalData))
writefile<-PCA_AmyUnfilteredSignalData $x[,1:4]
write.table(writefile, "PCA_1_4.txt", sep="\t")

PC1<-PCA_AmyUnfilteredSignalData $x[,1] #note: ask Megan why she originally named these variables with a "noOutliers" designation
PC2<-PCA_AmyUnfilteredSignalData $x[,2]
PC3<-PCA_AmyUnfilteredSignalData $x[,3]
PC4<-PCA_AmyUnfilteredSignalData $x[,4]

png("PCA Scree Plot - proportion variance explained.png")
plot(summary(PCA_AmyUnfilteredSignalData)$importance[2,]~(c(1:length(summary(PCA_AmyUnfilteredSignalData)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC #", ylab="Proportion of Variance Explained", col=2)
dev.off()

png("PCA Scree Plot - cumulative proportion of variance explained.png")
plot(summary(PCA_AmyUnfilteredSignalData)$importance[3,]~(c(1:length(summary(PCA_AmyUnfilteredSignalData)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC #", ylab="Cumulative Proportion of Variance Explained", col=3)
dev.off()
#ask Megan: what is the difference here?

png("PC1 vs PC2.png")
plot(PC1~PC2, main="Principal Components Analysis - PC1 vs PC2")
points(PC1[Diagnosis=="Control"]~PC2[Diagnosis =="Control"], col="black")
points(PC1[Diagnosis =="MD"]~PC2[Diagnosis =="MD"], col="magenta")
points(PC1[Diagnosis =="BP"]~PC2[Diagnosis =="BP"], col="blue")
points(PC1[Diagnosis =="Schiz"]~PC2[Diagnosis =="Schiz"], col="green")
#legend(min(PC1noOutliers), max(PC2noOutliers)+10, c("Control", "MDD", "BP", "Schiz"), text.col=c(3, 2, 4, 5), pch=19, col=c(3, 2, 4, 5)) #this line isn't working and I don't feel like figuring it out
dev.off()

# #Output a scatterplot illustrating the relationship between Principal components 3 & 4 (PC3 & PC4) and diagnosis:
png("PC3 vs PC3.png")
plot(PC3~PC4, main="Principal Components Analysis - PC3 vs PC4")
points(PC3[Diagnosis=="Control"]~PC4[Diagnosis =="Control"], col="black")
points(PC3[Diagnosis =="MD"]~PC4[Diagnosis =="MD"], col="magenta")
points(PC3[Diagnosis =="BP"]~PC4[Diagnosis =="BP"], col="blue")
points(PC3[Diagnosis =="Schiz"]~PC4[Diagnosis =="Schiz"], col="green")
#legend(min(PC1noOutliers), max(PC2noOutliers)+10, c("Control", "MDD", "BP", "Schiz"), text.col=c(3, 2, 4, 5), pch=19, col=c(3, 2, 4, 5)) #this line isn't working and I don't feel like figuring it out
dev.off()


SubjectPCA<-cbind(PC1, PC2, PC3, PC4)

PCAoutput<-cbind(FactorVariables, ContinuousVariables, SubjectPCA)
write.csv(PCAoutput, "PCAoutput.csv")

#********************#8. Characterizing potential sources of noise by and examining relationships between subject variables and principal components of variation (PCA) in the microarray signal data**************************
#note: as before, code mostly stolen/adapted from Megan

Amy_pca_pcofv <- paste(Amy_pca_unf, "/Characterizing Principal Components of Variation", sep = "")
dir.create(Amy_pca_pcofv)
setwd(Amy_pca_pcofv)


#Using a scatterplot with best fit line to visually examine the relationships between the continuous subject variables and SubjectPCA:

for (i in 1:length(SubjectPCA[1,])){
for(j in 1:length(ContinuousVariables[1,])){
png(paste(paste(colnames(SubjectPCA)[i], "vs", colnames(ContinuousVariables)[j], sep="  "), "png", sep="."))	
plot(SubjectPCA[,i]~ContinuousVariables[,j], main=paste(colnames(SubjectPCA)[i], "vs", colnames(ContinuousVariables)[j], sep="  "), xlab=colnames(ContinuousVariables)[j], ylab=colnames(SubjectPCA)[i])
RegressionLine<-lm(SubjectPCA[,i]~ContinuousVariables[,j])
abline(RegressionLine, col=2)
mtext(paste("p-value=", round(summary.lm(RegressionLine)$coefficients[8], digits=4)))
dev.off()		
}		
}

#Using boxplots to visually examine the relationships between the PCA and categorical subject variables:
for (i in 1:length(SubjectPCA[1,])){
for(j in 1:length(FactorVariables[1,])){
png(paste(paste(colnames(SubjectPCA)[i], "vs", colnames(FactorVariables)[j], sep="  "), "png", sep="."))	
boxplot(SubjectPCA[,i]~FactorVariables[,j], main=paste(colnames(SubjectPCA)[i], "vs", colnames(FactorVariables)[j], sep="  "), xlab=colnames(FactorVariables)[j], ylab=colnames(SubjectPCA)[i])
mtext(paste("p-value=", round(summary(aov(SubjectPCA[,i]~FactorVariables[,j]))[[1]][["Pr(>F)"]][1], digits=4)))
dev.off()		
}		
}


#Outputting a text file containing the statistical relationships between all of the subject variables and PCA:

StatisticalRelationshipsIVvsPCA<-file("Statistical Relationships between Subject Variables and PCA.txt")
out<-c(

capture.output(
summary.lm(lm(PC1~BrainPH + AgonalFactor +PMI+Diagnosis+Gender+Age+Suicide))
),

capture.output(
summary.lm(lm(PC2~BrainPH + AgonalFactor +PMI+Diagnosis+Gender+Age+Suicide))
),

capture.output(
summary.lm(lm(PC3~BrainPH + AgonalFactor +PMI+Diagnosis+Gender+Age+Suicide))
),

capture.output(
summary.lm(lm(PC4~BrainPH + AgonalFactor +PMI+Diagnosis+Gender+Age+Suicide))
),




#Using linear regression to examine the statistical relationships between PCA and the continuous subject variables:

capture.output(
for (i in 1:length(SubjectPCA[1,])){
for(j in 1:length(ContinuousVariables[1,])){
print(paste(colnames(SubjectPCA)[i], "vs", colnames(ContinuousVariables)[j], sep="  "))
print(summary.lm(lm(SubjectPCA[,i]~ContinuousVariables[,j])))
}		
}
),

#Using anova to examine the statistical relationships between PCA and categorical subject variables:

capture.output(
for (i in 1:length(SubjectPCA[1,])){
for(j in 1:length(FactorVariables[1,])){
print(paste(colnames(SubjectPCA)[i], "vs", colnames(FactorVariables)[j], sep="  "))
print(summary(aov(SubjectPCA[,i]~FactorVariables[,j])))		
}		
}
)

)
cat(out, file="Statistical Relationships between Subject Variables and PCA.txt", sep="\n", append=TRUE)
close(StatisticalRelationshipsIVvsPCA)
rm(out)

#Flagging PC that are collinear with subject variables (p<0.05 in simple bivariate model):

FlaggedRelationshipsBetweenIVandPCA<-file("Flagged Relationships Between Subject Variables and PCA.txt")
out<-c(

#Using linear regression to examine the statistical relationships between the continuous subject variables:
capture.output(
for (i in 1:length(SubjectPCA[1,])){
for(j in 1:length(ContinuousVariables[1,])){
if(summary.lm(lm(SubjectPCA[,i]~ContinuousVariables[,j]))$coefficient[8]<0.05){
print(paste(colnames(SubjectPCA)[i], "vs", colnames(ContinuousVariables)[j], "p-value=", summary.lm(lm(SubjectPCA[,i]~ContinuousVariables[,j]))$coefficient[8], sep="  "))}else{}
}		
}
),

#Using anova to examine the statistical relationships between the continuous subject variables and categorical subject variables:
capture.output(
for (i in 1:length(SubjectPCA[1,])){
for(j in 1:length(FactorVariables[1,])){
if(summary(aov(SubjectPCA[,i]~FactorVariables[,j]))[[1]][["Pr(>F)"]][1]<0.05){
print(paste(colnames(SubjectPCA)[i], "vs", colnames(FactorVariables)[j], "p-value=", summary(aov(SubjectPCA[,i]~FactorVariables[,j]))[[1]][["Pr(>F)"]][1], sep="  "))	
}else{}		
}		
}
)
)
cat(out, file="Flagged Relationships Between Subject Variables and PCA.txt", sep="\n", append=TRUE)
close(FlaggedRelationshipsBetweenIVandPCA)
rm(out)

setwd(Amy_pca)


#********************FILTERED by DETECTION P-VAL Characterizing the principal components of variation (PCA) in the microarray signal data********
#note: again, code largely stolen from Megan

Amy_pca_det <- paste(Amy_pca, "/Detected", sep = "")
dir.create(Amy_pca_det)
setwd(Amy_pca_det)


# # #Visualize the sample-sample correlations using a heatmap:
png("Sample Sample Correlations Heatmap.png")
image(cor(AmyDetectedSignalData), main="Visualizing the correlations between entire samples (by index #)", xlab="Red=Less correlated, Light yellow=Highly correlated")
dev.off()
# #Note that the heatmap can be tailored to focus on a certain level of correlation by using the command zlim=c(lower limit, upper limit)

#note that this generates a basically unreadable boxplot bc everything is on top of itself...
png("Boxplot Sample Sample Correlations.png")
boxplot(data.frame(cor(AmyDetectedSignalData)), cex=0.25, las=3, par(cex.axis=0.75))
Median10thQuantile<-median(apply((cor(AmyDetectedSignalData)), 1, quantile, 0.1))
abline(a=Median10thQuantile, b=0, col=2)
dev.off()

plot(AmyDetectedSignalData[,1]~ AmyDetectedSignalData[,2])
#Note that the normalized signal data (which is all median centered around 0 for each gene) has made it so that the sample-sample correlations are sometimes even negative!

#run the PCA itself
PCA_AmyDetectedSignalData<-prcomp(t(AmyDetectedSignalData))
writefile<-PCA_AmyDetectedSignalData $x[,1:4]
write.table(writefile, "PCA_1_4.txt", sep="\t")

PC1<-PCA_AmyDetectedSignalData $x[,1] #note: ask Megan why she originally named these variables with a "noOutliers" designation
PC2<-PCA_AmyDetectedSignalData $x[,2]
PC3<-PCA_AmyDetectedSignalData $x[,3]
PC4<-PCA_AmyDetectedSignalData $x[,4]

png("PCA Scree Plot - proportion variance explained.png")
plot(summary(PCA_AmyDetectedSignalData)$importance[2,]~(c(1:length(summary(PCA_AmyDetectedSignalData)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC #", ylab="Proportion of Variance Explained", col=2)
dev.off()

png("PCA Scree Plot - cumulative proportion of variance explained.png")
plot(summary(PCA_AmyDetectedSignalData)$importance[3,]~(c(1:length(summary(PCA_AmyDetectedSignalData)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC #", ylab="Cumulative Proportion of Variance Explained", col=3)
dev.off()
#ask Megan: what is the difference here?

png("PC1 vs PC2.png")
plot(PC1~PC2, main="Principal Components Analysis - PC1 vs PC2")
points(PC1[Diagnosis=="Control"]~PC2[Diagnosis =="Control"], col="black")
points(PC1[Diagnosis =="MD"]~PC2[Diagnosis =="MD"], col="magenta")
points(PC1[Diagnosis =="BP"]~PC2[Diagnosis =="BP"], col="blue")
points(PC1[Diagnosis =="Schiz"]~PC2[Diagnosis =="Schiz"], col="green")
#legend(min(PC1noOutliers), max(PC2noOutliers)+10, c("Control", "MDD", "BP", "Schiz"), text.col=c(3, 2, 4, 5), pch=19, col=c(3, 2, 4, 5)) #this line isn't working and I don't feel like figuring it out
dev.off()

# #Output a scatterplot illustrating the relationship between Principal components 3 & 4 (PC3 & PC4) and diagnosis:
png("PC3 vs PC3.png")
plot(PC3~PC4, main="Principal Components Analysis - PC3 vs PC4")
points(PC3[Diagnosis=="Control"]~PC4[Diagnosis =="Control"], col="black")
points(PC3[Diagnosis =="MD"]~PC4[Diagnosis =="MD"], col="magenta")
points(PC3[Diagnosis =="BP"]~PC4[Diagnosis =="BP"], col="blue")
points(PC3[Diagnosis =="Schiz"]~PC4[Diagnosis =="Schiz"], col="green")
#legend(min(PC1noOutliers), max(PC2noOutliers)+10, c("Control", "MDD", "BP", "Schiz"), text.col=c(3, 2, 4, 5), pch=19, col=c(3, 2, 4, 5)) #this line isn't working and I don't feel like figuring it out
dev.off()


SubjectPCA<-cbind(PC1, PC2, PC3, PC4)

PCAoutput<-cbind(FactorVariables, ContinuousVariables, SubjectPCA)
write.csv(PCAoutput, "PCAoutput.csv")

#********************#8. Characterizing potential sources of noise by and examining relationships between subject variables and principal components of variation (PCA) in the microarray signal data**************************
#note: as before, code mostly stolen/adapted from Megan

Amy_pca_pcofv <- paste(Amy_pca_det, "/Characterizing Principal Components of Variation", sep = "")
dir.create(Amy_pca_pcofv)
setwd(Amy_pca_pcofv)


#Using a scatterplot with best fit line to visually examine the relationships between the continuous subject variables and SubjectPCA:

for (i in 1:length(SubjectPCA[1,])){
for(j in 1:length(ContinuousVariables[1,])){
png(paste(paste(colnames(SubjectPCA)[i], "vs", colnames(ContinuousVariables)[j], sep="  "), "png", sep="."))	
plot(SubjectPCA[,i]~ContinuousVariables[,j], main=paste(colnames(SubjectPCA)[i], "vs", colnames(ContinuousVariables)[j], sep="  "), xlab=colnames(ContinuousVariables)[j], ylab=colnames(SubjectPCA)[i])
RegressionLine<-lm(SubjectPCA[,i]~ContinuousVariables[,j])
abline(RegressionLine, col=2)
mtext(paste("p-value=", round(summary.lm(RegressionLine)$coefficients[8], digits=4)))
dev.off()		
}		
}

#Using boxplots to visually examine the relationships between the PCA and categorical subject variables:
for (i in 1:length(SubjectPCA[1,])){
for(j in 1:length(FactorVariables[1,])){
png(paste(paste(colnames(SubjectPCA)[i], "vs", colnames(FactorVariables)[j], sep="  "), "png", sep="."))	
boxplot(SubjectPCA[,i]~FactorVariables[,j], main=paste(colnames(SubjectPCA)[i], "vs", colnames(FactorVariables)[j], sep="  "), xlab=colnames(FactorVariables)[j], ylab=colnames(SubjectPCA)[i])
mtext(paste("p-value=", round(summary(aov(SubjectPCA[,i]~FactorVariables[,j]))[[1]][["Pr(>F)"]][1], digits=4)))
dev.off()		
}		
}


#Outputting a text file containing the statistical relationships between all of the subject variables and PCA:

StatisticalRelationshipsIVvsPCA<-file("Statistical Relationships between Subject Variables and PCA.txt")
out<-c(

capture.output(
summary.lm(lm(PC1~BrainPH + AgonalFactor +PMI+Diagnosis+Gender+Age+Suicide))
),

capture.output(
summary.lm(lm(PC2~BrainPH + AgonalFactor +PMI+Diagnosis+Gender+Age+Suicide))
),

capture.output(
summary.lm(lm(PC3~BrainPH + AgonalFactor +PMI+Diagnosis+Gender+Age+Suicide))
),

capture.output(
summary.lm(lm(PC4~BrainPH + AgonalFactor +PMI+Diagnosis+Gender+Age+Suicide))
),




#Using linear regression to examine the statistical relationships between PCA and the continuous subject variables:

capture.output(
for (i in 1:length(SubjectPCA[1,])){
for(j in 1:length(ContinuousVariables[1,])){
print(paste(colnames(SubjectPCA)[i], "vs", colnames(ContinuousVariables)[j], sep="  "))
print(summary.lm(lm(SubjectPCA[,i]~ContinuousVariables[,j])))
}		
}
),

#Using anova to examine the statistical relationships between PCA and categorical subject variables:

capture.output(
for (i in 1:length(SubjectPCA[1,])){
for(j in 1:length(FactorVariables[1,])){
print(paste(colnames(SubjectPCA)[i], "vs", colnames(FactorVariables)[j], sep="  "))
print(summary(aov(SubjectPCA[,i]~FactorVariables[,j])))		
}		
}
)

)
cat(out, file="Statistical Relationships between Subject Variables and PCA.txt", sep="\n", append=TRUE)
close(StatisticalRelationshipsIVvsPCA)
rm(out)

#Flagging PC that are collinear with subject variables (p<0.05 in simple bivariate model):

FlaggedRelationshipsBetweenIVandPCA<-file("Flagged Relationships Between Subject Variables and PCA.txt")
out<-c(

#Using linear regression to examine the statistical relationships between the continuous subject variables:
capture.output(
for (i in 1:length(SubjectPCA[1,])){
for(j in 1:length(ContinuousVariables[1,])){
if(summary.lm(lm(SubjectPCA[,i]~ContinuousVariables[,j]))$coefficient[8]<0.05){
print(paste(colnames(SubjectPCA)[i], "vs", colnames(ContinuousVariables)[j], "p-value=", summary.lm(lm(SubjectPCA[,i]~ContinuousVariables[,j]))$coefficient[8], sep="  "))}else{}
}		
}
),

#Using anova to examine the statistical relationships between the continuous subject variables and categorical subject variables:
capture.output(
for (i in 1:length(SubjectPCA[1,])){
for(j in 1:length(FactorVariables[1,])){
if(summary(aov(SubjectPCA[,i]~FactorVariables[,j]))[[1]][["Pr(>F)"]][1]<0.05){
print(paste(colnames(SubjectPCA)[i], "vs", colnames(FactorVariables)[j], "p-value=", summary(aov(SubjectPCA[,i]~FactorVariables[,j]))[[1]][["Pr(>F)"]][1], sep="  "))	
}else{}		
}		
}
)
)
cat(out, file="Flagged Relationships Between Subject Variables and PCA.txt", sep="\n", append=TRUE)
close(FlaggedRelationshipsBetweenIVandPCA)
rm(out)

setwd(Amy_mainWD)

#********************#9. Running a Linear Model to examine the effects of diagnosis on the expression of all probes (including the FGF family) while accounting for important confounding variables (BrainPH, AgonalFactor, HoursFinal, Gender, Age)****

LinearModel_MeganBasic<-function(i){lm(AmyUnfilteredSignalData[i,]~Diagnosis+BrainPH+AgonalFactor +PMI+Gender+Age)}

NumberofXs<-length(summary.lm(LinearModel_MeganBasic(1))$coefficients[,1])
NameofXs<-dimnames(summary.lm(LinearModel_MeganBasic(1))$coefficients)[1][[1]]

ModelBetas<-matrix(0, length(AmyUnfilteredSignalData[,1]), NumberofXs)
colnames(ModelBetas)<-NameofXs
row.names(ModelBetas)<-row.names(AmyUnfilteredSignalData)

Modelpvalues<-matrix(0, length(AmyUnfilteredSignalData[,1]), NumberofXs)
colnames(Modelpvalues)<-NameofXs
row.names(Modelpvalues)<-row.names(AmyUnfilteredSignalData)

ModelSE<-matrix(0, length(AmyUnfilteredSignalData[,1]), NumberofXs)
colnames(ModelSE)<-NameofXs
row.names(ModelSE)<-row.names(AmyUnfilteredSignalData)

ModelTstat<-matrix(0, length(AmyUnfilteredSignalData[,1]), NumberofXs)
colnames(ModelTstat)<-NameofXs
row.names(ModelTstat)<-row.names(AmyUnfilteredSignalData)

#Running the linear model:
for(i in 1:length(AmyUnfilteredSignalData[,1])){
	TempModel<-LinearModel_MeganBasic(i)
	ModelBetas[i,]<-summary.lm(TempModel)$coefficients[,1]
	ModelSE[i,]<-summary.lm(TempModel)$coefficients[,2]
	ModelTstat[i,]<-summary.lm(TempModel)$coefficients[,3]
	Modelpvalues[i,]<-summary.lm(TempModel)$coefficients[,4]
}

#Outputting a histogram illustrating the raw p-values in comparison to what would be expected by chance:
Amy_pval_hist <- paste(Amy_mainWD, "/Histograms of p-values", sep = "")
dir.create(Amy_pval_hist)
setwd(Amy_pval_hist)

for (i in 1: NumberofXs){
png(paste(paste("Histogram of Raw Pvalues Using Megan's Linear Model for", NameofXs[i], sep="  "), "png", sep="."))	
hist(Modelpvalues[,i], breaks=100, col=i, main=paste("Raw P-values using Megan's Linear Model for", NameofXs[i], sep="  "), xlab="Unadjusted p-value", ylab="Count")
abline(a=(length(AmyUnfilteredSignalData[,1])/100), b=0)
dev.off()		
}		#megan says this is mostly to look at data to think through Benjamini Hotchberg mult comp correction, because small changes in many genes (population level) can artificially enhance effects

setwd(Amy_mainWD)

Amy_model <- paste(Amy_mainWD, "/Linear Model", sep = "")
dir.create(Amy_model)
setwd(Amy_model)

ModelBetasOutput<-cbind(row.names(GeneNamesAmyunfiltered),as.character(GeneNamesAmyunfiltered$Symbol),ModelBetas) 
ModelSEOutput<-cbind(row.names(GeneNamesAmyunfiltered),as.character(GeneNamesAmyunfiltered$Symbol),ModelSE)
ModelTstatOutput<-cbind(row.names(GeneNamesAmyunfiltered),as.character(GeneNamesAmyunfiltered$Symbol),ModelTstat)
ModelpvaluesOutput<-cbind(row.names(GeneNamesAmyunfiltered),as.character(GeneNamesAmyunfiltered$Symbol),Modelpvalues)

write.csv(ModelBetasOutput, "ModelBetas.csv")
write.csv(ModelSEOutput, "ModelSE.csv")
write.csv(ModelTstatOutput, "ModelTstat.csv")
write.csv(ModelpvaluesOutput, "ModelpvaluesRAW.csv")

#output the data for MD, BP, and Schiz results
#BP
TempPvalAdjBP<-mt.rawp2adjp(Modelpvalues[,2], proc=c("BH", "BY"))
ModelPvalAdjBP<-TempPvalAdjBP$adjp[order(TempPvalAdjBP$index),]
row.names(ModelPvalAdjBP)<-row.names(AmyUnfilteredSignalData)
colnames(ModelPvalAdjBP)<-c("BP rawp","BP BH-Corrected","BP BY-Corrected")
#MD
TempPvalAdjMD<-mt.rawp2adjp(Modelpvalues[,3], proc=c("BH", "BY"))
ModelPvalAdjMD<-TempPvalAdjMD$adjp[order(TempPvalAdjMD$index),]
row.names(ModelPvalAdjMD)<-row.names(AmyUnfilteredSignalData)
colnames(ModelPvalAdjMD)<-c("MD rawp","MD BH-Corrected","MD BY-Corrected")
#Schiz
TempPvalAdjSchiz<-mt.rawp2adjp(Modelpvalues[,4], proc=c("BH", "BY"))
ModelPvalAdjSchiz<-TempPvalAdjSchiz$adjp[order(TempPvalAdjSchiz$index),]
row.names(ModelPvalAdjSchiz)<-row.names(AmyUnfilteredSignalData)
colnames(ModelPvalAdjSchiz)<-c("Schiz rawp","Schiz BH-Corrected","Schiz BY-Corrected")


ModelOutput<-cbind(row.names(GeneNamesAmyunfiltered),as.character(GeneNamesAmyunfiltered$Symbol),ModelBetas[,2], ModelPvalAdjBP,ModelBetas[,3], ModelPvalAdjMD,ModelBetas[,4], ModelPvalAdjSchiz)
row.names(ModelOutput)<-row.names(AmyUnfilteredSignalData)
colnames(ModelOutput)<-c("Probe/Array Address ID","Gene Symbol","BP Betas","BP unadjusted p-value","BP BH-adjusted p-value","BP BY-adjusted p-value","MD Betas","MD unadjusted p-value","MD BH-adjusted p-value","MD BY-adjusted p-value","Schiz Betas","Schiz unadjusted p-value","Schiz BH-adjusted p-value","Schiz BY-adjusted p-value")

write.csv(ModelOutput, "Amy_Freeze3Unfiltered_MegansModel_DiagnosisOutput_WithMultCompCorrections.csv")

###E, 9/1 need to re-run this code filtered by detection p-val. So I will pull out the approproate probes from the row names of AmyDetectedSignalData and extract the appropriate subset to re-run multiple comparisons correction and output

#Extract only the detected probes for Modelpvalues and ModelBetas in order to rerun multiple comparisons correction:
ModelBetasDetected<-ModelBetas[row.names(AmyUnfilteredSignalData)%in%row.names(AmyDetectedSignalData),]
ModelpvaluesDetected<-Modelpvalues[row.names(AmyUnfilteredSignalData)%in%row.names(AmyDetectedSignalData),]

#output the data for MD, BP, and Schiz results
#BP
TempPvalAdjBPDet<-mt.rawp2adjp(ModelpvaluesDetected[,2], proc=c("BH", "BY"))
ModelPvalAdjBPDet<-TempPvalAdjBPDet$adjp[order(TempPvalAdjBPDet$index),]
row.names(ModelPvalAdjBPDet)<-row.names(AmyDetectedSignalData)
colnames(ModelPvalAdjBPDet)<-c("BP rawp","BP BH-Corrected","BP BY-Corrected")
#MD
TempPvalAdjMDDet<-mt.rawp2adjp(ModelpvaluesDetected[,3], proc=c("BH", "BY"))
ModelPvalAdjMDDet<-TempPvalAdjMDDet$adjp[order(TempPvalAdjMDDet$index),]
row.names(ModelPvalAdjMDDet)<-row.names(AmyDetectedSignalData)
colnames(ModelPvalAdjMDDet)<-c("MD rawp","MD BH-Corrected","MD BY-Corrected")
#Schiz
TempPvalAdjSchizDet<-mt.rawp2adjp(ModelpvaluesDetected[,4], proc=c("BH", "BY"))
ModelPvalAdjSchizDet<-TempPvalAdjSchizDet$adjp[order(TempPvalAdjSchizDet$index),]
row.names(ModelPvalAdjSchizDet)<-row.names(AmyDetectedSignalData)
colnames(ModelPvalAdjSchizDet)<-c("Schiz rawp","Schiz BH-Corrected","Schiz BY-Corrected")


ModelOutputDetected<-cbind(row.names(AmyDetectedSignalData),as.character(GeneNamesAmydetected$Symbol), ModelBetasDetected[,2], ModelPvalAdjBPDet, ModelBetasDetected[,3], ModelPvalAdjMDDet, ModelBetasDetected[,4], ModelPvalAdjSchizDet)
row.names(ModelOutputDetected)<-row.names(AmyDetectedSignalData)
colnames(ModelOutputDetected)<-c("Probe/Array Address ID","Gene Symbol","BP Betas","BP unadjusted p-value","BP BH-adjusted p-value","BP BY-adjusted p-value","MD Betas","MD unadjusted p-value","MD BH-adjusted p-value","MD BY-adjusted p-value","Schiz Betas","Schiz unadjusted p-value","Schiz BH-adjusted p-value","Schiz BY-adjusted p-value")

write.csv(ModelOutputDetected, "Amy_Freeze3Detected_MegansModel_DiagnosisOutput_WithMultCompCorrections.csv")


#************#10. Cleaning the signal data of the effects of confounds before making figures *********************

median(BrainPH)
#[1] 6.955

median(PMI)
#[1] 22.85

median(Age)
#[1] 49

median(AgonalFactor)
[1] 0

# In addition to median Brain pH, HoursFinal, and Age, an average subject was also defined as having male gender and an agonal factor of 0.

SignalCleanedofConfounds<-matrix(0, length(AmyUnfilteredSignalData[,1]), length(AmyUnfilteredSignalData[1,]))
row.names(SignalCleanedofConfounds)<-row.names(AmyUnfilteredSignalData)
colnames(SignalCleanedofConfounds)<-colnames(AmyUnfilteredSignalData)

GenderNumeric<-as.numeric(Gender)-1

for (i in 1:length(AmyUnfilteredSignalData[, 1])) {
	for (j in 1:length(AmyUnfilteredSignalData[1, ])) {
		SignalCleanedofConfounds[i, j] <- AmyUnfilteredSignalData[i, j] - ModelBetas[i, 5] * (BrainPH[j] - 6.955) - ModelBetas[i, 6] * AgonalFactor[j] - ModelBetas[i, 7] * (PMI[j] - 22.85) - ModelBetas[i, 8] * GenderNumeric[j] - ModelBetas[i, 9] * (Age[j] - 49)
	}
}

SignalCleanedofConfoundsOutput<-cbind(row.names(GeneNamesAmyunfiltered),as.character(GeneNamesAmyunfiltered$Symbol),SignalCleanedofConfounds)
colnames(SignalCleanedofConfoundsOutput)[1:2]<-c("Array_Address_Id","Gene Symbol")
write.csv(SignalCleanedofConfoundsOutput,"Amy_Unfiltered_SignalCleanedofConfounds.csv")

SignalCleanedofConfoundsDetected<-SignalCleanedofConfounds[row.names(AmyUnfilteredSignalData)%in%row.names(AmyDetectedSignalData),]
SignalCleanedofConfoundsDetectedOutput<-cbind(row.names(GeneNamesAmydetected),as.character(GeneNamesAmydetected $Symbol),SignalCleanedofConfoundsDetected)
colnames(SignalCleanedofConfoundsDetectedOutput)[1:2]<-c("Array_Address_Id","Gene Symbol")
write.csv(SignalCleanedofConfoundsDetectedOutput,"Amy_Detected_SignalCleanedofConfounds.csv")

setwd(Amy_mainWD)

Amy_fgf <- paste(Amy_model, "/FGF Family Results", sep = "")
dir.create(Amy_fgf)
setwd(Amy_fgf)


#Working below on more generaliz-able code that will analyze and plot gene-gene interactions with diagnosis for *all* the FGF-related probes that I've pulled from other lists. May go back and re-apply this approach to the Affy data to see if I can pull something out of the larger approach even if the FGF9-R1 interaction with diagnosis isn't replicating.


setwd("/Users/eaurbach/Dropbox/University of Michigan/Akil Lab/Postmortem Human Microarray/Freeze3/Freeze3 Data Reworked/FGF-associated molecules detection p-vals")
FGFAssocMoleculeList<-as.matrix(read.csv("fgf-associated molecules.csv",header=F,row.names=NULL))
colnames(FGFAssocMoleculeList)[1]="Gene Symbol"
length(FGFAssocMoleculeList)
# [1] 291
setwd(Amy_fgf)


#Now need to extract the idices of all the probes corresponding to these gene symbols, which is made more complicated by the presence of more than one probe for some of these genes (in other words, the length of the list of probes will be different from the length of the FGFAssocMoleculeList object)
#4) run gene-gene correlation stats for all the molecules against all the other molecules
#5) plot the significant ones
#6) for the list of genes for which there are multiple probes, plot probe-probe correlations for the sake of completeness to understand how much the probes are telling the same gene expression story

#first, get rid of the genes for which there are not probes in the Freeze3 datatset:
for (i in 1:length(FGFAssocMoleculeList)){
	Tmp<-which(GeneNamesAmyunfiltered$Symbol==FGFAssocMoleculeList[i])
	if (length(Tmp)==0){FGFAssocMoleculeList<-FGFAssocMoleculeList[-i]}
}
length(FGFAssocMoleculeList)
# [1] 278

#now pull all the indices for all the genes on the list (note, there will be dupliate probes for many of the genes)
FGFAssocMoleculeProbeIndices<-which(GeneNamesAmyunfiltered$Symbol==FGFAssocMoleculeList[1])
if(length(which(GeneNamesAmyunfiltered$Symbol==FGFAssocMoleculeList[1]))>1){
		FGFAssocMoleculesDuplicateProbeIndicesTF<-TRUE
	} else {
		FGFAssocMoleculesDuplicateProbeIndicesTF<-FALSE
}

for (i in 2:length(FGFAssocMoleculeList)){
	FGFAssocMoleculeProbeIndicesTmp<-which(GeneNamesAmyunfiltered$Symbol==FGFAssocMoleculeList[i])
	FGFAssocMoleculeProbeIndices<-c(FGFAssocMoleculeProbeIndices,FGFAssocMoleculeProbeIndicesTmp)
	
	if(length(which(GeneNamesAmyunfiltered$Symbol==FGFAssocMoleculeList[i]))>1){
		FGFAssocMoleculesDuplicateProbeIndicesTF<-c(FGFAssocMoleculesDuplicateProbeIndicesTF,TRUE)
	} else {
		FGFAssocMoleculesDuplicateProbeIndicesTF<-c(FGFAssocMoleculesDuplicateProbeIndicesTF,FALSE)
	}
	
}
length(FGFAssocMoleculeProbeIndices)
# [1] 394
length(FGFAssocMoleculesDuplicateProbeIndicesTF)
FGFAssocGenesWithDuplicateProbes<-FGFAssocMoleculeList[FGFAssocMoleculesDuplicateProbeIndicesTF==TRUE]
length(FGFAssocGenesWithDuplicateProbes)
# [1] 82

FGFAssocProbes<-row.names(SignalCleanedofConfounds[FGFAssocMoleculeProbeIndices,])

Labels<-paste(GeneNamesAmyunfiltered$Symbol[FGFAssocMoleculeProbeIndices], row.names(GeneNamesAmyunfiltered)[FGFAssocMoleculeProbeIndices],sep="-")

# ModelOutput #unfiltered
FGFsLinearModelUnfiltered<-ModelOutput[row.names(ModelOutput)%in%FGFAssocProbes,]
write.csv(FGFsLinearModelUnfiltered,"Amy_FGFsLinearModelUnfiltered.csv")
# ModelOutputDetected #detected
FGFsLinearModelDetected<-ModelOutputDetected[row.names(ModelOutputDetected)%in%FGFAssocProbes,]
write.csv(FGFsLinearModelDetected,"Amy_FGFsLinearModelDetected.csv")


LinearModel_GeneGeneIntDiagnosis<-function(i,j){lm(SignalCleanedofConfounds[FGFAssocMoleculeProbeIndices[i],] ~ SignalCleanedofConfounds[FGFAssocMoleculeProbeIndices[j],] + Diagnosis + Diagnosis*SignalCleanedofConfounds[FGFAssocMoleculeProbeIndices[j],])}

IntModelBetas<-matrix(0, length(FGFAssocMoleculeProbeIndices), length(FGFAssocMoleculeProbeIndices))
colnames(IntModelBetas)<-Labels
row.names(IntModelBetas)<-Labels

IntModelpvalues<-matrix(NA, length(FGFAssocMoleculeProbeIndices), length(FGFAssocMoleculeProbeIndices))
colnames(IntModelpvalues)<-Labels
row.names(IntModelpvalues)<-Labels

IntModelSE<-matrix(0, length(FGFAssocMoleculeProbeIndices), length(FGFAssocMoleculeProbeIndices))
colnames(IntModelSE)<-Labels
row.names(IntModelSE)<-Labels

IntModelTstat<-matrix(0, length(FGFAssocMoleculeProbeIndices), length(FGFAssocMoleculeProbeIndices))
colnames(IntModelTstat)<-Labels
row.names(IntModelTstat)<-Labels

#Running the linear model to assess the gene-gene relationship & interaction with diagnosis:

for (i in 1:length(FGFAssocMoleculeProbeIndices)) {
	for (j in (i):length(FGFAssocMoleculeProbeIndices)) {
		if (i == j) {
		} else {
			TempIntModel <- LinearModel_GeneGeneIntDiagnosis(i, j)
			IntModelBetas[i, j] <- summary.lm(TempIntModel)$coefficients[7, 1]
			IntModelSE[i, j] <- summary.lm(TempIntModel)$coefficients[7, 2]
			IntModelTstat[i, j] <- summary.lm(TempIntModel)$coefficients[7, 3]
			IntModelpvalues[i, j] <- summary.lm(TempIntModel)$coefficients[7, 4]

		}
	}
}
RawPList<-as.vector(IntModelpvalues)

TempIntAdjPList<-mt.rawp2adjp(RawPList, proc="BH",na.rm = TRUE)
IntAdjPList<-TempIntAdjPList$adjp[order(TempIntAdjPList$index),]
colnames(IntAdjPList)<-c("MD Int rawp","MD Int BH-Corrected")
#megan's code for permutation-based p-values (permute the various subject variables in the dataset to create a null distribution - should only apply to diagnosis) (package lmPerm). Come back to this after Copenhagen.
# Output<- lmp(formula=ValidationGeneSignalNoNA3[i,]~PsychiatricNoNA3+ Astrocyte + Microglia +  Neuron_Interneuron + Neuron_Projection + Oligodendrocyte, data=data.frame(ValidationGeneSignalNoNA3[i,], BrainPHCentered, AgonalFactorNoNA3, HoursFinalCorrectedNoNA3, AgeCentered, GenderNoNA3, PsychiatricNoNA3, DiagnosisNoNA3, SuicideNoNA3, PsychosisNoNA3, MoodNoNA3, CellTypeIndices_NoNA3_NormBestNoPrimaryOverlap), perm="Exact")

# summary(Output)


IntAdjPMatrix<-IntAdjPList[,2]
dim(IntAdjPMatrix)<-c(length(FGFAssocMoleculeProbeIndices),length(FGFAssocMoleculeProbeIndices))
colnames(IntAdjPMatrix)<-Labels
row.names(IntAdjPMatrix)<-Labels

write.csv(IntModelpvalues,"Amy_Gene-GeneInteractionWithDiagnosis_rawp.csv",row.names=TRUE,col.names=TRUE)
write.csv(IntAdjPMatrix,"Amy_Gene-GeneInteractionWithDiagnosis_BHcorrp.csv",row.names=TRUE,col.names=TRUE)

SortedIntAdjPList<-TempIntAdjPList$adjp #six entries survive multiple comparisons correction at a 0.05 threshold; 12 survive 0.1

#Extract significant data (i.e., what survives multiple comparisons) to graph
SigFGFsIntWithDiagnosis<-which(IntAdjPMatrix<0.1,arr.ind=TRUE) #note change the threshold value and rerun automated graphing code below in order to select a different subset of the gene pairs whose expression is significantly altered with diagnosis

Amy_fgfints <- paste(Amy_fgf, "/FGF Gene-Gene Interactions with Diagnosis", sep = "")
dir.create(Amy_fgfints)
setwd(Amy_fgfints)

for (k in 1:length(SigFGFsIntWithDiagnosis[, 1])) {

	TmpIntModelGraphing <- LinearModel_GeneGeneIntDiagnosis(SigFGFsIntWithDiagnosis[k, 1], SigFGFsIntWithDiagnosis[k, 2])

	pdf(paste(Labels[SigFGFsIntWithDiagnosis[k, 1]], " vs. ", Labels[SigFGFsIntWithDiagnosis[k, 2]], ".pdf", sep = ""))
	PlotTitle <- paste(Labels[SigFGFsIntWithDiagnosis[k, 1]], " and ", Labels[SigFGFsIntWithDiagnosis[k, 2]], " by Diagnosis", 
		sep = "")
	PlotXLab <- paste("Relative ", Labels[SigFGFsIntWithDiagnosis[k, 1]], "signal, cleaned of confounds")
	PlotYLab <- paste("Relative ", Labels[SigFGFsIntWithDiagnosis[k, 2]], "signal, cleaned of confounds")
	plot(SignalCleanedofConfounds[which(row.names(SignalCleanedofConfounds) == FGFAssocProbes[SigFGFsIntWithDiagnosis[k, 1]]), 
		] ~ SignalCleanedofConfounds[which(row.names(SignalCleanedofConfounds) == FGFAssocProbes[SigFGFsIntWithDiagnosis[k, 
		2]]), ], xlab = PlotXLab, ylab = PlotYLab, pch = ".", cex = 8)
	RegressionLine <- lm(SignalCleanedofConfounds[which(row.names(SignalCleanedofConfounds) == FGFAssocProbes[SigFGFsIntWithDiagnosis[k, 
		1]]), ] ~ SignalCleanedofConfounds[which(row.names(SignalCleanedofConfounds) == FGFAssocProbes[SigFGFsIntWithDiagnosis[k, 
		2]]), ])
	RegressionLineMD <- lm(SignalCleanedofConfounds[which(row.names(SignalCleanedofConfounds) == FGFAssocProbes[SigFGFsIntWithDiagnosis[k, 
		1]]), Diagnosis == "MD"] ~ SignalCleanedofConfounds[which(row.names(SignalCleanedofConfounds) == FGFAssocProbes[SigFGFsIntWithDiagnosis[k, 
		2]]), Diagnosis == "MD"])
	RegressionLineBP <- lm(SignalCleanedofConfounds[which(row.names(SignalCleanedofConfounds) == FGFAssocProbes[SigFGFsIntWithDiagnosis[k, 
		1]]), Diagnosis == "BP"] ~ SignalCleanedofConfounds[which(row.names(SignalCleanedofConfounds) == FGFAssocProbes[SigFGFsIntWithDiagnosis[k, 
		2]]), Diagnosis == "BP"])
	RegressionLineSchiz <- lm(SignalCleanedofConfounds[which(row.names(SignalCleanedofConfounds) == FGFAssocProbes[SigFGFsIntWithDiagnosis[k, 
		1]]), Diagnosis == "SC"] ~ SignalCleanedofConfounds[which(row.names(SignalCleanedofConfounds) == FGFAssocProbes[SigFGFsIntWithDiagnosis[k, 
		2]]), Diagnosis == "SC"])
	RegressionLineCtrl <- lm(SignalCleanedofConfounds[which(row.names(SignalCleanedofConfounds) == FGFAssocProbes[SigFGFsIntWithDiagnosis[k, 
		1]]), Diagnosis == "C"] ~ SignalCleanedofConfounds[which(row.names(SignalCleanedofConfounds) == FGFAssocProbes[SigFGFsIntWithDiagnosis[k, 
		2]]), Diagnosis == "C"])
	abline(RegressionLine, col = 2)
	abline(RegressionLineMD, col = "magenta")
	abline(RegressionLineBP, col = "blue")
	abline(RegressionLineSchiz, col = "green")
	abline(RegressionLineCtrl, col = "black")
	points(SignalCleanedofConfounds[which(row.names(SignalCleanedofConfounds) == FGFAssocProbes[SigFGFsIntWithDiagnosis[k, 1]]), 
		Diagnosis == "MD"] ~ SignalCleanedofConfounds[which(row.names(SignalCleanedofConfounds) == FGFAssocProbes[SigFGFsIntWithDiagnosis[k, 
		2]]), Diagnosis == "MD"], col = "magenta", pch = ".", cex = 8)
	points(SignalCleanedofConfounds[which(row.names(SignalCleanedofConfounds) == FGFAssocProbes[SigFGFsIntWithDiagnosis[k, 1]]), 
		Diagnosis == "BP"] ~ SignalCleanedofConfounds[which(row.names(SignalCleanedofConfounds) == FGFAssocProbes[SigFGFsIntWithDiagnosis[k, 
		2]]), Diagnosis == "BP"], col = "blue", pch = ".", cex = 8)
	points(SignalCleanedofConfounds[which(row.names(SignalCleanedofConfounds) == FGFAssocProbes[SigFGFsIntWithDiagnosis[k, 1]]), 
		Diagnosis == "SC"] ~ SignalCleanedofConfounds[which(row.names(SignalCleanedofConfounds) == FGFAssocProbes[SigFGFsIntWithDiagnosis[k, 
		2]]), Diagnosis == "SC"], col = "green", pch = ".", cex = 8)
	mtext(paste("All Data: p-value=", round(summary.lm(RegressionLine)$coefficients[8], digits = 6), ", r-squared=", round(summary.lm(RegressionLine)$r.squared, 
		digits = 3), sep = " "), line = 0)
	mtext(paste("MDD Data: p-value=", round(summary.lm(TmpIntModelGraphing)$coefficients[7, 4], digits = 6)), line = 2)
	mtext(paste("BP Data: p-value=", round(summary.lm(TmpIntModelGraphing)$coefficients[6, 4], digits = 6)), line = 3)
	mtext(paste("Schiz Data: p-value=", round(summary.lm(TmpIntModelGraphing)$coefficients[8, 4], digits = 6)), line = 1)
	dev.off()

}
#since the results are being driven highly by a couple of outliers, Megan suggests trying robust regression (function rlm), which will minimize the impact of outliers on the regression line


FGF2vsFGFR1<-summary.lm(lm(SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol=="FGF2",]~SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol[17237],]+Diagnosis+SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol[17237],]*Diagnosis))
pdf("FGF2vsFGFR1Cleaned_general_nolegend.pdf")
plot(SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol=="FGF2",]~SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol[17237],], xlab="Relative FGFR1 Signal: Cleaned of Confounds", ylab="Relative FGF2 Signal: Cleaned of Confounds",pch=".",cex=8)
RegressionLine<-lm(SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol=="FGF2",]~SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol[17237],])
RegressionLineMD<-lm(SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol=="FGF2", Diagnosis=="MD"]~SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol[17237], Diagnosis=="MD"])
RegressionLineBP<-lm(SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol=="FGF2", Diagnosis=="BP"]~SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol[17237], Diagnosis=="BP"])
RegressionLineSchiz<-lm(SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol=="FGF2", Diagnosis=="SC"]~SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol[17237], Diagnosis=="SC"])
RegressionLineCtrl<-lm(SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol=="FGF2", Diagnosis=="C"]~SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol[17237], Diagnosis=="C"])
abline(RegressionLine, col=2)
abline(RegressionLineMD,col="magenta")
abline(RegressionLineBP,col="blue")
abline(RegressionLineSchiz,col="green")
abline(RegressionLineCtrl,col="black")
points(SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol=="FGF2", Diagnosis=="MD"]~SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol[17237], Diagnosis=="MD"], col="magenta",pch=".",cex=8)
points(SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol=="FGF2", Diagnosis=="BP"]~SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol[17237], Diagnosis=="BP"], col="blue",pch=".",cex=8)
points(SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol=="FGF2", Diagnosis=="SC"]~SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol[17237], Diagnosis=="SC"], col="green",pch=".",cex=8)
mtext(paste("All Data: p-value=", round(summary.lm(RegressionLine)$coefficients[8], digits=6), ", r-squared=", round(summary.lm(RegressionLine)$r.squared, digits=3), sep=" "),line=0)
mtext(paste("MDD Data: p-value=", round(FGF2vsFGFR1$coefficients[7,4], digits=6)),line=2)
mtext(paste("BP Data: p-value=", round(FGF2vsFGFR1$coefficients[6,4], digits=6)),line=3)
mtext(paste("Schiz Data: p-value=", round(FGF2vsFGFR1$coefficients[8,4], digits=6)),line=1)
# legend(min(SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol[17237],])+0.1, max(SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol=="FGF2",])-0.025, c("Control", "MDD", "BP", "Schiz"), text.col=c("black","magenta","blue","green"), pch=19, col=c("black","magenta","blue","green"))
dev.off()


#****initial glimpses at cleaned data for FGF2 and FGF9


FGFGeneGeneInteractionsWithDiagnosis<-file("FGF Gene-Gene Interactions with Diagnosis.txt")
out<-c(

capture.output(
for (i in 1:length(FGFAssocMoleculeProbeIndices)){
for(j in (i):length(FGFAssocMoleculeProbeIndices)){
	if(i==j){}
	else {
# print(paste(GeneNamesAmyunfiltered$Symbol[FGFAssocMoleculeProbeIndices[i]],"-", GeneNamesAmyunfiltered$Symbol[FGFAssocMoleculeProbeIndices[j]],sep=""))

# print(summary.lm(lm(SignalCleanedofConfounds[FGFAssocMoleculeProbeIndices[i],] ~ SignalCleanedofConfounds[FGFAssocMoleculeProbeIndices[j],] + Diagnosis + Diagnosis*SignalCleanedofConfounds[FGFAssocMoleculeProbeIndices[j],])))

if(summary.lm(lm(SignalCleanedofConfounds[FGFAssocMoleculeProbeIndices[i],] ~ SignalCleanedofConfounds[FGFAssocMoleculeProbeIndices[j],] + Diagnosis + Diagnosis*SignalCleanedofConfounds[FGFAssocMoleculeProbeIndices[j],]))$coefficient[7,4] <= 0.05) {
	print(paste(GeneNamesAmyunfiltered$Symbol[FGFAssocMoleculeProbeIndices[i]],"-", GeneNamesAmyunfiltered$Symbol[FGFAssocMoleculeProbeIndices[j]],":  MDD interaction pval =", summary.lm(lm(SignalCleanedofConfounds[FGFAssocMoleculeProbeIndices[i],] ~ SignalCleanedofConfounds[FGFAssocMoleculeProbeIndices[j],] + Diagnosis + Diagnosis*SignalCleanedofConfounds[FGFAssocMoleculeProbeIndices[j],]))$coefficient[7,4],sep=" "))
}
}
}		
}
)
)
cat(out, file="FGF Gene-Gene Interactions with Diagnosis.txt", sep="\n", append=TRUE)
close(FGFGeneGeneInteractionsWithDiagnosis)
rm(out)


#*****Graphing******

png("FGF9vsDiagnosisCleaned.png")
boxplot(SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol=="FGF9",]~Diagnosis, ylab="Normalized Signal, Cleaned of Confounds")
mtext(paste("BP p-value=", round(Modelpvalues[GeneNamesAmyunfiltered$Symbol=="FGF9",2], digits=6), ",",
"MDD p-value=", round(Modelpvalues[GeneNamesAmyunfiltered$Symbol=="FGF9",3], digits=6), ",",
"Schiz p-value=", round(Modelpvalues[GeneNamesAmyunfiltered$Symbol=="FGF9",4], digits=6)
)) 
dev.off()

png("FGF2vsDiagnosisCleaned.png")
boxplot(SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol=="FGF2",]~Diagnosis, ylab="Normalized Signal, Cleaned of Confounds")
mtext(paste("BP p-value=", round(Modelpvalues[GeneNamesAmyunfiltered$Symbol=="FGF2",2], digits=6), ",",
"MDD p-value=", round(Modelpvalues[GeneNamesAmyunfiltered$Symbol=="FGF2",3], digits=6), ",",
"Schiz p-value=", round(Modelpvalues[GeneNamesAmyunfiltered$Symbol=="FGF2",4], digits=6)
)) 
dev.off()

GeneGeneCorsByDiagnosis<-file("Gene-Gene Relationships by Diagnosis.txt")
out<-c(
print("FGF2-FGF9"),
capture.output(
summary.lm(lm(SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol=="FGF9",]~SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol=="FGF2",]+Diagnosis+SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol=="FGF2",]*Diagnosis))
),print("FGF9-FGFR1"),
capture.output(
summary.lm(lm(SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol=="FGF9",]~SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol[17237],]+Diagnosis+SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol[17237],]*Diagnosis))
),print("FGF9-FGFR2"),
capture.output(
summary.lm(lm(SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol=="FGF9",]~SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol=="FGFR2",]+Diagnosis+SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol=="FGFR2",]*Diagnosis))
),print("FGF9-FGFR3"),
capture.output(
summary.lm(lm(SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol=="FGF9",]~SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol[4360],]+Diagnosis+SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol[4360],]*Diagnosis))
),print("FGF2-FGFR1"),
capture.output(
summary.lm(lm(SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol=="FGF2",]~SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol[17237],]+Diagnosis+SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol[17237],]*Diagnosis))
),print("FGF2-FGFR2"),
capture.output(
summary.lm(lm(SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol=="FGF2",]~SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol=="FGFR2",]+Diagnosis+SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol=="FGFR2",]*Diagnosis))
),print("FGF2-FGFR3"),
capture.output(
summary.lm(lm(SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol=="FGF2",]~SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol[4360],]+Diagnosis+SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol[4360],]*Diagnosis))
)
)
cat(out, file="Gene-Gene Relationships by Diagnosis.txt", sep="\n", append=TRUE)
close(GeneGeneCorsByDiagnosis)
rm(out)
#Note: some of these got replaced with index values because FGFR1 and FGFR3 both had multiple probes in the Freeze3 dataset


FGF2vsFGFR1<-summary.lm(lm(SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol=="FGF2",]~SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol[17237],]+Diagnosis+SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol[17237],]*Diagnosis))
pdf("FGF2vsFGFR1Cleaned_general_nolegend.pdf")
plot(SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol=="FGF2",]~SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol[17237],], xlab="Relative FGFR1 Signal: Cleaned of Confounds", ylab="Relative FGF2 Signal: Cleaned of Confounds",pch=".",cex=8)
RegressionLine<-lm(SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol=="FGF2",]~SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol[17237],])
RegressionLineMD<-lm(SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol=="FGF2", Diagnosis=="MD"]~SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol[17237], Diagnosis=="MD"])
RegressionLineBP<-lm(SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol=="FGF2", Diagnosis=="BP"]~SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol[17237], Diagnosis=="BP"])
RegressionLineSchiz<-lm(SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol=="FGF2", Diagnosis=="SC"]~SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol[17237], Diagnosis=="SC"])
RegressionLineCtrl<-lm(SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol=="FGF2", Diagnosis=="C"]~SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol[17237], Diagnosis=="C"])
abline(RegressionLine, col=2)
abline(RegressionLineMD,col="magenta")
abline(RegressionLineBP,col="blue")
abline(RegressionLineSchiz,col="green")
abline(RegressionLineCtrl,col="black")
points(SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol=="FGF2", Diagnosis=="MD"]~SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol[17237], Diagnosis=="MD"], col="magenta",pch=".",cex=8)
points(SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol=="FGF2", Diagnosis=="BP"]~SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol[17237], Diagnosis=="BP"], col="blue",pch=".",cex=8)
points(SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol=="FGF2", Diagnosis=="SC"]~SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol[17237], Diagnosis=="SC"], col="green",pch=".",cex=8)
mtext(paste("All Data: p-value=", round(summary.lm(RegressionLine)$coefficients[8], digits=6), ", r-squared=", round(summary.lm(RegressionLine)$r.squared, digits=3), sep=" "),line=0)
mtext(paste("MDD Data: p-value=", round(FGF2vsFGFR1$coefficients[7,4], digits=6)),line=2)
mtext(paste("BP Data: p-value=", round(FGF2vsFGFR1$coefficients[6,4], digits=6)),line=3)
mtext(paste("Schiz Data: p-value=", round(FGF2vsFGFR1$coefficients[8,4], digits=6)),line=1)
# legend(min(SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol[17237],])+0.1, max(SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol=="FGF2",])-0.025, c("Control", "MDD", "BP", "Schiz"), text.col=c("black","magenta","blue","green"), pch=19, col=c("black","magenta","blue","green"))
dev.off()


FGF9vsFGFR1<-summary.lm(lm(SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol=="FGF9",]~SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol[17237],]+Diagnosis+SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol[17237],]*Diagnosis))
pdf("FGF9vsFGFR1Cleaned_general_nolegend.pdf")
plot(SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol=="FGF9",]~SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol[17237],], xlab="Relative FGFR1 Signal: Cleaned of Confounds", ylab="Relative FGF9 Signal: Cleaned of Confounds",pch=".",cex=8)
RegressionLine<-lm(SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol=="FGF9",]~SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol[17237],])
RegressionLineMD<-lm(SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol=="FGF9", Diagnosis=="MD"]~SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol[17237], Diagnosis=="MD"])
RegressionLineBP<-lm(SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol=="FGF9", Diagnosis=="BP"]~SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol[17237], Diagnosis=="BP"])
RegressionLineSchiz<-lm(SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol=="FGF9", Diagnosis=="SC"]~SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol[17237], Diagnosis=="SC"])
RegressionLineCtrl<-lm(SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol=="FGF9", Diagnosis=="C"]~SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol[17237], Diagnosis=="C"])
abline(RegressionLine, col=2)
abline(RegressionLineMD,col="magenta")
abline(RegressionLineBP,col="blue")
abline(RegressionLineSchiz,col="green")
abline(RegressionLineCtrl,col="black")
points(SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol=="FGF9", Diagnosis=="MD"]~SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol[17237], Diagnosis=="MD"], col="magenta",pch=".",cex=8)
points(SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol=="FGF9", Diagnosis=="BP"]~SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol[17237], Diagnosis=="BP"], col="blue",pch=".",cex=8)
points(SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol=="FGF9", Diagnosis=="SC"]~SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol[17237], Diagnosis=="SC"], col="green",pch=".",cex=8)
mtext(paste("All Data: p-value=", round(summary.lm(RegressionLine)$coefficients[8], digits=6), ", r-squared=", round(summary.lm(RegressionLine)$r.squared, digits=3), sep=" "),line=0)
mtext(paste("MDD Data: p-value=", round(FGF9vsFGFR1$coefficients[7,4], digits=6)),line=2)
mtext(paste("BP Data: p-value=", round(FGF9vsFGFR1$coefficients[6,4], digits=6)),line=3)
mtext(paste("Schiz Data: p-value=", round(FGF9vsFGFR1$coefficients[8,4], digits=6)),line=1)
# legend(min(SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol[17237],])+0.1, max(SignalCleanedofConfounds[GeneNamesAmyunfiltered$Symbol=="FGF9",])-0.025, c("Control", "MDD", "BP", "Schiz"), text.col=c("black","magenta","blue","green"), pch=19, col=c("black","magenta","blue","green"))
dev.off()





























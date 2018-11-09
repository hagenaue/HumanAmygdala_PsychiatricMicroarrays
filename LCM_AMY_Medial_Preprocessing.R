##Generalized code for processing Vikram's amygdala Illumina microarray data:
##By Megan Hagenauer
##First attempt: Lateral Amygdala (6/24/2013, 6/25/2013, 6/27/2013)
##Generalized: 6/28/2013... 7/12/2013.


##R version 3.0.1 (2013-05-16) -- "Good Sport"


#******************************************************************

#STEP 1: Set the current working directory and workspace file:

##1. Click on the R Console 
##2. Go to File->change dir
##3. Z:\Beadchip files\Data\real\Genexp\  [...brain region folder...]



#*******************************************************************

#STEP 2: Install and load the necessary packages of R code:

##1. Go to Packages-> set CRAN mirror
##2. Choose USA(MI)
##3. Go to Packages->Install Packages
##4. Click on the packages of interest: gdata, fields, stats, car

##4. Go to Packages-> Set repositories
##5. Set the repositories to include "BioC software"
##6. Go to Packages->Install Packages
##7. Click on the packages of interest: preprocessCore, affy, multtest

##8. Load the packages by highlighting this code and hitting run (Ctrl+R):

#***CODE TO RUN FOR STEP 2****

library(gdata)
library(fields)
library(stats)
library(car)
library(affy)
library(preprocessCore)
library(multtest)
library(plyr)


#***END OF CODE FOR STEP 2****


#*******************************************************************

#STEP 3: Make sure the Illumina data is ready to be inputted into R:

##1. The data file should be outputted from GenomeStudio using "Export displayed data to a file."
##2. The data is a matrix of the raw gene expression data for each probe (by row) for each sample (by column).
##3. Save the file as "genexp_sample probe profile.csv" and place it in the appropriate folder for your brain region.
##4. Using Excel, open up the file and double check that it looks correct before inputting the data into R. The file should include:
	##A) A "header" of column names which identify the sample by the chip # and location on the chip.
	##B) The probe names + gene symbols should fill columns 1 and 2.
	##C) The final columns (farthest right) should be a slew of gene identification information (full name, chromosome, etc)
	##D) All columns in between should alternate between the Average Expression and Detection P-values for each sample.
	##E) There should be no additional information above the matrix of data (e.g. Illumina processing information).   
##5. If the data is outputted from GenomeStudio using "Export to Gene Spring Gx" it will have incorrect detection pvalues.



#*******************************************************************

#STEP 4: Read the Illumina data into R and convert the data into a workable format:

##1. Read in the data matrix by highlighting this code and hitting run (Ctrl+R) ***NOTE THAT THIS MAY TAKE A FEW MINUTES***


#***CODE TO RUN FOR STEP 4****

#Reads in the raw data:
RawData<-as.matrix(read.csv("Sample Probe Profile_Medial_genexp_rescan.csv", header=T, row.names=1))

#Determines how many probes there are for this raw data set:
TotalNumberOfProbes<-length(RawData[,1])

#Skipped - there are 142 columns of data now because there is not only AVGSignal and Detection P-value, but also BEADStError, and everything after 128 is gene identification information.
#Determines if there are the expected number of columns in the raw data set:
if(length(RawData[1,])==!90){print("ERROR! The data matrix is the wrong size! Fix now or everything else will be wrong!")}else{print("Input matrix is the correct size")}

#Extracts out the average signal values for each subject and places them into a numeric matrix with appropriate row and column names:
AVG_Signal<-matrix(as.numeric(RawData[, seq(from=2, to=127, by=3)]), nrow=TotalNumberOfProbes, ncol=42)
colnames(AVG_Signal)<-colnames(RawData[, seq(from=2, to=127, by=3)])
row.names(AVG_Signal)<-row.names(RawData)

#Double checks whether the average signal matrix is the correct size and numeric:
if(length(AVG_Signal[1,])==!42){print("ERROR! The AVG_Signal matrix is the wrong size! Fix now or everything else will be wrong!")}else{print("AVG_Signal matrix is the correct size")}
if(is.numeric(AVG_Signal)==F){print("ERROR! The AVG_Signal matrix is not numeric! Fix now or everything else will be wrong!")}else{print("AVG_Signal is properly numeric. Hurray!")}

#Extracts out the detection p-values for each subject and places them into a numeric matrix with appropriate row and column names:
Detection_Pvalue<-matrix(as.numeric(RawData[,seq(from=3, to=127, by=3)]), nrow=TotalNumberOfProbes, ncol=42)
colnames(Detection_Pvalue)<-colnames(RawData[, seq(from=3, to=127, by=3)])
row.names(Detection_Pvalue)<-row.names(RawData)

#Double checks whether the detection p-value matrix is the correct size and numeric:
if(length(Detection_Pvalue[1,])==!42){print("ERROR! The Detection_Pvalue matrix is the wrong size! Fix now or everything else will be wrong!")}else{print("Detection_Pvalue matrix is the correct size")}
if(is.numeric(Detection_Pvalue)==F){print("ERROR! The Detection_Pvalue matrix is not numeric! Fix now or everything else will be wrong!")}else{print("Detection_Pvalue is properly numeric. Hurray!")}

#Makes a matrix where the subject identifiers for the average signal and detection p-value matrices can be easily compared to make sure they are the same, and then outputs that data to a format that can be accessed by Excel:
DoubleCheckColumnNames<-cbind(colnames(AVG_Signal), colnames(Detection_Pvalue))
write.csv(DoubleCheckColumnNames, "04 DoubleCheckColumnNames.csv")

#Outputs the gene information for each of the probes (symbol, chromosome, full name, etc.)
GeneIdentifiers<-RawData[,c(1, c(128:142))]

#Outputs a text-file summarizing the process of reading in the Illumina data and any error messages that may have appeared during the process:
ReadingInIllumina<-file("04 Reading in Illumina Raw Data.txt")
out<-c(
print("
The total number of probes in the raw gene expression data set for this brain region:"),
capture.output(TotalNumberOfProbes),

print("
The dimensions of the raw data input matrix:"),
capture.output(dim(RawData)),
capture.output(if(length(RawData[1,])==!90){print("ERROR! The data matrix is the wrong size! Fix now or everything else will be wrong!")}else{print("Input matrix is the correct size")}),
capture.output(RawData[c(1:3), c(1:3)]),

print("
The dimensions of the raw AVG_Signal matrix:"),
capture.output(dim(AVG_Signal)),
capture.output(if(length(AVG_Signal[1,])==!42){print("ERROR! The AVG_Signal matrix is the wrong size! Fix now or everything else will be wrong!")}else{print("AVG_Signal matrix is the correct size")}),
capture.output(if(is.numeric(AVG_Signal)==F){print("ERROR! The AVG_Signal matrix is not numeric! Fix now or everything else will be wrong!")}else{print("AVG_Signal is properly numeric. Hurray!")}),
capture.output(AVG_Signal[c(1:3), c(1:2)]),

print("
The dimensions of the raw Detection_Pvalue matrix:"),
capture.output(dim(Detection_Pvalue)),
capture.output(if(length(Detection_Pvalue[1,])==!42){print("ERROR! The Detection_Pvalue matrix is the wrong size! Fix now or everything else will be wrong!")}else{print("Detection_Pvalue matrix is the correct size")}),
capture.output(if(is.numeric(Detection_Pvalue)==F){print("ERROR! The Detection_Pvalue matrix is not numeric! Fix now or everything else will be wrong!")}else{print("Detection_Pvalue is properly numeric. Hurray!")}),
capture.output(Detection_Pvalue[c(1:3), c(1:2)]),

print("
The column names for AVG_Signal and Detection_Pvalue... double check to make sure they are correct!"),
capture.output(DoubleCheckColumnNames),

print("
Gene Identifiers:"),
capture.output(GeneIdentifiers[c(1:3),])
)
cat(out, file="04 Reading in Illumina Raw Data.txt", sep="\n", append=TRUE)
close(ReadingInIllumina)
rm(out)


#***END OF CODE FOR STEP 4****

##2.  Take a look at "Reading in Illumina Raw Data.txt" and make sure everything looks o.k.  What to look for:
###A. Are the total number of probes what you would expect? (between 44,000-48,000)?
###B. Are there any error messages?
###C. Do the initial pieces of information in each matrix look correct?  (e.g., Do signal values look like signal values? Do p-values look like p-values? Do the gene identfiers look like gene information for the probes?)
###D. Do the column names for the AVG_Signal and Detection_Pvalue show the same subject identifiers? (chip+location on chip?)



#*******************************************************************

#STEP 5: Make sure the Subject Identification information is ready to be inputted into R:

##1. We need to convert the column names for AVG_Signal and Detection_Pvalue into Subject ID numbers.
##2. The columns are currently named by the chip# and the sample's location on the chip (outputted into file "DoubleCheckColumnNames.csv").
##3. Open the file "DoubleCheckColumnNames.csv" using Excel.
##4. Double check that the chip# and sample locations are the same for AVG_Signal and Detection_Pvalue.
##5. Copy and paste the column names for AVG_Signal into column 1 of a new excel file.
##6. Label the headers for columns 2 and 3: "Sample.Group" "Subject.Number".
##7. Decode (i.e. fill in the sample.group and subject.number for each AVG_Signal ID...but do not change the order!).
##8. Save the new excel file as "Sample ID file.csv" and place it in the appropriate folder for your brain region.


#*******************************************************************

#STEP 6: Read the Subject Identification file into R and convert the data into a workable format:

##1.  Read in the Subject Identification file by highlighting this code and hitting run (Ctrl+R):


#***CODE TO RUN FOR STEP 6****

#Reads in the subject identification information to decode the current chip and chip location information associated with each average signal and detection p-value column:
SubjectID<-as.matrix(read.csv("Sample ID file.csv", header=T))

#Double checks that the subject identification matrix is the correct size:
if(length(SubjectID[1,])==!3){print("ERROR! The data matrix is the wrong size! Fix now or everything else will be wrong!")}else{print("Subject ID matrix has the correct # of columns")}
if(length(SubjectID[,1])==!42){print("ERROR! The data matrix is the wrong size! Fix now or everything else will be wrong!")}else{print("Subject ID matrix has the correct # of rows")}

#Replaces the column names for the average signal and detection p-value matrices with the actual subject ID numbers:
colnames(AVG_Signal)<-as.numeric(SubjectID[,3])
colnames(Detection_Pvalue)<-as.numeric(SubjectID[,3])

#Makes diagnosis a factor with control as its reference level so that we can easily use it in later analyses:
Diagnosis<-as.factor(SubjectID[,2])
Diagnosis<-relevel(Diagnosis, ref="CTRL")

#Outputs a text file that contains a summary of reading in the subject ID information:
ReadingInSubjectID<-file("06 Reading in Subject ID Data.txt")
out<-c(
print("
The size of the Subject ID matrix:"),
capture.output(dim(SubjectID)),
capture.output(if(length(SubjectID[1,])==!3){print("ERROR! The data matrix is the wrong size! Fix now or everything else will be wrong!")}else{print("Subject ID matrix has the correct # of columns")}),
capture.output(if(length(SubjectID[,1])==!42){print("ERROR! The data matrix is the wrong size! Fix now or everything else will be wrong!")}else{print("Subject ID matrix has the correct # of rows")}),
capture.output(SubjectID[c(1:3), c(1:3)]),

print("
The column names for the data matrices have been changed to match the subject ID's:"),
capture.output(AVG_Signal[c(1:3), c(1:3)]),
capture.output(Detection_Pvalue[c(1:3), c(1:3)]),

print("
Diagnosis is now a factor with 2 levels:"),
capture.output(str(Diagnosis))
)
cat(out, file="06 Reading in Subject ID Data.txt", sep="\n", append=TRUE)
close(ReadingInSubjectID)
rm(out)

#***END OF CODE FOR STEP 6****

##2.  Take a look at "06 Reading in Subject ID Data.txt" and make sure everything looks o.k.  What to look for:
###A. Are there any error messages?
###B. Do the initial pieces of information in each matrix look correct?  (e.g., Does the subject ID matrix contain the AVG signal column names, "Sample.Group" and "Subject.Number", Does Diagnosis now have "Control" as its reference level?)




#*******************************************************************

#STEP 7: Log transforming the gene expression data.  
#The reason for log-transforming the data is because the variability present in the measurements with a large amount of signal tends to be greater than the variability present in the measurements with a low amount of signal (heteroskedasticity).  This heteroskedasticity can introduce error into some statistical tests. Log transforming the data reduces the relative variability in the measurements with a large amount of signal, and thus reduces the problem of heteroskedasticity (although it should be noted that our microarray data is still relatively heteroskedastic even after the log transform!).


##1. Highlight this code and hit run (Ctrl+R). ***Note that this code may take a few minutes!***


#***CODE TO RUN FOR STEP 7****

#Outputting the raw detection pvalue and average signal matrices before transforming them:
write.csv(Detection_Pvalue, "Detection_Pvalue.csv")
write.csv(AVG_Signal, "AVG_Signal.csv")

#Determining whether the minimum signal is below 0 and therefore can't be log-transformed:
MinimumSignals<-apply(AVG_Signal, 2, min)

png("07 MinimumSignalsLateral.png")
plot(MinimumSignals, main="Minimum Signal per Sample", xlab="Sample Index#", ylab="Minimum Signal")
dev.off()

if(min(MinimumSignals)<0|min(MinimumSignals)==0){
apply(AVG_Signal, 2, function(y) y+MinimumSignals)}else{print("Data can be log transformed: All minimum signals are above 0")}

#Log transforming the average signal (using log base 2)
Log.AVG_Signal<-log(AVG_Signal, 2)
colnames(Log.AVG_Signal)<-colnames(AVG_Signal)
row.names(Log.AVG_Signal)<-row.names(AVG_Signal)

#Double checking that the log transformed data matrix is the correct size:
if(length(Log.AVG_Signal[,1])==!TotalNumberOfProbes){print("ERROR: Log transformed data matrix is the wrong size.")}else{print("Log transformed data matrix has the correct # of rows.")}
if(length(Log.AVG_Signal[1,])==!42){print("ERROR: Log transformed data matrix is the wrong size.")}else{print("Log transformed data matrix has the correct # of columns.")}

#Double checking that one of the pieces of log-transformed data matches what we would expect:
if(2^AVG_Signal[1,2]==!Log.AVG_Signal[1,2]){print("ERROR: Something is wrong with the log transformation!")}else{print("The log transformation completed properly.")}

#Outputting a text file that summarizes the process of log-transforming the data:
LogTransformation<-file("07 LogTransformation.txt")
out<-c(
print("Log Transforming the Data: Matrix Dimensions:"),
capture.output(dim(Log.AVG_Signal)),
capture.output(if(length(Log.AVG_Signal[,1])==!TotalNumberOfProbes){print("ERROR: Log transformed data matrix is the wrong size.")}else{print("Log transformed data matrix has the correct # of rows.")}),
capture.output(if(length(Log.AVG_Signal[1,])==!42){print("ERROR: Log transformed data matrix is the wrong size.")}else{print("Log transformed data matrix has the correct # of columns.")}),
capture.output(Log.AVG_Signal[c(1:3), c(1:3)]),
capture.output(AVG_Signal[c(1:3), c(1:3)]),
capture.output(if(2^AVG_Signal[1,2]==!Log.AVG_Signal[1,2]){print("ERROR: Something is wrong with the log transformation!")}else{print("The log transformation completed properly.")})
)
cat(out, file="07 LogTransformation.txt", sep="\n", append=TRUE)
close(LogTransformation)
rm(out)

#Outputting a boxplot of the log-transformed values:
png("07 Boxplot logAVGSignal Lateral.png")
boxplot(data.frame(Log.AVG_Signal), cex=0.25, las=3, par(cex.axis=0.75), main="Boxplot of log signal values per sample (1 box=all 47k probes)", xlab="Sample ID", ylab="Log Signal")
dev.off()

#***END OF CODE FOR STEP 7****


##2.  Take a look at "07 MinimumSignalsLateral.png", "07 LogTransformation.txt", "07 Boxplot logAVGSignal Lateral.png" and make sure everything looks o.k.
### What to look for:
###A. Do any of the samples have a minimum signal that is close to 0 (<2) or that dramatically differs from the other samples' minimum signal?
###B. Are there any error messages in LogTransformation.txt?
###C. In LogTransformation.txt, do the initial values in the log transformed signal matrix look like they are the log transform (base 2) of the average signal matrix?
###D. Does the Boxplot showing the log-transformed average signal values for each sample indicate that any of the samples have overall a dramatically different amount of signal? (e.g. Do the interquartile ranges for all of the boxplots overlap?)  If not, you do not necessarily need to throw the questionable sample out, but you should mentally "flag" it as suspicious.

#*******************************************************************


#STEP 8: Filtering probes by whether they surpass a threshold for significant detection (i.e. accurate measurement). 
#The reason for doing this is two-fold: 1) those probes that do not surpass a threshold for significant detection for the vast majority of samples are less likely to be truely expressed in the brain region, and therefore of less theoretical interest, 2) probes that have signal values that are near the threshold for detection are likely to have restricted variability (a floor effect), meaning that the few samples that *do* have levels that surpass the threshold of detection for that probe are likely to have a disproportionate influence on statistical tests. This means that when running statistical comparisons based on diagnosis we are likely to have a larger number of genes that are falsely associated with diagnosis (Type 1 Error). 

##1. Highlight this code and hit run (Ctrl+R):


#***CODE TO RUN FOR STEP 8****

#Determines whether each probe measurement surpass the threshold of detection (with an alpha of 0.05) for each sample:
detect05.gene<-apply(Detection_Pvalue, 1, function(x) sum(x<0.05))

#Outputs a plot indicating the number of samples for which a probe surpasses the threshold of detection:
png("08 Detection pvalues Lateral.png")
plot(sort(detect05.gene), main="Number of probes surpassing detection threshold in 20% of samples", xlab="Probe Index #, sorted", ylab="Number of Samples with Significant Detection")
abline(a=8, b=0)
dev.off()

#Filters out all probes that do not surpass a threshold of detection in at least ~20% of the samples (8 subjects):
log.avg_signal_filtered<-Log.AVG_Signal[detect05.gene>8, ]
GeneIdentifiers_filtered<-GeneIdentifiers[detect05.gene>8, ]

#Calculates the number of probes left after filtering out the probes that weren't reaching the threshold of detection in 20% of the samples:
NumberofFilteredProbes<-length(log.avg_signal_filtered[,1])

#Double checks that the matrix of filtered probes still has the correct number of subjects:
if(length(log.avg_signal_filtered[1,])==!42){print("ERROR: Filtered log transformed data matrix is the wrong size.")}else{print("Filtered log transformed data matrix has the correct # of columns.")}

#Outputs a text file summarizing the process of filtering probes:
FilteringBySigDetection<-file("08 Filtering By Sig Detection.txt")
out<-c(
print("Probes were filtered by whether at least 20% of the samples showed significant detection (p-value<0.05). The number of probes meeting this threshold was:"),
print(NumberofFilteredProbes),
capture.output(if(length(log.avg_signal_filtered[1,])==!42){print("ERROR: Filtered log transformed data matrix is the wrong size.")}else{print("Filtered log transformed data matrix has the correct # of columns.")}),
capture.output(dim(log.avg_signal_filtered)),
capture.output(log.avg_signal_filtered[c(1:3), c(1:3)])
)
cat(out, file="08 Filtering By Sig Detection.txt", sep="\n", append=TRUE)
close(FilteringBySigDetection)
rm(out)

#Outputs a boxplot illustrating the distribution of log-transformed signal for each subject for only the filtered probes:
png("08 Boxplot logAVGSignal Lateral Filtered.png")
boxplot(data.frame(log.avg_signal_filtered), cex=0.25, las=3, par(cex.axis=0.75), main="Boxplot of log signal values per sample (1 box=all filtered probes)", xlab="Sample ID", ylab="Log Signal")
dev.off()

#***END OF CODE FOR STEP 8****


##2. Take a look at "08 Detection pvalues Lateral.png", "08 Filtering By Sig Detection.txt", and "08 Boxplot logAVGSignal Lateral Filtered.png" and make sure everything looks o.k. What to look for:
###A. When looking at the "08 Detection pvalues Lateral.png" plot, does it seem like there is a relatively steep slope? (i.e. are most of the probes either not detected in most samples or detected in most samples)
###B. Are there any error messages in "08 Filtering by Sig Detection.txt"?
###C. When looking at "08 Filtering by Sig Detection.txt", how many probes are left after filtering out the probes that weren't detected in at least 20% of the samples? Does this seem similar to your other brain regions?  
###D. Does the Boxplot showing the (now filtered) log-transformed average signal values for each sample indicate that any of the samples have overall a dramatically different amount of signal? (e.g. Do the interquartile ranges for all of the boxplots overlap?)  If not, you do not necessarily need to throw the questionable sample out, but you should mentally "flag" it as suspicious.


#*******************************************************************


#STEP 9: Quantile Normalization and visualizing sample-sample correlations. 
# Quantile normalization is a rank-based method of normalizing the data. The reason for doing it is to reduce the overall variability in the  signal levels for each sample, which is generally assumed to be due to technical variability. What quantile normalization does is it calculates the average distribution for the signals across all of the samples (sort of like determining the average sample box plot shape), and then forces ranked data for each sample to take on that distribution.
# The sample-sample correlations are a way of measuring how similar the overall gene expression profile is for each sample to all of the other samples in the analysis.  It is performed following quantile normalization so that the correlation values don't simply reflect the generall overall signal for each sample. These correlations can be analyzed using a heatmap, boxplots, or principal components analysis (PCA).
 

##1. Highlight this code and hit run (Ctrl+R):

#***CODE TO RUN FOR STEP 9****

#Quantile normalize the filtered log-transformed signal data:
NormFiltered<-normalize.quantiles(log.avg_signal_filtered)
row.names(NormFiltered)<-row.names(log.avg_signal_filtered)
colnames(NormFiltered)<-colnames(log.avg_signal_filtered)

#Output a boxplot illustrating the (now identical) signal distributions for each sample:
png("09 Boxplot NormSignal Lateral Filtered.png")
boxplot(data.frame(NormFiltered), cex=0.25, las=3, par(cex.axis=0.75), main="Boxplot of normalized signal values per sample (1 box= all filtered probes)", xlab="Sample ID", ylab="Quantile Normalized Signal")
dev.off()

#Output the quantile normalized data:
write.table(NormFiltered, "Quantile Normalized filtered data.txt", sep="\t")

#Visualize the sample-sample correlations using a heatmap:
png("09 Sample Sample Correlations Heatmap.png")
image(cor(NormFiltered), main="Visualizing the correlations between entire samples (by index#)", xlab="Red=Less correlated, Light yellow=Highly correlated")
dev.off()
#Note that the heatmap can be tailored to focus on a certain level of correlation by using the command zlim=c(lower limit, upper limit)

#Visualize the sample-sample correlations using a boxplot:
png("09 Boxplot Sample Sample Correlations.png")
boxplot(data.frame(cor(NormFiltered)), cex=0.25, las=3, par(cex.axis=0.75), main="Boxplot of sample-sample correlations", xlab="Subject", ylab="Sample-Sample Correlations")
Median10thQuantile<-median(apply((cor(NormFiltered)), 1, quantile, 0.1))
MedianQuantile<-median(apply((cor(NormFiltered)), 1, quantile, 0.5))
abline(a=Median10thQuantile, b=0, col=2)
abline(a=MedianQuantile, b=0, col=3)
mtext(paste("Median Sample-Sample Correlation=", round(MedianQuantile, digits=3), sep=" ")) 
dev.off()

#Run principal components analysis (PCA) to determine which major gradients of sample-sample correlations exist in the data (i.e. who is similar to whom):
pcaNormFiltered<-prcomp(t(NormFiltered))
tmp<-pcaNormFiltered$x[,1:4]
write.table(tmp, "PCA_1_4.txt", sep="\t")

PC1<-pcaNormFiltered$x[,1]
PC2<-pcaNormFiltered$x[,2]

PC3<-pcaNormFiltered$x[,3]
PC4<-pcaNormFiltered$x[,4]

#Output a scree plot for the PCA:
png("09 PCA Scree Plot1.png")
plot(summary(pcaNormFiltered)$importance[2,]~(c(1:42)), main="Variance Explained by Each Principal Component", xlab="PC #", ylab="Proportion of Variance Explained", col=2)
dev.off()

png("09 PCA Scree Plot2.png")
plot(summary(pcaNormFiltered)$importance[3,]~(c(1:42)), main="Variance Explained by Each Principal Component", xlab="PC #", ylab="Cumulative Proportion of Variance Explained", col=3)
dev.off()

#Output a scatterplot illustrating the relationship between Principal components 1 & 2 (PC1 & PC2):
png("09 PC1 vs PC2.png")
plot(PC1~PC2, main="Principal Components Analysis of Normalized Filtered Data")
points(PC1[Diagnosis=="CTRL"]~PC2[Diagnosis=="CTRL"], col=3)
points(PC1[Diagnosis=="MDD"]~PC2[Diagnosis=="MDD"], col=2)
legend(min(PC1), max(PC2)+10, c("Control", "MDD"), text.col=c(3, 2), pch=19, col=c(3, 2))
dev.off()

#Output a scatterplot illustrating the relationship between Principal components 3 & 4 (PC3 & PC4):
png("09 PC3 vs PC4.png")
plot(PC3~PC4, main="Principal Components Analysis of Normalized Filtered Data")
points(PC3[Diagnosis=="CTRL"]~PC4[Diagnosis=="CTRL"], col=3)
points(PC3[Diagnosis=="MDD"]~PC4[Diagnosis=="MDD"], col=2)
legend(min(PC3), max(PC4)+10, c("Control", "MDD"), text.col=c(3, 2), pch=19, col=c(3, 2))
dev.off()

#***END OF CODE FOR STEP 9****

##2. Take a look at "09 Boxplot NormSignal Lateral Filtered.png", "09 Sample Sample Correlations Heatmap.png", "09 Boxplot Sample Sample Correlations.png", "09 PCA Scree Plot1.png", "09 PCA Scree Plot2.png", "09 PC1 vs PC2.png", "09 PC3 vs PC4.png". What to look for:
###A. In the "Boxplot NormSignal Lateral Filtered.png" are all of the boxplots now identical following quantile normalization?
###B. Is there any obvious pattern of sample-sample correlations that pops out when looking at the heatmap?  Are there a few samples that stand out as less correlated with the rest (dark red stripes)? (if so, these might be considered outliers)
###C. Likewise, looking at the "Boxplot Sample Sample Correlations.png", are there any samples that have a distribution of sample-sample correlation coefficients (box) that is much lower than the rest?  Do any of them surpass the threshold of having a distribution of correlation coefficients (box) that falls completely beneath the median 10th quantile for all samples? (if so, these might be considered outliers)
###D. Looking at "Eigenvalues vs PC.png", how many PCs do you need to include in your analysis to account for 80% of the variation in your sample-sample correlations? Are PC1-4 enough?
###E. Looking at "PC1 vs PC2.png" and "PC3 vs PC4.png", are there any samples that stand out as being quite separate from the rest? (if so, these might be considered outliers)  Do there appear to be any clusters of samples?  

#*******************************************************************

#STEP 10: Outlier removal: At this point it would probably be good to stop and check in with me!

##1.  In order to remove the outlier samples, you will need to know their column index in R.  To do this, you can enter their subject number into this code (replacing 3686) and then hit ctrl+R to run (the index # will appear in the console): 
##Example:
which("4232"==colnames(NormFiltered))
which("3204"==colnames(NormFiltered))
which("4424"==colnames(NormFiltered))

##2.  After coming up with a list of indices for the outliers, you can remove them by inserting their index numbers into this code, preceded by a minus sign (replacing -4, -6, and -13). If there aren't any outliers to remove, we will instead rename the original quantile normalized data to indicate that it has made it past the step of outlier removal:

###If there are outliers (example):
OutlierRemovalIndices<-c(-28, -29, -35) 

###If there aren't outliers, just highlight this code and hit Ctrl+R to run. Skip the rest of this section (Step 10).
normNoOutliers<-NormFiltered
DiagnosisNoOutliers<-Diagnosis
PC1noOutliers<-PC1
PC2noOutliers<-PC2
PC3noOutliers<-PC3
PC4noOutliers<-PC4


##3. If there were outliers, we need to re-quantile normalize the data set following their removal. To do this, highlight this code and hit Ctrl+R to run:


#***CODE TO RUN FOR STEP 10****

#Filtering the diagnosis information so that it doesn't include outlier subjects:
DiagnosisNoOutliers<-Diagnosis[OutlierRemovalIndices]

#Repeating quantile normalization with the outliers removed:
normNoOutliers<-normalize.quantiles(log.avg_signal_filtered[,OutlierRemovalIndices])
rownames(normNoOutliers)<-rownames(log.avg_signal_filtered)
colnames(normNoOutliers)<-colnames(log.avg_signal_filtered [ ,OutlierRemovalIndices])

#Outputting a summary of the outlier removal process:
OutlierRemoval<-file("10 Outlier Removal.txt")
out<-c(
print("These sample indices were removed as outliers:"),
print(OutlierRemovalIndices),
print("...which correspond to the sample ID numbers:"),
capture.output(
for(i in 1:length(OutlierRemovalIndices)){
	print(colnames(log.avg_signal_filtered[ ,OutlierRemovalIndices[i]]))
}
),
print("Therefore, the dimensions of the new data matrix are now:"),
capture.output(dim(normNoOutliers)),
print("Check column names for quantile normalized data with outliers removed"),
print(colnames(normNoOutliers))
)
cat(out, file="10 Outlier Removal.txt", sep="\n", append=TRUE)
close(OutlierRemoval)
rm(out)

#Output a (new) boxplot illustrating the (now identical) signal distributions for each sample:
png("10 Boxplot NormSignal Lateral Filtered No Outliers.png")
boxplot(data.frame(normNoOutliers), cex=0.25, las=3, par(cex.axis=0.75), main="No outliers: Boxplot of normalized signal values per sample (1 box= all filtered probes)", xlab="Sample ID", ylab="Quantile Normalized Signal")
dev.off()

#Output the quantile normalized data again:
write.table(normNoOutliers, "Quantile Normalized filtered data No Outliers.txt", sep="\t")

#Visualize the sample-sample correlations using a heatmap:
png("10 Sample Sample Correlations Heatmap No Outliers.png")
image(cor(normNoOutliers), main="No Outliers: Visualizing the correlations between entire samples (by index#)", xlab="Red=Less correlated, Light yellow=Highly correlated")
dev.off()
#Note that the heatmap can be tailored to focus on a certain level of correlation by using the command zlim=c(lower limit, upper limit)

#Visualize the sample-sample correlations using a boxplot:
png("10 Boxplot Sample Sample Correlations No Outliers.png")
boxplot(data.frame(cor(normNoOutliers)), cex=0.25, las=3, par(cex.axis=0.75))
Median10thQuantile<-median(apply((cor(normNoOutliers)), 1, quantile, 0.1))
abline(a=Median10thQuantile, b=0, col=2)
dev.off()

#Run principal components analysis (PCA) to determine which major gradients of sample-sample correlations exist in the data (i.e. who is similar to whom):
pcaNormFilterednoOutliers<-prcomp(t(normNoOutliers))
tmp<-pcaNormFilterednoOutliers$x[,1:4]
write.table(tmp, "PCA_1_4 No Outliers.txt", sep="\t")

PC1noOutliers<-pcaNormFilterednoOutliers$x[,1]
PC2noOutliers<-pcaNormFilterednoOutliers$x[,2]

PC3noOutliers<-pcaNormFilterednoOutliers$x[,3]
PC4noOutliers<-pcaNormFilterednoOutliers$x[,4]

#Output a scree plot for the PCA (no outliers):
png("10 PCA Scree Plot1 no outliers.png")
plot(summary(pcaNormFilterednoOutliers)$importance[2,]~(c(1:length(summary(pcaNormFilterednoOutliers)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC #", ylab="Proportion of Variance Explained", col=2)
dev.off()

png("10 PCA Scree Plot2 no outliers.png")
plot(summary(pcaNormFilterednoOutliers)$importance[3,]~(c(1:length(summary(pcaNormFilterednoOutliers)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC #", ylab="Cumulative Proportion of Variance Explained", col=3)
dev.off()

#Output a scatterplot illustrating the relationship between Principal components 1 & 2 (PC1 & PC2):
png("10 PC1 vs PC2 no outliers.png")
plot(PC1noOutliers~PC2noOutliers, main="Principal Components Analysis of Normalized Filtered Data No Outliers")
points(PC1noOutliers[DiagnosisNoOutliers=="CTRL"]~PC2noOutliers[DiagnosisNoOutliers=="CTRL"], col=3)
points(PC1noOutliers[DiagnosisNoOutliers=="MDD"]~PC2noOutliers[DiagnosisNoOutliers=="MDD"], col=2)
legend(min(PC1noOutliers), max(PC2noOutliers)+10, c("Control", "MDD"), text.col=c(3, 2), pch=19, col=c(3, 2))
dev.off()

#Output a scatterplot illustrating the relationship between Principal components 3 & 4 (PC3 & PC4):
png("10 PC3 vs PC4 no outliers.png")
plot(PC3noOutliers~PC4noOutliers, main="Principal Components Analysis of Normalized Filtered Data No Outliers")
points(PC3noOutliers[DiagnosisNoOutliers=="CTRL"]~PC4noOutliers[DiagnosisNoOutliers=="CTRL"], col=3)
points(PC3noOutliers[DiagnosisNoOutliers=="MDD"]~PC4noOutliers[DiagnosisNoOutliers=="MDD"], col=2)
legend(min(PC3noOutliers), max(PC4noOutliers)+10, c("Control", "MDD"), text.col=c(3, 2), pch=19, col=c(3, 2))
dev.off()

#Outputting the new quantile normalized data (no outliers):
write.table(normNoOutliers,"Quantile normalized filtered data no outliers.txt",sep="\t")



#Making prettier plots for the paper

pdf("10 PC1 vs PC2 no outliers.pdf", height=4.5, width=4)
plot(PC1noOutliers~PC2noOutliers, xlab="PC2", ylab="PC1", lwd=1.5, cex.lab=1.3, col=as.numeric(DiagnosisNoOutliers)*-1+4)
mtext("Medial", line=1, cex=1.5)
#legend(min(PC1noOutliers), max(PC2noOutliers)+10, c("Control", "MDD"), text.col=c(3, 2), pch=19, col=c(3, 2))
dev.off()

pdf("10 PCA Scree Plot1 no outliers.pdf", height=4.5, width=4)
plot(summary(pcaNormFilterednoOutliers)$importance[2,]~(c(1:length(summary(pcaNormFilterednoOutliers)$importance[2,]))), xlab="PC #", ylab="Proportion of Variance Explained", cex.lab=1.3, col=1)
mtext("Medial", line=1, cex=1.5)
dev.off()

pdf("PC2_ByDiagnosis.pdf", height=4.5, width=4)
boxplot(PC2noOutliers~DiagnosisNoOutliers, xlab="Diagnosis", ylab="PC2", lwd=1.5, cex.lab=1.3)
stripchart(PC2noOutliers~DiagnosisNoOutliers, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 20, cex=1.5, col = c(3,2))
mtext("Medial", line=1, cex=1.5)
#legend(min(PC1noOutliers), max(PC2noOutliers)+10, c("Control", "MDD"), text.col=c(3, 2), pch=19, col=c(3, 2))
dev.off()





#***END OF CODE FOR STEP 10**** 


##4. Take a look at "10 Outlier Removal.txt", "10 Boxplot NormSignal Lateral Filtered No Outliers.png", "10 Sample Sample Correlations Heatmap No Outliers.png", "10 Boxplot Sample Sample Correlations No Outliers.png", "10 Eigenvalues vs PC No Outliers.png", "10 PC1 vs PC2 no outliers.png", "10 PC3 vs PC4 no outliers.png". What to look for:
###A. In "10 Outlier Removal.txt", does it look like the final column names are correct? (i.e. that the correct outlier subjects were truly removed?)
###B. In the "10 Boxplot NormSignal Lateral Filtered No Outliers.png" are all of the boxplots now identical following quantile normalization?
###C. Is there any obvious pattern of sample-sample correlations that pops out when looking at the new heatmap?  Are there still a few samples that stand out as less correlated with the rest (dark red stripes)? 
###D. Likewise, looking at the "10 Boxplot Sample Sample Correlations No Outliers.png", are there any samples that have a distribution of sample-sample correlation coefficients (box) that is much lower than the rest?  Do any of them still surpass the threshold of having a distribution of correlation coefficients (box) that falls completely beneath the median 10th quantile for all samples? 
###E. Looking at "10 Eigenvalues vs PC no Outliers.png", how many PCs do you need to include in your analysis to account for 80% of the variation in your sample-sample correlations? Are PC1-4 enough?
###F. Looking at "10 PC1 vs PC2 no outliers.png" and "10 PC3 vs PC4 no outliers.png", are there any samples that still stand out as being quite separate from the rest?  Do there appear to be any clusters of samples?  

 
#*******************************************************************

#STEP 11:  Placing the subject information in the correct format to import into R.  At this point, you need to talk with me again. 

##1. Open the "Template_Sample ID file extended for R.xls" from the Genexp folder.
##2. Open the "Sample ID file.csv" from the brain regions folder.
##3. Highlight the information in "Sample ID file.csv" and paste it (paste special: values) into a blank sheet in "Template_Sample ID file extended for R.xls"
##4. Insert two new columns to the right of the column containing the subject indentification by the sample's chip and location on chip.
##5. Highlight the column containing the subject identification by the sample's chip and location on chip. This data will need to be split into one column for the chip information and one column containing the location on the chip using "Data"->"Text to Columns" (delimited by "_"). The rest can be deleted.
##6. Insert a new column to the left (column A) called Sample Processing Order. Simply number the series 1, 2, 3, etc...
##7. At this point, there should be 5 columns: 1) Sample Processing Order, 2) Chip, 3) Location on Chip, 4) Diagnosis, 5) Subject Number. If there are any empty columns mixed in there, delete them.
##8. Copy and paste the 5 columns into the first 5 columns of the "New Brain Region" worksheet. At this point *most* of the subject information should update, except for...
##9. RNA concentration and RNA integrity need their formulas updated for each brain region.  Simply change the column # in the VLOOKUP function to reflect the brain region appropriate column in the "ConcInteg" worksheet.
##10. Pull the new formula for RNA concentration and RNA Integrity down to fill the rest of the 42 subjects.
##11. Copy Worksheet 1 and Paste Special -> Paste Values into a new file.  Name this file "Sample ID file extended for R.csv".

Summary of new columns of information:
###Column 6 (F): Cohort
###Column 7: Gender
###Column 8: Age
###Column 9: RNAAge
###Column 10: TOD
###Column 11: Suicide
###Column 12: Overdose
###Column 13: Brain pH
###Column 14: RNA concentration for the brain region of interest
###Column 15: RNA integrity for the brain region of interest
###Column 16: Hours Cold
###Column 17: Hours Ice
###Column 18 (R): Post Mortem Interval (Hours Final)

###Agonal Factor was left out of the data read in, because it is missing a lot of values and all of the values that we do have are identical (0)
###Column 19: Agonal Factor


************************************************************

#***CODE TO RUN FOR STEP 12****

#Read in the detailed subject information:
SampleInfo<-as.matrix(read.csv("Sample ID file extended for R.csv", header=T))
row.names(SampleInfo)<-SampleInfo[,5]
SampleInfo[c(1:3), ]

#Double check that the row names for the sample information match the original column names for the signal information:
IsMisMatched<-matrix(0, 42, 1)
for(i in 1:42){
if(row.names(SampleInfo)[i]==colnames(NormFiltered)[i]){IsMisMatched<<-0}else{IsMisMatched<<-1}	
}
if(sum(IsMisMatched)>0){print("ERROR: Subjects are in a different order in the subject information and signal database")}else{print("All good: Subject order matches in both the subject information and signal databases")}


#Separate out the columns into individual variables:

SampleProcessingOrder<-as.numeric(SampleInfo[,1])
Chip<-as.factor(SampleInfo[,2])
LocationOnChip<-as.factor(SampleInfo[,3])

Cohort<-as.factor(SampleInfo[,6])
Gender<-as.factor(SampleInfo[,7])
Gender<-relevel(Gender, ref="M")
Age<-as.numeric(SampleInfo[,8])
RNAAge<-as.numeric(SampleInfo[,9])
TOD<-as.numeric(SampleInfo[,10])
Suicide<-as.factor(SampleInfo[,11])
Overdose<-as.factor(SampleInfo[,12])
BrainpH<-as.numeric(SampleInfo[,13])
RNAconc<-as.numeric(SampleInfo[,14])
RNAintegrity<-as.numeric(SampleInfo[,15])
HoursCold<-as.numeric(SampleInfo[,16])
HoursIce<-as.numeric(SampleInfo[,17])
HoursFinal<-as.numeric(SampleInfo[,18])
HoursColdIce<-as.numeric(SampleInfo[,19])
HoursFinalCorrected<-as.numeric(SampleInfo[,20])

SubjectReadIn<-file("12 Reading in Subject Information.txt")
out<-c(
print("The Dimensions of the Detailed Subject Information Matrix:"),
capture.output(dim(SampleInfo)),
capture.output(if(sum(IsMisMatched)>0){print("ERROR: Subjects are in a different order in the subject information and signal database")}else{print("All good: Subject order matches in both the subject information and signal databases")}),
print("Here are each of the new subject variables with information about their output:"),
print("Chip:"),
capture.output(str(Chip)),
print("LocationOnChip:"),
capture.output(str(LocationOnChip)),
print("Cohort:"),
capture.output(str(Cohort)),
print("Gender:"),
capture.output(str(Gender)),
print("Age:"),
capture.output(str(Age)),
print("RNAAge:"),
capture.output(str(RNAAge)),
print("TOD:"),
capture.output(str(TOD)),
print("Suicide:"),
capture.output(str(Suicide)),
print("Overdose:"),
capture.output(str(Overdose)),
print("BrainpH:"),
capture.output(str(BrainpH)),
print("RNAconc:"),
capture.output(str(RNAconc)),
print("RNAintegrity:"),
capture.output(str(RNAintegrity)),
print("HoursCold:"),
capture.output(str(HoursCold)),
print("HoursIce:"),
capture.output(str(HoursIce)),
print("HoursFinal:"),
capture.output(str(HoursFinal)),
print("HoursColdIce:"),
capture.output(str(HoursColdIce)),
print("HoursFinalCorrected:"),
capture.output(str(HoursFinalCorrected))

)
cat(out, file="12 Reading in Subject Information.txt", sep="\n", append=TRUE)
close(SubjectReadIn)
rm(out)

#***END OF CODE FOR STEP 12****

##2. Take a look at "12 Reading in Subject Information.txt" and see if everything looks alright.  What to look for:
###A. Are there any error messages?
###B. Do each of the new variables have values that look appropriate?


#*******************************************************************

#STEP 13: Removing outliers from the detailed subject information:

##1. If there weren't any outliers, then just run this code by highlighting it and hitting ctrl+R:

#***CODE TO RUN FOR STEP 13 IF THERE AREN'T OUTLIERS****

SampleProcessingOrderNoOutliers<-SampleProcessingOrder
ChipNoOutliers<-Chip
LocationOnChipNoOutliers<-LocationOnChip
CohortNoOutliers<-Cohort
GenderNoOutliers<-Gender
AgeNoOutliers<-Age
RNAAgeNoOutliers<-RNAAge
TODNoOutliers<-TOD
SuicideNoOutliers<-Suicide
OverdoseNoOutliers<-Overdose
BrainpHNoOutliers<-BrainpH
RNAconcNoOutliers<-RNAconc
RNAintegrityNoOutliers<-RNAintegrity
HoursColdNoOutliers<-HoursCold
HoursIceNoOutliers<-HoursIce
HoursFinalNoOutliers<-HoursFinal
HoursColdIceNoOutliers<-HoursColdIce
HoursFinalCorrectedNoOutliers<-HoursFinalCorrected


#***END OF CODE TO RUN FOR STEP 13 IF THERE AREN'T OUTLIERS****




##2. If there were outliers, then run this code by highlighting it and hitting ctrl+R:

#***CODE TO RUN FOR STEP 13 IF THERE ARE OUTLIERS****

SampleInfoNoOutliers<-SampleInfo[OutlierRemovalIndices,]
cbind(colnames(normNoOutliers), SampleInfoNoOutliers[,5])
write.csv(SampleInfoNoOutliers, "SampleInfoNoOutliers.csv")


SampleProcessingOrderNoOutliers<-SampleProcessingOrder[OutlierRemovalIndices]
ChipNoOutliers<-Chip[OutlierRemovalIndices]
LocationOnChipNoOutliers<-LocationOnChip[OutlierRemovalIndices]
CohortNoOutliers<-Cohort[OutlierRemovalIndices]
GenderNoOutliers<-Gender[OutlierRemovalIndices]
AgeNoOutliers<-Age[OutlierRemovalIndices]
RNAAgeNoOutliers<-RNAAge[OutlierRemovalIndices]
TODNoOutliers<-TOD[OutlierRemovalIndices]
SuicideNoOutliers<-Suicide[OutlierRemovalIndices]
OverdoseNoOutliers<-Overdose[OutlierRemovalIndices]
BrainpHNoOutliers<-BrainpH[OutlierRemovalIndices]
RNAconcNoOutliers<-RNAconc[OutlierRemovalIndices]
RNAintegrityNoOutliers<-RNAintegrity[OutlierRemovalIndices]
HoursColdNoOutliers<-HoursCold[OutlierRemovalIndices]
HoursIceNoOutliers<-HoursIce[OutlierRemovalIndices]
HoursFinalNoOutliers<-HoursFinal[OutlierRemovalIndices]
HoursColdIceNoOutliers<-HoursColdIce[OutlierRemovalIndices]
HoursFinalCorrectedNoOutliers<-HoursFinalCorrected[OutlierRemovalIndices]

#***END OF CODE TO RUN FOR STEP 13 IF THERE ARE OUTLIERS****


#*******************************************************************

STEP 14:  Looking at the relationships between the independent variables

##1. To examine both visual and statistical relationships between the independent variables (with the outlier subjects removed, if there were any!), just highlight this code and click Ctrl+R:


#***CODE TO RUN FOR STEP 14****

SubjectFactorVariables<-cbind(DiagnosisNoOutliers, ChipNoOutliers, LocationOnChipNoOutliers, CohortNoOutliers, GenderNoOutliers, SuicideNoOutliers, OverdoseNoOutliers)

SubjectContinuousVariables<-cbind(AgeNoOutliers, RNAAgeNoOutliers, TODNoOutliers, BrainpHNoOutliers, RNAconcNoOutliers, RNAintegrityNoOutliers, HoursColdNoOutliers, HoursIceNoOutliers, HoursFinalNoOutliers, HoursColdIceNoOutliers, HoursFinalCorrectedNoOutliers)

#Outputting histograms of the continuous variables

for (i in 1:length(SubjectContinuousVariables[1,])){
png(paste(paste("14 Histogram of", colnames(SubjectContinuousVariables)[i], sep="  "), "png", sep="."))	
hist(SubjectContinuousVariables[,i], col=i+1)
dev.off()			
}


#Using a scatterplot with best fit line to visually examine the relationships between the continuous subject variables:
for (i in 1:length(SubjectContinuousVariables[1,])){
for(j in 1:length(SubjectContinuousVariables[1,])){
png(paste("14", paste(colnames(SubjectContinuousVariables)[i], "vs", colnames(SubjectContinuousVariables)[j], sep="  "), "png", sep="."))	
plot(SubjectContinuousVariables[,i]~SubjectContinuousVariables[,j], main=paste(colnames(SubjectContinuousVariables)[i], "vs", colnames(SubjectContinuousVariables)[j], sep="  "), xlab=colnames(SubjectContinuousVariables)[j], ylab=colnames(SubjectContinuousVariables)[i])
RegressionLine<-lm(SubjectContinuousVariables[,i]~SubjectContinuousVariables[,j])
abline(RegressionLine, col=2)
mtext(paste("p-value=", round(summary.lm(RegressionLine)$coefficients[8], digits=4)))
dev.off()		
}		
}



#Using boxplots to visually examine the relationships between the continuous subject variables and categorical subject variables:
for (i in 1:length(SubjectContinuousVariables[1,])){
for(j in 1:length(SubjectFactorVariables[1,])){
png(paste("14", paste(colnames(SubjectContinuousVariables)[i], "vs", colnames(SubjectFactorVariables)[j], sep="  "), "png", sep="."))	
boxplot(SubjectContinuousVariables[,i]~SubjectFactorVariables[,j], main=paste(colnames(SubjectContinuousVariables)[i], "vs", colnames(SubjectFactorVariables)[j], sep="  "), xlab=colnames(SubjectFactorVariables)[j], ylab=colnames(SubjectContinuousVariables)[i])
mtext(paste("p-value=", round(summary(aov(SubjectContinuousVariables[,i]~SubjectFactorVariables[,j]))[[1]][["Pr(>F)"]][1], digits=4)))
dev.off()		
}		
}





#Creating a text file of contingency tables to visually examine the relationships between categorical subject variables:
CrossTabsIV<-file("14 Cross Tabs Between Subject Factors.txt")
out<-c(
capture.output(
for (i in 1:length(SubjectFactorVariables[1,])){
for(j in 1:length(SubjectFactorVariables[1,])){
ContingencyTable<-table(SubjectFactorVariables[,i],SubjectFactorVariables[,j])
print(paste(colnames(SubjectFactorVariables)[i], "vs", colnames(SubjectFactorVariables)[j], sep="  "))
print(ContingencyTable)
print(paste("p-value=", chisq.test(ContingencyTable)$p.value))	
}		
}
)
)
cat(out, file="14 Cross Tabs Between Subject Factors.txt", sep="\n", append=TRUE)
close(CrossTabsIV)
rm(out)


#Outputting a text file containing the statistical relationships between all of the subject variables:

StatisticalRelationshipsIV<-file("14 Statistical Relationships between Subject Variables.txt")
out<-c(

capture.output(

summary.lm(lm(RNAconcNoOutliers~AgeNoOutliers+BrainpHNoOutliers+DiagnosisNoOutliers+GenderNoOutliers+RNAintegrityNoOutliers+RNAAgeNoOutliers+TODNoOutliers+HoursColdNoOutliers+HoursIceNoOutliers+HoursFinalNoOutliers+HoursFinalCorrectedNoOutliers))

),


capture.output(

summary.lm(lm(RNAintegrityNoOutliers~AgeNoOutliers+BrainpHNoOutliers+DiagnosisNoOutliers+GenderNoOutliers+RNAconcNoOutliers+RNAAgeNoOutliers+TODNoOutliers+HoursColdNoOutliers+HoursIceNoOutliers+HoursFinalNoOutliers+HoursFinalCorrectedNoOutliers))

),


capture.output(
#Calculating the variance inflation factor (vif) to determine which subject variables are highly related to other subject variables in the data set. Most important, of course, is whether any of the subject variables strongly correlate with Diagnosis. Note that "Location on Chip" has been removed as a variable because it is partially redundant with gender. 
vif(lm(normNoOutliers[1,]~AgeNoOutliers+BrainpHNoOutliers+ChipNoOutliers+DiagnosisNoOutliers+GenderNoOutliers+RNAconcNoOutliers+RNAintegrityNoOutliers+RNAAgeNoOutliers+TODNoOutliers+HoursColdNoOutliers+HoursIceNoOutliers+HoursFinalNoOutliers+HoursFinalCorrectedNoOutliers))

#It might make sense to remove RNAAge(because of the big outlier), or HoursCold/HoursIce since they seem to (almost) add up to HoursFinal... or agonal factor, depending on whether there is actually any subjects with agonal factor>0

#It might be nice to run a hierarchical cluster amongst these variables also...
),

#Using linear regression to examine the statistical relationships between the continuous subject variables:

capture.output(
for (i in 1:length(SubjectContinuousVariables[1,])){
for(j in 1:length(SubjectContinuousVariables[1,])){
print(paste(colnames(SubjectContinuousVariables)[i], "vs", colnames(SubjectContinuousVariables)[j], sep="  "))
print(summary.lm(lm(SubjectContinuousVariables[,i]~SubjectContinuousVariables[,j])))
}		
}
),
#Using anova to examine the statistical relationships between the continuous subject variables and categorical subject variables:

capture.output(
for (i in 1:length(SubjectContinuousVariables[1,])){
for(j in 1:length(SubjectFactorVariables[1,])){
print(paste(colnames(SubjectContinuousVariables)[i], "vs", colnames(SubjectFactorVariables)[j], sep="  "))
print(summary(aov(SubjectContinuousVariables[,i]~SubjectFactorVariables[,j])))		
}		
}
),

#Using chi-square to examine the statistical relationships between the categorical subject variables:

capture.output(
for (i in 1:length(SubjectFactorVariables[1,])){
for(j in 1:length(SubjectFactorVariables[1,])){
print(paste(colnames(SubjectFactorVariables)[i], "vs", colnames(SubjectFactorVariables)[j], sep="  "))
print(chisq.test(ContingencyTable))		
}		
}
)

)
cat(out, file="14 Statistical Relationships between Subject Variables.txt", sep="\n", append=TRUE)
close(StatisticalRelationshipsIV)
rm(out)


#Flagging variables that are collinear with other subject variables:
FlaggedRelationshipsBetweenIV<-file("14 Flagged Relationships Between Subject Variables.txt")
out<-c(

#Using linear regression to examine the statistical relationships between the continuous subject variables:
capture.output(
for (i in 1:length(SubjectContinuousVariables[1,])){
for(j in 1:length(SubjectContinuousVariables[1,])){
if(summary.lm(lm(SubjectContinuousVariables[,i]~SubjectContinuousVariables[,j]))$coefficient[8]<0.05){
print(paste(colnames(SubjectContinuousVariables)[i], "vs", colnames(SubjectContinuousVariables)[j], "p-value=", summary.lm(lm(SubjectContinuousVariables[,i]~SubjectContinuousVariables[,j]))$coefficient[8], sep="  "))}else{}
}		
}
),

#Using anova to examine the statistical relationships between the continuous subject variables and categorical subject variables:
capture.output(
for (i in 1:length(SubjectContinuousVariables[1,])){
for(j in 1:length(SubjectFactorVariables[1,])){
if(summary(aov(SubjectContinuousVariables[,i]~SubjectFactorVariables[,j]))[[1]][["Pr(>F)"]][1]<0.05){
print(paste(colnames(SubjectContinuousVariables)[i], "vs", colnames(SubjectFactorVariables)[j], "p-value=", summary(aov(SubjectContinuousVariables[,i]~SubjectFactorVariables[,j]))[[1]][["Pr(>F)"]][1], sep="  "))	
}else{}		
}		
}
),

#Using chi-square to examine the statistical relationships between the categorical subject variables:
capture.output(
for (i in 1:length(SubjectFactorVariables[1,])){
for(j in 1:length(SubjectFactorVariables[1,])){
ContingencyTable<-table(SubjectFactorVariables[,i], SubjectFactorVariables[,j])
if(chisq.test(ContingencyTable)$p.value<0.05){
print(paste(colnames(SubjectFactorVariables)[i], "vs", colnames(SubjectFactorVariables)[j], "p-value=", chisq.test(ContingencyTable)$p.value, sep="  "))	
}else{}		
}		
}
)
)
cat(out, file="14 Flagged Relationships Between Subject Variables.txt", sep="\n", append=TRUE)
close(FlaggedRelationshipsBetweenIV)
rm(out)




#***END OF CODE FOR STEP 14****


##2. Examine all of the step 14 output.  What to look for:
###A. Briefly look at the first input on the "14 Statistical Relationships between Subject Variables.txt" file. This is the Variance Inflation Factor and reflects the overall multicollinearity of the variables.  What is the VIF for diagnosis?
###B. Start with "14 Flagged Relationships Between Subject Variables.txt" - which subject variables are collinear with other subject variables? Is diagnosis collinear with anything?
###C. For those variables that were flagged, take a look at their output figures and statistics: Do you think the relationship is real or is it driven by outlier samples?
###D. Note that TOD has been treated as a linear variable for this analysis even though it is actually angular. Look at the figures for TOD to determine if there is an angular relationship, and if so come back and run the appropriate statistics.

#*******************************************************************

#STEP 15: Looking at the relationships between independent variables and principal components:

##1. To examine both visual and statistical relationships between the independent variables (with the outlier subjects removed, if there were any!) and the principal components of variation in the samples, just highlight this code and click Ctrl+R:


#***CODE TO RUN FOR STEP 15****

SubjectPCA<-cbind(PC1noOutliers, PC2noOutliers, PC3noOutliers, PC4noOutliers)



#Using a scatterplot with best fit line to visually examine the relationships between the continuous subject variables and SubjectPCA:
for (i in 1:length(SubjectPCA[1,])){
for(j in 1:length(SubjectContinuousVariables[1,])){
png(paste("15", paste(colnames(SubjectPCA)[i], "vs", colnames(SubjectContinuousVariables)[j], sep="  "), "png", sep="."))	
plot(SubjectPCA[,i]~SubjectContinuousVariables[,j], main=paste(colnames(SubjectPCA)[i], "vs", colnames(SubjectContinuousVariables)[j], sep="  "), xlab=colnames(SubjectContinuousVariables)[j], ylab=colnames(SubjectPCA)[i])
RegressionLine<-lm(SubjectPCA[,i]~SubjectContinuousVariables[,j])
abline(RegressionLine, col=2)
mtext(paste("p-value=", round(summary.lm(RegressionLine)$coefficients[8], digits=4)))
dev.off()		
}		
}



#Using boxplots to visually examine the relationships between the PCA and categorical subject variables:
for (i in 1:length(SubjectPCA[1,])){
for(j in 1:length(SubjectFactorVariables[1,])){
png(paste("15", paste(colnames(SubjectPCA)[i], "vs", colnames(SubjectFactorVariables)[j], sep="  "), "png", sep="."))	
boxplot(SubjectPCA[,i]~SubjectFactorVariables[,j], main=paste(colnames(SubjectPCA)[i], "vs", colnames(SubjectFactorVariables)[j], sep="  "), xlab=colnames(SubjectFactorVariables)[j], ylab=colnames(SubjectPCA)[i])
mtext(paste("p-value=", round(summary(aov(SubjectPCA[,i]~SubjectFactorVariables[,j]))[[1]][["Pr(>F)"]][1], digits=4)))
dev.off()		
}		
}






#Outputting a text file containing the statistical relationships between all of the subject variables and PCA:

StatisticalRelationshipsIVvsPCA<-file("15 Statistical Relationships between Subject Variables and PCA.txt")
out<-c(

capture.output(
#Calculating the variance inflation factor (vif) to determine which subject variables are highly related to other subject variables in the data set. Most important, of course, is whether any of the subject variables strongly correlate with Diagnosis. Note that "Location on Chip" has been removed as a variable because it is partially redundant with gender. 
summary.lm(lm(PC1noOutliers~AgeNoOutliers+BrainpHNoOutliers+ChipNoOutliers+DiagnosisNoOutliers+GenderNoOutliers+RNAconcNoOutliers+RNAintegrityNoOutliers+RNAAgeNoOutliers+TODNoOutliers+HoursColdNoOutliers+HoursIceNoOutliers+HoursFinalNoOutliers+HoursFinalCorrectedNoOutliers))
),

capture.output(
summary.lm(lm(PC2noOutliers~AgeNoOutliers+BrainpHNoOutliers+ChipNoOutliers+DiagnosisNoOutliers+GenderNoOutliers+RNAconcNoOutliers+RNAintegrityNoOutliers+RNAAgeNoOutliers+TODNoOutliers+HoursColdNoOutliers+HoursIceNoOutliers+HoursFinalNoOutliers+HoursFinalCorrectedNoOutliers))
),

capture.output(
summary.lm(lm(PC3noOutliers~AgeNoOutliers+BrainpHNoOutliers+ChipNoOutliers+DiagnosisNoOutliers+GenderNoOutliers+RNAconcNoOutliers+RNAintegrityNoOutliers+RNAAgeNoOutliers+TODNoOutliers+HoursColdNoOutliers+HoursIceNoOutliers+HoursFinalNoOutliers+HoursFinalCorrectedNoOutliers))
),

capture.output(
summary.lm(lm(PC4noOutliers~AgeNoOutliers+BrainpHNoOutliers+ChipNoOutliers+DiagnosisNoOutliers+GenderNoOutliers+RNAconcNoOutliers+RNAintegrityNoOutliers+RNAAgeNoOutliers+TODNoOutliers+HoursColdNoOutliers+HoursIceNoOutliers+HoursFinalNoOutliers+HoursFinalCorrectedNoOutliers))
),




capture.output(
#Calculating the variance inflation factor (vif) to determine which subject variables are highly related to other subject variables in the data set. Most important, of course, is whether any of the subject variables strongly correlate with Diagnosis. Note that "Location on Chip" has been removed as a variable because it is partially redundant with gender. 
vif(lm(normNoOutliers[1,]~AgeNoOutliers+BrainpHNoOutliers+ChipNoOutliers+DiagnosisNoOutliers+GenderNoOutliers+RNAconcNoOutliers+RNAintegrityNoOutliers+RNAAgeNoOutliers+TODNoOutliers+HoursColdNoOutliers+HoursIceNoOutliers+HoursFinalNoOutliers+PC1noOutliers+HoursFinalCorrectedNoOutliers))
),

capture.output(
vif(lm(normNoOutliers[1,]~AgeNoOutliers+BrainpHNoOutliers+ChipNoOutliers+DiagnosisNoOutliers+GenderNoOutliers+RNAconcNoOutliers+RNAintegrityNoOutliers+RNAAgeNoOutliers+TODNoOutliers+HoursColdNoOutliers+HoursIceNoOutliers+HoursFinalNoOutliers+PC2noOutliers+HoursFinalCorrectedNoOutliers))
),

capture.output(
vif(lm(normNoOutliers[1,]~AgeNoOutliers+BrainpHNoOutliers+ChipNoOutliers+DiagnosisNoOutliers+GenderNoOutliers+RNAconcNoOutliers+RNAintegrityNoOutliers+RNAAgeNoOutliers+TODNoOutliers+HoursColdNoOutliers+HoursIceNoOutliers+HoursFinalNoOutliers+PC3noOutliers+HoursFinalCorrectedNoOutliers))
),

capture.output(
vif(lm(normNoOutliers[1,]~AgeNoOutliers+BrainpHNoOutliers+ChipNoOutliers+DiagnosisNoOutliers+GenderNoOutliers+RNAconcNoOutliers+RNAintegrityNoOutliers+RNAAgeNoOutliers+TODNoOutliers+HoursColdNoOutliers+HoursIceNoOutliers+HoursFinalNoOutliers+PC4noOutliers+HoursFinalCorrectedNoOutliers))
),

#It might make sense to remove RNAAge(because of the big outlier), or HoursCold/HoursIce since they seem to (almost) add up to HoursFinal... or agonal factor, depending on whether there is actually any subjects with agonal factor>0
#It might be nice to run a hierarchical cluster amongst these variables also...


#Using linear regression to examine the statistical relationships between PCA and the continuous subject variables:

capture.output(
for (i in 1:length(SubjectPCA[1,])){
for(j in 1:length(SubjectContinuousVariables[1,])){
print(paste(colnames(SubjectPCA)[i], "vs", colnames(SubjectContinuousVariables)[j], sep="  "))
print(summary.lm(lm(SubjectPCA[,i]~SubjectContinuousVariables[,j])))
}		
}
),

#Using anova to examine the statistical relationships between PCA and categorical subject variables:

capture.output(
for (i in 1:length(SubjectPCA[1,])){
for(j in 1:length(SubjectFactorVariables[1,])){
print(paste(colnames(SubjectPCA)[i], "vs", colnames(SubjectFactorVariables)[j], sep="  "))
print(summary(aov(SubjectPCA[,i]~SubjectFactorVariables[,j])))		
}		
}
)

)
cat(out, file="15 Statistical Relationships between Subject Variables and PCA.txt", sep="\n", append=TRUE)
close(StatisticalRelationshipsIVvsPCA)
rm(out)



StatisticalRelationshipsIVvsPCA_2<-file("15 More Statistical Relationships between Subject Variables and PCA.txt")
out<-c(
  
  capture.output(
    
    summary.lm(lm(PC1noOutliers~DiagnosisNoOutliers))
  ),
  
  capture.output(
    summary.lm(lm(PC2noOutliers~DiagnosisNoOutliers))
  ),
  
  capture.output(
    summary.lm(lm(PC3noOutliers~DiagnosisNoOutliers))
  ),
  
  capture.output(
    summary.lm(lm(PC4noOutliers~DiagnosisNoOutliers))
  ),
  
  capture.output(
    summary.lm(lm(PC1noOutliers~DiagnosisNoOutliers+AgeNoOutliers+BrainpHNoOutliers+GenderNoOutliers+RNAintegrityNoOutliers+HoursFinalCorrectedNoOutliers))
  ),
  
  capture.output(
    summary.lm(lm(PC2noOutliers~DiagnosisNoOutliers+AgeNoOutliers+BrainpHNoOutliers+GenderNoOutliers+RNAintegrityNoOutliers+HoursFinalCorrectedNoOutliers))
  ),
  
  capture.output(
    summary.lm(lm(PC3noOutliers~DiagnosisNoOutliers+AgeNoOutliers+BrainpHNoOutliers+GenderNoOutliers+RNAintegrityNoOutliers+HoursFinalCorrectedNoOutliers))
  ),
  
  capture.output(
    summary.lm(lm(PC4noOutliers~DiagnosisNoOutliers+AgeNoOutliers+BrainpHNoOutliers+GenderNoOutliers+RNAintegrityNoOutliers+HoursFinalCorrectedNoOutliers))
  )
)
cat(out, file="15 More Statistical Relationships between Subject Variables and PCA.txt", sep="\n", append=TRUE)
close(StatisticalRelationshipsIVvsPCA_2)
rm(out)





#Flagging variables that are collinear with other subject variables:
FlaggedRelationshipsBetweenIVandPCA<-file("15 Flagged Relationships Between Subject Variables and PCA.txt")
out<-c(

#Using linear regression to examine the statistical relationships between the continuous subject variables:
capture.output(
for (i in 1:length(SubjectPCA[1,])){
for(j in 1:length(SubjectContinuousVariables[1,])){
if(summary.lm(lm(SubjectPCA[,i]~SubjectContinuousVariables[,j]))$coefficient[8]<0.05){
print(paste(colnames(SubjectPCA)[i], "vs", colnames(SubjectContinuousVariables)[j], "p-value=", summary.lm(lm(SubjectPCA[,i]~SubjectContinuousVariables[,j]))$coefficient[8], sep="  "))}else{}
}		
}
),

#Using anova to examine the statistical relationships between the continuous subject variables and categorical subject variables:
capture.output(
for (i in 1:length(SubjectPCA[1,])){
for(j in 1:length(SubjectFactorVariables[1,])){
if(summary(aov(SubjectPCA[,i]~SubjectFactorVariables[,j]))[[1]][["Pr(>F)"]][1]<0.05){
print(paste(colnames(SubjectPCA)[i], "vs", colnames(SubjectFactorVariables)[j], "p-value=", summary(aov(SubjectPCA[,i]~SubjectFactorVariables[,j]))[[1]][["Pr(>F)"]][1], sep="  "))	
}else{}		
}		
}
)
)
cat(out, file="15 Flagged Relationships Between Subject Variables and PCA.txt", sep="\n", append=TRUE)
close(FlaggedRelationshipsBetweenIVandPCA)
rm(out)


#***END OF CODE FOR STEP 15****


##2. Examine all of the step 15 output.  What to look for:
###A. Start with "15 Flagged Relationships Between Subject Variables and PCA.txt" - which subject variables are collinear with the principal components of variation in the sample? Is diagnosis collinear with any of the PC?
###B. For those variables that were flagged, take a look at their output figures and statistics: Do you think the relationship is real or is it driven by outlier samples?
###C. Note that TOD has been treated as a linear variable for this analysis even though it is actually angular. Look at the figures for TOD to determine if there is an angular relationship, and if so come back and run the appropriate statistics.
###D. Those variables which appear to have real relationships with the principal components should be included in STEP 17 of the analysis ("Calculating the effects of diagnosis in a model that controls for potential confounding variables").

#*************************************************

#Before starting diagnosis models:

##Outputting files combining p-values with basic stats about the probes:

Detection_Pvalue_filtered<-Detection_Pvalue[detect05.gene>8,]
Detection_Pvalue_filteredNoOutliers<-Detection_Pvalue_filtered
Detection_Pvalue_filteredNoOutliers<-Detection_Pvalue_filtered[,OutlierRemovalIndices]

Detection_PvalueMedianMDD<-apply(Detection_Pvalue_filteredNoOutliers[,DiagnosisNoOutliers=="MDD"], 1, median)
Detection_PvalueMedianCTRL<-apply(Detection_Pvalue_filteredNoOutliers[,DiagnosisNoOutliers=="CTRL"], 1, median)
plot(log(Detection_PvalueMedianMDD)~log(Detection_PvalueMedianCTRL))

NormNoOutliersMedianMDD<-apply(normNoOutliers[,DiagnosisNoOutliers=="MDD"], 1, median)
NormNoOutliersMedianCTRL<-apply(normNoOutliers[,DiagnosisNoOutliers=="CTRL"], 1, median)

plot(NormNoOutliersMedianMDD~NormNoOutliersMedianCTRL)

ProbeInfoOutput<-read.csv("ProbeInfoOutput.csv")
colnames(ProbeInfoOutput)

#*******************************************************************

#STEP 16:  Running Simple T-Tests by Diagnosis and Calculating Fold-Change:

##1. To compare gene expression between MDD and Control Subjects (*outliers removed) using a simple t-test (i.e. without controlling for confounding variables) just highlight this code and click Ctrl+R:


#***CODE TO RUN FOR STEP 16****


#I chose to run a Welch's t-test because it allows for unequal variance between groups. This was important because there are more MDD subjects (26) than Control Subjects (16):

SimpleTstat<-matrix(0, length(normNoOutliers[,1]), 1)
row.names(SimpleTstat)<-row.names(normNoOutliers)

SimpleTDf<-matrix(0, length(normNoOutliers[,1]), 2)
row.names(SimpleTDf)<-row.names(normNoOutliers)

SimplePval<-matrix(0, length(normNoOutliers[,1]), 1)
row.names(SimplePval)<-row.names(normNoOutliers)

for(i in 1:length(normNoOutliers[,1])){
	TempT<-t.test(normNoOutliers[i,]~DiagnosisNoOutliers)
	SimpleTstat[i,1]<-TempT$statistic
	SimpleTDf[i,]<-TempT$parameter
	SimplePval[i,1]<-TempT$p.value
}

#Outputting a figure illustrating the raw p-values vs. what would be expected by random chance:
png("16 Histogram of Unadjusted Pvalues Using Simple TTest.png")
hist(SimplePval, breaks=100, col=2, main="Unadjusted P-values Using Simple MDD vs. Control T-Test", xlab="Unadjusted p-value", ylab="Count")
abline(a=(length(normNoOutliers[,1])/100), b=0)
dev.off()

#At the moment, I skipped doing q-q plots because with an uneven sample size, Welch's test should produce different df for each gene.

#Applying two different types of multiple-comparison corrections to the raw p-values (Benjamini-Hochberg and Benjamini-Yekutieli):
TempPvalAdj<-mt.rawp2adjp(SimplePval[,1], proc=c("BH", "BY"))
SimplePvalAdj<-TempPvalAdj$adjp[order(TempPvalAdj$index),]
row.names(SimplePvalAdj)<-row.names(normNoOutliers)

#Calculating fold change:
FoldChange<-apply(NormFiltered[, Diagnosis=="MDD"], 1, mean)-apply(NormFiltered[, Diagnosis=="CTRL"], 1, mean)

#Outputting the p-values, fold change, and gene identifiers:
SimpleTtestOutput<-cbind(SimplePvalAdj, FoldChange, GeneIdentifiers_filtered)
row.names(SimpleTtestOutput)<-row.names(normNoOutliers)
colnames(SimpleTtestOutput)<-c("Unadjusted P-value", "BH Adjusted Pvalue", "BY Adjusted Pvalue", "Fold Change", "Symbol", "Search Key", "ILMN_Gene", "Chromosome", "Definition", "Synonyms", "ProbeID", "REFSEQ_ID", "UNIGENE_ID", "ENTREZ_GENE_ID", "PROBE_TYPE", "PROBE_START", "PROBE_SEQUENCE","ONTOLOGY_COMPONENT","ONTOLOGY_PROCESS", "ONTOLOGY_FUNCTION")

write.csv(SimpleTtestOutput, "SimpleTtestOutput.csv")


#***END OF CODE FOR STEP 16****


##2. Take a look at the "16 Histogram of Unadjusted Pvalues Using Simple Ttest". Does it look like there are a disproportionate number of genes with low pvalues?
##3. The output for the analysis using a simple t-test is in "SimpleTtestOutput.csv." When you open it in Excel, resave it as an Excel file. Then highlight the full data set using the corner triangle and sort by BH pvalue.



#*******************************************************************

#STEP 16:  Running a linear model to examine the effects of diagnosis while accounting for potential confounding factors: 

##1. The confounding factors should be chosen based on previous literature or based on which variables were flagged as being related to the PC's or diagnosis.
###For example, in the lateral amygdala, these variables are related to diagnosis:
####A. RNAconcNoOutliers
####B. SuicidenoOutliers (duh!)
####C. HoursFinalNoOutliers (almost significant)

###And these variables are related to the PCs:
####A. BrainpHNoOutliers (PC1)
####B. RNAconcNoOutliers  (PC2) (also correlates with diagnosis)
####C. DiagnosisNoOutliers (PC2)
####D. SuicideNoOutliers(PC2)
####E. SuicideNoOutliers (PC3)
####F. RNAconcNoOutliers (PC3)
####G. RNAAgeNoOutliers (PC4)


##Note that not all of these variables are worth including in the model. 
##For example, it is probably worth running suicide or overdose as a separate analysis since I believe all suicides and overdoses are MDD subjects.
##HoursCold and RNAAge both have extreme values that make them prone to producing false positives in statistical tests. They are also highly correlated with Cohort.


##2. Insert those variables into this code, placing diagnosis as the first variable.

LinearModel1<-function(i){lm(normNoOutliers[i,]~DiagnosisNoOutliers+AgeNoOutliers+RNAintegrityNoOutliers+RNAconcNoOutliers+BrainpHNoOutliers+GenderNoOutliers+HoursFinalCorrectedNoOutliers)}

##3. Then run this code (Ctrl+R) to get the p-values:


#***CODE TO RUN FOR STEP 17****

#Determining how many variables are in the linear model (since LM breaks factors down into dummy variables)
NumberofXsLM1<-length(summary.lm(LinearModel1(1))$coefficients[,1])
NameofXsLM1<-dimnames(summary.lm(LinearModel1(1))$coefficients)[1][[1]]

#Running the linear model:
LM1Betas<-matrix(0, length(normNoOutliers[,1]), NumberofXsLM1)
colnames(LM1Betas)<-NameofXsLM1
row.names(LM1Betas)<-row.names(normNoOutliers)

LM1pvalues<-matrix(0, length(normNoOutliers[,1]), NumberofXsLM1)
colnames(LM1pvalues)<-NameofXsLM1
row.names(LM1pvalues)<-row.names(normNoOutliers)

LM1SE<-matrix(0, length(normNoOutliers[,1]), NumberofXsLM1)
colnames(LM1SE)<-NameofXsLM1
row.names(LM1SE)<-row.names(normNoOutliers)

LM1Tstat<-matrix(0, length(normNoOutliers[,1]), NumberofXsLM1)
colnames(LM1Tstat)<-NameofXsLM1
row.names(LM1Tstat)<-row.names(normNoOutliers)

for(i in 1:length(normNoOutliers[,1])){
	TempLM1<-LinearModel1(i)
	LM1Betas[i,]<-summary.lm(TempLM1)$coefficients[,1]
	LM1SE[i,]<-summary.lm(TempLM1)$coefficients[,2]
	LM1Tstat[i,]<-summary.lm(TempLM1)$coefficients[,3]
	LM1pvalues[i,]<-summary.lm(TempLM1)$coefficients[,4]
}


#Outputting a histogram illustrating the raw p-values in comparison to what would be expected by chance:
for (i in 1:NumberofXsLM1){
png(paste(paste("17 Histogram of Unadjusted Pvalues Using LM1 for", NameofXsLM1[i], sep="  "), "png", sep="."))	
hist(LM1pvalues[,i], breaks=100, col=i, main=paste("Unadjusted P-values using LM1 for", NameofXsLM1[i], sep="  "), xlab="Unadjusted p-value", ylab="Count")
abline(a=(length(normNoOutliers[,1])/100), b=0)
dev.off()		
}		

#Outputting the raw pvalues and LM related statistics:
LM1BetasOutput<-cbind(LM1Betas, GeneIdentifiers_filtered) 
LM1SEOutput<-cbind(LM1Betas, GeneIdentifiers_filtered)
LM1TstatOutput<-cbind(LM1Tstat, GeneIdentifiers_filtered)
LM1pvaluesOutput<-cbind(LM1pvalues, GeneIdentifiers_filtered)

write.csv(LM1BetasOutput, "LM1Betas.csv")
write.csv(LM1SEOutput, "LM1SE.csv")
write.csv(LM1TstatOutput, "LM1Tstat.csv")
write.csv(LM1pvaluesOutput, "LM1pvaluesRAW.csv")


#Applying two different types of multiple-comparison corrections to the raw p-values (Benjamini-Hochberg and Benjamini-Yekutieli):
TempPvalAdjLM1<-mt.rawp2adjp(LM1pvalues[,2], proc=c("BH", "BY"))
LM1PvalAdj<-TempPvalAdjLM1$adjp[order(TempPvalAdjLM1$index),]
row.names(LM1PvalAdj)<-row.names(normNoOutliers)

#Outputting the p-values, betas, and gene identifiers:
LM1testOutput<-cbind(LM1PvalAdj, LM1Betas[,2], GeneIdentifiers_filtered)
row.names(LM1testOutput)<-row.names(normNoOutliers)
colnames(LM1testOutput)<-c("Unadjusted P-value", "BH Adjusted Pvalue", "BY Adjusted Pvalue", "Beta", "Symbol", "Search Key", "ILMN_Gene", "Chromosome", "Definition", "Synonyms", "ProbeID", "REFSEQ_ID", "UNIGENE_ID", "ENTREZ_GENE_ID", "PROBE_TYPE", "PROBE_START", "PROBE_SEQUENCE","ONTOLOGY_COMPONENT","ONTOLOGY_PROCESS", "ONTOLOGY_FUNCTION")


write.csv(LM1testOutput, "LM1testOutput.csv")

#***END OF CODE FOR STEP 17****


##4. The results are found in "LM1testOutput.csv."

#************************************************************************************************
##2. Insert those variables into this code, placing diagnosis as the first variable.

LinearModel2<-function(i){lm(normNoOutliers[i,]~DiagnosisNoOutliers+AgeNoOutliers+RNAintegrityNoOutliers+BrainpHNoOutliers+GenderNoOutliers+HoursFinalCorrectedNoOutliers)}

##3. Then run this code (Ctrl+R) to get the p-values:


#***CODE TO RUN FOR STEP 17****

#Determining how many variables are in the linear model (since LM breaks factors down into dummy variables)
NumberofXsLM2<-length(summary.lm(LinearModel2(1))$coefficients[,1])
NameofXsLM2<-dimnames(summary.lm(LinearModel2(1))$coefficients)[1][[1]]

#Running the linear model:
LM2Betas<-matrix(0, length(normNoOutliers[,1]), NumberofXsLM2)
colnames(LM2Betas)<-NameofXsLM2
row.names(LM2Betas)<-row.names(normNoOutliers)

LM2pvalues<-matrix(0, length(normNoOutliers[,1]), NumberofXsLM2)
colnames(LM2pvalues)<-NameofXsLM2
row.names(LM2pvalues)<-row.names(normNoOutliers)

LM2SE<-matrix(0, length(normNoOutliers[,1]), NumberofXsLM2)
colnames(LM2SE)<-NameofXsLM2
row.names(LM2SE)<-row.names(normNoOutliers)

LM2Tstat<-matrix(0, length(normNoOutliers[,1]), NumberofXsLM2)
colnames(LM2Tstat)<-NameofXsLM2
row.names(LM2Tstat)<-row.names(normNoOutliers)

for(i in 1:length(normNoOutliers[,1])){
	TempLM2<-LinearModel2(i)
	LM2Betas[i,]<-summary.lm(TempLM2)$coefficients[,1]
	LM2SE[i,]<-summary.lm(TempLM2)$coefficients[,2]
	LM2Tstat[i,]<-summary.lm(TempLM2)$coefficients[,3]
	LM2pvalues[i,]<-summary.lm(TempLM2)$coefficients[,4]
}


#Outputting a histogram illustrating the raw p-values in comparison to what would be expected by chance:
for (i in 1:NumberofXsLM2){
png(paste(paste("17 Histogram of Unadjusted Pvalues Using LM2 for", NameofXsLM2[i], sep="  "), "png", sep="."))	
hist(LM2pvalues[,i], breaks=100, col=i, main=paste("Unadjusted P-values using LM2 for", NameofXsLM2[i], sep="  "), xlab="Unadjusted p-value", ylab="Count")
abline(a=(length(normNoOutliers[,1])/100), b=0)
dev.off()		
}		

#Outputting the raw pvalues and LM related statistics:
LM2BetasOutput<-cbind(LM2Betas, GeneIdentifiers_filtered) 
LM2SEOutput<-cbind(LM2Betas, GeneIdentifiers_filtered)
LM2TstatOutput<-cbind(LM2Tstat, GeneIdentifiers_filtered)
LM2pvaluesOutput<-cbind(LM2pvalues, GeneIdentifiers_filtered)

write.csv(LM2BetasOutput, "LM2Betas.csv")
write.csv(LM2SEOutput, "LM2SE.csv")
write.csv(LM2TstatOutput, "LM2Tstat.csv")
write.csv(LM2pvaluesOutput, "LM2pvaluesRAW.csv")


#Output for MDD:

#Calculating Fold Change using the two methods most commonly used in our lab:
FC<-2^LM2Betas[,2]

FCv2<-FC

#NOTE: Update - this code is INCORRECT - IGNORE OUTPUT
for(i in 1:length(FC)){
if(FC[i]>1){FCv2[i]<-FC[i]}else{FCv2[i]<-(-1/FC[i])}
}

#Applying two different types of multiple-comparison corrections to the raw p-values (Benjamini-Hochberg and Benjamini-Yekutieli):
TempPvalAdjLM2<-mt.rawp2adjp(LM2pvalues[,2], proc=c("BH", "BY"))
LM2PvalAdj<-TempPvalAdjLM2$adjp[order(TempPvalAdjLM2$index),]
row.names(LM2PvalAdj)<-row.names(normNoOutliers)

#Outputting the p-values, betas, and gene identifiers:
LM2testOutputMDD<-cbind(row.names(normNoOutliers), LM2PvalAdj, LM2Betas[,2], FC, FCv2, NormNoOutliersMedianCTRL, NormNoOutliersMedianMDD, Detection_PvalueMedianCTRL, Detection_PvalueMedianMDD,GeneIdentifiers_filtered)
row.names(LM2testOutputMDD)<-row.names(normNoOutliers)
colnames(LM2testOutputMDD)[c(1:11)]<-c("X", "Unadjusted P-value", "BH Adjusted Pvalue", "BY Adjusted Pvalue", "Beta", "FoldChange", "FoldChangev2", "ControlMedianLog2NormSignal", "MDDMedianLog2NormSignal", "ControlMedianDetectionPval", "MDDMedianDetectionPval")

LM2testOutputMDDv2<- join(as.data.frame(LM2testOutputMDD), ProbeInfoOutput, by = "X")
write.csv(LM2testOutputMDDv2, "LM2testOutputMDDv2.csv")


#Output for Age:

#Applying two different types of multiple-comparison corrections to the raw p-values (Benjamini-Hochberg and Benjamini-Yekutieli):
TempPvalAdjLM2<-mt.rawp2adjp(LM2pvalues[,3], proc=c("BH", "BY"))
LM2PvalAdj<-TempPvalAdjLM2$adjp[order(TempPvalAdjLM2$index),]
row.names(LM2PvalAdj)<-row.names(normNoOutliers)

#Outputting the p-values, betas, and gene identifiers:
LM2testOutputAge<-cbind(row.names(normNoOutliers), LM2PvalAdj, LM2Betas[,3],NormNoOutliersMedianCTRL, Detection_PvalueMedianCTRL, GeneIdentifiers_filtered)
row.names(LM2testOutputAge)<-row.names(normNoOutliers)
colnames(LM2testOutputAge)[c(1:7)]<-c("X", "Unadjusted P-value", "BH Adjusted Pvalue", "BY Adjusted Pvalue", "Beta", "ControlMedianLog2NormSignal", "ControlMedianDetectionPval")

LM2testOutputAgev2<- join(as.data.frame(LM2testOutputAge), ProbeInfoOutput, by = "X")
write.csv(LM2testOutputAgev2, "LM2testOutputAgev2.csv")


#Output for RNA Integrity:

#Applying two different types of multiple-comparison corrections to the raw p-values (Benjamini-Hochberg and Benjamini-Yekutieli):
TempPvalAdjLM2<-mt.rawp2adjp(LM2pvalues[,4], proc=c("BH", "BY"))
LM2PvalAdj<-TempPvalAdjLM2$adjp[order(TempPvalAdjLM2$index),]
row.names(LM2PvalAdj)<-row.names(normNoOutliers)

#Outputting the p-values, betas, and gene identifiers:
LM2testOutputRNAIntegrity<-cbind(row.names(normNoOutliers), LM2PvalAdj, LM2Betas[,4], NormNoOutliersMedianCTRL, Detection_PvalueMedianCTRL, GeneIdentifiers_filtered)
row.names(LM2testOutputRNAIntegrity)<-row.names(normNoOutliers)
colnames(LM2testOutputRNAIntegrity)[c(1:7)]<-c("X", "Unadjusted P-value", "BH Adjusted Pvalue", "BY Adjusted Pvalue", "Beta", "ControlMedianLog2NormSignal", "ControlMedianDetectionPval")

LM2testOutputRNAIntegrityv2<- join(as.data.frame(LM2testOutputRNAIntegrity), ProbeInfoOutput, by = "X")
write.csv(LM2testOutputRNAIntegrityv2, "LM2testOutputRNAIntegrityv2.csv")


#Output for BrainPH:

#Applying two different types of multiple-comparison corrections to the raw p-values (Benjamini-Hochberg and Benjamini-Yekutieli):
TempPvalAdjLM2<-mt.rawp2adjp(LM2pvalues[,5], proc=c("BH", "BY"))
LM2PvalAdj<-TempPvalAdjLM2$adjp[order(TempPvalAdjLM2$index),]
row.names(LM2PvalAdj)<-row.names(normNoOutliers)

#Outputting the p-values, betas, and gene identifiers:
LM2testOutputBrainPH<-cbind(row.names(normNoOutliers), LM2PvalAdj, LM2Betas[,5], NormNoOutliersMedianCTRL, Detection_PvalueMedianCTRL, GeneIdentifiers_filtered)
row.names(LM2testOutputBrainPH)<-row.names(normNoOutliers)
colnames(LM2testOutputBrainPH)[c(1:7)]<-c("X", "Unadjusted P-value", "BH Adjusted Pvalue", "BY Adjusted Pvalue", "Beta", "ControlMedianLog2NormSignal", "ControlMedianDetectionPval")

LM2testOutputBrainPHv2<- join(as.data.frame(LM2testOutputBrainPH), ProbeInfoOutput, by = "X")
write.csv(LM2testOutputBrainPHv2, "LM2testOutputBrainPHv2.csv")


#Output for Gender:

#Applying two different types of multiple-comparison corrections to the raw p-values (Benjamini-Hochberg and Benjamini-Yekutieli):
TempPvalAdjLM2<-mt.rawp2adjp(LM2pvalues[,6], proc=c("BH", "BY"))
LM2PvalAdj<-TempPvalAdjLM2$adjp[order(TempPvalAdjLM2$index),]
row.names(LM2PvalAdj)<-row.names(normNoOutliers)

#Outputting the p-values, betas, and gene identifiers:
LM2testOutputGender<-cbind(row.names(normNoOutliers), LM2PvalAdj, LM2Betas[,6], NormNoOutliersMedianCTRL, Detection_PvalueMedianCTRL, GeneIdentifiers_filtered)
row.names(LM2testOutputGender)<-row.names(normNoOutliers)
colnames(LM2testOutputGender)[c(1:7)]<-c("X", "Unadjusted P-value", "BH Adjusted Pvalue", "BY Adjusted Pvalue", "Beta", "ControlMedianLog2NormSignal", "ControlMedianDetectionPval")

LM2testOutputGenderv2<- join(as.data.frame(LM2testOutputGender), ProbeInfoOutput, by = "X")
write.csv(LM2testOutputGenderv2, "LM2testOutputGenderv2.csv")



#Output for HoursFinal:

#Applying two different types of multiple-comparison corrections to the raw p-values (Benjamini-Hochberg and Benjamini-Yekutieli):
TempPvalAdjLM2<-mt.rawp2adjp(LM2pvalues[,7], proc=c("BH", "BY"))
LM2PvalAdj<-TempPvalAdjLM2$adjp[order(TempPvalAdjLM2$index),]
row.names(LM2PvalAdj)<-row.names(normNoOutliers)

#Outputting the p-values, betas, and gene identifiers:
LM2testOutputHoursFinal<-cbind(row.names(normNoOutliers), LM2PvalAdj, LM2Betas[,7], NormNoOutliersMedianCTRL, Detection_PvalueMedianCTRL, GeneIdentifiers_filtered)
row.names(LM2testOutputHoursFinal)<-row.names(normNoOutliers)
colnames(LM2testOutputHoursFinal)[c(1:7)]<-c("X", "Unadjusted P-value", "BH Adjusted Pvalue", "BY Adjusted Pvalue", "Beta", "ControlMedianLog2NormSignal", "ControlMedianDetectionPval")

LM2testOutputHoursFinalv2<- join(as.data.frame(LM2testOutputHoursFinal), ProbeInfoOutput, by = "X")
write.csv(LM2testOutputHoursFinalv2, "LM2testOutputHoursFinalv2.csv")


#Creating a normNoOutliers output that is corrected for each of the confounds in the dataset:

DiagnosisNumeric<-as.numeric(DiagnosisNoOutliers)-1
GenderNumeric<-as.numeric(GenderNoOutliers)-1

NormNoOutliersMedian<-apply(normNoOutliers, 1, median)
NormNoOutliersMean<-apply(normNoOutliers, 1, mean)

NormNoOutliersMeanCTRL<-apply(normNoOutliers[,DiagnosisNoOutliers=="CTRL"], 1, mean)

normNoOutliersCorrected<-matrix(0, length(normNoOutliers[,1]),length(normNoOutliers[1,]))

for(i in 1:length(normNoOutliers[,1])){
for(j in 1:length(normNoOutliers[1,])){

normNoOutliersCorrected[i,j]<-
(normNoOutliers[i,j]-(
LM2Betas[i,1]+
LM2Betas[i,3]*AgeNoOutliers[j]+
LM2Betas[i,4]*RNAintegrityNoOutliers[j]+
LM2Betas[i,5]*BrainpHNoOutliers[j]+
LM2Betas[i,6]*GenderNumeric[j]+
LM2Betas[i,7]*HoursFinalCorrectedNoOutliers[j]
)
+ NormNoOutliersMeanCTRL[i])

}
}



NormNoOutliersCorrectedMedian<-apply(normNoOutliersCorrected, 1, median)
NormNoOutliersCorrectedMean<-apply(normNoOutliersCorrected, 1, mean)
NormNoOutliersCorrectednoDiagnosisMean<-apply(normNoOutliersCorrectednoDiagnosis, 1, mean)

plot(NormNoOutliersCorrectedMedian~NormNoOutliersMedian)
plot(NormNoOutliersCorrectedMean~NormNoOutliersMean)
plot(NormNoOutliersCorrectednoDiagnosisMean~NormNoOutliersMean)


#This code cleans out *all* known effects:

normNoOutliersCorrectednoDiagnosis<-matrix(0, length(normNoOutliers[,1]),length(normNoOutliers[1,]))

for(i in 1:length(normNoOutliers[,1])){
for(j in 1:length(normNoOutliers[1,])){

normNoOutliersCorrectednoDiagnosis[i,j]<-
(normNoOutliers[i,j]-(
LM2Betas[i,1]+
LM2Betas[i,2]*DiagnosisNumeric[j]+
LM2Betas[i,3]*AgeNoOutliers[j]+
LM2Betas[i,4]*RNAintegrityNoOutliers[j]+
LM2Betas[i,5]*BrainpHNoOutliers[j]+
LM2Betas[i,6]*GenderNumeric[j]+
LM2Betas[i,7]*HoursFinalCorrectedNoOutliers[j]
)
+ NormNoOutliersMean[i])

}
}



#This code cleans out *all* known effects... and then sets the mean gene expression at what would be expected for 
#mean(Age)=51.92857
#Gender=Male
#Diagnosis=Control
#mean(RNAintegrity)+sd(RNAintegrity)=5.035714+0.9330834= 5.968798
#mean(BrainpH)+sd(BrainpH)=6.809048+0.156732= 6.96578
#mean(HoursFinalCorrected)-sd(HoursFinalCorrected)=21.40833-8.007861=13.40047


AgeNoOutliersCentered<-AgeNoOutliers-52
plot(AgeNoOutliersCentered~AgeNoOutliers)

RNAintegrityNoOutliersCentered<-RNAintegrityNoOutliers-6
plot(RNAintegrityNoOutliersCentered~RNAintegrityNoOutliers)

BrainpHNoOutliersCentered<-BrainpHNoOutliers-7
plot(BrainpHNoOutliersCentered~BrainpHNoOutliers)

HoursFinalCorrectedNoOutliersCentered<-HoursFinalCorrectedNoOutliers-13.4
plot(HoursFinalCorrectedNoOutliersCentered~HoursFinalCorrectedNoOutliers)


LinearModel4<-function(i){lm(normNoOutliers[i,]~DiagnosisNumeric+AgeNoOutliersCentered+RNAintegrityNoOutliersCentered+BrainpHNoOutliersCentered+GenderNumeric+HoursFinalCorrectedNoOutliersCentered)}


##3. Then run this code (Ctrl+R) to get the p-values:


#***CODE TO RUN FOR STEP 17****

#Determining how many variables are in the linear model (since LM breaks factors down into dummy variables)
NumberofXsLM4<-length(summary.lm(LinearModel4(1))$coefficients[,1])
NameofXsLM4<-dimnames(summary.lm(LinearModel4(1))$coefficients)[1][[1]]

#Running the linear model:
LM4Betas<-matrix(0, length(normNoOutliers[,1]), NumberofXsLM4)
colnames(LM4Betas)<-NameofXsLM4
row.names(LM4Betas)<-row.names(normNoOutliers)

LM4pvalues<-matrix(0, length(normNoOutliers[,1]), NumberofXsLM4)
colnames(LM4pvalues)<-NameofXsLM4
row.names(LM4pvalues)<-row.names(normNoOutliers)

LM4SE<-matrix(0, length(normNoOutliers[,1]), NumberofXsLM4)
colnames(LM4SE)<-NameofXsLM4
row.names(LM4SE)<-row.names(normNoOutliers)

LM4Tstat<-matrix(0, length(normNoOutliers[,1]), NumberofXsLM4)
colnames(LM4Tstat)<-NameofXsLM4
row.names(LM4Tstat)<-row.names(normNoOutliers)

LM4Residuals<-matrix(0, length(normNoOutliers[,1]), length(normNoOutliers[1,]))
colnames(LM4Residuals)<-colnames(normNoOutliers)
row.names(LM4Residuals)<-row.names(normNoOutliers)


for(i in 1:length(normNoOutliers[,1])){
	TempLM4<-LinearModel4(i)
	LM4Betas[i,]<-summary.lm(TempLM4)$coefficients[,1]
	LM4SE[i,]<-summary.lm(TempLM4)$coefficients[,2]
	LM4Tstat[i,]<-summary.lm(TempLM4)$coefficients[,3]
	LM4pvalues[i,]<-summary.lm(TempLM4)$coefficients[,4]
	LM4Residuals[i,]<-TempLM4$residuals+summary.lm(TempLM4)$coefficients[1,1]
}

#I think *hypothetically* the intercept (+/- SE) should represent the expected value for the probe in 52-year old control males with pH 7, PMI 13, RNAIntegrity of 6
#Although the centering of the variables away from the actual median values will mean that the SE for the Intercept may be optimistic

plot(LM2Betas[,1]~LM4Betas[,1])

write.csv(LM4Residuals, "LM4ResidualsCleanedSemiOptimal.csv")

LM4ResidualsRank<-LM4Residuals

for(i in 1:length(LM4Residuals[1,])){
LM4ResidualsRank[,i]<-round((rank(LM4Residuals[,i]))/(length(LM4Residuals[,1])), digits=6)
}

write.csv(LM4ResidualsRank, "LM4ResidualsRank.csv")

LM4rankIntercept<- round((rank(LM4Betas[,1]))/(length(LM4Betas[,1])), digits=6)


LM4InterceptOutput<-cbind(row.names(LM4Residuals), LM4Betas[,1], LM4SE[,1], LM4rankIntercept, GeneIdentifiers_filtered)
colnames(LM4InterceptOutput)[c(1:3)]<-c("X", "Intercept", "InterceptSE")

LM4InterceptOutputv2<- join(as.data.frame(LM4InterceptOutput), ProbeInfoOutput, by = "X")
write.csv(LM4InterceptOutputv2, "LM4InterceptOutputv2.csv")




#This code cleans out *all* known effects... and then sets the mean gene expression at what would be expected for 
#median(Age)=50
#Gender=Male
#Diagnosis=Control
#median(RNAintegrity)=5
#median(BrainpH)=6.805
#median(HoursFinalCorrected)=22.625


AgeNoOutliersCentered2<-AgeNoOutliers-50
plot(AgeNoOutliersCentered~AgeNoOutliers)

RNAintegrityNoOutliersCentered2<-RNAintegrityNoOutliers-5
plot(RNAintegrityNoOutliersCentered2~RNAintegrityNoOutliers)

BrainpHNoOutliersCentered2<-BrainpHNoOutliers-6.8
plot(BrainpHNoOutliersCentered2~BrainpHNoOutliers)

HoursFinalCorrectedNoOutliersCentered2<-HoursFinalCorrectedNoOutliers-22.6
plot(HoursFinalCorrectedNoOutliersCentered2~HoursFinalCorrectedNoOutliers)


LinearModel5<-function(i){lm(normNoOutliers[i,]~DiagnosisNumeric+AgeNoOutliersCentered2+RNAintegrityNoOutliersCentered2+BrainpHNoOutliersCentered2+GenderNumeric+HoursFinalCorrectedNoOutliersCentered2)}


##3. Then run this code (Ctrl+R) to get the p-values:


#***CODE TO RUN FOR STEP 17****

#Determining how many variables are in the linear model (since LM breaks factors down into dummy variables)
NumberofXsLM5<-length(summary.lm(LinearModel5(1))$coefficients[,1])
NameofXsLM5<-dimnames(summary.lm(LinearModel5(1))$coefficients)[1][[1]]

#Running the linear model:
LM5Betas<-matrix(0, length(normNoOutliers[,1]), NumberofXsLM5)
colnames(LM5Betas)<-NameofXsLM5
row.names(LM5Betas)<-row.names(normNoOutliers)

LM5pvalues<-matrix(0, length(normNoOutliers[,1]), NumberofXsLM5)
colnames(LM5pvalues)<-NameofXsLM5
row.names(LM5pvalues)<-row.names(normNoOutliers)

LM5SE<-matrix(0, length(normNoOutliers[,1]), NumberofXsLM5)
colnames(LM5SE)<-NameofXsLM5
row.names(LM5SE)<-row.names(normNoOutliers)

LM5Tstat<-matrix(0, length(normNoOutliers[,1]), NumberofXsLM5)
colnames(LM5Tstat)<-NameofXsLM5
row.names(LM5Tstat)<-row.names(normNoOutliers)

LM5Residuals<-matrix(0, length(normNoOutliers[,1]), length(normNoOutliers[1,]))

for(i in 1:length(normNoOutliers[,1])){
	TempLM5<-LinearModel5(i)
	LM5Betas[i,]<-summary.lm(TempLM5)$coefficients[,1]
	LM5SE[i,]<-summary.lm(TempLM5)$coefficients[,2]
	LM5Tstat[i,]<-summary.lm(TempLM5)$coefficients[,3]
	LM5pvalues[i,]<-summary.lm(TempLM5)$coefficients[,4]
	LM5Residuals[i,]<-TempLM5$residuals+summary.lm(TempLM5)$coefficients[1,1]
}

#I think *hypothetically* the intercept (+/- SE) should represent the expected value for the probe in 52-year old control males with pH 7, PMI 13, RNAIntegrity of 6
#Although the centering of the variables away from the actual median values will mean that the SE for the Intercept may be optimistic

plot(LM2Betas[,1]~LM5Betas[,1])


write.csv(LM5Residuals, "LM5ResidualsCleanedMedian.csv")


LM5ResidualsRank<-LM5Residuals

for(i in 1:length(LM5Residuals[1,])){
LM5ResidualsRank[,i]<-round((rank(LM5Residuals[,i]))/(length(LM5Residuals[,1])), digits=6)
}

write.csv(LM5ResidualsRank, "LM5ResidualsRank.csv")

LM5rankIntercept<- round((rank(LM5Betas[,1]))/(length(LM5Betas[,1])), digits=6)


LM5InterceptOutput<-cbind(row.names(LM5Residuals), LM5Betas[,1], LM5SE[,1], LM5rankIntercept, GeneIdentifiers_filtered)
colnames(LM5InterceptOutput)[c(1:3)]<-c("X", "Intercept", "InterceptSE")

LM5InterceptOutputv2<- join(as.data.frame(LM5InterceptOutput), ProbeInfoOutput, by = "X")
write.csv(LM5InterceptOutputv2, "LM5InterceptOutputv2.csv")











PvalueOutputSummary<-matrix(0, 6 , 5)
colnames(PvalueOutputSummary)<-c("p-val<0.001", "Random: p-val<0.001", "p-val<0.01", "Random: p-val<0.01", "BH p-val<0.10")
row.names(PvalueOutputSummary)<-c("MDD", "Age", "RNAIntegrity", "BrainPH", "Gender", "HoursFinal")

for(i in 2:7){
PvalueOutputSummary[(i-1), 1]<-sum(LM2pvalues[,i]<0.001)
PvalueOutputSummary[(i-1), 2]<-round(length(LM2pvalues[,i])*0.001, digits=0)
PvalueOutputSummary[(i-1), 3]<-sum(LM2pvalues[,i]<0.01)
PvalueOutputSummary[(i-1), 4]<-round(length(LM2pvalues[,i])*0.01, digits=0)
}

PvalueOutputSummary[1, 5]<-sum(LM2testOutputMDD[,3]<0.10)
PvalueOutputSummary[2, 5]<-sum(LM2testOutputAge[,3]<0.10)
PvalueOutputSummary[3, 5]<-sum(LM2testOutputRNAIntegrity[,3]<0.10)
PvalueOutputSummary[4, 5]<-sum(LM2testOutputBrainPH[,3]<0.10)
PvalueOutputSummary[5, 5]<-sum(LM2testOutputGender[,3]<0.10)
PvalueOutputSummary[6, 5]<-sum(LM2testOutputHoursFinal[,3]<0.10)

write.csv(PvalueOutputSummary, "PvalueOutputSummary.csv")



#***END OF CODE FOR STEP 17****


#*****To be prudent, I also included a model that includes more traditional confounds*************

##2. Insert those variables into this code, placing diagnosis as the first variable.

LinearModel2<-function(i){lm(normNoOutliers[i,]~DiagnosisNoOutliers+AgeNoOutliers+RNAintegrityNoOutliers+GenderNoOutliers+HoursFinalNoOutliers+BrainpHNoOutliers+RNAconcNoOutliers+RNAAgeNoOutliers)}


##3. Then run this code (Ctrl+R) to get the p-values:



#***CODE TO RUN FOR STEP 17****

#Determining how many variables are in the linear model (since LM breaks factors down into dummy variables)
NumberofXsLM2<-length(summary.lm(LinearModel2(1))$coefficients[,1])
NameofXsLM2<-dimnames(summary.lm(LinearModel2(1))$coefficients)[1][[1]]

#Running the linear model:
LM2Betas<-matrix(0, length(normNoOutliers[,1]), NumberofXsLM2)
colnames(LM2Betas)<-NameofXsLM2
row.names(LM2Betas)<-row.names(normNoOutliers)

LM2pvalues<-matrix(0, length(normNoOutliers[,1]), NumberofXsLM2)
colnames(LM2pvalues)<-NameofXsLM2
row.names(LM2pvalues)<-row.names(normNoOutliers)

LM2SE<-matrix(0, length(normNoOutliers[,1]), NumberofXsLM2)
colnames(LM2SE)<-NameofXsLM2
row.names(LM2SE)<-row.names(normNoOutliers)

LM2Tstat<-matrix(0, length(normNoOutliers[,1]), NumberofXsLM2)
colnames(LM2Tstat)<-NameofXsLM2
row.names(LM2Tstat)<-row.names(normNoOutliers)

for(i in 1:length(normNoOutliers[,1])){
	TempLM2<-LinearModel2(i)
	LM2Betas[i,]<-summary.lm(TempLM2)$coefficients[,1]
	LM2SE[i,]<-summary.lm(TempLM2)$coefficients[,2]
	LM2Tstat[i,]<-summary.lm(TempLM2)$coefficients[,3]
	LM2pvalues[i,]<-summary.lm(TempLM2)$coefficients[,4]
}


#Outputting a histogram illustrating the raw p-values in comparison to what would be expected by chance:
for (i in 1:NumberofXsLM2){
png(paste(paste("17 Histogram of Unadjusted Pvalues Using LM2 for", NameofXsLM2[i], sep="  "), "png", sep="."))	
hist(LM2pvalues[,i], breaks=100, col=i, main=paste("Unadjusted P-values using LM2 for", NameofXsLM2[i], sep="  "), xlab="Unadjusted p-value", ylab="Count")
abline(a=(length(normNoOutliers[,1])/100), b=0)
dev.off()		
}		

#Outputting the raw pvalues and LM related statistics:
LM2BetasOutput<-cbind(LM2Betas, GeneIdentifiers_filtered) 
LM2SEOutput<-cbind(LM2Betas, GeneIdentifiers_filtered)
LM2TstatOutput<-cbind(LM2Tstat, GeneIdentifiers_filtered)
LM2pvaluesOutput<-cbind(LM2pvalues, GeneIdentifiers_filtered)

write.csv(LM2BetasOutput, "LM2Betas.csv")
write.csv(LM2SEOutput, "LM2SE.csv")
write.csv(LM2TstatOutput, "LM2Tstat.csv")
write.csv(LM2pvaluesOutput, "LM2pvaluesRAW.csv")


#Applying two different types of multiple-comparison corrections to the raw p-values (Benjamini-Hochberg and Benjamini-Yekutieli):
TempPvalAdjLM2<-mt.rawp2adjp(LM2pvalues[,2], proc=c("BH", "BY"))
LM2PvalAdj<-TempPvalAdjLM2$adjp[order(TempPvalAdjLM2$index),]
row.names(LM2PvalAdj)<-row.names(normNoOutliers)

#Outputting the p-values, betas, and gene identifiers:
LM2testOutput<-cbind(LM2PvalAdj, LM2Betas[,2], GeneIdentifiers_filtered)
row.names(LM2testOutput)<-row.names(normNoOutliers)
colnames(LM2testOutput)<-c("Unadjusted P-value", "BH Adjusted Pvalue", "BY Adjusted Pvalue", "Beta", "Symbol", "Search Key", "ILMN_Gene", "Chromosome", "Definition", "Synonyms", "ProbeID", "REFSEQ_ID", "UNIGENE_ID", "ENTREZ_GENE_ID", "PROBE_TYPE", "PROBE_START", "PROBE_SEQUENCE","ONTOLOGY_COMPONENT","ONTOLOGY_PROCESS", "ONTOLOGY_FUNCTION")


write.csv(LM2testOutput, "LM2testOutput.csv")

#***END OF CODE FOR STEP 17****


##4. The results are found in "LM2testOutput.csv."


#************************************************************************************************

#STEP 18:  Double checking the data for expected effects:

##Gender:

png("18 Gender Check.png")
plot(normNoOutliers[GeneIdentifiers_filtered[,1]=="XIST",]~normNoOutliers[GeneIdentifiers_filtered[,1]=="RPS4Y1",], xlab="Male Gene: RPS4Y1", ylab="Female Gene: XIST", main = "Does subject gender match gendered gene expression?")
points(normNoOutliers[GeneIdentifiers_filtered[,1]=="XIST", GenderNoOutliers=="M"]~normNoOutliers[GeneIdentifiers_filtered[,1]=="RPS4Y1", GenderNoOutliers=="M"], col=2)
points(normNoOutliers[GeneIdentifiers_filtered[,1]=="XIST", GenderNoOutliers=="F"]~normNoOutliers[GeneIdentifiers_filtered[,1]=="RPS4Y1", GenderNoOutliers=="F"], col=3)
legend(min(normNoOutliers[GeneIdentifiers_filtered[,1]=="RPS4Y1",])+0.5, min(normNoOutliers[GeneIdentifiers_filtered[,1]=="XIST",])+0.5, c("M", "F"), text.col=c(2,3), pch=19, col=c(2,3))
dev.off()

GenderCheck<-rbind(normNoOutliers[GeneIdentifiers_filtered[,1]=="XIST",],normNoOutliers[GeneIdentifiers_filtered[,1]=="RPS4Y1",])

write.csv(GenderCheck, "18 GenderCheck.csv")

##PH: Using genes identified as pH-related in the Affy Amygdala data set:


PHGenesNoOutliers<-rbind(normNoOutliers[GeneIdentifiers_filtered[,1]=="NDUFA8",], normNoOutliers[GeneIdentifiers_filtered[,1]=="SFRP1",], normNoOutliers[GeneIdentifiers_filtered[,1]=="PTS",], normNoOutliers[GeneIdentifiers_filtered[,1]=="SFRP1",], normNoOutliers[GeneIdentifiers_filtered[,1]=="PSMD4",])

PredictedPHNoOutliers<-apply(PHGenesNoOutliers, 2, mean)

png("18 PH Check_Predicted vs Real.png")
plot(PredictedPHNoOutliers~BrainpHNoOutliers, xlab="Brain pH", ylab="Predicted pH", main="Predicting pH using gene expression: NDUFA8, SFRP1,PTS, SFRP1,PSMD4" )
dev.off()

png("18 PH Check_Predicted vs PC1.png")
plot(PredictedPHNoOutliers~PC1noOutliers, xlab="PC1", ylab="Predicted pH", main="Predicting pH using gene expression: NDUFA8, SFRP1,PTS, SFRP1,PSMD4" )
dev.off()

png("18 PH Check_Predicted vs PC2.png")
plot(PredictedPHNoOutliers~PC2noOutliers, xlab="PC2", ylab="Predicted pH", main="Predicting pH using gene expression: NDUFA8, SFRP1,PTS, SFRP1,PSMD4" )
dev.off()

PC1minusPC2noOutliers<-PC1noOutliers-PC2noOutliers

png("18 PH Check_Predicted vs PC1 and PC2.png")
plot(PredictedPHNoOutliers~PC1minusPC2noOutliers, xlab="PC1-PC2", ylab="Predicted pH", main="Predicting pH using gene expression: NDUFA8, SFRP1,PTS, SFRP1,PSMD4" )
dev.off()


##PH: Using genes identified as pH-related in the Affy data set, averaged across 4 brain regions (HC, DLPFC, ACG, NACC):

PHGenesNoOutliersPositive<-rbind(normNoOutliers[GeneIdentifiers_filtered[,1]=="ATP5O",], normNoOutliers[GeneIdentifiers_filtered[,1]=="NIT2",], normNoOutliers[GeneIdentifiers_filtered[,1]=="LRPPRC",], normNoOutliers[GeneIdentifiers_filtered[,1]=="USO1",])

PHGenesNoOutliersNegative<-rbind(normNoOutliers[GeneIdentifiers_filtered[,1]=="TUBB2B",], normNoOutliers[GeneIdentifiers_filtered[,1]=="ATP1B2",])

PredictedPHNoOutliers2<-apply(PHGenesNoOutliersPositive, 2, mean)-apply(PHGenesNoOutliersNegative, 2, mean)


png("18 PH Check_Predicted2 vs Real.png")
plot(PredictedPHNoOutliers2~BrainpHNoOutliers, xlab="Brain pH", ylab="Predicted pH", main="Predicting pH using gene expression: ATP5O, NIT2,LRPPRC, USO1,TUBB2B, ATP1B2" )
dev.off()

png("18 PH Check_Predicted2 vs PC1.png")
plot(PredictedPHNoOutliers2~PC1noOutliers, xlab="PC1", ylab="Predicted pH", main="Predicting pH using gene expression: ATP5O, NIT2,LRPPRC, USO1,TUBB2B, ATP1B2" )
dev.off()

png("18 PH Check_Predicted2 vs PC2.png")
plot(PredictedPHNoOutliers2~PC2noOutliers, xlab="PC2", ylab="Predicted pH", main="Predicting pH using gene expression: ATP5O, NIT2,LRPPRC, USO1,TUBB2B, ATP1B2" )
dev.off()


png("18 PH Check_Predicted2 vs PC1 and PC2.png")
plot(PredictedPHNoOutliers2~PC1minusPC2noOutliers, xlab="PC1-PC2", ylab="Predicted pH", main="Predicting pH using gene expression: ATP5O, NIT2,LRPPRC, USO1,TUBB2B, ATP1B2" )
dev.off()




FINALLY:  REMEMBER TO CLICK ON THE CONSOLE AND SAVE YOUR R WORKSPACE BEFORE EXITING!!!!!





########################

#SCRAP:


At the moment I am skipping doing qqplots:

Tstat<-mt.teststat(NormFiltered, (as.numeric(Diagnosis)-1))
names(Tstat)<-row.names(NormFiltered)

Tdf<-mt.teststat.num.denum(NormFiltered, (as.numeric(Diagnosis)-1))

tmp<-cbind(rt(18913, df=40), rt(18913, df=40), rt(18913, df=40), rt(18913, df=40), rt(18913, df=40))
RandomTdistribution<-apply(tmp, 1, mean)
length(RandomTdistribution)

png("qqplot for MDD vs Control Lateral Amygdala.png")
qqplot(RandomTdistribution, Tstat, xlab="Expected T", ylab="Observed T")
abline(a=0, b=1, col=2)
dev.off()


##Outputting files combining p-values with basic stats about the probes:

Detection_Pvalue_filtered<-Detection_Pvalue[detect05.gene>8,]

Detection_Pvalue_filteredNoOutliers<-Detection_Pvalue_filtered

Detection_Pvalue_filteredNoOutliers<-Detection_Pvalue_filtered[, OutlierRemovalIndices]

Detection_PvalueMedianMDD<-apply(Detection_Pvalue_filteredNoOutliers[,DiagnosisNoOutliers=="MDD"], 1, median)
Detection_PvalueMedianCTRL<-apply(Detection_Pvalue_filteredNoOutliers[,DiagnosisNoOutliers=="CTRL"], 1, median)
plot(log(Detection_PvalueMedianMDD)~log(Detection_PvalueMedianCTRL))

NormNoOutliersMedianMDD<-apply(normNoOutliers[,DiagnosisNoOutliers=="MDD"], 1, median)
NormNoOutliersMedianCTRL<-apply(normNoOutliers[,DiagnosisNoOutliers=="CTRL"], 1, median)

plot(NormNoOutliersMedianMDD~NormNoOutliersMedianCTRL)


#Outputting the p-values, betas, and gene identifiers:
LM1testOutput<-cbind(LM1PvalAdj, LM1Betas[,2], NormNoOutliersMedianCTRL, Detection_PvalueMedianCTRL, NormNoOutliersMedianMDD, Detection_PvalueMedianMDD, GeneIdentifiers_filtered)
row.names(LM1testOutput)<-row.names(normNoOutliers)
colnames(LM1testOutput)<-c("Unadjusted P-value", "BH Adjusted Pvalue", "BY Adjusted Pvalue", "Beta", "NormNoOutliersMedianCTRL", "Detection_PvalueMedianCTRL", "NormNoOutliersMedianMDD", "Detection_PvalueMedianMDD", "Symbol", "Search Key", "ILMN_Gene", "Chromosome", "Definition", "Synonyms", "ProbeID", "REFSEQ_ID", "UNIGENE_ID", "ENTREZ_GENE_ID", "PROBE_TYPE", "PROBE_START", "PROBE_SEQUENCE","ONTOLOGY_COMPONENT","ONTOLOGY_PROCESS", "ONTOLOGY_FUNCTION")

write.csv(LM1testOutput, "LM1testOutput.csv")


#Outputting the p-values, betas, and gene identifiers:
LM2testOutput<-cbind(LM2PvalAdj, LM2Betas[,2], NormNoOutliersMedianCTRL, Detection_PvalueMedianCTRL, NormNoOutliersMedianMDD, Detection_PvalueMedianMDD, GeneIdentifiers_filtered)
row.names(LM2testOutput)<-row.names(normNoOutliers)
colnames(LM2testOutput)<-c("Unadjusted P-value", "BH Adjusted Pvalue", "BY Adjusted Pvalue", "Beta", "NormNoOutliersMedianCTRL", "Detection_PvalueMedianCTRL", "NormNoOutliersMedianMDD", "Detection_PvalueMedianMDD", "Symbol", "Search Key", "ILMN_Gene", "Chromosome", "Definition", "Synonyms", "ProbeID", "REFSEQ_ID", "UNIGENE_ID", "ENTREZ_GENE_ID", "PROBE_TYPE", "PROBE_START", "PROBE_SEQUENCE","ONTOLOGY_COMPONENT","ONTOLOGY_PROCESS", "ONTOLOGY_FUNCTION")

write.csv(LM2testOutput, "LM2testOutput.csv")


####################################


library("BrainInABlender")

setwd("~/Documents/AMY LCM 10nuclei/Medial/13 BrainInABlender")

str(GeneIdentifiers_filtered)
head(GeneIdentifiers_filtered)
colnames(GeneIdentifiers_filtered)

sum(LM2testOutputMDDv2$X==row.names(GeneIdentifiers_filtered))
#So its in the same order

OriginalSymbolVsReannotated<-(as.character(LM2testOutputMDDv2$SYMBOLREANNOTATED)==GeneIdentifiers_filtered[,1])
table(OriginalSymbolVsReannotated)


#sum(is.na(GeneIdentifiers_filtered[,1]))
#[0]

sum(is.na(as.character(LM2testOutputMDDv2$SYMBOLREANNOTATED)))
head(normNoOutliers)

#temp<-data.frame(as.character(GeneIdentifiers_filtered[,1]), normNoOutliers, stringsAsFactors=F)
temp<-data.frame(as.character(LM2testOutputMDDv2$SYMBOLREANNOTATED[is.na(as.character(LM2testOutputMDDv2$SYMBOLREANNOTATED))==F]), normNoOutliers[is.na(as.character(LM2testOutputMDDv2$SYMBOLREANNOTATED))==F,], stringsAsFactors=F)
#Note - I changed this code to get rid of NAs, but double-checked and it turns out it doesn't matter.

write.csv(temp, "UserInput.csv")
str(temp)
#table(table(GeneIdentifiers_filtered[,1]))
#sort(table(GeneIdentifiers_filtered[,1]),decreasing=TRUE)[1:3]
#Ah - there are a whole bunch of blank gene ids.  That's o.k.

SirUnMixALotOutput<-Sir_UnMixALot(userInput=temp, dataColumns=c(2:40), geneColumn=1, species="human")

PublicationSpecific_CellTypeIndex<-SirUnMixALotOutput$PublicationSpecific_CellTypeIndex
AveragePrimary_CellTypeIndex<-SirUnMixALotOutput$AveragePrimary_CellTypeIndex

str(PublicationSpecific_CellTypeIndex)

#interesting - many of the support cell categories are not well differentiated (endothelial, mural, astrocyte, microglia) but all cluster together. All of the mature oligodendrocytes and all of the neurons are in their own cluster.


temp<-cbind(as.matrix(SubjectContinuousVariables), t(PublicationSpecific_CellTypeIndex), t(AveragePrimary_CellTypeIndex), PC1noOutliers, PC2noOutliers, PC3noOutliers, PC4noOutliers) 
str(temp)
SubjectContinuousVariables<-temp

#Jump to line 790 and re-run



row.names(AveragePrimary_CellTypeIndex)
for(i in c(1:10)){
  print(" 
        ")
  print(row.names(AveragePrimary_CellTypeIndex)[i])
  print(summary.lm(lm(AveragePrimary_CellTypeIndex[i,]~DiagnosisNoOutliers+AgeNoOutliers+RNAintegrityNoOutliers+BrainpHNoOutliers+GenderNoOutliers+HoursFinalCorrectedNoOutliers)))
}


#I need to re-run the PCA comparisons too.


write.csv(cor(cbind(as.numeric(SubjectFactorVariables[,2]), SubjectContinuousVariables), SubjectPCA), "Cor_SubjVar_vs_PCA.csv")


#Making prettier plots for the paper
row.names(AveragePrimary_CellTypeIndex)

pdf("MuralvsEndothelial_ByDiagnosis.pdf", height=4.5, width=4)
plot(AveragePrimary_CellTypeIndex[2,]~AveragePrimary_CellTypeIndex[4,], xlab="Mural Index", ylab="Endothelial Index", lwd=1.5, cex.lab=1.3, col=as.numeric(DiagnosisNoOutliers)*-1+4)
mtext("Medial", line=1, cex=1.5)
#legend(min(PC1noOutliers), max(PC2noOutliers)+10, c("Control", "MDD"), text.col=c(3, 2), pch=19, col=c(3, 2))
dev.off()


pdf("Endothelial_ByDiagnosis.pdf", height=4.5, width=4)
boxplot(AveragePrimary_CellTypeIndex[2,]~DiagnosisNoOutliers, xlab="Diagnosis", ylab="Endothelial Index", lwd=1.5, cex.lab=1.3)
stripchart(AveragePrimary_CellTypeIndex[2,]~DiagnosisNoOutliers, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 20, cex=1.5, col = c(3,2))
mtext("Medial", line=1, cex=1.5)
dev.off()


pdf("Mural_ByDiagnosis.pdf", height=4.5, width=4)
boxplot(AveragePrimary_CellTypeIndex[4,]~DiagnosisNoOutliers, xlab="Diagnosis", ylab="Mural Index", lwd=1.5, cex.lab=1.3)
stripchart(AveragePrimary_CellTypeIndex[4,]~DiagnosisNoOutliers, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 20, cex=1.5, col = c(3,2))
mtext("Medial", line=1, cex=1.5)
dev.off()


#############################

#Adding a volcano plot
#http://www.gettinggeneticsdone.com/2014/05/r-volcano-plots-to-visualize-rnaseq-microarray.html

#res <- read.table("results.txt", header=TRUE)

str(LM2Betas)
str(GeneIdentifiers_filtered)
str(LM2pvalues)
str(LM2PvalAdj)#This looks like it was overwritten at some point - wrong variables
str(LM2testOutputMDDv2) #Everything is stored as a factor. Goddamn my stupid early coding.
hist(as.numeric(as.character(LM2testOutputMDDv2$`BH Adjusted Pvalue`)))
padj<-as.numeric(as.character(LM2testOutputMDDv2$`BH Adjusted Pvalue`))
res<-data.frame(GeneIdentifiers_filtered[,1], LM2Betas[,2], LM2pvalues[,2], padj, stringsAsFactors = F) 
colnames(res)<-c("gene", "log2FoldChange", "pvalue", "padj")

head(res)

setwd("~/Documents/AMY LCM 10nuclei/Medial/09 LM2 Output_Best")

pdf("VolcanoPlot.pdf", height=4, width=4)
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot: Medial", xlim=c(-3,3), cex.lab=1.3, cex=0.6))

# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", cex=0.6))
with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", cex=0.6))
with(subset(res, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="green", cex=0.6))

# Label points with the textxy function from the calibrate plot
#I skipped this because there are too ficken many of them.
#library(calibrate)
#with(subset(res, padj<.05 & abs(log2FoldChange)>1), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=.8))
dev.off()

#Double checking the distribution of SE's for Diagnosis:

str(LM2SE)
colnames(LM2SE)

pdf("Histogram_SE.pdf", height=4, width=4)
hist(LM2SE[,2], breaks=100, col=2, xlab="SE", main="Medial", cex.lab=1.3)
dev.off()

#Maybe it would be better to read in all of the SE's for LM2 from each region and then make a violin plot illustrating them.
#Note: apparently the betas weren't outputted properly. I'm going to output it again now.
write.csv(LM2SE, "LM2SE.csv")



###################################################
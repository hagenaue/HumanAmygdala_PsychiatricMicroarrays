#Comparing qRT_PCR to Microarray for Vikram's Amygdala data
#Megan Hagenauer, April 2, 2018

#Brief methodological notes:

#For the qRT_PCR data, I have three values. 
#The first is the -DeltaDeltaCT or Log(2)FoldChange. To calculate this value, for each sample I subtracted the average of two PCR replicates for a housekeeping gene (AVE_Housekeeping_CT) from the average of two PCR replicates for the gene of interest (AVE_Gene_CT) to calculate the deltaCT.  Then I subtracted the AVE_deltaCT from the 10 MDD subjects from the AVE_deltaCT for the 10 CTRL subjects to calculate the -DeltaDeltaCT.  
#To calculate linear Fold Change, I calculated 2^Log(2)FoldChange.
#Finally, because the housekeeping gene sometimes appeared to potentially also show group differences, I did one version of the calculation that only compared the AVE_Gene_CT in the MDD and CTRL subjects (without first subtracting out the values for the reference housekeeping gene)

#For the microarray data, I have two values:
# Microarray_AveLogFC or Beta.  This is the beta from a larger linear model that controlled for other traditional subject variables (pH, PMI, age, gender, RNA Integrity). If there were two different probes representing the same gene, I averaged the betas from each of them.
# Microarray_AveFC(2^Beta) - in case someone wanted the fold change on a linear scale.


setwd("~/Documents/AMY LCM 10nuclei/10 Regional Comparison/MDD_RegionCompare/qRTPCR")

Data<-read.csv("qRT_PCR_forR.csv", header=T, stringsAsFactors = F)
str(Data)
# 'data.frame':	28 obs. of  7 variables:
#   $ Region                            : chr  "AAA" "AB" "AB" "AB" ...
# $ Gene                              : chr  "CTGF" "AOC2" "ARG2" "CTGF" ...
# $ logFC..DeltaDeltaCT.              : num  0.0468 1.6056 1.7515 1.0781 1.4509 ...
# $ FC                                : num  1.03 3.04 3.37 2.11 2.73 ...
# $ logFC..DeltaCT._noHousekeepingCTRL: num  -0.112 1.346 1.487 0.921 1.345 ...
# $ Microarray_AVELogFC.Beta.         : num  0.36 0.279 0.775 0.834 0.484 ...
# $ Microarray_AVEFC.2.Beta.          : num  1.28 1.21 1.71 1.78 1.4 ...


pdf("PCR_vs_microarray_logFC.pdf", height=5.5, width=5)
plot(Data$logFC..DeltaDeltaCT.~Data$Microarray_AVELogFC.Beta., ylab="Log(2) Fold Change (qRT-PCR)", xlab="Log(2) Fold Change (Microarray)", cex.axis=1, cex.lab=1.3, col=as.factor(Data$Region), pch=18, cex=1.3, xlim=c(min(Data$logFC..DeltaDeltaCT.), max(Data$logFC..DeltaDeltaCT.)))
TrendLine<-lm(Data$logFC..DeltaDeltaCT.~Data$Microarray_AVELogFC.Beta.)
abline(TrendLine, lwd=2)
TrendLineR2<-summary.lm(TrendLine)$r.squared
TrendLinePval<-summary.lm(TrendLine)$coefficients[2,4]
mtext(paste("R-squared=", signif(TrendLineR2, 2)), cex=1.3)
legend(x=-1.7, y=1.7, legend=levels(as.factor(Data$Region)), text.col=c(1:6), fill=c(1:6), cex=1.1, bty="n")
dev.off()

summary.lm(lm(Data$logFC..DeltaDeltaCT.~Data$Microarray_AVELogFC.Beta.))

Call:
  lm(formula = Data$logFC..DeltaDeltaCT. ~ Data$Microarray_AVELogFC.Beta.)

# Residuals:
#   Min      1Q  Median      3Q     Max 
# -1.0458 -0.4952 -0.1331  0.4987  1.0067 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                      0.1511     0.1143   1.322    0.198    
# Data$Microarray_AVELogFC.Beta.   1.6830     0.2223   7.569 4.91e-08 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.5969 on 26 degrees of freedom
# Multiple R-squared:  0.6879,	Adjusted R-squared:  0.6758 
# F-statistic: 57.29 on 1 and 26 DF,  p-value: 4.913e-08

#So the intercept is very close to zero, but the effects tend to be much larger in the qRT-PCR data.

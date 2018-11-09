
setwd("~/Documents/AMY LCM 10nuclei/10 Regional Comparison/MDD_RegionCompare/qRTPCR")

PCR_SST_NPY<-read.csv("SST_NPY_forGraphs.csv", header=T, stringsAsFactors = F)
str(PCR_SST_NPY)

PCR_SST_NPY$SST_AB_NegativeDeltaDeltaCT<-((PCR_SST_NPY$SST_AB_DeltaCT)-mean(PCR_SST_NPY$SST_AB_DeltaCT[PCR_SST_NPY$Diagnosis=="CTRL"]))*-1
PCR_SST_NPY$NPY_B_NegativeDeltaDeltaCT<-((PCR_SST_NPY$NPY_B_DeltaCT)-mean(PCR_SST_NPY$NPY_B_DeltaCT[PCR_SST_NPY$Diagnosis=="CTRL"]))*-1



AB_Microarray<-read.csv("AB_Microarray_ForPCRgraphs.csv", header=T, stringsAsFactors = F)
str(AB_Microarray)
AB_Microarray<-AB_Microarray[c(1:42),]


Basal_Microarray<-read.csv("Basal_Microarray_ForPCRgraphs.csv", header=T, stringsAsFactors = F)
str(Basal_Microarray)




pdf("PCR_AB_SSTvsDiagnosis.pdf", height=5, width=4)
boxplot(PCR_SST_NPY$SST_AB_NegativeDeltaDeltaCT~as.factor(PCR_SST_NPY$Diagnosis), xlab="Diagnosis", ylab="SST (DeltaDeltaCT)", lwd=1.5, cex.lab=1.3, main="AB")
stripchart(PCR_SST_NPY$SST_AB_NegativeDeltaDeltaCT~as.factor(PCR_SST_NPY$Diagnosis), vertical = TRUE, method = "jitter", add = TRUE, pch = 20, cex=1.5, col = c(3,2))
dev.off()


pdf("Microarray_AB_SSTvsDiagnosis.pdf", height=5, width=4)
boxplot(AB_Microarray$SST_ILMN_1812824~as.factor(AB_Microarray$Sample.Group), xlab="Diagnosis", ylab="SST (normalized log(2) signal)", lwd=1.5, cex.lab=1.3, main="AB")
stripchart(AB_Microarray$SST_ILMN_1812824~as.factor(AB_Microarray$Sample.Group), vertical = TRUE, method = "jitter", add = TRUE, pch = 20, cex=1.5, col = c(3,2))
dev.off()


pdf("Microarray_AB_SSTvsDiagnosis_Gender.pdf", height=5, width=5)
boxplot(AB_Microarray$SST_ILMN_1812824~as.factor(AB_Microarray$Sample.Group)*as.factor(AB_Microarray$Gender), xlab="Diagnosis", ylab="SST (normalized log(2) signal)", lwd=1.5, cex.lab=1.3, main="AB")

stripchart(AB_Microarray$SST_ILMN_1812824~as.factor(AB_Microarray$Sample.Group)*as.factor(AB_Microarray$Gender), vertical = TRUE, method = "jitter", add = TRUE, pch = 20, cex=1.5, col = c(3,2))
dev.off()

summary.lm(lm(AB_Microarray$SST_ILMN_1812824~as.factor(AB_Microarray$Sample.Group)))
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                               10.6873     0.1759  60.755  < 2e-16 ***
#   as.factor(AB_Microarray$Sample.Group)MDD  -0.9641     0.2236  -4.312 0.000103 ***
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.7036 on 40 degrees of freedom
# Multiple R-squared:  0.3173,	Adjusted R-squared:  0.3003 
# F-statistic:  18.6 on 1 and 40 DF,  p-value: 0.0001027

summary.lm(lm(AB_Microarray$SST_ILMN_1812824~as.factor(AB_Microarray$Sample.Group)*as.factor(AB_Microarray$Gender)))

# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                                               10.76018    0.50520  21.299   <2e-16 ***
#   as.factor(AB_Microarray$Sample.Group)MDD                                  -0.86769    0.55852  -1.554    0.129    
# as.factor(AB_Microarray$Gender)M                                          -0.08334    0.54008  -0.154    0.878    
# as.factor(AB_Microarray$Sample.Group)MDD:as.factor(AB_Microarray$Gender)M -0.17564    0.61516  -0.286    0.777    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.7145 on 38 degrees of freedom
# Multiple R-squared:  0.3314,	Adjusted R-squared:  0.2786 
# F-statistic: 6.278 on 3 and 38 DF,  p-value: 0.001443


pdf("Microarray_Basal_SSTvsDiagnosis_Gender.pdf", height=5, width=5)
boxplot(Basal_Microarray$SST_ILMN_1812824~as.factor(Basal_Microarray$Sample.Group)*as.factor(Basal_Microarray$Gender), xlab="Diagnosis", ylab="SST (normalized log(2) signal)", lwd=1.5, cex.lab=1.3, main="Basal")

stripchart(Basal_Microarray$SST_ILMN_1812824~as.factor(Basal_Microarray$Sample.Group)*as.factor(Basal_Microarray$Gender), vertical = TRUE, method = "jitter", add = TRUE, pch = 20, cex=1.5, col = c(3,2))
dev.off()

summary.lm(lm(Basal_Microarray$SST_ILMN_1812824~as.factor(Basal_Microarray$Sample.Group)*as.factor(Basal_Microarray$Gender)))
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                                                      11.5909     0.5470  21.190   <2e-16 ***
#   as.factor(Basal_Microarray$Sample.Group)MDD                                      -1.3689     0.6047  -2.264   0.0297 *  
#   as.factor(Basal_Microarray$Gender)M                                              -1.0598     0.5848  -1.812   0.0783 .  
# as.factor(Basal_Microarray$Sample.Group)MDD:as.factor(Basal_Microarray$Gender)M   1.6949     0.6696   2.531   0.0159 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.7736 on 36 degrees of freedom
# Multiple R-squared:  0.1649,	Adjusted R-squared:  0.0953 
# F-statistic: 2.369 on 3 and 36 DF,  p-value: 0.08676


pdf("PCR_Basal_NPYvsDiagnosis.pdf", height=5, width=4)
boxplot(PCR_SST_NPY$NPY_B_NegativeDeltaDeltaCT~as.factor(PCR_SST_NPY$Diagnosis), xlab="Diagnosis", ylab="NPY (DeltaDeltaCT)", lwd=1.5, cex.lab=1.3, main="Basal")
stripchart(PCR_SST_NPY$NPY_B_NegativeDeltaDeltaCT~as.factor(PCR_SST_NPY$Diagnosis), vertical = TRUE, method = "jitter", add = TRUE, pch = 20, cex=1.5, col = c(3,2))
dev.off()



pdf("Microarray_AB_NPYvsDiagnosis_Gender.pdf", height=5, width=5)
boxplot(AB_Microarray$NPY_ILMN_1731062~as.factor(AB_Microarray$Sample.Group)*as.factor(AB_Microarray$Gender), xlab="Diagnosis", ylab="NPY (normalized log(2) signal)", lwd=1.5, cex.lab=1.3, main="AB")

stripchart(AB_Microarray$NPY_ILMN_173106~as.factor(AB_Microarray$Sample.Group)*as.factor(AB_Microarray$Gender), vertical = TRUE, method = "jitter", add = TRUE, pch = 20, cex=1.5, col = c(3,2))
dev.off()

summary.lm(lm(AB_Microarray$NPY_ILMN_173106~as.factor(AB_Microarray$Sample.Group)*as.factor(AB_Microarray$Gender)))
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                                               12.06793    0.47277  25.526   <2e-16 ***
#   as.factor(AB_Microarray$Sample.Group)MDD                                   0.03239    0.52267   0.062    0.951    
# as.factor(AB_Microarray$Gender)M                                          -0.05030    0.50541  -0.100    0.921    
# as.factor(AB_Microarray$Sample.Group)MDD:as.factor(AB_Microarray$Gender)M  0.26043    0.57568   0.452    0.654    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.6686 on 38 degrees of freedom
# Multiple R-squared:  0.0405,	Adjusted R-squared:  -0.03525 
# F-statistic: 0.5346 on 3 and 38 DF,  p-value: 0.6613


pdf("Microarray_Basal_NPYvsDiagnosis_Gender.pdf", height=5, width=5)
boxplot(Basal_Microarray$NPY_ILMN_173106~as.factor(Basal_Microarray$Sample.Group)*as.factor(Basal_Microarray$Gender), xlab="Diagnosis", ylab="NPY (normalized log(2) signal)", lwd=1.5, cex.lab=1.3, main="Basal")

stripchart(Basal_Microarray$NPY_ILMN_173106~as.factor(Basal_Microarray$Sample.Group)*as.factor(Basal_Microarray$Gender), vertical = TRUE, method = "jitter", add = TRUE, pch = 20, cex=1.5, col = c(3,2))
dev.off()

summary.lm(lm(Basal_Microarray$NPY_ILMN_173106~as.factor(Basal_Microarray$Sample.Group)*as.factor(Basal_Microarray$Gender)))

# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                                                      13.2727     0.4320  30.721  < 2e-16 ***
#   as.factor(Basal_Microarray$Sample.Group)MDD                                      -1.5113     0.4776  -3.164  0.00316 ** 
#   as.factor(Basal_Microarray$Gender)M                                              -0.7055     0.4619  -1.527  0.13538    
# as.factor(Basal_Microarray$Sample.Group)MDD:as.factor(Basal_Microarray$Gender)M   1.2515     0.5289   2.366  0.02346 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.611 on 36 degrees of freedom
# Multiple R-squared:  0.2897,	Adjusted R-squared:  0.2305 
# F-statistic: 4.893 on 3 and 36 DF,  p-value: 0.005917


pdf("Microarray_AB_MBPvsDiagnosis_Gender.pdf", height=5, width=5)
boxplot(as.numeric(AB_Microarray$MBP_ILMN_2331544)~as.factor(AB_Microarray$Sample.Group)*as.factor(AB_Microarray$Gender), xlab="Diagnosis", ylab="MBP (normalized log(2) signal)", lwd=1.5, cex.lab=1.3, main="AB")

stripchart(as.numeric(AB_Microarray$MBP_ILMN_2331544)~as.factor(AB_Microarray$Sample.Group)*as.factor(AB_Microarray$Gender), vertical = TRUE, method = "jitter", add = TRUE, pch = 20, cex=1.5, col = c(3,2))
dev.off()

summary.lm(lm(as.numeric(AB_Microarray$MBP_ILMN_2331544)~as.factor(AB_Microarray$Sample.Group)*as.factor(AB_Microarray$Gender)))

# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                                                 7.4725     0.2044  36.566   <2e-16 ***
#   as.factor(AB_Microarray$Sample.Group)MDD                                   -0.1442     0.2259  -0.638   0.5271    
# as.factor(AB_Microarray$Gender)M                                            0.2764     0.2185   1.265   0.2136    
# as.factor(AB_Microarray$Sample.Group)MDD:as.factor(AB_Microarray$Gender)M  -0.5115     0.2488  -2.056   0.0467 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.289 on 38 degrees of freedom
# Multiple R-squared:  0.5131,	Adjusted R-squared:  0.4747 
# F-statistic: 13.35 on 3 and 38 DF,  p-value: 4.234e-06


pdf("Microarray_Basal_MBPvsDiagnosis_Gender.pdf", height=5, width=5)
boxplot(as.numeric(Basal_Microarray$MBP_ILMN_2331544)~as.factor(Basal_Microarray$Sample.Group)*as.factor(Basal_Microarray$Gender), xlab="Diagnosis", ylab="MBP (normalized log(2) signal)", lwd=1.5, cex.lab=1.3, main="Basal")

stripchart(as.numeric(Basal_Microarray$MBP_ILMN_2331544)~as.factor(Basal_Microarray$Sample.Group)*as.factor(Basal_Microarray$Gender), vertical = TRUE, method = "jitter", add = TRUE, pch = 20, cex=1.5, col = c(3,2))
dev.off()

summary.lm(lm(as.numeric(Basal_Microarray$MBP_ILMN_2331544)~as.factor(Basal_Microarray$Sample.Group)*as.factor(Basal_Microarray$Gender)))

# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                                                       7.1138     0.2135  33.321   <2e-16 ***
#   as.factor(Basal_Microarray$Sample.Group)MDD                                       0.2024     0.2360   0.857   0.3969    
# as.factor(Basal_Microarray$Gender)M                                               0.5682     0.2282   2.489   0.0176 *  
#   as.factor(Basal_Microarray$Sample.Group)MDD:as.factor(Basal_Microarray$Gender)M  -0.3894     0.2613  -1.490   0.1450    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.3019 on 36 degrees of freedom
# Multiple R-squared:  0.2452,	Adjusted R-squared:  0.1823 
# F-statistic: 3.899 on 3 and 36 DF,  p-value: 0.01645


pdf("Microarray_AB_MOBPvsDiagnosis_Gender.pdf", height=5, width=5)
boxplot(as.numeric(AB_Microarray$MOBP_ILMN_1768947)~as.factor(AB_Microarray$Sample.Group)*as.factor(AB_Microarray$Gender), xlab="Diagnosis", ylab="MOBP (normalized log(2) signal)", lwd=1.5, cex.lab=1.3, main="AB")

stripchart(as.numeric(AB_Microarray$MOBP_ILMN_1768947)~as.factor(AB_Microarray$Sample.Group)*as.factor(AB_Microarray$Gender), vertical = TRUE, method = "jitter", add = TRUE, pch = 20, cex=1.5, col = c(3,2))
dev.off()


pdf("Microarray_Basal_MOBPvsDiagnosis_Gender.pdf", height=5, width=5)
boxplot(as.numeric(Basal_Microarray$MOBP_ILMN_1768947)~as.factor(Basal_Microarray$Sample.Group)*as.factor(Basal_Microarray$Gender), xlab="Diagnosis", ylab="MOBP (normalized log(2) signal)", lwd=1.5, cex.lab=1.3, main="Basal")

stripchart(as.numeric(Basal_Microarray$MOBP_ILMN_1768947)~as.factor(Basal_Microarray$Sample.Group)*as.factor(Basal_Microarray$Gender), vertical = TRUE, method = "jitter", add = TRUE, pch = 20, cex=1.5, col = c(3,2))
dev.off()


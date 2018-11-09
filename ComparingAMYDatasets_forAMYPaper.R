#This was adapted from earlier code to compare amygdala results across datasets:
#Megan Hagenauer, 07-20-2018

setwd("~/Documents/R Code/MakingAMetaQueryDatabase/thefunction")

#First I inputted the functions from "IlluminaAndAffyIncBarnesNarayan.R" - Alek Pankonin's meta-query code for the Pritzker datasets.

#Example Function usage:
library(plyr)


setwd("~/Documents/R Code/MakingAMetaQueryDatabase/thefunction/LegendOfTheHiddenDataset_2017_04_07")

#3 necessary info files:

illuminaProbeInfo <- read.csv("IlluminaProbeInfo.csv", header = T)
ILMNDataSetModel <- read.csv("DatasetModel.csv", header = T)
AffyDataSetModel <- read.csv("AffyDataSetModel.csv", header = T)


setwd("~/Documents/R Code/MakingAMetaQueryDatabase/ComparingDatasets/AMY_UsingBetas")

All_Pritzker960AMY_Genes<-read.csv("AllDetectedGenes_Pritzker960_AMY.csv", header=T, stringsAsFactors = F)
str(All_Pritzker960AMY_Genes)
genes<-as.vector(All_Pritzker960AMY_Genes[,1])
str(genes)
#chr [1:19451] "A1BG" "A2BP1" "A2LD1" "A2M" "A2ML1" "A4GALT" "AAA1" "AAAS" "AACS" "AACSL" "AADACL1" "AADACL4" "AADAT" "AAGAB" "AAK1" "AAMP" "AARS" "AARS2" ...

testdatasets <- c("Pritzker960", "Freeze3", "LCM.AMY", "Affy6Region")
testtypeOfOutput <- c("Beta")
testvariablesInterest <- c("MDD")
testBrainRegion <- c("AMY", "AAA", "AB", "AHA", "Basal", "Central", "CO", "Lateral", "Medial", "PAC", "PL")

setwd("~/Documents/R Code/MakingAMetaQueryDatabase/AllDatasetsReformatted")

system.time(
  test <- AffyAndIllumina(ILMNDataModel = ILMNDataSetModel, AFFYDataModel = AffyDataSetModel, genesOfInterest = genes, DataSets = testdatasets, BrainRegion = testBrainRegion, variablesofInterest = testvariablesInterest, OutputsInterest = testtypeOfOutput)
)


TestOutputAffy<-test[[2]]
TestOutputIllumina<-test[[1]]
TestOutputJoined<-test[[3]]

str(TestOutputJoined)
str(TestOutputIllumina)

setwd("~/Documents/R Code/MakingAMetaQueryDatabase/ComparingDatasets/AMY_UsingBetas")

write.csv(TestOutputAffy, "TestOutputAffy.csv")
write.csv(TestOutputIllumina, "TestOutputIllumina.csv")
write.csv(TestOutputJoined, "TestOutputJoined.csv")

colnames(TestOutputIllumina)
temp<-cbind(TestOutputIllumina[,c(19:30)])

write.csv(cor(temp, use="complete.obs"), "MDD_Dataset_CorMatrix.csv")

colnames(TestOutputJoined)


testtypeOfOutput <- c("Beta", "NominalPvalue", "FDR")

setwd("~/Documents/R Code/MakingAMetaQueryDatabase/AllDatasetsReformatted")

system.time(
  test <- AffyAndIllumina(ILMNDataModel = ILMNDataSetModel, AFFYDataModel = AffyDataSetModel, genesOfInterest = genes, DataSets = testdatasets, BrainRegion = testBrainRegion, variablesofInterest = testvariablesInterest, OutputsInterest = testtypeOfOutput)
)


TestOutputAffy<-test[[2]]
TestOutputIllumina<-test[[1]]
TestOutputJoined<-test[[3]]

str(TestOutputJoined)
str(TestOutputIllumina)

setwd("~/Documents/R Code/MakingAMetaQueryDatabase/ComparingDatasets/AMY_UsingBetas")

write.csv(TestOutputAffy, "TestOutputAffy_MoreDetail.csv")
write.csv(TestOutputIllumina, "TestOutputIllumina_MoreDetail.csv")
write.csv(TestOutputJoined, "TestOutputJoined_MoreDetail.csv")


#************
#Just pulling out the results for the list of AMY genes associated with MDD and CMS in Sibille 2009:
str(genes)


setwd("~/Documents/R Code/MakingAMetaQueryDatabase/ComparingDatasets/AMY_UsingBetas")

Sibille_AMY_Genes<-read.csv("Sibille2009List_MDD_CMS_GenesAMY.csv", header=T, stringsAsFactors = F)
str(Sibille_AMY_Genes)
genes<-as.vector(Sibille_AMY_Genes[,1])
str(genes)

setwd("~/Documents/R Code/MakingAMetaQueryDatabase/AllDatasetsReformatted")

system.time(
  test <- AffyAndIllumina(ILMNDataModel = ILMNDataSetModel, AFFYDataModel = AffyDataSetModel, genesOfInterest = genes, DataSets = testdatasets, BrainRegion = testBrainRegion, variablesofInterest = testvariablesInterest, OutputsInterest = testtypeOfOutput)
)


SibilleOutputIllumina<-test[[1]]
SibilleOutputJoined<-test[[3]]


setwd("~/Documents/R Code/MakingAMetaQueryDatabase/ComparingDatasets/AMY_UsingBetas")

write.csv(SibilleOutputIllumina, "SibilleOutputIllumina_MoreDetail.csv")
write.csv(SibilleOutputJoined, "SibilleOutputJoined_MoreDetail.csv")

str(SibilleOutputIllumina)
colnames(SibilleOutputIllumina)
library(plyr)
colnames(Sibille_AMY_Genes)[1]<-"GeneSymbol"
colnames(SibilleOutputIllumina)[3]<-"GeneSymbol"

sum(is.na(SibilleOutputIllumina$GeneSymbol))
#[1] 5

SibilleOutputIllumina<-SibilleOutputIllumina[is.na(SibilleOutputIllumina$GeneSymbol)==F,]
dim(SibilleOutputIllumina)
#[1] 59 54

Sibille_AMY_Vs_OurAMYData<-join(Sibille_AMY_Genes, SibilleOutputIllumina[,c(3,19:30)], by="GeneSymbol", type="full")

dim(Sibille_AMY_Vs_OurAMYData)
#[1] 61 14

head(Sibille_AMY_Vs_OurAMYData)

Sibille_AMY_Vs_OurAMYData$GeneSymbol

CorMatrix_Sibille_AMY_Vs_OurAMYData<-cor(Sibille_AMY_Vs_OurAMYData[,-1], use="complete.obs")

write.csv(CorMatrix_Sibille_AMY_Vs_OurAMYData, "CorMatrix_Sibille_AMY_Vs_OurAMYData.csv")


pdf("LateralLogFC_VsSibilleDirectionOfEffect.pdf", height=5, width=4)
boxplot(Sibille_AMY_Vs_OurAMYData$Lateral_Illumina_LCM.AMY_ModelLM2_MDD_Beta~as.factor(Sibille_AMY_Vs_OurAMYData$Direction.of.Effect), xlab="Sibille (2009): Direction of Effect", ylab="Lateral: Log(2) Fold Change", lwd=1.5, cex.lab=1.3)

stripchart(Sibille_AMY_Vs_OurAMYData$Lateral_Illumina_LCM.AMY_ModelLM2_MDD_Beta~as.factor(Sibille_AMY_Vs_OurAMYData$Direction.of.Effect), vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 20, cex=1.5, col = c(4,2))
dev.off()

summary.lm(lm(Sibille_AMY_Vs_OurAMYData$Lateral_Illumina_LCM.AMY_ModelLM2_MDD_Beta~as.factor(Sibille_AMY_Vs_OurAMYData$Direction.of.Effect)))

# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                               -0.21991    0.04802  -4.579  4.9e-05 ***
#   as.factor(Sibille_AMY_Vs_OurAMYData$Direction.of.Effect)1  0.28940    0.08424   3.435  0.00145 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.2495 on 38 degrees of freedom
# (21 observations deleted due to missingness)
# Multiple R-squared:  0.237,	Adjusted R-squared:  0.2169 
# F-statistic:  11.8 on 1 and 38 DF,  p-value: 0.001446

summary.lm(lm(Sibille_AMY_Vs_OurAMYData$AB_Illumina_LCM.AMY_ModelLM2_MDD_Beta~as.factor(Sibille_AMY_Vs_OurAMYData$Direction.of.Effect)))

# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)                                                -0.1202     0.1072  -1.121    0.270
# as.factor(Sibille_AMY_Vs_OurAMYData$Direction.of.Effect)1   0.2629     0.1766   1.489    0.145
# 
# Residual standard error: 0.5252 on 36 degrees of freedom
# (23 observations deleted due to missingness)
# Multiple R-squared:  0.05799,	Adjusted R-squared:  0.03182 
# F-statistic: 2.216 on 1 and 36 DF,  p-value: 0.1453



pdf("BasalLogFC_VsSibilleDirectionOfEffect.pdf", height=5, width=4)
boxplot(Sibille_AMY_Vs_OurAMYData$Basal_Illumina_LCM.AMY_ModelLM2_MDD_Beta~as.factor(Sibille_AMY_Vs_OurAMYData$Direction.of.Effect), xlab="Sibille (2009): Direction of Effect", ylab="Basal: Log(2) Fold Change", lwd=1.5, cex.lab=1.3)

stripchart(Sibille_AMY_Vs_OurAMYData$Basal_Illumina_LCM.AMY_ModelLM2_MDD_Beta~as.factor(Sibille_AMY_Vs_OurAMYData$Direction.of.Effect), vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 20, cex=1.5, col = c(4,2))
dev.off()

pdf("ABLogFC_VsSibilleDirectionOfEffect.pdf", height=5, width=4)
boxplot(Sibille_AMY_Vs_OurAMYData$AB_Illumina_LCM.AMY_ModelLM2_MDD_Beta~as.factor(Sibille_AMY_Vs_OurAMYData$Direction.of.Effect), xlab="Sibille (2009): Direction of Effect", ylab="AB: Log(2) Fold Change", lwd=1.5, cex.lab=1.3)

stripchart(Sibille_AMY_Vs_OurAMYData$AB_Illumina_LCM.AMY_ModelLM2_MDD_Beta~as.factor(Sibille_AMY_Vs_OurAMYData$Direction.of.Effect), vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 20, cex=1.5, col = c(4,2))
dev.off()


pdf("AMY_Freeze3_LogFC_VsSibilleDirectionOfEffect.pdf", height=5, width=4)
boxplot(Sibille_AMY_Vs_OurAMYData$AMY_Illumina_Freeze3_ModelLM4_MDD_Beta~as.factor(Sibille_AMY_Vs_OurAMYData$Direction.of.Effect), xlab="Sibille (2009): Direction of Effect", ylab="Amygdala (Ref-8): Log(2) Fold Change", lwd=1.5, cex.lab=1.3)

stripchart(Sibille_AMY_Vs_OurAMYData$AMY_Illumina_Freeze3_ModelLM4_MDD_Beta~as.factor(Sibille_AMY_Vs_OurAMYData$Direction.of.Effect), vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 20, cex=1.5, col = c(4,2))
dev.off()

pdf("AMY_Pritzker960_LogFC_VsSibilleDirectionOfEffect.pdf", height=5, width=4)
boxplot(Sibille_AMY_Vs_OurAMYData$AMY_Illumina_Pritzker960_ModelLM4_MDD_Beta~as.factor(Sibille_AMY_Vs_OurAMYData$Direction.of.Effect), xlab="Sibille (2009): Direction of Effect", ylab="Amygdala (HT-12v3): Log(2) Fold Change", lwd=1.5, cex.lab=1.3)

stripchart(Sibille_AMY_Vs_OurAMYData$AMY_Illumina_Pritzker960_ModelLM4_MDD_Beta~as.factor(Sibille_AMY_Vs_OurAMYData$Direction.of.Effect), vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 20, cex=1.5, col = c(4,2))
dev.off()


#********************************************

#************
#Just pulling out the results for the list of AMY genes associated with MDD and CMS in Guillox 2012:
str(genes)

library(plyr)

setwd("~/Documents/R Code/MakingAMetaQueryDatabase/ComparingDatasets/AMY_UsingBetas")

Guillox_AMY_Genes<-read.csv("Guillox_2012_MDDAmygdalaMicroarray_Suppl_Table1.csv", header=T, stringsAsFactors = F)
str(Guillox_AMY_Genes)
# 'data.frame':	116 obs. of  3 variables:
#   $ PROBE         : chr  "ILMN_1812824 " "ILMN_1765966 " "ILMN_1679984" "ILMN_1697512 " ...
# $ GeneSymbol    : chr  "SST" "CHGB" "ZCCHC12" "SLC32A1" ...
# $ MDD.FoldChange: num  -0.84 -0.54 -0.53 -0.48 -0.39 -0.38 -0.38 -0.36 -0.34 -0.31 ...

genes<-as.vector(Guillox_AMY_Genes[,2])
str(genes)
# chr [1:116] "SST" "CHGB" "ZCCHC12" "SLC32A1" "AMPH" "TAC1" "ELMOD1" "NEFH" "TUSC3" "GNG3" "RGS7BP" "PFKM" "ACOT7" "KCNG1" "LONRF2" 
sum(is.na(genes))
#[1] 6

genes<-genes[is.na(genes)==F]
str(genes)
#chr [1:110] "SST" "CHGB" "ZCCHC12" "SLC32A1" "AMPH" "TAC1" "ELMOD1" "NEFH" "TUSC3" "GNG3" "RGS7BP" "PFKM" "ACOT7" "KCNG1" "LONRF2" ...

setwd("~/Documents/R Code/MakingAMetaQueryDatabase/AllDatasetsReformatted")

system.time(
  test <- AffyAndIllumina(ILMNDataModel = ILMNDataSetModel, AFFYDataModel = AffyDataSetModel, genesOfInterest = genes, DataSets = testdatasets, BrainRegion = testBrainRegion, variablesofInterest = testvariablesInterest, OutputsInterest = testtypeOfOutput)
)


GuilloxOutputIllumina<-test[[1]]
GuilloxOutputJoined<-test[[3]]


setwd("~/Documents/R Code/MakingAMetaQueryDatabase/ComparingDatasets/AMY_UsingBetas")

write.csv(GuilloxOutputIllumina, "GuilloxOutputIllumina_MoreDetail.csv")
write.csv(GuilloxOutputJoined, "GuilloxOutputJoined_MoreDetail.csv")

str(GuilloxOutputIllumina)
#'data.frame':	209 obs. of  54 variables:

colnames(GuilloxOutputIllumina)

colnames(Guillox_AMY_Genes)[2]<-"GeneSymbol"
colnames(GuilloxOutputIllumina)[3]<-"GeneSymbol"

sum(is.na(GuilloxOutputIllumina$GeneSymbol))
#[1] 24

GuilloxOutputIllumina<-GuilloxOutputIllumina[is.na(GuilloxOutputIllumina$GeneSymbol)==F,]
dim(GuilloxOutputIllumina)
#[1] 185  54

Guillox_AMY_Genes<-Guillox_AMY_Genes[is.na(Guillox_AMY_Genes$GeneSymbol)==F,]

Guillox_AMY_Vs_OurAMYDataDetailed<-join(Guillox_AMY_Genes, GuilloxOutputIllumina, by="GeneSymbol", type="full")
write.csv(Guillox_AMY_Vs_OurAMYDataDetailed, "Guillox_AMY_Vs_OurAMYDataDetailed.csv")


Guillox_AMY_Vs_OurAMYData<-join(Guillox_AMY_Genes, GuilloxOutputIllumina[,c(3,19:30)], by="GeneSymbol", type="full")

dim(Guillox_AMY_Vs_OurAMYData)
#[1] 205  15

head(Guillox_AMY_Vs_OurAMYData)

Guillox_AMY_Vs_OurAMYData$GeneSymbol

CorMatrix_Guillox_AMY_Vs_OurAMYData<-cor(Guillox_AMY_Vs_OurAMYData[,-c(1:2)], use="complete.obs")

write.csv(CorMatrix_Guillox_AMY_Vs_OurAMYData, "CorMatrix_Guillox_AMY_Vs_OurAMYData.csv")

#So the microarray results from all women look like our results from the medial nucleus instead of the basolateral nuclei???


pdf("MedialLogFC_VsGuilloxLogFoldChange.pdf", height=5, width=4)
plot(Guillox_AMY_Vs_OurAMYData$Medial_Illumina_LCM.AMY_ModelLM2_MDD_Beta~Guillox_AMY_Vs_OurAMYData$MDD.FoldChange, xlab="Guilloux et al. (2012): Log(2) FoldChange", ylab="Medial: Log(2) Fold Change", lwd=1.5, cex.lab=1.3)
BestfitLine<-lm(Guillox_AMY_Vs_OurAMYData$Medial_Illumina_LCM.AMY_ModelLM2_MDD_Beta~Guillox_AMY_Vs_OurAMYData$MDD.FoldChange)
abline(BestfitLine, col=2, lwd=1.5)
dev.off()

summary.lm(lm(Guillox_AMY_Vs_OurAMYData$Medial_Illumina_LCM.AMY_ModelLM2_MDD_Beta~Guillox_AMY_Vs_OurAMYData$MDD.FoldChange))

# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                              -0.06201    0.01976  -3.138  0.00209 ** 
#   Guillox_AMY_Vs_OurAMYData$MDD.FoldChange  0.42962    0.05525   7.776 1.71e-12 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.2286 on 135 degrees of freedom
# (68 observations deleted due to missingness)
# Multiple R-squared:  0.3094,	Adjusted R-squared:  0.3043 
# F-statistic: 60.47 on 1 and 135 DF,  p-value: 1.71e-12




pdf("Freeze3LogFC_VsGuilloxLogFoldChange.pdf", height=5, width=5)
plot(Guillox_AMY_Vs_OurAMYData$AMY_Illumina_Freeze3_ModelLM4_MDD_Beta~Guillox_AMY_Vs_OurAMYData$MDD.FoldChange, xlab="Guilloux et al. (2012): Log(2) FoldChange", ylab="Amygdala (HT12v3): Log(2) Fold Change", lwd=1.5, cex.lab=1.3)
BestfitLine<-lm(Guillox_AMY_Vs_OurAMYData$AMY_Illumina_Freeze3_ModelLM4_MDD_Beta~Guillox_AMY_Vs_OurAMYData$MDD.FoldChange)
abline(BestfitLine, col=2, lwd=1.5)
dev.off()

summary.lm(lm(Guillox_AMY_Vs_OurAMYData$AMY_Illumina_Freeze3_ModelLM4_MDD_Beta~Guillox_AMY_Vs_OurAMYData$MDD.FoldChange))

# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                               0.006150   0.009239   0.666    0.507    
# Guillox_AMY_Vs_OurAMYData$MDD.FoldChange -0.187426   0.026523  -7.067 2.25e-10 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.09174 on 99 degrees of freedom
# (104 observations deleted due to missingness)
# Multiple R-squared:  0.3353,	Adjusted R-squared:  0.3286 
# F-statistic: 49.94 on 1 and 99 DF,  p-value: 2.252e-10


pdf("LateralLogFC_VsGuilloxLogFoldChange.pdf", height=5, width=5)
plot(Guillox_AMY_Vs_OurAMYData$Lateral_Illumina_LCM.AMY_ModelLM2_MDD_Beta~Guillox_AMY_Vs_OurAMYData$MDD.FoldChange, xlab="Guilloux et al. (2012): Log(2) FoldChange", ylab="Lateral: Log(2) Fold Change", lwd=1.5, cex.lab=1.3)
BestfitLine<-lm(Guillox_AMY_Vs_OurAMYData$Lateral_Illumina_LCM.AMY_ModelLM2_MDD_Beta~Guillox_AMY_Vs_OurAMYData$MDD.FoldChange)
abline(BestfitLine, col=2, lwd=1.5)
dev.off()

summary.lm(lm(Guillox_AMY_Vs_OurAMYData$Lateral_Illumina_LCM.AMY_ModelLM2_MDD_Beta~Guillox_AMY_Vs_OurAMYData$MDD.FoldChange))

# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                              -0.03409    0.01804  -1.889   0.0611 .  
# Guillox_AMY_Vs_OurAMYData$MDD.FoldChange -0.27606    0.05039  -5.479 2.19e-07 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.2024 on 128 degrees of freedom
# (75 observations deleted due to missingness)
# Multiple R-squared:   0.19,	Adjusted R-squared:  0.1836 
# F-statistic: 30.02 on 1 and 128 DF,  p-value: 2.185e-07


pdf("BasalLogFC_VsGuilloxLogFoldChange.pdf", height=5, width=5)
plot(Guillox_AMY_Vs_OurAMYData$Basal_Illumina_LCM.AMY_ModelLM2_MDD_Beta~Guillox_AMY_Vs_OurAMYData$MDD.FoldChange, xlab="Guilloux et al. (2012): Log(2) FoldChange", ylab="Basal: Log(2) Fold Change", lwd=1.5, cex.lab=1.3)
BestfitLine<-lm(Guillox_AMY_Vs_OurAMYData$Basal_Illumina_LCM.AMY_ModelLM2_MDD_Beta~Guillox_AMY_Vs_OurAMYData$MDD.FoldChange)
abline(BestfitLine, col=2, lwd=1.5)
dev.off()

summary.lm(lm(Guillox_AMY_Vs_OurAMYData$Basal_Illumina_LCM.AMY_ModelLM2_MDD_Beta~Guillox_AMY_Vs_OurAMYData$MDD.FoldChange))

# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                              -0.08122    0.03073  -2.643   0.0092 ** 
#   Guillox_AMY_Vs_OurAMYData$MDD.FoldChange -0.52509    0.08647  -6.072 1.26e-08 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.3518 on 132 degrees of freedom
# (71 observations deleted due to missingness)
# Multiple R-squared:  0.2183,	Adjusted R-squared:  0.2124 
# F-statistic: 36.87 on 1 and 132 DF,  p-value: 1.257e-08

#********************************************

pdf("MDD_Freeze3vsPritzker960.pdf", width=6, height=6)
plot(TestOutputIllumina$AMY_Illumina_Freeze3_ModelLM4_MDD_Beta~TestOutputIllumina$AMY_Illumina_Pritzker960_ModelLM4_MDD_Beta, ylab="Freeze3 LogFC", xlab="Pritzker960 LogFC", main="MDD", cex.axis=1.3, cex.lab=1.3)
abline(lm(TestOutputIllumina$AMY_Illumina_Freeze3_ModelLM4_MDD_Beta~TestOutputIllumina$AMY_Illumina_Pritzker960_ModelLM4_MDD_Beta), col=2, lwd=3)
mtext(paste("R-squared=", signif(summary.lm(lm(TestOutputIllumina$AMY_Illumina_Freeze3_ModelLM4_MDD_Beta~TestOutputIllumina$AMY_Illumina_Pritzker960_ModelLM4_MDD_Beta))$r.squared, 3)), cex=1.3)
dev.off()

#Hmm... since they are similar, what if we averaged the two macro-dissected datasets and then ran comparisons with the other nuclei?

library(jmvcore)

MeanLogFC_Freeze3_Pritzker960<-apply(cbind(TestOutputIllumina$AMY_Illumina_Freeze3_ModelLM4_MDD_Beta, TestOutputIllumina$AMY_Illumina_Pritzker960_ModelLM4_MDD_Beta), 1, function(y) if(isError(mean(y))==T){return(NA)}else{mean(y)})

str(MeanLogFC_Freeze3_Pritzker960)
hist(MeanLogFC_Freeze3_Pritzker960[is.na(MeanLogFC_Freeze3_Pritzker960)==F], breaks=100)
#Excellent...

colnames(TestOutputIllumina)
temp<-cbind(MeanLogFC_Freeze3_Pritzker960, TestOutputIllumina[,c(19:30)])

write.csv(cor(temp, use="complete.obs"), "MDD_Dataset_CorMatrix_wAveFreeze3Pritzker960.csv")


#############

#How much do the samples overlap?

setwd("~/Documents/Freeze3/Amy/Demographics Files")

FreezeDemographics<-read.csv("Amy_SubjectDemographics_noNAs_OrderedBySubject.csv", header=T)

setwd("~/Documents/R Code/Amy/2 Input")

Pritzker960Demographics<-read.csv("Pritzker960SubjectInfo.csv", header=T)

setwd("~/Documents/AMY LCM 10nuclei/AB/02 Input")

AmyLCMDemographics<-read.csv("Sample ID file extended for Rv2.csv", header=T)

str(FreezeDemographics)

FreezeDemographics$subject_id
Pritzker960Demographics$Brain.number
AmyLCMDemographics$Subject.Number

Pritzker960Outliers<-c(4656, 3452, 2292, 3588)

Pritzker960Demographics$Outliers<-Pritzker960Demographics$Brain.number%in%Pritzker960Outliers

length(FreezeDemographics$subject_id)
#[1] 94

length(Pritzker960Demographics$Brain.number[Pritzker960Demographics$Outliers==F])
#[1] 92

sum(FreezeDemographics$subject_id%in%Pritzker960Demographics$Brain.number[Pritzker960Demographics$Outliers==F])
#[1] 55

55/94
#59%
55/92
#60%


  
#I wonder how that breaks down by diagnosis:

sum(FreezeDemographics$Disease=="C")
#[1] 40

sum(Pritzker960Demographics$Cohort=="Control"& Pritzker960Demographics$Outliers==F)
#[1] 22

sum(FreezeDemographics$subject_id[FreezeDemographics$Disease=="C"]%in%Pritzker960Demographics$Brain.number[Pritzker960Demographics$Cohort=="Control"& Pritzker960Demographics$Outliers==F])
#[1] 17
17/22
#So almost all of the Pritzker960 control subjects are found in Freeze 3 (77%)
17/40
#and 42% of the Freeze3 control subjects are in Pritzker960

sum(FreezeDemographics$Disease=="MD")
#[1] 32

sum(Pritzker960Demographics$Cohort=="MD"& Pritzker960Demographics$Outliers==F)
#[1] 24

sum(FreezeDemographics$subject_id[FreezeDemographics$Disease=="MD"]%in%Pritzker960Demographics$Brain.number[Pritzker960Demographics$Cohort=="MD"& Pritzker960Demographics$Outliers==F])
#[1] 19
19/24
#So almost all of the Pritzker960 MDD subjects are found in Freeze 3 (79%)
19/32
#... and 59% of the Freeze3 MDD subjects are in Pritzker960


sum(FreezeDemographics$Disease=="SC")
#[1] 10
sum(FreezeDemographics$Disease=="BP")
#[1] 12
#The Freeze3 sample sizes for SC and BP are really dinky.

sum(Pritzker960Demographics$Cohort=="SCHZ"& Pritzker960Demographics$Outliers==F)
#23
sum(Pritzker960Demographics$Cohort=="BP"& Pritzker960Demographics$Outliers==F)


#Overlap with AMY-LCM:

length(AmyLCMDemographics$Subject.Number)
#[1] 42

sum(AmyLCMDemographics$Subject.Number%in%Pritzker960Demographics$Brain.number[Pritzker960Demographics$Outliers==F])
#[1] 4

#Almost no overlap - awesome!

AmyLCMDemographics[AmyLCMDemographics$Subject.Number%in%Pritzker960Demographics$Brain.number,]

# Sample.Processing.Order GeneralChip LocationOnChip Sample.Group Subject.Number Cohort Gender Age Age.of.RNA..years. TOD..hrs. Suicide Overdose   pH RNAConc
# 5                        5           1              E         CTRL           4754     10      F  55           5.186301       5.0       0        0 6.79    10.6
# 30                      30           2              L          MDD           4867     10      M  48           4.523288      18.0       1        0 6.65    10.2
# 33                      33           3              C          MDD           4949     10      M  56           3.816438      14.5       1        0 6.96    11.7
# 34                      34           3              D          MDD           4982     10      M  49           3.539726       7.0       1        0 6.83     9.4
# RNAIntegrity Hours.Cold Hours.Ice  PMI HoursColdPlusIce PMICorrected Amplification.Efficiency
# 5           5.1          3        10  8.2               13         13.0                 39.62511
# 30          3.8          3        16 18.5               19         19.0                 53.94009
# 33          4.9          4        19 26.5               23         26.5                 42.11840
# 34          6.0          3        10 16.0               13         16.0                 42.74702

#1 CNTRL, 3 MDD

sum(AmyLCMDemographics$Sample.Group=="CTRL")
#[1] 16

sum(AmyLCMDemographics$Sample.Group=="MDD")
#[1] 26

1/16
#6.25% of control subjects found in Amy-LCM also in Pritzker960
1/22
#4.55% of control subjects found in Pritzker960 also found in Amy-LCM

3/26
#11.54% of MDD subjects found in Amy-LCM also in Pritzker960

3/24
#12.5% of MDD subjects found in Pritzker960 also found in Amy-LCM 

sum(AmyLCMDemographics$Subject.Number%in%FreezeDemographics$subject_id)
#[1] 0



#############


#Making some volcano plots for the macro-dissected Amygdala datasets:

setwd("~/Documents/R Code/MakingAMetaQueryDatabase/AllDatasetsReformatted")

AMY_Illumina_Freeze3_DEResults<-read.csv("AMY_Illumina_Freeze3_ModelLM4.csv", header=T)

AMY_Illumina_Pritzker960_DEResults<-read.csv("AMY_Illumina_Pritzker960_ModelLM4.csv", header=T)

#Adding a volcano plot
#http://www.gettinggeneticsdone.com/2014/05/r-volcano-plots-to-visualize-rnaseq-microarray.html


res<-data.frame(AMY_Illumina_Freeze3_DEResults$SYMBOLREANNOTATED, AMY_Illumina_Freeze3_DEResults$AMY_Illumina_Freeze3_ModelLM4_MDD_Beta, AMY_Illumina_Freeze3_DEResults$AMY_Illumina_Freeze3_ModelLM4_MDD_NominalPvalue, AMY_Illumina_Freeze3_DEResults$AMY_Illumina_Freeze3_ModelLM4_MDD_FDR, stringsAsFactors = F) 
colnames(res)<-c("gene", "log2FoldChange", "pvalue", "padj")

head(res)

setwd("~/Documents/R Code/MakingAMetaQueryDatabase/ComparingDatasets/AMY_UsingBetas")

pdf("VolcanoPlot_Freeze3_AMY.pdf", height=4, width=4)
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot: Freeze3 AMY", xlim=c(-3,3), cex.lab=1.3, cex=0.6))

# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", cex=0.6))
with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", cex=0.6))
with(subset(res, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="green", cex=0.6))
dev.off()

res<-data.frame(AMY_Illumina_Pritzker960_DEResults$SYMBOLREANNOTATED, AMY_Illumina_Pritzker960_DEResults$AMY_Illumina_Pritzker960_ModelLM4_MDD_Beta, AMY_Illumina_Pritzker960_DEResults$AMY_Illumina_Pritzker960_ModelLM4_MDD_NominalPvalue, AMY_Illumina_Pritzker960_DEResults$AMY_Illumina_Pritzker960_ModelLM4_MDD_FDR, stringsAsFactors = F) 
colnames(res)<-c("gene", "log2FoldChange", "pvalue", "padj")
str(res)
head(res)

pdf("VolcanoPlot_Pritzker960_AMY.pdf", height=4, width=4)
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot: Pritzker960 AMY", xlim=c(-3,3), cex.lab=1.3, cex=0.6))

# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", cex=0.6))
with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", cex=0.6))
with(subset(res, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="green", cex=0.6))
dev.off()


#############

#What happens if I run GSEA on just the average output from the macro-dissected datasets?

str(TestOutputIllumina)
str(MeanLogFC_Freeze3_Pritzker960)

library(fgsea)

setwd("~/Documents/Affy/NoPC1correct/DLPFC circadian/fGSEA/GMTfiles_Human")
temp<-gmtPathways("CombinedGMT_C2_C5_JustTraditional_CellType_20170928.gmt.txt")
names(temp)
temp[[1]]


MDD_betas_forGSEA_asVector<-MeanLogFC_Freeze3_Pritzker960
names(MDD_betas_forGSEA_asVector)<-TestOutputIllumina$SYMBOLREANNOTATED
str(MDD_betas_forGSEA_asVector)

#Remove NAs:
sum(is.na(TestOutputIllumina$SYMBOLREANNOTATED))
#[1] 2713
sum(is.na(MeanLogFC_Freeze3_Pritzker960))
#[1] 15936



MDD_betas_forGSEA_asVector_NoNA<-MDD_betas_forGSEA_asVector[is.na(TestOutputIllumina$SYMBOLREANNOTATED)==F & is.na(MeanLogFC_Freeze3_Pritzker960)==F]
str(MDD_betas_forGSEA_asVector_NoNA)

# Named num [1:11891] -0.1311 0.0659 0.0041 0.0161 -0.0686 ...
# - attr(*, "names")= chr [1:11891] "A2M" "A4GALT" "AAAS" "AACS" ...

NoDuplicatedGeneNames<-tapply(MDD_betas_forGSEA_asVector_NoNA, names(MDD_betas_forGSEA_asVector_NoNA), FUN=mean)
str(NoDuplicatedGeneNames)
# num [1:10574(1d)] -0.1311 0.0659 -0.0553 0.0041 0.0161 ...
# - attr(*, "dimnames")=List of 1
# ..$ : chr [1:10574] "A2M" "A4GALT" "AAA1" "AAAS" ...

names(NoDuplicatedGeneNames)[c(1:100)]
length(NoDuplicatedGeneNames)

setwd("~/Documents/R Code/MakingAMetaQueryDatabase/ComparingDatasets/AMY_UsingBetas")

MDD_betas_forGSEA_noDuplicatesRanked<-NoDuplicatedGeneNames[order(NoDuplicatedGeneNames)]
head(MDD_betas_forGSEA_noDuplicatesRanked)
write.csv(MDD_betas_forGSEA_noDuplicatesRanked, "MDD_AverageBetasFreeze3Pritzker960_forGSEA_noDuplicatesRanked.csv")


temp1<-fgsea(temp, MDD_betas_forGSEA_noDuplicatesRanked, nperm=10000, minSize = 15, maxSize = 500)
str(temp1)


temp1$leadingEdge<-vapply(temp1$leadingEdge, paste, collapse= ",", character(1L))

write.csv(temp1, "AverageBetasFreeze3Pritzker960_FGSEAResultswCellType2.csv")


##################

# Making a nice plot for the stereology results:


setwd("~/Documents/AMY LCM 10nuclei/Manuscript/Figures")

StereologyResults<-read.csv("StereologyResults.csv", header=T)
str(StereologyResults)

print(levels(StereologyResults$Region))
[1] "AB"      "Basal"   "Central" "Lateral" "Medial" 

StereologyResults$Region<-factor(StereologyResults$Region,levels(StereologyResults$Region)[c(4,2,1,3,5)])

print(levels(StereologyResults$Region))
[1] "Lateral" "Basal"   "AB"      "Central" "Medial" 

pdf("StereologyResults.pdf", height=8, width=4)
boxplot(StereologyResults$Volume.mm3.~StereologyResults$Diagnosis*StereologyResults$Region, ylab="Volume (mm3)", lwd=1, cex.lab=1.3, las=2)
stripchart(StereologyResults$Volume.mm3.~StereologyResults$Diagnosis*StereologyResults$Region, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 20, cex=1, col = c(3,2))
legend(x=7, y=400, legend=c("CTRL", "MDD"), col=c("green","red"), pch=20)
dev.off()





# 
# 
# #I wonder how MDD would compare to other diagnoses
# 
# 


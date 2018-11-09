##Generalized code for processing Freeze 3 Illumina microarray data:
##By Peter Blandino and  Megan Hagenauer
##First attempt: Median Centered Data for 6 regions (10/21/2014)
##Generalized: 10/21/2014...............

library(gdata)
library(fields)
library(stats)
library(car)
library(affy)
library(preprocessCore)
library(multtest)
library(gtools)
library(doBy)
library(outliers)
library(ggplot2)



##R version 3.1.1 (2014-07-10) -- "Sock it to Me"
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~Objects Generated and short description~~~~~~~~~~~~~~~~~~~~~~~~~~~

							"Original Array Data for 6 Regions"
##1		MedianCenteredData----Median Centered Data from Jun
##2		MedianCenteredDataMatrix--File denotes the above file except as a matrix

								"Chip IDs/Disease/Region"
##3		ChipIDs--this file has the subject connected to chip port/location, disease and region
##4		SampleIDsNoOutliersTF--Has the IDs for only samples in which there is chip data. Outliers removed based on Jun's notation of remove from the ChipID file.
##5		ChipIDsNoOutliers--Keeps all the ChipIDs that are true when the SampleIDs match the Chip sample IDs

				"Parceling Out Good Subjects by Region--Remove vs Keep"
##6		AMYChipIDNoOutliers --Pulled data from ChipIDsNoOutliers only for Amy
##7		ANCGChipIDNoOutliers--Pulled data from ChipIDsNoOutliers only for ANCG
##8		CBChipIDNoOutliers--Pulled data from ChipIDsNoOutliers only for CB
##9		DlpfcChipIDNoOutliers--Pulled data from ChipIDsNoOutliers only for DLPFC
##10	HCChipIDNoOutliers--Pulled data from ChipIDsNoOutliers only for HC
##11	NaccChipIDNoOutliers--Pulled data from ChipIDsNoOutliers only for Nacc

		"Parceling Out Median Centered Data by Region From Original Array Data 6 Regions"

##12	AMYMedianCenteredData--Pulled data from MedianCenteredDataMatrix only for Amy
##13	ANCGMedianCenteredData--Pulled data from MedianCenteredDataMatrix only for ANCG
##14	CBMedianCenteredData--Pulled data from MedianCenteredDataMatrix only for CB
##15	DlpfcMedianCenteredData--Pulled data from MedianCenteredDataMatrix only for Dlpfc
##16	HCMedianCenteredData--Pulled data from MedianCenteredDataMatrix only for HC
##17	NaccMedianCenteredData--Pulled data from MedianCenteredDataMatrix only for Nacc

					"Detectable P-values All Regions in One Object"

##18	DetectablePval is the DetectableP Table Updated.csv read in.
##19	DetectablePvalsorted is the DetectablePval file Sorted

			"Parceling Out Detectable P-values by Region--All 10 Regions"

##20	ANCGDetPval is the ANCG pulled out from the DetectablePvalsorted
##21	AmyDetPval is the Amy pulled out from the DetectablePvalsorted
##22	aThalDetPval is the aThal pulled out from the DetectablePvalsorted
##23	CBDetPval is the CB pulled out from the DetectablePvalsorted
##24	DlpfcDetPval is the Dlpfc pulled out from the DetectablePvalsorted
##25	HCDetPval is the HC pulled out from the DetectablePvalsorted
##26	mThalDetPval is the mThal pulled out from the DetectablePvalsorted
##27	NaccDetPval is the Nacc pulled out from the DetectablePvalsorted
##28	PCgDetPval is the PCg pulled out from the DetectablePvalsorted
##29	SCgDetPval is the SCg pulled out from the DetectablePvalsorted

									"Illumina Decode"
##30	IlluminaDecode is the HumanRef-8_V2_0_R4_11223162_A Master Decode.csv file

								"Probe List not Further Used"
##31	OurProbes is created from AMYMedianCenteredData row names.  Note these are a string.
##32	OurProbesSorted is OurProbes sorted.  This is a string.

			"Median Centered Data 2 for 6 Regions Ordered by Probes as.character
					--Probe Names are the Row Names for These Objects"
##33	AMYMedianCenteredData2 is the AMYMedianCenteredData ordered by row names (i.e. probes)
##34	CBMedianCenteredData2 is the CBMedianCenteredData ordered by row names (i.e. probes)
##35	DlpfcMedianCenteredData2 is DlpfcMedianCenteredData ordered by row names (i.e. probes)
##36	HCMedianCenteredData2 is the HCMedianCenteredData ordered by row names (i.e. probes)
##37	NaccMedianCenteredData2 is the NaccMedianCenteredData2 orded by row names (i.e. probes)
##38	ANCGMedianCenteredData2 was not accounted for and is ordered by row names (i.e. probes)
 
"P Value Objects by Region Ordered by Probe as a String, Column 2 in These Objects--10 Regions"
##39	HCDetPvalProbeOrdered is a result of changing HCDetPval data into a string.
##40	ACgDetPvalProbeOrdered is a result of changing ANCGDetPval data into a string.
##41	AmyDetPvalProbeOrdered is a result of changing AmyDetPval data into a string.
##42	aThalDetPvalProbeOrdered is a result of changing aThalDetPval data into a string.
##43	CBDetPvalProbeOrdered is a result of changing CBDetPval data into a string.
##44	DlpfcDetPvalProbeOrdered is a result of changing DlpfcDetPval data into a string.
##45	mThalDetPvalProbeOrdered is a result of changing mThalDetPval data into a string.
##46	NaccDetPvalProbeOrdered is a result of changing NaccDetPval data into a string.
##47	PCgDetPvalProbeOrdered is a result of changing PCgDetPval data into a string.
##48	SCgDetPvalProbeOrdered is a result of changing SCgDetPval data into a string.

			"Temp Files 1-5 Determining if Region_Median Centered Data 2 
					match Column 2 of Region_DetPvalProbeOrdered"
##49	Temp1 determines if probe names match up between ACgDetPvalProbeOrdered and ANCGMedianCenteredData2.  The object has "TRUE" values.
##50	Temp2 determines if probe names match up between AmyDetPvalProbeOrdered and AMYMedianCenteredData2.  The object has "TRUE" values.
##51	Temp3 determines if probe names match up between CBDetPvalProbeOrdered and CBMedianCenteredData2.  The object has "TRUE" values.
##52	Temp4 determines if probe names match up between DlpfcDetPvalProbeOrdered and DlpfcMedianCenteredData2.  The object has "TRUE" values.
##53	Temp5 determines if probe names match up between NaccDetPvalProbeOrdered and NaccMedianCenteredData2.  The object has "TRUE" values.

		"Illumina Decode 22177 Probes Ordered by Array_Address as a String"
		--NOTE:Array Address is equivalent to probe in all other objects
##54	IlluminaDecodeSorted is by making IlluminaDecode Array_Address_Id (i.e. probes, Column 15) a string and ordering it by that string.
##55	IlluminaDecodeSortedTF checks to see if IlluminaDecodeSorted Array_Address_Id matches the probe names of AMYMedianCenteredData2
##56	IlluminaDecodeSorted22177 is all the "True" from the IlluminaDecodeSorted--Removed the extra probes not on the chip

										"Experiment by Megan"
##57	HCMedianOfMedianCenteredData was median centered from HCMedianCenteredData2 however, the data was already median centered thus not useful.

		"Parceling Out Array Data w/Detectable P-Value < 0.055---For Each 6 Regions"
##58 HCMedianCenteredDataDetected is a result of taking a p value < 0.055 from column 3 of HCDetPvalProbeOrdered matched against HCMedianCenteredData2.  Rememeber these are in the order of the probes as a string.
##59	ACgMedianCenteredDataDetected is a result of taking a p value < 0.055 from column 3 of ACgDetPvalProbeOrdered matched against ACgMedianCenteredData2.  Rememeber these are in the order of the probes as a string.
##60	AmyMedianCenteredDataDetected is a result of taking a p value < 0.055 from column 3 of AmyDetPvalProbeOrdered matched against AMYMedianCenteredData2.  Rememeber these are in the order of the probes as a string.
##61	CBMedianCenteredDataDetected<-CBMedianCenteredData2[CBDetPvalProbeOrdered[,3]<0.055,] creates a new object which contains array data only when column 3 of the previous object had a detectable p-value of less than 0.055
##62	DlpfcMedianCenteredDataDetected<-DlpfcMedianCenteredData2[DlpfcDetPvalProbeOrdered[,3]<0.055,] creates a new object which is a subset of the Dlpfc array data that has a detectable p value from the object DlpfcDetPvalProbeOrdered 
##63	NaccMedianCenteredDataDetected<-NaccMedianCenteredData2[NaccDetPvalProbeOrdered[,3]<0.055,] creates a new object as a subset of array data from NaccMedianCenteredData2 that has a detectable p value less than 0.055 as criteria of the object NaccDetPvalProbeOrdered.

					"Detectable P-value <0.055 by Region--6 Regions"
##64	HCDetPvalDetected is the result of taking HCDetPvalProbeOrdered[HCDetPvalProbeOrdered[,3]<0.055,]
##65	ACgDetPvalDetected is the result of filtering ACgDetPvalProbeOrdered column 3 by p<0.055 and taking the remaining columns and rows and saving.
##66	AmyDetPvalDetected is the result of filtering AmyDetPvalProbeOrdered column 3 by p<0.055 and taking the remaining columns and rows and saving.
##67	CBDetPvalDetected<-CBDetPvalProbeOrdered[CBDetPvalProbeOrdered[,3]<0.055,] creates a new object with only probes that have a detectable pvalue of less than 0.055
##68	DlpfcDetPvalDetected<-DlpfcDetPvalProbeOrdered[DlpfcDetPvalProbeOrdered[,3]<0.055,] creates a sub list of probes wtih a detectable p value less that 0.055
##69	NaccDetPvalDetected<-NaccDetPvalProbeOrdered[NaccDetPvalProbeOrdered[,3]<0.055,] places only probes with a detectable p value of less than 0.055 from column 3 of the larger Nacc detecable p value object

			"Illumina Decode 22177 Relevant to Each Region w/Det P-value < 0.055"
	 **Abandoned because the far Right Columns are too large for Exproting into Excel**
##70	HCdetectedIlluminaDecodeSorted22177 is a file that combines IlluminaDecodeSorted22177 and HCDetPvalProbeOrdered with a pvalue <0.055
##71	ACgdetecedIlluminaDecodeSorted22177 pulls out the probes with a p value less than 0.055 from ACgDetPvalProbeOrdered and also matches with IlluminaDecodeSorted22177.  NOTE:  This has the REALLY large columns and is redone below.
##72	AmydetectedIlluminaDecodeSorted22177 pulls out the probes with a p value less than 0.055 from AmyDetPvalProbeOrdered and also matches with IlluminaDecodeSorted22177.  NOTE:  This has the REALLY large columns and is redone below.

	"Problematic Because the Content of the Last Few Columns is too Large for Excel"
##73	HC_Annotated_DetectedProbeData is a combination of the 3 files above which are: HCMedianCenteredDataDetected (median centered data that has a detectable p-value less than 0.055), HCDetPvalDetected (probes with a detected p-value less than 0.055) and HCdetectedIlluminaDecodeSorted22177 (which has all the probe information).  NOTE:  THIS FILE HAS INFORMATION TOO LARGE TO FIT INTO A CELL IN EXCEL THUS WILL INDUCE EXTRA ROWS IN EXCEL AS COMPARED TO DIMENSIONS SEEN IN R.  FILES BELOW WILL FIX THIS ISSUE BY REMOVING THE COLUMNS WITH THE ENORMOUS INFORMATION, SPECIFICALLY ONTOLOGICAL_COMPONENT, ONTOLOGICAL_PROCESS, AND ONTOLOGICAL_FUNCTION and two others.
##74	ACg_Annotated_detectedProbeData is a combination of the 3 files above which are: ACgMedianCenteredDataDetected (median centered data that has a detectable p-value less than 0.055), ACgDetPvalDetected (probes with a detected p-value less than 0.055) and ACgdetectedIlluminaDecodeSorted22177 (which has all the probe information).  NOTE:  THIS FILE HAS INFORMATION TOO LARGE TO FIT INTO A CELL IN EXCEL THUS WILL INDUCE EXTRA ROWS IN EXCEL AS COMPARED TO DIMENSIONS SEEN IN R.  FILES BELOW WILL FIX THIS ISSUE BY REMOVING THE COLUMNS WITH THE ENORMOUS INFORMATION, SPECIFICALLY ONTOLOGICAL_COMPONENT, ONTOLOGICAL_PROCESS, AND ONTOLOGICAL_FUNCTION and two others.
##75	Amy_Annotated_DetectedProbeData is a combination of the 3 files above which are: AmyMedianCenteredDataDetected (median centered data that has a detectable p-value less than 0.055), AmyDetPvalDetected (probes with a detected p-value less than 0.055) and AmydetectedIlluminaDecodeSorted22177 (which has all the probe information).  NOTE:  THIS FILE HAS INFORMATION TOO LARGE TO FIT INTO A CELL IN EXCEL THUS WILL INDUCE EXTRA ROWS IN EXCEL AS COMPARED TO DIMENSIONS SEEN IN R.  FILES BELOW WILL FIX THIS ISSUE BY REMOVING THE COLUMNS WITH THE ENORMOUS INFORMATION, SPECIFICALLY ONTOLOGICAL_COMPONENT, ONTOLOGICAL_PROCESS, AND ONTOLOGICAL_FUNCTION and two others.

			"Failed Attempts to Remove the Last 5 Columns From Annotation"
##76	ACg_Annotated_DetectedProbeDataMinus is ACg_Annotated_DetectedProbeData with the last 5 columns removed.  Not sure if this is really correct and will redo ACg with other regions below.
##77	ACg_Annotated_DetectedProbeDataAbrev is the same as ACg_Annotated_DetectedProbeDataMinus
##78	Amytemp is AmydetectedIlluminaDecodeSorted22177 with the last 5 columns removed. But this is not correct.

						"Illumina Decode with Last 5 Columns Removed"
##79	IlluminaDecodeSorted22177Abbreviated<-IlluminaDecodeSorted22177[1:23]- this removes the problematic columns from the Illumina decode and allows it to be applied to all data files.

						"Abbreviated Illumina Decode by Region-6 Regions--
			The Illumina Probes are in the Same Order as the Region_DetPvalProbeOrdered"
##80	ACgdetectedIlluminaDecodeSorted22177Abbreviated<-IlluminaDecodeSorted22177Abbreviated[ACgDetPvalProbeOrdered[,3]<0.055,] same as above but uses an object with an abbreviated Illumina Decode (i.e. has the last 5 columns removed)
##81	HCdetectedIlluminaDecodeSorted22177Abbreviated<-IlluminaDecodeSorted22177Abbreviated[HCDetPvalProbeOrdered[,3]<0.055,]  creates a new object by taking all rows in HCDetPvalProbeOrdered column 3 with the criteria of 0.055
##82	AmydetectedIlluminaDecodeSorted22177Abbreviated<-IlluminaDecodeSorted22177Abbreviated[AmyDetPvalProbeOrdered[,3]<0.055,] creates the object AmydetectedIlluminaDecodeSorted22177Abbreviated which is the abbreviated Illumina decode with only the probes that had a detectable pvalue of less than 0.055
##83	CBIlluminaDecodeSorted22177Abbreviated<-IlluminaDecodeSorted22177Abbreviated[CBDetPvalProbeOrdered[,3]<0.055,] created a new object which has Illumina decode probes that meet the criteria in column 3 of CBDetPvalProbeOrdered as determined across each row.
##84	DlpfcIlluminaDecodeSorted22177Abbreviated<-IlluminaDecodeSorted22177Abbreviated[DlpfcDetPvalProbeOrdered[,3]<0.055,] creates the object DlpfcIlluminaDecodeSorted22177Abbreviated which results from the Illumina decode object filtered by the Dlpfc detectable p-values. 
##85	NaccIlluminaDecodeSorted22177Abbreviated<-IlluminaDecodeSorted22177Abbreviated[NaccDetPvalProbeOrdered[,3]<0.055,] creates an object representing the probes of the Nacc with an detectable p value of less than 0.055 and pulls them from the Illumina decode that has the last 5 columns removed since those were too large to fit into excel.

				"Combining Regional Data, Det p-values and Illumina Decode Objects" 
			*NOTE:  The linking factor for all Objects are the Probe IDs(Row names for 							Region_MedianCenteredDataDetected, Column 2 of Region_DetPvalDetected 
				and Column 15 of Region_detectedIlluminaDecodeSorted22177Abbreviated* 
##86	ACg_Annotated_DetectedProbeData<-cbind(ACgMedianCenteredDataDetected,ACgDetPvalDetected,ACgdetectedIlluminaDecodeSorted22177Abbreviated) used cbind to combine ACgMedianCenteredDataDetected,ACgDetPvalDetected and ACgdetectedIlluminaDecodeSorted22177Abbreviated into a single object by columns maintaining all rows
##87	HC_Annotated_DetectedProbeData<-cbind(HCDetPvalDetected,HCMedianCenteredDataDetected, HCdetectedIlluminaDecodeSorted22177Abbreviated) created the object HC_Annotated_DetectedProbeData by combining 3 different files.
##88	Amy_Annotated_DetectedProbeData<-cbind(AmyMedianCenteredDataDetected,AmyDetPvalDetected,AmydetectedIlluminaDecodeSorted22177Abbreviated) creates a new object combining 3 amygdala files that contain detectable pvalues, Illumina decode and the array data.
##89	CB_Annotated_DectededProbeData<-cbind(CBMedianCenteredDataDetected,CBDetPvalDetected,CBIlluminaDecodeSorted22177Abbreviated) creates a new object combining the columns of all 3 files.
##90	Dlpfc_Annotated_DetectedProbeData<-cbind(DlpfcMedianCenteredDataDetected,DlpfcDetPvalDetected,DlpfcIlluminaDecodeSorted22177Abbreviated) creates a new object by combining 3 files.  This new object has Array data for Dlpfc, that has a detectable p value of 0.055 or less and annotated with the Illumina Decode.
##91	Nacc_Annotated_DetectedProbeData<-cbind(NaccMedianCenteredDataDetected,NaccDetPvalDetected,NaccIlluminaDecodeSorted22177Abbreviated) creates a new object with annotated Nacc array data with a detectable p value less than 0.055.

								"Subject Demographic Data"
##92	SubjectDemographics<-read.csv("All_Subject_Information_for_R.csv", header=T) created new object with all subject information.  SubjectDemographics has basic information I believe to be relevant for the overall analyses and my neuroimmune analyses.  DOES NOT COVER EVERYTHING FROM DAVID WALSH.
##93	SubjectDemographicsOrdered<-SubjectDemographics[order(SubjectDemographics$Sample.ID),] created a new object by ordering the column subject demographics by sample ID not subject ID.  Thus, the subject IDs are ordered within this column and made into a new object.
##94	SubjectDemographicsOrdered is the object in #98 ordered by Sample.ID
##95	SubjectDemographicsOrderedTF created a True/False object by matching all rows of column 4 of ChipIDsNoOutliers with all rows of column 6 of SubjectDemographicsOrdered and assigning a "True" value when the identifiers matched and "False" when they did not.
##96	SubjectDemographicsOrdered644<-SubjectDemographicsOrderedTF SubjectDemographicsOrdered644 takes SubjectDemographicsOrdered and kept all the "TRUE" values from SubjectDemographicsOrderedTF as criteria.  However, this is not in order by Region.
##97	SubjectDemographicsNoOutliers<-SubjectDemographicsOrdered[SampleIDsNoOutliersTF==TRUE,] is an object created by taking all "TRUE" values from SampleIDsNoOutliersTF and all columns of SubjectDemographicsOrdered and placing into a new object.  This object however still had "remove" subjects "outliers" included.

					"Subject Demographics Parceled by Region-6 Regions"
##98	AmySubjectsDemographicsOrdered644<-SubjectsDemographicsOrdered644[c(1:98),] created a new object by pulling rows 1 to 98 from SubjectsDemographicsOrdered644 as Amy.  This however did not work when tabling.  It only pulled 18 instead of 98 meaning SubjectsDemographicsOrdered644 is not ordered by region.
##99	AMYSubjectsDemographicsOrdered<-SubjectsDemographicsOrdered644[SubjectsDemographicsOrdered644[,2]=="AMY", ] created the object wanted by pulling all rows in column two of SubjectsDemographicsOrdered644 that met the criteria of AMY.  HOWEVER!!!!!!!!!!  The Sample IDs did not match in order between AMYSubjectsDemographicsOrdered rows with colnames of AMYMedianCenteredData2.  	
##100	CBSubjectsDemographicsOrdered<-SubjectsDemographicsOrdered644[SubjectsDemographicsOrdered644[,2]=="CB", ] created a new object by taking all rows where "CB" was noted in column 2 of SubjectsDemographicsOrdered644	
##101	DlpfcSubjectsDemographicsOrdered<-SubjectsDemographicsOrdered644[SubjectsDemographicsOrdered644[,2]=="DLPFC", ] creates a new object of subject demographics  when rows withing column 2 of SubjectsDemographicsOrdered644 equals "DLPFC"	
##102	HCSubjectsDemographicsOrder<-SubjectsDemographicsOrdered644[SubjectsDemographicsOrdered644[,2]=="HC",] creates an object by pulling out all columns when the rows in column 2 equal "HC"
##103	NaccSubjectsDemographicsOrdered
##104	ACgSubjectsDemographicsOrdered

			"Region_MedianCenteredData3 Ordered by Columns as Characters"
					*NOTE: Columns are Named After the Chip Position*
##105	AMYMedianCenteredData3<-AMYMedianCenteredData2[ , order(as.character(colnames(AMYMedianCenteredData2)))] created a new object by taking all rows and ordering the column names as characters.  This worked however, it contains 22177 rows of data, thus have all probes without detectable pvalue accounted for.
##106	CBMedianCenteredData3<-CBMedianCenteredData2[ , order(as.character(colnames(CBMedianCenteredData2)))] created a new object by taking column names from CBMedianCenteredData2, treated them as characters and ordered them so they would match CBSubjectsDemographicsOrdered created 2 in #106.  NOTE: this contains all data regardless of detectable p values.
##107	DlpfcMedianCenteredData3<-DlpfcMedianCenteredData2[ , order(as.character(colnames(DlpfcMedianCenteredData2)))] creates a new object which simply saves DlpfcMedianCenteredData2 order by column names when treated as characters.
##108	HCMedianCenteredData3
##109	ACgMedianCenteredData3
##110	NaccMedianCenteredData3

			"Region_MedianCenteredDataDetected2 p<0.055 and col ordered as characters"
##111	HCMedianCenteredDataDetected2
##112	ACgMedianCenteredDataDetected2
##113	AmyMedianCenteredDataDetected2 
##114	CBMedianCenteredDataDetected2
##115	DlpfcMedianCenteredDataDetected2  
##116	NaccMedianCenteredDataDetected2

		"Region_SubjectsDemographicsOrdered2 copy of Region_SubjectsDemographicsOrdered"
##	CBSubjectsDemographicsOrdered2<-CBSubjectsDemographicsOrdered
##	DlpfcSubjectsDemographicsOrdered2<-DlpfcSubjectsDemographicsOrdered
##	NaccSubjectsDemographicsOrdered2<-NaccSubjectsDemographicsOrdered
##	AMYSubjectsDemographicsOrdered2<-AMYSubjectsDemographicsOrdered
##	HCSubjectsDemographicsOrdered2<-HCSubjectsDemographicsOrdered
##	ACgSubjectsDemographicsOrdered2<-ACgSubjectsDemographicsOrdered

		"Region_MedianCenteredData4 renaming of Region_MedianCenteredData3"   							*Region_MedianCenteredData4, column names reassigned by subject ID 
		from Region_SubjectsDemographicsOrdered2. Each Contains 22177 probes*
##	CBMedianCenteredData4<-CBMedianCenteredData3
##	DlpfcMedianCenteredData4<-DlpfcMedianCenteredData3
##	NaccMedianCenteredData4<-NaccMedianCenteredData3
##	AMYMedianCenteredData4<-AMYMedianCenteredData3
##	HCMedianCenteredData4<-HCMedianCenteredData3
##	ACgMedianCenteredData4<-ACgMedianCenteredData3

	"Region_MedianCenteredDataDetected3 renaming of Region_MedianCenteredDataDetected2"
	 *Region_MedianCenteredDataDetected3 column names reassigned by subject IDs 
	 from Region_SubjectsDemographicsOrdered2.  Contains only detectable pvalue probes*
##	CBMedianCenteredDataDetected3<-CBMedianCenteredDataDetected2
##	DlpfcMedianCenteredDataDetected3<-DlpfcMedianCenteredDataDetected2
##	NaccMedianCenteredDataDetected3<-NaccMedianCenteredDataDetected2
##	AmyMedianCenteredDataDetected3<-AmyMedianCenteredDataDetected2
##	HCMedianCenteredDataDetected3<-HCMedianCenteredDataDetected2
##	ACgMedianCenteredDataDetected3<-ACgMedianCenteredDataDetected2

		"Region_SubjectsDemographicsOrdered3 ordered by subject_id 
				from Region_SubjectsDemographicsOrdered2"
##	CBSubjectsDemographicsOrdered3<-CBSubjectsDemographicsOrdered2[order(as.character(CBSubjectsDemographicsOrdered2[,6])),]
##	DlpfcSubjectsDemographicsOrdered3<-DlpfcSubjectsDemographicsOrdered2[order(as.character(DlpfcSubjectsDemographicsOrdered2[,6])),]
##	NaccSubjectsDemographicsOrdered3<-NaccSubjectsDemographicsOrdered2[order(as.character(NaccSubjectsDemographicsOrdered2[,6])),]
##	AmySubjectsDemographicsOrdered3<-AMYSubjectsDemographicsOrdered2[order(as.character(AMYSubjectsDemographicsOrdered2[,6])),]
##	HCSubjectsDemographicsOrdered3<-HCSubjectsDemographicsOrdered2[order(as.character(HCSubjectsDemographicsOrdered2[,6])),]
##	ACgSubjectsDemographicsOrdered3<-ACgSubjectsDemographicsOrdered2[order(as.character(ACgSubjectsDemographicsOrdered2[,6])),]

							"Region_MedianCenteredDataDetected4 
	is Region_MedianCenteredDataDetected3 ordered by column names i.e. subject ids"
		*NOTE:This is only Data with a Significant Detectable P value< 0.055*
##	CBMedianCenteredDataDetected4<-CBMedianCenteredDataDetected3[,order(colnames(CBMedianCenteredDataDetected3))]
##	DlpfcMedianCenteredDataDetected4<-DlpfcMedianCenteredDataDetected3[,order(colnames(DlpfcMedianCenteredDataDetected3))]
##	NaccMedianCenteredDataDetected4<-NaccMedianCenteredDataDetected3[,order(colnames(NaccMedianCenteredDataDetected3))]
##	AmyMedianCenteredDataDetected4<-AmyMedianCenteredDataDetected3[,order(colnames(AmyMedianCenteredDataDetected3))]
##	HCMedianCenteredDataDetected4<-HCMedianCenteredDataDetected3[,order(colnames(HCMedianCenteredDataDetected3))]
##	ACgMedianCenteredDataDetected4<-ACgMedianCenteredDataDetected3[,order(colnames(ACgMedianCenteredDataDetected3))]

									"Region_MedianCenteredData5 
			is Region_MedianCenteredData4 ordered by column names i.e. subject ids"
			  *NOTE: This is all the Data Irregardless of Detectable P-value*
##	CBMedianCenteredData5<-CBMedianCenteredData4[,order(colnames(CBMedianCenteredData4))]
##	DlpfcMedianCenteredData5<-DlpfcMedianCenteredData4[,order(colnames(DlpfcMedianCenteredData4))]
##	NaccMedianCenteredData5<-NaccMedianCenteredData4[,order(colnames(NaccMedianCenteredData4))]
##	AmyMedianCenteredData5<-AMYMedianCenteredData4[,order(colnames(AMYMedianCenteredData4))]
##	HCMedianCenteredData5<-HCMedianCenteredData4[,order(colnames(HCMedianCenteredData4))]
##	ACgMedianCenteredData5<-ACgMedianCenteredData4[,order(colnames(ACgMedianCenteredData4))]


"4 more regions readining in of microarray data corrected, normalized and outliers removed"  				*NOTE:  These also have the column name changed to Subject IDs and the probes 
						set as row names. *These did not work although created*
##	athal74readin is the athal74_corrected.txt microarray data 
##	athal74readin2 reiteration of athal74readin with rownames assigned
#	athal74readin3 has the coumn names adjusted to reflect just Subject IDs---data.frame
##	colnamesathal74readin3 just the column names of athal74readin3
##	mthal77readin is the mthal77_corrected3.txt microarray data
##	mthal77readin2 has column names adjusted to reflect Subject IDs
##	mthal77readin3 is just an updated name to keep consistent with other regions---data.frame
##	PCg85readin is the pcg85_corrected2.txt microarray data
##	PCg85readin2 has the column names adjusted to reflect the Subject IDs
##	PCg85readin3 is just an updated name to keep consistent with other regions-----data.frame
##	SCg86readin is the scg86_corrected_w_subject_IDs.txt microarray data
##	SCg86readin2 has additional columns removed originally, then written over by the row names of SCg86readin2
##	SCg86readin3 is an updated version of SCg86readin2 with row names and Subject IDs as column names with "X" removed-------data.frame
##	TestofgSub failed attempt at removing X from subject IDs of aThal74readin3
##	colnamesaThal74readin3 is just he column names
##	columnnames is a temporary file for pulling out column names from 1 of the 4 regions
##	columnnames1 is a temporary file to substitute the X in subject ID with a blank space
##	SCg86readin3rownamesfixed is the first column of SCg86readin3, which should be rownames
##	SCg86readin3a is an iteration of SCg86readin3rownamesfixed
## scg86_sampleIDs is read in data from SCg86_sampleEdittedfromSCg88.txt
##	scg86_sampleIDs2 iteration of scg86_sampleIDs
##	Temp6 examines if the row names of SCg86readin3a are equal to the rows of column two of SCgDetPvalProbeOrdered
##	SCg86readin2a is SCg86readin2 with row names ordered as characters, then re-written using row names as characters from SCg86readin3a
##	SCgDetPvalProbeOrdered has probes ordered as character from SCgDetPval
##	SCg86readin3b is an object created by reading SCg86readin3a.csv as a matrix
##	SCg86readin3c is an iteration of SCg86readin3b, however, row names are taken from column 1 of SCg86readin3b and ordered as numeric
##	SCg86readin3d is SCg86readin3c as a numeric matrix
##	SCgpvalprobesdatarownamesTF object comparing column 2 of SCgDetPvalProbeOrdered and the row.names of SCg86readin3c
##	junk iteration of SCgDetPvalProbeOrdered
##	Temp3 is a True/False object comparing CBDetPvalProbeOrdered column 2 and row.names  of CBMedianCenteredData2
##	CBDetPvalDetected contains only detectable p values < 0.055 from CBDetPvalProbeOrdered
##	PCgpvalprobesdatarownamesTF is a T/F object comparing column 2 of PCgDetPvalProbeOrdered to row.names(PCg85readin3)
##	PCg85readin4 iteration of PCg85readin3 with the row names ordered
##	SCg86readin4 iteration of SCg86readin3  with the row names ordered
##	mThalpvalprobesdatarownamesTF is a T/F object comparing column 2 of PCgDetPvalProbeOrdered to the row.names of mthal77readin3
## mthal77readin4 iteration of mthal77readin4  with the row names ordered
##	mThalpvalprobesdatarownamesTF2 is a T/F object comparing row.names of mthal77readin4 to column 2 of mThalDetPvalProbeOrdered
##	Junk3 combines row.names of mthal77readin4 with coumn 2 of mThalDetPvalProbeOrdered
##	Junk4 is <-mthal77readin4 with the row names ordered as characters
##	Junk5 is mthal77readin4 with the row names ordered as numeric 

								"4 More Regions Correctly Read IN"
##	athal74DataImported read in from athal74_corrected.txt and has 2 statistical columns at the end of the file
##	athal74DataImported2 iteration of athal74DataImported and has the last 2 statistical columns removed from the data and the row names are ordered.
##	aThalpvalprobesdatarownamesTF a T/F object comparing aThalDetPvalProbeOrdered column 2 and the row.names of athal74DataImported2
##	Junk9 has column 2 of aThalDetPvalProbeOrdered treated as numeric
##	Junk10 has the row names treated as numeric from athal74DataImported2))
##	mthal77DataImported read in from mthal77_corrected3.txt and has 1 statistical column at the end
##	mthal77DataImported2 iteration of mthal77DataImported with the statistical column removed and ordered by row names
##	mThalpvalprobesdatarownamesTF T/F object comparing mThalDetPvalProbeOrdered column 2 to the row.names in mthal77DataImported2
##	Junk7 has column 2 of mThalDetPvalProbeOrdered treated as numeric
##	Junk8 has the row names treated as numeric from mthal77DataImported2
##	SCg86DataImported read in from table scg86_corrected_w_subject_IDs.txt with 1 statistical column at the end
##	SCg86DataImported2 iteration of SCg86DataImported with the last statistical column removed to retain only data with row names ordered
##	SCgpvalprobesdatarownamesTF T/F object comparing columns of SCgDetPvalProbeOrdered with row.names in SCg86DataImported2
##	Junk11 has column 2 of SCgDetPvalProbeOrdered treated as numeric
##	Junk12 has the row names treated as numeric SCg86DataImported2
##__Remember the PCg was read in correctly however to make sure the nomenclature stays the same some new objects were created.

				"Removal of "X" from Subject IDs for 4 More Regions"

##	PCg85DataImported2 iteration of PCg85readin4 just to get the name correct
##	SCg86DataImported3 has the column names (i.e. Subject IDs) without X infront
##	mthal77DataImported3 has the column names (i.e. Subject IDs) without X infront
##	athal74DataImported3 has the column names (i.e. Subject IDs) without X infront
##	PCg85DataImported3 has the column names (i.e. Subject IDs) without X infront

		"4 More Regions-Median Centered Data Filtered by Detectable P Value < 0.055"
##	SCg86DataImported3MedCtrDataDetected<-SCg86DataImported3[SCgDetPvalProbeOrdered[,3]<0.055,]
##	PCg85DataImported3MedCtrDataDetected
##	aThal74DataImported3MedCtrDataDetected
##	mThal77DataImported3MedCtrDataDetected

				"4 More Regions-P Values Significantly Detected & Probes"
##	SCgDetPvalDetected
##	PCgDetPvalDetected
##	aThalDetPvalDetected
##	mThalDetPvalDetected

				"4 More Regions-Abbreviated Illumina Decode Relevant Probes"
##	SCgdetectedIlluminaDecodeSorted22177Abbreviated is Illumina Decode with last 5 columns removed and reduced to contain only probes with sig detection for this region
##	PCgdetectedIlluminaDecodeSorted22177Abbreviated is Illumina Decode with last 5 columns removed and reduced to contain only probes with sig detection for this region
##	aThaldetectedIlluminaDecodeSorted22177Abbreviated is Illumina Decode with last 5 columns removed and reduced to contain only probes with sig detection for this region
##	mThaldetectedIlluminaDecodeSorted22177Abbreviated is Illumina Decode with last 5 columns removed and reduced to contain only probes with sig detection for this region

		"4 More Regions-Detectable Data, Detection P Value < 0.055 & Illumina Annotation"
##	SCg_Annotated_DetectedProbeData is region specific median centered data that is significantly detected with Illumina Annotation
##	PCg_Annotated_DetectedProbeData is region specific median centered data that is significantly detected with Illumina Annotation
##	aThal_Annotated_DetectedProbeData is region specific median centered data that is significantly detected with Illumina Annotation
##	mThal_Annotated_DetectedProbeData is region specific median centered data that is significantly detected with Illumina Annotation

							"Subjects With Illness-Categorical"
##	SubjectsWithIllness Read in from Subjects with Illness.csv
##	SubjectsWithIllnessIDs<-SubjectsWithIllness[,1]
##	SubjectsWithIllnessIDsAsRowNames<-SubjectsWithIllness
##	row.names(SubjectsWithIllnessIDsAsRowNames)<SubjectsWithIllnessIDs
##	In doing this did I just order the one column and not the whole object?

							"Region_MedianCenteredDataDetected5"
##	aThalMedianCenteredDataDetected5 iteration of aThal74DataImported3MedCtrDataDetected
##	mThalMedianCenteredDataDetected5 iteration of mThal77DataImported3MedCtrDataDetected
##	PCgMedianCenteredDataDetected5 iteration of PCg85DataImported3MedCtrDataDetected
##	SCgMedianCenteredDataDetected5 iteration of SCg86DataImported3MedCtrDataDetected
##	HCMedianCenteredDataDetected5 iteration of HCMedianCenteredDataDetected4
##	ACgMedianCenteredDataDetected5 iteration of ACgMedianCenteredDataDetected4
##	AmyMedianCenteredDataDetected5 iteration of AmyMedianCenteredDataDetected4
##	CBMedianCenteredDataDetected5 iteration of CBMedianCenteredDataDetected4
##	DlpfcMedianCenteredDataDetected5 iteration of DlpfcMedianCenteredDataDetected4
##	NaccMedianCenteredDataDetected5 iteration of NaccMedianCenteredDataDetected4

							"Region_Subjects-4 More Regions"
##	aThalSubjects is column names from appropriate Region_MedianCenteredDataDetected5 
##	mThalSubjects is column names from appropriate Region_MedianCenteredDataDetected5
##	PCgSubjects is column names from appropriate Region_MedianCenteredDataDetected5
##	SCgSubjects is column names from appropriate Region_MedianCenteredDataDetected5

								"Region_SubjectIllnesTF"
##	AThalSubjectIllnessTF T/F Object comparing Col 1 SubjectsWithIllness and aThalSubjects
##	mThalSubjectIllnessTF T/F Object comparing Col 1 SubjectsWithIllness and mThalSubjects
##	PCgSubjectIllnessTF T/F Object comparing Col 1 SubjectsWithIllness and PCgSubjects
##	SCgSubjectIllnessTF T/F Object comparing Col 1 SubjectsWithIllness and SCgSubjects



							"Region_SubjectsIllnessCategories"
##	aThalSubjectsIllnessCategories keeps only subjects from SubjectsWithIllness that are "TRUE"
##	mThalSubjectsIllnessCategories keeps only subjects from SubjectsWithIllness that are "TRUE"
##	SCgSubjectsIllnessCategories keeps only subjects from SubjectsWithIllness that are "TRUE"
##	PCgSubjectsIllnessCategories keeps only subjects from SubjectsWithIllness that are "TRUE"
"TRUE" from "Region_SubjectIllnesTF"


Illness 6 more regions
NaccSubjectsIllnessCategories




									"Logic Flowchart"

## Data		Median Centered Data for 6 Regions -> Region_MedianCenteredData, Parcelled out by 6 different Regions (Amy, ACg, CB, Dlpfc, HC, Nacc) -> Region_MedianCentered Data 2, order by probes as characters for 6 different Regions (Amy, ACg, CB, Dlpfc, HC, Nacc)

## Detectable Pval is the read in -> "Region_DetPval" which parcelled out ALL 10 Regioins -> Region_DetPvalProbeOrdered -> Region_DetPvalDetected when p< 0.055 

## Merge	Region_MedianCentered Data 2 and Region_DetPvalProbeOrdered -> Region_MedianCenteredDataDetected


##_______________
##	Illumina Decode-> Array Address ID which is equivalent to Probe ID ordered as character -> IlluminaDecodeSorted22177Abbreviated (has last 5 columns removed)-> 

##	Merge	IlluminaDecodeSorted22177Abbreviated and Region_DetPvalProbeOrdered ->Region_detectedIlluminaDecodeSorted22177Abbreviated (both regional and det pval < 0.055)


##______________Used to Share with Consortium
##	Merge	Region_MedianCenteredDataDetected, Region_DetPvalDetected, Region_detectedIlluminaDecodeSorted22177Abbreviated TO CREATE:
 ***********************Annotated Median Centered Data w/Detection < 0.055 also reported for each of 6 Regions (Amy, ACg, CB, Dlpfc, HC, Nacc) **********************
 :NOW KNOWN AS "Region_Annotated_DetectedProbeData"


##_______________
## ChipIDsNoOutliers-> Parcelled out by 6 different Regioins (Amy, ACg, CB, Dlpfc, HC, Nacc).  This is used to translate Chip_IDs to Subject_IDs.

##	SubjectDemographics (666 Subjects)-> then ordered by Sample ID (SubjectDemographicsOrdered)

##	Merge SubjectDemographicsOrderedTF combines SubjectDemographicsOrdered Column 6 with ChipIDsNoOutliers Column 4 into a T/F object to create SubjectDemographicsOrdered644

##	SubjectDemographicsOrdered644 -> broken out into each region. . . -> Region_SubjectsDemographicsOrdered (Amy, ACg, CB, Dlpfc, HC, Nacc) plus criteria such as pH, agonal factor, etc. -> Region_SubjectDemographicsOrdered2 is simply a copy of the previous object -> Region_SubjectDemographicsOrdered3 is the previous object ordered by Subject IDs--This is for 644 Total Subjects and 6 Regions


##_______________
## Region_MedianCenteredData2 -> Region_MedianCenteredData3 has ordred column names aka Chip_IDs as characters -> Region_MedianCenteredData4 renames columns by Subject_IDs -> Region_MedianCenteredData5 column names ordered by Subject IDs-----Thus, Median Centered Data ordered by Subjects (columns) and has ALL 22177 probes (Amy, ACg, CB, Dlpfc, HC, Nacc)

##	Region_MedianCenteredDataDetected are the Region_MedianCenteredData which have a detectable p<0.055 and probes ordered as a string ->Region_MedianCenteredDataDetected2 has the columns ordered as characters -> Region_MedianCenteredDataDetected3 is the previous object with the columns names reassigned by Subject IDs -> Region_MedianCenteredDataDetected4 has the Subject IDs in order--------Thus, Median Centered Data ordered by Subjects (columns) and has ONLY Probes with a Detectable Pval < 0.055 (Amy, ACg, CB, Dlpfc, HC, Nacc)
*NOTE: Above Object (Region_MedianCenteredDataDetected4) should be used for Statistics*


##_______________4 More Regions
##	Region_DataImported -> trimmed unncessary columns and ordered row names Region_DataImported2 -> removed "X" from column names Region_DataImported3 -> filtered by detectable pval ojects farther above Region_DataImported3MedCtrDataDetected

## Region_DetPvalDetected are just the probes with detection of p<0.055 -same as done farther above but done with the 4 more region data prep

## Region_detectedIlluminaDecodeSorted22177Abbreviated are Illumina Decode with the last 5 columns removed, detectable pval only for each region.

Merge	Region_DataImported3MedCtrDataDetected, Region_DetPvalDetected, Region_detectedIlluminaDecodeSorted22177Abbreviated TO CREATE: ************************ Annotated Median Centered Data w/Detection < 0.055 for 4 Regions (athal, mThal, SCg, PCg) ************************* :NOW KNOWN AS "Region_Annotated_DetectedProbeData"


##_____________Illness Files
##	SubjectsWithIllness -> SubjectsWithIllnessIDs




##	Files to Use for Stats
SCgStats<-SCg86DataImported3MedCtrDataDetected
PCgStats<-PCg85DataImported3MedCtrDataDetected
aThalStats<-aThal74DataImported3MedCtrDataDetected
mThalStats<-mThal77DataImported3MedCtrDataDetected
HCStats<-HCMedianCenteredDataDetected4
AmyStats<-AmyMedianCenteredDataDetected4
ACgStats<-ACgMedianCenteredDataDetected4
CBStats<-CBMedianCenteredDataDetected4
DlpfcStats<-DlpfcMedianCenteredDataDetected4
NaccStats<-NaccMedianCenteredDataDetected4





#********************************USEFUL Summaries************************************************
# sum(HCDetPvalProbeOrdered[,3]<0.055)
#[1] 13334
# sum(ACgDetPvalProbeOrdered[,3]<0.055)
[1] 13408
# sum(AmyDetPvalProbeOrdered[,3]<0.055)
[1] 13319
# sum(aThalDetPvalProbeOrdered[,3]<0.055)
[1] 13771
# sum(CBDetPvalProbeOrdered[,3]<0.055)
[1] 12931
# sum(DlpfcDetPvalProbeOrdered[,3]<0.055)
[1] 13370
# sum(mThalDetPvalProbeOrdered[,3]<0.055)
[1] 13690
# sum(NaccDetPvalProbeOrdered[,3]<0.055)
[1] 13530
# sum(PCgDetPvalProbeOrdered[,3]<0.055)
[1] 13884
# sum(SCgDetPvalProbeOrdered[,3]<0.055)
[1] 13919



#When working across multiple days:
#1) Open R
#2) Set the working directory
#3) Load the work space (i.e. where you will store the files)
#4) open the code document
#4) reload packages---see Code to run for step 2 below

#******************************************************************

#STEP 1: Set the current working directory and workspace file:

##1. Click on the R Console 
##2. Go to File->change dir
##3. File on Desktop C:/Users/blandino/Desktop/Freeze 3 Analysis with Megan
##4.  NOTE: working directory is the place on the computer where a particular file of interest is located.  The directory must change if the file is not in the same location.



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

##8. Load the packages by highlighting this code and hitting run (Command and Enter):

#***CODE TO RUN FOR STEP 2****

library(gdata)
library(fields)
library(stats)
library(car)
library(affy)
library(preprocessCore)
library(multtest)
library(gtools)
# highlight and hit command and enter to run code which is loading the above packages.
#***END OF CODE FOR STEP 2****


#*******************************************************************

#STEP 3: Make sure the Illumina data is ready to be inputted into R:

##1. The data file should be outputted from GenomeStudio using "Export displayed data to a file."
##2. The data is a matrix of the raw gene expression data for each probe (by row) for each sample (by column).
##3. Save the file as "genexp_sample probe profile.csv" and place it in the appropriate folder for your brain region.  This does not apply to my data but good to know.
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
RawData<-as.matrix(read.csv("genexp_sample probe profile.csv", header=T, row.names=1))  #This does not apply to my data but good to know.
MedianCenteredData<-read.delim("medctr_region_chip.txt", header=T, sep = "\t")
#The above line does apply and coding really starts with this line above.  "\t" is used to denote the file is tab delimited.
dim(MedianCenteredData)
#[1] 22177   644   # Dimensions of the file MedianCenteredData
is.numeric(MedianCenteredData)
#[1] FALSE  #File is not numeric, thus probably be a string.
MedianCenteredDataMatrix<-as.matrix(MedianCenteredData)
#File read in as a Matrix
is.numeric(MedianCenteredDataMatrix)
#[1] TRUE  #File is now in a numeric format.

head(MedianCenteredDataMatrix)
# Lists the headers.  NOTE: What appears as the first column (Jun's Probe IDs) is actually the row headers.
        X1692264001_B X1692264001_H X1692264002_B X1692264005_A X1692264005_B
6960451    0.26341304    0.31664048   -0.08795335   -0.19959960   -0.33447153
2600731    0.02274852    0.09482501   -0.06437061   -0.02931029   -0.20257240
2120309    0.00000000   -0.14257006    0.02153450   -0.15255303   -0.01331088

colnames(MedianCenteredDataMatrix)
# Retrieves the column names
row.names(MedianCenteredDataMatrix)[c(1:5)]
#Retrieves the row names--these two functions are important to see how the data is organized.
# NOTE: "c" stands for concantonate which means connect.

ChipIDs<-read.delim("sample_ordered3.txt", header=T, sep = "\t")
dim(ChipIDs)
[1] 666  14
# This gives the dimensions of the file.  666 rows by 14 columns.
head(ChipIDs)

     Sample.ID      Name Region Individual Disease        chip remove order   label
1 X1692264001_B TZP008453    AMY       3281       C X1692264001   keep     1 5/31/07
2 X1692264001_E TZP008449    AMY       3452      BP X1692264001      r     2 5/31/07
3 X1692264001_H TZP008488    AMY       3878       C X1692264001   keep     3 5/31/07
4 X1692264002_B TZP008447    AMY       3426      MD X1692264002   keep     4 5/22/07
5 X1692264005_A TZP008476    AMY       3671      MD X1692264005   keep     5 5/11/07
6 X1692264005_B TZP008499    AMY       4236       C X1692264005   keep     6 5/11/07
      hyb Block chip.order label.order hyb.order
1 7/16/07     8        280         280       307
2 7/16/07     8        280         280       307
3 7/16/07     8        280         280       307
4 7/11/07     7        160         160       168
5  7/9/07     4         70          69        73
6  7/9/07     4         70          69        73


SampleIDsNoOutliersTF<-ChipIDs[,1]%in%colnames(MedianCenteredDataMatrix)
# Takes ChipIDs in the MedianCenteredDataMatrix an merges/aligns these two files and the resultant is titled SampleIDsNoOutliersTF

length(SampleIDsNoOutliersTF)
#[1] 666  Is the number of rows

SampleIDsNoOutliersTF[c(1:5)]
#[1]  TRUE FALSE  TRUE  TRUE  TRUE
# This . . . . . . 

sum(SampleIDsNoOutliersTF)
#[1] 644
#This simple gives a count of all the number of Trues that are in the file.

table(SampleIDsNoOutliersTF, ChipIDs[,7])
# Creates a table with column 7 of ChipIDs as criteria from SampleIDsNoOutliers


SampleIDsNoOutliersTF keep   r
                FALSE    0  22
                TRUE   644   0
                
#So it looks like the "Remove" column denotes outlier samples that were removed while Jun pre-preprocessed the Freeze3 dataset to remove batch effects.

ChipIDsNoOutliers<-ChipIDs[SampleIDsNoOutliersTF==TRUE, ]
dim(ChipIDsNoOutliers)
#[1] 644  14
# This pulls all the "true" and creates a new file.  The dimensions are checked.
head(ChipIDsNoOutliers)

      Sample.ID      Name Region Individual Disease        chip remove order   label
1 X1692264001_B TZP008453    AMY       3281       C X1692264001   keep     1 5/31/07
3 X1692264001_H TZP008488    AMY       3878       C X1692264001   keep     3 5/31/07
4 X1692264002_B TZP008447    AMY       3426      MD X1692264002   keep     4 5/22/07
5 X1692264005_A TZP008476    AMY       3671      MD X1692264005   keep     5 5/11/07
6 X1692264005_B TZP008499    AMY       4236       C X1692264005   keep     6 5/11/07
7 X1692264005_G TZP008435    AMY       3031      MD X1692264005   keep     7 5/11/07
      hyb Block chip.order label.order hyb.order
1 7/16/07     8        280         280       307
3 7/16/07     8        280         280       307
4 7/11/07     7        160         160       168
5  7/9/07     4         70          69        73
6  7/9/07     4         70          69        73
7  7/9/07     4         70          69        73



table(ChipIDsNoOutliers[,3])
#creates a table of column 3 and tells us there are 98 AMY counted in column 3, 109 ANCG counted and so on.
  AMY  ANCG    CB DLPFC    HC  NACC 
   98   109   113   109   109   106 

#Below is done by hand to determine where each region begins and ends so they can be extracted.

ChipIDsNoOutliers
AMY Rows 1:98 
ANCG Rows 99:207
CB Rows 208:320
DLPFC Rows 321:429 
HC Rows 430:538
NACC Rows 539:644 

AMYChipIDNoOutliers<-ChipIDsNoOutliers[c(1:98),] #defined this file
table(AMYChipIDNoOutliers[,3])  #checked to make sure the count is correct
  AMY  ANCG    CB DLPFC    HC  NACC 
   98     0     0     0     0     0 
   
   
ANCGChipIDNoOutliers<-ChipIDsNoOutliers[c(99:207),]  #defined this file
table(ANCGChipIDNoOutliers[,3])  #checked to make sure the count is correct
AMY  ANCG    CB DLPFC    HC  NACC 
    0   109     0     0     0     0 

CBChipIDNoOutliers<-ChipIDsNoOutliers[c(208:320),]  #defined this file
table(CBChipIDNoOutliers[,3])  #checked to make sure the count is correct
  AMY  ANCG    CB DLPFC    HC  NACC 
    0     0   113     0     0     0 

DlpfcChipIDNoOutliers<-ChipIDsNoOutliers[c(321:429),]  #defined this file
table(DlpfcChipIDNoOutliers[,3])  #checked to make sure the count is correct
AMY  ANCG    CB DLPFC    HC  NACC 
    0     0     0   109     0     0 

HCChipIDNoOutliers<-ChipIDsNoOutliers[c(430:538),]  #defined this file
table(HCChipIDNoOutliers[,3])  #checked to make sure the count is correct
AMY  ANCG    CB DLPFC    HC  NACC 
    0     0     0     0   109     0 

NaccChipIDNoOutliers<-ChipIDsNoOutliers[c(539:644),]  #defined this file
table(NaccChipIDNoOutliers[,3])  #checked to make sure the count is correct
AMY  ANCG    CB DLPFC    HC  NACC 
    0     0     0     0     0   106 

#Extracting Median Centered Data from the large Jun file into each separate region
AMYMedianCenteredData<-MedianCenteredDataMatrix[,c(1:98)]
ANCGMedianCenteredData<-MedianCenteredDataMatrix[,c(99:207)]
CBMedianCenteredData<-MedianCenteredDataMatrix[,c(208:320)]
DlpfcMedianCenteredData<-MedianCenteredDataMatrix[,c(321:429)]
HCMedianCenteredData<-MedianCenteredDataMatrix[,c(430:538)]
NaccMedianCenteredData<-MedianCenteredDataMatrix[,c(539:644)]


#Everything after this point needs to be re-run and saved to the workspace because R crashed.  We left off in the middle of splitting up the detectable p-value file by region. After that point, we were going to sort those region-specific detectable p-value files by probe#. We also needed to draw out the probes for our data from the detailed Probe information from Illumina.


#Doublechecking that the ChipID and MedianCentered data files have subjects in the same order:
colnames(AMYMedianCenteredData)==AMYChipIDNoOutliers[,1]
# command verifies that the column names in each file match each other.  The ouput will be used as an example for the other 5 regions below.
#colnames(AMYMedianCenteredData)==AMYChipIDNoOutliers[,1]
 [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
[21] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
[41] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
[61] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
[81] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE

colnames(ANCGMedianCenteredData)==ANCGChipIDNoOutliers[,1]
colnames(CBMedianCenteredData)==CBChipIDNoOutliers[,1]
colnames(DlpfcMedianCenteredData)==DlpfcChipIDNoOutliers[,1]
colnames(HCMedianCenteredData)==HCChipIDNoOutliers[,1]
colnames(NaccMedianCenteredData)==NaccChipIDNoOutliers[,1]

#I changed the working directory to access the detection p-values
DetectablePval<-read.csv("DetectableP Table Updated.csv", header=T)
dim(DetectablePval)
#[1] 221770      3
head(DetectablePval)
#  Region Jun.Probe Detectable.Pval
1    Amy   6040050               0
2    Amy   2640291               0
3    Amy   3130008               0
4    Amy   6180408               0
5    Amy   6560228               0
6    Amy   4390201               0

table(DetectablePval[,1])

  ACg   Amy aThal    CB Dlpfc    HC mThal  Nacc   PCg   SCg 
22177 22177 22177 22177 22177 22177 22177 22177 22177 22177

#ProbeIDs by Region
ACG Rows 1:22177
Amy Rows 22178:44354
aThal Rows 44355:66531
CB Rows 66532:88708
Dlpfc Rows 88709:110885
HC Rows 110886:133062
mThal Rows 133063:155239
Nacc Rows 155240:177416
PCg Rows 177417:199593
SCg Rows 199594:221770

DetectablePvalsorted<-DetectablePval[order(DetectablePval$Region),]
head(DetectablePvalsorted)

      Region Jun.Probe Detectable.Pval
43483    ACg   2600731               0
43484    ACg   7510608               0
43485    ACg   5130373               0
43486    ACg   1030215               0
43487    ACg   6770646               0
43488    ACg   6350386               0

ANCGDetPval<-DetectablePvalsorted[c(1:22177),]

head(ANCGDetPval)
      Region Jun.Probe Detectable.Pval
43483    ACg   2600731               0
43484    ACg   7510608               0
43485    ACg   5130373               0
43486    ACg   1030215               0
43487    ACg   6770646               0
43488    ACg   6350386               0


table(ANCGDetPval[,1])
  ACg   Amy aThal    CB Dlpfc    HC mThal  Nacc   PCg   SCg 
22177     0     0     0     0     0     0     0     0     0 

AmyDetPval<-DetectablePvalsorted[c(22178:44354),]
table(AmyDetPval[,1])
#   ACg   Amy aThal    CB Dlpfc    HC mThal  Nacc   PCg   SCg 
    0 22177     0     0     0     0     0     0     0     0 

aThalDetPval<-DetectablePvalsorted[c(44355:66531),]
table(aThalDetPval[,1])
#ACg   Amy aThal    CB Dlpfc    HC mThal  Nacc   PCg   SCg 
    0     0 22177     0     0     0     0     0     0     0 

CBDetPval<-DetectablePvalsorted[c(66532:88708),]
table(CBDetPval[,1])
#ACg   Amy aThal    CB Dlpfc    HC mThal  Nacc   PCg   SCg 
    0     0     0 22177     0     0     0     0     0     0 

DlpfcDetPval<-DetectablePvalsorted[c(88709:110885),]
table(DlpfcDetPval[,1])
#ACg   Amy aThal    CB Dlpfc    HC mThal  Nacc   PCg   SCg 
    0     0     0     0 22177     0     0     0     0     0 

HCDetPval<-DetectablePvalsorted[c(110886:133062),]
table(HCDetPval[,1])
#ACg   Amy aThal    CB Dlpfc    HC mThal  Nacc   PCg   SCg 
    0     0     0     0     0 22177     0     0     0     0

mThalDetPval<-DetectablePvalsorted[c(133063:155239),]
table(mThalDetPval[,1])
#ACg   Amy aThal    CB Dlpfc    HC mThal  Nacc   PCg   SCg 
    0     0     0     0     0     0 22177     0     0     0 

NaccDetPval<-DetectablePvalsorted[c(155240:177416),]
table(NaccDetPval[,1])
ACg   Amy aThal    CB Dlpfc    HC mThal  Nacc   PCg   SCg 
    0     0     0     0     0     0     0 22177     0     0  

PCgDetPval<-DetectablePvalsorted[c(177417:199593),]
table(PCgDetPval[,1])
#ACg   Amy aThal    CB Dlpfc    HC mThal  Nacc   PCg   SCg 
    0     0     0     0     0     0     0     0 22177     0 
SCgDetPval<-DetectablePvalsorted[c(199594:221770),]
table(SCgDetPval[,1])
#ACg   Amy aThal    CB Dlpfc    HC mThal  Nacc   PCg   SCg 
    0     0     0     0     0     0     0     0     0 22177 
    
#I changede the working directory again
IlluminaDecode<-read.csv("HumanRef-8_V2_0_R4_11223162_A Master Decode.csv", header=T)

dim(IlluminaDecode)
#[1] 22184    28
#The decode has more probes than the chip does!

str(IlluminaDecode)


OurProbes<-row.names(AMYMedianCenteredData)
#Takes the Amy file and just the rows and makes them a file in-and-of itself to act as our list of probes for the data.
length(OurProbes)
[1] 22177
head(OurProbes)
[1] "6960451" "2600731" "2120309" "7510608" "1570494" "6520451"
#Note: Probes are being treated as string variable.  We know this because they are in quotations.


OurProbesSorted<-sort(OurProbes)
head(OurProbesSorted)
[1] "10008" "10017" "10019" "10021" "10025" "10048"


AMYMedianCenteredData2<-AMYMedianCenteredData[order(row.names(AMYMedianCenteredData)), ]
      
(row.names(AMYMedianCenteredData))[c(1:5)]      
      
      X1692264001_B X1692264001_H X1692264002_B X1692264005_A X1692264005_B X1692264005_G
10008    0.24758237  -0.007965018  0.0140427694  -0.406526374     0.1796235    0.10860487
10017    0.00000000   0.031518659 -0.0658551030  -0.177372384     0.1128019    0.05115629
10019    0.56004238  -0.117464051  0.0963222356   0.003268713    -0.1012926   -0.09958782
10021    0.08886211  -0.047386695 -0.1176407854   0.099289904    -0.5693162   -0.06567642
10025    0.38587837   0.098676779 -0.0093789861  -0.115623502     0.1704954    0.14005634
10048   -0.04969145   0.223485571 -0.0004024703  -0.043854675    -0.3421210    0.04346414

ANCGMedianCenteredData2<-ANCGMedianCenteredData[order(row.names(ANCGMedianCenteredData)), ]
CBMedianCenteredData2<-CBMedianCenteredData[order(row.names(CBMedianCenteredData)), ]
DlpfcMedianCenteredData2<-DlpfcMedianCenteredData[order(row.names(DlpfcMedianCenteredData)), ]
HCMedianCenteredData2<-HCMedianCenteredData[order(row.names(HCMedianCenteredData)), ]
NaccMedianCenteredData2<-NaccMedianCenteredData[order(row.names(NaccMedianCenteredData)), ]

is.numeric(DlpfcDetPval[,2]) #yes but the data was a string.  we changed the numeric to a string using the as.character command.  Note Illumina decode is also treating the array address ID as numeric and will need to be changed.

head(HCDetPval)
HCDetPvalProbeOrdered<-HCDetPval[order(as.character(HCDetPval[,2])),]
head(HCDetPvalProbeOrdered)
ACgDetPvalProbeOrdered<-ANCGDetPval[order(as.character(ANCGDetPval[,2])),]
AmyDetPvalProbeOrdered<-AmyDetPval[order(as.character(AmyDetPval[,2])),]
aThalDetPvalProbeOrdered<-aThalDetPval[order(as.character(aThalDetPval[,2])),]
CBDetPvalProbeOrdered<-CBDetPval[order(as.character(CBDetPval[,2])),]
DlpfcDetPvalProbeOrdered<-DlpfcDetPval[order(as.character(DlpfcDetPval[,2])),]
mThalDetPvalProbeOrdered<-mThalDetPval[order(as.character(mThalDetPval[,2])),]
NaccDetPvalProbeOrdered<-NaccDetPval[order(as.character(NaccDetPval[,2])),]
PCgDetPvalProbeOrdered<-PCgDetPval[order(as.character(PCgDetPval[,2])),]
SCgDetPvalProbeOrdered<-SCgDetPval[order(as.character(SCgDetPval[,2])),]

head(HCDetPvalProbeOrdered)
Temp<-(HCDetPvalProbeOrdered[,2]==row.names(HCMedianCenteredData2))
length(Temp)
Temp[c(1:5)]
sum(Temp)

Temp1<-(ACgDetPvalProbeOrdered[,2]==row.names(ANCGMedianCenteredData2))
sum(Temp1)
#[1] 22177  True statements

Temp2<-(AmyDetPvalProbeOrdered[,2]==row.names(AMYMedianCenteredData2))
sum(Temp2)

Temp3<-(CBDetPvalProbeOrdered[,2]==row.names(CBMedianCenteredData2))
sum(Temp3)

Temp4<-(DlpfcDetPvalProbeOrdered[,2]==row.names(DlpfcMedianCenteredData2))
sum(Temp4)

Temp5<-(NaccDetPvalProbeOrdered[,2]==row.names(NaccMedianCenteredData2))
length(NaccDetPvalProbeOrdered[,2])
length(NaccMedianCenteredData2[,1])
sum(Temp5)


str(IlluminaDecode)

length(IlluminaDecode$Array_Address_Id)

IlluminaDecodeSorted<-IlluminaDecode[order(as.character(IlluminaDecode$Array_Address_Id)),] 

head(IlluminaDecodeSorted)

IlluminaDecodeSortedTF<-IlluminaDecodeSorted$Array_Address_Id%in%row.names(AMYMedianCenteredData2)

length(IlluminaDecodeSortedTF)
sum(IlluminaDecodeSortedTF)

IlluminaDecodeSortedTF[c(1:30)]

IlluminaDecodeSorted22177<-IlluminaDecodeSorted[IlluminaDecodeSortedTF==T,]
dim(IlluminaDecodeSorted22177)

length(IlluminaDecodeSorted22177$Array_Address_Id)

plot(as.numeric(IlluminaDecodeSorted22177$Array_Address_Id)~as.numeric(row.names(AMYMedianCenteredData2)))


HCDetPvalProbeOrdered[c(1:20),]
HCMedianCenteredData2[c(1:20),]

dim(HCDetPvalProbeOrdered)
dim(HCMedianCenteredData2)

plot(as.numeric(HCDetPvalProbeOrdered[,2])~as.numeric(row.names(HCMedianCenteredData2)))


HCDetPvalProbeOrdered[c(1:20),]
HCMedianCenteredData2[c(1:20),]

dim(HCDetPvalProbeOrdered)
dim(HCMedianCenteredData2)

plot(as.numeric(HCDetPvalProbeOrdered[,2])~as.numeric(row.names(HCMedianCenteredData2)))


ACgDetPvalProbeOrdered[c(1:20),]
ANCGMedianCenteredData2[c(1:20),]

dim(ACgDetPvalProbeOrdered)
dim(ANCGMedianCenteredData2)

plot(as.numeric(ACgDetPvalProbeOrdered[,2])~as.numeric(row.names(ANCGMedianCenteredData2)))


AmyDetPvalProbeOrdered[c(1:20),]
AMYMedianCenteredData2[c(1:20),]

dim(AmyDetPvalProbeOrdered)
dim(AMYMedianCenteredData2)

plot(as.numeric(AmyDetPvalProbeOrdered[,2])~as.numeric(row.names(AMYMedianCenteredData2)))

DlpfcDetPvalProbeOrdered[c(1:20),]
DlpfcMedianCenteredData2[c(1:20),]

dim(DlpfcDetPvalProbeOrdered)
dim(DlpfcMedianCenteredData2)

plot(as.numeric(DlpfcDetPvalProbeOrdered[,2])~as.numeric(row.names(DlpfcMedianCenteredData2)))

CBDetPvalProbeOrdered[c(1:20),]
CBMedianCenteredData2[c(1:20),]

dim(CBDetPvalProbeOrdered)
dim(CBMedianCenteredData2)

plot(as.numeric(CBDetPvalProbeOrdered[,2])~as.numeric(row.names(CBMedianCenteredData2)))


#*******Pulling out the data based on detection p-value p=0.05 or p<0.05***********

sum(HCDetPvalProbeOrdered[,3]<0.055)
#[1] 13334

png("HC_DistributionOfDetPval.png")
plot(sort(HCDetPvalProbeOrdered[,3]), xlab="Probe Index: Ordered by Pvalue", ylab="Detection P-value for Probe", main="# of Probes at each level of detection")
dev.off()

png("HC_DistributionOfDetPval2.png")
hist(HCDetPvalProbeOrdered[,3], xlab="Detection P-value for Probes", ylab="# of Probes at each detection p-value", main="# of Probes at each level of detection", col="blue")
dev.off()


#Not useful because data is already median centered...
HCMedianOfMedianCenteredData<-apply(HCMedianCenteredData2, 1, median)
plot(HCMedianOfMedianCenteredData~HCDetPvalProbeOrdered[,3])


head(HCDetPvalProbeOrdered)

HCMedianCenteredDataDetected<-HCMedianCenteredData2[HCDetPvalProbeOrdered[,3]<0.055,]
head(HCMedianCenteredDataDetected)

dim(HCMedianCenteredDataDetected)

HCDetPvalDetected<-HCDetPvalProbeOrdered[HCDetPvalProbeOrdered[,3]<0.055,]
dim(HCDetPvalDetected)

hist(HCDetPvalDetected[,3])

HCdetectedIlluminaDecodeSorted22177<-IlluminaDecodeSorted22177[HCDetPvalProbeOrdered[,3]<0.055,]
dim(HCdetectedIlluminaDecodeSorted22177)

str(HCDetPvalDetected)
HC_Annotated_DetectedProbeData<-cbind(HCDetPvalDetected,HCMedianCenteredDataDetected, HCdetectedIlluminaDecodeSorted22177) 
dim(HC_Annotated_DetectedProbeData)

write.csv(HC_Annotated_DetectedProbeData, "HC_Annotated_DetectedProbeData.csv")


## Catch up here for the other 5 regions


*********************Problems ensued due to Illumina Decode having too much information**********************

ACgMedianCenteredDataDetected<-ANCGMedianCenteredData2[ACgDetPvalProbeOrdered[,3]<0.055,]
##dim (ACgMedianCenteredDataDetected)
[1] 13408   109

ACgDetPvalDetected<-ACgDetPvalProbeOrdered[ACgDetPvalProbeOrdered[,3]<0.055,]
dim(ACgDetPvalDetected)
[1] 13408     3

ACgdetectedIlluminaDecodeSorted22177<-IlluminaDecodeSorted22177[ACgDetPvalProbeOrdered[,3]<0.055,]
dim(ACgdetectedIlluminaDecodeSorted22177)
[1] 13408    28

ACg_Annotated_DetectedProbeData<-cbind(ACgMedianCenteredDataDetected,ACgDetPvalDetected,ACgdetectedIlluminaDecodeSorted22177)
dim(ACg_Annotated_DetectedProbeData)
[1] 13408   140

write.csv(ACg_Annotated_DetectedProbeData,"ACg_Annotated_DetectedProbeData")

ACg_Annotated_DetectedProbeData["2690537:2690541",]

print(ACg_Annotated_DetectedProbeData[3510,140])

head(ACg_Annotated_DetectedProbeData)
ACg_Annotated_DetectedProbeDataMinus<-ACg_Annotated_DetectedProbeData(,c(-136:-138))
Note: the last column of the file ACg_Annotated_DetectedProbeData has cells so large they spill over into another row.  Thus the last, 3 columns are removed ontology_component, ontology_function, ontology_process.  See code below. 
dim(ACg_Annotated_DetectedProbeData)
ACg_Annotated_DetectedProbeDataMinus<-ACg_Annotated_DetectedProbeData[1:135]
dim(ACg_Annotated_DetectedProbeDataMinus)
[1] 13408   135
ACg_Annotated_DetectedProbeDataAbrev<-ACg_Annotated_DetectedProbeDataMinus
write.csv(ACg_Annotated_DetectedProbeDataAbrev, "ACg_Annotated_DetectedProbeDataAbrev.csv")


HC_Annotated_DetectedProbeDataAbbrev<-HC_Annotated_DetectedProbeData[1:135]
write.csv(HC_Annotated_DetectedProbeDataAbbrev, "HC_Annotated_DetectedProbeDataAbbrev.csv")

AmyMedianCenteredDataDetected<-AMYMedianCenteredData2[AmyDetPvalProbeOrdered[,3]<0.055,]
dim (AmyMedianCenteredDataDetected)
[1] 13319    98

AmyDetPvalDetected<-AmyDetPvalProbeOrdered[AmyDetPvalProbeOrdered[,3]<0.055,]
dim(AmyDetPvalDetected)
[1] 13319     3

AmydetectedIlluminaDecodeSorted22177<-IlluminaDecodeSorted22177[AmyDetPvalProbeOrdered[,3]<0.055,]
dim(AmydetectedIlluminaDecodeSorted22177)
[1] 13319    28

Amy_Annotated_DetectedProbeData<-cbind(AmyMedianCenteredDataDetected,AmyDetPvalDetected,AmydetectedIlluminaDecodeSorted22177)
dim(Amy_Annotated_DetectedProbeData)
[1] 13319   129

Amy_Annotated_DetectedProbeDataAbbrev<-Amy_Annotated_DetectedProbeData[1:135]
dim(Amy_Annotated_DetectedProbeDataAbbrev)

Amytemp<-AmydetectedIlluminaDecodeSorted22177[1:23]
head(Amytemp)

dim(HC_Annotated_DetectedProbeData)
dim(HCdetectedIlluminaDecodeSorted22177)
HCtemp<-


dim(IlluminaDecodeSorted22177)
[1] 22177    28

head(IlluminaDecodeSorted22177)

IlluminaDecodeSorted22177Abbreviated<-IlluminaDecodeSorted22177[1:23]

dim(IlluminaDecodeSorted22177Abbreviated)
[1] 22177    23

head(IlluminaDecodeSorted22177Abbreviated)


**********************************************Fixes****************************************
##ACg cbinding
 
dim (ACgMedianCenteredDataDetected)
##[1] 13408   109

dim(ACgDetPvalDetected)
##[1] 13408     3

ACgdetectedIlluminaDecodeSorted22177<-IlluminaDecodeSorted22177[ACgDetPvalProbeOrdered[,3]<0.055,]
dim(ACgdetectedIlluminaDecodeSorted22177)
##[1] 13408    28

ACgdetectedIlluminaDecodeSorted22177Abbreviated<-IlluminaDecodeSorted22177Abbreviated[ACgDetPvalProbeOrdered[,3]<0.055,]
dim(ACgdetectedIlluminaDecodeSorted22177Abbreviated)
##[1] 13408    23

ACg_Annotated_DetectedProbeData<-cbind(ACgMedianCenteredDataDetected,ACgDetPvalDetected,ACgdetectedIlluminaDecodeSorted22177Abbreviated)
dim(ACg_Annotated_DetectedProbeData)
##[1] 13408   135

write.csv(ACg_Annotated_DetectedProbeData,"ACg_Annotated_DetectedProbeData.csv")

##HC cbinding

dim(HCMedianCenteredDataDetected)
##[1] 13334   109

dim(HCDetPvalDetected)
##[1] 13334     3

HCdetectedIlluminaDecodeSorted22177Abbreviated<-IlluminaDecodeSorted22177Abbreviated[HCDetPvalProbeOrdered[,3]<0.055,]
dim(HCdetectedIlluminaDecodeSorted22177Abbreviated)
##[1] 13334    23

HC_Annotated_DetectedProbeData<-cbind(HCDetPvalDetected,HCMedianCenteredDataDetected, HCdetectedIlluminaDecodeSorted22177Abbreviated) 
dim(HC_Annotated_DetectedProbeData)
##[1] 13334   135

write.csv(HC_Annotated_DetectedProbeData,"HC_Annotated_DetectedProbeData.csv")

##Amy cbinding

AmyMedianCenteredDataDetected<-AMYMedianCenteredData2[AmyDetPvalProbeOrdered[,3]<0.055,]
dim (AmyMedianCenteredDataDetected)
[1] 13319    98

AmyDetPvalDetected<-AmyDetPvalProbeOrdered[AmyDetPvalProbeOrdered[,3]<0.055,]
dim(AmyDetPvalDetected)
[1] 13319     3

AmydetectedIlluminaDecodeSorted22177Abbreviated<-IlluminaDecodeSorted22177Abbreviated[AmyDetPvalProbeOrdered[,3]<0.055,]
dim(AmydetectedIlluminaDecodeSorted22177Abbreviated)
[1] 13319    23

Amy_Annotated_DetectedProbeData<-cbind(AmyMedianCenteredDataDetected,AmyDetPvalDetected,AmydetectedIlluminaDecodeSorted22177Abbreviated)
dim(Amy_Annotated_DetectedProbeData)
[1] 13319   124

write.csv(Amy_Annotated_DetectedProbeData,"Amy_Annotated_DetectedProbeData.csv")

###CB cbinding

CBMedianCenteredDataDetected<-CBMedianCenteredData2[CBDetPvalProbeOrdered[,3]<0.055,]
dim(CBMedianCenteredDataDetected)
[1] 12931   113

CBDetPvalDetected<-CBDetPvalProbeOrdered[CBDetPvalProbeOrdered[,3]<0.055,]
dim(CBDetPvalDetected)
##[1] 12931     3

CBIlluminaDecodeSorted22177Abbreviated<-IlluminaDecodeSorted22177Abbreviated[CBDetPvalProbeOrdered[,3]<0.055,]
dim(CBIlluminaDecodeSorted22177Abbreviated)
[1] 12931    23

CB_Annotated_DectededProbeData<-cbind(CBMedianCenteredDataDetected,CBDetPvalDetected,CBIlluminaDecodeSorted22177Abbreviated)
dim(CB_Annotated_DectededProbeData)
[1] 12931   139

write.csv(CB_Annotated_DectededProbeData,"CB_Annotated_DectededProbeData.csv")

##Dlpfc cbinding

DlpfcMedianCenteredDataDetected<-DlpfcMedianCenteredData2[DlpfcDetPvalProbeOrdered[,3]<0.055,]
dim(DlpfcMedianCenteredDataDetected)
[1] 13370   109

DlpfcDetPvalDetected<-DlpfcDetPvalProbeOrdered[DlpfcDetPvalProbeOrdered[,3]<0.055,]
dim(DlpfcDetPvalDetected)
[1] 13370     3

DlpfcIlluminaDecodeSorted22177Abbreviated<-IlluminaDecodeSorted22177Abbreviated[DlpfcDetPvalProbeOrdered[,3]<0.055,]
dim(DlpfcIlluminaDecodeSorted22177Abbreviated)
[1] 13370    23

Dlpfc_Annotated_DetectedProbeData<-cbind(DlpfcMedianCenteredDataDetected,DlpfcDetPvalDetected,DlpfcIlluminaDecodeSorted22177Abbreviated)
dim(Dlpfc_Annotated_DetectedProbeData)
[1] 13370   135

write.csv(Dlpfc_Annotated_DetectedProbeData,"Dlpfc_Annotated_DetectedProbeData.csv")

##Nacc cbinding

NaccMedianCenteredDataDetected<-NaccMedianCenteredData2[NaccDetPvalProbeOrdered[,3]<0.055,]
dim(NaccMedianCenteredData)
[1] 13530   106

NaccDetPvalDetected<-NaccDetPvalProbeOrdered[NaccDetPvalProbeOrdered[,3]<0.055,]
dim(NaccDetPvalDetected)
[1] 13530     3

NaccIlluminaDecodeSorted22177Abbreviated<-IlluminaDecodeSorted22177Abbreviated[NaccDetPvalProbeOrdered[,3]<0.055,]
dim(NaccIlluminaDecodeSorted22177Abbreviated)
[1] 13530    23

Nacc_Annotated_DetectedProbeData<-cbind(NaccMedianCenteredDataDetected,NaccDetPvalDetected,NaccIlluminaDecodeSorted22177Abbreviated)
dim(Nacc_Annotated_DetectedProbeData)
[1] 13530   132

write.csv(Nacc_Annotated_DetectedProbeData,"Nacc_Annotated_DetectedProbeData.csv")

***********************Subject Demographics*****************************************

SubjectDemographics<-read.csv("All_Subject_Information_for_R.csv", header=T)


#ToDo:
#Pull out Subject Demographics by Region
#Remove Outliers (or just match subject ids with our microarray data) - %in%
#Put subjects in the same order as microarray data - order subject IDs (as.character)

 
SubjectDemographicsOrdered<-SubjectDemographics[order(SubjectDemographics$Sample.ID),]
# dim(SubjectDemographicsOrdered)
[1] 666  32
##There are data for 644 subjects, thus with 666 subjects, this has the outliers in it.

head(SampleIDsNoOutliersTF)
head(ChipIDsNoOutliers)
dim(ChipIDsNoOutliers)
[1] 644  14

head(SubjectDemographicsOrdered)
dim(SubjectDemographicsOrdered)

SubjectDemographicsOrderedTF<-ChipIDsNoOutliers[,4]%in%SubjectDemographicsOrdered[,6]
##Column 4 of ChipIDsNoOutliers has the subject ID and column 6 of SubjectDemographicsOrdered has the subject IDs.  
length(SubjectDemographicsOrderedTF)
[1] 644
print(SubjectDemographicsOrderedTF)
##Created a file names SubjectDemographics644Ordered but the looks like 6440 and is confusing.  Thus abandoned it.
SubjectDemographicsOrdered644<-SubjectDemographicsOrderedTF

head(SubjectDemographicsOrdered644)

This process created a file with all true false.  Trying again without changint to TF
SubjectDemographicsOrdered644<-ChipIDsNoOutliers[,4]%in%SubjectDemographicsOrdered[,6]
head(SubjectDemographicsOrdered644)
##Still a true/false file

head(ChipIDsNoOutliers)
##Subject IDs are in column 4
dim(ChipIDsNoOutliers)
##[1] 644  14
sort(ChipIDsNoOutliers [,4])
##This does sort the subject IDs
head(ChipIDsNoOutliers)
##however, the subject IDs do not stay in order when looking at the first 6 subjects in the header format


head(SubjectDemographicsOrdered)
##Subject IDs are in column 6
dim(SubjectDemographicsOrdered)
[1] 666  32

dim(SubjectDemographicsOrdered644)

sort(SubjectDemographicsOrdered[,6])
sort(ChipIDsNoOutliers [,4])
SubjectDemographicsOrdered644<-ChipIDsNoOutliers[,4]%in%SubjectDemographicsOrdered[,6]
dim(SubjectDemographicsOrdered644)
head(SubjectDemographicsOrdered644)
table(SubjectDemographicsOrdered644, ChipIDsNoOutliers [,3])
##SubjectDemographicsOrdered644 AMY ANCG  CB DLPFC  HC NACC
                         TRUE  98  109 113   109 109  106
sum(SubjectDemographicsOrdered644)
[1] 644

match(ChipIDsNoOutliers [,4],SubjectDemographicsOrdered[,6])
MergreTest<-merge(ChipIDsNoOutliers [,4],(SubjectDemographicsOrdered[,6])
dim(MergreTest)
head(MergreTest)
##     x    y
1 3281 3038
2 3878 3038
3 3426 3038
4 3671 3038
5 4236 3038
6 3031 3038

Matchtest<-Match(ChipIDsNoOutliers [,1],(SubjectDemographicsOrdered[,1])
MergeTest<-merge(ChipIDsNoOutliers [,1],(SubjectDemographicsOrdered[,1])
dim(MergeTest)

Match(ChipIDsNoOutliers,SubjectDemographicsOrdered)

##**************Didn't work**********
ChipIDsNoOutliersOrdered<-ChipIDsNoOutliers[order(ChipIDsNoOutliers[,4],)]
ChipIDsNoOutliersOrdered<-ChipIDsNoOutliers[order(ChipIDsNoOutliers[,4])]
##************************************

SubjectDemographicsNoOutliers<-SubjectDemographicsOrdered[SampleIDsNoOutliersTF==TRUE,]
dim(SubjectDemographicsNoOutliers)
##[1] 644  32
head(SubjectDemographicsNoOutliers)
head(SampleIDsNoOutliersTF)
length(SampleIDsNoOutliersTF)
## [1] 666
sum(SampleIDsNoOutliersTF)
##[1] 644sort()


SubjectDemographicsOrdered

SubjectDemographicsOrderedTF<-SubjectDemographicsOrdered[,6]%in%ChipIDsNoOutliers[,4]
print(SubjectDemographicsOrderedTF)
head(SubjectDemographicsOrdered)
str(SubjectDemographicsOrdered)

head(ChipIDsNoOutliers)
SubjectDemographicsOrderedTF<-SubjectDemographicsOrdered[,1]%in%ChipIDsNoOutliers[,1]
print(SubjectDemographicsOrderedTF)

table(SubjectDemographicsOrderedTF)

SubjectsDemographicsOrdered644<-SubjectDemographicsOrdered [SubjectDemographicsOrderedTF==TRUE,]
dim(SubjectsDemographicsOrdered644)
##[1] 644  32
head(SubjectsDemographicsOrdered644)
table(SubjectsDemographicsOrdered644 [,2])
##AMY  ANCG    CB DLPFC    HC  NACC 
   98   109   113   109   109   106

AmySubjectsDemographicsOrdered644<-SubjectsDemographicsOrdered644[c(1:98),]
table(AmySubjectsDemographicsOrdered644[,2])
##AMY  ANCG    CB DLPFC    HC  NACC 
   18    19    18    12    17    14 

##Order by column 2.


AMYSubjectsDemographicsOrdered<-SubjectsDemographicsOrdered644[SubjectsDemographicsOrdered644[,2]=="AMY", ]
dim(AMYSubjectsDemographicsOrdered)

head(AMYSubjectsDemographicsOrdered)

as.character(AMYSubjectsDemographicsOrdered[,2])==colnames(AMYMedianCenteredData2)

cbind(as.character(AMYSubjectsDemographicsOrdered[,1]),colnames(AMYMedianCenteredData2))


SubjectsDemographicsOrdered644RegionOrdered<-order(SubjectsDemographicsOrdered644()[,2],))
dim(SubjectsDemographicsOrdered644RegionOrdered)
table(SubjectsDemographicsOrdered644RegionOrdered)

head(AMYMedianCenteredData2)

#1. Get Median Centered Data in same order as Demographics. Probably the easiest way to do this is to order both data sets by ChipID (column 1 as character in Demographics and colnames Median Centered Data), then make sure that the order matches.

#2. Then assign the SubjectIDs from the demographcis file as the new colnames for a new MedianCentered data file.

#3. (optional) Then reorder both the demographics file and the Median Centered Data file using the subject ID so that the order will vaguely match the order for other brain regions.

order(as.character(colnames(AMYMedianCenteredData2)))
AMYMedianCenteredData3<-AMYMedianCenteredData2[ , order(as.character(colnames(AMYMedianCenteredData2)))]
cbind(as.character(AMYSubjectsDemographicsOrdered[,1]),colnames(AMYMedianCenteredData3))
as.character(AMYSubjectsDemographicsOrdered[,1])==colnames(AMYMedianCenteredData3)

CBSubjectsDemographicsOrdered<-SubjectsDemographicsOrdered644[SubjectsDemographicsOrdered644[,2]=="CB", ]
CBMedianCenteredData3<-CBMedianCenteredData2[ , order(as.character(colnames(CBMedianCenteredData2)))]
cbind(as.character(CBSubjectsDemographicsOrdered[,1]),colnames(CBMedianCenteredData3))
as.character(CBSubjectsDemographicsOrdered[,1])==colnames(CBMedianCenteredData3)

DlpfcSubjectsDemographicsOrdered<-SubjectsDemographicsOrdered644[SubjectsDemographicsOrdered644[,2]=="DLPFC", ]
DlpfcMedianCenteredData3<-DlpfcMedianCenteredData2[ , order(as.character(colnames(DlpfcMedianCenteredData2)))]
cbind(as.character(DlpfcSubjectsDemographicsOrdered[,1]),colnames(DlpfcMedianCenteredData3))
as.character(DlpfcSubjectsDemographicsOrdered[,1])==colnames(DlpfcMedianCenteredData3)

HCSubjectsDemographicsOrder<-SubjectsDemographicsOrdered644[SubjectsDemographicsOrdered644[,2]=="HC",]
HCMedianCenteredData3<-HCMedianCenteredData2[,order(as.character(colnames(HCMedianCenteredData2)))]
cbind(as.character(HCSubjectsDemographicsOrder[,1]),colnames(HCMedianCenteredData3))
as.character(HCSubjectsDemographicsOrder[,1])==colnames(HCMedianCenteredData3)

ACgSubjectsDemographicsOrdered<-SubjectsDemographicsOrdered644[SubjectDemographicsOrdered644[,2]=="ANCG",]
##s was missing in Subject (S) demographics ordered
ACgSubjectsDemographicsOrder<-SubjectsDemographicsOrdered644[SubjectsDemographicsOrdered644[,2]=="ANCG",]
#mispelling redone
ACgSubjectsDemographicsOrdered<-SubjectsDemographicsOrdered644[SubjectsDemographicsOrdered644[,2]=="ANCG",]
table(SubjectsDemographicsOrdered644 [,2])

HCSubjectsDemographicsOrder<-SubjectsDemographicsOrdered644[SubjectsDemographicsOrdered644[,2]=="HC",]
##mispelling redone
HCSubjectsDemographicsOrdered<-SubjectsDemographicsOrdered644[SubjectsDemographicsOrdered644[,2]=="HC",]
HCMedianCenteredData3<-HCMedianCenteredData2[order(as.character(HCMedianCenteredData2))]

NaccSubjectsDemographicsOrdered<-SubjectsDemographicsOrdered644[SubjectsDemographicsOrdered644[,2]=="NACC",]
#Megan Messing around, not really necessary:
NaccSubjectsDemographicsOrdered2<-NaccSubjectsDemographicsOrdered[order(as.character(NaccSubjectsDemographicsOrdered[,1])),]


HCMedianCenteredData3<-HCMedianCenteredData2[ ,order(as.character(colnames(HCMedianCenteredData2)))]
ACgMedianCenteredData3<-ANCGMedianCenteredData2[ ,order(as.character(colnames(ANCGMedianCenteredData2)))]
NaccMedianCenteredData3<-NaccMedianCenteredData2[ ,order(as.character(colnames(NaccMedianCenteredData2)))]
cbind(as.character(ACgSubjectsDemographicsOrder[,1]),colnames(ACgMedianCenteredData3))
as.character(ACgSubjectsDemographicsOrder[,1])==colnames(ACgMedianCenteredData3)

cbind(as.character(NaccSubjectsDemographicsOrdered[,1]),colnames(NaccMedianCenteredData3))
as.character(NaccSubjectsDemographicsOrdered[,1])==colnames(NaccMedianCenteredData3)
dim(NaccSubjectsDemographicsOrdered)
dim(NaccMedianCenteredData3)
order(as.character(NaccSubjectsDemographicsOrdered[,1]))
head(NaccSubjectsDemographicsOrdered)

colnames(ACg_Annotated_DetectedProbeData)
ACg_Annotated_DetectedProbeData2<-ACg_Annotated_DetectedProbeData[ ,order(as.character(colnames(ACg_Annotated_DetectedProbeData)))]
HC_Annotated_DetectedProbeData2<-HC_Annotated_DetectedProbeData[ ,order(as.character(colnames(HC_Annotated_DetectedProbeData)))]
Amy_Annotated_DetectedProbeData2<-Amy_Annotated_DetectedProbeData[ ,order(as.character(colnames(Amy_Annotated_DetectedProbeData)))]
CB_Annotated_DectededProbeData2<-CB_Annotated_DectededProbeData[ ,order(as.character(colnames(CB_Annotated_DectededProbeData)))]
Dlpfc_Annotated_DetectedProbeData2<-Dlpfc_Annotated_DetectedProbeData[ ,order(as.character(colnames(Dlpfc_Annotated_DetectedProbeData)))]
Nacc_Annotated_DetectedProbeData2<-Nacc_Annotated_DetectedProbeData[ ,order(as.character(colnames(Nacc_Annotated_DetectedProbeData)))]

cbind(as.character(ACgSubjectsDemographicsOrder[,1]),colnames(ACg_Annotated_DetectedProbeData2))


#How to change colnames:
MatrixOfInterest2<-MatrixOfInterest
colnames(MatrixOfInterest2)<-newcolnames_whateveryoulike
row.names(MatrixOfInterest2)<-row.names(MatrixOfInterest)


##Notes for further processing:
#At this point the subjects in demographics info is in the same order as in MedianCenteredData3, we originally intended to change the colnames in MedianCenteredData3 to the SubjectIDs instead of using the Chip IDs.
#It would be nice to do this for the filtered median centered data too (i.e. the data that has been filtered by detection p-value)

cbind(as.character(AMYSubjectsDemographicsOrdered[,1]),colnames(AMYMedianCenteredData3))
cbind(as.character(AMYSubjectsDemographicsOrdered[,1]),colnames(AmyMedianCenteredDataDetected2))
##created MedianCenteredDataDetected2 (p<0.055) ordered Chip ID
AmyMedianCenteredDataDetected2<-AmyMedianCenteredDataDetected[,order(as.character(colnames(AmyMedianCenteredDataDetected)))]
## used cbind(as.character(AMYSubjectsDemographicsOrdered[,1]),colnames(AmyMedianCenteredDataDetected2)) to check order and it did work
## as.character(AMYSubjectsDemographicsOrdered[,1])==colnames(AmyMedianCenteredDataDetected2)
## All True

ACgMedianCenteredDataDetected2<-ACgMedianCenteredDataDetected[,order(as.character(colnames(ACgMedianCenteredDataDetected)))]
## check as.character(ACgSubjectsDemographicsOrder[,1])==colnames(ACgMedianCenteredDataDetected2)  did not work only first 2 were true.
cbind(as.character(ACgSubjectsDemographicsOrder[,1]),colnames(ACgMedianCenteredDataDetected2))
order(as.character(colnames(ACgMedianCenteredDataDetected2)))
TemporaryObject<-ACgMedianCenteredDataDetected2
fix(TemporaryObject)
str(TemporaryObject)
TemporaryObject2<-data.frame(TemporaryObject)
rm(TemporaryObject)
ACgMedianCenteredDataDetected3<-data.frame(ACgMedianCenteredDataDetected)
(ACgSubjectsDemographicsOrder[,1])==colnames(ACgMedianCenteredDataDetected3)
order(as.character(colnames(ACgMedianCenteredDataDetected3)))
as.numeric(as.character(colnames(ACgMedianCenteredDataDetected3)))

CBMedianCenteredDataDetected2<-CBMedianCenteredDataDetected[,ordered(as.character(colnames(CBMedianCenteredDataDetected)))]
as.character(CBSubjectsDemographicsOrdered[,1])==colnames(CBMedianCenteredDataDetected2)
##All True

DlpfcMedianCenteredDataDetected2<-DlpfcMedianCenteredDataDetected[,ordered(as.character(colnames(DlpfcMedianCenteredDataDetected)))]
as.character(DlpfcSubjectsDemographicsOrdered[,1])==colnames(DlpfcMedianCenteredDataDetected2)
##False

NaccMedianCenteredDataDetected2<-NaccMedianCenteredDataDetected[,ordered(as.character(colnames(NaccMedianCenteredDataDetected)))]
as.character(NaccSubjectsDemographicsOrdered[,1])==colnames(NaccMedianCenteredDataDetected2)
##False
HCMedianCenteredDataDetected2<-HCMedianCenteredDataDetected[,ordered(as.character(colnames(HCMedianCenteredDataDetected)))]
##To rectify the descrepancy and since I can not figure out how to do so in R I will export the files, change the formating in Excel and read the files back in as new objects
Changed directory to get the files manipulated
write.csv(NaccMedianCenteredDataDetected2,"NaccMedianCenteredDataDetected2.csv")
write.csv(DlpfcMedianCenteredDataDetected2,"DlpfcMedianCenteredDataDetected2.csv")
write.csv(ACgMedianCenteredDataDetected2,"ACgMedianCenteredDataDetected2.csv")
write.csv(HCMedianCenteredDataDetected2,"HCMedianCenteredDataDetected2.csv")
ACgMedianCenteredDataDetected3<-read.csv("ACgMedianCenteredDataDetected2.csv",header=T,row.names=1)
ordered(as.character(colnames(ACgMedianCenteredDataDetected3)))
## check as.character(ACgSubjectsDemographicsOrder[,1])==colnames(ACgMedianCenteredDataDetected2)  did not work only first 2 were true.
cbind(as.character(ACgSubjectsDemographicsOrder[,1]),colnames(ACgMedianCenteredDataDetected3))
##The above DID NOT WORK
mixedsort(colnames(ACgMedianCenteredDataDetected3))
##New function from gtools
as.character(ACgSubjectsDemographicsOrder[,1])==colnames(ACgMedianCenteredDataDetected3)
##Did not align
mixedsort(row.names(ACgSubjectsDemographicsOrder))
mixedsort(ACgSubjectsDemographicsOrder[,1])==colnames(ACgMedianCenteredDataDetected3)
cbind(as.character(ACgSubjectsDemographicsOrder[,1]),colnames(ACgMedianCenteredDataDetected3))
cbind(as.character(ACgSubjectsDemographicsOrder[,1]),mixedsort(colnames(ACgMedianCenteredDataDetected3)))
as.character(ACgSubjectsDemographicsOrder[,1])==mixedsort(colnames(ACgMedianCenteredDataDetected3))
ACgMedianCenteredDataDetected4<-ACgMedianCenteredDataDetected3 [,mixedsort(colnames(ACgMedianCenteredDataDetected3))]
as.character(ACgSubjectsDemographicsOrder[,1])==(colnames(ACgMedianCenteredDataDetected4)
as.character(ACgSubjectsDemographicsOrder[,1])==as.character(colnames(ACgMedianCenteredDataDetected4))
dim(ACgSubjectsDemographicsOrder)
dim(ACgMedianCenteredDataDetected4)
##ACgMedianCenteredDataDetected3 was reassigned when changing colnames to subject IDs.  It no longer supports the object ACgMedianCenteredDataDetected4

##Median Centered Data Detected 2 had ordered instead of order fixed this as.  THIS WORKED!  Note: the mixedsort function above worked also.
CBMedianCenteredDataDetected2<-CBMedianCenteredDataDetected[,order(as.character(colnames(CBMedianCenteredDataDetected)))]
DlpfcMedianCenteredDataDetected2<-DlpfcMedianCenteredDataDetected[,order(as.character(colnames(DlpfcMedianCenteredDataDetected)))]
NaccMedianCenteredDataDetected2<-NaccMedianCenteredDataDetected[,order(as.character(colnames(NaccMedianCenteredDataDetected)))]
AmyMedianCenteredDataDetected2<-AmyMedianCenteredDataDetected[,order(as.character(colnames(AmyMedianCenteredDataDetected)))]
HCMedianCenteredDataDetected2<-HCMedianCenteredDataDetected[,order(as.character(colnames(HCMedianCenteredDataDetected)))]     
as.character(CBSubjectsDemographicsOrdered[,1])==colnames(CBMedianCenteredDataDetected2)
as.character(DlpfcSubjectsDemographicsOrdered[,1])==colnames(DlpfcMedianCenteredDataDetected2)
as.character(NaccSubjectsDemographicsOrdered[,1])==colnames(NaccMedianCenteredDataDetected2)
as.character(AMYSubjectsDemographicsOrdered[,1])==colnames(AmyMedianCenteredDataDetected2)
as.character(HCSubjectsDemographicsOrdered[,1])==colnames(HCMedianCenteredDataDetected2)
##made but no longer needed ACgMedianCenteredDataDetected3 and ACgMedianCenteredDataDetected4

dim(CBSubjectsDemographicsOrdered)
##[1] 113  32
dim(CBMedianCenteredDataDetected2)
##[1] 12931   113

TempSubj<-CBSubjectsDemographicsOrdered
TempMedCtrDataDetPval<-CBMedianCenteredDataDetected2
dim(TempSubj)
TempSubj[,1]==colnames(TempMedCtrDataDetPval)

#How to change colnames:
MatrixOfInterest2<-MatrixOfInterest
colnames(MatrixOfInterest2)<-newcolnames_whateveryoulike
row.names(MatrixOfInterest2)<-row.names(MatrixOfInterest)
colnames(TempMedCtrDataDetPval)<-row.names(TempSubj)
head(TempSubj)
head(TempMedCtrDataDetPval)
##Above did not work because the row names are not the subject_id.  That column is 6 in the TempSubj object.
colnames(TempMedCtrDataDetPval)<-TempSubj[,6]
head(TempMedCtrDataDetPval)
TempSubj[,6]==colnames(TempMedCtrDataDetPval)
##Yay Works!!! :)

##Creating Region_SubjectsDemographicsOrdered2 and Region_MedianCenteredDataDetected3 for 6 Regions followed by assigning subject_id from Region_SubjectsDemographicsOrdered2, i.e. column 6 to the column names of the corresponding Region_MedianCenteredDataDetected3.  Lastly, a TRUE/FALSE equation to determine if the subject_id of Region_SubjectsDemographicsOrdered2 is equivalent to the subject ids as column names of Region_MedianCenteredDataDetected3
CBSubjectsDemographicsOrdered2<-CBSubjectsDemographicsOrdered
CBMedianCenteredDataDetected3<-CBMedianCenteredDataDetected2
DlpfcSubjectsDemographicsOrdered2<-DlpfcSubjectsDemographicsOrdered
DlpfcMedianCenteredDataDetected3<-DlpfcMedianCenteredDataDetected2
NaccSubjectsDemographicsOrdered2<-NaccSubjectsDemographicsOrdered
NaccMedianCenteredDataDetected3<-NaccMedianCenteredDataDetected2
AMYSubjectsDemographicsOrdered2<-AMYSubjectsDemographicsOrdered
AmyMedianCenteredDataDetected3<-AmyMedianCenteredDataDetected2
HCSubjectsDemographicsOrdered2<-HCSubjectsDemographicsOrdered
HCMedianCenteredDataDetected3<-HCMedianCenteredDataDetected2
ACgSubjectsDemographicsOrdered2<-ACgSubjectsDemographicsOrdered
ACgMedianCenteredDataDetected3<-ACgMedianCenteredDataDetected2

colnames(CBMedianCenteredDataDetected3)<-CBSubjectsDemographicsOrdered2[,6]
colnames(CBMedianCenteredDataDetected3)==CBSubjectsDemographicsOrdered2[,6]
colnames(DlpfcMedianCenteredDataDetected3)<-DlpfcSubjectsDemographicsOrdered2[,6]
colnames(DlpfcMedianCenteredDataDetected3)==DlpfcSubjectsDemographicsOrdered2[,6]
colnames(NaccMedianCenteredDataDetected3)<-NaccSubjectsDemographicsOrdered2[,6]
colnames(NaccMedianCenteredDataDetected3)==NaccSubjectsDemographicsOrdered2[,6]
colnames(AmyMedianCenteredDataDetected3)<-AMYSubjectsDemographicsOrdered2[,6]
colnames(AmyMedianCenteredDataDetected3)==AMYSubjectsDemographicsOrdered2[,6]
colnames(HCMedianCenteredDataDetected3)<-HCSubjectsDemographicsOrdered2[,6]
colnames(HCMedianCenteredDataDetected3)==HCSubjectsDemographicsOrdered2[,6]
colnames(ACgMedianCenteredDataDetected3)<-ACgSubjectsDemographicsOrdered2[,6]
colnames(ACgMedianCenteredDataDetected3)==ACgSubjectsDemographicsOrdered2[,6]

CBMedianCenteredData4<-CBMedianCenteredData3
DlpfcMedianCenteredData4<-DlpfcMedianCenteredData3
NaccMedianCenteredData4<-NaccMedianCenteredData3
AMYMedianCenteredData4<-AMYMedianCenteredData3
HCMedianCenteredData4<-HCMedianCenteredData3
ACgMedianCenteredData4<-ACgMedianCenteredData3

colnames(CBMedianCenteredData4)<-CBSubjectsDemographicsOrdered2[,6]
colnames(CBMedianCenteredData4)==CBSubjectsDemographicsOrdered2[,6]
colnames(DlpfcMedianCenteredData4)<-DlpfcSubjectsDemographicsOrdered2[,6]
colnames(DlpfcMedianCenteredData4)==DlpfcSubjectsDemographicsOrdered2[,6]
colnames(NaccMedianCenteredData4)<-NaccSubjectsDemographicsOrdered2[,6]
colnames(NaccMedianCenteredData4)==NaccSubjectsDemographicsOrdered2[,6]
colnames(AMYMedianCenteredData4)<-AMYSubjectsDemographicsOrdered2[,6]
colnames(AMYMedianCenteredData4)==AMYSubjectsDemographicsOrdered2[,6]
colnames(HCMedianCenteredData4)<-HCSubjectsDemographicsOrdered2[,6]
colnames(HCMedianCenteredData4)==HCSubjectsDemographicsOrdered2[,6]
colnames(ACgMedianCenteredData4)<-ACgSubjectsDemographicsOrdered2[,6]
colnames(ACgMedianCenteredData4)==ACgSubjectsDemographicsOrdered2[,6]

CBSubjectsDemographicsOrdered2[1:5,1:6]
str(CBSubjectsDemographicsOrdered2)
str(CBMedianCenteredDataDetected3)
dim(CBMedianCenteredDataDetected3)
dim(CBSubjectsDemographicsOrdered2)

CBSubjectsDemographicsOrdered3<-CBSubjectsDemographicsOrdered2[order(as.character(CBSubjectsDemographicsOrdered2[,6])),]
CBMedianCenteredDataDetected4<-CBMedianCenteredDataDetected3[,order(colnames(CBMedianCenteredDataDetected3))]
CBMedianCenteredData5<-CBMedianCenteredData4[,order(colnames(CBMedianCenteredData4))]
colnames(CBMedianCenteredDataDetected4)==CBSubjectsDemographicsOrdered3[,6]
colnames(CBMedianCenteredData5)==CBSubjectsDemographicsOrdered3[,6]

DlpfcSubjectsDemographicsOrdered3<-DlpfcSubjectsDemographicsOrdered2[order(as.character(DlpfcSubjectsDemographicsOrdered2[,6])),]
DlpfcMedianCenteredDataDetected4<-DlpfcMedianCenteredDataDetected3[,order(colnames(DlpfcMedianCenteredDataDetected3))]
DlpfcMedianCenteredData5<-DlpfcMedianCenteredData4[,order(colnames(DlpfcMedianCenteredData4))]
colnames(DlpfcMedianCenteredDataDetected4)==DlpfcSubjectsDemographicsOrdered3[,6]
colnames(DlpfcMedianCenteredData5)==DlpfcSubjectsDemographicsOrdered3[,6]

NaccSubjectsDemographicsOrdered3<-NaccSubjectsDemographicsOrdered2[order(as.character(NaccSubjectsDemographicsOrdered2[,6])),]
NaccMedianCenteredDataDetected4<-NaccMedianCenteredDataDetected3[,order(colnames(NaccMedianCenteredDataDetected3))]
NaccMedianCenteredData5<-NaccMedianCenteredData4[,order(colnames(NaccMedianCenteredData4))]
colnames(NaccMedianCenteredDataDetected4)==NaccSubjectsDemographicsOrdered3[,6]
colnames(NaccMedianCenteredData5)==NaccSubjectsDemographicsOrdered3[,6]

AmySubjectsDemographicsOrdered3<-AMYSubjectsDemographicsOrdered2[order(as.character(AMYSubjectsDemographicsOrdered2[,6])),]
AmyMedianCenteredDataDetected4<-AmyMedianCenteredDataDetected3[,order(colnames(AmyMedianCenteredDataDetected3))]
AmyMedianCenteredData5<-AMYMedianCenteredData4[,order(colnames(AMYMedianCenteredData4))]
colnames(AmyMedianCenteredDataDetected4)==AmySubjectsDemographicsOrdered3[,6]
colnames(AmyMedianCenteredData5)==AmySubjectsDemographicsOrdered3[,6]

HCSubjectsDemographicsOrdered3<-HCSubjectsDemographicsOrdered2[order(as.character(HCSubjectsDemographicsOrdered2[,6])),]
HCMedianCenteredDataDetected4<-HCMedianCenteredDataDetected3[,order(colnames(HCMedianCenteredDataDetected3))]
HCMedianCenteredData5<-HCMedianCenteredData4[,order(colnames(HCMedianCenteredData4))]
colnames(HCMedianCenteredDataDetected4)==HCSubjectsDemographicsOrdered3[,6]
colnames(HCMedianCenteredData5)==HCSubjectsDemographicsOrdered3[,6]

ACgSubjectsDemographicsOrdered3<-ACgSubjectsDemographicsOrdered2[order(as.character(ACgSubjectsDemographicsOrdered2[,6])),]
ACgMedianCenteredDataDetected4<-ACgMedianCenteredDataDetected3[,order(colnames(ACgMedianCenteredDataDetected3))]
ACgMedianCenteredData5<-ACgMedianCenteredData4[,order(colnames(ACgMedianCenteredData4))]
colnames(ACgMedianCenteredDataDetected4)==ACgSubjectsDemographicsOrdered3[,6]
colnames(ACgMedianCenteredData5)==ACgSubjectsDemographicsOrdered3[,6]


athal74readin<-read.delim("athal74_corrected.txt", header=T, sep="\t")
dim(athal74readin)
##[1] 22177    78
str(athal74readin)
head(athal74readin)
colnames(athal74readin[c(1:5),1:3])
row.names(athal74readin)[c(1:5)]
athal74readin2<-athal74readin[,c(-1, -76, -77, -78)]
row.names(athal74readin2)<-athal74readin[,1]
dim(athal74readin2)
head(athal74readin2)

athal74readin3<-athal74readin2
TestofgSub<-gsub("[:alnum:]","[^[:alnum]]",athal74readin3)
colnamesaThal74readin3<-colnames(athal74readin3)
dim(colnamesaThal74readin3)
colnames(colnamesaThal74readin3)
length(colnamesaThal74readin3)
row.names(athal74readin3)
colnames(athal74readin3)
columnnames<-colnames(athal74readin3)
columnnames1<-gsub("X","", columnnames)
colnames(athal74readin3)<-columnnames1

mthal77readin<-read.table("mthal77_corrected3.txt", header=T, sep="")

dim(mthal77readin)
str(mthal77readin)
row.names(mthal77readin)
columnnames<-colnames(mthal77readin)
columnnames1<-gsub("X","", columnnames)
colnames(mthal77readin2)<-columnnames1
str(mthal77readin2)
head(mthal77readin2[1:5,1:5])
mthal77readin3<-mthal77readin2

PCg85readin<-read.table(as.character("pcg85_corrected2.txt",header=T))
row.names(PCg85readin)
str(PCg85readin)
PCg85readin2<-PCg85readin
columnnames<-colnames(PCg85readin)
columnnames1<-gsub("X","", columnnames)
colnames(PCg85readin2)<-columnnames1
head(PCg85readin2[1:5,1:5])

colnames(PCg85readin2)==colnames(mthal77readin2)
cbind(colnames(PCg85readin2),colnames(mthal77readin2))
PCg85readin3<-PCg85readin2

SCg86readin<-read.delim("scg86_corrected_w_subject_IDs.txt",header=T)
str(SCg86readin)
dim(SCg86readin)
SCg86readin2<-SCg86readin[,c(-88,-89)]
str(SCg86readin2)
dim(SCg86readin2)
colnames(SCg86readin2)
SCg86readin2[1:5,1]
SCg86readin3<-SCg86readin2
row.names(SCg86readin3)<-SCg86readin2[,1]
str(SCg86readin3)
columnnames<-colnames(SCg86readin2)
columnnames1<-gsub("X","", columnnames)
colnames(SCg86readin3)<-columnnames1
head(SCg86readin3[1:5,1:5])
SCg86readin3rownamesfixed<-SCg86readin3[,-1]
head(SCg86readin3rownamesfixed[1:5,1:5])
SCg86readin3a<-SCg86readin3rownamesfixed

##Not used because I could not get the code below to accept numeric names for the column headers.
scg86_sampleIDs<-read.delim("SCg86_sampleEdittedfromSCg88.txt", header=T)
str(scg86_sampleIDs)
row.names(scg86_sampleIDs)
dim(scg86_sampleIDs)
## read as a data frame 86 observations with 6 variables, column 1 are the subject IDs, no row names.
scg86_sampleIDs[1:5,1]
##[1] 1834 1881 2169 2208 2248
scg86_sampleIDs2<-scg86_sampleIDs
scg86_sampleIDs2[,1]
colnames(scg86_sampleIDs2)==scg86_sampleIDs2[,1]
##There are no colnames yet.
SCg86readin2<-row.names(SCg86readin2[,1])
colnames(SCg86readin2)==row.names(as.character(scg86_sampleIDs))
colnames(SCg86readin2)
##Above useless
##Fixed by using
##SCg86readin3rownamesfixed<-SCg86readin3[,-1]
##Then SCg86readin3a<-SCg86readin3rownamesfixed


##Fix for Region_Annotated_DetectedProbeData and Region_Annotated_DetectedProbeDataWithSubjectID

colnames(ACgMedianCenteredDataDetected4)==ACgSubjectsDemographicsOrdered3[,6]
colnames(ACgMedianCenteredData5)==ACgSubjectsDemographicsOrdered3[,6]

write.csv(Dlpfc_Annotated_DetectedProbeData,"Dlpfc_Annotated_DetectedProbeData.csv")
write.csv(CB_Annotated_DectededProbeData,"CB_Annotated_DectededProbeData.csv")
write.csv(HC_Annotated_DetectedProbeData,"HC_Annotated_DetectedProbeData.csv")
write.csv(Amy_Annotated_DetectedProbeData,"Amy_Annotated_DetectedProbeData.csv")
write.csv(ACg_Annotated_DetectedProbeData,"ACg_Annotated_DetectedProbeData.csv")
write.csv(Nacc_Annotated_DetectedProbeData, "Nacc_Annotated_DetectedProbeData.csv")

write.csv(ACg_Annotated_DetectedProbeDataWithSubjectID,"ACg_Annotated_DetectedProbeDataWithSubjectID.csv")
write.csv(HC_Annotated_DetectedProbeDataWithSubjectID,"HC_Annotated_DetectedProbeDataWithSubjectID.csv")
write.csv(Amy_Annotated_DetectedProbeDataWithSubjectID,"Amy_Annotated_DetectedProbeDataWithSubjectID.csv")
write.csv(CB_Annotated_DectededProbeDataWithSubjectID,"CB_Annotated_DectededProbeDataWithSubjectID.csv")
write.csv(Dlpfc_Annotated_DetectedProbeDataWithSubjectID,"Dlpfc_Annotated_DetectedProbeDataWithSubjectID.csv")
write.csv(Nacc_Annotated_DetectedProbeDataWithSubjectID,"Nacc_Annotated_DetectedProbeDataWithSubjectID.csv")

## 4 more regions pval and getting things in order
SCg86readin2
SCg86readin3a
dim(SCgDetPvalProbeOrdered)
dim(SCg86readin2)
str(SCg86readin2)
str(SCg86readin3a)
row.names(SCg86readin3a)

Temp6<-(SCgDetPvalProbeOrdered[,2]==row.names(SCg86readin3a))
sum(Temp6)
dim(Temp6)
str(Temp6)
str(SCgDetPvalProbeOrdered)
SCg86readin2a<-SCg86readin2[order(as.character(row.names(SCg86readin2))),]
SCgDetPvalProbeOrdered<-SCgDetPval[order(as.character(SCgDetPval[,2])),]
SCgDetPvalProbeOrdered[,2]==row.names(as.character(SCg86readin3a))
row.names(SCg86readin3a)
SCg86readin2a<-SCg86readin2[order(as.character(row.names(SCg86readin3a))),]
##Row names of SCg86readin2 and SCg86readin3a are not lining up with SCgDetPvalProbeOrdered[,2]

plot(as.numeric(row.names(SCg86readin2))~as.numeric(row.names(SCg86readin3a)))
row.names(SCg86readin3a)
min(row.names(SCg86readin3a))
min(row.names(SCg86readin2[,1]))
min(SCg86readin2[,1])
plot(order(as.numeric(SCg86readin2[,2]))~order(as.numeric(row.names(SCg86readin3a))))
plot(as.character(SCg86readin2[,2])~as.character(row.names(SCg86readin3a)))
SCgDetPvalProbeOrdered[,2]
range(SCgDetPvalProbeOrdered[,2])
##[1]   10008 7650767
range(row.names(SCg86readin3a))
[1] "10008"  "990768"
##Note the string values in SCg86readin3a
row.names(as.numeric(SCg86readin3a))
str(SCgDetPvalProbeOrdered[,2])
row.names(as.numeric(unlist(SCg86readin3a)))
SCgDetPvalProbeOrdered[,2]==row.names(SCg86readin3a)
str(SCg86readin3a)
order(row.names(SCg86readin3a))
order(SCgDetPvalProbeOrdered[,2])
plot(as.numeric(SCgDetPvalProbeOrdered[,2])~as.numeric(row.names(SCg86readin3a)))
order(as.numeric(row.names(SCg86readin3a)))
write.csv(SCg86readin3a,"SCg86readin3a.csv")


SCg86readin3b<-read.csv(as.matrix("SCg86readin3a.csv",header=T))
str(SCg86readin3b)
SCg86readin3c<-SCg86readin3b
row.names(as.numeric(SCg86readin3c))<-SCg86readin3b[,1]
order(row.names(SCg86readin3c))
plot(as.numeric(SCgDetPvalProbeOrdered[,2])~as.numeric(row.names(SCg86readin3c)))
range(SCgDetPvalProbeOrdered[,2])
range(row.names(SCg86readin3c))
SCg86readin3d<-as.matrix(as.numeric(SCg86readin3c))
str(SCg86readin3d)
SCgpvalprobesdatarownamesTF<-SCgDetPvalProbeOrdered[,2]%in%row.names(SCg86readin3c)
sum(SCgpvalprobesdatarownamesTF)
##[1] 22177
junk<-SCgDetPvalProbeOrdered
junk[,2]==row.names(SCg86readin3c)
sort(row.names(SCg86readin3d))
sort(SCgDetPvalProbeOrdered[,2])
plot(as.numeric(SCgDetPvalProbeOrdered[,2])~as.character(row.names(SCg86readin3c)))


Temp3<-(CBDetPvalProbeOrdered[,2]==row.names(CBMedianCenteredData2))
CBDetPvalDetected<-CBDetPvalProbeOrdered[CBDetPvalProbeOrdered[,3]<0.055,]


PCg85readin
str(PCg85readin3)
row.names(PCg85readin3)
PCgDetPvalProbeOrdered[,2]
PCgpvalprobesdatarownamesTF<-PCgDetPvalProbeOrdered[,2]%in%row.names(PCg85readin3)
sum(PCgpvalprobesdatarownamesTF)
plot(as.numeric(PCgDetPvalProbeOrdered[,2])~as.numeric(row.names(PCg85readin3)))
order(as.numeric(PCgDetPvalProbeOrdered[,2]))
PCg85readin4<-PCg85readin3[order(row.names(PCg85readin3)),]
plot(as.numeric(PCgDetPvalProbeOrdered[,2])~as.numeric(row.names(PCg85readin4)))
PCgDetPvalProbeOrdered[,2]==row.names(PCg85readin4)


SCg86readin4<-SCg86readin3b[order(as.numeric(row.names(SCg86readin3b))),]
plot(as.numeric(SCgDetPvalProbeOrdered[,2])~as.numeric(row.names(SCg86readin4)))


mthal77readin
str(mthal77readin3)
row.names(mthal77readin3)
mThalpvalprobesdatarownamesTF<-PCgDetPvalProbeOrdered[,2]%in%row.names(mthal77readin3)
sum(mThalpvalprobesdatarownamesTF)
plot(as.numeric(mThalDetPvalProbeOrdered[,2])~as.numeric(row.names(mthal77readin3)))
order(as.numeric(mThalDetPvalProbeOrdered[,2]))
mthal77readin4<-mthal77readin3[order(row.names(mthal77readin3)),]
plot(as.numeric(mThalDetPvalProbeOrdered[,2])~as.numeric(row.names(mthal77readin4)))
mThalDetPvalProbeOrdered[,2]==row.names(mthal77readin4)
range(row.names(mthal77readin4))
##[1] "10008"  "990768"
range(mThalDetPvalProbeOrdered[,2])
##[1]   10008 7650767
length(mThalDetPvalProbeOrdered[,2])
##[1] 22177
length(row.names(mthal77readin4))
##[1] 22177
mThalpvalprobesdatarownamesTF<-mThalDetPvalProbeOrdered[,2]%in%row.names(mthal77readin4)
sum(mThalpvalprobesdatarownamesTF)
## [1] 22177
is.numeric(row.names(mthal77readin4))
##[1] FALSE
mthal77readin4<-mthal77readin3[order(as.numeric(row.names(mthal77readin3))),]
is.numeric(row.names(mthal77readin4))
##[1] FALSE
mThalpvalprobesdatarownamesTF2<-row.names(mthal77readin4)%in%mThalDetPvalProbeOrdered[,2]
sum(mThalpvalprobesdatarownamesTF2)

row.names(mthal77readin4)==mThalDetPvalProbeOrdered[,2]
Junk3<-cbind(row.names(mthal77readin4), mThalDetPvalProbeOrdered[,2])
is.numeric(mThalDetPvalProbeOrdered[,2])
##[1] TRUE
Junk4<-mthal77readin4[as.character(row.names(mthal77readin4)),]
Junk5<-mthal77readin4[as.numeric(row.names(mthal77readin4)),]
row.names(Junk4)==row.names(Junk5)
row.names(Junk4)
row.names(Junk5)
str(Junk5)



athal74readin2
athal74readin3
##~~~~~~~~~~~~~~~~~~~~Trying to read in the 4 regions again, cannot seem to get the row names to convert to numeric and align with Region_DetPvalProbeOrdered[,2] objects~~~~~~~~~~~~~~~~~~~~~~~~~~

##as an example for the other files:  PCg85readin<-read.table(as.character("pcg85_corrected2.txt",header=T))
athal74DataImported<-read.table(as.character("athal74_corrected.txt", header=T))
athal74DataImported2<-athal74DataImported
str(athal74DataImported2)##'data.frame':	22177 obs. of  76 variables: This removed 2 statistical columns from the end.
athal74DataImported2<-athal74DataImported2[,c(-75,-76)]
str(athal74DataImported2)
##'data.frame':	22177 obs. of  74 variables:
row.names(athal74DataImported2)
is.numeric(row.names(athal74DataImported2))
##False
aThalpvalprobesdatarownamesTF<-aThalDetPvalProbeOrdered[,2]%in%row.names(athal74DataImported2)
sum(aThalpvalprobesdatarownamesTF)
plot(as.numeric(aThalDetPvalProbeOrdered[,2])~as.numeric(row.names(athal74DataImported2)))
aThalDetPvalProbeOrdered(order(as.numeric(aThalDetPvalProbeOrdered[,2]))
athal74DataImported2<-athal74DataImported2[order(row.names(athal74DataImported2)),]
plot(as.numeric(aThalDetPvalProbeOrdered[,2])~as.numeric(row.names(athal74DataImported2)))
Junk9<-as.numeric(aThalDetPvalProbeOrdered[,2])
Junk10<-as.numeric(row.names(athal74DataImported2))
Junk9==Junk10
##True

mthal77DataImported<-read.table(as.character("mthal77_corrected3.txt", header=T))
mthal77DataImported2<-mthal77DataImported
str(mthal77DataImported2)##'data.frame':	22177 obs. of  78 variables: This removed 1 statistical columns from the end.
mthal77DataImported2<-mthal77DataImported2[,-78]
str(mthal77DataImported2)
##'data.frame':	22177 obs. of  77 variables:
row.names(mthal77DataImported2)
is.numeric(row.names(mthal77DataImported2))
##False
mThalpvalprobesdatarownamesTF<-mThalDetPvalProbeOrdered[,2]%in%row.names(mthal77DataImported2)
sum(mThalpvalprobesdatarownamesTF)
plot(as.numeric(mThalDetPvalProbeOrdered[,2])~as.numeric(row.names(mthal77DataImported2)))
mThalDetPvalProbeOrdered(order(as.numeric(mThalDetPvalProbeOrdered[,2]))
mthal77DataImported2<-mthal77DataImported2[order(row.names(mthal77DataImported2)),]
plot(as.numeric(mThalDetPvalProbeOrdered[,2])~as.numeric(row.names(mthal77DataImported2)))
Junk7<-as.numeric(mThalDetPvalProbeOrdered[,2])
Junk8<-as.numeric(row.names(mthal77DataImported2))
Junk7==Junk8
##True

SCg86DataImported<-read.table(as.character("scg86_corrected_w_subject_IDs.txt", header=T))
SCg86DataImported2<-SCg86DataImported
str(SCg86DataImported2)##'data.frame':	22177 obs. of  87 variables: This removed 1 statistical columns from the end.
SCg86DataImported2<-SCg86DataImported2[,-87]
str(SCg86DataImported2)
##'data.frame':	22177 obs. of  86 variables:
row.names(SCg86DataImported2)
is.numeric(row.names(SCg86DataImported2))
##False
SCgpvalprobesdatarownamesTF<-SCgDetPvalProbeOrdered[,2]%in%row.names(SCg86DataImported2)
sum(SCgpvalprobesdatarownamesTF)
plot(as.numeric(SCgDetPvalProbeOrdered[,2])~as.numeric(row.names(SCg86DataImported2)))
SCgDetPvalProbeOrdered(order(as.numeric(SCgDetPvalProbeOrdered[,2]))
SCg86DataImported2<-SCg86DataImported2[order(row.names(SCg86DataImported2)),]
plot(as.numeric(SCgDetPvalProbeOrdered[,2])~as.numeric(row.names(SCg86DataImported2)))
Junk11<-as.numeric(SCgDetPvalProbeOrdered[,2])
Junk12<-as.numeric(row.names(SCg86DataImported2))
Junk11==Junk12
##True
##Checking up where I left off to make sure the PCg files are in the correcte format and aligned properly.
PCg85DataImported2<-PCg85readin4
plot(as.numeric(PCgDetPvalProbeOrdered[,2])~as.numeric(row.names(PCg85readin4)))
PCgDetPvalProbeOrdered[,2]==row.names(PCg85readin4)

##files for gsub function to fix column names~~~~~SCg86DataImported2, mthal77DataImported2, athal74DataImported2,PCg85DataImported2
fix(athal74DataImported2)
##Removed 3933.1 to 3933
colnames(athal74DataImported2)
##column names are fixed for athal74DataImported2
SCg86DataImported3<-SCg86DataImported2
columnnames<-colnames(SCg86DataImported2)
columnnames1<-gsub("X","", columnnames)
colnames(SCg86DataImported3)<-columnnames1
colnames(SCg86DataImported3)

mthal77DataImported3<-mthal77DataImported2
columnnames<-colnames(mthal77DataImported2)
columnnames1<-gsub("X","", columnnames)
colnames(mthal77DataImported3)<-columnnames1
colnames(mthal77DataImported3)

athal74DataImported3<-athal74DataImported2
columnnames<-colnames(athal74DataImported2)
columnnames1<-gsub("X","", columnnames)
colnames(athal74DataImported3)<-columnnames1
colnames(athal74DataImported3)

PCg85DataImported3<-PCg85DataImported2
columnnames<-colnames(PCg85DataImported2)
columnnames1<-gsub("X","", columnnames)
colnames(PCg85DataImported3)<-columnnames1
colnames(PCg85DataImported3)
row.names(PCg85DataImported3[1:5,1:5])
row.names(PCg85DataImported3[1:5])
row.names(PCg85DataImported3[c(1:5),])
row.names(PCg85DataImported3)[c(1:5)]
row.names(athal74DataImported3[1:5,1:5])
row.names(mthal77DataImported3[1:5,1:5])
row.names(SCg86DataImported3[1:5,1:5])

SCgDetPvalProbeOrdered[,2]==row.names(SCg86DataImported3)
##TRUE
SCg86DataImported3MedCtrDataDetected<-SCg86DataImported3[SCgDetPvalProbeOrdered[,3]<0.055,]
dim(SCg86DataImported3MedCtrDataDetected)
##[1] 13919    86
SCgDetPvalDetected<-SCgDetPvalProbeOrdered[SCgDetPvalProbeOrdered[,3]<0.055,]
dim(SCgDetPvalDetected)
##[1] 13919     3
plot(row.names(SCg86DataImported3)~IlluminaDecodeSorted22177Abbreviated[,15])
row.names(SCg86DataImported3)==IlluminaDecodeSorted22177Abbreviated[,15]
SCgdetectedIlluminaDecodeSorted22177Abbreviated<-IlluminaDecodeSorted22177Abbreviated[SCgDetPvalProbeOrdered[,3]<0.055,]
dim(SCgdetectedIlluminaDecodeSorted22177Abbreviated)
##[1] 13919    23
SCg_Annotated_DetectedProbeData<-cbind(SCg86DataImported3MedCtrDataDetected,SCgDetPvalDetected, SCgdetectedIlluminaDecodeSorted22177Abbreviated)
dim(SCg_Annotated_DetectedProbeData)
##[1] 13919   112

##PCg85DataImported3
plot((PCgDetPvalProbeOrdered[,2])~row.names(PCg85DataImported3))
PCg85DataImported3MedCtrDataDetected<-PCg85DataImported3[PCgDetPvalProbeOrdered[,3]<0.055,]
dim(PCg85DataImported3MedCtrDataDetected)
##[1] 13884    85
PCgDetPvalDetected<-PCgDetPvalProbeOrdered[PCgDetPvalProbeOrdered[,3]<0.055,]
dim(PCgDetPvalDetected)
##[1] 13884     3
plot(row.names(PCg85DataImported3)~IlluminaDecodeSorted22177Abbreviated[,15])
PCgdetectedIlluminaDecodeSorted22177Abbreviated<-IlluminaDecodeSorted22177Abbreviated[PCgDetPvalProbeOrdered[,3]<0.055,]
dim(PCgdetectedIlluminaDecodeSorted22177Abbreviated)
##[1] 13884    23
PCg_Annotated_DetectedProbeData<-cbind(PCg85DataImported3MedCtrDataDetected,PCgDetPvalDetected, PCgdetectedIlluminaDecodeSorted22177Abbreviated)
dim(PCg_Annotated_DetectedProbeData)
##[1] 13884   111

##athal74DataImported3
plot((aThalDetPvalProbeOrdered[,2])~row.names(athal74DataImported3))
aThal74DataImported3MedCtrDataDetected<-athal74DataImported3[aThalDetPvalProbeOrdered[,3]<0.055,]
dim(aThal74DataImported3MedCtrDataDetected)
##[1] 13771    74
aThalDetPvalDetected<-aThalDetPvalProbeOrdered[aThalDetPvalProbeOrdered[,3]<0.055,]
dim(aThalDetPvalDetected)
##[1] 13771     3
plot(row.names(athal74DataImported3)~IlluminaDecodeSorted22177Abbreviated[,15])
aThaldetectedIlluminaDecodeSorted22177Abbreviated<-IlluminaDecodeSorted22177Abbreviated[aThalDetPvalProbeOrdered[,3]<0.055,]
dim(aThaldetectedIlluminaDecodeSorted22177Abbreviated)
##[1] 13771    23
aThal_Annotated_DetectedProbeData<-cbind(aThal74DataImported3MedCtrDataDetected, aThalDetPvalDetected, aThaldetectedIlluminaDecodeSorted22177Abbreviated)
dim(aThal_Annotated_DetectedProbeData)
##[1] 13771   100

##mthal77DataImported3
plot((mThalDetPvalProbeOrdered[,2])~row.names(mthal77DataImported3))
mThal77DataImported3MedCtrDataDetected<-mthal77DataImported3[mThalDetPvalProbeOrdered[,3]<0.055,]
dim(mThal77DataImported3MedCtrDataDetected)
##[1] 13690    77
mThalDetPvalDetected<-mThalDetPvalProbeOrdered[mThalDetPvalProbeOrdered[,3]<0.055,]
dim(mThalDetPvalDetected)
##[1] 13690     3
plot(row.names(mthal77DataImported3)~IlluminaDecodeSorted22177Abbreviated[,15])
mThaldetectedIlluminaDecodeSorted22177Abbreviated<-IlluminaDecodeSorted22177Abbreviated[mThalDetPvalProbeOrdered[,3]<0.055,]
dim(mThaldetectedIlluminaDecodeSorted22177Abbreviated)
##[1] 13690    23
mThal_Annotated_DetectedProbeData<-cbind(mThal77DataImported3MedCtrDataDetected, mThalDetPvalDetected, mThaldetectedIlluminaDecodeSorted22177Abbreviated)
dim(mThal_Annotated_DetectedProbeData)
##[1] 13690   103

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##Changing directory and reading in Subjects with Illness
SubjectsWithIllness<-read.csv("Subjects with Illness.csv", header=T)
str(SubjectsWithIllness)
##Subject IDs are in column 1
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
aThalMedianCenteredDataDetected5<-aThal74DataImported3MedCtrDataDetected
mThalMedianCenteredDataDetected5<-mThal77DataImported3MedCtrDataDetected
PCgMedianCenteredDataDetected5<-PCg85DataImported3MedCtrDataDetected
SCgMedianCenteredDataDetected5<-SCg86DataImported3MedCtrDataDetected
HCMedianCenteredDataDetected5<-HCMedianCenteredDataDetected4
ACgMedianCenteredDataDetected5<-ACgMedianCenteredDataDetected4
AmyMedianCenteredDataDetected5<-AmyMedianCenteredDataDetected4
CBMedianCenteredDataDetected5<-CBMedianCenteredDataDetected4
DlpfcMedianCenteredDataDetected5<-DlpfcMedianCenteredDataDetected4
NaccMedianCenteredDataDetected5<-NaccMedianCenteredDataDetected4

AThalSubjectIllnessTF<-SubjectsWithIllness[,1])%in%row.names(aThalMedianCenteredDataDetected5)
AThalSubjectIllnessTF
sum(AThalSubjectIllnessTF)
AThalSubjectIllnessTF<-as.numeric(SubjectsWithIllness[,1])%in%row.names(as.numeric(aThalMedianCenteredDataDetected5))
str(SubjectsWithIllness)
str(aThalMedianCenteredDataDetected5)

aThalSubjects<-colnames(aThalMedianCenteredDataDetected5)
AThalSubjectIllnessTF<-SubjectsWithIllness[,1]%in%aThalSubjects
sum(AThalSubjectIllnessTF)
##[1] 74

aThalSubjectsIllnessCategories<-SubjectsWithIllness[AThalSubjectIllnessTF==TRUE,]
dim(aThalSubjectsIllnessCategories)
[1] 74  9

mThalSubjects<-colnames(mThalMedianCenteredDataDetected5)
mThalSubjectIllnessTF<-SubjectsWithIllness[,1]%in%mThalSubjects
sum(mThalSubjectIllnessTF)
##[1] 77
mThalSubjectIllnessTF

mThalSubjectsIllnessCategories<-SubjectsWithIllness[mThalSubjectIllnessTF==TRUE,]
dim(mThalSubjectsIllnessCategories)
## [1] 77  9

PCgSubjects<-colnames(PCgMedianCenteredDataDetected5)
PCgSubjectIllnessTF<-SubjectsWithIllness[,1]%in%PCgSubjects
sum(PCgSubjectIllnessTF)
PCgSubjectIllnessTF
colnames(PCgMedianCenteredDataDetected5)==SubjectsWithIllness[,1]

temp<-PCgSubjects%in%SubjectsWithIllness[,1]
sum(temp)
PCgSubjects[temp==F]
[1] "3015" "3956"
cbind(SubjectsWithIllness[,1], PCgSubjectIllnessTF)
##Looking Back at the Pritzker Subject Info, 3015 is highly likely to be subject 4015 and subject 3956 is highly likely to be subject 3965 based on diagnoses, sex and cohort.
fix(PCgMedianCenteredDataDetected5)
colnames(PCgMedianCenteredDataDetected5)
PCgMedianCenteredDataDetected5With2SubjFixed<-PCgMedianCenteredDataDetected5[,sort(colnames(PCgMedianCenteredDataDetected5))]
colnames(PCgMedianCenteredDataDetected5With2SubjFixed)

head(PCgMedianCenteredDataDetected5[,c(67:69)])
head(PCgMedianCenteredDataDetected5With2SubjFixed[,c(67:69)])

plot(PCgMedianCenteredDataDetected5[,69]~PCgMedianCenteredDataDetected5With2SubjFixed[,69])
plot(PCgMedianCenteredDataDetected5[,69]~PCgMedianCenteredDataDetected5With2SubjFixed[,68])


PCgSubjects<-colnames(PCgMedianCenteredDataDetected5With2SubjFixed)
PCGSubjectIllnessTF<-SubjectsWithIllness[,1]%in% PCgSubjects
sum(PCGSubjectIllnessTF)
[1] 85
PCgSubjectsIllnessCategories<-SubjectsWithIllness[PCGSubjectIllnessTF==TRUE,]
dim(PCgSubjectsIllnessCategories)
[1] 85  9


SCgSubjects<-colnames(SCgMedianCenteredDataDetected5)
SCgSubjectIllnessTF<-SubjectsWithIllness[,1]%in%SCgSubjects
sum(SCgSubjectIllnessTF)
##[1] 86

SCgSubjectsIllnessCategories<-SubjectsWithIllness[SCgSubjectIllnessTF==TRUE,]
dim(SCgSubjectsIllnessCategories)
[1] 86  9


CBSubjects<-CBSubjectsDemographicsOrdered3[,6]
CBSubjectIllnessTF<-SubjectsWithIllness[,1]%in%CBSubjects
sum(CBSubjectIllnessTF)
[1] 111
CBSubjectsIllnessCategories<-SubjectsWithIllness[CBSubjectIllnessTF==TRUE,]
dim(CBSubjectsIllnessCategories)
[1] 111   9

DlpfcSubjects<-DlpfcSubjectsDemographicsOrdered3[,6]
DlpfcSubjectIllnessTF<-SubjectsWithIllness[,1]%in%DlpfcSubjects
sum(DlpfcSubjectIllnessTF)
[1] 108
DlpfcSubjectsIllnessCategories<-SubjectsWithIllness[DlpfcSubjectIllnessTF==TRUE,]
dim(DlpfcSubjectsIllnessCategories)
[1] 108   9

NaccSubjects<-NaccSubjectsDemographicsOrdered3[,6]
NaccSubjectIllnessTF<-SubjectsWithIllness[,1]%in%NaccSubjects
sum(NaccSubjectIllnessTF)
[1] 105
NaccSubjectsIllnessCategories<-SubjectsWithIllness[NaccSubjectIllnessTF==TRUE,]
dim(NaccSubjectsIllnessCategories)
[1] 105   9

AmySubjects<-AmySubjectsDemographicsOrdered3[,6]
AmySubjectIllnessTF<-SubjectsWithIllness[,1]%in%AmySubjects
sum(AmySubjectIllnessTF)
[1] 97
AmySubjectsIllnessCategories<-SubjectsWithIllness[AmySubjectIllnessTF==TRUE,]
dim(AmySubjectsIllnessCategories)
[1] 97  9

ACgSubjects<-ACgSubjectsDemographicsOrdered3[,6]
ACgSubjectIllnessTF<-SubjectsWithIllness[,1]%in%ACgSubjects
sum(ACgSubjectIllnessTF)
[1] 109
ACgSubjectsIllnessCategories<-SubjectsWithIllness[ACgSubjectIllnessTF==TRUE,]
dim(ACgSubjectsIllnessCategories)
[1] 109   9

HCSubjects<-HCSubjectsDemographicsOrdered3[,6]
HCSubjectIllnessTF<-SubjectsWithIllness[,1]%in%HCSubjects
sum(HCSubjectIllnessTF)
[1] 109
HCSubjectsIllnessCategories<-SubjectsWithIllness[HCSubjectIllnessTF==TRUE,]
dim(HCSubjectsIllnessCategories)
[1] 109   9

str(HC_Annotated_DetectedProbeData)
dim(SCgDetPvalDetected)
str(SCgDetPvalDetected)
SCgDetPvalDetected[1:1000,1]
table(SCgDetPvalDetected[,1])

str(SCg86DataImported3MedCtrDataDetected)
row.names(SCg86DataImported3MedCtrDataDetected)

table(PCgDetPvalDetected[,1])
  ACg   Amy aThal    CB Dlpfc    HC mThal  Nacc   PCg   SCg 
    0     0     0     0     0     0     0     0 13884     0 
table(aThalDetPvalDetected[,1])
  ACg   Amy aThal    CB Dlpfc    HC mThal  Nacc   PCg   SCg 
    0     0 13771     0     0     0     0     0     0     0 
table(mThalDetPvalDetected[,1])
  ACg   Amy aThal    CB Dlpfc    HC mThal  Nacc   PCg   SCg 
    0     0     0     0     0     0 13690     0     0     0 
table(SCgDetPvalDetected[,1])
  ACg   Amy aThal    CB Dlpfc    HC mThal  Nacc   PCg   SCg 
    0     0     0     0     0     0     0     0     0 13919 

HCtest<-HCDetPvalProbeOrdered[,3]<0.055

PCgtest<-(PCgDetPvalProbeOrdered[,3]<0.055)
PCgtest2<-PCgDetPvalDetected[,1]
PCgtest==PCgtest3
length(PCgtest)
length(PCgtest2)
sum(PCgDetPvalProbeOrdered[,3]<0.055)
PCgtest3<-PCgtest==TRUE

str(HCSubjectsDemographicsOrder)
str(HCSubjectsIllnessCategories)
HCSubjectsDemographicsOrder[,6]== HCSubjectsIllnessCategories[,1]
str(HCSubjectsDemographicsOrdered3)

HCSubjectIllnessCategories2<-HCSubjectsIllnessCategories[order(HCSubjectsIllnessCategories[,1]),]
HCSubjectsDemographicsOrdered3[,6]==HCSubjectIllnessCategories2[,1]
HCSubjectsDemographicsOrdered3[,6]
HCSubjectsIllnessCategories[,1]
colnames(HCSubjectsDemographicsOrdered3)
dim(HCSubjectsDemographicsOrdered3)
##[1] 109  32
dim(HCSubjectIllnessCategories2)
##[1] 109   9

DlpfcSubjectsDemographicsOrdered3[,6]== DlpfcSubjectsIllnessCategories[,1]
DlpfcSubjectIllnessCategories2<-DlpfcSubjectsIllnessCategories[order(DlpfcSubjectsIllnessCategories[,1]),]
DlpfcSubjectsDemographicsOrdered3[,6]==DlpfcSubjectIllnessCategories2[,1]
dim(DlpfcSubjectsDemographicsOrdered3)
##	[1] 109  32
dim(DlpfcSubjectIllnessCategories2)
##	[1] 108   9

NaccSubjectsDemographicsOrdered3[,6]== NaccSubjectsIllnessCategories[,1]
NaccSubjectsIllnessCategories2<-NaccSubjectsIllnessCategories[order(NaccSubjectsIllnessCategories[,1]),]
dim(NaccSubjectsDemographicsOrdered3)
##	[1] 106  32
dim(NaccSubjectsIllnessCategories2)
##	[1] 105   9

AmySubjectsDemographicsOrdered3[,6]== AmySubjectsIllnessCategories[,1]
AmySubjectsIllnessCategories2<-AmySubjectsIllnessCategories[order(AmySubjectsIllnessCategories[,1]),]
dim(AmySubjectsDemographicsOrdered3)
##	[1] 98 32
dim(AmySubjectsIllnessCategories)
##	[1] 97  9



CBSubjectsDemographicsOrdered3[,6]== CBSubjectsIllnessCategories[,1]
CBSubjectsIllnessCategories2<-CBSubjectsIllnessCategories[order(CBSubjectsIllnessCategories[,1]),]
dim(CBSubjectsDemographicsOrdered3)
##	[1] 113  32
dim(CBSubjectsIllnessCategories)
##	[1] 111   9




ACgSubjectsDemographicsOrdered3[,6]== ACgSubjectsIllnessCategories[,1]
##FALSE
dim(ACgSubjectsDemographicsOrdered3)
##	[1] 109  32
dim(ACgSubjectsIllnessCategories)
##	[1] 109   9

aThalSubjectsDemographicsOrdered3[,6]== aThalSubjectsIllnessCategories[,1]
dim(aThalSubjectsDemographicsOrdered3)
dim(aThalSubjectsIllnessCategories)




CBSubjectsDemographicsOrdered3[,6]== CBSubjectsIllnessCategories[,1]
dim(mThalSubjectsDemographicsOrdered3)
dim(mThalSubjectsIllnessCategories)





CBSubjectsDemographicsOrdered3[,6]== CBSubjectsIllnessCategories[,1]
dim(PCgSubjectsDemographicsOrdered3)
dim(PCgSubjectsIllnessCategories)


CBSubjectsDemographicsOrdered3[,6]== CBSubjectsIllnessCategories[,1]
dim(SCgSubjectsDemographicsOrdered3)
dim(SCgSubjectsIllnessCategories)


library(doBy)
PCgSubsinLargerSubsTF<-SubjectsDemographicsOrdered644[,6]%in%PCgSubjects
PCgSubjectsfromSubjectsDemographicsOrdered644<-SubjectsDemographicsOrdered644[PCgSubsinLargerSubsTF==TRUE,]
PCgUnique<-PCgSubjectsfromSubjectsDemographicsOrdered644[firstobs(PCgSubjectsfromSubjectsDemographicsOrdered644[,6]),]
dim(PCgUnique)
##[1] 85 32
SCgSubsinLargerSubsTF<-SubjectsDemographicsOrdered644[,6]%in%SCgSubjects
SCgSubjectsfromSubjectsDemographicsOrdered644<-SubjectsDemographicsOrdered644[SCgSubsinLargerSubsTF==TRUE,]
SCgUnique<-SCgSubjectsfromSubjectsDemographicsOrdered644[firstobs(SCgSubjectsfromSubjectsDemographicsOrdered644[,6]),]
dim(SCgUnique)
##[1] 86 32
aThalSubsinLargerSubsTF<-SubjectsDemographicsOrdered644[,6]%in%aThalSubjects
aThalSubjectsfromSubjectsDemographicsOrdered644<-SubjectsDemographicsOrdered644[aThalSubsinLargerSubsTF==TRUE,]
aThalUnique<-aThalSubjectsfromSubjectsDemographicsOrdered644[firstobs(aThalSubjectsfromSubjectsDemographicsOrdered644[,6]),]
dim(aThalUnique)
##[1] 74 32
mThalSubsinLargerSubsTF<-SubjectsDemographicsOrdered644[,6]%in%mThalSubjects
mThalSubjectsfromSubjectsDemographicsOrdered644<-SubjectsDemographicsOrdered644[mThalSubsinLargerSubsTF==TRUE,]
mThalUnique<-mThalSubjectsfromSubjectsDemographicsOrdered644[firstobs(mThalSubjectsfromSubjectsDemographicsOrdered644[,6]),]
dim(mThalUnique)
##[1] 77 32

PCgUnique[,6]==PCgSubjectsIllnessCategories[,1]
PCgSubjectIllnessCategories2<-PCgSubjectsIllnessCategories[order(PCgSubjectsIllnessCategories[,1]),]
PCgUnique[,6]==PCgSubjectIllnessCategories2[,1]
str(PCgUnique)
str(PCgSubjectIllnessCategories2)
PCgUniqueOrdered<-PCgUnique[order(PCgUnique[,6]),]
cbind(PCgUniqueOrdered[,6],PCgSubjectIllnessCategories2[,1])
PCgUniqueOrdered[,6]==PCgSubjectIllnessCategories2[,1]

SCgSubjectIllnessCategories2<-SCgSubjectsIllnessCategories[order(SCgSubjectsIllnessCategories[,1]),]
SCgUniqueOrdered<-SCgUnique[order(SCgUnique[,6]),]
SCgUniqueOrdered[,6]==SCgSubjectIllnessCategories2[,1]
dim(SCgUniqueOrdered)
##[1] 86 32

mThalSubjectIllnessCategories2<-mThalSubjectsIllnessCategories[order(mThalSubjectsIllnessCategories[,1]),]
mThalUniqueOrdered<-mThalUnique[order(mThalUnique[,6]),]
mThalUniqueOrdered[,6]==mThalSubjectIllnessCategories2[,1]
dim(mThalUniqueOrdered)
##[1] 77 32

aThalSubjectIllnessCategories2<-aThalSubjectsIllnessCategories[order(aThalSubjectsIllnessCategories[,1]),]
aThalUniqueOrdered<-aThalUnique[order(aThalUnique[,6]),]
aThalUniqueOrdered[,6]==aThalSubjectIllnessCategories2[,1]
dim(aThalUniqueOrdered)
##[1] 74 32

colnames(aThalUniqueOrdered)
colnames(aThalSubjectIllnessCategories2)
athalSubjectDemogrphaicsforAnalysis<-cbind(aThalUniqueOrdered,aThalSubjectIllnessCategories2[c(4,5,6,8,1)])
dim(athalSubjectDemogrphaicsforAnalysis)
##[1] 74 37
colnames(athalSubjectDemogrphaicsforAnalysis)

mthalSubjectDemogrphaicsforAnalysis<-cbind(mThalUniqueOrdered,mThalSubjectIllnessCategories2[c(4,5,6,8,1)])
dim(mthalSubjectDemogrphaicsforAnalysis)
##[1] 77 37
colnames(mthalSubjectDemogrphaicsforAnalysis)
colnames(mthalSubjectDemogrphaicsforAnalysis)==colnames(athalSubjectDemogrphaicsforAnalysis)

PCgSubjectDemogrphaicsforAnalysis<-cbind(PCgUniqueOrdered,PCgSubjectIllnessCategories2[c(4,5,6,8,1)])
dim(PCgSubjectDemogrphaicsforAnalysis)
##[1] 85 37
colnames(PCgSubjectDemogrphaicsforAnalysis)
colnames(PCgSubjectDemogrphaicsforAnalysis)==colnames(athalSubjectDemogrphaicsforAnalysis)

SCgSubjectDemogrphaicsforAnalysis<-cbind(SCgUniqueOrdered,SCgSubjectIllnessCategories2[c(4,5,6,8,1)])
dim(SCgSubjectDemogrphaicsforAnalysis)
##[1] 86 37
colnames(SCgSubjectDemogrphaicsforAnalysis)
colnames(SCgSubjectDemogrphaicsforAnalysis)==colnames(athalSubjectDemogrphaicsforAnalysis)




##	Creating files for Stats easily named

SCgStats<-SCg86DataImported3MedCtrDataDetected
PCgStats<-PCg85DataImported3MedCtrDataDetected
aThalStats<-aThal74DataImported3MedCtrDataDetected
mThalStats<-mThal77DataImported3MedCtrDataDetected
HCStats<-HCMedianCenteredDataDetected4
AmyStats<-AmyMedianCenteredDataDetected4
ACgStats<-ACgMedianCenteredDataDetected4
CBStats<-CBMedianCenteredDataDetected4
DlpfcStats<-DlpfcMedianCenteredDataDetected4
NaccStats<-NaccMedianCenteredDataDetected4

#Fun stuff: Characterizing demographics for each region


png("18 Gender Check.png")

plot(normNoOutliers[GeneIdentifiers_filtered[,1]=="XIST",]~normNoOutliers[GeneIdentifiers_filtered[,1]=="RPS4Y1",], xlab="Male Gene: RPS4Y1", ylab="Female Gene: XIST", main = "Does subject gender match gendered gene expression?")
points(normNoOutliers[GeneIdentifiers_filtered[,1]=="XIST", GenderNoOutliers=="M"]~normNoOutliers[GeneIdentifiers_filtered[,1]=="RPS4Y1", GenderNoOutliers=="M"], col=2)
points(normNoOutliers[GeneIdentifiers_filtered[,1]=="XIST", GenderNoOutliers=="F"]~normNoOutliers[GeneIdentifiers_filtered[,1]=="RPS4Y1", GenderNoOutliers=="F"], col=3)
egend(min(normNoOutliers[GeneIdentifiers_filtered[,1]=="RPS4Y1",])+0.5, min(normNoOutliers[GeneIdentifiers_filtered[,1]=="XIST",])+0.5, c("M", "F"), text.col=c(2,3), pch=19, col=c(2,3))
dev.off()

str(PCg85DataImported3)
str(PCgStats)
row.names(PCgStats)
row.names(PCg85DataImported3)

plot(as.numeric(PCgStats[PCgdetectedIlluminaDecodeSorted22177Abbreviated[,12]=="RPS4Y1",])~PCgSubjectDemogrphaicsforAnalysis[,7])
sum(PCgdetectedIlluminaDecodeSorted22177Abbreviated[,12]=="RPS4Y1")
sum(PCgdetectedIlluminaDecodeSorted22177Abbreviated[,12]=="CD99")
#######################################################################################
#compiling a list of *strongly* y chromosome related genes:
Ychromosomegenes<-c("RPS4Y1", "DDX3Y", "EIF1AY", "USP9Y", "TTTY15", "UTY", "NLGN4Y")

#extracting signal for y-chromosome related genes by creating a true-false statement for which probes on the decode list are found in the y-chromosome gene list. 
PCG_Ychromosomesignal<-PCgStats[PCgdetectedIlluminaDecodeSorted22177Abbreviated[,12]%in%Ychromosomegenes,]

#extracting the gene names for y-chromosome related genes by creating a true-false statement for which probes on the decode list are found in the y-chromosome gene list.
PCG_Ychromosomegenenames<-PCgdetectedIlluminaDecodeSorted22177Abbreviated[PCgdetectedIlluminaDecodeSorted22177Abbreviated[,12]%in%Ychromosomegenes, 12]

#made the gene names for the y chromosome genes the row-names for the y-chromosome signal data:
row.names(PCG_Ychromosomesignal)<-PCG_Ychromosomegenenames

PCGychromosomeIndex<-apply(PCG_Ychromosomesignal, 2, mean)

png("PCG_Histogram_YChromosomeIndex.png")
hist(PCGychromosomeIndex, col=2)
dev.off()

#Number of genetic females - cut-off determined by separation in bimodal histogram:
sum(PCGychromosomeIndex<(-1))
[1] 18

sum(PCgSubjectDemogrphaicsforAnalysis[,7]=="F")
[1] 18

#The number of females in the sample matches the number of genetic females.

png("PCG_Boxplot_YChromosomeIndex.png")
boxplot(PCGychromosomeIndex~PCgSubjectDemogrphaicsforAnalysis[,7])
dev.off()
#... but it looks like two of the samples are switched. And there is one super male.



#bound together the gender information with the y-chromosome signal data just to make comparisons easy.
temp<-rbind(PCgSubjectDemogrphaicsforAnalysis[,7], PCG_Ychromosomesignal)

#wrote the gender infor + y-chromosome data out to an excel file where we could easily visualize the data using conditional formatting (i.e. pretty colors)
write.csv(temp, "PCG_GenderCheck.csv")

table(PCgSubjectDemogrphaicsforAnalysis[,7])
levels(PCgSubjectDemogrphaicsforAnalysis[,7])
[1] "F" "M"
#########SCg  No problem with Gender

Ychromosomegenes<-c("RPS4Y1", "DDX3Y", "EIF1AY", "USP9Y", "TTTY15", "UTY", "NLGN4Y")

SCG_Ychromosomesignal<-SCgStats[SCgdetectedIlluminaDecodeSorted22177Abbreviated[,12]%in%Ychromosomegenes,]

SCG_Ychromosomegenenames<-SCgdetectedIlluminaDecodeSorted22177Abbreviated[SCgdetectedIlluminaDecodeSorted22177Abbreviated[,12]%in%Ychromosomegenes, 12]

row.names(SCG_Ychromosomesignal)<-SCG_Ychromosomegenenames

temp<-rbind(SCgSubjectDemogrphaicsforAnalysis[,7], SCG_Ychromosomesignal)
write.csv(temp, "SCG_GenderCheck.csv")

table(SCgSubjectDemogrphaicsforAnalysis[,7])
# F  M 
18 68 
levels(SCgSubjectDemogrphaicsforAnalysis[,7])

###########mThal--Only 1 off with Gender

Ychromosomegenes<-c("RPS4Y1", "DDX3Y", "EIF1AY", "USP9Y", "TTTY15", "UTY", "NLGN4Y")

mThal_Ychromosomesignal<-mThalStats[mThaldetectedIlluminaDecodeSorted22177Abbreviated[,12]%in%Ychromosomegenes,]

mThal_Ychromosomegenenames<-mThaldetectedIlluminaDecodeSorted22177Abbreviated[mThaldetectedIlluminaDecodeSorted22177Abbreviated[,12]%in%Ychromosomegenes, 12]

row.names(mThal_Ychromosomesignal)<-mThal_Ychromosomegenenames

temp<-rbind(mthalSubjectDemogrphaicsforAnalysis[,7], mThal_Ychromosomesignal)
write.csv(temp, "mThal_GenderCheck.csv")

table(mthalSubjectDemogrphaicsforAnalysis[,7])
## F  M 
16 61 
levels(mthalSubjectDemogrphaicsforAnalysis[,7])

###########aThal  No problem with Gender

Ychromosomegenes<-c("RPS4Y1", "DDX3Y", "EIF1AY", "USP9Y", "TTTY15", "UTY", "NLGN4Y")

aThal_Ychromosomesignal<-aThalStats[aThaldetectedIlluminaDecodeSorted22177Abbreviated[,12]%in%Ychromosomegenes,]

aThal_Ychromosomegenenames<-aThaldetectedIlluminaDecodeSorted22177Abbreviated[aThaldetectedIlluminaDecodeSorted22177Abbreviated[,12]%in%Ychromosomegenes, 12]

row.names(aThal_Ychromosomesignal)<-aThal_Ychromosomegenenames

temp<-rbind(athalSubjectDemogrphaicsforAnalysis[,7], aThal_Ychromosomesignal)
write.csv(temp, "aThal_GenderCheck.csv")

table(athalSubjectDemogrphaicsforAnalysis[,7])
levels(athalSubjectDemogrphaicsforAnalysis[,7])

##ACg
Ychromosomegenes<-c("RPS4Y1", "DDX3Y", "EIF1AY", "USP9Y", "TTTY15", "UTY", "NLGN4Y")

ACg_Ychromosomesignal<-ACgStats[ACgdetectedIlluminaDecodeSorted22177Abbreviated[,12]%in%Ychromosomegenes,]

ACg_Ychromosomegenenames<-ACgdetectedIlluminaDecodeSorted22177Abbreviated[ACgdetectedIlluminaDecodeSorted22177Abbreviated[,12]%in%Ychromosomegenes, 12]

row.names(ACg_Ychromosomesignal)<-ACg_Ychromosomegenenames

temp<-rbind(ACgSubjectsDemographicsOrdered3[,7], ACg_Ychromosomesignal)
write.csv(temp, "ACg_GenderCheck.csv")

table(ACgSubjectsDemographicsOrdered3[,7])
## F  M 
25 84 
levels(ACgSubjectsDemographicsOrdered3[,7])

##Amy
Ychromosomegenes<-c("RPS4Y1", "DDX3Y", "EIF1AY", "USP9Y", "TTTY15", "UTY", "NLGN4Y")

AMY_Ychromosomesignal<-AmyStats[AmydetectedIlluminaDecodeSorted22177Abbreviated[,12]%in%Ychromosomegenes,]

AMY_Ychromosomegenenames<-AmydetectedIlluminaDecodeSorted22177Abbreviated[AmydetectedIlluminaDecodeSorted22177Abbreviated[,12]%in%Ychromosomegenes, 12]

row.names(AMY_Ychromosomesignal)<-AMY_Ychromosomegenenames

temp<-rbind(AmySubjectsDemographicsOrdered3[,7], AMY_Ychromosomesignal)
write.csv(temp, "AMY_GenderCheck.csv")

table(AmySubjectsDemographicsOrdered3[,7])
## F  M 
23 75 

levels(AmySubjectsDemographicsOrdered3[,7])

##HC
Ychromosomegenes<-c("RPS4Y1", "DDX3Y", "EIF1AY", "USP9Y", "TTTY15", "UTY", "NLGN4Y")

HC_Ychromosomesignal<-HCStats[HCdetectedIlluminaDecodeSorted22177Abbreviated[,12]%in%Ychromosomegenes,]

HC_Ychromosomegenenames<-HCdetectedIlluminaDecodeSorted22177Abbreviated[HCdetectedIlluminaDecodeSorted22177Abbreviated[,12]%in%Ychromosomegenes, 12]

row.names(HC_Ychromosomesignal)<-HC_Ychromosomegenenames

temp<-rbind(HCSubjectsDemographicsOrdered3[,7], HC_Ychromosomesignal)
write.csv(temp, "HC_GenderCheck.csv")

table(HCSubjectsDemographicsOrdered3[,7])
## F  M 
25 84 
levels(HCSubjectsDemographicsOrdered3[,7])

##CB
Ychromosomegenes<-c("RPS4Y1", "DDX3Y", "EIF1AY", "USP9Y", "TTTY15", "UTY", "NLGN4Y")

CB_Ychromosomesignal<-CBStats[CBIlluminaDecodeSorted22177Abbreviated[,12]%in%Ychromosomegenes,]

CB_Ychromosomegenenames<-CBIlluminaDecodeSorted22177Abbreviated[CBIlluminaDecodeSorted22177Abbreviated[,12]%in%Ychromosomegenes, 12]

row.names(CB_Ychromosomesignal)<-CB_Ychromosomegenenames

temp<-rbind(CBSubjectsDemographicsOrdered3[,7], CB_Ychromosomesignal)
write.csv(temp, "CB_GenderCheck.csv")

table(CBSubjectsDemographicsOrdered3[,7])
## F  M 
25 88 
levels(CBSubjectsDemographicsOrdered3[,7])

##Dlpfc  These are WHACKY!
Ychromosomegenes<-c("RPS4Y1", "DDX3Y", "EIF1AY", "USP9Y", "TTTY15", "UTY", "NLGN4Y")

Dlpfc_Ychromosomesignal<-DlpfcStats[DlpfcIlluminaDecodeSorted22177Abbreviated[,12]%in%Ychromosomegenes,]

Dlpfc_Ychromosomegenenames<-DlpfcIlluminaDecodeSorted22177Abbreviated[DlpfcIlluminaDecodeSorted22177Abbreviated[,12]%in%Ychromosomegenes, 12]

row.names(Dlpfc_Ychromosomesignal)<-Dlpfc_Ychromosomegenenames

temp<-rbind(DlpfcSubjectsDemographicsOrdered3[,7], Dlpfc_Ychromosomesignal)
write.csv(temp, "Dlpfc_GenderCheck.csv")

table(DlpfcSubjectsDemographicsOrdered3[,7])
## F  M 
24 85 
levels(DlpfcSubjectsDemographicsOrdered3[,7])

##Nacc  Looks WHACKY!
Ychromosomegenes<-c("RPS4Y1", "DDX3Y", "EIF1AY", "USP9Y", "TTTY15", "UTY", "NLGN4Y")

Nacc_Ychromosomesignal<-NaccStats[NaccIlluminaDecodeSorted22177Abbreviated[,12]%in%Ychromosomegenes,]

Nacc_Ychromosomegenenames<-NaccIlluminaDecodeSorted22177Abbreviated[NaccIlluminaDecodeSorted22177Abbreviated[,12]%in%Ychromosomegenes, 12]

row.names(Nacc_Ychromosomesignal)<-Nacc_Ychromosomegenenames

temp<-rbind(NaccSubjectsDemographicsOrdered3[,7], Nacc_Ychromosomesignal)
write.csv(temp, "Nacc_GenderCheck.csv")

table(NaccSubjectsDemographicsOrdered3[,7])
## F  M 
24 82 
levels(NaccSubjectsDemographicsOrdered3[,7])































plot(as.numeric(PCgStats[PCgdetectedIlluminaDecodeSorted22177Abbreviated[,12]=="CD99",])~PCgSubjectDemogrphaicsforAnalysis[,7])
##Y SPecific
plot(as.numeric(PCgStats[PCgdetectedIlluminaDecodeSorted22177Abbreviated[,12]=="MSN",])~PCgSubjectDemogrphaicsforAnalysis[,7])

plot(as.numeric(PCgStats[PCgdetectedIlluminaDecodeSorted22177Abbreviated[,12]=="NLGN4Y",])~PCgSubjectDemogrphaicsforAnalysis[,7])

plot(as.numeric(PCgStats[PCgdetectedIlluminaDecodeSorted22177Abbreviated[,12]=="NLGN4X",])~PCgSubjectDemogrphaicsforAnalysis[,7])

plot(as.numeric(PCgStats[PCgdetectedIlluminaDecodeSorted22177Abbreviated[,12]=="DAXX",])~PCgSubjectDemogrphaicsforAnalysis[,7])

plot(as.numeric(PCgStats[PCgdetectedIlluminaDecodeSorted22177Abbreviated[,12]=="EIF1B",])~PCgSubjectDemogrphaicsforAnalysis[,7])

plot(as.numeric(PCgStats[PCgdetectedIlluminaDecodeSorted22177Abbreviated[,12]=="GTPBP6",])~PCgSubjectDemogrphaicsforAnalysis[,7])

plot(as.numeric(PCgStats[PCgdetectedIlluminaDecodeSorted22177Abbreviated[,12]=="IL3RA",])~PCgSubjectDemogrphaicsforAnalysis[,7])

plot(as.numeric(PCgStats[PCgdetectedIlluminaDecodeSorted22177Abbreviated[,12]=="JARID1D",])~PCgSubjectDemogrphaicsforAnalysis[,7])

plot(as.numeric(PCgStats[PCgdetectedIlluminaDecodeSorted22177Abbreviated[,12]=="NLK",])~PCgSubjectDemogrphaicsforAnalysis[,7])

plot(as.numeric(PCgStats[PCgdetectedIlluminaDecodeSorted22177Abbreviated[,12]=="PRKY",])~PCgSubjectDemogrphaicsforAnalysis[,7])

plot(as.numeric(PCgStats[PCgdetectedIlluminaDecodeSorted22177Abbreviated[,12]=="RPS4X",])~PCgSubjectDemogrphaicsforAnalysis[,7])

plot(as.numeric(PCgStats[PCgdetectedIlluminaDecodeSorted22177Abbreviated[,12]=="SPRY3",])~PCgSubjectDemogrphaicsforAnalysis[,7])

plot(as.numeric(PCgStats[PCgdetectedIlluminaDecodeSorted22177Abbreviated[,12]=="USP9Y",])~PCgSubjectDemogrphaicsforAnalysis[,7])

plot(as.numeric(PCgStats[PCgdetectedIlluminaDecodeSorted22177Abbreviated[,12]=="UTX",])~PCgSubjectDemogrphaicsforAnalysis[,7])

plot(as.numeric(PCgStats[PCgdetectedIlluminaDecodeSorted22177Abbreviated[,12]=="UTY",])~PCgSubjectDemogrphaicsforAnalysis[,7])

plot(as.numeric(SCgStats[SCgdetectedIlluminaDecodeSorted22177Abbreviated[,12]=="UTX",])~SCgSubjectDemogrphaicsforAnalysis[,7])

plot(as.numeric(SCgStats[SCgdetectedIlluminaDecodeSorted22177Abbreviated[,12]=="UTY",])~SCgSubjectDemogrphaicsforAnalysis[,7])

plot(as.numeric(mThalStats[mThaldetectedIlluminaDecodeSorted22177Abbreviated[,12]=="UTX",])~ mthalSubjectDemogrphaicsforAnalysis[,7])

plot(as.numeric(mThalStats[mThaldetectedIlluminaDecodeSorted22177Abbreviated[,12]=="UTY",])~ mthalSubjectDemogrphaicsforAnalysis[,7])

plot(as.numeric(aThalStats[aThaldetectedIlluminaDecodeSorted22177Abbreviated[,12]=="UTX",])~ athalSubjectDemogrphaicsforAnalysis[,7])

plot(as.numeric(aThalStats[aThaldetectedIlluminaDecodeSorted22177Abbreviated[,12]=="UTY",])~ athalSubjectDemogrphaicsforAnalysis[,7])

plot(as.numeric(AmyStats[AmydetectedIlluminaDecodeSorted22177Abbreviated[,12]=="UTY",])~ AmySubjectsDemographicsOrdered3[,7])

plot(as.numeric(AmyStats[AmydetectedIlluminaDecodeSorted22177Abbreviated[,12]=="UTX",])~ AmySubjectsDemographicsOrdered3[,7])

plot(as.numeric(HCStats[HCdetectedIlluminaDecodeSorted22177Abbreviated[,12]=="UTY",])~ HCSubjectsDemographicsOrdered3[,7])

plot(as.numeric(HCStats[HCdetectedIlluminaDecodeSorted22177Abbreviated[,12]=="UTX",])~ HCSubjectsDemographicsOrdered3[,7])

plot(as.numeric(ACgStats[ACgdetectedIlluminaDecodeSorted22177Abbreviated[,12]=="UTY",])~ ACgSubjectsDemographicsOrdered3[,7])

plot(as.numeric(ACgStats[ACgdetectedIlluminaDecodeSorted22177Abbreviated[,12]=="UTX",])~ ACgSubjectsDemographicsOrdered3[,7])

plot(as.numeric(DlpfcStats[DlpfcIlluminaDecodeSorted22177Abbreviated[,12]=="UTY",])~ DlpfcSubjectsDemographicsOrdered3[,7])

plot(as.numeric(DlpfcStats[DlpfcIlluminaDecodeSorted22177Abbreviated[,12]=="UTX",])~ DlpfcSubjectsDemographicsOrdered3[,7])

plot(as.numeric(CBStats[CBIlluminaDecodeSorted22177Abbreviated[,12]=="UTY",])~ CBSubjectsDemographicsOrdered3[,7])

plot(as.numeric(CBStats[CBIlluminaDecodeSorted22177Abbreviated[,12]=="UTX",])~ CBSubjectsDemographicsOrdered3[,7])

plot(as.numeric(NaccStats[NaccIlluminaDecodeSorted22177Abbreviated[,12]=="UTY",])~ NaccSubjectsDemographicsOrdered3[,7])

plot(as.numeric(NaccStats[NaccIlluminaDecodeSorted22177Abbreviated[,12]=="UTX",])~ NaccSubjectsDemographicsOrdered3[,7])


plot(as.numeric(SCgStats[SCgdetectedIlluminaDecodeSorted22177Abbreviated[,12]=="HLA-C",])~SCgSubjectDemogrphaicsforAnalysis[,7])

table(PCgSubjectDemogrphaicsforAnalysis[,7])
## F  M 
##18 67 

######################################################
##FARTING AROUND TRYING TO MAKE SENSE OF IT ALL
str(SCgSubjectDemogrphaicsforAnalysis)
dim(SCgSubjectDemogrphaicsforAnalysis)
SCgSubjectDemogrphaicsforAnalysis[,c(6,7,8,3)]
##Column 6= Subject ID, 7=Sex, 8=Age, 3=Disease
table(SCgSubjectDemogrphaicsforAnalysis[,7])
## F  M 
##18 68 

SCgCheckM<-SCgSubjectDemogrphaicsforAnalysis[SCgSubjectDemogrphaicsforAnalysis[,7]=="M",]
SCgCheckM[,c(6,7,8,3)]

SCgCheckF<-SCgSubjectDemogrphaicsforAnalysis[SCgSubjectDemogrphaicsforAnalysis[,7]=="F",]
SCgCheckF[,c(6,7,8,3)]

PCgSubjectDemogrphaicsforAnalysis[,c(6,7,8,3)]
table(PCgSubjectDemogrphaicsforAnalysis[,7])
## F  M 
##18 67 
table(SCgSubjectDemogrphaicsforAnalysis[,7])
## F  M 
##18 68 

str(PCgSubjectDemogrphaicsforAnalysis)
PCgCheckM<-PCgSubjectDemogrphaicsforAnalysis[PCgSubjectDemogrphaicsforAnalysis[,7]=="M",]
PCgCheckM[,c(6,7,8,3)]

PCgCheckF<-PCgSubjectDemogrphaicsforAnalysis[PCgSubjectDemogrphaicsforAnalysis[,7]=="F",]
PCgCheckF[,c(6,7,8,3)]

SCgCheckF[,6]==PCgCheckF[,6]

plot(as.numeric(SCgStats[SCgdetectedIlluminaDecodeSorted22177Abbreviated[,12]=="UTY",])~SCgSubjectDemogrphaicsforAnalysis[,7])

max(SCgStats[SCgdetectedIlluminaDecodeSorted22177Abbreviated[,12]=="UTX",])
##[1] 0.4831073
min(SCgStats[SCgdetectedIlluminaDecodeSorted22177Abbreviated[,12]=="UTX",])
##[1] -0.502373

doesthiswork<-SCgStats[SCgdetectedIlluminaDecodeSorted22177Abbreviated[,12]=="UTX",]
doesthisworkm<-doesthiswork[SCgSubjectDemogrphaicsforAnalysis[,7]=="M"]
doesthisworkm

######################################################

plot(as.numeric(PCgStats[PCgdetectedIlluminaDecodeSorted22177Abbreviated[,12]=="CYorf15A",])~PCgSubjectDemogrphaicsforAnalysis[,7])



plot(as.numeric(PCg85DataImported3[PCgdetectedIlluminaDecodeSorted22177Abbreviated[,12]=="SRY",])~PCgSubjectDemogrphaicsforAnalysis[,7])

hist(as.numeric(PCgStats[PCgdetectedIlluminaDecodeSorted22177Abbreviated[,12]=="RPS4Y1",]))
hist(as.numeric(PCg85DataImported3[PCgdetectedIlluminaDecodeSorted22177Abbreviated[,12]=="SRY",]))
dim(PCg85DataImported3)

png("18 Gender Check.png")

plot(PCgStats[PCgdetectedIlluminaDecodeSorted22177Abbreviated[,12]=="UTX",]~ PCgStats[PCgdetectedIlluminaDecodeSorted22177Abbreviated[,12]=="UTY",], xlab="Male Gene: UTY", ylab="Female Gene: UTX", main = "Does subject gender match gendered gene expression?")

points(PCgStats[PCgdetectedIlluminaDecodeSorted22177Abbreviated[,12]=="UTX", PCgSubjectDemogrphaicsforAnalysis[,7]=="M"]~ PCgStats[PCgdetectedIlluminaDecodeSorted22177Abbreviated[,12]=="UTY", PCgSubjectDemogrphaicsforAnalysis[,7]=="M"], col=2)

points(PCgStats[PCgdetectedIlluminaDecodeSorted22177Abbreviated[,12]=="UTX", PCgSubjectDemogrphaicsforAnalysis[,7]=="F"]~ PCgStats[PCgdetectedIlluminaDecodeSorted22177Abbreviated[,12]=="UTY", PCgSubjectDemogrphaicsforAnalysis[,7]=="F"], col=3)

legend(min(PCgStats[PCgdetectedIlluminaDecodeSorted22177Abbreviated[,12]=="UTY",])+0.5, min(PCgStats[PCgdetectedIlluminaDecodeSorted22177Abbreviated[,12]=="UTX",])+0.5, c("M", "F"), text.col=c(2,3), pch=19, col=c(2,3))
dev.off()






















#Make a demographics folder for each brain region:

head(DlpfcSubjectsDemographicsOrdered)

str(DlpfcSubjectsDemographicsOrdered)

DlpfcSubjectsDemographicsOrderedCategory<-DlpfcSubjectsDemographicsOrdered[,c(3,7, 9, 11,12,13, 15,16,17,19,20,21,22,24,25,26, 28, 29, 30, 31, 32)]

colnames(DlpfcSubjectsDemographicsOrderedCategory)

DlpfcSubjectsDemographicsOrderedContinuous<-DlpfcSubjectsDemographicsOrdered[,c(8, 10,14, 18)]

colnames(DlpfcSubjectsDemographicsOrderedContinuous)

GeneralDemographicsCategorical<-file("GeneralDemographicsCategorical.txt")
out<-c(
capture.output(
for(i in 1:length(colnames(DlpfcSubjectsDemographicsOrderedCategory))){
print(colnames(DlpfcSubjectsDemographicsOrderedCategory)[i])
print(table(DlpfcSubjectsDemographicsOrderedCategory[,i]))
}
)
)
cat(out, file="GeneralDemographicsCategorical.txt", sep="\n", append=TRUE)
close(GeneralDemographicsCategorical)
rm(out)


DiagnosisDemographicsCategorical<-file("DiagnosisGeneralDemographicsCategorical.txt")
out<-c(
capture.output(
for(i in 1:length(colnames(DlpfcSubjectsDemographicsOrderedCategory))){
print(colnames(DlpfcSubjectsDemographicsOrderedCategory)[i])
print(table(DlpfcSubjectsDemographicsOrderedCategory[,i], DlpfcSubjectsDemographicsOrderedCategory[,1]))
}
)
)
cat(out, file="DiagnosisDemographicsCategorical.txt", sep="\n", append=TRUE)
close(DiagnosisDemographicsCategorical)
rm(out)




GeneralDemographicsContinuous<-file("GeneralDemographicsContinuous.txt")
out<-c(
capture.output(
for(i in 1:length(colnames(DlpfcSubjectsDemographicsOrderedContinuous))){
print(colnames(DlpfcSubjectsDemographicsOrderedContinuous)[i])
print(summary(DlpfcSubjectsDemographicsOrderedContinuous[,i]))
}
)
)
cat(out, file="GeneralDemographicsContinuous.txt", sep="\n", append=TRUE)
close(GeneralDemographicsContinuous)
rm(out)





i<-1

DiagnosisDemographicsContinuous<-file("DiagnosisDemographicsContinuous.txt")
out<-c(
capture.output(
for(i in 1:length(colnames(DlpfcSubjectsDemographicsOrderedContinuous))){
print(colnames(DlpfcSubjectsDemographicsOrderedContinuous)[i])
#We need to put the continuous variable and diagnosis in the same dataframe for the aggregate function to be happy.
temp<-data.frame(DlpfcSubjectsDemographicsOrderedContinuous[,i], DlpfcSubjectsDemographicsOrderedCategory[,1])
print("Mean by Diagnosis")
print(aggregate(temp[,1]~temp[,2], temp, mean))
print("StDev by Diagnosis")
print(aggregate(temp[,1]~temp[,2], temp, sd))
}
)
)
cat(out, file="DiagnosisDemographicsContinuous.txt", sep="\n", append=TRUE)
close(DiagnosisDemographicsContinuous)
rm(out)

`


SubjectContinuousVariables<-DlpfcSubjectsDemographicsOrderedContinuous
SubjectFactorVariables<-DlpfcSubjectsDemographicsOrderedCategory[,-12]

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


i<-2
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
`ContingencyTable<-table(SubjectFactorVariables[,i],SubjectFactorVariables[,j])
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


























































































































































#####Vikrim Coding as example below**********************************************************



#Extracts out the detection p-values for each subject and places them into a numeric matrix with appropriate row and column names:
Detection_Pvalue<-matrix(as.numeric(RawData[,seq(from=3, to=85, by=2)]), nrow=TotalNumberOfProbes, ncol=42)
colnames(Detection_Pvalue)<-colnames(RawData[, seq(from=3, to=85, by=2)])
row.names(Detection_Pvalue)<-row.names(RawData)

#Double checks whether the detection p-value matrix is the correct size and numeric:
if(length(Detection_Pvalue[1,])==!42){print("ERROR! The Detection_Pvalue matrix is the wrong size! Fix now or everything else will be wrong!")}else{print("Detection_Pvalue matrix is the correct size")}
if(is.numeric(Detection_Pvalue)==F){print("ERROR! The Detection_Pvalue matrix is not numeric! Fix now or everything else will be wrong!")}else{print("Detection_Pvalue is properly numeric. Hurray!")}

#Makes a matrix where the subject identifiers for the average signal and detection p-value matrices can be easily compared to make sure they are the same, and then outputs that data to a format that can be accessed by Excel:
DoubleCheckColumnNames<-cbind(colnames(AVG_Signal), colnames(Detection_Pvalue))
write.csv(DoubleCheckColumnNames, "04 DoubleCheckColumnNames.csv")

#Outputs the gene information for each of the probes (symbol, chromosome, full name, etc.)
GeneIdentifiers<-RawData[,c(1, c(86:91))]

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
Diagnosis<-relevel(Diagnosis, ref="Control")

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
points(PC1[Diagnosis=="Control"]~PC2[Diagnosis=="Control"], col=3)
points(PC1[Diagnosis=="MDD"]~PC2[Diagnosis=="MDD"], col=2)
legend(min(PC1), max(PC2)+10, c("Control", "MDD"), text.col=c(3, 2), pch=19, col=c(3, 2))
dev.off()

#Output a scatterplot illustrating the relationship between Principal components 3 & 4 (PC3 & PC4):
png("09 PC3 vs PC4.png")
plot(PC3~PC4, main="Principal Components Analysis of Normalized Filtered Data")
points(PC3[Diagnosis=="Control"]~PC4[Diagnosis=="Control"], col=3)
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
which("3686"==colnames(NormFiltered))

##2.  After coming up with a list of indices for the outliers, you can remove them by inserting their index numbers into this code, preceded by a minus sign (replacing -4, -6, and -13). If there aren't any outliers to remove, we will instead rename the original quantile normalized data to indicate that it has made it past the step of outlier removal:

###If there are outliers (example):
OutlierRemovalIndices<-c(-4, -6, -13) 

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
rownames(normNoOutliers)<<-rownames(log.avg_signal_filtered)
colnames(normNoOutliers)<<-colnames(log.avg_signal_filtered [ ,OutlierRemovalIndices])

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
plot(summary(pcaNormFilterednoOutliers)$importance[3,]~(c(1:length(summary(pcaNormFilterednoOutliers)$importance[2,]))))), main="Variance Explained by Each Principal Component", xlab="PC #", ylab="Cumulative Proportion of Variance Explained", col=3)
dev.off()

#Output a scatterplot illustrating the relationship between Principal components 1 & 2 (PC1 & PC2):
png("10 PC1 vs PC2 no outliers.png")
plot(PC1noOutliers~PC2noOutliers, main="Principal Components Analysis of Normalized Filtered Data No Outliers")
points(PC1noOutliers[DiagnosisNoOutliers=="Control"]~PC2noOutliers[DiagnosisNoOutliers=="Control"], col=3)
points(PC1noOutliers[DiagnosisNoOutliers=="MDD"]~PC2noOutliers[DiagnosisNoOutliers=="MDD"], col=2)
legend(min(PC1noOutliers), max(PC2noOutliers)+10, c("Control", "MDD"), text.col=c(3, 2), pch=19, col=c(3, 2))
dev.off()

#Output a scatterplot illustrating the relationship between Principal components 3 & 4 (PC3 & PC4):
png("10 PC3 vs PC4 no outliers.png")
plot(PC3noOutliers~PC4noOutliers, main="Principal Components Analysis of Normalized Filtered Data No Outliers")
points(PC3noOutliers[DiagnosisnoOutliers=="Control"]~PC4noOutliers[DiagnosisnoOutliers=="Control"], col=3)
points(PC3noOutliers[DiagnosisnoOutliers=="MDD"]~PC4noOutliers[DiagnosisnoOutliers=="MDD"], col=2)
legend(min(PC3noOutliers), max(PC4noOutliers)+10, c("Control", "MDD"), text.col=c(3, 2), pch=19, col=c(3, 2))
dev.off()

#Outputting the new quantile normalized data (no outliers):
write.table(normNoOutliers,"Quantile normalized filtered data no outliers.txt",sep="\t")

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

#STEP 12:  Importing the detailed subject information into R and putting it into the correct format:

##1. Read in the detailed subject information by highlighting this code and hitting run (Ctrl+R):


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
capture.output(str(HoursFinal))
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

#***END OF CODE TO RUN FOR STEP 13 IF THERE AREN'T OUTLIERS****




##2. If there were outliers, then run this code by highlighting it and hitting ctrl+R:

#***CODE TO RUN FOR STEP 13 IF THERE ARE OUTLIERS****

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

#***END OF CODE TO RUN FOR STEP 13 IF THERE ARE OUTLIERS****


#*******************************************************************

STEP 14:  Looking at the relationships between the independent variables

##1. To examine both visual and statistical relationships between the independent variables (with the outlier subjects removed, if there were any!), just highlight this code and click Ctrl+R:


#***CODE TO RUN FOR STEP 14****

SubjectFactorVariables<-cbind(DiagnosisNoOutliers, ChipNoOutliers, LocationOnChipNoOutliers, CohortNoOutliers, GenderNoOutliers, SuicideNoOutliers, OverdoseNoOutliers)

SubjectContinuousVariables<-cbind(AgeNoOutliers, RNAAgeNoOutliers, TODNoOutliers, BrainpHNoOutliers, RNAconcNoOutliers, RNAintegrityNoOutliers, HoursColdNoOutliers, HoursIceNoOutliers, HoursFinalNoOutliers)


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
#Calculating the variance inflation factor (vif) to determine which subject variables are highly related to other subject variables in the data set. Most important, of course, is whether any of the subject variables strongly correlate with Diagnosis. Note that "Location on Chip" has been removed as a variable because it is partially redundant with gender. 
vif(lm(normNoOutliers[1,]~AgeNoOutliers+BrainpHNoOutliers+ChipNoOutliers+DiagnosisNoOutliers+GenderNoOutliers+RNAconcNoOutliers+RNAintegrityNoOutliers+RNAAgeNoOutliers+TODNoOutliers+HoursColdNoOutliers+HoursIceNoOutliers+HoursFinalNoOutliers))

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
vif(lm(normNoOutliers[1,]~AgeNoOutliers+BrainpHNoOutliers+ChipNoOutliers+DiagnosisNoOutliers+GenderNoOutliers+RNAconcNoOutliers+RNAintegrityNoOutliers+RNAAgeNoOutliers+TODNoOutliers+HoursColdNoOutliers+HoursIceNoOutliers+HoursFinalNoOutliers+PC1noOutliers))
),

capture.output(
vif(lm(normNoOutliers[1,]~AgeNoOutliers+BrainpHNoOutliers+ChipNoOutliers+DiagnosisNoOutliers+GenderNoOutliers+RNAconcNoOutliers+RNAintegrityNoOutliers+RNAAgeNoOutliers+TODNoOutliers+HoursColdNoOutliers+HoursIceNoOutliers+HoursFinalNoOutliers+PC2noOutliers))
),

capture.output(
vif(lm(normNoOutliers[1,]~AgeNoOutliers+BrainpHNoOutliers+ChipNoOutliers+DiagnosisNoOutliers+GenderNoOutliers+RNAconcNoOutliers+RNAintegrityNoOutliers+RNAAgeNoOutliers+TODNoOutliers+HoursColdNoOutliers+HoursIceNoOutliers+HoursFinalNoOutliers+PC3noOutliers))
),

capture.output(
vif(lm(normNoOutliers[1,]~AgeNoOutliers+BrainpHNoOutliers+ChipNoOutliers+DiagnosisNoOutliers+GenderNoOutliers+RNAconcNoOutliers+RNAintegrityNoOutliers+RNAAgeNoOutliers+TODNoOutliers+HoursColdNoOutliers+HoursIceNoOutliers+HoursFinalNoOutliers+PC4noOutliers))
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
FoldChange<-apply(NormFiltered[, Diagnosis=="MDD"], 1, mean)-apply(NormFiltered[, Diagnosis=="Control"], 1, mean)

#Outputting the p-values, fold change, and gene identifiers:
SimpleTtestOutput<-cbind(SimplePvalAdj, FoldChange, GeneIdentifiers_filtered)
row.names(SimpleTtestOutput)<-row.names(normNoOutliers)
colnames(SimpleTtestOutput)<-c("Unadjusted P-value", "BH Adjusted Pvalue", "BY Adjusted Pvalue", "Fold Change", "Symbol", "Search Key", "ILMN_Gene", "Chromosome", "Definition", "Synonyms", "Entrez ID")

write.csv(SimpleTtestOutput, "SimpleTtestOutput.csv")


#***END OF CODE FOR STEP 16****


##2. Take a look at the "16 Histogram of Unadjusted Pvalues Using Simple Ttest". Does it look like there are a disproportionate number of genes with low pvalues?
##3. The output for the analysis using a simple t-test is in "SimpleTtestOutput.csv." When you open it in Excel, resave it as an Excel file. Then highlight the full data set using the corner triangle and sort by BH pvalue.



#*******************************************************************

#STEP 16:  Running a linear model to examine the effects of diagnosis while accounting for potential confounding factors: 

##1. The confounding factors should be chosen based on previous literature or based on which variables were flagged as being related to the PC's or diagnosis.
###For example, in the lateral amygdala, these variables are related to diagnosis:
####A. HoursFinalNoOutliers
####B. SuicidenoOutliers (duh!)
###And these variables are related to the PCs:
####A. DiagnosisNoOutliers(PC1,4)
####B. RNAconcNoOutliers (PC2)
####C. RNAintegrityNoOutliers (PC2)
####D. GenderNoOutliers (PC2)
####E. RNAAgeNoOutliers (PC3)
####F. HoursColdNoOutliers (PC3)
####F. CohortNoOutliers (PC3)
####G. ChipNoOutliers (PC3)
####H. AgeNoOutliers (PC4)	


##Note that not all of these variables are worth including in the model. 
##For example, it is probably worth running suicide as a separate analysis since all suicides are MDD subjects.
##HoursCold and RNAAge both have extreme values that make them prone to producing false positives in statistical tests. They are also highly correlated with Cohort.
##There was also a pattern of association between PC1 to location on chip that might be worth going back and looking at.


##2. Insert those variables into this code, placing diagnosis as the first variable.

LinearModel1<-function(i){lm(normNoOutliers[i,]~DiagnosisNoOutliers+RNAconcNoOutliers+RNAintegrityNoOutliers+GenderNoOutliers+HoursFinalNoOutliers+CohortNoOutliers+ChipNoOutliers+AgeNoOutliers)}


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
#***END OF CODE FOR STEP 17****


##4. The results are found in "LM1testOutput.csv."


#************************************************************************************************

#STEP 18:  Double checking the data for expected effects:

##Gender:

png("18 Gender Check.png")
plot(normNoOutliers[GeneIdentifiers_filtered[,1]=="XIST",]~normNoOutliers[GeneIdentifiers_filtered[,1]=="RPS4Y1",], xlab="Male Gene: RPS4Y1", ylab="Female Gene: XIST", main = "Does subject gender match gendered gene expression?")
points(normNoOutliers[GeneIdentifiers_filtered[,1]=="XIST", GenderNoOutliers=="M"]~normNoOutliers[GeneIdentifiers_filtered[,1]=="RPS4Y1", GenderNoOutliers=="M"], col=2)
points(normNoOutliers[GeneIdentifiers_filtered[,1]=="XIST", GenderNoOutliers=="F"]~normNoOutliers[GeneIdentifiers_filtered[,1]=="RPS4Y1", GenderNoOutliers=="F"], col=3)
legend(min(normNoOutliers[GeneIdentifiers_filtered[,1]=="RPS4Y1",])+0.5, min(normNoOutliers[GeneIdentifiers_filtered[,1]=="XIST",])+0.5, c("M", "F"), text.col=c(2,3), pch=19, col=c(2,3))
dev.off()


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









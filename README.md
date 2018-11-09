# HumanAmygdala_PsychiatricMicroarrays

This is the code release to accompany the submission of the manuscript "Transcriptional profiling of ten subnuclei of the human amygdala in Major Depressive Disorder" by Sharma, Hagenauer et al.

Although it is released under my github account, the code was written by several individuals - Dr. Jun Z. Li, Dr. Peter Blandino, Dr. Elyse Aurbach, Alek Pankonin, and myself. The analysis of the different datasets was performed between the years 2008-2018, and the code used to perform the analysis was written with single usage in mind - I apologize pre-emptively that it is neither well-annotated nor structured to easily allow for rapid reproducibility (the earliest code documentation is actually in .docx and .pptx format which also includes example output instead of a true code document!). There are also additional analyses included in the code documentation that were not included in the manuscript due to space/focus (e.g., additional brain regions, subject variables, and follow-up and control analyses). If you have questions, please ask - I am happy to provide more information!

Brief Description of the files:

1) The Illumina Ref-8 macrodissected tissue analysis code files:

    A) IlluminaRef8_step1_QualityControl_Preprocessing_Freeze3.docx - written by Dr. Jun Z. Li circa 2008 - only some of this output was used for this paper
  
    B) IlluminaRef8_step2_QualityControl_011215.R - written by Dr. Peter Blandino circa 2014-2015 - only some of this output was used for this paper
  
    C) IlluminaRef8_step3_AMY_VariableDiagnostics_Regression_20161013.R - written by Dr. Elyse Aurbach circa 2016

2) The Illumina HT12v3 macrodissected amygdala analysis code files:

    A) IlluminaHT12v3_AMY_Preprocessing.pptx - written by Dr. Megan Hagenauer circa 2012-2013

3) The Illumina LCM-dissected amygdala analysis code files: - written by Dr. Megan Hagenauer circa 2013-2014, some additional analyses 2015-2017. You may note that the code is annotated with instructions - I originally wrote a generalized version of the code for all 10 subnuclei with the intention that someone else would run the analyses, but then ended up just running all of the analysis myself anyway.
  
    A) LCM_AMY_AAA_Preprocessing.R
  
    B) LCM_AMY_AB_Preprocessing.R - Note: starting after line 2000, there is some code for analyzing results from all 10 subnuclei (not just the AB)
  
    C) LCM_AMY_AHA_Preprocessing.R
  
    D) LCM_AMY_Basal_Preprocessing.R
  
    E) LCM_AMY_Central_Preprocessing.R
  
    F) LCM_AMY_CO_Preprocessing.R
  
    G) LCM_AMY_Lateral_Preprocessing.R
  
    H) LCM_AMY_Medial_Preprocessing.R
  
    I) LCM_AMY_PAC_Preprocessing.R
  
    J) LCM_AMY_PL_Preprocessing.R
  
    K) RegionalComparison_July24.R - Note: starting around line 2700 there is code for a failed project to compare relative gene expression levels for genes across subnuclei - despite the fact that the subnulei were in different microarray processing batches! (stupid idea, I know - but it was worth a try, right?) 
  
    L) FiguresForPaper_20180820.R
 
4) Code files for comparing across the major datasets:

    A) ComparingAMY_qRTPCR_vs_Microarray.R - written by Dr. Megan Hagenauer circa 2018
  
    B) ComparingAMYDatasets_forAMYPaper.R - written by Alek Pankonin and Dr. Megan Hagenauer circa 2016-2018
    
    C) IlluminaAndAffyIncBarnesNarayan.R - written by Alek Pankonin circa 2016 - this code is not really analysis code for amygdala datasets in particular, just a general meta-query to draw results from different Pritzker datasets and join them together. It is referenced by ComparingAMYDatasets_forAMYPaper.R.

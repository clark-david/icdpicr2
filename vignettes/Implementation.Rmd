Implementation Details

ICDPIC -- International Classification of Diseases Programs for Injury Categorization was originally developed using ICD Version 9 Clinical Modification (ICD-9-CM) diagnosis codes and Stata statistical software (Statacorp, College Station, Texas).  After the introduction of ICD-10-CM to US hospitals in 2015, an update to accommodate this change was developed using R statistical software (R Project, Vienna, Austria). The context for ICDPIC and ICDPICR, along with a general history of injury severity scoring, has been presented in a previous publication.[1]

Initial development of the ICDPIC Stata programs occurred as part of research projects funded by the National Center for Injury Prevention and Control through the Harvard Injury Control Research Center (CDC R49/CCR 115279) and by the Maine Medical Center (MMC) Research Strategic Plan. The translation of ICDPIC to R was initially supported by funding from the MMC Division of Trauma and Surgical Critical Care and MMC Center for Outcomes Research and Evaluation. The authors are grateful for this support.

ICDPICR was developed using the American College of Surgeons (ACS) Trauma Quality Program (TQP) Participant Use File (PUF) and the Agency for Healthcare Research and Quality (AHRQ) Healthcare Cost and Utilization Project (HCUP) National Inpatient Sample (NIS).  The TQP PUF is the successor to the ACS Trauma Quality Improvement Program (TQIP) Research Data File and the ACS National Trauma Data Bank (NTDB) Research Data Set; NIS was previously called the Nationwide Inpatient Sample.  The original data are not provided as part of ICDPICR, but can be obtained by others following the Data Use Agreements of the sponsoring organizations.  Content reproduced from the TQP PUF remains the full and exclusive copyrighted property of the ACS, which is not responsible for any claims arising from works based on the original data.  Content reproduced from the NIS does not constitute the findings, policies, or recommendations of the U.S. Government, the U.S. Department of Health, or AHRQ.

ICDPICR Version 0.1.0 was provided on github.com, and Version 1.0.0 was added to the Comprehensive R Archive Network (CRAN) in January 2021.  Version 2.0.0 is a further update in response to numerous inquiries and suggestions. The most important changes are as follows:
   ICDPICR initially used TQIP and NIS data from 2017 and 2016.  Version 2.0.0 uses data from 2020.  Version 2.0.0 requires that data be in ICD-10 format.  
   ICDPICR Version 0.1.0 was designed to use only data coded with ICD-10-CM (US Clinical Modification), which limited its value for international users.[2,3]  ICDPICR Version 1.0.0 allowed the user to specify whether data are in ICD-10-CM or a basic ICD-10 format.  Eskesen and colleagues [4] have pointed out that some valid ICD-10 codes were still missing from Version 1.0.0.  ICDPICR Version 2.0.0 attempts to remedy this and has all ICD-10-CM or basic ICD-10 codes (including those previously omitted) in the same table.
   The default “ROCmax” option for calculating Abbreviated Injury Scores [5] in ICDPICR Version 0.1.0 was based upon mortality data in the NTDB, using an ad hoc algorithm to quantify the relative severity of each individual diagnosis code. The ROCmax option in ICDPICR Version 1.0.0 replaced the original ad hoc algorithm with the well-established methodology of ridge regression to estimate the independent effect of each injury diagnosis. Version 1.0.0 also allowed the user to choose either the TQIP PUF or the NIS as the reference database. Version 2.0.0 computes an average of the AIS obtained from TQIP and the AIS obtained from NIS.  In response to user feedback, Version 2.0.0 also modifies the results to make them more consistent (e.g., a right-sided injury has the same AIS as an otherwise identical left-sided injury).  The lookup table of injury codes now also includes valid codes not in either reference database (including basic or truncated ICD-10 codes) and AIS values are assigned to these codes using a weighted average of similar codes. 
   Version 2.0.0 updates the CDC mechanism and intent categories to the most recent versions, and will eventually add the CDC categories of body area and nature of injury (corresponding to the “Barell Matrix” that was developed using ICD-9).[8,9]  
   Version 2.0.0 has also been simplified in several ways so that it will run faster than previous versions. In particular, the user may omit the determination of injury mechanism codes, which has been the most time-consuming step.

Programs used to derive severity scores in Version 2.0.0
   icd10aisA – Reads in raw data from the 2020 TQP PUF or the 2020 NIS.  Identifies cases with a principal injury diagnosis specified by an ICD-10-CM code.  The National Trauma Data Standard used by TQP considers valid ICD-10-CM injury codes to be those in the ranges S00-S99, T07, T14, T20-T28, and T30-T32, so icdaisA recognizes only these codes in the calculation of injury severity. ICDPICR also requires that ICD-10-CM injury codes starting with the letter “S” conclude with the letter “A” (indicating an initial encounter), except for codes indicating a fracture, where codes concluding with the letters “B” or “C” indicate an initial encounter with an open fracture. 
   icd10aisB – Reads in each of the data sets prepared by icdaisA, transforms them into matrices, and performs logistic ridge regression with death as an outcome, using R package glmnet, which is described in detail in the documentation for that package. For each reference dataset (TQP PUF or NIS), the logistic ridge regression results in an independent estimate of effect (log odds ratio) for each diagnosis code. These are tabulated and can be combined with the estimated model intercept to estimate the probability of mortality for individual subjects.  Body regions defined for the original ISS [6] are determined for each ICD-10-CM code.
   icd10aisC – Reads in the tabulated effect estimates for each diagnosis code produced by icdaisB and determines the largest effect estimate in each ISS body region for each subject, which will subsequently be stratified into Abbreviated Injury Scores (AIS) [5] and used to estimate ISS.  For each reference dataset, icdaisC initializes cutpoints categorizing the largest effect estimate for each body region into an AIS score of 1, 2, 3, 4, or 5.  The program then uses an adaptive algorithm that randomly varies the cutpoints to determine the combination of cutpoints for which the c-statistic (area under a Receiver Operator Characteristic curve) for ISS to predict mortality is maximized.  For each diagnosis and reference dataset, the program tabulates the optimal AIS estimates along with the effect estimates and intercepts from ridge regression. 
   icd10aisD (New in Version 2.0.0) – Modifies the AIS scores obtained by icdaisC in the following ways:
1)	Creates a common AIS score for each ICD-10-CM code as a weighted average of the results from TQP PUF and NIS.
2)	Creates a common AIS score for ICD-10-CM codes S06.##2A … S06.##9A, which specify lengths of coma and in some cases whether a subject lived or died after intracranial injury.  Retaining outcome information in the diagnosis code would result in overfitting of the ISS model.
3)	Creates a common AIS score for otherwise identical ICD-10-CM codes that specify whether an injury was right-sided, left-sided, or unspecified with respect to laterality.
4)	Assigns an AIS score to truncated ICD-10-CM codes (3-6 digits long) as a weighted average of the corresponding 7-digit ICD-10-CM codes. This includes most basic ICD-10 codes, but those not produced by truncation are added and similarly assigned an AIS score.
5)	Calculates ISS and NISS.
   icd10aisE (New in Version 2.0.0) – Produces a table of injury mechanisms and intents based upon codes starting with U, V, W, X, or Y in ICD-10 data.  Assigns a cell in the CDC Framework for ICD-10 data (extending it to truncated versions of ICD-10-CM codes).  The TQP PUF and NIS mortality for subjects with an ICD-10-CM diagnosis in each given cell is calculated, and the maximum cell mortality for a given subject is provided as a rough estimate of severity.[10]
   
Programs icd10aisA, icd10aisB, icd10aisC, icd10aisD, and icd10aisE can be found at 
https://github.com/clark-david/icdpicr2/


The following tables are available in R after downloading icdpicr2
•	i10_map_sev: This table is new in Version 2.0.0, replacing i10_map_roc.  It uses the results of icdaisC and idsaisD to produce an approximate AIS score and body region for any ICD-10 code and effect estimates from the ridge regression for any ICD-10-CM code contained either in TQP or NIS.
•	i10_map_mech:  This table is new in Version 2.0.0, replacing i10_ecode.  It uses the results of icdaisE to assign a mechanism category for any ICD-10 external cause of injury code (starting with letters U, V, W, X, or Y), using a CDC proposed table.
•	testdata:  Sample data that can be used to demonstrate the functioning of the program.


Program cat_trauma2
cat_trauma2 – Reads in user data in the specified format.  Calculates AIS, ISS, NISS, mortality predictions from the TQP and NIS models, and injury mechanisms.  Further details about the options available in cat_trauma2 are provided in the help file for this function.


References
1.	Clark DE, Black AW, Skavdahl DH, Hallagan LD. Open-access programs for injury categorization using ICD-9 or ICD-10. Inj Epidemiol 2018; 5:11.
2.	Airaksinen NK, Heinänen MT, Handolin LE. The reliability of the ICD-AIS map in identifying serious road traffic injuries from the Helsinki Trauma Registry. Injury 2019; 50:1545-1551.
3.	Niemann M, Märdian S, Niemann P, Tetteh L, Tsitsilonis S, Braun KF, Stöckle U, Graef F.  Transforming the German ICD-10 (ICD-10-GM) into Injury Severity Score (ISS) – Introducing a new method for automated re-coding.  Plos One 2021: 16(9):e0257183.
4.	Eskesen TO, Sillesen M, Rasmussen LS, Steinmetz J.  Agreement between standard and ICD-10-based Injury Severity Scores.  Clin Epidemiol 2022; 14:201-210.
5.	Committee on Medical Aspects of Automotive Safety, AMA. Rating the severity of tissue damage. I. The abbreviated scale. JAMA 1971; 215:277-280.
6.	Baker SP, O’Neill B, Haddon W Jr., Long WB. The injury severity score: A method for describing patients with multiple injuries and evaluating emergency care. J Trauma 1974; 14:187-196.
7.	Osler T, Baker SP, Long WA. Modification of the injury severity score that both improves accuracy and simplifies scoring. J Trauma 1997; 43:922-925.
8.	Annest JL, Hedegaard H, Chen L, Warner M, Smalls E. Proposed framework for presenting injury data using ICD-10-CM external cause of injury codes. 2014. https://www.cdc.gov/injury/wisqars/pdf/ICD-10-CM_External_Cause_Injury_Codes-a.pdf.
9.	Hedegaard H, Johnson RL, Garnett MF, Thomas KE.  The 2020 International Classification of Diseases, 10th Revision, Clinical Modification injury diagnosis framework for categorizing injuries by body region and nature of injury.  Nat Health Stat Reports 2020; 150:1-26.
10.	Clark DE, Ahmad S.  Estimating injury severity using the Barell matrix.  Inj Prev 2006; 12:111-116.
11.	Copes WS, Champion HR, Sacco WJ, Lawnick MM, Keast SL, Bain LW. The injury severity score revisited. J Trauma 1988; 28:69-77.
12.	Wan V, Reddy S, Thomas A, Issa N, Posluszny J, Schwulst S, Shapiro M, Alam H, Bilimoria KY, Stey AM.  How does Injury Severity Score derived from ICDPIC utilizing ICD-10-CM codes perform compared to Injury Severity Score derived from TQIP?  J Trauma Acute Care Surg 2023; 94:141-147.
13.	Sebastião YV, Metzger GA, Chisolm DJ, Xiang H, Cooper JN. Impact of ICD-9-CM to ICD-10-CM coding transition on trauma hospitalization trends among young adults in 12 states. Inj Epidemiol 2021; 8:4


Contents
•	Programs to derive severity scores
•	Tables used in calculating severity scores
•	Program cat_trauma
•	References

%export

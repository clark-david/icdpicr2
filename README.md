
The purpose of this package is to provide injury researchers with easy
methods for categorization of injuries described using ICD-10 diagnosis
codes and for deriving measures of injury severity.

ICDPIC (International Classification of Diseases Programs for Injury
Categorization) was originally developed for ICD Version 9 Clinical
Modification (ICD-9-CM) diagnosis codes using Stata, providing an easy
way to convert these codes to standard injury severity scores and
categories. After the introduction of ICD-10-CM to US hospitals in 2015,
an update to accommodate this change was developed using R, and was
eventually published on CRAN as package “icdpicr”.

Package “icdpicr2” is a successor to “icdpicr”. The “Implementation
Details” vignette provided with this package describes the history of
injury severity scoring, ICDPIC, and “icdpicr”“; it discusses issues
with the original version of”icdpicr”” and the aims for “icdpicr2”.

Package “icdpicr2”” consists primarily of the function “cat_trauma2”.
This function first reads in user data in a specified format, then
calculates AIS, ISS, NISS, mortality predictions, and injury mechanisms.
It returns the original data file with these fields added. Further
details about the options available in “cat_trauma2” are provided in the
help file for this function.

Other functions in the package derive an “iciss” (International
Classification of Diseases-related Injury Severity Score) and a
“framework” (categorization by injury type and location developed by the
US CDC). Further information is provided in the Implementation Details
and the help files for each function.

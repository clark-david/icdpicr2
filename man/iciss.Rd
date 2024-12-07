% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iciss.R
\name{iciss}
\alias{iciss}
\title{Compute International Classification of Diseases-Based Injury Severity Score (ICISS)}
\usage{
iciss(df, dx_pre, conservative = TRUE, messages = TRUE)
}
\arguments{
\item{df}{A dataframe in wide format containing ICD-10 diagnosis codes with a common column name prefix.
Diagnosis codes should be character strings and may have a decimal or not.}

\item{dx_pre}{Prefix for diagnosis code column names (example: dx1, dx2, etc.)}

\item{conservative}{Should the program exclude DSP based on fewer than 5 observations? Must be TRUE (default) or FALSE.
\itemize{
\item TRUE - DSP based on fewer than 5 observations in the reference database will be excluded.
\item FALSE - DSP based on fewer than 5 observations in the reference database will be included.
}}

\item{messages}{Should the program report completion of each step? Must be TRUE or FALSE (default).
\itemize{
\item TRUE - Messages will report completion of each step (may be helpful for large data sets).
\item FALSE - Messages will not be reported.
}}
}
\value{
A dataframe identical to the dataframe passed to the function with the following additional variables
         added:
         \itemize{
         \item DSP_1-DSP_n: DSP for diagnosis codes 1..n
         \item totaln_1-totaln_n: Number of cases from which DSP was calculated for diagnosis codes 1..n
         \item PS_iciss_prod: ICISS calculated as the product of DSP
         \item PS_iciss_min: ICISS calculated as the minimum of DSP
         }
}
\description{
This function adds Diagnosis-specific Survival Probabilities (DSP) to a dataframe, based on the table
   provided by Gedeborg and colleagues (J Trauma Acute Care Surgery 2014).  Codes longer than four digits
   are treated as if they were four-digit ICD-10 codes as published by the World Health Organization.
   Thus, an ICD-10-CM code like S00552A is considered the same as S005.
For each observation this function will
\enumerate{
   \item assign a severity (DSP) to each valid ICD-10 injury diagnosis code,
   \item calculate one version of ICISS as the product of these DSP, and
   \item calculate another version of ICISS as the minimum of these DSP.
}
}
\details{
Data should be in wide format, as in the example below
}
\examples{
df_in <- read.table(header = TRUE, text = "
    ident   dx1       dx2       dx3
    31416   S32110A   S3251     NA
    31417   S72141A   T07XXXA   D62
")
df_out <- iciss(df_in, "dx", TRUE, FALSE)

}

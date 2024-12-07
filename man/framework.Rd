% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/framework.R
\name{framework}
\alias{framework}
\title{Categorize ICD-10 injury diagnosis codes similar to the "Barell Matrix" used for ICD-9}
\usage{
framework(df, dx_pre, severity = FALSE, messages = FALSE)
}
\arguments{
\item{df}{A dataframe in wide format containing ICD-10 diagnosis codes with a common column name prefix.
Diagnosis codes should be character strings and may have a decimal or not.}

\item{dx_pre}{Prefix for diagnosis code column names (example: dx1, dx2, etc.)}

\item{severity}{Should the program calculate a severity score? Must be TRUE or FALSE (default).
\itemize{
\item TRUE - Program will calculate the estimated survival for each diagnosis and the minimum overall.
\item FALSE - No severity scores will be calculated.
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
         \item cell_1-cell_n: Cell assigned for diagnosis codes 1..n
         \item PsCell_1-PsCell_n: Survival in TQP and NIS for patients with a diagnosis in this cell
         \item bPSmin: The minimum of PsCell_1-PsCell_n for this subject
         }
}
\description{
This function categorizes ICD-10 diagnosis codes according to the table given in
         Hedegaard H, Johnson RL, Garnett MF, Thomas KE. The 2020 International Classification of Diseases,
         10th Revision, Clinical Modification injury diagnosis framework for categorizing injuries by
         body region and nature of injury. Nat Health Stat Reports 2020;150:1-26,
   and (if option selected) predicts mortality for each subject as described in
         Clark DE, Ahmad S. Estimating injury severity using the Barell matrix. Inj Prev 2006;12:111-116.
}
\details{
Data should be in wide format, as in the example below:
}
\examples{
df_in <- read.table(header = TRUE, text = "
    ident   dx1       dx2       dx3
    31416   S32110A   S3251     NA
    31417   S72141A   T07XXXA   D62
")
df_out <- framework(df_in, "dx", TRUE, FALSE)

}

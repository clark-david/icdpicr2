% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cat_trauma2.R
\name{cat_trauma2}
\alias{cat_trauma2}
\title{Categorize trauma data and calculate scores}
\usage{
cat_trauma2(df, dx_pre, messages = TRUE)
}
\arguments{
\item{df}{A dataframe in wide format containing ICD-10 diagnosis codes with a common column name prefix.
Diagnosis codes should be character strings and may have a decimal or not.}

\item{dx_pre}{Prefix for diagnosis code column names (example: dx1, dx2, etc.)}

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
         \item sev_1-sev_n: AIS severity for diagnosis codes 1..n
         \item issbr_1-issbr_n: ISS body region for diagnosis codes 1..n
         \item mxaisbr1-mxaisbr6: maximum AIS severity for each of the 6 ISS body regions
         \item maxais: maximum AIS severity over all ISS body regions
         \item riss: computed injury severity score
         \item niss: computed "new injury severity score"
         \item PmortTQP: TQP model predicted probability of mortality
         \item PmortNIS: NIS model predicted probability of mortality
         \item mechcode_1-mechcode_4: first 4 mechanism codes found in each row of data
         \item mech_1-mech_4: CDC external cause of injury major mechanism for each mechanism code captured
         \item intent_1-intent_4: intent for each mechanism code captured
         }
}
\description{
This function adds Abbreviated Injury Scores (AIS), Injury Severity Scores (ISS), and other descriptors of injury to a dataframe.
For each observation this function will
\enumerate{
   \item assign a severity (AIS) and ISS body region values to each valid ICD-10 injury diagnosis code,
   \item add variables for maximum severity of each body region,
   \item calculate ISS, "New ISS", maximum AIS, and regression-based mortality predictions,
   \item select the first 4 mechanism (external cause) codes and categorize mechanism and intent following CDC guidelines
}
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
df_out <- cat_trauma2(df_in, "dx", TRUE)

}

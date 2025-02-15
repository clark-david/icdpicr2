% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{i10_map_sev}
\alias{i10_map_sev}
\title{Table i10_map_sev}
\format{
A data frame with 20,551 rows and 8 variables:
\describe{
  \item{dx}{A valid ICD-10 diagnosis code.}
  \item{severity}{The associated Abbreviated Injury Score for this diagnosis.}
  \item{issbr}{The associated ISS body region for this diagnosis.}
  \item{TQIPeffect}{The coefficient for this diagnosis in the TQP regression model.}
  \item{TQIPint}{The intercept in the TQP regression model.}
  \item{NISeffect}{The coefficient for this diagnosis in the NIS regression model.}
  \item{NISint}{The intercept in the TQP regression model.}
  \item{version}{The most recent date when this dataset was revised.  Encoded vyymmdd.}
  ...
}
}
\usage{
i10_map_sev
}
\description{
A dataset containing ICD-10 diagnosis codes and their properties
for determination of ISS and other injury severity scores.
}
\keyword{datasets}

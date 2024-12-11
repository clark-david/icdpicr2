#' ICD10CM injury codes
#'
#' A dataset containing ICD-10 diagnosis codes and their properties
#' for determination of ISS and other injury severity scores.
#'
#' @format A data frame with 20,551 rows and 8 variables:
#' \describe{
#'   \item{dx}{A valid ICD-10 diagnosis code.}
#'   \item{severity}{The associated Abbreviated Injury Score for this diagnosis.}
#'   \item{issbr}{The associated ISS body region for this diagnosis.}
#'   \item{TQIPeffect}{The coefficient for this diagnosis in the TQP regression model.}
#'   \item{TQIPint}{The intercept in the TQP regression model.}
#'   \item{NISeffect}{The coefficient for this diagnosis in the NIS regression model.}
#'   \item{NISint}{The intercept in the TQP regression model.}
#'   \item{version}{The most recent date when this dataset was revised.  Encoded vyymmdd.}
#'   ...
#' }
"i10_map_sev"


#' ICD10CM injury codes
#'
#' A dataset containing ICD-10 diagnosis codes and their properties
#' for determination of injury mechanism.
#'
#' @format A data frame with 8,117 rows and 5 variables:
#' \describe{
#'   \item{dx}{A valid ICD-10 diagnosis code.}
#'   \item{mechmaj}{The CDC major mechanism category.}
#'   \item{intent}{The CDC intent category.}
#'   \item{mechmin}{The CDC minor mechanism category, if any.}
#'   \item{version}{The most recent date when this dataset was revised.  Encoded vyymmdd.}
#'   ...
#' }
"i10_map_mech"


#' ICD10CM injury codes
#'
#' A dataset containing ICD-10 diagnosis codes and their properties
#' for calculation of ICD Injury Severity Score (ICISS).
#'
#' @format A data frame with 20,551 rows and 5 variables:
#' \describe{
#'   \item{dx}{A valid ICD-10 diagnosis code.}
#'   \item{totaln}{The number of cases from which the observed survival was calculated.}
#'   \item{dsp_noncons}{The observed survival for subjects with this diagnosis.}
#'   \item{dsp_cons}{Same as above, but NA if totaln < 5.}
#'   \item{version}{The most recent date when this dataset was revised.  Encoded vyymmdd.}
#'   ...
#' }
"i10_map_iciss"


#' ICD10CM injury codes
#'
#' A dataset containing ICD-10 diagnosis codes and their properties
#' for assignment to the CDC Framework.
#'
#' @format A data frame with 20,499 rows and 5 variables:
#' \describe{
#'   \item{dx}{A valid ICD-10 diagnosis code.}
#'   \item{cell}{The corresponding cell in the CDC Framework.}
#'   \item{PmCell}{The observed mortality for subjects with a diagnosis in this cell.}
#'   \item{PmCell}{The observed survival for subjects with a diagnosis in this cell.}
#'   \item{version}{The most recent date when this dataset was revised.  Encoded vyymmdd.}
#'   ...
#' }
"i10_map_frame"


#' ICD10CM injury codes
#'
#' A dataset containing a sample of trauma registry data for use in examples and tests.
#'
#' @format A data frame with 20,000 rows and 13 variables:
#' \describe{
#'   \item{temp_id}{A sequential number to identify individual subjects.}
#'   \item{died}{A binary indicator variable for death. 1 = died. 0 = survived.}
#'   \item{dx1}{1st ICD-10-CM injury code recorded on an encounter.}
#'   \item{dx2}{2nd ICD-10-CM injury code recorded on an encounter.}
#'   \item{dx3}{3rd ICD-10-CM injury code recorded on an encounter.}
#'   \item{dx4}{4th ICD-10-CM injury code recorded on an encounter.}
#'   \item{dx5}{5th ICD-10-CM injury code recorded on an encounter.}
#'   \item{dx6}{6th ICD-10-CM injury code recorded on an encounter.}
#'   \item{dx7}{7th ICD-10-CM injury code recorded on an encounter.}
#'   \item{dx8}{8th ICD-10-CM injury code recorded on an encounter.}
#'   \item{dx9}{9th ICD-10-CM injury code recorded on an encounter.}
#'   \item{dx10}{10th ICD-10-CM injury code recorded on an encounter.}
#'   \item{dx11}{11th ICD-10-CM injury code recorded on an encounter.}
#'   \item{dx12}{12th ICD-10-CM injury code recorded on an encounter.}
#'   \item{dx13}{13th ICD-10-CM injury code recorded on an encounter.}
#'   \item{dx14}{14th ICD-10-CM injury code recorded on an encounter.}
#'   \item{dx15}{15th ICD-10-CM injury code recorded on an encounter.}
#'   \item{dx16}{16th ICD-10-CM injury code recorded on an encounter.}
#'   \item{dx17}{17th ICD-10-CM injury code recorded on an encounter.}
#'   \item{dx18}{18th ICD-10-CM injury code recorded on an encounter.}
#'   \item{dx19}{19th ICD-10-CM injury code recorded on an encounter.}
#'   \item{dx20}{20th ICD-10-CM injury code recorded on an encounter.}
#'   ...
#' }
"testdata"

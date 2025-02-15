#' Table i10_map_sev
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


#' Table i10_map_mech
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


#' Table i10_map_iciss
#'
#' A dataset containing ICD-10 diagnosis codes and their properties
#' for calculation of several versions of the ICD Injury Severity Score (ICISS).
#'
#' @format A data frame with 20,551 rows and 12 variables:
#' \describe{
#'   \item{dx}{A valid ICD-10 diagnosis code.}
#'   \item{digits1234}{The 4-digit version of dx.}
#'   \item{tot_int}{The number of cases in international data from which the observed survival was calculated.}
#'   \item{dsp_int}{The observed survival in international data for subjects with this diagnosis.}
#'   \item{dsp_int_c}{Same as above, but NA if tot_int < 5.}
#'   \item{tot_TQP}{The number of cases in TQP data from which the observed survival was calculated.}
#'   \item{dsp_TQP}{The observed survival in TQP data for subjects with this diagnosis.}
#'   \item{dsp_TQP_c}{Same as above, but NA if tot_TQP < 5.}
#'   \item{tot_NIS}{The number of cases in NIS data from which the observed survival was calculated.}
#'   \item{dsp_NIS}{The observed survival in NIS data for subjects with this diagnosis.}
#'   \item{dsp_NIS_c}{Same as above, but NA if tot_NIS < 5.}
#'   \item{Version}{The most recent date when this dataset was revised.  Encoded vyymmdd.}
#'   ...
#' }
"i10_map_iciss"


#' Table i10_map_frame
#'
#' A dataset containing ICD-10 diagnosis codes and their properties
#' for assignment to the CDC Framework.
#'
#' @format A data frame with 20,499 rows and 5 variables:
#' \describe{
#'   \item{dx}{A valid ICD-10 diagnosis code.}
#'   \item{cell}{The corresponding cell in the CDC Framework.}
#'   \item{PmCell}{The observed mortality for subjects with a diagnosis in this cell.}
#'   \item{PsCell}{The observed survival for subjects with a diagnosis in this cell.}
#'   \item{version}{The most recent date when this dataset was revised.  Encoded vyymmdd.}
#'   ...
#' }
"i10_map_frame"


#' Test Data
#'
#' A dataset containing a sample of trauma registry data for use in examples and tests.
#'
#' @format A data frame with 20,000 rows and 13 variables:
#' \describe{
#'   \item{temp_id}{A sequential number to identify individual subjects.}
#'   \item{died}{A binary indicator variable for death. 1 = died. 0 = survived.}
#'   \item{I10_DX1}{1st ICD-10-CM injury code recorded on an encounter.}
#'   \item{I10_DX2}{2nd ICD-10-CM injury code recorded on an encounter.}
#'   \item{I10_DX3}{3rd ICD-10-CM injury code recorded on an encounter.}
#'   \item{I10_DX4}{4th ICD-10-CM injury code recorded on an encounter.}
#'   \item{I10_DX5}{5th ICD-10-CM injury code recorded on an encounter.}
#'   \item{I10_DX6}{6th ICD-10-CM injury code recorded on an encounter.}
#'   \item{I10_DX7}{7th ICD-10-CM injury code recorded on an encounter.}
#'   \item{I10_DX8}{8th ICD-10-CM injury code recorded on an encounter.}
#'   \item{I10_DX9}{9th ICD-10-CM injury code recorded on an encounter.}
#'   \item{I10_DX10}{10th ICD-10-CM injury code recorded on an encounter.}
#'   \item{Version}{The most recent date when this dataset was revised. Encoded vyymmdd.}
#'   ...
#' }
"testdata"

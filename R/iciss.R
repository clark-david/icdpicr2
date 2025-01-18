
#' Compute International Classification of Diseases-Based Injury Severity Score (ICISS)
#'
#' This function adds Diagnosis-specific Survival Probabilities (DSP) to a dataframe, based on the table provided by
#'    Gedeborg and colleagues (J Trauma Acute Care Surgery 2014) and also TQP and NIS.  ICD-10-CM codes longer than
#'    four digits are treated as if they were four-digit ICD-10 codes as published by the World Health Organization.
#'    Thus, an ICD-10-CM code like S00552A is considered the same as S005.
#' For each observation this function will
#' \enumerate{
#'    \item assign a severity (DSP) to each valid ICD-10 injury diagnosis code,
#'    \item calculate one version of ICISS as the product of these DSP, and
#'    \item calculate another version of ICISS as the minimum of these DSP.
#'    \item It repeats the above using international data, TQP, and NIS as reference data.
#'}
#'
#'
#' @param df A dataframe in wide format containing ICD-10 diagnosis codes with a common column name prefix.
#'           Diagnosis codes should be character strings and may have a decimal or not.
#'
#' @param dx_pre Prefix for diagnosis code column names (example: dx1, dx2, etc.)
#'
#' @param conservative Should the program exclude DSP based on fewer than 5 observations? Must be TRUE (default) or FALSE.
#'          \itemize{
#'          \item TRUE - DSP based on fewer than 5 observations in the reference database will be excluded.
#'          \item FALSE - DSP based on fewer than 5 observations in the reference database will be included.
#'          }
#'
#' @param messages Should the program report completion of each step? Must be TRUE or FALSE (default).
#'          \itemize{
#'          \item TRUE - Messages will report completion of each step (may be helpful for large data sets).
#'          \item FALSE - Messages will not be reported.
#'          }
#'
#' @return A dataframe identical to the dataframe passed to the function with the following additional variables
#'          added:
#'          \itemize{
#'          \item dsp_int1-dsp_intn: DSP for diagnosis codes 1..n, from international data
#'          \item dsp_TQP1-dsp_TQPn: DSP for diagnosis codes 1..n, from TQP data
#'          \item dsp_TQP1-dsp_NISn: DSP for diagnosis codes 1..n, from NIS data
#'          \item PS_int_prod: ICISS calculated as the product of dsp_int1-dsp_intn
#'          \item PS_int_min: ICISS calculated as the minimum of dsp_int1-dsp_intn
#'          \item PS_TQP_prod: ICISS calculated as the product of dsp_TQP1-dsp_TQPn
#'          \item PS_TQP_min: ICISS calculated as the minimum of dsp_TQP1-dsp_TQPn
#'          \item PS_NIS_prod: ICISS calculated as the product of dsp_NIS1-dsp_NISn
#'          \item PS_NIS_min: ICISS calculated as the minimum of dsp_NIS1-dsp_NISn
#'          }
#'
#' @details  Data should be in wide format, as in the example below
#'
#' @examples
#' df_in <- read.table(header = TRUE, text = "
#'     ident   dx1       dx2       dx3
#'     31416   S32110A   S3251     NA
#'     31417   S72141A   T07XXXA   D62
#' ")
#' df_out <- iciss(df_in, "dx", TRUE, FALSE)
#'
#' @importFrom stats na.omit
#' @export


iciss <- function(df, dx_pre, conservative=TRUE, messages=TRUE) {

  #Version 250111

  starttime=Sys.time()

  # Verify input
  if(!is.data.frame(df)) stop("First argument must be a dataframe")
  if(NROW(df) == 0) stop("Data file contains no observations. It must contain at least one row")
  if(!is.character(dx_pre)) stop("Second argument must be a character string")
  # Ensure dx_pre is a valid variable name
  if(make.names(dx_pre) != dx_pre) stop("Second argument must be a valid variable name in R")
  # Check if user entered a correct prefix for the diagnosis code variables in the input file
  # Determine how many diagnosis code variables there are in the data
  regex_dx <- paste0("^", dx_pre, "([0-9]+)$")
  dx_colnames <- grep(regex_dx, names(df), value = TRUE)
  # Replace full column name with first capture group and convert to number
  dx_nums <- as.numeric(sub(regex_dx, "\\1", dx_colnames))
  num_dx <- length(dx_nums)
  if(num_dx == 0) stop("No variables with prefix found in data")

  # Make sure df is not a tibble and if it is convert back to regular dataframe
  df <- data.frame(df)

  itab <- i10_map_iciss

  if(conservative==TRUE) {
    itab <- itab[ , c("dx","dsp_int_c","dsp_TQP_c","dsp_NIS_c")]
    itab <- dplyr::rename(itab,dsp_int=dsp_int_c,dsp_TQP=dsp_TQP_c,dsp_NIS=dsp_NIS_c)
  }
  if(conservative==FALSE) {
    itab <- itab[ , c("dx","dsp_int","dsp_TQP","dsp_NIS")]
  }


  #--------------------------------------------------------------------#
  #  Merge diagnosis code variables with reference table to obtain DSP #
  #     for each diagnosis code and add them to the data               #
  #--------------------------------------------------------------------#
  for(i in dx_nums){

    if( (messages==TRUE) && (i%%5==0 || i==num_dx) ){
      message("Determining DSPs for Diagnosis ", i, " of ", num_dx)
    }

    # Create column name
    dx_name <- paste0(dx_pre, i)

    # Pull just the diagnosis code column of interest
    df_ss <- df[ , dx_name, drop = FALSE]

    # Add row variable for sorting back to original order
    df_ss$n <- 1:NROW(df_ss)

    # Strip out decimal in all codes
    df_ss[ , dx_name] <- sub("\\.", "", df_ss[ , dx_name])

    # Get rid of codes that do not start with a valid character
    i10_valid <- c("S","T")
    df_ss[ , dx_name] <- ifelse(substr(df_ss[,dx_name], 1, 1) %in% c(i10_valid), df_ss[,dx_name], NA)

    # Merge with lookup table for severity
    temp <- merge(df_ss, itab, by.x = dx_name, by.y = "dx", all.x = TRUE, all.y = FALSE, sort = FALSE)

    # Reorder rows after merge
    temp <- temp[order(temp$n), ]

    # Reorder columns and drop dx and n
    temp <- temp[ , c("dsp_int","dsp_TQP","dsp_NIS")]
    # Rename columns
    names(temp) <- paste0(c("dsp_int","dsp_TQP","dsp_NIS"), i)

    # Add temp columns to dataframe
    df <- .insert_columns(df, dx_name, temp)

  }   #END FOR LOOP (i in dxnums)


  #-----------------------------------------------------------#
  # Add survival predictions using product or minimum of DSP  #
  #-----------------------------------------------------------#

  if(messages==TRUE){
    mindiff=round(as.double(difftime(Sys.time(),starttime,units="secs"))/60)
    message("Time elapsed ", mindiff, " minutes")
    message("Calculating survival predictions based on international data")
  }

  # Duplicate table of diagnoses and convert to long form
  df <- dplyr::mutate(df,RowID=dplyr::row_number())
  df_calc <- df
  df_calc <- dplyr::select(df_calc,RowID,dplyr::starts_with("dsp_int"))
  df_calc <- tidyr::pivot_longer(df_calc,cols=tidyr::starts_with("dsp_int"),names_to="ColName")
  df_calc <- dplyr::group_by(df_calc,RowID)
  df_calc <- dplyr::mutate(df_calc,ColID1=dplyr::row_number())
  df_calc <- dplyr::ungroup(df_calc)
  df_calc <- dplyr::rename(df_calc,dsp_int=value)
  df_calc <- dplyr::select(df_calc,-ColName)

  # Calculate minimum and product for each individual
  # Without modification, prod() returns 1 for empty set,
  #    and min() returns Inf or -Inf for empty set, with warnings
  # The following code returns NA if no diagnosis has a DSP from the lookup table
  df_calc2 <- dplyr::group_by(df_calc,RowID)
  df_calc2 <- dplyr::mutate(df_calc2,all_na=all(is.na(dsp_int)))
  df_calc2 <- dplyr::mutate(df_calc2,PS_int_prod=prod(dsp_int,na.rm=TRUE))
  df_calc2 <- dplyr::mutate(df_calc2,PS_int_prod=dplyr::if_else(all_na==TRUE,NA,PS_int_prod))
  df_calc2 <- suppressWarnings(dplyr::mutate(df_calc2,PS_int_min=min(dsp_int,na.rm=TRUE)))
  df_calc2 <- dplyr::mutate(df_calc2,PS_int_min=dplyr::if_else(all_na==TRUE,NA,PS_int_min))
  df_calc2 <- dplyr::ungroup(df_calc2)

  # Keep one set of results for each individual and add to original dataframe
  df_results <- dplyr::filter(df_calc2,ColID1==1)
  df_results <- dplyr::select(df_results,RowID,PS_int_prod,PS_int_min)
  df_results <- dplyr::rename(df_results,RowID2=RowID)
  df_results <- dplyr::arrange(df_results,RowID2)

  df <- dplyr::arrange(df,RowID)
  df <- dplyr::bind_cols(df,df_results)
  df <- dplyr::select(df,-dplyr::starts_with("RowID"))


  if(messages==TRUE){
    mindiff=round(as.double(difftime(Sys.time(),starttime,units="secs"))/60)
    message("Time elapsed ", mindiff, " minutes")
    message("Calculating survival predictions based on TQP data")
  }

  # Duplicate table of diagnoses and convert to long form
  df <- dplyr::mutate(df,RowID=dplyr::row_number())
  df_calc <- df
  df_calc <- dplyr::select(df_calc,RowID,dplyr::starts_with("dsp_TQP",ignore.case=FALSE))
  df_calc <- tidyr::pivot_longer(df_calc,cols=tidyr::starts_with("dsp_TQP",ignore.case=FALSE),names_to="ColName")
  df_calc <- dplyr::group_by(df_calc,RowID)
  df_calc <- dplyr::mutate(df_calc,ColID1=dplyr::row_number())
  df_calc <- dplyr::ungroup(df_calc)
  df_calc <- dplyr::rename(df_calc,dsp_TQP=value)
  df_calc <- dplyr::select(df_calc,-ColName)

  # Calculate minimum and product for each individual
  # Without modification, prod() returns 1 for empty set,
  #    and min() returns Inf or -Inf for empty set, with warnings
  # The following code returns NA if no diagnosis has a DSP from the lookup table
  df_calc2 <- dplyr::group_by(df_calc,RowID)
  df_calc2 <- dplyr::mutate(df_calc2,all_na=all(is.na(dsp_TQP)))
  df_calc2 <- dplyr::mutate(df_calc2,PS_TQP_prod=prod(dsp_TQP,na.rm=TRUE))
  df_calc2 <- dplyr::mutate(df_calc2,PS_TQP_prod=dplyr::if_else(all_na==TRUE,NA,PS_TQP_prod))
  df_calc2 <- suppressWarnings(dplyr::mutate(df_calc2,PS_TQP_min=min(dsp_TQP,na.rm=TRUE)))
  df_calc2 <- dplyr::mutate(df_calc2,PS_TQP_min=dplyr::if_else(all_na==TRUE,NA,PS_TQP_min))
  df_calc2 <- dplyr::ungroup(df_calc2)

  # Keep one set of results for each individual and add to original dataframe
  df_results <- dplyr::filter(df_calc2,ColID1==1)
  df_results <- dplyr::select(df_results,RowID,PS_TQP_prod,PS_TQP_min)
  df_results <- dplyr::rename(df_results,RowID2=RowID)
  df_results <- dplyr::arrange(df_results,RowID2)

  df <- dplyr::arrange(df,RowID)
  df <- dplyr::bind_cols(df,df_results)
  df <- dplyr::select(df,-dplyr::starts_with("RowID"))


  if(messages==TRUE){
    mindiff=round(as.double(difftime(Sys.time(),starttime,units="secs"))/60)
    message("Time elapsed ", mindiff, " minutes")
    message("Calculating survival predictions based on NIS data")
  }

  # Duplicate table of diagnoses and convert to long form
  df <- dplyr::mutate(df,RowID=dplyr::row_number())
  df_calc <- df
  df_calc <- dplyr::select(df_calc,RowID,dplyr::starts_with("dsp_NIS",ignore.case=FALSE))
  df_calc <- tidyr::pivot_longer(df_calc,cols=tidyr::starts_with("dsp_NIS",ignore.case=FALSE),names_to="ColName")
  df_calc <- dplyr::group_by(df_calc,RowID)
  df_calc <- dplyr::mutate(df_calc,ColID1=dplyr::row_number())
  df_calc <- dplyr::ungroup(df_calc)
  df_calc <- dplyr::rename(df_calc,dsp_NIS=value)
  df_calc <- dplyr::select(df_calc,-ColName)

  # Calculate minimum and product for each individual
  # Without modification, prod() returns 1 for empty set,
  #    and min() returns Inf or -Inf for empty set, with warnings
  # The following code returns NA if no diagnosis has a DSP from the lookup table
  df_calc2 <- dplyr::group_by(df_calc,RowID)
  df_calc2 <- dplyr::mutate(df_calc2,all_na=all(is.na(dsp_NIS)))
  df_calc2 <- dplyr::mutate(df_calc2,PS_NIS_prod=prod(dsp_NIS,na.rm=TRUE))
  df_calc2 <- dplyr::mutate(df_calc2,PS_NIS_prod=dplyr::if_else(all_na==TRUE,NA,PS_NIS_prod))
  df_calc2 <- suppressWarnings(dplyr::mutate(df_calc2,PS_NIS_min=min(dsp_NIS,na.rm=TRUE)))
  df_calc2 <- dplyr::mutate(df_calc2,PS_NIS_min=dplyr::if_else(all_na==TRUE,NA,PS_NIS_min))
  df_calc2 <- dplyr::ungroup(df_calc2)

  # Keep one set of results for each individual and add to original dataframe
  df_results <- dplyr::filter(df_calc2,ColID1==1)
  df_results <- dplyr::select(df_results,RowID,PS_NIS_prod,PS_NIS_min)
  df_results <- dplyr::rename(df_results,RowID2=RowID)
  df_results <- dplyr::arrange(df_results,RowID2)

  df <- dplyr::arrange(df,RowID)
  df <- dplyr::bind_cols(df,df_results)
  df <- dplyr::select(df,-dplyr::starts_with("RowID"))


  if(messages==TRUE){
    mindiff=round(as.double(difftime(Sys.time(),starttime,units="secs"))/60)
    message("Time elapsed ", mindiff, " minutes")
  }

  message("=============================================")
  message("REMINDER")
  message("ICDPICR Version 2.0.5 IS BEING TESTED")
  message("Major bugs and flaws may still exist")
  message("Please report issues to david.clark@tufts.edu")
  message("or at github/clark-david/icdpicr2/issues")
  message("==============================================")

  # Return dataframe
  df

} #END iciss


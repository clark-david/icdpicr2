
#' Compute International Classification of Diseases-Based Injury Severity Score (ICISS)
#'
#' This function adds Diagnosis-specific Survival Probabilities (DSP) to a dataframe, based on the table
#'    provided by Gedeborg and colleagues (J Trauma Acute Care Surgery 2014).  Codes longer than four digits
#'    are treated as if they were four-digit ICD-10 codes as published by the World Health Organization.
#'    Thus, an ICD-10-CM code like S00552A is considered the same as S005.
#' For each observation this function will
#' \enumerate{
#'    \item assign a severity (DSP) to each valid ICD-10 injury diagnosis code,
#'    \item calculate one version of ICISS as the product of these DSP, and
#'    \item calculate another version of ICISS as the minimum of these DSP.
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
#'          \item DSP_1-DSP_n: DSP for diagnosis codes 1..n
#'          \item totaln_1-totaln_n: Number of cases from which DSP was calculated for diagnosis codes 1..n
#'          \item PS_iciss_prod: ICISS calculated as the product of DSP
#'          \item PS_iciss_min: ICISS calculated as the minimum of DSP
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

  #Version 241216

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
    itab <- itab[ , c("dx","totaln","dsp_cons")]
    itab <- dplyr::rename(itab,dsp=dsp_cons)
  }
  if(conservative==FALSE) {
    itab <- itab[ , c("dx","totaln","dsp_noncons")]
    itab <- dplyr::rename(itab,dsp=dsp_noncons)
  }


  #--------------------------------------------------------------------#
  #  Merge diagnosis code variables with reference table to obtain DSP #
  #     for each diagnosis code and add them to the data               #
  #--------------------------------------------------------------------#
  for(i in dx_nums){

    if( (messages==TRUE) && (i%%5==0 || i==num_dx) ){
      message("Determining DSP for Diagnosis ", i, " of ", num_dx)
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
    temp <- temp[ , c("totaln","dsp")]

    # Rename columns
    names(temp) <- paste0(c("totaln","dsp"), i)

    # Add temp columns to dataframe
    df <- .insert_columns(df, dx_name, temp)

  }   #END FOR LOOP (i in dxnums)


  #-----------------------------------------------------------#
  # Add survival predictions using product or minimum of DSP  #
  #-----------------------------------------------------------#

  if(messages==TRUE){
    mindiff=round(as.double(difftime(Sys.time(),starttime,units="secs"))/60)
    message("Time elapsed ", mindiff, " minutes")
    message("Calculating survival predictions")
  }

  # Duplicate table of diagnoses and convert to long form
  df <- dplyr::mutate(df,RowID=dplyr::row_number())
  df_calc <- df
  df_calc <- dplyr::select(df_calc,RowID,dplyr::starts_with("dsp"))
  df_calc <- tidyr::pivot_longer(df_calc,cols=tidyr::starts_with("dsp"),names_to="ColName")
  df_calc <- dplyr::group_by(df_calc,RowID)
  df_calc <- dplyr::mutate(df_calc,ColID1=dplyr::row_number())
  df_calc <- dplyr::ungroup(df_calc)
  df_calc <- dplyr::rename(df_calc,dsp=value)
  df_calc <- dplyr::select(df_calc,-ColName)

  # Calculate minimum and product for each individual
  # Without modification, prod() returns 1 for empty set,
  #    and min() returns Inf or -Inf for empty set, with warnings
  # The following code returns NA if no diagnosis has a DSP from the lookup table
  df_calc2 <- dplyr::group_by(df_calc,RowID)
  df_calc2 <- dplyr::mutate(df_calc2,all_na=all(is.na(dsp)))
  df_calc2 <- dplyr::mutate(df_calc2,PS_iciss_prod=prod(dsp,na.rm=TRUE))
  df_calc2 <- dplyr::mutate(df_calc2,PS_iciss_prod=dplyr::if_else(all_na==TRUE,NA,PS_iciss_prod))
  df_calc2 <- suppressWarnings(dplyr::mutate(df_calc2,PS_iciss_min=min(dsp,na.rm=TRUE)))
  df_calc2 <- dplyr::mutate(df_calc2,PS_iciss_min=dplyr::if_else(all_na==TRUE,NA,PS_iciss_min))
  df_calc2 <- dplyr::ungroup(df_calc2)

  # Keep one set of results for each individual and add to original dataframe
  df_results <- dplyr::filter(df_calc2,ColID1==1)
  df_results <- dplyr::select(df_results,RowID,PS_iciss_prod,PS_iciss_min)
  df_results <- dplyr::rename(df_results,RowID2=RowID)
  df_results <- dplyr::arrange(df_results,RowID2)

  df <- dplyr::select(df,-dplyr::starts_with("totaln"))
  df <- dplyr::arrange(df,RowID)
  df <- dplyr::bind_cols(df,df_results)
  df <- dplyr::select(df,-dplyr::starts_with("RowID"))

  if(messages==TRUE){
    mindiff=round(as.double(difftime(Sys.time(),starttime,units="secs"))/60)
    message("Time elapsed ", mindiff, " minutes")
  }

  message("=============================================")
  message("REMINDER")
  message("ICDPICR Version 2.0.4 IS BEING TESTED")
  message("Major bugs and flaws may still exist")
  message("Please report issues to david.clark@tufts.edu")
  message("or at github/clark-david/icdpicr2/issues")
  message("==============================================")

  # Return dataframe
  df

} #END iciss


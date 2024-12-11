

#' Categorize ICD-10 injury diagnosis codes similar to the "Barell Matrix" used for ICD-9
#'
#'
#' This function categorizes ICD-10 diagnosis codes according to the table given in
#'          Hedegaard H, Johnson RL, Garnett MF, Thomas KE. The 2020 International Classification of Diseases,
#'          10th Revision, Clinical Modification injury diagnosis framework for categorizing injuries by
#'          body region and nature of injury. Nat Health Stat Reports 2020;150:1-26,
#'    and (if option selected) predicts mortality for each subject as described in
#'          Clark DE, Ahmad S. Estimating injury severity using the Barell matrix. Inj Prev 2006;12:111-116.
#'
#' @param df A dataframe in wide format containing ICD-10 diagnosis codes with a common column name prefix.
#'           Diagnosis codes should be character strings and may have a decimal or not.
#'
#' @param dx_pre Prefix for diagnosis code column names (example: dx1, dx2, etc.)
#'
#' @param severity Should the program calculate a severity score? Must be TRUE or FALSE (default).
#'          \itemize{
#'          \item TRUE - Program will calculate the estimated survival for each diagnosis and the minimum overall.
#'          \item FALSE - No severity scores will be calculated.
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
#'          \item cell_1-cell_n: Cell assigned for diagnosis codes 1..n
#'          \item PsCell_1-PsCell_n: Survival in TQP and NIS for patients with a diagnosis in this cell
#'          \item bPSmin: The minimum of PsCell_1-PsCell_n for this subject
#'          }
#'
#' @details  Data should be in wide format, as in the example below:
#'
#'
#' @examples
#' df_in <- read.table(header = TRUE, text = "
#'     ident   dx1       dx2       dx3
#'     31416   S32110A   S3251     NA
#'     31417   S72141A   T07XXXA   D62
#' ")
#' df_out <- framework(df_in, "dx", TRUE, FALSE)
#'
#' @importFrom stats na.omit
#' @export


framework <- function(df, dx_pre, severity=FALSE, messages=FALSE) {

  #Version 241211

  require(dplyr)
  require(readr)
  require(tidyr)
  require(stringr)
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

  ftab <- i10_map_frame
  ftab <- ftab[ , c("dx","cell","PsCell")]

  #---------------------------------------------------------------------------------------#
  #  Merge diagnosis code variables with reference table to obtain cell name              #
  #  and cell survival probability for each diagnosis code and add them to the data       #
  #---------------------------------------------------------------------------------------#
  for(i in dx_nums){

    if( (messages==TRUE) && (i%%5==0 || i==num_dx) ){
      message("Determining cell for Diagnosis ", i, " of ", num_dx)
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
    temp <- merge(df_ss, ftab, by.x = dx_name, by.y = "dx", all.x = TRUE, all.y = FALSE, sort = FALSE)

    # Reorder rows after merge
    temp <- temp[order(temp$n), ]

    # Reorder columns and drop dx and n
    temp <- temp[ , c("cell","PsCell")]

    # Rename columns
    names(temp) <- paste0(c("cell_","PsCell_"), i)

    # Add temp columns to dataframe
    df <- .insert_columns(df, dx_name, temp)

  }   #END FOR LOOP (i in dxnums)


  if (severity==TRUE) {
    #----------------------------------------------------------#
    # Add mortality prediction for ICD-10-cm codes from cells  #
    #----------------------------------------------------------#


    if(messages==TRUE){
      mindiff=round(as.double(difftime(Sys.time(),starttime,units="secs"))/60)
      message("Time elapsed ", mindiff, " minutes")
      message("Calculating survival prediction")
    }

    # Duplicate table of diagnoses and convert to long form
    df <- dplyr::mutate(df,RowID=row_number())
    df_calc <- df
    df_calc <- dplyr::select(df_calc,RowID,starts_with("PsCell"))
    df_calc <- tidyr::pivot_longer(df_calc,cols=starts_with("PsCell"),names_to="ColName")
    df_calc <- dplyr::group_by(df_calc,RowID)
    df_calc <- dplyr::mutate(df_calc,ColID1=row_number())
    df_calc <- dplyr::ungroup(df_calc)
    df_calc <- dplyr::rename(df_calc,PsCell=value)
    df_calc <- dplyr::select(df_calc,-ColName)

    # Calculate minimum DSP for each individual
    # Without modification, min() returns Inf or -Inf for empty set, with warnings
    # The following code returns NA if no diagnosis has a DSP in the lookup table
    df_calc2 <- dplyr::group_by(df_calc,RowID)
    df_calc2 <- dplyr::mutate(df_calc2,all_na=all(is.na(PsCell)))
    df_calc2 <- suppressWarnings(mutate(df_calc2,PS_cell_min=min(PsCell,na.rm=TRUE)))
    df_calc2 <- dplyr::mutate(df_calc2,PS_cell_min=if_else(all_na==TRUE,NA,PS_cell_min))
    df_calc2 <- dplyr::ungroup(df_calc2)

    # Keep one set of results for each individual and add to original dataframe
    df_results <- dplyr::filter(df_calc2,ColID1==1)
    df_results <- dplyr::select(df_results,RowID,PS_cell_min)
    df_results <- dplyr::rename(df_results,RowID2=RowID)
    df_results <- dplyr::arrange(df_results,RowID2)

    df <- dplyr::arrange(df,RowID)
    df <- dplyr::bind_cols(df,df_results)
    df <- dplyr::select(df,-starts_with("RowID"))

  } #END if severity==TRUE

  # Set rownames
  rownames(df) <- 1:nrow(df)

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

} #END framework

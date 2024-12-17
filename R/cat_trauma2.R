
#' Categorize trauma data and calculate scores
#'
#' This function adds Abbreviated Injury Scores (AIS), Injury Severity Scores (ISS), and other descriptors of injury to a dataframe.
#' For each observation this function will
#' \enumerate{
#'    \item assign a severity (AIS) and ISS body region values to each valid ICD-10 injury diagnosis code,
#'    \item add variables for maximum severity of each body region,
#'    \item calculate ISS, "New ISS", maximum AIS, and regression-based mortality predictions,
#'    \item select the first 4 mechanism (external cause) codes and categorize mechanism and intent following CDC guidelines
#'}
#'
#' @param df A dataframe in wide format containing ICD-10 diagnosis codes with a common column name prefix.
#'           Diagnosis codes should be character strings and may have a decimal or not.
#'
#' @param dx_pre Prefix for diagnosis code column names (example: dx1, dx2, etc.)
#'
#' @param messages Should the program report completion of each step? Must be TRUE or FALSE (default).
#'          \itemize{
#'          \item TRUE - Messages will report completion of each step (may be helpful for large data sets).
#'          \item FALSE - Messages will not be reported.
#'          }
#'
#'
#' @return A dataframe identical to the dataframe passed to the function with the following additional variables
#'          added:
#'          \itemize{
#'          \item sev_1-sev_n: AIS severity for diagnosis codes 1..n
#'          \item issbr_1-issbr_n: ISS body region for diagnosis codes 1..n
#'          \item mxaisbr1-mxaisbr6: maximum AIS severity for each of the 6 ISS body regions
#'          \item maxais: maximum AIS severity over all ISS body regions
#'          \item riss: computed injury severity score
#'          \item niss: computed "new injury severity score"
#'          \item PmortTQP: TQP model predicted probability of mortality
#'          \item PmortNIS: NIS model predicted probability of mortality
#'          \item mechcode_1-mechcode_4: first 4 mechanism codes found in each row of data
#'          \item mech_1-mech_4: CDC external cause of injury major mechanism for each mechanism code captured
#'          \item intent_1-intent_4: intent for each mechanism code captured
#'          }
#'
#' @details  Data should be in wide format, as in the example below:
#'
#' @examples
#' df_in <- read.table(header = TRUE, text = "
#'     ident   dx1       dx2       dx3
#'     31416   S32110A   S3251     NA
#'     31417   S72141A   T07XXXA   D62
#' ")
#' df_out <- cat_trauma2(df_in, "dx", TRUE)
#'
#' @importFrom stats na.omit
#' @export


cat_trauma2 <- function(df, dx_pre, messages = TRUE) {

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

  ntab <- i10_map_sev[ , c("dx","severity","issbr")]
  rtab <- i10_map_sev[ , c("dx","TQIPeffect","TQIPint","NISeffect","NISint")]

  #---------------------------------------------------------------------------------#
  #  Merge diagnosis code variables with N-Code reference table to obtain severity  #
  #  and ISS body region variables for each diagnosis code and add them to the data #
  #---------------------------------------------------------------------------------#
  for(i in dx_nums){

    if( (messages==TRUE) && (i%%5==0 || i==num_dx) ){
      message("Determining severity of Diagnosis ", i, " of ", num_dx)
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
    temp <- merge(df_ss, ntab, by.x = dx_name, by.y = "dx", all.x = TRUE, all.y = FALSE, sort = FALSE)

    # Reorder rows after merge
    temp <- temp[order(temp$n), ]

    # Reorder columns and drop dx and n
    temp <- temp[ , c("severity","issbr")]

    # Rename columns
    names(temp) <- paste0(c("sev_","issbr_"), i)

    # Add temp columns to dataframe
    df <- .insert_columns(df, dx_name, temp)

  }   #END FOR LOOP (i in dxnums)


  #----------------------------------------------------#
  # Create variables for maximum AIS/ISS body region.  #
  #----------------------------------------------------#
  if(messages==TRUE){
    mindiff=round(as.double(difftime(Sys.time(),starttime,units="secs"))/60)
    message("Time elapsed ", mindiff, " minutes")
    message("Calculating Maximum AIS for each body region")
  }
  # Body regions are coded as text
  body_regions <- unique(i10_map_sev$issbr)
  # Make usable for column names
  issbr_names <- gsub("/", "_", body_regions)

  # For each of the 6 body regions loop through dx codes and get max ais for that body region
  for(i in body_regions){

    # Get severity columns and multiply by 1 if they are for body region i and 0 otherwise
    # This uses element-wise multiplication of matrices
    # All severity columns as a matrix * Indicator matrix of body region columns (entries 1 or 0)
    temp <- df[ , grepl("sev_", names(df)), drop = FALSE] * (1*(df[ , grepl("issbr_", names(df))] == i))

    # Take max (excluding 9) and assign to mxaisbr_i
    # Severity score of 9 implies unknown severity
    # Thus we want to exclude these as long as there is at least one known severity for the body region
    # However if all severity scores for the body region are 9 then we will assign maxaisbr a value of 9
    # by defining a function that converts all zeros to NA,
    # where zeros represent severity values not associated with body region i
    df[ , paste0("mxaisbr_", gsub("/","",i))] <- apply(temp, 1, function(row){
      row <- ifelse(row == 0, NA, row)
      if(all(is.na(row))){
        maxaisbr <- 0
      } else if(all(row == 9, na.rm = TRUE)){
        maxaisbr <- 9
      } else {
        maxaisbr <- max(c(0, row[row != 9]), na.rm = TRUE)
      }
      return(maxaisbr)
    })
  } #END FOR i in body_regions


  #----------------------------------------------------------------------#
  #  Calculate maximum severity over all ISS body regions, excluding 9s  #
  #----------------------------------------------------------------------#
  # Define function to convert 9 to 0
  c9to0 <- function(x) ifelse(x == 9, 0, x)

  # Assign maxais
  # This is a bit complicated
  #   if all maxbr_i are in (9,0,NA) then the max should be 9
  #   if there is at least one positive maxbr_i that is not 9 then we need to exclude the 9s
  # For each row in df...
  df$maxais <- apply(df, 1, function(row){
    # Select mxaisbr columns
    row <- row[grepl("mxaisbr", names(row))]
    # If the max excluding 9 is zero then include 9 so that if there is a 9 then the max will be 9
    if(all(is.na(row))){
      maxais <- as.numeric(NA)
    } else if(max(c9to0(row), na.rm = TRUE) == 0){
      maxais <- max(row, na.rm = TRUE)
    } else {
      maxais <- max(c9to0(row), na.rm = TRUE)
    }
    return(maxais)
  })
  df$maxais <- as.numeric(df$maxais)


  #-----------------------#
  #  Calculate ISS value  #
  #-----------------------#
  if(messages==TRUE){
    message("Calculating Injury Severity Score ")
  }

  # ISS is calculated as the sum of the squared three highest maxaisbr varaiables for a given person.
  # We need to exclude 9s in this calculation.
  df$riss <- apply(df, 1, function(row){
    # Select the max ais variables for a given row
    temp <- row[grepl("^mxaisbr", names(row))]
    # For some reason apply is converting these to char. We will convert them back to numeric
    # Also convert 9s to 0
    temp <- as.numeric(c9to0(temp))
    # Take the three highest, square them, and sum the result
    sum(temp[order(-temp)[1:3]]^2)
  })

  # Replace ISS value with 75 if maximum severity is 6. This implies that the person is dead.
  # 75 is the max ISS score 3(5^2) = 75
  df[df$maxais == 6,"riss"] <- 75

  # Replace ISS value with NA if maximum severity is 9.
  # If maxais is 9, this implies that there were only injuries of unknown severity
  df[df$maxais == 9, "riss"] <- NA


  #------------------------------------------------#
  # Calculate the New Injury Severity Score (NISS) #
  #------------------------------------------------#
  if(messages==TRUE){
    message("Calculating 'New Injury Severity Score' ")
  }
  # NISS is calculated as the sum of the squared three highest AIS variables for a given person.
  # We need to exclude 9s in this calculation.
  df$niss <- apply(df, 1, function(row){
    # select the max ais variables for a given row
    temp <- row[grepl("^sev_", names(row))]
    # convert NA to 0
    temp <- as.numeric(temp)
    temp <- ifelse(is.na(temp) | temp == 9, 0, temp)
    # Take the three highest, square them, and sum the result
    sum(temp[order(-temp)[1:3]]^2)
  })

  # Replace NISS value with 75 if maximum severity is 6. This implies that the person is dead.
  # 75 is the max NISS score 3(5^2) = 75
  df[df$maxais == 6,"niss"] <- 75

  # Replace NISS value with NA if maximum severity is 9.
  # If maxais is 9 this implies that there were only injuries of unknown severity
  df[df$maxais == 9, "niss"] <- NA

  if(messages==TRUE){
    mindiff=round(as.double(difftime(Sys.time(),starttime,units="secs"))/60)
    message("Time elapsed ", mindiff, " minutes")
    message("Calculating predicted mortality from regression models")
  }

  #--------------------------------------------------#
  # Add mortality predictions from regression models #
  #--------------------------------------------------#

  # Duplicate table of diagnoses and convert to long form
  df <- dplyr::mutate(df,RowID=dplyr::row_number())
  df_calc <- df
  df_calc <- dplyr::select(df_calc,RowID,dplyr::starts_with(dx_pre))
  df_calc <- tidyr::pivot_longer(df_calc,cols=tidyr::starts_with(dx_pre),names_to="ColName")
  df_calc <- dplyr::group_by(df_calc,RowID)
  df_calc <- dplyr::mutate(df_calc,ColID1=dplyr::row_number())
  df_calc <- dplyr::ungroup(df_calc)
  df_calc <- dplyr::rename(df_calc,dx=value)
  df_calc <- dplyr::select(df_calc,-ColName)
  # Strip out decimal in all codes, if present
  df_calc <- dplyr::mutate(df_calc,dx=stringr::str_replace(dx,"\\.",""))

  # Merge tables and calculate regression prediction for each individual
  df_merged <- dplyr::left_join(df_calc,rtab,by="dx",relationship="many-to-many")
  df_merged <- dplyr::group_by(df_merged,RowID)
  df_merged <- dplyr::mutate(df_merged,all_na=all(is.na(TQIPeffect)))
  df_merged <- dplyr::mutate(df_merged,TQIPsum=sum(TQIPeffect,na.rm=TRUE))
  df_merged <- dplyr::mutate(df_merged,xTQP=TQIPsum+TQIPint)
  df_merged <- dplyr::mutate(df_merged,PmortTQP=(1/(1+exp(-xTQP))))
  df_merged <- dplyr::mutate(df_merged,NISsum=sum(NISeffect,na.rm=TRUE))
  df_merged <- dplyr::mutate(df_merged,xNIS=NISsum+NISint)
  df_merged <- dplyr::mutate(df_merged,PmortNIS=(1/(1+exp(-xNIS))))
  df_merged <- dplyr::ungroup(df_merged)

  # Keep one set of results for each individual and add to original dataframe
  df_results <- dplyr::filter(df_merged,ColID1==1)
  df_results <- dplyr::select(df_results,RowID,PmortTQP,PmortNIS)
  df_results <- dplyr::rename(df_results,RowID2=RowID)
  df_results <- dplyr::arrange(df_results,RowID2)

  df <- dplyr::arrange(df,RowID)
  df <- dplyr::bind_cols(df,df_results)
  df <- dplyr::select(df,-dplyr::starts_with("RowID"))


  #--------------------------------------------------#
  # Extract mechanism codes and add them to the data #
  #--------------------------------------------------#

  if(messages==TRUE) {
    mindiff=round(as.double(difftime(Sys.time(),starttime,units="secs"))/60)
    message("Time elapsed ", mindiff, " minutes")
    message("Extracting mechanism codes")
  }

  # Obtain table of valid mechanism codes
  etab <- i10_map_mech
  etab <- etab[ , c("dx","mechmaj","intent")]

  # Duplicate table of diagnoses and convert to long form
  df_mech <- df
  df_mech <- dplyr::mutate(df_mech,RowID=dplyr::row_number())
  df_mech <- dplyr::select(df_mech,RowID,dplyr::starts_with(dx_pre))
  df_mech <- tidyr::pivot_longer(df_mech,cols=tidyr::starts_with(dx_pre),names_to="ColName")
  df_mech <- dplyr::group_by(df_mech,RowID)
  df_mech <- dplyr::mutate(df_mech,ColID1=dplyr::row_number())
  df_mech <- dplyr::ungroup(df_mech)
  df_mech <- dplyr::rename(df_mech,dx=value)
  df_mech <- dplyr::select(df_mech,-ColName)
  # Strip out decimal in all codes, if present
  df_mech <- dplyr::mutate(df_mech,dx=stringr::str_replace(dx,"\\.",""))

  # Merge tables and select first four "columns" with mechanism data
  df_merged <- dplyr::left_join(df_mech,etab,by="dx",relationship="many-to-many")
  df_merged <- dplyr::mutate(df_merged,ColID1=dplyr::if_else(is.na(mechmaj),99,ColID1))
  df_merged <- dplyr::group_by(df_merged,RowID)
  df_merged <- dplyr::arrange(df_merged,ColID1)
  df_merged <- dplyr::mutate(df_merged,ColID2=dplyr::row_number())
  df_merged <- dplyr::ungroup(df_merged)
  df_merged <- dplyr::filter(df_merged,ColID2<=4)
  df_merged <- dplyr::mutate(df_merged,dx=dplyr::if_else(is.na(mechmaj),NA,dx))
  df_merged <- dplyr::select(df_merged,-ColID1)
  df_merged <- dplyr::rename(df_merged,mechcode=dx,mech=mechmaj)

  # Convert back to wide form and merge with original dataframe
  df_merged_wide <- tidyr::pivot_wider(df_merged,id_cols=RowID,values_from=c(mechcode,mech,intent),
                                       names_from=ColID2)
  df_merged_wide <- dplyr::arrange(df_merged_wide,RowID)
  df_merged_wide <- dplyr::select(df_merged_wide,-RowID)
  df <- dplyr::bind_cols(df,df_merged_wide)


  if(messages==TRUE){
    mindiff=round(as.double(difftime(Sys.time(),starttime,units="secs"))/60)
    message("Time elapsed ", mindiff, " minutes")
  }

  # Set rownames
  rownames(df) <- 1:nrow(df)

  message("=============================================")
  message("REMINDER")
  message("ICDPICR Version 2.0.4 IS BEING TESTED")
  message("Major bugs and flaws may still exist")
  message("Please report issues to david.clark@tufts.edu")
  message("or at github/clark-david/icdpicr2/issues")
  message("==============================================")

  # Return dataframe
  df

} #END cat_trauma2



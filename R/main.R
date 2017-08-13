#' Return of a list of discretized likelihood curves from a data frame
#'
#' Returns a list of discretized likelihood curves for a specific column of a given dataframe. Determines for each value in that column a specific standard deviation if no general standard deviation is available.
#'
#' The main.R function determines the individual standard deviation through the function sd_computation_from_the_non_linear_regression.R for each value and delivers both to the function discretiziation.R in which the discretisation list is created.
#' The complete discretization list is returned by the main.R function and contains for each value:\cr
#' min = minimal value of the discretiziation\cr
#' max = maximal value of the discretiziation\cr
#' h = step length of the discretiziation\cr
#' probs = dicretized distribution function of the value\cr
#' @param input_table Dataset as a dataframe which should be used for the computations.
#' @param column_For_Discr Numeric column number from which the likelihood curves should be computed.
#' @param start_Row Numeric value which indicates the start of the computation.
#' @param end_Row Numeric value which indicates the end of the computation.
#' @param number_Of_Discretizations Numeric value which defines the number of discretizations for each value.
#' @param h Optional numeric value which defines a step length h for the discretizations of each value.
#' @param halfdist Optional numeric value for the minimum and maximum value of the discretizations. \cr
#'minimal value: value - halfdist \cr
#'maximal value: value + halfdist
#' @param min Optional numeric value for the minimum of the discretizations.
#' @param max Optional numeric value for the maximum of the discretizations.
#' @param general_Sd Optional numeric value which defines a general standard deviation for the discretizations.
#' @examples main (input_table, column_For_Discr =  3, start_Row = 1,
#'       end_Row =  nrow(input_table), number_Of_Discretizations = 200)
#' @export


main <- function(input_table,
                 column_For_Discr,
                 start_Row,
                 end_Row,
                 number_Of_Discretizations,
                 h = NA,
                 halfdist = NA,
                 min = NA,
                 max = NA,
                 general_Sd = NA){

  source ("R/discretization.R")
  source("R/sd_computation_from_the_non_linear_regression.R")

  selected_Rows_For_Computation <- seq (start_Row, end_Row, 1)

  discr_list <- list()

  for (i in seq (start_Row, end_Row, 1)){
    table_value <- input_table[[column_For_Discr]][i]
    # sd from the sd function, if peptide value is out of range the min/max value is used
    if (is.na (input_table[[column_For_Discr]][i]) | is.infinite(input_table[[column_For_Discr]][i])){
      # TODO Falls kein vorhanden ersetzen durch NA, 0 od. LOD 
      table_value <- NA
    } else {
        if (!is.na (general_Sd)) sd <- general_Sd
        else sd <- sd_computation_from_the_non_linear_regression (table_value)
    }

    discr_list[[i]] <- discretization (table_value,
                                       sd,
                                       number_Of_Discretizations,
                                       h,
                                       halfdist,
                                       min,
                                       max)
  }


  #Normalize the probabilities
  for (i in seq (1, length (discr_list), 1)){
    if (!is.null (discr_list[[i]][1])){
      if (!is.na (discr_list[[i]][1]) ){
        #divide probabilities by the increment h
        discr_list[[i]]$probs <- sapply (discr_list[[i]]$probs, function(x) x / discr_list[[i]]$h)
      }
    }
  }


  # name_Of_Returned_Object <- paste ("discr_list", colnames (input_table [column_For_Discr]), sep = "_")
  name_Of_Returned_Object <- paste ("discr_list", column_For_Discr, sep = "_")
  assign (name_Of_Returned_Object, discr_list, envir = .GlobalEnv)

  #return(discr_list)
}

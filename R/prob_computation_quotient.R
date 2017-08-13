#' Computation of the discretized probabilities of each dataset value
#'
#' Used by the function discretization.R. Returns for each value a list with the discretized likelihood curve based on the range, step length h and the standard deviation.
#' @param input_table_value Value for the discretization.
#' @param sd Dynamic or general standard deviation.
#' @param h1 Chosen or computed step length h.
#' @param halfdist Optional numeric value for the minimum and maximum value of the discretizations.
#' @param min Optional numeric value for the minimum of the discretizations.
#' @param max Optional numeric value for the maximum of the discretizations.


prob_computation_quotient <- function (input_table,
                                       rowName,
                                       prot_Col = NA, # first Column
                                       prot_Col2 = NA, # second Column
                                aggregation_Mode = 2, # multiplicative
                                steps_For_Aggregation = NA,
                                ){
  
  #uncertainty_curve from uncertainty_curve_X
  name <- paste ("uncertainty_curve", rowName, sep = "_")
  uncertainty_curve <- eval (as.symbol (name))
  
  
  #number of the column used as legend
  col_Index <- which (colnames (input_table) == prot_Col)
  col_Index2 <- which (colnames (input_table) == prot_Col2)

  steps <- c()
  probs <- c()
  
  for (i in seq (1, length (col_Index), 1)){
    for(j in steps_For_Aggregation){
      # quotient
      steps[j] <- uncertainty_curve[col_Index2[i]]$steps[j] - uncertainty_curve[col_Index[i]]$steps[j]
      
      # additional
      if (aggregation_Mode == 1){
        
      }
      # multiplicative
      if (aggregation_Mode == 2){
        probs[j] <- uncertainty_curve[col_Index2[i]]$probsList[j] * uncertainty_curve[col_Index[i]]$probsList[j]
      }
    }
  }
  
  probsList <- list(probs)
  
  return(probsList)
}


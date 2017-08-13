#' Plotting of the uncertainty of dataset values
      #'
      #' The function can plot the uncertainty of the chosen dataset values in two possible ways. First the superior element like aggregation and second the row like aggregation.
      #' The superior element like aggregation plots all values from the chosen columns related to a given superior element. Whereas the row like aggregation plots all
      #' values from the chosen columns related to a given subordinate element in a specific row. Furthermore a list is returned containing the discretized values of the uncertainty curves. This list is used for the uncertainty aggregation.
      #' @param steps_Col1 Dataset as a dataframe which should be used for the computations.
      #' @param steps_Col2 Name of the superior element from which the values should be plotted (superior element like aggregation).
      #' @param probs_Col1 Name of the subordinate element from which the values should be plotted (Row like aggregation).
      #' @param probs_Col2 Column name of the superior element column (superior element like aggregation).
      #' @param col_Legend Column name of the column for legend / column with subordinate elements.
      #'If rowNames of the col_Legend are not unique, you have to add an additional
      #'colum to the input_table with the row numbers:
      #' likelihood_curves (input_table, rowName = "7",col_Legend = "Row_Numbers",
      #'                    prot_Col2 = c("3min", "6min"), aggregation_Type = 2)
      #'
      #' @export

quotient_utils <- function (uncertainty_curve1, 
                            uncertainty_curve2){
      
      number_Of_Columns <- length(uncertainty_curve1$steps)
      number_Of_Columns2 <- length(uncertainty_curve2$steps)

      increment <- uncertainty_curve1$h*uncertainty_curve2$h
      
      steps_List <- list()
      prob_List <- list()
      
      k <- 1
      for (i in seq(1, number_Of_Columns2, 1)) {
            for (j in seq(1, number_Of_Columns, 1)) {
              steps_List[k] <- round(uncertainty_curve2$steps[i] - uncertainty_curve1$steps[j], 2)
              prob_List[k] <- uncertainty_curve2$prob_List_Normalized[i]*uncertainty_curve1$prob_List_Normalized[[j]]
              k <- k + 1
            }
      }
    
      
      DT <- data.table(steps = unlist(steps_List), probs = unlist(prob_List))[order(-steps, decreasing = TRUE)]
      # sum up over the same steps values
      DT <- DT[, lapply(.SD, sum), by = steps]
      
      prob_List_Normalization_Factor <- 1 / sum (unlist (DT$probs)) # sum should be 1
      prob_List_Normalized0 <- lapply (DT$probs, function (x) x * prob_List_Normalization_Factor)
      prob_List_Normalized <- lapply (prob_List_Normalized0, function (x) x / increment)

      x <- unlist(prob_List_Normalized)

      uncertainty_quotient <- list(steps = DT$steps, prob_List_Normalized = x, h = increment)
}
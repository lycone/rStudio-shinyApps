#' Return of the range of a dataset value
#'
#' The uncertainty_evaluation of each value happens in a specific range related to the value peak. This function returns this
#' range in steps related to the step length h of the value with a maximum and minimum of approximately 1/100 of the peak (from the uncertainty curve of that specific value).
#' @param input_table Dataset as a dataframe which should be used for the computations.
#' @param discr_list Discretization list of the used values (see main.R for more information).
#' @param column_Number Column_Number of the values.
#' @param peptide_Row Row number of the values.
#' @param protein_Rows Position of the values in the generated discretization list (only important for the row like aggregation).
#' @param increment Step length h.
#' @param general_Sd Optional general standard deviation.


range_of_uncertainty <- function (input_table,
                                  discr_list,
                                  column_Number,
                                  peptide_Row,
                                  protein_Rows,
                                  increment,
                                  general_Sd){


#peptide value from the table for the peak value
peptide_Value <- input_table[[column_Number]][peptide_Row]

#steps h from peptide
steps <- seq (discr_list[[protein_Rows]]$min, discr_list[[protein_Rows]]$max, by = discr_list[[protein_Rows]]$h)

#sd of the mid value
if (!is.na(general_Sd))
  sd_Step <- general_Sd
if (is.na(general_Sd)){
  sd_Step <- sd_computation_from_the_non_linear_regression (steps[length(steps)/2])
}

#peak value of the standard distribution
peak_value <- pnorm (input_table[[column_Number]][peptide_Row] + increment, peptide_Value, sd_Step) - pnorm(input_table[[column_Number]][peptide_Row], peptide_Value, sd_Step)

#min / max range of the farthest right / left point ---- ~ 1/100 of the peak
peak_value_min <- peak_value / 120
peak_value_max <- peak_value / 98

#list with prob from the min and max step
prob_List_Determine_Min_Max <- list()
steps_Min_Max <- c()

a <- T
b <- T
k <- 0
while (a == T || b == T){

  k <- k + 1 #terminates the loop, if x times pass (no endless loop, because of a or b)

  steps_Min_Max[1] <- steps[1]
  steps_Min_Max[2] <- steps[length(steps)]

  #first and last prob_value
  for ( i in seq (1, 2, 1)){
    if (i == 1){
      if (!is.na(general_Sd))
        sd_Step <- general_Sd
      if (is.na(general_Sd)){
        sd_Step <- sd_computation_from_the_non_linear_regression (steps[1]) # sd for first step
      }

    }
    if (i == 2){
      if (!is.na(general_Sd))
        sd_Step <- general_Sd
      if (is.na(general_Sd)){
        sd_Step <- sd_computation_from_the_non_linear_regression (steps[length(steps)]) # sd for last step
      }

    }
    # if one of these NULL -> value out of the regression range
    prob_List_Determine_Min_Max[i] <- (pnorm (input_table[[column_Number]][peptide_Row] + increment, steps_Min_Max[i], sd_Step) - pnorm(input_table[[column_Number]][peptide_Row], steps_Min_Max[i], sd_Step))
  }


  # MIN
  if (prob_List_Determine_Min_Max[1] >= peak_value_min && prob_List_Determine_Min_Max[1] <= peak_value_max){
    a <- F
  }
   else if (prob_List_Determine_Min_Max[1] >= peak_value_max){
    step_1 <- steps[1] - increment
    steps <- append (step_1, steps)
  }
  else {
    steps <- steps[-1]
  }


  # MAX
  if (prob_List_Determine_Min_Max[2] >= peak_value_min && prob_List_Determine_Min_Max[2] <= peak_value_max){
    b <- F
  }
  else if (prob_List_Determine_Min_Max[2] >= peak_value_max){
    steps <- append (steps, steps[length(steps)] + increment)
  }
  else{
    steps <- steps[-length(steps)]
  }

  # terminates the loop, if x times pass
  if (k == 1000){
    a <- F
    b <- F
  }
}


return (steps)
}

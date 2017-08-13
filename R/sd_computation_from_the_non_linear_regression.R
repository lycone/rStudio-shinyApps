#' Determination of the dynamic standard deviation for one dataset value
#'
#' Based on the dynamic standard deviation regression curve a standard deviation is determined for the dataset values if no general standard deviation is chosen.
#' @param input_table_value Value from which the standard deviation is determined.


sd_computation_from_the_non_linear_regression <- function (input_table_value){


  if (input_table_value >= max(median_List)){
    sd <- predict (lo, max(median_List))
  }
  else if (input_table_value <= min(median_List)){
    sd <- predict (lo, min(median_List))
  }
  else sd <- predict (lo, input_table_value)

  return (sd)
}

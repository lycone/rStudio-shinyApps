#' Creation of the Discretization List
#'
#' Creates the list with the discretized likelihood curves for each given dataset value. Uses the function h_determination.R for the computation of the steplength h if no step length is chosen and the function prob_computation.R for the computation of the specific discretization values.
#'
#'If no range is chosen via min/max or halfdist, the halfdist is determined by 4 times the standard deviation.
#'
#' Each dataset value list contains:	\cr
#' min = minimal value of the discretization \cr
#' max = maximal value of the discretization \cr
#' h = step length of the discretization \cr
#' probs = dicretized distribution function of the value
#'
#' @param input_table_value Value for the discretization.
#' @param sd Standard deviation for the computation of the likelihood function.
#' @param number_Of_Discretizations Numeric value which defines the number of discretizations for each value.
#' @param h Optional numeric value which defines a step length h for the descitizations of each value.
#' @param halfdist Optional numeric value for the minimum and maximum value of the descritizations.
#' @param min Optional numeric value for the minimum of the descritizations.
#' @param max Optional numeric value for the maximum of the descritizations.



discretization <- function (input_table_value,
                            sd,
                            number_Of_Discretizations,
                            h = NULL,
                            halfdist = NULL,
                            min = NULL,
                            max = NULL){

  source ("R/h_determination.R")
  source ("R/prob_computation.R")

  # TODO Inf-Werte sind sonder zu behandeln?
  if (is.na(input_table_value) | is.infinite(input_table_value)){# skip NA values
    return (NA)
  }
  
  if (is.na (halfdist))
    if (is.na (min) || is.na (max))
      halfdist <- 4*sd    #default halfdist -> 4*sd

  #halfdist available
  if (!is.na (halfdist)){

    min0 <- round (input_table_value - halfdist, digits = 2)
    minround <- round (min0, digits = 1)
    if (abs (minround - min0) < 0.03){
      if (minround < min0) min1 <- minround + 0.05
      if (minround >= min0) min1 <- minround - 0.05
    }
    else min1 <- minround

    max1 <- input_table_value + halfdist

    h1 <- h_determination (h, halfdist, min1, max1, number_Of_Discretizations)

  }

  #min/max available
  if (!is.na (min) && !is.na (max)){
    min1 <- min
    max1 <- max
    h1 <- h_determination (h, halfdist, min, max, number_Of_Discretizations)
  }


#prob computation
  probsList <- prob_computation (input_table_value, sd, h1, halfdist, min1, max1)
  discr_list <- as.list( c (min1,max1,h1) )
  discr_list[4] <- probsList
  names(discr_list) <- c ("min","max","h","probs")

  return (discr_list)
}

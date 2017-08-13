#' Computation of the steplength h
#'
#' Used by the function discretization.R. Computes the step length h for the discretized likelihood curves, if no step length is chosen.
#'
#' The computation is based on range and number_Of_Discretizations.
#' A rounding is performed after the first computation to ensure consistent discretization values for the performed aggregation during the uncertainty evaluation.
#' @param h Optional numeric value which defines a step length h for the discretizations of each value.
#' @param halfdist Optional numeric value for the minimum and maximum value of the discretizations.
#' @param min Optional numeric value for the minimum of the discretizations.
#' @param max Optional numeric value for the maximum of the discretizations.
#' @param number_Of_Discretizations Numeric value which defines the number of discretizations for each dataset value.

h_determination <- function (h,
                             halfdist,
                             min,
                             max,
                             number_Of_Discretizations){


  if (!is.na (h)){
    h1 <- h
  }
  # TODO Wurde für Weiterberechnung mit generierten Daten gesetzt. Gehört später weg!
  h1bearb <- 5
  
  if (is.na (h)){
    if (!is.na (halfdist)){
      h1 <- halfdist *2 / number_Of_Discretizations
    }
    if (!is.na (min) && !is.na (max)){
      h1 <- abs (max-min)/number_Of_Discretizations
    }
    if (h1 < 0.0025 && h1 > 0.00001){
      h1bearb <- 0.001
    }
    if (h1 < 0.0075 && h1 >= 0.0025){
      h1bearb <- 0.005
    }
    if (h1 < 0.025 && h1 >= 0.0075){
      h1bearb <- 0.01
    }
    if (h1 < 0.075 && h1 >= 0.025){
      h1bearb <- 0.05
    }
    if (h1 < 0.25 && h1 >= 0.075){
      h1bearb <- 0.1
    }
    if (h1 < 0.75 && h1 >= 0.25){
      h1bearb <- 0.5
    }
    if (h1 < 2.5 && h1 >= 0.75){
      h1bearb <- 1
    }
    if (h1 < 100 && h1 >= 2.5){
      h1bearb <- 5
    }

    h1 <- h1bearb
  }

  return (h1)
}

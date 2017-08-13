#' Computation of the discretized probabilities of each dataset value
#'
#' Used by the function discretization.R. Returns for each value a list with the discretized likelihood curve based on the range, step length h and the standard deviation.
#' @param input_table_value Value for the discretization.
#' @param sd Dynamic or general standard deviation.
#' @param h1 Chosen or computed step length h.
#' @param halfdist Optional numeric value for the minimum and maximum value of the discretizations.
#' @param min Optional numeric value for the minimum of the discretizations.
#' @param max Optional numeric value for the maximum of the discretizations.


prob_computation <- function (input_table_value,
                              sd,
                              h1,
                              halfdist,
                              min,
                              max){

  if (!is.na (halfdist)){
    steps <- seq ( round (input_table_value - halfdist, digits = 2), input_table_value + halfdist, by = h1) # define intervals
  }
  if (!is.null (min) && !is.null (max)){
    steps <- seq (min, max, by = h1)
  }

  probs <- c()
  probs[1] <- pnorm(steps[1], input_table_value, sd) # calculation 1. prob

  for (i in seq (2, length (steps) - 1, 1)){ # calculation probs,  2. to number of probs-1 in single Steps
    probs[i] <- (pnorm (steps[i+1], input_table_value, sd) - pnorm(steps[i], input_table_value, sd))
  }

  probs[length(steps)] <- 1 - pnorm(steps[length(steps)-1], input_table_value, sd) # calculation last prob 1-penultimate

  probsList <- list (probs)

  return (probsList)
}
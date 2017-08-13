#' Dynamic standard deviation
#'
#' Generates a dynamic standard deviation based on specified columns of the input table. The dynamic standard deviation is represented by the means of the standard deviation of several areas. The areas contain the means of the rows from the chosen columns and are determined by the parameter n (Amount of values per SD area). The function returns a list with the standard deviations of the different sd areas and the related regression curve as an object of class "loess". Also creates a plot with the dynamic standard deviation curve.
#' @param input_table Dataset as a dataframe which should be used for the computations.
#' @param n Amount of peptides per SD area.
#' @param columns_For_Sd A numeric vector which contains the column numbers from which the sd should be computed.
#' @param use.mad Logical value. TRUE returns the estimation of the standard deviation based on the MAD. FALSE returns the estimation of the standard deviation based on the standard deviation of the single dataset values. (default = TRUE).
#' @param xmin Minimal range of x axis.
#' @param xmax Maximal range of x axis.
#' @param ymin Minimal range of y axis.
#' @param ymax Maximal range of y axis.
#' @examples dynamic_standard_deviation (input_table, n = 100, columns_For_Sd = c (3, 4, 5), use.mad = F,
#'                             xmin = -1, xmax = 1.5, ymin = 0, ymax = 1.5)
#' @export


dynamic_standard_deviation <- function (input_table,
                                        n,
                                        columns_For_Sd,
                                        columns_For_Sd2,
                                        use.mad = T,
                                        xlab = "Intensity Value",
                                        ylab = "Standard Deviation",
                                        xmin = NA,
                                        xmax = NA,
                                        ymin = NA,
                                        ymax = NA){

  source("R/sd_sample_MAD.R")

  #new dataframe consisting of the choosen columns for the sd computation
  # table1 <- input_table[, columns_For_Sd]
  
  tableTemp <- rbindlist(list(input_table[, columns_For_Sd, with = FALSE], 
                           input_table[, columns_For_Sd2, with = FALSE]))
  
  number_Of_Columns_For_Sd <- seq (1,length(columns_For_Sd), by = 1)

  # Data Preparation
  # replace none finite values with NA
  tableTemp[mapply(is.infinite, tableTemp)] <- NA
  table1 <- tableTemp[complete.cases(tableTemp)]
  
  #mean computation for rows
  # table1_means <- rowMeans(table1, na.rm = T)
  table1_means <- rowMeans(table1)
  
  x <- sort(table1_means, decreasing = F) # means sortiert
  
  #Number of different sd computations (number of rows(peptides) / number of choosen peptides per sd area )
  cycle_TIME <- as.integer (length(x)/n)
  
  remainder <- length(x) %% n
  
  # Mean breakpoints for sd area classification
  breakpoints = NULL
  for (i in seq (1, cycle_TIME, 1)){
    if (i < as.integer (cycle_TIME/2)){
      breakpoints <- append(breakpoints,x[n*i])
    }
    else{
      breakpoints <- append(breakpoints,x[n*i+remainder])
    }
  }

 # Data Clustering
 # divide the dataset into different tables depending on the means of the peptide intensities/regulation factors
 for (i in seq (1, cycle_TIME, 1)){
  name <- paste ("intensity_Window", i, sep = "")

  if (i == 1){
    assign (name, subset(table1, rowMeans(table1, na.rm = T) <= breakpoints[i]))
  }

  if (i > 1 && i <= cycle_TIME){
    assign (name, subset(table1, rowMeans(table1, na.rm = T) <= breakpoints[i] & rowMeans(table1, na.rm = T) > breakpoints[i-1]))
  }
 }


 #intensity_Window1
 # compute different sd's depending on the means of the peptide intensities
 sd <- list()
 for (i in seq (1, cycle_TIME, 1)){
  name <- paste ("intensity_Window", i, sep = "")
  sd[i] <- sd_sample_MAD(eval(as.symbol(name)), number_Of_Columns_For_Sd, use.mad)
 }


 # compute median of the different tables for plotting against their sd
 list1 <- c()
 median_List <- c()
 for (i in seq (1,cycle_TIME, 1)){
  name <- paste ("intensity_Window", i, sep = "")

  for (j in seq (1, length(number_Of_Columns_For_Sd), 1)){
    list1 <- append (list1, eval(as.symbol(name))[[j]])
  }

  median_List[i] <- median (list1, na.rm = T)
  list1 <- c()
 }


 #median_List
 sd_List <- unlist(sd)
  
 if(is.na(xmin) | is.null(xmin)){
  xmin <- min(median_List)
 }
 if(is.na(xmax) | is.null(xmax)){
  xmax <- max(median_List)
 }
 
 if(is.na(ymin) | is.null(ymin)){
   ymin <- min(sd_List)
 }
 if(is.na(ymax) | is.null(ymax)){
   ymax <- 2*max(sd_List)
 }
 
 sd_plot <<- plot(median_List, sd_List, xlim = c(xmin, xmax), ylim = c(ymin, ymax), xlab = xlab, ylab = ylab)

 lo <<- loess(sd_List~median_List) #loess function to global environment (used by other functions)
 xl <- seq(min(median_List),max(median_List), (max(median_List) - min(median_List))/1000)
 lines(xl, predict(lo,xl), col='red', lwd=2)

 median_List <<- median_List
}

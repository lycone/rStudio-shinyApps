#' Plotting of the uncertainty of dataset values
#'
#' The function can plot the uncertainty of the chosen dataset values in two possible ways. First the superior element like aggregation and second the row like aggregation.
#' The superior element like aggregation plots all values from the chosen columns related to a given superior element. Whereas the row like aggregation plots all
#' values from the chosen columns related to a given subordinate element in a specific row. Furthermore a list is returned containing the discretized values of the uncertainty curves. This list is used for the uncertainty aggregation.
#' @param input_table Dataset as a dataframe which should be used for the computations.
#' @param Peptide_Name Name of the subordinate element from which the values should be plotted (Row like aggregation).
#' @param Peptide_Name_List Column name of the superior element column (superior element like aggregation).
#' @param col_Legend Column name of the column for legend / column with subordinate elements.
#' @param prot_Col Column names which contains the values that should be plotted (superior element like aggregation).
#' @param prot_Col Column names which contain the values that should be plotted (Row like aggregation).
#' @param title Optional character for the headline of the plot.
#' @param general_Sd Optional numeric value which defines a general standard deviation for the uncertainty evaluation.
#' @param title Optional character for the headline of the plot.
#' @param ymax Maximal range of y axis.
#' @param xlabel Label of the x axis.
#' @param ylabel Label of the y axis.
#' @param aggregation_Type 1 = superior element like aggregation, 2 = Row like aggregation.
#' @param plotF Optional logical value. TRUE = plot returned, FALSE = no plot returned  (default = TRUE).
#' @examples superior element like aggregation:
#' uncertainty_eval (input_table, proteinName = "Protein_X", Peptide_Name_List = "Protein",
#' col_Legend = "Peptide", prot_Col = c("3min", "6min"), aggregation_Type = 1)
#'
#' @export


uncertainty_eval <- function (input_table,
                              Peptide_Name = NA,
                              Peptide_Name_List = NA,
                              col_Legend,
                              prot_Col = NA,
                              general_Sd = NA,
                              title = "",
                              ymax = 1, # max probability doesn't surpace 1
                              xLabel = "x",
                              yLabel = "y",
                              aggregation_Type,
                              plotF = T){

  source ("R/range_of_uncertainty_evaluation.R")

  discr_list <- list()
  legend <- c()
  
  name_Of_Returned_Object <- c()
  
  #colum numbers from choosen replicates
  protein_Column <- c()

  if (aggregation_Type == 2){ # Row like Aggregation
    name_Of_Returned_Object <- paste ("uncertainty_curve", Peptide_Name, sep = "_")
    
    # how much disccritizations for loop and for the ifs
    protein_Rows <- seq (1, length (prot_Col), 1)
    
    #Row from the selected Peptide
    peptide_Row <- which(input_table[[1]] == Peptide_Name)
    peptide_Row <- peptide_Row[1]
    
    #colum numbers from choosen replicates
    
    for (i in seq (1, length(prot_Col), 1)){
      protein_Column[i] <- which (colnames (input_table) == prot_Col[i])
    }
    
    
    column_Number <- c()

    for (i in seq (1, length (prot_Col), 1)){
      name <- paste ("discr_list", prot_Col[i], sep = "_")
      
      if (!is.null (eval(as.symbol(name)) [[peptide_Row]] [1] )){
        if (!is.na(eval(as.symbol(name))[peptide_Row] [[1]] [1])){

          discr_list[i] <- eval(as.name(name))[peptide_Row]
          legend <- append(legend, as.character(prot_Col[i]))
        } else{

          discr_list[i] <- NA
        }
      
        column_Number <- append (column_Number, protein_Column[i])
      }
    }
    
    #title for plot
    if (title == ""){
      title <- Peptide_Name
    }
  }
  
  if(aggregation_Type == 1){

    protein_Rows <- NULL
    
    i <- 1
    
    for(proteinName in Peptide_Name_List){
      #row number in input_table / discr_list
      protein_Row <- which(input_table[[1]] == proteinName)
      protein_Rows[i] <- protein_Row[1]

      i <- i + 1
    }

    
    peptide_Row <- c() # for determine_range_of_values...()
    column_Number <- c() # for determine_range_of_values...()
    
    for (i in seq (1, length(prot_Col), 1)){
      protein_Column[i] <- which (colnames (input_table) == prot_Col[i])
      
      name <- paste ("discr_list", prot_Col[i], sep = "_") #"V1_3min"

      for (j in seq (1, length (protein_Rows), 1)){
        
        if (protein_Rows [j] <= length (eval(as.symbol(name))) ){
          if (!is.null (eval(as.symbol(name)) [[protein_Rows[j] ]] [1]) &&
              !is.na(eval(as.symbol(name))[protein_Rows[j]] [[1]] [1])){
              
            discr_list[j] <- eval(as.symbol(name))[protein_Rows[j]]
            legend <- append (legend, as.character(input_table[[col_Legend]][protein_Rows[j]]))
            name_Of_Returned_Object <- paste ("uncertainty_curve", Peptide_Name_List[i], sep = "_")
          }
          else{

            discr_list[j] <- NA
          }
          
          peptide_Row <- append (peptide_Row, protein_Rows[j])
          column_Number <- append (column_Number, protein_Column[i])
        }
      }
    }

    # TODO Fehler bei length(discr_list) = 0 (zB. Einziger discr_list = NA) 
    # numbers for new discr_list
    protein_Rows <- seq (1, length(discr_list), 1)

    #title for plot
    if (title == ""){
      title <- prot_Col # Column name
    }
  }


  #--------------------------------------------------------
  # aggregation_Type 1 and 2

  # first time create new plot
  ploti <- 0
  #number of non NA's for colors
  color_count <- length(legend)
  #corlornumber from colorcount
  color <- 0


  steps <- vector("list", length (discr_list))
  xmin <- NULL
  xmax <- NULL
  ymax <- NULL
  
  for (i in seq (1, length(discr_list), 1)){
    #multiple column_Numbers
    if(is.na(discr_list[i])){
      next
    }

    if (aggregation_Type == 2){
      increment <- discr_list[[i]]$h
      #Min, Max Determination of the useful Steps
      steps[[i]] <- range_of_uncertainty (input_table,
                                          discr_list,
                                          column_Number[i], #unterschiedliche columns for the Max value
                                          peptide_Row, # always the same peptide_Row
                                          protein_Rows[i], # but different "protein_Rows" in discr_list
                                          increment,  #peptide Row != position in Replicate Aggregation
                                          general_Sd)
    }
    
    #one column_Number
    if (aggregation_Type == 1){
      increment <- discr_list[[i]]$h
      #Min, Max Determination of the useful Steps
      steps[[i]] <- range_of_uncertainty (input_table,
                                          discr_list,
                                          column_Number[i], #unterschiedliche columns for the Max value
                                          peptide_Row[i],  # same as protein_Rows
                                          protein_Rows[i], # unterschiedliche protein Rows
                                          increment,
                                          general_Sd)
    }
  }

  uncertainty_curve <- list()

  for (i in seq (1, length(discr_list), 1)){
    
    if(is.na(discr_list[i])){
      
      uncertainty_curve[[i]] <- NA
      next
    }

    increment <- discr_list[[i]]$h
    sd_Step <- NULL
    prob_List <- list()

      #calculation probs,  2. to number of probs-1 in single Steps
    if (aggregation_Type == 1)
      for (j in seq (1, length (steps[[i]]), 1)){
        if (!is.na(general_Sd))
          sd_Step <- general_Sd
        if (is.na(general_Sd)){
            sd_Step <- sd_computation_from_the_non_linear_regression (steps[[i]][j])
          }
        prob_List[j] <- (pnorm (input_table[[column_Number[i]]][peptide_Row[i]] + increment, steps[[i]][j], sd_Step) - pnorm(input_table[[column_Number[i]]][peptide_Row[i]], steps[[i]][j], sd_Step))
      }

    if (aggregation_Type == 2)
      for (j in seq (1, length (steps[[i]]), 1)){
        if (!is.na(general_Sd))
          sd_Step <- general_Sd
        if (is.na(general_Sd)){
          sd_Step <- sd_computation_from_the_non_linear_regression (steps[[i]][j])
        }

        prob_List[j] <- (pnorm (input_table[[column_Number[i]]][peptide_Row] + increment, steps[[i]][j], sd_Step) - pnorm(input_table[[column_Number[i]]][peptide_Row], steps[[i]][j], sd_Step))
      }

    #Normaliaztion
    #change the sum of the probs to 1
    # sum_prob <- sum (unlist (prob_List)) = 1
    prob_List_Normalization_Factor <- 1 / sum (unlist (prob_List)) # sum should be 1
    prob_List_Normalized0 <- lapply (prob_List, function (x) x * prob_List_Normalization_Factor)

    prob_List_Normalized <- lapply (prob_List_Normalized0, function (x) x / increment)

    #save in list for Aggregation
    x <- unlist(prob_List_Normalized)
    # x <- unlist(prob_List_Normalized0)

    uncertainty_curve[[i]] <- list (steps = steps[[i]], prob_List_Normalized = x, h = increment)
    
    if(is.null(xmax)){
      xmax <- max(uncertainty_curve[[i]]$steps)
    }
    if(is.null(xmin)){
      xmin <- min(uncertainty_curve[[i]]$steps)
    }
    if(is.null(ymax)){
      ymax <- max(uncertainty_curve[[i]]$prob_List_Normalized)
    }

    xmin <- min(xmin, min(uncertainty_curve[[i]]$steps))
    xmax <- max(xmax, max(uncertainty_curve[[i]]$steps))
    ymax <- max(ymax, max(uncertainty_curve[[i]]$prob_List_Normalized))
}

if(plotF){
  ploti <- 0
  #number of non NA's for colors
  color_count <- length(legend)
  #corlornumber from colorcount
  color <- 0
  
  for (i in seq(1, length(uncertainty_curve), 1)) {
    if(is.na(uncertainty_curve[[i]][1])){
      
      next # skip uncertainty_quotient
    }
    
    ploti <- ploti + 1
    color <- color + 1
    
    if(ploti == 1){ # first time create new plot
      plot (uncertainty_curve[[i]]$steps, uncertainty_curve[[i]]$prob_List_Normalized, type = "l", main = title, xlim = c(xmin, xmax), ylim = c(0, ymax), xlab = xLabel, ylab = yLabel, col = rich.colors(color_count)[color])
    }
    if (ploti > 1){
      lines (uncertainty_curve[[i]]$steps, uncertainty_curve[[i]]$prob_List_Normalized, type = "l", col = rich.colors(color_count)[color])
    }
  }
  
  legend("topright", legend = legend, col = rich.colors(color_count), pch = 20, cex = 0.7)  
}

if(aggregation_Type == 2) # only save if aggregation_Type because complete
  assign (name_Of_Returned_Object, uncertainty_curve, envir = .GlobalEnv)
}
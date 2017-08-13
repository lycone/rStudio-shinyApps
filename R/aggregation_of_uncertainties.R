#' Aggregation of the uncertainty curves
#'
#' The uncertainty curves of different dataset values can be aggregated to a single curve. This can be realized in an additive or multiplicative way. For a subsequent filtering, the
#' additive aggregation is recommended, since the filtering is based on multimodal curves. Unimodal curves, which are generally created by the multiplicative aggregation, won't pass the filter.
#' @param input_table Dataset as a dataframe which should be used for the computations.
#' @param proteinName Name of the superior element from which the values should be plotted (superior element like aggregation).
#' @param rowName Name of the subordinate element from which the values should be plotted (Row like aggregation).
#' @param name_Of_SupEl_Col Column name of the superior element column (superior element like aggregation).
#' @param col_Legend Column name of the column for legend / column with subordinate elements.
#' @param prot_Col Column names which contains the values that should be plotted (superior element like aggregation).
#' @param prot_Col2 Column names which contain the values that should be plotted (Row like aggregation).
#' @param title Optional character for the headline of the plot.
#' @param xlabel Label of the x axis.
#' @param ylabel Label of the y axis.
#' @param aggregation_Mode Optional numeric number. 1 = additively or 2 = multiplicative aggregation of the likelihood curves (default = 1).
#' @param aggregation_Type Numeric number. 1 = superior element like aggregation, 2 = Row like aggregation.
#' @param plotA Optional logical value. T = plot aggregation curve, F = return aggregation_curve as numeric vector for filtering with the function filter_aggregated_values.R (default = TRUE).
#' @examples superior element like aggregation, plot return:
#' aggregation_of_uncertainties (input_table, proteinName = "Protein_X", name_Of_SupEl_Col = "Protein",
#'         col_Legend = "Peptide", prot_Col = c("3min", "6min"), aggregation_Type = 1)
#'
#' superior element like aggregation, aggregation curve return:
#' aggregation_curve <- aggregation_of_uncertainties (input_table, proteinName = "Protein_X",
#'         name_Of_SupEl_Col = "Protein", col_Legend = "Peptide", prot_Col = c("3min", "6min"), aggregation_Type = 1, plotA = F)
#'
#' row like aggregation:
#' uncertainty_aggr(rowName = "Peptide_X", col_Legend = "Peptide", prot_Col2 = c ("3min", "6min"), aggregation_Type = 2)
#'
#' If rowNames of the col_Legend are not unique, you have to add an additional
#' colum to the input_table with the row numbers:
#' uncertainty_aggr(rowName = "7", col_Legend = "Row_Numbers", prot_Col2 = c ("3min", "6min"), aggregation_Type = 2)
#' @export


aggregation_of_uncertainties <- function (input_table,
                              Peptide_Name,
                              prot_Col = NA,
                              prot_Col_First = NA,
                              prot_Col_Second = NA,
                              title = "",
                              xLabel = "x",
                              yLabel = "y",
                              aggregation_Type = 1, # elements or quotients
                              aggregation_Mode = 1, # Additiv - Multiplicatikv
                              plotA = T){
  
  source("R/uncertainty_evaluation.R")
  source("R/uncertainty_quotient.R")
  
  uncertainty_aggragation <- list()
  name_Of_Returned_Object <- c()

  legend <- c()

  peptide_Column <- c()

  peptide_Column_First <- c()
  peptide_Column_Second <- c()
  
  xmin <- NULL
  xmax <- NULL
  ymax <- NULL

  steps_List <- list()
  probs_List <- list()

  # ------------- Elements
  # -------------
  if(aggregation_Type == 1){ # aggregation for elements
    name_Of_Returned_Object <- paste("intensity_aggragation", Peptide_Name, sep = "_")

    #Row from the selected Peptide
    peptide_Rows <- which(input_table[[1]] == Peptide_Name)
    peptide_Row <- peptide_Rows[1] # if multiple rows with same name pick the first one (should not appen!)
    
    name <- paste("uncertainty_curve", Peptide_Name, sep = "_")
    
    uncertainty_curve <- NULL
    
    tryCatch({
      message("aggregation_of_uncertainties: loading ", name, " (aggregation_Type == 1)")
      
      uncertainty_curve <- eval(as.symbol(name))
    },
    
    error = function(name){
      
      message("aggregation_of_uncertainties: Error while loading ", name, " (aggregation_Type == 1)")
    },
    
    warning = function(name){
      
      message("aggregation_of_uncertainties: Warning while loading ", name, " (aggregation_Type == 1)")
    },
    
    finally = {
      if(is.null(uncertainty_curve)){
        message("aggregation_of_uncertainties: calculating ", name, " (aggregation_Type == 1)")
        
        uncertainty_eval( input_table = input_table,
                          Peptide_Name = Peptide_Name,
                          col_Legend = names(input_table[1,1]),
                          Peptide_Name_List = NA,
                          prot_Col = prot_Col,
                          general_Sd = NA,
                          aggregation_Type = 2,
                          plotF = F)
        
        message("aggregation_of_uncertainties: reloading ", name, " (aggregation_Type == 1)")
        uncertainty_curve <- eval(as.symbol(name))
      }
    }
    )

    increment <- 0
    
    for(i in seq(1, length(prot_Col), 1)){
      peptide_Column <- which(colnames(input_table) == prot_Col[i])

      if(!is.null(uncertainty_curve[[i]][1])){ # uncertainty_curve exist?
        # -1 to skip the Sequence index of Peptids
        if(!is.na(uncertainty_curve[[peptide_Column - 1]][1])){

          steps_List <- append(steps_List, round(uncertainty_curve[[peptide_Column - 1]]$steps, 2))
          probs_List <- append(probs_List, uncertainty_curve[[peptide_Column - 1]]$prob_List_Normalized)
          increment <- increment + uncertainty_curve[[peptide_Column - 1]]$h
          
          legend <- paste(legend, as.character(prot_Col[i]))
        } 
      } 
    }

    DT <- data.table(steps = unlist(steps_List), probs = unlist(probs_List))[order(-steps, decreasing = TRUE)]
    
    if(aggregation_Mode == 1){
      # sum up over the same steps values
      DT <- DT[, lapply(.SD, sum), by = steps]
    }
    else if(aggregation_Mode == 2){
      # mult up over the same steps values
      DT <- DT[, lapply(.SD, prod), by = steps]
    }
    else{
      # Error!
    }
    
    probs_List <- unlist(DT$probs) 
    #Normaliaztion
    #change the sum of the probs to 1
    # sum_prob <- sum (unlist (prob_List)) = 1
    prob_List_Normalization_Factor <- 1 / sum (probs_List) # sum should be 1
    prob_List_Normalized0 <- lapply (probs_List, function (x) x * prob_List_Normalization_Factor)
    
    prob_List_Normalized <- lapply (prob_List_Normalized0, function (x) x / increment)
    
    #save in list for Aggregation
    x <- unlist(prob_List_Normalized)
    # x <- unlist(prob_List_Normalized0)
    
    uncertainty_aggragation[[1]] <- list(steps = DT$steps, probs = x)

    #title for plot
    if (title == ""){
    
      title <- paste("Intensitiy Aggregation", Peptide_Name, sep = "-")
    }
  }

  
  # ------------- Quotient
  # -------------
  if(aggregation_Type == 2){ # aggregation for quotients
    name_Of_Returned_Object <- paste("quotient_aggragation", Peptide_Name, sep = "_")

    #Row from the selected Peptide
    peptide_Rows <- which(input_table[[1]] == Peptide_Name)
    peptide_Row <- peptide_Rows[1] # if multiple rows with same name pick the first one (should not appen!)

    name <- paste("uncertainty_quotient", Peptide_Name, sep = "_")

    uncertainty_quotient <- NULL
    
    tryCatch({
        message("aggregation_of_uncertainties: loading ", name, " (aggregation_Type == 2)")
      
        uncertainty_quotient <- eval(as.symbol(name))
      },
      
      error = function(name){
        
        message("aggregation_of_uncertainties: Error while loading ", name, " (aggregation_Type == 2)")
      },
      
      warning = function(name){

        message("aggregation_of_uncertainties: Warning while loading ", name, " (aggregation_Type == 2)")
      },
      
      finally = {
        if(is.null(uncertainty_quotient)){
          message("aggregation_of_uncertainties: calculating ", name, " (aggregation_Type == 2)")
          
          uncertainty_quot(input_table = input_table,
                           Peptide_Name = Peptide_Name,
                           prot_Col_First = prot_Col_First,
                           prot_Col_Second = prot_Col_Second,
                           defaul_Col_First = prot_Col_First,
                           defaul_Col_Second = prot_Col_Second,
                           general_Sd = NA,
                           aggregation_Type = 2,
                           plotQ = F)
          
          message("aggregation_of_uncertainties: reloading ", name, " (aggregation_Type == 2)")
          uncertainty_quotient <- eval(as.symbol(name))
        }
      }
    )
    
    increment <- 0

    for(i in seq(1, length(prot_Col_First), 1)){
      
      if(!is.null(uncertainty_quotient[[i]][1])){ # uncertainty_quotient exist?
        # -1 to skip the Sequence index of Peptids
        if(!is.na(uncertainty_quotient[[i]][1])){
          
          steps_List <- append(steps_List, round(uncertainty_quotient[[i]]$steps, 2))
          probs_List <- append(probs_List, uncertainty_quotient[[i]]$prob_List_Normalized)
          increment <- increment + uncertainty_quotient[[i]]$h
          
          legend <- paste(legend, paste(as.character(prot_Col_Second[i]), as.character(prot_Col_First[i]), sep = "/")) # F5..Area/B5..Area
        }
      }
    }

    DT <- data.table(steps = unlist(steps_List), probs = unlist(probs_List))[order(-steps, decreasing = TRUE)]

    if(aggregation_Mode == 1){
      # sum up over the same steps values
      DT <- DT[, lapply(.SD, sum), by = steps]
    }
    else if(aggregation_Mode == 2){
      # mult up over the same steps values
      DT <- DT[, lapply(.SD, prod), by = steps]
    }
    else{
      # Error!
    }

    probs_List <- unlist(DT$probs)
    #Normaliaztion
    #change the sum of the probs to 1
    prob_List_Normalization_Factor <- 1 / sum (probs_List) # sum should be 1
    prob_List_Normalized0 <- lapply (probs_List, function (x) x * prob_List_Normalization_Factor)

    prob_List_Normalized <- lapply (prob_List_Normalized0, function (x) x / increment)

    #save in list for Aggregation
    x <- unlist(prob_List_Normalized)
    # x <- unlist(prob_List_Normalized0)

    uncertainty_aggragation[[1]] <- list(steps = DT$steps, probs = x)

    #title for plot
    if (title == ""){

      title <- paste("Quotient Aggregation", Peptide_Name, sep = "-")
    }
  }
  

  if(plotA){
    ploti <- 0
    #number of non NA's for colors
    color_count <- length(legend)
    #corlornumber from colorcount
    color <- 0

    for (i in seq(1, length(uncertainty_aggragation), 1)) {
      if(is.na(uncertainty_aggragation[[i]][1])){

        next # skip uncertainty_aggragation
      }

      ploti <- ploti + 1
      color <- color + 1

      if(ploti == 1){ # first time create new plot
        if(is.null(xmax)){
          xmax <- max(uncertainty_aggragation[[i]]$steps)
        }
        if(is.null(xmin)){
          xmin <- min(uncertainty_aggragation[[i]]$steps)
        }
        if(is.null(ymax)){
          ymax <- max(uncertainty_aggragation[[i]]$probs)
        }
          
        xmin <- min(xmin, min(uncertainty_aggragation[[i]]$steps))
        xmax <- max(xmax, max(uncertainty_aggragation[[i]]$steps))
        ymax <- max(ymax, max(uncertainty_aggragation[[i]]$probs))
        
        plot (uncertainty_aggragation[[i]]$steps, uncertainty_aggragation[[i]]$probs, type = "l", main = title, xlim = c(xmin, xmax), ylim = c(0, ymax), xlab = xLabel, ylab = yLabel, col = rich.colors(color_count)[color])
        shadowtext(x = (xmax + xmin)/2, y = 0, legend, col='lightsteelblue')
      }
      if (ploti > 1){
        lines (uncertainty_aggragation$steps, uncertainty_aggragation$probs, type = "l", col = rich.colors(color_count)[color])
      }
    }

    assign (name_Of_Returned_Object, uncertainty_aggragation, envir = .GlobalEnv)
    
    # legend("topright", legend = legend, col = rich.colors(color_count), pch = 20, cex = 0.7)  
  }
}

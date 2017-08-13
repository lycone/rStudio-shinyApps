      #' Plotting of the uncertainty of dataset values
      #'
      #' The function can plot the uncertainty of the chosen dataset values in two possible ways. First the superior element like aggregation and second the row like aggregation.
      #' The superior element like aggregation plots all values from the chosen columns related to a given superior element. Whereas the row like aggregation plots all
      #' values from the chosen columns related to a given subordinate element in a specific row. Furthermore a list is returned containing the discretized values of the uncertainty curves. This list is used for the uncertainty aggregation.
      #' @param input_table Dataset as a dataframe which should be used for the computations.
      #' @param proteinName Name of the superior element from which the values should be plotted (superior element like aggregation).
      #' @param rowName Name of the subordinate element from which the values should be plotted (Row like aggregation).
      #' @param name_Of_SupEl_Col Column name of the superior element column (superior element like aggregation).
      #' @param col_Legend Column name of the column for legend / column with subordinate elements.
      #' @param prot_Col Column names which contains the values that should be plotted (superior element like aggregation).
      #' @param prot_Col2 Column names which contain the values that should be plotted (Row like aggregation).
      #' @param title Optional character for the headline of the plot.
      #' @param general_Sd Optional numeric value which defines a general standard deviation for the uncertainty evaluation.
      #' @param title Optional character for the headline of the plot.
      #' @param ymax Maximal range of y axis.
      #' @param xlabel Label of the x axis.
      #' @param ylabel Label of the y axis.
      #' @param aggregation_Type 1 = superior element like aggregation, 2 = Row like aggregation.
      #' @param plotF Optional logical value. TRUE = plot returned, FALSE = no plot returned  (default = TRUE).
      #' @examples superior element like aggregation:
      #' uncertainty_eval (input_table, proteinName = "Protein_X", name_Of_SupEl_Col = "Protein",
      #' col_Legend = "Peptide", prot_Col = c("3min", "6min"), aggregation_Type = 1)
      #'
      #'row like aggregation:
      #' likelihood_curves (input_table, rowName = "Peptide_X",col_Legend = "Peptide",
      #'                    prot_Col2 = c("3min", "6min"), aggregation_Type = 2)
      #'
      #'If rowNames of the col_Legend are not unique, you have to add an additional
      #'colum to the input_table with the row numbers:
      #' likelihood_curves (input_table, rowName = "7",col_Legend = "Row_Numbers",
      #'                    prot_Col2 = c("3min", "6min"), aggregation_Type = 2)
      #'
      #' @export

uncertainty_quot <- function (input_table,
            Peptide_Name = NA,
            Peptide_Name_List = NA, # Peptids list
            prot_Col_First = NA, # 1st list of selected columns
            prot_Col_Second = NA, # 2nd list of selected columns
            defaul_Col_First = NA, # 1st list default
            defaul_Col_Second = NA, # 2nd list default
            general_Sd = NA,
            title = "",
            ymax = 1, # max probability doesn't surpace 1
            xLabel = "x",
            yLabel = "y",
            aggregation_Type,
            plotQ = T){

  source("R/uncertainty_evaluation.R")
  source ("R/quotient_utils.R")

  uncertainty_quotient <- list()
  name_Of_Returned_Object <- c()

  legend <- c()
  
  #colum numbers from choosen replicates
  peptide_Column_First <- c()
  peptide_Column_Second <- c()
  
  xmin <- NULL
  xmax <- NULL
  ymax <- NULL

  if(aggregation_Type == 2){
    name_Of_Returned_Object <- paste("uncertainty_quotient", Peptide_Name, sep = "_")
    
    #Row from the selected Peptide
    peptide_Rows <- which(input_table[[1]] == Peptide_Name)
    peptide_Row <- peptide_Rows[1] # if multiple rows with same name pick the first one (should not appen!)

    name <- paste("uncertainty_curve", Peptide_Name, sep = "_")
    
    uncertainty_curve <- NULL
    
    tryCatch({
      message("uncertainty_quot: loading ", name, " (aggregation_Type == 2)")
      
      uncertainty_curve <- eval(as.symbol(name))
      },
      
      error = function(name){
        
        message("uncertainty_quot: Error while loading", name, " (aggregation_Type == 2)")
      },
      
      warning = function(name){

        message("uncertainty_quot: Warning while loading ", name, " (aggregation_Type == 2)")
      },
      
      finally = {
        if(is.null(uncertainty_curve)){
          message("uncertainty_quot: calculating ", name, " (aggregation_Type == 2)")
          uncertainty_eval( input_table = input_table,
                            Peptide_Name = Peptide_Name,
                            col_Legend = names(input_table[1,1]),
                            Peptide_Name_List = NA,
                            prot_Col = c(defaul_Col_First, defaul_Col_Second),
                            general_Sd = NA,
                            aggregation_Type = 2,
                            plotF = F)
          
          message("uncertainty_quot: reloading ", name, " (aggregation_Type == 2)")
          uncertainty_curve <- eval(as.symbol(name))
        }
      }
    )
  
    for (i in seq(1, length(prot_Col_First), 1)) {
      peptide_Column_First <- which(colnames(input_table) == prot_Col_First[i])
      peptide_Column_Second <- which(colnames(input_table) == prot_Col_Second[i])
      
      if(!is.null(uncertainty_curve[[i]][1])){ # uncertainty_curve exist?
        # -1 to skip the Sequence index of Peptids
        if(!is.na(uncertainty_curve[[peptide_Column_First - 1]][1]) && 
          !is.na(uncertainty_curve[[peptide_Column_Second - 1]][1])){
          
          uncertainty_quotient[[i]] <- quotient_utils(uncertainty_curve[[peptide_Column_First - 1]],
                                                      uncertainty_curve[[peptide_Column_Second - 1]])
          
          legend <- append(legend, paste(as.character(prot_Col_Second[i]), as.character(prot_Col_First[i]), sep = "/")) # F5..Area/B5..Area
          
          if(is.null(xmax)){
            xmax <- max(uncertainty_quotient[[i]]$steps)
          }
          if(is.null(xmin)){
            xmin <- min(uncertainty_quotient[[i]]$steps)
          }
          if(is.null(ymax)){
            ymax <- max(uncertainty_quotient[[i]]$prob_List_Normalized)
          }
          
          xmin <- min(xmin, min(uncertainty_quotient[[i]]$steps))
          xmax <- max(xmax, max(uncertainty_quotient[[i]]$steps))
          ymax <- max(ymax, max(uncertainty_quotient[[i]]$prob_List_Normalized))
        } else{ # uncertainty_curve does not exist
  
          uncertainty_quotient[[i]] <- NA
        }
      } 
    }

    #title for plot
    if (title == ""){
    
      title <- Peptide_Name
    }
    
    # only save object from aggregation_Type == 2 (complete objects)
    assign (name_Of_Returned_Object, uncertainty_quotient, envir = .GlobalEnv)
  }

  if(aggregation_Type == 1){

    i <- 1

    for (Peptide_Name in Peptide_Name_List) {
      name_Of_Returned_Object <- paste("uncertainty_quotient", Peptide_Name, sep = "_")
      #Row from the selected Peptide
      peptide_Rows <- which(input_table[[1]] == Peptide_Name)
      peptide_Row <- peptide_Rows[1] # if multiple rows with same name pick the first one (should not appen!)


      name <- paste("uncertainty_curve", Peptide_Name, sep = "_")
      # uncertainty_curve <- eval(as.symbol(name))
      uncertainty_curve <- NULL
      
      tryCatch({
        message("uncertainty_quot: loading ", name, " (aggregation_Type == 1)")
        
        uncertainty_curve <- eval(as.symbol(name))
      },
      
      error = function(name){
        
        message("uncertainty_quot: Error while loading ", name, " (aggregation_Type == 1)")
      },
      
      warning = function(name){
        
        message("uncertainty_quot: Warning while loading ", name, " (aggregation_Type == 1)")
      },
      
      finally = {
        if(is.null(uncertainty_curve)){
          message("uncertainty_quot: calculating ", name, " (aggregation_Type == 1)")
          uncertainty_eval( input_table = input_table,
                            Peptide_Name = Peptide_Name,
                            col_Legend = names(input_table[1,1]),
                            Peptide_Name_List = NA,
                            prot_Col = c(defaul_Col_First, defaul_Col_Second),
                            general_Sd = NA,
                            aggregation_Type = 2,
                            plotF = F)
          
          message("uncertainty_quot: reloading ", name, " (aggregation_Type == 1)")
          uncertainty_curve <- eval(as.symbol(name))
        }
      }
      )
      
      # for (k in seq(1, length(prot_Col_First), 1)){
        peptide_Column_First <- which(colnames(input_table) == prot_Col_First[1])
        peptide_Column_Second <- which(colnames(input_table) == prot_Col_Second[1])
        
        if(!is.null(uncertainty_curve[[1]][1])){ # uncertainty_curve exist?
          # -1 to skip the Sequence index of Peptids
          if(!is.na(uncertainty_curve[[peptide_Column_First - 1]]) && 
            !is.na(uncertainty_curve[[peptide_Column_Second - 1]])){
            
            uncertainty_quotient[[i]] <- quotient_utils(uncertainty_curve[[peptide_Column_First - 1]],
                                                        uncertainty_curve[[peptide_Column_Second - 1]])
            
            legend <- append(legend, as.character(Peptide_Name))
            
            if(is.null(xmax)){
              xmax <- max(uncertainty_quotient[[i]]$steps)
            }
            if(is.null(xmin)){
              xmin <- min(uncertainty_quotient[[i]]$steps)
            }
            if(is.null(ymax)){
              ymax <- max(uncertainty_quotient[[i]]$prob_List_Normalized)
            }
            
            xmin <- min(xmin, min(uncertainty_quotient[[i]]$steps))
            xmax <- max(xmax, max(uncertainty_quotient[[i]]$steps))
            ymax <- max(ymax, max(uncertainty_quotient[[i]]$prob_List_Normalized))
          } else{ # uncertainty_curve does not exist
    
            uncertainty_quotient[[i]] <- NA
          }
        } 
      # }
      
      i <- i + 1
    }
    #title for plot
    if (title == ""){
    
      title <- paste(prot_Col_Second, prot_Col_First, sep = "/")
    }
  }

  if(plotQ){
    ploti <- 0
    #number of non NA's for colors
    color_count <- length(legend)
    #corlornumber from colorcount
    color <- 0

    for (i in seq(1, length(uncertainty_quotient), 1)) {
      if(is.na(uncertainty_quotient[[i]][1])){

        next # skip uncertainty_quotient
      }

      ploti <- ploti + 1
      color <- color + 1

      if(ploti == 1){ # first time create new plot
        plot (uncertainty_quotient[[i]]$steps, uncertainty_quotient[[i]]$prob_List_Normalized, type = "l", main = title, xlim = c(xmin, xmax), ylim = c(0, ymax), xlab = xLabel, ylab = yLabel, col = rich.colors(color_count)[color])
      }
      if (ploti > 1){
        lines (uncertainty_quotient[[i]]$steps, uncertainty_quotient[[i]]$prob_List_Normalized, type = "l", col = rich.colors(color_count)[color])
      }
    }
    
    legend("topright", legend = legend, col = rich.colors(color_count), pch = 20, cex = 0.7)  
  }
}
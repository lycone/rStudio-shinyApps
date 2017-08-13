# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
# Autor: Boris Nguema Bekale

library(shiny)

library("XLConnectJars")
library("XLConnect")
library(data.table)
library(gplots)
library(TeachingDemos)

library(stats)

################ Wrapper function for datapath
fixUploadedFilename <- function(fileX){
  dir = fileX$datapath
  name = file.path(dirname(fileX$datapath), fileX$name)
  file.rename(from = dir, to = name)
  fileX$datapath <- name
  
  fileX
}


############### Server
############### 
############### 
############### 

shinyServer(function(input, output, session) {
  
  ############### database
  getElementData <- reactive({
    validate(need(!is.null(input$fileInput), "please select file"))
    
    filename <- fixUploadedFilename(input$fileInput)
    
    datEx <- NULL
    
    if(input$fileInputTyp == 1){ # Excel
      datEx <- readWorksheet(loadWorkbook(filename$datapath), 1)
    }else if(input$fileInputTyp == 2){ # text/csv
      datEx <- read.csv2(filename$datapath)
    }else{
      #obviously an error has occurred
    }
    
    if(is.null(datEx)) return(NULL)
    
    element_data <<- data.table(datEx)
    
    return(element_data)
  })
  
  # on second element_data calculated
  observeEvent(getElementData(), {
    if(!is.null(element_data)){
      ############################## Database  ##############################
      updateSelectInput(session, "Columns_For_Discr", choices = colnames(element_data[1, -1]), selected = F)
      updateSelectInput(session, "Columns_For_Discr2", choices = colnames(element_data[1, -1]), selected = F)
      
      ############################## Frenquency Distribution ##############################
      updateSelectInput(session, "element_Columns_Freq_Dist", choices = colnames(element_data[1, -1]), selected = colnames(element_data[1, 1]))
    }
  })
  
  ############### output element dataset
  output$dataSet <- renderDataTable({
    return(getElementData()[, input$dataSetStartColumn:input$dataSetEndColumn])
  })
  

    
  ############### quotient dataset
  getQuotientData <- reactive({
    quotient_data <<- NULL
    
    if(!is.null(input$Columns_For_Discr) && !is.null(input$Columns_For_Discr2) && 
       length(input$Columns_For_Discr) == length(input$Columns_For_Discr2) &&
       length(input$Columns_For_Discr) > 1 &&
       !is.null(input$RB_Replace)){
      
      replace = function(value, replaceBy){
        
        if(value < 1){
          if(replaceBy == TRUE)
            value <- 1
        }
        
        return(value)
      }
      
      sequence <- colnames(element_data[1, 1])
      quotient_data <<- element_data[, lapply(.SD, function(x) log2(x)),  by = sequence, .SDcols = c(input$Columns_For_Discr, input$Columns_For_Discr2)]
      # quotient_data <<- element_data[, lapply(.SD, function(x) replace(x, input$RB_Replace)),  by = Sequence, .SDcols = c(input$Columns_For_Discr, input$Columns_For_Discr2)]
    }

    return(quotient_data)
  })
  
  
  # on second quotient_data calculated
  observeEvent(getQuotientData(), {
    if(!is.null(quotient_data)){
      ############################## Uncertainty  Computation ##############################
      updateSelectInput(session, "element_Column_For_Discr", choices = c(input$Columns_For_Discr, input$Columns_For_Discr2), selected = input$Columns_For_Discr[1])
      updateSelectInput(session, "element_Rows_For_Discr", choices = quotient_data[ , 1], selected = quotient_data[1, 1])
      updateSelectInput(session, "element_Columns_For_Discr", choices = c(input$Columns_For_Discr, input$Columns_For_Discr2), selected = c(input$Columns_For_Discr, input$Columns_For_Discr2))
      updateSelectInput(session, "element_Row_For_Disc", choices = quotient_data[ , 1], selected = quotient_data[1, 1])

      updateSelectInput(session, "first_Element_Column_Selection", choices = input$Columns_For_Discr, selected = input$Columns_For_Discr[1])
      updateSelectInput(session, "second_Element_Column_Selection", choices = input$Columns_For_Discr2, selected = input$Columns_For_Discr2[1])
      updateSelectInput(session, "quotient_Rows_For_Uncertainty", choices = quotient_data[ , 1], selected = quotient_data[1, 1])
      
      updateSelectInput(session, "first_Element_Columns_Selection", choices = input$Columns_For_Discr, selected = input$Columns_For_Discr)
      updateSelectInput(session, "second_Element_Columns_Selection", choices = input$Columns_For_Discr2, selected = input$Columns_For_Discr2)
      updateSelectInput(session, "quotient_Row_For_Uncertainty", choices = quotient_data[ , 1], selected = quotient_data[1, 1])
      
      updateSelectInput(session, "element_Row_For_Aggregation", choices = quotient_data[ , 1], selected = quotient_data[1, 1])
      updateSelectInput(session, "quotient_Row_For_Aggregation", choices = quotient_data[ , 1], selected = quotient_data[1, 1])
    }
  })
  
  # Render of DataTable of QuotientData
  output$quotient <- renderDataTable({
    getQuotientData()
  })
  
  
  ############### sd computation
  getDynamicSD <- reactive({
    source("R/dynamic_standard_deviation.R")
    source("R/main.R")
    
    if(!is.null(quotient_data) &&
       !is.null(input$Columns_For_Discr) &&
       !is.null(input$Columns_For_Discr2) &&
       !is.null(input$column_Rows_For_Discr_Default)){
       # !is.null(input$column_Rows_For_Discr_Custom_min) &&
       # !is.null(input$column_Rows_For_Discr_Custom_max)){
      
      Columns_For_Discr <- input$Columns_For_Discr
      Columns_For_Discr2 <- input$Columns_For_Discr2
      
      general_Sd <- dynamic_standard_deviation(input_table = quotient_data,
                                 n = 100,
                                 columns_For_Sd = Columns_For_Discr,
                                 columns_For_Sd2 = Columns_For_Discr2,
                                 use.mad = T)
      
      # compute main for all data      
      if(input$column_Rows_For_Discr_Default == 1){
        start_Row <- 1
        end_Row <- nrow(quotient_data)
      }else if(input$column_Rows_For_Discr_Default == 2){
        
        start_Row <- input$column_Rows_For_Discr_Custom_min
        end_Row <- input$column_Rows_For_Discr_Custom_max
      }
      
      withProgress(message = "Data discretisation",
                   detail = "This may take a while...", value = 0, {
                  for (column_Name in c(Columns_For_Discr, Columns_For_Discr2)){
                    incProgress(1/length(c(Columns_For_Discr,Columns_For_Discr2)))
                      
                    main(input_table = quotient_data,
                            column_For_Discr = column_Name,
                            start_Row = start_Row,
                            end_Row = end_Row,
                            number_Of_Discretizations = 100)
                  }
      })
    }
    
    return(general_Sd)
  })
  
  
  output$element_SD_plot <- renderPlot({
    getDynamicSD()    
  })


  
  ############################## uncertainty computation  ##############################
  
  ############################## element
  
  # show row wise uncertainty Plots
  output$element_row_wise_plot <- renderPlot({
    source("R/uncertainty_evaluation.R")
    
    if(getDynamicSD() &&
       !is.null(input$element_Rows_For_Discr) &&
       !is.null(input$element_Column_For_Discr)){
      
      uncertainty_eval(input_table = quotient_data,
                       Peptide_Name_List = input$element_Rows_For_Discr,
                       col_Legend = names(quotient_data[1,1]),
                       # col_Legend = "Sequence",
                       prot_Col = input$element_Column_For_Discr,
                       title = input$element_Column_For_Discr,
                       xLabel = "Intensity Value (log2(x))",
                       yLabel = "Likelihood",
                       aggregation_Type = 1,
                       plotF = T)
    }
  })
  
  
  # show column wise uncertainty Plots
  output$element_column_wise_plot <- renderPlot({
    source("R/uncertainty_evaluation.R")
    
    if(!is.null(quotient_data) && 
       !is.null(input$element_Row_For_Disc) &&
       !is.null(input$element_Columns_For_Discr)){

      uncertainty_eval(input_table = quotient_data,
                      Peptide_Name = input$element_Row_For_Disc, # Proteinname/Peptidename
                      # col_Legend = "Sequence",
                      col_Legend = colnames(quotient_data[1, 1]),
                      prot_Col = input$element_Columns_For_Discr,
                      general_Sd = NA,
                      title = "",
                      xLabel = "Intensity Value (log2(x))",
                      yLabel = "Likelihood",
                      aggregation_Type = 2,
                      plotF = T)
    }
  })
  
  ############################## quotient
  
  # show row wise uncertainty Plots
  output$quotient_row_wise_plot <- renderPlot({
    source("R/uncertainty_quotient.R")
    
    if(!is.null(quotient_data) && 
          !is.null(input$first_Element_Column_Selection) &&
          !is.null(input$second_Element_Column_Selection) &&
          !is.null(input$quotient_Rows_For_Uncertainty)){
      
      uncertainty_quot(input_table = quotient_data,
                       Peptide_Name_List = input$quotient_Rows_For_Uncertainty,
                       prot_Col_First = input$first_Element_Column_Selection,
                       prot_Col_Second = input$second_Element_Column_Selection,
                       defaul_Col_First = input$Columns_For_Discr,
                       defaul_Col_Second = input$Columns_For_Discr2,
                       general_Sd = NA,
                       title = "",
                       ymax = NA,
                       xLabel = "Intensity Value (log2(x))",
                       yLabel = "Likelihood",
                       aggregation_Type = 1,
                       plotQ = T)
    }
  })
  
  
  # show column wise uncertainty Plots
  output$quotient_column_wise_plot <- renderPlot({
    source("R/uncertainty_quotient.R")
    
    if(!is.null(getQuotientData()) && 
       !is.null(input$first_Element_Columns_Selection) &&
       !is.null(input$second_Element_Columns_Selection) &&
       length(input$first_Element_Columns_Selection) == length(input$second_Element_Columns_Selection) &&
       !is.null(input$quotient_Row_For_Uncertainty)){
      
      uncertainty_quot(input_table = quotient_data,
                      Peptide_Name = input$quotient_Row_For_Uncertainty,
                      prot_Col_First = input$first_Element_Columns_Selection,
                      prot_Col_Second = input$second_Element_Columns_Selection,
                      defaul_Col_First = input$Columns_For_Discr,
                      defaul_Col_Second = input$Columns_For_Discr2,
                      general_Sd = NA,
                      title = "",
                      ymax = NA,
                      xLabel = "Regulation Factor (log2(x))",
                      yLabel = "Likelihood",
                      aggregation_Type = 2,
                      plotQ = T)
    }
  })
  
  
  ############################## aggregation computation  ##############################
  
  # ############################## element
  
  output$element_aggregation_row_wise_plot <- renderPlot({
    source("R/aggregation_of_uncertainties.R")
    
    if(!is.null(getQuotientData()) &&
       !is.null(input$element_Row_For_Aggregation) &&
       !is.null(input$element_Columns_For_Discr) &&
       !is.null(input$element_Aggregation_Mode)){ 
      
      aggregation_of_uncertainties(input_table = quotient_data,
                                  Peptide_Name = input$element_Row_For_Aggregation,
                                  prot_Col = input$element_Columns_For_Discr,
                                  prot_Col_First = NA,
                                  prot_Col_Second = NA,
                                  title = "",
                                  xLabel = "Intensity Value (log2(x))",
                                  yLabel = "Likelihood",
                                  aggregation_Type = 1, # 1 = elements or 2 = quotients
                                  aggregation_Mode = input$element_Aggregation_Mode, # Additiv - Multiplicatikv
                                  plotA = T)
    }
  })
  
  
  # ############################## quotient
  # 
  # show row wise aggregation Plots
  output$quotient_aggregation_row_wise_plot <- renderPlot({
    source("R/aggregation_of_uncertainties.R")

    if(!is.null(getElementData()) &&
       !is.null(input$quotient_Row_For_Aggregation) &&
       !is.null(input$first_Element_Columns_Selection) &&
       !is.null(input$second_Element_Columns_Selection) &&
       !is.null(input$quotient_Aggregation_Mode)){

      aggregation_of_uncertainties(input_table = quotient_data,
                                   Peptide_Name = input$quotient_Row_For_Aggregation,
                                   prot_Col = NA,
                                   prot_Col_First = input$Columns_For_Discr,
                                   prot_Col_Second = input$Columns_For_Discr2,
                                   title = "",
                                   xLabel = "Regulation Factor (log2(x))",
                                   yLabel = "Likelihood",
                                   aggregation_Type = 2, # 1 = elements or 2 = quotients
                                   aggregation_Mode = input$quotient_Aggregation_Mode, # Additiv - Multiplicatikv
                                   plotA = T)
    }
  })
  
  # clean up the enviroment
  session$onSessionEnded(function(){
    rm(list = ls())
  })
})
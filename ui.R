
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
# Autor: Boris Nguema Bekale

library(shiny)

shinyUI(
  navbarPage(id = "NavbarPageId", "Uncertainty Modeling for Omics Data",
             # Tab: Database
             tabPanel("Database",   titlePanel("Data Configuration"),
                      sidebarLayout(
                        sidebarPanel(
                          radioButtons("fileInputTyp", h5(strong("1. Set Datatyp")), c("Excel" = 1, "csv" = 2), selected = 1, inline = TRUE),
                          fileInput("fileInput", "1. Load Data File", multiple = FALSE,
                                    accept=c('application/vnd.ms-excel','application/vnd.openxmlformats-officedocument.spreadsheetml.sheet','.xls','.xlsx', # For Excel-Files
                                             "text/csv","text/comma-separated-values,text/plain",".csv") # For CSV-Files
                          ),
                          hr(),
                          
                          h5(strong("2. Data Presentation")),
                          # checkboxInput("showRowNumber", "Show Row Number", value = T),
                          # checkboxInput("skipEmptyResult", "Skip Empty Hits", value = T),
                          # radioButtons("RB_Replace", "Replace Empty Hits With:",
                          #              c("1" = 1,
                          #                "NA" = NA,
                          #                "LoD" = "LoD"), selected = NA, inline = TRUE),
                          checkboxInput("RB_Replace"
                                        , "Replace Empty Hits With [1]", value = TRUE, width = NULL),
                          
                          numericInput("dataSetStartColumn", "Start Column", min = 0, value = 1, step = 1),
                          numericInput("dataSetEndColumn", "End Column", min = 3, value = 5, step = 2),
                          hr(),
                          
                          h5(strong("3.Setting of Replicates for Computation")),
                          selectInput("Columns_For_Discr", "Select 1st. Replicates (Columns)", choices = "", multiple = T, 
                                      helpText("select 1st. columns for futher calculations")),
                          selectInput("Columns_For_Discr2", "Select 2nd. Replicates (Columns)", choices = "", multiple = T, 
                                      helpText("select 2nd. columns for futher calculations")),
                          hr(),
                          
                          radioButtons("column_Rows_For_Discr_Default", h5(strong("4. Set Number of Rows for Computation")), c("All" = 1, "Custom" = 2), selected = 2),
                          conditionalPanel(
                            condition = "input.column_Rows_For_Discr_Default == 2",
                            numericInput("column_Rows_For_Discr_Custom_min", "Start Row-Index", value = 1, min = 1, step = 100),
                            numericInput("column_Rows_For_Discr_Custom_max", "End Row-Index", value = 100, min = 2, step = 100))
                        ),
                        
                        # Show data
                        mainPanel(
                          tags$style(type="text/css",".shiny-output-error {visibility: hidden;}",
                                     ".shiny-output-error:before {visibility: hidden;}"),
                          tabsetPanel(id ="dataSet", type = "tab", #position = "above",
                                      # element factors overview
                                      tabPanel("Element",
                                               h3(strong("Element Factors Overview")),
                                               dataTableOutput("dataSet")
                                      ),
                                      
                                      # quotient factors overview
                                      tabPanel("Quotient",
                                               h3(strong("Quotient Factors Overview")),
                                               dataTableOutput("quotient")
                                      )
                          )
                          )
                        )
                      ),
             
          # SD Computation
          tabPanel("SD Computation",   titlePanel("SD Settings"),
                   sidebarLayout(
                     sidebarPanel(
                       h5(strong("Regulation-Factors")),
                       
                       checkboxInput("dataSetSD", "Use SD from Dataset", value = T),
                       # br(),
                       numericInput("selectedSD", "SD for Uncertainty of Elements", min = 0, value = 1, step = 10),
                       
                       br(),
                       
                       tags$hr(),
                       h5(strong("Quotient-Factors")),
                       
                       checkboxInput("quotientDataSD", "Use SD from Dataset", value = T),
                       # br(),
                       numericInput("selectedQuotientDataSD", "SD for Uncertainty of Quotients", min = 0, value = 1, step = 10)
                     ),
                     
                     # Show data
                     # aggregation curves
                     mainPanel(tags$style(type="text/css",".shiny-output-error {visibility: hidden;}",
                                          ".shiny-output-error:before {visibility: hidden;}"),
                               tabPanel("element_SD_plot",
                                        # h3(strong("Normal Distributions of selected Elements")),
                                        h3(strong("SD Distribution of Elements")),
                                        plotOutput("element_SD_plot")
                               )
                     )
                   )
          ),
          

          # Uncertainty Computation
          navbarMenu(title = "Uncertainty Computation",
             # Uncertainty of Elements
             tabPanel("uncertainty of elements", titlePanel("Uncertainty of Elements"), value = "UncertaintyElement",
                      sidebarLayout(
                        sidebarPanel(
                          tabsetPanel(
                            tabPanel("Column_Wise",
                                      br(),
                                      selectInput("element_Columns_For_Discr", h4(strong("1. Select Column(s) for Uncertainty Computation")), choices = "", multiple = T),
                                               
                                      hr(),
                                      selectInput("element_Row_For_Disc", h4(strong("2. Select an Element for Uncertainty Computation")), choices = "", multiple = F)
                                    ),
                            
                            tabPanel("Row_Wise",
                                      br(),
                                      selectInput("element_Column_For_Discr", h4(strong("1. Select a Column for Uncertainty Computation")), choices = "", multiple = F),
                                     
                                      hr(),
                                      selectInput("element_Rows_For_Discr", h4(strong("2. Select Element(s) for Uncertainty Computation")), choices = "", multiple = T)
                                    )
                          )
                        ),
                        
                        # uncertainty curves
                        mainPanel(tags$style(type="text/css",".shiny-output-error {visibility: hidden;}",
                                             ".shiny-output-error:before {visibility: hidden;}"),
                                  tabsetPanel(id ="UncertaintyElement", type = "tab", #position = "above",
                                              # Uncertainty column wise
                                              tabPanel("Column_Wise",
                                                       h3(strong("Normal Distributions of selected Elements")),
                                                       plotOutput("element_column_wise_plot")
                                              ),
                                              
                                              # Uncertainty row wise
                                              tabPanel("Row_Wise",
                                                       h3(strong("Normal Distributions of selected Element")),
                                                       plotOutput("element_row_wise_plot")
                                              )
                                  )
                        )
                      )
             ),
             
             
             # Uncertainty of Quotients
             tabPanel("uncertainty of quotients", titlePanel("Uncertainty of Quotients"), value = "UncertaintyQuotient",
              sidebarLayout(
                sidebarPanel(
                  tabsetPanel(
                    tabPanel("Column_Wise",
                             br(),
                             
                             h4(strong("1. Setting of Column-Factors for Uncertainty Computation")),
                             selectInput("first_Element_Columns_Selection", "Select First Factors", choices = "", multiple = T, 
                                         helpText("First parameter for the calculation of Regulation-Factor")),
                             selectInput("second_Element_Columns_Selection", "Select Second Factors", choices = "", multiple = T,
                                         helpText("Second parameter for the calculation of Regulation-Factor")),
                             
                             br(),
                             hr(),
                             selectInput("quotient_Row_For_Uncertainty", h4(strong("2. Select an Element for Uncertainty Computation")), choices = "", multiple = F)
                          ),
                    
                    tabPanel("Row_Wise",
                             br(),
                             
                             h4(strong("1. Setting of Column-Factors for Uncertainty Computation")),
                             selectInput("first_Element_Column_Selection", "Select First Factor", choices = "", multiple = F, 
                                         helpText("First parameter for the calculation of Regulation-Factor")),
                             selectInput("second_Element_Column_Selection", "Select Second Factor", choices = "", multiple = F,
                                         helpText("Second parameter for the calculation of Regulation-Factor")),
                             
                             br(),
                             hr(),
                             selectInput("quotient_Rows_For_Uncertainty", h4(strong("2. Select Element(s) for Uncertainty Computation")), choices = "", multiple = T)
                    )
                  )
              ),
                        
            # uncertainty curves
            mainPanel(tags$style(type="text/css",".shiny-output-error {visibility: hidden;}",
                      ".shiny-output-error:before {visibility: hidden;}"),
              tabsetPanel(id ="UncertaintyQuotient", type = "tab", #position = "above",
                # uncertainty column wise
                tabPanel("Column_Wise",
                          h3(strong("Normal Distributions of selected Quotients")),
                          plotOutput("quotient_column_wise_plot")
                ),
                
                # uncertainty row wise
                tabPanel("Row_Wise",
                         h3(strong("Normal Distributions of selected Quotient")),
                         plotOutput("quotient_row_wise_plot")
                )
              )
            )
          )
        )
      ),
     


    # Aggregation Computation
    navbarMenu(title = "Aggregation Computation",
      # Aggregation of Elements
      tabPanel("aggregation of elements", titlePanel("Aggregation of Elements"), value = "AggregationElement",
        sidebarLayout(
          sidebarPanel(
            tabPanel("Row_Wise",
                br(),
                selectInput("element_Row_For_Aggregation", h4(strong("1. Select Element(s) for Aggregation Computation")), choices = "", multiple = F),
    
                hr(),
                radioButtons("element_Aggregation_Mode", h4(strong("2. Choose an Aggregation Mode")), c("Additiv" = 1, "Multiplicativ" = 2), selected = 1),
                 
                hr(),
                radioButtons("element_Aggregation_Type", h4(strong("3. Choose an Aggregation Type")), c("All" = 1, "Individual" = 2), selected = 1)
            )
          ),

            # aggregation curves
          mainPanel(tags$style(type="text/css",".shiny-output-error {visibility: hidden;}",
                                  ".shiny-output-error:before {visibility: hidden;}"),
                    tabPanel("Row_Wise",
                      h3(strong("Normal Distributions of selected Elements")),
                      plotOutput("element_aggregation_row_wise_plot")
                    )
          )
        )
      ),
             
             
      # Aggregation of Quotients
      tabPanel("aggregation of quotients", titlePanel("Aggregation of Quotients"), value = "AggregationQuotient",
        sidebarLayout(
          sidebarPanel(
            tabPanel("Row_Wise",
                br(),
                selectInput("quotient_Row_For_Aggregation", h4(strong("1. Select Element(s) for Aggregation Computation")), choices = "", multiple = F),

                hr(),
                radioButtons("quotient_Aggregation_Mode", h4(strong("2. Choose an Aggregation Mode")), c("Additiv" = 1, "Multiplicativ" = 2), selected = 1),

                hr(),
                radioButtons("quotient_Aggregation_Type", h4(strong("3. Choose an Aggregation Type")), c("All" = 1, "Individual" = 2), selected = 1)
            )
          ),

            # aggregation curves
          mainPanel(tags$style(type="text/css",".shiny-output-error {visibility: hidden;}",
                                  ".shiny-output-error:before {visibility: hidden;}"),
                    tabPanel("Row_Wise",
                      h3(strong("Normal Distributions of selected Elements")),
                      plotOutput("quotient_aggregation_row_wise_plot")
                    )
          )
        )
      )
    )
  )
)
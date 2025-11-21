library(shiny)
library(shinyFiles)
library(plotly)
library(glue)
library(rsconnect)
library(DT)

# Define UI ----
# Define UI for random distribution app ----


# TAB LAYOUT

intro_tab_ui <- function(){
  tabPanel(
    "Intro",
    sidebarLayout(
      
      sidebarPanel(
        "", width=3,
        
        h3(strong("Variable requirements"), style = "font-size:170%"),
        p("Please read the following variable input requirements before proceeding with the app", 
          style = "font-size:150%"),
        
        # Dilutions
        h3(tags$u("Dilutions")),
        p("- dilutions should be formatted as 1/100, 1/200, etc.", style = "font-size:130%"),
        p("- dilutions should be listed from highest to lowest.", style = "font-size:130%"),
        p("- Correct example format: 1/3200, 1/1600...", style = "font-size:130%"),
        
        # Standard curve
        h3(tags$u("Standard Curve")),
        
        # Background samples
        h3(tags$u("Background samples"))
      ),
      
      mainPanel(
        fluidRow(
          column(12, align="center",
                 img(src='logo.png', width=200, height=200))
        ),
        
        fluidRow(
          column(12, align="center",
                 h1(strong("Luminex Data Analysis (LDA)"))
          )
        ),
        
        # Intro to the app
        p("Welcome to the Luminex Data Analysis app. This app is designed to help you pre-process your raw MBA Luminex dataset into an analysable format. Please read the variable input requirements before proceeding with the app.",
          style = "font-size:140%"),
        
        hr(),
        
        # About
        h2(strong("About")),
        p("This R Shiny app was created by researchers in the Drakeley group at the London School of Hygiene and Tropical Medicine to standardise the quality control steps of analysing raw multiplex datasets from a Luminex machine.", 
          style = "font-size:140%"),
        
        div(
          strong("Composition"),
          tags$ul(
            tags$li(strong("Bead count:"), " Visualise wells with low bead counts on each plate. Default threshold is 30 beads."),
            tags$li(strong("Standard curve:"), " Plots MFI values against log dilution values."),
            tags$li(strong("Coefficient of variation plots:"), " Measures variation between samples in a plate."),
            tags$li(strong("Levey-Jennings plots:"), " Tracks control samples over time."),
            tags$li(strong("Normalisation:"), " Seeks differences between units or batches."),
            tags$li(strong("Loess:"), " Visualise differences between raw and normalised plate."),
            tags$li(strong("Download:"), " Data can be cleaned and downloaded for downstream analysis.")
          ), style = "font-size:140%"
        ),
        
        p("No coding experience is necessary to interpret and use this app.", 
          style = "font-size:140%"),
        
        hr(),
        
        # Contact information
        h2(strong("Contact"))
      )
    )
  )
}


upload_tab_ui<-function(){
  tabPanel(
        "Upload", 
        sidebarLayout(
            sidebarPanel("", width=3,
                         
                         # File Upload
                         h5(strong("File upload"), style = "font-size:120%"),
                           
                           # Adds text
                           p("Upload your data using the file upload button below.", style = "font-size:120%"),
                           p(
                             fileInput(
                               inputId = "fileUpload",
                               label = "",
                               multiple = TRUE,
                               buttonLabel = "Browse...",
                               placeholder = "No file selected"
                               ), 
                             style = "font-size:120%"),
                           
                           # Inserts a horizontal line (a divider in the UI)
                           hr(),
                         
                         # Dilutions
                         h5(strong("Dilutions"), style = "font-size:120%"),
                         
                           # Adds a paragraph of text under the heading
                           p(
                             textOutput("dilutions"),
                             textInput(inputId = "dilutionInput",label=""),
                             actionButton("dilutionButton","Submit dilutions (comma seperated)"),
                             style = "font-size:120%"
                           ),
                           
                           # Inserts a horizontal line (a divider in the UI)
                           hr(),
      
                         # Standard curves
                         h5(strong("Standard curve"), style = "font-size:120%"),
                         
                           # Adds a paragraph of text under the heading
                           p(
                             textOutput("std_label"),
                             textInput(inputId = "std_labelInput",label=""),
                             actionButton("std_labelButton","Submit standard curves"),
                             style = "font-size:120%"
                           ),
                           
                           # Inserts a horizontal line (a divider in the UI)
                           hr(),
                         
                         # Background
                         h5(strong("Background samples"), style = "font-size:120%"),
                         
                           # Adds a paragraph of text under the heading
                           p(
                             textOutput("bkg_label"),
                             textInput(inputId = "bkg_labelInput",label=""),
                             actionButton("bkg_labelButton","Submit Background samples"),
                             style = "font-size:120%"
                           )
            ),
      
                         
            ### SEE IF YOU NEED THIS   
                         # # Add level 5 header for Background adjustment
                         # h5(strong("Background adjustment"), style = "font-size:120%"),
                         # 
                         #   # Adds a paragraph of text under the heading
                         #   p(
                         #     "Check this box if you want to adjust for background. This takes",
                         #     "the mean of all background readings for all plates for each antigen",
                         #     "in turn, and subtracts this from all other readings for that antigen"
                         #     , style = "font-size:120%"),
                         #   
                         #   checkboxInput("background_adjustment", "Adjust for background", value = FALSE), # subtraction. do the same with division
                         #  
                         # ),
            mainPanel(
              
              fluidRow(
                column(12, align="center",
                       h1("Luminex Data Analysis"))
              ),
              
              p("Please upload your data using the file upload button on the left. You can upload multiple files at once. Once you have uploaded your data, you can use the tabs on the top of the app to visualise and process your data.",
                style = "font-size:140%"),
              
              hr(),
              
              h5("Upload summary", style = "font-size:150%"),
              
              # Adds text under the heading
              p(tableOutput("upload_summary"), style = "font-size:150%"),
              
              fluidRow(
                column(3, tableOutput("plate_table")),
                column(3, tableOutput("antigen_table"))
            )
        )
  )
  )
}

bead_count_ui <- function() {
  
  tabPanel("Bead Count",
           sidebarLayout(
             sidebarPanel(width=3,
                          
                          # Add level 5 header for Bead threshold
                          h5(strong("Bead count threshold"), style = "font-size:120%"),
                          
                          p(
                            "Low bead counts may not ensure statistically valid results.",
                            "Please add the minimum number beads required to proceed with analysis (usually 30 or 50).", style = "font-size:120%"),
                          
                          # User input for bead threshold
                          numericInput(inputId = "bead_threshold",
                                       label = "Please enter a minimum bead count: ",
                                       value = 30, # default
                                       min = 0),
                          
                          # Plate selection
                          selectInput(inputId = "plate_selection",
                                      label = "Select Plate:",
                                      choices = NULL)  # will populate dynamically in server

             ),
             mainPanel(plotOutput("bead_count_plots", height = "1200px")))
           )
}


std_curve_plots_ui <- function(){
    
    tabPanel("Standard Curves", 
        sidebarLayout(
            sidebarPanel("", width=3,
                
                         p("Standard curves show MFI plotted against the log of the dilution. ", 
                           "These plots can be used to determine the MFI of a sample at a ",
                           "given dilution."
                           )
            ),
            mainPanel(
              DT::dataTableOutput("std_table_test")
              #uiOutput("std_curve_plot")   
              
            )
        )
    )

}

levey_jennings_plots_ui <- function(){

    tabPanel("Variation", 
        sidebarLayout(
            sidebarPanel("", width=3,
                p(
                    "The coefficient of variation is a measure of the variation between samples",
                    "in a plate. It is calculated as the standard deviation divided by the mean."
                ),
                p(
                    "Levey-Jennings plots tracks how control samples behave over time across different plates.",
                    " The Levey-Jennings plot is a plot shows the median MFI for each antigen and",
                    " plate. The lines indicate 1 and 2 standard deviations from the mean. ",
                    "This plot can be used to identify antigens that are not performing well."
                ),
                hr(),
                h5("Levey-Jennings controls"),
                uiOutput("levey_jennings_select_points")
            ),
            mainPanel(
                h4("Coefficient of variation plots"),
                uiOutput("coefv_plots"),
                h4("Levey-Jennings plots"),
                uiOutput("levey_jennings_plots")
            )
        )
    )
}

normalisation_ui <- function(){
    tabPanel("Normalisation", 
        sidebarLayout(
            sidebarPanel("", width=3,
                
                         p("A curve is fitted to the control samples for each antigen. The curve",
                           "is then used to estimate 100 MFI values at 100 incremental dilutions",
                           "within curve's original range.", style = "font-size:130%"),
                         
                         p("Different ways of fitting the curve have been implemented. Please",
                           "choose one and wait a few seconds for the visualisation to appear:", style = "font-size:130%"),
                         
                         p(strong("poly method"),": fits a cubic polynomial curve to the data", style = "font-size:130%"),
                         p(strong("nls method"), ": uses non-linear least squares fitting using the nls2 library", style = "font-size:130%"),
                         p(strong("kernsmooth method"), ": estimates a regression function using local polynomials using the KernSmooth library", style = "font-size:130%"),
                         p(strong("cgam method"), ": fits a generalised additive model using the cgam library", style = "font-size:130%"),
                         uiOutput("normControls")
                         ),
            mainPanel(
                uiOutput("normalisation_plots"),
                # uiOutput("normalised_ag_plot")
            )
        )
    )
}

loess_ui <- function(){
    tabPanel("Loess", 
        sidebarLayout(
            sidebarPanel("", width=3,
                p(
                    "MFIs are subtracted between the chosen reference and each",
                    "target plate curve at 100 dilutions for each antigen,",
                    "producing Δ MFI values for the dilutions. A loess",
                    "curve is then fitted to the Δ MFI values for each",
                    "plate/antigen. The loess curve is then used to adjust the",
                    "raw MFI values for each plate/antigen.",
                    ),
                p(paste(
                    "Select a plate to use as a reference and use the tabs to",
                    "visualise the outputs"
                )),
                uiOutput("loessControls1"),
                uiOutput("loessControls2"),
                uiOutput("download_button")
            ),
            mainPanel(
                tabsetPanel(
                    id="loess_tabs",
                    type = "tabs",
                    tabPanel("loess_plots",
                        br(),
                        p(paste(
                            "These plots visualise the difference between",
                            "the reference raw data and the normalised plate"
                        )),
                        uiOutput("loess_plots")
                    ),
                    tabPanel("raw_comparison",
                        br(),
                        p(paste(
                            "These plots visualise the loess adjusted MFI",
                            "value for each plate and each antigen. Each",
                            "line represents an individual plate"
                        )),
                        uiOutput("loess_ag_plot"),
                        uiOutput("loess_data")
                    )
                )

            )
        )
    )
}


download_data_ui <- function() {
  
  tabPanel("Download data",
           sidebarLayout(
             sidebarPanel(width=3,
                         
                          # Add description that this page is to adjust for data if happy. After all the steps if it looks good then download
                          
                          # Adjusting low bead count
                          h5(strong("Adjusting low bead count"), style = "font-size:120%"),
                            
                            # User input for bead threshold
                            numericInput(inputId = "bead_threshold",
                                         label = "Please enter a minimum bead count: ",
                                         value = 30, # default
                                         min = 0),
                            
                            # Adds a paragraph of text under the heading
                            p("Check this box if you want to remove samples with a low bead count for the cleaned dataset",
                              style = "font-size:120%"),
                            
                            checkboxInput("bead_adjustment", "Remove low beads", value = FALSE),
                          
                          # Inserts a horizontal line (a divider in the UI)
                          hr(),
                          
                          # Background adjustment
                          h5(strong("Remove background"), style = "font-size:120%"),
                          
                            # Adds text
                            p("Please select how you would like to remove the background MFI", style = "font-size:120%"),
                            
                            # Adds a paragraph of text under the heading
                            p(
                              radioButtons(
                                inputId = "backgroundremoval",
                                label = "Choose one option:",
                                choices = c("Subtract", "Divide"),
                                selected = "Subtract"  # default choice
                              )
                            )
                          
             ),
             mainPanel()
             )
  )
}


ui <- fluidPage(
  
  # Increase the size of the tab panel
    tags$head(tags$style(HTML("
      /* Tab headers */
      .nav-tabs > li > a {
        font-size: 18px !important;       /* text size */
        padding: 15px 20px !important;    /* increase height & horizontal padding */
      }
  
      /* Active tab */
      .nav-tabs > li.active > a,
      .nav-tabs > li.active > a:focus,
      .nav-tabs > li.active > a:hover {
        font-size: 18px !important;
        padding: 15px 20px !important;
      }
  
      /* Space between tabs and content */
      .tab-content {
        padding-top: 20px !important;
      }
    "))),
  
  # Tabset panel
    tabsetPanel(
      id = "tabs_id",
      type = "tabs",
      intro_tab_ui(),
      upload_tab_ui(),
      bead_count_ui(),
      std_curve_plots_ui(),
      levey_jennings_plots_ui(),
      normalisation_ui(),
      loess_ui(),
      download_data_ui()
    )
)


# Define server logic ----

parse_dilution_string <- function(input){
    trimws(unlist(strsplit(trimws(input),",")))
}
  
server <- function(input, output, session) {
  
  # Read in functions
    source("luminex_dx_functions_v11.R")
    source("normalisation_functions.R")
  
  
  # Set up dilutions
    # Defines a default set of dilutions
    dilutions <- c("1/31250", "1/6250", "1/1250", "1/1000", "1/250", "1/50", "1/10")
    
    # Add text to home button on what the value the dilutions are set to
    output$dilutions <- renderText(paste("The default dilutions are set to", paste(dilutions,collapse=", "),". If you have your own dilutions, then you can manually update this here. Please separate the dilutions with a comma."))
  
    # Reactive dilutions (default or user input)
    user_dilutions <- reactiveVal(dilutions)
    
    observeEvent(input$dilutionButton, {
      parsed <- parse_dilution_string(input$dilutionInput)
      if(length(parsed) > 0){
        user_dilutions(parsed)
        output$dilutions <- renderText(paste("The dilutions are set to", paste(parsed, collapse=", "),". You can update this here."))
      }
    })
  
      # # Defines a default set of dilutions
      # dilutions <- c("1/31250", "1/6250", "1/1250", "1/1000", "1/250", "1/50", "1/10")
      # 
      # output$dilutions <- renderText(paste("The default dilutions are set to", paste(dilutions,collapse=", "),". If you have your own dilutions, then you can manually update this here. Please separate the dilutions with a comma."))
      # 
      # # Update dilutions interactively (flexibility)
      # observeEvent(input$dilutionButton, {
      #     
      #     # Reads user input (dilutionInput) and parses it into a proper format
      #     dilutions <- parse_dilution_string(input$dilutionInput)
      #     output$dilutions <- renderText(paste("The dilutions are set to", paste(dilutions,collapse=", "),". You can update this here."))
      # })
  
  # Set up standard curves
    std_label <- reactiveVal("(CP3|Std Curve|WHO)")
    
    output$std_label <- renderText({paste("The default Standard Curves are set to variables that contain CP3, Std Curve, WHO. If you have your own standard curves, then you can manually update this here. Please seperate the standard curves with a comma.")})
    
    # Update std_label interactively (flexibility)
    observeEvent(input$std_labelButton, {
      
      # Parse user input: split by comma, trim whitespace
      labels <- trimws(unlist(strsplit(input$std_labelInput, ",")))
      
      if(length(labels) > 0){
        # Combine into regex pattern
        regex_pattern <- paste(labels, collapse = ",")
        std_label(regex_pattern)  # update reactiveVal
        
        # Update displayed text
        output$std_label <- renderText({
          paste("The standard curve labels are now set to: ", paste(labels, collapse = ", "))
        })
      }
    })
    #   # Reads user input (std_labelInput) and parses it into a proper format
    #   std_label <- parse_dilution_string(input$std_labelInput)
    #   output$std_label <- renderText(paste("The standard curve labels are set to ", paste(std_label,collapse=", "),". You can update this here."))
    # })
    
    cleaned_data <- reactive({
      req(input$fileUpload)  # make sure a file is uploaded
      files <- input$fileUpload$datapath
      
      # Call your read.batch.std function
      read.batch.std(path = dirname(files[1]), std_label = std_label)
    })
   
    
  # Set up background values 
    bkg_label <- ("(blank|background)") # contains any cell that has blank or background
    
    output$bkg_label <- renderText({paste("The default background samples are set to variables that contain <i>blank</i> or <i>background</i>. If your background samples are labelled differently, then you can manually update this here. Please seperate each background sample name with a comma.")})
    
    # Update bkg_label interactively (flexibility)
    observeEvent(input$bkg_labelButton, {
      
      # Parse user input: split by comma, trim whitespace
      labels <- trimws(unlist(strsplit(input$bkg_labelInput, ",")))
      
      if(length(labels) > 0){
        # Combine into regex pattern
        regex_pattern <- paste(labels, collapse = ",")
        bkg_label(regex_pattern)  # update reactiveVal
        
        # Update displayed text
        output$bkg_label <- renderText({
          paste("The background samples are now set to: ", paste(labels, collapse = ", "))
        })
      }
    })
    
  # File upload and data preparation
    d <- reactive({
      
      # Make sure a file has been uploaded
      req(input$fileUpload)
      
      # Extract the folder path and the file label, drop .csv extension
      raw_data_path <- dirname(input$fileUpload[1,"datapath"])
      plate_lab <- substr(input$fileUpload$name,1,nchar(input$fileUpload$name)-4)
      
      # Reading uploaded data
      d = list()
      d$rawdatapath <- raw_data_path
      d$uploadedfilenumber <- nrow(input$fileUpload)
      
        ## Load batch data from the uploaded files
        d$plates <- read.batch(path=raw_data_path)
        
        ## Assign user-friendly names to the plates
        names(d$plates)<-plate_lab
      
      # Join plates
        ## Store metadata: plate names, number of plates
        d$platenames <- names(d$plates)
        d$numplates <- length(d$platenames)
      
        ## Combine plates into one dataset
        d$combinedplates <- join.plates(d$plates)
        
          ## Skip NAs
          d$combinedplates <- na.omit(d$combinedplates)
          for (i in seq_along(d$plates)){
            d$plates[[i]] <- na.omit(d$plates[[i]])
          }
        
        ## Extracts antigen names (columns 3 onward = data columns)
        d$ags <- colnames(d$combinedplates)[3:ncol(d$combinedplates)]
        
        ## Extracts Sample names (columns 3 onward = data columns)
          d$sample_type <- colnames(d$combinedplates)[2]
          
          # Replace any cell containing "Unknown" (case-insensitive) with "Sample"
          d$sample_type <- ifelse(grepl("Unknown", d$sample_type, ignore.case = TRUE), "Sample", d$sample_type)[2]
      
        ## Metadata with dates - include date & plate info
        d$datesplates <- read.batch(path=raw_data_path,inc_date=T,inc_plate=T)
        d$combineddatesplates <- join.plates(d$datesplates)
        d
    })
    
    # Data preparation beads - Reactive for bead count
    bead_data <- reactive({
        # Make sure a file has been uploaded
        req(input$fileUpload)
        
        # Extract the folder path and the file label, drop .csv extension
        raw_data_path <- dirname(input$fileUpload[1,"datapath"])
        plate_lab <- substr(input$fileUpload$name,1,nchar(input$fileUpload$name)-4)
        
        # Reading uploaded data
        bead_data = list()
        bead_data$rawdatapath <- raw_data_path
        bead_data$uploadedfilenumber <- nrow(input$fileUpload)
        
        ## Load batch data from the uploaded files
        bead_data$plates <- read.batch.beads(path = raw_data_path)
        
        ## Assign user-friendly names to the plates
        names(bead_data$plates) <- plate_lab
        
        # Join plates
        ## Store metadata: plate names, number of plates
        bead_data$platenames <- names(bead_data$plates)
        bead_data$numplates <- length(bead_data$platenames)
        
        ## Combine plates into one dataset
        bead_data$combinedplates <- join.plates(bead_data$plates)
        
        ## Skip NAs
        bead_data$combinedplates <- na.omit(bead_data$combinedplates)
        for (i in seq_along(bead_data$plates)){
          bead_data$plates[[i]] <- na.omit(bead_data$plates[[i]])
        }
        
        ## Extracts antigen names (columns 3 onward = data columns)
        bead_data$ags <- colnames(bead_data$combinedplates)[3:ncol(bead_data$combinedplates)]
        
        ## Metadata with dates - include date & plate info
        bead_data$datesplates <- read.batch(path=raw_data_path,inc_date=T,inc_plate=T)
        bead_data$combineddatesplates <- join.plates(bead_data$datesplates)
        bead_data
        })
    
    # Data preparation - standard curve
    std_curve_data <- reactive({
      req(input$std_labelInput)   # ensures input is not NULL
      
      # Parse and clean user input
      std_labels <- trimws(unlist(strsplit(input$std_labelInput, ",")))
      
      # Filter standard rows dynamically
      std_rows <- d()$combinedplates$Sample %in% std_labels
      standard_data <- d()$combinedplates[std_rows, ]
      
      # Pass filtered data to get.standard
      get.standard(
        data = standard_data,
        std_label = NULL,          # no need, data is already filtered
        dilutions = user_dilutions(),
        n_points = length(std_labels)
      )
    })
## MIKA NOTES - MAYBE MAKE DATASET WHERE FIRST COLUMN IS THE STD CURVE VARIABLE, SECOND IS THE DILUTION AND THEN MERGE THAT WITH THE MEDIAN MFI TABLE FOR ALL TESTED ANTIGENS

  
    # Create a reactive that extracts antigen column names from the uploaded data for app use
    ag <- reactive({
        colnames(d()$plates[[1]])[3:ncol(d()$plates[[1]])]
    })

    # Normalisation
    nd <- reactive({
        # ag_list <- colnames(d()$plates[[1]])[3:ncol(d()$plates[[1]])]
        progress <- Progress$new(session, min=1, max=d()$numplates)
        on.exit(progress$close())

        progress$set(message = 'Calculation in progress',
                    detail = 'This may take a while...')
        
        # Use normalise_plates function
        normalise_plates(
            all_plates = d()$plates,
            dilutions = user_dilutions(),
            ag_list = ag(),
            fit_type = input$fitting_function,
            progress = progress )
    })

    # Loess adjustment - apply loess regression across plates to smooth out systematic plate-to-plate variation
    ld <- reactive({
        
        # reference plate selected by the user
        req(input$loess_ref_plate)
      
        progress <- Progress$new(session, min=1, max=d()$numplates)
        on.exit(progress$close())

        progress$set(message = 'Calculation in progress',
                    detail = 'This may take a while...')

        d <- loess_adjustment(
            all_plates = d()$plates,
            fitted_data = nd()$data,
            ref_plate = input$loess_ref_plate,
            ag_list = ag(),
            progress = progress
        )
        d$data <- lapply(d$data, as.data.frame)
        d
    })
    
    # Bead threshold adjustment
      ## Set default beads
      bead_control <- reactiveValues(threshold = 30)
      
      ## Update bead threshold from 30 if the user defines an alternative amount
      observeEvent(input$bead_threshold, {
        # Only update if numeric and not NA
        if (!is.null(input$bead_threshold) && !is.na(as.numeric(input$bead_threshold))) {
          bead_control$threshold <- as.numeric(input$bead_threshold)
        }
      })
      
      # Update plate selection dropdown
      observe({
        req(bead_data())
        updateSelectInput(
          session,
          "plate_selection",
          choices = names(bead_data()$plates),
          selected = names(bead_data()$plates)[1]
        )
      })

    # Upload summary: number of files, number of plates, antigens
          output$upload_summary<- renderText({
              if (!is.null(input$fileUpload)){
                  paste(d()$uploadedfilenumber," files were uploaded, leading to ", d()$numplates," plates being loaded with the following plates and antigens:")
              } else {
                  "Please upload a file using the upload button on the left sidebar."
              }
          })
      
          # Tables of antigen names
          output$antigen_table<-renderTable({
              req(input$fileUpload)
              data.frame(Antigens=d()$ags)
          })
          
          # Tables of Sample names
          output$sample_type_table<-renderTable({
            req(input$fileUpload)
            data.frame(Sample=d()$sample_type)
          })
          
          # Tables of plate names
          output$plate_table<-renderTable({
              req(input$fileUpload)
              data.frame(Plates=d()$platenames)
          })

    output$download_button <- renderUI({
        req(ld())
        downloadButton("download_combined", "Download normalised dataset")
    })

    # After normalisation and adjustment, allows the user to download the final dataset as normdata.csv
    output$download_combined <- downloadHandler(
        
        filename = function() {
        paste0("normdata.csv")
        },
        content = function(file) {
        write.csv(ld()$data, file)
        }
    )

    output$raw_data_files <- renderTable({
        req(input$fileUpload)
        as.data.frame(names(d()$plates))

    })


    output$moreControls <- renderUI({
        req(input$fileUpload)
        tagList(
            selectInput("ref_plate", "Reference plate", d()$platenames),
            selectInput("standard_plate", "Reference plate", d()$platenames),
            selectInput("selected_ag", "Select antigen", d()$ags) )
    })

    output$normControls <- renderUI({
        req(input$fileUpload)
        tagList(
            selectInput("fitting_function","Select fitting function",c("poly","nls","kernsmooth","cgam")),
            renderText("fitting_function_text"),
            selectInput("norm_plot_plate", "Select plate to visualise", d()$platenames), )
    })

    output$fitting_function_text <- renderText({
        req(input$fitting_function)
        if (input$fitting_function=="poly"){
            "The poly method fits a cubic polynomial curve to the data"
        } else if (input$fitting_function=="nls"){
            "The nls method uses non-linear least squares fitting using the nls2 library"
        } else if (input$fitting_function=="kernsmooth"){
            "The kernsmooth method estimates a regression function using local polynomials using the KernSmooth library"
        } else if (input$fitting_function=="cgam"){
            "The cgam method fits a generalised additive model using the cgam library"
        }
    })
  
    output$loessControls1 <- renderUI({
        req(input$fileUpload)   
        tagList(
            selectizeInput(
                inputId = "loess_ref_plate",
                label = "Select reference plate",
                choices = d()$platenames,
                options = list(
                    onInitialize = I('function() { this.setValue(""); }') )
            ) )
    })

    output$loessControls2 <- renderUI({
        req(input$loess_ref_plate)   
        tagList(
            selectizeInput(
                inputId = "loess_plot_plate",
                label = "Select plate to visualise",
                choices = d()$platenames) )
    })

    # Standard Curves
    output$cleaned_table <- DT::renderDataTable({
      req(all_std_data())
      DT::datatable(
        all_std_data(),
        options = list(pageLength = 10, scrollX = TRUE),
        rownames = FALSE
      )
    })
    
    
    output$std_curve_plot <- renderPlotly({
      req(std_curve_data())
      
      plot.std.curve3_interactive(
        data = list(std_curve_data()),
        antigen = ag(),        
        dilutions = user_dilutions(),
        plate_labels = input$std_curve_input
      )
    })
          # ORIGINAL CODE
          # output$std_curves <- renderUI({
          #     req(input$fileUpload)
          #     fluidRow(
          #         lapply(ag(),function(antigen){
          #             id <- paste0("plot_std_", antigen)
          #             plotOutput(outputId = id)
          #             
          #             std_curves <- lapply(d()$platenames, function(x) {
          #               
          #               # Generates interactive standard curve per antigen and plate  
          #               get.standard(
          #                     data = d()$plates[[x]],
          #                     std_label = std_label(), # standard curves are defined as those that contain CP3, the word Std Curve, or WHO (for WHO)
          #                     dilutions = dilutions,
          #                     n_points = length(dilutions))
          #               })
          #             
          #             names(std_curves) <- d()$platenames
          # 
          #             column(6,
          #                 renderPlotly({
          #                     # Generates interactive standard curve plots per antigen and plate
          #                     plot.std.curve3_interactive(
          #                         data = std_curves,
          #                         antigen = antigen,
          #                         dilutions = dilutions,
          #                         plate_labels = input$std_curve_input )
          #                 }) )
          #         }) )
          # })
      
          # Allow user to select up to 5 plates for comparison (compare standard curves visually) 
          output$std_curve_plate_selection <- renderUI({
              selectizeInput(
                  inputId="std_curve_input",
                  label="Select plates (max 5)",
                  choices=d()$platenames,
                  selected = NULL,
                  multiple=TRUE,
                  options = list(maxItems = 5) ) # Can remove this if we want unlimited plate selection
          })

    # QC
    ## Levey-Jennings controls for tracking high, mid, and low dilution points across plates
    output$levey_jennings_select_points <- renderUI({
        tagList(
            selectInput("levey_high", "Select high dilution", user_dilutions(), selected="1/10"),
            selectInput("levey_mid", "Select mid dilution", user_dilutions(), selected="1/50"),
            selectInput("levey_low", "Select low dilution", user_dilutions(),selected="1/250") )
    })

    ## Plots coefficient of variation across dilutions to see consistency of replicates
    output$coefv_plots <- renderUI({
        req(input$fileUpload)
        fluidRow(
            column(12,
                renderPlotly({
                    plot.coefv(
                        data = d()$plates,
                        ag(),
                        std_label = std_label(),
                        dilutions=user_dilutions(),
                        target_dilution = "1/250" ) 
                })  ),
            column(12,
                renderPlotly({
                    plot.coefv(
                        data = d()$plates,
                        ag(),
                        std_label = std_label(),
                        dilutions=user_dilutions(),
                        target_dilution = "1/1250" )
                }) )
        ) })

    ## Generate Levey-Jennings plots for each antigen
    output$levey_jennings_plots <- renderUI({
        req(input$fileUpload)
        fluidRow(
            lapply(ag(), function(antigen) {
                column(6,
                    renderPlot({
                        levey.jennings(
                            data = d()$combineddatesplates,
                            std_label = std_label(),
                            blank_label = bkg_label, # "Background0",
                            dil_high = input$levey_high,
                            dil_mid = input$levey_mid,
                            dil_low = input$levey_low,
                            by_var = "plate",
                            ag_list = c(antigen) )
                }) )
            }) )
    })
    
    # Generate low bead plot for each antigen
    output$bead_count_plots <- renderPlot({
      req(bead_data())
      req(input$plate_selection)
      
      plate_selected <- bead_data()$plates[[input$plate_selection]]
      
      plot_beads_by_well(
        plate_selected,
        plate_name = input$plate_selection,
        threshold = bead_control$threshold
      )
    })
    
    # Normalised Data & Plots
    ## Shows normalised matrix
    output$normalised_matrix<-renderTable({
        nd()$data[["MP01"]]
    })
    
    ## Produces normalisation diagnostic plots per antigen and plate
    output$normalisation_plots <- renderUI({
        req(input$fileUpload)
        fluidRow(
            lapply(ag(), function(a) {
                column(4,
                    renderPlot({
                        nd()$plots[[input$norm_plot_plate]][[a]]

                    }) )
            }) )
    })

    # Loess Adjustment Visualisation - Compares raw vs loess-adjusted plots for each antigen and plate
    output$loess_ag_plot <- renderUI({
        req(input$fileUpload)
        fluidRow(
            lapply(ag(), function(a) {
                column(6,
                    renderPlot({
                        par(mfrow=c(1,2))

                        # Compares raw vs loess-adjusted antigen plots
                        plot.ag.plates.raw(data=d()$plates, ag=a, dilutions=user_dilutions())
                        plot.ag.plates.raw(data=ld()$data, ag=a, dilutions=user_dilutions())
                        par(mfrow=c(1,1))
                    }) )
            }) )
    })

    # Loess-adjusted diagnostic plots
    output$loess_plots <- renderUI({
        req(input$fileUpload)
        req(input$loess_ref_plate)
        fluidRow(
            lapply(ag(), function(a) {
                column(6,
                    renderPlot({
                        ld()$plots[[input$loess_plot_plate]][[a]]

                    }) )
            }) )
    })

    # Display loess-adjusted dataset as a table
    output$loess_data <- renderTable({
        ld()$data[[1]]
    })}


    # # Downloading cleaned data
    # output$cleaned_data <- renderDataTable({
    #   
    #   # Remove background
    #   
    #     ## Create variable that is the selected background removal
    #     bkg_select <- reactive({
    #       input$backgroundremoval
    #     })

          ## Create cleaned dataset

    # })
  


# Run the app ----
shinyApp(ui = ui, server = server)

# rsconnect::deployApp(
#   appName = "Luminex_Data_Analysis"
# )
library(shiny)
library(shinyFiles)
library(plotly)
library(glue)
# Define UI ----
# Define UI for random distribution app ----

upload_tab_ui<-function(){
  tabPanel(
        "Upload", 
        sidebarLayout(
            sidebarPanel("", width=3,
                h5("Dilutions"),
                p(
                    textOutput("dilutions"),
                    textInput(inputId = "dilutionInput",label=""),
                    actionButton("dilutionButton","Submit dilutions")
                ),
                hr(),
                h5("File upload"),
                p("Upload your data using the file upload button below."),
                p(
                    fileInput(
                        inputId = "fileUpload",
                        label = "",
                        multiple = TRUE,
                        buttonLabel = "Browse...",
                        placeholder = "No file selected"
                    )
                ),
                hr(),
                h5("Background adjustment"),
                p(
                    "Check this box if you want to adjust for background. This takes",
                    "the mean of all background readings for all plates for each antigen",
                    "in turn, and subtracts this from all other readings for that antigen"
                ),
                checkboxInput("background_adjustment", "Adjust for background", value = FALSE)

            ),
            mainPanel(
                fluidRow(
                    column(12, align="center",
                        img(src='logo.png', width=200, height=200)  
                    ),
                ),
                fluidRow(
                    column(12, align="center",
                        h1(
                            "Luminex Data analysis"
                        ),
                    ),
                ),
                p(
                    "Welcome to the Luminex data analysis app. This app is designed",
                    "to help you analyse your Luminex data. Please upload your data",
                    "using the file upload button on the left. You can upload multiple",
                    "files at once. Once you have uploaded your data, you can use the",
                    "tabs on the top of the app to visualise and process your data."
                ),
                hr(),
                h5("Upload summary"),
                p(
                    tableOutput("upload_summary")
                ),
                fluidRow(
                    column(6,
                        tableOutput("plate_table")
                    ),
                    column(6,
                        tableOutput("antigen_table")
                    )
                )
                
                
            )
          
        )
      )
}

std_curve_plots_ui<-function(){

    
    tabPanel("Standard Curves", 
        sidebarLayout(
            sidebarPanel("", width=3,
                p(
                    "Standard curves show MFI plotted against the log of the dilution. ",
                    "These plots can be used to determine the MFI of a sample at a ",
                    "given dilution."

                )
            ),
            mainPanel(
                uiOutput("std_curves")
            )
        )
    )

}

levy_jennings_plots_ui<-function()(

    tabPanel("Variation", 
        sidebarLayout(
            sidebarPanel("", width=3,
                p(
                    "The coefficient of variation is a measure of the variation between samples",
                    "in a plate. It is calculated as the standard deviation divided by the mean."
                ),
                p(
                    "The Levy Jennings plot is a plot shows the median MFI for each antigen and",
                    "plate. The lines indicate 1 and 2 standard deviations from the mean.",
                    "This plot can be used to identify antigens that are not performing well."
                ),
                hr(),
                h5("Levy Jennings controls"),
                uiOutput("levy_jennings_select_points")
            ),
            mainPanel(
                h4("Coefficient of variation plots"),
                uiOutput("coefv_plots"),
                h4("Levy Jennings plots"),
                uiOutput("levy_jennings_plots")
            )
        )
    )
)

normalisation_ui<-function(){
    tabPanel("Normalisation", 
        sidebarLayout(
            sidebarPanel("", width=3,
                p(
                    "A curve is fitted to the control samples for each antigen. The curve",
                    "is then used to estimate 100 MFI values at 100 incremental dilutions",
                    "within curve's original range."),
                p(
                    "Different ways of fitting the curve have been implemented. Please",
                    "choose one and wait a few seconds for the visualisation to appear"
                ),
                p(
                    "The poly method fits a cubic polynomial curve to the data"
                ),
                p(
                    "The nls method uses non-linear least squares fitting using the nls2 library"
                ),
                p(
                    "The kernsmooth method estimates a regression function using local polynomials using the KernSmooth library"
                ),
                p(
                    "The cgam method fits a generalised additive model using the cgam library"
                ),
                uiOutput("normControls")
            ),
            mainPanel(
                uiOutput("normalisation_plots"),
                # uiOutput("normalised_ag_plot")
            )
        )
    )
}

loess_ui<-function(){
    tabPanel("Loess", 
        sidebarLayout(
            sidebarPanel("", width=3,
                p(
                    "MFIs are subtracted between the chosen reference and each",
                    "target plate curve at 100 dilutions for each antigen,",
                    "producing delta MFI values for the dilutions. A loess",
                    "curve is then fitted to the delta MFI values for each",
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

ui <- fluidPage(
  
  # App title ----

  

    
   
    
    # Main panel for displaying outputs ----

    
    # Output: Tabset w/ plot, summary, and table ----
    tabsetPanel(
      id="tabs_id",
      type = "tabs",
      upload_tab_ui(),
      std_curve_plots_ui(),
      levy_jennings_plots_ui(),
      normalisation_ui(),
      loess_ui()
    )
    


)

# Define server logic ----

parse_dilution_string<-function(input){
    trimws(unlist(strsplit(trimws(input),",")))
}
  
server <- function(input, output, session) {
    source("luminex_dx_functions_v11.R")
    source("normalisation_functions.R")
    dilutions <- c("1/31250", "1/6250", "1/1250", "1/250", "1/50", "1/10")
    output$dilutions<-renderText(paste("The dilutions are set to",paste(dilutions,collapse=", "),". You can update this here."))
    observeEvent(input$dilutionButton, {
        dilutions <- parse_dilution_string(input$dilutionInput)
        output$dilutions<-renderText(paste("The dilutions are set to",paste(dilutions,collapse=", "),". You can update this here."))
    })
    
    d<-reactive({
        raw_data_path<-dirname(input$fileUpload[1,"datapath"])

        plate_lab <- substr(input$fileUpload$name,1,nchar(input$fileUpload$name)-4)
        d = list()
        d$rawdatapath<-raw_data_path
        d$uploadedfilenumber<-nrow(input$fileUpload)
        d$plates<-read.batch(path=raw_data_path)
        names(d$plates)<-plate_lab
        d$platenames<-names(d$plates)
        d$numplates<-length(d$platenames)
        d$combinedplates<-join.plates(d$plates)
        d$ags<-colnames(d$combinedplates)[3:ncol(d$combinedplates)]
        d$datesplates <- read.batch(path=raw_data_path,inc_date=T,inc_plate=T)
        d$combineddatesplates <- join.plates(d$datesplates)
        d
    })

    ag<-reactive({
        colnames(d()$plates[[1]])[3:ncol(d()$plates[[1]])]
    })

    nd<-reactive({
        # ag_list <- colnames(d()$plates[[1]])[3:ncol(d()$plates[[1]])]
        progress <- Progress$new(session, min=1, max=d()$numplates)
        on.exit(progress$close())

        progress$set(message = 'Calculation in progress',
                    detail = 'This may take a while...')
        normalise_plates(
            all_plates = d()$plates,
            dilutions = dilutions,
            ag_list = ag(),
            fit_type = input$fitting_function,
            progress = progress
        )
    })

    ld<-reactive({
        req(input$loess_ref_plate)
        progress <- Progress$new(session, min=1, max=d()$numplates)
        on.exit(progress$close())

        progress$set(message = 'Calculation in progress',
                    detail = 'This may take a while...')

        d<-loess_adjustment(
            all_plates = d()$plates,
            fitted_data = nd()$data,
            ref_plate = input$loess_ref_plate,
            ag_list = ag(),
            progress = progress
        )
        d$data<-lapply(d$data, as.data.frame)
        d
    })

    output$upload_summary<- renderText({
        if (!is.null(input$fileUpload)){
            paste(d()$uploadedfilenumber,"files were uploaded, leading to", d()$numplates,"plates being loaded with the following plates and antigens:")
        } else {
            "Please upload a file using the upload button on the left sidebar."
        }
    })

    output$antigen_table<-renderTable({
        req(input$fileUpload)
        data.frame(Antigens=d()$ags)
    })
    output$plate_table<-renderTable({
        req(input$fileUpload)
        data.frame(Plates=d()$platenames)
    })

    output$download_button <- renderUI({
        req(ld())
        downloadButton("download_combined", "Download normalised dataset")
    })

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
            selectInput("selected_ag", "Select antigen", d()$ags)
        )
    })

    output$normControls <- renderUI({
        req(input$fileUpload)
        tagList(
            selectInput("fitting_function","Select fitting function",c("poly","nls","kernsmooth","cgam")),
            renderText("fitting_function_text"),
            selectInput("norm_plot_plate", "Select plate to visualise", d()$platenames),
        )
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
                    onInitialize = I('function() { this.setValue(""); }')
                )
            )
        )
    })

    output$loessControls2 <- renderUI({
        req(input$loess_ref_plate)   
        tagList(
            selectizeInput(
                inputId = "loess_plot_plate",
                label = "Select plate to visualise",
                choices = d()$platenames
            )
        )
    })

    output$std_curves<-renderUI({
        req(input$fileUpload)
        fluidRow(
            lapply(ag(),function(antigen){
                id <- paste0("plot_std_", antigen)
                plotOutput(outputId = id)
                std_curves <- lapply(d()$platenames, function(x) {
                    get.standard(
                        data = d()$plates[[x]],
                        std_label = "CP3",
                        dilutions = dilutions,
                        n_points = 6
                    )
                })
                names(std_curves) <- d()$platenames

                column(6,
                    renderPlotly({
                        plot.std.curve3_interactive(
                            data = std_curves,
                            antigen = antigen,
                            dilutions = dilutions,
                            plate_labels = input$std_curve_input
                        )
                    })
                )
            })
        )
    })

    output$std_curve_plate_selection<-renderUI({
        selectizeInput(
            inputId="std_curve_input",
            label="Select plates (max 5)",
            choices=d()$platenames,
            selected = NULL,
            multiple=TRUE,
            options = list(maxItems = 5)
        )
    })

    output$levy_jennings_select_points <- renderUI({
        tagList(
            selectInput("levy_high", "Select high dilution", dilutions, selected="1/10"),
            selectInput("levy_mid", "Select mid dilution", dilutions, selected="1/50"),
            selectInput("levy_low", "Select low dilution", dilutions,selected="1/250")
        )
    })

    output$coefv_plots <- renderUI({
        req(input$fileUpload)
        fluidRow(
            column(12,
                renderPlotly({
                    plot.coefv(
                        data = d()$plates,
                        ag(),
                        std_label = "CP3",
                        dilutions=dilutions,
                        target_dilution = "1/250"
                    )
                })
            ),
            column(12,
                renderPlotly({
                    plot.coefv(
                        data = d()$plates,
                        ag(),
                        std_label = "CP3",
                        dilutions=dilutions,
                        target_dilution = "1/1250"
                    )
                })
            )
        )
    })

    output$levy_jennings_plots <- renderUI({
        req(input$fileUpload)
        fluidRow(
            lapply(ag(), function(antigen) {
                column(6,
                    renderPlot({
                        levy.jennings(
                            data = d()$combineddatesplates,
                            std_label = "CP3",
                            blank_label = "Background0",
                            dil_high = input$levy_high,
                            dil_mid = input$levy_mid,
                            dil_low = input$levy_low,
                            by_var = "plate",
                            ag_list = c(antigen)
                        )
                })
                )
            })
        )
    })

    output$plate_plan<-renderPlot({
        plate_data <- view.plate.data(data=d()$plates[[input$ref_plate]],antigen=input$selected_ag)
        plot.plate(data1=plate_data$sample.data,antigen=input$selected_ag)


    })
    
    output$normalised_matrix<-renderTable({
        nd()$data[["MP01"]]
    })

    output$normalisation_plots <- renderUI({
        req(input$fileUpload)
        fluidRow(
            lapply(ag(), function(a) {
                column(4,
                    renderPlot({
                        nd()$plots[[input$norm_plot_plate]][[a]]

                    })
                )
            })
        )
    })

    output$loess_ag_plot <- renderUI({
        req(input$fileUpload)
        fluidRow(
            lapply(ag(), function(a) {
                column(6,
                    renderPlot({
                        par(mfrow=c(1,2))
                        plot.ag.plates.raw(data=d()$plates, ag=a, dilutions=dilutions)
                        plot.ag.plates.raw(data=ld()$data, ag=a, dilutions=dilutions)
                        par(mfrow=c(1,1))
                    })
                )
            })
        )
    })

    output$loess_plots <- renderUI({
        req(input$fileUpload)
        req(input$loess_ref_plate)
        fluidRow(
            lapply(ag(), function(a) {
                column(6,
                    renderPlot({
                        ld()$plots[[input$loess_plot_plate]][[a]]

                    })
                )
            })
        )
    })


    output$loess_data<-renderTable({
        ld()$data[[1]]
    })




}
  


# Run the app ----
shinyApp(ui = ui, server = server)
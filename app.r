#//TODO implement sample ordering check
#//TODO better error message when km can't be plotted when no points are significant
#//TODO n_sig as a percentage
#//NOTE if deployment is failing, try reinstalling survivALL via install_github()

library(shiny)
library(Biobase)
library(survivALL)
library(survcomp)
library(survival)
library(ggplot2)
library(crukCIMisc)
library(GGally)
library(ggthemes)

ui = fluidPage(theme = "bootstrap.css",
               titlePanel("survivAPP"),
               sidebarLayout(
                             sidebarPanel(
                                          helpText("Welcome to survivAPP!"),
                                          helpText("Data should be formatted as two tab-delimited .txt 
                                                   files - these are the measurement data
                                                   (such as a matrix of gene expression) and the 
                                                   corresponding clinical data which includes event and
                                                   time-to-event information. The two inputs must be in 
                                                   the same sample ordering. Once loaded, select
                                                   the relevant measure ID and clinical 
                                                   columns containing the survival information
                                                   and click Run. Be sure to click Run for each new calculation"),
                                          helpText("survivALL, optimal Kaplan-meier and 
                                                   summary outputs will then be calculated."),
                                          helpText("Example data can be downloaded
                                                   from the following link - http://bit.ly/2z4r8YH"),
                                          helpText("Details of the full survivALL package can be found at -
                                            https://cran.rstudio.com/web/packages/survivALL/index.html &                                             https://www.biorxiv.org/content/early/2017/10/25/208660"),
                                          fileInput('file1', 'Select measure data (.txt)',
                                                    accept=c('text/csv',
                                                             'text/comma-separated-values,text/plain',
                                                             '.csv')),
                                          fileInput('file2', 'Select clinical data (.txt)',
                                                    accept=c('text/csv',
                                                             'text/comma-separated-values,text/plain',
                                                             '.csv')),

                                          uiOutput("select_measure", placeholder = "asdf"),
                                          uiOutput("select_event"),
                                          uiOutput("select_time"),
                                          actionButton("go", "Run")
                                          ),
                             mainPanel(
                                       tabsetPanel(
                                                   tabPanel("survivALL", plotOutput('p_all')),
                                                   tabPanel("Kaplan-Meier", plotOutput('p_km', 
                                                                                       width = "50%")),
                                                   tabPanel("Summary", tableOutput('t_sum'))
                                                   )

                                       )
               ),
               hr(),
               print("Users are advised that whilst no user data is stored persistantly, shinyapps.io hosts via HIPPA & PHI non-compliant infrastructure. For complete security when analysing proprietary datasets we reccommend running the survivALL R package on a local machine")
                                          
)

server = function(input, output, session){
    myData <- reactive({
            inFile <- input$file1
                if (is.null(inFile)) return(NULL)
                data <- read.delim(inFile$datapath, header = TRUE, row.names = 1)
                data
                  
    })
    myData2 <- reactive({
            inFile <- input$file2
                if (is.null(inFile)) return(NULL)
                data <- read.delim(inFile$datapath, header = TRUE, row.names = 1)
                data
                  
    })

    output$select_measure <- renderUI({
        selectInput("select_measure", 
                    label = "Select measure ID",
                    choices = row.names(myData()))
    })
    output$select_time <- renderUI({
        selectInput("select_time", 
                    label = "Select time-to-event column",
                    choices = colnames(myData2()))
    })
    output$select_event <- renderUI({
        selectInput("select_event", 
                    label = "Select event column",
                    choices = colnames(myData2()))
    })

    observeEvent(input$go, {
        measure <- myData()
        clinical <- myData2()
        srvall <- survivALL(
                measure = as.numeric(measure[input$select_measure,]), 
                srv <- clinical,
                time = input$select_time,
                event = input$select_event,
                bs = c(),
                measure_name = input$select_measure)
        contsig <- checkContSig(srvall$measure,
                                srvall$event_time,
                                srvall$event)

        output$p_all <- renderPlot({
            plotALL(
                    measure = measure[input$select_measure,], 
                    srv <- clinical,
                    time = input$select_time,
                    event = input$select_event,
                    bs = c(),
                    title = input$select_measure)
        })
        output$p_km <- renderPlot({
            Class <- ifelse(srvall$clsf == 0, "low", "high") 
            #stats
            srv_obj <- Surv(as.numeric(srvall$event_time), srvall$event)
            broom::tidy(coxph(srv_obj ~ Class))
            #km-plot
            data_fit <- survfit(srv_obj ~ Class)
            ggsurv(data_fit, surv.col = c("#525252", "#bdbdbd")) +
                theme_pander() + 
                ggtitle(input$select_measure) +
                theme(plot.title = element_text(hjust = 0.5),
                      axis.line = element_line(size = 0.1), 
                      axis.ticks = element_line(size = 0.85),
                      legend.position = 'bottom')
        }) 
        output$t_sum <- renderTable({
            if(sum(!is.na(srvall$log10_p)) == 0){
                sum_dfr <- data.frame(single_test_sig = contsig,
                                      measure = input$select_measure,
                                      n_significant = 0,
                                      optimal_p = NA,
                                      optimal_position = NA)
            } else {
                optimal_pos <- which(diff(srvall$clsf) == 1) + 1
                sum_dfr <- data.frame(single_test_sig = contsig, 
                                      measure = input$select_measure,
                                      n_significant = sum(!is.na(srvall$log10_p)),
                                      optimal_p = srvall$p[optimal_pos],
                                      optimal_position = paste0(round(optimal_pos / 
                                                                      nrow(srvall) * 100, 3),
                                                                "%"))
            }
        })
    })

}

shinyApp(ui, server)

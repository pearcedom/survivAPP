library(shiny)
library(survivALL)
library(Biobase)
library(survival)
library(survcomp)
library(crukCIMisc)
library(ggthemes)
library(GGally)
pal <- sapply(c("Blue", "Magenta", "Teal", "Green", "Orange", "Yellow", "Purple", "Pink", "LightBlue", "Turquoise"), CRUKcol)


ui <- fluidPage(
                #title
                titlePanel("survivAPP"),
                sidebarLayout(
                              sidebarPanel(
                                           fileInput("measure", "select .csv",
                                                     accept = c("text/csv", 
                                                                "text/comma-separtaed-values, 
                                                                text/plain", ".csv"))
                                          textInput("affy_probe", 
                                                      "Select Affymetrix probe ID:", 
                                                      value = "NM_004448", 
                                                      placeholder = "e.g. NM_004448 (ERBB2)"),
                                          ),
                              mainPanel(
                                        plotOutput(outputId = "p_all"), 
                                        plotOutput(outputId = "p_km")
                                        )
                              )
                )

server <- function(input, output) {
    data(nki_subset)

    output$contents <- renderTable({
            # input$file1 will be NULL initially. After the user selects
            # and uploads a file, it will be a data frame with 'name',
            # 'size', 'type', and 'datapath' columns. The 'datapath'
            # column will contain the local filenames where the data can
            # be found.
            inFile <- input$measure
                if (is.null(inFile))
                          return(NULL)
                read.csv(inFile$datapath, header = input$header)
    })

    output$p_all <- renderPlot({
        plotALL(measure = exprs(nki_subset)[input$affy_probe,], 
                srv <- pData(nki_subset),
                time = "t.dmfs",
                event = "e.dmfs",
                bs = c(),
                measure_name = input$affy_probe)
    }) 
    output$p_km <- renderPlot({
        srvall <- survivALL(
                            measure = exprs(nki_subset)[input$affy_probe,], 
                            srv <- pData(nki_subset),
                            time = "t.dmfs",
                            event = "e.dmfs",
                            bs = c(),
                            measure_name = input$affy_probe)

        Class <- ifelse(srvall$clsf == 0, "low", "high") 
        #stats
        srv_obj <- survival::Surv(as.numeric(srvall$event_time), srvall$event)
        broom::tidy(survival::coxph(srv_obj ~ Class))
        #km-plot
        data_fit <- survival::survfit(srv_obj ~ Class)
       ggsurv(data_fit, surv.col = c("#525252", "#bdbdbd")) +
            ylim(0.4, 1) +
            theme_pander() + 
            theme(plot.title = element_text(hjust = 0.5),
                  axis.line = element_line(size = 0.1), 
                  axis.ticks = element_line(size = 0.85),
                  legend.position = 'none')
    }) 
}
shinyApp(ui, server)

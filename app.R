#
# Web application for visualisation of in vitro decidualisation timecourse.
##

library(shiny)
source("functions.R")

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Decidualisation Timecourse"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(

            # Choose gene name
            textInput("genename", label = h3("Gene Name"), value = NULL, placeholder="Enter gene name..."),

            #hr(),
            # choose a genomic region
            textInput("genome_coord", label = h3("Genomic region"), value = NULL, placeholder="Enter chromosomal location..."),


            #fluidRow(column(3, verbatimTextOutput("value"))),

            # select Biopsies
            checkboxGroupInput(inputId="which", label="Select samples:",
                               choices = c("Sample 1"="S169","Sample 2"="S170","Sample 3"="S506", "Sample 4"="S508"),
                               selected = c("S169","S170","S506","S508"))#,

        ),

        # Show a plot of the generated distribution
        mainPanel(
            h3(textOutput("sampleChoice", container = span)),
            plotOutput("distPlot")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$sampleChoice <- renderText(input$which, sep=", ")

    output$distPlot <- renderPlot({
        # draw the line plot of RNA expression based on genename/ensemblID and selected biopsies
        plotRNA.FUN(mydata = AH_EL_RNA_ALLREPS, geneName=input$genename, ensemblID=NULL, biopsies = input$which)

    })
}

# Run the application
shinyApp(ui = ui, server = server)

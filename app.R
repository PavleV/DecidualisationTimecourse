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

            # Choose how to deal with gene names and genomic locations
            radioButtons("link", label = h3("Link genes and genomic coordinates"),
                         choices = list("By gene name" = 1, "By genomic location" = 2, "Separate genes and genomic coordinates" = 3),
                         selected = 1),

            # select Biopsies
            checkboxGroupInput(inputId="which", label="Select samples for RNA plot:",
                               choices = c("Sample 1"="S169","Sample 2"="S170","Sample 3"="S506", "Sample 4"="S508"),
                               selected = c("S169","S170","S506","S508")),

            # Choose gene name
            conditionalPanel(
                condition = "input.link != '2'",
            textInput("genename", label = h3("Gene Name"), value = NULL, placeholder="Enter one or more gene names..."),
            ),

            conditionalPanel(
                condition = "input.link == '1'",
                #numeric input of up and downstream
                numericInput("add_up", label = h3("Add upstream bases"), value = 0),
                numericInput("add_up", label = h3("Add downstream bases"), value = 0)
            ),

            # choose a genomic region

            conditionalPanel(
                condition = "input.link != '1'",
                textInput("genome_coord", label = h3("Genomic region"), value = NULL, placeholder="Enter chromosomal location...")

            )


        ),

        # Show the generated plots
        mainPanel(
            h3(textOutput("sampleChoice", container = span)),
            plotOutput("rnaPlot"),
            h3(textOutput("genome_coord", container = span)),
            plotOutput("atacPlot")
        )
    )
)

# Define server logic required to draw plots
server <- function(input, output) {

    #output$sampleChoice <- renderText(input$which, sep=",")

    output$rnaPlot <- renderPlot({
        req(input$genename)
        # split multiple entries into character vector
        geneNames.to.plot <- str_split_1(input$genename, pattern=regex("[:,_-[:space:]]"))
        # draw the line plot of RNA expression based on genename/ensemblID and selected biopsies
        plotRNA.FUN(mydata = AH_EL_RNA_ALLREPS, geneName=geneNames.to.plot, ensemblID=NULL, biopsies = input$which)
    })

    output$genome_coord <- renderText(input$genome_coord, sep=", ")

    output$atacPlot <- renderPlot({
        req(input$genome_coord)
        # draw the line plot based on peaks within the specified genomic coordinates
        plotATAC.FUN(mydata=ATAC_countsmatrix_cleaned, coordinate.key = AllPeaks.granges, coordinates = input$genome_coord)
    })
}

# Run the application
shinyApp(ui = ui, server = server)

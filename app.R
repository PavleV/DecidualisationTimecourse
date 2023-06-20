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
            textInput("genename", label = h3("Gene Name(s)"), value = "SCARA5", placeholder="Enter one or more gene names..."),
            ),

            conditionalPanel(
                condition = "input.link == '1'",
                #numeric input of up and downstream
                numericInput("add_up", label = h3("Add upstream bases"), value = 0),
                numericInput("add_down", label = h3("Add downstream bases"), value = 0)
            ),

            # choose a genomic region

            conditionalPanel(
                condition = "input.link != '1'",
                textInput("genome_coord", label = h3("Genomic region (hg19)"), value = "chr8:27696656-27881112", placeholder="Enter chromosomal location...")

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

        if(!is.null(input$genename) & input$link != '2')
        {
            # split multiple entries into character vector
            geneNames.to.plot <- str_split_1(input$genename, pattern=regex("[:,_-[:space:]]"))
        }
        if(!is.null(input$genome_coord) & input$link == '2')
        {
            genomic_coordinates <- extractCoord(input$genome_coord)
            geneNames.to.plot <- geneKey.ranges$mcols.GeneName[subjectHits(findOverlaps(genomic_coordinates,geneKey.ranges))]
        }

        req(geneNames.to.plot)
        # draw the line plot of RNA expression based on genename/ensemblID and selected biopsies
        plotRNA.FUN(mydata = AH_EL_RNA_ALLREPS, geneName=geneNames.to.plot, ensemblID=NULL, biopsies = input$which)
    })

    output$genome_coord <- renderText(input$genome_coord, sep=", ")


    output$atacPlot <- renderPlot({

            if(!is.null(input$genename) & input$link == '1')
            {
                gene_names <- str_split_1(input$genename, pattern=regex("[:,_-[:space:]]"))
                genomic_coordinates <- extractCoordfromGene(genes=gene_names, geneCoordKey=geneKey.ranges, add.upstream = input$add_up, add.downstream = input$add_down)

            }
            if(!is.null(input$genome_coord) & input$link != '1')
            {
                genomic_coordinates <- extractCoord(input$genome_coord)
            }

        req(genomic_coordinates)
        # draw the line plot based on peaks within the specified genomic coordinates

        plotATAC.FUN(mydata=ATAC_TPMmatrix, coordinate.key = AllPeaks.granges, coordinates = genomic_coordinates)
    })
}

# Run the application
shinyApp(ui = ui, server = server)

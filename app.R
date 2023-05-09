#
# Web application for visualisation of in vitro decidualisation timecourse.
##

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Decidualisation Timecourse"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            
            # Choose gene name
            textInput("genename", label = h3("Gene Name"), value = "Enter gene name..."),
            
            hr(),
            fluidRow(column(3, verbatimTextOutput("value"))),
            
            # select Biopsies
            checkboxGroupInput(inputId="which", label="Select:", 
                               choices = c("Choice 1"="Choice1","Choice 2"="Choice2","Choice 3"="Choice3"),
                               selected = c("Choice1","Choice2","Choice3")),
            hr(),
            sliderInput("bins",
                        "Number of bins:",
                        min = 1,
                        max = 50,
                        value = 30)
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
    
    output$sampleChoice <- renderPrint(print(input$which))

    output$distPlot <- renderPlot({
        # generate bins based on input$bins from ui.R
        x    <- faithful[, 2]
        bins <- seq(min(x), max(x), length.out = input$bins + 1)

        # draw the histogram with the specified number of bins
        hist(x, breaks = bins, col = 'darkgray', border = 'white',
             xlab = 'Waiting time to next eruption (in mins)',
             main = 'Histogram of waiting times')
    })
}

# Run the application 
shinyApp(ui = ui, server = server)

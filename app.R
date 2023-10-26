#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)
library(readxl)
library(openxlsx)
# PG2

# Define UI for application
ui <- fluidPage(

    titlePanel("Analiza danych z qPCR"),

    sidebarLayout(
      sidebarPanel(
        fileInput("plik", label = "Dane:"),
        selectInput("zakladka", "Wybierz zakładkę:",
                    choices = NULL),
        checkboxGroupInput("refer1", "Wybierz gen referencyjny:",
                    choices = NULL),
        checkboxGroupInput("refer2", "Inny gen referencyjny:",
                    choices = NULL),
        selectInput("kontrol","Wybierz próbę kontrolną:",
                    choices = NULL),
      ),
     
      mainPanel(
        tabsetPanel(
        tabPanel("Tabela",
          tableOutput("tabela")
        ),
        tabPanel("Tabela Cq",
                 tableOutput("tabela_Cq")
        ),
        tabPanel("Tabela Eff",
                 tableOutput("tabela_Eff")
        ),
        tabPanel("Wykres",
          plotOutput("wykres_qpcr"))
    )
    )
)
)





# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  wb <- reactiveVal(NULL)
  
  wybor_gen <- reactive ({
    if (!is.null(input$plik)) {
      read_excel(input$plik$datapath, sheet = 1)
    }
  })
  
  observeEvent(input$plik, {
  wb(input$plik$datapath)
  strony <- getSheetNames(wb())
  updateSelectInput(session, "zakladka", choices = strony)
})


  
observeEvent(wybor_gen(), {
  if (!is.null(wybor_gen())) {
    updateCheckboxGroupInput(session, "refer1", choices = colnames(wybor_gen()))
    updateCheckboxGroupInput(session, "refer2", choices = colnames(wybor_gen()))
    updateSelectInput(session, "kontrol", choices = select(wybor_gen(), "group"))
  }
})
  output$tabela <- renderTable({
    if (!is.null(input$zakladka) && !is.null(wb())) {
      dane <- read.xlsx(input$plik$datapath, sheet = input$zakladka)
      dane
    }

  })
  
  

}

# Run the application 
shinyApp(ui = ui, server = server)

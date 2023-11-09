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
library(psych)

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
        selectInput("kontrol","Wybierz próbę kontrolną:",
                    choices = NULL),
      ),
     
      mainPanel(
        tabsetPanel(
          id = "panele",
        tabPanel("Tabela",
          tableOutput("tabela")
         ),
        tabPanel("Wykres",
          plotOutput("wykres_qpcr")
          ),
        tabPanel("Tabela z wynikami qPCR",
                 actionButton("przycisk", "qpcr"),
                 uiOutput("wyniki")
                 )
        
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
    updateCheckboxGroupInput(session, "refer1", choices = colnames(wybor_gen())[-c(1,2)])
    updateSelectInput(session, "kontrol", choices = select(wybor_gen(), "group"))
  }
})
  output$tabela <- renderTable({
    if (!is.null(input$zakladka) && !is.null(wb())) {
      dane <- read.xlsx(input$plik$datapath, sheet = input$zakladka)
      dane
    }

  })
  
  my_reads <- reactive ({
    if (!is.null(input$plik)) {
      read_excel(input$plik$datapath, sheet = 1)
    }
  })
  linreg_effic <- reactive ({
    if (!is.null(input$plik)) {
      read_excel(input$plik$datapath, sheet = 2)
    }
  })
  normalize <- reactive({input$kontrol})
  reference <- reactive({c(input$refer1)})  

  qPCR.expression <- function(my_reads, linreg_effic, normalize, reference) {
    require(dplyr); require(psych); require(data.table)
    my_reads$sample <- factor(my_reads$sample)
    my_reads$group <- factor(my_reads$group)
    genes <- colnames(my_reads[,3:length(my_reads)])
    for (i in 3:ncol(my_reads)) {
      y <- my_reads[,c(1,2,i)]
      x <- linreg_effic %>% filter(gene == genes[i-2])
      y <- cbind(y, effic=x$value)
      g_eff <- y[,4]^y[,3]
      y[ , ncol(y) + 1] <- g_eff
      colnames(y)[ncol(y)] <- "g_eff"
      colnames(y)[3] <- "gene"
      assign(genes[i-2], y)
    }
    
    m1 <- data.table(get(reference()[1])[1])
    for (i in 1:length(reference())) {
      x <- get(reference()[i])[5]
      names(x)[1] <- reference()[i]
      m1 <- cbind(m1, x)
    }
    m2 <- as.matrix(m1[,2:ncol(m1)])
    rownames(m2) <- m1$sample
    geo_ref <- exp(rowMeans(log(m2)))
    geo_ref2 <- cbind(m1[,1], geo_ref)
    
    genes_list <- list()
    for (i in 1:length(genes)) {
      x <- get(genes[i])
      genes_list[[genes[i]]] <- x
    }
    
    names <- c("sample", names(genes_list))
    exp_summary <- data.frame(levels(my_reads$group))
    exp_sd <- data.frame(levels(my_reads$group))
    for (g in genes_list) {
      g <- g %>% left_join(geo_ref2, by = "sample") %>% group_by(sample) %>% 
        mutate(georef_eg = geo_ref/g_eff)
      n <- g %>% group_by(group) %>% 
        summarise(geo_reps = geometric.mean(georef_eg))
      g <- merge(x = g, y = n, by = "group", all.x = T)
      normal <- g %>% filter(group == normalize())
      normal <- normal[1,c(1,8)]
      g <- g %>% mutate(expression = georef_eg/normal$geo_reps)
      e <- g %>% group_by(group) %>% 
        summarise(express_mean = geometric.mean(expression), sd = sd(expression))
      exp_summary <- cbind(exp_summary, e$express_mean)
      exp_sd <- cbind(exp_sd, e$sd)
    }
    colnames(exp_summary) <- names
    colnames(exp_sd) <- names
    results <- list(expression = exp_summary, sd = exp_sd)
    return(results)
  }

  observeEvent(input$przycisk, {output$wyniki <- renderUI({
    if (length(reference()) >= 2) {
      return(tableOutput("qpcr_table"))
    } else {
      return(textOutput("wymogi"))
      }
    })
  })
  
  output$qpcr_table <- renderTable({
    (qPCR.expression(my_reads = my_reads(), linreg_effic = linreg_effic(), normalize = normalize(), reference = reference()))
  }, digits = 12)
  
  output$wymogi <- renderText({
    paste("Wybierz więcej genów referencyjnych!")
  })
  
 

}

# Run the application 
shinyApp(ui = ui, server = server)

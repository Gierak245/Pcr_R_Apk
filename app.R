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
library(plotly)

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
        downloadButton('downloadData', 'Pobierz dane:')
      ),
     
      mainPanel(
        tabsetPanel(
          id = "panele",
        tabPanel("Tabela",
          tableOutput("tabela")
         ),
        tabPanel("Wykres",
          plotlyOutput("wykres_qpcr")
          ),
        tabPanel("Tabela z wynikami qPCR",
                 actionButton("przycisk", "qPCR"),
                 br(),
                 h3(uiOutput("Expression")),
                 uiOutput("wyniki"),
                 h3(uiOutput("Sd")),
                 uiOutput("wyniki2")
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
  }

  qpcr_results <- reactive(qPCR.expression(my_reads = my_reads(), linreg_effic = linreg_effic(), normalize = normalize(), reference = reference()))
  
  observeEvent(input$przycisk, {output$wyniki <- renderUI({
    if (length(reference()) >= 1) {
      return(tableOutput("qpcr_table1",
                         caption = "Expression"))
    } else {
      return(textOutput("wymogi"))
      }
    })
  })
  
wide_gen <- reactive({
  df <- as.data.frame(qpcr_results()[1:2]) |> 
    rename("sample" = "expression.sample") |> select(-sd.sample)
  
  result <- pivot_longer(df, cols = -sample, names_to = c(".value", "Gene"), names_sep =  "\\.")
  return(result)
})

output$wykres_qpcr <- renderPlotly({
  ggplotly(
    ggplot(wide_gen(), aes(x = sample, y = expression, fill = Gene)) +
      geom_col(position = position_dodge(width = 0.8)) +
      geom_errorbar(aes(ymin = expression - sd, ymax = expression + sd), position = position_dodge(width = 0.8),linewidth = 3) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  )
})
  
  
  
  observeEvent(input$przycisk, {output$wyniki2 <- renderUI({
    if (length(reference()) >= 1) {
      return(tableOutput("qpcr_table2"))
    } else {
      return(textOutput("wymogi"))
    }
    })
  })
  
  
  output$qpcr_table1 <- renderTable({
    df <- as.data.frame(qpcr_results()[1])
    colnames(df) <- substring(colnames(df), 12)
    return(df)
  }, digits = 3)
  
  output$qpcr_table2 <- renderTable({
    df <- as.data.frame(qpcr_results()[2])
    colnames(df) <- substring(colnames(df), 4)
    return(df)
  }, digits = 3)
  
  output$wymogi <- renderText({
    paste("Wybierz więcej genów referencyjnych!")
  })
  
  output$expr <- renderText({
    "Expression"
  })
  
  output$sds <-  renderText({
    "SD"
  })
  observeEvent(input$przycisk, {output$Expression <- renderUI({
    if (length(reference()) >= 1) {
      return(textOutput("expr"))
    } else {
      return(NULL)
    }
    })
  })
  observeEvent(input$przycisk, {output$Sd <- renderUI({
    if (length(reference()) >= 1) {
      return(textOutput("sds"))
    } else {
      return(NULL)
    }
  })
  })
  
  output$downloadData <- downloadHandler(
    filename = function() { 
      "Results.txt"
    },
    content = function(file) {
      if (!is.null(wide_gen())) {
        write.csv(wide_gen(), file)
      }
    }
  )
  
}

# Run the application 
shinyApp(ui = ui, server = server)

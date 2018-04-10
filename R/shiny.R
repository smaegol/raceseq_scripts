#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggplot2)
library(dplyr)
library(ggsignif)
library(plotly)
library(shinycssloaders)

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("Raceseq data"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
        helpText("Choose dataset"),
         selectInput(inputId = "trans", label = "transcript",
                     choices = levels(data_processed$transcript),selected="REPORTERL1_overexp"),
      selectInput(inputId = "exp_type", label = "experiment type",
                  choices = levels(data_processed$exp_type),selected="OVR"),
      uiOutput("condition") %>% withSpinner(),
      uiOutput("project") %>% withSpinner(),
      uiOutput("min_pos")%>% withSpinner(),
      uiOutput("max_pos")%>% withSpinner(),
      checkboxInput("facetProjects","separate projects (sequencing runs)",value = FALSE),
      actionButton("go", "Plot")
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
         plotOutput("distPlot",height=500) %>% withSpinner(),
         plotOutput("mappingPos") %>% withSpinner()
      )
   )
)





# Define server logic required to draw a histogram
server <- function(input, output) {
   
  v <- reactiveValues(doPlot = FALSE)
  
  observeEvent(input$go, {
    # 0 will be coerced to FALSE
    # 1+ will be coerced to TRUE
    v$doPlot <- input$go
  })
  
  observeEvent(input$trans, {
    v$doPlot <- FALSE
  })  
  
  observeEvent(input$exp_type, {
    v$doPlot <- FALSE
  })  
  
  observeEvent(input$max_pos, {
    v$doPlot <- FALSE
  })  
  
  observeEvent(input$min_pos, {
    v$doPlot <- FALSE
  })  
  
  observeEvent(input$exp_type, {
    v$doPlot <- FALSE
  })  
  
  observeEvent(input$facetProjects, {
    v$doPlot <- FALSE
  })  

  observeEvent(input$project, {
    v$doPlot <- FALSE
  })
    
   output$distPlot <- renderPlot({
     
     if (v$doPlot == FALSE) return()
     
      # generate bins based on input$bins from ui.R
      transcript = input$trans
      condition = input$conditions
      exp_type = input$exp_type
      projects = input$project
   
      # draw the histogram with the specified number of bins
      if(length(condition)>=3) {
        tt<-analyze_uridylation(data_processed, transcript,exp_type,conditions=condition,facet_projects = input$facetProjects,mapping_position_min = input$min_mapping_pos,mapping_position_max = input$max_mapping_pos,project=projects)
        print(tt$plot)
      } else {
        print("Please select at least 3 conditions")
      }
   })
   output$mappingPos <- renderPlot({
     if (v$doPlot == FALSE) return()
     
     # generate bins based on input$bins from ui.R
     transcript = input$trans
     condition = input$conditions
     exp_type = input$exp_type
     projects=input$project
     # draw the histogram with the specified number of bins
     map_plot<-plot_mapping_positions(data_processed, transcript,exp_type,conditions=condition,facet_projects=TRUE,mapping_position_min = input$min_mapping_pos,mapping_position_max = input$max_mapping_pos,project=projects)
     print(map_plot)
   })
   
   output$condition <- renderUI({
     transcript_sel = input$trans
     exp2 = input$exp_type
     get_conditions <- data_processed %>% filter(transcript == transcript_sel,exp_type == exp2,) %>% ungroup() %>% select(condition) %>% mutate_all(as.character) %>% mutate_all(as.factor)
     conditions2 <- levels(get_conditions$condition) 
     checkboxGroupInput("conditions","available conditions",conditions2)
     })
   
   output$project <- renderUI({
     transcript_sel = input$trans
     exp2 = input$exp_type
     get_projects <- data_processed %>% filter(transcript == transcript_sel,exp_type == exp2,) %>% ungroup() %>% select(project_name) %>% mutate_all(as.character) %>% mutate_all(as.factor)
     projects2 <- levels(get_projects$project_name) 
     checkboxGroupInput("project","sequencing runs",projects2)
   })
   
   output$min_pos <- renderUI({
     transcript_sel = input$trans
     exp2 = input$exp_type
     get_max_mapping_pos <- data_processed %>% filter(transcript == transcript_sel,exp_type == exp2) %>% ungroup() %>% select(mapping_position) %>% max()
     sliderInput("min_mapping_pos","Minimum mapping position",0,get_max_mapping_pos,0,step=1)
   })
   output$max_pos <- renderUI({
     transcript_sel = input$trans
     exp2 = input$exp_type
     min_mapping_pos <- input$min_mapping_pos
     get_max_mapping_pos <- data_processed %>% filter(transcript == transcript_sel,exp_type == exp2) %>% ungroup() %>% select(mapping_position) %>% max()
     sliderInput("max_mapping_pos","Maximum mapping position",min_mapping_pos+10,get_max_mapping_pos,get_max_mapping_pos,step=1)
   })
   
   
}

# Run the application 
shinyApp(ui = ui, server = server)


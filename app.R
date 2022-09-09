#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
# install.packages("shiny")
# install.packages("shinythemes")
# install.packages('scattermore')
# install.packages('shinydashboard')
# install.packages('feather')
# install.packages('promises')
# install.packages('future')
# devtools::install_github("Roche/ggtips")

library(shiny)
library(shinythemes)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(scattermore)
library(shinydashboard)
library(ggtext)
library(feather)
library(promises)
library(future)
library(viridis)
library(shinyWidgets)
library(ggtips) 

# Load data
dataset <- read_feather("~/Desktop/shiny_data.feather")
reduced <- read_feather("~/Desktop/reduced.feather")
mapping <- read_feather("~/Desktop/mapping.feather")
index <- read_feather("~/Desktop/index.feather")
dot <- read_feather("~/Desktop/dot.feather")



#-----------------------------------------------------------------------------

col_list <- c("Gene Name", "Gene ID", "Counts", "Cluster", "Annotation", "Treatment")

# Load color palette
dakota <-c("#cd4c42", "#5c8492", "#b25757", "#fe906a", "#6f636b", "#6a9491", "#82ac92", "#a26f6a", "#184459", "#596c7f", "#d97f64", "#263946", "#bebab6", "#7a7f84", "#cab6b2", "#fae2af", "#f3933b","#65838d", "#82aca7", "#a0b4ac", "#b5b9b0", "#fbc1c1", "#e89690", "#d76660", "#cac6b9", "#878787", "#cb8034", "#7f93a2", "#ac8287", "#c1d6d3")

clusters <- read.csv("~/GitHub/BmSC_shiny/fig2a_newlabels.csv")


# Define UI for application
ui <- fluidPage(theme = shinytheme("cerulean"),
                navbarPage("B. malayi SC",
                           
 # About Page            
        tabPanel("About",
                 mainPanel(title = "About",
                           p(h2("Welcome")),
                           p("Single-cell atlas of B. malayi microfilariae. Come and explore"),
                           p(tags$a(href="https://www.biorxiv.org/content/10.1101/2022.08.30.505865v1","bioRxiv")),
                           p(tags$a(href="https://github.com/zamanianlab/Bmsinglecell-ms","Github")))),
                  
 
 #Data Table
        tabPanel("Dataset Table",
                  sidebarLayout(
                    sidebarPanel(
                      checkboxGroupInput("show_vars", "Columns:", names(reduced), selected=col_list),
                      width = 2),
                  mainPanel(title = "BmaSC Table", value = "reduced", DT::dataTableOutput("table")))),
 
 
 #UMAP Exploration  
       tabPanel("UMAP Exploration",
                sidebarLayout(
                  sidebarPanel(
                     searchInput(inputId = "Search", label = NULL, placeholder = "Bma-myo-3", btnSearch = icon("search")),
                     width = 3),
              mainPanel("Exploring the single-cell atlas by searching a WBGene ID or gene name.",
                        plotOutput(outputId = "global_umap"), 
                        plotOutput(outputId = "dotplot"),
                        p(""),
                        p("MS=Muscle, MD=Mesoderm, C=Coelomocyte, S=Secretory, CA=Canal-associated, IB=Inner body"),
                        width = 9))),

 # Treatment Comparison
 tabPanel("Treatment",
     sidebarLayout(
       sidebarPanel(
         selectInput(inputId = "Condition",
             label = "Treatment",
             choices = c("Untreated", "Ivermectin (1µM)"),
             selected = NULL)),
      mainPanel("Under Construction"))),
 
  # Downloads
        tabPanel("Downloads",
                 mainPanel("Under Construction"))
   
))






# Define server logic
server <- function(input, output, session) {
  
  # render datatable
  output$table = DT::renderDataTable(reduced)
  

   # # render reactive umap
   #  filteredGeneIDs <- reactive({
   #    mapping %>%
   #    {if(input$Search != "") filter(Gene_ID == input$Search)
   #    else .}
   #  })

  
  # Custom content for the tooltip label of the UMAP
      # customContentFunction <- function(mapping) {
      #   Annotation <- as.character(mapping$Annotation)
      #   Cluster <- as.character(mapping$Cluster)
      #   paste(Annotation, Cluster)
      # }
     
  # make umap based on filtered data
      output$global_umap <- renderPlot({
      plotdata <- subset(mapping, mapping$`Gene ID` == input$Search | mapping$`Gene Name` == input$Search) 
      
      ggplot(data = plotdata, aes(x = UMAP_1, y = UMAP_2))+
          geom_point(data = index, color = "grey", size = 0.5)+
          geom_point(aes(color = Counts), size = 1)+
          geom_text(data = clusters, aes(x = x, y = y, label = str_wrap(text, width = 8)), size = 4, fontface = "plain")+
          scale_color_viridis(guide = "colourbar")+
          labs(color = "Counts")+
          theme(axis.text = element_blank(),
                axis.title= element_blank(),
                axis.ticks = element_blank(),
                legend.text = element_text(size = 12, face = "plain"),
                legend.title = element_text(size = 12, face = "plain"),
                panel.background = element_blank(),
                axis.line = element_blank(),
                legend.background=element_blank(),
                legend.key = element_blank())
    
    })


      
      
    # dot plot based on filtered input (reactive)
      output$dotplot <- renderPlot({
        dotdata <- subset(dot, dot$gene_id == req(input$Search) | dot$gene_name == req(input$Search)) 
      
        ggplot(data = dotdata, aes(y = id, x = gene_name))+
                 geom_point(aes(size = pct.exp, color = avg.exp.scaled))+
                 scale_size("Proportion (%)", range = c(-1, 5))+
                 scale_color_viridis()+
                 labs(x = "Genes", y = "Cluster", size = "Proportion (%)", color = "Avg. Exp.")+
                 facet_grid(cols = vars(ID), rows = vars(gene_name), space = "free", scales = "free", drop = TRUE)+
                 theme(#text=element_text(family="Helvetica"),
                       panel.background = element_blank(),
                       axis.line = element_line (colour = "black"),
                       legend.background=element_blank(),
                       legend.text = element_text(size = 12),
                       legend.title = element_text(size = 12, vjust = 1),
                       legend.key = element_blank(),
                       axis.text.x = ggplot2::element_text(size = 12, angle = 90, vjust = 0.5),
                       axis.text.y = ggplot2::element_text(size = 12, hjust = 1, face = "italic"),
                       axis.title.x = ggplot2::element_text(size = 12, vjust = -1),
                       axis.title.y = ggplot2::element_text(size = 12), 
                       strip.text.x = element_text(size = 12),
                       strip.text.y = element_blank(),
                       strip.background = element_blank(),
                       panel.spacing.x = unit(0.5, "lines"), 
                       #legend.key.width = unit(0.35, "cm"),
                       #legend.key.height = unit(0.25, "cm"),
                       #legend.key.size = unit(0.25, "cm"), 
                       legend.position = "right",
                       panel.grid = element_line(color = "#ededed", size = 0.1))+
                 coord_flip()
               
      })
      
    
    }
    
    


# Run the application 
shinyApp(ui = ui, server = server)












#scratchpad

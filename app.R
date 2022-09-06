#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
install.packages("shiny")
install.packages("shinythemes")
install.packages('scattermore')
install.packages('shinydashboard')
install.packages('feather')
install.packages('promises')
install.packages('future')

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


# Load data
#dataset <- readRDS("~/GitHub/BmaSC_shiny/shiny_data.RDS")
dataset <- read_feather("~/Desktop/shiny_data.feather")
reduced <- read_feather("~/Desktop/reduced.feather")
mapping <- read_feather("~/Desktop/mapping.feather")




#saveRDS(sc_data, "~/Desktop/shiny_data.RDS")
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
                           p(tags$a(href="https://www.biorxiv.org/content/10.1101/2022.08.30.505865v1",
                                    "bioRxiv")))),
 
 #Data Panel
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
                     sidebarSearchForm(label = "Search Gene_id", "Search","searchgene"),
                     width = 3),
              mainPanel("Under Construction", value = "mapping", plotOutput(outputId = "global_umap")))),


 
  # Downloads
        tabPanel("Downloads",
                 mainPanel("Under Construction"))
                
    
   
))






# Define server logic
server <- function(input, output) {
  
  # render datatable for plot
  output$table = DT::renderDataTable(reduced)
  
  
  output$global_umap <- renderPlot({
    ggplot(data = mapping, aes(x = UMAP_1, y = UMAP_2))+
      geom_scattermore(aes(color = Cluster), size = 0.5, show.legend = FALSE)+
      geom_text(data = clusters, aes(x = x, y = y, label = str_wrap(text, width = 8)), size = 3, fontface = "plain")+
      scale_size_area(max_size = 15)+
      scale_color_manual(values = c("#c1d6d3", "#5c8492", "#b25757", "#6a9491", "#7a7f84", "#cab6b2", "#fae2af", "#f3933b","#ac8287", "#65838d", "#82aca7", "#fe906a", "#e3e2e1", "#e89690","#cd4c42", "#6f636b", "#82ac92", "#a26f6a", "#184459", "#596c7f","#263946", "#d97f64", "#a0b4ac", "#e3e2e1", "#fbc1c1", "#7f93a2", "#d76660", "#cac6b9", "#e3e2e1", "#cb8034"), labels = function(color) str_wrap(color, width = 8))+
      labs( color = "Cell Type")+
      theme(text=element_text(family="Helvetica"),
            axis.text = element_blank(),
            axis.title= element_blank(),
            axis.ticks = element_blank(),
            legend.text = element_markdown(size = 10, face = "plain"),
            legend.key.size = unit(0.2, "cm"),
            panel.background = element_blank(),
            legend.margin = margin(0, 0, 0, 0.5, "cm"),
            axis.line = element_blank(),
            legend.background=element_blank(),
            legend.key = element_blank(),
            plot.margin = margin(0.5, 0.5, 0.5, 0.25, unit = "cm"))+
      guides(color = guide_legend(override.aes = list(size=3), ncol = 1))+
      NULL
    })
  
  
  
  

  # render interactive umap
    # interactive search
      #make a data frame that changes when the user searches a WBGeneID
      # filteredGeneIDs = reactive({
      #   if(input$Search == "") {
      #     df = dataset
      #   }
      #   else {
      #     df = dataset %>% dplyr::filter(gene_id == toupper(input$Search))
      #   }
      #   df
      # })
      # 
      # 
      # filteredGeneIDs = reactive({ #make umap based on filtered data
      #   ggplot(data = filteredGeneIDs(), aes(x = UMAP_1, y = UMAP_2))+
      #     geom_scattermore(aes(color = cluster), size = 0.1, show.legend = FALSE)+
      #     geom_text(data = clusters, aes(x = x, y = y, label = str_wrap(text, width = 8)), size = 3, fontface = "plain")+
      #     scale_size_area(max_size = 15)+
      #     scale_color_manual(values = c("#c1d6d3", "#5c8492", "#b25757", "#6a9491", "#7a7f84", "#cab6b2", "#fae2af", "#f3933b","#ac8287", "#65838d", "#82aca7", "#fe906a", "#e3e2e1", "#e89690","#cd4c42", "#6f636b", "#82ac92", "#a26f6a", "#184459", "#596c7f","#263946", "#d97f64", "#a0b4ac", "#e3e2e1", "#fbc1c1", "#7f93a2", "#d76660", "#cac6b9", "#e3e2e1", "#cb8034"), labels = function(color) str_wrap(color, width = 8))+
      #     labs( color = "Cell Type")+
      #     theme(text=element_text(family="Helvetica"),
      #           axis.text = element_blank(),
      #           axis.title= element_blank(),
      #           legend.text = element_markdown(size = 10, face = "plain"),
      #           legend.key.size = unit(0.2, "cm"),
      #           panel.background = element_blank(),
      #           legend.margin = margin(0, 0, 0, 0.5, "cm"),
      #           axis.line = element_blank(),
      #           legend.background=element_blank(),
      #           legend.key = element_blank(),
      #           plot.margin = margin(0.5, 0.5, 0.5, 0.25, unit = "cm"))+
      #     guides(color = guide_legend(override.aes = list(size=3), ncol = 1))+
      #     NULL
      #   })


      # output$global_umap <- renderPlot({ #render the plot with an error if there is no data
      #   validate(
      #     need(nrow(filteredGeneIDs()) > 0, "GeneID not found"))
      #   })
      

    }
    
    

# Run the application 
shinyApp(ui = ui, server = server)

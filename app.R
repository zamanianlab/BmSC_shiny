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

library(shiny)
library(shinythemes)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(scattermore)
library(shinydashboard)
library(ggtext)
library(fst)

# Load data
#dataset <- readRDS("~/GitHub/BmaSC_shiny/shiny_data.RDS")
dataset<- read_fst("~/Desktop/tmp.fst", from = 1, to = NULL, as.data.table = TRUE)
#-----------------------------------------------------------------------------
# new_combined <- readRDS("~/Library/CloudStorage/Box-Box/ZamanianLab/LabMembers/Clair/project_singlecell/BmSinglecell-ms/2_10XGenomics/scmulti_integrated.RDS")
# data <- as_tibble(new_combined@reductions$umap@cell.embeddings, rownames = 'index') # UMAP coordinates for each cell
# 
# md <- as_tibble(new_combined@meta.data, rownames = 'index') # metadata detailing ut/t identity and cluster information
# 
# counts <- as_tibble(new_combined@assays[["RNA"]]@data, rownames = "gene_id") %>%  # gene expression matrix of normalized counts
#   pivot_longer(!gene_id, names_to = "index", values_to = "counts") 
# 
# 
# sc_data <- data %>% 
#   left_join(md) %>% 
#   left_join(counts) %>%
#   filter(counts > 0)
# 
# 
# sc_data <- sc_data %>% 
#   mutate(annotation = case_when(
#     integrated_snn_res.0.5 == "0" ~ "Unannotated",
#     integrated_snn_res.0.5 == "1" ~ "Body Wall Muscle",
#     integrated_snn_res.0.5 == "2" ~ "Unannotated",
#     integrated_snn_res.0.5 == "3" ~ "Unannotated",
#     integrated_snn_res.0.5 == "4" ~ "Unannotated",
#     integrated_snn_res.0.5 == "5" ~ "Coelomocyte",
#     integrated_snn_res.0.5 == "6" ~ "Unannotated",
#     integrated_snn_res.0.5 == "7" ~ "Unannotated",
#     integrated_snn_res.0.5 == "8" ~ "Mesoderm",
#     integrated_snn_res.0.5 == "9" ~ "Unannotated",
#     integrated_snn_res.0.5 == "10" ~ "Motor Neuron",
#     integrated_snn_res.0.5 == "11" ~ "Neuron",
#     integrated_snn_res.0.5 == "12" ~ "Neuron",
#     integrated_snn_res.0.5 == "13" ~ "Canal-associated",
#     integrated_snn_res.0.5 == "14" ~ "Secretory",
#     integrated_snn_res.0.5 == "15" ~ "Unannotated",
#     integrated_snn_res.0.5 == "16" ~ "Mesoderm",
#     integrated_snn_res.0.5 == "17" ~ "Neuron",
#     integrated_snn_res.0.5 == "18" ~ "Body Wall Muscle",
#     integrated_snn_res.0.5 == "19" ~ "Unannotated",
#     integrated_snn_res.0.5 == "20" ~ "Unannotated",
#     integrated_snn_res.0.5 == "21" ~ "Inner Body",
#     integrated_snn_res.0.5 == "22" ~ "Neuron",
#     integrated_snn_res.0.5 == "23" ~ "Neuron",
#     integrated_snn_res.0.5 == "24" ~ "Neuron",
#     integrated_snn_res.0.5 == "25" ~ "Neuron",
#     integrated_snn_res.0.5 == "26" ~ "Neuron"
#   )) %>% 
#   mutate(cluster = case_when(
#     integrated_snn_res.0.5 == "0" ~ 1,
#     integrated_snn_res.0.5 == "1" ~ 2,
#     integrated_snn_res.0.5 == "2" ~ 3,
#     integrated_snn_res.0.5 == "3" ~ 4,
#     integrated_snn_res.0.5 == "4" ~ 5,
#     integrated_snn_res.0.5 == "5" ~ 6,
#     integrated_snn_res.0.5 == "6" ~ 7,
#     integrated_snn_res.0.5 == "7" ~ 8,
#     integrated_snn_res.0.5 == "8" ~ 9,
#     integrated_snn_res.0.5 == "9" ~ 10,
#     integrated_snn_res.0.5 == "10" ~ 11,
#     integrated_snn_res.0.5 == "11" ~ 12,
#     integrated_snn_res.0.5 == "12" ~ 13,
#     integrated_snn_res.0.5 == "13" ~ 14,
#     integrated_snn_res.0.5 == "14" ~ 15,
#     integrated_snn_res.0.5 == "15" ~ 16,
#     integrated_snn_res.0.5 == "16" ~ 17,
#     integrated_snn_res.0.5 == "17" ~ 18,
#     integrated_snn_res.0.5 == "18" ~ 19,
#     integrated_snn_res.0.5 == "19" ~ 20,
#     integrated_snn_res.0.5 == "20" ~ 22,
#     integrated_snn_res.0.5 == "21" ~ 22,
#     integrated_snn_res.0.5 == "22" ~ 23,
#     integrated_snn_res.0.5 == "23" ~ 24,
#     integrated_snn_res.0.5 == "24" ~ 25,
#     integrated_snn_res.0.5 == "25" ~ 26,
#     integrated_snn_res.0.5 == "26" ~ 27
#   ))
# 
# sc_data <- sc_data %>% 
#   select(-"percent.mt", -"Percent.Largest.Gene", -"integrated_snn_res.0.5")
# col_order <- c("gene_id", "counts", "cluster", "annotation", "orig.ident", "index", "nCount_RNA", "nFeature_RNA", "UMAP_1", "UMAP_2")
#dataset <- dataset[,col_order]


# Round the values in columns
dataset$counts <- round(dataset$counts, digits = 3)
dataset$nCount_RNA <- round(dataset$nCount_RNA, digits = 3)

#saveRDS(sc_data, "~/Desktop/shiny_data.RDS")
#-----------------------------------------------------------------------------

col_list <- c("gene_id", "counts", "cluster", "annotation", "orig.ident", "index", "nCount_RNA", "nFeature_RNA", "UMAP_1", "UMAP_2")

# Load color palette
dakota <-c("#cd4c42", "#5c8492", "#b25757", "#fe906a", "#6f636b", "#6a9491", "#82ac92", "#a26f6a", "#184459", "#596c7f", "#d97f64", "#263946", "#bebab6", "#7a7f84", "#cab6b2", "#fae2af", "#f3933b","#65838d", "#82aca7", "#a0b4ac", "#b5b9b0", "#fbc1c1", "#e89690", "#d76660", "#cac6b9", "#878787", "#cb8034", "#7f93a2", "#ac8287", "#c1d6d3")

clusters <- read.csv("~/GitHub/BmaSC_shiny/fig2a_newlabels.csv")


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
                      checkboxGroupInput("show_vars", "Columns:", names(dataset), selected=col_list)),
                  mainPanel(title = "BmaSC Table", value = "dataset", DT::dataTableOutput("table")))),
 
 
 
 #UMAP Exploration  
       tabPanel("UMAP Exploration",
                sidebarLayout(
                  sidebarPanel(
                     sidebarSearchForm(label = "Search Gene_id", "Search","searchgene")),
              mainPanel("Under Construction", value = "dataset", plotOutput(outputId = "global_umap")))),


 
  # Downloads
        tabPanel("Downloads",
                 mainPanel("Under Construction"))
                
    
   
))






# Define server logic
server <- function(input, output) {
  
  # render datatable for plot
  output$table = DT::renderDataTable(dataset)
  
  
  output$global_umap <- renderPlot({
    ggplot(data = dataset, aes(x = UMAP_1, y = UMAP_2))+
      geom_scattermore(aes(color = cluster), size = 0.5, show.legend = FALSE)+
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

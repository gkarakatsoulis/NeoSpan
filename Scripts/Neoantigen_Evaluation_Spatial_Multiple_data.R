# Load necessary libraries
library(shiny)
library(shinydashboard)
library(shinycssloaders)
library(shinyjs)
library(DT)
library(ggplot2)
library(rstatix)
library(dplyr)
library(Seurat)
library(SPOTlight)
library(BayesSpace)
library(spdep)
library(mistyR)
library(plotly)
library(forestplot)
library(nicheDE)
library(SingleCellExperiment)
library(scran)
library(NMF)
# --------------------------------------------------------------------
#                             User Interface
# --------------------------------------------------------------------
options(shiny.maxRequestSize = 12 * 1024^3) # To accept large file sizes

ui <- dashboardPage(
  
  
  dashboardHeader(title = "NeoSpan"),
  
  dashboardSidebar(
    sidebarMenu(
      id = 'tabs',
      menuItem("About", tabName = "about", icon = icon("home")),
      menuItem("Import Data", tabName = "data", icon = icon("database")),
      menuItem("Preprocessing", tabName = "preprocessing", icon = icon("filter")),
      menuItem("Spatial Visualization", tabName = "spatial_visual", icon = icon("map")),
      menuItem("Deconvolution", tabName = "deconv", icon = icon("project-diagram")),
      menuItem("Differential Expression Analysis", tabName = "dea", icon = icon("chart-bar")),
      menuItem("Spatial Autocorrelation", tabName = "spatial_autocorrelation", icon = icon("globe-americas")),
      # menuItem("Spatial Context", tabName = "spatial_context", icon = icon("globe-americas")),
      menuItem("Evidence Synthesis", tabName = "ev_synth", icon = icon("book"))
      
    )
  ),
  
  dashboardBody(
    
    tags$head(
      tags$link(rel = "stylesheet",
                href = "https://cdnjs.cloudflare.com/ajax/libs/animate.css/4.1.1/animate.min.css"),
      tags$style(HTML("
      
      .sleek-grid {
      display: flex;
      flex-wrap: wrap;
      justify-content: center;
      gap: 10px;
      }
    
      .sleek-block {
      min-height: 200px;
      flex: 0 0 calc(50% - 5px); /* 2 columns */
      display: flex;
      align-items: center;
      gap: 0px;
      margin-bottom: 5px;
      padding: 5px 5px;
      border-radius: 12px;
      background-color: #6dabe8;
      box-shadow: 0 4px 12px rgba(0,0,0,0.05);
      transition: transform 0.3s ease;
      }
    
      .sleek-block:hover {
        transform: translateY(-3px);
      }
      
      .sleek-half {
        width: 50%;
        padding: 1px;
      }
      
      .sleek-image-container {
        display: flex;
        width: 100%;
        height: 100%;
        padding: 0;
        margin: 0;
        align-items: center;
      }

      .sleek-image {
        width: 200px;
        height: 200px;
        object-fit: contain;
      }
  
      .sleek-text {
        font-size: 18px;
        font-weight: 500;
        text-align: left;
      }
      
      .right-align {
      flex-direction: row-reverse;
      text-align: center;
      }
    
    @media (max-width: 768px) {
      .sleek-block {
        flex: 0 0 100%;
        flex-direction: column !important;
        text-align: center;
      }
      .sleek-half {
      width: 100%;
      
      }
    .sleek-text { 
      justify-content: center;
      }
    }
    
    "
      ))
      
    ),
    

    tabItems(
      
      # ------------------
      # About
      # ------------------
      tabItem(tabName = "about",
              
              tabsetPanel(id = 'about_tabs',
                
                tabPanel('Dashboard',
                         fluidRow(
                           box(
                             # title = 'About',
                             status = 'primary',
                             # solidHeader = TRUE,
                             width = 12, 
                             
                             {
                               
                               sleek_items <- c(
                                 "Complete data preprocessing",
                                 "Spatial visualization of the spots",
                                 "Spot Deconvolution",
                                 "Spatial Co-Localization of Immune Cells and Neoantigens",
                                 "Differential Expression Analysis - Neoantigen-Positive vs Negative spots",
                                 "Evidence synthesis (in case of multiple individuals/samples)"
                               )
                               
                               sleek_divs <- lapply(seq_along(sleek_items), function(i) {
                                 direction_class <- if (i %% 2 == 0) "right-align" else "left-align"
                                 
                                 div(class = paste("sleek-block", direction_class, "animate__animated", 
                                                   ifelse(i %% 2 == 0, "animate__fadeInRight", "animate__fadeInLeft")),
                                     div(
                                       class = "sleek-half sleek-image-wrapper",
                                       img(src = paste0("img", i, ".png"), class = "sleek-image")
                                     ),
                                     div(
                                       class = "sleek-half sleek-text-wrapper",
                                       span(sleek_items[[i]], class = "sleek-text")
                                     )
                                 )
                               })
                               
                               
                               tagList(
                                 div(
                                   
                                   h1("Welcome to the NeoSpan Dashboard", align = "center",
                                      style = "margin-top: -80px; padding-top: 0px;"),
                                   
                                   p(
                                     "A user-friendly environment for neoantigen evaluation using spatial and non-spatial statistics.",
                                     align = "center"
                                   ),
                                   
                                   div(class = "sleek-grid", sleek_divs), 
                                   
                                   style = "text-align:center; font-size: 20px; margin-top: 100px;"
                                   
                                 ),
                                 column(
                                   width = 4,
                                   align = 'center',
                                   actionButton("neospan_project_btn", 'About NeoSpan project',
                                                icon("info-circle"),
                                                style = "color: white; background-color: #6e5494; border-color: #6e5494")
                                 ),
                                 column(
                                   width = 4,
                                   align = 'center',
                                   actionButton("guidebook_btn", 'Show Guidebook',
                                                icon("question-circle"),
                                                style = "color: white; background-color: #6e5494; border-color: #6e5494")
                                 ),
                                 column(
                                   width = 4,
                                   align = 'center',
                                   actionButton("github_btn", "Visit GitHub", 
                                                icon = icon("github"),
                                                style = "color: white; background-color: #6e5494; border-color: #6e5494")
                                 )
                               )
                               
                             }
                             
                           )
                         )
                    ),
                
                tabPanel('Dashboard Guidebook',
                  
                         tabsetPanel(id = 'guide_tabs',
                                     tabPanel('Preprocessing',
                                              uiOutput("PreprocessingContent")
                                     ),
                                     tabPanel('Deconvolution',
                                              uiOutput("DeconvolutionContent")),
                                     tabPanel('Spatial Autocorrelation',
                                              uiOutput("SpatialAutoContent")),
                                     tabPanel('Differential Expression Analysis',
                                              uiOutput('DEAContent')),
                                     tabPanel('Evidence Synthesis',
                                              uiOutput('EvSyntContent'))
                         )
                  
                  
                ),
                
                tabPanel('NeoSpan Project Pipeline',
                         
                         uiOutput('neospan_project')
                         
                    
                         
                  )
                
                
              )

      ),
      
      
      # Data Viewer Tab
      tabItem(tabName = "data",
              fluidRow(
                useShinyjs(),
                box(
                  title = "Import Dataset(s)", status = "primary", solidHeader = TRUE,
                  width = 12,
                  radioButtons('object_upload', label = 'Select the type of your data', choices = c('10x Visium'), inline = T),
                  # Upload a file
                  div(id = 'hide_inputs',
                      fluidRow(
                        
                        conditionalPanel(
                          
                          "condition = input.object_upload == 'Seurat (.rds)'",
                          column(
                            width = 6,
                            fileInput("rds_file", 'Browse your Seurat object file', accept = '.rds')
                          )
                          
                        ),
                        
                        conditionalPanel(
                          
                          "condition = input.object_upload == '10x Visium'",
                          column(
                            width = 3,
                            fileInput("h5_file", 'filtered_feature_bc_matrix.h5', accept = '.h5')
                          ),
                          column(
                            width = 3,
                            fileInput('spatial_folder', 'Spatial folder (.tar file)', accept = c(".zip", ".tar.gz", ".tar", ".gz"),)
                          )
                          
                        ),

                        column(
                          width = 3,
                          uiOutput('neo_antigen_file_reset')
                        ),
                        column(
                          width = 3,
                          uiOutput('region_ann_file_reset')
                        ),
                        column(
                          width = 4,
                          textInput('dataset_name', 'Dataset ID', value = 'Individual_1')
                        ),
                        column(
                          width = 8,
                          div(style = "margin-top: 25px;",
                              actionButton('load_visium', 'Load Data')
                          )
                        )
                        
                        
                      ),

                      
                    ),
                  div(id = 'hide_add_remove_buttons',
                      column(
                        width = 6,
                        uiOutput('add_individual')
                      ),
                      column(
                        width = 6,
                        uiOutput('remove_individual')
                      )
                  ),
                  div(id = 'hide_remove_ind',
                      column(
                        width = 6,
                        uiOutput('remove_individual_sp')
                      ),
                      column(
                        width = 6,
                        div(style = "margin-top: 25px;",
                            uiOutput('remove_individual_confirm')
                        )
                      )
                    )
                ),
                
                uiOutput('box_summary')

                
              ),
              
      ),
      
      # ------------------
      # Preprocessing
      # ------------------
      tabItem(tabName = "preprocessing",
              
              tabsetPanel(
                
                tabPanel('Quality Control',
                         
                         fluidRow(
                           
                           box(
                             title = "Quality Control Visualization", status = "primary", solidHeader = TRUE,
                             width = 9,
                             # uiOutput('qc_vis'),
                             column(
                               width = 6,
                               selectInput('qc_plot_type_data', 'Select the Individual', choices = NULL, multiple = F)
                             ),
                             column(
                               width = 6,
                               selectInput("qc_plot_type", "Select the plot type", choices = c('Violin plots', 'Scatterplot', 'Spatial plot'), multiple = F)
                             ),
                             actionButton("qc_vis_Button", "Run"),
                             plotlyOutput('qc_plot'),
                             column(
                               width = 12,
                               align = 'right',
                               # uiOutput('get_download_qc_plot')
                             )
                           ),
                           
                           box(
                             title = "Quality Control Filters", status = "primary", solidHeader = TRUE,
                             width = 3,
                             uiOutput('qc'),
                             radioButtons("qc_filter_data", "Select the Individual(s)", choices = c('Any', 'Specific'), inline = T),
                             conditionalPanel(
                               condition = "input.qc_filter_data == 'Specific'",
                               selectInput('qc_filter_data_sp', 'Specify the Individual(s)', choices = NULL, multiple = T)
                             ),
                             selectInput("n_genes_filter", "Number of Genes Filtering:",
                                          choices = c("No gene filtering" = "none",
                                                      "Minimum" = "min",
                                                      "Maximum" = "max",
                                                      "Between" = "both"),
                                          selected = "none"),
                             
                             conditionalPanel(
                               condition = "input.n_genes_filter == 'min' || input.n_genes_filter == 'both'",
                               numericInput("min_genes", "Min Genes per Spot:", 
                                            min = 0, max = 10000, value = 500)
                             ),
                             
                             conditionalPanel(
                               condition = "input.n_genes_filter == 'max' || input.n_genes_filter == 'both'",
                               numericInput("max_genes", "Max Genes per Spot:", 
                                            min = 0, max = 10000, value = 700)
                             ),
                             
                             
                             selectInput("n_counts_filter", "Number of Counts Filtering:",
                                         choices = c("No count filtering" = "none",
                                                     "Minimum" = "min",
                                                     "Maximum" = "max",
                                                     "Between" = "both"),
                                         selected = "none"),
                             
                             conditionalPanel(
                               condition = "input.n_counts_filter == 'min' || input.n_counts_filter == 'both'",
                               numericInput("min_counts", "Min Counts per Spot:", 
                                            min = 0, max = 10000, value = 500)
                             ),
                             
                             conditionalPanel(
                               condition = "input.n_counts_filter == 'max' || input.n_counts_filter == 'both'",
                               numericInput("max_counts", "Max Counts per Spot:", 
                                            min = 0, max = 10000, value = 700)
                             ),
                             
                             numericInput("max_mito", "Max % Mitochondrial Genes:", min = 0, max = 100, value = 10),
                             column(
                               width = 6,
                               align = 'left',
                               actionButton("qc_Button", "Run")
                             ),
                             column(
                               width = 6,
                               align = 'right',
                               actionButton("qc_reset_Button", 'Reset')
                             )
                           ),
                         
                         
                         )
                  
                  
                  
                ),
                
                tabPanel('Normalization and Feature Selection (Primary)',
                  
                  fluidRow(
                    
                    box(
                      title = "Data Normalization", 
                      status = "primary", 
                      solidHeader = TRUE,
                      width = 3,
                      column(
                        width = 12,
                        radioButtons("normal_data", "Select the Individual", choices = c('Any', 'Specific'), inline = T)
                      ),
                      column(
                        width = 12,
                        conditionalPanel(
                          condition = "input.normal_data == 'Specific'",
                          selectInput('normal_data_sp', 'Specify the Individual(s)', choices = NULL, multiple = T)
                        )
                      ),
                      column(
                        width = 12,
                        selectInput("normal", "Select the normalization method", 
                                    choices = c('SCTransform','Log-normalization'), 
                                    multiple = FALSE)
                      ),
                      
                      column(
                        width = 12,
                        conditionalPanel(
                          condition = "input.normal == 'Log-normalization'",
                          uiOutput('hvg') 
                        )
                        
                      ),
                      
                      column(
                        width = 12,
                        conditionalPanel(
                          condition = "input.normal == 'Log-normalization'",
                          uiOutput('hvg_selection_method')
                        )
                      ),
                      
                      column(
                        width = 12,  
                        div(style = "margin-top: 25px;",
                            actionButton("normal_Button", "Run")
                        )
                      )
                      
                    ),
                    
                    uiOutput('box_normal_plot')
                    

                    
                    
                  )
                  
                  
                ),
                
                # tabPanel('Feature Selection (Primary)',
                #          
                #          fluidRow(
                #            
                #            box(
                #              title = "Highly Variable Genes", status = "primary", solidHeader = TRUE,
                #              width = 12,
                #              # column(
                #              #   width = 3,
                #              #   radioButtons("hvg_data", "Select the Individual(s)",
                #              #                choices = c('Any', 'Specific'),
                #              #                inline = T)
                #              # ),
                #              # column(
                #              #   width = 3,
                #              #   conditionalPanel(
                #              #     condition = "input.hvg_data == 'Specific'",
                #              #     selectInput('hvg_data_sp', 'Specify the Individual(s)', choices = NULL, multiple = T)
                #              #   )
                #              # ),
                #              # column(
                #              #   width = 3,
                #              #   uiOutput('hvg') 
                #              #   
                #              # ),
                #              # 
                #              # column(
                #              #   width = 3,
                #              #   uiOutput('hvg_selection_method')
                #              #   
                #              # ),
                #              # column(
                #              #   width = 12,
                #              #   div(style = "margin-top: 25px;",
                #              #       actionButton("hvg_Button", "Run")
                #              #   )
                #              # ),
                #              
                #              fluidRow(
                #                
                #                column(
                #                  width = 4,
                #                  align = 'left',
                #                  uiOutput('hvg_plot_data')
                #                ),
                #                
                #                column(
                #                 
                #                  width = 8,
                #                  align = 'left',
                #                  div(style = "margin-top: 25px;",
                #                      uiOutput('hvg_plot_Button')
                #                  )
                #                  
                #                  
                #                ),
                #                
                #                column(
                #                  
                #                  width = 12,
                #                  align = 'center',
                #                  plotOutput('hvg_plot')
                #                  
                #                ),
                #                
                #                column(
                #                  
                #                  width = 12,
                #                  align = 'right',
                #                  uiOutput('get_download_hvg_plot')
                #                  
                #                )
                #                
                #                
                #              )
                #            )
                #            
                #          )
                #          
                #       ),
                
                tabPanel('Dimensionality Reduction',
                         
                         fluidRow(
                           
                           box(
                             title = "Dimensionality Reduction", status = "primary", solidHeader = TRUE,
                             width = 12,
                             
                             fluidRow(
                               
                               column(
                                 width = 4,
                                 selectInput('dim_red_data_sp', 'Specify the Individual', choices = NULL, multiple = F)
                               ),
                               column(
                                 
                                 width = 4,
                                 uiOutput('dim_red')
                                 
                               ),
                               
                               column(
                                 
                                 width = 4,
                                 numericInput('number_dims', 'Number of dimensions:', value = 30, min = 1, max = 100)
                                 
                               ),
                               
                               
                             ),
                             
                             fluidRow(
                               
                               column(
                                 
                                 width = 12,
                                 align = 'left',
                                 actionButton("dim_red_Button", "Run")
                                 
                               )
                               
                             ),
                             
                             
                             fluidRow(
                               
                               box(
                                 
                                 title = 'Dimensionality Reduction Plot', status = "primary", solidHeader = TRUE,
                                 width = 6,
                                 fluidRow(
                                   
                                   column(
                                     
                                     width = 12,
                                     plotOutput('dim_red_plot'),
                                     
                                   ),
                                   
                                   column(
                                     
                                     width = 12,
                                     align = 'right',
                                     uiOutput('get_download_dim_red_plot')
                                     
                                   ),
                                   
                                   
                                 )
                                 
                               ),
                               
                               box(
                                 
                                 title = 'Elbow Plot', status = "primary", solidHeader = TRUE,
                                 width = 6,
                                 
                                 fluidRow(
                                   
                                   column(
                                     
                                     width = 12,
                                     align = 'left',
                                     uiOutput('elbow_Button')
                                     
                                   ),
                                   
                                   column(
                                     
                                     width = 12,
                                     align = 'left',
                                     plotOutput('elbow_plot')
                                     
                                   )
                                   
                                   
                                 )
                                 
                               )
                               
                               
                             )

                             
                           ),
                           
                           
                         )
                         
                         
                      ),
                
                tabPanel('Clustering',
                         
                         fluidRow(
                           
                           box(
                             title = "Nearest-Neighbor Graph and Clustering", status = "primary", solidHeader = TRUE,
                             width = 12,
                             column(
                               width = 4,
                               selectInput('clustering_data_sp', 'Specify the Individual', choices = NULL, multiple = F)
                             ),
                             column(
                               width = 8,
                               numericInput('dimens_num', 'Select number of dimensions', value = 2)
                             ),
                             actionButton("neighbors_Button", "Find Neighbors"),
                             column(
                               width = 12,
                               uiOutput('dimens_num_ui')
                             ),
                             uiOutput('clustering_algorithm'),
                             uiOutput('clustering_Button'),
                             uiOutput('clustering_plot_Button'),
                             plotOutput('clustering_plot')
                           ),
                           
                         )
                         
                         
                      ),
                tabPanel('Preprocessing History',
                         
                         box(
                           title = "Preprocessing Steps Implemented", status = "primary", solidHeader = TRUE,
                           width = 12,
                           column(
                             width = 12,
                             dataTableOutput("preprocess_history_table")
                           )
                         ),
                    )
                         
                
              ),
              
      ),
      
      # ------------------
      # Deconvolution
      # ------------------
      
      tabItem(tabName = 'deconv',
              
              fluidRow(
                box(
                  title = "Deconvolution", status = "primary", solidHeader = TRUE,
                  width = 12,
                  
                  radioButtons("deconv_data", "Select the Individual(s)", choices = c('Any', 'Specific'), inline = T),
                  
                  conditionalPanel(
                    condition = "input.deconv_data == 'Specific'",
                    selectInput('deconv_data_sp', 'Specify the Individual(s)', choices = NULL, multiple = T)
                  ),
                  
                  radioButtons("deconv_method", "Choose Method:",
                               choices = c("SPOTlight (Cell Type Deconvolution)" = "spotlight",
                                           "BayesSpace (Spatial Clustering)" = "bayespace")),
                  
                  # Conditional panel for SPOTlight (only shown if SPOTlight is selected)
                  conditionalPanel(
                    condition = "input.deconv_method == 'spotlight'",
                    fluidRow(
                      column(width = 6,
                             fileInput("scrna_upload", "Upload scRNA-seq Reference (RDS/Seurat)", 
                                       buttonLabel = "Browse...")),
                      
                      column(width = 4,
                             textInput("celltype_column", "Cell Type Column Name", 
                                       value = "cell_type")),
                      column(
                        width = 2,
                        div(style = "margin-top: 25px;",
                            actionButton("scrna_upload_Button", "Upload")
                        )
                      )
                      
                    )

                   
                  ),
                  

                  
                  # BayesSpace parameters (only shown if BayesSpace is selected)
                  conditionalPanel(
                    condition = "input.deconv_method == 'bayespace'",
                    sliderInput("n_clusters_deconv", "Number of Clusters:", 
                                min = 2, max = 10, value = 5),
                    selectInput("bayespace_model", "Model:",
                                choices = c("t", "normal"), selected = "t")
                  ),
                  
                  uiOutput('deconv_Button_ui'),
                  hr(),
                  # conditionalPanel(
                  #   condition = "input.deconv_method == 'spotlight'",
                  #   uiOutput('results_validate_ui'),
                  #   uiOutput('results_report_ui')
                  # ),
                  uiOutput('deconv_data_sp_show_ui'),
                  uiOutput('deconv_show_Button_ui'),
                  plotOutput("deconv_plot"),
                  dataTableOutput("deconv_table"),
                  # uiOutput("download_deconv_ui")
                )
                
              )
                
      ),
      
      # ------------------
      # Spatial Visualization
      # ------------------
      
      tabItem(tabName = 'spatial_visual',
              
              fluidRow(
                box(
                  title = "Spatial Visualization of data", status = "primary", solidHeader = TRUE,
                  width = 12,
                  column(
                    
                    width = 4,
                    align = 'center',
                    selectInput('spatial_vis_data_ID','Select the Individual',
                                choices = NULL, multiple = F)
                    
                  ),
                  column(
                    
                    width = 4,
                    align = 'center',
                    uiOutput('spatial_vis_group')
                    
                  ),
                  
                  column(
                    
                    width = 4,
                    align = 'center',
                    uiOutput('spatial_vis_data')
                    
                  ),
                  
                  # conditionalPanel(
                  #   condition = "input.spatial_vis_group == 'Neoantigen Status'",
                  #   radioButtons("neoantigens_specific", "Specific Neoantigens", choices = c('Any', 'Specific'))
                  # ),
                  # 
                  # conditionalPanel(
                  #   condition = "input.spatial_vis_group == 'Neoantigen Status' & input.neoantigens_specific == 'Specific'",
                  #   selectInput("neoantigens_define", "Select the desired Neoantigens", choices = c('N1', 'N2'))
                  # ),
                  div(style = "margin-top: 25px;",  # Adjust this value as needed
                      actionButton("spatial_vis_Button", "Run")
                  ),
                  column(
                    width = 12,
                    align = 'center',
                    plotOutput('spatial_vis_plot',
                               width = "800px", 
                               height = "640px")
                  ),
                  column(
                    width = 12,
                    align = 'right',
                    uiOutput("download_spatial_plot_ui")
                    
                  )
                )
                
              )
              
      ),
      
      # ------------------
      # Differential Expression Analysis
      # ------------------
      tabItem(tabName = "dea",
              fluidRow(
                uiOutput('dea_appear'),
                uiOutput('dea_results_box')
                
              )
      ),
      
      # ------------------
      # Spatial Autocorrelation
      # ------------------
      tabItem(tabName = "spatial_autocorrelation",
              fluidRow(
                box(
                  title = "Spatial Autocorrelation", status = "primary", solidHeader = TRUE,
                  width = 12,
                  column(
                    width = 4,
                    selectInput('spatial_auto_data', 'Select the Individual*',
                                choices = NULL, multiple = F)
                  ),
                  
                  column(
                    
                    width = 4,
                    radioButtons('morans_genes_specific', 'Desired Genes**',
                                 choices = c('Any', 'Specific'), inline = T)
                    
                  ),
                  
                  column(
                    width = 4,
                    conditionalPanel(
                      
                      condition = "input.morans_genes_specific == 'Specific'",
                      selectizeInput("morans_genes_define", 
                                     "Select the desired Gene(s)",
                                     choices = NULL,
                                     multiple = T,
                                     options = list(
                                       placeholder = 'Type to search genes...',
                                       maxOptions = 20, 
                                       server = TRUE    
                                     ))
                      
                    )
                    
                  ),
                  
                  column(width = 12,
                         textOutput('spatial_auto_preproc_message')),
                  column(width = 12,
                         textOutput('spatial_auto_genes_message')),
                  
                  column(
                    
                    width = 12,
                    actionButton("morans_I", "Run")
                    
                  ),
                  

                  box(
                    
                    title = "Moran's I table", status = "primary", solidHeader = TRUE,
                    width = 6,
                    uiOutput('moran_table_Button'),
                    withSpinner(DTOutput("moran_table"), type = 6)
                    
                    
                  ),
                  
                  
                  box(
                    
                    title = "Moran's I plot", status = "primary", solidHeader = TRUE,
                    width = 6,
                    uiOutput('moran_plot_Button'),
                    plotlyOutput("moran_plot"),
                    uiOutput('get_download_moran_plot')
                    
                    
                  ),
                  
                  box(
                    
                    title = "Spatial plot", status = "primary", solidHeader = TRUE,
                    width = 12,
                    uiOutput('spatial_plot_gene'),
                    uiOutput('spatial_plot_Button'),
                    column(
                      width = 12,
                      align = 'center',
                      plotOutput("spatial_plot",
                                 width = "800px",
                                 height = "640px")
                    )
                    
                    
                  ),
                    
                )
                
              )
      ),
      # ------------------
      # Spatial Context
      # ------------------
      tabItem(tabName = "spatial_context",
              fluidRow(
                box(
                  title = "Spatial Modelling", status = "primary", solidHeader = TRUE,
                  width = 12,
                  selectInput('spatial_context_ind_sp',
                              'Select individual',
                              choice = NULL, multiple = F),
                  actionButton('spatial_context_Button', 'Run')
                  
                )
                
              )
      ),
      # ------------------
      # Evidence Synthesis
      # ------------------
      tabItem(tabName = "ev_synth",
              fluidRow(
                box(
                  title = "Evidence Synthesis", status = "primary", solidHeader = TRUE,
                  width = 12,
                  column(
                    width = 12,
                    verbatimTextOutput('ev_synth_summary')
                  ),
                  column(
                    width = 4,
                    selectizeInput('ev_synth_gene', 'Select the Gene',
                                choices = NULL, multiple = T)
                  ),
                  column(
                    width = 4,
                    uiOutput('ev_synth_data_ui')
                  ),
                  column(
                    width = 4,
                    conditionalPanel(
                      condition = "input.ev_synth_data == 'Specific'",
                      selectInput('ev_synth_data_sp', 'Specify the Individual(s)', choices = NULL, multiple = T)
                    )
                  ),
                  column(
                    width = 12,
                    uiOutput('ev_synth_forest_Button_ui')
                  ),
                  column(
                    width = 12,
                    plotOutput('ev_synth_forest_plot')
                  ),
                  column(
                    width = 12,
                    align = 'right',
                    uiOutput('get_download_ev_synth_forestplot')
                  )
                )
              )
      )
    )
  )
)

# Server Definition
server <- function(input, output, session) {
  
  # --------------------------------------------------------------
  #                           HOME PAGE
  # --------------------------------------------------------------
  observeEvent(input$guidebook_btn, {
    updateTabsetPanel(session, "about_tabs", selected = "Dashboard Guidebook")
    updateTabsetPanel(session, "guide_tabs", selected = "Preprocessing")
  })
  
  observeEvent(input$neospan_project_btn, {
    updateTabsetPanel(session, "about_tabs", selected = "NeoSpan Project Pipeline")
  })
  
  observeEvent(input$github_btn, {
    browseURL("https://github.com/gkarakatsoulis/NeoSpan")
  })
  
  
  output$neospan_project <- renderUI({
    tags$iframe(
      src = "Neospan_Project.html",
      width = "100%",
      height = "950px",
      style = "border: none;"
    )
    
    # includeHTML("www/Neospan_Project.html")
  })
  
  output$SpatialAutoContent <- renderUI({
    tags$iframe(
      src = "Spatial_Autocorrelation.html",
      width = "100%",
      height = "950px",
      style = "border: none"
    )
    # includeHTML("www/Spatial_Autocorrelation.html")
  })
  
  output$DeconvolutionContent <- renderUI({
    # includeHTML("www/Deconvolution.html")
    tags$iframe(
      src = "Deconvolution.html",
      width = "100%",
      height = "950px",
      style = "border: none"
    )
  })
  
  output$DEAContent <- renderUI({
    # includeHTML("www/DEA.html")
    tags$iframe(
      src = "DEA.html",
      width = "100%",
      height = "950px",
      style = "border: none"
    )
  })
  
  output$PreprocessingContent <- renderUI({
    # includeHTML("www/Preprocessing.html")
    tags$iframe(
      src = "Preprocessing.html",
      width = "100%",
      height = "950px",
      style = "border: none;"
    )
    
  })
  
  

  
  
  
  # --------------------------------------------------------------
  #                   DATASETS AS LISTS (PER INDIVIDUAL)
  # --------------------------------------------------------------
  all_visium = reactiveVal(list())
  all_original = reactiveVal(list())
  all_neo_status_exists = reactiveVal(list())
  all_region_ann_exists = reactiveVal(list())
  
  # --------------------------------------------------------------
  #                   IMPORT DATA
  # --------------------------------------------------------------
  
  # Dataset ID defeault
  n_of_inds = reactiveVal(1)
  dataset_id_defeault = reactiveVal(NULL)
  individual_ids = reactiveVal(NULL)

  neo_file_input_id <- reactiveVal("dataset_neoantigen_status")
  region_file_input_id <- reactiveVal("dataset_region_ann")
  
  output$neo_antigen_file_reset <- renderUI({
    fileInput(neo_file_input_id(), 'Neoantigen Status (optional)', accept = c('.csv', '.txt', '.xlsx'))
  })
  
  output$region_ann_file_reset <- renderUI({
    fileInput(region_file_input_id(), 'Region annotation (optional)', accept = c('.csv', '.txt', '.xlsx'))
  })
  
  outputOptions(output, "neo_antigen_file_reset", suspendWhenHidden = FALSE)
  outputOptions(output, "region_ann_file_reset", suspendWhenHidden = FALSE)
  
  
  observeEvent(input$add_individual_Button, {
    
    shinyjs::hide("hide_remove_ind")
    shinyjs::show("hide_inputs")
    # shinyjs::hide("hide_add_remove_buttons")
    
    enable('h5_file')
    enable("spatial_folder")
    enable("dataset_neoantigen_status")
    enable("dataset_region_ann")
    
    individual_number = n_of_inds() + 1
    
    dataset_id_defeault(paste0('Individual_', individual_number))
    updateTextInput(session, 'dataset_name', 'Dataset ID',
                    value = dataset_id_defeault())
    
    n_of_inds(individual_number)
    
    new_id_file_input = paste0('dataset_neoantigen_status_', individual_number)
    neo_file_input_id(new_id_file_input)
    
    new_id_file_input_reg = paste0('dataset_region_ann_', individual_number)
    region_file_input_id(new_id_file_input_reg)
    
    
    
  })
  
  observeEvent(input$remove_individual_Button, {
    
    shinyjs::show("hide_remove_ind")
    
    output$remove_individual_sp = renderUI({
      
      selectInput('remove_individual_sp',
                  'Select Individual(s) to remove',
                  choices = individual_ids(),
                  multiple = T)
      
    })
    
    output$remove_individual_confirm = renderUI({
      
      
      actionButton('remove_individual_confirm', 'Remove')
      
    })
    
    
  })
  
  observeEvent(input$remove_individual_confirm, {
    
    req(input$remove_individual_sp)
    
    myids = names(all_visium())
    tmp_all_visium = all_visium()
    tmp_all_original = all_original()
    
    for (i in input$remove_individual_sp){
      
      tmp_all_visium[[i]] = NULL
      tmp_all_original[[i]] = NULL
      
      
    }
    
    all_visium(tmp_all_visium)
    all_original(tmp_all_original)
    
    individual_ids(names(all_visium()))
    
    if (length(all_visium()) == 0){
      
      output$box_summary = renderUI({
        
        box(
          title = 'Dataset(s) Summary', status = "primary", solidHeader = TRUE,
          width = 12
          
        )
        
      })
      
    } else {
      
      output$box_summary = renderUI({
        box(
          
          title = 'Dataset(s) Summary', status = "primary", solidHeader = TRUE,
          width = 12,
          column(
            width = 3,
            DTOutput("uploaded_datasets_table") 
          ),
          column(
            width = 4,
            selectInput('data_summary_select', 'Select data to summarise', choices = individual_ids(), multiple = F)
            
          ),
          column(
            width = 1,
            div(style = "margin-top: 25px;",
                actionButton('data_summary_select_Button', 'Show')
            )
            
          ),
          column(
            width = 4,
            verbatimTextOutput("data_visium_summary")
          )
          
        )
        
      })
      
    }
    
    updateSelectInput(session,
                      'qc_plot_type_data', 'Select the Individual',
                      choices = individual_ids())
    
    updateSelectInput(session,
                      'qc_filter_data_sp', 'Specify the Individual(s)',
                      choices = individual_ids())
    
    updateSelectInput(session,
                      'normal_data_sp', 'Specify the Individual(s)',
                      choices = individual_ids())
    
    updateSelectInput(session,
                      'hvg_data_sp', 'Specify the Individual(s)',
                      choices = individual_ids())
    
    updateSelectInput(session,
                      'dim_red_data_sp', 'Specify the Individual(s)',
                      choices = individual_ids())
    
    updateSelectInput(session,
                      'clustering_data_sp', 'Specify the Individual(s)',
                      choices = individual_ids())
    
    updateSelectInput(session,
                      'spatial_vis_data_ID', 'Select the Individual',
                      choices = individual_ids())
    
    updateSelectInput(session,
                      'deconv_data_sp', 'Specify the Individual(s)',
                      choices = individual_ids())
    
    
    updateSelectInput(session,
                      'dea_data_sp', 'Specify the Individual(s)',
                      choices = individual_ids())
    
    updateSelectInput(session,
                      'spatial_auto_data', 'Select the Individual*',
                      choices = individual_ids())
    
    updateSelectInput(session,
                      'spatial_context_ind_sp', 'Select individual',
                      choices = individual_ids())
    
    updateSelectInput(session,
                      'ev_synth_data', 'Select the Individual',
                      choices = individual_ids())
    
    
  })
  
  
  
  # Dataset read
  dataset_visium = reactiveVal(NULL)
  original_visium = reactiveVal(NULL) # This will be the original file without changes
  dataset_neo_status = reactiveVal(NULL)
  dataset_region_ann = reactiveVal(NULL)
  

  # 1. Load the data when file is uploaded
  observeEvent(input$load_visium, {
    
    if (input$object_upload == '10x Visium') {
      
      req(is.null(input$h5_file) | is.null(input$spatial_folder))
      
      showNotification('Error: Please upload both the filtered_feature_bc_matrix and the spatial_folder')
      
    }
    
    
  })
  
  
  observeEvent(input$load_visium, {
    
    if (input$object_upload == '10x Visium') {
      
      req(input$h5_file, input$spatial_folder)
      
      
      showModal(modalDialog(
        title = "Please wait",
        "Uploading...",
        footer = NULL,
        easyClose = FALSE
      ))
      
      shinyjs::hide("hide_inputs")
      shinyjs::show("hide_add_remove_buttons")
      
      
      disable('h5_file')
      disable("spatial_folder")
      disable("dataset_neoantigen_status")
      disable("dataset_region_ann")
      
      # Create a temporary directory for processing
      temp_dir <- tempdir()
      
      # Copy .h5 file to temp_dir
      h5_path <- file.path(temp_dir, input$h5_file$name)
      file.copy(input$h5_file$datapath, h5_path)
      
      # Extract spatial.tar to temp_dir/spatial/
      untar(input$spatial_folder$datapath, exdir = temp_dir)
      
      # Load data using Seurat
      tmp_seurat_obj = Load10X_Spatial(
        data.dir = temp_dir,
        filename = input$h5_file$name
      )
      
      
    } else {
      
      req(input$rds_file)
      
      tmp_seurat_obj <- readRDS(input$rds_file)
      
    }
    
    
    
    dataset_visium(tmp_seurat_obj)
    original_visium(tmp_seurat_obj)

    
    # Load neoantigen status
    if (!is.null(input[[neo_file_input_id()]])){
      
        z = read.csv(input[[neo_file_input_id()]]$datapath, sep = ';')


        if (ncol(z) == 3){

          colnames(z) = c('Barcode', 'Neoantigen_Status', 'Neoantigen_Counts')

        } else {

          colnames(z) = c('Barcode', 'Neoantigen_Status')

        }

        dataset_neo_status(z)
        
        
        current_visium <- dataset_visium()
        current_original = original_visium()
        
        
        # Prepare annotation data
        status_data <- z
        rownames(status_data) <- status_data$Barcode
        
        # Add metadata only for barcodes that exist in the Seurat object
        common_barcodes <- intersect(rownames(status_data), colnames(current_visium))
        if(length(common_barcodes) == 0) {
          showNotification("No matching barcodes found between dataset and annotations!",
                           type = "error")
          return()
        }
        
        status_data = status_data[common_barcodes, ]
        
        # Add metadata
        current_visium <- AddMetaData(
          object = current_visium,
          metadata = status_data["Neoantigen_Status"]
        )
        
        current_original <- AddMetaData(
          object = current_original,
          metadata = status_data["Neoantigen_Status"]
        )
        
        # Update the reactive value
        dataset_visium(current_visium)
        original_visium(current_original)
        
        tmp_neo_exists = all_neo_status_exists()
        tmp_neo_exists[[input$dataset_name]] = T
        all_neo_status_exists(tmp_neo_exists)

        
    }
    
    
    
    
    # Load the region annotation info
    if (!is.null(input[[region_file_input_id()]])){
      
      z = read.csv(input[[region_file_input_id()]]$datapath)
      
      colnames(z) = c('Barcode', 'Region')
      
      dataset_region_ann(z)
      
      
      # Get current Seurat object
      current_visium <- dataset_visium()
      current_original = original_visium()
      
      
      # Prepare annotation data
      region_data <- dataset_region_ann()
      rownames(region_data) <- region_data$Barcode
      
      # Add metadata only for barcodes that exist in the Seurat object
      common_barcodes <- intersect(rownames(region_data), colnames(current_visium))
      if(length(common_barcodes) == 0) {
        showNotification("No matching barcodes found between dataset and annotations!",
                         type = "error")
        return()
      }
      
      region_data = region_data[common_barcodes, ]
      
      # Add metadata
      current_visium <- AddMetaData(
        object = current_visium,
        metadata = region_data["Region"]
      )
      
      current_original <- AddMetaData(
        object = current_original,
        metadata = region_data["Region"]
      )
      
      # Update the reactive value
      dataset_visium(current_visium)
      original_visium(current_original)

      
      tmp_region_exists = all_region_ann_exists()
      tmp_region_exists[[input$dataset_name]] = T
      all_region_ann_exists(tmp_region_exists)
      
    }
    
    mydata_visium = all_visium()
    mydata_visium[[input$dataset_name]] = dataset_visium()
    
    all_visium(mydata_visium)
    all_original(mydata_visium)
    
    myids = individual_ids()
    myids = c(myids, input$dataset_name)
    individual_ids(myids)
    
    updateSelectInput(session,
      'qc_plot_type_data', 'Select the Individual',
      choices = individual_ids())
    
    updateSelectInput(session,
                      'qc_filter_data_sp', 'Specify the Individual(s)',
                      choices = individual_ids())
    
    updateSelectInput(session,
                      'normal_data_sp', 'Specify the Individual(s)',
                      choices = individual_ids())
    
    updateSelectInput(session,
                      'hvg_data_sp', 'Specify the Individual(s)',
                      choices = individual_ids())
    
    updateSelectInput(session,
                      'dim_red_data_sp', 'Specify the Individual(s)',
                      choices = individual_ids())
    
    updateSelectInput(session,
                      'clustering_data_sp', 'Specify the Individual(s)',
                      choices = individual_ids())
    
    output$uploaded_datasets_table <- renderDT({
      req(length(individual_ids()) > 0)
      myresult = cbind('Datasets Uploaded' = individual_ids()) |> as.data.frame()
    })
    
    updateSelectInput(session,
                      'spatial_vis_data_ID', 'Select the Individual',
                      choices = individual_ids())
    
    output$add_individual = renderUI({
      
      
      actionButton('add_individual_Button', 'Add new individual')
      
    })
    
    output$remove_individual = renderUI({
      
      
      actionButton('remove_individual_Button', 'Remove existing individual(s)')
      
    })
    

    
    updateSelectInput(session,
                      'dea_data_sp', 'Specify the Individual(s)',
                      choices = individual_ids())
    
    updateSelectInput(session,
                      'spatial_auto_data', 'Select the Individual',
                      choices = individual_ids())
    
    updateSelectInput(session,
                      'spatial_context_ind_sp', 'Select individual',
                      choices = individual_ids())
    
    updateSelectInput(session,
                      'ev_synth_data', 'Select the Individual',
                      choices = individual_ids())
    
    reset("h5_file")
    reset("spatial_folder")
    reset("dataset_neoantigen_status")
    reset("dataset_region_ann")
    
    removeModal()
    
  })
  

  
  observeEvent(input$load_visium,{
    
    req(input$h5_file, input$spatial_folder)
    
    output$box_summary = renderUI({
      box(
        
        title = 'Dataset(s) Summary', status = "primary", solidHeader = TRUE,
        width = 12,
        column(
          width = 3,
          DTOutput("uploaded_datasets_table") 
        ),
        column(
          width = 4,
          selectInput('data_summary_select', 'Select data to summarise', choices = individual_ids(), multiple = F)
          
        ),
        column(
          width = 1,
          div(style = "margin-top: 25px;",  # Adjust this value as needed
              actionButton('data_summary_select_Button', 'Show')
          )
          
        ),
        column(
          width = 4,
          verbatimTextOutput("data_visium_summary")
        )
        
      )
      
    })
    
  })
  
  
  observeEvent(input$data_summary_select_Button,{
    
    req(all_visium(), input$data_summary_select)
    mydata = all_visium()[[input$data_summary_select]]
    output$data_visium_summary <- renderPrint({
      
      cat("=== Summary ===\n")
      cat("Number of spots (barcodes):", ncol(mydata), "\n")
      cat("Number of genes/features:", nrow(mydata), "\n")
      
    })
  })
  


  
  
  
  # --------------------------------------------------------------
  # --------------------------------------------------------------
  # --------------------------------------------------------------
  
  
  # --------------------------------------------------------------
  #                          PREPROCESSING
  # --------------------------------------------------------------
  
  # History track
  preproc_history_track = reactiveVal(list())
  
  # Quality Control Plot
  qc_plot = reactiveVal(NULL)
  
  observeEvent(input$qc_vis_Button, {
    req(all_visium(),input$qc_plot_type_data %in% names(all_visium()))
    data = all_visium()[[input$qc_plot_type_data]]
    
    data[["percent.mt"]] = PercentageFeatureSet(data, pattern = "^MT-")
    
    tmp_qc_data <- data@meta.data %>%
      select(any_of(c('nFeature_Spatial', 'nCount_Spatial', 'percent.mt'))) %>%
      mutate(cell = rownames(.))
    
    if (!(any(colnames(tmp_qc_data) == 'percent.mt'))){
      
      tmp_qc_data['percent.mt'] = 0
    }
    
    if (input$qc_plot_type == 'Violin plots'){
      
      

      # Define the metrics you want to facet
      metrics <- c("nFeature_Spatial", "nCount_Spatial", 'percent.mt')

      # Generate one violin plot per metric
      plots <- lapply(metrics, function(metric) {
        plot_ly(
          data = tmp_qc_data,
          y = as.formula(paste0("~", metric)),
          type = 'violin',
          box = list(visible = TRUE),
          meanline = list(visible = TRUE),
          points = 'all',
          jitter = 0.3,
          marker = list(
            color = 'black',
            size = 3,
            opacity = 0.5),
          name = metric,
          hoverinfo = 'y'
        ) %>%
          layout(title = 'Violin Plots', yaxis = list(title = metric), showlegend = FALSE)
      })

      # Combine using subplot
      p = subplot(plots, nrows = 1, shareX = TRUE, titleX = TRUE, titleY = F)

      
    } else if (input$qc_plot_type == 'Scatterplot') {
      
      
      p = plot_ly(
        data = tmp_qc_data,
        x = ~nCount_Spatial,
        y = ~nFeature_Spatial,
        type = "scatter",
        mode = "markers",
        color = ~percent.mt,
        text = ~paste(
          "Counts:", nCount_Spatial,
          "<br>Genes:", nFeature_Spatial,
          "<br>Mitochondrial:", percent.mt
        ),
        hoverinfo = "text",
        marker = list(size = 6, opacity = 0.7)
      ) %>%
        layout(
          title = list(text = "Scatterplot", x = 0.5),
          xaxis = list(title = "nCount_Spatial"),
          yaxis = list(tile = 'nFeature_Spatial'),
          legend = list(
            orientation = "h",
            x = 0.5,
            y = -0.25,
            xanchor = "center")
        ) %>%
        
        config(
          displaylogo = FALSE,
          modeBarButtonsToRemove = list('hoverCompareCartesian')
        )
      
      
    } else {
      
      spatial_data <- GetTissueCoordinates(data)
      spatial_data = spatial_data[, c('x', 'y')]
      expression_data <- FetchData(data, vars = "nCount_Spatial")
      plot_data <- cbind(spatial_data, expression_data)
      
      p = plot_ly(
        data = plot_data,
        x = ~y, 
        y = ~-x,  
        type = "scatter",
        mode = "markers",
        color = ~nCount_Spatial,
        colors = viridis::viridis(100),  # Using viridis color scale
        marker = list(
          size = 5,  
          opacity = 0.8,
          line = list(width = 0)
        ),
        text = ~paste("nCount_Spatial:", round(nCount_Spatial, 2)),
        hoverinfo = "text"
      ) %>%
        layout(
          title = "Spatial Plot coloured by number of counts per Spots",
          xaxis = list(
            title = "",
            showgrid = FALSE,
            zeroline = FALSE,
            showticklabels = FALSE
          ),
          yaxis = list(
            title = "",
            scaleanchor = "x",  # Maintain aspect ratio
            showgrid = FALSE,
            zeroline = FALSE,
            showticklabels = FALSE
          ),
          legend = list(
            title = list(text = "<b>nCount_Spatial</b>"),
            orientation = "v",
            x = 1.05,
            xanchor = "left"
          ),
          plot_bgcolor = 'rgba(0,0,0,0)',  # Transparent background
          paper_bgcolor = 'rgba(0,0,0,0)'
        ) %>%
        colorbar(
          title = "nCount",
          len = 0.5,
          yanchor = "middle"
        ) %>%
        config(
          displayModeBar = TRUE,
          scrollZoom = TRUE
        )
      
      
    }
    
    qc_plot(p)
    

    output$qc_plot <- renderPlotly({
      qc_plot()
    })
    
    output$get_download_qc_plot = renderUI({
      
      downloadButton("download_qc_plot_button", "Download Plot")
      
    })
    
    
  })
  
  output$download_qc_plot_button <- downloadHandler(
    filename = function() {
      paste("qc_plot_", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      # Ensure the plot exists before downloading
      req(qc_plot())
      
      # Save the plot
      png(file, width = 8, height = 6, units = "in", res = 300)
      print(qc_plot())
      dev.off()
    }
  )
  

  
  
  # Quality Control Filters
  qc_filter_exists = reactiveVal(list())
  
  observeEvent(input$qc_Button, {
    req(all_visium(), input$qc_filter_data)
    
    mydata = all_visium()
    
    if (input$qc_filter_data == 'Any'){
      
      myids = names(mydata)
      
    } else {
      
      myids = input$qc_filter_data_sp
      
    }
    
    tmp_qc_filter_exists = qc_filter_exists()
    
    mydata_history_track = preproc_history_track()
    
    for (i in myids){
      
      tmp_qc_filter_exists[[i]] = T
      
      filtered_data = mydata[[i]]
      
      filtered_data[["percent.mt"]] <- PercentageFeatureSet(filtered_data, pattern = "^MT-")
      
      filtered_data <- subset(filtered_data, 
                              subset = percent.mt < input$max_mito)
      
      ## Number of genes
      
      if(input$n_genes_filter == "min") {
        filtered_data <- subset(filtered_data, 
                                subset = nFeature_Spatial > input$min_genes)
        
        if ('QC_Filter' %in% names(mydata_history_track[[i]])){
          
          mydata_history_track[[i]]['QC_Filter'] = paste0(mydata_history_track[[i]]['QC_Filter'], ' & Min genes: ', input$min_genes)
          
        } else {
          
          mydata_history_track[[i]] = c(mydata_history_track[[i]], 'QC_Filter' = paste0('Min genes: ', input$min_genes))
        }
        
      } else if (input$n_genes_filter == "max") {
        
        filtered_data <- subset(filtered_data, 
                                subset = nFeature_Spatial < input$max_genes)
        
        if ('QC_Filter' %in% names(mydata_history_track[[i]])){
          
          mydata_history_track[[i]]['QC_Filter'] = paste0(mydata_history_track[[i]]['QC_Filter'], ' & Max genes: ', input$max_genes)
          
        } else {
          
          mydata_history_track[[i]] = c(mydata_history_track[[i]], 'QC_Filter' = paste0('Max genes: ', input$max_genes))
        }
        
        
      } else if (input$n_genes_filter == "both") {
        
        filtered_data <- subset(filtered_data, 
                                subset = nFeature_Spatial > input$min_genes & 
                                  nFeature_Spatial < input$max_genes)
        
        if ('QC_Filter' %in% names(mydata_history_track[[i]])){
          
          mydata_history_track[[i]]['QC_Filter'] = paste0(mydata_history_track[[i]]['QC_Filter'], ' & Genes within: (', input$min_genes, ' - ', input$max_genes, ')')
          
        } else {
          
          mydata_history_track[[i]] = c(mydata_history_track[[i]], 'QC_Filter' = paste0('Genes within: (', input$min_genes, ' - ', input$max_genes, ')'))
        }
        
      }
      
      
      ## Number of Counts
      
      if(input$n_counts_filter == "min") {
        filtered_data <- subset(filtered_data, 
                                subset = nCount_Spatial > input$min_counts)
        
        if ('QC_Filter' %in% names(mydata_history_track[[i]])){
          
          mydata_history_track[[i]]['QC_Filter'] = paste0(mydata_history_track[[i]]['QC_Filter'], ' & Min counts: ', input$min_counts)
          
        } else {
          
          mydata_history_track[[i]] = c(mydata_history_track[[i]], 'QC_Filter' = paste0('Min counts: ', input$min_counts))
        }
        
      } 
      else if (input$n_counts_filter == "max") {
        
        filtered_data <- subset(filtered_data, 
                                subset = nCount_Spatial < input$max_counts)
        
        if ('QC_Filter' %in% names(mydata_history_track[[i]])){
          
          mydata_history_track[[i]]['QC_Filter'] = paste0(mydata_history_track[[i]]['QC_Filter'], ' & Max counts: ', input$max_counts)
          
        } else {
          
          mydata_history_track[[i]] = c(mydata_history_track[[i]], 'QC_Filter' = paste0('Max counts: ', input$max_counts))
        }
        
        
      } else if (input$n_counts_filter == "both") {
        
        filtered_data <- subset(filtered_data, 
                                subset = nCount_Spatial > input$min_counts & 
                                  nCount_Spatial < input$max_counts)
        
        if ('QC_Filter' %in% names(mydata_history_track[[i]])){
          
          mydata_history_track[[i]]['QC_Filter'] = paste0(mydata_history_track[[i]]['QC_Filter'], ' & counts within: (', input$min_counts, ' - ', input$max_counts, ')')
          
        } else {
          
          mydata_history_track[[i]] = c(mydata_history_track[[i]], 'QC_Filter' = paste0('counts within: (', input$min_counts, ' - ', input$max_counts, ')'))
        }
        
      }
      
      mydata[[i]] = filtered_data
      
      
    }
    
    all_visium(mydata)
    qc_filter_exists(tmp_qc_filter_exists)

    preproc_history_track(mydata_history_track)
    
    showNotification("Filters applied successfully!")
  })
  
  observeEvent(input$qc_reset_Button, {
    req(input$qc_filter_data)
    
    showModal(modalDialog(
      title = "Are you sure?",
      ifelse(input$qc_filter_data == 'Any',
             paste0("This will remove all filters and restore all the original datasets."),
             paste0("This will remove all filters from ", paste(input$qc_filter_data_sp, collapse = ', '), " and restore their original datasets")),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("qc_reset_Button_confirm", "Yes, Reset")
      )
    ))
    
    
  })
  
  
  observeEvent(input$qc_reset_Button_confirm, {
    
    req(all_visium(), all_original(), input$qc_filter_data)
    
    mydata = all_visium()
    mydata_original = all_original()
    
    myhistory = preproc_history_track()
    
    
    if (input$qc_filter_data == 'Any'){
      
      myids = names(mydata)
      
    } else {
      
      myids = input$qc_filter_data_sp
      
    }
    
    tmp_qc_filter_exists = qc_filter_exists()
    tmp_normal_exists = normal_exists()
    tmp_hvg_exists = hvg_exists()
    
    
    for (i in myids){
      
      mydata[[i]] = mydata_original[[i]]
      tmp_qc_filter_exists[[i]] = F
      tmp_normal_exists[[i]] = F
      tmp_hvg_exists[[i]] = F
      
      
      myhistory[[i]] = c()
      
    }
    
    all_visium(mydata)
    qc_filter_exists(tmp_qc_filter_exists)
    normal_exists(tmp_normal_exists)
    hvg_exists(tmp_hvg_exists)
    preproc_history_track(myhistory)
    
    removeModal()
  })
  
  
  
  
  # Normalization
  output$normal_plot_data <- renderUI({
    selectInput("normal_plot_data", "Select the Individual to plot", choices = individual_ids(), multiple = F)
  })
  
  output$normal_plot_Button <- renderUI({
    actionButton("normal_plot_Button", "Display Plot")
  })
  
  
  normal_exists = reactiveVal()
  
  
  observeEvent(input$normal_Button, {
    req(all_visium())
    
    mydata = all_visium()
    
    if (input$normal_data == 'Any'){
      
      myids = names(mydata)
      
    } else {
      
      myids = input$normal_data_sp
      
    }
    
    
    tmp_normal_exists = normal_exists()
    
    mydata_history_track = preproc_history_track()
    
    for (i in myids){
      
      # Update the normalization history
      if ('Normalization' %in% names(mydata_history_track[[i]])){
        
        mydata_history_track[[i]]['Normalization'] = paste0(mydata_history_track[[i]]['Normalization'], ' & ', input$normal)
        
      } else {
        
        mydata_history_track[[i]] = c(mydata_history_track[[i]], 'Normalization' = input$normal)
      }
      
     
      tryCatch({
        showNotification(paste0("Normalizing ", i, '...'), duration = NULL, id = paste0('normalizing_', i))
        
        options(future.globals.maxSize = 1000 * 1024^2)
        
        # Apply Normalization
        if (input$normal == 'Log-normalization'){
          
          tmp_hvg_exists = hvg_exists()
          
          mydata[[i]] = NormalizeData(mydata[[i]], normalization.method = "LogNormalize", scale.factor = 10000)
          
          mydata[[i]] <- FindVariableFeatures(
            mydata[[i]],
            selection.method = input$hvg_selection_method,
            nfeatures = input$hvg_number
          )
          
          # Update the Feaure Selection
          if ('Feature_Selection' %in% names(mydata_history_track[[i]])){
            
            mydata_history_track[[i]]['Feature_Selection'] = paste0(mydata_history_track[[i]]['Feature_Selection'],
                                                                    ' & HVG:',
                                                                    input$hvg_selection_method)
            
          } else {
            
            mydata_history_track[[i]] = c(mydata_history_track[[i]], 'Feature_Selection' = input$hvg_selection_method)
          }
          
        } else {
          
          mydata[[i]] <- SCTransform(mydata[[i]], assay = "Spatial")
          
        }
        
        
        removeNotification(id = paste0('normalizing_', i))
        showNotification(paste("Normalization of", i, "complete!"))
        
      }, error = function(e) {
        showNotification(paste("Normalization of ", i, "failed:", e$message), type = "error")
      }) 
      
    }
    
    all_visium(mydata)
    
    normal_exists(tmp_normal_exists)
    
    if (input$normal == 'Log-normalization'){hvg_exists(tmp_hvg_exists)}
    
    preproc_history_track(mydata_history_track)
    
    
    # output$hvg_plot_data <- renderUI({
    #   selectInput("hvg_plot_data", "Select the Individual to plot", choices = individual_ids(), multiple = F)
    # })
    
    
    output$box_normal_plot = renderUI({
      
      box(
        
        title = "Mean-Variance Relationship", 
        status = "primary", 
        solidHeader = TRUE,
        width = 9,
        
        fluidRow(
          column(
            width = 4,
            uiOutput('normal_plot_data')
          ),
          
          column(
            width = 4,  # Plot controls
            div(style = "margin-top: 25px;",  # Adjust this value as needed
                uiOutput('normal_plot_Button')
            )
          )
        ),
        
        fluidRow(
          
          column(
            
            width = 12,
            align = 'center',
            plotOutput('normal_plot')
            
          )
          
        ),
        
        
        fluidRow(
          
          column(
            
            width = 12,
            align = 'right',
            uiOutput('get_download_normal_plot')
            
            
          )
          
          
        )
        
      )
      
      
    })

  })
  
  
  normal_plot = reactiveVal(NULL)
  
  observeEvent(input$normal_plot_Button,{
    
    req(all_visium(),input$normal_plot_data %in% names(all_visium()))
    mydata <- all_visium()[[input$normal_plot_data]]
    
    p = VariableFeaturePlot(mydata) +
      
      ggtitle("Mean-Variance Relationship")
    
    # p = VlnPlot(mydata,
    #             cols = 'steelblue',
    #             features = c("nFeature_Spatial", "nCount_Spatial"), 
    #             pt.size = 0.1, ncol = 2) &
    #   
    #   theme(
    #     plot.title = element_text(size = 10),
    #     axis.title.x = element_blank(), 
    #     axis.text.x = element_blank(),  
    #     axis.ticks.x = element_blank()
    #   )
    
    normal_plot(p)
    
    output$normal_plot <- renderPlot({
      
      normal_plot()

    })
    
    
    output$get_download_normal_plot = renderUI({
      
      downloadButton("download_normal_plot_button", "Download Plot")
      
    })
    
    
    
  })
  
  
  output$download_normal_plot_button <- downloadHandler(
    
    filename = function() {
      paste("normal_plot_", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      # Ensure the plot exists before downloading
      req(normal_plot())
      
      # Save the plot
      png(file, width = 8, height = 6, units = "in", res = 300)
      print(normal_plot())
      dev.off()
    }
    
    )
  
  


  
  
  
  # Feature Selection (Primary)
  hvg_exists = reactiveVal()
  
  
  output$hvg <- renderUI({
    numericInput("hvg_number", "Number of highly variable genes", value = 2000, min = 500, max = 5000)
  })
  
  
  output$hvg_selection_method <- renderUI({
    selectInput("hvg_selection_method", "Select the HVG method", choices = c('vst', 'mean.var.plot', 'dispersion'), multiple = F)
  })
  
  
  # observeEvent(input$hvg_Button, {
  #   req(all_visium(), input$hvg_data, input$hvg_number, input$hvg_selection_method)  # Require data and user input
  #   
  #   mydata = all_visium()
  #   
  #   if (input$hvg_data == 'Any'){
  #     
  #     myids = names(mydata)
  #     
  #   } else {
  #     
  #     myids = input$hvg_data_sp
  #     
  #   }
  #   
  #   tmp_hvg_exists = hvg_exists()
  #   
  #   mydata_history_track = preproc_history_track()
  #   
  #   for (i in myids){
  #     
  #     tryCatch({
  #       # Show processing message
  #       showNotification(paste0("Selecting HVGs for" ,i, "..."), duration = NULL, id = paste0("hvg_msg_", i))
  #       
  #       # Update Seurat object with HVGs
  #       mydata[[i]] <- FindVariableFeatures(
  #         mydata[[i]],
  #         selection.method = input$hvg_selection_method,
  #         nfeatures = input$hvg_number  # User-defined number of features
  #       )
  #       
  #       if ('Feature_Selection' %in% names(mydata_history_track[[i]])){
  #         
  #         mydata_history_track[[i]]['Feature_Selection'] = paste0(mydata_history_track[[i]]['Feature_Selection'],
  #                                                                 ' & HVG:',
  #                                                                 input$hvg_selection_method)
  #         
  #       } else {
  #         
  #         mydata_history_track[[i]] = c(mydata_history_track[[i]], 'Feature_Selection' = input$hvg_selection_method)
  #       }
  #       
  #       
  #       # Status message
  #       output$hvg_status <- renderText({
  #         paste("Success! Selected", input$hvg_number, "HVGs.")
  #       })
  #       
  #       removeNotification(id = paste0("hvg_msg_", i))
  #       showNotification(paste0("HVG selection of ", i, " complete!"), type = "message")
  #       
  #     }, error = function(e) {
  #       showNotification(paste("Error:", e$message), type = "error")
  #       output$hvg_status <- renderText("HVG selection failed.")
  #     })
  #     
  #   }
  #   
  #   # Update reactive object
  #   all_visium(mydata)
  #   
  #   hvg_exists(tmp_hvg_exists)
  #   
  #   preproc_history_track(mydata_history_track)
  #   
  #   output$hvg_plot_data <- renderUI({
  #     selectInput("hvg_plot_data", "Select the Individual to plot", choices = individual_ids(), multiple = F)
  #   })
  #   
  #   
  #   output$hvg_plot_Button <- renderUI({
  #     actionButton("hvg_plot_Button", "Display Plot")
  #   })
  # 
  # })
  

  
  
  
  hvg_plot = reactiveVal(NULL)
  
  observeEvent(input$hvg_plot_Button,{
    
    req(all_visium(),input$hvg_plot_data %in% names(all_visium()))
    mydata <- all_visium()[[input$hvg_plot_data]]
    
    
    p = VariableFeaturePlot(mydata) +
      
      ggtitle("Mean-Variance Relationship After Normalization")
    
    hvg_plot(p)
    
    output$hvg_plot <- renderPlot({
      
      hvg_plot()
      
    })
    
    output$get_download_hvg_plot = renderUI({
      
      downloadButton("download_hvg_plot_button", "Download Plot")
      
    })
    
    
  })
  
  output$download_hvg_plot_button <- downloadHandler(
    filename = function() {
      paste("hvg_plot_", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      # Ensure the plot exists before downloading
      req(hvg_plot())
      
      # Save the plot
      png(file, width = 8, height = 6, units = "in", res = 300)
      print(hvg_plot())
      dev.off()
    }
  )
  

  # Dimensionality Reduction
  output$dim_red <- renderUI({
    selectInput("dim_red_method", "Select the method", choices = c('PCA', 'UMAP'), multiple = F)
  })
  
  dim_red_plot = reactiveVal(NULL)
  
  observeEvent(input$dim_red_Button, {
    req(all_visium(), input$dim_red_data_sp, input$dim_red_method, input$number_dims)
    
    mydata = all_visium()
    
    tryCatch({
      showNotification(paste0("Running ", input$dim_red_method, " analysis for ", input$dim_red_data_sp, "..."),
                       duration = NULL, id = paste0("running_msg_", input$dim_red_method, '_', input$dim_red_data_sp))
      
      obj <- mydata[[input$dim_red_data_sp]]
      
      # SCALE THE DATA (CRITICAL STEP)
      obj <- ScaleData(
        obj,
        features = VariableFeatures(obj),  # Use variable features if defined
        verbose = FALSE
      )
      
      if (input$dim_red_method == "PCA") {
        # Run PCA
        obj <- RunPCA(
          obj,
          npcs = input$number_dims,
          verbose = FALSE
        )
        reduction <- "pca"
        title <- "PCA Visualization"
      } else {
        # Run UMAP (requires PCA first)
        if (!"pca" %in% names(obj@reductions)) {
          obj <- RunPCA(obj, npcs = input$number_dims, verbose = FALSE)
        }
        obj <- RunUMAP(
          obj,
          dims = 1:input$number_dims,
          verbose = FALSE
        )
        reduction <- "umap"
        title <- "UMAP Visualization"
      }
      
      
      # Update the reactive object
      mydata[[input$dim_red_data_sp]] = obj
      all_visium(mydata)
      
      
      p = DimPlot(obj,
                  reduction = reduction,
                  pt.size = 1.5
                  ) +
        ggtitle(title) +
        theme_minimal() +
        theme(legend.position = 'none')
      
      dim_red_plot(p)
      
      
      # Plot results
      output$dim_red_plot <- renderPlot({
        
        dim_red_plot()
      })
      
      output$get_download_dim_red_plot = renderUI({
        
        downloadButton("download_dim_red_plot_button", "Download Plot")
        
      })
      
      
      output$elbow_Button = renderUI({
        
        actionButton("elbow_Button", "Display Elbow Plot")
        
      })
      
      
      
      output$status <- renderText({
        paste("Success! Applied", input$dim_red_method, "with", input$number_dims, "dimensions.")
      })
      
      removeNotification(id = paste0("running_msg_", input$dim_red_method, '_', input$dim_red_data_sp))
      showNotification(paste0(input$dim_red_method, " analysis for ", input$dim_red_data_sp, " complete!"), type = "message")
      
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
      output$status <- renderText(paste("Failed:", e$message))
    })
  })
  
  output$download_dim_red_plot_button <- downloadHandler(
    filename = function() {
      paste("dim_red_plot_", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      # Ensure the plot exists before downloading
      req(output$dim_red_plot)
      
      # Save the plot
      png(file, width = 8, height = 6, units = "in", res = 300)
      print(output$dim_red_plot)
      dev.off()
    }
  )
  
  elbow_plot = reactiveVal()
  
  observeEvent(input$elbow_Button,{
    
    req(all_visium(), input$dim_red_data_sp)
    
    mydata = all_visium()[[input$dim_red_data_sp]]
    
    p = ElbowPlot(mydata)
    
    elbow_plot(p)
    
    output$elbow_plot = renderPlot({
      
      elbow_plot()
      
    })
    
    
  })
  
  # Clustering
  num_dimensions = reactiveVal(NULL)
  neighbors_algorithm = reactiveVal(NULL)
  
  
  observeEvent(input$neighbors_Button, {
    
    req(all_visium(), input$clustering_data_sp)
    
    mydata = all_visium()
    
    mydata[[input$clustering_data_sp]] = FindNeighbors(mydata[[input$clustering_data_sp]], dims = 1:input$dimens_num)
    
    num_dimensions(input$dimens_num)
    
    all_visium(mydata)
    
    output$clustering_algorithm = renderUI({
      
      selectInput('clustering_algorithm_select',
                  'Select the clustering algorithm',
                  choices = c('Louvain', 'Leiden'), multiple = F)
      
    })
    
    neighbors_algorithm(input$clustering_algorithm_select)
    
    output$clustering_Button = renderUI({
      
      actionButton('clustering_Button', 'Run Clustering')
      
    })
    
  })
  
  
  observeEvent(input$clustering_Button, {
    
    req(all_visium(), input$clustering_algorithm_select, input$clustering_data_sp)
    
    mydata = all_visium()
    
    algorithm = ifelse(input$clustering_algorithm_select == 'Louvain', 1, 2)
    
    mydata[[input$clustering_data_sp]] = FindClusters(
      mydata[[input$clustering_data_sp]],
      resolution = 0.5,    # Start with 0.5 (adjust as needed)
      algorithm = algorithm,       # 1 = Louvain (default), 2 = Leiden
      random.seed = 42      # For reproducibility
    )
    
    all_visium(mydata)
    
    
    output$clustering_plot_Button = renderUI({
      
      actionButton('clustering_plot_Button', 'Plot Clusters')
      
    })
    
    
  })
  
  
  clustering_plot = reactiveVal(NULL)
  
  
  observeEvent(input$clustering_plot_Button,{
    
    req(all_visium(), input$clustering_data_sp)
    
    mydata = all_visium()[[input$clustering_data_sp]]
    
    p = DimPlot(mydata, group.by = "seurat_clusters", label = TRUE) +
      ggtitle("UMAP: Clusters")
    
    output$clustering_plot = renderPlot({
      
      p
      
    })
    
    clustering_plot(p)
    
  })
  
  
  output$preprocess_history_table = renderDataTable({
    
    myhistory = preproc_history_track()
    
    myhistory = lapply(myhistory, function(x){return(as.data.frame(t(x)))})
    
    myhistory = data.table::rbindlist(myhistory, idcol = 'Individual_ID', fill = T)
    
    myhistory
    
  })
  
  

  

  # --------------------------------------------------------------
  # --------------------------------------------------------------
  # --------------------------------------------------------------
  
  # --------------------------------------------------------------
  #               Deconvolution
  # --------------------------------------------------------------
  single_cell_data = reactiveVal()
  
  observeEvent(input$scrna_upload_Button, {
    req(input$scrna_upload)
    
    showModal(modalDialog(
      title = "Please wait",
      "Single-cell data Uploading...",
      footer = NULL,
      easyClose = FALSE
    ))
    
    scrna <- anndata::read_h5ad(input$scrna_upload$datapath)
    
    single_cell_data(scrna)
    
    print(scrna)
    
    removeModal()
    
  })
  
  output$deconv_Button_ui = renderUI({
    
    if (input$deconv_method == "bayespace"){
      
      actionButton("deconv_Button", "Run Deconvolution")
      
    } else if (input$deconv_method == "spotlight" & !is.null(single_cell_data())){
      
      actionButton("deconv_Button", "Run Deconvolution")
      
    } else {
      
      return(NULL)
      
    }
    
  })
  
  
  
  output$results_validate_ui = renderUI({
    
    box(
      title = 'Model Validation Plots', status = "primary", solidHeader = TRUE,
      width = 12,
      
      column(
        width = 6,
        plotOutput('topic_plot')
      ),
      column(
        width = 6,
        plotOutput('topic_facet_plot')
      )
      
      
    )
    
  })

  
  
  deconv_results <- eventReactive(input$deconv_Button, {
    req(all_visium(), input$deconv_data)
    
    mydata = all_visium()
    
    myresults = list()
    
    if (input$deconv_data == 'Any'){
      
      myids = names(mydata)
      
    } else {
      
      myids = input$deconv_data_sp
      
    }
    
    if (input$deconv_method == "spotlight") {
      
      showModal(modalDialog(
        title = "Please wait",
        "Deconvolution with SPOTlight Running...",
        footer = NULL,
        easyClose = FALSE
      ))
      
      req(single_cell_data())
      
      scrna = single_cell_data()
      
      scrna_counts <- Matrix::Matrix(t(as.matrix(scrna$X)), sparse = TRUE)
      
      print(head(scrna_counts))
      
      cell_types <- scrna$obs[['cell_type']] |> as.character()
      
      print(unique(cell_types))
      
      for (i in myids){
        
        mydata_i = mydata[[i]]
        
        visium_counts = mydata_i[['Spatial']]$counts
        
        
        common_genes = intersect(rownames(visium_counts), rownames(scrna_counts)) # No common genes: They use different names!!
        
        print(common_genes)
        
        # Map gene names to each other
        print('Start mapping gene names to each other')
        ensembl <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
        
        gene_symbols <- rownames(visium_counts)
        
        mapping <- biomaRt::getBM(
          attributes = c("ensembl_gene_id", "hgnc_symbol"),  # Columns to retrieve
          filters = "hgnc_symbol",                           # Filter by gene symbol
          values = gene_symbols,                             # Your gene list
          mart = ensembl                                    # Ensembl connection
        )
        
        mapping = mapping %>% filter(ensembl_gene_id %in% rownames(scrna_counts))
        
        print('Mapping done!')
        
        print(head(mapping))
        
        scrna_counts = scrna_counts[mapping$ensembl_gene_id, ]
        
        visium_counts = visium_counts[mapping$hgnc_symbol, ]
        
        
        rownames(scrna_counts) = rownames(visium_counts)
        
        
        # ==============================================================================
        # 4. Convert to SC and ST Experiments
        # ==============================================================================
        sce <- SingleCellExperiment(
          assays = list(counts = scrna_counts,
                        logcounts = as.matrix(scrna_counts))
        )
        
        
        colData(sce)$cell_type <- cell_types
        
        colLabels(sce) = colData(sce)$cell_type
        
        coords <- GetTissueCoordinates(mydata_i)
        coords = coords[, c('x', 'y')]
        
        spe <- SpatialExperiment::SpatialExperiment(
          assays = list(counts = visium_counts),
          colData = coords
        )
        
        
        mgs <- scoreMarkers(sce, groups = sce$cell_type)
        
        mgs_fil = list()
        
        for (j in names(mgs)) {
          
          mgs_fil[[j]] = mgs[[j]] %>%
            
            as.data.frame() %>%
            
            mutate(Gene = rownames(.)) %>%
            
            filter(mean.AUC > 0.4) %>%
            
            arrange(desc(mean.AUC))
          
          
        }
        
        mgs_df = data.table::rbindlist(mgs_fil, idcol = 'Cluster')
        
        print('MGS derived!')
        
        # split cell indices by identity
        idx <- split(seq(ncol(sce)), sce$cell_type)
        
        # downsample to at most 20 per identity & subset
        # We are using 5 here to speed up the process but set to 75-100 for your real
        # life analysis
        
        n_cells <- 5
        
        cs_keep <- lapply(idx, function(i) {
          n <- length(i)
          if (n < n_cells)
            n_cells <- n
          sample(i, n_cells)
        })
        
        sce <- sce[, unlist(cs_keep)]
        
        print('Training the model...')
        
        res <- SPOTlight(
          x = sce,
          y = spe,
          groups = as.character(sce$cell_type),
          mgs = mgs_df,
          weight_id = "mean.AUC",
          group_id = "Cluster",
          gene_id = "Gene")
      
        print(head(res$mat))
        
        mydata_i <- AddMetaData(mydata_i, metadata = res$mat)
        
        myresults[[i]] = list(
          method = "SPOTlight",
          visium = mydata_i,
          results = res$mat,
          mod = res$NMF
        )
        
        
      }
      
      removeModal()
      
      print('Finished!')
      
      return(myresults)
      
      
    } else if (input$deconv_method == "bayespace") {
      
      showModal(modalDialog(
        title = "Please wait",
        "Deconvolution with BayesSpace Running...",
        footer = NULL,
        easyClose = FALSE
      ))
      
      
      for (i in myids){
        
        mydata_i = mydata[[i]]
        
        mydata_i = ScaleData(mydata_i)
        
        # Run BayesSpace clustering
        sce <- as.SingleCellExperiment(mydata_i)
        
        coords <- GetTissueCoordinates(mydata_i)
        coords = coords[, c('x', 'y')]
        
        # Ensure spatial coordinates are in metadata
        colData(sce)$row <- coords[, 'y']  
        colData(sce)$col <- coords[, 'x'] 
        
        visium <- BayesSpace::spatialPreprocess(sce,
                                                platform = "Visium",
                                                skip.PCA = F,
                                                n.PCs = 15)
        
        colData(visium)$array_row <- colData(visium)$row
        colData(visium)$array_col <- colData(visium)$col
        colData(visium)$spot.idx <- seq_len(ncol(visium))
        
        
        visium_cluster <- spatialCluster(
          visium,              
          q = input$n_clusters_deconv,             
          platform = "Visium",
          use.dimred = "PCA", 
          nrep = 1000,       
          burn.in = 100,
          model = input$bayespace_model
        )
        
        myresults[[i]] = list(
          method = "BayesSpace",
          visium = visium_cluster,
          results = colData(visium_cluster)[, "spatial.cluster", drop = FALSE]
        )
        
        # tryCatch({
        #   
        #   
        #   
        #   
        #   return(myresults[[i]])
        #   
        #   
        # }, error = function(e) {
        #   showNotification(paste("Error analyzing", i, ":", e$message), 
        #                    type = "error")
        # })

      }
      
      removeModal()
      
      print('Finished!')
      
      return(myresults)
      

    }
    
  })
  
  
  deconv_data_sp_show_options = reactiveVal()
  
  observeEvent(input$deconv_Button, {
    
    req(deconv_results())
    
    myids = names(deconv_results())
    deconv_data_sp_show_options(myids)
    
    output$deconv_data_sp_show_ui = renderUI({
      
      selectInput('deconv_data_sp_show', 'Specify the Individual(s)',
                  choices = deconv_data_sp_show_options(), multiple = F)
      
    })
    
    
    output$deconv_show_Button_ui = renderUI({
      
      actionButton('deconv_show_Button_ui', 'Show results')
      
    })
    
  })
  
  observeEvent(input$deconv_show_Button_ui,{
    
    req(deconv_results(), all_visium())
    
    res <- deconv_results()[[input$deconv_data_sp_show]]
    
    mydata <- all_visium()[[input$deconv_data_sp_show]]
    
    # Plot spatial results
    output$deconv_plot <- renderPlot({
      
      if (res$method == "SPOTlight") {
        
        # Plot the first cell type (can be modified to allow selection)
        SpatialFeaturePlot(visium, features = colnames(res$results)[1])
        
      } else if (res$method == 'BayesSpace') {
        
        # Add to Seurat metadata (with identical barcode order)
        mydata$bayespace_clusters <- res$visium$spatial.cluster
        mydata$bayespace_clusters = as.factor(mydata$bayespace_clusters)
        
        print(SpatialDimPlot(mydata,
                             group.by = 'bayespace_clusters',
                             image.alpha = 0.5, 
                             pt.size.factor = 1.5) +
                theme(legend.position = 'bottom') +
                
                labs(fill = 'BayesSpace Cluster'))
        
      }
    })
    
    # # Show results table
    # output$deconv_table <- renderDataTable({
    #   
    #   res$results
    # })
    # 
    # # Download results
    # output$download_deconv <- downloadHandler(
    #   filename = function() {
    #     paste0(res$method, "_results.csv")
    #   },
    #   content = function(file) {
    #     write.csv(res$results, file)
    #   }
    # )
    
    
  })

  
  
  
  
  # Show download button
  output$download_deconv_ui <- renderUI({
    req(deconv_results())
    div(style = "margin-top: 20px;",
        downloadButton("download_deconv", "Download Results")
    )
  })
  
  
  # --------------------------------------------------------------
  #               Spatial Visualization
  # --------------------------------------------------------------
  selected_spatial_vis_group <- reactiveVal("None")
  
  spatial_vis_group_options = reactiveVal('None')
  
  # Reactive to generate visualization options based on current ID and filters
  potential_spatial_vis_data <- reactive({
    req(input$spatial_vis_data_ID)
    
    if (isTRUE(qc_filter_exists()[[input$spatial_vis_data_ID]])) {
      c('Original - Unfiltered', 'Filtered (from Quality Control)')
    } else {
      'Original - Unfiltered'
    }
  })
  
  
  observeEvent(input$spatial_vis_data_ID,{
    
    myoptions = 'None'
    
    if (isTRUE(all_neo_status_exists()[[input$spatial_vis_data_ID]])){
      myoptions = c(myoptions, 'Neoantigen Status') |> unique()


    }

    if (isTRUE(all_region_ann_exists()[[input$spatial_vis_data_ID]])){

      myoptions = c(myoptions, 'Region') |> unique()

    }
    
    spatial_vis_group_options(myoptions)
    
    output$spatial_vis_group <- renderUI({
      selectInput("spatial_vis_group",
                  "Select the grouping variable", 
                  choices = spatial_vis_group_options(),
                  multiple = FALSE)
    })
    
    
    output$spatial_vis_data <- renderUI({
      selectInput("spatial_vis_data",
                  "Select the dataset for visualization", 
                  choices = potential_spatial_vis_data(),
                  multiple = FALSE)
    })
    
    
  })

  
  
  
  
  
  # Store the last generated plot
  last_spatial_plot <- reactiveVal(NULL)
  
  # Update selected_group and generate plot when button is clicked
  observeEvent(input$spatial_vis_Button, {
    req(all_visium(), input$spatial_vis_group, input$spatial_vis_data_ID)
    
    selected_spatial_vis_group(gsub(' ', '_', input$spatial_vis_group))
    
    if (input$spatial_vis_data == 'Original - Unfiltered'){
      
      mydata = all_original()[[input$spatial_vis_data_ID]]
      
    } else {
      
      mydata = all_visium()[[input$spatial_vis_data_ID]]
      
    }
    
    
    # Generate the plot
    current_plot = if (selected_spatial_vis_group() == 'None') {
      
      Idents(mydata) <- "all_spots"
      
      SpatialDimPlot(mydata,
                     image.alpha = 0.5, 
                     pt.size.factor = 1.5) +
        theme(legend.position = 'none')
      
    } else {
      
      SpatialDimPlot(mydata,
                     group.by = selected_spatial_vis_group(),
                     image.alpha = 0.5, 
                     pt.size.factor = 1.5) +
        theme(legend.position = 'bottom',
              legend.text = element_text(size = 14),
              legend.title = element_text(size = 16, face = "bold"))
    }
    
    # Store the plot for downloading
    last_spatial_plot(current_plot)
    
    # Render the plot
    output$spatial_vis_plot = renderPlot({ current_plot })
    
    # Show download button
    output$download_spatial_plot_ui <- renderUI({
      div(style = "margin-top: 20px;",
          downloadButton("download_spatial_plot", "Download Plot")
      )
    })
  })
  
  # Download handler
  output$download_spatial_plot <- downloadHandler(
    filename = function() {
      paste("spatial_plot_", selected_spatial_vis_group(), "_", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      req(last_spatial_plot())
      ggsave(file, 
             plot = last_spatial_plot(),
             device = "png",
             width = 8, 
             height = 6,
             dpi = 300)
    }
  )
  
  # --------------------------------------------------------------
  #               DIFFERENTIAL EXPRESSION ANALYSIS
  # --------------------------------------------------------------
  
  output$dea_appear_message <- renderPrint({
    
    cat('Please import data through "Import data" tab to activate this functionality')
    
  })
  
  output$dea_preproc_message <- renderPrint({
    
    cat('*Note: Differential Expression Analysis can only be conducted in normalized data. Please verify that you only select datasets under normalization. You can check the normalization status in the "Preprocessing History" within the "Preprocessing" tab.  \n')

    
  })
  
  
  output$dea_appear = renderUI({
    
    if (length(all_visium()) == 0){
      
      box(
        title = "Differential Expression Analysis between Neoantigen-positive and negative spots",
        status = "primary", solidHeader = TRUE,
        width = 12,
        column(
          width = 12,
          align = 'center',
          verbatimTextOutput('dea_appear_message') %>%
            tagAppendAttributes(style = "color: red; font-weight: bold;")
        )
        
      )


      
    } else {
      
      box(
        title = "Differential Expression Analysis between Neoantigen-positive and negative spots",
        status = "primary", solidHeader = TRUE,
        width = 12,
        column(
          width = 4,
          selectInput('dea_test', 'Select the DEA method',
                      choices = c('Wilcoxon Rank Sum' = 'wilcox',
                                  'Wilcoxon Rank Sum (limma)' = 'wilcox_limma',
                                  'Likelihood-ratio' = 'bimod',
                                  # 'ROC analysis' = 'roc',
                                  'Student t-test' = 't',
                                  'Negative binomial' = 'negbinom',
                                  'Poisson' = 'poisson',
                                  'Logistic Regression' = 'LR',
                                  'MAST', 'DESeq2'),
                      multiple = F),
          
        ),
        column(
          width = 4,
          radioButtons("dea_data", "Select the Individual(s)*", choices = c('Any', 'Specific'), inline = T)
        ),
        column(
          width = 4,
          conditionalPanel(
            condition = "input.dea_data == 'Specific'",
            selectInput('dea_data_sp', 'Specify the Individual(s)', choices = NULL, multiple = T)
          )
        ),
        column(width = 12, actionButton("deaButton", "Run")),
        column(width = 12,
               textOutput('dea_preproc_message'))
      )
      
    }
    
    
    
  })
  
  
  
  potential_dea_individuals = reactiveVal(NULL)
  
  
  # Differential Expression Analysis
  results_dea = reactiveVal()
  

  observeEvent(input$deaButton,{
    req(all_visium(), all_neo_status_exists(), input$dea_data)
    
    mydata = all_visium()
    
    if (input$dea_data == 'Any'){
      
      myids = names(mydata)
      
    } else {
      
      myids = input$dea_data_sp
      
    }
    
    tmp_all_neo_status_exists = all_neo_status_exists()
    
    dea_results_list = list()
    
    showModal(modalDialog(
      title = "Please wait",
      "DEA Running...",
      footer = NULL,
      easyClose = FALSE
    ))
    
    withProgress(message = 'Running Differential Expression Analysis...', value = 0.5, {
      
      for (i in myids){
        
        if (isTRUE(tmp_all_neo_status_exists[[i]])){
          
          # Set active identity for the analysis
          
          tmp_obj = mydata[[i]]
          
          Idents(tmp_obj) <- "Neoantigen_Status"
          
          tryCatch({
            markers <- FindMarkers(
              tmp_obj,
              ident.1 = "Positive",
              ident.2 = "Negative",
              min.pct = 0.1,
              test.use = input$dea_test
              # logfc.threshold = 0.25
            )
            
            # Prepare results table
            markers$Gene <- rownames(markers)
            rownames(markers) <- NULL
            
            # Calculate SE based on test type
            test <- input$dea_test
            
            if (test %in% c("wilcox", "wilcox_limma", "t", "bimod", "LR", "poisson", "negbinom", 'MAST')) {
              # Approximate Z and SE
              markers$Z <- qnorm(1 - markers$p_val_adj / 2)
              markers$SE <- abs(markers$avg_log2FC / markers$Z)
              
            } else if (test == "DESeq2" && "lfcSE" %in% colnames(markers)) {
              markers$SE <- markers$lfcSE
              
            } else {
              # Fallback: NA with a warning
              markers$SE <- NA
              warning("SE could not be estimated for the selected test.")
            }
            
            dea_results_list[[i]] <- markers %>%
              select(Gene, `pct.1`, `pct.2`, avg_log2FC, SE, p_val, p_val_adj) %>%
              arrange(p_val_adj)
            
          }, error = function(e) {
            showNotification(paste("Error analyzing", i, ":", e$message), 
                             type = "error")
            dea_results_list[[i]] <- NULL
          })
          
        } else {
          
          showNotification(paste0(i, ' has no neoantigen status information. Analysis omitted'), type = "warning")
          
        }
        
      }
      
      incProgress(1, detail = "Done!")
      
      dea_results_list = Filter(Negate(is.null), dea_results_list)
      
      potential_dea_individuals(names(dea_results_list))
      
      if (length(dea_results_list) == 0){
        
        showNotification("No valid DEA results generated", type = "warning")
        return(NULL)
      }
      
      results_dea(dea_results_list)
      
    })
    
    removeModal()

    output$dea_results_box = renderUI({
      
      box(
        title = "DEA Results",
        status = "primary", solidHeader = TRUE,
        width = 12,
        uiOutput('dea_results_data_show'),
        uiOutput('dea_results_data_show_Button'),
        
        dataTableOutput('dea_results'),
        uiOutput('get_download_results_dea'),
        uiOutput("volcano_controls"),
        plotlyOutput('volcano_plot'),
        uiOutput('get_download_results_dea_plot')
      )
      
    })
    
    output$dea_results_data_show = renderUI({

      selectInput('dea_results_data_show', 'Select the Individual',
                  choice = potential_dea_individuals(),
                  multiple = F)

    })

    output$dea_results_data_show_Button = renderUI({

      actionButton('dea_results_data_show_Button', 'Show results')

    })

  })
  
  # Display results table
  observeEvent(input$dea_results_data_show_Button, {
    
    
    req(results_dea(), input$dea_results_data_show)
    
    myresults = results_dea()
    

    myresults = myresults[[input$dea_results_data_show]]

    output$dea_results <- renderDataTable({
      
      myresults[, 2:7] = sapply(myresults[, 2:7], function(x){x = round(x, 4)})
      
      myresults$GeneCard = paste0('https://www.genecards.org/cgi-bin/carddisp.pl?gene=', myresults$Gene)
      
      myresults$GeneCard = sprintf(
        '<a href="%s" target="_blank">Browse</a>', 
        myresults$GeneCard
      )
      
      myresults = myresults[, c(1,8,2:7)]

      DT::datatable(
        myresults,
        escape = FALSE
      )
      


    })
    
    # Download DEA results
    output$get_download_results_dea <- renderUI({
      req(results_dea())
      downloadButton("download_results_dea", "Download Table")
    })
    
  })
  
  
  output$download_results_dea <- downloadHandler(
    filename = function() {
      paste("DEA_Results_", Sys.Date(), ".xlsx", sep = "")
    },
    content = function(file) {
      openxlsx::write.xlsx(results_dea(), file)
    }
  )
  
  
  # Track if DEA has been run
  dea_has_run <- reactiveVal(FALSE)
  
  # Update flag when DEA runs
  observeEvent(input$deaButton, {
    dea_has_run(TRUE)
  })
  
  # Render volcano controls conditionally
  output$volcano_controls <- renderUI({
    req(dea_has_run())  # Only show after DEA has run
    
    tagList(
      fluidRow(
        column(
          width = 6,
          selectInput("p_value_type", "P-value Type:",
                      choices = c("Adjusted p-value" = "p_val_adj",
                                  "Raw p-value" = "p_val"),
                      selected = "p_val_adj")
        ),
        column(
          width = 3,
          numericInput("p_cutoff", "Significance Cutoff:",
                       value = 0.05, min = 0.0001, max = 1, step = 0.01)
        ),
        column(
          width = 3,
          numericInput("fc_cutoff", "Log2 FC Threshold:",
                       value = 0.5, min = 0, step = 0.1)
        )
      ),
      
      fluidRow(
        column(
          width = 4,
          selectInput("sig_color", "Significant Color:",
                      choices = c("Red" = "red", "Green" = "green3",
                                  "Purple" = "purple", "Orange" = "orange"),
                      selected = "red")
        ),
        column(
          width = 4,
          selectInput("ns_color", "Non-significant Color:",
                      choices = c("Gray" = "gray", "Blue" = "blue",
                                  "Black" = "black", "Light blue" = "lightblue"),
                      selected = "gray")
        ),

      ),
      
      actionButton('volcano_create', 'Create the volcano plot')
      )
      
    
    
  })
  
  # Volcano plot (Dynamic)
  results_dea_plot <- eventReactive(input$volcano_create, {
    req(results_dea(), input$dea_results_data_show, input$p_value_type, input$p_cutoff,
        input$fc_cutoff, input$sig_color, input$ns_color)
    
    dea_df <- results_dea()[[input$dea_results_data_show]]
    
    # Create significance column based on user parameters
    dea_df$significant <- dea_df[[input$p_value_type]] < input$p_cutoff & 
      abs(dea_df$avg_log2FC) > input$fc_cutoff
    
    plot_ly(
      data = dea_df,
      x = ~avg_log2FC,
      y = ~-log10(.data[[input$p_value_type]]),
      type = "scatter",
      mode = "markers",
      color = ~factor(significant, 
                      levels = c(FALSE, TRUE),
                      labels = c("Not significant", "Significant")),
      colors = c("Not significant" = input$ns_color,
                 "Significant" = input$sig_color),
      text = ~paste(
        "Gene:", Gene,
        "<br>Log2FC:", round(avg_log2FC, 2),
        "<br>P-value:", signif(.data[[input$p_value_type]], 3)
      ),
      hoverinfo = "text",
      marker = list(size = 6, opacity = 0.7)
    ) %>%
      layout(
        title = list(text = "Volcano Plot", x = 0.5),
        xaxis = list(title = "Log2 Fold Change"),
        yaxis = list(
          title = if (input$p_value_type == "p_val_adj") {
            "-Log10 Adjusted P-value"
          } else {
            "-Log10 P-value"
          }
        ),
        legend = list(
          orientation = "h",
          x = 0.5,
          y = -0.25,
          xanchor = "center"),
        shapes = list(
          # Horizontal p-value cutoff line
          list(
            type = "line",
            x0 = min(dea_df$avg_log2FC),
            x1 = max(dea_df$avg_log2FC),
            y0 = -log10(input$p_cutoff),
            y1 = -log10(input$p_cutoff),
            line = list(dash = "dash", color = "black")
          ),
          # Vertical logFC cutoff lines
          list(
            type = "line",
            x0 = input$fc_cutoff,
            x1 = input$fc_cutoff,
            y0 = 0,
            y1 = max(-log10(dea_df[[input$p_value_type]]), na.rm = TRUE),
            line = list(dash = "dash", color = "black")
          ),
          list(
            type = "line",
            x0 = -input$fc_cutoff,
            x1 = -input$fc_cutoff,
            y0 = 0,
            y1 = max(-log10(dea_df[[input$p_value_type]]), na.rm = TRUE),
            line = list(dash = "dash", color = "black")
          )
        )
      ) %>%
      
      config(
        displaylogo = FALSE,
        modeBarButtonsToRemove = list("toImage", 'hoverCompareCartesian')
      )

    
  })
  
  results_dea_plot_download <- eventReactive(input$volcano_create, {
    req(results_dea(), input$dea_results_data_show, input$p_value_type, input$p_cutoff,
        input$fc_cutoff, input$sig_color, input$ns_color)
    
    dea_df <- results_dea()[[input$dea_results_data_show]]
    
    # Create significance column based on user parameters
    dea_df$significant <- dea_df[[input$p_value_type]] < input$p_cutoff & 
      abs(dea_df$avg_log2FC) > input$fc_cutoff

    
    ggplot(dea_df, aes(x = avg_log2FC, y = -log10(.data[[input$p_value_type]]))) +
      geom_point(aes(color = significant), alpha = 0.7, size = 2.5) +
      scale_color_manual(
        values = c("TRUE" = input$sig_color, "FALSE" = input$ns_color),
        labels = c("TRUE" = "Significant", "FALSE" = "Not significant"),
        name = paste0("p < ", input$p_cutoff, " & |FC| > ", input$fc_cutoff)
      ) +
      geom_hline(yintercept = -log10(input$p_cutoff),
                 linetype = "dashed", color = "black") +
      geom_vline(xintercept = c(-input$fc_cutoff, input$fc_cutoff),
                 linetype = "dashed", color = "black") +
      theme_minimal(base_size = 14) +
      labs(
        title = "Volcano Plot",
        x = "Log2 Fold Change",
        y = ifelse(input$p_value_type == "p_val_adj",
                   "-Log10 Adjusted P-value",
                   "-Log10 P-value")
      ) +
      theme(
        legend.position = "bottom",
        plot.title = element_text(face = "bold", hjust = 0.5),
        panel.grid.major = element_line(color = "gray90"),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)
      )
  })
  
  output$volcano_plot <- renderPlotly({
    results_dea_plot()
  })
  
  output$get_download_results_dea_plot <- renderUI({
    req(results_dea_plot_download())
    downloadButton("download_results_dea_plot", "Download Plot")
  })
  
  output$download_results_dea_plot <- downloadHandler(
    filename = function() {
      paste("Volcano_Plot_", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      ggsave(file, plot = results_dea_plot_download(),
             width = 8, height = 6, dpi = 300)
    }
  )

  
  # --------------------------------------------------------------
  #               SPATIAL AUTOCORRELATION
  # --------------------------------------------------------------
  output$spatial_auto_preproc_message <- renderPrint({
    
    cat('*Note: Spatial Autocorrelation Analysis can only be conducted in normalized data. Please verify that you only select datasets under normalization. You can check the normalization status in the "Preprocessing History" within the "Preprocessing" tab.  \n')
    
    
  })
  
  output$spatial_auto_genes_message <- renderPrint({
    
    cat('**Note: A typical dataset may contain thousands of genes. If all of them are selected, the analysis may take too long to run.  \n')
    
    
  })
  
  genes_available = reactive({

    req(all_visium(), input$tabs == "spatial_autocorrelation", input$spatial_auto_data)

    data = all_visium()[[input$spatial_auto_data]]

    mygens = rownames(data)
    mygens = mygens[!is.na(mygens)]
    mygens = sort(mygens)

  })
  
  observe({
    
    
    req(genes_available())  # Wait until genes are available
    updateSelectizeInput(
      session,
      "morans_genes_define",
      choices = genes_available(),
      server = T
    )
  })
  
  genes_selected_for_moran = reactiveVal(NULL)
  

  
  
  moran_results = eventReactive(input$morans_I, {
    
    showModal(modalDialog(
      title = "Please wait",
      "Running...",
      footer = NULL,
      easyClose = FALSE
    ))
    
    
    req(all_visium(), input$spatial_auto_data)
    
    mydata = all_visium()[[input$spatial_auto_data]]

    # Get coordinates and data
    coords <- GetTissueCoordinates(mydata)
    coords = coords[, c('x', 'y')]
    
    data <- GetAssayData(mydata, slot = "data")
    
    
    
    spatial_points <- SpatialPoints(coords)
    knn_graph <- knn2nb(knearneigh(spatial_points, k = 6))
    weights <- nb2listw(knn_graph)
    
    results = list()
    
    mygenes = rownames(data)
    
    if (input$morans_genes_specific == 'Specific'){
      
      mygenes = input$morans_genes_define
      
    }
    
    mygenes = mygenes[!is.na(mygenes)]
    
    genes_selected_for_moran(mygenes)
    
    for (i in mygenes){
      
      if (any(data[i,] != 0)){
        
        mytest = moran.test(data[i,], listw = weights)
        
        results[[i]] = cbind(
          Moran_I = mytest$estimate["Moran I statistic"],
          p_value = mytest$p.value,
          Expected = mytest$estimate["Expectation"],
          Variance = mytest$estimate["Variance"]
        ) %>% as.data.frame()
        
      } else {
        
        results[[i]] = cbind(
          Moran_I = NA,
          p_value = NA,
          Expected = NA,
          Variance = NA
        ) %>% as.data.frame()
        
        
      }
      
    }
    
    
    results = data.table::rbindlist(results, idcol = 'Gene')
    
    results$p_adj = p.adjust(results$p_value, method = 'fdr')
    
    removeModal()
    
    results
    
  })
  
  output$moran_table = renderDT({
    req(moran_results())
    datatable(
      moran_results(),
      options = list(pageLength = 10, scrollX = TRUE),
      rownames = TRUE
    ) %>% formatRound(columns = 2:6, digits = 4)
  })
  
  
  observeEvent(input$morans_I,{
    
    output$spatial_plot_gene <- renderUI({
      selectInput("spatial_plot_gene", 'Select Gene', choices = genes_selected_for_moran(), multiple = F)
    })
    
    output$spatial_plot_Button <- renderUI({
      actionButton("spatial_plot_Button", "Display Plot")
    })
    
    output$moran_plot_Button <- renderUI({
      actionButton("moran_plot_Button", "Display Plot")
    })
    
    
  })
  
  
  observeEvent(input$spatial_plot_Button, {
    
    output$spatial_plot <- renderPlot({
      req(all_visium(), input$spatial_auto_data, input$spatial_plot_gene)
      
      SpatialFeaturePlot(
        all_visium()[[input$spatial_auto_data]],
        features = input$spatial_plot_gene,
        pt.size.factor = 1.5,
        alpha = c(0.8, 1)
      ) + ggtitle(paste("Spatial Expression:", input$spatial_plot_gene)) +
        theme(legend.position = 'right')
    })
    
    
  })
  
  moran_plot = reactiveVal(NULL)
  moran_plot_download = reactiveVal(NULL)
  
  observeEvent(input$moran_plot_Button, {
    
    output$moran_plot <- renderPlotly({
      req(moran_results())
      
      data = moran_results() %>%
        
        mutate(Significant = ifelse(p_value < 0.05, 'Yes', 'No'))
      
      
      p = plot_ly(
        data = data,
        x = ~Moran_I,
        y = ~-log10(p_value),
        type = 'scatter',
        mode = 'markers',
        color = ~Significant,
        colors = c("grey", "red"),
        marker = list(size = 6, opacity = 0.7),
        text = ~paste(
          "Gene:", Gene,
          "<br>Moran's I:", round(Moran_I, 3),
          "<br>p-value:", signif(p_value, 3)
        ),
        hoverinfo = "text"
      ) %>%
        layout(
          title = list(text = "Moran's I Volcano Plot", x = 0.5),
          xaxis = list(title = "Moran's I Statistic"),
          yaxis = list(title = "-log10(p-value)"),
          shapes = list(
            list(
              type = "line",
              x0 = min(data$Moran_I, na.rm = TRUE),
              x1 = max(data$Moran_I, na.rm = TRUE),
              y0 = -log10(0.05),
              y1 = -log10(0.05),
              line = list(dash = "dash", color = "black")
            )
          ),
          legend = list(
            orientation = "h",
            x = 0.5,
            y = -0.25,
            xanchor = "center"
          )
        ) %>%
        config(displaylogo = FALSE,
               modeBarButtonsToRemove = list("toImage", 'hoverCompareCartesian'))
      

      

      p_download = ggplot(data %>% filter(!is.na(Significant)),
                 aes(x = Moran_I, y = -log10(p_value), color = Significant)) +

        geom_point(size = 3, alpha = 0.7) +

        scale_color_manual(values = c("grey", "red")) +

        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +

        labs(title = "Moran's I Volcano Plot",
             x = "Moran's I Statistic",
             y = "-log10(p-value)") +

        theme_minimal() +
        
        theme(legend.position = 'bottom')


      
      moran_plot(p)
      moran_plot_download(p_download)
      
      p
        
    })
    
    output$get_download_moran_plot = renderUI({
      
      downloadButton("download_moran_plot_button", "Download Plot")
      
    })
    
    
  })
  
  output$download_moran_plot_button <- downloadHandler(
    filename = function() {
      paste("moran_plot_", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      # Ensure the plot exists before downloading
      req(moran_plot())
      
      # Save the plot
      png(file, width = 8, height = 6, units = "in", res = 300)
      print(moran_plot_download())
      dev.off()
    }
  )
  
  
  

  # --------------------------------------------------------------
  # --------------------------------------------------------------
  # --------------------------------------------------------------
  
  # --------------------------------------------------------------
  #  ------------        SPATIAL CONTEXT               -----------
  # --------------------------------------------------------------
  
  observeEvent(input$spatial_context_Button,{
    
    req(all_visium(), input$spatial_context_ind_sp)
    
    mydata = all_visium()[[input$spatial_context_ind_sp]]

    # Prepare data for MISTy
    expr_data <- GetAssayData(mydata, slot = "data")
    coords <- GetTissueCoordinates(mydata)
    coords = coords[, c('x', 'y')]
    
    misty = mistyR::misty(
      data = as.data.frame(t(expr_data)),
      spatial_coordinates = coords,
      views = c("intrinsic", "juxtaview", "paraview"),
      group_variable = "Neoantigen_Status"
    )
    
    
  })

  # --------------------------------------------------------------
  #  ------------        EVIDENCE SYNTHESIS            -----------
  # --------------------------------------------------------------
  
  common_genes_dea = reactiveVal()
  
  observe({
    
    req(results_dea())
    
    output$ev_synth_summary <- renderPrint({
      
      cat("=== Summary ===\n")
      if (length(results_dea()) == 0){
        
        cat('There are no datasets with DEA results')
        
      } else{
        
        
      }
      cat('Number of datasets with DEA results: ', length(results_dea()),  '\n')
      
      myresults = data.table::rbindlist(results_dea(), idcol = 'Individual_ID')
      total_genes_DEA = myresults$Gene |> unique() |> length()
      cat('Number of genes included in DEA results:', total_genes_DEA, '\n')
      
      n_common_genes = myresults$Gene |> table()
      n_common_genes = n_common_genes[n_common_genes == length(results_dea())] |> length()
      common_genes_dea(n_common_genes)
      cat("Number of common genes in all DEA results:", common_genes_dea(), "\n")
      
    })
    
    
  })
  

  genes_available_dea = reactive({
    
    req(results_dea())
    
    data = data.table::rbindlist(results_dea(), idcol = 'Individual_ID')
    
    mygens = data$Gene |> unique()
    mygens = mygens[!is.na(mygens)]
    mygens = sort(mygens)
    
  })
  
  observe({
    
    
    req(genes_available_dea())  # Wait until genes are available
    updateSelectizeInput(
      session,
      "ev_synth_gene",
      choices = genes_available_dea(),
      server = T
    )
  })
  
  genes_selected_for_ev_synth = reactiveVal(NULL)
    
  inds_selected_for_ev_synth = reactiveVal(NULL)
  
  observeEvent(input$deaButton, {
    req(results_dea())
    
    myresults = data.table::rbindlist(results_dea(), idcol = 'Individual_ID')
    
    myids = myresults$Individual_ID |> unique()
    
    inds_selected_for_ev_synth(myids)
    
    output$ev_synth_data_ui = renderUI({
      
      radioButtons("ev_synth_data", "Select the Individual(s) to include",
                   choices = c('Any', 'Specific'), inline = T)
      
    })
    
    updateSelectInput(session,
                      'ev_synth_data_sp', 'Select the Individual(s)',
                      choices = inds_selected_for_ev_synth())
    
    output$ev_synth_forest_Button_ui = renderUI({
      
      actionButton('ev_synth_forest_Button', 'Display Forest plot')
      
    })
    
  })
  
  ev_synth_forestplot = reactiveVal()
  
  observeEvent(input$ev_synth_forest_Button,{
    
    req(results_dea(), input$ev_synth_gene, input$ev_synth_data)
    
    myresults = data.table::rbindlist(results_dea(), idcol = 'Individual_ID')
    
    myids = myresults$Individual_ID |> unique()
    
    if (input$ev_synth_data == 'Any') {
      myids <- myresults$Individual_ID |> unique()
    } else {
      myids <- input$ev_synth_data_sp
    }

    
    df = myresults %>% 
      
      filter(Individual_ID %in% myids, Gene %in% input$ev_synth_gene) %>%
      
      mutate(SE = ifelse(is.na(SE), 0, SE))
    
    # 95% CI
    df$lower <- df$avg_log2FC - 1.96 * df$SE
    df$upper <- df$avg_log2FC + 1.96 * df$SE
    

    p = ggplot(df, aes(y = avg_log2FC, x = Individual_ID, colour = Gene)) +
      
      geom_point(size = 2,
                 position = position_dodge(width = 0.2)) +
      
      geom_errorbar(aes(ymin = lower, ymax = upper),
                    width = 0,
                    linewidth = 0.7,
                    position = position_dodge(width = 0.2)) +
      
      geom_hline(yintercept = 0, linetype = "dashed", color = "blue", linewidth = 0.5) +
      
      labs(
        y = "log2 Fold Change",
        x = "Sample"
      ) +
      theme_classic() +
      
      theme(legend.position = 'bottom')
    
    ev_synth_forestplot(p)
    
    output$ev_synth_forest_plot = renderPlot({
      
      ev_synth_forestplot()
      
    })

    output$get_download_ev_synth_forestplot = renderUI({
      
      downloadButton("download_ev_synth_forest_plot_button", "Download Plot")
      
    })
    
    
  })
  
  output$download_ev_synth_forest_plot_button <- downloadHandler(
    filename = function() {
      paste("forest_plot_", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      # Ensure the plot exists before downloading
      req(ev_synth_forestplot())
      
      # Save the plot
      png(file, width = 8, height = 6, units = "in", res = 300)
      print(ev_synth_forestplot())
      dev.off()
    }
  )
    
    
  
}



# Run the application
shinyApp(ui = ui, server = server)


pkgs.list <- c("shiny", "shinyjs","shinydashboard", "reticulate", "plotly", "dplyr", "ggplot2", "BiocParallel", "Biostrings", "ggplot2", "tidyr", "DT")#, "circlize", "BiocManager", 
# , "viridis", "ggrepel", "forcats", "shinycssloaders",
#, "Seurat", "reshape2", "qs")

# setwd("~/6_TransposableElements/TEammo_app-dev/")
for (i in 1:length(pkgs.list)){
  if(!require(pkgs.list[i], character.only = T)){
    install.packages(pkgs.list[i])
    require(pkgs.list[i], character.only = TRUE)
  }else{require(pkgs.list[1],character.only = TRUE)}
}
teClassification <- read.table("orozco_classification-2024_mchelper.tsv.csv", sep = ";", header = TRUE)
teClassification <- teClassification %>%
  mutate(AppName = paste(Class, Order, Superfamily, sep ="/") %>% sub("/$", "", .) %>% sub("/$", "", .))
teClassification <- rbind(
  teClassification
  , c("", "Unclassified", "", "Unclassified")
)
use_condaenv("MCHelper")

source("TEammo_functions.R")
options(shiny.maxRequestSize = 100*1024^2) # 100 MB, change according to genome sizes to upload
# #####-
# # Debug----
# # #####-
setwd("~/6_TransposableElements/TEammo_app_v2.1/")
# input <- list()
# input$teBuscoLib <- "BUSCO_libs/nematoda_odb10_ALL.hmm"
# input$leftTrim <- 2
# input$rightTrim <- 70
# fileNames <- list()
# fileNames$teGenomes <- list.files(path = "0_raw", pattern = glob2rx("*a"), recursive = TRUE, full.names = TRUE)
# teLibs <- list()
# clickData <- list()
# teTables <- list()
# # input$te80InputLib <- "./MCHelper_automatic/C.elegans/N2-test/1_MI_MCH/complete_models.fa"



FL_thresh <- 0.9 #set by default by TE-Aid, MCHelper states it at 0.94, optimum would be to take the input of MCHelper
os <- 400 #min ORF length
mchelper_ui <- dashboardPage(
  title = "MCHelper", skin = "yellow"
  , dashboardHeader(title = span(tagList(icon("yammer", lib = "font-awesome"),"TEammo")))
  , dashboardSidebar(
    sidebarMenu(
      menuItem("Genome assembly", tabName = "teGenomeAssembly", badgeColor = "red", badgeLabel = "dev")
      , menuItem("Repeat modelling", tabName = "teRepeatModelling", badgeColor = "red", badgeLabel = "dev")
      , menuItem("MCHelper", tabName = "teMchelper", selected = TRUE)
      , menuItem("HELIANO", tabName = "teHeliano", selected = TRUE)
    )
  )# end dashboardSidebar
  , dashboardBody(
    tabItems(
      tabItem(tabName = "teGenomeAssembly"
              , h1("Space to assemble genomes, either uploading self-generated data (long reads) or by SRA ids")
              , tabsetPanel(
                tabPanel(
                  title = "SRA ID"
                ) #end of tabPanel for SRA ID
                , tabPanel(
                  title = "Upload data"
                )
              ) #end of tabsetPanel for teGenomeAssembly
      )# end of tabItem for teGenomeAssembly
      , tabItem(tabName = "teRepeatModelling"
                , h1("Space to run either RM2 or any other repeat modeller")
      ) # end of tabItem for teRepeatModelling
      , tabItem(
        tabName = "teMchelper"
        , tabsetPanel(
          id = "mchTabs"
          # #############-
          # ## MCHelper Input UI----
          # #############-
          , tabPanel(
            title = "Input data"
            , DTOutput("inputData")
            , fluidRow(
              column(6,
                     fileInput("teGenomeFile", "Genome fasta:", accept = c(".fna", ".fa", ".fasta", ".zip", ".tar.gz")),
                     fileInput("teRm2File", "Raw library:", accept = c(".fna", ".fa", ".fasta", ".zip", ".tar.gz"))
              ),
              column(6,
                     textInput("teSpecies", "Species:", placeholder = "Enter species name"),
                     textInput("teStrain", "Strain:", placeholder = "Enter strain name")
              )
            ),
            actionButton("teSubmit", "Submit Data"),
            verbatimTextOutput("teUploadSummary", placeholder = TRUE)
            , verbatimTextOutput("debug")
          )
          # #############-
          # ## MCHelper curated UI----
          # #############-
          , tabPanel(
            title = "Curated seqs", value = "1234"
            # , h1("You have to run automatic module to start with annotation")
            , sidebarLayout(
              sidebarPanel = sidebarPanel(
                width = 4
                , fluidRow(
                  column(6, numericInput("teAutoModuleCores", label = "Number of cores", min = 1, max = parallel::detectCores(), value = round(parallel::detectCores()/1.25), width = "80%"))
                  , column(6, selectizeInput("teBuscoLib", label = "Busco library", choices = list.files("./BUSCO_libs", pattern = "ALL.hmm", recursive = TRUE)))#if more libs are needed they could be added to the BUSCO_libs folder, either manually or download throuhg app (this has to be coded)
                )
                , actionButton("teAutoModuleRun", label = "Run or load\nautomatic module")
              )#end of sidebarPanel for automatic tabPanel
              , mainPanel = mainPanel(
                fluidRow(
                  plotlyOutput("teNonCurated")
                )
              )#end of mainPanel for automatic tabPanle
            )#end of sidebarLayout for automatic tabpanel
            , h3("Manual reannotation")
            , sidebarLayout(
              sidebarPanel = sidebarPanel(
                width = 4
                , h5("TE library under inspection")
                , verbatimTextOutput("teManualInputLib", placeholder = TRUE)
                , selectizeInput("teManualSeq", label = "Select sequence to inspect", choices = NULL)
                , fluidRow(
                  column(3, numericInput("leftTrim", "Left trim", value = 0))
                  , column(3, numericInput("rightTrim", "Right trim", value = 0))
                  , column(6, selectizeInput("teDecision", label = "TE annotation", choices = teClassification$AppName))#Server side loading as we want to show current annotation # numericInput("teDecision", "Decision (1-42):", value = NA, min = 1, max = 42))
                )
                , fluidRow(
                  column(6, actionButton("teTrim", label = "Trim seq"))
                  , column(6, actionButton("teAnnotate", label = "Re-annotate"))
                  
                )
                , verbatimTextOutput("teTrimError")
                # , verbatimTextOutput("debug", placeholder = TRUE)
                , actionButton("teDiscard", label = "Discard seq")
                , tableOutput("teManualClass")
                , h4("TE coding")
                , tableOutput("teManualCoding")
                , h4("TE structure")
                , tableOutput("teManualStructure")
                , tableOutput("teManualStats")
              ) #end sidebarPanel for MCHelper_test
              , mainPanel = mainPanel(
                width = 8
                # , tableOutput("teManualCuratedClass")
                , fluidRow(
                  column(6, plotOutput("teConsDivergence"))
                  , column(6, plotlyOutput("teConsCoverage"))
                )
                , fluidRow(
                  column(6, plotOutput("teConsSelf"))
                  , column(6, plotOutput("teConsStructure"))
                )
              )# end of mainPanel for MCHelper
            )#end sidebarLayout for MCHelper
          ) #end tabPanel for Curated UI
          # #############-
          # ## MCHelper 80-80-80 UI----
          # #############-
          , tabPanel(
            title = "80-80-80 models", value = "56"
            , fluidRow(
              sidebarLayout(
                sidebarPanel = sidebarPanel(
                  width = 4
                  , h5("TE library to 80-80-80 blast")
                  , verbatimTextOutput("te80InputLib", placeholder = TRUE)
                  # , selectizeInput(inputId = "te80InputLib", label = "Complete sequences Input"
                  #                  , choices = list.files(pattern = "^complete_models.fa", recursive = TRUE, full.names = TRUE)
                  # )
                  , actionButton("te80ModuleRun", label = "Run 80-80-80 model")
                )# end of sidebarPanel for 80-80-80 models
                , mainPanel = mainPanel(
                  plotlyOutput("te80InputPlot")
                  # column(6, plotOutput("te80InputPlot"))
                  # , column(6, verbatimTextOutput("te80ModelLog", placeholder = TRUE))
                )#end of mainPanel for 80-80-80 models input
              )#end of sidebarLayout for 80-80-80 input
            )#end of fluidRow for contain run 808080 models input
            , fluidRow(
              sidebarLayout(
                sidebarPanel = sidebarPanel(
                  width = 4
                  , h3("Manual reannotation")
                  , h5("TE library under inspection")
                  , verbatimTextOutput("teNon80Lib", placeholder = TRUE)
                  , selectizeInput("te80Seq", label = "Select sequence to inspect", choices = NULL)
                  , textOutput("te80Message")
                  , fluidRow(
                    column(6, selectizeInput("te80Decision", label = "TE annotation", choices = teClassification$AppName))#Server side loading as we want to show current annotation # numericInput("teDecision", "Decision (1-42):", value = NA, min = 1, max = 42))
                    , column(6, actionButton("te80Annotate", label = "Re-annotate and Trim TE"))
                  )
                  , fluidRow(
                    column(3, numericInput("leftTrim80", "Left trim", value = 0))
                    , column(3, numericInput("rightTrim80", "Right trim", value = 0))
                    , column(6, actionButton("te80Stand", label = "Stand by seq"))
                  )
                  , verbatimTextOutput("te80TrimError")
                  , column(3, actionButton("te80Discard", label = "Discard seq")) #deprecated since trimmed library joins curated library
                  , tableOutput("te80Stats")
                ) #end sidebarPanel for MCHelper_test
                , mainPanel = mainPanel(
                  width = 8
                  , h1("TE Aid plots for curation")
                  , tableOutput("te80CuratedClass")
                  , fluidRow(
                    column(6, plotOutput("te80ConsDivergence"))
                    , column(6, plotlyOutput("te80ConsCoverage"))
                  )
                  , fluidRow(
                    column(6, plotOutput("te80ConsSelf"))
                    , column(6, plotOutput("te80ConsStructure"))
                  )
                )# end of mainPanel for MCHelper
              )#end sidebarLayout for MCHelper
            )#end of fluidRow for non808080 Manual inspection
          )# end of tabPanel for 80-80-80
          # #############-
          # ## MCHelper 70-70-70 UI----
          # #############-
          , tabPanel(
            title = "70-70-70 models", value = "7"
            , fluidRow(
              sidebarLayout(
                sidebarPanel = sidebarPanel(
                  width = 4
                  , h5("TE library under inspection")
                  , verbatimTextOutput("te70InputLib", placeholder = TRUE)
                  # , selectizeInput(inputId = "te70InputLib", label = "Complete sequences Input"
                  #                  , choices = list.files(pattern = "standby_complete.fa", recursive = TRUE, full.names = TRUE)
                  # )
                  , actionButton("te70ModuleRun", label = "Run 70-70-70 model")
                )# end of sidebarPanel for 70-70-70 models
                , mainPanel = mainPanel(
                  plotlyOutput("te70InputPlot")
                  # column(6, plotOutput("te70InputPlot"))
                  # , column(6, verbatimTextOutput("te70ModelLog", placeholder = TRUE))
                )#end of mainPanel for 70-70-70 models input
              )#end of sidebarLayout for 70-70-70 input
            )#end of fluidRow for contain run 707070 models
          )# end of tabPanel for 70-70-70
          # #############-
          # ## MCHelper Incomplete UI----
          # #############-
          , tabPanel(
            title = "Incomplete models", value = "8"
            , fluidRow(
              sidebarLayout(
                sidebarPanel = sidebarPanel(
                  width = 4
                  , selectizeInput(inputId = "teIncompleteInputLib", label = "Inomplete sequences Input"
                                   , choices = list.files(pattern = "^incomplete_models.fa", recursive = TRUE, full.names = TRUE)
                  )
                  , actionButton("teIncompleteModuleRun", label = "Run Incomplete module")
                )# end of sidebarPanel for Incomplete-Incomplete-Incomplete models
                , mainPanel = mainPanel(
                  # column(6, plotOutput("teIncompleteInputPlot"))
                  # , column(6, verbatimTextOutput("teIncompleteModelLog", placeholder = TRUE))
                  plotlyOutput("teIncompleteInputPlot")
                )#end of mainPanel for Incomplete-Incomplete-Incomplete models input
              )#end of sidebarLayout for Incomplete-Incomplete-Incomplete input
            )#end of fluidRow for contain run IncompleteIncompleteIncomplete models
            , fluidRow(
              sidebarLayout(
                sidebarPanel = sidebarPanel(
                  width = 4
                  , h3("Manual reannotation")
                  , h5("TE library under inspection")
                  , verbatimTextOutput("teIncompleteNon80InputLib", placeholder = TRUE)
                  , selectizeInput("teIncompleteSeq", label = "Select sequence to inspect", choices = NULL)
                  , textOutput("teIncompleteMessage")
                  , fluidRow(
                    column(6, selectizeInput("teIncompleteDecision", label = "TE annotation", choices = teClassification$AppName))#Server side loading as we want to show current annotation # numericInput("teDecision", "Decision (1-42):", value = NA, min = 1, max = 42))
                    # , column(3, numericInput("leftTrimIncomplete", "Left trim", value = 0))
                    # , column(3, numericInput("rightTrimIncomplete", "Right trim", value = 0))
                    , column(6, checkboxGroupInput(inputId = "teIncompleteType", label = "Elements to review", choices = c("MAVERICK" = "MAVERICK", "TIR/MITE" = "TIR|MITE")))
                  )
                  , fluidRow(
                    column(6, actionButton("teIncompleteAnnotate", label = "Re-annotate and Trim TE"))
                    # , column(3, actionButton("teIncompleteStand", label = "Stand by seq"))
                    , column(3, actionButton("teIncompleteDiscard", label = "Discard seq"))
                  )
                  , verbatimTextOutput("teIncompleteTrimError")
                  , tableOutput("teIncompleteStats")
                ) #end sidebarPanel for MCHelper_test
                , mainPanel = mainPanel(
                  width = 8
                  , h1("TE Aid plots for curation")
                  , tableOutput("teIncompleteCuratedClass")
                  , fluidRow(
                    column(6, plotOutput("teIncompleteConsDivergence"))
                    , column(6, plotlyOutput("teIncompleteConsCoverage"))
                  )
                  , fluidRow(
                    column(6, plotOutput("teIncompleteConsSelf"))
                    , column(6, plotOutput("teIncompleteConsStructure"))
                  )
                )# end of mainPanel for MCHelper
              )#end sidebarLayout for MCHelper
            )#end of fluidRow for nonIncompleteIncompleteIncomplete Manual inspection
            
          )# end of tabPanel for Incomplete
          # #############-
          # ## MCHelper Final UI----
          # #############-
          , tabPanel(
            title = "Final library", value = 9
            # , h1("You have to run automatic module to start with annotation")
            , sidebarLayout(
              sidebarPanel = sidebarPanel(
                width = 3
                , fluidRow(
                  column(6, numericInput("teFinalCores", label = "Number of cores", min = 1, max = parallel::detectCores(), value = round(parallel::detectCores()/1.25), width = "80%"))
                  , column(6, actionButton("teAddLibToDb", label = "Add final library to DB"))
                )
                , fluidRow(
                  column(6, actionButton("teFinalModuleRun", label = "Run or load Final module"))
                  , column(6, checkboxInput(inputId = "teFinalShowLog", label = "Show log file", value = FALSE))
                )
              )#end of sidebarPanel for automatic tabPanel
              , mainPanel = mainPanel(
                plotlyOutput("teFinal")
                # uiOutput("debug")
                # , h2("MCHelper automatic log file")
              )#end of mainPanel for automatic tabPanle
            )#end of sidebarLayout for automatic tabpanel
            , fluidRow(
              h1("TE Structure information for review")
              , sidebarLayout(
                sidebarPanel = sidebarPanel(
                  width = 3
                  , h3("Manual Review")
                  , h4("Final TE library")
                  , verbatimTextOutput("teFinalLib", placeholder = TRUE)
                  , fluidRow(
                    column(7, selectizeInput("teFinalSeq", label = "Select sequence to inspect", choices = NULL))
                    , column(5, selectizeInput("teFinalAnnotation", label = "TE annotation", choices = teClassification$AppName))#Server side loading as we want to show current annotation # numericInput("teDecision", "Decision (1-42):", value = NA, min = 1, max = 42))
                  )
                  , fluidRow(
                    column(4, actionButton("teFinalPrev", label = "Previous seq"))
                    , column(3, actionButton("teFinalNext", label = "Next seq"))
                    , column(3, actionButton("teFinalDiscard", label = "Discard"))
                    , column(2, actionButton("teFinalSave", label = "Save"))
                  )
                  , tableOutput("teFinalClass")
                  , h4("TE coding")
                  , tableOutput("teFinalCoding")
                  , h4("TE structure")
                  , tableOutput("teFinalStructure")
                  , tableOutput("teFinalStats")
                ), mainPanel = mainPanel(
                  width = 9
                  , tableOutput("teFinalCuratedClass")
                  , fluidRow(
                    column(6, plotOutput("teFinalConsDivergence"))
                    , column(6, plotlyOutput("teFinalConsCoverage"))
                  )
                  , fluidRow(
                    column(6, plotOutput("teFinalConsSelf"))
                    , column(6, plotOutput("teFinalConsStructure"))
                  )
                )
              )
            )#end TE Aid plots for the final library
          )# end of tabPanel for Final lib construction and inspection
        )# end MCHelper tabsetPanel
      )#end of tabItem for MCHelper
      , tabItem(tabName = "teHeliano"
                , h1("Select a genome or library to inspect for Helitron-like TE!")
                , tabsetPanel(
                  id = "helianoTabs"
                  # #############-
                  # ## Heliano Input UI----
                  # #############-
                  , tabPanel(
                    title = "Input data"
                    , DTOutput("helitronData")
                  )
                  # #############-
                  # ## Heliano Run UI----
                  # #############-
                  , tabPanel(
                    title = "HELIANO inspection", value = "helianoRun"
                    , fluidRow(
                      sidebarLayout(
                        sidebarPanel = sidebarPanel(
                          numericInput("helCores", label = "Number of cores", min = 1, max = parallel::detectCores(), value = round(parallel::detectCores()/1.25))
                          , actionButton("helianoRun", label = "Run or load HELIANO")
                        )
                        , mainPanel = mainPanel(
                          plotlyOutput("helianoRunLibs")
                        )
                      )
                    )
                    , sidebarLayout(
                      sidebarPanel = sidebarPanel(
                        width = 3
                        , selectizeInput("helManualSeq", label = "Select sequence to inspect", choices = NULL)
                      )
                      , mainPanel = mainPanel(
                        width = 9
                        # , tableOutput("helManualCuratedClass")
                        , fluidRow(
                          column(6, plotOutput("helConsDivergence"))
                          , column(6, plotlyOutput("helConsCoverage"))
                        )
                        , fluidRow(
                          column(6, plotOutput("helConsSelf"))
                          , column(6, plotOutput("helConsStructure"))
                        )
                      )
                    )
                  )
                )
      ) # end of tabItem for teRepeatModelling
      )#end of tabItems
    ) #end dashboardBody
)#end of dashboardPage


mchelper_server <- function(input, output, session) {
  #Initialize internal objects----
  updateSelectizeInput(session, inputId = "teBuscoLib", choices = list.files("./BUSCO_libs", pattern = "ALL.hmm", recursive = TRUE))
  fileNames <- reactiveValues(
    teGenomes = list.files(path = "0_raw", pattern = glob2rx("*a"), recursive = TRUE, full.names = TRUE)
    , teRm2Files = list.files("./RM2_output", pattern = "families.fa", recursive = TRUE, full.names = TRUE)
  ) #initialize the genome file paths
  # fileNames <- list(teGenomes = list.files(path = "0_raw", pattern = glob2rx("*a"), recursive = TRUE, full.names = TRUE), teRm2Files = list.files("./RM2_output", pattern = "families.fa", recursive = TRUE, full.names = TRUE)) #debugging
  teLibs <- reactiveValues()
  # teLibs <- list() #debugging
  teTables <- reactiveValues()
  observeEvent(input$mchTabs, {
    if(input$mchTabs == "Input data"){
      #restart the objects
      fileNames <- reactiveValues(
        teGenomes = list.files(path = "0_raw", pattern = glob2rx("*a"), recursive = TRUE, full.names = TRUE)
        , teRm2Files = list.files("./RM2_output", pattern = "families.fa", recursive = TRUE, full.names = TRUE)
      )#initialize the genome file paths
      
      teLibs <- reactiveValues()
      
      output$inputData <- renderDT({
        data_with_progress <- teData()
        data_with_progress$Status <- lapply(data_with_progress$Status, create_progress_bar)
        datatable(
          data_with_progress %>% dplyr::select(Species, Strain, Genome_file, Raw_lib, Status),
          escape = FALSE, # Permite HTML en la tabla
          callback = JS(
            "table.on('click', 'td', function() {",
            "  var row = table.row(this).index();", # Obtener índice de la fila
            "  Shiny.setInputValue('row_click', row + 1);", # Enviar fila al servidor (1-indexed)
            "});"
          ),
          options = list(autoWidth = TRUE, dom = "t"), 
          selection = "single" # Permitir selección de una sola celda
          , rownames = FALSE
        ) %>%
          formatStyle(
            "Status",
            target = "row",
            backgroundColor = styleColorBar(c(0, 100), "#4caf50")
          )
      })
    }
  })
  teData <- reactive({
    # rv$teDataRefresh  # React to the trigger value
    # fileNames$teRm2Files <- list.files("./RM2_output", pattern = "families.fa", recursive = TRUE, full.names = TRUE) #debugging
    data <- data.frame(file_name = fileNames$teRm2Files) %>% 
      mutate(
        Species = file_name %>% strsplit("/", fixed = TRUE) %>% lapply("[[", 3) %>% unlist()
        , Strain = file_name %>% strsplit("/", fixed = TRUE) %>% lapply("[[", 4) %>% unlist()
        , Raw_lib = file_name %>% strsplit("/", fixed = TRUE) %>% lapply("[[", 5) %>% unlist()
        , Genome_file = sapply(Strain, function(s) grep(s, x = fileNames$teGenomes, value = TRUE))
        # , Status = lapply(1:nrow(data), function(i) getMchStatus(data$Species[i], data$Strain[i], data$Raw_lib[i], outDir = "./MCHelper")) %>% unlist()
      )
    data$Status <- lapply(1:nrow(data), function(i) getMchStatus(data$Species[i], data$Strain[i], data$Raw_lib[i], outDir = "./MCHelper")) %>% unlist()
    
    data
  })
  
  # Select TE lib to continue analyzing----
  observeEvent(input$row_click, {
    # req(input$inputData_cells_selected)
    # fileNames <- list() #debugging
    updateSelectizeInput(session, inputId = "teManualSeq", selected = NULL, choices = NULL)
    selectedRow <- input$inputData_rows_selected
    species <- teData()[selectedRow, "Species"]
    strain <- teData()[selectedRow, "Strain"]
    status <- teData()[selectedRow, "Status"]
    fileNames$teGenome <- teData()[selectedRow, "Genome_file"]
    fileNames$teRepeatLib <- teData()[selectedRow, "file_name"]
    # selectedRow <- 1 #debugging
    # species <- data[selectedRow, "Species"] #debugging
    # strain <- data[selectedRow, "Strain"] #debugging
    # status <- data[selectedRow, "Status"]#debugging
    # fileNames$teGenome <- data[selectedRow, "Genome_file"] #debugging
    # fileNames$teRepeatLib <- data[selectedRow, "file_name"] #debugging
    # # # #
    
    fileNames$teOutDir <- paste("./MCHelper", species, strain, sep = "/")
    # fileNames$teRepeatLib <- paste("./RM2_output", species, strain, "N2_test-families.fa", sep = "/") #debugging
    teLibs$teRepeatLib <- teReadLib(fileNames$teRepeatLib, libIdentifier = "Raw lib")
    fileNames$teCleanLib <- paste0(fileNames$teOutDir, "/", strain, "-clean_families.fa")
    teLibs$teCleanLib <- teReadLib(fileNames$teCleanLib, libIdentifier = "Clean lib")
    
    fileNames$teAutoCuratedLib <- paste0(fileNames$teOutDir, "/curated_sequences_NR.fa")
    teLibs$teAutoCuratedLib <- teReadLib(fileNames$teAutoCuratedLib, libIdentifier = "Curated sequences")
    
    fileNames$teManualOutDir <- paste0(fileNames$teOutDir,"/1_MI_MCH")
    
    fileNames$teCompleteModelsLib <- paste0(fileNames$teManualOutDir, "/complete_models.fa")
    teLibs$teCompleteModelsLib <- teReadLib(fileNames$teCompleteModelsLib, libIdentifier = "Complete models")
    fileNames$teIncompleteModelsLib <- paste0(fileNames$teManualOutDir, "/incomplete_models.fa")
    teLibs$teIncompleteModelsLib <- teReadLib(fileNames$teIncompleteModelsLib, libIdentifier = "Incomplete models")
    # fileNames$teManualSplittedLib <- paste0(fileNames$teManualOutDir, "/splitted_curated_sequences_NR.fa")
    fileNames$teAutoCuratedTmpLib <- paste0(fileNames$teManualOutDir, "/tmp_curated_sequences_NR.fa")
    teLibs$teAutoCuratedTmp <- teReadLib(fileNames$teAutoCuratedTmpLib)
    fileNames$teManualDiscardedLib <- paste0(fileNames$teManualOutDir, "/discarded_sequences.fa")
    
    fileNames$te80Lib <- paste0(fileNames$teManualOutDir, "/complete_808080.fa")
    
    output$teNonCurated <- renderPlotly({
      tePlotLib(teLibList = list(teLibs$teRepeatLib, teLibs$teCleanLib, teLibs$teAutoCuratedLib))#
    })
    # fileNames$te80Lib <- gsub(fileNames$teCompleteModelsLib, pattern = "complete_models.fa", replacement = "complete_808080.fa")
    teLibs$te80Lib <- teReadLib(fileNames$te80Lib, libIdentifier = "808080 lib")
    
    fileNames$teNon80Lib <- gsub(fileNames$teCompleteModelsLib, pattern = "complete_models.fa", replacement = "complete_non808080.fa")
    fileNames$te80OutDir <- gsub(fileNames$teNon80Lib, pattern = "complete_non808080.fa", replacement = "2_MI_MCH")
    fileNames$teNon80TmpLib <- paste0(fileNames$te80OutDir, "/tmp_non808080.fa")
    fileNames$te80RecovLib <- gsub(fileNames$teCompleteModelsLib, pattern = "complete_models.fa", replacement = "complete_808080_recovered.fa")
    teLibs$te80RecovLib <- teReadLib(fileNames$te80RecovLib, libIdentifier = "Complete 808080 recovered")
    fileNames$teStandbyLib <- paste0(fileNames$te80OutDir, "/standby_sequences.fa")
    teLibs$teStandbyLib <- teReadLib(fileNames$teStandbyLib, libIdentifier = "Stand by lib")
    
    fileNames$te70InputDb <- gsub(pattern = "standby_sequences.fa", replacement = "allDatabases.clustered_merged.fa", x = fileNames$teStandbyLib)
    fileNames$te70Lib <- gsub(pattern = "standby_sequences.fa", replacement = "standby_707070.fa", x = fileNames$teStandbyLib)
    teLibs$te70Lib <- teReadLib(fileNames$te70Lib, libIdentifier = "707070 lib")
    fileNames$te70FilterLog <- gsub(fileNames$te70Lib, pattern = "standby_707070.fa", replacement = "707070_filter.log")

    fileNames$teIncompleteNon80 <- gsub(fileNames$teIncompleteModelsLib, pattern = "incomplete_models.fa", replacement = "incomplete_non808080.fa")
    fileNames$teIncompleteNewfam <- gsub(x = fileNames$teIncompleteNon80, pattern = "incomplete_non808080.fa", replacement = "3_MI_MCH/incomplete_newfam.fa")
    fileNames$teIncomopleteRecovLib <- gsub(x = fileNames$teIncompleteNon80, pattern = "incomplete_non808080.fa", replacement = "incomplete_808080_recovered.fa")
    teLibs$teIncompleteRecovLib <- teReadLib(fileNames$teIncomopleteRecovLib, libIdentifier = "Incomplete 808080 recovered")
    fileNames$teIncompleteNon80Tmp <- gsub(pattern = "incomplete_newfam.fa", replacement = "tmp_incomplete_non808080.fa", fileNames$teIncompleteNewfam)
    
    fileNames$te70FiltLib <- gsub(pattern = "standby_707070.fa", replacement = "standby_707070_filtered.fa",fileNames$te70Lib)
    teLibs$te70FiltLib <- teReadLib(fileNames$te70FiltLib, libIdentifier = "Standby 707070 filtered")
    
    fileNames$teIncompleteOutDir <- paste0(fileNames$teManualOutDir, "/3_MI_MCH")
    
    fileNames$teFinalMchDir <- paste0(fileNames$teOutDir, "/MCH_final")
    
    updateTabsetPanel(session, inputId = "mchTabs"
                      , selected = ifelse(grepl(as.character(status), x = "1234"), yes = "1234"
                                          , no = ifelse(grepl(as.character(status), x = "56"), yes = "56", no = as.character(status)))
    )
  }, once = FALSE)
  
  # Upload TE lib data----
  observeEvent(input$teSubmit, {
    req(input$teGenomeFile)
    req(input$teRm2File)
    req(input$teSpecies != "")
    req(input$teStrain != "")
    genomeDir <- paste("./0_raw", input$teSpecies, input$teStrain, sep = "/")
    if(dir.exists(genomeDir)){
      showModal(modalDialog(title = "Genome directory already exists, please check species and strain"))
    } else{
      dir.create(genomeDir, recursive = TRUE)
      genomeSrc <- input$teGenomeFile$datapath # Path temporal del archivo subido
      genomeDest <- file.path(genomeDir, input$teGenomeFile$name) # Nuevo path con el nombre original
      file.copy(genomeSrc, genomeDest) # Copiar el archivo al destino
      file.remove(genomeSrc)
      if(tools::file_ext(genomeDest) == "zip") {unzip(genomeDest, exdir = genomeDir)}
    }
    
    rawTeDir <- paste("./RM2_output", input$teSpecies, input$teStrain, sep = "/")
    if(dir.exists(rawTeDir)){
      showModal(modalDialog(title = "Raw library already exists, please check species and strain"))
    } else{
      dir.create(rawTeDir, recursive = TRUE)
      rawTeSrc <- input$teRm2File$datapath # Path temporal del archivo subido
      rawTeDest <- file.path(rawTeDir, input$teRm2File$name) # Nuevo path con el nombre original
      file.copy(rawTeSrc, rawTeDest) # Copiar el archivo al destino
      file.remove(rawTeSrc)
    }
    output$teUploadSummary <- renderText({
      paste0("Raw TE library uploaded ", rawTeSrc, " and moved to ", rawTeDest
             , "\n\n"
             , "Genome file uploaded ", genomeSrc, " and moved to ", genomeDest)
    })
  })
  
  
  # 
  #############-
  # Main for MCHelper pipeline----
  #############-
  
  ## 1 Automatic module against clean families  from RM2----
  observeEvent(input$teAutoModuleRun, {
    cat("Automatic module executing ----\n")
    fileNames$teAutoMchLog <- paste(fileNames$teOutDir, "mchelper_automatic.log", sep = "/")
    fileNames$teManualMchLog <- gsub(pattern = "automatic", replacement = "manual", x = fileNames$teAutoMchLog)
    if(file.exists(fileNames$teAutoMchLog)){# looks like the automatic module has been run
      cat("Automatic module ran ----\n")
      if(file.exists(fileNames$teAutoCuratedLib)){#if curated sequences library exists
        teLibs$teAutoCuratedLib <- teReadLib(fileNames$teAutoCuratedLib, libIdentifier = "Curated sequences")
        output$teNonCurated <- renderPlotly({
          tePlotLib(teLibList = list(teLibs$teRepeatLib, teLibs$teCleanLib, teLibs$teAutoCuratedLib))
        })
        if(file.exists(fileNames$teManualMchLog)){# if manual module has been executed
          logFile <- readLines(fileNames$teManualMchLog)
          
          #if MCHelper manual has ended, log file will have a particular string in the last line
          if(logFile[length(logFile)] == "MCHelper with TE-Aid successfully run"){
            output$teManualInputLib <- renderText({fileNames$teAutoCuratedLib})#automatic curated lib under inspection
            if(file.exists(fileNames$teCompleteModelsLib)){
              cat("Esta ejecutando como si existiera la libreria completa----\n")
              teLibs$teAutoCuratedTmp <- Biostrings::readDNAStringSet(fileNames$teAutoCuratedTmpLib)
              updateSelectizeInput(session, inputId = "teManualSeq", choices = names(teLibs$teAutoCuratedTmp))
                    if(length(teLibs$teAutoCuratedTmp) == 0) {
                      showModal(modalDialog(
                        title = "Looks like you already reviewed the automatic curated TE models!!",
                        # div(textOutput(outputId = "teMessage"))
                        paste0("There are no sequences to analyze in selected library: "
                               , gsub(pattern = "curated_sequences_NR.fa", replacement = "1_MI_MCH/tmp_curated_sequences_NR.fa", fileNames$teAutoCuratedLib)
                               , "\nPlease change your input library for Manual Inspection models"
                        )
                      ))
                      updateTabsetPanel(session, inputId = "mchTabs", selected = "56")
                    }
            } else{#if no manual inspection has been done yet
              cat("La libreria temporal automatica se carga por primera vez ----\n")
              teLibs$teAutoCuratedTmp <- Biostrings::readDNAStringSet(fileNames$teAutoCuratedLib)
              # teLibs$teCompleteModelsLib <- Biostrings::DNAStringSet() #initialize complete models set
              # teLibs$teManualSplittedLib <- Biostrings::DNAStringSet()
              # Biostrings::writeXStringSet(teLibs$teManualSplittedLib, filepath = fileNames$teManualSplittedLib, append = FALSE) #append false to ensure a new lib is created
              updateSelectizeInput(session, inputId = "teManualSeq", choices = names(teLibs$teAutoCuratedTmp))
              # input$teManualSeq <- names(teLibs$teAutoCuratedTmp)[1]#for debugging
            }
            
            # teTables$manualCuratedClassif <- read.table(paste0(fileNames$teManualOutDir, "/denovoLibTEs_PC.classif"), sep = "\t", header = TRUE)
            # fileNames$te_aidGenomeBnFiles <- list.files(fileNames$teManualOutDir, pattern = "genome.blastn", recursive = TRUE, full.names = TRUE)
            # fileNames$te_aidSelfBnFiles <- list.files(fileNames$teManualOutDir, pattern = "self.blastn", recursive = TRUE, full.names = TRUE)
          } else{#if MCHelper manual has not ended, show a message
            cat("MCHelper automatic module still running---\n", fileNames$teAutoCuratedLib, "\n")#show modeldialog
            showModal(modalDialog(
              title = "Looks like MCHelper Manual module is still running!!",
              paste0("Wait until MCHelper processes: "
                     , fileNames$teAutoCuratedLib
                     , "\n or please change your input library for Manual Inspection models"
              )
            ))
            # rv$teDataRefresh <- rv$teDataRefresh + 1
          }
        }
        paste(c("MCHelper automatic seems to have been run\n" #it may be running too
                , readLines(fileNames$teAutoMchLog)), collapse = "\n")
      } else{#automatic module has not ended
        showModal(modalDialog(
          title = "MCHelper automatic seems to be still running"
          , paste("MCHelper automatic log exists but curated sequences library is not found, looks like automatic module has not ended yet. LogFile:", fileNames$teAutoMchLog)
        ))
        paste(c("MCHelper automatic seems to be still running \n"
                , readLines(fileNames$teAutoMchLog)), collapse = "\n")
      }
      #     
    } else{#run automatic module
      cat("########\n## Indexing genome ", fileNames$teGenome, "\n########")
      # cat("Creating database index----\n")
      system2(
        "conda"
        , args = c(
          "run -n MCHelper makeblastdb -in"
          , fileNames$teGenome
          , "-dbtype nucl"
        )
      )
      cat("########\n##Running automatic module for ", fileNames$teCleanLib, "\n########")
      system2(
        "conda",
        args = c(
          "run",
          "-n", "MCHelper",
          "--no-capture-output",
          "python3",
          "mchelper-ats/MCHelper.py",
          "-r", "A",
          "-a", "F",
          "--input_type", "fasta",
          "-l", fileNames$teCleanLib, 
          "-g", fileNames$teGenome #fileNames$teGenomes[grepl(pattern = strain, x = fileNames$teGenomes)]
          ,"-o", fileNames$teOutDir,
          "-c", "1",
          "-t", input$teAutoModuleCores,
          "-b", paste0("BUSCO_libs/", input$teBuscoLib)
          , "-v", "Y >", fileNames$teAutoMchLog #verbose output to be able to monitor the script
          , "&& conda run -n MCHelper --no-capture-output python3",
          "mchelper-ats/MCHelper.py"
          ,"-r", "T" #use T module as it only runs TE-Aid that is what we want in this part
          , "--input_type", "fasta",
          "-l", fileNames$teAutoCuratedLib,
          "-g", fileNames$teGenome
          , "-o", fileNames$teManualOutDir,
          "-t", input$teAutoModuleCores
          , "-v", "Y > ", fileNames$teManualMchLog
        )
        # , stdout = fileNames$teAutoMchLog
        # , stderr = fileNames$teAutoMchLog
        , wait = FALSE
      )## esto genera el archivo curated_sequences_NR.fa que usaremos para revisar manualmente
    }
    output$teManualStats <- renderTable({
      # teList <- teLibs #debugging
      teList <- reactiveValuesToList(teLibs)
      teCurationStats(teLibList = teList[grep(names(teList), pattern = "Repeat|Clean|Auto", value = TRUE)])#list(teLibs$teRepeatLib, teLibs$teCleanLib, teLibs$teAutoCuratedLib, teLibs$teAutoCuratedTmp)
    }, rownames = TRUE)
  })
  
  observeEvent(input$teManualSeq, {
    # input$teManualSeq <- names(teLibs$teAutoCuratedTmp)[217] #debugging
    req(input$teManualSeq)
    clickData$click1 <- NULL
    clickData$click2 <- NULL
    teAnnotation <- input$teManualSeq %>% strsplit("#") %>% lapply("[[", 2) %>% unlist()
    updateSelectInput(session, inputId = "teDecision"
                      # , choices = c("Hola", "Hola3", "Hols3", teAnnotation)#this should be changed to the full table of annotations (42 entries?)
                      , selected = teAnnotation)
    seqName <- input$teManualSeq %>% strsplit("#", fixed = TRUE) %>% unlist() %>% head(1)
    teSeq <- teLibs$teAutoCuratedTmp[input$teManualSeq]
    
    # input$teDecision <- teAnnotation #for debugging
    # fileNames$te_aidGenomeBnFiles <- list.files(fileNames$teManualOutDir, pattern = glob2rx(paste0(seqName, ".fa*genome.blastn*")), recursive = TRUE, full.names = TRUE) #%>% tail(1)
    # fileNames$te_aidSelfBnFiles <- list.files(fileNames$teManualOutDir, pattern = glob2rx(paste0(seqName, ".fa*self.blastn*")), recursive = TRUE, full.names = TRUE) #%>% tail(1)
    # teTable <- teTables$manualCuratedClassif
    teTable <- read.table(paste0(fileNames$teManualOutDir, "/denovoLibTEs_PC.classif"), sep = "\t", header = TRUE)
    teBlastGenome <- list.files(fileNames$teManualOutDir, pattern = glob2rx(paste0(seqName, ".fa*genome.blastn*")), recursive = TRUE, full.names = TRUE) #%>% tail(1)
    teBlastSelf <- list.files(fileNames$teManualOutDir, pattern = glob2rx(paste0(seqName, ".fa*self.blastn*")), recursive = TRUE, full.names = TRUE) #%>% tail(1)
    
    
    updateNumericInput(inputId = "leftTrim", value = 1, max = Biostrings::width(teSeq)) # right trim cannot be higher than the length of the sequence
    updateNumericInput(inputId = "rightTrim", value = Biostrings::width(teSeq), max = Biostrings::width(teSeq)) # right trim cannot be higher than the length of the sequence
    
    output$teManualCuratedClass <- renderTable({
      teTable %>% subset(Seq_name == seqName) %>% select(Seq_name, length, class, order, sFamily, struct, coding, other, strand)
    })
    output$teManualClass <- renderTable(
      teClassOutput(teTable, seqName)$teInfo
    )
    output$teManualCoding <- renderTable(
      teClassOutput(teTable, seqName)$teCoding, colnames = FALSE
    )
    output$teManualStructure <- renderTable(
      teClassOutput(teTable, seqName)$teStruc
    )
    
    
    output$teConsDivergence <- renderPlot({
      teConsDivergencePlot(seqName, blastFiles = teBlastGenome, leftTrim = input$leftTrim, rightTrim = input$rightTrim)
    })
    output$teConsCoverage <- renderPlotly({
      teConsCoveragePlot(seqName, blastFiles = teBlastGenome, leftTrim = input$leftTrim, rightTrim = input$rightTrim, source = "A")
    })
    output$teConsSelf <- renderPlot({
      teConsSelfPlot(seqName, blastFiles = teBlastSelf, leftTrim = input$leftTrim, rightTrim = input$rightTrim)
    })
    output$teConsStructure <- renderPlot({
      teConsStructurePlot(seqName, blastFiles = teBlastSelf, leftTrim = input$leftTrim, rightTrim = input$rightTrim)
    })
    
  })
  
  ######-
  #### Capture for trimming sequences----
  ######-
  clickData <- reactiveValues(
    click1 = NULL, #this should be x=0
    click2 = NULL # this should be x=length(consensus)
  )
    # clickData$click1 <- list(x = 1, y = 2) #debugging
    # clickData$click2 <- list(x = 20, y = 2)
  # Capturar clics del grafico
  observeEvent(event_data("plotly_click", source = "A"), {
    click <- event_data("plotly_click", source = "A")

    if (is.null(clickData$click1)) {
      # Primer clic
      clickData$click1 <- list(x = click$x, y = click$y)
      # Actualiza los valores de los numericInput
      if(is.null(clickData$click2)){#this is if it is the first click
        updateNumericInput(session, "rightTrim", value = max(clickData$click1$x, clickData$click2$x))
      } else{
        updateNumericInput(session, "leftTrim", value = min(clickData$click1$x, clickData$click2$x))
        updateNumericInput(session, "rightTrim", value = max(clickData$click1$x, clickData$click2$x))
      }

      # nulificamos el segundo click para que escuche de nuevo el click
      clickData$click2 <- NULL
      # Actualiza los valores de los numericInput
    } else if (is.null(clickData$click2)) {
      # Segundo clic
      clickData$click2 <- list(x = click$x, y = click$y)
      # nulificamos el primer click para que escuche de nuevo el click
      # Actualiza los valores de los numericInput
      updateNumericInput(session, "rightTrim", value = max(clickData$click1$x, clickData$click2$x))
      updateNumericInput(session, "leftTrim", value = min(clickData$click1$x, clickData$click2$x))
      clickData$click1 <- NULL
    }

  })
  # output$debug <- renderText({
  #   # event <- event_data("plotly_relayout", source = "A")
  #   # if (!is.null(event$xaxis.range[1]) && !is.null(event$xaxis.range[2])) {
  #   #   paste("Gráfico 1 - Límites actuales: [",
  #   #         event$xaxis.range[1], ", ", event$xaxis.range[2], "]")
  #   #   updateNumericInput(session, "rightTrim", value = event$xaxis.range[2])
  #   #   updateNumericInput(session, "leftTrim", value = event$xaxis.range[1])
  #   #   # input$leftTrim <- event$xaxis.range[1]
  #   #   # input$rightTrim <- event$xaxis.range[2]
  #   # } else {
  #   #   "No hay interacción aún en Gráfico 1."
  #   # }
  #   event <- event_data("plotly_relayout", source = "A")
  #   if (!is.null(event[["xaxis.range[0]"]]) && !is.null(event[["xaxis.range[1]"]])) {
  #     paste("Gráfico 1 - Límites actuales: [",
  #           event[["xaxis.range[0]"]], ", ", event[["xaxis.range[1]"]], "]")
  #   } else {
  #     # "No hay interacción aún en Gráfico 1."
  #     paste("Límites actuales: [",
  #           input$teConsCoverage_brush$x, ", ", input$teConsCoverage_brush$x, "]")
  #   }
  # })
  
  #### Reannotating TE ------
  observeEvent(input$teAnnotate, {#if no change annotation, new annotation == old annotation
    req(length(teLibs$teAutoCuratedTmp) != 0)
    output$teTrimError <- NULL #to remove message of trimming
    seqName <- input$teManualSeq
    teSeq <- teLibs$teAutoCuratedTmp[seqName]
    newAnnotation <- paste(
      seqName %>% strsplit("#", fixed = TRUE) %>% lapply("[[", 1) %>% unlist()
      , input$teDecision
      , sep = "#"
    )
    #re-naming in manual.v2
    names(teSeq) <- newAnnotation
    if(grepl(pattern = "_inc|_unconfirmed", x = names(teSeq))) {
      #it is an incomplete or unconfirmed model
      cat("Writting incomplete model curated library----\n")
      Biostrings::writeXStringSet(teSeq, file = fileNames$teIncompleteModelsLib, append = TRUE)
      # teLibs$teCompleteModelsLib[newAnnotation] <- NULL
    }else{
      #it is a complete model
      cat("Writting complete model curated library----\n")
      Biostrings::writeXStringSet(teSeq, file = fileNames$teCompleteModelsLib, append = TRUE)
    }
    
    #remove original seq from automatic annotated lib
    teLibs$teAutoCuratedTmp[input$teManualSeq] <- NULL
    cat("Writting temporary analyzed library----\n")
    Biostrings::writeXStringSet(teLibs$teAutoCuratedTmp, file = fileNames$teAutoCuratedTmpLib)
    clickData$click1 <- NULL
    clickData$click2 <- NULL
    
    #update the TE seq to analyze
    updateSelectizeInput(session, inputId = "teManualSeq", choices = names(teLibs$teAutoCuratedTmp))
    
    if(length(teLibs$teAutoCuratedTmp) == 0) {
      updateSelectizeInput(session, inputId = "te80InputLib"
                           , choices = list.files(pattern = "^complete_models.fa", recursive = TRUE, full.names = TRUE)
                           , server = TRUE)
      showModal(modalDialog(
        title = "You finished reviewing the automatic curated TE models!!",
        # div(textOutput(outputId = "teMessage"))
        paste0("There are no more sequences to analyze in selected library: "
               , fileNames$teAutoCuratedTmpLib
               , "\nPlease change your input library for Manual Inspection models"
        )
      ))
      #update to 808080 models tab
      updateTabsetPanel(session, inputId = "mchTabs", selected = "5")
      # file.remove(fileNames$teAutoCuratedTmpLib)
    }
  }, once = FALSE) #eliminate the observe for the re-annotation button
  
  #### Trim TE ------
  observeEvent(input$teTrim, {
    fileNames$teManualTrimDir <- paste0(fileNames$teManualOutDir, "/trimming")
    if(dir.exists(fileNames$teManualTrimDir)){
      cat(fileNames$teManualTrimDir, " already existing ----\n")
    } else{cat(fileNames$teManualTrimDir, " created directory ---\n"); dir.create(fileNames$teManualTrimDir)}
    
    # trimming
    # teLibs$teAutoCuratedTmp <- Biostrings::readDNAStringSet(fileNames$teAutoCuratedTmpLib) #debugging
    # input$teManualSeq <- names(teLibs$teAutoCuratedTmp[1])#debugging
    seqName <- input$teManualSeq %>% strsplit("#", fixed = TRUE) %>% unlist() %>% head(1)
    # teSeq <- teLibs$teAutoCuratedTmp[input$teManualSeq] #we switch to original sequence for being able to trim it again.
    teSeq <- teLibs$teAutoCuratedLib[input$teManualSeq]
    if((input$leftTrim>=1) & (input$rightTrim > input$leftTrim) & (input$rightTrim<=Biostrings::width(teSeq))){#if trimming parameters are ok
      req(length(teLibs$teAutoCuratedTmp) != 0)
      cat("Trimming sequence-----\n")
      teSeq <- subseq(teSeq, start = input$leftTrim, end = input$rightTrim)
      seqFa <- paste0(fileNames$teManualTrimDir, "/", seqName,".fa")
      if(file.exists(seqFa)){
        cat(seqFa, " already trimmed ----\n")
        system(paste0("rm -rd ", seqFa, "*"))}
      Biostrings::writeXStringSet(teSeq, filepath = seqFa, append = FALSE)
      cat("Running MCHelper TE-aid for seq", seqName, " ----\n")
      system2(
        "conda"
        , args = c(
          "run -n MCHelper --no-capture-output python3 mchelper-ats/MCHelper.py"
          ,"-r", "T" #use T module as it only runs TE-Aid that is what we want in this part
          , "--input_type", "fasta",
          "-l", seqFa,
          "-g", fileNames$teGenome
          , "-o", paste0(seqFa, "_MCH"),
          "-t", 1
        )
        , wait = FALSE
      )
      write.table(data.frame(leftTrim = input$leftTrim, rightTrim = input$rightTrim), row.names = FALSE
                  , file = paste0(seqFa, "_trimPositions.tsv"), sep = "\t")
      showModal(modalDialog(
        title = paste("Running MCHelper TE-aid for seq", seqName, " ----\n")
        , "Please continue with other sequences while MCHelper is running for the trimmed sequence"
      ))
      #remove original seq from automatic annotated lib
      teLibs$teAutoCuratedTmp[input$teManualSeq] <- NULL
      cat("Writting temporary analyzed library----\n")
      Biostrings::writeXStringSet(teLibs$teAutoCuratedTmp, file = fileNames$teAutoCuratedTmpLib)
      
      #append the trimmed seq to the temporary library
      Biostrings::writeXStringSet(teSeq, filepath = fileNames$teAutoCuratedTmpLib, append = TRUE)
      teLibs$teAutoCuratedTmp <- c(teLibs$teAutoCuratedTmp, teSeq)
      clickData$click1 <- NULL
      clickData$click2 <- NULL
      #update the TE seq to analyze
      updateSelectizeInput(session, inputId = "teManualSeq", choices = names(teLibs$teAutoCuratedTmp))
      # paste("conda run -n MCHelper --no-capture-output python3 mchelper-ats/MCHelper.py -r T --input_type fasta -l"
      #       , seqFa, "-g", fileNames$teGenome, "-o")
    } else{
      cat("Trimming parameters not met----\n")
      output$teTrimError <- renderText("Trimming parameters are not correct. Please check them")
      showModal(modalDialog(
        title = "¡Trimming parameters are not met!",
        "You are trying to trim a sequence larger thant the sequence itself. Please review your trimming parameters and re-trim"
      ))
    }
  })
  
  #### Discard TE model --------
  observeEvent(input$teDiscard, {
    req(length(teLibs$teAutoCuratedTmp) != 0)
    teSeq <- teLibs$teAutoCuratedTmp[input$teManualSeq]
    cat("Writting discarded library----\n")
    Biostrings::writeXStringSet(teSeq, file = fileNames$teManualDiscardedLib, append = TRUE)
    #remove original seq from automatic annotated lib
    teLibs$teAutoCuratedTmp[input$teManualSeq] <- NULL
    cat("Writting temporary analyzed library----\n")
    Biostrings::writeXStringSet(teLibs$teAutoCuratedTmp, file = fileNames$teAutoCuratedTmpLib)
    updateSelectizeInput(session, inputId = "teManualSeq", choices = names(teLibs$teAutoCuratedTmp))
    if(length(teLibs$teAutoCuratedTmp) == 0) {
      updateSelectizeInput(session, inputId = "te80InputLib"
                           , choices = list.files(pattern = "^complete_models.fa", recursive = TRUE, full.names = TRUE)
                           , server = TRUE)
      showModal(modalDialog(
        title = "You finished reviewing the automatic curated TE models!!",
        # div(textOutput(outputId = "teMessage"))
        paste0("There are no more sequences to analyze in selected library: "
               , fileNames$teAutoCuratedTmpLib
               , "\nPlease change your input library for Manual Inspection models"
        )
      ))
      #update to 808080 models tab
      updateTabsetPanel(session, inputId = "mchTabs", selected = 5)
      # file.remove(fileNames$teAutoCuratedTmpLib)
    }
  })
  
  output$te80InputLib <- renderText({fileNames$teCompleteModelsLib})
  output$te80InputPlot <- renderPlotly({
    tePlotLib(teLibList = list(teLibs$teAutoCuratedLib, teLibs$teCompleteModelsLib, teLibs$te80Lib))
  })
  
  # 3 Blastn 80-80-80 module for complete models ----
  observeEvent(input$te80ModuleRun,{
    # req(teLibs$teAutoCuratedTmp)
    if(length(teLibs$teAutoCuratedTmp) != 0){
      showModal(modalDialog(
        title = "Please finish with the Manual Inspection module before proceeding!"
        , paste0("You must have reviewed all the sequences before going to 80-80-80 model. There are still sequences in the curated library to review manually: \n"
                 , fileNames$teAutoCuratedTmpLib
        )
      ))
    } else{#if the curated library has been completely reviewed
      cat("The temporary curated library is empty, executing the 808080 module\n")
  # #     
      strain <- fileNames$teCompleteModelsLib %>% strsplit("/", fixed = TRUE) %>% lapply("[[", 4) %>% unlist() #"N2" #input$teStrain
  # #     # species <- fileNames$te80InputLib %>% strsplit("/", fixed = TRUE) %>% lapply("[[", 3) %>% unlist() #"C. elegans"#input$teSpecies
      fileNames$te80BlastLog <- gsub(fileNames$teCompleteModelsLib, pattern = "complete_models.fa", replacement = "blast808080.log")
      fileNames$te80Blast <- gsub(fileNames$teCompleteModelsLib, pattern = "complete_models.fa", replacement = "complete_models_vs_allDBs.blast")
      if(file.exists(fileNames$te80BlastLog)){
        #if blast 808080 outputs are present, raise a warning
        warning(paste0("Looks like blast agains TE database has been already executed for library ", fileNames$teCompleteModelsLib, ", please check\n"))
  # #       logFileBlast <- readLines(fileNames$te80BlastLog)
      }else{#if blast 808080 outputs are not present, do the blast
        cat("Performing the blast --------------- \n")
        system2(
          "conda",
          args = c(
            "run"
            , "-n", "MCHelper"
            , "--no-capture-output"
            , "blastn"
            ,"-query", fileNames$teCompleteModelsLib
            ,"-db", "mchelper-ats/db/allDatabases.clustered_rename.fa"
            , "-num_threads", input$teAutoModuleCores
            , "-outfmt 6 -qcov_hsp_perc 80 -perc_identity 80 -max_hsps 1"
            , "-out", fileNames$te80Blast
          )
          , stdout = fileNames$te80BlastLog
          , stderr = fileNames$te80BlastLog
          , wait = TRUE
        )
  # #       # conda run -n MCHelper --no-capture-output blastn -query MCHelper_automatic/C.elegans/N2-test/1_MI_MCH/complete_models.fa -db mchelper-ats/db/allDatabases.clustered_rename.fa -num_threads 32 -outfmt 6 -qcov_hsp_perc 80 -perc_identity 80 -max_hsps 1 -out MCHelper_automatic/C.elegans/N2-test/1_MI_MCH/N2-test_vs_allDBs.blast
        cat("Blast successfully run --------------- \n")
  # #       logFileBlast <- fileNames$te80BlastLog
      }
  # #     
      fileNames$te80AssignLog <- gsub(fileNames$teCompleteModelsLib, pattern = "complete_models.fa", replacement = "assign_families.log")
      fileNames$te80Assign <- gsub(fileNames$teCompleteModelsLib, pattern = "complete_models.fa", replacement = "")
      if(file.exists(fileNames$te80AssignLog)){
  # #       #if assign_families outputs are present, raise a warning
        warning(paste0("Looks like assign_families has been already executed for library", fileNames$teCompleteModelsLib, ", please check\n"))
  # #       logFileAssign <- readLines(fileNames$te80AssignLog)
      }else{#if assign_families has not been executed, execute it
        cat("Running assign_families --------------- \n")
        system2(
          "conda"
          , args = c(
            "run"
            , "-n", "MCHelper"
            , "--no-capture-output"
            , "python3", "mchelper-ats/tools/Curation_toolkit/curation_toolkit.py"
            , "assign_families" #curation toolkit module
            , strain #strainID
            , fileNames$teCompleteModelsLib #complete TE library
            , fileNames$te80Blast #blast result
            , gsub(fileNames$teCompleteModelsLib, pattern = "complete_models.fa", replacement = "complete") # output files prefix
          )
          , stdout = fileNames$te80AssignLog
          , stderr = fileNames$te80AssignLog
          , wait = TRUE
        )
  # #       # conda run -n MCHelper --no-capture-output python3 mchelper-ats/tools/Curation_toolkit/curation_toolkit.py assign_families N2-test MCHelper_automatic/C.elegans/N2-test/1_MI_MCH/complete_models.fa MCHelper_automatic/C.elegans/N2-test/1_MI_MCH/N2-test_vs_allDBs.blast MCHelper_automatic/C.elegans/N2-test/1_MI_MCH/N2-test_complete
        cat("Assign families successfully run --------------- \n")
  # #       # system("# copy 808080.fa file to the 808080_recovered.fa") #let's do it in the mchelper curation_toolkit side
  # #       logFileAssign <- readLines(fileNames$te80AssignLog)
      }
  # #     
  # #     ### Run MCHelper on non808080.fa
      fileNames$te80MchLog <- gsub(fileNames$teCompleteModelsLib, pattern = "complete_models.fa", replacement = "mchelper_manual.log")
      output$teNon80Lib <- renderText(fileNames$teNon80Lib)
      if(!file.exists(fileNames$te80MchLog)){# if manual module has not run (does the log file exist?), run it
        cat("Running MCHelper in T mode for non808080 sequence, please be patient...\nInput library:", fileNames$teNon80Lib)
        system2(
          "conda"
          , args = c(
            "run",
            "-n", "MCHelper",
            "--no-capture-output",
            "python3",
            "mchelper-ats/MCHelper.py"
            # "-a", "F"
            ,"-r", "T" #use T module as it only runs TE-Aid that is what we want in this part
            , "--input_type", "fasta",
            "-l", fileNames$teNon80Lib
            , "-g", fileNames$teGenome
            , "-o", fileNames$te80OutDir,
            "-t", input$teAutoModuleCores
            , "-v", "Y"
          )
          , stdout = fileNames$te80MchLog
          , stderr = fileNames$te80MchLog
          , wait = FALSE #avoid frozen app while running MCHelper (send process to background)
        )
        cat("MCHelper manual inspection module to background for ", fileNames$teNon80Lib, "----\n")
  # #       logFile80 <- "MCHelper manual inspection module sent to background"
      } else{#looks like the automatic module has been run
        #   #if manual curation stamp is present
        if(file.exists(fileNames$teNon80TmpLib)){
          cat("There is a temporal stamp for manual curation of incomplete models ----\n")
          teLibs$teNon80TmpLib <- Biostrings::readDNAStringSet(fileNames$teNon80TmpLib)
  # #         # # teLibs$non808080Manual <- Biostrings::readDNAStringSet(gsub(fileNames$te80InputLib, pattern = "complete_models.fa", replacement = "2_MI_MCH/complete_newfam.fa"))
  # #         # teLibs$non808080StandbyLib <- Biostrings::readDNAStringSet(gsub(fileNames$te80InputLib, pattern = "complete_models.fa", replacement = "2_MI_MCH/standby_complete.fa"))
          updateSelectizeInput(session, inputId = "te80Seq", choices = names(teLibs$teNon80TmpLib))
          # # input$te80Seq <- names(teLibs$non808080Tmp)[1] #debugging
          if(length(teLibs$teNon80TmpLib) == 0){
            showModal(modalDialog(
              title = "Looks like you already reviewed the non808080 TE models! Moving to 70-70-70 models",
              paste0("There are no sequences to analyze in selected library: "
                     , fileNames$teNon80TmpLib
                     , "\nPlease change your input library for 808080 models"
                     )
              ))
            updateTabsetPanel(session, inputId = "mchTabs", selected = "7")
          }
        }else{#if no manual review has been saved, initialize temporary DNAstringsets
          cat("There is NOT a temporal stamp for manual curation of incomplete models ----\n")
          teLibs$teNon80TmpLib <- Biostrings::readDNAStringSet(fileNames$teNon80Lib)
  # #         teLibs$teStandbyLib <- Biostrings::DNAStringSet()
  # #         Biostrings::writeXStringSet(teLibs$teStandbyLib, filepath = fileNames$teStandbyLib, append = FALSE) #append false to ensure a new file is created
          updateSelectizeInput(session, inputId = "te80Seq", choices = names(teLibs$teNon80TmpLib))
        }
        teTables$te80CuratedClassif <- read.table(paste0(fileNames$te80OutDir,"/denovoLibTEs_PC.classif"), sep = "\t", header = TRUE)
        # fileNames$te_aid80GenomeBnFiles <- list.files(fileNames$te80OutDir, pattern = "genome.blastn", recursive = TRUE, full.names = TRUE)
        # fileNames$te_aid80SelfBnFiles <- list.files(fileNames$te80OutDir, pattern = "self.blastn", recursive = TRUE, full.names = TRUE)
  # #       logFile80 <- readLines(fileNames$te80MchLog)
        #   if(input$teManualShowLog){paste(c("MCHelper manual TE-aid seems to have been run, please double check it\n", logFile), collapse = "\n") #it may be running too}
        if(length(teLibs$teNon80TmpLib) == 0){
          showModal(modalDialog(
            title = "There are no sequences in the non808080 library, maybe all of them have a good model, congrats!!",
            paste0("You can manually check the non 80-80-80 file:"
                   , fileNames$teNon80Lib
                   , "or even check the complete models library:", fileNames$teCompleteModelsLib
                   , "\nPlease change your complete sequences input library to manually inspect the non 80-80-80 models"
            )
          ))
          updateTabsetPanel(session, inputId = "mchTabs", selected = "7")
        }
      }
    }#end of execute 808080 model
    
    output$te80Stats <- renderTable({
      # teList <- teLibs #debugging
      teList <- reactiveValuesToList(teLibs)
      teCurationStats(teLibList = teList[grep(names(teList), pattern = "Complete|80Lib|80Tmp|Stand", value = TRUE)])#list(teLibs$teRepeatLib, teLibs$teCleanLib, teLibs$teAutoCuratedLib, teLibs$teAutoCuratedTmp)
    }, rownames = TRUE)
  # #   # 
  # #   # 
  # #   paste(c("Blast results --------"
  # #           , fileNames$te80BlastLog
  # #           , logFileBlast
  # #           ,"\n\n Assign family results -----"
  # #           , fileNames$te80AssignLog
  # #           , logFileAssign
  # #           , "\n\n MCHelper manual inspection module ----"
  # #           , logFile80
  # #   ), collapse = "\n"
  # #   )
  })#end of run te80Module

  ######-
  #### Plotting non808080 sequences-------
  ######-
  observeEvent(input$te80Seq, {
    # input$te80Seq <- names(teLibs$teNon80TmpLib)[1] #for debugging
    req(input$te80Seq)
    seqName <- input$te80Seq %>% strsplit("#") %>% unlist() %>% head(1)
    te80Annotation <- input$te80Seq %>% strsplit("#") %>% lapply("[[", 2) %>% unlist()
    updateSelectInput(session, inputId = "te80Decision"
                      # , choices = c("Hola", "Hola3", "Hols3", teAnnotation)#this should be changed to the full table of annotations (42 entries?)
                      , selected = te80Annotation)
    
    # input$te80Decision <- teAnnotation #for debugging
    fileNames$te_aid80GenomeBnFiles <- list.files(fileNames$te80OutDir, pattern = glob2rx(paste0(seqName, ".fa*genome.blastn*")), recursive = TRUE, full.names = TRUE) %>% tail(1)
    fileNames$te_aid80SelfBnFiles <- list.files(fileNames$te80OutDir, pattern = glob2rx(paste0(seqName, ".fa*self.blastn*")), recursive = TRUE, full.names = TRUE) %>% tail(1)
    teBlastSelf <- fileNames$te_aid80SelfBnFiles
    teBlastGenome <- fileNames$te_aid80GenomeBnFiles
    
    
    updateNumericInput(inputId = "leftTrim80", value = 1, max = Biostrings::width(teLibs$teNon80TmpLib[input$te80Seq])) # right trim cannot be higher than the length of the sequence
    updateNumericInput(inputId = "rightTrim80", value = Biostrings::width(teLibs$teNon80TmpLib[input$te80Seq]), max = Biostrings::width(teLibs$teNon80TmpLib[input$te80Seq])) # right trim cannot be higher than the length of the sequence
    output$te80CuratedClass <- renderTable({
      teTables$te80CuratedClassif %>% subset(Seq_name == seqName) %>% select(Seq_name, length, class, order, sFamily, struct, coding, other, strand)
    })
    output$te80ConsDivergence <- renderPlot({
      teConsDivergencePlot(seqName, blastFiles = teBlastGenome, leftTrim = input$leftTrim80, rightTrim = input$rightTrim80)
    })
    output$te80ConsCoverage <- renderPlotly({
      teConsCoveragePlot(seqName, blastFiles = teBlastGenome, leftTrim = input$leftTrim80, rightTrim = input$rightTrim80, source = "B")
    })
    output$te80ConsSelf <- renderPlot({
      teConsSelfPlot(seqName, blastFiles = teBlastSelf, leftTrim = input$leftTrim80, rightTrim = input$rightTrim80)
    })
    output$te80ConsStructure <- renderPlot({
      teConsStructurePlot(seqName, blastFiles = teBlastSelf, leftTrim = input$leftTrim80, rightTrim = input$rightTrim80)
    })
  })
  
    ######-
    #### Capture for trimming sequences-------
    ######-
    clickData2 <- reactiveValues(
      click1 = NULL,
      click2 = NULL
    )
    # clickData$click1 <- list(x = 1, y = 2) #debugging
    # clickData$click2 <- list(x = 20, y = 2) #debugging
    # Capturar clics del grafico
    observeEvent(event_data("plotly_click", source = "B"), {
      click <- event_data("plotly_click", source = "B")

      if (is.null(clickData2$click1)) {
        # Primer clic
        clickData2$click1 <- list(x = click$x, y = click$y)
        # Actualiza los valores de los numericInput
        if(is.null(clickData2$click2)){#this is if it is the first click
          #initially we assign the click to the right trim, just because we need to assign it to any position, if not the plot will not show the line.
          updateNumericInput(session, "rightTrim80", value = max(clickData2$click1$x, clickData2$click2$x))
        } else{
          updateNumericInput(session, "rightTrim80", value = max(clickData2$click1$x, clickData2$click2$x))
          updateNumericInput(session, "leftTrim80", value = min(clickData2$click1$x, clickData2$click2$x))
        }

        # nulificamos el segundo click para que escuche de nuevo el click
        clickData2$click2 <- NULL
        # Actualiza los valores de los numericInput
      } else if (is.null(clickData2$click2)) {
        # Segundo clic
        clickData2$click2 <- list(x = click$x, y = click$y)
        # nulificamos el primer click para que escuche de nuevo el click
        # Actualiza los valores de los numericInput
        updateNumericInput(session, "rightTrim80", value = max(clickData2$click1$x, clickData2$click2$x))
        updateNumericInput(session, "leftTrim80", value = min(clickData2$click1$x, clickData2$click2$x))
        clickData2$click1 <- NULL
      }
    })
    #### Reannotating and Trim TE ------
    observeEvent(input$te80Annotate, {#if no change annotation, new annotation == old annotation
      req(length(teLibs$teNon80TmpLib) > 0)
      seqName <- input$te80Seq
      teSeq <- teLibs$teNon80TmpLib[seqName]
  #     strain <- input$te80InputLib %>% strsplit("/", fixed = TRUE) %>% lapply("[[", 4) %>% unlist()
      # input$te80Decision <- "CLASSI/Prueba" # debugging
      newAnnotation <- paste(
        seqName %>% strsplit("#", fixed = TRUE) %>% lapply("[[", 1) %>% unlist()
        , input$te80Decision
        , sep = "#"
      )
  #     #re-naming in manual.v2
      names(teSeq) <- newAnnotation
      # trimming
      if((input$leftTrim80>=1) & (input$rightTrim80 > input$leftTrim80) & (input$rightTrim80<=Biostrings::width(teLibs$teNon80TmpLib[seqName]))){#if trimming parameters are ok
        cat("Trimming sequence-----\n")
        teSeq <- subseq(teSeq[newAnnotation], start = input$leftTrim80, end = input$rightTrim80)
        clickData$click1 <- NULL
        clickData$click2 <- NULL

        cat("Writting newfamilies model curated library----\n")
        Biostrings::writeXStringSet(teSeq, filepath = gsub(fileNames$teNon80TmpLib, pattern = "tmp_non808080.fa", replacement = "complete_newfam.fa"), append = TRUE)
  #       Biostrings::writeXStringSet(teLibs$non808080Manual, file = gsub(fileNames$te80InputLib, pattern = "complete_models.fa", replacement = "2_MI_MCH/complete_newfam.fa"))

        cat("Appending TE sequence to the recovered library")
        Biostrings::writeXStringSet(teSeq, filepath = fileNames$te80RecovLib, append = TRUE)#append the non80 sequence to the recovered library
        # Biostrings::writeXStringSet(teLibs$non808080Manual[newAnnotation], file = gsub(fileNames$te80InputLib, pattern = "complete_models.fa", replacement = paste0(strain, "_808080_recovered.fa")), append = TRUE)

        #remove original seq from automatic annotated lib
        teLibs$teNon80TmpLib[input$te80Seq] <- NULL
        cat("Writting temporary analyzed library----\n")
        Biostrings::writeXStringSet(teLibs$teNon80TmpLib, filepath = fileNames$teNon80TmpLib)

        #update the TE seq to analyze
        updateSelectizeInput(session, inputId = "te80Seq", choices = names(teLibs$teNon80TmpLib))
      } else{
        cat("Trimming parameters not met----\n")
        output$te80TrimError <- renderText("Trimming parameters are not correct. Please check them")
      }
      if(length(teLibs$teNon80TmpLib) == 0){
        showModal(modalDialog(
          title = "You finished reviewing the non 80-80-80 TE models!!",
          paste0("There are sequences to analyze in selected library: "
                 , fileNames$teNon80TmpLib
                 , "\nPlease change your input library for 80-80-80 models"
          )
        ))
        updateTabsetPanel(session, inputId = "mchTabs", selected = "7")
      }
      updateSelectizeInput(session, inputId = "te70InputLib"
                           , choices = list.files(pattern = "^standby_sequences.fa", recursive = TRUE, full.names = TRUE)
                           , server = TRUE)
    }, once = FALSE) #eliminate the observe for the re-annotation button

    #### Standby TE model --------
    observeEvent(input$te80Stand, {
      req(length(teLibs$teNon80TmpLib) > 0)
      seqName <- input$te80Seq
      teSeq <- teLibs$teNon80TmpLib[seqName]
      #     strain <- input$te80InputLib %>% strsplit("/", fixed = TRUE) %>% lapply("[[", 4) %>% unlist()
      # input$te80Decision <- "CLASSI/Prueba" # debugging
      newAnnotation <- paste(
        seqName %>% strsplit("#", fixed = TRUE) %>% lapply("[[", 1) %>% unlist()
        , input$te80Decision
        , sep = "#"
      )
      names(teSeq) <- newAnnotation
      #remove original seq from automatic annotated lib
      cat("Sending the sequence to stand by model library----- \n")
      Biostrings::writeXStringSet(teSeq, file = fileNames$teStandbyLib, append = TRUE)
      teLibs$teNon80TmpLib[seqName] <- NULL
      cat("Writting temporary analyzed library----\n")
      Biostrings::writeXStringSet(teLibs$teNon80TmpLib, file = fileNames$teNon80TmpLib)

      updateSelectizeInput(session, inputId = "te80Seq", choices = names(teLibs$teNon80TmpLib))
      if(length(teLibs$teNon80TmpLib) == 0){
        showModal(modalDialog(
          title = "You finished reviewing the non 80-80-80 TE models!!",
          paste0("There are sequences to analyze in selected library: "
                 , fileNames$teNon80TmpLib
                 , "\nPlease change your input library for 80-80-80 models"
          )
        ))
        updateTabsetPanel(session, inputId = "mchTabs", selected = "7")
      }
      updateSelectizeInput(session, inputId = "te70InputLib"
                           , choices = list.files(pattern = "^standby_sequences.fa", recursive = TRUE, full.names = TRUE)
                           , server = TRUE)
    })

    #### Discard TE model --------
    observeEvent(input$te80Discard, {
      req(length(teLibs$teNon80TmpLib) > 0)
      #remove original seq from automatic annotated lib
      cat("Discarding the sequence and sending it to discarded fasta----- \n")
      Biostrings::writeXStringSet(teLibs$teNon80TmpLib[input$te80Seq], filepath = gsub(fileNames$teNon80TmpLib, pattern = "tmp_non808080.fa", replacement = "discarded_non808080.fa"), append = TRUE)
      teLibs$teNon80TmpLib[input$te80Seq] <- NULL
      cat("Writting temporary analyzed library----\n")
      Biostrings::writeXStringSet(teLibs$teNon80TmpLib, file = fileNames$teNon80TmpLib)

      updateSelectizeInput(session, inputId = "te80Seq", choices = names(teLibs$teNon80TmpLib))
      if(length(teLibs$teNon80TmpLib) == 0){
        showModal(modalDialog(
          title = "You finished reviewing the non 80-80-80 TE models!!",
          paste0("There are sequences to analyze in selected library: "
                 , fileNames$teNon80TmpLib
                 , "\nPlease change your input library for 80-80-80 models"
          )
        ))
      }
      # updateSelectizeInput(session, inputId = "te70InputLib"
      #                      , choices = list.files(pattern = "^standby_sequences.fa", recursive = TRUE, full.names = TRUE)
      #                      , server = TRUE)
      updateTabsetPanel(session, inputId = "mchTabs", selected = "7")
      })
  
  output$te70InputPlot <- renderPlotly({
    tePlotLib(teLibList = list(teLibs$teStandbyLib, teLibs$te70Lib, teLibs$te70FiltLib))
  })
  output$te70InputLib <- renderText(fileNames$teStandbyLib)
  
  observeEvent(input$te70ModuleRun,{
    if(length(teLibs$teNon80TmpLib) != 0){
      cat("Please finish with the manual inspection of non 808080 models", fileNames$teNon80TmpLib)
      showModal(modalDialog(
        title = "Please finish with the Manual Inspection of 808080 models before proceeding!"
        , paste0("You must have reviewed all the sequences before going to 70-70-70 model. There are still sequences in the non808080 library to review manually: \n"
                 , gsub(pattern = "standby_complete.fa", replacement = "tmp_non808080.fa", fileNames$te70InputLib)
        )
      ))
    } else{#if the curated library has been completely reviewed
      cat("The temporary curated library is empty, executing the 707070 module\n")
      # #       req(file.exists(fileNames$teStandbyLib))
      # 4 Blastn 70-70-70 module for complete models ----
      # # #       strain <- fileNames$te70InputLib %>% strsplit("/", fixed = TRUE) %>% lapply("[[", 4) %>% unlist() #"N2" #input$teStrain
      # # #       species <- fileNames$te70InputLib %>% strsplit("/", fixed = TRUE) %>% lapply("[[", 3) %>% unlist() #"C. elegans"#input$teSpecies
      fileNames$teStandBlast <- gsub(fileNames$teStandbyLib, pattern = "standby_sequences.fa", replacement = paste0("standby", "_vs_allDBs.blast"))
      fileNames$teStandBlastLog <- gsub(fileNames$teStandbyLib, pattern = "standby_sequences.fa", replacement = "blast707070.log")
      
      if(file.exists(fileNames$teStandBlastLog)){
        #if blast 707070 outputs are present, raise a warning
        warning(paste0("Looks like blast agains TE database has been already executed for library ", fileNames$teStandbyLib, ", please check\n"))
        # #         logFileBlast70 <- readLines(fileNames$teStandBlastLog)
      }else{#if blast 707070 outputs are not present, do the blast
        cat("Adding high confidence models to the original database----\n")
        system(
          paste0("cat mchelper-ats/db/allDatabases.clustered_rename.fa "
                 , fileNames$te80RecovLib #gsub(pattern = "2_MI_MCH/standby_complete.fa", replacement = "N2-test_808080_recovered.fa", x = fileNames$te70InputLib)#MCHelper_automatic/C.elegans/N2-test/1_MI_MCH/N2-test_808080_recovered.fa >,
                 , " > "
                 , fileNames$te70InputDb #"MCHelper_automatic/C.elegans/N2-test/1_MI_MCH/allDatabases.clustered_merged.fa"
          )
        )
        cat("Creating database index----\n")
        system2(
          "conda"
          , args = c(
            "run -n MCHelper makeblastdb -in"
            , fileNames$te70InputDb
            , "-dbtype nucl"
          )
        )
        cat("Performing the blast --------------- \n")
        system2(
          "conda",
          args = c(
            "run"
            , "-n", "MCHelper"
            , "--no-capture-output"
            , "blastn"
            ,"-query", fileNames$teStandbyLib
            ,"-db", fileNames$te70InputDb
            , "-num_threads", input$teAutoModuleCores
            , "-outfmt 6 -qcov_hsp_perc 70 -perc_identity 70 -max_hsps 1"
            , "-out", fileNames$teStandBlast
          )
          , stdout = fileNames$teStandBlastLog
          , stderr = fileNames$teStandBlastLog
          , wait = TRUE
        )
                # conda run -n MCHelper --no-capture-output blastn -query  MCHelper_automatic/C.elegans/N2-test/1_MI_MCH/2_MI_MCH/standby_complete.fa  -db ${species}/curation/complete_models/allDatabases.clustered_rename.fa -num_threads 32 -outfmt 6 -out ${species}/curation/complete_models/${species}_standBy_vs_allDBs.blast -qcov_hsp_perc 70 -perc_identity 70 -max_hsps 1
        #
        cat("Blast successfully run --------------- \n")
        # #         logFileBlast70 <- readLines(fileNames$teStandBlastLog)
      }#end of 707070 blast
      # # #       
      #write stand by models with 707070 hits to a fasta file for cross match test
      teTables$teStandBlast <- read.table(file = fileNames$teStandBlast, sep = "\t", comment.char = ""
                                          , col.names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
      )
      Biostrings::writeXStringSet(teLibs$teStandbyLib[teTables$teStandBlast$qseqid], filepath = fileNames$te70Lib)
      # # 
      # #       #run 707070 filtering
      if(file.exists(fileNames$te70FilterLog)){
        warning("Looks like 707070 filtering has been run")
        # #         te70FilterLog <- readLines(fileNames$te70FilterLog)
      } else{
        cat("Running 707070 filtering ------ \n")
        system2(
          "conda"
          , args = c(
            "run -n MCHelper"
            , "python3 mchelper-ats/tools/Curation_toolkit/curation_toolkit.py 707070_filtering "
            , fileNames$te70Lib #input lib path
            , fileNames$te70FiltLib#output lib path
          )
          , stdout = fileNames$te70FilterLog
          , stderr = fileNames$te70FilterLog
          , wait = TRUE
        )
        te70FilterLog <- readLines(fileNames$te70FilterLog)
      }
    }
    if(file.exists(fileNames$te70FilterLog)) {
      showModal(modalDialog(
        title = "You already run 70-70-70, moving to Incomplete models!"
      ))
      updateTabsetPanel(session, inputId = "mchTabs", selected = "8")
    }
    # # 
    # # 
    # #     paste(c("Blast results --------"
    # #             , fileNames$teStandBlastLog
    # #             , logFileBlast70
    # #             ,"\n\n Filter 707070 results -----"
    # #             , fileNames$te70FilterLog
    # #             , te70FilterLog
    # #     ), collapse = "\n"
    # #     )
  })#end of run te70Module
  
  output$teIncompleteInputPlot <- renderPlotly({
    tePlotLib(list(teLibs$teAutoCuratedLib, teLibs$teIncompleteModelsLib, teLibs$teIncompleteRecovLib))
  })
  output$teIncompleteNon80InputLib <- renderText(fileNames$teIncompleteNon80)
  observeEvent(input$teIncompleteModuleRun,{
      ### 5 Blastn Incomplete module for complete models ----
      strain <- fileNames$teIncompleteModelsLib %>% strsplit("/", fixed = TRUE) %>% lapply("[[", 4) %>% unlist() #"N2" #input$teStrain
      fileNames$teIncompleteBlast <- gsub(x = fileNames$teIncompleteModelsLib, pattern = "incomplete_models.fa", replacement = "incomplete_models_vs_allDBs.blast")
      fileNames$teIncompleteBlastLog <- gsub(fileNames$teIncompleteModelsLib, pattern = "incomplete_models.fa", replacement = "blast808080_incomplete.log")
      fileNames$teIncompleteNon80Log <- gsub(fileNames$teIncompleteNon80, pattern = "incomplete_non808080.fa", replacement = "mchelper_manual_incomplete.log")
      fileNames$teIncompleteAssignLog <- gsub(fileNames$teIncompleteModelsLib, pattern = "incomplete_models.fa", replacement = "assign_families_incomplete.log")
      if(length(teLibs$teAutoCuratedTmp) != 0){
        showModal(modalDialog(
          title = "Please finish with the Manual Inspection of Automatic module until before proceeding!"
          , paste0("You must have obtained a high confidence database (80-80-80 models) before going to Incomplete module. There are no sequences in the Incomplete library:\n"
                   , fileNames$teIncompleteModelsLib
          )
        ))
        if(!file.exists(fileNames$te70InputDb)){#check that te database 2 has been created and indexed
          showModal(modalDialog(
            title = "Please check that extended TE database has been created and indexed (808080 module)!"
            , paste0("You must have created an extended TE database with high confidence models (80-80-80 models) before going to Incomplete module. I could not find the database: \n"
                     , fileNames$te70InputDb
            )
          ))
        }
      } else{#if the curated library has been completely reviewed
        cat("The temporary curated library is empty, executing the Incomplete module\n")
  # #       req(file.exists(fileNames$teIncompleteModelsLib))
        if(file.exists(fileNames$teIncompleteBlastLog)){
          #if blast Incomplete outputs are present, raise a warning
          warning(paste0("Looks like blast agains TE database has been already executed for library ", fileNames$teIncompleteModelsLib, ", please check\n"))
  # #         logFileBlastIncomplete <- readLines(fileNames$teIncompleteBlastLog)
        }else{#if blast Incomplete outputs are not present, do the blast
          cat("Performing the blast --------------- \n")
          system2(
            "conda",
            args = c(
              "run"
              , "-n", "MCHelper"
              , "--no-capture-output"
              , "blastn"
              ,"-query", fileNames$teIncompleteModelsLib
              ,"-db", fileNames$te70InputDb #is the same input as db for 7070707 models
              , "-num_threads", input$teAutoModuleCores
              , "-outfmt 6 -qcov_hsp_perc 80 -perc_identity 80 -max_hsps 1"
              , "-out", fileNames$teIncompleteBlast
            )
            , stdout = fileNames$teIncompleteBlastLog
            , stderr = fileNames$teIncompleteBlastLog
            , wait = TRUE
          )
          # conda run -n MCHelper --no-capture-output blastn -query  MCHelper_automatic/C.elegans/N2-test/1_MI_MCH/2_MI_MCH/standby_complete.fa  -db ${species}/curation/complete_models/allDatabases.clustered_rename.fa -num_threads 32 -outfmt 6 -out ${species}/curation/complete_models/${species}_standBy_vs_allDBs.blast -qcov_hsp_perc Incomplete -perc_identity Incomplete -max_hsps 1
  # # 
          cat("Blast successfully run --------------- \n")
  # #         logFileBlastIncomplete <- readLines(fileNames$teIncompleteBlastLog)
        }#end of Incomplete blast
  # # 
      }
  # # 
  # #     
      if(file.exists(fileNames$teIncompleteAssignLog)){
        #if assign_families outputs are present, raise a warning
        warning(paste0("Looks like assign_families has been already executed for library", fileNames$te80InputLib, ", please check\n"))
  # #       logFileAssign <- readLines(fileNames$teIncompleteAssignLog)
      }else{#if assign_families has not been executed, execute it
        cat("Running assign_families --------------- \n")
        system2(
          "conda"
          , args = c(
            "run"
            , "-n", "MCHelper"
            , "--no-capture-output"
            , "python3", "mchelper-ats/tools/Curation_toolkit/curation_toolkit.py"
            , "assign_families" #curation toolkit module
            , strain #strainID
            , fileNames$teIncompleteModelsLib #complete TE library
            , fileNames$teIncompleteBlast #blast result
            , gsub(fileNames$teIncompleteModelsLib, pattern = "incomplete_models.fa", replacement = "incomplete") # output files prefix
          )
          , stdout = fileNames$teIncompleteAssignLog
          , stderr = fileNames$teIncompleteAssignLog
          , wait = TRUE
        )
        # conda run -n MCHelper --no-capture-output python3 mchelper-ats/tools/Curation_toolkit/curation_toolkit.py assign_families N2-test MCHelper_automatic/C.elegans/N2-test/1_MI_MCH/complete_models.fa MCHelper_automatic/C.elegans/N2-test/1_MI_MCH/N2-test_vs_allDBs.blast MCHelper_automatic/C.elegans/N2-test/1_MI_MCH/N2-test_complete
        cat("Assign families successfully run --------------- \n")
  # #       logFileAssign <- readLines(fileNames$teIncompleteAssignLog)
      }
  # # 
  # #     ### Run MCHelper on non808080.fa
  # # 
      if(!file.exists(fileNames$teIncompleteNon80Log)){# if manual module has not run (does the log file exist?), run it
        cat("Running MCHelper in T mode for non808080 sequence, please be patient...\nInput library:", fileNames$teIncompleteNon80)
        system2(
          "conda"
          , args = c(
            "run",
            "-n", "MCHelper",
            "--no-capture-output",
            "python3",
            "mchelper-ats/MCHelper.py"
            # "-a", "F"
            ,"-r", "T" #use T module as it only runs TE-Aid that is what we want in this part
            , "--input_type", "fasta",
            "-l", fileNames$teIncompleteNon80
            , "-g", fileNames$teGenome
            , "-o", gsub(fileNames$teIncompleteNon80, pattern = "incomplete_non808080.fa", replacement = "3_MI_MCH"),
            "-t", input$teAutoModuleCores
            , "-v", "Y"
          )
          , stdout = fileNames$teIncompleteNon80Log
          , stderr = fileNames$teIncompleteNon80Log
          , wait = FALSE #avoid frozen app while running MCHelper (send process to background)
        )
  # #       logFileIncomplete <- paste0("Running MCHelper in T mode for library "
  # #                           , fileNames$teIncompleteNon80
  # #                           , ", please be patient...\n")
      } else{#looks like the manual module has been run
  # #       logFileIncomplete <- readLines(fileNames$teIncompleteNon80Log)
  # #         #if manual curation stamp is present
        if(file.exists(fileNames$teIncompleteNon80Tmp)){
          cat("A manual curation stamp is present, loading it ---- \n")
          teLibs$teIncompleteNon80Tmp <- teReadLib(fileNames$teIncompleteNon80Tmp)
  # # #         teLibs$non808080IncompleteManual <- Biostrings::readDNAStringSet(fileNames$teIncompleteNewfam)
          updateSelectizeInput(session, inputId = "teIncompleteSeq", choices = names(teLibs$teIncompleteNon80Tmp))
  # # #         # input$teIncompleteSeq <- names(teLibs$non808080IncompleteTmp)[1] #debugging
          if(length(teLibs$teIncompleteNon80Tmp) == 0){
            showModal(modalDialog(
              title = "Looks like you already reviewed the non808080 Incomplete TE models!!",
              paste0("There are no sequences to analyze in selected library: "
                     , fileNames$teIncompleteNon80Tmp
                     , "\nPlease change your input library for 808080 models"
              )
            ))
          }
        } else{#if no manual review has been saved, initialize temporary DNAstringsets
          cat("A manual curation stamp is NOT present, creating files ---- \n")
          teLibs$teIncompleteNon80Tmp <- Biostrings::readDNAStringSet(fileNames$teIncompleteNon80)
  # #         teLibs$teIncomopleteRecovLib <- Biostrings::DNAStringSet()
  # #         # Biostrings::writeXStringSet(teLibs$non808080IncompleteManual, file = fileNames$teIncompleteNewfam, append = FALSE) # append false to ensure a new file is created
          updateSelectizeInput(session, inputId = "teIncompleteSeq", choices = names(teLibs$teIncompleteNon80Tmp))
        }
        teTables$teIncompleteCuratedClassif <- read.table(gsub(fileNames$teIncompleteNewfam, pattern = "incomplete_newfam.fa", replacement = "denovoLibTEs_PC.classif"), sep = "\t", header = TRUE)
  # #       # teLibs$non808080 <- Biostrings::readDNAStringSet(fileNames$te80InputLib)
        # fileNames$te_aidIncompleteGenomeBnFiles <- list.files(gsub(fileNames$teIncompleteNewfam, pattern = "incomplete_newfam.fa", replacement = ""), pattern = "genome.blastn", recursive = TRUE, full.names = TRUE)
        # fileNames$te_aidIncompleteSelfBnFiles <- list.files(gsub(fileNames$teIncompleteNewfam, pattern = "incomplete_newfam.fa", replacement = ""), pattern = "self.blastn", recursive = TRUE, full.names = TRUE)
  # #       
        if(length(teLibs$teIncompleteNon80Tmp) == 0){
          showModal(modalDialog(
            title = "There are no sequences in the non808080 temporal library, looks like you already reviewed this library!!",
            paste0("You can manually check the non80-80-80 file:"
                   , fileNames$teIncompleteNon80
                   , "or even check the incomplete models library:", fileNames$teIncompleteModelsLib
                   , "\nPlease change your complete sequences input library to manually inspect the non 80-80-80 models"
            )
          ))
          updateTabsetPanel(session, inputId = "mchTabs", selected = "9")
        }
      }#end of run 808080 models for incomplete
      
      output$teIncompleteStats <- renderTable({
        # teList <- teLibs #debugging
        teList <- reactiveValuesToList(teLibs)
        teCurationStats(teLibList = teList[grep(names(teList), pattern = "CuratedLib|Incomplete", value = TRUE)])
      }, rownames = TRUE)
  # # #   #   
  # # #   #   
  # #     paste(c("Blast results --------"
  # #             , fileNames$teIncompleteBlastLog
  # #             , logFileBlastIncomplete
  # #             ,"\n\n Assign families results -----"
  # #             , fileNames$teIncompleteAssignLog
  # #             , logFileAssign
  # #             , "\n\n MCHelper -T results ----"
  # #             , logFileIncomplete
  # #     ), collapse = "\n"
  # #     )
    })#end of run teIncompleteModule
  
    ######-
    #### Plotting non808080 sequences----
    ######-

    observeEvent(input$teIncompleteSeq, {
      # input$teIncompleteSeq <- names(teLibs$teIncompleteNon80Tmp)[1] #for debugging
      req(input$teIncompleteSeq)
      seqName <- input$teIncompleteSeq %>% strsplit("#") %>% unlist() %>% head(1)
      teIncompleteAnnotation <- input$teIncompleteSeq %>% strsplit("#") %>% lapply("[[", 2) %>% unlist()
      updateSelectInput(session, inputId = "teIncompleteDecision"
                        # , choices = c("Hola", "Hola3", "Hols3", teAnnotation)#this should be changed to the full table of annotations (42 entries?)
                        , selected = teIncompleteAnnotation)
      # input$teIncompleteDecision <- teAnnotation #for debugging
      # fileNames$te_aidIncompleteGenomeBnFiles <- list.files(gsub(fileNames$teIncompleteNewfam, pattern = "incomplete_newfam.fa", replacement = ""), pattern = "genome.blastn", recursive = TRUE, full.names = TRUE)
      # fileNames$te_aidIncompleteSelfBnFiles <- list.files(gsub(fileNames$teIncompleteNewfam, pattern = "incomplete_newfam.fa", replacement = ""), pattern = "self.blastn", recursive = TRUE, full.names = TRUE)
      
      teBlastGenome <- list.files(fileNames$teIncompleteOutDir, pattern = glob2rx(paste0(seqName, "*.fa*genome.blast*")), recursive = TRUE, full.names = TRUE) %>% tail(1)
      teBlastSelf <- list.files(fileNames$teIncompleteOutDir, pattern = glob2rx(paste0(seqName, "*.fa*self.blast*")), recursive = TRUE, full.names = TRUE) %>% tail(1)
      
      output$teIncompleteCuratedClass <- renderTable({
        teTables$teIncompleteCuratedClassif %>% subset(Seq_name == seqName) %>% select(Seq_name, length, class, order, sFamily, struct, coding, other, strand)
      })
      output$teIncompleteConsDivergence <- renderPlot({
        teConsDivergencePlot(seqName, blastFiles = teBlastGenome, leftTrim = NULL, rightTrim = NULL)
      })
      output$teIncompleteConsCoverage <- renderPlotly({
        teConsCoveragePlot(seqName, blastFiles = teBlastGenome, leftTrim = NULL, rightTrim = NULL, source = "C")
      })
      output$teIncompleteConsSelf <- renderPlot({
        teConsSelfPlot(seqName, blastFiles = teBlastSelf, leftTrim = NULL, rightTrim = NULL)
      })
      output$teIncompleteConsStructure <- renderPlot({
        teConsStructurePlot(seqName, blastFiles = teBlastSelf, leftTrim = NULL, rightTrim = NULL)
      })

    }) #end of observeEvent for Incomplete input sequence
    #### Reannotating and Trim TE ------
    observeEvent(input$teIncompleteAnnotate, {#if no change annotation, new annotation == old annotation
      req(length(teLibs$teIncompleteNon80Tmp) > 0)
      seqName <- input$teIncompleteSeq
      teSeq <- teLibs$teIncompleteNon80Tmp[seqName]
      # input$teIncompleteDecision <- "CLASSI/Prueba" # debugging
      newAnnotation <- paste(
        seqName %>% strsplit("#", fixed = TRUE) %>% lapply("[[", 1) %>% unlist()
        , input$teIncompleteDecision
        , sep = "#"
      )
      #     #re-naming in manual.v2
      names(teSeq) <- newAnnotation
      cat("Writting TE into the recovered library -----\n")
      Biostrings::writeXStringSet(teSeq, filepath = fileNames$teIncomopleteRecovLib, append = TRUE)

      cat("Writting TE into the newfamilies library -----\n")
      Biostrings::writeXStringSet(teSeq, filepath = fileNames$teIncompleteNewfam, append = TRUE)

      #remove original seq from automatic annotated lib
      teLibs$teIncompleteNon80Tmp[input$teIncompleteSeq] <- NULL
      cat("Writting temporary analyzed library----\n")
      Biostrings::writeXStringSet(teLibs$teIncompleteNon80Tmp, filepath = fileNames$teIncompleteNon80Tmp, append = FALSE) #append false to ensure we save this file in a destructive way

      #update the TE seq to analyze
      updateSelectizeInput(session, inputId = "teIncompleteSeq", choices = names(teLibs$teIncompleteNon80Tmp))

      # # trimming
      # if((input$leftTrimIncomplete>=1) & (input$rightTrimIncomplete > input$leftTrimIncomplete) & (input$rightTrimIncomplete<=Biostrings::width(teLibs$teNon80TmpLib[seqName]))){#if trimming parameters are ok
      #   cat("Trimming sequence-----\n")
      #   teSeq <- subseq(teSeq[newAnnotation], start = input$leftTrimIncomplete, end = input$rightTrimIncomplete)
      #   clickData$click1 <- NULL
      #   clickData$click2 <- NULL
      #
      #
      #
      # } else{
      #   cat("Trimming parameters not met----\n")
      #   output$teIncompleteTrimError <- renderText("Trimming parameters are not correct. Please check them")
      # }

      if(length(teLibs$teIncompleteNon80Tmp) == 0){
        showModal(modalDialog(
          title = "You finished reviewing the non 80-80-80 TE models!!",
          paste0("There are sequences to analyze in selected library: "
                 , fileNames$teIncompleteNon80Tmp
                 , "\nPlease change your input library for 80-80-80 models"
          )
        ))
        updateTabsetPanel(session, inputId = "mchTabs", selected = "9")
      }
    })

    #### Discard TE model --------
    observeEvent(input$teIncompleteDiscard, {
      req(length(teLibs$teIncompleteNon80Tmp) > 0)
      #remove original seq from automatic annotated lib
      cat("Discarding the sequence and sending it to discarded fasta----- \n")
      Biostrings::writeXStringSet(teLibs$teIncompleteNon80Tmp[input$teIncompleteSeq], filepath = gsub(fileNames$teIncompleteNon80Tmp, pattern = "tmp_incomplete_non808080.fa", replacement = "discarded_non808080.fa"), append = TRUE)
      teLibs$teIncompleteNon80Tmp[input$teIncompleteSeq] <- NULL
      cat("Writting temporary analyzed library----\n")
      Biostrings::writeXStringSet(teLibs$teIncompleteNon80Tmp, file = fileNames$teIncompleteNon80Tmp)

      updateSelectizeInput(session, inputId = "teIncompleteSeq", choices = names(teLibs$teIncompleteNon80Tmp))
      if(length(teLibs$teIncompleteNon80Tmp) == 0){
        showModal(modalDialog(
          title = "You finished reviewing the non 80-80-80 TE models!!",
          paste0("There are sequences to analyze in selected library: "
                 , fileNames$teIncompleteNon80Tmp
                 , "\nPlease change your input library for 80-80-80 models"
          )
        ))
        updateTabsetPanel(session, inputId = "mchTabs", selected = "9")
      }
    })
  
  #############-
  ## 6 Final MCHelper module----
  #############-
  ### Plots for automatic comparison----
  output$teFinal <- renderPlotly({
    tePlotLib(teLibList = list(teLibs$teRepeatLib, teLibs$te80RecovLib, teLibs$te70FiltLib, teLibs$teIncompleteRecovLib, teLibs$teFinalNR), libType = "Final TE libraries to join")
  })
  # #   observeEvent(input$teFinalShowLog, {
  # #     output$debug <- renderUI({
  # #       if(input$teFinalShowLog){
  # #         fluidRow(
  # #           column(6, plotlyOutput("teFinal", height = "120%"))
  # #           , column(6, verbatimTextOutput("teFinalLog"))
  # #         )
  # #       } else {
  # #         plotlyOutput("teFinal", height = "130%")
  # #       }
  # #     })
  # #   })
  # #   
    observeEvent(input$teFinalModuleRun, {
  # #     # req(length(teLibs$teRepeatLib) > 0)
  # #     req(fileNames$teRepeatLib)
  # #     strain <- fileNames$teRepeatLib %>% strsplit("/", fixed = TRUE) %>% lapply("[[", 4) %>% unlist()
      if (dir.exists(fileNames$teFinalMchDir)){
        cat(fileNames$teFinalMchDir, " already existing ---- \n")
      } else {dir.create(fileNames$teFinalMchDir); cat(fileNames$teFinalMchDir, " directory created ---- \n")}
      fileNames$teFinalJoinLog <- paste0(x = fileNames$teFinalMchDir, "/join_libraries.log")
      fileNames$teFinalJoin <- paste0(x = fileNames$teFinalMchDir, "/final_curated_R_tmp.fa") #gsub(x = fileNames$teRepeatLib, pattern = "clean-families.fa", replacement = "")
      fileNames$teFinalJoinRename <- paste0(fileNames$teFinalMchDir, "/final_curated_R.fa")
      fileNames$teFinalRenameLog <- paste0(fileNames$teFinalMchDir, "/rename_libraries.log")
      fileNames$teFinalNR <- paste0(fileNames$teFinalMchDir, "/final_curated_NR.fa")
      fileNames$teFinalCdhitLog <- paste0(fileNames$teFinalMchDir, "/cd-hit_final.log")
      fileNames$teFinalMchLog <- paste0(fileNames$teFinalMchDir, "/mchelper_manual_final.log")
      # fileNames$teFinalMchDir <- paste0(fileNames$teRepeatLib, "/MCH_final")
  # #     #### Join libraries ----
      if(file.exists(fileNames$teFinalJoinLog)){
        warning("Libraries already joined ----\n")
  # #       joinLog <- readLines(fileNames$teFinalJoinLog)
      }else{
        system2("cat"
                , args = c(fileNames$te80RecovLib, fileNames$te70FiltLib, fileNames$teIncomopleteRecovLib)
                , stdout = fileNames$teFinalJoin
                , stderr = fileNames$teFinalJoinLog
        )
        # paste("cat", fileNames$te80RecovLib, fileNames$te70FiltLib, fileNames$teIncomopleteRecovLib, "> final_curated_R_tmp.fa")
        cat("Join libraries executed ---- \n")
  # #       joinLog <- "Join libraries being executed ---- \n"
      }
  # #     
  # #     ### Rename duplicates----
      if(file.exists(fileNames$teFinalRenameLog)){
        warning("Rename duplicates already run ----\n")
  # #       renameLog <- readLines(fileNames$teFinalRenameLog)
      } else{
        system2(
          "conda"
          , args = c(
            "run -n MCHelper --no-capture-output python3 mchelper-ats/tools/Curation_toolkit/curation_toolkit.py"
            , "rename_duplicates"
            , fileNames$teFinalJoin
            , fileNames$teFinalJoinRename
          )
          , stdout = fileNames$teFinalRenameLog
          , stderr = fileNames$teFinalRenameLog
        )
        cat("Rename duplicates executed ---- \n")
  # #       renameLog <- "Rename duplicates being executed ----\n"
      }
  # #     
  # #     ### CD-hit reduction----
      if(file.exists(fileNames$teFinalCdhitLog)){
  # #       warning("CD hit already run ----\n")
  # #       cdHitLog <- readLines(fileNames$teFinalCdhitLog)
      } else {
        system2(
          "conda"
          , args = c(
            "run -n MCHelper --no-capture-output cd-hit-est"
            ,"-i", fileNames$teFinalJoinRename
            , "-o", fileNames$teFinalNR
            , "-c 0.95 -aS 0.98 -G 0 -g 1 -b 500 -M 0 -d 0"
            , "-T", input$teFinalCores
          )
          , stdout = fileNames$teFinalCdhitLog
          , stderr = fileNames$teFinalCdhitLog
        )
  # #       cdHitLog <- "CD-hit being executed ----\n"
        cat("CD-hit executed ---- \n")
  # #       # cd-hit-est -i ${species}/curation/final_curated_TE_lib_R.fa -o ${species}/curation/final_curated_TE_lib_NR.fa -c 0.95 -aS 0.98 -G 0 -g 1 -b 500 -T 32 -M 0 -d 0
      }
  # #     
  # #     ### MCHelper visualization----
      if(file.exists(fileNames$teFinalMchLog)){
        warning("MCHelper manual mode already executed for final TE library ----\n")
        mchLog <- readLines(fileNames$teFinalMchLog)
        if(mchLog[length(mchLog)] == "MCHelper with TE-Aid successfully run"){
          teLibs$teFinalNR <- teReadLib(fileNames$teFinalNR, libIdentifier = "Final curated TE library")
          updateSelectizeInput(session, inputId = "teFinalSeq", choices = names(teLibs$teFinalNR))
          teTables$teFinalCuratedClassif <- read.table(paste0(fileNames$teFinalMchDir, "/denovoLibTEs_PC.classif"), header = TRUE, sep = "\t", comment.char = "")
        } else{
          showModal(modalDialog(
            title = "MCHelper still running for final curated TE library"
            ,"Wait until it finishes"
          ))
        }

      } else {
        cat("Running MCHelper manual inspection for ", fileNames$teFinalNR, "----\n")
        system2(
          "conda",
          args = c(
            "run -n MCHelper --no-capture-output python3 mchelper-ats/MCHelper.py"
            # "-a", "F"
            ,"-r", "T" #use T module as it only runs TE-Aid that is what we want in this part
            , "--input_type", "fasta",
            "-l", fileNames$teFinalNR,
            "-g", fileNames$teGenome
            , "-o", fileNames$teFinalMchDir,
            "-t", input$teFinalCores
            , "-v", "Y"
          )
          , stdout = fileNames$teFinalMchLog
          , stderr = fileNames$teFinalMchLog
          , wait = FALSE #avoid frozen app while running MCHelper (send process to background)
        )
        cat("MCHelper manual inspection module to background for ", fileNames$teAutoCuratedLib, "----\n")
  # #       mchLog <- "MCHelper manual inspection module sent to background"
      }
  # #     
      output$teFinalLib <- renderText(fileNames$teFinalNR)
      
      
      output$teFinalStats <- renderTable({
        # teList <- teLibs #debugging
        teList <- reactiveValuesToList(teLibs)
        teCurationStats(teLibList = teList[grep(names(teList), pattern = "Repeat|Recov|Filt|Final", value = TRUE)])
      }, rownames = TRUE)
      last_idx <- tryCatch({read.table(paste0(fileNames$teFinalMchDir, "/last_index.txt"))$V1}, error = function(e) 1)
      teIndex(last_idx)
      updateSelectizeInput(session, inputId = "teFinalSeq", choices = names(teLibs$teFinalNR), selected = names(teLibs$teFinalNR)[last_idx])
  # #     output$teFinalLog <- renderText(
  # #       paste(c("Join libraries already run ---\n"
  # #               , joinLog
  # #               , "Rename duplicates already run ----\n"
  # #               , renameLog
  # #               , "CD-hit reduction already run ----\n"
  # #               , cdHitLog
  # #               , "MCHelper TE aid already run ----\n"
  # #               , mchLog
  # #       )
  # #       , collapse = "\n"
  # #       )
  # #     )
    })
  # #   
  # #   
  # #   
    observeEvent(input$teFinalSeq, {
      # input$teFinalSeq <- names(teLibs$teFinalNon80Tmp)[1] #for debugging
      # input$teFinalSeq <- "Mariner9_CB#CLASSII/MITE"
      idx <- match(input$teFinalSeq, names(teLibs$teFinalNR))
      if (!is.na(idx)) {
        teIndex(idx)
      }
      req(input$teFinalSeq)
      write.table(teIndex(), file = paste0(fileNames$teFinalMchDir, "/last_index.txt"), row.names = FALSE, col.names = FALSE)
      seqName <- names(teLibs$teFinalNR)[teIndex()] %>% strsplit("#") %>% unlist() %>% head(1)
      teFinalAnnotation <- input$teFinalSeq %>% strsplit("#") %>% lapply("[[", 2) %>% unlist()
      updateSelectInput(session, inputId = "teFinalAnnotation"
                        # , choices = c("Hola", "Hola3", "Hols3", teAnnotation)#this should be changed to the full table of annotations (42 entries?)
                        , selected = teFinalAnnotation)
      # input$teFinalDecision <- teAnnotation #for debugging
      teBlastGenome <- list.files(fileNames$teFinalMchDir, pattern = glob2rx(paste0(seqName, ".fa*genome.blast*")), recursive = TRUE, full.names = TRUE) %>% tail(1)
      teBlastSelf <- list.files(fileNames$teFinalMchDir, pattern = glob2rx(paste0(seqName, ".fa*self.blast*")), recursive = TRUE, full.names = TRUE) %>% tail(1)
      # output$teFinalCuratedClass <- renderTable({
      #   teTables$teFinalCuratedClassif %>% subset(Seq_name == seqName) %>% select(Seq_name, length, class, order, sFamily, struct, coding, other, strand)
      # })
      teTable <- teTables$teFinalCuratedClassif
      output$teFinalClass <- renderTable(
        teClassOutput(teTable, seqName)$teInfo
      )
      output$teFinalCoding <- renderTable(
        teClassOutput(teTable, seqName)$teCoding, colnames = FALSE
      )
      output$teFinalStructure <- renderTable(
        teClassOutput(teTable, seqName)$teStruc
      )
      
      output$teFinalConsDivergence <- renderPlot({
        teConsDivergencePlot(seqName, blastFiles = teBlastGenome, leftTrim = NULL, rightTrim = NULL)
      })
      output$teFinalConsCoverage <- renderPlotly({
        teConsCoveragePlot(seqName, blastFiles = teBlastGenome, leftTrim = NULL, rightTrim = NULL, source = "D")
      })
      output$teFinalConsSelf <- renderPlot({
        teConsSelfPlot(seqName, blastFiles = teBlastSelf, leftTrim = NULL, rightTrim = NULL)
      })
      output$teFinalConsStructure <- renderPlot({
        teConsStructurePlot(seqName, blastFiles = teBlastSelf, leftTrim = NULL, rightTrim = NULL)
      })
    }) #end of observeEvent for Final input sequence
  
  teIndex <- reactiveVal(1)
  # "Next" button
  observeEvent(input$teFinalNext, {
    new_idx <- teIndex() + 1
    if (new_idx >= 1 & new_idx <= length(teLibs$teFinalNR)) {
      teIndex(new_idx)
      updateSelectizeInput(session, "teFinalSeq", selected = names(teLibs$teFinalNR)[new_idx])
    }
    write.table(teIndex(), file = paste0(fileNames$teFinalMchDir, "/last_index.txt"), row.names = FALSE, col.names = FALSE)
  })
  # "Previous" button
  observeEvent(input$teFinalPrev, {
    new_idx <- teIndex() - 1
    if (new_idx >= 1) {
      teIndex(new_idx)
      updateSelectizeInput(session, "teFinalSeq", selected = names(teLibs$teFinalNR)[new_idx])
    }
    write.table(teIndex(), file = paste0(fileNames$teFinalMchDir, "/last_index.txt"), row.names = FALSE, col.names = FALSE)
  })
  observeEvent(input$teFinalDiscard, {
    Biostrings::writeXStringSet(teLibs$teFinalNR[input$teFinalSeq], filepath = paste0(fileNames$teFinalMchDir, "/discarded_sequences.fa"), append = TRUE)
    teLibs$teFinalNR[input$teFinalSeq] <- NULL
    Biostrings::writeXStringSet(teLibs$teFinalNR, filepath = fileNames$teFinalNR)
    updateSelectizeInput(session, inputId = "teFinalSeq", selected = names(teLibs$teFinalNR)[teIndex()])
  })
  observeEvent(input$teFinalSave, {
    seqName <- names(teLibs$teFinalNR)[teIndex()] %>% strsplit("#") %>% unlist() %>% head(1)
    newAnnotation <- paste(
      seqName
      , input$teFinalAnnotation
      , sep = "#"
    )
    names(teLibs$teFinalNR)[teIndex()] <- newAnnotation
    Biostrings::writeXStringSet(teLibs$teFinalNR, filepath = fileNames$teFinalNR)
    updateSelectizeInput(session, inputId = "teFinalSeq", selected = newAnnotation, choices = names(teLibs$teFinalNR))
  })
  
  observeEvent(input$teAddLibToDb, {
    req(teLibs$teFinalNR != 0)
    # if log add db exists, show a modal to prevent adding twice the library
    # backup the previous db
    # cat the final library to the previous library
    # makeblastdb of the merged library
    # generate a log file to prevent adding the same library twice
  })
  # #############-
  # ## End of MCHelper pipeline----
  # #############-
  
  # #############-
  # ## HELIANO ----
  # #############-
  helData <- reactive({
    # rv$teDataRefresh  # React to the trigger value
    fileNames$helFiles <- list.files("./HELIANO", pattern = "NR-clean.fa", recursive = TRUE, full.names = TRUE) #debugging
    data <- data.frame(file_name = fileNames$helFiles) %>% 
      mutate(
        Species = file_name %>% strsplit("/", fixed = TRUE) %>% lapply("[[", 3) %>% unlist()
        , Strain = file_name %>% strsplit("/", fixed = TRUE) %>% lapply("[[", 4) %>% unlist()
        , Raw_lib = file_name %>% strsplit("/", fixed = TRUE) %>% lapply("[[", 5) %>% unlist()
        , Genome_file = grep(Strain, fileNames$teGenomes, value = TRUE)
      )
    # data$Status <- lapply(1:nrow(data), function(i) getMchStatus(data$Species[i], data$Strain[i], data$Raw_lib[i], outDir = "./MCHelper")) %>% unlist()
    data
  })
  output$helitronData <- renderDT({
    data_with_progress <- helData()
    # data_with_progress$Status <- lapply(data_with_progress$Status, create_progress_bar)
    datatable(
      data_with_progress #%>% dplyr::select(Species, Strain, Genome_file, Raw_lib, Status)
      ,
      escape = FALSE, # Permite HTML en la tabla
      callback = JS(
        "table.on('click', 'td', function() {",
        "  var row = table.row(this).index();", # Obtener índice de la fila
        "  Shiny.setInputValue('helRowClick', row + 1);", # Enviar fila al servidor (1-indexed)
        "});"
      ),
      options = list(autoWidth = TRUE, dom = "t"), 
      selection = "single" # Permitir selección de una sola celda
      , rownames = FALSE
    ) %>%
      formatStyle(
        "Species",
        target = "row",
        backgroundColor = styleColorBar(c(0, 100), "#4caf50")
      )
  })
  helLibs <- reactiveValues()
  output$helianoRunLibs <- renderPlotly({
    tePlotLib(teLibList = list(helLibs$helRawLib))
  })
  observeEvent(input$helRowClick, {
    # req(input$inputData_cells_selected)
    # fileNames <- list() #debugging
    # selectedRow <- input$helitronData_rows_selected
    selectedRow <- input$helRowClick
    species <- helData()[selectedRow, "Species"]
    strain <- helData()[selectedRow, "Strain"]
    status <- helData()[selectedRow, "Status"]
    fileNames$teGenome <- helData()[selectedRow, "Genome_file"]
    fileNames$helRawLib <- helData()[selectedRow, "file_name"]
    # showModal(modalDialog(title = "You selected the row", paste(species, strain, fileNames$helRawLib)))
    # selectedRow <- 1 #debugging
    # species <- data[selectedRow, "Species"] #debugging
    # strain <- data[selectedRow, "Strain"] #debugging
    # status <- data[selectedRow, "Status"]#debugging
    # fileNames$teGenome <- data[selectedRow, "Genome_file"] #debugging
    # fileNames$helRawLib <- data[selectedRow, "file_name"] #debugging
    # # # #
    
    fileNames$helOutDir <- paste("./HELIANO", species, strain, sep = "/")
    # # fileNames$teRepeatLib <- paste("./RM2_output", species, strain, "N2_test-families.fa", sep = "/") #debugging
    helLibs$helRawLib <- teReadLib(fileNames$helRawLib, libIdentifier = "Raw lib")
    # fileNames$teCleanLib <- paste0(fileNames$teOutDir, "/", strain, "-clean_families.fa")
    # teLibs$teCleanLib <- teReadLib(fileNames$teCleanLib, libIdentifier = "Clean lib")
    # 
    # fileNames$teAutoCuratedLib <- paste0(fileNames$teOutDir, "/curated_sequences_NR.fa")
    # teLibs$teAutoCuratedLib <- teReadLib(fileNames$teAutoCuratedLib, libIdentifier = "Curated sequences")
    # 
    fileNames$helManualOutDir <- paste0(fileNames$helOutDir,"/MCH_manual")
    # 
    # fileNames$teCompleteModelsLib <- paste0(fileNames$teManualOutDir, "/complete_models.fa")
    # teLibs$teCompleteModelsLib <- teReadLib(fileNames$teCompleteModelsLib, libIdentifier = "Complete models")
    # fileNames$teIncompleteModelsLib <- paste0(fileNames$teManualOutDir, "/incomplete_models.fa")
    # teLibs$teIncompleteModelsLib <- teReadLib(fileNames$teIncompleteModelsLib, libIdentifier = "Incomplete models")
    # # fileNames$teManualSplittedLib <- paste0(fileNames$teManualOutDir, "/splitted_curated_sequences_NR.fa")
    # fileNames$teAutoCuratedTmpLib <- paste0(fileNames$teManualOutDir, "/tmp_curated_sequences_NR.fa")
    # teLibs$teAutoCuratedTmp <- teReadLib(fileNames$teAutoCuratedTmpLib)
    # fileNames$teManualDiscardedLib <- paste0(fileNames$teManualOutDir, "/discarded_sequences.fa")
    # 
    # fileNames$te80Lib <- paste0(fileNames$teManualOutDir, "/complete_808080.fa")
    # 
    # output$teNonCurated <- renderPlotly({
    #   tePlotLib(teLibList = list(teLibs$teRepeatLib, teLibs$teCleanLib, teLibs$teAutoCuratedLib))#
    # })
    # # fileNames$te80Lib <- gsub(fileNames$teCompleteModelsLib, pattern = "complete_models.fa", replacement = "complete_808080.fa")
    # teLibs$te80Lib <- teReadLib(fileNames$te80Lib, libIdentifier = "808080 lib")
    # 
    # fileNames$teNon80Lib <- gsub(fileNames$teCompleteModelsLib, pattern = "complete_models.fa", replacement = "complete_non808080.fa")
    # fileNames$te80OutDir <- gsub(fileNames$teNon80Lib, pattern = "complete_non808080.fa", replacement = "2_MI_MCH")
    # fileNames$teNon80TmpLib <- paste0(fileNames$te80OutDir, "/tmp_non808080.fa")
    # fileNames$te80RecovLib <- gsub(fileNames$teCompleteModelsLib, pattern = "complete_models.fa", replacement = "complete_808080_recovered.fa")
    # teLibs$te80RecovLib <- teReadLib(fileNames$te80RecovLib, libIdentifier = "Complete 808080 recovered")
    # fileNames$teStandbyLib <- paste0(fileNames$te80OutDir, "/standby_sequences.fa")
    # teLibs$teStandbyLib <- teReadLib(fileNames$teStandbyLib, libIdentifier = "Stand by lib")
    # 
    # fileNames$te70InputDb <- gsub(pattern = "standby_sequences.fa", replacement = "allDatabases.clustered_merged.fa", x = fileNames$teStandbyLib)
    # fileNames$te70Lib <- gsub(pattern = "standby_sequences.fa", replacement = "standby_707070.fa", x = fileNames$teStandbyLib)
    # teLibs$te70Lib <- teReadLib(fileNames$te70Lib, libIdentifier = "707070 lib")
    # fileNames$te70FilterLog <- gsub(fileNames$te70Lib, pattern = "standby_707070.fa", replacement = "707070_filter.log")
    # 
    # fileNames$teIncompleteNon80 <- gsub(fileNames$teIncompleteModelsLib, pattern = "incomplete_models.fa", replacement = "incomplete_non808080.fa")
    # fileNames$teIncompleteNewfam <- gsub(x = fileNames$teIncompleteNon80, pattern = "incomplete_non808080.fa", replacement = "3_MI_MCH/incomplete_newfam.fa")
    # fileNames$teIncomopleteRecovLib <- gsub(x = fileNames$teIncompleteNon80, pattern = "incomplete_non808080.fa", replacement = "incomplete_808080_recovered.fa")
    # teLibs$teIncompleteRecovLib <- teReadLib(fileNames$teIncomopleteRecovLib, libIdentifier = "Incomplete 808080 recovered")
    # fileNames$teIncompleteNon80Tmp <- gsub(pattern = "incomplete_newfam.fa", replacement = "tmp_incomplete_non808080.fa", fileNames$teIncompleteNewfam)
    # 
    # fileNames$te70FiltLib <- gsub(pattern = "standby_707070.fa", replacement = "standby_707070_filtered.fa",fileNames$te70Lib)
    # teLibs$te70FiltLib <- teReadLib(fileNames$te70FiltLib, libIdentifier = "Standby 707070 filtered")
    # 
    # fileNames$teIncompleteOutDir <- paste0(fileNames$teManualOutDir, "/3_MI_MCH")
    # 
    # fileNames$teFinalMchDir <- paste0(fileNames$teOutDir, "/MCH_final")
    
    updateTabsetPanel(session, inputId = "helianoTabs"
                      , selected = "helianoRun"#ifelse(grepl(as.character(status), x = "1234"), yes = "1234"
                                          #, no = ifelse(grepl(as.character(status), x = "56"), yes = "56", no = as.character(status)))
    )
  })
  observeEvent(input$helianoRun, {
    tryCatch({
      mchLogFile <- readLines(paste0(fileNames$helOutDir, "_mchelper-manual.log"))
      if(mchLogFile[length(mchLogFile)] == "MCHelper with TE-Aid successfully run") {
        #MCH_TE-aid has finished
        updateSelectizeInput(session, inputId = "helManualSeq", choices = names(helLibs$helRawLib))#fileNames$helRawLib
      }
    }
    , error = function(e) {}
    )
  })
  
  observeEvent(input$helManualSeq, {
    # input$helManualSeq <- names(helLibs$helRawLib)[2] # debugging
    req(input$helManualSeq)
    seqName <- input$helManualSeq %>% strsplit("#", fixed = TRUE) %>% unlist() %>% head(1)
    helAnnot <- input$helManualSeq %>% strsplit("#", fixed = TRUE) %>% lapply("[[", 2) %>% unlist()
    helBlastGenome <- list.files(paste0(fileNames$helOutDir, "/MCH_manual"), pattern = paste0(seqName,".fa.genome.blastn"), recursive = TRUE, full.names = TRUE)
    helBlastSelf <- list.files(paste0(fileNames$helOutDir, "/MCH_manual"), pattern = paste0(seqName,".fa.self.blastn"), recursive = TRUE, full.names = TRUE)
    
    
    # output$teManualCuratedClass <- renderTable({
    #   teTable %>% subset(Seq_name == seqName) %>% select(Seq_name, length, class, order, sFamily, struct, coding, other, strand)
    # })
    # output$teManualClass <- renderTable(
    #   teClassOutput(teTable, seqName)$teInfo
    # )
    # output$teManualCoding <- renderTable(
    #   teClassOutput(teTable, seqName)$teCoding, colnames = FALSE
    # )
    # output$teManualStructure <- renderTable(
    #   teClassOutput(teTable, seqName)$teStruc
    # )
    
    
    output$helConsDivergence <- renderPlot({
      teConsDivergencePlot(seqName, blastFiles = helBlastGenome)#, leftTrim = input$leftTrim, rightTrim = input$rightTrim
    })
    output$helConsCoverage <- renderPlotly({
      teConsCoveragePlot(seqName, blastFiles = helBlastGenome, leftTrim = 1, rightTrim = 1, source = "helCons")
    })
    output$helConsSelf <- renderPlot({
      teConsSelfPlot(seqName, blastFiles = helBlastSelf)#, leftTrim = input$leftTrim, rightTrim = input$rightTrim
    })
    output$helConsStructure <- renderPlot({
      teConsStructurePlot(seqName, blastFiles = helBlastSelf)# , leftTrim = input$leftTrim, rightTrim = input$rightTrim
    })
    
  })
  # mkdir -p HELIANO/C.elegans
  # conda run -n HELIANO --no-capture-output heliano -g 0_raw/C.elegans/N2/GCF_000002985.6_WBcel235_genomic.fna -o HELIANO/C.elegans/N2/ &> HELIANO/C.elegans/N2_heliano.log &
  # conda run -n HELIANO --no-capture-output heliano -g 0_raw/C.briggsae/VX34/GCA_022453885-chromosomes.fasta -o HELIANO/C.briggsae/VX34/ &> HELIANO/C.briggsae/VX34_heliano.log &
  # conda run -n MCHelper --no-capture-output cd-hit-est -i HELIANO/C.elegans/N2/RC.representative.fa -o HELIANO/C.elegans/N2/RC.representative_NR.fa -c 0.95 -aS 0.98 -G 0 -g 1 -b 500 -M 0 -d 0 &> HELIANO/C.elegans/N2/cd-hit-est.log
  # cp MCHelper/C.elegans/N2/MCH_final/final_curated_NR.fa HELIANO/C.elegans/N2/
  # "conda run -n MCHelper makeblastdb -in HELIANO/C.elegans/N2/final_curated_NR.fa -dbtype nucl"
  # conda run -n MCHelper --no-capture-output blastn -query HELIANO/C.elegans/N2/RC.representative_NR.fa -db HELIANO/C.elegans/N2/final_curated_NR.fa -num_threads 32 -outfmt 6 -qcov_hsp_perc 80 -perc_identity 80 -max_hsps 1 -out HELIANO/C.elegans/N2/hel-NR_vs_final-curated-NR.blast
  # conda run -n MCHelper --no-capture-output blastn -query HELIANO/C.elegans/N2/RC.representative_NR.fa -db MCHelper/C.elegans/N2/1_MI_MCH/2_MI_MCH/allDatabases.clustered_merged.fa -num_threads 32 -outfmt 6 -qcov_hsp_perc 80 -perc_identity 80 -max_hsps 1 -out HELIANO/C.elegans/N2/hel-NR_vs_allDB-merged.blast
  # cat mchelper-ats/db/allDatabases.clustered_rename.fa HELIANO/C.elegans/N2/final_curated_NR.fa > HELIANO/C.elegans/N2/allDatabases_final-NR_merged.fa
  # conda run -n MCHelper --no-capture-output makeblastdb -in HELIANO/C.elegans/N2/allDatabases_final-NR_merged.fa -dbtype nucl
  # conda run -n MCHelper --no-capture-output blastn -query HELIANO/C.elegans/N2/RC.representative_NR.fa -db HELIANO/C.elegans/N2/allDatabases_final-NR_merged.fa -num_threads 32 -outfmt 6 -qcov_hsp_perc 80 -perc_identity 80 -max_hsps 1 -out HELIANO/C.elegans/N2/hel-NR_vs_allDB-final.blast 
  # conda run -n MCHelper --no-capture-output blastn -query HELIANO/C.elegans/N2/RC.representative_NR.fa -db mchelper-ats/db/allDatabases.clustered_rename.fa -num_threads 32 -outfmt 6 -qcov_hsp_perc 80 -perc_identity 80 -max_hsps 1 -out HELIANO/C.elegans/N2/hel-NR_vs_allDB.blast 
  # sed '/^>/ s/$/#Unclassified/' HELIANO/C.elegans/N2/RC.representative.fa > HELIANO/C.elegans/N2/RC.representative-clean.fa
  # conda run -n MCHelper --no-capture-output python3 mchelper-ats/MCHelper.py -r A -a F --input_type fasta -l HELIANO/C.elegans/N2/RC.representative-clean.fa -g 0_raw/C.elegans/N2/GCF_000002985.6_WBcel235_genomic.fna -o HELIANO/C.elegans/N2/MCH_auto -c 1 -t 20 -b BUSCO_libs/nematoda_odb10_ALL.hmm -v Y &> HELIANO/C.elegans/N2_mchelper-automatic.log &
  # previous command results in 13 sequences curated automatically, 11 classified as HELITRON, 1 as LTR and 1 as TIR. I am not sure if MCHelper automatic is a good aproach as it removes too many sequences
  # conda run -n MCHelper --no-capture-output python3 mchelper-ats/MCHelper.py -r A -a F --input_type fasta -l HELIANO/C.elegans/N2/RC.representative-clean.fa -g 0_raw/C.elegans/N2/GCF_000002985.6_WBcel235_genomic.fna -o HELIANO/C.elegans/N2/MCH_auto -c 1 -t 20 -b BUSCO_libs/nematoda_odb10_ALL.hmm -v Y &> HELIANO/C.elegans/N2_mchelper-automatic.log &
  # sed '/^>/ s/$/#Unclassified/' HELIANO/C.elegans/N2/RC.representative_NR.fa > HELIANO/C.elegans/N2/RC.representative_NR-clean.fa
  # conda run -n MCHelper --no-capture-output python3 mchelper-ats/MCHelper.py -r T --input_type fasta -l HELIANO/C.elegans/N2/RC.representative_NR-clean.fa -g 0_raw/C.elegans/N2/GCF_000002985.6_WBcel235_genomic.fna -o HELIANO/C.elegans/N2/MCH_manual -c 1 -t 20 -b BUSCO_libs/nematoda_odb10_ALL.hmm -v Y &> HELIANO/C.elegans/N2_mchelper-manual.log &
}

shinyApp(
  ui = mchelper_ui,
  server = mchelper_server, options = list(host="0.0.0.0", port = 3001)
)



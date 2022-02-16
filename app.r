list.of.packages <- c("shiny", "shinyalert", "shinyFiles", "shinydashboard", "DT", "tools")
# TODO: 4D, FIT i PREDICT
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
# shiny
library(shiny)
library(shinyalert)
library(shinydashboard)
library(DT)
library(tools)
library(shinyFiles)
source('mripredict.r')
source('mripredict_cv.R')
source('mripredict_fit.r')
source('mripredict_predict.r')
source('mripredict_library.r')
source('mripredict_cv.R')
source('mripredict_library.r')
source('mripredict.r')
source('rotation3d.r')
source('cv_functions.r')
source('predict.r')
source('fit.r')
source('combat_quim.R')
source('combat_utils.R')
source('glmnet.utils.R')
source('extra.R')
suppressWarnings(options(warn=-1))

mri_files = data.frame(files=NULL);
families = c('Binomial', 'Gaussian', 'Cox')
families_short = c('binomial', 'gaussian', 'cox')
modulations = c('1 image (GM or WM)', '2 images (GM/WM + GM/WM modulated)', '4 images (GM + GM modulated + WM + WM modulated)', 'Only clinical')
modulations_short = c('un','fu','all','clinical')
response <- predictors <- modulation <- c('None')

####################################################################################################
####################################################################################################
version=1.0
version_online=tryCatch(
  {
    read.csv('https://drive.google.com/uc?export=download&id=1LAKstd-hfj39WzGTGCpACpxI8_UMmOHA')
  },
  error = function(e) version
)

cat('\nCurrent version: MRIPredict', version)
cat('\nLatest version available: MRIPredict', version_online$versio)
response_family = 'binomial' 
folder = 'output'
modulation='un'
if(!dir.exists(folder)) dir.create(folder)
covariates = c('None')
####################################################################################################
####################################################################################################

# # Define server logic to summarize and view selected dataset ----
server <- function(input, output, session) {
  global <- reactiveValues(datapath = getwd())
  if(version_online>version){
    showModal(modalDialog(
      title = sprintf("New version available. Please update the software."),
      HTML(sprintf("New version: %s <br>
                    You have:    %s. <br>You can find the new version <a href='https://mripredict.com'>here</a>"))
    ))
  }
  
  out_table <- reactiveValues()
  volumes = getVolumes()
  
  # select nifti files
  # shinyFileChoose(input, 'files', roots=volumes, defaultPath = getwd(), filetype=c('nii','gz'),
  #                 session = session)
  # mri_files <- reactive({
  #   as.character(parseFilePaths(roots = volumes, input$files)$datapath)
  # })
  mri_files <- reactiveValues(paths=NULL)
  observe({
    if(!is.null(input$gm_un)){
      mri_files$paths_gm_un <- as.character(parseFilePaths(roots = volumes, input$gm_un)$datapath)
      output$gm_selected <- renderText({sprintf('Selected: %s images.', length(mri_files$paths_gm_un))})
      #out_table$GM_unmodulated <- mri_files$paths_gm_un
    }
    if(!is.null(input$gm_fu)){
      mri_files$paths_gm_fu <- as.character(parseFilePaths(roots = volumes, input$gm_fu)$datapath)
      output$gmfu_selected <- renderText({sprintf('Selected: %s images.', length(mri_files$paths_gm_fu))})
      #out_table$GM_modulated <- mri_files$paths_gm_fu
    }
    if(!is.null(input$wm_un)){
      mri_files$paths_wm_un <- as.character(parseFilePaths(roots = volumes, input$wm_un)$datapath)
      output$wm_selected <- renderText({sprintf('Selected: %s images.', length(mri_files$paths_wm_un))})
      #out_table$WM_unmodulated <- mri_files$paths_wm_un
    }
    if(!is.null(input$wm_fu)){
      mri_files$paths_wm_fu <- as.character(parseFilePaths(roots = volumes, input$wm_fu)$datapath)
      output$wmfu_selected <- renderText({sprintf('Selected: %s images.', length(mri_files$paths_wm_fu))})
      #out_table$WM_modulated <- mri_files$paths_wm_fu
    }
    if(!is.null(input$modulation))
      mri_files$modulation <- modulations_short[which(input$modulation==modulations)]
    if(!is.null(input$modality)){
      mri_files$response_family <- families_short[which(input$modality==families)]
      if(mri_files$response_family == 'cox'){
        mri_files$time = input$response_time
        mri_files$status = input$response_status
      }
    }
    if(!is.null(input$response)){
      mri_files$response <- input$response
      
    }
    if(!is.null(input$covariates))
      mri_files$covariates <- input$covariates
    if(!is.null(input$predictors))
      mri_files$predictors <- input$predictors
    if(!is.null(input$ensemble_learning))
      mri_files$ensemble_learning <- input$ensemble_learning
    if(!is.null(input$model))
      mp <- input$model
    if(!is.null(input$clinical)){
      output$txt_file <- renderText({as.character(mri_files$path_clinical)})
    }
    if(!is.null(input$load_model)) {
      output$txt_model <- renderText({as.character(mri_files$path_model)})
    }
    
  })
  observeEvent(input$folds, {
    mri_files$folds <- input$folds
  })
  
  # select nifti directory
  shinyDirChoose(input, 'dirs', roots=volumes) #volumes)
  
  # clinical file
  shinyFileChoose(input, 'clinical', roots=volumes,  filetype=c('txt','xlsx','xls','csv'),
                  session = session)
  
  
  shinyFileChoose(input, 'gm_un', roots=volumes, filetype=c('nii','gz'),
                  session = session)
  shinyFileChoose(input, 'gm_fu', roots=volumes,  filetype=c('nii','gz'),
                  session = session)
  shinyFileChoose(input, 'wm_un', roots=volumes, filetype=c('nii','gz'),
                  session = session)
  shinyFileChoose(input, 'wm_fu', roots=volumes,  filetype=c('nii','gz'),
                  session = session)
  
  
  
  shinyFileChoose(input, 'wm_un', roots=volumes,  filetype=c('nii','gz'),
                  session = session)
  shinyFileChoose(input, 'wm_fu', roots=volumes,  filetype=c('nii','gz'),
                  session = session)
  observe({
    shinyFileChoose(input, 'load_model', roots=volumes(),  filetype=c('rds'), session= session)
    req(input$load_model)
    if(!is.null(input$load_model)){
      mri_files$path_model <- as.character(parseFilePaths(volumes, input$load_model)$datapath)
      if(length(mri_files$path_model)>0) {
        mri_files$mp <- readRDS(mri_files$path_model)
        mp <- mri_files$mp
        
        updateSelectInput(session, "modulation",
                          label = "Input images",
                          choices = modulations,
                          selected = modulations[which(mp$modulation==modulations_short)])
        updateSelectInput(session, "modality",
                          label = "Family",
                          choices = families,
                          selected = families[which(mp$response_family==families_short)])
        updateSelectInput(session, "response",
                          label = "Response",
                          choices = mp$response_var,
                          selected = mp$response_var)
        updateSelectInput(session, "response_time",
                          label = "Time",
                          choices = mp$response_time,
                          selected = mp$response_time)
        updateSelectInput(session, "response_status",
                          label = "Status",
                          choices = mp$response_status,
                          selected = mp$response_status)
        updateCheckboxGroupInput(session, 'covariates', 
                                 label='Covariates',
                                 choices=mp$covX_var,
                                 selected = mp$covX_var)
        
        updateCheckboxGroupInput(session, 'predictors', 
                                 label='Predictors',
                                 choices=mp$pred_var,
                                 selected = mp$pred_var)
        updateCheckboxInput(session,inputId = 'ensemble_learning',value = mp$n_iterations>1)
        
        # update model representation
        variables_used<-.most_frequent_variables(model_list = mp$models,
                                                 mp = mp)
        mri_files$variables_used <- variables_used
        
      }
    }
  })
  
  
  # EVENTS
  # quan s'apreti el botó de CV s'executarà això
  observeEvent(input$run_cv, {
    showModal(modalDialog(
      tagList(
        textInput("saveName", label = "Save name", placeholder = "test1")
      ), 
      title="Select a name for your project",
      footer = tagList(actionButton("confirmCV", "Run"),
                       modalButton("Cancel")
      )
    ))
    
  })
  observeEvent(input$run_fit, {
    
    showModal(modalDialog(
      tagList(
        textInput("saveName", label = "Save name", placeholder = "test1")
      ), 
      title="Select a name for your project",
      footer = tagList(actionButton("confirmFit", "Run"),
                       modalButton("Cancel")
      )
    ))
    
  })
  observeEvent(input$run_predict, {
    
    showModal(modalDialog(
      tagList(
        textInput("saveName", label = "Save name", placeholder = "test1")
      ), 
      title="Select a name for your project",
      footer = tagList(actionButton("confirmPredict", "Run"),
                       modalButton("Cancel")
      )
    ))
    
  })
  
  # CALL: CV  
  observeEvent(input$confirmCV, {
    removeModal()
    save_name <- sprintf('output/%s',input$saveName)
    if(!dir.exists(save_name)) dir.create(save_name)
    save_name <- sprintf('%s/%s',save_name, input$saveName)
    if(file.exists(sprintf("%s_list_folds_n%s.txt", save_name,mri_files$folds))){
      showModal(modalDialog(
        tagList(
          selectInput("nameExists", label = "A project with the same name exists", choices = list.files(pattern="*.txt"))
        ), 
        title="Change name or continue",
        footer = tagList(actionButton("run_cv", "Rename"),
                         modalButton("Continue")
        )
      ))
    }
    print(modulation)
    print(mri_files$covariates)
    if(is.null(mri_files$covariates) || mri_files$covariates == 'None') mri_files$covariates <- ''
    if(is.null(mri_files$predictors) || mri_files$predictors == 'None') mri_files$predictors <- ''
    if(is.null(mri_files$ensemble_learning)) mri_files$ensemble_learning <- F
    information_variables = c(mri_files$covariates, mri_files$predictors)
    print(information_variables)
    print(sprintf("%s_list_folds_n%s.txt", input$saveName,mri_files$folds))
    if(mri_files$response_family=='cox'){
      mri_files$response=c("time","status")#c(mri_files$time,mri_files$status)
    }
    paths = cbind(mri_files$paths_gm_un, mri_files$paths_gm_fu, mri_files$paths_wm_un, mri_files$paths_wm_fu)
    mp <- mripredict(paths, mri_files$clinical, mri_files$response, mri_files$covariates, mri_files$predictors, mri_files$response_family, mri_files$modulation, information_variables = information_variables)
    withProgress(message = "Performing Cross-validation", value = 0, {
      
      mp <- mripredict_cv(mp, space = "NO_CHECK", save_name = save_name, folds_file = "", n_cores = 1, n_folds = mri_files$folds,
                          use_significant_voxels = FALSE, use_ensemble_learning = mri_files$ensemble_learning, use_ensemble_voxels = mri_files$ensemble_learning, use_ensemble_subjects = FALSE, ide_shiny = TRUE)
    })
    mri_files$results <- mp$cv_results
    write.csv(mp$cv_results, file = sprintf('output/%s/cv_results.csv', input$saveName), row.names = F)
    mri_files$metrics <- mp$metrics
    write.csv(mp$metrics, file=sprintf('output/%s/cv_metrics.txt', input$saveName), row.names = F)
    variables_used<-.most_frequent_variables(model_list = mp$models,
                                             mp = mp,
                                             file = sprintf("%s_betas_summary.csv", save_name))
    mri_files$variables_used <- variables_used
    updateTabsetPanel(session, "allTabs",
                      selected = "Results")
  })
  
  # CALL: FIT 
  observeEvent(input$confirmFit, {
    removeModal()
    save_name <- sprintf('output/%s',input$saveName)
    if(!dir.exists(save_name)) dir.create(save_name)
    save_name <- sprintf('%s/%s',save_name, input$saveName)
    if(file.exists(sprintf("%s_list_folds_n%s.txt", input$saveName,mri_files$folds))){
      showModal(modalDialog(
        tagList(
          selectInput("nameExists", label = "A project with the same name exists", choices = list.files(pattern="*.rds"))
        ), 
        title="Change name or continue",
        footer = tagList(actionButton("run_fit", "Rename"),
                         modalButton("Continue")
        )
      ))
    }   
    print(mri_files$covariates)
    if(is.null(mri_files$covariates) || mri_files$covariates == 'None') mri_files$covariates <- ''
    if(is.null(mri_files$predictors) || mri_files$predictors == 'None') mri_files$predictors <- ''
    if(is.null(mri_files$ensemble_learning)) mri_files$ensemble_learning = F
    information_variables = c(mri_files$covariates, mri_files$predictors)
    print(information_variables)
    print(sprintf("%s_list_folds_n%s.txt", input$saveName,mri_files$folds))
    
    paths = cbind(mri_files$paths_gm_un, mri_files$paths_gm_fu, mri_files$paths_wm_un, mri_files$paths_wm_fu)
    mp <- mripredict(paths, mri_files$clinical, mri_files$response, mri_files$covariates, mri_files$predictors, mri_files$response_family, mri_files$modulation, information_variables = information_variables)
    
    withProgress(message = "Performing the training of the model", value = 0, {
      mp <- mripredict_fit(mp = mp, space = "NO_CHECK", n_cores=1, use_ensemble_learning = mri_files$ensemble_learning, use_ensemble_voxels = mri_files$ensemble_learning)
    })
    mri_files$results <- mp$cv_results
    saveRDS(mp, sprintf("%s_model.rds",save_name))
    variables_used<-.most_frequent_variables(model_list = mp$models,
                                             mp = mp,
                                             file = sprintf("%s_betas_summary.csv", save_name))
    mri_files$variables_used <- variables_used
    mri_files$mp <- mp
  })
  
  # CALL: PREDICT
  observeEvent(input$confirmPredict, {
    removeModal()
    if(file.exists(sprintf("%s_list_folds_n%s.txt", input$saveName,mri_files$folds))){
      showModal(modalDialog(
        tagList(
          selectInput("nameExists", label = "A project with the same name exists", choices = list.files(pattern="*.rds"))
        ), 
        title="Change name or continue",
        footer = tagList(actionButton("run_predict", "Rename"),
                         modalButton("Continue")
        )
      ))
    }
    save_name <- sprintf('output/%s',input$saveName)
    if(!dir.exists(save_name)) dir.create(save_name)
    print(mri_files$covariates)
    if(is.null(mri_files$covariates) || mri_files$covariates == 'None') mri_files$covariates <- ''
    if(is.null(mri_files$predictors) || mri_files$predictors == 'None') mri_files$predictors <- ''
    information_variables = c(mri_files$covariates, mri_files$predictors)
    print(information_variables)
    print(sprintf("%s_list_folds_n%s.txt", input$saveName,mri_files$folds))
    
    paths = cbind(mri_files$paths_gm_un, mri_files$paths_gm_fu, mri_files$paths_wm_un, mri_files$paths_wm_fu)
    #mp <- mripredict(paths, mri_files$clinical, mri_files$response, mri_files$covariates, mri_files$predictors, mri_files$response_family, mri_files$modulation, information_variables = information_variables)
    mp <- mri_files$mp
    withProgress(message = "Calculating predictions", value = 0, {
      if(mp$modulation == 'fu'){
        mp <- mripredict_predict(mp = mp, mri_paths_file = mri_files$paths_gm_un, mri_fu_paths_file = mri_files$paths_gm_fu, data_table_file = mri_files$clinical, space='NO_CHECK', n_cores=1)
      } else {
        mp <- mripredict_predict(mp = mp, mri_paths_file = paths, data_table_file = mri_files$clinical, space='NO_CHECK', n_cores=1)
      }
    })
    mri_files$results <- cbind(mp$prediction,paths,mri_files$clinical)
    write.csv(x = mri_files$results, file = sprintf('%s_pred.csv',save_name))
    if(mp$response_family == 'binomial') {
      colnames(mri_files$results)[1] <- sprintf("Prediction (0 = %s; 1 = %s)",mp$response_ref, mp$response_event)
    } else {
      colnames(mri_files$results)[1] <- 'Prediction'
    }
    updateTabsetPanel(session, "allTabs",
                      selected = "Results")
  })
  
  
  # reactive UI values
  observeEvent(input$clinical, {
    inFile <- parseFilePaths(roots=volumes,input$clinical)
    mri_files$path_clinical <- as.character(inFile$datapath)
    if(NROW(inFile)){
      if(file_ext(as.character(inFile$datapath))=='csv'){
        mri_files$clinical <- read.table(as.character(inFile$datapath), header=TRUE, sep=',')
      }
      else{
        mri_files$clinical <- read.table(as.character(inFile$datapath), header=TRUE)
      }
      updateSelectInput(session, 'response', choices = colnames(mri_files$clinical))
      updateSelectInput(session, 'response_time', choices = colnames(mri_files$clinical))
      updateSelectInput(session, 'response_status', choices = colnames(mri_files$clinical))
      updateCheckboxGroupInput(session, 'covariates', choices = colnames(mri_files$clinical))
      updateCheckboxGroupInput(session, 'predictors', choices = colnames(mri_files$clinical))
    }
  })
  # hide or show response time/status or response 
  observeEvent(input$response_family, {
    if(mri_files$response_family=='cox'){
      
    } else {
      
    }
  })
  modulation <- eventReactive(input$modulation, {
    modulations_short[which(input$modulation==modulations)]
  })
  # mri_files$response_family <- eventReactive(input$modality, {
  #   browser()
  #   families_short[which(input$modality==families)]
  # })
  # outputs
  output$files <- renderPrint({mri_files})
  output$clinical_table <- DT::renderDataTable(mri_files$clinical)
  output$table <- DT::renderDataTable(cbind(mri_files$paths_gm_un, mri_files$paths_gm_fu, mri_files$paths_wm_un, mri_files$paths_wm_fu))
  output$results <- DT::renderDataTable(mri_files$results)
  output$metrics <- DT::renderDataTable(mri_files$metrics)
  output$txt_model <- renderPrint(mri_files$model_path)
  output$model <- DT::renderDataTable(mri_files$variables_used)
}

ui <- dashboardPage(
  skin="red",
  #useShinyalert(),
  dashboardHeader(title=span(img(src="icon_back_w.png",height=50),"MRIPredict"),
                  
                  tags$li(
                    a(
                      span(icon('question-circle'),strong("ABOUT MRIPredict")),
                      height = 40,
                      href = "https://mripredict.com"
                    ),
                    class = "dropdown"
                  ),
                  tags$li(
                    a(
                      span(icon('envelope'),strong("Contact")),
                      height = 40,
                      href = "mailto:solanes@clinic.cat"
                    ),
                    class = "dropdown"
                  )
  ),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Configuration", tabName = "Configuration", icon = icon("flask")),
      menuItem("Data", tabName = "Data", icon = icon("table")),
      menuItem("Results", tabName = "Results", icon = icon("poll")),
      menuItem("Model", tabName = "Model", icon = icon("square-root-alt")),
      menuItem("Who are we?", tabName = "who", icon = icon("users"))
    )
  ),
  dashboardBody(
    tags$head(tags$style(HTML('
      .main-header .logo {
        font-family: "Georgia", Times, "Times New Roman", serif;
        font-weight: bold;
        font-size: 20px;
      }
    '))),
    tabItems(
      # 1st tab
      tabItem(tabName = "Configuration", 
              fluidRow(
                box(
                  h4("Configuration"),
                  selectInput('modulation','Input images',modulations),
                  selectInput('modality','Family',families),
                  conditionalPanel(
                    condition = "input.modality == 'Cox'",
                    selectInput('response_time','Time',response),
                    selectInput('response_status','Status',response)
                  ),
                  conditionalPanel(
                    condition = "input.modality != 'Cox'",
                    selectInput('response','Response',response)
                  ),
                  #selectInput('response','Response',response),
                  #selectInput('response_time','Time',response),
                  #selectInput('response_status','Status',response),
                  checkboxInput(
                    inputId = 'ensemble_learning',
                    label = 'Perform ensemble learning',
                    value = FALSE
                  ),
                  checkboxGroupInput(inputId = 'covariates',label = 'Covariates',choices = covariates),
                  checkboxGroupInput(inputId = 'predictors',label = 'Predictors',choices = predictors)
                )
                ,
                box(
                  h4("Input images"),
                  conditionalPanel(
                    condition = "input.modulation != 'Only clinical'",
                    shinyFilesButton('gm_un', label='Gray matter', title='Please select GM unmodulated files', multiple=TRUE),
                    h6(textOutput("gm_selected"))
                  ),
                  conditionalPanel(
                    condition = "input.modulation != 'Only clinical' && input.modulation != '1 image (GM or WM)'",
                    shinyFilesButton('gm_fu', label='Gray matter modulated', title='Please select GM modulated files', multiple=TRUE),
                    h6(textOutput("gmfu_selected"))
                  ),
                  conditionalPanel(
                    condition = "input.modulation == '4 images (GM + GM modulated + WM + WM modulated)'",
                    shinyFilesButton('wm_un', label='White matter', title='Please select WM unmodulated files', multiple=TRUE),
                    h6(textOutput("wm_selected")),
                    shinyFilesButton('wm_fu', label='White matter modulated', title='Please select WM modulated files', multiple=TRUE),
                    h6(textOutput("wmfu_selected"))
                  )
                  # , shinyFilesButton('gm_fu', label='Gray matter modulated', title='Please select GM modulated files', multiple=TRUE),
                  # h6(textOutput("gmfu_selected")),
                  # shinyFilesButton('wm_un', label='White matter', title='Please select WM unmodulated files', multiple=TRUE),
                  # h6(textOutput("wm_selected")),
                  # shinyFilesButton('wm_fu', label='White matter modulated', title='Please select WM modulated files', multiple=TRUE),
                  # h6(textOutput("wmfu_selected"))
                ),
                #actionButton("files", "File select", class = "btn-lg btn-success"),
                #shinyDirButton('dirs',label='Folder containing images', title='Please select the folder where the images are stored'),
                box(
                  h4("Clinical information"),
                  shinyFilesButton('clinical', label='Clinical file', title='Please select the clinical file', multiple=FALSE),
                  h6('File selected:'),
                  textOutput("txt_file")
                ),
                br(),br(),
                box(
                  sliderInput("folds", "Number of folds:", min = 3, max = 10, value = 10),
                  actionButton("run_cv", "Cross-validation", class = "btn-lg btn-success")
                ),
                box(
                  #actionButton("run_cv", "Cross-validation", class = "btn-lg btn-success"),
                  actionButton("run_fit", "Train model", class = "btn-lg btn-success"),
                  #actionButton("load_model", "Load model", class = "btn-lg btn-success"),
                  shinyFilesButton('load_model', label='Load model', title='Select a .rds model', multiple=FALSE),
                  textOutput("txt_model"),
                  actionButton("run_predict", "Apply model", class = "btn-lg btn-success")
                )
                
                
                
              ))
      ,
      # 2nd tab
      tabItem(tabName = "Data",
              fluidRow(
                box(
                  
                  title = "Images",
                  DT::dataTableOutput('table')),
                box(
                  title = "Clinical",
                  DT::dataTableOutput('clinical_table'))
              )),
      # 3rd tab
      tabItem(tabName = "Results",
              fluidRow(
                column(12,
                       h4('Cross-validation metrics'),
                       DT::dataTableOutput('metrics')),
                column(3,
                       h4('Results'),
                       DT::dataTableOutput('results'))
              )),
      # 4th tab
      tabItem(tabName = "Model",
              fluidRow(
                column(12,
                       h4('Model information'),
                       DT::dataTableOutput('model'))
              )),
      tabItem(tabName = "who",
              fluidRow(
                box(span("This software has been developed by Joaquim Radua and Aleix Solanes from ",strong(a("IMARD Group", href="https:imardgroup.com")), " - IDIBAPS, Hospital Clinic de Barcelona")),
                box(
                  strong(a("Joaquim Radua", href='https://pubmed.ncbi.nlm.nih.gov/?term=Radua+J&cauthor_id=32454268'))
                ),
                box(
                  strong(a("Aleix Solanes", href='https://pubmed.ncbi.nlm.nih.gov/?term=Solanes+A&cauthor_id=33831461')),
                  br(),
                  "Mail:", a("solanes@clinic.cat", href="mailto:solanes@clinic.cat")
                )
              ))
    )
  )
)






# 
# 
# 
# tabsetPanel(id = "allTabs", type = "tabs",
#             tabPanel("Configuration", 
#                      fluidRow(
#                        column(3,
#                               h4("Configuration"),
#                               selectInput('modulation','Input images',modulations),
#                               selectInput('modality','Family',families),
#                               selectInput('response','Response',response),
#                               checkboxInput(
#                                 inputId = 'ensemble_learning',
#                                 label = 'Perform ensemble learning',
#                                 value = TRUE
#                               ),
#                               checkboxGroupInput(inputId = 'covariates',label = 'Covariates',choices = covariates),
#                               checkboxGroupInput(inputId = 'predictors',label = 'Predictors',choices = predictors)
#                               
#                        ),
#                        column(4,
#                               h4("Input images"),
#                               shinyFilesButton('gm_un', label='Select GM unmodulated files', title='Please select GM unmodulated files', multiple=TRUE),
#                               h6(textOutput("gm_selected")),
#                               shinyFilesButton('gm_fu', label='Select GM modulated files', title='Please select GM modulated files', multiple=TRUE),
#                               h6(textOutput("gmfu_selected")),
#                               shinyFilesButton('wm_un', label='Select WM unmodulated files', title='Please select WM unmodulated files', multiple=TRUE),
#                               h6(textOutput("wm_selected")),
#                               shinyFilesButton('wm_fu', label='Select WM modulated files', title='Please select WM modulated files', multiple=TRUE),
#                               h6(textOutput("wmfu_selected")),
#                               #actionButton("files", "File select", class = "btn-lg btn-success"),
#                               #shinyDirButton('dirs',label='Folder containing images', title='Please select the folder where the images are stored'),
#                               h4("Clinical information"),
#                               shinyFilesButton('clinical', label='Clinical file', title='Please select the clinical file', multiple=FALSE),
#                               h6('File selected:'),
#                               textOutput("txt_file"),   
#                               br(),br(),
#                               sliderInput("folds", "Number of folds:", min = 3, max = 10, value = 10),
#                               actionButton("run_cv", "Cross-validation", class = "btn-lg btn-success"),
#                               actionButton("run_fit", "Train model", class = "btn-lg btn-success"),
#                               #actionButton("load_model", "Load model", class = "btn-lg btn-success"),
#                               shinyFilesButton('load_model', label='Load model', title='Select a .rds model', multiple=FALSE),
#                               textOutput("txt_model"),
#                               actionButton("run_predict", "Apply model", class = "btn-lg btn-success")  
#                               
#                        ),
#                        column(3,
#                               img(src = "icon_back_w.png"))
#                      )),
#             tabPanel("Data", 
#                      fluidRow(
#                        column(12,
#                               tableOutput('table')),
#                        column(3,
#                               tableOutput('clinical_table'))
#                      )
#             ),
#             tabPanel("Results",
#                      fluidRow(
#                        column(12,
#                               h4('Cross-validation metrics'),
#                               tableOutput('metrics')),
#                        column(3,
#                               h4('Results'),
#                               tableOutput('results'))
#                      )
#             ),
#             tabPanel("Model",
#                      fluidRow(
#                        column(12,
#                               h4('Model information'),
#                               tableOutput('model'))
#                      )))
# ))

shinyApp(ui = ui, server = server)





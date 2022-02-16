# FIX MANIA
# S'HA EXECUTAT EN DUES TONGADES, I EL QUE HA FET AMB ELS RDS ÉS:
# 1era vegada: feia bé, anava guardant models un a continuació de l'altre a ritme de 360 per fold (20 imputacions * 18 talls)
# 2ona vegada: els índex de model recomençaven a 1
# TIME FROM MRI
# RELAPSE DEPRESSIO, RELAPSE MANIA, TOTS RELAPSES
# TIME FROM REMISSION
# RELAPSE DEPRESSIO, RELAPSE MANIA, TOTS RELAPSES
base_names = c('/home/aleix/Documents/papers_aleix/MRIPredict/MRIPredict/R/output_mania/mri_and_data/mania_totes_timefromMRI_relapsedepressio',
               '/home/aleix/Documents/papers_aleix/MRIPredict/MRIPredict/R/output_mania/mri_and_data/mania_totes_timefromMRI_relapsemania',
               '/home/aleix/Documents/papers_aleix/MRIPredict/MRIPredict/R/output_mania/mri_and_data/mania_totes_timefromMRI_totsrelapses',
               '/home/aleix/Documents/papers_aleix/MRIPredict/MRIPredict/R/output_mania/mri_and_data/mania_totes_timefromremission_relapsedepressio',
               '/home/aleix/Documents/papers_aleix/MRIPredict/MRIPredict/R/output_mania/mri_and_data/mania_totes_timefromremission_relapsemania',
               '/home/aleix/Documents/papers_aleix/MRIPredict/MRIPredict/R/output_mania/mri_and_data/mania_totes_timefromremission_totsrelapses'
)
for(base_name in base_names){
  #if(!file.exists(sprintf("%s_models.rds", base_name))){
    m=list()
    n_models_previous = 0
    out = list()
    
    for(i in 1:10) { # primer s'agafen els últims 360 models
      m_temp = readRDS(sprintf('%s_model_FOLD_%s_model_list.rds',base_name,i))
      n_models = length(m_temp)
      print(n_models)
      if(n_models_previous == n_models) break
      m[[i]]=m_temp[(n_models-359):n_models] 
      n_models_previous = n_models
    }
    
    k = 1 # indicador de fold
    for(j in i:10) { # aquí toca agafar els k*360:(k+1)*360 models
      m_temp = readRDS(sprintf('%s_model_FOLD_%s_model_list.rds',base_name,j))
      n_models = length(m_temp)
      m[[j]]=m_temp[(k*360-359):(k*360)] 
      k = k + 1
    }
    
    
    ind_out = 1
    for(i in 1:10){
      for(j in 1:length(m[[i]])){
        out[[ind_out]] = m[[i]][[j]]
        ind_out = ind_out + 1
      }
    }
    cat("\n Model size:",length(out))
    #saveRDS(m, sprintf("%s_models_folds",base_name))
    saveRDS(out, sprintf("%s_models.rds", base_name))
 # }
}

rds_names = paste0(base_names,'_models.rds')









#::::::::::::::::::::::::::::: IMPORTS :::::::::::::::::::::::::::::::::
# totes les proves
list.of.packages <- c("sys", "methods", "readxl", "tcltk")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, library, character.only=TRUE)
library(sys)
args = commandArgs(trailingOnly=TRUE)
options(warn=-1)

source('mripredict_cv_parallel.r')
source('mripredict_library.r')
source('mripredict.r')
source('rotation3d.r')
source('cv_functions.r')
source('predict.r')
source('fit.r')
source('combat_quim.R')
source('combat_utils.R')
source('glmnet-utils.R')
#::::::::::::::::::::::::::: FI IMPORTS :::::::::::::::::::::::::::::::

#::::::::::::::::::::::::::: COMMON PARAMETERS ::::::::::::::::::::::::
response_family = 'cox' 
modulation <- 'all'
covars=c()
if(!dir.exists('tmp')) dir.create(folder, recursive = TRUE)
# information_variables are the variables that contain information (only interesting when there are NAs, used in multiple imputation)


mri_paths_file = '~/Documents/dades/mania_totes/gm/subsamp/subsamp/subsamp/s4/mri_files.txt'
mri_fu_paths_file = '~/Documents/dades/mania_totes/gm_mod/subsamp/subsamp/subsamp/s4/mri_files.txt'
mri_wm_paths_file = '~/Documents/dades/mania_totes/wm/subsamp/subsamp/subsamp/s4/mri_files.txt'
mri_wm_fu_paths_file = '~/Documents/dades/mania_totes/wm_mod/subsamp/subsamp/subsamp/s4/mri_files.txt'
paths = cbind(read.table(mri_paths_file, stringsAsFactors = FALSE)$V1,
              read.table(mri_fu_paths_file, stringsAsFactors = FALSE)$V1,
              read.table(mri_wm_paths_file, stringsAsFactors = FALSE)$V1,
              read.table(mri_wm_fu_paths_file, stringsAsFactors = FALSE)$V1)
data_table_file = read.table("~/Documents/dades/mania_totes/data_mania_final.csv", sep=',', header = TRUE, na.strings = c('','NA','na','N/A'))
data_table_file$status_relapse_mania <- 0
data_table_file$status_relapse_mania[data_table_file$polaritat_relapse=='mania'] <- 1
data_table_file$status_relapse_depressio <- 0
data_table_file$status_relapse_depressio[data_table_file$polaritat_relapse=='depressio'] <- 1
data_table_file$time_from_MRI[data_table_file$time_from_MRI==0] = data_table_file$time_from_MRI+0.5
data_table_file$time_from_remission[data_table_file$time_from_remission==0] = data_table_file$time_from_remission+0.5

########################################################
########## PROVES DE MRI + DATA ########################
########################################################
folder = 'output_mania/mri_and_data'
if(!dir.exists(folder)) dir.create(folder, recursive = TRUE)


# mirem quines variables tenen massa pocs valors uniques
ab=apply(data_table_file,2,table)
# descartem les que no tinguin un minim de 2 per categoria
# Ham_3, Ham_19, Ham_21, PANNS_N2, PANNS_N7, PANNS_PG16

###################### 1. TIME FROM REMISSION, tots els relapses
cat('DOING MRI+DATA: 1 b. MANIA, TIME FROM REMISSION, tots els relapses \n')
predictor = c('TAP','EdatDebut','EdatPrimerIngres','duration_of_illness','AP.EQP.total_Andreasen','TEA',
              'HDRS17','PANSS_P_Total','PANSS_N_Total','PANSS_PG_Total','PANSS_Total','YMRS',
              'Ham_1','Ham_2',
              'Ham_3',# low 
              'Ham_4','Ham_5','Ham_6','Ham_7','Ham_8','Ham_9','Ham_10','Ham_11','Ham_12','Ham_13','Ham_14','Ham_15','Ham_16','Ham_17','Ham_18',
              'Ham_19',#
              'Ham_20',
              'Ham_21',#
              'PANNS_P1','PANNS_P2','PANNS_P3','PANNS_P4','PANNS_P5','PANNS_P6','PANNS_P7','PANNS_N1',
              #'PANNS_N2',
              #'PANNS_N6',#
              'PANNS_N3','PANNS_N4','PANNS_N5',
              'PANNS_N7',#
              'PANNS_PG1','PANNS_PG2','PANNS_PG3','PANNS_PG4','PANNS_PG5','PANNS_PG6','PANNS_PG7','PANNS_PG8','PANNS_PG9','PANNS_PG10','PANNS_PG11','PANNS_PG12','PANNS_PG13','PANNS_PG14','PANNS_PG15',
              'PANNS_PG16',#
              'Young_1','Young_2','Young_3','Young_4','Young_5','Young_6','Young_7','Young_8','Young_9','Young_10','Young_11')
response_var   = c('time_from_remission', 'status')
covariates = c('age','Dona')
information_variables = c(covariates, predictor)
save_name <- paste(folder, 'mania_totes_timefromremission_totsrelapses', sep="/")
tmp_file <- 'tmp/mania1b.txt'


    mp <- mripredict(paths, data_table_file, response_var, covariates, predictor, response_family, modulation, information_variables = information_variables)
    mp$models <- readRDS(rds_names[6])
    variables_used<-.most_frequent_variables(model_list = mp$models,
                                             mp = mp,
                                             file = sprintf("%s_betas_summary.csv", save_name))

###################### 2. TIME FROM REMISSION, relapse mania
cat('DOING MRI+DATA: 2 b. MANIA, TIME FROM REMISSION, relapse mania \n')
response_var   = c('time_from_remission', 'status_relapse_mania')
covariates = c('age','Dona')
information_variables = c(covariates, predictor)
save_name <- paste(folder, 'mania_totes_timefromremission_relapsemania', sep="/")

    mp <- mripredict(paths, data_table_file, response_var, covariates, predictor, response_family, modulation, information_variables = information_variables)
    mp$models <- readRDS(rds_names[5])
    variables_used<-.most_frequent_variables(model_list = mp$models,
                                             mp = mp,
                                             file = sprintf("%s_betas_summary.csv", save_name))

###################### 3. TIME FROM REMISSION, relapse depressiu
cat('DOING MRI+DATA: 3 b. MANIA, TIME FROM REMISSION, relapse depressio \n')
response_var   = c('time_from_remission', 'status_relapse_depressio')
covariates = c('age','Dona')
information_variables = c(covariates, predictor)
save_name <- paste(folder, 'mania_totes_timefromremission_relapsedepressio', sep="/")

    mp <- mripredict(paths, data_table_file, response_var, covariates, predictor, response_family, modulation, information_variables = information_variables)
    mp$models <- readRDS(rds_names[4])
    variables_used<-.most_frequent_variables(model_list = mp$models,
                                             mp = mp,
                                             file = sprintf("%s_betas_summary.csv", save_name))

###################### 4. TIME FROM MRI, tots els relapses
cat('DOING MRI+DATA: 4 b. MANIA, TIME FROM MRI, tots els relapses \n')
response_var   = c('time_from_MRI', 'status')
covariates = c('age','Dona')
information_variables = c(covariates, predictor)
save_name <- paste(folder, 'mania_totes_timefromMRI_totsrelapses', sep="/")

    mp <- mripredict(paths, data_table_file, response_var, covariates, predictor, response_family, modulation, information_variables = information_variables)
    mp$models <- readRDS(rds_names[3])
    variables_used<-.most_frequent_variables(model_list = mp$models,
                                             mp = mp,
                                             file = sprintf("%s_betas_summary.csv", save_name))

###################### 5. TIME FROM MRI, relapse mania
cat('DOING MRI+DATA: 5 b. MANIA, TIME FROM MRI, relapse mania \n')
response_var   = c('time_from_MRI', 'status_relapse_mania')
covariates = c('age','Dona')
information_variables = c(covariates, predictor)
save_name <- paste(folder, 'mania_totes_timefromMRI_relapsemania', sep="/")

    mp <- mripredict(paths, data_table_file, response_var, covariates, predictor, response_family, modulation, information_variables = information_variables)
    mp$models <- readRDS(rds_names[2])
    variables_used<-.most_frequent_variables(model_list = mp$models,
                                             mp = mp,
                                             file = sprintf("%s_betas_summary.csv", save_name))

###################### 6. TIME FROM MRI, relapse depressiu
cat('DOING MRI+DATA: 6 b. MANIA, TIME FROM MRI, relapse depressio \n')
response_var   = c('time_from_MRI', 'status_relapse_depressio')
covariates = c('age','Dona')
information_variables = c(covariates, predictor)
save_name <- paste(folder, 'mania_totes_timefromMRI_relapsedepressio', sep="/")

    mp <- mripredict(paths, data_table_file, response_var, covariates, predictor, response_family, modulation, information_variables = information_variables)
    mp$models <- readRDS(rds_names[1])
    variables_used<-.most_frequent_variables(model_list = mp$models,
                                             mp = mp,
                                             file = sprintf("%s_betas_summary.csv", save_name))

#!/usr/bin/env Rscript

#::::::::::::::::::::::::::::: IMPORTS :::::::::::::::::::::::::::::::::
# totes les proves
list.of.packages <- c("sys", "methods", "readxl", "tcltk")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, library, character.only=TRUE)
library(sys)
args = commandArgs(trailingOnly=TRUE)
options(warn=-1)

source('mripredict_cv.R')
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


mri_paths_file = '~/Documents/dades/feps_totes/gm/subsamp/subsamp/subsamp/s4/mri_files.txt'
mri_fu_paths_file = '~/Documents/dades/feps_totes/gm_mod/subsamp/subsamp/subsamp/s4/mri_files.txt'
mri_wm_paths_file = '~/Documents/dades/feps_totes/wm/subsamp/subsamp/subsamp/s4/mri_files.txt'
mri_wm_fu_paths_file = '~/Documents/dades/feps_totes/wm_mod/subsamp/subsamp/subsamp/s4/mri_files.txt'
paths = cbind(read.table(mri_paths_file, stringsAsFactors = FALSE)$V1,
              read.table(mri_fu_paths_file, stringsAsFactors = FALSE)$V1,
              read.table(mri_wm_paths_file, stringsAsFactors = FALSE)$V1,
              read.table(mri_wm_fu_paths_file, stringsAsFactors = FALSE)$V1)
data_table_file = read.table("~/Documents/dades/feps_totes/data_feps_final.csv", sep=',', header = TRUE)
data_table_file$time_from_MRI[data_table_file$time_from_MRI==0] = data_table_file$time_from_MRI+0.5
data_table_file$time_from_remission[data_table_file$time_from_remission==0] = data_table_file$time_from_remission+0.5
########################################################
########## PROVES DE NOMÉS MRI #########################
########################################################
folder = 'output_feps/only_mri'
if(!dir.exists(folder)) dir.create(folder, recursive = TRUE)

###################### MRIPREDICT
predictor = c()
response_var   = c('time_from_remission', 'status')
covariates = c('age','Dona')
information_variables = c(covariates, predictor)
save_name <- paste(folder, 'feps_totes_timefromremission_totsrelapses', sep="/")
tmp_file <- 'tmp/feps1.txt'
if(!file.exists(tmp_file)){
  file.create(tmp_file)
  if (!file.exists(sprintf("%s_results_cv.csv",save_name))){
    mp <- mripredict(paths, data_table_file, response_var, covariates, predictor, response_family, modulation, information_variables = information_variables)
    
    mp <- mripredict_cv(mp, space = "NO_CHECK", save_name = save_name, folds_file = "", n_cores = 1,
                        use_significant_voxels = FALSE, use_ensemble_learning = TRUE, use_ensemble_voxels = TRUE, use_ensemble_subjects = FALSE, name_combat = 'feps1')
    
    write.csv(mp$cv_results,
              file=sprintf("%s_results_cv.csv",save_name),
              row.names = FALSE)
    variables_used<-.most_frequent_variables(model_list = mp$models,
                                             mp = mp,
                                             file = sprintf("%s_betas_summary.csv", save_name))
  }
  file.remove(tmp_file)
}
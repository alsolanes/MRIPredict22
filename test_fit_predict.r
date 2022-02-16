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
source('mripredict_fit.r')
source('mripredict_predict.r')
source('rotation3d.r')
source('cv_functions.r')
source('predict.r')
source('fit.r')
source('combat_quim.R')
source('combat_utils.R')
source('glmnet-utils.R')
#::::::::::::::::::::::::::: FI IMPORTS :::::::::::::::::::::::::::::::

#::::::::::::::::::::::::::: COMMON PARAMETERS ::::::::::::::::::::::::

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
data_table_file = read.table("~/Documents/dades/feps_totes/data_feps_final.csv", sep=',', header = TRUE, na.strings = c('','NA','na','N/A'))
data_table_file$time_from_MRI[data_table_file$time_from_MRI==0] = data_table_file$time_from_MRI+0.5
data_table_file$time_from_remission[data_table_file$time_from_remission==0] = data_table_file$time_from_remission+0.5
########################################################
########## PROVES DE NOMÉS MRI #########################
########################################################
folder = 'output_feps/only_mri'
if(!dir.exists(folder)) dir.create(folder, recursive = TRUE)

###################### 1. TIME FROM REMISSION, tots els relapses
cat('DOING MRI: 1. feps, TIME FROM REMISSION, tots els relapses \n')

predictors = c()
response_family = 'gaussian' 
modulation <- 'all'
response_var   = c('age')
covariates = c('Dona')
information_variables = c(covariates, predictor)
save_name <- paste(folder, 'feps_totes_timefromremission_totsrelapses', sep="/")
tmp_file <- 'tmp/feps1.txt'















######################################

#paths = cbind(as.character(mri_files$paths_gm_un), as.character(mri_files$paths_gm_fu), as.character(mri_files$paths_wm_un), as.character(mri_files$paths_wm_fu))
paths = paths

mp <- mripredict(mri_paths_file = paths, data_table_file = data_table_file, response_var = response_var, covariates = covariates, predictor = predictors, response_family = response_family, modulation = modulation, information_variables = information_variables)

mp <- mripredict_fit(mp = mp, space = "NO_CHECK", n_cores=1, ENSEMBLE_LEARNING = FALSE, name_combat = '')

mpres <- mripredict_predict(mp, mri_paths_file = paths, data_table_file = data_table_file, space='NO_CHECK')

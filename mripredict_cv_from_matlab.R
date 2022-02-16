#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
options(warn=-1)


path_r <- ifelse(dir.exists('R'),'R/','../R/')
source(paste(path_r,'mripredict_cv.R',sep=""))
source(paste(path_r,'mripredict_library.r',sep=""))
.require('sys')
source(paste(path_r,'mripredict.r',sep=""))
source(paste(path_r,'rotation3d.r',sep=""))
source(paste(path_r,'cv_functions.r',sep=""))
source(paste(path_r,'combat_quim.R',sep=""))
source(paste(path_r,'combat_utils.R',sep=""))
source(paste(path_r,'predict.r',sep=""))
source(paste(path_r,'fit.r',sep=""))
.require('glmnet')
source(paste(path_r,'glmnet-utils.R',sep=""))
.require("methods")
.require('readxl')
.require('tcltk')
.require("glmnet")
.require("oro.nifti")
.require("survival")
#.require("doParallel")
# .require("parallel")
.require("logistf")

# què necessitem per córrer mripredict_cv
# mripredict_cv = function(mp, space = "MNI", save_name = "results_cv", preloaded_covB_path = NULL, preloaded_covB_fu_path = NULL, folds_file = "", n_cores = 1,
#                   use_significant_voxels = FALSE, use_ensemble_learning = FALSE, use_ensemble_voxels = FALSE, use_ensemble_subjects = FALSE, name_combat = '')

# path: a fitxer on hi haurà columnes per cada modalitat seguint GM GM-MOD WM WM-MOD
# data_table_file: fitxer amb les dades clíniques
# response_var:
# covariates:
# predictor:
# response_family:
# modulation:
# information_variables:
##
# mp :: ve de sobre
# n_cores :: ho agafarà de la GUI (crear selector)
# resta de paràmetres per defecte

if(FALSE){
  
# TEST 1: 4 IMATGES AMB COVARS, BINOMIAL  
  mri_paths_file = 'C:/Users/alsol/Documents/MRIPredict/matlab/output/ff/mri_files.txt'
  data_table_file = 'C:/Users/alsol/Downloads/test10sujetos.csv'
  covars_file = 'C:/Users/alsol/Documents/MRIPredict/matlab/output/ff/covs.txt'
# TEST 2: 1 IMATGE AMB COVARS, BINOMIAL
  folder= 'C:/Users/alsol/Documents/MRIPredict/matlab/output/ff'
  

  covars = .read_covars_file(covars_file)
  response_var          <- covars[[1]]
  covariates            <- covars[[2]]
  predictor             <- covars[[3]]
  information_variables <- covars[[4]]
  if(length(covars)<4 || covars[[4]]=='') {
    information_variables = c(covariates, predictor)
  } else {
    information_variables = covars[[4]]
  }
  response_family = 'binomial' 
  modulation = 'un'
  #folder= '/home/aleix/Documents/papers_aleix/MRIPredict/MRIPredict/matlab/output'
  #folder = dirname(mri_paths_file)
  save_name = paste(folder, 'study',sep="/")
} else {
  # mp :: provinent de mripredict = function (mri_paths_file, data_table_file, response_var, covariates, predictor, response_family, modulation, information_variables = "")
  # mri_paths_file
  mri_paths_file        <- args[1] # mri files
  data_table_file       <- args[2] # clinical data
  covars_file           <- args[3] # file containing in line 1: outcome vars; 2: covars; 3: predictor vars; 4: informative vars.
  response_family       <- args[4] # 'binomial', 'gaussian', 'cox'
  modulation            <- args[5] # 'un', 'fu', 'op', 'all' (default mode with 4 images)
  folder                <- args[6] # folder where output will be saved
  save_name             <- args[7]
  if(is.na(save_name) || save_name == 'NA')
    save_name <-  paste(folder,'study',sep="/")
  else 
    save_name <- paste(folder, save_name, sep = "/")
  
  save_i = 1
  while(file.exists(sprintf("%s%s_results_cv_%s.csv",save_name,save_i,save_i))){
    save_i = save_i + 1   
  }
  save_name <- paste(save_name,save_i, sep="")
  print(save_name)
  # carregar covariables
  covars = .read_covars_file(covars_file)
  response_var          <- covars[[1]]
  covariates            <- covars[[2]]
  predictor             <- covars[[3]]
  information_variables <- covars[[4]]
  if(length(covars)<4 || covars[[4]]=='') {
    information_variables = c(covariates, predictor)
  } else {
    information_variables = covars[[4]]
  }
}


mp <- mripredict(mri_paths_file, data_table_file, response_var, covariates, predictor, response_family, modulation, information_variables = information_variables)

mp <- mripredict_cv(mp, space = "NO_CHECK", save_name = save_name, preloaded_covB_path = NULL, preloaded_covB_fu_path = NULL, folds_file = "", n_cores = 1,
                    use_significant_voxels = FALSE, use_ensemble_learning = FALSE, use_ensemble_voxels = FALSE, use_ensemble_subjects = FALSE, name_combat = '',
                    n_folds=10)
save_i = 1
while(file.exists(sprintf("%s_results_cv_%s.csv",save_name,save_i)))
  save_i = save_i + 1  
write.csv(mp$cv_results, file=sprintf("%s_results_cv%s.csv",save_name,ifelse(save_i>1,save_i,"")), row.names = FALSE)
write.csv(mp$cv_results, file=paste(folder,"results_cv.csv",sep="/"), row.names = FALSE)
variables_used<-.most_frequent_variables(model_list = mp$models, mp = mp, file=paste(folder,"betas_summary.csv",sep="/"))

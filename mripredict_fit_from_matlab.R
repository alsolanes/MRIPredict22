#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
options(warn=-1) # disable warnings

# NOTICE THAT MRIPREDICT_SAVE REQUIRES  "sudo apt-get install r-cran-xml" IN UBUNTU
#install.packages("XML", repos = "https://cran.r-project.org/")

source('../R/mripredict_fit.r')
source('../R/mripredict_library.r')
source('../R/mripredict.r')
source('../R/fit.r')
source('../R/imputation.R')
source('../R/cv_functions.r')
.require("methods")

space = "NO_CHECK"
# mri_file és llista de fitxers a llegir

DEBUG = TRUE

DEBUG_family = 'gaussian'
USE_PRELOADED_COVB = FALSE
TAKE_SIGNIFICANT = FALSE
ENSEMBLE_LEARNING = FALSE
if (!DEBUG) {
  
  mri_file = args[1]
  data_file = args[2]
  folds_file = args[3]
  response_family = args[4]
  covars_file = args[5]
  modulation = args[6]
  
} else {
  
  switch(DEBUG_family,
         binomial = {
           ## un i fu
           mri_file = "data/IXI/train/mri_files.txt" # subsampling
           data_file = "data/IXI/train/table_extra_field.txt"
           response_family = "binomial"
           #covars_file = "data_test/tmp_nomesPreds2.txt"
           covars_file = "data/IXI/covs.txt"
           #save_name = "tmp_test_binomial.xml"
           #folds_file = "data/IXI_folds.txt"
           modulation = "un"
         },
         gaussian = {
           ## un i fu
           mri_file = "data/IXI/train/mri_files.txt" # subsampling
           #mri_file=""
           data_file = "data/IXI/train/table_extra_field.txt"
           response_family = "gaussian"
           covars_file = "data/IXI/covs_gauss.txt"
           #covars_file = "data_test/tmp_covs.txt"
           save_name = "tmp_test_gaussian.xml"
           #folds_file = "data/IXI//IXI_folds.txt"
           modulation = "un"
         },
         cox = {
           # dis
           mri_file="/home/asolanesf/Documents/data/data_survival/mania_unitat_pol/rec/SUBSAMP/SUBSAMP/mri_files_mania_rec_subsubsamp.txt"
           data_file = "/home/asolanesf/Documents/data/data_survival/mania_unitat_pol/rec/clinical_rec.txt"
           response_family = "cox"
           covars_file = "/home/asolanesf/Documents/data/data_survival/mania_unitat_pol/rec/covs.txt"
           save_name = "tmp_cox_mania_unitat_rec.xml"
           folds_file = ""
           modulation = "un"
         },
         rec_sub = {
           # rec
           mri_file="/home/aleix/Documents/data/data_survival/mania_unitat_pol/rec/SUBSAMP/SUBSAMP/mri_files_mania_rec_subsubsamp.txt"
           data_file = "/home/aleix/Documents/data/data_survival/mania_unitat_pol/rec/clinical_rec_full.txt"
           response_family = "cox"
           covars_file = "/home/aleix/Documents/data/data_survival/mania_unitat_pol/rec/covs.txt"
           #covars_file = "/home/aleix/Documents/data/data_survival/mania_unitat_pol/rec/covs_full.txt" # falta crear
           save_name = sprintf("%stmp_cox_mania_unitat_rec_sub", folder)
           folds_file = ""
           modulation = "un"
         }
  )
}


mri_file = "/media/BACKUP_750/SCH_VS_CNT/gm/unmod/s2/subsamp/subsamp/mri_files.txt"
data_file ="/media/BACKUP_750/SCH_VS_CNT/clinical.txt"
#covars_file = "/media/BACKUP_750/SCH_VS_CNT/covs_clinical_age.txt"
covars_file = "/media/BACKUP_750/SCH_VS_CNT/covs_clinical_nocov.txt"
family_measure = "gaussian"
modulation = "un"
matter = "gm"


if (USE_PRELOADED_COVB) {
  
  if (SUBSAMP){
    preloaded_covB_path = c(sprintf('/media/BACKUP_750/effects/subsamp/covB_controls_%s_unmod_%s_beta1_subsamp.nii.gz', matter, smoothing),
                            sprintf('/media/BACKUP_750/effects/subsamp/covB_controls_%s_unmod_%s_beta2_subsamp.nii.gz', matter, smoothing),
                            sprintf('/media/BACKUP_750/effects/subsamp/covB_controls_%s_unmod_%s_beta3_subsamp.nii.gz', matter, smoothing))
    preloaded_covB_fu_path = c(sprintf('/media/BACKUP_750/effects/subsamp/covB_controls_%s_fumod_%s_beta1_subsamp.nii.gz', matter, smoothing),
                               sprintf('/media/BACKUP_750/effects/subsamp/covB_controls_%s_fumod_%s_beta2_subsamp.nii.gz', matter, smoothing),
                               sprintf('/media/BACKUP_750/effects/subsamp/covB_controls_%s_fumod_%s_beta3_subsamp.nii.gz', matter, smoothing))
  } else {
    preloaded_covB_path = c(sprintf('/media/BACKUP_750/effects/covB_controls_%s_unmod_%s_beta1_subsamp.nii.gz', matter, smoothing),
                            sprintf('/media/BACKUP_750/effects/covB_controls_%s_unmod_%s_beta2_subsamp.nii.gz', matter, smoothing),
                            sprintf('/media/BACKUP_750/effects/covB_controls_%s_unmod_%s_beta3_subsamp.nii.gz', matter, smoothing))
    preloaded_covB_fu_path = c(sprintf('/media/BACKUP_750/effects/covB_controls_%s_fumod_%s_beta1_subsamp.nii.gz', matter, smoothing),
                               sprintf('/media/BACKUP_750/effects/covB_controls_%s_fumod_%s_beta2_subsamp.nii.gz', matter, smoothing),
                               sprintf('/media/BACKUP_750/effects/covB_controls_%s_fumod_%s_beta3_subsamp.nii.gz', matter, smoothing))
  }
} else {
  preloaded_covB = ""
  preloaded_covB_path = ""
}
covars = .read_covars_file(covars_file)

response   = covars[[1]]
covariates = covars[[2]]
if(length(covars) < 3 || covars[[3]] == ''){
  predictors = c()
}else{
  predictors = covars[[3]]
}
if(length(covars)<4 || covars[[4]]==''){
  information_variables = c(covariates, predictors)
}else{
  information_variables = covars[[4]]
}
#mp = mripredict(mri_file, data_file, response, covariates, predictors, response_family, modulation)
mp = mripredict(mri_paths_file = mri_file, data_table_file = data_file, response_var = response,covariates =  covariates, predictor = predictors, response_family = response_family, modulation = modulation, information_variables = information_variables)

#fit
mp = mripredict_fit(mp = mp, space = space, preloaded_covB_path =  preloaded_covB_path, preloaded_covB_fu_path = preloaded_covB_fu_path, SIGNIFICANT = TAKE_SIGNIFICANT, n_cores = 1, ENSEMBLE_LEARNING = ENSEMBLE_LEARNING)
#print(save_name)
#mripredict_save(mp, save_name)
save(mp,file='mp_model.RData')
cat(sprintf('Model saved in mp_model.RData. '))

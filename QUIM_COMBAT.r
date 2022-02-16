# debug file for Quim - COMBAT
#!/usr/bin/env Rscript

# totes les proves
library(sys)
args = commandArgs(trailingOnly=TRUE)
options(warn=-1)

source('../R/mripredict_cv.R')
source('../R/mripredict_library.r')
source('../R/mripredict.r')
source('../R/imputation.R')
source('../R/rotation3d.r')
source('../R/cv_functions.r')
source('../R/predict.r')
source('../R/fit.r')
source('../R/combat_quim.R')
source('../R/combat_utils.R')
.require("methods")
.require('readxl')
.require('tcltk')

##########################
### DEFAULT PARAMETERS ###
##########################
space = "NO_CHECK" 

##########################
### PARAMETERS TO FILL ###
##########################
data_file = "/media/BACKUP1/SCH_VS_CNT/clinical_invented.txt"
n_cores = 1
# Where files will be saved
folder = '/home/aleix/Documents/MRIPredict/R/QUIM_COMBAT'
if (!file.exists(folder))
  system(sprintf("mkdir %s", folder))

family_measure = "binomial" # 'binomial', 'gaussian' or 'cox' 
folds_file = "/home/aleix/Documents/papers_aleix/MRIPredict/MRIPredict/R/results/sch/bin_rai_list_folds.txt"

########################## PROVA 7 ###############################################
########################## GM un S4 subsamp ######################################
cat('DOING: GM un S4 subsamp\n')
modulation = 'un'
smoothing = 's4'
matter = 'gm'
SIGNIFICANT = FALSE
ENSEMBLE_LEARNING = FALSE
USE_PRELOADED_COVB = FALSE
# carregar covariables
covars=c()
response   = 'grupos'
covariates = c('Sexe','Edat')
predictors = c()#c('TAP')
# information_variables are the variables that contain information (only interesting when there are NAs, used in multiple imputation)
# at least, must include covariables and predictors used in the model.
if(length(covars)<4 || covars[[4]]=='') {
  information_variables = c(covariates, predictors)
} else {
  information_variables = covars[[4]]
}
save_name = sprintf("%s/schVScnt_%s_%s_%s_%s_%s", folder, matter, modulation, smoothing, family_measure, ifelse(USE_PRELOADED_COVB,'usingPreloadedCov','noPreloadedCov'))
mri_file = "/media/BACKUP1/SCH_VS_CNT/gm/unmod/s4/subsamp/subsamp/mri_files.txt"
paths_fumod = "/media/BACKUP1/SCH_VS_CNT/gm/fumod/s4/subsamp/subsamp/mri_files.txt"
preloaded_covB_path = c(sprintf('/media/BACKUP1/effects/subsamp/subsamp/covB_controls_%s_unmod_%s_beta1_subsamp_subsamp.nii.gz', matter, smoothing), # Intercept
                        sprintf('/media/BACKUP1/effects/subsamp/subsamp/covB_controls_%s_unmod_%s_beta2_subsamp_subsamp.nii.gz', matter, smoothing), # Sexe
                        sprintf('/media/BACKUP1/effects/subsamp/subsamp/covB_controls_%s_unmod_%s_beta3_subsamp_subsamp.nii.gz', matter, smoothing)) # Edat
preloaded_covB_fu_path = c(sprintf('/media/BACKUP1/effects/subsamp/subsamp/covB_controls_%s_fumod_%s_beta1_subsamp_subsamp.nii.gz', matter, smoothing), # Intercept
                           sprintf('/media/BACKUP1/effects/subsamp/subsamp/covB_controls_%s_fumod_%s_beta2_subsamp_subsamp.nii.gz', matter, smoothing), # Sexe
                           sprintf('/media/BACKUP1/effects/subsamp/subsamp/covB_controls_%s_fumod_%s_beta3_subsamp_subsamp.nii.gz', matter, smoothing)) # Edat
paths = read.table(mri_file)
if (modulation != "un") {
  paths_fu = read.table(paths_fumod)
  paths = cbind(paths,paths_fu)
}
mp = mripredict(mri_paths_file = paths, 
                data_table_file = data_file, 
                response_var = response, 
                covariates = covariates, 
                predictor = predictors, 
                response_family = family_measure, 
                modulation = modulation, 
                information_variables = information_variables)
mp = mripredict_cv(mp=mp, 
                   space=space, 
                   save_name=save_name, 
                   n_cores=n_cores, 
                   preloaded_covB_path = ifelse(USE_PRELOADED_COVB,preloaded_covB_path,""), 
                   preloaded_covB_fu_path = ifelse(USE_PRELOADED_COVB,preloaded_covB_fu_path,""), 
                   folds_file=folds_file, 
                   SIGNIFICANT = SIGNIFICANT, 
                   ENSEMBLE_LEARNING = ENSEMBLE_LEARNING)

save_i = 1
while(file.exists(sprintf("%s_results_gm_%s_cv_%s_SUBSAMP.csv",save_name,smoothing,save_i)))
  save_i = save_i + 1  
write.csv(mp$cv_results, file=sprintf("%s_results_gm_%s_cv_%s_SUBSAMP.csv",save_name,smoothing,save_i), row.names = FALSE)

########################## FI: GM un S4 subsamp ##################################
##################################################################################


options(warn=1) # reenable warnings
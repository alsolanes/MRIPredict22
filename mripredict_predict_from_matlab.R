
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
options(warn=-1) # disable warnings

#source('../R/mripredict_fit.r')
source('../R/mripredict_library.r')
source('../R/mripredict.r')
source('../R/mripredict_predict.r')
source('../R/mripredict_load.r')
source('../R/predict.r')
.require("methods")

space = "NO_CHECK"
# mri_file és llista de fitxers a llegir
#mri_file = "/home/asolanesf/Documents/data/toy_ds/forR/optmod_subsamp_files.txt"
xml_path = args[1]
#data file ha d'estar separat per tabs (copiant a excel i passant a txt per exemple)
mri_file = args[2]
data_file = args[3]
results_file_name = args[4]


# xml_path = "../R/current_model.xml"
# mri_file = "../matlab/R_tmp/tmp_mri_files.txt"
# #data_file = "/home/asolanesf/Documents/MRIPredict/matlab/data/cv/oasis_clinic.txt"
# data_file = "/home/asolanesf/Documents/MRIPredict/matlab/data/predict/oasis_predict.txt"
# results_file_name = "output_file_results.txt"

# mri_file = "data/IXI/test/tmp_mri_files.txt" # subsampling
# #mri_file = "/home/asolanesf/Documents/data/optmod/list_unfu.txt" # dades no subsamplejades
# data_file = "data/IXI/test/table_extra_field.txt"
# response_family = "binomial"
# covars_file = "data/IXI/test/tmp_covs.txt"
# save_name = "tmp_test.xml"
# folds_file = "data/IXI/test/IXI_folds.txt"
# modulation = "op"

#mp = mripredict_load(xml_path)
# load model
load('mp_model.RData')

DEBUG = TRUE
DEBUG_family = 'gaussian'

if (!DEBUG) {
  
  xml_path = args[1]
  mri_file = args[2]
  data_file = args[3]
  results_file_name = args[4]
  
} else {
  
  switch(DEBUG_family,
         binomial = {
           ## un i fu
           mri_file="data/IXI/test/tmp_mri_files.txt"
           data_file = "data/IXI/test/table_extra_field.txt"
           response_family = "binomial"
           #covars_file = "data/IXI/test/tmp_nomesPreds2.txt"
           covars_file = "data/IXI/covs.txt"
           #save_name = sprintf("%stmp_test_binomial.xml",folder)
           #xml_path = "tmp_test_binomial.xml"
           #folds_file = "data/IXI/test/IXI_folds.txt"
           modulation = "un"
         },
         gaussian = {
           ## un i fu
           mri_file = "data/IXI/test/tmp_mri_files.txt" # subsampling

           data_file = "data/IXI/test/table_extra_field.txt"
           response_family = "gaussian"
           covars_file = "data/IXI/covs_gauss.txt"
           covars_file = ""
           #covars_file = "data/IXI/test/tmp_covs.txt"
           xml_path = "tmp_test_gaussian.xml"
           #folds_file = "data/IXI/test/IXI_folds.txt"
           modulation = "un"
           
           
           
           mri_file = "/media/BACKUP_750/SCH_VS_CNT/gm/unmod/s2/subsamp/subsamp/mri_files.txt"
           data_file = "/media/BACKUP_750/SCH_VS_CNT/clinical.txt"
           data_file = ""
           covars_file = "/media/BACKUP_750/SCH_VS_CNT/covs_clinical_age.txt"
           family_measure = "gaussian"
           modulation = "un"
           matter = "gm"
           
           
         },
         cox = {
           # 
           xml_path = "/home/asolanesf/Documents/MRIPredict/R/tmp_cox_mania_unitat_rec.xml"
           
           # rec
           mri_file="/home/asolanesf/Documents/data/data_survival/mania_unitat_pol/rec/SUBSAMP/SUBSAMP/mri_files_mania_rec_subsubsamp.txt"

           data_file = "/home/asolanesf/Documents/data/data_survival/mania_unitat_pol/rec/clinical_rec.txt"
           response_family = "cox"
           covars_file = "/home/asolanesf/Documents/data/data_survival/mania_unitat_pol/rec/covs.txt"
           save_name = "tmp_cox_mania_unitat_rec.xml"
           folds_file = ""
           modulation = "un"
         }
  )
  
}
results_file_name = "/home/aleix/Documents/MRIPredict_2019/MRIPredict/matlab/output/test_fit/predicted_results.txt"


results = mripredict_predict(mp = mp, mri_paths_file = mri_file, data_table_file = data_file, space = space, n_cores=1)
cat('Predictions:\n')
print(rowMeans(results))
write.csv(rowMeans(results),file = 'results.csv')
cat(sprintf('Results saved in: results.csv'))

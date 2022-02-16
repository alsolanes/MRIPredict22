Fit a model:
mp = mripredict(mri_paths_file, data_table_file, response_var, covariates, predictor, response_family, modulation, information_variables = "")
mripredict_fit(mp, space = "MNI", save_name = "results_cv", preloaded_covB_path = NULL, preloaded_covB_fu_path = NULL, folds_file = "", n_cores = 1,
                          use_significant_voxels = FALSE, use_ensemble_learning = FALSE, use_ensemble_voxels = FALSE, use_ensemble_subjects = FALSE,
                          ide_shiny=FALSE, standardize_images=FALSE) 

Use a previously fitted model on new data:
mp = load('previous_model.rds')
mripredict_predict(mp, mri_paths_file = new_files, mri_fu_paths_file = new_fully_modulated_files, data_table_file=NULL, space, n_cores=1)

Perform a cross-validation:
mp = mripredict(mri_paths_file, data_table_file, response_var, covariates, predictor, response_family, modulation, information_variables = "")
mripredict_cv(mp, space = "MNI", save_name = "results_cv", preloaded_covB_path = NULL, preloaded_covB_fu_path = NULL, folds_file = NULL, n_cores = 1,
                         use_significant_voxels = FALSE, use_ensemble_learning = FALSE, use_ensemble_voxels = FALSE, use_ensemble_subjects = FALSE, n_folds = 10,
                         ide_shiny=FALSE, standardize_images=FALSE)

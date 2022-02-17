mripredict_predict = function(mp, mri_paths_file = "", mri_fu_paths_file = "", data_table_file=NULL, space, n_cores=1) {
  .require("glmnet")
  .require("oro.nifti")
  .require("survival")
  .require("doParallel")
  .require("parallel")
  .require("logistf")
  SIGNIFICANCE_THRESHOLD = qnorm(0.975)
  
  N_ITERATIONS = mp$n_iterations
  N_M_IMPUTATIONS = mp$n_imputations
  
  #####################################################################################
  ### parallel definition
  if(n_cores=='auto') {
    n_cores = 2
  }
  #####################################################################################
  
  .print_action("Setting new MRIPredict model")
  mri_paths = ""
  mri_fu_paths = ""
  if (length(mri_paths_file)<=2){
    if (mri_paths_file != '' && dim(mri_paths_file)[1]!=0) {
      if(file.exists(mri_paths_file)) {
        mri_paths = .read_1col_file(mri_paths_file)
        mri_fu_paths = .read_2col_file(mri_paths_file)
        if(file.exists(mri_fu_paths_file)){
          mri_fu_paths = .read_1col_file(mri_fu_paths_file)
        }
      } else {
        stop("MRI paths file does not exist.")
      }
    }
  } else{
    mri_paths = mri_paths_file
  }
  if (data_table_file!="" && dim(data_table_file)[1]!=0){
    if (length(data_table_file)==1){
      data_table = .read_data_table(data_table_file)
    } else {
      data_table = data_table_file
    }
  }
  else {
    data_table = ""
  }
  #####################################################################################
  
  ## load MRI data, covariates and folds
  .print_action("Loading MRI data")
  a = Sys.time()
  mri = NULL
  if (mri_paths != "")
    mri    = .load_mri(mri_paths, space = space)
  # load fu data
  mri_fu = NULL
  if (!mp$modulation=='un' && mri_fu_paths != "") {
    mri_fu = .load_mri(mri_fu_paths, space=space)
  }
  mri_wm = NULL
  mri_wm_fu = NULL
  if (!is.null(mp$mri_wmmod_paths) && mp$mri_wmmod_paths != "") {
    mri_fu = .load_mri(mp$mri_fu_paths, mp$mask_fu$data, space = space)
    mri_wm = .load_mri(mp$mri_wm_paths, mp$mask_wm$data, space = space)
    mri_wm_fu = .load_mri(mp$mri_wmmod_paths, mp$mask_wmmod$data, space = space)
  }
  
  .print_action("Creating response vector and covariate matrix")
  ## define covariates matrix
  # if(!is.null(mp$covX_transf)){
  #   covX = .create_covX(data_table, mp$covX_transf, mri$n) # mri$n is the number of subjects, which is not present in data_table if there are no covariates
  # }  else{
  #   covX = matrix(1,nrow(data_table_file))
  # }
  
  n_subjects = nrow(data_table)
  ## define covariates matrix
  if(!is.null(mp$covX_transf))
    covX = .create_covX(data_table, mp$covX_transf)
  #data_informative_table = .create_covX(mp$data_table, mp$data_table_transf)

  sites = NULL
  if(!is.null(mp$data_table_transf)){
    
    
    data_informative_table = data.frame2glmnet.matrix(mp$data_table_transf, data_table)
    if ('site' %in% colnames(data_table)) {
      sites = factor(x = data_table[, "site"])
    } 
  }
  else{
    data_informative_table = NULL
    data_table_imputed_test = NULL
  }
  
  
  
  .print_ok()
  
  ####################################################################################
  linPreds = c()
  i = 1
  for (iter in 1:N_ITERATIONS){
    for (iter_imputation in 1:N_M_IMPUTATIONS) {
      model_list = mp$models[[i]] # load models stored
      i = i+1
      lasso_ijk = model_list$lasso_ijk
      signif_indx = model_list$signif_indx
      lasso_covB = model_list$lasso_covB
      lasso_covB_fu = model_list$lasso_fu_covB
      mask = model_list$mask
      scale_clinical = model_list$predictors_scale
      
      ###################################### MULTIPLE IMPUTATION #################################
      if(!is.null(mp$impute_obj)){
        imp.data_informative_table_test = impute.glmnet.matrix(m = mp$impute_obj, x = data_informative_table, nimp = N_M_IMPUTATIONS)
      }
      
      ###########################################################################################
      if(!is.null(data_table_imputed_test)){
        if(!is.null(mp$covX_var)){
          covX_test = cbind(1,data_table_imputed_test[, .indx_equivalent_columns_from_transf(mp$data_table_transf, mp$covX_var)])
        } else {
          covX_test = data_table_imputed_test
        }
      } else {
        if(is.null(data_table_file)){
          covX_test = matrix(1,nrow=mri$n)
        } else {
          covX_test = matrix(1, nrow=nrow(data_table))
        }
      }
      if(!is.null(mp$pred_var)){
        predX_test = matrix(data_table_imputed_test[, mp$pred_var],nrow=nrow(data_table_imputed_test))
      }
      tipett_take_un = model_list$tipett_take_un

      
      ### TESTING ##
      preds = apply_model(mp = mp, 
                          mri_data = mri$data, 
                          mri_fu_data = mri_fu$data, 
                          mri_wm_data = mri_wm$data,
                          mri_wm_fu_data = mri_wm_fu$data,
                          covX_test = covX_test, 
                          signif_indx = model_list$signif_indx,
                          lasso_covB = model_list$lasso_covB,
                          lasso_covB_fu = model_list$lasso_covB_fu, 
                          mask = model_list$mask,
                          masks_3d = model_list$masks_3d,
                          predX_test = predX_test,
                          scale_clinical = model_list$scale_predictors,
                          scale_mri = model_list$scale_mri,
                          lasso = model_list$lasso,
                          lasso_predX_indx = model_list$lasso_predX_indx,
                          covX_site = sites,
                          n_voxels_mask = model_list$n_voxels_mask,
                          combat = model_list$combat
      )
      linPreds = cbind(linPreds,preds)
    }
  }
  linPreds
}

if(F) {
  mp = readRDS('results_cv_mp.rds')
  a = mripredict_predict(mp, mri_paths_file = mp$mri_paths, mri_fu_paths_file = "", data_table_file = mp$data_table, space = "NO_CHECK", n_cores=1)
}

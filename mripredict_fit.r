mripredict_fit = function(mp, space = "MNI", save_name = "results_cv", preloaded_covB_path = NULL, preloaded_covB_fu_path = NULL, folds_file = "", n_cores = 1,
                          use_significant_voxels = FALSE, use_ensemble_learning = FALSE, use_ensemble_voxels = FALSE, use_ensemble_subjects = FALSE,
                          ide_shiny=FALSE, standardize_images=FALSE) {
  # space refers to the space that the data should be
  .require("glmnet")
  .require("oro.nifti")
  .require("survival")
  .require("doParallel")
  # .require("parallel")
  .require("logistf")
  data_informative_table = NULL
  imp.data_informative_table_train = NULL
  imp.data_informative_table_test= NULL
  #.require("doSNOW") # dubtós 
  SIGNIFICANCE_THRESHOLD = qnorm(0.975)
  #n_folds = 10
  if (use_ensemble_learning) {
    N_ITERATIONS = 18
  } else {
    N_ITERATIONS = 1
  }
  
  n_multiple_imputations = 1
  N_M_IMPUTATIONS = 20 # multiple imputations
  #####################################################################################
  ### parallel definition
  if (n_cores == 'auto') {
    n_cores = max(round(detectCores() / 2), 1)
  }
  #####################################################################################
  # load masks
  
  mp$mask$data <- mp$mask_fu$data <- mp$mask_wm$data <- mp$mask_wmmmod$data <- NULL
  if (!is.null(mp$mask_path) && mp$mask_path!="") mp$mask = .load_mri(mp$mask_path, space = 'NO_CHECK')
  if (!is.null(mp$mask_fu_path) && mp$mask_fu_path!="") mp$mask_fu = .load_mri(mp$mask_fu_path, space = 'NO_CHECK')
  if (!is.null(mp$mask_wm_path) && mp$mask_wm_path!="") mp$mask_wm = .load_mri(mp$mask_wm_path, space = 'NO_CHECK')
  if (!is.null(mp$mask_wmmod_path) && mp$mask_wmmod_path!="") mp$mask_wmmod = .load_mri(mp$mask_wmmod_path, space = 'NO_CHECK')
  ## load MRI data, covariates and folds
  .print_action("Loading MRI data")
  a = Sys.time()
  mri = NULL
  
  if (!is.null(mp$mri_paths) && mp$mri_paths != ""){
    mri = .load_mri(mp$mri_paths, mask=mp$mask$data, space = space)
  }
  
  # load fully modulated data
  mri_fu = NULL
  if (!mp$modulation == 'un' && !is.null(mp$mri_fu_paths)) {
    mri_fu = .load_mri(mp$mri_fu_paths, mp$mask_fu$data, space = space)
    
  }
  mri_wm = NULL
  mri_wm_fu = NULL
  if (!is.null(mp$mri_wmmod_paths) && mp$mri_wmmod_paths != "") {
    mri_fu = .load_mri(mp$mri_fu_paths, mp$mask_fu$data, space = space)
    mri_wm = .load_mri(mp$mri_wm_paths, mp$mask_wm$data, space = space)
    mri_wm_fu = .load_mri(mp$mri_wmmod_paths, mp$mask_wmmod$data, space = space)
  }
  cat("Time loading data:", difftime(Sys.time(), a, units = "mins"), "mins.\n")
  .print_ok()
  
  
  .print_action("Creating response vector and covariate matrix")
  ## define Y depending if cox or gaussian/binomial
  
  Y = c()
  if (mp$response_family == "cox") {
    Y = mp$data_table[, match(mp$response_var, colnames(mp$data_table))]
    Y = Surv(time = as.numeric(Y[, 1]), event = as.numeric(Y[, 2]))
  } else {
    Y = .create_Y(mp$data_table, mp$response_var, mp$response_event)
  }
  
  
  n_subjects = nrow(mp$data_table)
  ## define covariates matrix
  if(!is.null(mp$covX_transf))
    covX = .create_covX(mp$data_table, mp$covX_transf)
  #data_informative_table = .create_covX(mp$data_table, mp$data_table_transf)
  if(!is.null(mp$data_table_transf))
    data_informative_table = data.frame2glmnet.matrix(mp$data_table_transf, mp$data_table)
  
  # as.matrix(data_informative_table[, .indx_equivalent_columns_from_transf(mp$data_table_transf, mp$pred_var)])
  #mp$pred_var_dummies = colnames(data_informative_table[, .indx_equivalent_columns_from_transf(mp$data_table_transf, mp$pred_var)])
  .print_ok()
  ## get list of folds or create it
  
  if ('site' %in% colnames(mp$data_table)) {
    sites = factor(x = mp$data_table[, "site"])
  } else {
    sites = NULL
  }
  
  
  
  
  n_subjects = min(length(Y), nrow(Y)) # -number of subjects used
  subjects_used = matrix(nrow = n_subjects, ncol = N_ITERATIONS)
  subjects_used = list()
  z_cuts = c()
  
  
  train_thresholds_all = c()
  
  # if (file.exists(sprintf("%s_train_thresholds.csv", save_name))) {
  #   train_thresholds_all = read.csv(sprintf("%s_train_thresholds.csv", save_name))
  # }
  ##################################################################################
  ##################################################################################
  # N_ITERATIONS: one iteration splitting the brain in different angles and axes.
  mp$models = list()
  model_counter = 1
  time_points = c(30, 60, 180, 360, 720)
  
  
  
  
  # Define train and test
  
  
  
  
  # multiple imputation
  #data_table_imputed_train = matrix(data_informative_table[training,], ncol=ncol(data_informative_table))
  #data_table_imputed_test = matrix(data_informative_table[test,], ncol=ncol(data_informative_table))
  if(!is.null(data_informative_table)) {
    data_table_imputed_train = data_informative_table
    
    if ( any(is.na(data_informative_table))) {
      n_multiple_imputations = N_M_IMPUTATIONS #20
      impute_obj = impute.glmnet.matrix_fit(x = data_table_imputed_train, n_cores = 4)
      imp.data_informative_table_train = impute.glmnet.matrix(m = impute_obj, x = data_table_imputed_train, nimp = N_M_IMPUTATIONS)
      mp$impute_obj = impute_obj
      
    } else {
      n_multiple_imputations = 1
      mp$impute_obj = NULL
    }
  }
  else{
    data_table_imputed_train = NULL
    imp.data_informative_table_train = NULL
  }
  start_fold.time <- Sys.time()
  
  mp$n_iterations = N_ITERATIONS
  mp$n_imputations = n_multiple_imputations
  # for (i_time_point in 1:length(time_points)) {
  #   train_linear_predictor_to_find_the_threshold_all[[i_time_point]] = c()
  # }
  
  #####################################################################################
  #####################################################################################
  ### CROSSVALIDATION ###
  for (iter in 1:N_ITERATIONS) {
    #sprintf("%s_iteration%d_%d.txt", save_name, iter, fold, iter)
    
    if (N_ITERATIONS > 1 && use_ensemble_subjects) {
      subjects_used[[iter]] = training[sample(1:length(training), replace = TRUE)] ### a canviar!
      training = subjects_used[[iter]]
    }
    if (mp$response_family == "cox") {
      trainY = Y
    } else {
      trainY = Y
    }
    
    linPreds = c()
    if (mp$response_family == "cox") {
      train_thresholds = NULL
    }
    start.time <- Sys.time()
    
    if (N_ITERATIONS > 1)
      cat("Iteration ", iter, "of", N_ITERATIONS, "\n", save_name)
    
    sprintf('/nIteration: %d', iter)
    ##########################################################################
    # IF MULTIPLE IMPUTATION. PRED = MEAN ( PRED_PER_FOLD )
    ##########################################################################
    preds = c()
    
    for (iter_imputation in 1:n_multiple_imputations) {    
      if (n_multiple_imputations > 1 & !is.null(imp.data_informative_table_train)) {
        data_table_imputed_train = imp.data_informative_table_train[[iter_imputation]]
        #data_table_imputed_test = as.matrix(imp.data_informative_table_test[[iter_imputation]])
      }
      if (n_multiple_imputations > 1)
        .print_action(paste("Imputation", iter_imputation, "of", n_multiple_imputations, "\n", save_name))
      
      # CALCULATE MODEL
      internal_folds = c()
      
      internal_folds = assign.folds(y = trainY, family = mp$response_family, nfolds = 10, site = sites)
      
      
      
      model_list = fit_model(mp = mp,
                             data_informative_table = data_table_imputed_train,
                             Y = trainY,
                             mri = mri$data,
                             mri_fu = mri_fu$data,
                             mri_wm = mri_wm$data,
                             mri_fu_wm = mri_wm_fu$data,
                             iter = ifelse(use_ensemble_voxels == FALSE, -1, iter),
                             internal_folds = internal_folds,
                             n_cores = n_cores,
                             use_significant_voxels = use_significant_voxels,
                             covX_site = sites,
                             #name_combat = name_combat,
                             standardize_images = standardize_images
      )
      #mp$combat = model_list$combat
      if (!mp$modulation %in% c('cl','clinical')) {
        # per passar a mni la ijk =>  mp$mni = (mri$sto.xyz %*% mp$ijk)[1:3, ]          
        lasso_ijk = matrix(model_list$lasso_ijk, nrow = 3)
        # lasso_ijk està amb totes les imatges juntes
        ##############################################################
        dim_x = dim(mri$data)[1]
        if (mp$modulation == 'fu') {
          idx_gm_mod = which(lasso_ijk[1,] > dim_x * 0 & lasso_ijk[1,] <= dim_x * 1)
          idx_gm = which(lasso_ijk[1,] > dim_x * 1 & lasso_ijk[1,] <= dim_x * 2)
          
          ijk_gm     = lasso_ijk[, idx_gm] - c(dim_x * 0, 0, 0)
          model_list$lasso_mni$gm     = ijk2mni(rbind(matrix(ijk_gm, nrow = 3), 1), mri$sto.ijk)[, 1:length(idx_gm)]
          model_list$lasso_mni$gm_betas = model_list$lasso$beta[idx_gm]
          
          ijk_gm_mod = lasso_ijk[, idx_gm_mod] - c(dim_x * 1, 0, 0)
          model_list$lasso_mni$gm_mod = ijk2mni(rbind(matrix(ijk_gm_mod, nrow = 3), 1), mri$sto.ijk)[, 1:length(idx_gm_mod)]
          model_list$lasso_mni$gm_mod_betas = model_list$lasso$beta[idx_gm_mod]
        } else {
          idx_gm = which(lasso_ijk[1,] > dim_x * 0 & lasso_ijk[1,] <= dim_x * 1)
          idx_gm_mod = which(lasso_ijk[1,] > dim_x * 1 & lasso_ijk[1,] <= dim_x * 2)
          
          
          ijk_gm     = lasso_ijk[, idx_gm] - c(dim_x * 0, 0, 0)
          model_list$lasso_mni$gm     = ijk2mni(rbind(matrix(ijk_gm, nrow = 3), 1), mri$sto.ijk)[, 1:length(idx_gm)]
          model_list$lasso_mni$gm_betas = model_list$lasso$beta[idx_gm]
          if(length(idx_gm_mod)>0){
            ijk_gm_mod = lasso_ijk[, idx_gm_mod] - c(dim_x * 1, 0, 0)
            model_list$lasso_mni$gm_mod = ijk2mni(rbind(matrix(ijk_gm_mod, nrow = 3), 1), mri$sto.ijk)[, 1:length(idx_gm_mod)]
            model_list$lasso_mni$gm_mod_betas = model_list$lasso$beta[idx_gm_mod]
          }
        }
        idx_wm = which(lasso_ijk[1,] > dim_x * 2 & lasso_ijk[1,] <= dim_x * 3)
        if(length(idx_wm)>0){
          ijk_wm     = lasso_ijk[, idx_wm] - c(dim_x * 2, 0, 0)
          model_list$lasso_mni$wm     = ijk2mni(rbind(matrix(ijk_wm, nrow=3),1), mri$sto.ijk)[, 1:length(idx_wm)]
          model_list$lasso_mni$wm_betas = model_list$lasso$beta[idx_wm]
        }
        idx_wm_mod = which(lasso_ijk[1,] > dim_x * 3 & lasso_ijk[1,] <= dim_x * 4)
        if(length(idx_wm_mod)>0){
          ijk_wm_mod = lasso_ijk[, idx_wm_mod] - c(dim_x * 3, 0, 0)
          model_list$lasso_mni$wm_mod =  ijk2mni(rbind(matrix(ijk_wm_mod, nrow=3),1), mri$sto.ijk)[, 1:length(idx_wm_mod)]
          model_list$lasso_mni$wm_mod_betas = model_list$lasso$beta[idx_wm_mod]
        }
        ##############################################################
        
        signif_indx = model_list$signif_indx 
        lasso_covB = model_list$lasso_covB
        lasso_covB_fu = model_list$lasso_fu_covB
        
        n_voxels_mask = model_list$n_voxels_mask
        
        #################################################################################################################################################
      } else {
        n_voxels_mask = 0
      }
      
      tipett_take_un = model_list$tipett_take_un
      # clinical vars to scale
      scale_clinical = model_list$scale_predictors
      scale_mri = model_list$scale_mri
      ###########################################################################################################################################################
      # END CALCULATING MODEL
      
      
      
      
      mp$models[[model_counter]] <- model_list 
      model_counter <- model_counter + 1
      ###########################################################################################################################################################
      .print_ok()
      ### END TRAINING ###
      
      if(ide_shiny)
        incProgress(1/(n_folds*N_ITERATIONS*n_multiple_imputations), detail = paste("Doing crossvalidation. Fold",fold, ". Iteration:", iter, ". Imputation:",iter_imputation))
    }
    # fi for n_multiple_imputations  
    
    
    
    # average of linear predictors in multiple_imputation
    
    
    .print_ok()
    
    ##########################################################################
    # FI FOR MULTIPLE IMPUTATION. PRED = MEAN ( PRED_PER_FOLD )
    ##########################################################################
    
    
    
    end.time <- Sys.time()
    
    
    #######################
    if (iter > 1) message("Time per iteration:", difftime(end.time, start.time, units = 'mins'), ' mins.')
    
    
    # fi ITERATIONS (BRAIN CUTS)
    
    
    
  }
  
  
  # save rds
  # cat('\n[Saving] - Saving fold model to',sprintf("%s_model", save_name))
  # name_base <- sprintf("%s_model", save_name)
  # saveRDS(object = mp, file = sprintf("%s_mp.rds", name_base))
  #saveRDS(object = mp$models, file = sprintf("%s_model_list.rds", name_base))
  
  
  
  
  
  # fi cv  
  
  cat('\n[Saving] - Saving MRIPredict object to:', sprintf('%s_mp.rds',save_name))
  
  saveRDS(mp,file=sprintf('%s_mp.rds',save_name))
  mp
  
}
if(F){
  load('debug.RData')
  mripredict_fit(mp, space = "NO_CHECK", save_name = "results_cv", preloaded_covB_path = NULL, preloaded_covB_fu_path = NULL, folds_file = "", n_cores = 1,
                 use_significant_voxels = FALSE, use_ensemble_learning = T, use_ensemble_voxels = T, use_ensemble_subjects = F,
                 ide_shiny=FALSE, standardize_images=FALSE)
}

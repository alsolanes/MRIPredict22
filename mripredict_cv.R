mripredict_cv = function(mp, space = "MNI", save_name = "results_cv", preloaded_covB_path = NULL, preloaded_covB_fu_path = NULL, folds_file = NULL, n_cores = 1,
                         use_significant_voxels = FALSE, use_ensemble_learning = FALSE, use_ensemble_voxels = FALSE, use_ensemble_subjects = FALSE, n_folds = 10,
                         ide_shiny=FALSE, standardize_images=FALSE) {
  # space refers to the space that the data should be
  
  .require("glmnet")
  .require("oro.nifti")
  .require("survival")
  .require("doParallel")
  # .require("parallel")
  .require("logistf")
  #.require("doSNOW") # dubtós 
  SIGNIFICANCE_THRESHOLD = qnorm(0.975)
  #n_folds = 10
  if (use_ensemble_learning) {
    N_ITERATIONS = 18
  } else {
    N_ITERATIONS = 1
  }
  N_M_IMPUTATIONS = 20 # multiple imputations
  # if save_name does not contain any folder, let's save it in output folder
  if(!grepl("/", save_name, fixed = TRUE)){
    if(!dir.exists(paste0("output/",save_name)))
      dir.create(paste0('output/', save_name),recursive = T)
    save_name = paste0('output/', save_name, "/",save_name)
  }
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
  cat("Time loading data:", round(difftime(Sys.time(), a, units = "mins"),2), "mins.\n")
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
  
  if (mp$response_family == "binomial" && min(table(Y))<n_folds) {
    message('Warning: Too few samples in one outcome group. The program may be unable to perform the cross-validation.')
  }  
  n_subjects = nrow(mp$data_table)
  ## define covariates matrix
  
  covX = .create_covX(mp$data_table, mp$covX_transf)
  
  if(!is.null(mp$data_table_transf))
    data_informative_table = data.frame2glmnet.matrix(mp$data_table_transf, mp$data_table)
  else
    data_informative_table = NULL
  # as.matrix(data_informative_table[, .indx_equivalent_columns_from_transf(mp$data_table_transf, mp$pred_var)])
  #mp$pred_var_dummies = colnames(data_informative_table[, .indx_equivalent_columns_from_transf(mp$data_table_transf, mp$pred_var)])
  .print_ok()
  ## get list of folds or create it
  
  if ('site' %in% colnames(mp$data_table)) {
    sites = factor(x = mp$data_table[, "site"])
  } else {
    sites = NULL
  }
  
  if (is.null(folds_file))
    folds_file = sprintf("%s_list_folds_n%s.txt", save_name,n_folds)
  
  if (file.exists(folds_file)) {
    # if folds file has been previously created
    list_folds = .read_folds_file(folds_file)
  } else {
    # if empty folds file name
    list_folds = NULL
  }
  
  if (is.null(list_folds)) {
    message('No folds file preloaded. New folds distribution will be created.')
    assigned_fold = assign.folds(y = Y, family = mp$response_family, nfolds = n_folds, site = sites)
    ### CREATE FILE 
    write(assigned_fold, file = sprintf("%s_list_folds_n%s.txt", save_name,n_folds), sep = ",", ncolumns = length(assigned_fold))
    is_folds_loaded = FALSE
  } else if (n_folds != length(unique(list_folds[[1]]))) {
    message('The folds file specified has different number of folds. New folds distribution will be created.')
    assigned_fold = assign.folds(y = Y, family = mp$response_family, nfolds = n_folds, site = sites)
    ### CREATE FILE 
    write(assigned_fold, file = sprintf("%s_list_folds_n%s.txt", save_name,n_folds), sep = ",", ncolumns = length(assigned_fold))
    is_folds_loaded = FALSE
  } else {
    assigned_fold = list_folds[[1]]
    if (length(list_folds) < (n_folds + 1)) {
      #write(assigned_fold, file=sprintf("%s_list_folds.txt",save_name), sep=",", ncolumns = length(assigned_fold))
      is_folds_loaded = FALSE
    } else {
      is_folds_loaded = TRUE
    }
    if (length(assigned_fold) != n_subjects) {
      stop('The number of subjects in the fold file, and the subjects selected are not the same. Please check that you selected the correct fold file, or the correct subjects.')
    }
  }
  if (!is.null(preloaded_covB_path) && preloaded_covB_path != "") {
    preloaded_covB = .load_mri(preloaded_covB_path, space = 'NO_CHECK')
    if (!mp$modulation == 'un' && mp$mri_fu_paths != "") {
      preloaded_covB_fu = .load_mri(preloaded_covB_fu_path, space = 'NO_CHECK')
    }
    if (!(all(dim(preloaded_covB$data[, , , 1]) == dim(mri$data[, , , 1])))) {
      stop('Preloaded effects do not have the same dimensions that the images.')
    }
  } else {
    preloaded_covB = ""
  }
  
  cv_table = matrix(nrow = nrow(mp$data_table), ncol = 2)
  cv_table[, 1] = assigned_fold
  
  cv_accuracy = c()
  cv_betas = c()
  cv_table_predictions = data.frame()
  
  results = c()
  
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
  for (fold in 1:n_folds) {
    
    if (!file.exists(sprintf('%s_fold%s.Rdata', save_name, fold))){
      
      cat(paste("Starting fold", fold, "of", n_folds),"...\n")
      
      
      
      # Define train and test
      training = which(assigned_fold != fold)
      test = which(assigned_fold == fold)
      if (N_ITERATIONS > 1 && use_ensemble_subjects) {
        subjects_used[[iter]] = training[sample(1:length(training), replace = TRUE)] ### a canviar!
        training = subjects_used[[iter]]
      }
      if (mp$response_family == "cox") {
        trainY = Y[training,]
        testY = Y[test,]
      } else {
        trainY = Y[training]
        testY = Y[test]
      }
      sites_training = sites[training]
      sites_test = sites[test]
      # CALCULATE MODEL
      internal_folds = c()
      
      if (is_folds_loaded) {            
        internal_folds = list_folds[[fold + 1]]
      } else {
        internal_folds = assign.folds(y = trainY, family = mp$response_family, nfolds = n_folds, site = sites_training)
        ### APPEND TO FILE FOLDS
        write(internal_folds, file = sprintf("%s_list_folds_n%s.txt", save_name,n_folds), sep = ",", ncolumns = length(internal_folds), append = TRUE)
      }
      
      
      
      # multiple imputation
      #data_table_imputed_train = matrix(data_informative_table[training,], ncol=ncol(data_informative_table))
      #data_table_imputed_test = matrix(data_informative_table[test,], ncol=ncol(data_informative_table))
      n_multiple_imputations = 1
      if(!is.null(data_informative_table)){
        data_table_imputed_train = data_informative_table[training,,drop = F]
        data_table_imputed_test = data_informative_table[test,,drop=F]
        if (any(is.na(data_informative_table))) {
          n_multiple_imputations = N_M_IMPUTATIONS #20
          if (length(Sys.glob(sprintf("%s_fold%d_iteration*.txt", save_name, fold, iter)))!=N_ITERATIONS) {
            impute_obj = impute.glmnet.matrix_fit(x = data_table_imputed_train, n_cores = n_cores)
            imp.data_informative_table_train = impute.glmnet.matrix(m = impute_obj, x = data_table_imputed_train, nimp = N_M_IMPUTATIONS)
            imp.data_informative_table_test = impute.glmnet.matrix(m = impute_obj, x = data_table_imputed_test, nimp = N_M_IMPUTATIONS)
          }
        }
      } else{
        data_table_imputed_train = NULL
        imp.data_informative_table_train = NULL
        data_table_imputed_test = NULL
        imp.data_informative_table_test = NULL
      }
      start_fold.time <- Sys.time()
      
      train_linear_predictor_to_find_the_threshold_all = vector(mode = 'list', length = length(time_points))
      
      # for (i_time_point in 1:length(time_points)) {
      #   train_linear_predictor_to_find_the_threshold_all[[i_time_point]] = c()
      # }
      
      #####################################################################################
      #####################################################################################
      ### CROSSVALIDATION ###
      for (iter in 1:N_ITERATIONS) {
        #sprintf("%s_iteration%d_%d.txt", save_name, iter, fold, iter)
        if (!file.exists(sprintf("%s_fold%d_iteration%d.txt", save_name, fold, iter))) {
          
          
          
          linPreds = c()
          if (mp$response_family == "cox") {
            train_thresholds = NULL
          }
          start.time <- Sys.time()
          
          if (N_ITERATIONS > 1)
            cat("Fold:", fold, ". Iteration ", iter, "of", N_ITERATIONS, "\n", save_name)
          
          sprintf('/nIteration: %d', iter)
          ##########################################################################
          # IF MULTIPLE IMPUTATION. PRED = MEAN ( PRED_PER_FOLD )
          ##########################################################################
          preds = c()
          
          for (iter_imputation in 1:n_multiple_imputations) {    
            
            if (n_multiple_imputations > 1) {
              data_table_imputed_train = imp.data_informative_table_train[[iter_imputation]]
              #data_table_imputed_test = as.matrix(imp.data_informative_table_test[[iter_imputation]])
            }
            if (n_multiple_imputations > 1)
              .print_action(paste("Fold:",fold,"Imputation", iter_imputation, "of", n_multiple_imputations, "\n", save_name))
            
  
            model_list = fit_model(mp = mp,
                                   data_informative_table = data_table_imputed_train,
                                   Y = trainY,
                                   mri = mri$data[,,, training],
                                   mri_fu = mri_fu$data[,,, training],
                                   mri_wm = mri_wm$data[,,, training],
                                   mri_fu_wm = mri_wm_fu$data[,,, training],
                                   preloaded_covB = preloaded_covB,
                                   preloaded_covB_fu = preloaded_covB_fu,
                                   iter = ifelse(use_ensemble_voxels == FALSE, -1, iter),
                                   internal_folds = internal_folds,
                                   n_cores = n_cores,
                                   use_significant_voxels = use_significant_voxels,
                                   covX_site = sites_training,
                                   #name_combat = name_combat,
                                   standardize_images = standardize_images
            )
            mp$combat = model_list$combat
            if (!mp$modulation %in% c('cl','clinical')) {
              # per passar a mni la ijk =>  mp$mni = (mri$sto.xyz %*% mp$ijk)[1:3, ]          
              lasso_ijk = matrix(model_list$lasso_ijk, nrow = 3)
              # lasso_ijk està amb totes les imatges juntes
              ##############################################################
              dim_x = dim(mri$data)[1]
              if (mp$modulation == 'fu') {
                idx_gm_mod = which(lasso_ijk[1,] > dim_x * 0 & lasso_ijk[1,] <= dim_x * 1)
                idx_gm = which(lasso_ijk[1,] > dim_x * 1 & lasso_ijk[1,] <= dim_x * 2)
                if (length(idx_gm) > 0) {
                  ijk_gm     = lasso_ijk[, idx_gm] - c(dim_x * 0, 0, 0)
                  model_list$lasso_mni$gm     = ijk2mni(rbind(matrix(ijk_gm, nrow = 3), 1), mri$sto.ijk)[, 1:length(idx_gm)]
                  model_list$lasso_mni$gm_betas = model_list$lasso$beta[idx_gm]
                } else {
                  model_list$lasso_mni$gm     = NULL
                  model_list$lasso_mni$gm_betas = NULL
                }
                ijk_gm_mod = lasso_ijk[, idx_gm_mod] - c(dim_x * 1, 0, 0)
                if(length(idx_gm_mod)>0){
                  model_list$lasso_mni$gm_mod = ijk2mni(rbind(matrix(ijk_gm_mod, nrow = 3), 1), mri$sto.ijk)[, 1:length(idx_gm_mod)]
                  model_list$lasso_mni$gm_mod_betas = model_list$lasso$beta[idx_gm_mod]
                } else {
                  model_list$lasso_mni$gm_mod = NULL
                  model_list$lasso_mni$gm_mod_betas = NULL
                }
              } else {
                idx_gm = which(lasso_ijk[1,] > dim_x * 0 & lasso_ijk[1,] <= dim_x * 1)
                idx_gm_mod = which(lasso_ijk[1,] > dim_x * 1 & lasso_ijk[1,] <= dim_x * 2)
                
                
                ijk_gm     = lasso_ijk[, idx_gm] - c(dim_x * 0, 0, 0)
                if (length(idx_gm) > 0) {
                  model_list$lasso_mni$gm     = ijk2mni(rbind(matrix(ijk_gm, nrow = 3), 1), mri$sto.ijk)[, 1:length(idx_gm)]
                  model_list$lasso_mni$gm_betas = model_list$lasso$beta[idx_gm]
                } else {
                  model_list$lasso_mni$gm     = NULL
                  model_list$lasso_mni$gm_betas = NULL
                }
                if(length(idx_gm_mod)>0){
                  ijk_gm_mod = lasso_ijk[, idx_gm_mod] - c(dim_x * 1, 0, 0)
                  model_list$lasso_mni$gm_mod = ijk2mni(rbind(matrix(ijk_gm_mod, nrow = 3), 1), mri$sto.ijk)[, 1:length(idx_gm_mod)]
                  model_list$lasso_mni$gm_mod_betas = model_list$lasso$beta[idx_gm_mod]
                } else {
                  model_list$lasso_mni$gm_mod = NULL
                  model_list$lasso_mni$gm_mod_betas = NULL
                }
              }
              idx_wm = which(lasso_ijk[1,] > dim_x * 2 & lasso_ijk[1,] <= dim_x * 3)
              if(length(idx_wm)>0){
                ijk_wm     = lasso_ijk[, idx_wm] - c(dim_x * 2, 0, 0)
                model_list$lasso_mni$wm     = ijk2mni(rbind(matrix(ijk_wm, nrow=3),1), mri$sto.ijk)[, 1:length(idx_wm)]
                model_list$lasso_mni$wm_betas = model_list$lasso$beta[idx_wm]
              } else {
                model_list$lasso_mni$wm     = NULL
                model_list$lasso_mni$wm_betas = NULL
              }
              idx_wm_mod = which(lasso_ijk[1,] > dim_x * 3 & lasso_ijk[1,] <= dim_x * 4)
              if(length(idx_wm_mod)>0){
                ijk_wm_mod = lasso_ijk[, idx_wm_mod] - c(dim_x * 3, 0, 0)
                model_list$lasso_mni$wm_mod =  ijk2mni(rbind(matrix(ijk_wm_mod, nrow=3),1), mri$sto.ijk)[, 1:length(idx_wm_mod)]
                model_list$lasso_mni$wm_mod_betas = model_list$lasso$beta[idx_wm_mod]
              } else {
                model_list$lasso_mni$wm_mod =  NULL
                model_list$lasso_mni$wm_mod_betas = NULL
              }
              ##############################################################
              
              signif_indx = model_list$signif_indx 
              lasso_covB = model_list$lasso_covB
              lasso_covB_fu = model_list$lasso_fu_covB
              
              n_voxels_mask = model_list$n_voxels_mask
              #################################################################################################################################################
              #data_table_imputed_train = model_list$data_table_imputed_train
              #data_table_imputed_test = .impute_Xtest(data_table_imputed_train,data_informative_table[test,], n_iter=20)
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
            
            if (mp$response_family == "cox") {
              #predX_training= matrix(data_table_imputed_train[, .indx_equivalent_columns_from_transf(mp$data_table_transf, mp$pred_var)],nrow=nrow(data_table_imputed_train))
              predX_training= matrix(data_table_imputed_train[, mp$pred_var],nrow=nrow(data_table_imputed_train))
              # Apply the model to the training, and then find the optimal threshold to separate high vs low HR
              #covX_training = cbind(1, as.matrix(data_table_imputed_train[, .indx_equivalent_columns_from_transf(mp$data_table_transf, mp$covX_var)]))
              covX_training = cbind(1, as.matrix(data_table_imputed_train[, mp$covX_var]))
              # XXX PASSAR LA VARIABLE TIME_POINTS --- ATENCIO, SI NO HI HA PROU EVENTS AL TIME POINT EL THREHSOLD ES CALCULA MALAMENT
              for (i_time_point in 1:length(time_points)) {
                status0_at_time_point = which((trainY[, 1] == time_points[i_time_point] & trainY[, 2] == 0) | (trainY[, 1] > time_points[i_time_point]))
                status1_at_time_point = which(trainY[, 1] <= time_points[i_time_point] & trainY[, 2] == 1)
                if (length(status0_at_time_point) > 0 && length(status1_at_time_point) > 0) {
                  
                  any_status_at_time_point = c(status0_at_time_point, status1_at_time_point)
                  training_for_threshold = training[any_status_at_time_point]
                  
                  train_linear_predictor_to_find_the_threshold = apply_model(mp = mp,
                                                                             mri_data = mri$data[,,, training_for_threshold],
                                                                             mri_fu_data = mri_fu$data[,,, training_for_threshold],
                                                                             mri_wm_data = mri_wm$data[,,, training_for_threshold],
                                                                             mri_wm_fu_data = mri_wm_fu$data[,,, training_for_threshold],
                                                                             covX_test = covX_training[any_status_at_time_point,],
                                                                             signif_indx = signif_indx,
                                                                             lasso_covB = lasso_covB,
                                                                             lasso_covB_fu = lasso_covB_fu,
                                                                             mask = model_list$mask,
                                                                             predX_test = predX_training[any_status_at_time_point,],
                                                                             scale_clinical = scale_clinical,
                                                                             scale_mri = scale_mri,
                                                                             lasso = model_list$lasso,
                                                                             lasso_predX_indx = model_list$lasso_predX_indx,
                                                                             tipett_take_un = tipett_take_un,
                                                                             img_kappa = NULL,
                                                                             use_significant_voxels = use_significant_voxels,
                                                                             #covX_site = sites_training[[any_status_at_time_point]],
                                                                             covX_site = sites_training[any_status_at_time_point],
                                                                             masks_3d = model_list$masks_3d, # ALEIX, AIXÒ ÉS PEL COMBAT, POTSER NO HO NECESSITAREM,
                                                                             #name_combat = name_combat,
                                                                             n_voxels_mask = n_voxels_mask,
                                                                             combat = model_list$combat,
                                                                             standardize_images = standardize_images
                  )
                  train_linear_predictor_to_find_the_threshold_all[[i_time_point]] = cbind(train_linear_predictor_to_find_the_threshold_all[[i_time_point]], train_linear_predictor_to_find_the_threshold)
                  ################################################## copiat a fora els loops
                  ##################################################
                }
              }
            }
            
            
            mp$models[[model_counter]] <- model_list 
            model_counter <- model_counter + 1
            ###########################################################################################################################################################
            .print_ok()
            ### END TRAINING ###
            ### TESTING ###
            # WE APPLY EVERY IMPUTATION TRAINING TO EVERY IMPUTATION TEST, 20 TRAIN MODELS * 20 TEST IMPUTATIONS         
            .print_action("Test sample: applying the model")
            test_preds = c()
            for (iter_imputation_test in 1:n_multiple_imputations) {
              if (n_multiple_imputations > 1) {
                data_table_imputed_test = as.matrix(imp.data_informative_table_test[[iter_imputation_test]])
              }
              #covX_test = cbind(1, matrix(data_table_imputed_test[, .indx_equivalent_columns_from_transf(mp$data_table_transf, mp$covX_var)],nrow=nrow(data_table_imputed_test)))
              
              if(!is.null(mp$covX_var)) {
                covX_test = cbind(1, matrix(data_table_imputed_test[, mp$covX_var],nrow=nrow(data_table_imputed_test)))
              } else {
                covX_test = matrix(1, nrow = length(test))
              }
              #predX_test = matrix(data_table_imputed_test[, .indx_equivalent_columns_from_transf(mp$data_table_transf, mp$pred_var)],nrow=nrow(data_table_imputed_test))
              if(!is.null(mp$pred_var)) {
                predX_test = matrix(data_table_imputed_test[, mp$pred_var],nrow=nrow(data_table_imputed_test))
              } else {
                predX_test = NULL
              }
              preds = apply_model(mp = mp,
                                  mri_data = mri$data[,,, test,drop=FALSE],
                                  mri_fu_data = mri_fu$data[,,, test,drop=FALSE],
                                  mri_wm_data = mri_wm$data[,,, test,drop=FALSE],
                                  mri_wm_fu_data = mri_wm_fu$data[,,, test,drop=FALSE],
                                  covX_test = covX_test,
                                  signif_indx = signif_indx,
                                  lasso_covB = lasso_covB,
                                  lasso_covB_fu = lasso_covB_fu,
                                  mask = model_list$mask,
                                  predX_test = predX_test,
                                  scale_clinical = scale_clinical,
                                  scale_mri = scale_mri,
                                  lasso = model_list$lasso,
                                  lasso_predX_indx = model_list$lasso_predX_indx,
                                  tipett_take_un = tipett_take_un,
                                  img_kappa = NULL,
                                  use_significant_voxels = use_significant_voxels,
                                  covX_site = sites_test,
                                  masks_3d = model_list$masks_3d, # ALEIX, AIXÒ ÉS PEL COMBAT, POTSER NO HO NECESSITAREM,
                                  #name_combat = name_combat,
                                  n_voxels_mask = n_voxels_mask,
                                  standardize_images = standardize_images,
                                  combat = model_list$combat
              )
              test_preds = cbind(test_preds, preds)
            }
            preds = rowMeans(test_preds)
            linPreds = cbind(linPreds, preds)
            if(ide_shiny)
              incProgress(1/(n_folds*N_ITERATIONS*n_multiple_imputations), detail = paste("Doing crossvalidation. Fold",fold, ". Iteration:", iter, ". Imputation:",iter_imputation))
          }
          # fi for n_multiple_imputations  
          
         
          
          # average of linear predictors in multiple_imputation
          linPred = matrix(rowMeans(linPreds))
          
          .print_ok()
          .print_action("Saving the predictions")
          cat('\n')
          pred = c()
          if (mp$response_family == "binomial") {
            prob = 1 / (1 + exp(-linPred))
            pred = prob
            bac = .metrics_binary(testY, pred > (sum(trainY == 1) / length(trainY)))$bac
            cv_accuracy[fold] = bac # accuracy
            cv_table[test, 2] = prob
          }
          else if (mp$response_family == "gaussian") {
            pred = as.matrix(linPred)
            cv_accuracy[fold] = sqrt(mean((pred - testY) ^ 2))
            cv_table[test, 2] = pred
          }
          else if (mp$response_family == "cox") {
            pred = linPred
            #rd2 = .coxph_RD2(pred, testY[,1], testY[,2])
            
            cv_accuracy[fold] = NA
            cv_betas[fold] = NA
            cv_table[test, 2] = pred
            
            results = data.frame(id = test, linear_predictor = pred[, 1], time = testY[, 1], status = testY[, 2])
            if (n_multiple_imputations > 1)
              print('Mean results for multiple imputations:')
            
          }
          ##########################################################################
          # FI FOR MULTIPLE IMPUTATION. PRED = MEAN ( PRED_PER_FOLD )
          ##########################################################################
          
          if (mp$response_family == "cox") {
            
            cv_table_predictions = rbind(cv_table_predictions, cbind(
              data.frame(id = test, linear_predictor = pred[, 1], rd2_beta = NA, time = testY[, 1], status = testY[, 2], fold = fold, iteration = iter)
              #,preds_status_at_time_points
            )) #, z_cuts = z_cuts[iter]))
          } else if (mp$response_family == "gaussian") {
            cv_table_predictions = rbind(cv_table_predictions, data.frame(id = test, linear_predictor = pred[, 1], response = testY, fold = fold, iteration = iter)) #, z_cuts = z_cuts[iter]))
          } else {
            cv_table_predictions = rbind(cv_table_predictions, data.frame(id = test, linear_predictor = linPred[, 1], prob = pred[, 1], response = testY, fold = fold, iteration = iter)) #, z_cuts = z_cuts[iter]))
          }
          .print_ok()
          switch(mp$response_family,
                 binomial = message("BAC: ", cv_accuracy[fold]),
                 gaussian = message("RMSE: ", cv_accuracy[fold]),
                 cox = { }) # message("RD2: ", cv_accuracy[fold]))
          end.time <- Sys.time()
          
          #### fold cox results
          if (mp$response_family == 'cox') {
            saveRDS(object = train_linear_predictor_to_find_the_threshold_all, file = sprintf("%s_train_linpred_to_find_threhshold_FOLD_%d_ITER_%d_IMP_%d.rds", save_name, fold, iter, iter_imputation))
            frontier_time = .find_best_time(testY[, 1], testY[, 2])
            metrics_cox = .metrics_cox(results = data.frame(id = test, linear_predictor = pred[, 1], rd2_beta = NA, time = testY[, 1], status = testY[, 2], fold = fold), frontier_time = frontier_time, iteration = iter, folder = sprintf("%s_iteration%d", save_name, iter), save = TRUE)
            write.csv(metrics_cox, sprintf("%s_fold%d_iteration%d.txt", save_name, fold, iter), row.names = FALSE)
            cat("Fold", fold, "performance\n")
            
            print(as.data.frame(metrics_cox))
          }
          #######################
          if (iter > 1) message("Time per iteration: ", round(difftime(end.time, start.time, units = 'mins'),2), ' mins.')
          
          
          # fi ITERATIONS (BRAIN CUTS)
          #write.csv(t(matrix(unlist(subjects_used),nrow=N_ITERATIONS)),file=sprintf('subjects_sampled_fold_%s',fold), col.names = c(1:18))
          #####################################################################################
          #####################################################################################
          
          
          # if (mp$response_family == "cox") {
          #   preds_status_at_time_points = NULL
          #   
          #   for (i in 1:ncol(train_thresholds)) {
          #     preds_status_at_time_points = cbind(preds_status_at_time_points, linPred > mean(train_thresholds[,i]))
          #   }
          # }
        }
        else {
          cat('Skipping fold: ', fold, ' Iteration: ', iter, '\n') # aquí hauriem de carregar a cv_table... el fitxer .csv que ha trobat guardat
          res_csv = read.csv(sprintf("%s_iteration%d_%s_results_fold%d.csv", save_name, iter, mp$response_family, fold))
          cv_table_predictions = rbind(cv_table_predictions, res_csv)
          name_base <- sprintf("%s_model_FOLD_%d", save_name,fold)
          #mp$models <- readRDS(file = sprintf("%s_model_list.rds", name_base))
          mp <- readRDS(file=sprintf("%s_mp.rds",name_base))
        }
        
      }
      
      # save rds
      cat('\n[Saving] - Saving fold model to',sprintf("%s_model_FOLD_%d", save_name,fold))
      name_base <- sprintf("%s_model_FOLD_%d", save_name,fold)
      if(fold>0) name_base_previous <- sprintf("%s_model_FOLD_%d", save_name,fold-1)
      saveRDS(object = mp, file = sprintf("%s_mp.rds", name_base))
      #saveRDS(object = mp$models, file = sprintf("%s_model_list.rds", name_base))
      #file.remove(sprintf("%s_mp.rds", name_base_previous))
      if (mp$response_family == "cox") {
        trainY = Y[training,]
        train_thresholds_row = c()
        for (i_time_point in 1:length(time_points)) {
          status0_at_time_point = which((trainY[, 1] == time_points[i_time_point] & trainY[, 2] == 0) | (trainY[, 1] > time_points[i_time_point]))
          status1_at_time_point = which(trainY[, 1] <= time_points[i_time_point] & trainY[, 2] == 1)
          if (!is.null(train_linear_predictor_to_find_the_threshold_all[[i_time_point]])) {
            train_linear_predictor_to_find_the_threshold = rowMeans(train_linear_predictor_to_find_the_threshold_all[[i_time_point]])
          } else {
            train_linear_predictor_to_find_the_threshold = NA
          }
          sorted_train_linear_predictor_to_find_the_threshold = sort(unique(train_linear_predictor_to_find_the_threshold))
          thresholds = (sorted_train_linear_predictor_to_find_the_threshold[-length(sorted_train_linear_predictor_to_find_the_threshold)] + sorted_train_linear_predictor_to_find_the_threshold[-1]) / 2
          bacs = c()
          for (threshold in thresholds) {
            true_status = c(rep(0, length(status0_at_time_point)), rep(1, length(status1_at_time_point)))
            predicted_status = 1 * (train_linear_predictor_to_find_the_threshold > threshold)
            sensitivity = sum(predicted_status == 1 & true_status == 1) / sum(true_status == 1)
            specificity = sum(predicted_status == 0 & true_status == 0) / sum(true_status == 0)
            bac = (sensitivity + specificity) / 2
            bacs = c(bacs, bac)
          }
          train_thresholds_row = c(train_thresholds_row, thresholds[which.max(bacs)])
        }
        train_thresholds_all = rbind(train_thresholds_all, train_thresholds_row)
        
      }
      message("\n[Fold ended] - Time per fold:", difftime(Sys.time(), start_fold.time, units = 'mins'), ' mins.')
      save(list=ls(), 
           file =sprintf('%s_fold%s.Rdata', save_name, fold))
    }
    else{
      load(sprintf('%s_fold%s.Rdata', save_name, fold))
    }
  }
  # fi cv  
  
  
  if(!is.null(train_thresholds_all))
    write.csv(train_thresholds_all, sprintf("%s_train_thresholds.csv", save_name), row.names = FALSE)
  
  if (mp$response_family == "cox") {
    mean_pred_subjs = data.frame(id = sort(unique(cv_table_predictions$id)),
                                 lin_pred = tapply(cv_table_predictions$linear_predictor, cv_table_predictions$id, mean), # returns means sorted by $id
                                 times = tapply(cv_table_predictions$time, cv_table_predictions$id, unique), # times sorted by $id
                                 status = tapply(cv_table_predictions$status, cv_table_predictions$id, unique), # status sorted by $id
                                 fold = tapply(cv_table_predictions$fold, cv_table_predictions$id, unique) # fold sorted by $id
    )
    
    #train_thresholds_means = colMeans(train_thresholds_all)
    n_pred_subjs_cols = ncol(mean_pred_subjs)
    th_i = 1
    for (time_point in time_points) {
      for (fold in 1:n_folds) {
        mean_pred_subjs[mean_pred_subjs$fold == fold, n_pred_subjs_cols + th_i] <- NA
        if (ncol(train_thresholds_all) >= th_i) {
          mean_pred_subjs[mean_pred_subjs$fold == fold, n_pred_subjs_cols + th_i] <- mean_pred_subjs$lin_pred[mean_pred_subjs$fold == fold] > train_thresholds_all[fold, th_i]
        }
        colnames(mean_pred_subjs)[n_pred_subjs_cols + th_i] <- sprintf("threshold_%s", time_point)
      }
      th_i = th_i + 1
    }
    mp$cv_results = mean_pred_subjs
  } else {
    mp$cv_results = cv_table_predictions
  }
  # for(time_point in time_points){
  #   mean_pred_subjs[,length(colnames(mean_pred_subjs))+1] <- mean_pred_subjs$lin_pred>train
  #   colnames(mean_pred_subjs)[length(colnames(mean_pred_subjs))] <- sprintf("threshold_%s",time_point)
  # }
  ### END CROSSVALIDATION ###
  mp$cv_table = cv_table;
  mp$cv_accuracy = cv_accuracy;
  #;

  mp$cv_betas = cv_betas;
  mp$subjects_used = subjects_used;
  cat("\n[End] - FOLDS performance:", cv_accuracy, "\n")
  flush.console()
  switch(mp$response_family,
         binomial = {
           bin_threshold = sum(cv_table_predictions$response == 1) / length(cv_table_predictions$response)
           message("Mean BAC: ", .metrics_binary(cv_table_predictions$response, cv_table_predictions$prob > bin_threshold)$bac)
         },
         gaussian = message("Mean RMSE: ", sqrt(mean((cv_table_predictions$linear_predictor - cv_table_predictions$response) ^ 2))),
         cox = {
           #mean_rd2 = .mean_RD2(cv_table_predictions$rd2_beta)
           final_rd2 = .coxph_RD2(predictor = mp$cv_results$lin_pred, stime = mp$cv_results$time, sevent = mp$cv_results$status)
           #message("Mean RD2: ", mean_rd2, "\n")
           message("Final RD2: ", final_rd2$RD2, " Beta:", final_rd2$b)
         }
  )
  
  
  if (mp$response_family == "cox") {
    mp$frontier_time = .find_best_time(trainY[, 1], trainY[, 2])
    mp$metrics = .metrics_cox(cv_table_predictions, mp$frontier_time, save = FALSE)
    write.csv(mp$metrics, sprintf("%s_custom_metric.txt", save_name))
    
    write.csv(.coxph_RD2(predictor = mp$cv_results$lin_pred, stime = mp$cv_results$time, sevent = mp$cv_results$status), sprintf("%s_rd2.txt", save_name))
    #write(cat("Mean RD2: ", mean_rd2, "\n"), sprintf("%s_rd2.txt",save_name))
    print(mp$metrics)
    #colnames(mp$cv_results) = c('id','pred','time','status','fold','iteration')
  } else if (mp$response_family == "gaussian") {
    colnames(mp$cv_results) = c('id', 'pred', 'real', 'fold', 'iteration')
    mp$cv_results = aggregate(mp$cv_results[, 2:4], list(mp$cv_results$id), mean)
    colnames(mp$cv_results) = c('id', 'prediction', 'real', 'fold')
    mp$metrics = sqrt(mean((mp$cv_results$prediction - mp$cv_results$real) ^ 2))
    cat(sprintf("Mean RMSE: %f", mp$metrics), file = sprintf("%s_gaussian.txt", save_name))
  } else {
    write.csv(mp$cv_results, sprintf("%s_bin_raw.csv",save_name))
    mp$cv_results = aggregate(mp$cv_results[, 2:5], list(mp$cv_results$id), mean)
    colnames(mp$cv_results) = c('id', 'linear_predictor', 'probability', 'real', 'fold')
    mp$metrics <- .metrics_binary(mp$cv_results$real, mp$cv_results$probability > bin_threshold)
    colnames(mp$cv_results)[3] <- sprintf("Prediction (0 = %s; 1 = %s)",mp$response_ref, mp$response_event)
    colnames(mp$cv_results)[4] <- sprintf("Real (0 = %s; 1 = %s)",mp$response_ref, mp$response_event)
    cat(sprintf("Mean BAC: %f", mp$metrics$BAC, file = sprintf("%s_bin.txt", save_name)))
    #write.csv(mp$cv_results, sprintf("%s_bin.txt",save_name))
  }
  cat('\n[Saving] - Saving MRIPredict object to:', sprintf('%s_mp.rds',save_name))
  mp$cv_results = mp$cv_results[order(mp$cv_results$id),]
  saveRDS(mp,file=sprintf('%s_mp.rds',save_name))
  #file.remove(sprintf("%s_model_FOLD_%s_model_list.rds", save_name,n_folds))
  file.remove(sprintf("%s_model_FOLD_%s_mp.rds",save_name,n_folds-1))
  mp
}


if(F) {
  source('mripredict.r')
  source('mripredict_fit.r')
  source('mripredict_predict.r')
  source('mripredict_library.r')
  source('mripredict_library.r')
  source('mripredict.r')
  source('rotation3d.r')
  source('cv_functions.r')
  source('predict.r')
  source('fit.r')
  source('combat_quim.R')
  source('combat_utils.R')
  source('glmnet.utils.R')
  source('extra.R')
  mp <- readRDS("mp_debugfit.rds")
  mp <- mripredict_cv(mp, space = "NO_CHECK", 
                      save_name = "testing", 
                      folds_file = NULL, 
                      n_cores = 1, 
                      n_folds = 10,
                      use_significant_voxels = FALSE, 
                      use_ensemble_learning = F, 
                      use_ensemble_voxels = F, 
                      use_ensemble_subjects = FALSE, 
                      ide_shiny = F)
}

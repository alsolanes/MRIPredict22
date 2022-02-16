fit_model = function(mp, data_informative_table , Y, mri, mri_fu, preloaded_covB=NULL, preloaded_covB_fu=NULL, 
                     SIGNIFICANCE_THRESHOLD=qnorm(0.975), iter=-1, internal_folds = c(), n_cores=1, 
                     use_significant_voxels=FALSE, covX_site = c(), mri_wm=NULL,mri_fu_wm=NULL,
                     name_combat='_combat', mask=NULL, standardize_images = FALSE) {
  tipett_take_un=c()
  
  # if (any(is.na(data_informative_table))){
  #   
  #   # FER AIXO UN SOL COP PER TOITES LES IMPUTACIONS:
  #   impute_obj = .impute_fit(data_informative_table_train)
  #   imp.data_informative_table_train = .impute_apply(impute_obj, data_informative_table_train)
  #   imp.data_informative_table_test = .impute_apply(impute_obj, data_informative_table_test)
  #   # FER AIXO PER CADA IMPUTATION SET (canviant [[1]] per [[2]], etc:
  #   data_table_imputed_train = imp.data_informative_table[[1]]
  #   data_table_imputed_test = imp.data_informative_table_tes[[1]]
  #   
  #   warning('Check the imputation, still not tested.')
  # } else {
  #   data_table_imputed_train = data_informative_table
  # }
  #

  if(!is.null(data_informative_table)){
    if (ncol(data_informative_table)>1) { # if there are more than one variable
      covX_training = as.matrix(cbind(1,data_informative_table[,mp$covX_var]))
      predX_train = as.matrix(data_informative_table[, mp$pred_var])
    } else {
      # ALEIX, IMPUTAR SI HI HA MISSINGS (?)
      covX_training = as.matrix(data_informative_table)
    }
    data_informative_table = as.matrix(data_informative_table)
  } else {
    covX_training = matrix(1, nrow=nrow(as.matrix(Y)))
  }
  
  .print_ok()
  ### TRAINING ###
  
  #.print_action("Training sample: masking\n")
  # Initialize masks
  masked_X = c()
  masked_X_sign = c()
  ### PREPROCESS MRI VOXELS
  
  # Scale mri images
  # STANDARDIZE IMAGES
  scale_mri = list()
  scale_mri$sd_mri = NA
  scale_mri$sd_mri_fu = NA
  scale_mri$sd_mri_wm = NA
  scale_mri$sd_mri_wm_fu = NA
  scale_mri$mean_mri = NA
  scale_mri$mean_mri_fu = NA
  scale_mri$mean_mri_wm = NA
  scale_mri$mean_mri_wm_fu = NA
  if (standardize_images){
    scale_mri$sd_mri = sd(mri)
    scale_mri$mean_mri = mean(mri)
    mri = mri - scale_mri$mean_mri
    if(!is.null(mri_fu)){
      scale_mri$sd_mri_fu = sd(mri_fu)
      scale_mri$mean_mri_fu = mean(mri_fu)
      mri_fu = (mri_fu - scale_mri$mean_mri_fu)/scale_mri$sd_mri_fu * scale_mri$sd_mri
    }
    if(!is.null(mri_wm)){
      scale_mri$sd_mri_wm = sd(mri_wm)
      scale_mri$mean_mri_wm = mean(mri_wm)
      mri_wm = (mri_wm - scale_mri$mean_mri_wm)/scale_mri$sd_mri_wm * scale_mri$sd_mri
    }
    if(!is.null(mri_fu_wm)){
      scale_mri$sd_mri_wm_fu = sd(mri_wm_fu)
      scale_mri$mean_mri_wm_fu = mean(mri_wm_fu)
      mri_fu_wm = (mri_fu_wm - scale_mri$mean_mri_wm_fu)/scale_mri$sd_mri_wm_fu * scale_mri$sd_mri
    }
  }
  
  
  signif_indx=c()
  
  masks_3d = list()
  combat = list()
  if (!is.null(mp$mri_paths) && mp$mri_paths != "") {
    # Create mask
    if (!is.null(mask)) {
      mask = mask
    } else {
      mask = c()
      mask = array(TRUE, dim = dim(mri[,,,1]))
      mask_o = mask
    }
    # FIRST SINGLE MASK USED TO MAKE 'CUTS' TO THE BRAIN
    mask = mask & (apply(mri, 1:3, var) > 0)
    mask_ijk = which(mask, arr.ind = TRUE) # cooridnates 3d of the mask 
    if (iter != -1) {
      .print_action('Training sample: cutting brain')
      mask = rotate_coordinates(list_coordinates = mask_ijk, iter = iter, return_3d=TRUE, img_3d = mask_o)>0
    }
    
    # COMBAT AND MASK (sd>0) FOR EACH modality
    n_image_modalities = 0
    if (!is.null(mri)){
      n_image_modalities = n_image_modalities + 1
      masks_3d = list()
      if (!is.null(covX_site)) {
        .print_action('Applying combat')#'Removing effects of covariates (for UN+FU-GM+WM) - False, only applying combat')
      }
      list_masked = remove_effects(mask_3d = mask, mri = mri, covX_training = covX_training, trainY = Y, response_family=mp$response_family, REMOVE_EFFECT_TO_X = FALSE,
                                   n_cores = n_cores, SIGNIFICANT = (use_significant_voxels || mp$modulation=='op'), covX_site = covX_site, modality = paste0("_un_gm",name_combat))
      combat$gm = list_masked$combat
      img_3d = mri[,,,1]
      img_3d = img_3d * 0
      for (img_i in 1:dim(mri)[4]){ #for every image
        img_3d[which(list_masked$mask_3d)] = list_masked$X[img_i,]
        mri[,,,img_i] = img_3d
      }
      masks_3d$un_gm = which(list_masked$mask_3d) # This mask is the intersection between 'mask' and voxels sd > 0 from 'remove_effects'
    }
    if (!is.null(mri_fu)) {
      n_image_modalities = n_image_modalities + 1
      list_masked = remove_effects(mask_3d = mask, mri = mri_fu, covX_training = covX_training, trainY = Y, response_family=mp$response_family, REMOVE_EFFECT_TO_X = FALSE,
                                   n_cores = n_cores, SIGNIFICANT = (use_significant_voxels || mp$modulation=='op'), covX_site = covX_site, modality = paste0("_fu_gm",name_combat))
      combat$gm_mod = list_masked$combat
      img_3d = mri_fu[,,,1]
      img_3d = img_3d * 0
      for (img_i in 1:dim(mri_fu)[4]){ #for every image
        img_3d[which(list_masked$mask_3d)] = list_masked$X[img_i,]
        mri_fu[,,,img_i] = img_3d
      }
      masks_3d$fu_gm = which(list_masked$mask_3d)
    }
    if (!is.null(mri_wm)) {
      n_image_modalities = n_image_modalities + 1
      list_masked = remove_effects(mask_3d = mask, mri = mri_wm, covX_training = covX_training, trainY = Y, response_family=mp$response_family, REMOVE_EFFECT_TO_X = FALSE,
                                   n_cores = n_cores, SIGNIFICANT = (use_significant_voxels || mp$modulation=='op'), covX_site = covX_site, modality = paste0("_un_wm",name_combat))
      combat$wm = list_masked$combat
      img_3d = mri_wm[,,,1]
      img_3d = img_3d * 0
      for (img_i in 1:dim(mri_wm)[4]){ #for every image
        img_3d[which(list_masked$mask_3d)] = list_masked$X[img_i,]
        mri_wm[,,,img_i] = img_3d
      }
      masks_3d$un_wm = which(list_masked$mask_3d)
    }
    if (!is.null(mri_fu_wm)) {
      n_image_modalities = n_image_modalities + 1
      list_masked = remove_effects(mask_3d = mask, mri = mri_fu_wm, covX_training = covX_training, trainY = Y, response_family=mp$response_family, REMOVE_EFFECT_TO_X = FALSE,
                                   n_cores = n_cores, SIGNIFICANT = (use_significant_voxels || mp$modulation=='op'), covX_site = covX_site, modality = paste0("_fu_wm",name_combat))
      combat$wm_mod = list_masked$combat
      img_3d = mri_fu_wm[,,,1]
      img_3d = img_3d * 0
      for (img_i in 1:dim(mri_fu_wm)[4]){ #for every image
        img_3d[which(list_masked$mask_3d)] = list_masked$X[img_i,]
        mri_fu_wm[,,,img_i] = img_3d
      }
      masks_3d$fu_wm = which(list_masked$mask_3d)
    }
    
    
    mask_single_image <- mask
    # tindrem 4 mascares, les creem a new_mask
    #writeNIfTI(nifti(mask*mri[,,,1],sform=1,datatype = 16), filename=sprintf("mask_%s",iter), onefile = TRUE, gzipped = TRUE, verbose = FALSE, warn = -1, compression = 6)
    ################################################################################################################################      
    
    # WE APPEND THE DIFFERENT IMAGES (GM, GM_MOD, WM, WM_MOD)
    # THE SAME WITH THE DIFFERENT MASKS
    if (mp$modulation != 'un' && mp$modulation!='fu'){ # if we have mri_Wm means we are not in 'un' or 'fu' modulation
      orig_dims = dim(mask)
      new_dims = orig_dims * c(n_image_modalities,1,1)
      new_mask = array(NA, dim = new_dims)
      for (img_i in 1:n_image_modalities) { # for every image (wm, gm, ...)
        mask_i <- switch(as.character(img_i),
                         '1'=masks_3d$un_gm,
                         '2'=masks_3d$fu_gm,
                         '3'=masks_3d$un_wm,
                         '4'=masks_3d$fu_wm)
        first_index_x = (((img_i-1)*(orig_dims[1]))+1)
        new_mask[first_index_x:(first_index_x+orig_dims[1]-1),,] = .mask_from_which(mask_i, orig_dims[1])
      }
      mask = new_mask
      new_mask <- NULL
      ##
      new_mri = array(NA, dim = dim(mri) * c(n_image_modalities,1,1,1))
      if(n_image_modalities>0){
        img_j = 1
        first_index_x = (((img_j-1)*(orig_dims[1]))+1)
        new_mri[first_index_x:(first_index_x+orig_dims[1]-1),,,] = mri
      }
      if(n_image_modalities>1){
        img_j = 2
        first_index_x = (((img_j-1)*(orig_dims[1]))+1) 
        new_mri[first_index_x:(first_index_x+orig_dims[1]-1),,,] = mri_fu
      }
      if(n_image_modalities>2){
        img_j = 3
        first_index_x = (((img_j-1)*(orig_dims[1]))+1)
        new_mri[first_index_x:(first_index_x+orig_dims[1]-1),,,] = mri_wm
      }
      if(n_image_modalities>3){
        img_j = 4
        first_index_x = (((img_j-1)*(orig_dims[1]))+1) 
        new_mri[first_index_x:(first_index_x+orig_dims[1]-1),,,] = mri_fu_wm
      }
      mri <- new_mri
      mri_wm <- NULL
      mri_fu_wm <- NULL
      new_mri <- NULL
    }
    ###############################################################################################################################   
    mask_ijk = cbind(which(mask, arr.ind = TRUE),1)
    
    mask_ijk = t(mask_ijk[which(mask_ijk[,4]==1),])
    .print_ok()
    ### UN-MODULATED ###
    .print_action('Training sample: removing effects of covariates\n')
    if (is.null(preloaded_covB) || preloaded_covB == "") {
      # ALEIX, HE TRET D'AQUÍ covX_site perquè NO faci combat, ja que ja l'ha fet abans
      # QUAN FEM EL PAQUET FINAL, SERIA MÉS RÀPID TREURE L'EFECTE DE LES COVARIABLES UN SOL COP, I NO PER CADA TALL.
      list_masked = remove_effects(mask_3d = mask, mri = mri, covX_training = covX_training, trainY = Y, response_family=mp$response_family, n_cores = n_cores, SIGNIFICANT = (use_significant_voxels || mp$modulation=='op'))
    } else {
      list_masked = remove_effects_precalculated(mask = mask, preloaded_covB = preloaded_covB, mri = mri, covX = covX_training, trainY = Y, response_family = mp$response_family, n_cores = n_cores, SIGNIFICANCE = (use_significant_voxels || mp$modulation=='op'), covX_site = covX_site)
    }
    
    mri<-NULL
    
    mask_ijk = list_masked$mask_ijk
    
    mask = list_masked$mask_3d # OJUT!! TESTING
    
    
    .print_ok()
    
    ### END UN-MODULATED ###
    
    
    ### FULLY-MODULATED ###
    if (mp$modulation == 'fu' || mp$modulation == 'op') {
      .print_action('Training sample: Removing effect of covariates on fully modulated data')
      if (preloaded_covB == "")
        list_fu_masked = remove_effects(mask_3d = mask, mri = mri_fu, covX_training = covX_training, trainY = Y, n_cores = n_cores, response_family=mp$response_family)
      else
        list_fu_masked = remove_effects_precalculated(mask = mask, preloaded_covB = preloaded_covB_fu, mri = mri_fu, covX = covX_training, trainY = Y, response_family = mp$response_family, n_cores = n_cores, SIGNIFICANCE = (use_significant_voxels || mp$modulation=='op'))
      
      
      
      .print_ok()
    }
    mri_fu <- NULL
    ### END FULLY-MODULATED ###
    
    ### OPTMOD ###
    if (mp$modulation == 'op') {
      # save images, paths and call optimal_modulation script, load the resulting img_kappa
      
      .print_action('Training sample: optimal modulation')
      list_op_masked = tipett_modulation(mask_X = list_masked$X, t_or_z_vals = list_fu_masked$t_or_z_vals, mask_fu_X = list_fu_masked$X, t_or_z_fu_vals=list_fu_masked$t_or_z_vals, covX_training = covX_training, trainY = Y, response_family = mp$response_family, SIGNIFICANCE = use_significant_voxels)
      # code optimal modulation call
      #list_op_masked = optimal_modulation(mask_X = mask_X, mask_fu_X = mask_fu_X, mask_ijk = mask_ijk, covX_training = covX_training)
      
      #mask_op_signif_indx = list_op_masked$mask_op_signif_indx
      tipett_take_un = list_op_masked$tipett_take_un
      #####################################################################################################
      .print_ok()
    }
    ### FI OPTMOD ###
    
    if (mp$modulation == 'un') {
      masked_X = list_masked$X
      mask_X = ''
      tvals = list_masked$t_or_z_vals
      masked_covB = list_masked$covB
      masked_covB_fu = c()
    } else if (mp$modulation == 'fu') {
      masked_X = list_fu_masked$X
      mask_fu_X = ""
      tvals = list_fu_masked$t_or_z_vals
      masked_covB = list_fu_masked$covB
      masked_covB_fu = list_fu_masked$covB
    } else if (mp$modulation == 'op') {
      masked_X = list_op_masked$mask_op_X
      mask_op_X = ""
      tvals = list_op_masked$t_or_z_op_vals
      masked_covB = list_masked$covB
      masked_covB_fu = list_fu_masked$covB
    } else {
      masked_X = list_masked$X
      mask_X = ''
      tvals = list_masked$t_or_z_vals
      masked_covB = list_masked$covB
      masked_covB_fu = c()
    }
    
    signif_indx = c()
    if (use_significant_voxels) {
      signif_indx = abs(tvals) > SIGNIFICANCE_THRESHOLD # 
      if(sum(signif_indx) < 2)
        signif_indx = which(abs(tvals)==sort(abs(tvals),decreasing = TRUE)[1:2])
      mask_ijk = mask_ijk[, signif_indx]
      if (ncol(as.matrix(masked_covB))==1) {
        masked_covB = masked_covB[signif_indx]
      } else {
        masked_covB = masked_covB[,signif_indx]
      }
      if (mp$modulation == 'op' || mp$modulation == 'fu') {
        if (ncol(as.matrix(masked_covB_fu))==1){
          masked_covB_fu = masked_covB_fu[signif_indx]
        } else {
          masked_covB_fu = masked_covB_fu[,signif_indx]
        }
      }
      masked_X = masked_X[,signif_indx]
    }
  }
  n_voxels_mask = ifelse(is.null(masked_X),0,ncol(masked_X))
  
  ### END MRI VOXELS
  ### PREDICTOR VARIABLES ###
  scale_clinical = list()
  scale_clinical$MODE_SD = NA
  scale_clinical$MEANS = NA
  scale_clinical$SDS = NA
  
  if(!is.null(mp$pred_transf)) {
    if(!is.null(masked_X)){
      # ESCALAR AQUÍ
      scale_clinical = list(MODE_SD=c(),MEANS=c(),SDS=c())
      # Estimate the mode SD of voxels after discarding voxels with SD near zero
      DF = nrow(masked_X) - 1
      SD = apply(masked_X, 2, function (x) {sqrt(sum(x^2) / DF)}) # Voxel-wise standard deviation, after removing the effects of age and sex
      HIST = hist(SD, plot = FALSE, breaks =  100) # Histogram of SD - there is a mixture of ~constant 0 and a normal with mean = ~0.06
      LOWESS = lowess(HIST$mids, HIST$density) # Lowess to remove irregularities that would prevent us from finding the peak of the normal distribution
      scale_clinical$MODE_SD = LOWESS$x[which(LOWESS$y == max(LOWESS$y[-1:-10]))]
      
      # Scale clinical variables (e.g. covX) in the training dataset
      scale_clinical$MEANS = apply(predX_train, 2, mean) # Calculate the mean of each covariate / as.matrix used because errors when one single predictor var selected
      predX_train = as.matrix(apply(predX_train, 1, function(x) {x - scale_clinical$MEANS})) # Center the covariates
      if(length(scale_clinical$MEANS)>1)
        predX_train = t(predX_train)
      #scale_clinical$SDS = apply(predX_train, 2, sd) # Calculate the sd of each covariate
      #predX_train = apply(predX_train, 1, function(x) {x / scale_clinical$SDS * scale_clinical$MODE_SD}) # Scale the covariates to MODE_SD
      scale_clinical$SDS = apply(predX_train, 2, sd) # Calculate the sd of each covariate
      predX_train = as.matrix(apply(predX_train, 1, function(x) {x / scale_clinical$SDS * scale_clinical$MODE_SD})) # Scale the covariates to MODE_SD
      
      if(length(scale_clinical$MEANS)>1)
        predX_train = t(predX_train)
      # replace NaNs by 0
      columnsNaN <- which(apply(predX_train,2,function(x) {all(is.na(x))}))
      predX_train[,columnsNaN]<-0
      masked_X = cbind(masked_X, predX_train)
    } else {
      # replace NaNs by 0
      columnsNaN <- which(apply(predX_train,2,function(x) {all(is.na(x))}))
      predX_train[,columnsNaN]<-0
      masked_X = predX_train
    }
    
  }
  ### END PREDICTOR ###
  .print_ok()
  ### TRAINING ###
  .print_action("Training sample: fitting lasso regression")
  
  lasso = glmnet_fit(
    x= masked_X,
    y = Y,
    family = mp$response_family,
    foldid = internal_folds,
    nfolds = length(unique(internal_folds)),
    standardize = FALSE
  )
  lassoB0 = lasso$a0
  
  
  cat('\nVariables used in the model: \n')
  voxels_used = lasso$i <= n_voxels_mask
  predictors_used = lasso$i > n_voxels_mask
  
  if (!all(voxels_used == FALSE))
    cat(sprintf('\nNumber of voxels: %s; Lasso index (voxel variable): %s', n_voxels_mask, lasso$i[voxels_used]))
  if (!all(predictors_used==FALSE))
    cat(sprintf('\nNumber of voxels: %s; Lasso index: %s; Predictor variable: %s', n_voxels_mask, lasso$i[predictors_used], names(predictors_used[predictors_used])))
  
  cat('\n')
  if (length(lasso$beta)==1) {
    lassoB = t(lasso$beta)
  } else {
    lassoB = lasso$beta
  }
  lasso_covB = c()
  lasso_covB_fu = c()
  lasso_ijk = c()
  
  if (!is.null(mp$mri_paths) && mp$mri_paths != "") {
    if (ncol(as.matrix(masked_covB))==1) {
      lasso_covB = as.matrix(masked_covB)[lasso$i[which(lasso$i<=n_voxels_mask)]] # save covariate betas only the ones selected by lasso
      if (mp$modulation %in% c('fu','op'))
        lasso_covB_fu = as.matrix(masked_covB_fu)[lasso$i[which(lasso$i<=n_voxels_mask)]]
    } else {
      lasso_covB = as.matrix(masked_covB)[, lasso$i[which(lasso$i<=n_voxels_mask)]] # save covariate betas only the ones selected by lasso
      if (mp$modulation %in% c('fu','op'))
        lasso_covB_fu = as.matrix(masked_covB_fu[, lasso$i[which(lasso$i<=n_voxels_mask)]])
    }
    lasso_ijk = matrix(mask_ijk[, lasso$i[which(lasso$i<=n_voxels_mask)]], nrow=4)[1:3,]
    
    lasso_predX_indx = lasso$i[which(lasso$i>n_voxels_mask)] - n_voxels_mask
  } else {
    lasso_predX_indx = lasso$i
  }
  
  # exportar: 
  #      data_table_imputed_train
  #      signif_indx # index of significant variables
  #      mask_covB # betas of the voxels
  #      mask_fu_covB # betas of the fully modulated voxels
  #      img_kappa ??
  #      lasso_predX_indx # index of significant predictor variables
  #      lasso
  #      masked_X_sign # significant voxels
  
  # DEBUG:
  #  linPred = lasso$a0 + as.matrix(masked_X[,lasso$i]) %*% lasso$beta
  #  pred = 1 / (1 + exp(-linPred))
  #  .metrics_binary(Y, pred > (sum(Y==1) / length(Y)))$bac
  
  list(data_table_imputed_train=data_informative_table, 
       signif_indx=signif_indx, 
       lasso_ijk = lasso_ijk,
       lasso_covB= lasso_covB, 
       lasso_fu_covB=lasso_covB_fu,
       mask = mask,
       lasso=lasso, 
       lasso_predX_indx=lasso_predX_indx, 
       scale_predictors = scale_clinical,
       scale_mri = scale_mri,
       #masked_X_sign=masked_X,
       tipett_take_un=tipett_take_un,
       take_significant = use_significant_voxels, # index tipett corresponds to index where X_un abs(t_value) is greater than X_fu calculated from significant over mask
       masks_3d = masks_3d,
       n_voxels_mask= n_voxels_mask,
       combat = combat
  )
  
}

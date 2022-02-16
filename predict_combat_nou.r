apply_model_combatindependent = function(mp, mri_data, mri_fu_data, mri_wm_data=NULL, mri_wm_fu_data=NULL, covX_test, signif_indx, lasso_covB, lasso_covB_fu = NULL, mask, predX_test, scale_clinical, lasso, lasso_predX_indx, tipett_take_un=NULL,img_kappa=NULL,use_significant_voxels=FALSE, covX_site=c(),
                       masks_3d, name_combat='',n_voxels_mask, combat=NULL, n_cores = 1 # ALEIX, ARA LA NECESSITEM PER COMBAT AMB GM+WM, NO SÉ SI LA NECESSITAREM DESPRES
) {
  X = NULL
  X_signif = c()
  if (mp$modulation=='op')
    TIPETT=TRUE
  else
    TIPETT=FALSE
  
  # remove effect of covariates 
  # passar mri_test i mri_fu_test a 2D
  # ALEIX, ARA FEM QUE NOMÉS APLIQUI EL COMBAT EN GM+WM, PERÒ EN REALITAT SERIA SEMPRE
  # 
  if (!is.null(n_voxels_mask) && any(lasso$i <= n_voxels_mask)) { # if the model contains at least one voxel
    if (!is.null(mri_data) ){
      
      n_subj = nrow(covX_test)
      mri_test = lapply(seq_len(n_subj),function(i){mri_data[,,,i][which(masks_3d$un_gm)]})
      mri_test = matrix(unlist(mri_test), nrow = n_subj, byrow = TRUE)
      # Apply combat
      if(!is.null(covX_site)) {
        combat = combat_fit(dat = mri_test, batch = covX_site, mod = covX_test, verbose = FALSE)
        mri_test = combat_apply(combat, mri_test, covX_site, mod=covX_test, verbose = FALSE)$dat.combat # ALEIX, AIXÒ ÉS SI NO PRECALCULAT. PER PRECALCULAT, NO S'HA DE POSAR "mod"
      }
      # FI combat
      img_3d = mri_data[,,,1]
      img_3d = img_3d * 0
      for (img_i in 1:dim(mri_data)[4]){ #for every image
        img_3d[which(masks_3d$un_gm)] = mri_test[img_i,]
        mri_data[,,,img_i] = img_3d
      }
    }
    if (!is.null(mri_fu_data)) {
      mri_test = lapply(seq_len(n_subj),function(i){mri_fu_data[,,,i][which(masks_3d$fu_gm)]})
      mri_test = matrix(unlist(mri_test), nrow = n_subj, byrow = TRUE)
      if(!is.null(covX_site)) {
        combat = combat_fit(dat = mri_test, batch = covX_site, mod = covX_test, verbose = FALSE)
        mri_test = combat_apply(combat, mri_test, covX_site, mod=covX_test, verbose = FALSE)$dat.combat # ALEIX, AIXÒ ÉS SI NO PRECALCULAT. PER PRECALCULAT, NO S'HA DE POSAR "mod"
      }
      img_3d = mri_fu_data[,,,1]
      img_3d = img_3d * 0
      for (img_i in 1:dim(mri_fu_data)[4]){ #for every image
        img_3d[which(masks_3d$fu_gm)] = mri_test[img_i,]
        mri_fu_data[,,,img_i] = img_3d
      }
    }
    if (!is.null(mri_wm_data)) {
      mri_test = lapply(seq_len(n_subj),function(i){mri_wm_data[,,,i][which(masks_3d$un_wm)]})
      mri_test = matrix(unlist(mri_test), nrow = n_subj, byrow = TRUE)
      if(!is.null(covX_site)) {
        combat = combat_fit(dat = mri_test, batch = covX_site, mod = covX_test, verbose = FALSE)
        mri_test = combat_apply(combat, mri_test, covX_site, mod=covX_test, verbose = FALSE)$dat.combat # ALEIX, AIXÒ ÉS SI NO PRECALCULAT. PER PRECALCULAT, NO S'HA DE POSAR "mod"
      }
      img_3d = mri_wm_data[,,,1]
      img_3d = img_3d * 0
      for (img_i in 1:dim(mri_wm_data)[4]){ #for every image
        img_3d[which(masks_3d$un_wm)] = mri_test[img_i,]
        mri_wm_data[,,,img_i] = img_3d
      }
    }
    if (!is.null(mri_wm_fu_data)) {
      mri_test = lapply(seq_len(n_subj),function(i){mri_wm_fu_data[,,,i][which(masks_3d$fu_wm)]})
      mri_test = matrix(unlist(mri_test), nrow = n_subj, byrow = TRUE)
      if(!is.null(covX_site)) {
        combat = combat_fit(dat = mri_test, batch = covX_site, mod = covX_test, verbose = FALSE)
        mri_test = combat_apply(combat, mri_test, covX_site,mod=covX_test, verbose = FALSE)$dat.combat # ALEIX, AIXÒ ÉS SI NO PRECALCULAT. PER PRECALCULAT, NO S'HA DE POSAR "mod"
      }
      img_3d = mri_wm_fu_data[,,,1]
      img_3d = img_3d * 0
      for (img_i in 1:dim(mri_wm_fu_data)[4]){ #for every image
        img_3d[which(masks_3d$fu_wm)] = mri_test[img_i,]
        mri_wm_fu_data[,,,img_i] = img_3d
      }
    }
    
    if (!is.null(mri_data)) {
      n_subj = dim(mri_data)[4]
      
      #################################################################################################################
      if (!is.null(mri_wm_fu_data) ){
        # orig_dims = dim(mask)
        # new_dims = orig_dims * c(4,1,1)
        # new_mask = array(NA, dim = new_dims)
        # for (img_i in 1:4) { # for every image (wm, gm, ...)
        #   first_index_x = (((img_i-1)*(orig_dims[1]))+1)
        #   new_mask[first_index_x:(first_index_x+orig_dims[1]-1),,] = mask
        # }
        # mask = new_mask
        
        orig_dims = dim(mri_data)
        new_mri = array(NA, dim = dim(mri_data) * c(4,1,1,1))
        for (img_i in 1:dim(mri_data)[4]){ #for every image
          img_j = 1
          first_index_x = (((img_j-1)*(orig_dims[1]))+1)
          new_mri[first_index_x:(first_index_x+orig_dims[1]-1),,,img_i] = mri_data[,,,img_i]
          img_j = 2
          first_index_x = (((img_j-1)*(orig_dims[1]))+1) 
          new_mri[first_index_x:(first_index_x+orig_dims[1]-1),,,img_i] = mri_fu_data[,,,img_i]
          img_j = 3
          first_index_x = (((img_j-1)*(orig_dims[1]))+1) 
          new_mri[first_index_x:(first_index_x+orig_dims[1]-1),,,img_i] = mri_wm_data[,,,img_i]
          img_j = 4
          first_index_x = (((img_j-1)*(orig_dims[1]))+1) 
          new_mri[first_index_x:(first_index_x+orig_dims[1]-1),,,img_i] = mri_wm_fu_data[,,,img_i]
        }
        mri_data = new_mri
      }
      #################################################################################################################
      
      mri_test = lapply(seq_len(n_subj),function(i){mri_data[,,,i][which(mask)]})
      mri_test = matrix(unlist(mri_test), nrow = n_subj, byrow = TRUE)
      
      n_voxels = sum(mask)
      if(use_significant_voxels) {
        mri_test = mri_test[,signif_indx] # signif_indx respect mask
        n_voxels = sum(signif_indx)
      }
      mri_test = mri_test[, lasso$i[which(lasso$i <= n_voxels)]] # lasso$i respect signif_indx
      
      if (mp$modulation %in% c('fu','op')) {
        mri_fu_test = lapply(seq_len(n_subj),function(i){mri_fu_data[,,,i][which(mask)]})
        mri_fu_test = matrix(unlist(mri_fu_test), nrow = n_subj, byrow = TRUE)
        if(use_significant_voxels) {
          mri_fu_test = mri_fu_test[,signif_indx] # afegit 
        }
        mri_fu_test = mri_fu_test[, lasso$i[which(lasso$i <= n_voxels)]]
      }
      # 2. remove effect sex and age
      #X_un = abs(mri_test - covX_test %*% mask_covB)
      
      if(covX_test != "") {
        X_un = mri_test - covX_test %*% lasso_covB # X_un will be lasso$i indices over significant voxels
      } else {
        X_un = mri_test
      }
      #X_un = cbind(X_un, abs(X_un))
      #X_un = X_un[,lasso$i]
      if (!is.null(lasso_covB_fu)) {
        if (covX_test != "") {
          X_fu = mri_fu_test - covX_test %*% lasso_covB_fu
        } else {
          X_fu = mri_fu_test
        }
      } else {
        X_fu = NULL
      }
      # X_un = matrix(NA, nrow(covX_test), ncol(mask_ijk))
      # X_fu = X_un
      # X_op = X_un
      ########################################################
      
      if (TIPETT) {
        #tipett_take_un = tipett_take_un[which(mask)] # NOT NECESSARY IF THERE IS NO MASK !!!
        if(use_significant_voxels) {
          tipett_take_un = tipett_take_un[which(signif_indx)]
        }
        take_X_un_val = tipett_take_un[lasso$i] # select the significant tippet tvals of the mask
        X_op = matrix(0, nrow = nrow(X_un), ncol = ncol(X_un))
        X_op[,  take_X_un_val] = X_un[,take_X_un_val] # which tipett abs(tval_un)>abs(tval_fu)
        X_op[, !take_X_un_val] = X_fu[, !take_X_un_val]
        
      } else if (mp$modulation=='op') {
        for (i in 1:ncol(lasso_covB)) {
          #signif_i = signif_indx[i]
          # if (mp$modulation=='un' || mp$modulation == 'op') {
          #   if (ncol(mask_covB)==1)
          #     covB_temp = mask_covB[i]
          #   else
          #     covB_temp = mask_covB[,i]
          #   
          #   X_un[,i] = list_mri_test[[i]] - as.matrix(covX_test) %*% covB_temp
          # }
          # if (mp$modulation=='fu' || mp$modulation == 'op') {
          #   if (ncol(as.matrix(mask_covB))==1)
          #     covB_temp = mask_covB_fu[i]
          #   else
          #     covB_temp = mask_covB_fu[,i]
          #   
          #   X_fu[,i] = list_mri_fu_test[[i]] - as.matrix(covX_test) %*% covB_temp
          # }
          
          X_op[,i] = (1 - img_kappa[i]) * X_un[,i] + img_kappa[i] * X_fu[,i]
        }
      }
      
      
      X_signif = switch(mp$modulation,
                        un=X_un,
                        fu=X_fu,
                        op=X_op,
                        all=X_un)
      #X_signif = X_signif[,lasso$i[which(lasso$i<=ncol())]]
    }
    ### TEST: PREDICTOR VARIABLES ###
  }
  if(!is.null(mp$pred_transf)) {
    predX_test = as.matrix(predX_test)
    #predX_test = predX_test[,lasso_predX_indx]
    # Scale clinical variables in the test dataset
    if(!is.na(scale_clinical$MEANS)){
      predX_test = t(apply(predX_test, 1, function(x) {x - scale_clinical$MEANS})) # Center the covariates
      predX_test = t(apply(predX_test, 1, function(x) {x / scale_clinical$SDS * scale_clinical$MODE_SD})) # Scale the covariates to MODE_SD
    }
    predX_test = predX_test[,lasso_predX_indx]
    X = cbind(X_signif, predX_test)
  } else {
    X = X_signif
  }
  ### TEST: END PREDICTOR ###
  linPred = as.matrix(X) %*% lasso$beta
  if (mp$response_family != 'cox') {
    linPred = lasso$a0 + linPred
  }
  #.print_metrics(linPred, mp$response_family, .find_best_time(trainY[,1], trainY[,2]), testY, fold, test)
  
  linPred
}

remove_effects_precalculated = function(mask_3d, preloaded_covB, mri, covX_training="", trainY="", response_family="", SIGNIFICANCE_THRESHOLD=qnorm(0.975), SIGNIFICANCE=FALSE, n_cores = 1, covX_site = c()) {
  # if precalculated effects are loaded, then this function is run
  # S'HA DE CANVIAR EL NOM DE LA FUNCIÓ, AQUESTA CORRE QUAN HI HA EFECTES PRECARREGATS, remove_effects_and_mask ARA S'EXECUTA QUAN NO N'HI HA, PERÒ NO FA MASCARA, TAMBÉ S'HA DE CANVIAR EL NOM
  if (response_family == "gaussian") {
    SIGNIFICANCE_THRESHOLD = qt(0.975, nrow(covX_training) - ncol(covX_training))
  }
  trainY = switch(response_family,
                  "cox" = Surv(trainY[, 1], trainY[, 2]),
                  trainY)
  
  a = Sys.time()
  cat('Removing effects...\n')
  n_subj = length(mri[1,1,1,])
  n_covs = length(preloaded_covB$data[1,1,1,])
  
  # Mask constant voxels
  # ALEIX: AIXO TAMBÉ A LA FUNCIO QUE NO USA EDAT SEXE PRECALCULATS !!!
  #mask_3d = mask_3d & (apply(mri, 1:3, var) > 0)
  
  # CALCULATE RESIDUALS
  mri_matrix  = lapply(seq_len(n_subj),function(i){mri[,,,i][which(mask_3d)]})
  mri_matrix  = matrix(unlist(mri_matrix), nrow = n_subj, byrow = TRUE)
  covB_matrix = lapply(seq_len(n_covs),function(i){preloaded_covB$data[,,,i][which(mask_3d)]})
  covB_matrix = matrix(unlist(covB_matrix) ,ncol = length(covB_matrix[[1]]), byrow = TRUE)
  
  X = mri_matrix - covX_training %*% covB_matrix # X is masked
  
  
  combat = c('')
  if (!is.null(covX_site)) {
    ######################################################################################################
    ##################################     COMBAT    #####################################################
    ######################################################################################################
    # COMBAT
    # ALEIX: he posat que es fa servir la variable "covX_site", que és un factor amb els nivells sempre iguals
    ###################################################
    combat = combat_fit(dat = X, batch = covX_site, verbose = FALSE) ### ALEIX, AQUí NO POSAR PRECALCULATS !!!
    X = combat_apply(tmp = combat, dat = X, batch = covX_site, verbose = FALSE)$dat.combat
    
    # ALEIX: Això ho podem posar a l'estructura general MRIPredict
    save(combat, file = "combat.Rdata")
    
  }    
  
  ######################################################################################################
  ######################################################################################################
  ######################################################################################################
  
  
  n_voxels = dim(mri_matrix)[2] # number of voxels
  message('Running on ',n_cores, ' thread(s).')
  if (n_cores>1) {
    time1=Sys.time()
    
    if (SIGNIFICANCE){
      
      pb <- txtProgressBar(max=n_voxels, style=3)
      progress <- function(n) setTxtProgressBar(pb, n)
      opts <- list(progress=progress)
      
      cl <- makeCluster(n_cores, setup_strategy = "sequential", timeout=0.5)
      registerDoParallel(cl)
      
      sig <- foreach (i=1:n_voxels, .combine='c', .options.snow=opts) %dopar% {
        if(length(unique(X[,i]))>1){
          if (response_family == "cox") {
            sig= summary(coxph(trainY ~ X[,i]))$coefficients[1, 4] # > SIGNIFICANCE_THRESHOLD
            if (is.na(as.numeric(sig)))
              sig = 0
          } else if (response_family == "binomial") {
            sig = summary(glm(trainY ~ X[,i], family = binomial))$coefficients[2, 3]
          } else {
            sig = summary(lm(trainY ~ X[,i]))$coefficients[2, 3]
          }
        } else {
          sig = 0
        }
        sig
      }
      close(pb)
      stopCluster(cl)
      cat('Time removing effects (run with ', n_cores,' cores):',difftime(Sys.time(),a,units="mins"),"mins.\n")
    } else {
      sig <- rep(NA, n_voxels)
    }
    
  } else {
    
    sig = rep(NA, n_voxels)
    
    if (SIGNIFICANCE) {
      .require('pbapply')
      pbo = pboptions(type="timer")
      time1=Sys.time()
      
      sig=pbapply::pbsapply(1:n_voxels, function(i){
        if (length(unique(X[,i]))>1){
          if (response_family == "cox") {
            sig = summary(coxph(trainY ~ X[,i]))$coefficients[1, 4] # > SIGNIFICANCE_THRESHOLD
            if (is.na(as.numeric(sig)))
              sig = 0
          } else if (response_family == "binomial") {
            sig = summary(glm(trainY ~ X[,i], family = binomial))$coefficients[2, 3]
          } else {
            sig = summary(lm(trainY ~ X[,i]))$coefficients[2, 3]
          }
        } else {
          sig = 0
        }
        sig
      })
    } else {
      sig <- rep(NA, n_voxels)
    }
  }
  
  out_X = X#[, out_signif_indx]
  out_covB = covB_matrix#[,out_signif_indx]
  mask_ijk = which(mask_3d, arr.ind = TRUE)
  out_mask_ijk = mask_ijk#[out_signif_indx,]
  cat('Time removing effects:',difftime(Sys.time(),a,units="mins"),"mins.\n")
  list(X = out_X, covB = out_covB, t_or_z_vals = sig, mask_ijk = t(out_mask_ijk), mask_3d = mask_3d, combat = combat)
}

remove_effects = function(mask_3d, mri, covX_training, trainY=NULL, n_cores = 1, response_family, SIGNIFICANCE_THRESHOLD = qnorm(0.975), SIGNIFICANT = TRUE, covX_site = c(),
                          modality = "", REMOVE_EFFECT_TO_X = TRUE # label to save the combat file
) {
  # AQUESTA CORRE QUAN NO HI HA EFECTES PRECARREGATS, remove_effects s'EXECUTA QUAN N'HI HA, PERÒ NO FA MASCARA, TAMBÉ S'HA DE CANVIAR EL NOM
  # if no precalculated effects, then instead of remove_effects, this function is run
  if (response_family == "gaussian") {
    SIGNIFICANCE_THRESHOLD = qt(0.975, nrow(covX_training) - ncol(covX_training))
  }
  if(length(dim(mri))==4){
    # Mask constant voxels
    mask_constant = mri[,,,1]>-1
    
    mask_constant = mask_constant & (apply(mri, 1:3, var) > 0)
    mask_3d = mask_3d & mask_constant
    mask_ijk = cbind(which(mask_3d, arr.ind = TRUE),1)
    mask_ijk = t(mask_ijk[which(mask_ijk[,4]==1),])
    #~mask_ijk_2 = which(mask_3d, arr.ind = TRUE) # s'han d'agafar index q estiguin a mask_ijk_2 i a mask_ijk
    
    # S'HAN D'AFEGIR LES COVARIABLES QUAN CALCULEM COMBAT
    
    list_mri_training = lapply(seq_len(ncol(mask_ijk)), function(i) {
      mri[as.numeric(mask_ijk[1, i]), as.numeric(mask_ijk[2, i]), as.numeric(mask_ijk[3, i]), ]
    }) # masked mriw
    
  } else if (length(dim(mri)==2)){
    mask_constant = mri[,1]>-1 & (apply(mri, 2, var) > 0)
    mask_3d = mask_3d & mask_constant
    list_mri_training = lapply(seq_len(ncol(mask_3d)), function(i) {
      mri[,i]
    }) # masked mriw
    mask_ijk = NULL
  }
  mri <- NULL
  n_voxels = length(list_mri_training)
  message('n_cores:',n_cores)
  
  #cl <- makeCluster(n_cores)
  #registerDoParallel(cl)
  # find tvalues per voxel
  time.start = Sys.time()
  covX_training = as.matrix(covX_training)
  
  
  ######################################################################################################
  ##################################     COMBAT 2: S'HAURA D'ADAPTAR   #####################################################
  ######################################################################################################
  # COMBAT
  # ALEIX: he posat que es fa servir la variable "covX_site", que és un factor amb els nivells sempre iguals
  # ALEIX: list_mri_training és una llista o una matriu?
  ###################################################
  # ALEIX: Funcions per posar a una altra banda si cal:
  
  
  ######################################################################################################
  ######################################################################################################
  ######################################################################################################
  if(!is.null(trainY)) {
    #mask_all = parLapply(cl=cl,list_mri_training, function(voxel_training, trainY, covX_training, SIGNIFICANT) {
    mask_all = lapply(list_mri_training, function(voxel_training, trainY, covX_training, SIGNIFICANT) {
      # ALEIX, QUAN HI HA COMBAT NO S'HAURIA DE TREURE L'EFECTE ABANS SINÓ DESPRÉS
      if (!is.null(covX_site) || !REMOVE_EFFECT_TO_X) {
        covB = rep(0, ncol(covX_training))
        X = voxel_training
      } else {
        covm = lm.fit(covX_training, voxel_training)
        covB = coefficients(covm)
        X = residuals(covm)
      }
      if (SIGNIFICANT){
        if(length(unique(X))>1){
          if (response_family == "cox") {
            sig = summary(coxph(trainY ~ X))$coefficients[1, 4] # > SIGNIFICANCE_THRESHOLD
            if (is.na(as.numeric(sig))) {
              sig = 0
            }
          } else if (response_family == "binomial") {
            sig = summary(glm(trainY ~ X, family = binomial))$coefficients[2, 3]
          } else {
            sig = summary(lm(trainY ~ X))$coefficients[2, 3]
          }
        }else{
          sig = 0
        }
      } else {
        sig = NA # TO FIX
      }
      c(covB, X, sig)
    }, switch(response_family,
              "cox" = Surv(trainY[, 1], trainY[, 2]),
              trainY
    ), covX_training, SIGNIFICANT)
    #################################
    trainY = switch(response_family,
                    "cox" = Surv(trainY[, 1], trainY[, 2]),
                    trainY)
  }
  cat(sprintf('Time:%s mins.',difftime(time1 = Sys.time(), time2 = time.start, units="mins")))#, 
  
  
  #stopCluster(cl)
  list_mri_training = NULL # free space
  mask_all = matrix(unlist(mask_all), ncol = n_voxels)
  
  mask_covB = mask_all[1:ncol(covX_training),]
  X = mask_all[-c(1:ncol(covX_training), nrow(mask_all)),]
  combat = c('')

  if (!is.null(covX_site)) {
    ######################################################################################################
    ##################################     COMBAT    #####################################################
    ######################################################################################################
    combat = combat_fit(dat = X, batch = covX_site, mod = covX_training, verbose = FALSE)
    X = combat_apply(tmp = combat, dat = X, batch = covX_site, mod = covX_training, verbose = FALSE)
    mask_all = apply(X$dat.combat, 2, function(voxel_training, covX_training) {
      # ALEIX, QUAN HI HA COMBAT S'HA DE TREURE L'EFECTE DESPRÉS
      covB = rep(0, ncol(covX_training))
      X = voxel_training
      c(covB, X, NA)
    }, covX_training)
    mask_covB = mask_all[1:ncol(covX_training),]
    X = mask_all[-c(1:ncol(covX_training), nrow(mask_all)),]
  }
  t_or_z_vals = mask_all[nrow(mask_all),]
  mask_all = ''
  
  mask_signif_indx = which(abs(t_or_z_vals[1:n_voxels]) > SIGNIFICANCE_THRESHOLD)
  
  list(covB = mask_covB, 
       X = X, 
       t_or_z_vals = t_or_z_vals, 
       signif_indx = mask_signif_indx, 
       mask_ijk = mask_ijk,
       mask_3d = mask_3d,
      combat = combat)
}

calculate_effects_controls = function(mri, covX, path) {
  B = solve(t(covX) %*% covX) %*% t(covX)
  beta = apply(mri, 1:3, function (x) {B %*% x}) # ATENCIO: LA IMATGE BETA0 és beta[1,,,], i no beta[,,,1] !!!!!!!
  paths = c()
  for (i in 1:ncol(covX)) {
    mri_example[] = nifti(beta[i,,,])
    writeNIfTI(mri_example, paste(path, "_beta", i, sep = ""))
    paths = c(paths, paste(path, "_beta", i, sep = ""))
  }
  paths
}

load_preloaded_covB = function(path) {
  preloaded_covB = read.csv(file=path)
  ijk = preloaded_covB[,1:3]
  covB = preloaded_covB[,4:ncol(preloaded_covB)]
  list(covB=covB, ijk=ijk)
}

tipett_modulation = function(mask_X, t_or_z_vals, mask_fu_X, t_or_z_fu_vals, mask_ijk, covX_training, trainY, response_family, SIGNIFICANCE=FALSE, SIGNIFICANCE_THRESHOLD=qnorm(0.975)) {
  indx_max_abs_un = abs(t_or_z_vals)>abs(t_or_z_fu_vals) # which abs tvals are bigger in unmodulated data
  t_or_z_op_vals = t_or_z_fu_vals
  t_or_z_op_vals[indx_max_abs_un] = t_or_z_vals[indx_max_abs_un]
  mask_op_X = mask_fu_X
  mask_op_X[,indx_max_abs_un] = mask_X[,indx_max_abs_un]
  # n_voxels = ncol(mask_op_X)
  
  # if (SIGNIFICANCE){
  # cl <- parallel::makeCluster(n_cores, type=ifelse(Sys.info()[[1]]=="Windows",'PSOCK','FORK'))
  # registerDoParallel(cl)
  # t_or_z_op_vals <- foreach (i=1:n_voxels, .combine='c') %dopar% {
  #   if (length(unique(mask_op_X[,i]))>1){
  #     if (response_family == "cox") {
  #       sig= summary(coxph(trainY ~ mask_op_X[,i]))$coefficients[1, 4] # > SIGNIFICANCE_THRESHOLD
  #       if (is.na(as.numeric(sig)))
  #         sig = 0
  #     } else if (response_family == "binomial") {
  #       sig = summary(glm(trainY ~ mask_op_X[,i], family = binomial))$coefficients[2, 3]
  #     } else {
  #       sig = summary(lm(trainY ~ mask_op_X[,i]))$coefficients[2, 3]
  #     }
  #   } else {
  #     sig = 0
  #   }
  #   sig
  # }
  # stopCluster(cl)
  #mask_op_signif_indx = which(abs(t_or_z_op_vals) > SIGNIFICANCE_THRESHOLD)
  #} else {
  # t_or_z_op_vals <- rep(NA, n_voxels)
  #mask_op_signif_indx = seq_len(n_voxels)
  #}
  
  #mask_op_all = matrix(unlist(list_mask_all), ncol = n_voxels)
  #mask_op_covB = mask_op_all[1:ncol(covX_training),]
  #mask_op_X = mask_op_all[-c(1:ncol(covX_training), nrow(mask_all)),]
  #t_or_z_op_vals = mask_op_all[nrow(mask_op_all),]
  
  # mask_op_signif_indx = which(abs(t_or_z_op_vals) > SIGNIFICANCE_THRESHOLD)
  
  list(mask_op_X = mask_op_X, t_or_z_op_vals = t_or_z_op_vals,
       # mask_op_signif_indx = mask_op_signif_indx,
       tipett_take_un = indx_max_abs_un)
  
}

optimal_modulation = function(mask_X, mask_fu_X, mask_ijk, covX_training, trainY, response_family) {
  
  nx = ceiling(sqrt(ncol(mask_X)))
  # mask
  
  tmp_mask = mask_X[1,]
  length(tmp_mask) <- prod(dim(matrix(mask_X[1,], ncol = nx)))
  m_tmp_mask = matrix(!is.na(tmp_mask), ncol = nx, byrow = FALSE)
  nim = nifti(array(m_tmp_mask, dim = c(1, nx, nx)), datatype = 16)
  nim@scl_slope = 1
  writeNIfTI(nim, 'temp_mask', onefile = TRUE, gzipped = TRUE, verbose = FALSE, warn = -1, compression = 6)
  
  path_un = c()
  path_fu = c()
  #save mask_X i mask_fu_X
  for (subj_i in 1:nrow(covX_training)) {
    # save the different masked images into nifti. No recycling values, and replace NA with 0
    subj_x = mask_X[subj_i,]
    length(subj_x) <- prod(dim(matrix(subj_x, ncol = nx)))
    m_subj_x = matrix(subj_x, ncol = nx, byrow = FALSE)
    m_subj_x[is.na(m_subj_x)] = 0
    nim = nifti(array(m_subj_x, dim = c(1, nx, nx)), datatype = 16)
    nim@scl_slope = 1
    
    path_un_tmp = paste(getwd(), "/R_tmp/un_", subj_i, sep = "")
    path_un = c(path_un, paste(path_un_tmp, ".nii.gz", sep = ""))
    writeNIfTI(nim, path_un_tmp, onefile = TRUE, gzipped = TRUE, verbose = FALSE, warn = -1, compression = 6)
    
    subj_x = mask_fu_X[subj_i,]
    length(subj_x) <- prod(dim(matrix(subj_x, ncol = nx)))
    m_subj_x = matrix(subj_x, ncol = nx, byrow = FALSE)
    m_subj_x[is.na(m_subj_x)] = 0
    nim = nifti(array(m_subj_x, dim = c(1, nx, nx)), datatype = 16)
    nim@scl_slope = 1
    
    path_fu_tmp = paste(getwd(), "/R_tmp/fu_", subj_i, sep = "")
    path_fu = c(path_fu, paste(path_fu_tmp, ".nii.gz", sep = ""))
    writeNIfTI(nim, path_fu_tmp, onefile = TRUE, gzipped = TRUE, verbose = FALSE, warn = -1, compression = 6)
  }
  write(path_un, file = "temp_un.txt")
  write(path_fu, file = "temp_fu.txt")
  # call to script optimal_modulation
  
  # crear design_mat.txt amb c(1s, Y, covars)
  if (ncol(as.matrix(covX_training)) == 1)
    design_mat = cbind(1, trainY)
  else
    design_mat = cbind(1, trainY, covX_training[, 2:ncol(covX_training)])
  write.table(design_mat, "design_mat.txt", col.names = FALSE, row.names = FALSE)
  system('../matlab/optimal_modulation_linux64 -d design_mat.txt -i temp_un.txt temp_fu.txt -o temp -m temp_mask.nii.gz')
  img_kappa = readNIfTI('temp_kappa.nii.gz')
  
  img_kappa = as.vector(img_kappa@.Data[1,,])
  
  #mri_op_training
  #####################################################################################################
  #####################################################################################################
  
  ######### pendent d'arreglar      
  list_mri_op_training = lapply(seq_len(ncol(mask_ijk)), function(i) {
    (1 - img_kappa[i]) * mri$data[as.numeric(mask_ijk[1, i]), as.numeric(mask_ijk[2, i]), as.numeric(mask_ijk[3, i]), training] + img_kappa[i] * mri_fu$data[as.numeric(mask_ijk[1, i]), as.numeric(mask_ijk[2, i]), as.numeric(mask_ijk[3, i]), test]
  })
  
  #########
  ###find tvalues for optimal modulated data
  list_mask_all = pbmclapply(list_mri_op_training, function(voxel_training, trainY, covX_training) {
    covm = lm.fit(covX_training, voxel_training)
    covB = coefficients(covm)
    X = residuals(covm)
    if (response_family == "binomial") {
      sig = summary(glm(trainY ~ X, family = binomial()))$coefficients[2, 3]
    } else if (response_family == "cox") {
      surv = Surv(trainY[, 1], trainY[, 2])
      sig = summary(coxph(surv ~ X))$coefficients[4] # > 1.96
    } else {
      sig = summary(lm(trainY ~ X))$coefficients[2, 3]
    }
    c(covB, X, sig)
  }, trainY, covX_training, mc.cores = n_cores)
  list_mri_op_training = "" # free space
  mask_op_all = matrix(unlist(list_mask_all), ncol = n_voxels)
  mask_op_covB = mask_op_all[1:ncol(covX_training),]
  mask_op_X = mask_op_all[-c(1:ncol(covX_training), nrow(mask_all)),]
  t_or_z_op_vals = mask_op_all[nrow(mask_op_all),]
  mask_op_signif_indx = which(abs(t_or_z_op_vals[1:n_voxels]) > SIGNIFICANCE_THRESHOLD)
  
  list(mask_op_covB = mask_op_covB, mask_op_X = mask_op_X, t_or_z_op_vals = t_or_z_op_vals, mask_op_signif_indx = mask_op_signif_indx)
}

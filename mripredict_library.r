#source("glmnet-utils.R")
.check_response_var = function (data_table, response_var, response_settings) {
  response_col = match(response_var, colnames(data_table))
  if (anyNA(response_col)) {
    stop(paste("Response variable", response_var, "not found"))
  }
  data_table[, response_col]
}
.create_covX_old = function (data_table, covX_transf, nrows = nrow(data_table)) {
  covX = matrix(1, nrows)
  colnames_covX = c()
  if (length(covX_transf)>0) {
    for (i in 1:length(covX_transf)) {
      covX_transf_i = covX_transf[[i]]
      data_col = data_table[, match(covX_transf_i[1], colnames(data_table))]
      covX = cbind(covX, as.numeric(switch(
        covX_transf_i[2],
        "factor" = data_col == covX_transf_i[3],
        "numeric" = data_col))
      )
      colnames_covX = c(colnames_covX, switch(covX_transf_i[2], "factor" = sprintf("%s_%s",covX_transf_i[1],covX_transf_i[3]), "numeric"=covX_transf_i[1]))
    }
  }
  colnames(covX)<- c('intercept',colnames_covX)
  covX
}
.create_covX = function (data_table, covX_transf, nrows = nrow(data_table)) {
  covX = matrix(1, nrows)
  colnames_covX = c()
  if (length(covX_transf)>0) {
    for (i in 1:length(covX_transf)) {
      covX_transf_i = covX_transf[[i]]
      data_col = data_table[, match(covX_transf_i[1], colnames(data_table))]
      covX = cbind(covX, as.numeric(switch(
        covX_transf_i[2],
        "factor" = data_col == covX_transf_i[3],
        "numeric" = data_col))
      )
      colnames_covX = c(colnames_covX, switch(covX_transf_i[2], "factor" = sprintf("%s_%s",covX_transf_i[1],covX_transf_i[3]), "numeric"=covX_transf_i[1]))
    }
  }
  colnames(covX)<- c('intercept',colnames_covX)
  covX
}
.add_predictor = function (mri_data, data_table, pred_transf) {
  out = mri_data
  for (i in 1:length(pred_transf)) {
    pred_transf_i = pred_transf[[i]]
    data_col = data_table[, match(pred_transf_i[1], colnames(data_table))]
    out[[length(out)+1]] = as.numeric(switch(
      pred_transf_i[2],
      "factor" = data_col == pred_transf_i[3],
      "numeric" = data_col))
  }
  out
}
.get_predictor_matrix = function (data_table, pred_transf) {
  out = c()
  for (i in 1:length(pred_transf)) {
    pred_transf_i = pred_transf[[i]]
    data_col = data_table[, match(pred_transf_i[1], colnames(data_table))]
    out = c(out,as.numeric(switch(
      pred_transf_i[2],
      "factor" = data_col == pred_transf_i[3],
      "numeric" = data_col)))
  }
  matrix(out, ncol = length(pred_transf))
}
.create_Y = function (data_table, response_var, response_event) {
  out = c();
  if (length(unique(data_table[,response_var])) > 2){
    out = as.numeric(data_table[, response_var])
  } else {
    out = as.numeric(data_table[, match(response_var, colnames(data_table))] == response_event)
  } 
  out
}
.load_mri = function (mri_paths, mask=c() , space = "MNI") {
  if(!is.null(mask) && mask=="") mask=NULL
  n = length(mri_paths);
  if (!n) {
    stop("No MRI data have been specified")
  }
  mri <- .read_mri(as.character(mri_paths[1]), read_data= 0)
  dim <- mri$dim
  mri$data <- array(numeric(), c(dim,n))
  # if (length(mask)==0) {
  #   
  #   mri$data <- array(numeric(), c(dim,n))
  # } else {
  #   mri$data <- array(numeric(), c(sum(mask),n))
  # }
  mat <- mri$sto.xyz
  mri$n = n
  for (i in 1:n) {
    mri_i = .read_mri(as.character(mri_paths[i]))
    if (length(mri_i$dim) != 3) {
      stop(paste(mri_paths[i], "should be 3D"))
    }
    if (any(mri_i$dim != dim)) {
      stop(paste(mri_paths[i], "should have the same dimensions than other MRI data"))
    }
    if ((mri_i$sform.code != "NIFTI.XFORM.MNI.152") && (space != "NO_CHECK")) {
      stop(paste(mri_paths[i], "should be in MNI space"))
    }
    if (any(mri_i$sto.xyz != mat)) {
      stop(paste(mri_paths[i], "should have the same MNI transformation matrix than other MRI data"))
    }
    if (mri_i$scl.slope != 0) {
      if (length(mask)==0){
        mri$data[,,,i] = mri_i$scl.inter + mri_i$scl.slope * mri_i$data
      } else {
        mri_tmp = mri_i$scl.inter + mri_i$scl.slope * mri_i$data
        mri$data[,,,i] = mri_tmp * mask[,,,1]#mri_tmp[which(mask==1 | mask==TRUE)]
      }
    } else {
      if (length(mask)==0){
        mri$data[,,,i] = mri_i$data
      }else {
        mri_tmp = mri_i$scl.inter + mri_i$scl.slope * mri_i$data
        mri$data[,,,i] = mri_tmp * mask[,,,1]#mri_tmp[which(mask==1 | mask == TRUE)] 
      }
    }
  }
  #mri$data <- ff::as.ff(x=mri$data) # DEBUG FF!
  mri
}
.read_mri = function(mri_paths, read_data = TRUE) {
  mri_tmp <- readNIfTI(mri_paths, read_data = read_data)
  dim <- mri_tmp@dim_[2:(2+mri_tmp@dim_[1]-1)]
  sform_code = ""
  
  # from NIFtI documentation
  #   NIFTI_XFORM_UNKNOWN      0 /! Arbitrary coordinates (Method 1). /
  #   NIFTI_XFORM_SCANNER_ANAT 1 /! Scanner-based anatomical coordinates /
  #   NIFTI_XFORM_ALIGNED_ANAT 2 /! Coordinates aligned to another file's, or to anatomical "truth".            /
  #   NIFTI_XFORM_TALAIRACH    3 /! Coordinates aligned to Talairach-Tournoux Atlas; (0,0,0)=AC, etc. /
  #   NIFTI_XFORM_MNI_152      4 /! MNI 152 normalized coordinates. /
  
  if (mri_tmp@sform_code==4) {
    sform_code = "NIFTI.XFORM.MNI.152"
  }
  
  mri <- list(dim = dim,
              sto.xyz = matrix(c(mri_tmp@srow_x, mri_tmp@srow_y, mri_tmp@srow_z, c(0,0,0,1)), nrow = 4, ncol = 4, byrow = TRUE),
              scl.slope = mri_tmp@scl_slope,
              scl.inter = mri_tmp@scl_inter,
              sform.code = sform_code,
              data = mri_tmp@.Data
  )
  if (sum(mri$sto.xyz)==1){
    mri$sto.ijk <- solve(matrix(c(1,0,0,0,
                                  0,1,0,0,
                                  0,0,1,0,
                                  1,1,1,1), nrow = 4))
  } else {
    mri$sto.ijk <- (matrix(c(1,0,0,0,
                             0,1,0,0,
                             0,0,1,0,
                             1,1,1,1), nrow = 4) %*% solve(mri$sto.xyz))
  }
  mri
}
.read_folds_file = function(folds_file_path) {
  con=file(folds_file_path, open="r")
  lines=readLines(con)
  # number of lines (should be 11, 1 for general distribution and 10 for each internal CV in a 10-fold)
  linn = list()
  long=length(lines)
  for (i in 1:long){
    linn[i] = list(as.numeric(unlist(strsplit(lines[i],split=","))))
  }
  close(con)
  linn
}
.read_covars_file = function(path_file, split_txt = ",") {
  # 1st line should be response, 2nd covariates, 3rd predictors, 4th information variables
  con = file(path_file, open="r")
  linn=list()
  lines = readLines(con)
  for (i in 1:4){
    linn[i] = list(unlist(strsplit(lines[i],split=split_txt)))
    if (is.na(linn[i])) linn[i] = ""
  }
  close(con)
  linn
}
.print_action = function (message) {
  cat(message, "... ", sep = "")
  flush.console()
}
.print_ok = function (prefix="") {
  cat(prefix, "Ok\n", sep = "")
  flush.console()
}
.read_1col_file = function (mri_paths_file) {
  read.table(mri_paths_file, as.is = TRUE)$V1
}
.read_2col_file = function(mri_paths_file) {
  read.table(mri_paths_file, as.is = TRUE)$V2
}
.read_3col_file = function(mri_paths_file) {
  read.table(mri_paths_file, as.is = TRUE)$V3
}
.read_4col_file = function(mri_paths_file) {
  read.table(mri_paths_file, as.is = TRUE)$V4
}
.read_data_table = function (data_table_file) {
  text = gsub(",", "\t", readLines(data_table_file))
  text = read.csv(textConnection(text), sep = "\t", check.names = FALSE)
  as.matrix(text)
}
.require = function (package) {
  if (!require(package, character.only = TRUE, quietly = TRUE)) {
    install.packages(package, repos="http://cran.uk.r-project.org")
    library(package, character.only = TRUE, quietly = TRUE)
  }
}
.check_stop = function(message) {
  stop(paste("\n",message, sep=""))
}
.mask_from_which = function(which_array, mask_dim) {
  mask = array(F, mask_dim)
  mask[which_array] = T
}
.indx_equivalent_columns_from_transf = function(data_table_transf, columns, has_intercept=TRUE) {
  #given transf mat from data_table, the function returns the equivalent indxs that correspond to "columns" selected
  out = c()
  for (i in 1:length(data_table_transf)){
    if(data_table_transf[[i]][2]=='factor'){
      if ( !is.null(data_table_transf) && sprintf("%s_%s",data_table_transf[[i]][1],data_table_transf[[i]][3]) %in% columns ) {
        out = c(out, i)
      }
    }else {
      if ( !is.null(data_table_transf) && data_table_transf[[i]][1] %in% columns ) {
        out = c(out, i)
      }
    }
  }
  if (has_intercept){
    out = out + 1  
  }
  out
}
.most_frequent_variables <- function(model_list, mp = NULL, file=NULL) {

  mnis = list(gm=NULL, gm_mod = NULL, wm = NULL, wm_mod = NULL)
  out = data.frame(mni_or_variable = NULL, ijk=NULL,modality = NULL, beta = NULL)
  # iterate through 
  for(i in seq_len(length(model_list))) {
    tmp <- data.frame(mni_or_variable=NULL, modality=NULL)
    if(length(model_list[[i]]$lasso_mni$gm)>0){
      model_list[[i]]$lasso_mni$gm <- matrix(model_list[[i]]$lasso_mni$gm, nrow = 3)
      gm <- paste(model_list[[i]]$lasso_mni$gm[1,], model_list[[i]]$lasso_mni$gm[2,], model_list[[i]]$lasso_mni$gm[3,], sep="_")
      mnis$gm <- c(mnis$gm, gm)
      tmp <- data.frame(mni_or_variable=gm, modality='gm')
    }
    if(length(model_list[[i]]$lasso_mni$gm_mod)>0){
      model_list[[i]]$lasso_mni$gm_mod <- matrix(model_list[[i]]$lasso_mni$gm_mod, nrow = 3)
      gm_mod <- paste(model_list[[i]]$lasso_mni$gm_mod[1,], model_list[[i]]$lasso_mni$gm_mod[2,], model_list[[i]]$lasso_mni$gm_mod[3,], sep="_")
      mnis$gm_mod <- c(mnis$gm_mod, gm_mod)
      tmp <- rbind(tmp, data.frame(mni_or_variable=gm_mod, modality='gm_mod'))
    }
    if(length(model_list[[i]]$lasso_mni$wm)>0){
      model_list[[i]]$lasso_mni$wm <- matrix(model_list[[i]]$lasso_mni$wm, nrow = 3)
      wm <- paste(model_list[[i]]$lasso_mni$wm[1,], model_list[[i]]$lasso_mni$wm[2,], model_list[[i]]$lasso_mni$wm[3,], sep="_")
      mnis$wm <- c(mnis$wm, wm)
      tmp <- rbind(tmp, data.frame(mni_or_variable=wm, modality='wm'))
    }
    if(length(model_list[[i]]$lasso_mni$wm_mod)>0){
      model_list[[i]]$lasso_mni$wm_mod <- matrix(model_list[[i]]$lasso_mni$wm_mod, nrow = 3)
      wm_mod <- paste(model_list[[i]]$lasso_mni$wm_mod[1,], model_list[[i]]$lasso_mni$wm_mod[2,], model_list[[i]]$lasso_mni$wm_mod[3,], sep="_")
      mnis$wm_mod <- c(mnis$wm_mod, wm_mod)
      tmp <- rbind(tmp, data.frame(mni_or_variable=wm_mod, modality='wm_mod'))
    }
    # mni, modality, betas, fold
    
    if(nrow(tmp)>0)
      tmp <- cbind(tmp,data.frame(betas = model_list[[i]]$lasso$beta[c(1:nrow(tmp))], fold=i))
    
    if(length(model_list[[i]]$lasso_predX_indx)>0) {# when multiple predictors fails
      #predictor_vars = sapply(mp$pred_transf, function(x) x[1]) # take all variable names
      predictor_vars = mp$pred_var
      betas_i = model_list[[i]]$lasso$i > model_list[[i]]$n_voxels_mask
      tmp <- rbind(tmp,data.frame(mni_or_variable=predictor_vars[model_list[[i]]$lasso_predX_indx], modality="predictor_variable", betas = model_list[[i]]$lasso$beta[betas_i], fold=i))
    }
    if(mp$response_family!='cox')
      tmp <- rbind(data.frame(mni_or_variable="intercept",modality="intercept",betas=model_list[[i]]$lasso$a0,fold=i),tmp)
    # sample
    
    out <- rbind(out, tmp)
  }
  if(!is.null(file))
    write.csv(out, file=file)
  out
}
.most_frequent_variables_ijk <- function(model_list, mp = NULL, file=NULL) {
  
  ijks = c()
  out = data.frame(mni_or_variable = NULL, ijk=NULL, modality = NULL, beta = NULL,covB_intercept = NULL, covB_age = NULL, covB_sex = NULL)
  
  for(i in seq_len(length(model_list))) {
    tmp <- data.frame(mni_or_variable=NULL, modality=NULL)
    if(length(model_list[[i]]$lasso_ijk)>0){
      model_list[[i]]$lasso_ijk <- matrix(model_list[[i]]$lasso_ijk, nrow = 3)
      coords <- paste(model_list[[i]]$lasso_ijk[1,], model_list[[i]]$lasso_ijk[2,], model_list[[i]]$lasso_ijk[3,], sep="_")
      ijks <- c(ijks, coords)
      tmp <- data.frame(mni_or_variable=coords, modality='voxel')
    }
    # mni, modality, betas, fold
    
    if(nrow(tmp)>0)
      model_list[[i]]$lasso_covB = as.matrix(model_list[[i]]$lasso_covB)
      tmp <- cbind(tmp, data.frame(betas = model_list[[i]]$lasso$beta[c(1:nrow(tmp))], 
                                   fold=i,
                                   covB_intercept = model_list[[i]]$lasso_covB[1,],
                                   covB_age = model_list[[i]]$lasso_covB[2,],
                                   covB_sex = model_list[[i]]$lasso_covB[3,]))
    
    # predictor variables
    if(length(model_list[[i]]$lasso_predX_indx)>0) {# when multiple predictors fails
      #predictor_vars = sapply(mp$pred_transf, function(x) x[1]) # take all variable names
      predictor_vars = mp$pred_var
      betas_i = model_list[[i]]$lasso$i > model_list[[i]]$n_voxels_mask
      tmp <- rbind(tmp, data.frame(mni_or_variable=predictor_vars[model_list[[i]]$lasso_predX_indx], 
                                   modality="predictor_variable", 
                                   betas = model_list[[i]]$lasso$beta[betas_i], 
                                   fold=i,
                                   covB_intercept = NA,
                                   covB_age = NA,
                                   covB_sex = NA))
    }
    if(mp$response_family!='cox')
      tmp <- rbind(data.frame(mni_or_variable="intercept", modality="intercept", betas=model_list[[i]]$lasso$a0,fold=i), tmp)
    # sample
    
    out <- rbind(out, tmp)
  }
  if(!is.null(file))
    write.csv(out, file=file)
  out
}
.find_best_time = function (time, status) {
  # Function that tries to give a number that splits sample keeping the proportion between
  # (time<x & status == 1, high risk) and (time >= x, low risk)
  best.x = NA;
  best.err2 = Inf;
  for (x in unique(sort(time))) {
    err2 = (sum(time < x & status == 1) - sum(time >= x))^2;
    if (err2 < best.err2) {
      best.x = x;
      best.err2 = err2;
    }
  }
  best.x;
}
.metrics_binary = function (real, predictions, folder="", save=FALSE) {
  classlabels = c(1,0)
  if (length(classlabels) > 2) {
    stop("Only binary classification supported")
  }
  class1 = classlabels[1]
  class2 = classlabels[2]
  tp = sum(predictions == class1 & real == class1)
  fp = sum(predictions == class1 & real == class2)
  fn = sum(predictions == class2 & real == class1)
  tn = sum(predictions == class2 & real == class2)
  acc_class1 = tp/(tp+fn)
  acc_class2 = tn/(tn+fp)
  bac = 0.5*(acc_class1 + acc_class2)
  metric <- data.frame(class1 = class1,
                       class2 = class2,
                       n_class1 = sum(real==class1),
                       n_class2 = sum(real==class2),
                       tp = tp,
                       fp = fp,
                       fn = fn,
                       tn = tn,
                       ppv = tp/(tp+fp),
                       sensitivity = acc_class1,
                       specificity = acc_class2,
                       bac = bac)
  if (save)
    write.csv(metric,sprintf("%s/binary_results_fold%d.csv",folder))
  metric
}
.metrics_cox = function (results, frontier_time, iteration = 1, folder="", save=TRUE) {

  sorted_results = results[order(results$linear_predictor),] # sort by second col, linear predictor
  sorted_results_desc = results[order(results$linear_predictor, decreasing = TRUE),] # sort by second col, linear predictor
  idx_g1 = which(results$time < frontier_time & !is.na(results$linear_predictor) & results$status == 1)
  idx_g2 = which(results$time >= frontier_time & !is.na(results$linear_predictor))
  if (length(idx_g1) > 0 && length(idx_g2) > 0) {
    g1 = cbind(results[idx_g1,],risk=1)
    g2 = cbind(results[idx_g2,],risk=0)
    g1_th_linPred = sorted_results_desc$linear_predictor[dim(g1)[1]]
    g2_th_linPred = sorted_results$linear_predictor[dim(g2)[1]]
    th_linPred = (g1_th_linPred + g2_th_linPred) / 2

    #th_linPred = sorted_results_desc$linear_predictor[ceiling(dim(g2)[1] / 2)] # this 2 can be changed by a parameter to priorize sensitiviy or specificity
    g1_predicted = g1$linear_predictor >= th_linPred # these should have increased risk
    g2_predicted = !(g2$linear_predictor < th_linPred) # reduced risk

    g12 = rbind(g1,g2)
    g12_predicted = c(g1_predicted, g2_predicted)
    perf = .metrics_binary(g12[,ncol(g12)],g12_predicted)
  } else {
    message("Warning: not enough samples with status equal to 1 to calculate the performance of the model.")
    perf = NULL
  }
  results$iteration = iteration
  if (save){
    if (!is.null(results$fold))
      write.csv(results,sprintf("%s_cox_results_fold%d.csv",folder,results$fold[1]), row.names = FALSE)
    else
      write.csv(results,sprintf("%s_cox_results.csv",folder), row.names = FALSE)
  }

  perf
}
ijk2mni <- function(ijk,sto_ijk){
  ijk <- round(ijk)
  mni <- solve(sto_ijk) %*% rbind(matrix(ijk[1:3,], nrow=3), 1)
  as.matrix(mni[1:3,])
}

mni2ijk <- function(mni, sto_ijk){
  coordinate <- round(sto_ijk %*% rbind(matrix(mni[1:3,], nrow=3), 1)) # sto_ijk = solve(sto_xyz)
  as.matrix(coordinate[1:3,])
}

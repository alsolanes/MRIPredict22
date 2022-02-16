# FOLDS ########################################################################
assign.folds = function (y, family = c("binomial", "cox", "gaussian"), site = NULL, nfolds = 10) {

  # Check that family and y are correct ########################################
  if (!(
    (is.vector(family) || is.factor(family))
    && length(family) == 1
    && family %in% c("binomial", "cox", "gaussian")
  )) {
    stop('family must be "binomial", "cox", or "gaussian"')
  }
  n = switch(family,
             "binomial" = {
               if (!(is.vector(y) && all(y %in% 0:1))) {
                 stop('for "binomial", y must be a binary vector')
               }
               length(y)
             },
             "cox" = {
               if (class(y) != "Surv") {
                 stop('for "cox", y must be a "Surv" object')
               }
               nrow(y)
             },
             "gaussian" = {
               if (!(is.vector(y) && is.numeric(y))) {
                 stop('for "gaussian", y must be a numeric vector')
               }
               length(y)
             }
  )
  
  # Check that site is correct #################################################
  if (!(is.null(site) || (
    (is.vector(site) || is.factor(site))
    && length(site) == n
  ))) {
    stop("site must be a vector with the same length as y, or NULL")
  }
  
  # Check that nfolds is correct ###############################################
  if (!(is.vector(nfolds) && is.numeric(nfolds) && length(nfolds) == 1 && nfolds > 0)) {
    stop("nfolds must be a positive number")
  }
  
  # Assign folds ###############################################################
  if (is.null(site)) {
    return(.assign.folds_one_site(y, family, nfolds))
  }
  folds = rep(NA, length(y))
  for (site_i in site) {
    i = which(site == site_i)
    folds[i] = .assign.folds_one_site(y[i], family, nfolds)
  }
  folds
}
.assign.folds_one_site = function (y, family, nfolds) {
  folds = rep(NA, length(y));
  folds = switch(family,
                 "binomial" = {
                   for (i in 0:1) {
                     indexs_i = which(y == i);
                     folds[indexs_i] = .assign.folds_simple(indexs_i, nfolds);
                   }
                   folds
                 },
                 "cox" = {
                   sorted_idx = sort.int(y[,1],index.return=TRUE)$ix # Same as "order" ??
                   y_w_indx = cbind(y, sorted_idx)
                   idx_g0 = y_w_indx[y[,2] == 0, 3];
                   idx_g1 = y_w_indx[y[,2] == 1, 3];
                   n_g0 = length(idx_g0)
                   if (n_g0>0) {
                     for (i in seq(0, n_g0-1, by=nfolds)) {
                       n_i = min(c(nfolds, n_g0-i))
                       indexs_i = idx_g0[i + (1:n_i)]
                       folds[indexs_i] = sample(n_i)
                     }
                   }
                   n_g1 = length(idx_g1)
                   if (n_g1>0) {
                     for (i in seq(0, n_g1-1, by=nfolds)){
                       n_i = min(c(nfolds, n_g1-i))
                       indexs_i = idx_g1[i + (1:n_i)]
                       folds[indexs_i] = sample(n_i)
                     }
                   }
                   y_w_indx = data.frame(time=y_w_indx[,1],status=y_w_indx[,2],sorted_idx=y_w_indx[,3],folds=folds);
                   y_w_indx_original_order = y_w_indx[order(sorted_idx),]#apply(y_w_indx,3);
                   folds = y_w_indx_original_order$folds;
                   folds
                 },
                 "gaussian" = {
                   sorted_idx = order(y)
                   n = length(y)
                   for (i in seq(0, n - 1, by=nfolds)) {
                     n_i = min(c(nfolds, n - i))
                     indexs_i = sorted_idx[i + (1:n_i)]
                     # assign randomly one to each folder from 1 to nfolds
                     folds[indexs_i] = .assign.folds_simple(indexs_i, nfolds)
                   }
                   # If at least 2 y are different from the rest, we can detect constant training samples and solve this problem
                   table_y = table(y)
                   if (length(table_y) > 2 || (length(table_y) == 2 && all(table_y) >= 2)) {
                     fold_with_constant_training_sample = NA
                     for (i in 1:nfolds) {
                       if (var(y[which(folds != i)]) == 0) {
                         fold_with_constant_training_sample = i # There might be only one fold with constant training sample
                       }
                     }
                     if (!is.na(fold_with_constant_training_sample)) {
                       # The repeated y is in all of  folds but the fold with constant training sample
                       index_of_one_repeated_y_from_other_folds = sample(which(folds != fold_with_constant_training_sample), 1)
                       index_of_one_different_y_from_the_fold   = sample(which(y != y[index_of_one_repeated_y_from_other_folds]), 1)
                       tmp = folds[index_of_one_different_y_from_the_fold] # Swap the folds of the indices
                       folds[index_of_one_different_y_from_the_fold] = folds[index_of_one_repeated_y_from_other_folds]
                       folds[index_of_one_repeated_y_from_other_folds] = tmp
                     }
                   }
                   folds
                 }
  )
  # Add a random integer so that if a fold has to be smaller or larger, make it a different fold each time
  1 + (folds + sample(1:nfolds, 1)) %% nfolds
}
.assign.folds_simple = function (y, nfolds) {
  n = length(y)
  if (n < nfolds) {
    folds = sample(1:nfolds, n)
  } else {
    folds = sample(1 + (1:n %% nfolds))
  }
  # Add a random integer so that if a fold has to be smaller or larger, make it a different fold each time
  1 + (folds + sample(1:nfolds, 1)) %% nfolds
}


# DATA.FRAME TO GLMNET.MATRIX ##################################################
data.frame2glmnet.matrix_fit = function (x) {
  
  # Check that x is correct ####################################################
  if (!is.data.frame(x)) {
    stop("x must be a data.frame")
  }
  
  # Estimate the conversion ####################################################
  m = list()
  if (ncol(x) > 0) {
    for (j in 1:ncol(x)) {
      xj_name = colnames(x)[j]
      xj_char = as.character(x[, j])
      xj_not_na_num = suppressWarnings(as.numeric(xj_char[which(!is.na(xj_char))]))
      xj_levels = sort(unique(xj_char))
      if (length(xj_levels) < 2) {
        stop(paste("variable", xj_name, "has no different values"))
      }
      if (any(is.na(xj_not_na_num))) {
        m[[length(m) + 1]] = c(xj_name, "factor", xj_levels)
      } else {
        m[[length(m) + 1]] = c(xj_name, "numeric")
      }
    }
  }
  class(m) = "data.frame2glmnet.matrix_fit"
  m
}
data.frame2glmnet.matrix = function (m, x) {
  # Check that x and m are correct ##########################################
  if (!is.data.frame(x)) {
    stop("x must be a data.frame")
  }
  if (class(m) != "data.frame2glmnet.matrix_fit"){
    stop ('m must be a "data.frame2glmnet.matrix" object')
  } 
  
  # Apply the conversion #######################################################
  xp = NULL
  if (length(m) > 0) {
    for (i in 1:length(m)) {
      transf_i = m[[i]]
      xj = x[, match(transf_i[1], colnames(x))]
      xp = cbind(xp, switch(
        transf_i[2],
        "factor" = {
          if (length(transf_i) == 4) {
            xpj = matrix(as.numeric(xj == transf_i[4]))
            colnames(xpj) = paste0(transf_i[1], ":" , transf_i[4])
          } else {
            xpj = NULL
            for (k in 3:length(transf_i)) {
              xpj = cbind(xpj, as.numeric(xj == transf_i[k]))
            }
            colnames(xpj) = paste0(transf_i[1], ":" , transf_i[3:length(transf_i)])
          }   
          xpj
        },
        "numeric" = {
          xpj = matrix(xj)
          colnames(xpj) = transf_i[1]
          xpj
        }
      ))
    }
  }
  xp
}

# GLMNET #######################################################################
glmnet_fit = function (x, y, family = "binomial", foldid = NULL,
                       nfolds = 10, standardize = TRUE, min.beta = 1e-12) {

  # Check that family and y are correct ########################################
  if (!(
    (is.vector(family) || is.factor(family))
    && length(family) == 1
    && family %in% c("binomial", "cox", "gaussian")
  )) {
    stop('family must be "binomial", "cox", or "gaussian"')
  } 
  if (family == "binomial") {
    if (!(is.vector(y) && all(y %in% 0:1))) {
      stop('for "binomial", y must be a binary vector')
    }
    n = length(y)
  }
  if (family == "cox") {
    if (class(y) != "Surv") {
      stop('for "cox", y must be a Surv object')
    }
    n = nrow(y)
  }
  if (family == "gaussian") {
    if (!(is.vector(y) && is.numeric(y))) {
      stop('for "gaussian", y must be a numeric vector')
    }
    n = length(y)
  }
  
  # Check that x is correct ####################################################
  if (!(is.matrix(x) && nrow(x) == n)) {
    stop("x must be a matrix with the same height as y")
  } 
  
  # Check that foldid is correct ###############################################
  if (!(is.null(foldid)
        || (is.vector(foldid) && is.numeric(foldid) && length(foldid) == n)
  )) {
    stop("foldid must be a numeric vector with the same length as y, or NULL")
  }
  
  # Check that nfolds, standardize, and min.beta are correct ###################
  if (!(is.vector(nfolds) && is.numeric(nfolds) && length(nfolds) == 1 && nfolds > 0)) {
    stop("nfolds must be a positive number")
  }
  if (!(is.vector(standardize) && is.logical(standardize) && length(standardize) == 1)) {
    stop("standardize must be TRUE or FALSE")
  }
  if (!(is.vector(min.beta) && is.numeric(min.beta) && length(min.beta) == 1 && min.beta > 0)) {
    stop("min.beta must be a positive number")
  }
  
  if (ncol(x) == 1) {
    # x has only one column ####################################################
    coef = switch(family,
                  "binomial" = coef(glm(y ~ x, family = binomial)),
                  "cox" = coef(coxph(Surv(time = y[,1], event = y[,2]) ~ x)),
                  "gaussian" = coef(lm(y ~ x))
    )
    if (family == "cox") {
      a0 = NULL
      betas = coef
    } else {
      a0 = coef[1]
      betas = coef[2]
    }
    i = 1
  } else {
    # x has more than one column ###############################################
    type_measure = switch(family,
                          "binomial" = "class",
                          "cox" = "deviance",
                          "gaussian" = "mse"
    )
    if (family == "cox") {
      colnames(y) <- c("time", "status")
    }
    # Cross validation to estimate best lambda
    if (is.null(foldid)) {
      foldid = assign.folds(y, family, nfolds = nfolds)
    }
    if (family == "binomial" && min(table(y)) < 3){
      stop("too few subjects of one group")
    }
    cv = cv.glmnet(x, y, type.measure = type_measure, family = family,
                   foldid = foldid, nfolds = nfolds, standardize = standardize)
    idx_lambda = match(cv$lambda.min, cv$lambda)
    # if lambda.min corresponds to 0 degrees of freedom, we will use the closest lambda which results in at least 1 beta value
    if (cv$glmnet.fit$df[idx_lambda] == 0) {
      idx_lambda = which(cv$glmnet.fit$df > 0)[1]
    }
    glmnet.control(fdev = 0)
    lasso = glmnet(x, y, family, lambda = cv$lambda, standardize = standardize)
    a0 = lasso$a0[idx_lambda]
    betas = lasso$beta[, idx_lambda]
    betas[which(abs(betas) < min.beta)] = 0 # Remove negligible betas
    i = which(betas != 0)
  }
  m = list(family = family, a0 = a0, i = i, beta = betas[i])
  class(m) = "glmnet_fit"
  m
}
glmnet_predict = function (m, x) {
  # Check m and x ##############################################################    
  if (class(m) != "glmnet_fit") {
    stop('m must be a "glmnet_fit" object')
  } 
  if (!is.matrix(x)) {
    stop("x must be a matrix")
  } 
  
  # Predict ####################################################################
  y = matrix(x[, m$i], ncol = length(m$i)) %*% m$beta
  if (m$family != "cox") {
    y = m$a0 + y
  }
  if (m$family == "binomial") {
    y = 1 / (1 + exp(-y))
  }
  y
}

# SELECT #######################################################################
glmnet_select = function (x) {
  if (!(is.list(x) && length(x) > 0 && class(x[[1]]) == "glmnet_fit")) {
    stop('x must be a list of objects of class "glmnet_fit" objects')
  }
  if (length(x) == 1) {
    warning("The list contains only one model")
    return(x[[1]])
  }
  # Create a models matrix
  vars = c()
  for (i in 1:length(x)) {
    vars = c(vars, x[[i]]$i)
  }
  vars = unique(sort(vars))
  models = matrix(0, ncol = length(vars), nrow = length(x))
  for (i in 1:length(x)) {
    models[i, match(x[[i]]$i, vars)] = 1
  }
  # Select the best combination most correlated with the other best combination
  cat("[glmnet_select] - Calculating Dice coefficients...\n")
  dice = matrix(NA, ncol = nrow(models), nrow = nrow(models))
  for (i1 in 1:(nrow(models) - 1)) {
    for (i2 in (i1 + 1):nrow(models)) {
      dice_ij = 2 * sum(models[i1,] * models[i2,]) / (sum(models[i1,]) + sum(models[i2,]))
      dice[i1, i2] = dice_ij
      dice[i2, i1] = dice_ij
    }
  }
  mean_dice = apply(dice, 1, mean, na.rm = TRUE)
  cat("[glmnet_select] - Selecting the model with the highest Dice coefficient...\n")
  selected = which(mean_dice == max(mean_dice))
  # Return the model
  y = x[[selected[1]]]
  if (length(selected) > 1) {
    for (i in 2:length(selected)) {
      y$a0 = c(y$a0, x[[selected[i]]]$a0)
      y$beta = rbind(y$beta, x[[selected[i]]]$beta)
    }
    if (!is.null(y$a0)) { # a0 is null for Cox
      y$a0 = mean(y$a0)
    }
    y$beta = apply(y$beta, 2, mean)
  }
  y
}

# IMPUTE #######################################################################
impute.glmnet.matrix_fit = function (x, n_cores = 1) {
  
  # Check that x is correct ####################################################
  if (!is.matrix(x)) {
    stop("x must be a matrix")
  }
  ###### PARALLEL
  library(doParallel)
  cl<-makeCluster(n_cores)
  registerDoParallel(cl)
  ######################
  # Fit imputation models ###################################################### 
  start.time <- Sys.time()
  X_na = is.na(x)
  m = list()

  cat("[impute.glmnet.matrix_fit] Estimating imputation models...\n")
  pb <- txtProgressBar(min=0, max=ncol(x),style=3)
  m = foreach(j=1:ncol(x), .export=c('assign.folds','.assign.folds_one_site','cv.glmnet','glmnet.control','glmnet','glmnet_predict','glmnet_fit','.assign.folds_simple')) %dopar% {
  
  #for (j in 1:ncol(x)) {
    X_na_j = X_na[, j]; 
    x.x = matrix(x[which(!X_na_j), -j], ncol = ncol(x) - 1)
    x.y = x[which(!X_na_j), j]
    family = ifelse(setequal(x.y, c(0, 1)), "binomial", "gaussian")
    # Find all subsets of complete data
    completeCols = apply(x.x, 1, function (tmp) {
      complete = which(!is.na(tmp))
      ifelse(length(complete) > 0, paste(complete, collapse = ","), NA)
    })
    completeCols = setdiff(unique(completeCols), NA)
    # Find a lasso model for each subset, and its error
    imp.models = list()
    errors = c()
    for (k in seq_len(length(completeCols))) {
      name = paste("col", j, "-set", k, sep = "")
      cols.complete = as.numeric(strsplit(completeCols[k], ",", fixed = TRUE)[[1]])
      rows.complete = which(apply(x.x, 1, function (tmp) {all(cols.complete %in% which(!is.na(tmp)))}))
      x.x.complete = matrix(x.x[rows.complete, cols.complete], ncol = length(cols.complete))
      x.y.complete = x.y[rows.complete]
      if ( length(unique(x.y.complete))>1 && 
           ((family == "binomial" && length(table(x.y.complete, exclude=NULL)) == 2 && 
             min(table(x.y.complete, exclude=NULL)) > 2
           )
           || (family == "gaussian" && 
               (length(table(x.y.complete, exclude=NULL)) > 2 ||
                (length(table(x.y.complete, exclude=NULL)) == 2 && min(table(x.y.complete, exclude = NULL)) > 2))
           )
           )
      ) {
        imp.model = glmnet_fit(x.x.complete, x.y.complete, family)
        imp.model$name = name
        x.y.complete.pred = glmnet_predict(imp.model, x.x.complete) # This must be done BEFORE change imp.model$i
        imp.model$i = cols.complete[imp.model$i]
        error = ifelse(family == "binomial",
                       mean(((x.y.complete.pred > 0.5) != x.y.complete)),
                       sqrt(mean((x.y.complete.pred - x.y.complete)^2))
        )
        names(error) = name
      } else {
        imp.model = NULL
        error = Inf
      }
      
      imp.models[[k]] = imp.model
      errors = c(errors, error)
    }
    setTxtProgressBar(pb, j)
    
    # m[[j]] = list(
    #   family = family,
    #   data = x.y,
    #   imp.models = imp.models[order(errors)],
    #   errors = errors[order(errors)]
    # )
    list(
      family = family,
      data = x.y,
      imp.models = imp.models[order(errors)],
      errors = errors[order(errors)]
    )
  }
  close(pb)
  stopCluster(cl)
  cat("[impute.glmnet.matrix_fit] Running time:", difftime(Sys.time(), start.time, units = "secs"), "\n")
  class(m) = "impute.glmnet.matrix_fit"
  m
}
impute.glmnet.matrix = function (m, x, nimp = 20) {

  # Check that m, x and nimp are correct #######################################
  if (class(m) != "impute.glmnet.matrix_fit"){
    stop ('m must be a "impute.glmnet.matrix_fit" object')
  } 
  if (!is.matrix(x)) {
    stop("x must be a matrix")
  }
  if (!(is.vector(nimp) && is.numeric(nimp) && length(nimp) == 1 && nimp > 0)) {
    stop("nimp must be a positive number")
  }
  
  # Impute #####################################################################
    start.time <- Sys.time()
  X_na = is.na(x)
  x.imp = list()
  for (imp in seq_len(nimp)) {
    x.imp[[imp]] = x
  }
  for (j in seq_len(ncol(x))) {
    X_na_j = X_na[,j]
    mj = m[[j]]
    family = mj$family
    # Apply all lasso models (starting with the smallest error)
    for (i in 1:length(mj$imp.models)) {
      if (any(X_na_j)) {
        imp.model = mj$imp.models[[i]]
        if (!is.null(imp.model)) {
          x.x = matrix(x[,-j], ncol = ncol(x) - 1)
          cols.used = imp.model$i
          rows.to_predict.complete = which(X_na_j & apply(x.x, 1, function (tmp) {all(cols.used %in% which(!is.na(tmp)))}))
          x.x.complete = matrix(x.x[rows.to_predict.complete,], ncol = ncol(x) - 1)
          x.y.to_predict = glmnet_predict(imp.model, x.x.complete)
          for (imp in 1:nimp) {
            if (family == "binomial") {
              Ximp = as.numeric(x.y.to_predict > runif(length(rows.to_predict.complete)))
            } else {
              Ximp = x.y.to_predict + rnorm(length(rows.to_predict.complete), 0, mj$errors[i])
            }
            x.imp[[imp]][rows.to_predict.complete, j] = Ximp
          }
          X_na_j[rows.to_predict.complete][which(!is.na(x.y.to_predict))] = FALSE
        }
      }
    }
    # Fill with random values if there are still missings
    if (any(X_na_j)) {
      for (imp in 1:nimp) {
        Ximp = sample(mj$data, sum(X_na_j), replace = TRUE)
        x.imp[[imp]][which(X_na_j), j] = Ximp
      }
    }
  }
  cat("[impute.glmnet.matrix] Running time:", difftime(Sys.time(), start.time, units = "secs"), "\n")
  x.imp
}

# CV ###########################################################################
cv = function (x, y, family = c("binomial", "cox", "gaussian"), fit_fun, predict_fun, site = NULL, covar = NULL, nfolds = 10, ...) {
  if (!is.null(site)) {
    if (!is.null(covar)) {
      cat("[cv] Cross-validation with sites and covariates\n")
      type = "site+covar"
    } else {
      cat("[cv] Cross-validation with sites\n")
      type = "site"
    }
  } else {
    if (!is.null(covar)) {
      cat("[cv] Cross-validation with covariates\n")
      type = "covar"
    } else {
      cat("[cv] Simple cross-validation\n")
      type = "simple"
    }
  }
  folds = assign.folds(y, family, site = site, nfolds = nfolds)
  models = list()
  output = data.frame(fold = folds, y, y.pred = NA)
  for (fold in 1:nfolds) {
    # Training
    cat("[cv] Fold", fold, " - Training\n")
    training = which(folds != fold)
    x_training = x[training,]
    y_training = switch(family,
                        "cox" = y[training,],
                        y[training]
    )
    model = switch(type,
                   "simple" =     fit_fun(x_training, y_training, NULL,           NULL,             ...),
                   "site" =       fit_fun(x_training, y_training, site[training], NULL,             ...),
                   "covar" =      fit_fun(x_training, y_training, NULL,           covar[training,], ...),
                   "site+covar" = fit_fun(x_training, y_training, site[training], covar[training,], ...)
    )
    models[[fold]] =  model
    # Test
    cat("[cv] Fold", fold, " - Test\n")
    test = which(folds == fold)
    x_test = matrix(x[test,], ncol = ncol(x)) # Ensure matrix if length(test) = 1
    output$y.pred[test] = switch(type,
                                 "simple" =     predict_fun(model, x_test, NULL,       NULL,         ...),
                                 "site" =       predict_fun(model, x_test, site[test], NULL,         ...),
                                 "covar" =      predict_fun(model, x_test, NULL,       covar[test,], ...),
                                 "site+covar" = predict_fun(model, x_test, site[test], covar[test,], ...),
    )
  }
  list(predictions = output, models = models)
}

seleccio_covariables_data.frame2glmnet.matrix_fit = function (x, covariates) {
  j = match(covariates, colnames(x))
  if (any(is.na(j))){
    stop(paste("Covariate", covariates[which(is.na(j))], "not found"))
  } 
  x = x[,j]
  for (j in 1:ncol(x)) {
    x[,j] = factor(x[,j])
  }
  x
}


.create_undummy_table <- function(x, transf) {
  dummy_colnames <- sapply(strsplit(colnames(x), ':'), function (tmp) {tmp[1]})
  unique_dummy_colnames <- unique(dummy_colnames)
  xp <- setNames(data.frame(matrix(data = NA,ncol = length(unique_dummy_colnames), nrow = nrow(x))), unique_dummy_colnames)
  dummy_values <- sapply(strsplit(colnames(x), ':'), function(tmp) { tmp[2] })
  for (i in 1:length(dummy_values)) {
    if (is.na(dummy_values[i])) { # not a dummy, copy directly
      xp[,dummy_colnames[i]] <- x[,i]
    } else {
      xp[which(x[,i]==1),dummy_colnames[i]] <- dummy_values[i]
    }
  }
  for (i in 1:nrow(transf$factor_ref)) { # set reference values
    xp[which(is.na(xp[,transf$factor_ref[i,1]])),transf$factor_ref[i,1]] <- transf$factor_ref[i,2]
  }
  xp
}
.calculate_ensemble_accuracy = function(cv_results, family) {
  acc_per_fold = c()
  for (i in 1:18) {
    fold_results = mp$cv_results[mp$cv_results$iteration==i,]
    if (family == "binomial"){
      print(.metrics_binary(real = fold_results$real, predictions = fold_results$pred>0.5))
    } else if (family == "gaussian"){
      print(sqrt(mean((fold_results$pred-fold_results$real)^2)))
    } else { # COX
      
    }
  }
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
.correc = function(i, n){
  c1 = c(9.5,28.7,1.9,0.,-7.,-6.2,-1.6)
  c2 = c(-6195.,-9569.,-6728.,-17614.,-8278.,-3570., 1075.)
  c3 = c(93380.,175160.,410400.,2157600.,2.376e6,2.065e6,2.065e6)
  
  mic = 1e-6
  c14 = 1.9e-5
  
  if (i*n == 4) {return(c14)}
  if (i<1 || i > 7){return(0)}
  if (i!=4 && n>20){return(0)}
  if (i==4 && n>40){return(0)}
  
  n = 1./n^2
  i = i-1
  out = (c1[i]+n*(c2[i] + n*c3[i]))*mic
}
.normOrder = function (N) {
  
  if (length(N)>1) {
    n = length(N)
  } else {
    n = N
  }
  
  s = N/2
  
  n2 = floor(s)
  
  eps = c(.419885,.450536, .456936, .468488)
  dl1 = c(.112063,.12177,  .239299, .215159)
  dl2 = c(.080122,.111348,-.211867,-.115049)
  gam = c(.474798,.469051, .208597, .259784)
  lam = c(.282765,.304856,.407708,.414093)
  bb = -.283833
  d = -.106136
  b1 = .5641896
  
  s[1] = b1
  
  # calculate normal tail areas for first 3 order statistics
  k = 3
  len_s = k
  if (n2 < k){
    k = n2
  }
  
  s = vector(mode="numeric",length = n2)
  for (i in 1:k){
    e1 = (i - eps[i]) / (n + gam[i])
    e2 = e1^(lam[i])
    s[i] = e1 + e2 * (dl1[i] + e2 * dl2[i]) / n - .correc(i+1, n)
  }
  
  if (n2 > k){
    for (i in 4:n2){
      e1 = (i - eps[4]) / (n + gam[4])
      e2 = e1^(lam[4] + bb / (i+d))
      s[i] = e1 + e2 * (dl1[4] + e2 * dl2[4]) / n - .correc(i+1, n)
    }
  }
  
  for (i in 1:n2) {
    s[i] = -qnorm(s[i],0,1)
  }
  
  if (n%%2 == 0) {
    out = c(-s, s[length(s):1])
  } else {
    out = c(-s, 0,s[length(s):1])
  }
  out
}
.print_metrics = function(linPred, family, frontier_time=0, testY, id) {
  cat("\n")
  pred = c()
  if (family == "binomial") {
    prob = 1 / (1 + exp(-linPred))
    pred = prob
    bac = .metrics_binary(testY,pred>mean(testY))$bac
    results = data.frame(id = id, linear_predictor = pred[,1], real = testY[,1])
    print(results)
  }
  else if (family == "gaussian") {
    pred = as.matrix(linPred)
    
    results = data.frame(id = id, linear_predictor = pred[,1], real = testY[,1])
    print(results)
  }
  else if (family == "cox") {
    pred = linPred
    
    rd2 = .coxph_RD2(pred, testY[,1], testY[,2])
    
    results = data.frame(id = id, linear_predictor = pred[,1], time = testY[,1], status = testY[,2])
    print(results)
    print(.metrics_cox(results, frontier_time, rd2$b, save=FALSE))
  }
}
.mean_RD2 = function(betes) {
  betes_mean = mean(betes)
  D = exp(betes_mean)
  kappa = sqrt(8 / pi)
  sigma2 = pi^2 / 6;
  RD2 = (D^2 / kappa^2) / (sigma2 + D^2 / kappa^2)
  RD2
}
.coxph_RD2 = function(predictor, stime, sevent) {
  surv_obj <- Surv(stime, sevent, type="right")
  coxphfit <- summary(coxph(surv_obj ~ predictor))
  b <- coxphfit$coefficients[1]
  
  PI = sort.int((predictor - mean(predictor)) * b, index.return=TRUE) # calculate linear predictors
  
  rankits = .normOrder(length(PI$x))
  
  kappa = sqrt(8 / pi)
  
  rankits = rankits / kappa
  
  surv_obj <- Surv(stime[PI$ix],sevent[PI$ix],type="right")
  coxphfit <- summary(coxph(surv_obj ~ rankits))
  b <- coxphfit$coefficients[1]
  D = exp(b)
  sigma2 = pi^2 / 6
  
  RD2 = (D^2 / kappa^2) / (sigma2 + D^2 / kappa^2)
  list(RD2=RD2, b=b)
}


# PROBABLEMENT DESCARTADES
# .find_best_time = function (time, status) {
#   # Function that tries to give a number that splits sample keeping the proportion between
#   # (time<x & status == 1, high risk) and (time >= x, low risk)
#   best.x = NA;
#   best.err2 = Inf;
#   for (x in unique(sort(time))) {
#     err2 = (sum(time < x & status == 1) - sum(time >= x))^2;
#     if (err2 < best.err2) {
#       best.x = x;
#       best.err2 = err2;
#     }
#   }
#   best.x;
# }
# .metrics_cox = function (results, frontier_time, iteration = 1, folder="", save=TRUE) {
#   
#   sorted_results = results[order(results$linear_predictor),] # sort by second col, linear predictor
#   sorted_results_desc = results[order(results$linear_predictor, decreasing = TRUE),] # sort by second col, linear predictor
#   idx_g1 = which(results$time < frontier_time & !is.na(results$linear_predictor) & results$status == 1)
#   idx_g2 = which(results$time >= frontier_time & !is.na(results$linear_predictor))
#   if (length(idx_g1) > 0 && length(idx_g2) > 0) {
#     g1 = cbind(results[idx_g1,],risk=1)
#     g2 = cbind(results[idx_g2,],risk=0)
#     g1_th_linPred = sorted_results_desc$linear_predictor[dim(g1)[1]]
#     g2_th_linPred = sorted_results$linear_predictor[dim(g2)[1]]
#     th_linPred = (g1_th_linPred + g2_th_linPred) / 2
#     
#     #th_linPred = sorted_results_desc$linear_predictor[ceiling(dim(g2)[1] / 2)] # this 2 can be changed by a parameter to priorize sensitiviy or specificity
#     g1_predicted = g1$linear_predictor >= th_linPred # these should have increased risk
#     g2_predicted = !(g2$linear_predictor < th_linPred) # reduced risk
#     
#     g12 = rbind(g1,g2)
#     g12_predicted = c(g1_predicted, g2_predicted)
#     perf = .metrics_binary(g12[,ncol(g12)],g12_predicted)
#   } else {
#     message("Warning: not enough samples with status equal to 1 to calculate the performance of the model.")
#     perf = NULL
#   }
#   results$iteration = iteration
#   if (save){
#     if (!is.null(results$fold))
#       write.csv(results,sprintf("%s_cox_results_fold%d.csv",folder,results$fold[1]), row.names = FALSE)
#     else
#       write.csv(results,sprintf("%s_cox_results.csv",folder), row.names = FALSE)
#   }
#   
#   perf
# }

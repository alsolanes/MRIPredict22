mripredict_load = function (xml_file) {

  .require("XML")

  .print_action("Reading the MRIPredict model")
  xml = try(xmlParse(xml_file))
  if (class(xml)[1] == "try-error") {
    stop(paste("Could not open or parse", xml_file))
  }
  xml_mripredict = xmlRoot(xml)
  if (xmlName(xml_mripredict) != "mripredict") {
    stop(paste(xml_file, "is not an MRIPredict file"))
  }

  # DATA ###
  # Get mri_paths and data_table
  xml_subj = xml_mripredict[["data"]]["subject"]

  xml_var = xml_subj[[1]]["variable"]
  data_colnames = c()
  for (j in 1:length(xml_var)) {
    # The order of the variables in the table is not specified in the XML
    data_colnames = c(data_colnames, xmlAttrs(xml_var[[j]])[["name"]])
  }

  mri_paths = c()
  mri_fu_paths = c()
  data_table = matrix(nrow = length(xml_subj), ncol = length(data_colnames))
  colnames(data_table) = data_colnames
  for (i in 1:length(xml_subj)) {
    xml_subj_i = xml_subj[[i]]
    subj_i_id = as.numeric(xmlAttrs(xml_subj_i)[["order"]])

    xml_subj_i_mri = xml_subj_i["mri"]
    for (j in 1:length(xml_subj_i_mri)){
      attrs_subj_i_var_j = xmlAttrs(xml_subj_i_mri[[j]])
      if (attrs_subj_i_var_j[["modality"]] == "un-modulatedGM"){
        mri_paths[subj_i_id] = attrs_subj_i_var_j[["path"]]
      } else if (attrs_subj_i_var_j[["modality"]] == "fully-modulatedGM"){
        mri_fu_paths[subj_i_id] = attrs_subj_i_var_j[["path"]]
      }
    }
    
    xml_subj_i_var = xml_subj_i["variable"]
    for (j in 1:length(xml_subj_i_var)) {
      attrs_subj_i_var_j = xmlAttrs(xml_subj_i_var[[j]])
      data_table[subj_i_id, match(attrs_subj_i_var_j[["name"]], colnames(data_table))] =
                    attrs_subj_i_var_j[["value"]]
    }
    modulation = switch
  }

  # MODEL ###
  # Get response_var, response_ref, response_event, covX_factor_ref, covX_transf, lassoB0, mni, covB and lassoB
  xml_model = xml_mripredict[["model"]]

  attrs_response = xmlAttrs(xml_model[["response"]])
  response_family = attrs_response[["family"]]

  switch(response_family,
         binomial = {
            response_var = attrs_response[["variable"]]
            response_ref = attrs_response[["reference"]]
            response_event = attrs_response[["event"]]
            },
         gaussian = {
           response_var = attrs_response[["variable"]]
           response_ref = attrs_response[["reference"]]
           response_event = attrs_response[["event"]]
         },
         cox = {
           response_ref = ""
           response_event = ""
           response_var = c(attrs_response[["time"]], attrs_response[["status"]])
         })
         
         
  
  
  xml_var = xml_model[["covariates"]]["variable"]
  covX_factor_ref = NULL
  covX_transf = list()
  j2 = 1
  if (length(xml_var)>0){
    for (j in 1:length(xml_var)) {
      xml_var_j = xml_var[[j]]
      attrs_var_j = xmlAttrs(xml_var_j)
      variable_j_name = attrs_var_j[["name"]]
      if (attrs_var_j[["type"]] == "factor") {
        xml_var_j_level = xml_var_j["level"]
        covX_factor_ref = rbind(covX_factor_ref, c(variable_j_name, xmlAttrs(xml_var_j[["reference"]])[["label"]]))
        j2 = j2 + 1
        for (k in 1:length(xml_var_j_level)) {
          attrs_var_j_level_k = xmlAttrs(xml_var_j_level[[k]])
          covX_transf[[as.numeric(attrs_var_j_level_k[["column"]])]] =
            c(variable_j_name, "factor", attrs_var_j_level_k[["label"]])
        }
      } else {
        covX_transf[[as.numeric(attrs_var_j[["column"]])]] = c(variable_j_name, "numeric")
      }
    }
  }


##################  
  

  # CROSSVALIDATION ###
  xml_cv = xml_model[["crossvalidation"]]
  if (!is.null(xml_cv)) {
    xml_fold = xml_cv['fold']
    cv_table = matrix(nrow = length(xml_subj), ncol = 2)
    cv_accuracy = c()
    for (j in 1:length(xml_fold)) {
      xml_fold_j = xml_fold[[j]]
      fold_j_id = as.numeric(xmlAttrs(xml_fold_j)[["order"]])
      xml_fold_j_subj = xml_fold_j['subject']
      for (i in 1:length(xml_fold_j_subj)) {
        attrs_fold_j_subj_i = xmlAttrs(xml_fold_j_subj[[i]])
        cv_table[as.numeric(attrs_fold_j_subj_i[["id"]]), ] =
          c(fold_j_id, as.numeric(attrs_fold_j_subj_i[["prediction"]]))
      }
      cv_accuracy[fold_j_id] = as.numeric(xmlAttrs(xml_fold_j[["summary"]])[["accuracy"]])
    }
  }

  # COEFFICIENTS ###

  xml_coef = xml_model[["coefficients"]]
  lassoB = c()
  mni = c()
  covB = c()
  kappa = c()

  if (!is.null(xml_coef)) {
    modality = xmlAttrs(xml_coef[['mri']])[['modality']]
    lassoB0 = as.numeric(xmlAttrs(xml_coef[["intercept"]][["lasso"]])[["beta"]])
    xml_voxel = xml_coef["voxel"]
    
    mni = matrix(nrow = 3, ncol = length(xml_voxel))
    covB = matrix(nrow = length(lassoB0) + length(covX_transf), ncol = length(xml_voxel))
    covB_fu = matrix(nrow = length(lassoB0) + length(covX_transf), ncol = length(xml_voxel))
    for (j in 1:length(xml_voxel)) {
      xml_voxel_j = xml_voxel[[j]]
      voxel_j_id = as.numeric(xmlAttrs(xml_voxel_j)[["id"]])
      attrs_mni = xmlAttrs(xml_voxel_j[["mni"]])
      mni[, voxel_j_id] =  as.numeric(c(attrs_mni[["x"]], attrs_mni[["y"]], attrs_mni[["z"]]))
      covB[, voxel_j_id] = as.numeric(strsplit(xmlAttrs(xml_voxel_j["covariates"][[1]])[["beta"]], ",")[[1]])

      if (modality != 'un-modulated')
        covB_fu[, voxel_j_id] = as.numeric(strsplit(xmlAttrs(xml_voxel_j["covariates"][[2]])[["beta"]], ",")[[1]])
      lassoB[voxel_j_id] = as.numeric(xmlAttrs(xml_voxel_j[["lasso"]])[["beta"]])
      if (modality == 'optimally-modulated'){
        kappa[voxel_j_id] = as.numeric(xmlAttrs(xml_voxel_j[["optimal_modulation"]])[['kappa']])
      }
    }
    
  }

  
  ### PREDICTORS ###    
  xml_var = xml_model[["predictors"]]["variable"]
  pred_factor_ref = NULL
  pred_transf = list()
  j2 = 1
  if (length(xml_var)>0){
    for (j in 1:length(xml_var)) {
      xml_var_j = xml_var[[j]]
      attrs_var_j = xmlAttrs(xml_var_j)
      variable_j_name = attrs_var_j[["name"]]
      if (attrs_var_j[["type"]] == "factor") {
        xml_var_j_level = xml_var_j["level"]
        pred_factor_ref = rbind(pred_factor_ref, c(variable_j_name, xmlAttrs(xml_var_j[["reference"]])[["label"]]))
        j2 = j2 + 1
        for (k in 1:length(xml_var_j_level)) {
          attrs_var_j_level_k = xmlAttrs(xml_var_j_level[[k]])
          pred_transf[[as.numeric(attrs_var_j_level_k[["column"]])]] =
            c(variable_j_name, "factor", attrs_var_j_level_k[["label"]])
          lassoB = c(lassoB, as.numeric(attrs_var_j_level_k[["lasso_beta"]]))
        }
      } else {
        pred_transf[[as.numeric(attrs_var_j[["column"]])]] = c(variable_j_name, "numeric")
        lassoB = c(lassoB, as.numeric(attrs_var_j[["lasso_beta"]]))
      }
    }
  }
  
  mp = list(
    mri_paths = mri_paths,
    mri_fu_paths = mri_fu_paths,
    data_table = data_table,
    response_var = response_var,
    response_ref = response_ref,
    response_event = response_event,
    response_family = response_family,
    covX_factor_ref = covX_factor_ref,
    covX_transf = covX_transf,
    pred_factor_ref = pred_factor_ref,
    pred_transf = pred_transf,
    mni = mni,
    lassoB = lassoB,
    lassoB0 = lassoB0,
    kappa = kappa,
    covB = covB,
    covB_fu = covB_fu,
    modality = modality
  )
  if (!is.null(xml_cv)){
    mp$cv_table = cv_table
    mp$cv_accuracy = cv_accuracy
  }
  class(mp) = "mripredict"
  
  free(xml)
  .print_ok()

  mp

}
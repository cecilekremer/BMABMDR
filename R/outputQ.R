
#' output function for estimated model parameter from each model
#'
#' @param par_obj list containing parameter draws from either Laplace or MCMC fit
#' @param pvec probability vector to compute credible interval for the BMD.
#' @param max.dose maximum dose tested
#' @param rho should be set to TRUE if data is clustered
#'
#' @return dataframe containing model estimates
#'
#'
#' @export outLPQ
#'
outLPQ <- function(par_obj, pvec, max.dose, rho = FALSE) {
  ret <- t(data.frame(
    # # transformed parameters
    # par.at = quantile(par_obj$p1, c(0.025,0.5,0.975), na.rm = TRUE),
    par.dt = quantile(par_obj$p3, c(0.025,0.5,0.975), na.rm = TRUE),
    par.k = quantile(par_obj$p2, pvec, na.rm = TRUE),
    # parameters on original scale
    par.a = quantile(par_obj$a, c(0.025,0.5,0.975), na.rm = TRUE),
    par.b = quantile(par_obj$b, c(0.025,0.5,0.975), na.rm = TRUE),
    par.d = quantile(par_obj$d, c(0.025,0.5,0.975), na.rm = TRUE),
    # natural parameters
    BMD = quantile(par_obj$BMD, pvec, na.rm = TRUE)*max.dose,
    par.rho = ifelse(rep(rho==TRUE, 3), quantile(par_obj$rho, c(0.025,0.5,0.975), na.rm = TRUE), 0)
    #min.resp = quantile(par_obj$min_response, c(0.025,0.5,0.975), na.rm = TRUE)
  ))

  return(ret)
}


#### function to extract the sampled values for the model parameters and the BMD
#' function to extract the sampled values for the model parameters and the BMD
#'
#' @param mod_obj model object either of class stanfit or stanfitOptim
#' @param pars parameters to be extracted.
#' @param model_name name of the model to be extracted
#' @param rho should be set to TRUE if data is clustered
#'
#' @return dataframe containing model estimates
#'
#' @export parq_extract
#'
parq_extract <- function(mod_obj, pars = c(letters[c(1,2,4)], "BMD",
                                           paste0("par",1:3), "min_response"),
                         model_name, rho = FALSE) {

  if(is.stanfit(mod_obj)){

    pars_samples <- rstan::extract(mod_obj, pars)
    pars_samplesd <- data.frame(ModelName = model_name,
                                a = pars_samples$a,
                                b = pars_samples$b,
                                d = pars_samples$d,
                                BMD = pars_samples$BMD,
                                rho = ifelse(rep(rho==TRUE, length(pars_samples$a)), pars_samples$rho, 0),
                                p1 = pars_samples$par1,
                                p2 = pars_samples$par2,
                                p3 = pars_samples$par3
    )

  } else if(is.stanfitOptim(mod_obj)) {

    pars_samples <- mod_obj$theta_tilde[,pars]
    pars_samplesd <- data.frame(ModelName = model_name,
                                a = pars_samples[,"a"],
                                b = pars_samples[,"b"],
                                d = pars_samples[,"d"],
                                BMD = pars_samples[,"BMD"],
                                rho = ifelse(rep(rho==TRUE, nrow(pars_samples)), pars_samples[,"rho[1]"], 0),
                                p1 = pars_samples[,"par1"],
                                p2 = pars_samples[,"par2"],
                                p3 = pars_samples[,"par3"]
    )

    return(pars_samplesd)

  } else stop("please check if input is of class stanfit or stanfitOptim")



}

# function to turn the vector of weights into dataframe
#' function to extract the sampled values for the model parameters and the BMD
#'
#' @param mod_obj model object either of class stanfit or stanfitOptim
#' @param type type of weight to be extracted. It can be either of BS (bridge sampling), LP (Laplace) and
#'             both to extract bridge sampling and laplace approximation.
#' @param model_name name of the model to be extracted
#'
#' @return dataframe containing model estimates
#'
#' @export weightsQ_extract
#'
weightsQ_extract <- function(mod.obj, type = c("BS", "LP", "both")) {
  type <- match.arg(type)
  if(is.BMADRQ2(mod.obj)[3] == 2 & type == 'BS'){

    return(data.frame(Model = names(mod.obj$weights_bridge_sampling),
                      BS_Weights = mod.obj$weights_bridge_sampling))

  } else if(is.BMADRQ2(mod.obj)[3] == 2 & type == 'LP') {

    return(data.frame(Model = names(mod.obj$weights_bridge_sampling),
                      LP_Weights = mod.obj$weights_laplace))

  } else if(is.BMADRQ2(mod.obj)[3] == 2 & type == "both") {

    return(data.frame(Model = names(mod.obj$weights_bridge_sampling),
                      BS_Weights = mod.obj$weights_bridge_sampling,
                      LP_Weights = mod.obj$weights_laplace))

  } else if(is.BMADRQ2(mod.obj)[2] == 2){

    return(data.frame(Model = names(mod.obj$weights),
                      LP_Weights = mod.obj$weights))

  }else stop('please check mod.obj or provide type')

}

# function to extract mixture
#' function to extract BMD mixture
#' @param mod_obj BMDBMA model object
#' @param conv logical to indicate if mixture should be only based on converged
#'             models or all the models regardless of the convergence status
#' @return dataframe containing Mixture per model
#' @export BMDQmixture_extract
#'
BMDQmixture_extract <- function(mod.obj, conv = FALSE){

  if(is.BMADRQ2(mod.obj)[3] == 2 & conv == TRUE) {
    BMDMixture <- data.frame(Model = "Model Averaged",
                             BMDMixture = mod.obj$BMDMixture.conv
    )
  } else if(is.BMADRQ2(mod.obj)[3] == 2 & conv == FALSE) {
    BMDMixture <- data.frame(Model = "Model Averaged",
                             BMDMixture = mod.obj$BMDMixture
    )
  } else if(is.BMADRQ2(mod.obj)[2] == 2) {

    BMDMixture <- data.frame(Model = "Model Averaged",
                             BMDMixture = mod.obj$BMDMixture
    )

  } else stop("please check mod.obj or supply conv")

  return(BMDMixture)

}

#' function to extract model averaged BMD
#'
#' @param mod_obj BMDBMA model object
#' @param conv logical to indicate if mixture should be only based on converged
#'             models or all the models regardless of the convergence status
#' @return dataframe containing BMDL, BMD and BMDU
#' @export BMDMAQ_extract
#'
BMDMAQ_extract <- function(mod.obj, conv = FALSE) {

  if(is.BMADRQ2(mod.obj)[3] == 2 & conv == TRUE) {
    BMDMA <- data.frame(Model = "Model Averaged",
                        Type = c("BS", "LP"),
                        BMDL = c(mod.obj$MA_bs_conv[1], mod.obj$MA_ls_conv[1]),
                        BMD = c(mod.obj$MA_bs_conv[2], mod.obj$MA_ls_conv[2]),
                        BMDU = c(mod.obj$MA_bs_conv[3], mod.obj$MA_ls_conv[3])
    )
  } else if(is.BMADRQ2(mod.obj)[3] == 2 & conv == FALSE) {
    BMDMA <- data.frame(Model = "Model Averaged",
                        Type = c("BS", "LP"),
                        BMDL = c(mod.obj$MA_bridge_sampling[1], mod.obj$MA_laplace[1]),
                        BMD = c(mod.obj$MA_bridge_sampling[2], mod.obj$MA_laplace[2]),
                        BMDU = c(mod.obj$MA_bridge_sampling[3], mod.obj$MA_laplace[3])
    )
  } else if(is.BMADRQ2(mod.obj)[2] == 2) {

    BMDMA <- data.frame(Model = "Model Averaged",
                        Type = "LP",
                        BMDL = mod.obj$MA[1],
                        BMD = mod.obj$MA[2],
                        BMDU = mod.obj$MA[3]
    )

  } else stop("check mod.obj or supply conv")

  return(BMDMA)

}

#' function to extract the BMDL,BMD and BMDU
#'
#' @param mod_obj BMDBMA model object
#' @return datatframe containing BMDLs, BMDs and BMDUs per model
#' @export BMDLU
#'
BMDLU <- function(mod.obj, pvec = c(0.05, 0.5, 0.95), type = c('continuous', 'quantal')) {

  type <- match.arg(type)
  if(type == 'continuous') {

    model.listN <- mod.obj$parsN
    model.listLN <- mod.obj$parsLN

    md_N <- vapply(model.listN[!is.na(model.listN)], function(x){
      y <- quantile(x$BMD, probs = pvec)
      yy <- data.frame(t(y) * mod.obj$max.dose)
      return(yy)
    }, data.frame(BMDL=1,BMD=2,BMDU=3))

    mdnames_N <- vapply(model.listN[!is.na(model.listN)], function(x){
      unique(x[,"ModelName"])
    }, character(1))
    md_N <- cbind(Model = mdnames_N, as.data.frame(t(md_N)))

    md_LN <- vapply(model.listLN[!is.na(model.listLN)], function(x){
      y <- quantile(x$BMD, probs = pvec)
      yy <- data.frame(t(y)* mod.obj$max.dose)
      return(yy)
    }, data.frame(BMDL=1,BMD=2,BMDU=3))

    mdnames_LN <- vapply(model.listLN[!is.na(model.listLN)], function(x){
      unique(x[,"ModelName"])
    }, character(1))

    md_LN <- cbind(Model = mdnames_LN, as.data.frame(t(md_LN)))
    ret <- rbind(md_N, md_LN)

  } else if(type == 'quantal') {

    model.listQ <- mod.obj$parsQ

    md_Q <- vapply(model.listQ[!is.na(model.listQ)], function(x){
      y <- quantile(x$BMD, probs = pvec)
      yy <- data.frame(t(y) * mod.obj$max.dose)
      return(yy)
    }, data.frame(BMDL=1,BMD=2,BMDU=3))

    mdnames_Q <- vapply(model.listQ[!is.na(model.listQ)], function(x){
      unique(x[,"ModelName"])
    }, character(1))
    md_Q <- cbind(Model = mdnames_Q, as.data.frame(t(md_Q)))
    ret <- md_Q

  } else stop('supply data type to be quantal or continuous')

  ret2 <- do.call(cbind.data.frame,
                  apply(ret, 2, function(x)unlist(x),
                        simplify = FALSE))
  return(ret2)
}


#' function to get predicted values from the DRMS
#'
#' @param mod.obj BMDBMA model object
#' @param dose vector of dose values to use for prediction
#' @param what prediction of interest. It can either be "predicted" or "response_at_BMD"
#' @param model_averaged option to compute model averaged prediction. Defaults to FALSE
#' @param weight_type type of weight to be used. It can either be "BS" for bridge sampling or
#'              "LP" for Laplace approximation
#' @return dataframe of model prediction per model or
#'         list with dataframe of model preditions and model-averaged predictions
#' @export predict.BMADRQ
#'
predict.BMADRQ <- function(mod.obj, dose,
                           what = c("predicted", "resp_at_BMD"),
                           model_averaged = FALSE,
                           weight_type = c("BS", "LP")) {

  type <- 'quantal'
  pars <- par_med(mod.obj, type)
  what <- match.arg(what)
  q <- mod.obj$q


  if(what == "predicted") {

    preds <- vector("list", nrow(pars))
    for(i in seq_len(nrow(pars))) {
      DRMP <- DRM(pars$Model[i], type)
      ps <-  unlist(pars[i,paste0("p", 1:3)])
      rpd <- data.frame(Model = pars$Model[i],
                        Dose = dose,
                        predicted = DRMP(ps, dose, q))
      rpd <- dplyr::mutate(rpd,
                           predicted = replace(predicted,
                                               stringr::str_detect(Model, "_LN"),
                                               exp(predicted))
      )
      preds[[i]] <- rpd
    }

    preds <- do.call(rbind.data.frame,preds)
    if(model_averaged == TRUE & is.BMADRQ2(mod.obj)[3] == 2 & weight_type == "BS") {

      wts <- weightsQ_extract(mod.obj, type = "both")
      pred_wts <- merge(preds, wts, by = "Model")
      MA_pred <- dplyr::summarise(dplyr::group_by(pred_wts, Dose),
                                  model_averaged = sum(predicted*BS_Weights)
      )
      return(list(predicted = preds,
                  model_averaged = MA_pred)
      )

    } else if(model_averaged == TRUE & is.BMADRQ2(mod.obj)[3] == 2 & weight_type == "LP") {

      wts <- weightsQ_extract(mod.obj, type = "both")
      pred_wts <- merge(preds, wts, by = "Model")
      MA_pred <- dplyr::summarise(dplyr::group_by(pred_wts, Dose),
                                  model_averaged = sum(predicted*LP_Weights)
      )
      return(list(predicted = preds,
                  model_averaged = MA_pred)
      )

    } else if(model_averaged == TRUE & is.BMADRQ2(mod.obj)[2] == 2 & weight_type == "LP") {

      wts <- weightsQ_extract(mod.obj)
      pred_wts <- merge(preds, wts, by = "Model")
      MA_pred <- dplyr::summarise(dplyr::group_by(pred_wts, Dose),
                                  model_averaged = sum(predicted*LP_Weights)
      )
      return(list(predicted = preds,
                  model_averaged = MA_pred)
      )


    } else return(preds)

  } else if(what == "resp_at_BMD") {

    preds <- vector("list", nrow(pars))
    for(i in seq_len(nrow(pars))){
      DRMP <- DRM(pars$Model[i], type)
      ps <-  unlist(pars[i,paste0("p", 1:3)])
      rpd <- data.frame(Model = pars$Model[i],
                        BMD = pars$BMD[i]*mod.obj$max.dose,
                        resp_at_BMD = DRMP(ps, pars$BMD[i], q))
      rpd <- dplyr::mutate(rpd,
                           resp_at_BMD = replace(resp_at_BMD,
                                                 stringr::str_detect(Model, "_LN"),
                                                 exp(resp_at_BMD))
      )
      preds[[i]] <- rpd
    }

    preds <- do.call(rbind.data.frame,preds)

    BMDMA <- BMDMAQ_extract(mod.obj, conv = FALSE)

    if(model_averaged == TRUE & is.BMADRQ2(mod.obj)[3] == 2 & weight_type == "BS") {

      wts <- weightsQ_extract(mod.obj, type = "both")
      pred_wts <- merge(preds, wts, by = "Model")

      MA_pred <- data.frame(dplyr::filter(BMDMA, Type == "BS"),
                            model_averaged_response = sum(pred_wts$resp_at_BMD*pred_wts$BS_Weights))
      return(list(resp_at_BMD = preds,
                  model_averaged = MA_pred)
      )

    } else if(model_averaged == TRUE & is.BMADRQ2(mod.obj)[3] == 2 & weight_type == "LP") {

      wts <- weightsQ_extract(mod.obj, type = "both")
      pred_wts <- merge(preds, wts, by = "Model")
      MA_pred <- data.frame(dplyr::filter(BMDMA, Type == "LP"),
                            model_averaged_response = sum(pred_wts$resp_at_BMD*pred_wts$LP_Weights))
      return(list(resp_at_BMD = preds,
                  model_averaged = MA_pred)
      )

    } else if(model_averaged == TRUE & is.BMADRQ2(mod.obj)[2] == 2 & weight_type == "LP") {

      wts <- weightsQ_extract(mod.obj)
      pred_wts <- merge(preds, wts, by = "Model")
      MA_pred <- data.frame(dplyr::filter(BMDMA, Type == "LP"),
                            model_averaged_response = sum(pred_wts$resp_at_BMD*pred_wts$LP_Weights))
      return(list(resp_at_BMD = preds,
                  model_averaged = MA_pred)
      )

    } else return(preds)

  } else stop("Dose can only be NULL if what = 'resp_at_BMD'") ## error message seems to be in the wrong place?

}


#' function to get BMDs and weights
#'
#' @param mod_obj BMDBMA model object
#' @return datatframe containing BMDs and weights per model
#' @export BMDWeights
#'
BMDWeights <- function(mod.obj, type = c('continuous', 'quantal')) {
  type <-  match.arg(type)
  BMDs <- BMDLU(mod.obj, type = type)
  if(type == 'continuous') {

    if(is.BMADR2(mod.obj)[3] == 2 ) {
      wts_BS <- weights_extract(mod.obj, type = "BS")
      wts_LP <- weights_extract(mod.obj, type = "LP")

      return(merge(merge(BMDs, wts_BS, by = "Model", sort = FALSE),
                   wts_LP, by = "Model", sort = FALSE))
    } else if(is.BMADR2(mod.obj)[2] == 2) {

      wts_LP <- weights_extract(mod.obj, type = "LP")
      return(merge(BMDs, wts_LP, by = "Model", sort = FALSE))
    }

  } else if(type == 'quantal') {

    if(is.BMADRQ2(mod.obj)[3] == 2 ) {
      wts_BS <- weightsQ_extract(mod.obj, type = "BS")
      wts_LP <- weightsQ_extract(mod.obj, type = "LP")

      return(merge(merge(BMDs, wts_BS, by = "Model", sort = FALSE),
                   wts_LP, by = "Model", sort = FALSE))
    } else if(is.BMADRQ2(mod.obj)[2] == 2) {

      wts_LP <- weightsQ_extract(mod.obj, type = "LP")
      return(merge(BMDs, wts_LP, by = "Model", sort = FALSE))

    }


  } else stop("please check mod.obj")


}

#' print function for BMABMDR objects
#'
#' @param mod_obj BMDBMA model object
#' @return datatframe containing BMDs and weights per model
#'
#' @export print.BMADRQ
#'
print.BMADRQ <- function(mod.obj) {

  if(is.BMADRQ(mod.obj)){
    type <- 'quantal'
    print(knitr::kable(BMDWeights(mod.obj, type),
                       row.names = FALSE))
  } else stop('wrong object type supplied')

}

#' summary function for BMABMDR objects
#'
#' @param mod_obj BMDBMA model object
#' @param conv logical indicating if only converged models should be summarized.
#' @return datatframe containing BMDs, weights and parameters per model
#'
#' @export summary.BMADRQ
#'
summary.BMADRQ <- function(mod.obj, conv = FALSE) {

  BMDWeights(mod.obj, type = 'quantal')
  par_med(mod.obj, type = 'quantal')

  return(list(BMDWeights = BMDWeights(mod.obj, type = 'quantal'),
              ModelAverage = BMDMAQ_extract(mod.obj, conv = FALSE),
              ModelParameters = par_med(mod.obj, type = 'quantal')

  ))
}

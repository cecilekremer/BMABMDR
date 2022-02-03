#' Function for internal use
#'
#' @param par_obj value
#' @param pvec value
#' @param max.dose logical
#'
#' @return .
#'
outLP <- function(par_obj, pvec, max.dose){
  ret <- t(data.frame(
    # model parameters
    # transformed
    par.dt = quantile(par_obj$p4, c(0.025,0.5,0.975), na.rm=T),
    par.is2t = quantile(par_obj$is2t, c(0.025,0.5,0.975), na.rm=T),
    par.k = quantile(par_obj$k, pvec, na.rm=T),
    # original scale
    par.a = quantile(par_obj$a, c(0.025,0.5,0.975), na.rm=T),
    par.b = quantile(par_obj$b, c(0.025,0.5,0.975), na.rm=T),
    par.c = quantile(par_obj$c, c(0.025,0.5,0.975), na.rm=T),
    par.d = quantile(par_obj$d, c(0.025,0.5,0.975), na.rm=T),
    par.s2 = quantile(par_obj$s2, c(0.025,0.5,0.975), na.rm=T),
    # natural parameters
    BMD = quantile(par_obj$BMD*max.dose, pvec, na.rm=T),
    min.resp = quantile(par_obj$min_response, c(0.025,0.5,0.975), na.rm=T),
    max.resp = quantile(par_obj$max_response, c(0.025,0.5,0.975), na.rm=T),
    fold.change = quantile(par_obj$p3, c(0.025, 0.5, 0.975), na.rm=T)
  ))

  return(ret)
}

#' Function for internal use
#'
#' @param par_obj value
#' @param pvec value
#' @param max.dose logical
#' @param rho logical
#'
#' @return .
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


#' Function for internal use
#'
#' @param mod_obj value
#' @param pars value
#' @param model_name model name
#'
#' @return .
#'
par_extract <- function(mod_obj, pars = c(letters[1:4],"k",
                                          paste0("par",1:5,""),
                                          "invsigma2","mu_inf","mu_0"), model_name) {

  if(is.stanfit(mod_obj)){

    pars_samples <- rstan::extract(mod_obj, pars)
    pars_samplesd <- data.frame(ModelName = model_name,
                                a = pars_samples$a,
                                b = pars_samples$b,
                                c = pars_samples$c,
                                d = pars_samples$d,
                                k = pars_samples$k,
                                BMD = pars_samples$`par2`,
                                p1 = pars_samples$`par1`, # mu(0)
                                p2 = pars_samples$`par2`, # BMD
                                p3 = pars_samples$`par3`, # fold change
                                p4 = pars_samples$`par4`, # log(d)
                                is2t = pars_samples$`par5`,
                                s2 = 1/pars_samples$invsigma2,
                                min_response = pars_samples$`mu_0`,
                                max_response = pars_samples$`mu_inf`
    )

  } else if(is.stanfitOptim(mod_obj)) {

    pars_samples <- mod_obj$theta_tilde[,pars]
    pars_samplesd <- data.frame(ModelName = model_name,
                                a = pars_samples[,"a"],
                                b = pars_samples[,"b"],
                                c = pars_samples[,"c"],
                                d = pars_samples[,"d"],
                                k = pars_samples[,"k"],
                                BMD = pars_samples[,"par2"],
                                p1 = pars_samples[,"par1"],
                                p2 = pars_samples[,"par2"],
                                p3 = pars_samples[,"par3"],
                                p4 = pars_samples[,"par4"],
                                is2t = pars_samples[,"par5"],
                                s2 = 1/pars_samples[,"invsigma2"],
                                min_response = pars_samples[,"mu_0"],
                                max_response = pars_samples[,"mu_inf"]
    )

    return(pars_samplesd)

  } else stop("please check if input is of class stanfit or stanfitOptim")



}

#' Function for internal use
#'
#' @param mod_obj value
#' @param pars value
#' @param model_name model name
#' @param rho logical
#'
#' @return .
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


#' Function for internal use
#'
#' @param mod.obj value
#' @param type value
#'
#' @return .
#'
weights_extract <- function(mod.obj, type = c("BS", "LP", "both")) {
  type <- match.arg(type)
  if(is.BMADR2(mod.obj)[3] == 2 & type == 'BS'){

    return(data.frame(Model = names(mod.obj$weights_bridge_sampling),
                      BS_Weights = mod.obj$weights_bridge_sampling))

  } else if(is.BMADR2(mod.obj)[3] == 2 & type == 'LP') {

    return(data.frame(Model = names(mod.obj$weights_bridge_sampling),
                      LP_Weights = mod.obj$weights_laplace))

  } else if(is.BMADR2(mod.obj)[3] == 2 & type == "both"){

    return(data.frame(Model = names(mod.obj$weights_bridge_sampling),
                      BS_Weights = mod.obj$weights_bridge_sampling,
                      LP_Weights = mod.obj$weights_laplace))

  } else if(is.BMADR2(mod.obj)[2] == 2){

    return(data.frame(Model = names(mod.obj$weights),
                      LP_Weights = mod.obj$weights))

  }else stop('please check mod.obj or provide type')

}

#' Function for internal use
#'
#' @param mod.obj value
#' @param type value
#'
#' @return .
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


#' Function for internal use
#'
#' @param mod.obj value
#' @param conv value
#'
#' @return .
#'
BMDmixture_extract <- function(mod.obj, conv = FALSE){

  if(is.BMADR2(mod.obj)[3] == 2 & conv == TRUE) {
    BMDMixture <- data.frame(Model = "Model Averaged",
                             BMDMixture = mod.obj$BMDMixture.conv
    )
  } else if(is.BMADR2(mod.obj)[3] == 2 & conv == FALSE) {
    BMDMixture <- data.frame(Model = "Model Averaged",
                             BMDMixture = mod.obj$BMDMixture
    )
  } else if(is.BMADR2(mod.obj)[2] == 2) {

    BMDMixture <- data.frame(Model = "Model Averaged",
                             BMDMixture = mod.obj$BMDMixture
    )

  } else stop("please check mod.obj or supply conv")

  return(BMDMixture)

}

#' Function for internal use
#'
#' @param mod.obj value
#' @param conv value
#'
#' @return .
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

#' Function for internal use
#'
#' @param mod.obj value
#' @param conv value
#'
#' @return .
#'
BMDMA_extract <- function(mod.obj, conv = FALSE) {

  if(is.BMADR2(mod.obj)[3] == 2 & conv == TRUE) {
    BMDMA <- data.frame(Model = "Model Averaged",
                        Type = c("BS", "LP"),
                        BMDL = c(mod.obj$MA_bs_conv[1], mod.obj$MA_ls_conv[1]),
                        BMD = c(mod.obj$MA_bs_conv[2], mod.obj$MA_ls_conv[2]),
                        BMDU = c(mod.obj$MA_bs_conv[3], mod.obj$MA_ls_conv[3])
    )
  } else if(is.BMADR2(mod.obj)[3] == 2 & conv == FALSE) {
    BMDMA <- data.frame(Model = "Model Averaged",
                        Type = c("BS", "LP"),
                        BMDL = c(mod.obj$MA_bridge_sampling[1], mod.obj$MA_laplace[1]),
                        BMD = c(mod.obj$MA_bridge_sampling[2], mod.obj$MA_laplace[2]),
                        BMDU = c(mod.obj$MA_bridge_sampling[3], mod.obj$MA_laplace[3])
    )
  } else if(is.BMADR2(mod.obj)[2] == 2) {

    BMDMA <- data.frame(Model = "Model Averaged",
                        Type = "LP",
                        BMDL = mod.obj$MA[1],
                        BMD = mod.obj$MA[2],
                        BMDU = mod.obj$MA[3]
    )

  } else stop("check mod.obj or supply conv")

  return(BMDMA)

}

#' Function for internal use
#'
#' @param mod.obj value
#' @param conv value
#'
#' @return .
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


#' Function for internal use
#'
#' @param mod.obj value
#' @param type value
#'
#' @return .
#'
par_med <- function(mod.obj, type = c('continuous','quantal')) {

  type <- match.arg(type)
  if(type == 'continuous') {

    model.listN <- mod.obj$parsN
    model.listLN <- mod.obj$parsLN

    md_N <- vapply(model.listN[!is.na(model.listN)], function(x){
      y <- apply(x[,-c(1)], 2, median, na.rm = TRUE)
      yy <- data.frame(t(y))
      names(yy) <- names(x)[-1]
      return(yy)
    }, data.frame(a=1,b=2,c=3,d=4,k=5,BMD=6,p1=7,p2=8,p3=9,p4=10,
                  is2t=11, s2=12, min_response=13, max_response=14))

    mdnames_N <- vapply(model.listN[!is.na(model.listN)], function(x){
      unique(x[,"ModelName"])
    }, character(1))

    md_N <- cbind(Model = mdnames_N, as.data.frame(t(md_N)))

    md_LN <- vapply(model.listLN[!is.na(model.listLN)], function(x){
      y <- apply(x[,-c(1)], 2, median, na.rm = TRUE)
      yy <- data.frame(t(y))
      names(yy) <- names(x)[-1]
      return(yy)
    }, data.frame(a=1,b=2,c=3,d=4,k=5,BMD=6,p1=7,p2=8,p3=9,p4=10,
                  is2t=11, s2=12, min_response=13, max_response=14))

    mdnames_LN <- vapply(model.listLN[!is.na(model.listLN)], function(x){
      unique(x[,"ModelName"])
    }, character(1))


    md_LN <- cbind(Model = mdnames_LN, as.data.frame(t(md_LN)))
    ret <- rbind(md_N, md_LN)

  }else if(type == 'quantal') {

    model.listQ <- mod.obj$parsQ

    md_Q <- vapply(model.listQ[!is.na(model.listQ)], function(x){
      y <- apply(x[,-c(1)], 2, median, na.rm = TRUE)
      yy <- data.frame(t(y))
      names(yy) <- names(x)[-1]
      return(yy)
    }, data.frame(a=1,b=2,d=3,BMD=4,rho=0,p1=5,p2=6,p3=7))

    mdnames_Q <- vapply(model.listQ[!is.na(model.listQ)], function(x){
      unique(x[,"ModelName"])
    }, character(1))

    md_Q <- cbind(Model = mdnames_Q, as.data.frame(t(md_Q)))
    ret <- md_Q
  }else stop('supply data type to be quantal or continuous')

  ret2 <- do.call(cbind.data.frame,
                  apply(ret, 2, function(x)unlist(x),
                        simplify = FALSE))
  return(ret2)
}

#' Function for internal use
#'
#' @param mod.obj value
#' @param pvec value
#' @param type value
#'
#' @return .
#'
BMDLU <- function(mod.obj, pvec = c(0.05, 0.5, 0.95), type = c('continuous','quantal')) {

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
  }else if(type == 'quantal') {

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

  }else stop('supply data type to be quantal or continuous')

  ret2 <- do.call(cbind.data.frame,
                  apply(ret, 2, function(x)unlist(x),
                        simplify = FALSE))
  return(ret2)
}

#' Function to get predicted values from the DRMS
#'
#' @param mod.obj is the result of pars_med
#' @param dose dose
#' @param type dose response model type
#' @param what prediction of interest
#' @param model_averaged option to compute model averaged prediction
#' @param weight_type type of weight to be used
#'
#' @return .
#'
#' @export predict.BMADR
#'
predict.BMADR <- function(mod.obj, dose,
                          type = c("increasing", "decreasing"),
                          what = c("predicted", "resp_at_BMD"),
                          model_averaged = FALSE,
                          weight_type = c("BS", "LP")) {

  pars <- par_med(mod.obj) # model parameters (median), for each model
  type <- match.arg(type)
  what <- match.arg(what)
  q <- mod.obj$q


  if(what == "predicted") {

    preds <- vector("list", nrow(pars))
    for(i in seq_len(nrow(pars))) {
      DRMP <- DRM(pars$Model[i], type)
      ps <- unlist(pars[i,paste0("p", 1:4)])
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
    if(model_averaged == TRUE & is.BMADR2(mod.obj)[3] == 2 & weight_type == "BS") {

      wts <- weights_extract(mod.obj, type = "both")
      pred_wts <- merge(preds, wts, by = "Model")
      MA_pred <- dplyr::summarise(dplyr::group_by(pred_wts, Dose),
                                  model_averaged = sum(predicted*BS_Weights)
      )
      return(list(predicted = preds,
                  model_averaged = MA_pred)
      )

    } else if(model_averaged == TRUE & is.BMADR2(mod.obj)[3] == 2 & weight_type == "LP") {

      wts <- weights_extract(mod.obj, type = "both")
      pred_wts <- merge(preds, wts, by = "Model")
      MA_pred <- dplyr::summarise(dplyr::group_by(pred_wts, Dose),
                                  model_averaged = sum(predicted*LP_Weights)
      )
      return(list(predicted = preds,
                  model_averaged = MA_pred)
      )

    } else if(model_averaged == TRUE & is.BMADR2(mod.obj)[2] == 2 & weight_type == "LP") {

      wts <- weights_extract(mod.obj)
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
      ps <- unlist(pars[i,paste0("p", 1:4)])
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

    BMDMA <- BMDMA_extract(mod.obj, conv = FALSE)

    if(model_averaged == TRUE & is.BMADR2(mod.obj)[3] == 2 & weight_type == "BS") {

      wts <- weights_extract(mod.obj, type = "both")
      pred_wts <- merge(preds, wts, by = "Model")

      MA_pred <- data.frame(dplyr::filter(BMDMA, Type == "BS"),
                            model_averaged_response = sum(pred_wts$resp_at_BMD*pred_wts$BS_Weights))
      return(list(resp_at_BMD = preds,
                  model_averaged = MA_pred)
      )

    } else if(model_averaged == TRUE & is.BMADR2(mod.obj)[3] == 2 & weight_type == "LP") {

      wts <- weights_extract(mod.obj, type = "both")
      pred_wts <- merge(preds, wts, by = "Model")
      MA_pred <- data.frame(dplyr::filter(BMDMA, Type == "LP"),
                            model_averaged_response = sum(pred_wts$resp_at_BMD*pred_wts$LP_Weights))
      return(list(resp_at_BMD = preds,
                  model_averaged = MA_pred)
      )

    } else if(model_averaged == TRUE & is.BMADR2(mod.obj)[2] == 2 & weight_type == "LP") {

      wts <- weights_extract(mod.obj)
      pred_wts <- merge(preds, wts, by = "Model")
      MA_pred <- data.frame(dplyr::filter(BMDMA, Type == "LP"),
                            model_averaged_response = sum(pred_wts$resp_at_BMD*pred_wts$LP_Weights))
      return(list(resp_at_BMD = preds,
                  model_averaged = MA_pred)
      )

    } else return(preds)

  } else stop("Dose can only be NULL if what = 'resp_at_BMD'") ## error message seems to be in the wrong place?

}

#' Function to get predicted values from the DRMS
#'
#' @param mod.obj is the result of pars_med
#' @param type dose response model type
#' @param what prediction of interest
#' @param model_averaged option to compute model averaged prediction
#' @param weight_type type of weight to be used
#'
#' @return .
#'
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


#' Function for internal use
#'
#' @param mod.obj is the result of pars_med
#' @param type value
#'
#' @return .
#'
BMDWeights <- function(mod.obj, type = c('continuous','quantal')) {

  type <-  match.arg(type)
  BMDs <- BMDLU(mod.obj, type = type)

  if(type == 'continuous') {

    if(is.BMADR2(mod.obj)[3] == 2) {
      wts_BS <- weights_extract(mod.obj, type = "BS")
      wts_LP <- weights_extract(mod.obj, type = "LP")

      return(merge(merge(BMDs, wts_BS, by = "Model", sort = FALSE),
                   wts_LP, by = "Model", sort = FALSE))
    } else if(is.BMADR2(mod.obj)[2] == 2) {

      wts_LP <- weights_extract(mod.obj, type = "LP")
      return(merge(BMDs, wts_LP, by = "Model", sort = FALSE))
    }

  }else if(type == 'quantal') {

    if(is.BMADRQ2(mod.obj)[3] == 2 ) {
      wts_BS <- weightsQ_extract(mod.obj, type = "BS")
      wts_LP <- weightsQ_extract(mod.obj, type = "LP")

      return(merge(merge(BMDs, wts_BS, by = "Model", sort = FALSE),
                   wts_LP, by = "Model", sort = FALSE))
    } else if(is.BMADRQ2(mod.obj)[2] == 2) {

      wts_LP <- weightsQ_extract(mod.obj, type = "LP")
      return(merge(BMDs, wts_LP, by = "Model", sort = FALSE))

    }


  }else stop("please check mod.obj")


}


#' Function for internal use
#'
#' @param mod.obj is the result of pars_med
#'
#' @return .
#'
#' @export
#'
print.BMADR <- function(mod.obj) {
  message("Here are the estimated BMDs per model along with their weights")
  print(knitr::kable(BMDWeights(mod.obj),
                     row.names = FALSE))
}

#' Function for internal use
#'
#' @param mod.obj is the result of pars_med
#'
#' @return .
#'
#' @export
#'
print.BMADRQ <- function(mod.obj) {

  if(is.BMADRQ(mod.obj)){
    type <- 'quantal'
    print(knitr::kable(BMDWeights(mod.obj, type),
                       row.names = FALSE))
  } else stop('wrong object type supplied')

}

#' Function for internal use
#'
#' @param mod.obj is the result of pars_med
#' @param conv
#'
#' @return .
#'
#' @export
#'
summary.BMADR <- function(mod.obj, conv = FALSE) {
  BMDWeights(mod.obj)
  par_med(mod.obj)
  BMDMA_extract(mod.obj, conv = FALSE)
  return(list(BMDWeights = BMDWeights(mod.obj),
              ModelAverage = BMDMA_extract(mod.obj, conv = FALSE),
              ModelParameter = par_med(mod.obj)

  ))
}

#' Function for internal use
#'
#' @param mod.obj is the result of pars_med
#' @param conv
#'
#' @return .
#'
#' @export
#'
summary.BMADRQ <- function(mod.obj, conv = FALSE) {

  BMDWeights(mod.obj, type = 'quantal')
  par_med(mod.obj, type = 'quantal')

  return(list(BMDWeights = BMDWeights(mod.obj, type = 'quantal'),
              ModelAverage = BMDMAQ_extract(mod.obj, conv = FALSE),
              ModelParameters = par_med(mod.obj, type = 'quantal')

  ))
}


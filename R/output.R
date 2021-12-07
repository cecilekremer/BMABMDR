#' Function for internal use
#'
#' @param par_obj value
#' @param pvec value
#' @param max.dose logical
#'
#' @return .
#'
outLP <- function(par_obj, pvec, max.dose) {
  ret <- t(data.frame(
    # transformed parameters
    par.at = quantile(par_obj$p1, c(0.025,0.5,0.975), na.rm=T),
    par.ct = quantile(par_obj$p3, c(0.025,0.5,0.975), na.rm=T),
    par.dt = quantile(par_obj$p4, c(0.025,0.5,0.975), na.rm=T),
    par.is2t = quantile(par_obj$is2t, c(0.025,0.5,0.975), na.rm=T),
    par.k = quantile(par_obj$p2, pvec, na.rm=T),
    # parameters on original scale
    par.a = quantile(par_obj$a, c(0.025,0.5,0.975), na.rm=T),
    par.b = quantile(par_obj$b, c(0.025,0.5,0.975), na.rm=T),
    par.c = quantile(par_obj$c, c(0.025,0.5,0.975), na.rm=T),
    par.d = quantile(par_obj$d, c(0.025,0.5,0.975), na.rm=T),
    par.s2 = quantile(par_obj$s2, c(0.025,0.5,0.975), na.rm=T),
    # natural parameters
    BMD = quantile(par_obj$BMD, pvec, na.rm=T)*max.dose,
    min.resp = quantile(par_obj$min_response, c(0.025,0.5,0.975), na.rm=T),
    max.resp = quantile(par_obj$max_response, c(0.025,0.5,0.975), na.rm=T)
  ))

  return(ret)
}


#' Function for internal use
#'
#' @param mod_obj value
#' @param pars value
#'
#' @return .
#'
par_extract <- function(mod_obj, pars = c(letters[1:4], "BMD",
                                          paste0("par[",1:5,"]"),
                                          "invsigma2", "min_response",
                                          "max_response"), model_name) {

  if(is.stanfit(mod_obj)){

    pars_samples <- rstan::extract(mod_obj, pars)
    pars_samplesd <- data.frame(ModelName = model_name,
                                a = pars_samples$a,
                                b = pars_samples$b,
                                c = pars_samples$c,
                                d = pars_samples$d,
                                BMD = pars_samples$BMD,
                                p1 = pars_samples$`par[1]`,
                                p2 = pars_samples$`par[2]`,
                                p3 = pars_samples$`par[3]`,
                                p4 = pars_samples$`par[4]`,
                                is2t = pars_samples$`par[5]`,
                                s2 = 1/pars_samples$invsigma2,
                                min_response = pars_samples$min_response,
                                max_response = pars_samples$max_response
    )

  } else if(is.stanfitOptim(mod_obj)) {

    pars_samples <- mod_obj$theta_tilde[,pars]
    pars_samplesd <- data.frame(ModelName = model_name,
                                a = pars_samples[,"a"],
                                b = pars_samples[,"b"],
                                c = pars_samples[,"c"],
                                d = pars_samples[,"d"],
                                BMD = pars_samples[,"BMD"],
                                p1 = pars_samples[,"par[1]"],
                                p2 = pars_samples[,"par[2]"],
                                p3 = pars_samples[,"par[3]"],
                                p4 = pars_samples[,"par[4]"],
                                is2t = pars_samples[,"par[5]"],
                                s2 = 1/pars_samples[,"invsigma2"],
                                min_response = pars_samples[,"min_response"],
                                max_response = pars_samples[,"max_response"]
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
#'
#' @return .
#'
par_med <- function(mod.obj) {

  model.listN <- mod.obj$parsN
  model.listLN <- mod.obj$parsLN

  md_N <- vapply(model.listN[!is.na(model.listN)], function(x){
    y <- apply(x[,-c(1)], 2, median, na.rm = TRUE)
    yy <- data.frame(t(y))
    names(yy) <- names(x)[-1]
    return(yy)
  }, data.frame(a=1,b=2,c=3,d=4,BMD=5,p1=6,p2=7,p3=8,p4=9,
                is2t=10, s2=11, min_response=12, max_response=13))

  mdnames_N <- vapply(model.listN[!is.na(model.listN)], function(x){
    unique(x[,"ModelName"])
  }, character(1))

  md_N <- cbind(Model = mdnames_N, as.data.frame(t(md_N)))

  md_LN <- vapply(model.listLN[!is.na(model.listLN)], function(x){
    y <- apply(x[,-c(1)], 2, median, na.rm = TRUE)
    yy <- data.frame(t(y))
    names(yy) <- names(x)[-1]
    return(yy)
  }, data.frame(a=1,b=2,c=3,d=4,BMD=5,p1=6,p2=7,p3=8,p4=9,
                is2t=10, s2=11, min_response=12, max_response=13))

  mdnames_LN <- vapply(model.listLN[!is.na(model.listLN)], function(x){
    unique(x[,"ModelName"])
  }, character(1))


  md_LN <- cbind(Model = mdnames_LN, as.data.frame(t(md_LN)))
  ret <- rbind(md_N, md_LN)
  ret2 <- do.call(cbind.data.frame,
                  apply(ret, 2, function(x)unlist(x),
                        simplify = FALSE))
  return(ret2)
}

#' Function for internal use
#'
#' @param mod.obj value
#' @param pvec value
#'
#' @return .
#'
BMDLU <- function(mod.obj, pvec = c(0.05, 0.5, 0.95)) {

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
  ret2 <- do.call(cbind.data.frame,
                  apply(ret, 2, function(x)unlist(x),
                        simplify = FALSE))
  return(ret2)
}

#' Function for internal use
#'
#' @param mod.obj is the result of pars_med
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

#' Function for internal use
#'
#' @param mod.obj is the result of pars_med
#'
#' @return .
#'
BMDWeights <- function(mod.obj) {

  BMDs <- BMDLU(mod.obj)
  if(is.BMADR2(mod.obj)[3] == 2) {
    wts_BS <- weights_extract(mod.obj, type = "BS")
    wts_LP <- weights_extract(mod.obj, type = "LP")

    return(merge(merge(BMDs, wts_BS, by = "Model", sort = FALSE),
                 wts_LP, by = "Model", sort = FALSE))
  } else if(is.BMADR2(mod.obj)[2] == 2) {

    wts_LP <- weights_extract(mod.obj, type = "LP")
    return(merge(BMDs, wts_LP, by = "Model", sort = FALSE))
  } else stop("please check mod.obj")


}

#' Function for internal use
#'
#' @param mod.obj is the result of pars_med
#'
#' @return .
#'
print.BMADR <- function(mod.obj) {
  message("Here are the estimated BMDs per model along with their weights")
  print(knitr::kable(BMDWeights(mod.obj),
                     row.names = FALSE))
}

#' Function for internal use
#'
#' @param mod.obj is the result of pars_med
#' @param conv
#'
#' @return .
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

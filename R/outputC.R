
#' output function for estimated model parameter from each model
#'
#' @param par_obj list containing parameter draws from either Laplace or MCMC fit
#' @param pvec probability vector to compute credible interval for the BMD.
#' @param max.dose maximum dose tested
#'
#' @return dataframe containing model estimates
#'
#'
#' @export outLPC
#'
outLPC <- function(par_obj, pvec, max.dose){
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
    fold.change = quantile(par_obj$p3, c(0.025, 0.5, 0.975), na.rm=T),
    rho = quantile(par_obj$rho, c(0.025, 0.5, 0.975), na.rm=T)
  ))

  return(ret)
}


#### function to extract the sampled values for the model parameters and the BMD
#' function to extract the sampled values for the model parameters and the BMD
#'
#' @param mod_obj model object either of class stanfit or stanfitOptim
#' @param pars parameters to be extracted.
#' @param model_name name of the model to be extracted
#'
#' @return dataframe containing model estimates
#'
#' @export par_extract
#'
par_extractC <- function(mod_obj, pars = c(letters[1:4],"k",
                                           paste0("par",1:6,""),
                                           "invsigma2","mu_inf","mu_0","rho_cluster"), model_name) {

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
                                max_response = pars_samples$`mu_inf`,
                                rho = pars_samples$`rho_cluster`
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
                                max_response = pars_samples[,"mu_inf"],
                                rho = pars_samples[,"rho_cluster"]
    )

    return(pars_samplesd)

  } else stop("please check if input is of class stanfit or stanfitOptim")



}

#' function to get the median of the needed parameters for predicted values. This is done over all
#' fitted models
#'
#' @param mod_obj BMDBMA model object
#' @return dataframe containing parameter estimates per model
#' @export par_med
#'
par_med <- function(mod.obj, type = c('continuous', 'quantal'), clustered = F) {

  type <- match.arg(type)
  if(type == 'continuous') {

    model.listN <- mod.obj$parsN
    model.listLN <- mod.obj$parsLN

    if(clustered == T){
      md_N <- vapply(model.listN[!is.na(model.listN)], function(x){
        y <- apply(x[,-c(1)], 2, median, na.rm = TRUE)
        yy <- data.frame(t(y))
        names(yy) <- names(x)[-1]
        return(yy)
      }, data.frame(a=1,b=2,c=3,d=4,k=5,BMD=6,p1=7,p2=8,p3=9,p4=10,
                    is2t=11, s2=12, min_response=13, max_response=14, rho=15))
    }else if(clustered == F){
      md_N <- vapply(model.listN[!is.na(model.listN)], function(x){
        y <- apply(x[,-c(1)], 2, median, na.rm = TRUE)
        yy <- data.frame(t(y))
        names(yy) <- names(x)[-1]
        return(yy)
      }, data.frame(a=1,b=2,c=3,d=4,k=5,BMD=6,p1=7,p2=8,p3=9,p4=10,
                    is2t=11, s2=12, min_response=13, max_response=14))
    }


    mdnames_N <- vapply(model.listN[!is.na(model.listN)], function(x){
      unique(x[,"ModelName"])
    }, character(1))

    md_N <- cbind(Model = mdnames_N, as.data.frame(t(md_N)))

    if(clustered == T){
      md_LN <- vapply(model.listLN[!is.na(model.listLN)], function(x){
        y <- apply(x[,-c(1)], 2, median, na.rm = TRUE)
        yy <- data.frame(t(y))
        names(yy) <- names(x)[-1]
        return(yy)
      }, data.frame(a=1,b=2,c=3,d=4,k=5,BMD=6,p1=7,p2=8,p3=9,p4=10,
                    is2t=11, s2=12, min_response=13, max_response=14, rho=15))
    }else if(clustered == F){
      md_LN <- vapply(model.listLN[!is.na(model.listLN)], function(x){
        y <- apply(x[,-c(1)], 2, median, na.rm = TRUE)
        yy <- data.frame(t(y))
        names(yy) <- names(x)[-1]
        return(yy)
      }, data.frame(a=1,b=2,c=3,d=4,k=5,BMD=6,p1=7,p2=8,p3=9,p4=10,
                    is2t=11, s2=12, min_response=13, max_response=14))
    }


    mdnames_LN <- vapply(model.listLN[!is.na(model.listLN)], function(x){
      unique(x[,"ModelName"])
    }, character(1))


    md_LN <- cbind(Model = mdnames_LN, as.data.frame(t(md_LN)))
    ret <- rbind(md_N, md_LN)

  } else if(type == 'quantal') {

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
  } else stop('supply data type to be quantal or continuous')

  ret2 <- do.call(cbind.data.frame,
                  apply(ret, 2, function(x)unlist(x),
                        simplify = FALSE))
  return(ret2)
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
#' @param type dose response model type
#' @param what prediction of interest. It can either be "predicted" or "response_at_BMD"
#' @param model_averaged option to compute model averaged prediction. Defaults to FALSE
#' @param weight_type type of weight to be used. It can either be "BS" for bridge sampling or
#'              "LP" for Laplace approximation
#' @return dataframe of model prediction per model or
#'         list with dataframe of model preditions and model-averaged predictions
#' @export predict.BMADR
#'
predict.BMADR <- function(mod.obj, dose,
                          type = c("increasing", "decreasing"),
                          what = c("predicted", "resp_at_BMD"),
                          model_averaged = FALSE,
                          clustered = FALSE,
                          weight_type = c("BS", "LP")){

  if(clustered == F){
    pars <- par_med(mod.obj) # model parameters (median), for each model
  }else if(clustered == T){
    pars <- par_med(mod.obj, clustered = T)
  }
  type <- match.arg(type)
  what <- match.arg(what)
  q <- mod.obj$q
  shift <- mod.obj$shift


  if(what == "predicted") {

    preds <- vector("list", nrow(pars))
    for(i in seq_len(nrow(pars))) {
      DRMP <- DRM(pars$Model[i], type)
      ps <- unlist(pars[i,paste0("p", 1:4)])
      rpd <- data.frame(Model = pars$Model[i],
                        Dose = dose,
                        predicted = DRMP(ps, dose, q, shift))
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
                        resp_at_BMD = DRMP(ps, pars$BMD[i], q, shift))
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

#' function to get BMDs and weights
#'
#' @param mod_obj BMDBMA model object
#'
#' @return datatframe containing BMDs and weights per model
#'
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
#' @export print.BMADR
#'
print.BMADR <- function(mod.obj) {
  message("Here are the estimated BMDs per model along with their weights")
  print(knitr::kable(BMDWeights(mod.obj),
                     row.names = FALSE))
}

#' summary function for BMABMDR objects
#'
#' @param mod_obj BMDBMA model object
#' @param conv logical indicating if only converged models should be summarized.
#' @return datatframe containing BMDs, weights and parameters per model
#'
#' @export summary.BMADR
#'
summary.BMADR <- function(mod.obj, conv = FALSE, clustered = FALSE) {
  BMDWeights(mod.obj)
  if(clustered == F){
    par_med(mod.obj)
  }else if(clustered == T){
    par_med(mod.obj, clustered = T)
  }
  BMDMA_extract(mod.obj, conv = FALSE)
  return(list(BMDWeights = BMDWeights(mod.obj),
              ModelAverage = BMDMA_extract(mod.obj, conv = FALSE),
              ModelParameter = ifelse(clustered == F, par_med(mod.obj), par_med(mod.obj, clustered= T))

  ))
}

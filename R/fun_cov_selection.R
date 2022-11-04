
#' Function to perform model selection when including a covariate effect (used internally)
#'
#' @param optMod model
#' @param lld loglikelihood
#' @param min.ll minimum loglikelihood among submodels
#' @param nlevels number of covariate levels
#' @param dataMod data for optMod
#' @param covar which parameters depend on covariate
#'
#' @description This function computes Laplace weights for model selection
#'
#' @importFrom truncnorm dtruncnorm
#'
#' @return .
#'
#' @export fun.w
#'
fun.w <- function(optMod, lld, min.ll, nlevels,
                  dataMod,
                  covar = c('a_sigma2', 'BMD_d', 'all', 'none')) {

  covar = match.arg(covar, c('a_sigma2', 'BMD_d', 'all', 'none'))
  pars <- optMod$par
  nmpar <- names(pars)

  if(covar != 'none'){

    a <- pars[grep("par1\\[", nmpar)]
    #b <- pars[grep("b\\[", nmpar)]
    bmd <- pars[grep("par2\\[", nmpar)]
    # c <- pars[grep("c\\[", nmpar)]
    c <- pars[nmpar == "par3"]
    d <- pars[grep("par4\\[", nmpar)]
    s2 <- pars[grep("par5\\[", nmpar)]

  } else {
    a <- pars[nmpar == "par1"]
    #b <- pars[nmpar == "b"]
    bmd <- pars[nmpar == "par2"]
    c <- pars[nmpar == "par3"]
    d <- pars[nmpar == "par4"]
    s2 <- pars[nmpar == "par5"]
  }

  DIH <- det(-solve(optMod$hessian))
  DIH <- ifelse(DIH<0, 0, DIH)

  npar <- sum(c(length(a), length(bmd), length(c), length(d), length(s2)))

  if(covar == 'a_sigma2') {

    w1 <- sum((2*pi)^(npar/2)*sqrt(DIH)*exp(lld-min.ll)*

                # sigma2
                dnorm(s2, dataMod$data$priormu[5,], dataMod$data$priorSigma[5,5])*
                # d
                truncnorm::dtruncnorm(d, b = dataMod$data$truncd, mean = dataMod$data$priormu[4,1],
                                      sd = dataMod$data$priorSigma[4,4])*
                # a
                mc2d::dpert(a, min = dataMod$data$priorlb[1,], max = dataMod$data$priorub[1,],
                            mode = dataMod$data$priormu[1,], shape = dataMod$data$shape.a)*
                # BMD
                mc2d::dpert(bmd, min = dataMod$data$priorlb[2,1], max = dataMod$data$priorub[2,1],
                            mode = dataMod$data$priormu[2,1], shape = dataMod$data$shape.BMD)*
                # c
                mc2d::dpert(c, min = dataMod$data$priorlb[3,], max = dataMod$data$priorub[3,],
                            mode = dataMod$data$priormu[3,], shape = dataMod$data$shape.c))

  } else if(covar == 'BMD_d') {

    w1 <- sum((2*pi)^(npar/2)*sqrt(DIH)*exp(lld-min.ll)*

                # sigma2
                dnorm(s2, dataMod$data$priormu[5,1], dataMod$data$priorSigma[5,5])*
                # d
                truncnorm::dtruncnorm(d, b = dataMod$data$truncd, mean = dataMod$data$priormu[4,],
                                      sd = dataMod$data$priorSigma[4,4])*
                # a
                mc2d::dpert(a, min = dataMod$data$priorlb[1,1], max = dataMod$data$priorub[1,1],
                            mode = dataMod$data$priormu[1,1], shape = dataMod$data$shape.a)*
                # BMD
                mc2d::dpert(bmd, min = dataMod$data$priorlb[2,], max = dataMod$data$priorub[2,],
                            mode = dataMod$data$priormu[2,], shape = dataMod$data$shape.BMD)*
                # c
                mc2d::dpert(c, min = dataMod$data$priorlb[3,1], max = dataMod$data$priorub[3,1],
                            mode = dataMod$data$priormu[3,1], shape = dataMod$data$shape.c))

  } else if(covar == 'all') {

    w1 <- sum((2*pi)^(npar/2)*sqrt(DIH)*exp(lld-min.ll)*

                # sigma2
                dnorm(s2, dataMod$data$priormu[5,], dataMod$data$priorSigma[5,5])*
                # d
                truncnorm::dtruncnorm(d, b = dataMod$data$truncd, mean = dataMod$data$priormu[4,],
                                      sd = dataMod$data$priorSigma[4,4])*
                # a
                mc2d::dpert(a, min = dataMod$data$priorlb[1,], max = dataMod$data$priorub[1,],
                            mode = dataMod$data$priormu[1,], shape = dataMod$data$shape.a)*
                # BMD
                mc2d::dpert(bmd, min = dataMod$data$priorlb[2,], max = dataMod$data$priorub[2,],
                            mode = dataMod$data$priormu[2,], shape = dataMod$data$shape.BMD)*
                # c
                mc2d::dpert(c, min = dataMod$data$priorlb[3,], max = dataMod$data$priorub[3,],
                            mode = dataMod$data$priormu[3,], shape = dataMod$data$shape.c))

  } else {

    npar <- 5

    w1 <- (2*pi)^(npar/2)*sqrt(DIH)*exp(lld-min.ll)*

      # sigma2
      dnorm(s2, dataMod$data$priormu[5], dataMod$data$priorSigma[5,5])*
      # d
      truncnorm::dtruncnorm(d, b = dataMod$data$truncd, mean = dataMod$data$priormu[4],
                            sd = dataMod$data$priorSigma[4,4])*
      # a
      mc2d::dpert(a, min = dataMod$data$priorlb[1], max = dataMod$data$priorub[1],
                  mode = dataMod$data$priormu[1], shape = dataMod$data$shape.a)*
      # BMD
      mc2d::dpert(bmd, min = dataMod$data$priorlb[2], max = dataMod$data$priorub[2],
                  mode = dataMod$data$priormu[2], shape = dataMod$data$shape.BMD)*
      # c
      mc2d::dpert(c, min = dataMod$data$priorlb[3], max = dataMod$data$priorub[3],
                  mode = dataMod$data$priormu[3], shape = dataMod$data$shape.c)
  }

  return(w1)
}

#' @rdname fun.w
#' @export
fun.w.QE4 <- function(optMod, lld, min.ll, nlevels,
                      dataMod,
                      covar = c('a_sigma2', 'BMD_d', 'all', 'none')) {

  covar = match.arg(covar, c('a_sigma2', 'BMD_d', 'all', 'none'))
  pars <- optMod$par
  nmpar <- names(pars)

  if(covar != 'none'){

    a <- pars[grep("par1\\[", nmpar)]
    #b <- pars[grep("b\\[", nmpar)]
    bmd <- pars[grep("par2\\[", nmpar)]
    # c <- pars[grep("c\\[", nmpar)]
    c <- pars[nmpar == "par3"]
    d <- pars[grep("par4\\[", nmpar)]
    s2 <- pars[grep("par5\\[", nmpar)]

  } else {
    a <- pars[nmpar == "par1"]
    #b <- pars[nmpar == "b"]
    bmd <- pars[nmpar == "par2"]
    c <- pars[nmpar == "par3"]
    d <- pars[nmpar == "par4"]
    s2 <- pars[nmpar == "par5"]
  }

  DIH <- det(-solve(optMod$hessian))
  DIH <- ifelse(DIH<0, 0, DIH)

  npar <- sum(c(length(a), length(bmd), length(c), length(d), length(s2)))

  if(covar == 'a_sigma2') {

    w1 <- sum((2*pi)^(npar/2)*sqrt(DIH)*exp(lld-min.ll)*

                # sigma2
                dnorm(s2, dataMod$data$priormuQ[5,], dataMod$data$priorSigmaQ[5,5])*
                # d
                truncnorm::dtruncnorm(d, b = dataMod$data$truncdQ, mean = dataMod$data$priormuQ[4,1],
                                      sd = dataMod$data$priorSigmaQ[4,4])*
                # a
                mc2d::dpert(a, min = dataMod$data$priorlb[1,], max = dataMod$data$priorub[1,],
                            mode = dataMod$data$priormuQ[1,], shape = dataMod$data$shape.a)*
                # BMD
                mc2d::dpert(bmd, min = dataMod$data$priorlb[2,1], max = dataMod$data$priorub[2,1],
                            mode = dataMod$data$priormuQ[2,1], shape = dataMod$data$shape.BMD)*
                # c
                mc2d::dpert(c, min = dataMod$data$priorlb[3,], max = dataMod$data$priorub[3,],
                            mode = dataMod$data$priormuQ[3,], shape = dataMod$data$shape.c))

  } else if(covar == 'BMD_d') {

    w1 <- sum((2*pi)^(npar/2)*sqrt(DIH)*exp(lld-min.ll)*

                # sigma2
                dnorm(s2, dataMod$data$priormuQ[5,1], dataMod$data$priorSigmaQ[5,5])*
                # d
                truncnorm::dtruncnorm(d, b = dataMod$data$truncdQ, mean = dataMod$data$priormuQ[4,],
                                      sd = dataMod$data$priorSigmaQ[4,4])*
                # a
                mc2d::dpert(a, min = dataMod$data$priorlb[1,1], max = dataMod$data$priorub[1,1],
                            mode = dataMod$data$priormuQ[1,1], shape = dataMod$data$shape.a)*
                # BMD
                mc2d::dpert(bmd, min = dataMod$data$priorlb[2,], max = dataMod$data$priorub[2,],
                            mode = dataMod$data$priormuQ[2,], shape = dataMod$data$shape.BMD)*
                # c
                mc2d::dpert(c, min = dataMod$data$priorlb[3,1], max = dataMod$data$priorub[3,1],
                            mode = dataMod$data$priormuQ[3,1], shape = dataMod$data$shape.c))

  } else if(covar == 'all') {

    w1 <- sum((2*pi)^(npar/2)*sqrt(DIH)*exp(lld-min.ll)*

                # sigma2
                dnorm(s2, dataMod$data$priormuQ[5,], dataMod$data$priorSigmaQ[5,5])*
                # d
                truncnorm::dtruncnorm(d, b = dataMod$data$truncdQ, mean = dataMod$data$priormuQ[4,],
                                      sd = dataMod$data$priorSigmaQ[4,4])*
                # a
                mc2d::dpert(a, min = dataMod$data$priorlb[1,], max = dataMod$data$priorub[1,],
                            mode = dataMod$data$priormuQ[1,], shape = dataMod$data$shape.a)*
                # BMD
                mc2d::dpert(bmd, min = dataMod$data$priorlb[2,], max = dataMod$data$priorub[2,],
                            mode = dataMod$data$priormuQ[2,], shape = dataMod$data$shape.BMD)*
                # c
                mc2d::dpert(c, min = dataMod$data$priorlb[3,], max = dataMod$data$priorub[3,],
                            mode = dataMod$data$priormuQ[3,], shape = dataMod$data$shape.c))

  } else {

    npar <- 5

    w1 <- (2*pi)^(npar/2)*sqrt(DIH)*exp(lld-min.ll)*

      # sigma2
      dnorm(s2, dataMod$data$priormuQ[5], dataMod$data$priorSigmaQ[5,5])*
      # d
      truncnorm::dtruncnorm(d, b = dataMod$data$truncdQ, mean = dataMod$data$priormuQ[4],
                            sd = dataMod$data$priorSigmaQ[4,4])*
      # a
      mc2d::dpert(a, min = dataMod$data$priorlb[1], max = dataMod$data$priorub[1],
                  mode = dataMod$data$priormuQ[1], shape = dataMod$data$shape.a)*
      # BMD
      mc2d::dpert(bmd, min = dataMod$data$priorlb[2], max = dataMod$data$priorub[2],
                  mode = dataMod$data$priormuQ[2], shape = dataMod$data$shape.BMD)*
      # c
      mc2d::dpert(c, min = dataMod$data$priorlb[3], max = dataMod$data$priorub[3],
                  mode = dataMod$data$priormuQ[3], shape = dataMod$data$shape.c)
  }

  return(w1)
}



#' Function to perform model selection when including a covariate effect (used internally)
#'
#' @param model stan model including covariates
#' @param model_name model name
#' @param model.none stan model without covariates
#' @param loglik function to compute log-likelihood
#' @param data_asigma2 data when covariate is a_sigma2
#' @param data_dBMD data when covariate is BMD_d
#' @param data_all data when covariate is all
#' @param data_none data when covariate is none
#' @param prior.weightsCov vector of prior weights for the 4 submodels (defaults to rep(1,4), equal weight)
#' @param pvec vector of probabilities for BMD credible interval
#' @param ndraws ndraws
#' @param td truncation for d in QE4 model
#' @param seed random seed for reproducibility
#'
#' @description This function selects the best-fitting submodel
#'
#' @return .
#'
#' @export fun_cov_selection
#'
fun_cov_selection <- function(model, model_name, model.none, loglik, data_asigma2, data_dBMD,
                              data_all, data_none, prior.weightsCov = rep(1, 4),
                              pvec, ndraws, td, seed){

  # fit all four submodels
  optMod_asigma2 <- fun_optimCov(model, data_asigma2$data,
                                 data_asigma2$start, ndraws, seed, pvec)

  if((ifelse(is.na(optMod_asigma2[[3]]),TRUE,(optMod_asigma2[[3]]!=0)) | length(optMod_asigma2)!=9)){
    prior.weightsCov[1] <- 0
    # warning('Difficulties fitting the Probit (Lognormal) model; prior weight was set to 0 and the model is not included in model averaging')
    ll_asigma2 <- NA
  }else{
    ll_asigma2 <- loglik(pars = optMod_asigma2$par, x = data_asigma2$data$x, n = data_asigma2$data$n,
                         m = data_asigma2$data$m, s2 = data_asigma2$data$s2,
                         qval = data_asigma2$data$q, shift = data_asigma2$data$shift, covar = 'a_sigma2',
                         nlevels = data_asigma2$data$nlevels,
                         trt_ind = data_asigma2$data$trt_ind)
  }


  optMod_dBMD <- fun_optimCov(model, data_dBMD$data,
                              data_dBMD$start, ndraws, seed, pvec)

  if((ifelse(is.na(optMod_dBMD[[3]]),TRUE,(optMod_dBMD[[3]]!=0)) | length(optMod_dBMD)!=9)){
    prior.weightsCov[2] <- 0
    # warning('Difficulties fitting the Probit (Lognormal) model; prior weight was set to 0 and the model is not included in model averaging')
    ll_dBMD <- NA
  }else{
    ll_dBMD <- loglik(pars = optMod_dBMD$par, x = data_dBMD$data$x, n = data_dBMD$data$n,
                      m = data_dBMD$data$m, s2 = data_dBMD$data$s2,
                      qval = data_dBMD$data$q, shift = data_dBMD$data$shift, covar = 'BMD_d',
                      nlevels = data_dBMD$data$nlevels,
                      trt_ind = data_dBMD$data$trt_ind)
  }

  optMod_all <- fun_optimCov(model, data_all$data,
                             data_all$start, ndraws, seed, pvec)

  if((ifelse(is.na(optMod_all[[3]]),TRUE,(optMod_all[[3]]!=0)) | length(optMod_all)!=9)){
    prior.weightsCov[3] <- 0
    # warning('Difficulties fitting the Probit (Lognormal) model; prior weight was set to 0 and the model is not included in model averaging')
    ll_all <- NA
  }else{
    ll_all <- loglik(pars = optMod_all$par, x = data_all$data$x, n = data_all$data$n,
                     m = data_all$data$m, s2 = data_all$data$s2,
                     qval = data_all$data$q, shift = data_all$data$shift, covar = 'all',
                     nlevels = data_all$data$nlevels,
                     trt_ind = data_all$data$trt_ind)
  }

  optMod_none <- fun_optim(model.none, data_none$data,
                           data_none$start, ndraws, seed, pvec)

  if((ifelse(is.na(optMod_none[[3]]),TRUE,(optMod_none[[3]]!=0)) | length(optMod_none)!=9)){
    prior.weightsCov[4] <- 0
    # warning('Difficulties fitting the Probit (Lognormal) model; prior weight was set to 0 and the model is not included in model averaging')
    ll_none <- NA
  }else{
    ll_none <- loglik(pars = optMod_none$par, x = data_none$data$x, n = data_none$data$n,
                      m = data_none$data$m, s2 = data_none$data$s2,
                      qval = data_none$data$q, shift = data_none$data$shift, covar = 'none',
                      nlevels = 1,
                      trt_ind = data_none$data$trt_ind)
  }

  lls <- c(ll_asigma2, ll_dBMD, ll_all, ll_none)
  if(sum(is.na(lls)) == 4){
    # print(warning(paste0('Problems fitting model ', model_name, ', this model was excluded from the analysis.')))
    warning(paste0('Problems fitting model ', model_name, ', this model was excluded from the analysis.'))
    return(NULL)
  }else{
    max.ll = max(lls, na.rm = T)
    if(is.na(lls[which((max.ll-lls[!is.na(lls)]) < 709)][1])){
      lpw.sub <- rep(0, 4)
      lpw.sub[which(lls == max.ll)] <- 1
    }else{
      minll <- min(lls[which((max.ll-lls[!is.na(lls)]) < 709)], na.rm = T)

      # compute the weight of each model
      if(prior.weightsCov[1] > 0){
        if(grepl('QE4', model_name)){
          w_asigma2 <- fun.w.QE4(optMod_asigma2, lld = ll_asigma2, min.ll = minll, nlevels = data_asigma2$data$nlevels,
                                 dataMod = data_asigma2, covar = 'a_sigma2')
        }else{
          w_asigma2 <- fun.w(optMod_asigma2, lld = ll_asigma2, min.ll = minll, nlevels = data_asigma2$data$nlevels,
                             dataMod = data_asigma2, covar = 'a_sigma2')
        }
      }else{w_asigma2 <- 0}

      if(prior.weightsCov[2] > 0){
        if(grepl('QE4', model_name)){
          w_dBMD <- fun.w.QE4(optMod_dBMD, lld = ll_dBMD, min.ll = minll, nlevels = data_dBMD$data$nlevels,
                              dataMod = data_dBMD, covar = 'BMD_d')
        }
        w_dBMD <- fun.w(optMod_dBMD, lld = ll_dBMD, min.ll = minll, nlevels = data_dBMD$data$nlevels,
                        dataMod = data_dBMD, covar = 'BMD_d')
      }else{w_dBMD <- 0}

      if(prior.weightsCov[3] > 0){
        if(grepl('QE4', model_name)){
          w_all <- fun.w.QE4(optMod_all, lld = ll_all, min.ll = minll, nlevels = data_all$data$nlevels,
                             dataMod = data_all, covar = 'all')
        }
        w_all <- fun.w(optMod_all, lld = ll_all, min.ll = minll, nlevels = data_all$data$nlevels,
                       dataMod = data_all, covar = 'all')
      }else{w_all <- 0}

      if(prior.weightsCov[4] > 0){
        if(grepl('QE4', model_name)){
          w_none <- fun.w.QE4(optMod_none, lld = ll_none, min.ll = minll, nlevels = 1,
                              dataMod = data_none, covar = 'none')
        }
        w_none <- fun.w(optMod_none, lld = ll_none, min.ll = minll, nlevels = 1,
                        dataMod = data_none, covar = 'none')
      }else{w_none <- 0}

      w <- c(w_asigma2, w_dBMD, w_all, w_none)
      w <- ifelse(w == 'Inf' | is.na(w), 0, w)
      prior.weightsCov = prior.weightsCov/sum(prior.weightsCov==1)
      lpw.sub=(prior.weightsCov*w)/sum(prior.weightsCov*w)
    }

    if(max(lpw.sub) == lpw.sub[4]){
      best.submodel <- optMod_none
      best.submodel.name <- 'none'
    }else if(max(lpw.sub) == lpw.sub[1]){
      best.submodel <- optMod_asigma2
      best.submodel.name <- 'a_sigma2'
    }else if(max(lpw.sub) == lpw.sub[2]){
      best.submodel <- optMod_dBMD
      best.submodel.name <- 'BMD_d'
    }else if(max(lpw.sub) == lpw.sub[3]){
      best.submodel <- optMod_all
      best.submodel.name <- 'all'
    }

    return(list(Weights = data.frame(Covariate = c('a_sigma2', 'BMD_d', 'all', 'none'), Loglik = lls,
                                     # BIC = c(BIC.asigma2, BIC.dBMD, BIC.all, BIC.none)
                                     Weights = lpw.sub),
                a_sigma2 = optMod_asigma2,
                BMD_d = optMod_dBMD,
                all = optMod_all,
                none = optMod_none,
                best.submodel = best.submodel,
                best.submodel.name = best.submodel.name

    ))
  }
}

#' @rdname fun.w
#' @export
fun.wQ <- function(optMod, lld, min.ll, nlevels,
                   dataMod,
                   covar = c('background', 'BMD_d', 'all', 'none')) {

  covar = match.arg(covar, c('background', 'BMD_d', 'all', 'none'))
  pars <- optMod$par
  nmpar <- names(pars)

  if(covar != 'none'){

    a <- pars[grep("par1\\[", nmpar)]
    bmd <- pars[grep("par2\\[", nmpar)]
    d <- pars[grep("par3\\[", nmpar)]

  } else {
    a <- pars[nmpar == "par1"]
    bmd <- pars[nmpar == "par2"]
    d <- pars[nmpar == "par3"]
  }

  DIH <- det(-solve(optMod$hessian))
  DIH <- ifelse(DIH<0, 0, DIH)

  npar <- sum(c(length(a), length(bmd), length(d)))

  if(covar == 'background') {

    w1 <- sum((2*pi)^(npar/2)*sqrt(DIH)*exp(lld-min.ll)*
                # d
                truncnorm::dtruncnorm(d, b = dataMod$data$truncd, mean = dataMod$data$priormu[3,1],
                                      sd = dataMod$data$priorSigma[3,3])*
                # a
                mc2d::dpert(a, min = dataMod$data$priorlb[1,], max = dataMod$data$priorub[1,],
                            mode = dataMod$data$priormu[1,], shape = dataMod$data$priorgama[1, ])*
                # BMD
                mc2d::dpert(bmd, min = dataMod$data$priorlb[2,1], max = dataMod$data$priorub[2,1],
                            mode = dataMod$data$priormu[2,1], shape = dataMod$data$priorgama[2,1]))

  } else if(covar == 'BMD_d') {

    w1 <- sum((2*pi)^(npar/2)*sqrt(DIH)*exp(lld-min.ll)*
                # d
                truncnorm::dtruncnorm(d, b = dataMod$data$truncd, mean = dataMod$data$priormu[3,],
                                      sd = dataMod$data$priorSigma[3,3])*
                # a
                mc2d::dpert(a, min = dataMod$data$priorlb[1,1], max = dataMod$data$priorub[1,1],
                            mode = dataMod$data$priormu[1,1], shape = dataMod$data$priorgama[1,1])*
                # BMD
                mc2d::dpert(bmd, min = dataMod$data$priorlb[2,], max = dataMod$data$priorub[2,],
                            mode = dataMod$data$priormu[2,], shape = dataMod$data$priorgama[2, ]))

  } else if(covar == 'all') {

    w1 <- sum((2*pi)^(npar/2)*sqrt(DIH)*exp(lld-min.ll)*
                # d
                truncnorm::dtruncnorm(d, b = dataMod$data$truncd, mean = dataMod$data$priormu[3,],
                                      sd = dataMod$data$priorSigma[3,3])*
                # a
                mc2d::dpert(a, min = dataMod$data$priorlb[1,], max = dataMod$data$priorub[1,],
                            mode = dataMod$data$priormu[1,], shape = dataMod$data$priorgama[1,])*
                # BMD
                mc2d::dpert(bmd, min = dataMod$data$priorlb[2,], max = dataMod$data$priorub[2,],
                            mode = dataMod$data$priormu[2,], shape = dataMod$data$priorgama[2,]))

  } else {

    npar <- 3

    w1 <- (2*pi)^(npar/2)*sqrt(DIH)*exp(lld-min.ll)*
      # d
      truncnorm::dtruncnorm(d, b = dataMod$data$truncd, mean = dataMod$data$priormu[3],
                            sd = dataMod$data$priorSigma[3,3])*
      # a
      mc2d::dpert(a, min = dataMod$data$priorlb[1], max = dataMod$data$priorub[1],
                  mode = dataMod$data$priormu[1], shape = dataMod$data$priorgama[1])*
      # BMD
      mc2d::dpert(bmd, min = dataMod$data$priorlb[2], max = dataMod$data$priorub[2],
                  mode = dataMod$data$priormu[2], shape = dataMod$data$priorgama[2])
  }

  return(w1)
}

#' @rdname fun.w
#' @export
fun.wQ.QE4 <- function(optMod, lld, min.ll, nlevels,
                       dataMod,
                       covar = c('background', 'BMD_d', 'all', 'none')) {

  covar = match.arg(covar, c('background', 'BMD_d', 'all', 'none'))
  pars <- optMod$par
  nmpar <- names(pars)

  if(covar != 'none'){

    a <- pars[grep("par1\\[", nmpar)]
    bmd <- pars[grep("par2\\[", nmpar)]
    d <- pars[grep("par3\\[", nmpar)]

  } else {
    a <- pars[nmpar == "par1"]
    bmd <- pars[nmpar == "par2"]
    d <- pars[nmpar == "par3"]
  }

  DIH <- det(-solve(optMod$hessian))
  DIH <- ifelse(DIH<0, 0, DIH)

  npar <- sum(c(length(a), length(bmd), length(d)))

  if(covar == 'background') {

    w1 <- sum((2*pi)^(npar/2)*sqrt(DIH)*exp(lld-min.ll)*
                # d
                truncnorm::dtruncnorm(d, b = dataMod$data$truncdQ, mean = dataMod$data$priormuQ[3,1],
                                      sd = dataMod$data$priorSigmaQ[3,3])*
                # a
                mc2d::dpert(a, min = dataMod$data$priorlb[1,], max = dataMod$data$priorub[1,],
                            mode = dataMod$data$priormuQ[1,], shape = dataMod$data$priorgama[1,])*
                # BMD
                mc2d::dpert(bmd, min = dataMod$data$priorlb[2,1], max = dataMod$data$priorub[2,1],
                            mode = dataMod$data$priormuQ[2,1], shape = dataMod$data$priorgama[2,1]))

  } else if(covar == 'BMD_d') {

    w1 <- sum((2*pi)^(npar/2)*sqrt(DIH)*exp(lld-min.ll)*
                # d
                truncnorm::dtruncnorm(d, b = dataMod$data$truncdQ, mean = dataMod$data$priormuQ[3,],
                                      sd = dataMod$data$priorSigmaQ[3,3])*
                # a
                mc2d::dpert(a, min = dataMod$data$priorlb[1,1], max = dataMod$data$priorub[1,1],
                            mode = dataMod$data$priormuQ[1,1], shape = dataMod$data$priorgama[1,1])*
                # BMD
                mc2d::dpert(bmd, min = dataMod$data$priorlb[2,], max = dataMod$data$priorub[2,],
                            mode = dataMod$data$priormuQ[2,], shape = dataMod$data$priorgama[2,]))

  } else if(covar == 'all') {

    w1 <- sum((2*pi)^(npar/2)*sqrt(DIH)*exp(lld-min.ll)*
                # d
                truncnorm::dtruncnorm(d, b = dataMod$data$truncdQ, mean = dataMod$data$priormuQ[3,],
                                      sd = dataMod$data$priorSigmaQ[3,3])*
                # a
                mc2d::dpert(a, min = dataMod$data$priorlb[1,], max = dataMod$data$priorub[1,],
                            mode = dataMod$data$priormuQ[1,], shape = dataMod$data$priorgama[1,])*
                # BMD
                mc2d::dpert(bmd, min = dataMod$data$priorlb[2,], max = dataMod$data$priorub[2,],
                            mode = dataMod$data$priormuQ[2,], shape = dataMod$data$priorgama[2,]))

  } else {

    npar <- 3

    w1 <- (2*pi)^(npar/2)*sqrt(DIH)*exp(lld-min.ll)*
      # d
      truncnorm::dtruncnorm(d, b = dataMod$data$truncdQ, mean = dataMod$data$priormuQ[3],
                            sd = dataMod$data$priorSigmaQ[3,3])*
      # a
      mc2d::dpert(a, min = dataMod$data$priorlb[1], max = dataMod$data$priorub[1],
                  mode = dataMod$data$priormuQ[1], shape = dataMod$data$priorgama[1])*
      # BMD
      mc2d::dpert(bmd, min = dataMod$data$priorlb[2], max = dataMod$data$priorub[2],
                  mode = dataMod$data$priormuQ[2], shape = dataMod$data$priorgama[2])
  }

  return(w1)
}

#' @rdname fun_cov_selection
#' @export
fun_cov_selectionQ <- function(model, model_name, model.none, loglik, data_bkg, data_dBMD,
                               data_all, data_none, prior.weightsCov = rep(1, 4),
                               pvec, ndraws, td, seed){

  # fit all four submodels
  optMod_bkg <- fun_optimQCov(model, data_bkg$data,
                              data_bkg$start, ndraws, seed, pvec)

  if((ifelse(is.na(optMod_bkg[[3]]),TRUE,(optMod_bkg[[3]]!=0)) | length(optMod_bkg)!=9)){
    prior.weightsCov[1] <- 0
    # warning('Difficulties fitting the Probit (Lognormal) model; prior weight was set to 0 and the model is not included in model averaging')
    ll_bkg <- NA
  }else{
    ll_bkg <- loglik(pars = optMod_bkg$par, x = data_bkg$data$x,
                     n = data_bkg$data$n, y = data_bkg$data$y,
                     qval = data_bkg$data$q, covar = 'background',
                     nlevels = data_bkg$data$nlevels,
                     trt_ind = data_bkg$data$trt_ind)
  }

  optMod_dBMD <- fun_optimQCov(model, data_dBMD$data,
                               data_dBMD$start, ndraws, seed, pvec)

  if((ifelse(is.na(optMod_dBMD[[3]]),TRUE,(optMod_dBMD[[3]]!=0)) | length(optMod_dBMD)!=9)){
    prior.weightsCov[1] <- 0
    # warning('Difficulties fitting the Probit (Lognormal) model; prior weight was set to 0 and the model is not included in model averaging')
    ll_dBMD <- NA
  }else{
    ll_dBMD <- loglik(pars = optMod_dBMD$par, x = data_dBMD$data$x,
                      n = data_dBMD$data$n, y = data_dBMD$data$y,
                      qval = data_dBMD$data$q, covar = 'BMD_d',
                      nlevels = data_dBMD$data$nlevels,
                      trt_ind = data_dBMD$data$trt_ind)
  }

  optMod_all <- fun_optimQCov(model, data_all$data,
                              data_all$start, ndraws, seed, pvec)

  if((ifelse(is.na(optMod_all[[3]]),TRUE,(optMod_all[[3]]!=0)) | length(optMod_all)!=9)){
    prior.weightsCov[1] <- 0
    # warning('Difficulties fitting the Probit (Lognormal) model; prior weight was set to 0 and the model is not included in model averaging')
    ll_all <- NA
  }else{
    ll_all <- loglik(pars = optMod_all$par, x = data_all$data$x,
                     n = data_all$data$n, y = data_all$data$y,
                     qval = data_all$data$q, covar = 'all',
                     nlevels = data_all$data$nlevels,
                     trt_ind = data_all$data$trt_ind)
  }

  optMod_none <- fun_optimQ(model.none, data_none$data,
                            data_none$start, ndraws, seed, pvec)

  if((ifelse(is.na(optMod_none[[3]]),TRUE,(optMod_none[[3]]!=0)) | length(optMod_none)!=9)){
    prior.weightsCov[1] <- 0
    # warning('Difficulties fitting the Probit (Lognormal) model; prior weight was set to 0 and the model is not included in model averaging')
    ll_none <- NA
  }else{
    ll_none <- loglik(pars = optMod_none$par, x = data_none$data$x,
                      n = data_none$data$n, y = data_none$data$y,
                      qval = data_none$data$q, covar = 'none',
                      nlevels = 1,
                      trt_ind = data_none$data$trt_ind)
  }

  lls <- c(ll_bkg, ll_dBMD, ll_all, ll_none)
  if(sum(is.na(lls)) == 4){
    # print(warning(paste0('Problems fitting model ', model_name, ', this model was excluded from the analysis.')))
    warning(paste0('Problems fitting model ', model_name, ', this model was excluded from the analysis.'))
    return(NULL)
  }else{
    max.ll = max(lls, na.rm = T)
    if(is.na(lls[which((max.ll-lls[!is.na(lls)]) < 709)][1])){
      lpw.sub <- rep(0, 4)
      lpw.sub[which(lls == max.ll)] <- 1
    }else{
      minll <- min(lls[which((max.ll-lls[!is.na(lls)]) < 709)], na.rm = T)

      # compute the weight of each model
      if(prior.weightsCov[1] > 0){
        if(grepl('QE4', model_name)){
          w_bkg <- fun.wQ.QE4(optMod_bkg, lld = ll_bkg, min.ll = minll, nlevels = data_bkg$data$nlevels,
                              dataMod = data_bkg, covar = 'background')
        }else{
          w_bkg <- fun.wQ(optMod_bkg, lld = ll_bkg, min.ll = minll, nlevels = data_bkg$data$nlevels,
                          dataMod = data_bkg, covar = 'background')
        }
      }else{w_bkg <- 0}

      if(prior.weightsCov[2] > 0){
        if(grepl('QE4', model_name)){
          w_dBMD <- fun.wQ.QE4(optMod_dBMD, lld = ll_dBMD, min.ll = minll, nlevels = data_dBMD$data$nlevels,
                               dataMod = data_dBMD, covar = 'BMD_d')
        }
        w_dBMD <- fun.wQ(optMod_dBMD, lld = ll_dBMD, min.ll = minll, nlevels = data_dBMD$data$nlevels,
                         dataMod = data_dBMD, covar = 'BMD_d')
      }else{w_dBMD <- 0}

      if(prior.weightsCov[3] > 0){
        if(grepl('QE4', model_name)){
          w_all <- fun.wQ.QE4(optMod_all, lld = ll_all, min.ll = minll, nlevels = data_all$data$nlevels,
                              dataMod = data_all, covar = 'all')
        }
        w_all <- fun.wQ(optMod_all, lld = ll_all, min.ll = minll, nlevels = data_all$data$nlevels,
                        dataMod = data_all, covar = 'all')
      }else{w_all <- 0}

      if(prior.weightsCov[4] > 0){
        if(grepl('QE4', model_name)){
          w_none <- fun.wQ.QE4(optMod_none, lld = ll_none, min.ll = minll, nlevels = 1,
                               dataMod = data_none, covar = 'none')
        }
        w_none <- fun.wQ(optMod_none, lld = ll_none, min.ll = minll, nlevels = 1,
                         dataMod = data_none, covar = 'none')
      }else{w_none <- 0}

      w <- c(w_bkg, w_dBMD, w_all, w_none)
      w <- ifelse(w == 'Inf' | is.na(w), 0, w)
      prior.weightsCov = prior.weightsCov/sum(prior.weightsCov==1)
      lpw.sub=(prior.weightsCov*w)/sum(prior.weightsCov*w)
    }

    if(max(lpw.sub) == lpw.sub[4]){
      best.submodel <- optMod_none
      best.submodel.name <- 'none'
    }else if(max(lpw.sub) == lpw.sub[1]){
      best.submodel <- optMod_bkg
      best.submodel.name <- 'background'
    }else if(max(lpw.sub) == lpw.sub[2]){
      best.submodel <- optMod_dBMD
      best.submodel.name <- 'BMD_d'
    }else if(max(lpw.sub) == lpw.sub[3]){
      best.submodel <- optMod_all
      best.submodel.name <- 'all'
    }


    return(list(Weights = data.frame(Covariate = c('background', 'BMD_d', 'all', 'none'), Loglik = lls,
                                     # BIC = c(BIC.bkg, BIC.dBMD, BIC.all, BIC.none)
                                     Weights = lpw.sub),
                bkg = optMod_bkg,
                BMD_d = optMod_dBMD,
                all = optMod_all,
                none = optMod_none,
                best.submodel = best.submodel,
                best.submodel.name = best.submodel.name

    ))
  }


}


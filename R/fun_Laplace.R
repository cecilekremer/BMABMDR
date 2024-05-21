
#' Model averaging using Full Laplace method
#'
#' @param data.N the input data as returned by function PREP_DATA_N or PREP_DATA_N_C for clustered data
#' @param data.LN the input data as returned by function PREP_DATA_LN or PREP_DATA_LN_C for clustered data
#' @param data.Q the input data as returned by function PREP_DATA_Q
#' @param prior.weights a vector specifying which of the 16 (continuous) or 8 (quantal) models should be included (1 = include, 0 = exclude)
#' @param ndraws the number of draws, default 30000
#' @param seed random seed, default 123
#' @param pvec vector specifying the three BMD quantiles of interest (default 90% CrI)
#'
#' @examples
#'  data_N <- PREP_DATA_N(data = as.data.frame(immunotoxicityData[1:5,]),
#'                        sumstats = TRUE, sd = TRUE, q = 0.1)
#'  data_LN <- PREP_DATA_LN(data = as.data.frame(immunotoxicityData[1:5,]),
#'                          sumstats = TRUE, sd = TRUE, q = 0.1) #'
#'  FLBMD <- full.laplace_MA(data_N, data_LN, prior.weights = c(rep(1,4), rep(0,12)))
#'
#' @description This function performs model averaging using the full Laplace approximation. By default, all 16 models are included for continuous data, and all 8 models for quantal data.
#'              Models can be excluded by setting their respective weight to 0 in \code{prior.weights}. The order of models fitted can be obtained using the \code{\link{get_models()}} function.
#'
#' `full.laplace_MA` is used for continuous data
#'
#' `full.laplace_MAc` is used for clustered continuous data (i.e. with litter effect)
#'
#' `full.laplaceQ_MA` is used for quantal data (with or without litter effect)
#'
#' @importFrom truncnorm dtruncnorm
#'
#' @return The function returns a BMDBMA model object, which is a list containing the following objects:
#' @return `E4_N` parameter estimates from the (Exponential, Normal) model
#' @return `IE4_N` parameter estimates from the (Inverse Exponential, Normal) model
#' @return `H4_N` parameter estimates from the (Hill, Normal) model
#' @return `LN4_N` parameter estimates from the (Log-normal, Normal)
#' @return `G4_N` parameter estimates from the (Gamma, Normal) model
#' @return `QE4_N` parameter estimates from the (Quadratic Exponential, Normal) model
#' @return `P4_N` parameter estimates from the (Probit, Normal) model
#' @return `L4_N` parameter estimates from the (Logit, Normal) model
#' @return `E4_LN` parameter estimates from the (Exponential, Log-normal) model
#' @return `IE4_LN` parameter estimates from the (Inverse Exponential, Log-normal) model
#' @return `H4_LN` parameter estimates from the (Hill, Log-normal) model
#' @return `LN4_LN` parameter estimates from the (Log-normal, Log-normal)
#' @return `G4_LN` parameter estimates from the (Gamma, Log-normal) model
#' @return `QE4_LN` parameter estimates from the (Quadratic Exponential, Log-normal) model
#' @return `P4_LN` parameter estimates from the (Probit, Log-normal) model
#' @return `L4_LN` parameter estimates from the (Logit, Log-normal) model
#' @return `E4_Q` parameter estimates from the (Exponential, Quantal) model
#' @return `IE4_Q` parameter estimates from the (Inverse Exponential, Quantal) model
#' @return `H4_Q` parameter estimates from the (Hill, Quantal) model
#' @return `LN4_Q` parameter estimates from the (Lognormal, Quantal) model
#' @return `G4_Q` parameter estimates from the (Gamma, Quantal) model
#' @return `QE4_Q` parameter estimates from the (Quadratic Exponential, Quantal) model
#' @return `P4_Q` parameter estimates from the (Probit, Quantal) model
#' @return `L4_Q` parameter estimates from the (Logit, Quantal) model
#' @return `MA`  Model averaged BMD credible interval
#' @return `weights` Model weights used in the averaging
#' @return `bf` Bayes factor comparing the best model against saturated ANOVA model
#' @return `covs` matrix with covariances between parameters b-d and BMD-d
#' @return `corrs` matrix with correlation between parameters b-d and BMD-d
#' @return `p.msg` warning message if model averaged posterior has been truncated
#' @return `w.msg` warning message if Laplace weights could not be computed and one model gets all the weight
#' @return `shift` shift value for lognormal data
#' @return `BIC.SM` BIC value of saturated model, used to test for GOF
#' @return `BIC.bestfit` BIC value of best fitting model, used to test for GOF
#' @return `means.SM` mean response per dose level estimated from the saturated model, used to test for GOF
#' @return `gof_check` GOF message
#' @return `models_included_laplace` vector containing the names of models included in the model averaging
#' @return `q` BMR
#' @return `max.dose` maximum dose level (original scale)
#' @return `dataN` normal summary data used for analysis
#' @return `dataLN` lognormal summary data used for analysis
#' @return `BMDMixture` vector of length \code{ndraws} containing the draws from the model-averaged posterior
#' @return `MA_dr` vector containing the model-averaged response at each dose level
#' @return `MA_post` vector containing the 0.5%-percentiles of the model-averaged posterior
#' @return `llN` vector containing the loglikelihood values of the Normal models
#' @return `llLN` vector containing the loglikelihood values of the Lognormal models
#' @return `parsN` list containing the fitted Normal models
#' @return `parsLN` list containing the fitted Lognormal models
#' @return `is_bin` logical indicating whether binomial (no litter) model was used
#' @return `is_betabin` logical indicating whether betabinomial (with litter) model was used
#' @return `data` quantal summary data used for analysis
#' @return `parsQ` list containing the fitted Quantal models
#' @return `llQ` vector containing the loglikelihood values of the Quantal models
#'
#' @export full.laplace_MA
#'
full.laplace_MA=function(data.N, data.LN,
                         prior.weights = rep(1,16),
                         ndraws=30000,seed=123,
                         pvec=c(0.05,0.5,0.95),
                         plot=FALSE){

  # prior.weights = prior.weights/sum(prior.weights==1) #this is now done below when calculating weights

  out.stop = 'ok'

  data = data.N$data
  start = data.N$start
  startQ = data.N$startQ

  llN = c() # likelihoods
  BMDL = c()
  BMD = c()
  BMDU = c()

  ### Getting the posterior modes and the hessian to obtain the model specific posterior distributions
  # optimizing function: obtain point estimates by maximizing the joint posterior from the Stan model
  ### Additionally save BMD(L/U) estimates and loglikelihood values; set to NA if model has prior weight 0

  nModels = 16

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[1]))
  }
  if(prior.weights[1]>0){

    # print(1)

    optE4_NI <- fun_optim(stanmodels$mE4, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optE4_NI[[3]]),TRUE,(optE4_NI[[3]]!=0)) | length(optE4_NI)!=9)){
      prior.weights[1] <- 0
      warning('difficulties fitting the Exponential (Normal) model; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsE4N <- par_extract(optE4_NI, model_name = "E4_N")
      E4resNI <- quantile(parsE4N$BMD*data$maxD, pvec, na.rm = T)
      E4resNI <- c(E4resNI,apply(parsE4N[,c(paste0("p",1:4), "is2t")], 2, median))
      names(E4resNI) <- c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t")

      E4outNI <- outLP(parsE4N, pvec, data$maxD)

      # Covariance between b-d and between BMD-d
      E4covNI <- c(cov(parsE4N[,c("b","d")], use="na.or.complete")["b","d"],
                   cov(parsE4N[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      E4corrNI <- c(cor(parsE4N[,c("b","d")], use="na.or.complete")["b","d"],
                    cor(parsE4N[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      if(data$data_type == 1){
        DRM_E4_N <- DRM.E4_NI(E4resNI[4:7], data$x, data$q)
        # obtain loglikelihood
        llE4N=llfE4_NI(optE4_NI$par[c(1,2,9,4,5)],
                       nvec=data$n,
                       dvec=data$x,
                       mvec=data$m,
                       s2vec=data$s2,
                       qval=data$q)
      }else if(data$data_type == 3){
        DRM_E4_N <- DRM.E4_ND(E4resNI[4:7], data$x, data$q)
        # obtain loglikelihood
        llE4N=llfE4_ND(optE4_NI$par[c(1,2,9,4,5)],
                       nvec=data$n,
                       dvec=data$x,
                       mvec=data$m,
                       s2vec=data$s2,
                       qval=data$q)
      }
      # llE4fN=function(x) llfE4_NI(x,nvec=data$n,dvec=data$x,mvec=data$m,s2vec=data$s2,qval=data$q)+
      #   mvtnorm::dmvnorm(x[c(4,5)],mean=data$priormu[c(4,5)],
      #                    sigma=data$priorSigma[c(4,5),c(4,5)], log = T)+
      #   mc2d::dpert(x[1], min = data$priorlb[1], max = data$priorub[1],
      #               mode = data$priormu[1], shape = data$shape.a, log=T)+
      #   mc2d::dpert(x[2], min = data$priorlb[2], max = data$priorub[2],
      #               mode = data$priormu[2], shape = data$shape.BMD,log=T)+
      #   mc2d::dpert(x[3], min = data$priorlb[3], max = data$priorub[3],
      #               mode = data$priormu[3], shape = data$shape.c,log=T)


      llN = c(llN, llE4N)
      BMDL = c(BMDL, E4resNI[1])
      BMD = c(BMD, E4resNI[2])
      BMDU = c(BMDU, E4resNI[3])
    }
  }
  if(prior.weights[1] == 0){E4resNI=NULL; llN = c(llN,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA);
  E4covNI=rep(NA,2); E4corrNI=rep(NA,2); DRM_E4_N=rep(NA, length(data$x))
  parsE4N <- NA
  E4outNI <- t(data.frame(
    par.dt = rep(NA,3),
    par.is2t = rep(NA,3),
    par.k = rep(NA,3),
    par.a = rep(NA,3),
    par.b = rep(NA,3),
    par.c = rep(NA,3),
    par.d = rep(NA,3),
    par.s2 = rep(NA,3),
    BMD = rep(NA,3),
    min.resp = rep(NA,3),
    max.resp = rep(NA,3),
    fold.change = rep(NA,3)
  ))
  }

  #
  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[2]))
  }
  if(prior.weights[2]>0){

    # print(2)

    optIE4_NI <- fun_optim(stanmodels$mIE4, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optIE4_NI[[3]]),TRUE,(optIE4_NI[[3]]!=0)) | length(optIE4_NI)!=9)){
      prior.weights[2] <- 0
      warning('difficulties fitting the Inverse Exponential (Normal) model; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsIE4N <- par_extract(optIE4_NI, model_name = "IE4_N")
      IE4resNI <- quantile(parsIE4N$BMD*data$maxD, pvec, na.rm = T)
      IE4resNI <- c(IE4resNI,apply(parsIE4N[,c(paste0("p",1:4), "is2t")], 2, median))
      names(IE4resNI) <- c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t")

      IE4outNI <- outLP(parsIE4N, pvec, data$maxD)

      # Covariance between b-d and between BMD-d
      IE4covNI <- c(cov(parsIE4N[,c("b","d")], use="na.or.complete")["b","d"],
                    cov(parsIE4N[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      IE4corrNI <- c(cor(parsIE4N[,c("b","d")], use="na.or.complete")["b","d"],
                     cor(parsIE4N[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      if(data$data_type == 1){
        DRM_IE4_N <- DRM.IE4_NI(IE4resNI[4:7], data$x, data$q)
        # obtain loglikelihood
        llIE4N=llfIE4_NI(optIE4_NI$par[c(1,2,9,4,5)],
                         nvec=data$n,
                         dvec=data$x,
                         mvec=data$m,
                         s2vec=data$s2,
                         qval=data$q)
      }else if(data$data_type == 3){
        DRM_IE4_N <- DRM.IE4_ND(IE4resNI[4:7], data$x, data$q)
        # obtain loglikelihood
        llIE4N=llfIE4_ND(optIE4_NI$par[c(1,2,9,4,5)],
                         nvec=data$n,
                         dvec=data$x,
                         mvec=data$m,
                         s2vec=data$s2,
                         qval=data$q)
      }

      llN = c(llN, llIE4N)
      BMDL = c(BMDL, IE4resNI[1])
      BMD = c(BMD, IE4resNI[2])
      BMDU = c(BMDU, IE4resNI[3])
    }
  }
  if(prior.weights[2] == 0){
    IE4resNI=NULL; llN = c(llN,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA);
    IE4covNI=rep(NA,2); IE4corrNI=rep(NA,2); DRM_IE4_N=rep(NA,length(data$x))
    parsIE4N <- NA
    IE4outNI <- t(data.frame(
      par.dt = rep(NA,3),
      par.is2t = rep(NA,3),
      par.k = rep(NA,3),
      par.a = rep(NA,3),
      par.b = rep(NA,3),
      par.c = rep(NA,3),
      par.d = rep(NA,3),
      par.s2 = rep(NA,3),
      BMD = rep(NA,3),
      min.resp = rep(NA,3),
      max.resp = rep(NA,3),
      fold.change = rep(NA,3)
    ))
  }
  #
  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[3]))
  }
  if(prior.weights[3]>0){

    # print(3)

    optH4_NI <- fun_optim(stanmodels$mH4, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optH4_NI[[3]]),TRUE,(optH4_NI[[3]]!=0)) | length(optH4_NI)!=9)){
      prior.weights[3] <- 0
      warning('difficulties fitting the Hill (Normal) model; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsH4N <- par_extract(optH4_NI, model_name = "H4_N")
      H4resNI <- quantile(parsH4N$BMD*data$maxD, pvec, na.rm = T)  #exp(quantile(optE4_NI$theta_tilde[,"par[2]"],pvec))
      H4resNI <- c(H4resNI,apply(parsH4N[,c(paste0("p",1:4), "is2t")], 2, median))
      names(H4resNI) <- c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t")

      H4outNI <- outLP(parsH4N, pvec, data$maxD)

      # Covariance between b-d and between BMD-d
      H4covNI = c(cov(parsH4N[,c("b","d")], use="na.or.complete")["b","d"],
                  cov(parsH4N[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      H4corrNI = c(cor(parsH4N[,c("b","d")], use="na.or.complete")["b","d"],
                   cor(parsH4N[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      if(data$data_type == 1){
        DRM_H4_N <- DRM.H4_NI(H4resNI[4:7], data$x, data$q)
        # obtain loglikelihood
        llH4N=llfH4_NI(optH4_NI$par[c(1,2,9,4,5)],
                       nvec=data$n,
                       dvec=data$x,
                       mvec=data$m,
                       s2vec=data$s2,
                       qval=data$q)
      }else if(data$data_type == 3){
        DRM_H4_N <- DRM.H4_ND(H4resNI[4:7], data$x, data$q)
        # obtain loglikelihood
        llH4N=llfH4_ND(optH4_NI$par[c(1,2,9,4,5)],
                       nvec=data$n,
                       dvec=data$x,
                       mvec=data$m,
                       s2vec=data$s2,
                       qval=data$q)
      }

      llN = c(llN, llH4N)
      BMDL = c(BMDL, H4resNI[1])
      BMD = c(BMD, H4resNI[2])
      BMDU = c(BMDU, H4resNI[3])
    }
  }
  if(prior.weights[3] == 0){H4resNI=NULL; llN = c(llN,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA);
  BMDU=c(BMDU,NA); H4covNI=rep(NA,2); H4corrNI=rep(NA,2); DRM_H4_N=rep(NA,length(data$x))
  parsH4N <- NA
  H4outNI <- t(data.frame(
    par.dt = rep(NA,3),
    par.is2t = rep(NA,3),
    par.k = rep(NA,3),
    par.a = rep(NA,3),
    par.b = rep(NA,3),
    par.c = rep(NA,3),
    par.d = rep(NA,3),
    par.s2 = rep(NA,3),
    BMD = rep(NA,3),
    min.resp = rep(NA,3),
    max.resp = rep(NA,3),
    fold.change = rep(NA,3)
  ))
  }

  #
  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[4]))
  }
  if(prior.weights[4]>0){

    # print(4)

    optLN4_NI <- fun_optim(stanmodels$mLN4, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optLN4_NI[[3]]),TRUE,(optLN4_NI[[3]]!=0)) | length(optLN4_NI)!=9)){
      prior.weights[4] <- 0
      warning('difficulties fitting the Lognormal (Normal) model; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsLN4N <- par_extract(optLN4_NI, model_name = "LN4_N")
      LN4resNI <- quantile(parsLN4N$BMD*data$maxD, pvec, na.rm = T)
      LN4resNI <- c(LN4resNI,apply(parsLN4N[,c(paste0("p",1:4), "is2t")], 2, median))
      names(LN4resNI) <- c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t")

      LN4outNI <- outLP(parsLN4N, pvec, data$maxD)

      # Covariance between b-d and between BMD-d
      LN4covNI = c(cov(parsLN4N[,c("b","d")], use="na.or.complete")["b","d"],
                   cov(parsLN4N[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      LN4corrNI = c(cor(parsLN4N[,c("b","d")], use="na.or.complete")["b","d"],
                    cor(parsLN4N[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      if(data$data_type == 1){
        DRM_LN4_N <- DRM.LN4_NI(LN4resNI[4:7], data$x, data$q)
        # obtain loglikelihood
        llLN4N=llfLN4_NI(optLN4_NI$par[c(1,2,9,4,5)],
                         nvec=data$n,
                         dvec=data$x,
                         mvec=data$m,
                         s2vec=data$s2,
                         qval=data$q)
      }else if(data$data_type == 3){
        DRM_LN4_N <- DRM.LN4_ND(LN4resNI[4:7], data$x, data$q)
        # obtain loglikelihood
        llLN4N=llfLN4_ND(optLN4_NI$par[c(1,2,9,4,5)],
                         nvec=data$n,
                         dvec=data$x,
                         mvec=data$m,
                         s2vec=data$s2,
                         qval=data$q)
      }

      llN = c(llN, llLN4N)
      BMDL = c(BMDL, LN4resNI[1])
      BMD = c(BMD, LN4resNI[2])
      BMDU = c(BMDU, LN4resNI[3])
    }
  }
  if(prior.weights[4] == 0){LN4resNI=NULL; llN = c(llN,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA);
  LN4covNI=rep(NA,2); LN4corrNI=rep(NA,2); DRM_LN4_N=rep(NA,length(data$x))
  parsLN4N <- NA
  LN4outNI <- t(data.frame(
    par.dt = rep(NA,3),
    par.is2t = rep(NA,3),
    par.k = rep(NA,3),
    par.a = rep(NA,3),
    par.b = rep(NA,3),
    par.c = rep(NA,3),
    par.d = rep(NA,3),
    par.s2 = rep(NA,3),
    BMD = rep(NA,3),
    min.resp = rep(NA,3),
    max.resp = rep(NA,3),
    fold.change = rep(NA,3)
  ))
  }
  #
  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[5]))
  }
  if(prior.weights[5]>0){
    # print(5)

    # if(data$is_increasing == 1){
    #   data$init_b = qgamma(data$q/
    #                          (optE4_NI$par[8]-1), rate=1.0, shape=optE4_NI$par[10])/optE4_NI$par[2]
    # }else if(data$is_decreasing == 1){
    #   data$init_b = qgamma((-data$q)/
    #                          (optE4_NI$par[8]-1), rate=1.0, shape=optE4_NI$par[10])/optE4_NI$par[2]
    # }

    optG4_NI <- fun_optim(stanmodels$mG4, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optG4_NI[[3]]),TRUE,(optG4_NI[[3]]!=0)) | length(optG4_NI)!=9)){
      prior.weights[5] <- 0
      warning('difficulties fitting the Gamma (Normal) model; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsG4N <- par_extract(optG4_NI, model_name = "G4_N")
      G4resNI <- quantile(parsG4N$BMD*data$maxD, pvec, na.rm = T)
      G4resNI <- c(G4resNI,apply(parsG4N[,c(paste0("p",1:4), "is2t")], 2, median))
      names(G4resNI) <- c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t")

      G4outNI <- outLP(parsG4N, pvec, data$maxD)

      # Covariance between b-d and between BMD-d
      G4covNI = c(cov(parsG4N[,c("b","d")], use="na.or.complete")["b","d"],
                  cov(parsG4N[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      G4corrNI = c(cor(parsG4N[,c("b","d")], use="na.or.complete")["b","d"],
                   cor(parsG4N[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      if(data$data_type == 1){
        DRM_G4_N <- DRM.G4_NI(G4resNI[4:7], data$x, data$q)
        # obtain loglikelihood
        llG4N=llfG4_NI(optG4_NI$par[c(1,2,9,4,5)],
                       nvec=data$n,
                       dvec=data$x,
                       mvec=data$m,
                       s2vec=data$s2,
                       qval=data$q)
      }else if(data$data_type == 3){
        DRM_G4_N <- DRM.G4_ND(G4resNI[4:7], data$x, data$q)
        # obtain loglikelihood
        llG4N=llfG4_ND(optG4_NI$par[c(1,2,9,4,5)],
                       nvec=data$n,
                       dvec=data$x,
                       mvec=data$m,
                       s2vec=data$s2,
                       qval=data$q)
      }

      llN = c(llN, llG4N)
      BMDL = c(BMDL, G4resNI[1])
      BMD = c(BMD, G4resNI[2])
      BMDU = c(BMDU, G4resNI[3])
    }
  }
  if(prior.weights[5] == 0){G4resNI=NULL; llN = c(llN,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA);
  BMDU=c(BMDU,NA); G4covNI=rep(NA,2); G4corrNI=rep(NA,2); DRM_G4_N=rep(NA, length(data$x))
  parsG4N <-NA
  G4outNI <- t(data.frame(
    par.dt = rep(NA,3),
    par.is2t = rep(NA,3),
    par.k = rep(NA,3),
    par.a = rep(NA,3),
    par.b = rep(NA,3),
    par.c = rep(NA,3),
    par.d = rep(NA,3),
    par.s2 = rep(NA,3),
    BMD = rep(NA,3),
    min.resp = rep(NA,3),
    max.resp = rep(NA,3),
    fold.change = rep(NA,3)
  ))
  }
  #
  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[6]))
  }
  if(prior.weights[6]>0){
    # print(6)

    optQE4_NI <- fun_optim(stanmodels$mQE4, data, startQ, ndraws, 123, pvec)

    if((ifelse(is.na(optQE4_NI[[3]]),TRUE,(optQE4_NI[[3]]!=0)) | length(optQE4_NI)!=9)){
      prior.weights[6] <- 0
      warning('difficulties fitting the Quadratic Exponential (Normal) model; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsQE4N <- par_extract(optQE4_NI, model_name = "QE4_N")
      QE4resNI <- quantile(parsQE4N$BMD*data$maxD, pvec, na.rm = T)
      QE4resNI <- c(QE4resNI,apply(parsQE4N[,c(paste0("p",1:4), "is2t")], 2, median))
      names(QE4resNI) <- c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t")

      QE4outNI <- outLP(parsQE4N, pvec, data$maxD)

      # Covariance between b-d and between BMD-d
      QE4covNI = c(cov(parsQE4N[,c("b","d")], use="na.or.complete")["b","d"],
                   cov(parsQE4N[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      QE4corrNI = c(cor(parsQE4N[,c("b","d")], use="na.or.complete")["b","d"],
                    cor(parsQE4N[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      if(data$data_type == 1){
        DRM_QE4_N <- DRM.QE4_NI(QE4resNI[4:7], data$x, data$q)
        # obtain loglikelihood
        llQE4N=llfQE4_NI(optQE4_NI$par[c(1,2,9,4,5)],
                         nvec=data$n,
                         dvec=data$x,
                         mvec=data$m,
                         s2vec=data$s2,
                         qval=data$q)
      }else if(data$data_type == 3){
        DRM_QE4_N <- DRM.QE4_ND(QE4resNI[4:7], data$x, data$q)
        # obtain loglikelihood
        llQE4N=llfQE4_ND(optQE4_NI$par[c(1,2,9,4,5)],
                         nvec=data$n,
                         dvec=data$x,
                         mvec=data$m,
                         s2vec=data$s2,
                         qval=data$q)
      }

      llN = c(llN, llQE4N)
      BMDL = c(BMDL, QE4resNI[1])
      BMD = c(BMD, QE4resNI[2])
      BMDU = c(BMDU, QE4resNI[3])
    }
  }
  if(prior.weights[6] == 0){QE4resNI=NULL; llN = c(llN,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA);
  BMDU=c(BMDU,NA); QE4covNI=rep(NA,2); QE4corrNI=rep(NA,2); DRM_QE4_N=rep(NA,length(data$x))
  parsQE4N <- NA
  QE4outNI <- t(data.frame(
    par.dt = rep(NA,3),
    par.is2t = rep(NA,3),
    par.k = rep(NA,3),
    par.a = rep(NA,3),
    par.b = rep(NA,3),
    par.c = rep(NA,3),
    par.d = rep(NA,3),
    par.s2 = rep(NA,3),
    BMD = rep(NA,3),
    min.resp = rep(NA,3),
    max.resp = rep(NA,3),
    fold.change = rep(NA,3)
  ))
  }
  #
  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[7]))
  }
  if(prior.weights[7]>0){
    # print(7)

    optP4_NI <- fun_optim(stanmodels$mP4, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optP4_NI[[3]]),TRUE,(optP4_NI[[3]]!=0)) | length(optP4_NI)!=9)){
      prior.weights[7] <- 0
      warning('difficulties fitting the Probit (Normal) model; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsP4N <- par_extract(optP4_NI, model_name = "P4_N")
      P4resNI <- quantile(parsP4N$BMD*data$maxD, pvec, na.rm = T)
      P4resNI <- c(P4resNI,apply(parsP4N[,c(paste0("p",1:4), "is2t")], 2, median))
      names(P4resNI) <- c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t")

      P4outNI <- outLP(parsP4N, pvec, data$maxD)

      # Covariance between b-d and between BMD-d
      P4covNI = c(cov(parsP4N[,c("b","d")], use="na.or.complete")["b","d"],
                  cov(parsP4N[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      P4corrNI = c(cor(parsP4N[,c("b","d")], use="na.or.complete")["b","d"],
                   cor(parsP4N[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      if(data$data_type == 1){
        DRM_P4_N <- DRM.P4_NI(P4resNI[4:7], data$x, data$q)
        # obtain loglikelihood
        llP4N=llfP4_NI(optP4_NI$par[c(1,2,9,4,5)],
                       nvec=data$n,
                       dvec=data$x,
                       mvec=data$m,
                       s2vec=data$s2,
                       qval=data$q)
      }else if(data$data_type == 3){
        DRM_P4_N <- DRM.P4_ND(P4resNI[4:7], data$x, data$q)
        # obtain loglikelihood
        llP4N=llfP4_ND(optP4_NI$par[c(1,2,9,4,5)],
                       nvec=data$n,
                       dvec=data$x,
                       mvec=data$m,
                       s2vec=data$s2,
                       qval=data$q)
      }


      llN = c(llN, llP4N)
      BMDL = c(BMDL, P4resNI[1])
      BMD = c(BMD, P4resNI[2])
      BMDU = c(BMDU, P4resNI[3])
    }
  }
  if(prior.weights[7] == 0){P4resNI=NULL; llN = c(llN,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA);
  P4covNI=rep(NA,2); P4corrNI=rep(NA,2); DRM_P4_N=rep(NA,length(data$x)); parsP4N <- NA
  P4outNI <- t(data.frame(
    par.dt = rep(NA,3),
    par.is2t = rep(NA,3),
    par.k = rep(NA,3),
    par.a = rep(NA,3),
    par.b = rep(NA,3),
    par.c = rep(NA,3),
    par.d = rep(NA,3),
    par.s2 = rep(NA,3),
    BMD = rep(NA,3),
    min.resp = rep(NA,3),
    max.resp = rep(NA,3),
    fold.change = rep(NA,3)
  ))
  }
  #
  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[8]))
  }
  if(prior.weights[8]>0){
    # print(8)

    optL4_NI <- fun_optim(stanmodels$mL4, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optL4_NI[[3]]),TRUE,(optL4_NI[[3]]!=0)) | length(optL4_NI)!=9)){
      prior.weights[8] <- 0
      warning('difficulties fitting the Logit (Normal) model; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsL4N <- par_extract(optL4_NI, model_name = "L4_N")
      L4resNI <- quantile(parsL4N$BMD*data$maxD, pvec, na.rm = T)
      L4resNI <- c(L4resNI,apply(parsL4N[,c(paste0("p",1:4), "is2t")], 2, median))
      names(L4resNI) <- c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t")

      L4outNI <- outLP(parsL4N, pvec, data$maxD)


      # Covariance between b-d and between BMD-d
      L4covNI = c(cov(parsL4N[,c("b","d")], use="na.or.complete")["b","d"],
                  cov(parsL4N[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      L4corrNI = c(cor(parsL4N[,c("b","d")], use="na.or.complete")["b","d"],
                   cor(parsL4N[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      if(data$data_type == 1){
        DRM_L4_N <- DRM.L4_NI(L4resNI[4:7], data$x, data$q)
        # obtain loglikelihood
        llL4N=llfL4_NI(optL4_NI$par[c(1,2,9,4,5)],
                       nvec=data$n,
                       dvec=data$x,
                       mvec=data$m,
                       s2vec=data$s2,
                       qval=data$q)
      }else if(data$data_type == 3){
        DRM_L4_N <- DRM.L4_ND(L4resNI[4:7], data$x, data$q)
        # obtain loglikelihood
        llL4N=llfL4_ND(optL4_NI$par[c(1,2,9,4,5)],
                       nvec=data$n,
                       dvec=data$x,
                       mvec=data$m,
                       s2vec=data$s2,
                       qval=data$q)
      }


      llN = c(llN, llL4N)
      BMDL = c(BMDL, L4resNI[1])
      BMD = c(BMD, L4resNI[2])
      BMDU = c(BMDU, L4resNI[3])
    }
  }
  if(prior.weights[8] == 0){L4resNI=NULL; llN = c(llN,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA);
  L4covNI=rep(NA,2); L4corrNI=rep(NA,2); DRM_L4_N=rep(NA,length(data$x)); parsL4N <- NA
  L4outNI <- t(data.frame(
    par.dt = rep(NA,3),
    par.is2t = rep(NA,3),
    par.k = rep(NA,3),
    par.a = rep(NA,3),
    par.b = rep(NA,3),
    par.c = rep(NA,3),
    par.d = rep(NA,3),
    par.s2 = rep(NA,3),
    BMD = rep(NA,3),
    min.resp = rep(NA,3),
    max.resp = rep(NA,3),
    fold.change = rep(NA,3)
  ))
  }


  ## Data to use for Lognormal distribution

  data = data.LN$data
  start = data.LN$start
  startQ = data.LN$startQ

  llLN = c()

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[9]))
  }
  if(prior.weights[9]>0){
    # print(9)

    optE4_LNI <- fun_optim(stanmodels$mE4, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optE4_LNI[[3]]),TRUE,(optE4_LNI[[3]]!=0)) | length(optE4_LNI)!=9)){
      prior.weights[9] <- 0
      warning('difficulties fitting the Exponential (Lognormal) model; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsE4LN <- par_extract(optE4_LNI, model_name = "E4_LN")
      E4resLNI <- quantile(parsE4LN$BMD*data$maxD, pvec, na.rm = T)
      E4resLNI <- c(E4resLNI, apply(parsE4LN[,c(paste0("p",1:4), "is2t")], 2, median))
      names(E4resLNI) <- c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t")

      E4outLNI <- outLP(parsE4LN, pvec, data$maxD)


      # Covariance between b-d and between BMD-d
      E4covLNI = c(cov(parsE4LN[,c("b","d")], use="na.or.complete")["b","d"],
                   cov(parsE4LN[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      E4corrLNI = c(cor(parsE4LN[,c("b","d")], use="na.or.complete")["b","d"],
                    cor(parsE4LN[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      if(data$data_type == 2){
        DRM_E4_LN <- exp(DRM.E4_LNI(E4resLNI[4:7], data$x, data$q, data$shift))
        # obtain loglikelihood
        llE4LN=llfE4_LNI(optE4_LNI$par[c(1,2,9,4,5)],
                         nvec=data$n,
                         dvec=data$x,
                         mvec=data$m,
                         s2vec=data$s2,
                         qval=data$q,
                         shift=data$shift)
      }else if(data$data_type == 4){
        DRM_E4_LN <- exp(DRM.E4_LND(E4resLNI[4:7], data$x, data$q, data$shift))
        # obtain loglikelihood
        llE4LN=llfE4_LND(optE4_LNI$par[c(1,2,9,4,5)],
                         nvec=data$n,
                         dvec=data$x,
                         mvec=data$m,
                         s2vec=data$s2,
                         qval=data$q,
                         shift=data$shift)
      }
      # llE4fLN=function(x) llfE4_LNI(x,nvec=data$n,dvec=data$x,mvec=data$m,s2vec=data$s2,qval=data$q)+
      #   mvtnorm::dmvnorm(x[c(4,5)],mean=data$priormu[c(4,5)],
      #                    sigma=data$priorSigma[c(4,5),c(4,5)], log = T)+
      #   mc2d::dpert(x[1], min = data$priorlb[1], max = data$priorub[1],
      #               mode = data$priormu[1], shape = data$shape.a, log=T)+
      #   mc2d::dpert(x[2], min = data$priorlb[2], max = data$priorub[2],
      #               mode = data$priormu[2], shape = data$shape.BMD,log=T)+
      #   mc2d::dpert(x[3], min = data$priorlb[3], max = data$priorub[3],
      #               mode = data$priormu[3], shape = data$shape.c,log=T)


      llLN = c(llLN, llE4LN)
      BMDL = c(BMDL, E4resLNI[1])
      BMD = c(BMD, E4resLNI[2])
      BMDU = c(BMDU, E4resLNI[3])
    }
  }
  if(prior.weights[9] == 0){E4resLNI=NULL; llLN = c(llLN,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA);
  E4covLNI=rep(NA,2); E4corrLNI=rep(NA,2); DRM_E4_LN=rep(NA,length(data$x))
  parsE4LN <- NA
  E4outLNI <- t(data.frame(
    par.dt = rep(NA,3),
    par.is2t = rep(NA,3),
    par.k = rep(NA,3),
    par.a = rep(NA,3),
    par.b = rep(NA,3),
    par.c = rep(NA,3),
    par.d = rep(NA,3),
    par.s2 = rep(NA,3),
    BMD = rep(NA,3),
    min.resp = rep(NA,3),
    max.resp = rep(NA,3),
    fold.change = rep(NA,3)
  ))
  }
  #
  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[10]))
  }
  if(prior.weights[10]>0){
    # print(10)

    optIE4_LNI <- fun_optim(stanmodels$mIE4, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optIE4_LNI[[3]]),TRUE,(optIE4_LNI[[3]]!=0)) | length(optIE4_LNI)!=9)){
      prior.weights[10] <- 0
      warning('difficulties fitting the Inverse Exponential (Lognormal) model; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsIE4LN <- par_extract(optIE4_LNI, model_name = "IE4_LN")
      IE4resLNI <- quantile(parsIE4LN$BMD*data$maxD, pvec, na.rm = T)
      IE4resLNI <- c(IE4resLNI, apply(parsIE4LN[,c(paste0("p",1:4), "is2t")], 2, median))
      names(IE4resLNI) <- c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t")

      IE4outLNI <- outLP(parsIE4LN, pvec, data$maxD)

      # Covariance between b-d and between BMD-d
      IE4covLNI = c(cov(parsIE4LN[,c("b","d")], use="na.or.complete")["b","d"],
                    cov(parsIE4LN[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      IE4corrLNI = c(cor(parsIE4LN[,c("b","d")], use="na.or.complete")["b","d"],
                     cor(parsIE4LN[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      if(data$data_type == 2){
        DRM_IE4_LN <- exp(DRM.IE4_LNI(IE4resLNI[4:7], data$x, data$q, data$shift))
        # obtain loglikelihood
        llIE4LN=llfIE4_LNI(optIE4_LNI$par[c(1,2,9,4,5)],
                           nvec=data$n,
                           dvec=data$x,
                           mvec=data$m,
                           s2vec=data$s2,
                           qval=data$q,
                           shift=data$shift)
      }else if(data$data_type == 4){
        DRM_IE4_LN <- exp(DRM.IE4_LND(IE4resLNI[4:7], data$x, data$q, data$shift))
        # obtain loglikelihood
        llIE4LN=llfIE4_LND(optIE4_LNI$par[c(1,2,9,4,5)],
                           nvec=data$n,
                           dvec=data$x,
                           mvec=data$m,
                           s2vec=data$s2,
                           qval=data$q,
                           shift=data$shift)
      }

      llLN = c(llLN, llIE4LN)
      BMDL = c(BMDL, IE4resLNI[1])
      BMD = c(BMD, IE4resLNI[2])
      BMDU = c(BMDU, IE4resLNI[3])
    }
  }
  if(prior.weights[10] == 0){IE4resLNI=NULL; llLN = c(llLN,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA);
  IE4covLNI=rep(NA,2); IE4corrLNI=rep(NA,2); DRM_IE4_LN=rep(NA,length(data$x))
  parsIE4LN <- NA
  IE4outLNI <- t(data.frame(
    par.dt = rep(NA,3),
    par.is2t = rep(NA,3),
    par.k = rep(NA,3),
    par.a = rep(NA,3),
    par.b = rep(NA,3),
    par.c = rep(NA,3),
    par.d = rep(NA,3),
    par.s2 = rep(NA,3),
    BMD = rep(NA,3),
    min.resp = rep(NA,3),
    max.resp = rep(NA,3),
    fold.change = rep(NA,3)
  ))
  }
  #
  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[11]))
  }
  if(prior.weights[11]>0){
    # print(11)

    optH4_LNI <- fun_optim(stanmodels$mH4, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optH4_LNI[[3]]),TRUE,(optH4_LNI[[3]]!=0)) | length(optH4_LNI)!=9)){
      prior.weights[11] <- 0
      warning('difficulties fitting the Hill (Lognormal) model; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsH4LN <- par_extract(optH4_LNI, model_name = "H4_LN")
      H4resLNI <- quantile(parsH4LN$BMD*data$maxD, pvec, na.rm = T)
      H4resLNI <- c(H4resLNI, apply(parsH4LN[,c(paste0("p",1:4), "is2t")], 2, median))
      names(H4resLNI) <- c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t")

      H4outLNI <- outLP(parsH4LN, pvec, data$maxD)

      # Covariance between b-d and between BMD-d
      H4covLNI = c(cov(parsH4LN[,c("b","d")], use="na.or.complete")["b","d"],
                   cov(parsH4LN[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      H4corrLNI = c(cor(parsH4LN[,c("b","d")], use="na.or.complete")["b","d"],
                    cor(parsH4LN[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      if(data$data_type == 2){
        DRM_H4_LN <- exp(DRM.H4_LNI(H4resLNI[4:7], data$x, data$q, data$shift))
        # obtain loglikelihood
        llH4LN=llfH4_LNI(optH4_LNI$par[c(1,2,9,4,5)],
                         nvec=data$n,
                         dvec=data$x,
                         mvec=data$m,
                         s2vec=data$s2,
                         qval=data$q,
                         shift=data$shift)
      }else if(data$data_type == 4){
        DRM_H4_LN <- exp(DRM.H4_LND(H4resLNI[4:7], data$x, data$q, data$shift))
        # obtain loglikelihood
        llH4LN=llfH4_LND(optH4_LNI$par[c(1,2,9,4,5)],
                         nvec=data$n,
                         dvec=data$x,
                         mvec=data$m,
                         s2vec=data$s2,
                         qval=data$q,
                         shift=data$shift)
      }

      llLN = c(llLN, llH4LN)
      BMDL = c(BMDL, H4resLNI[1])
      BMD = c(BMD, H4resLNI[2])
      BMDU = c(BMDU, H4resLNI[3])
    }
  }
  if(prior.weights[11] == 0){H4resLNI=NULL; llLN = c(llLN,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA);
  H4covLNI=rep(NA,2); H4corrLNI=rep(NA,2); DRM_H4_LN=rep(NA,length(data$x))
  parsH4LN <- NA
  H4outLNI <- t(data.frame(
    par.dt = rep(NA,3),
    par.is2t = rep(NA,3),
    par.k = rep(NA,3),
    par.a = rep(NA,3),
    par.b = rep(NA,3),
    par.c = rep(NA,3),
    par.d = rep(NA,3),
    par.s2 = rep(NA,3),
    BMD = rep(NA,3),
    min.resp = rep(NA,3),
    max.resp = rep(NA,3),
    fold.change = rep(NA,3)
  ))
  }
  #
  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[12]))
  }
  if(prior.weights[12]>0){
    # print(12)

    optLN4_LNI <- fun_optim(stanmodels$mLN4, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optLN4_LNI[[3]]),TRUE,(optLN4_LNI[[3]]!=0)) | length(optLN4_LNI)!=9)){
      prior.weights[12] <- 0
      warning('difficulties fitting the Lognormal (Lognormal) model; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsLN4LN <- par_extract(optLN4_LNI, model_name = "LN4_LN")
      LN4resLNI <- quantile(parsLN4LN$BMD*data$maxD, pvec, na.rm = T)
      LN4resLNI <- c(LN4resLNI, apply(parsLN4LN[,c(paste0("p",1:4), "is2t")], 2, median))
      names(LN4resLNI) <- c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t")

      LN4outLNI <- outLP(parsLN4LN, pvec, data$maxD)

      # Covariance between b-d and between BMD-d
      LN4covLNI = c(cov(parsLN4LN[,c("b","d")], use="na.or.complete")["b","d"],
                    cov(parsLN4LN[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      LN4corrLNI = c(cor(parsLN4LN[,c("b","d")], use="na.or.complete")["b","d"],
                     cor(parsLN4LN[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      if(data$data_type == 2){
        DRM_LN4_LN <- exp(DRM.LN4_LNI(LN4resLNI[4:7], data$x, data$q, data$shift))
        # obtain loglikelihood
        llLN4LN=llfLN4_LNI(optLN4_LNI$par[c(1,2,9,4,5)],
                           nvec=data$n,
                           dvec=data$x,
                           mvec=data$m,
                           s2vec=data$s2,
                           qval=data$q,
                           shift=data$shift)
      }else if(data$data_type == 4){
        DRM_LN4_LN <- exp(DRM.LN4_LND(LN4resLNI[4:7], data$x, data$q, data$shift))
        # obtain loglikelihood
        llLN4LN=llfLN4_LND(optLN4_LNI$par[c(1,2,9,4,5)],
                           nvec=data$n,
                           dvec=data$x,
                           mvec=data$m,
                           s2vec=data$s2,
                           qval=data$q,
                           shift=data$shift)
      }


      llLN = c(llLN, llLN4LN)
      BMDL = c(BMDL, LN4resLNI[1])
      BMD = c(BMD, LN4resLNI[2])
      BMDU = c(BMDU, LN4resLNI[3])
    }
  }
  if(prior.weights[12] == 0){LN4resLNI=NULL; llLN = c(llLN,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA);
  LN4covLNI=rep(NA,2); LN4corrLNI=rep(NA,2); DRM_LN4_LN=rep(NA,length(data$x))
  parsLN4LN <- NA
  LN4outLNI <- t(data.frame(
    par.dt = rep(NA,3),
    par.is2t = rep(NA,3),
    par.k = rep(NA,3),
    par.a = rep(NA,3),
    par.b = rep(NA,3),
    par.c = rep(NA,3),
    par.d = rep(NA,3),
    par.s2 = rep(NA,3),
    BMD = rep(NA,3),
    min.resp = rep(NA,3),
    max.resp = rep(NA,3),
    fold.change = rep(NA,3)
  ))
  }
  #
  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[13]))
  }
  if(prior.weights[13]>0){
    # print(13)

    # if(data$is_increasing == 1){
    #   data$init_b = qgamma(log(1+data$q)/
    #                          (optE4_LNI$par[7]*(optE4_LNI$par[8]-1)), rate=1.0, shape=optE4_LNI$par[10])/optE4_LNI$par[2]
    # }else if(data$is_decreasing == 1){
    #   data$init_b = qgamma(log(1-data$q)/
    #                          (optE4_LNI$par[7]*(optE4_LNI$par[8]-1)), rate=1.0, shape=optE4_LNI$par[10])/optE4_LNI$par[2]
    # }

    optG4_LNI <- fun_optim(stanmodels$mG4, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optG4_LNI[[3]]),TRUE,(optG4_LNI[[3]]!=0)) | length(optG4_LNI)!=9)){
      prior.weights[13] <- 0
      warning('difficulties fitting the Gamma (Lognormal) model; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsG4LN <- par_extract(optG4_LNI, model_name = "G4_LN")
      G4resLNI <- quantile(parsG4LN$BMD*data$maxD, pvec, na.rm = T)
      G4resLNI <- c(G4resLNI, apply(parsG4LN[,c(paste0("p",1:4), "is2t")], 2, median))
      names(G4resLNI) <- c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t")

      G4outLNI <- outLP(parsG4LN, pvec, data$maxD)

      # Covariance between b-d and between BMD-d
      G4covLNI = c(cov(parsG4LN[,c("b","d")], use="na.or.complete")["b","d"],
                   cov(parsG4LN[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      G4corrLNI = c(cor(parsG4LN[,c("b","d")], use="na.or.complete")["b","d"],
                    cor(parsG4LN[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      if(data$data_type == 2){
        DRM_G4_LN <- DRM.G4_LNI(G4resLNI[4:7], data$x, data$q, data$shift)
        # obtain loglikelihood
        llG4LN=llfG4_LNI(optG4_LNI$par[c(1,2,9,4,5)],
                         nvec=data$n,
                         dvec=data$x,
                         mvec=data$m,
                         s2vec=data$s2,
                         qval=data$q,
                         shift=data$shift)
      }else if(data$data_type == 4){
        DRM_G4_LN <- DRM.G4_LND(G4resLNI[4:7], data$x, data$q, data$shift)
        # obtain loglikelihood
        llG4LN=llfG4_LND(optG4_LNI$par[c(1,2,9,4,5)],
                         nvec=data$n,
                         dvec=data$x,
                         mvec=data$m,
                         s2vec=data$s2,
                         qval=data$q,
                         shift=data$shift)
      }

      llLN = c(llLN, llG4LN)
      BMDL = c(BMDL, G4resLNI[1])
      BMD = c(BMD, G4resLNI[2])
      BMDU = c(BMDU, G4resLNI[3])
    }
  }
  if(prior.weights[13] == 0){G4resLNI=NULL; llLN = c(llLN,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA);
  BMDU=c(BMDU,NA); G4covLNI=rep(NA,2); G4corrLNI=rep(NA,2); DRM_G4_LN=rep(NA,length(data$x))
  parsG4LN <- NA
  G4outLNI <- t(data.frame(
    par.dt = rep(NA,3),
    par.is2t = rep(NA,3),
    par.k = rep(NA,3),
    par.a = rep(NA,3),
    par.b = rep(NA,3),
    par.c = rep(NA,3),
    par.d = rep(NA,3),
    par.s2 = rep(NA,3),
    BMD = rep(NA,3),
    min.resp = rep(NA,3),
    max.resp = rep(NA,3),
    fold.change = rep(NA,3)
  ))
  }
  #
  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[14]))
  }
  if(prior.weights[14]>0){
    # print(14)

    optQE4_LNI <- fun_optim(stanmodels$mQE4, data, startQ, ndraws, 123, pvec)

    if((ifelse(is.na(optQE4_LNI[[3]]),TRUE,(optQE4_LNI[[3]]!=0)) | length(optQE4_LNI)!=9)){
      prior.weights[14] <- 0
      warning('difficulties fitting the Quadratic Exponential (Lognormal) model; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsQE4LN <- par_extract(optQE4_LNI, model_name = "QE4_LN")
      QE4resLNI <- quantile(parsQE4LN$BMD*data$maxD, pvec, na.rm = T)
      QE4resLNI <- c(QE4resLNI, apply(parsQE4LN[,c(paste0("p",1:4), "is2t")], 2, median))
      names(QE4resLNI) <- c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t")

      QE4outLNI <- outLP(parsQE4LN, pvec, data$maxD)

      # Covariance between b-d and between BMD-d
      QE4covLNI = c(cov(parsQE4LN[,c("b","d")], use="na.or.complete")["b","d"],
                    cov(parsQE4LN[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      QE4corrLNI = c(cor(parsQE4LN[,c("b","d")], use="na.or.complete")["b","d"],
                     cor(parsQE4LN[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      if(data$data_type == 2){
        DRM_QE4_LN <- exp(DRM.QE4_LNI(QE4resLNI[4:7], data$x, data$q, data$shift))
        # obtain loglikelihood
        llQE4LN=llfQE4_LNI(optQE4_LNI$par[c(1,2,9,4,5)],
                           nvec=data$n,
                           dvec=data$x,
                           mvec=data$m,
                           s2vec=data$s2,
                           qval=data$q,
                           shift=data$shift)
      }else if(data$data_type == 4){
        DRM_QE4_LN <- exp(DRM.QE4_LND(QE4resLNI[4:7], data$x, data$q, data$shift))
        # obtain loglikelihood
        llQE4LN=llfQE4_LND(optQE4_LNI$par[c(1,2,9,4,5)],
                           nvec=data$n,
                           dvec=data$x,
                           mvec=data$m,
                           s2vec=data$s2,
                           qval=data$q,
                           shift=data$shift)
      }

      llLN = c(llLN, llQE4LN)
      BMDL = c(BMDL, QE4resLNI[1])
      BMD = c(BMD, QE4resLNI[2])
      BMDU = c(BMDU, QE4resLNI[3])
    }
  }
  if(prior.weights[14] == 0){QE4resLNI=NULL; llLN = c(llLN,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA);
  QE4covLNI=rep(NA,2); QE4corrLNI=rep(NA,2); DRM_QE4_LN=rep(NA,length(data$x))
  parsQE4LN <- NA
  QE4outLNI <- t(data.frame(
    par.dt = rep(NA,3),
    par.is2t = rep(NA,3),
    par.k = rep(NA,3),
    par.a = rep(NA,3),
    par.b = rep(NA,3),
    par.c = rep(NA,3),
    par.d = rep(NA,3),
    par.s2 = rep(NA,3),
    BMD = rep(NA,3),
    min.resp = rep(NA,3),
    max.resp = rep(NA,3),
    fold.change = rep(NA,3)
  ))
  }
  #
  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[15]))
  }
  if(prior.weights[15]>0){
    # print(15)

    optP4_LNI <- fun_optim(stanmodels$mP4, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optP4_LNI[[3]]),TRUE,(optP4_LNI[[3]]!=0)) | length(optP4_LNI)!=9)){
      prior.weights[15] <- 0
      warning('difficulties fitting the Probit (Lognormal) model; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsP4LN <- par_extract(optP4_LNI, model_name = "P4_LN")
      P4resLNI <- quantile(parsP4LN$BMD*data$maxD, pvec, na.rm = T)
      P4resLNI <- c(P4resLNI, apply(parsP4LN[,c(paste0("p",1:4), "is2t")], 2, median))
      names(P4resLNI) <- c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t")

      P4outLNI <- outLP(parsP4LN, pvec, data$maxD)


      # Covariance between b-d and between BMD-d
      P4covLNI = c(cov(parsP4LN[,c("b","d")], use="na.or.complete")["b","d"],
                   cov(parsP4LN[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      P4corrLNI = c(cor(parsP4LN[,c("b","d")], use="na.or.complete")["b","d"],
                    cor(parsP4LN[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      if(data$data_type == 2){
        DRM_P4_LN <- exp(DRM.P4_LNI(P4resLNI[4:7], data$x, data$q, data$shift))
        # obtain loglikelihood
        llP4LN=llfP4_LNI(optP4_LNI$par[c(1,2,9,4,5)],
                         nvec=data$n,
                         dvec=data$x,
                         mvec=data$m,
                         s2vec=data$s2,
                         qval=data$q,
                         shift=data$shift)
      }else if(data$data_type == 4){
        DRM_P4_LN <- exp(DRM.P4_LND(P4resLNI[4:7], data$x, data$q, data$shift))
        # obtain loglikelihood
        llP4LN=llfP4_LND(optP4_LNI$par[c(1,2,9,4,5)],
                         nvec=data$n,
                         dvec=data$x,
                         mvec=data$m,
                         s2vec=data$s2,
                         qval=data$q,
                         shift=data$shift)
      }

      llLN = c(llLN, llP4LN)
      BMDL = c(BMDL, P4resLNI[1])
      BMD = c(BMD, P4resLNI[2])
      BMDU = c(BMDU, P4resLNI[3])
    }
  }
  if(prior.weights[15] == 0){P4resLNI=NULL; llLN = c(llLN,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA);
  P4covLNI=rep(NA,2); P4corrLNI=rep(NA,2); DRM_P4_LN=rep(NA,length(data$x))
  parsP4LN <- NA
  P4outLNI <- t(data.frame(
    par.dt = rep(NA,3),
    par.is2t = rep(NA,3),
    par.k = rep(NA,3),
    par.a = rep(NA,3),
    par.b = rep(NA,3),
    par.c = rep(NA,3),
    par.d = rep(NA,3),
    par.s2 = rep(NA,3),
    BMD = rep(NA,3),
    min.resp = rep(NA,3),
    max.resp = rep(NA,3),
    fold.change = rep(NA,3)
  ))
  }
  #
  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[16]))
  }
  if(prior.weights[16]>0){
    # print(16)

    optL4_LNI <- fun_optim(stanmodels$mL4, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optL4_LNI[[3]]),TRUE,(optL4_LNI[[3]]!=0)) | length(optL4_LNI)!=9)){
      prior.weights[16] <- 0
      warning('difficulties fitting the Logit (Lognormal) model; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsL4LN <- par_extract(optL4_LNI, model_name = "L4_LN")
      L4resLNI <- quantile(parsL4LN$BMD*data$maxD, pvec, na.rm = T)
      L4resLNI <- c(L4resLNI, apply(parsL4LN[,c(paste0("p",1:4), "is2t")], 2, median))
      names(L4resLNI) <- c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t")

      L4outLNI <- outLP(parsL4LN, pvec, data$maxD)

      # Covariance between b-d and between BMD-d
      L4covLNI = c(cov(parsL4LN[,c("b","d")], use="na.or.complete")["b","d"],
                   cov(parsL4LN[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      L4corrLNI = c(cor(parsL4LN[,c("b","d")], use="na.or.complete")["b","d"],
                    cor(parsL4LN[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      if(data$data_type == 2){
        DRM_L4_LN <- exp(DRM.L4_LNI(L4resLNI[4:7], data$x, data$q, data$shift))
        # obtain loglikelihood
        llL4LN=llfL4_LNI(optL4_LNI$par[c(1,2,9,4,5)],
                         nvec=data$n,
                         dvec=data$x,
                         mvec=data$m,
                         s2vec=data$s2,
                         qval=data$q,
                         shift=data$shift)
      }else if(data$data_type == 4){
        DRM_L4_LN <- exp(DRM.L4_LND(L4resLNI[4:7], data$x, data$q, data$shift))
        # obtain loglikelihood
        llL4LN=llfL4_LND(optL4_LNI$par[c(1,2,9,4,5)],
                         nvec=data$n,
                         dvec=data$x,
                         mvec=data$m,
                         s2vec=data$s2,
                         qval=data$q,
                         shift=data$shift)
      }

      llLN = c(llLN, llL4LN)
      BMDL = c(BMDL, L4resLNI[1])
      BMD = c(BMD, L4resLNI[2])
      BMDU = c(BMDU, L4resLNI[3])
    }
  }
  if(prior.weights[16] == 0){L4resLNI=NULL; llLN = c(llLN,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA);
  L4covLNI=rep(NA,2); L4corrLNI=rep(NA,2); DRM_L4_LN=rep(NA,length(data$x))
  parsL4LN <- NA
  L4outLNI <- t(data.frame(
    par.dt = rep(NA,3),
    par.is2t = rep(NA,3),
    par.k = rep(NA,3),
    par.a = rep(NA,3),
    par.b = rep(NA,3),
    par.c = rep(NA,3),
    par.d = rep(NA,3),
    par.s2 = rep(NA,3),
    BMD = rep(NA,3),
    min.resp = rep(NA,3),
    max.resp = rep(NA,3),
    fold.change = rep(NA,3)
  ))
  }

  minll = min(llN,llLN,na.rm=T)


  # the weights

  fun.w <- function(DIH, ll, min.ll, opt, mu, sig, lb, ub, s1, s2, s3, td){

    w1 <- (2*pi)^(2.5)*sqrt(DIH)*exp(ll-min.ll)*
      # d, sigma2
      # mvtnorm::dmvnorm(opt$par[c(4,5)],mean=mu[c(4,5)],
      #                  sigma=sig[c(4,5),c(4,5)])*
      # sigma2
      dnorm(opt$par[5], mu[5], sig[5,5])*
      # d
      truncnorm::dtruncnorm(opt$par[4], b = td, mean = mu[4], sd = sig[4,4])*

      # a
      mc2d::dpert(opt$par[1], min = lb[1], max = ub[1],
                  mode = mu[1], shape = s1)*
      # BMD
      mc2d::dpert(opt$par[2], min = lb[2], max = ub[2],
                  mode = mu[2], shape = s2)*
      # c
      mc2d::dpert(opt$par[9], min = lb[3], max = ub[3],
                  mode = mu[3], shape = s3)
    return(w1)
  }

  w=c()

  # normal

  data = data.N$data
  start=data.N$start
  startQ=data.N$startQ

  lls <- c(llN, llLN)
  w.msg <- ''

  max.ll = max(lls, na.rm = T)
  if(!(1 %in% prior.weights)){
    w.msg <- 'Laplace approximation could not be performed, problem fitting all models'
    warning('Laplace approximation could not be performed, problem fitting all models')
    out.stop <- 'not.fit'
  }else if(is.na(lls[which((max.ll-lls < 709) & prior.weights>0)][1])){
    lpw <- rep(0, 16)
    lpw[which(lls == max.ll)] <- 1
    w.msg <- 'Laplace weights could not be computed and one model gets all the weight; using another prior for parameter d might help'
    warning('Laplace weights could not be computed and one model gets all the weight; using another prior for parameter d might help')
  }else{

    if(FALSE %in% ((max.ll-lls[!is.na(lls)]) < 709)){
      w.msg <- 'not all models were used in computation of Laplace weights, some models set to 0; using another prior for parameter d might help'
      warning('not all models were used in computation of Laplace weights, some models set to 0; using another prior for parameter d might help')
    }

    minll <- min(lls[which((max.ll-lls < 709) & prior.weights>0)], na.rm = T)

    if(prior.weights[1]>0){
      DIHE4h=try(det(-pracma::pinv(optE4_NI$hessian)), silent = T)
      DIHE4=ifelse(DIHE4h<0 | class(DIHE4h)[1] == 'try-error',0,DIHE4h)
      w = c(w, fun.w(DIHE4, llE4N, minll, optE4_NI, data$priormu, data$priorSigma,
                     data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c,
                     data$truncd))
      # DIHE4h=try(det(-pracma::pinv(hessian(func = llE4fN,x=optE4_NI$par[c(1,2,9,3,4)],method.args=list(eps=1e-4, d=0.1, zero.tol=sqrt(.Machine$double.eps/7e-7), r=8, v=2, show.details=FALSE))))
      # DIHE4=ifelse(DIHE4h<0,0,DIHE4h)
      ## Aproximation of marginal (i.e. integrated) likelihood (= 'model evidence')

    }else{w=c(w,0)}

    if(prior.weights[2]>0){
      DIHIE4h=try(det(-pracma::pinv(optIE4_NI$hessian)), silent = T)
      DIHIE4=ifelse(DIHIE4h<0 | class(DIHIE4h)[1] == 'try-error',0,DIHIE4h)
      w = c(w, fun.w(DIHIE4, llIE4N, minll, optIE4_NI, data$priormu, data$priorSigma,
                     data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c,
                     data$truncd))
    }else{w=c(w,0)}

    if(prior.weights[3]>0){
      DIHH4h=try(det(-pracma::pinv(optH4_NI$hessian)), silent = T)
      DIHH4=ifelse(DIHH4h<0 | class(DIHH4h)[1] == 'try-error',0,DIHH4h)
      w = c(w, fun.w(DIHH4, llH4N, minll, optH4_NI, data$priormu, data$priorSigma,
                     data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c,
                     data$truncd))
    }else{w=c(w,0)}

    if(prior.weights[4]>0){
      DIHLN4h=try(det(-pracma::pinv(optLN4_NI$hessian)), silent = T)
      DIHLN4=ifelse(DIHLN4h<0 | class(DIHLN4h)[1] == 'try-error',0,DIHLN4h)
      w = c(w, fun.w(DIHLN4, llLN4N, minll, optLN4_NI, data$priormu, data$priorSigma,
                     data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c,
                     data$truncd))
    }else{w=c(w,0)}

    if(prior.weights[5]>0){
      DIHG4h=try(det(-pracma::pinv(optG4_NI$hessian)), silent = T)
      DIHG4=ifelse(DIHG4h<0 | class(DIHG4h)[1] == 'try-error',0,DIHG4h)
      w = c(w, fun.w(DIHG4, llG4N, minll, optG4_NI, data$priormu, data$priorSigma,
                     data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c,
                     data$truncd))
    }else{w=c(w,0)}

    if(prior.weights[6]>0){
      DIHQE4h=try(det(-pracma::pinv(optQE4_NI$hessian)), silent = T)
      DIHQE4=ifelse(DIHQE4h<0 | class(DIHQE4h)[1] == 'try-error',0,DIHQE4h)
      # w = c(w, fun.w(DIHQE4, llQE4N, minll, optQE4_NI, data$priormu, data$priorSigma,
      #                data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c))
      w = c(w, fun.w(DIHQE4, llQE4N, minll, optQE4_NI, data$priormuQ, data$priorSigmaQ,
                     data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c,
                     data$truncdQ))
    }else{w=c(w,0)}

    if(prior.weights[7]>0){
      DIHP4h=try(det(-pracma::pinv(optP4_NI$hessian)), silent = T)
      DIHP4=ifelse(DIHP4h<0 | class(DIHP4h)[1] == 'try-error',0,DIHP4h)
      w = c(w, fun.w(DIHP4, llP4N, minll, optP4_NI, data$priormu, data$priorSigma,
                     data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c,
                     data$truncd))
    }else{w=c(w,0)}

    if(prior.weights[8]>0){
      DIHL4h=try(det(-pracma::pinv(optL4_NI$hessian)), silent = T)
      DIHL4=ifelse(DIHL4h<0 | class(DIHL4h)[1] == 'try-error',0,DIHL4h)
      w = c(w, fun.w(DIHL4, llL4N, minll, optL4_NI, data$priormu, data$priorSigma,
                     data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c,
                     data$truncd))
    }else{w=c(w,0)}


    # lognormal

    data=data.LN$data
    start=data.LN$start
    startQ=data.LN$startQ

    if(prior.weights[9]>0){
      DIHE4h=try(det(-pracma::pinv(optE4_LNI$hessian)), silent = T)
      DIHE4=ifelse(DIHE4h<0 | class(DIHE4h)[1] == 'try-error',0,DIHE4h)
      w = c(w, fun.w(DIHE4, llE4LN, minll, optE4_LNI, data$priormu, data$priorSigma,
                     data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c,
                     data$truncd))
      # DIHE4h=try(det(-pracma::pinv(hessian(func = llE4fLN,x=optE4_LNI$par[c(1,2,9,3,4)],method.args=list(eps=1e-4, d=0.1, zero.tol=sqrt(.Machine$double.eps/7e-7), r=8, v=2, show.details=FALSE))))
      # DIHE4=ifelse(DIHE4h<0,0,DIHE4h)

    }else{w=c(w,0)}

    if(prior.weights[10]>0){
      DIHIE4h=try(det(-pracma::pinv(optIE4_LNI$hessian)), silent = T)
      DIHIE4=ifelse(DIHIE4h<0 | class(DIHIE4h)[1] == 'try-error',0,DIHIE4h)
      w = c(w, fun.w(DIHIE4, llIE4LN, minll, optIE4_LNI, data$priormu, data$priorSigma,
                     data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c,
                     data$truncd))
    }else{w=c(w,0)}

    if(prior.weights[11]>0){
      DIHH4h=try(det(-pracma::pinv(optH4_LNI$hessian)), silent = T)
      DIHH4=ifelse(DIHH4h<0 | class(DIHH4h)[1] == 'try-error',0,DIHH4h)
      w = c(w, fun.w(DIHH4, llH4LN, minll, optH4_LNI, data$priormu, data$priorSigma,
                     data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c,
                     data$truncd))
    }else{w=c(w,0)}

    if(prior.weights[12]>0){
      DIHLN4h=try(det(-pracma::pinv(optLN4_LNI$hessian)), silent = T)
      DIHLN4=ifelse(DIHLN4h<0 | class(DIHLN4h)[1] == 'try-error',0,DIHLN4h)
      w = c(w, fun.w(DIHLN4, llLN4LN, minll, optLN4_LNI, data$priormu, data$priorSigma,
                     data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c,
                     data$truncd))
    }else{w=c(w,0)}

    if(prior.weights[13]>0){
      DIHG4h=try(det(-pracma::pinv(optG4_LNI$hessian)), silent = T)
      DIHG4=ifelse(DIHG4h<0 | class(DIHG4h)[1] == 'try-error',0,DIHG4h)
      w = c(w, fun.w(DIHG4, llG4LN, minll, optG4_LNI, data$priormu, data$priorSigma,
                     data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c,
                     data$truncd))
    }else{w=c(w,0)}

    if(prior.weights[14]>0){
      DIHQE4h=try(det(-pracma::pinv(optQE4_LNI$hessian)), silent = T)
      DIHQE4=ifelse(DIHQE4h<0 | class(DIHQE4h)[1] == 'try-error',0,DIHQE4h)
      # w = c(w, fun.w(DIHQE4, llQE4LN, minll, optQE4_LNI, data$priormu, data$priorSigma,
      #                data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c))
      w = c(w, fun.w(DIHQE4, llQE4LN, minll, optQE4_LNI, data$priormuQ, data$priorSigmaQ,
                     data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c,
                     data$truncdQ))
    }else{w=c(w,0)}

    if(prior.weights[15]>0){
      DIHP4h=try(det(-pracma::pinv(optP4_LNI$hessian)), silent = T)
      DIHP4=ifelse(DIHP4h<0 | class(DIHP4h)[1] == 'try-error',0,DIHP4h)
      w = c(w, fun.w(DIHP4, llP4LN, minll, optP4_LNI, data$priormu, data$priorSigma,
                     data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c,
                     data$truncd))
    }else{w=c(w,0)}

    if(prior.weights[16]>0){
      DIHL4h=try(det(-pracma::pinv(optL4_LNI$hessian)), silent = T)
      DIHL4=ifelse(DIHL4h<0 | class(DIHL4h)[1] == 'try-error',0,DIHL4h)
      w = c(w, fun.w(DIHL4, llL4LN, minll, optL4_LNI, data$priormu, data$priorSigma,
                     data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c,
                     data$truncd))
    }else{w=c(w,0)}

    w <- ifelse(w == 'Inf' | is.na(w), 0, w)
    prior.weights = prior.weights/sum(prior.weights==1)
    lpw=(prior.weights*w)/sum(prior.weights*w)

  }

  # the model average posterior as a mixture
  if(out.stop != 'not.fit'){
    count=round(lpw*ndraws)
    mabmd=(c(# normal
      if(prior.weights[1]>0) sample(optE4_NI$theta_tilde[,2],count[1],replace=T),
      if(prior.weights[2]>0) sample(optIE4_NI$theta_tilde[,2],count[2],replace=T),
      if(prior.weights[3]>0) sample(optH4_NI$theta_tilde[,2],count[3],replace=T),
      if(prior.weights[4]>0) sample(optLN4_NI$theta_tilde[,2],count[4],replace=T),
      if(prior.weights[5]>0) sample(optG4_NI$theta_tilde[,2],count[5],replace=T),
      if(prior.weights[6]>0) sample(optQE4_NI$theta_tilde[,2],count[6],replace=T),
      if(prior.weights[7]>0) sample(optP4_NI$theta_tilde[,2],count[7],replace=T),
      if(prior.weights[8]>0) sample(optL4_NI$theta_tilde[,2],count[8],replace=T),
      # lognormal
      if(prior.weights[9]>0) sample(optE4_LNI$theta_tilde[,2],count[9],replace=T),
      if(prior.weights[10]>0) sample(optIE4_LNI$theta_tilde[,2],count[10],replace=T),
      if(prior.weights[11]>0) sample(optH4_LNI$theta_tilde[,2],count[11],replace=T),
      if(prior.weights[12]>0) sample(optLN4_LNI$theta_tilde[,2],count[12],replace=T),
      if(prior.weights[13]>0) sample(optG4_LNI$theta_tilde[,2],count[13],replace=T),
      if(prior.weights[14]>0) sample(optQE4_LNI$theta_tilde[,2],count[14],replace=T),
      if(prior.weights[15]>0) sample(optP4_LNI$theta_tilde[,2],count[15],replace=T),
      if(prior.weights[16]>0) sample(optL4_LNI$theta_tilde[,2],count[16],replace=T)
    ))

    mabkg=(c(# normal
      if(prior.weights[1]>0) sample(optE4_NI$theta_tilde[,"mu_0"],count[1],replace=T),
      if(prior.weights[2]>0) sample(optIE4_NI$theta_tilde[,"mu_0"],count[2],replace=T),
      if(prior.weights[3]>0) sample(optH4_NI$theta_tilde[,"mu_0"],count[3],replace=T),
      if(prior.weights[4]>0) sample(optLN4_NI$theta_tilde[,"mu_0"],count[4],replace=T),
      if(prior.weights[5]>0) sample(optG4_NI$theta_tilde[,"mu_0"],count[5],replace=T),
      if(prior.weights[6]>0) sample(optQE4_NI$theta_tilde[,"mu_0"],count[6],replace=T),
      if(prior.weights[7]>0) sample(optP4_NI$theta_tilde[,"mu_0"],count[7],replace=T),
      if(prior.weights[8]>0) sample(optL4_NI$theta_tilde[,"mu_0"],count[8],replace=T),
      # lognormal
      if(prior.weights[9]>0) sample(optE4_LNI$theta_tilde[,"mu_0"],count[9],replace=T),
      if(prior.weights[10]>0) sample(optIE4_LNI$theta_tilde[,"mu_0"],count[10],replace=T),
      if(prior.weights[11]>0) sample(optH4_LNI$theta_tilde[,"mu_0"],count[11],replace=T),
      if(prior.weights[12]>0) sample(optLN4_LNI$theta_tilde[,"mu_0"],count[12],replace=T),
      if(prior.weights[13]>0) sample(optG4_LNI$theta_tilde[,"mu_0"],count[13],replace=T),
      if(prior.weights[14]>0) sample(optQE4_LNI$theta_tilde[,"mu_0"],count[14],replace=T),
      if(prior.weights[15]>0) sample(optP4_LNI$theta_tilde[,"mu_0"],count[15],replace=T),
      if(prior.weights[16]>0) sample(optL4_LNI$theta_tilde[,"mu_0"],count[16],replace=T)
    ))

    mamaxy=(c(# normal
      if(prior.weights[1]>0) sample(optE4_NI$theta_tilde[,"mu_inf"],count[1],replace=T),
      if(prior.weights[2]>0) sample(optIE4_NI$theta_tilde[,"mu_inf"],count[2],replace=T),
      if(prior.weights[3]>0) sample(optH4_NI$theta_tilde[,"mu_inf"],count[3],replace=T),
      if(prior.weights[4]>0) sample(optLN4_NI$theta_tilde[,"mu_inf"],count[4],replace=T),
      if(prior.weights[5]>0) sample(optG4_NI$theta_tilde[,"mu_inf"],count[5],replace=T),
      if(prior.weights[6]>0) sample(optQE4_NI$theta_tilde[,"mu_inf"],count[6],replace=T),
      if(prior.weights[7]>0) sample(optP4_NI$theta_tilde[,"mu_inf"],count[7],replace=T),
      if(prior.weights[8]>0) sample(optL4_NI$theta_tilde[,"mu_inf"],count[8],replace=T),
      # lognormal
      if(prior.weights[9]>0) sample(optE4_LNI$theta_tilde[,"mu_inf"],count[9],replace=T),
      if(prior.weights[10]>0) sample(optIE4_LNI$theta_tilde[,"mu_inf"],count[10],replace=T),
      if(prior.weights[11]>0) sample(optH4_LNI$theta_tilde[,"mu_inf"],count[11],replace=T),
      if(prior.weights[12]>0) sample(optLN4_LNI$theta_tilde[,"mu_inf"],count[12],replace=T),
      if(prior.weights[13]>0) sample(optG4_LNI$theta_tilde[,"mu_inf"],count[13],replace=T),
      if(prior.weights[14]>0) sample(optQE4_LNI$theta_tilde[,"mu_inf"],count[14],replace=T),
      if(prior.weights[15]>0) sample(optP4_LNI$theta_tilde[,"mu_inf"],count[15],replace=T),
      if(prior.weights[16]>0) sample(optL4_LNI$theta_tilde[,"mu_inf"],count[16],replace=T)
    ))


    maci=quantile(mabmd,pvec, na.rm = T)*data$maxD ## original scale
    names(maci)=c("BMDL","BMD","BMDU")
    if(maci[1] == 0){
      maci[1] = 0.000000001
    }
    if(maci[2] == 0){
      maci[2] = 0.000000001
    }


    if(TRUE %in% (mabmd > data$maxD) && data$maxD > 1){
      mabmd = ifelse(mabmd > data$maxD, data$maxD, mabmd)
      p.msg = 'The model averaged posterior distribution has been truncated at max(Dose)^2'
      warnings('The model averaged posterior distribution has been truncated at max(Dose)^2')
    }else{
      p.msg = ''
    }

    BMDq = quantile(mabmd, seq(0,1,0.005), na.rm = T)*data$maxD ## original scale

    BMDL = c(BMDL, maci[1]); BMD = c(BMD, maci[2]); BMDU = c(BMDU, maci[3])

    names(BMDL) <- c("E4_N","IE4_N","H4_N","LN4_N","G4_N","QE4_N","P4_N","L4_N","E4_LN","IE4_LN","H4_LN",
                     "LN4_LN","G4_LN","QE4_LN","P4_LN","L4_LN","MA")
    model = c("E4_N","IE4_N","H4_N","LN4_N","G4_N","QE4_N","P4_N","L4_N","E4_LN","IE4_LN","H4_LN",
              "LN4_LN","G4_LN","QE4_LN","P4_LN","L4_LN","MA")
    model = as.factor(model)
    weight = c(lpw[1], lpw[2], lpw[3], lpw[4], lpw[5], lpw[6], lpw[7], lpw[8],
               lpw[9], lpw[10], lpw[11], lpw[12], lpw[13], lpw[14], lpw[15], lpw[16], 1)


    names(lpw) = model[1:16]

    ### Model-averaged response per dose level
    dr.MA <- c()
    for(i in 1:length(data$x)){
      dr.MA[i] = weighted.mean(x = c(DRM_E4_N[i],DRM_IE4_N[i],DRM_H4_N[i],DRM_LN4_N[i],DRM_G4_N[i],DRM_QE4_N[i],
                                     DRM_P4_N[i], DRM_L4_N[i] ,
                                     DRM_E4_LN[i], DRM_IE4_LN[i], DRM_H4_LN[i], DRM_LN4_LN[i], DRM_G4_LN[i],
                                     DRM_QE4_LN[i], DRM_P4_LN[i],DRM_L4_LN[i]),
                               w = lpw,
                               na.rm = T)
    }


    ## Covariances
    covs = t(data.frame(
      E4_N = E4covNI,
      IE4_N = IE4covNI,
      H4_N = H4covNI,
      LN4_N = LN4covNI,
      G4_N = G4covNI,
      QE4_N = QE4covNI,
      P4_N = P4covNI,
      L4_N = L4covNI,
      E4_LN = E4covLNI,
      IE4_LN = IE4covLNI,
      H4_LN = H4covLNI,
      LN4_LN = LN4covLNI,
      G4_LN = G4covLNI,
      QE4_LN = QE4covLNI,
      P4_LN = P4covLNI,
      L4_LN = L4covLNI
    ))
    colnames(covs) = c("b-d", "BMD-d")

    corrs = t(data.frame(
      E4_N = E4corrNI,
      IE4_N = IE4corrNI,
      H4_N = H4corrNI,
      LN4_N = LN4corrNI,
      G4_N = G4corrNI,
      QE4_N = QE4corrNI,
      P4_N = P4corrNI,
      L4_N = L4corrNI,
      E4_LN = E4corrLNI,
      IE4_LN = IE4corrLNI,
      H4_LN = H4corrLNI,
      LN4_LN = LN4corrLNI,
      G4_LN = G4corrLNI,
      QE4_LN = QE4corrLNI,
      P4_LN = P4corrLNI,
      L4_LN = L4corrLNI
    ))
    colnames(corrs) = c("b-d", "BMD-d")

    modelnames = c("E4_N","IE4_N","H4_N","LN4_N","G4_N","QE4_N","P4_N","L4_N","E4_LN","IE4_LN","H4_LN","LN4_LN","G4_LN","QE4_LN","P4_LN","L4_LN")

    ### Some additional checks

    # print(maci)

    if(maci[2]/maci[1] > 20){
      warning('BMD/BMDL is larger than 20')
    }
    if(maci[3]/maci[1] > 50){
      warning('BMDU/BMDL is larger than 50')
    }
    if(maci[2] < (data.N$data$x[2]*data.N$data$maxD/10)){
      warning('BMD is 10 times lower than the lowest non-zero dose')
    }

    ### best fitting model vs saturated ANOVA model
    best.fit = modelnames[which(weight[1:16] == max(weight[1:16]))][1]
    nrchains = 3; nriterations = 3000; warmup = 1000; delta = 0.8; treedepth = 10
    bfTest <- modelTest(best.fit, data.N, data.LN, get(paste0('opt', best.fit, 'I')), type = 'Laplace',
                        seed, ndraws, nrchains, nriterations, warmup, delta, treedepth)

    warning(bfTest$warn.bf)

  }

  if(out.stop == 'not.fit'){
    ret_results <- list(MA = c(NA,NA,NA),
                        MA_post = NA,
                        bkg_post = NA,
                        maxy_post = NA)
  }else{

    ret_results <- list(
      # model parameters
      E4_N=E4outNI,IE4_N=IE4outNI,H4_N=H4outNI,LN4_N=LN4outNI,G4_N=G4outNI,QE4_N=QE4outNI,P4_N=P4outNI,L4_N=L4outNI,
      E4_LN=E4outLNI,IE4_LN=IE4outLNI,H4_LN=H4outLNI,LN4_LN=LN4outLNI,G4_LN=G4outLNI,QE4_LN=QE4outLNI,
      P4_LN=P4outLNI,L4_LN=L4outLNI,
      # covariances
      covs = covs, corrs = corrs,
      # weights and MA
      weights=lpw,
      MA=maci,
      llN=llN, llLN=llLN, MA_post = BMDq,
      bkg_post = mabkg,
      maxy_post = mamaxy,
      # dose-response
      MA_dr = dr.MA,
      parsN = list(parsE4N, parsIE4N, parsH4N, parsLN4N, parsG4N,
                   parsQE4N, parsP4N, parsL4N),
      parsLN = list(parsE4LN, parsIE4LN, parsH4LN, parsLN4LN, parsG4LN,
                    parsQE4LN, parsP4LN, parsL4LN),
      BMDMixture = (mabmd)*data$maxD,
      dataN = data.frame(
        dose = c(data.N$data$x),
        sd = sqrt(data.N$data$s2),
        m = data.N$data$m),
      dataLN = data.frame(
        dose = c(data.LN$data$x),
        sd = sqrt(data.LN$data$s2),
        m = data.LN$data$m.org),
      max.dose = data.N$data$maxD,
      q = data.N$data$q,
      # increasing = T,
      models_included_laplace = modelnames[prior.weights > 0],
      bf = bfTest$bayesFactor, gof_check = bfTest$warn.bf,
      # means.SM = bfTest$means.SM, parBestFit = bfTest$par.best,
      # BIC.bestfit = bfTest$BIC.bestfit, BIC.SM = bfTest$BIC.SM,
      shift = data.LN$data$shift,
      w.msg = w.msg, p.msg = p.msg
    )

    attr(ret_results, "class") <- c("BMADR", "LP")

  }

  return(ret_results)

}


#' @rdname full.laplace_MA
#' @export
full.laplace_MAc=function(data.N, data.LN,
                          prior.weights = rep(1,16),
                          ndraws=30000,seed=123,
                          pvec=c(0.05,0.5,0.95)){

  # prior.weights = prior.weights/sum(prior.weights==1) #this is now done below when calculating weights

  out.stop = 'ok'

  data = data.N$data
  start = data.N$start
  startQ = data.N$startQ

  llN = c() # likelihoods
  BMDL = c()
  BMD = c()
  BMDU = c()

  ### Getting the posterior modes and the hessian to obtain the model specific posterior distributions
  # optimizing function: obtain point estimates by maximizing the joint posterior from the Stan model
  ### Additionally save BMD(L/U) estimates and loglikelihood values; set to NA if model has prior weight 0

  nModels = 16

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[1]))
  }
  if(prior.weights[1]>0){

    # print(1)

    optE4_NI <- fun_optimC(stanmodels$mE4c, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optE4_NI[[3]]),TRUE,(optE4_NI[[3]]!=0)) | length(optE4_NI)!=9)){
      prior.weights[1] <- 0
      warning('difficulties fitting the Exponential (Normal) model; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsE4N <- par_extractC(optE4_NI, model_name = "E4_N")
      E4resNI <- quantile(parsE4N$BMD*data$maxD, pvec, na.rm = T)
      E4resNI <- c(E4resNI,apply(parsE4N[,c(paste0("p",1:4), "is2t", "rho")], 2, median))
      names(E4resNI) <- c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t","rho")

      E4outNI <- outLP(parsE4N, pvec, data$maxD, clustered = T)

      # Covariance between b-d and between BMD-d
      E4covNI <- c(cov(parsE4N[,c("b","d")], use="na.or.complete")["b","d"],
                   cov(parsE4N[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      E4corrNI <- c(cor(parsE4N[,c("b","d")], use="na.or.complete")["b","d"],
                    cor(parsE4N[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      if(data$data_type == 1){
        DRM_E4_N <- DRM.E4_NI(E4resNI[4:7], data$x, data$q)
        # obtain loglikelihood
        llE4N=llfE4_NIc(optE4_NI$par[c(1,2,10,4,5,6)],
                        d=data$x,
                        n=data$n,
                        nij=data$nij,
                        y=data$y,
                        qval=data$q)
      }else if(data$data_type == 3){
        DRM_E4_N <- DRM.E4_ND(E4resNI[4:7], data$x, data$q)
        # obtain loglikelihood
        llE4N=llfE4_NDc(optE4_NI$par[c(1,2,10,4,5,6)],
                        d=data$x,
                        n=data$n,
                        nij=data$nij,
                        y=data$y,
                        qval=data$q)
      }
      llN = c(llN, llE4N)
      BMDL = c(BMDL, E4resNI[1])
      BMD = c(BMD, E4resNI[2])
      BMDU = c(BMDU, E4resNI[3])
    }
  }
  if(prior.weights[1] == 0){E4resNI=NULL; llN = c(llN,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA);
  E4covNI=rep(NA,2); E4corrNI=rep(NA,2); DRM_E4_N=rep(NA, length(data$x))
  parsE4N <- NA
  E4outNI <- t(data.frame(
    par.dt = rep(NA,3),
    par.is2t = rep(NA,3),
    par.k = rep(NA,3),
    par.a = rep(NA,3),
    par.b = rep(NA,3),
    par.c = rep(NA,3),
    par.d = rep(NA,3),
    par.s2 = rep(NA,3),
    BMD = rep(NA,3),
    min.resp = rep(NA,3),
    max.resp = rep(NA,3),
    fold.change = rep(NA,3), rho = rep(NA,3)
  ))
  }

  #
  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[2]))
  }
  if(prior.weights[2]>0){

    # print(2)

    optIE4_NI <- fun_optimC(stanmodels$mIE4c, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optIE4_NI[[3]]),TRUE,(optIE4_NI[[3]]!=0)) | length(optIE4_NI)!=9)){
      prior.weights[2] <- 0
      warning('difficulties fitting the Inverse Exponential (Normal) model; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsIE4N <- par_extractC(optIE4_NI, model_name = "IE4_N")
      IE4resNI <- quantile(parsIE4N$BMD*data$maxD, pvec, na.rm = T)
      IE4resNI <- c(IE4resNI,apply(parsIE4N[,c(paste0("p",1:4), "is2t","rho")], 2, median))
      names(IE4resNI) <- c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t","rho")

      IE4outNI <- outLP(parsIE4N, pvec, data$maxD, clustered = T)

      # Covariance between b-d and between BMD-d
      IE4covNI <- c(cov(parsIE4N[,c("b","d")], use="na.or.complete")["b","d"],
                    cov(parsIE4N[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      IE4corrNI <- c(cor(parsIE4N[,c("b","d")], use="na.or.complete")["b","d"],
                     cor(parsIE4N[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      if(data$data_type == 1){
        DRM_IE4_N <- DRM.IE4_NI(IE4resNI[4:7], data$x, data$q)
        # obtain loglikelihood
        llIE4N=llfIE4_NIc(optIE4_NI$par[c(1,2,10,4,5,6)],
                          d=data$x,
                          n=data$n,
                          nij=data$nij,
                          y=data$y,
                          qval=data$q)
      }else if(data$data_type == 3){
        DRM_IE4_N <- DRM.IE4_ND(IE4resNI[4:7], data$x, data$q)
        # obtain loglikelihood
        llIE4N=llfIE4_NDc(optIE4_NI$par[c(1,2,10,4,5,6)],
                          d=data$x,
                          n=data$n,
                          nij=data$nij,
                          y=data$y,
                          qval=data$q)
      }

      llN = c(llN, llIE4N)
      BMDL = c(BMDL, IE4resNI[1])
      BMD = c(BMD, IE4resNI[2])
      BMDU = c(BMDU, IE4resNI[3])
    }}
  if(prior.weights[2] == 0){
    IE4resNI=NULL; llN = c(llN,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA);
    IE4covNI=rep(NA,2); IE4corrNI=rep(NA,2); DRM_IE4_N=rep(NA,length(data$x))
    parsIE4N <- NA
    IE4outNI <- t(data.frame(
      par.dt = rep(NA,3),
      par.is2t = rep(NA,3),
      par.k = rep(NA,3),
      par.a = rep(NA,3),
      par.b = rep(NA,3),
      par.c = rep(NA,3),
      par.d = rep(NA,3),
      par.s2 = rep(NA,3),
      BMD = rep(NA,3),
      min.resp = rep(NA,3),
      max.resp = rep(NA,3),
      fold.change = rep(NA,3), rho = rep(NA,3)
    ))
  }
  #
  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[3]))
  }
  if(prior.weights[3]>0){

    # print(3)

    optH4_NI <- fun_optimC(stanmodels$mH4c, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optH4_NI[[3]]),TRUE,(optH4_NI[[3]]!=0)) | length(optH4_NI)!=9)){
      prior.weights[3] <- 0
      warning('difficulties fitting the Hill (Normal) model; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsH4N <- par_extractC(optH4_NI, model_name = "H4_N")
      H4resNI <- quantile(parsH4N$BMD*data$maxD, pvec, na.rm = T)  #exp(quantile(optE4_NI$theta_tilde[,"par[2]"],pvec))
      H4resNI <- c(H4resNI,apply(parsH4N[,c(paste0("p",1:4), "is2t","rho")], 2, median))
      names(H4resNI) <- c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t","rho")

      H4outNI <- outLP(parsH4N, pvec, data$maxD, clustered = T)

      # Covariance between b-d and between BMD-d
      H4covNI = c(cov(parsH4N[,c("b","d")], use="na.or.complete")["b","d"],
                  cov(parsH4N[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      H4corrNI = c(cor(parsH4N[,c("b","d")], use="na.or.complete")["b","d"],
                   cor(parsH4N[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      if(data$data_type == 1){
        DRM_H4_N <- DRM.H4_NI(H4resNI[4:7], data$x, data$q)
        # obtain loglikelihood
        llH4N=llfH4_NIc(optH4_NI$par[c(1,2,10,4,5,6)],
                        d=data$x,
                        n=data$n,
                        nij=data$nij,
                        y=data$y,
                        qval=data$q)
      }else if(data$data_type == 3){
        DRM_H4_N <- DRM.H4_ND(H4resNI[4:7], data$x, data$q)
        # obtain loglikelihood
        llH4N=llfH4_NDc(optH4_NI$par[c(1,2,10,4,5,6)],
                        d=data$x,
                        n=data$n,
                        nij=data$nij,
                        y=data$y,
                        qval=data$q)
      }

      llN = c(llN, llH4N)
      BMDL = c(BMDL, H4resNI[1])
      BMD = c(BMD, H4resNI[2])
      BMDU = c(BMDU, H4resNI[3])
    }}
  if(prior.weights[3] == 0){H4resNI=NULL; llN = c(llN,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA);
  BMDU=c(BMDU,NA); H4covNI=rep(NA,2); H4corrNI=rep(NA,2); DRM_H4_N=rep(NA,length(data$x))
  parsH4N <- NA
  H4outNI <- t(data.frame(
    par.dt = rep(NA,3),
    par.is2t = rep(NA,3),
    par.k = rep(NA,3),
    par.a = rep(NA,3),
    par.b = rep(NA,3),
    par.c = rep(NA,3),
    par.d = rep(NA,3),
    par.s2 = rep(NA,3),
    BMD = rep(NA,3),
    min.resp = rep(NA,3),
    max.resp = rep(NA,3),
    fold.change = rep(NA,3), rho = rep(NA,3)
  ))
  }

  #
  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[4]))
  }
  if(prior.weights[4]>0){

    # print(4)

    optLN4_NI <- fun_optimC(stanmodels$mLN4c, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optLN4_NI[[3]]),TRUE,(optLN4_NI[[3]]!=0)) | length(optLN4_NI)!=9)){
      prior.weights[4] <- 0
      warning('difficulties fitting the Lognormal (Normal) model; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsLN4N <- par_extractC(optLN4_NI, model_name = "LN4_N")
      LN4resNI <- quantile(parsLN4N$BMD*data$maxD, pvec, na.rm = T)
      LN4resNI <- c(LN4resNI,apply(parsLN4N[,c(paste0("p",1:4), "is2t","rho")], 2, median))
      names(LN4resNI) <- c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t","rho")

      LN4outNI <- outLP(parsLN4N, pvec, data$maxD, clustered = T)


      # Covariance between b-d and between BMD-d
      LN4covNI = c(cov(parsLN4N[,c("b","d")], use="na.or.complete")["b","d"],
                   cov(parsLN4N[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      LN4corrNI = c(cor(parsLN4N[,c("b","d")], use="na.or.complete")["b","d"],
                    cor(parsLN4N[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      if(data$data_type == 1){
        DRM_LN4_N <- DRM.LN4_NI(LN4resNI[4:7], data$x, data$q)
        # obtain loglikelihood
        llLN4N=llfLN4_NIc(optLN4_NI$par[c(1,2,10,4,5,6)],
                          d=data$x,
                          n=data$n,
                          nij=data$nij,
                          y=data$y,
                          qval=data$q)
      }else if(data$data_type == 3){
        DRM_LN4_N <- DRM.LN4_ND(LN4resNI[4:7], data$x, data$q)
        # obtain loglikelihood
        llLN4N=llfLN4_NDc(optLN4_NI$par[c(1,2,10,4,5,6)],
                          d=data$x,
                          n=data$n,
                          nij=data$nij,
                          y=data$y,
                          qval=data$q)
      }

      llN = c(llN, llLN4N)
      BMDL = c(BMDL, LN4resNI[1])
      BMD = c(BMD, LN4resNI[2])
      BMDU = c(BMDU, LN4resNI[3])
    }}
  if(prior.weights[4] == 0){LN4resNI=NULL; llN = c(llN,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA);
  LN4covNI=rep(NA,2); LN4corrNI=rep(NA,2); DRM_LN4_N=rep(NA,length(data$x))
  parsLN4N <- NA
  LN4outNI <- t(data.frame(
    par.dt = rep(NA,3),
    par.is2t = rep(NA,3),
    par.k = rep(NA,3),
    par.a = rep(NA,3),
    par.b = rep(NA,3),
    par.c = rep(NA,3),
    par.d = rep(NA,3),
    par.s2 = rep(NA,3),
    BMD = rep(NA,3),
    min.resp = rep(NA,3),
    max.resp = rep(NA,3),
    fold.change = rep(NA,3), rho = rep(NA,3)
  ))
  }
  #
  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[5]))
  }
  if(prior.weights[5]>0){
    # print(5)

    optG4_NI <- fun_optimC(stanmodels$mG4c, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optG4_NI[[3]]),TRUE,(optG4_NI[[3]]!=0)) | length(optG4_NI)!=9)){
      prior.weights[5] <- 0
      warning('difficulties fitting the Gamma (Normal) model; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsG4N <- par_extractC(optG4_NI, model_name = "G4_N")
      G4resNI <- quantile(parsG4N$BMD*data$maxD, pvec, na.rm = T)
      G4resNI <- c(G4resNI,apply(parsG4N[,c(paste0("p",1:4), "is2t", "rho")], 2, median))
      names(G4resNI) <- c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t","rho")

      G4outNI <- outLP(parsG4N, pvec, data$maxD, clustered = T)

      # Covariance between b-d and between BMD-d
      G4covNI = c(cov(parsG4N[,c("b","d")], use="na.or.complete")["b","d"],
                  cov(parsG4N[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      G4corrNI = c(cor(parsG4N[,c("b","d")], use="na.or.complete")["b","d"],
                   cor(parsG4N[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      if(data$data_type == 1){
        DRM_G4_N <- DRM.G4_NI(G4resNI[4:7], data$x, data$q)
        # obtain loglikelihood
        llG4N=llfG4_NIc(optG4_NI$par[c(1,2,10,4,5,6)],
                        d=data$x,
                        n=data$n,
                        nij=data$nij,
                        y=data$y,
                        qval=data$q)
      }else if(data$data_type == 3){
        DRM_G4_N <- DRM.G4_ND(G4resNI[4:7], data$x, data$q)
        # obtain loglikelihood
        llG4N=llfG4_NDc(optG4_NI$par[c(1,2,10,4,5,6)],
                        d=data$x,
                        n=data$n,
                        nij=data$nij,
                        y=data$y,
                        qval=data$q)
      }

      llN = c(llN, llG4N)
      BMDL = c(BMDL, G4resNI[1])
      BMD = c(BMD, G4resNI[2])
      BMDU = c(BMDU, G4resNI[3])
    }}
  if(prior.weights[5] == 0){G4resNI=NULL; llN = c(llN,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA);
  BMDU=c(BMDU,NA); G4covNI=rep(NA,2); G4corrNI=rep(NA,2); DRM_G4_N=rep(NA, length(data$x))
  parsG4N <-NA
  G4outNI <- t(data.frame(
    par.dt = rep(NA,3),
    par.is2t = rep(NA,3),
    par.k = rep(NA,3),
    par.a = rep(NA,3),
    par.b = rep(NA,3),
    par.c = rep(NA,3),
    par.d = rep(NA,3),
    par.s2 = rep(NA,3),
    BMD = rep(NA,3),
    min.resp = rep(NA,3),
    max.resp = rep(NA,3),
    fold.change = rep(NA,3), rho = rep(NA,3)
  ))
  }
  #

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[6]))
  }
  if(prior.weights[6]>0){
    # print(6)

    optQE4_NI <- fun_optimC(stanmodels$mQE4c, data, startQ, ndraws, 123, pvec)

    if((ifelse(is.na(optQE4_NI[[3]]),TRUE,(optQE4_NI[[3]]!=0)) | length(optQE4_NI)!=9)){
      prior.weights[6] <- 0
      warning('difficulties fitting the Quadratic Exponential (Normal) model; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsQE4N <- par_extractC(optQE4_NI, model_name = "QE4_N")
      QE4resNI <- quantile(parsQE4N$BMD*data$maxD, pvec, na.rm = T)
      QE4resNI <- c(QE4resNI,apply(parsQE4N[,c(paste0("p",1:4), "is2t", "rho")], 2, median))
      names(QE4resNI) <- c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t","rho")

      QE4outNI <- outLP(parsQE4N, pvec, data$maxD, clustered = T)

      # Covariance between b-d and between BMD-d
      QE4covNI = c(cov(parsQE4N[,c("b","d")], use="na.or.complete")["b","d"],
                   cov(parsQE4N[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      QE4corrNI = c(cor(parsQE4N[,c("b","d")], use="na.or.complete")["b","d"],
                    cor(parsQE4N[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      if(data$data_type == 1){
        DRM_QE4_N <- DRM.QE4_NI(QE4resNI[4:7], data$x, data$q)
        # obtain loglikelihood
        llQE4N=llfQE4_NIc(optQE4_NI$par[c(1,2,10,4,5,6)],
                          d=data$x,
                          n=data$n,
                          nij=data$nij,
                          y=data$y,
                          qval=data$q)
      }else if(data$data_type == 3){
        DRM_QE4_N <- DRM.QE4_ND(QE4resNI[4:7], data$x, data$q)
        # obtain loglikelihood
        llQE4N=llfQE4_NDc(optQE4_NI$par[c(1,2,10,4,5,6)],
                          d=data$x,
                          n=data$n,
                          nij=data$nij,
                          y=data$y,
                          qval=data$q)
      }

      llN = c(llN, llQE4N)
      BMDL = c(BMDL, QE4resNI[1])
      BMD = c(BMD, QE4resNI[2])
      BMDU = c(BMDU, QE4resNI[3])
    }}
  if(prior.weights[6] == 0){QE4resNI=NULL; llN = c(llN,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA);
  BMDU=c(BMDU,NA); QE4covNI=rep(NA,2); QE4corrNI=rep(NA,2); DRM_QE4_N=rep(NA,length(data$x))
  parsQE4N <- NA
  QE4outNI <- t(data.frame(
    par.dt = rep(NA,3),
    par.is2t = rep(NA,3),
    par.k = rep(NA,3),
    par.a = rep(NA,3),
    par.b = rep(NA,3),
    par.c = rep(NA,3),
    par.d = rep(NA,3),
    par.s2 = rep(NA,3),
    BMD = rep(NA,3),
    min.resp = rep(NA,3),
    max.resp = rep(NA,3),
    fold.change = rep(NA,3), rho = rep(NA,3)
  ))
  }
  #

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[7]))
  }
  if(prior.weights[7]>0){
    # print(7)

    optP4_NI <- fun_optimC(stanmodels$mP4c, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optP4_NI[[3]]),TRUE,(optP4_NI[[3]]!=0)) | length(optP4_NI)!=9)){
      prior.weights[7] <- 0
      warning('difficulties fitting the Probit (Normal) model; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsP4N <- par_extractC(optP4_NI, model_name = "P4_N")
      P4resNI <- quantile(parsP4N$BMD*data$maxD, pvec, na.rm = T)
      P4resNI <- c(P4resNI,apply(parsP4N[,c(paste0("p",1:4), "is2t","rho")], 2, median))
      names(P4resNI) <- c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t","rho")

      P4outNI <- outLP(parsP4N, pvec, data$maxD, clustered = T)

      # Covariance between b-d and between BMD-d
      P4covNI = c(cov(parsP4N[,c("b","d")], use="na.or.complete")["b","d"],
                  cov(parsP4N[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      P4corrNI = c(cor(parsP4N[,c("b","d")], use="na.or.complete")["b","d"],
                   cor(parsP4N[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      if(data$data_type == 1){
        DRM_P4_N <- DRM.P4_NI(P4resNI[4:7], data$x, data$q)
        # obtain loglikelihood
        llP4N=llfP4_NIc(optP4_NI$par[c(1,2,10,4,5,6)],
                        d=data$x,
                        n=data$n,
                        nij=data$nij,
                        y=data$y,
                        qval=data$q)
      }else if(data$data_type == 3){
        DRM_P4_N <- DRM.P4_ND(P4resNI[4:7], data$x, data$q)
        # obtain loglikelihood
        llP4N=llfP4_NDc(optP4_NI$par[c(1,2,10,4,5,6)],
                        d=data$x,
                        n=data$n,
                        nij=data$nij,
                        y=data$y,
                        qval=data$q)
      }


      llN = c(llN, llP4N)
      BMDL = c(BMDL, P4resNI[1])
      BMD = c(BMD, P4resNI[2])
      BMDU = c(BMDU, P4resNI[3])
    }}
  if(prior.weights[7] == 0){P4resNI=NULL; llN = c(llN,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA);
  P4covNI=rep(NA,2); P4corrNI=rep(NA,2); DRM_P4_N=rep(NA,length(data$x)); parsP4N <- NA
  P4outNI <- t(data.frame(
    par.dt = rep(NA,3),
    par.is2t = rep(NA,3),
    par.k = rep(NA,3),
    par.a = rep(NA,3),
    par.b = rep(NA,3),
    par.c = rep(NA,3),
    par.d = rep(NA,3),
    par.s2 = rep(NA,3),
    BMD = rep(NA,3),
    min.resp = rep(NA,3),
    max.resp = rep(NA,3),
    fold.change = rep(NA,3), rho = rep(NA,3)
  ))
  }
  #
  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[8]))
  }
  if(prior.weights[8]>0){
    # print(8)

    optL4_NI <- fun_optimC(stanmodels$mL4c, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optL4_NI[[3]]),TRUE,(optL4_NI[[3]]!=0)) | length(optL4_NI)!=9)){
      prior.weights[8] <- 0
      warning('difficulties fitting the Logit (Normal) model; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsL4N <- par_extractC(optL4_NI, model_name = "L4_N")
      L4resNI <- quantile(parsL4N$BMD*data$maxD, pvec, na.rm = T)
      L4resNI <- c(L4resNI,apply(parsL4N[,c(paste0("p",1:4), "is2t","rho")], 2, median))
      names(L4resNI) <- c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t","rho")

      L4outNI <- outLP(parsL4N, pvec, data$maxD, clustered = T)


      # Covariance between b-d and between BMD-d
      L4covNI = c(cov(parsL4N[,c("b","d")], use="na.or.complete")["b","d"],
                  cov(parsL4N[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      L4corrNI = c(cor(parsL4N[,c("b","d")], use="na.or.complete")["b","d"],
                   cor(parsL4N[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      if(data$data_type == 1){
        DRM_L4_N <- DRM.L4_NI(L4resNI[4:7], data$x, data$q)
        # obtain loglikelihood
        llL4N=llfL4_NIc(optL4_NI$par[c(1,2,10,4,5,6)],
                        d=data$x,
                        n=data$n,
                        nij=data$nij,
                        y=data$y,
                        qval=data$q)
      }else if(data$data_type == 3){
        DRM_L4_N <- DRM.L4_ND(L4resNI[4:7], data$x, data$q)
        # obtain loglikelihood
        llL4N=llfL4_NDc(optL4_NI$par[c(1,2,10,4,5,6)],
                        d=data$x,
                        n=data$n,
                        nij=data$nij,
                        y=data$y,
                        qval=data$q)
      }


      llN = c(llN, llL4N)
      BMDL = c(BMDL, L4resNI[1])
      BMD = c(BMD, L4resNI[2])
      BMDU = c(BMDU, L4resNI[3])
    }}
  if(prior.weights[8] == 0){L4resNI=NULL; llN = c(llN,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA);
  L4covNI=rep(NA,2); L4corrNI=rep(NA,2); DRM_L4_N=rep(NA,length(data$x)); parsL4N <- NA
  L4outNI <- t(data.frame(
    par.dt = rep(NA,3),
    par.is2t = rep(NA,3),
    par.k = rep(NA,3),
    par.a = rep(NA,3),
    par.b = rep(NA,3),
    par.c = rep(NA,3),
    par.d = rep(NA,3),
    par.s2 = rep(NA,3),
    BMD = rep(NA,3),
    min.resp = rep(NA,3),
    max.resp = rep(NA,3),
    fold.change = rep(NA,3), rho = rep(NA,3)
  ))
  }


  ## Data to use for Lognormal distribution

  data = data.LN$data
  start = data.LN$start
  startQ = data.LN$startQ

  llLN = c()

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[9]))
  }
  if(prior.weights[9]>0){
    # print(9)

    optE4_LNI <- fun_optimC(stanmodels$mE4c, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optE4_LNI[[3]]),TRUE,(optE4_LNI[[3]]!=0)) | length(optE4_LNI)!=9)){
      prior.weights[9] <- 0
      warning('difficulties fitting the Exponential (Lognormal) model; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsE4LN <- par_extractC(optE4_LNI, model_name = "E4_LN")
      E4resLNI <- quantile(parsE4LN$BMD*data$maxD, pvec, na.rm = T)
      E4resLNI <- c(E4resLNI, apply(parsE4LN[,c(paste0("p",1:4), "is2t","rho")], 2, median))
      names(E4resLNI) <- c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t","rho")

      E4outLNI <- outLP(parsE4LN, pvec, data$maxD, clustered = T)


      # Covariance between b-d and between BMD-d
      E4covLNI = c(cov(parsE4LN[,c("b","d")], use="na.or.complete")["b","d"],
                   cov(parsE4LN[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      E4corrLNI = c(cor(parsE4LN[,c("b","d")], use="na.or.complete")["b","d"],
                    cor(parsE4LN[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      if(data$data_type == 2){
        DRM_E4_LN <- exp(DRM.E4_LNI(E4resLNI[4:7], data$x, data$q, data$shift))
        # obtain loglikelihood
        llE4LN=llfE4_LNIc(optE4_LNI$par[c(1,2,10,4,5,6)],
                          d=data$x,
                          n=data$n,
                          nij=data$nij,
                          y=data$y,
                          qval=data$q,
                          shift=data$shift)
      }else if(data$data_type == 4){
        DRM_E4_LN <- exp(DRM.E4_LND(E4resLNI[4:7], data$x, data$q, data$shift))
        # obtain loglikelihood
        llE4LN=llfE4_LNDc(optE4_LNI$par[c(1,2,10,4,5,6)],
                          d=data$x,
                          n=data$n,
                          nij=data$nij,
                          y=data$y,
                          qval=data$q,
                          shift=data$shift)
      }

      llLN = c(llLN, llE4LN)
      BMDL = c(BMDL, E4resLNI[1])
      BMD = c(BMD, E4resLNI[2])
      BMDU = c(BMDU, E4resLNI[3])
    }}
  if(prior.weights[9] == 0){E4resLNI=NULL; llLN = c(llLN,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA);
  E4covLNI=rep(NA,2); E4corrLNI=rep(NA,2); DRM_E4_LN=rep(NA,length(data$x))
  parsE4LN <- NA
  E4outLNI <- t(data.frame(
    par.dt = rep(NA,3),
    par.is2t = rep(NA,3),
    par.k = rep(NA,3),
    par.a = rep(NA,3),
    par.b = rep(NA,3),
    par.c = rep(NA,3),
    par.d = rep(NA,3),
    par.s2 = rep(NA,3),
    BMD = rep(NA,3),
    min.resp = rep(NA,3),
    max.resp = rep(NA,3),
    fold.change = rep(NA,3), rho = rep(NA,3)
  ))
  }
  #

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[10]))
  }
  if(prior.weights[10]>0){
    # print(10)

    optIE4_LNI <- fun_optimC(stanmodels$mIE4c, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optIE4_LNI[[3]]),TRUE,(optIE4_LNI[[3]]!=0)) | length(optIE4_LNI)!=9)){
      prior.weights[10] <- 0
      warning('difficulties fitting the Inverse Exponential (Lognormal) model; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsIE4LN <- par_extractC(optIE4_LNI, model_name = "IE4_LN")
      IE4resLNI <- quantile(parsIE4LN$BMD*data$maxD, pvec, na.rm = T)
      IE4resLNI <- c(IE4resLNI, apply(parsIE4LN[,c(paste0("p",1:4), "is2t","rho")], 2, median))
      names(IE4resLNI) <- c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t","rho")

      IE4outLNI <- outLP(parsIE4LN, pvec, data$maxD, clustered = T)

      # Covariance between b-d and between BMD-d
      IE4covLNI = c(cov(parsIE4LN[,c("b","d")], use="na.or.complete")["b","d"],
                    cov(parsIE4LN[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      IE4corrLNI = c(cor(parsIE4LN[,c("b","d")], use="na.or.complete")["b","d"],
                     cor(parsIE4LN[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      if(data$data_type == 2){
        DRM_IE4_LN <- exp(DRM.IE4_LNI(IE4resLNI[4:7], data$x, data$q, data$shift))
        # obtain loglikelihood
        llIE4LN=llfIE4_LNIc(optIE4_LNI$par[c(1,2,10,4,5,6)],
                            d=data$x,
                            n=data$n,
                            nij=data$nij,
                            y=data$y,
                            qval=data$q,
                            shift=data$shift)
      }else if(data$data_type == 4){
        DRM_IE4_LN <- exp(DRM.IE4_LND(IE4resLNI[4:7], data$x, data$q, data$shift))
        # obtain loglikelihood
        llIE4LN=llfIE4_LNDc(optIE4_LNI$par[c(1,2,10,4,5,6)],
                            d=data$x,
                            n=data$n,
                            nij=data$nij,
                            y=data$y,
                            qval=data$q,
                            shift=data$shift)
      }

      llLN = c(llLN, llIE4LN)
      BMDL = c(BMDL, IE4resLNI[1])
      BMD = c(BMD, IE4resLNI[2])
      BMDU = c(BMDU, IE4resLNI[3])
    }}
  if(prior.weights[10] == 0){IE4resLNI=NULL; llLN = c(llLN,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA);
  IE4covLNI=rep(NA,2); IE4corrLNI=rep(NA,2); DRM_IE4_LN=rep(NA,length(data$x))
  parsIE4LN <- NA
  IE4outLNI <- t(data.frame(
    par.dt = rep(NA,3),
    par.is2t = rep(NA,3),
    par.k = rep(NA,3),
    par.a = rep(NA,3),
    par.b = rep(NA,3),
    par.c = rep(NA,3),
    par.d = rep(NA,3),
    par.s2 = rep(NA,3),
    BMD = rep(NA,3),
    min.resp = rep(NA,3),
    max.resp = rep(NA,3),
    fold.change = rep(NA,3), rho = rep(NA,3)
  ))
  }
  #

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[11]))
  }
  if(prior.weights[11]>0){
    # print(11)

    optH4_LNI <- fun_optimC(stanmodels$mH4c, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optH4_LNI[[3]]),TRUE,(optH4_LNI[[3]]!=0)) | length(optH4_LNI)!=9)){
      prior.weights[11] <- 0
      warning('difficulties fitting the Hill (Lognormal) model; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsH4LN <- par_extractC(optH4_LNI, model_name = "H4_LN")
      H4resLNI <- quantile(parsH4LN$BMD*data$maxD, pvec, na.rm = T)
      H4resLNI <- c(H4resLNI, apply(parsH4LN[,c(paste0("p",1:4), "is2t","rho")], 2, median))
      names(H4resLNI) <- c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t","rho")

      H4outLNI <- outLP(parsH4LN, pvec, data$maxD, clustered = T)

      # Covariance between b-d and between BMD-d
      H4covLNI = c(cov(parsH4LN[,c("b","d")], use="na.or.complete")["b","d"],
                   cov(parsH4LN[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      H4corrLNI = c(cor(parsH4LN[,c("b","d")], use="na.or.complete")["b","d"],
                    cor(parsH4LN[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      if(data$data_type == 2){
        DRM_H4_LN <- exp(DRM.H4_LNI(H4resLNI[4:7], data$x, data$q, data$shift))
        # obtain loglikelihood
        llH4LN=llfH4_LNIc(optH4_LNI$par[c(1,2,10,4,5,6)],
                          d=data$x,
                          n=data$n,
                          nij=data$nij,
                          y=data$y,
                          qval=data$q,
                          shift=data$shift)
      }else if(data$data_type == 4){
        DRM_H4_LN <- exp(DRM.H4_LND(H4resLNI[4:7], data$x, data$q, data$shift))
        # obtain loglikelihood
        llH4LN=llfH4_LNDc(optH4_LNI$par[c(1,2,10,4,5,6)],
                          d=data$x,
                          n=data$n,
                          nij=data$nij,
                          y=data$y,
                          qval=data$q,
                          shift=data$shift)
      }

      llLN = c(llLN, llH4LN)
      BMDL = c(BMDL, H4resLNI[1])
      BMD = c(BMD, H4resLNI[2])
      BMDU = c(BMDU, H4resLNI[3])
    }}
  if(prior.weights[11] == 0){H4resLNI=NULL; llLN = c(llLN,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA);
  H4covLNI=rep(NA,2); H4corrLNI=rep(NA,2); DRM_H4_LN=rep(NA,length(data$x))
  parsH4LN <- NA
  H4outLNI <- t(data.frame(
    par.dt = rep(NA,3),
    par.is2t = rep(NA,3),
    par.k = rep(NA,3),
    par.a = rep(NA,3),
    par.b = rep(NA,3),
    par.c = rep(NA,3),
    par.d = rep(NA,3),
    par.s2 = rep(NA,3),
    BMD = rep(NA,3),
    min.resp = rep(NA,3),
    max.resp = rep(NA,3),
    fold.change = rep(NA,3), rho = rep(NA,3)
  ))
  }
  #

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[12]))
  }
  if(prior.weights[12]>0){
    # print(12)

    optLN4_LNI <- fun_optimC(stanmodels$mLN4c, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optLN4_LNI[[3]]),TRUE,(optLN4_LNI[[3]]!=0)) | length(optLN4_LNI)!=9)){
      prior.weights[12] <- 0
      warning('difficulties fitting the Lognormal (Lognormal) model; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsLN4LN <- par_extractC(optLN4_LNI, model_name = "LN4_LN")
      LN4resLNI <- quantile(parsLN4LN$BMD*data$maxD, pvec, na.rm = T)
      LN4resLNI <- c(LN4resLNI, apply(parsLN4LN[,c(paste0("p",1:4), "is2t","rho")], 2, median))
      names(LN4resLNI) <- c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t","rho")

      LN4outLNI <- outLP(parsLN4LN, pvec, data$maxD, clustered = T)

      # Covariance between b-d and between BMD-d
      LN4covLNI = c(cov(parsLN4LN[,c("b","d")], use="na.or.complete")["b","d"],
                    cov(parsLN4LN[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      LN4corrLNI = c(cor(parsLN4LN[,c("b","d")], use="na.or.complete")["b","d"],
                     cor(parsLN4LN[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      if(data$data_type == 2){
        DRM_LN4_LN <- exp(DRM.LN4_LNI(LN4resLNI[4:7], data$x, data$q, data$shift))
        # obtain loglikelihood
        llLN4LN=llfLN4_LNIc(optLN4_LNI$par[c(1,2,10,4,5,6)],
                            d=data$x,
                            n=data$n,
                            nij=data$nij,
                            y=data$y,
                            qval=data$q,
                            shift=data$shift)
      }else if(data$data_type == 4){
        DRM_LN4_LN <- exp(DRM.LN4_LND(LN4resLNI[4:7], data$x, data$q, data$shift))
        # obtain loglikelihood
        llLN4LN=llfLN4_LNDc(optLN4_LNI$par[c(1,2,10,4,5,6)],
                            d=data$x,
                            n=data$n,
                            nij=data$nij,
                            y=data$y,
                            qval=data$q,
                            shift=data$shift)
      }


      llLN = c(llLN, llLN4LN)
      BMDL = c(BMDL, LN4resLNI[1])
      BMD = c(BMD, LN4resLNI[2])
      BMDU = c(BMDU, LN4resLNI[3])
    }}
  if(prior.weights[12] == 0){LN4resLNI=NULL; llLN = c(llLN,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA);
  LN4covLNI=rep(NA,2); LN4corrLNI=rep(NA,2); DRM_LN4_LN=rep(NA,length(data$x))
  parsLN4LN <- NA
  LN4outLNI <- t(data.frame(
    par.dt = rep(NA,3),
    par.is2t = rep(NA,3),
    par.k = rep(NA,3),
    par.a = rep(NA,3),
    par.b = rep(NA,3),
    par.c = rep(NA,3),
    par.d = rep(NA,3),
    par.s2 = rep(NA,3),
    BMD = rep(NA,3),
    min.resp = rep(NA,3),
    max.resp = rep(NA,3),
    fold.change = rep(NA,3), rho = rep(NA,3)
  ))
  }
  #

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[13]))
  }
  if(prior.weights[13]>0){
    # print(13)

    optG4_LNI <- fun_optimC(stanmodels$mG4c, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optG4_LNI[[3]]),TRUE,(optG4_LNI[[3]]!=0)) | length(optG4_LNI)!=9)){
      prior.weights[13] <- 0
      warning('difficulties fitting the Gamma (Lognormal) model; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsG4LN <- par_extractC(optG4_LNI, model_name = "G4_LN")
      G4resLNI <- quantile(parsG4LN$BMD*data$maxD, pvec, na.rm = T)
      G4resLNI <- c(G4resLNI, apply(parsG4LN[,c(paste0("p",1:4), "is2t","rho")], 2, median))
      names(G4resLNI) <- c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t","rho")

      G4outLNI <- outLP(parsG4LN, pvec, data$maxD, clustered = T)

      # Covariance between b-d and between BMD-d
      G4covLNI = c(cov(parsG4LN[,c("b","d")], use="na.or.complete")["b","d"],
                   cov(parsG4LN[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      G4corrLNI = c(cor(parsG4LN[,c("b","d")], use="na.or.complete")["b","d"],
                    cor(parsG4LN[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      if(data$data_type == 2){
        DRM_G4_LN <- DRM.G4_LNI(G4resLNI[4:7], data$x, data$q, data$shift)
        # obtain loglikelihood
        llG4LN=llfG4_LNIc(optG4_LNI$par[c(1,2,10,4,5,6)],
                          d=data$x,
                          n=data$n,
                          nij=data$nij,
                          y=data$y,
                          qval=data$q,
                          shift=data$shift)
      }else if(data$data_type == 4){
        DRM_G4_LN <- DRM.G4_LND(G4resLNI[4:7], data$x, data$q, data$shift)
        # obtain loglikelihood
        llG4LN=llfG4_LNDc(optG4_LNI$par[c(1,2,10,4,5,6)],
                          d=data$x,
                          n=data$n,
                          nij=data$nij,
                          y=data$y,
                          qval=data$q,
                          shift=data$shift)
      }

      llLN = c(llLN, llG4LN)
      BMDL = c(BMDL, G4resLNI[1])
      BMD = c(BMD, G4resLNI[2])
      BMDU = c(BMDU, G4resLNI[3])
    }}
  if(prior.weights[13] == 0){G4resLNI=NULL; llLN = c(llLN,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA);
  BMDU=c(BMDU,NA); G4covLNI=rep(NA,2); G4corrLNI=rep(NA,2); DRM_G4_LN=rep(NA,length(data$x))
  parsG4LN <- NA
  G4outLNI <- t(data.frame(
    par.dt = rep(NA,3),
    par.is2t = rep(NA,3),
    par.k = rep(NA,3),
    par.a = rep(NA,3),
    par.b = rep(NA,3),
    par.c = rep(NA,3),
    par.d = rep(NA,3),
    par.s2 = rep(NA,3),
    BMD = rep(NA,3),
    min.resp = rep(NA,3),
    max.resp = rep(NA,3),
    fold.change = rep(NA,3), rho = rep(NA,3)
  ))
  }
  #

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[14]))
  }
  if(prior.weights[14]>0){
    # print(14)

    optQE4_LNI <- fun_optimC(stanmodels$mQE4c, data, startQ, ndraws, 123, pvec)

    if((ifelse(is.na(optQE4_LNI[[3]]),TRUE,(optQE4_LNI[[3]]!=0)) | length(optQE4_LNI)!=9)){
      prior.weights[14] <- 0
      warning('difficulties fitting the Quadratic Exponential (Lognormal) model; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsQE4LN <- par_extractC(optQE4_LNI, model_name = "QE4_LN")
      QE4resLNI <- quantile(parsQE4LN$BMD*data$maxD, pvec, na.rm = T)
      QE4resLNI <- c(QE4resLNI, apply(parsQE4LN[,c(paste0("p",1:4), "is2t","rho")], 2, median))
      names(QE4resLNI) <- c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t","rho")

      QE4outLNI <- outLP(parsQE4LN, pvec, data$maxD, clustered = T)

      # Covariance between b-d and between BMD-d
      QE4covLNI = c(cov(parsQE4LN[,c("b","d")], use="na.or.complete")["b","d"],
                    cov(parsQE4LN[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      QE4corrLNI = c(cor(parsQE4LN[,c("b","d")], use="na.or.complete")["b","d"],
                     cor(parsQE4LN[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      if(data$data_type == 2){
        DRM_QE4_LN <- exp(DRM.QE4_LNI(QE4resLNI[4:7], data$x, data$q, data$shift))
        # obtain loglikelihood
        llQE4LN=llfQE4_LNIc(optQE4_LNI$par[c(1,2,10,4,5,6)],
                            d=data$x,
                            n=data$n,
                            nij=data$nij,
                            y=data$y,
                            qval=data$q,
                            shift=data$shift)
      }else if(data$data_type == 4){
        DRM_QE4_LN <- exp(DRM.QE4_LND(QE4resLNI[4:7], data$x, data$q, data$shift))
        # obtain loglikelihood
        llQE4LN=llfQE4_LNDc(optQE4_LNI$par[c(1,2,10,4,5,6)],
                            d=data$x,
                            n=data$n,
                            nij=data$nij,
                            y=data$y,
                            qval=data$q,
                            shift=data$shift)
      }

      llLN = c(llLN, llQE4LN)
      BMDL = c(BMDL, QE4resLNI[1])
      BMD = c(BMD, QE4resLNI[2])
      BMDU = c(BMDU, QE4resLNI[3])
    }}
  if(prior.weights[14] == 0){QE4resLNI=NULL; llLN = c(llLN,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA);
  QE4covLNI=rep(NA,2); QE4corrLNI=rep(NA,2); DRM_QE4_LN=rep(NA,length(data$x))
  parsQE4LN <- NA
  QE4outLNI <- t(data.frame(
    par.dt = rep(NA,3),
    par.is2t = rep(NA,3),
    par.k = rep(NA,3),
    par.a = rep(NA,3),
    par.b = rep(NA,3),
    par.c = rep(NA,3),
    par.d = rep(NA,3),
    par.s2 = rep(NA,3),
    BMD = rep(NA,3),
    min.resp = rep(NA,3),
    max.resp = rep(NA,3),
    fold.change = rep(NA,3), rho = rep(NA,3)
  ))
  }
  #

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[15]))
  }
  if(prior.weights[15]>0){
    # print(15)

    optP4_LNI <- fun_optimC(stanmodels$mP4c, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optP4_LNI[[3]]),TRUE,(optP4_LNI[[3]]!=0)) | length(optP4_LNI)!=9)){
      prior.weights[15] <- 0
      warning('difficulties fitting the Probit (Lognormal) model; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsP4LN <- par_extractC(optP4_LNI, model_name = "P4_LN")
      P4resLNI <- quantile(parsP4LN$BMD*data$maxD, pvec, na.rm = T)
      P4resLNI <- c(P4resLNI, apply(parsP4LN[,c(paste0("p",1:4), "is2t","rho")], 2, median))
      names(P4resLNI) <- c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t","rho")

      P4outLNI <- outLP(parsP4LN, pvec, data$maxD, clustered = T)


      # Covariance between b-d and between BMD-d
      P4covLNI = c(cov(parsP4LN[,c("b","d")], use="na.or.complete")["b","d"],
                   cov(parsP4LN[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      P4corrLNI = c(cor(parsP4LN[,c("b","d")], use="na.or.complete")["b","d"],
                    cor(parsP4LN[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      if(data$data_type == 2){
        DRM_P4_LN <- exp(DRM.P4_LNI(P4resLNI[4:7], data$x, data$q, data$shift))
        # obtain loglikelihood
        llP4LN=llfP4_LNIc(optP4_LNI$par[c(1,2,10,4,5,6)],
                          d=data$x,
                          n=data$n,
                          nij=data$nij,
                          y=data$y,
                          qval=data$q,
                          shift=data$shift)
      }else if(data$data_type == 4){
        DRM_P4_LN <- exp(DRM.P4_LND(P4resLNI[4:7], data$x, data$q, data$shift))
        # obtain loglikelihood
        llP4LN=llfP4_LNDc(optP4_LNI$par[c(1,2,10,4,5,6)],
                          d=data$x,
                          n=data$n,
                          nij=data$nij,
                          y=data$y,
                          qval=data$q,
                          shift=data$shift)
      }

      llLN = c(llLN, llP4LN)
      BMDL = c(BMDL, P4resLNI[1])
      BMD = c(BMD, P4resLNI[2])
      BMDU = c(BMDU, P4resLNI[3])
    }}
  if(prior.weights[15] == 0){P4resLNI=NULL; llLN = c(llLN,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA);
  P4covLNI=rep(NA,2); P4corrLNI=rep(NA,2); DRM_P4_LN=rep(NA,length(data$x))
  parsP4LN <- NA
  P4outLNI <- t(data.frame(
    par.dt = rep(NA,3),
    par.is2t = rep(NA,3),
    par.k = rep(NA,3),
    par.a = rep(NA,3),
    par.b = rep(NA,3),
    par.c = rep(NA,3),
    par.d = rep(NA,3),
    par.s2 = rep(NA,3),
    BMD = rep(NA,3),
    min.resp = rep(NA,3),
    max.resp = rep(NA,3),
    fold.change = rep(NA,3), rho = rep(NA,3)
  ))
  }
  #

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[16]))
  }
  if(prior.weights[16]>0){
    # print(16)

    optL4_LNI <- fun_optimC(stanmodels$mL4c, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optL4_LNI[[3]]),TRUE,(optL4_LNI[[3]]!=0)) | length(optL4_LNI)!=9)){
      prior.weights[16] <- 0
      warning('difficulties fitting the Logit (Lognormal) model; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsL4LN <- par_extractC(optL4_LNI, model_name = "L4_LN")
      L4resLNI <- quantile(parsL4LN$BMD*data$maxD, pvec, na.rm = T)
      L4resLNI <- c(L4resLNI, apply(parsL4LN[,c(paste0("p",1:4), "is2t","rho")], 2, median))
      names(L4resLNI) <- c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t","rho")

      L4outLNI <- outLP(parsL4LN, pvec, data$maxD, clustered = T)

      # Covariance between b-d and between BMD-d
      L4covLNI = c(cov(parsL4LN[,c("b","d")], use="na.or.complete")["b","d"],
                   cov(parsL4LN[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      L4corrLNI = c(cor(parsL4LN[,c("b","d")], use="na.or.complete")["b","d"],
                    cor(parsL4LN[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      if(data$data_type == 2){
        DRM_L4_LN <- exp(DRM.L4_LNI(L4resLNI[4:7], data$x, data$q, data$shift))
        # obtain loglikelihood
        llL4LN=llfL4_LNIc(optL4_LNI$par[c(1,2,10,4,5,6)],
                          d=data$x,
                          n=data$n,
                          nij=data$nij,
                          y=data$y,
                          qval=data$q,
                          shift=data$shift)
      }else if(data$data_type == 4){
        DRM_L4_LN <- exp(DRM.L4_LND(L4resLNI[4:7], data$x, data$q, data$shift))
        # obtain loglikelihood
        llL4LN=llfL4_LNDc(optL4_LNI$par[c(1,2,10,4,5,6)],
                          d=data$x,
                          n=data$n,
                          nij=data$nij,
                          y=data$y,
                          qval=data$q,
                          shift=data$shift)
      }

      llLN = c(llLN, llL4LN)
      BMDL = c(BMDL, L4resLNI[1])
      BMD = c(BMD, L4resLNI[2])
      BMDU = c(BMDU, L4resLNI[3])
    }}
  if(prior.weights[16] == 0){L4resLNI=NULL; llLN = c(llLN,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA);
  L4covLNI=rep(NA,2); L4corrLNI=rep(NA,2); DRM_L4_LN=rep(NA,length(data$x))
  parsL4LN <- NA
  L4outLNI <- t(data.frame(
    par.dt = rep(NA,3),
    par.is2t = rep(NA,3),
    par.k = rep(NA,3),
    par.a = rep(NA,3),
    par.b = rep(NA,3),
    par.c = rep(NA,3),
    par.d = rep(NA,3),
    par.s2 = rep(NA,3),
    BMD = rep(NA,3),
    min.resp = rep(NA,3),
    max.resp = rep(NA,3),
    fold.change = rep(NA,3), rho = rep(NA,3)
  ))
  }

  minll = min(llN,llLN,na.rm=T)


  # the weights

  fun.w <- function(DIH, ll, min.ll, opt, mu, sig, lb, ub, s1, s2, s3, td){

    w1 <- (2*pi)^(2.5)*sqrt(DIH)*exp(ll-min.ll)*
      # d, sigma2
      # mvtnorm::dmvnorm(opt$par[c(4,5)],mean=mu[c(4,5)],
      # sigma=sig[c(4,5),c(4,5)])*
      # sigma2
      dnorm(opt$par[5], mu[5], sig[5,5])*
      # d
      truncnorm::dtruncnorm(opt$par[4], b = td, mean = mu[4], sd = sig[4,4])*
      # a
      mc2d::dpert(opt$par[1], min = lb[1], max = ub[1],
                  mode = mu[1], shape = s1)*
      # BMD
      mc2d::dpert(opt$par[2], min = lb[2], max = ub[2],
                  mode = mu[2], shape = s2)*
      # c
      mc2d::dpert(opt$par[10], min = lb[3], max = ub[3],
                  mode = mu[3], shape = s3)*
      # rho
      mc2d::dpert(opt$par[6], min = lb[6], max = ub[6],
                  mode = mu[6], shape = 0.0001)
    return(w1)
  }

  w=c()

  # normal

  data = data.N$data
  start=data.N$start
  startQ=data.N$startQ

  lls <- c(llN, llLN)
  w.msg <- ''

  max.ll = max(lls, na.rm = T)
  if(!(1 %in% prior.weights)){
    w.msg <- 'Laplace approximation could not be performed, problem fitting all models'
    warning('Laplace approximation could not be performed, problem fitting all models')
    out.stop <- 'not.fit'
  }else if(is.na(lls[which((max.ll-lls[!is.na(lls)]) < 709 & prior.weights>0)][1])){
    lpw <- rep(0, 16)
    lpw[which(lls == max.ll)] <- 1
    w.msg <- 'Laplace weights could not be computed and one model gets all the weight; using another prior for parameter d might help'
    warning('Laplace weights could not be computed and one model gets all the weight; using another prior for parameter d might help')
  }else{

    if(FALSE %in% ((max.ll-lls[!is.na(lls)]) < 709)){
      w.msg <- 'Not all models were used in computation of Laplace weights, some models set to 0; using another prior for parameter d might help'
      warning('Not all models were used in computation of Laplace weights, some models set to 0; using another prior for parameter d might help')
    }

    minll <- min(lls[which((max.ll-lls[!is.na(lls)]) < 709 & prior.weights>0)], na.rm = T)

    if(prior.weights[1]>0){
      DIHE4h=try(det(-pracma::pinv(optE4_NI$hessian)))
      DIHE4=ifelse(DIHE4h<0 | class(DIHE4h)[1] == 'try-error',0,DIHE4h)
      w = c(w, fun.w(DIHE4, llE4N, minll, optE4_NI, data$priormu, data$priorSigma,
                     data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c,
                     data$truncd))
      # DIHE4h=try(det(-pracma::pinv(hessian(func = llE4fN,x=optE4_NI$par[c(1,2,9,3,4)],method.args=list(eps=1e-4, d=0.1, zero.tol=sqrt(.Machine$double.eps/7e-7), r=8, v=2, show.details=FALSE))))
      # DIHE4=ifelse(DIHE4h<0,0,DIHE4h)
      ## Aproximation of marginal (i.e. integrated) likelihood (= 'model evidence')

    }else{w=c(w,0)}

    if(prior.weights[2]>0){
      DIHIE4h=try(det(-pracma::pinv(optIE4_NI$hessian)))
      DIHIE4=ifelse(DIHIE4h<0 | class(DIHIE4h)[1] == 'try-error',0,DIHIE4h)
      w = c(w, fun.w(DIHIE4, llIE4N, minll, optIE4_NI, data$priormu, data$priorSigma,
                     data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c,
                     data$truncd))
    }else{w=c(w,0)}

    if(prior.weights[3]>0){
      DIHH4h=try(det(-pracma::pinv(optH4_NI$hessian)))
      DIHH4=ifelse(DIHH4h<0 | class(DIHH4h)[1] == 'try-error',0,DIHH4h)
      w = c(w, fun.w(DIHH4, llH4N, minll, optH4_NI, data$priormu, data$priorSigma,
                     data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c,
                     data$truncd))
    }else{w=c(w,0)}

    if(prior.weights[4]>0){
      DIHLN4h=try(det(-pracma::pinv(optLN4_NI$hessian)))
      DIHLN4=ifelse(DIHLN4h<0 | class(DIHLN4h)[1] == 'try-error',0,DIHLN4h)
      w = c(w, fun.w(DIHLN4, llLN4N, minll, optLN4_NI, data$priormu, data$priorSigma,
                     data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c,
                     data$truncd))
    }else{w=c(w,0)}

    if(prior.weights[5]>0){
      DIHG4h=try(det(-pracma::pinv(optG4_NI$hessian)))
      DIHG4=ifelse(DIHG4h<0 | class(DIHG4h)[1] == 'try-error',0,DIHG4h)
      w = c(w, fun.w(DIHG4, llG4N, minll, optG4_NI, data$priormu, data$priorSigma,
                     data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c,
                     data$truncd))
    }else{w=c(w,0)}

    if(prior.weights[6]>0){
      DIHQE4h=try(det(-pracma::pinv(optQE4_NI$hessian)))
      DIHQE4=ifelse(DIHQE4h<0 | class(DIHQE4h)[1] == 'try-error',0,DIHQE4h)
      w = c(w, fun.w(DIHQE4, llQE4N, minll, optQE4_NI, data$priormuQ, data$priorSigmaQ,
                     data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c,
                     data$truncdQ))
    }else{w=c(w,0)}

    if(prior.weights[7]>0){
      DIHP4h=try(det(-pracma::pinv(optP4_NI$hessian)))
      DIHP4=ifelse(DIHP4h<0 | class(DIHP4h)[1] == 'try-error',0,DIHP4h)
      w = c(w, fun.w(DIHP4, llP4N, minll, optP4_NI, data$priormu, data$priorSigma,
                     data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c,
                     data$truncd))
    }else{w=c(w,0)}

    if(prior.weights[8]>0){
      DIHL4h=try(det(-pracma::pinv(optL4_NI$hessian)))
      DIHL4=ifelse(DIHL4h<0 | class(DIHL4h)[1] == 'try-error',0,DIHL4h)
      w = c(w, fun.w(DIHL4, llL4N, minll, optL4_NI, data$priormu, data$priorSigma,
                     data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c,
                     data$truncd))
    }else{w=c(w,0)}


    # lognormal

    data=data.LN$data
    start=data.LN$start
    startQ=data.LN$startQ

    if(prior.weights[9]>0){
      DIHE4h=try(det(-pracma::pinv(optE4_LNI$hessian)))
      DIHE4=ifelse(DIHE4h<0 | class(DIHE4h)[1] == 'try-error',0,DIHE4h)
      w = c(w, fun.w(DIHE4, llE4LN, minll, optE4_LNI, data$priormu, data$priorSigma,
                     data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c,
                     data$truncd))
      # DIHE4h=try(det(-pracma::pinv(hessian(func = llE4fLN,x=optE4_LNI$par[c(1,2,9,3,4)],method.args=list(eps=1e-4, d=0.1, zero.tol=sqrt(.Machine$double.eps/7e-7), r=8, v=2, show.details=FALSE))))
      # DIHE4=ifelse(DIHE4h<0,0,DIHE4h)

    }else{w=c(w,0)}

    if(prior.weights[10]>0){
      DIHIE4h=try(det(-pracma::pinv(optIE4_LNI$hessian)))
      DIHIE4=ifelse(DIHIE4h<0 | class(DIHIE4h)[1] == 'try-error',0,DIHIE4h)
      w = c(w, fun.w(DIHIE4, llIE4LN, minll, optIE4_LNI, data$priormu, data$priorSigma,
                     data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c,
                     data$truncd))
    }else{w=c(w,0)}

    if(prior.weights[11]>0){
      DIHH4h=try(det(-pracma::pinv(optH4_LNI$hessian)))
      DIHH4=ifelse(DIHH4h<0 | class(DIHH4h)[1] == 'try-error',0,DIHH4h)
      w = c(w, fun.w(DIHH4, llH4LN, minll, optH4_LNI, data$priormu, data$priorSigma,
                     data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c,
                     data$truncd))
    }else{w=c(w,0)}

    if(prior.weights[12]>0){
      DIHLN4h=try(det(-pracma::pinv(optLN4_LNI$hessian)))
      DIHLN4=ifelse(DIHLN4h<0 | class(DIHLN4h)[1] == 'try-error',0,DIHLN4h)
      w = c(w, fun.w(DIHLN4, llLN4LN, minll, optLN4_LNI, data$priormu, data$priorSigma,
                     data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c,
                     data$truncd))
    }else{w=c(w,0)}

    if(prior.weights[13]>0){
      DIHG4h=try(det(-pracma::pinv(optG4_LNI$hessian)))
      DIHG4=ifelse(DIHG4h<0 | class(DIHG4h)[1] == 'try-error',0,DIHG4h)
      w = c(w, fun.w(DIHG4, llG4LN, minll, optG4_LNI, data$priormu, data$priorSigma,
                     data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c,
                     data$truncd))
    }else{w=c(w,0)}

    if(prior.weights[14]>0){
      DIHQE4h=try(det(-pracma::pinv(optQE4_LNI$hessian)))
      DIHQE4=ifelse(DIHQE4h<0 | class(DIHQE4h)[1] == 'try-error',0,DIHQE4h)
      w = c(w, fun.w(DIHQE4, llQE4LN, minll, optQE4_LNI, data$priormuQ, data$priorSigmaQ,
                     data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c,
                     data$truncdQ))
    }else{w=c(w,0)}

    if(prior.weights[15]>0){
      DIHP4h=try(det(-pracma::pinv(optP4_LNI$hessian)))
      DIHP4=ifelse(DIHP4h<0 | class(DIHP4h)[1] == 'try-error',0,DIHP4h)
      w = c(w, fun.w(DIHP4, llP4LN, minll, optP4_LNI, data$priormu, data$priorSigma,
                     data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c,
                     data$truncd))
    }else{w=c(w,0)}

    if(prior.weights[16]>0){
      DIHL4h=try(det(-pracma::pinv(optL4_LNI$hessian)))
      DIHL4=ifelse(DIHL4h<0 | class(DIHL4h)[1] == 'try-error',0,DIHL4h)
      w = c(w, fun.w(DIHL4, llL4LN, minll, optL4_LNI, data$priormu, data$priorSigma,
                     data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c,
                     data$truncd))
    }else{w=c(w,0)}

    w <- ifelse(w == 'Inf' | is.na(w), 0, w)
    prior.weights = prior.weights/sum(prior.weights==1)
    lpw=(prior.weights*w)/sum(prior.weights*w)

  }


  if(out.stop != 'not.fit'){
    # the model average posterior as a mixture
    count=round(lpw*ndraws)
    mabmd=(c(# normal
      if(prior.weights[1]>0) sample(optE4_NI$theta_tilde[,2],count[1],replace=T),
      if(prior.weights[2]>0) sample(optIE4_NI$theta_tilde[,2],count[2],replace=T),
      if(prior.weights[3]>0) sample(optH4_NI$theta_tilde[,2],count[3],replace=T),
      if(prior.weights[4]>0) sample(optLN4_NI$theta_tilde[,2],count[4],replace=T),
      if(prior.weights[5]>0) sample(optG4_NI$theta_tilde[,2],count[5],replace=T),
      if(prior.weights[6]>0) sample(optQE4_NI$theta_tilde[,2],count[6],replace=T),
      if(prior.weights[7]>0) sample(optP4_NI$theta_tilde[,2],count[7],replace=T),
      if(prior.weights[8]>0) sample(optL4_NI$theta_tilde[,2],count[8],replace=T),
      # lognormal
      if(prior.weights[9]>0) sample(optE4_LNI$theta_tilde[,2],count[9],replace=T),
      if(prior.weights[10]>0) sample(optIE4_LNI$theta_tilde[,2],count[10],replace=T),
      if(prior.weights[11]>0) sample(optH4_LNI$theta_tilde[,2],count[11],replace=T),
      if(prior.weights[12]>0) sample(optLN4_LNI$theta_tilde[,2],count[12],replace=T),
      if(prior.weights[13]>0) sample(optG4_LNI$theta_tilde[,2],count[13],replace=T),
      if(prior.weights[14]>0) sample(optQE4_LNI$theta_tilde[,2],count[14],replace=T),
      if(prior.weights[15]>0) sample(optP4_LNI$theta_tilde[,2],count[15],replace=T),
      if(prior.weights[16]>0) sample(optL4_LNI$theta_tilde[,2],count[16],replace=T)
    ))
    maci=quantile(mabmd,pvec, na.rm = T)*data$maxD ## original scale
    names(maci)=c("BMDL","BMD","BMDU")

    if(TRUE %in% (mabmd > data$maxD)  && data$maxD > 1){
      mabmd = ifelse(mabmd > data$maxD, data$maxD, mabmd)
      p.msg = 'The model averaged posterior distribution has been truncated at max(Dose)^2'
      warnings('The model averaged posterior distribution has been truncated at max(Dose)^2')
    }else{
      p.msg = ''
    }

    BMDq = quantile(mabmd, seq(0,1,0.005), na.rm = T)*data$maxD ## original scale

    BMDL = c(BMDL, maci[1]); BMD = c(BMD, maci[2]); BMDU = c(BMDU, maci[3])

    names(BMDL) <- c("E4_N","IE4_N","H4_N","LN4_N","G4_N","QE4_N","P4_N","L4_N","E4_LN","IE4_LN","H4_LN",
                     "LN4_LN","G4_LN","QE4_LN","P4_LN","L4_LN","MA")
    model = c("E4_N","IE4_N","H4_N","LN4_N","G4_N","QE4_N","P4_N","L4_N","E4_LN","IE4_LN","H4_LN",
              "LN4_LN","G4_LN","QE4_LN","P4_LN","L4_LN","MA")
    model = as.factor(model)
    weight = c(lpw[1], lpw[2], lpw[3], lpw[4], lpw[5], lpw[6], lpw[7], lpw[8],
               lpw[9], lpw[10], lpw[11], lpw[12], lpw[13], lpw[14], lpw[15], lpw[16], 1)


    names(lpw) = model[1:16]

    ### Model-averaged response per dose level
    dr.MA <- c()
    for(i in 1:length(data$x)){
      dr.MA[i] = weighted.mean(x = c(DRM_E4_N[i],DRM_IE4_N[i],DRM_H4_N[i],DRM_LN4_N[i],DRM_G4_N[i],DRM_QE4_N[i],
                                     DRM_P4_N[i], DRM_L4_N[i] ,
                                     DRM_E4_LN[i], DRM_IE4_LN[i], DRM_H4_LN[i], DRM_LN4_LN[i], DRM_G4_LN[i],
                                     DRM_QE4_LN[i], DRM_P4_LN[i],DRM_L4_LN[i]),
                               w = lpw,
                               na.rm = T)
    }


    ## Covariances
    covs = t(data.frame(
      E4_N = E4covNI,
      IE4_N = IE4covNI,
      H4_N = H4covNI,
      LN4_N = LN4covNI,
      G4_N = G4covNI,
      QE4_N = QE4covNI,
      P4_N = P4covNI,
      L4_N = L4covNI,
      E4_LN = E4covLNI,
      IE4_LN = IE4covLNI,
      H4_LN = H4covLNI,
      LN4_LN = LN4covLNI,
      G4_LN = G4covLNI,
      QE4_LN = QE4covLNI,
      P4_LN = P4covLNI,
      L4_LN = L4covLNI
    ))
    colnames(covs) = c("b-d", "BMD-d")

    corrs = t(data.frame(
      E4_N = E4corrNI,
      IE4_N = IE4corrNI,
      H4_N = H4corrNI,
      LN4_N = LN4corrNI,
      G4_N = G4corrNI,
      QE4_N = QE4corrNI,
      P4_N = P4corrNI,
      L4_N = L4corrNI,
      E4_LN = E4corrLNI,
      IE4_LN = IE4corrLNI,
      H4_LN = H4corrLNI,
      LN4_LN = LN4corrLNI,
      G4_LN = G4corrLNI,
      QE4_LN = QE4corrLNI,
      P4_LN = P4corrLNI,
      L4_LN = L4corrLNI
    ))
    colnames(corrs) = c("b-d", "BMD-d")

    modelnames = c("E4_N","IE4_N","H4_N","LN4_N","G4_N","QE4_N","P4_N","L4_N","E4_LN","IE4_LN","H4_LN","LN4_LN","G4_LN","QE4_LN","P4_LN","L4_LN")

    ### Some additional checks

    if(maci[2]/maci[1] > 20){
      warning('BMD/BMDL is larger than 20')
    }
    if(maci[3]/maci[1] > 50){
      warning('BMDU/BMDL is larger than 50')
    }
    if(maci[2] < (data.N$data$x[2]*data.N$data$maxD/10)){
      warning('BMD is 10 times lower than the lowest non-zero dose')
    }

    ### best fitting model vs saturated ANOVA model
    best.fit = modelnames[which(weight[1:16] == max(weight[1:16]))][1]
    nrchains = 3; nriterations = 3000; warmup = 1000; delta = 0.8; treedepth = 10
    bfTest <- modelTestC(best.fit, data.N, data.LN, get(paste0('opt', best.fit, 'I')), type = 'Laplace',
                         seed, ndraws, nrchains, nriterations, warmup, delta, treedepth)
    print(warning(bfTest$warn.bf))
  }

  if(out.stop == 'not.fit'){
    ret_results <- list(MA = c(NA,NA,NA),
                        MA_post = NA,
                        bkg_post = NA,
                        maxy_post = NA)
  }else{

    ret_results <- list(
      # model parameters
      E4_N=E4outNI,IE4_N=IE4outNI,H4_N=H4outNI,LN4_N=LN4outNI,G4_N=G4outNI,QE4_N=QE4outNI,P4_N=P4outNI,L4_N=L4outNI,
      E4_LN=E4outLNI,IE4_LN=IE4outLNI,H4_LN=H4outLNI,LN4_LN=LN4outLNI,G4_LN=G4outLNI,QE4_LN=QE4outLNI,
      P4_LN=P4outLNI,L4_LN=L4outLNI,
      # covariances
      covs = covs, corrs = corrs,
      # weights and MA
      weights=lpw,
      MA=maci,
      llN=llN, llLN=llLN, MA_post = BMDq,
      # dose-response
      MA_dr = dr.MA,
      parsN = list(parsE4N, parsIE4N, parsH4N, parsLN4N, parsG4N,
                   parsQE4N, parsP4N, parsL4N),
      parsLN = list(parsE4LN, parsIE4LN, parsH4LN, parsLN4LN, parsG4LN,
                    parsQE4LN, parsP4LN, parsL4LN),
      BMDMixture = (mabmd)*data$maxD,
      # data = data.frame(
      #   dose = c(data.N$data$x),
      #   sd = sqrt(data.N$data$s2),
      #   m = data.N$data$m),
      data = data.N$data$data,
      max.dose = data.N$data$maxD,
      q = data.N$data$q,
      # increasing = T,
      models_included_laplace = modelnames[prior.weights > 0],
      bf = bfTest$bayesFactor, gof_check = bfTest$warn.bf,
      # means.SM = bfTest$means.SM, parBestFit = bfTest$par.best,
      # BIC.bestfit = bfTest$BIC.bestfit, BIC.SM = bfTest$BIC.SM,
      shift = data.LN$data$shift,
      w.msg = w.msg, p.msg = p.msg
    )

    attr(ret_results, "class") <- c("BMADR", "LP")
  }


  return(ret_results)

}


#' @rdname full.laplace_MA
#' @export
full.laplaceQ_MA=function(data.Q, prior.weights = rep(1, 8),
                          ndraws=30000,seed=123,pvec=c(0.05,0.5,0.95)){

  # prior.weights = prior.weights/sum(prior.weights==1) #this is now done below when calculating weights

  out.stop = 'ok'

  ## Data to use for Normal distribution
  data = data.Q$data
  start = data.Q$start
  startQ = data.Q$startQ

  llQ = c() # likelihoods
  BMDL = c()
  BMD = c()
  BMDU = c()

  ### Getting the posterior modes and the hessian to obtain the model specific posterior distributions
  # optimizing function: obtain point estimates by maximizing the joint posterior from the Stan model
  ### Additionally save BMD(L/U) estimates and loglikelihood values; set to NA if model has prior weight 0

  nModels = 8

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='quantal')[1]))
  }
  if(prior.weights[1]>0){

    # print(1)
    optE4_Q <- fun_optimQ(stanmodels$mE4_Q, data, start, ndraws, seed, pvec)

    if((ifelse(is.na(optE4_Q[[3]]),TRUE,(optE4_Q[[3]]!=0)) | length(optE4_Q)!=9)){
      prior.weights[1] <- 0
      warning('difficulties fitting the Exponential model; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      if(data$is_bin == 1) {
        parsE4Q <- parq_extract(optE4_Q, model_name = "E4_Q", pars = c('a', 'b', 'd', 'BMD',
                                                                       paste0('par',1:3)))
      } else {
        parsE4Q <- parq_extract(optE4_Q, model_name = "E4_Q", pars = c('a', 'b', 'd', 'BMD','rho[1]',
                                                                       paste0('par',1:3)),
                                rho = TRUE)
      }

      E4resQ <- quantile(parsE4Q$BMD, pvec, na.rm = T)*data$maxD
      if(data$is_bin == 1){
        E4resQ <- c(E4resQ, apply(parsE4Q[,c('a', 'b', 'd')], 2, median))
      } else {
        E4resQ <- c(E4resQ, apply(parsE4Q[,c('a', 'b', 'd','rho')], 2, median))
      }

      names(E4resQ) <- ifelse(rep(data$is_bin == 1, ifelse(data$is_bin == 1,6,7)),
                              c("BMDL","BMD","BMDU","a","b","d"),
                              c("BMDL","BMD","BMDU","a","b","d","rho")
      )

      if(data$is_bin == 1){
        E4outQ <- outLPQ(parsE4Q, pvec, data$maxD)
      } else {
        E4outQ <- outLPQ(parsE4Q, pvec, data$maxD, rho=TRUE)
      }

      # Covariance between b-d and between BMD-d
      E4covQ <- c(cov(parsE4Q[,c("b","d")], use="na.or.complete")["b","d"],
                  cov(parsE4Q[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      E4corrQ <- c(cor(parsE4Q[,c("b","d")], use="na.or.complete")["b","d"],
                   cor(parsE4Q[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      DRM_E4_Q <- DRM.E4_Q(optE4_Q$par[1:3], data$x, data$q)

      # obtain loglikelihood
      llE4Q <- ifelse(data$is_bin == 1,
                      llfE4_Q(optE4_Q$par[1:3],
                              nvec=data$n,
                              dvec=data$x,
                              yvec=data$y,
                              qval=data$q),
                      llfE42_Q(optE4_Q$par[1:3],
                               nvec=data$n,
                               dvec=data$x,
                               yvec=data$y,
                               qval=data$q,
                               rho = optE4_Q$par[stringr::str_detect(names(optE4_Q$par),'rho') &
                                                   !stringr::str_detect(names(optE4_Q$par),'eta')])
      )
      # save LL and BMD(L/U) estimate
      llQ = c(llQ, llE4Q)
      BMDL = c(BMDL, E4resQ[1])
      BMD = c(BMD, E4resQ[2])
      BMDU = c(BMDU, E4resQ[3])
    }
  }
  if(prior.weights[1] == 0){E4resQ=NULL; llQ = c(llQ,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA);
  E4covQ=rep(NA,2); E4corrQ=rep(NA,2); DRM_E4_Q=rep(NA, length(data$x))
  parsE4Q <- NA
  E4outQ <- t(data.frame(
    # transformed parameters
    par.at = rep(NA,3),
    par.dt = rep(NA,3),
    par.k = rep(NA,3),
    # parameters on original scale
    par.a = rep(NA,3),
    par.b = rep(NA,3),
    par.c = rep(NA,3),
    par.d = rep(NA,3),
    # natural parameters
    BMD = rep(NA,3),
    min.resp = rep(NA,3)
  ))
  }

  #
  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='quantal')[2]))
  }
  if(prior.weights[2]>0){

    # print(2)
    optIE4_Q <- fun_optimQ(stanmodels$mIE4_Q, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optIE4_Q[[3]]),TRUE,(optIE4_Q[[3]]!=0)) | length(optIE4_Q)!=9)){
      prior.weights[2] <- 0
      warning('difficulties fitting the Inverse Exponential model; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      if(data$is_bin == 1) {
        parsIE4Q <- parq_extract(optIE4_Q, model_name = "IE4_Q", pars = c('a', 'b', 'd', 'BMD',
                                                                          paste0('par',1:3)))
      } else {
        parsIE4Q <- parq_extract(optIE4_Q, model_name = "IE4_Q", pars = c('a', 'b', 'd', 'BMD','rho[1]',
                                                                          paste0('par',1:3)),
                                 rho = TRUE)
      }
      IE4resQ <- quantile(parsIE4Q$BMD, pvec, na.rm = T)*data$maxD  #exp(quantile(optE4_Q$theta_tilde[,"par[2]"],pvec))

      if(data$is_bin == 1){
        IE4resQ <- c(IE4resQ, apply(parsIE4Q[,c('a', 'b', 'd')], 2, median))
      } else {
        IE4resQ <- c(IE4resQ, apply(parsIE4Q[,c('a', 'b', 'd','rho')], 2, median))
      }

      names(IE4resQ) <- ifelse(rep(data$is_bin == 1, ifelse(data$is_bin == 1,6,7)),
                               c("BMDL","BMD","BMDU","a","b","d"),
                               c("BMDL","BMD","BMDU","a","b","d","rho")
      )

      if(data$is_bin == 1){
        IE4outQ <- outLPQ(parsIE4Q, pvec, data$maxD)
      } else {
        IE4outQ <- outLPQ(parsIE4Q, pvec, data$maxD, rho = TRUE)
      }

      # Covariance between b-d and between BMD-d
      IE4covQ <- c(cov(parsIE4Q[,c("b","d")], use="na.or.complete")["b","d"],
                   cov(parsIE4Q[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      IE4corrQ <- c(cor(parsIE4Q[,c("b","d")], use="na.or.complete")["b","d"],
                    cor(parsIE4Q[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      DRM_IE4_Q <- DRM.IE4_Q(optIE4_Q$par[1:3], data$x, data$q)

      llIE4Q <- ifelse(data$is_bin == 1,
                       llfIE4_Q(optIE4_Q$par[1:3],nvec=data$n,
                                dvec=data$x,
                                yvec=data$y,
                                qval=data$q),
                       llfIE42_Q(optIE4_Q$par[1:3],
                                 nvec=data$n,
                                 dvec=data$x,
                                 yvec=data$y,
                                 qval=data$q,
                                 rho = optIE4_Q$par[stringr::str_detect(names(optIE4_Q$par),'rho') &
                                                      !stringr::str_detect(names(optIE4_Q$par),'eta')])
      )

      llQ = c(llQ, llIE4Q)
      BMDL = c(BMDL, IE4resQ[1])
      BMD = c(BMD, IE4resQ[2])
      BMDU = c(BMDU, IE4resQ[3])
    }
  }
  if(prior.weights[2] == 0){
    IE4resQ=NULL; llQ = c(llQ,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA);
    IE4covQ=rep(NA,2); IE4corrQ=rep(NA,2); DRM_IE4_Q=rep(NA,length(data$x))
    parsIE4Q <- NA
    IE4outQ <- t(data.frame(
      # transformed parameters
      par.at = rep(NA,3),
      par.dt = rep(NA,3),
      par.k = rep(NA,3),
      # parameters on original scale
      par.a = rep(NA,3),
      par.b = rep(NA,3),
      par.c = rep(NA,3),
      par.d = rep(NA,3),
      # natural parameters
      BMD = rep(NA,3),
      min.resp = rep(NA,3)
    ))
  }
  #

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='quantal')[3]))
  }
  if(prior.weights[3]>0){

    # print(3)

    optH4_Q <- fun_optimQ(stanmodels$mH4_Q, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optH4_Q[[3]]),TRUE,(optH4_Q[[3]]!=0)) | length(optH4_Q)!=9)){
      prior.weights[3] <- 0
      warning('difficulties fitting the Hill model; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      if(data$is_bin == 1) {
        parsH4Q <- parq_extract(optH4_Q, model_name = "H4_Q", pars = c('a', 'b', 'd', 'BMD',
                                                                       paste0('par',1:3)))
      } else {
        parsH4Q <- parq_extract(optH4_Q, model_name = "H4_Q", pars = c('a', 'b', 'd', 'BMD','rho[1]',
                                                                       paste0('par',1:3)),
                                rho = TRUE)
      }

      H4resQ <- quantile(parsH4Q$BMD, pvec, na.rm = T)*data$maxD
      if(data$is_bin == 1){
        H4resQ <- c(H4resQ, apply(parsH4Q[,c('a', 'b', 'd')], 2, median))
      } else {
        H4resQ <- c(H4resQ, apply(parsH4Q[,c('a', 'b', 'd','rho')], 2, median))
      }

      names(H4resQ) <- ifelse(rep(data$is_bin == 1, ifelse(data$is_bin == 1,6,7)),
                              c("BMDL","BMD","BMDU","a","b","d"),
                              c("BMDL","BMD","BMDU","a","b","d","rho")
      )

      if(data$is_bin == 1){
        H4outQ <- outLPQ(parsH4Q, pvec, data$maxD)
      } else {
        H4outQ <- outLPQ(parsH4Q, pvec, data$maxD, rho = TRUE)
      }


      # Covariance between b-d and between BMD-d
      H4covQ = c(cov(parsH4Q[,c("b","d")], use="na.or.complete")["b","d"],
                 cov(parsH4Q[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      H4corrQ = c(cor(parsH4Q[,c("b","d")], use="na.or.complete")["b","d"],
                  cor(parsH4Q[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      DRM_H4_Q <- DRM.H4_Q(optH4_Q$par[1:3], data$x, data$q)

      llH4Q <- ifelse(data$is_bin == 1,
                      llfH4_Q(optH4_Q$par[1:3],nvec=data$n,
                              dvec=data$x,
                              yvec=data$y,
                              qval=data$q),
                      llfH42_Q(optH4_Q$par[1:3],
                               nvec=data$n,
                               dvec=data$x,
                               yvec=data$y,
                               qval=data$q,
                               rho = optH4_Q$par[stringr::str_detect(names(optH4_Q$par),'rho') &
                                                   !stringr::str_detect(names(optH4_Q$par),'eta')])
      )

      llQ = c(llQ, llH4Q)
      BMDL = c(BMDL, H4resQ[1])
      BMD = c(BMD, H4resQ[2])
      BMDU = c(BMDU, H4resQ[3])
    }
  }
  if(prior.weights[3] == 0){H4resQ=NULL; llQ = c(llQ,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA);
  BMDU=c(BMDU,NA); H4covQ=rep(NA,2); H4corrQ=rep(NA,2); DRM_H4_Q=rep(NA,length(data$x))
  parsH4Q <- NA
  H4outQ <- t(data.frame(
    # transformed parameters
    par.at = rep(NA,3),
    par.dt = rep(NA,3),
    par.k = rep(NA,3),
    # parameters on original scale
    par.a = rep(NA,3),
    par.b = rep(NA,3),
    par.c = rep(NA,3),
    par.d = rep(NA,3),
    # natural parameters
    BMD = rep(NA,3),
    min.resp = rep(NA,3)
  ))
  }

  #
  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='quantal')[4]))
  }
  if(prior.weights[4]>0){

    # print(4)
    optLN4_Q <- fun_optimQ(stanmodels$mLN4_Q, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optLN4_Q[[3]]),TRUE,(optLN4_Q[[3]]!=0)) | length(optLN4_Q)!=9)){
      prior.weights[4] <- 0
      warning('difficulties fitting the Lognormal model; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      if(data$is_bin == 1) {
        parsLN4Q <- parq_extract(optLN4_Q, model_name = "LN4_Q", pars = c('a', 'b', 'd', 'BMD',
                                                                          paste0('par',1:3)))
      } else {
        parsLN4Q <- parq_extract(optLN4_Q, model_name = "LN4_Q", pars = c('a', 'b', 'd', 'BMD','rho[1]',
                                                                          paste0('par',1:3)),
                                 rho = TRUE)
      }
      LN4resQ <- quantile(parsLN4Q$BMD, pvec, na.rm = T)*data$maxD  #exp(quantile(optE4_Q$theta_tilde[,"par[2]"],pvec))

      if(data$is_bin == 1){
        LN4resQ <- c(LN4resQ, apply(parsLN4Q[,c('a', 'b', 'd')], 2, median))
      } else {
        LN4resQ <- c(LN4resQ, apply(parsLN4Q[,c('a', 'b', 'd','rho')], 2, median))
      }

      names(LN4resQ) <- ifelse(rep(data$is_bin == 1, ifelse(data$is_bin == 1,6,7)),
                               c("BMDL","BMD","BMDU","a","b","d"),
                               c("BMDL","BMD","BMDU","a","b","d","rho")
      )

      if(data$is_bin == 1){
        LN4outQ <- outLPQ(parsLN4Q, pvec, data$maxD)
      } else {
        LN4outQ <- outLPQ(parsLN4Q, pvec, data$maxD, rho = TRUE)
      }

      # Covariance between b-d and between BMD-d
      LN4covQ = c(cov(parsLN4Q[,c("b","d")], use="na.or.complete")["b","d"],
                  cov(parsLN4Q[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      LN4corrQ = c(cor(parsLN4Q[,c("b","d")], use="na.or.complete")["b","d"],
                   cor(parsLN4Q[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      DRM_LN4_Q = DRM.LN4_Q(optLN4_Q$par[1:3], data$x, data$q)


      llLN4Q <- ifelse(data$is_bin == 1,
                       llfLN4_Q(optLN4_Q$par[1:3],nvec=data$n,
                                dvec=data$x,
                                yvec=data$y,
                                qval=data$q),
                       llfLN42_Q(optLN4_Q$par[1:3],
                                 nvec=data$n,
                                 dvec=data$x,
                                 yvec=data$y,
                                 qval=data$q,
                                 rho = optLN4_Q$par[stringr::str_detect(names(optLN4_Q$par),'rho') &
                                                      !stringr::str_detect(names(optLN4_Q$par),'eta')])
      )

      llQ = c(llQ, llLN4Q)
      BMDL = c(BMDL, LN4resQ[1])
      BMD = c(BMD, LN4resQ[2])
      BMDU = c(BMDU, LN4resQ[3])
    }
  }
  if(prior.weights[4] == 0){LN4resQ=NULL; llQ = c(llQ,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA);
  LN4covQ=rep(NA,2); LN4corrQ=rep(NA,2); DRM_LN4_Q=rep(NA,length(data$x))
  parsLN4Q <- NA
  LN4outQ <- t(data.frame(
    # transformed parameters
    par.at = rep(NA,3),
    par.dt = rep(NA,3),
    par.k = rep(NA,3),
    # parameters on original scale
    par.a = rep(NA,3),
    par.b = rep(NA,3),
    par.c = rep(NA,3),
    par.d = rep(NA,3),
    # natural parameters
    BMD = rep(NA,3),
    min.resp = rep(NA,3)
  ))
  }
  #

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='quantal')[5]))
  }
  if(prior.weights[5]>0){
    # print(5)

    # data$init_b <- qgamma(data$q, rate=1.0, shape=optE4_Q$par[6])/optE4_Q$par[2]

    optG4_Q <- fun_optimQ(stanmodels$mG4_Q, data, start, ndraws, 123, pvec)
    if((ifelse(is.na(optG4_Q[[3]]),TRUE,(optG4_Q[[3]]!=0)) | length(optG4_Q)!=9)){
      prior.weights[5] <- 0
      warning('difficulties fitting the Gamma model; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      if(data$is_bin == 1) {
        parsG4Q <- parq_extract(optG4_Q, model_name = "G4_Q", pars = c('a', 'b', 'd', 'BMD',
                                                                       paste0('par',1:3)))
      } else {
        parsG4Q <- parq_extract(optG4_Q, model_name = "G4_Q", pars = c('a', 'b', 'd', 'BMD','rho[1]',
                                                                       paste0('par',1:3)),
                                rho = TRUE)
      }

      G4resQ <- quantile(parsG4Q$BMD, pvec, na.rm = T)*data$maxD
      if(data$is_bin == 1){
        G4resQ <- c(G4resQ, apply(parsG4Q[,c('a', 'b', 'd')], 2, median))
      } else {
        G4resQ <- c(G4resQ, apply(parsG4Q[,c('a', 'b', 'd','rho')], 2, median))
      }

      names(G4resQ) <- ifelse(rep(data$is_bin == 1, ifelse(data$is_bin == 1,6,7)),
                              c("BMDL","BMD","BMDU","a","b","d"),
                              c("BMDL","BMD","BMDU","a","b","d","rho")
      )

      if(data$is_bin == 1){
        G4outQ <- outLPQ(parsG4Q, pvec, data$maxD)
      } else {
        G4outQ <- outLPQ(parsG4Q, pvec, data$maxD, rho = TRUE)
      }

      # Covariance between b-d and between BMD-d
      G4covQ = c(cov(parsG4Q[,c("b","d")], use="na.or.complete")["b","d"],
                 cov(parsG4Q[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      G4corrQ = c(cor(parsG4Q[,c("b","d")], use="na.or.complete")["b","d"],
                  cor(parsG4Q[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      DRM_G4_Q = DRM.G4_Q(optG4_Q$par[1:3], data$x, data$q)

      llG4Q <- ifelse(data$is_bin == 1,
                      llfG4_Q(optG4_Q$par[1:3],nvec=data$n,
                              dvec=data$x,
                              yvec=data$y,
                              qval=data$q),
                      llfG42_Q(optG4_Q$par[1:3],
                               nvec=data$n,
                               dvec=data$x,
                               yvec=data$y,
                               qval=data$q,
                               rho = optG4_Q$par[stringr::str_detect(names(optG4_Q$par),'rho') &
                                                   !stringr::str_detect(names(optG4_Q$par),'eta')])
      )

      llQ = c(llQ, llG4Q)
      BMDL = c(BMDL, G4resQ[1])
      BMD = c(BMD, G4resQ[2])
      BMDU = c(BMDU, G4resQ[3])
    }
  }
  if(prior.weights[5] == 0){G4resQ=NULL; llQ = c(llQ,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA);
  BMDU=c(BMDU,NA); G4covQ=rep(NA,2); G4corrQ=rep(NA,2); DRM_G4_Q=rep(NA, length(data$x))
  parsG4Q <-NA
  G4outQ <- t(data.frame(
    # transformed parameters
    par.at = rep(NA,3),
    par.dt = rep(NA,3),
    par.k = rep(NA,3),
    # parameters on original scale
    par.a = rep(NA,3),
    par.b = rep(NA,3),
    par.c = rep(NA,3),
    par.d = rep(NA,3),
    # natural parameters
    BMD = rep(NA,3),
    min.resp = rep(NA,3)
  ))
  }
  #

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='quantal')[6]))
  }
  if(prior.weights[6]>0){
    # print(6)
    optQE4_Q = fun_optimQ(stanmodels$mQE4_Q, data, startQ, ndraws, 123, pvec)

    if((ifelse(is.na(optQE4_Q[[3]]),TRUE,(optQE4_Q[[3]]!=0)) | length(optQE4_Q)!=9)){
      prior.weights[6] <- 0
      warning('difficulties fitting the Quadratic Exponential model; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      if(data$is_bin == 1) {
        parsQE4Q <- parq_extract(optQE4_Q, model_name = "QE4_Q", pars = c('a', 'b', 'd', 'BMD',
                                                                          paste0('par',1:3)))
      } else {
        parsQE4Q <- parq_extract(optQE4_Q, model_name = "QE4_Q", pars = c('a', 'b', 'd', 'BMD','rho[1]',
                                                                          paste0('par',1:3)),
                                 rho = TRUE)
      }
      QE4resQ <- quantile(parsQE4Q$BMD, pvec, na.rm = T)*data$maxD  #exp(quantile(optE4_Q$theta_tilde[,"par[2]"],pvec))

      if(data$is_bin == 1){
        QE4resQ <- c(QE4resQ, apply(parsQE4Q[,c('a', 'b', 'd')], 2, median))
      } else {
        QE4resQ <- c(QE4resQ, apply(parsQE4Q[,c('a', 'b', 'd','rho')], 2, median))
      }

      names(QE4resQ) <- ifelse(rep(data$is_bin == 1, ifelse(data$is_bin == 1,6,7)),
                               c("BMDL","BMD","BMDU","a","b","d"),
                               c("BMDL","BMD","BMDU","a","b","d","rho")
      )

      if(data$is_bin == 1){
        QE4outQ <- outLPQ(parsQE4Q, pvec, data$maxD)
      } else {
        QE4outQ <- outLPQ(parsQE4Q, pvec, data$maxD, rho = TRUE)
      }

      # Covariance between b-d and between BMD-d
      QE4covQ = c(cov(parsQE4Q[,c("b","d")], use="na.or.complete")["b","d"],
                  cov(parsQE4Q[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      QE4corrQ = c(cor(parsQE4Q[,c("b","d")], use="na.or.complete")["b","d"],
                   cor(parsQE4Q[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      DRM_QE4_Q = DRM.QE4_Q(optQE4_Q$par[1:3], data$x, data$q)


      llQE4Q = ifelse(data$is_bin == 1,
                      llfQE4_Q(optQE4_Q$par[1:3],nvec=data$n,
                               dvec=data$x,
                               yvec=data$y,
                               qval=data$q),
                      llfQE42_Q(optQE4_Q$par[1:3],nvec=data$n,
                                dvec=data$x,
                                yvec=data$y,
                                qval=data$q,
                                rho = optQE4_Q$par[stringr::str_detect(names(optQE4_Q$par),'rho') &
                                                     !stringr::str_detect(names(optQE4_Q$par),'eta')])
      )

      llQ = c(llQ, llQE4Q)
      BMDL = c(BMDL, QE4resQ[1])
      BMD = c(BMD, QE4resQ[2])
      BMDU = c(BMDU, QE4resQ[3])
    }
  }
  if(prior.weights[6] == 0){QE4resQ=NULL; llQ = c(llQ,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA);
  BMDU = c(BMDU,NA); QE4covQ = rep(NA,2); QE4corrQ=rep(NA,2); DRM_QE4_Q = rep(NA,length(data$x))
  parsQE4Q <- NA
  QE4outQ <- t(data.frame(
    # transformed parameters
    par.at = rep(NA,3),
    par.dt = rep(NA,3),
    par.k = rep(NA,3),
    # parameters on original scale
    par.a = rep(NA,3),
    par.b = rep(NA,3),
    par.c = rep(NA,3),
    par.d = rep(NA,3),
    # natural parameters
    BMD = rep(NA,3),
    min.resp = rep(NA,3)
  ))
  }
  #

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='quantal')[7]))
  }
  if(prior.weights[7]>0){
    # print(7)
    optP4_Q = fun_optimQ(stanmodels$mP4_Q, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optP4_Q[[3]]),TRUE,(optP4_Q[[3]]!=0)) | length(optP4_Q)!=9)){
      prior.weights[7] <- 0
      warning('difficulties fitting the Probit model; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      if(data$is_bin == 1) {
        parsP4Q <- parq_extract(optP4_Q, model_name = "P4_Q", pars = c('a', 'b', 'd', 'BMD',
                                                                       paste0('par',1:3)))
      } else {
        parsP4Q <- parq_extract(optP4_Q, model_name = "P4_Q", pars = c('a', 'b', 'd', 'BMD','rho[1]',
                                                                       paste0('par',1:3)),
                                rho = TRUE)
      }

      P4resQ <- quantile(parsP4Q$BMD, pvec, na.rm = T)*data$maxD
      if(data$is_bin == 1){
        P4resQ <- c(P4resQ, apply(parsP4Q[,c('a', 'b', 'd')], 2, median))
      } else {
        P4resQ <- c(P4resQ, apply(parsP4Q[,c('a', 'b', 'd','rho')], 2, median))
      }

      names(P4resQ) <- ifelse(rep(data$is_bin == 1, ifelse(data$is_bin == 1,6,7)),
                              c("BMDL","BMD","BMDU","a","b","d"),
                              c("BMDL","BMD","BMDU","a","b","d","rho")
      )

      if(data$is_bin == 1){
        P4outQ <- outLPQ(parsP4Q, pvec, data$maxD)
      } else {
        P4outQ <- outLPQ(parsP4Q, pvec, data$maxD, rho = TRUE)
      }


      # Covariance between b-d and between BMD-d
      P4covQ = c(cov(parsP4Q[,c("b","d")], use="na.or.complete")["b","d"],
                 cov(parsP4Q[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      P4corrQ = c(cor(parsP4Q[,c("b","d")], use="na.or.complete")["b","d"],
                  cor(parsP4Q[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      DRM_P4_Q = DRM.P4_Q(optP4_Q$par[1:3], data$x, data$q)

      llP4Q <- ifelse(data$is_bin == 1,
                      llfP4_Q(optP4_Q$par[1:3],nvec=data$n,
                              dvec=data$x,
                              yvec=data$y,
                              qval=data$q),
                      llfP42_Q(optP4_Q$par[1:3],
                               nvec=data$n,
                               dvec=data$x,
                               yvec=data$y,
                               qval=data$q,
                               rho = optP4_Q$par[stringr::str_detect(names(optP4_Q$par),'rho') &
                                                   !stringr::str_detect(names(optP4_Q$par),'eta')])
      )

      llQ = c(llQ, llP4Q)
      BMDL = c(BMDL, P4resQ[1])
      BMD = c(BMD, P4resQ[2])
      BMDU = c(BMDU, P4resQ[3])
    }
  }
  if(prior.weights[7] == 0){P4resQ=NULL; llQ = c(llQ,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA);
  P4covQ=rep(NA,2); P4corrQ=rep(NA,2); DRM_P4_Q=rep(NA,length(data$x)); parsP4Q <- NA
  P4outQ <- t(data.frame(
    # transformed parameters
    par.at = rep(NA,3),
    par.dt = rep(NA,3),
    par.k = rep(NA,3),
    # parameters on original scale
    par.a = rep(NA,3),
    par.b = rep(NA,3),
    par.c = rep(NA,3),
    par.d = rep(NA,3),
    # natural parameters
    BMD = rep(NA,3),
    min.resp = rep(NA,3)
  ))
  }
  #

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='quantal')[8]))
  }
  if(prior.weights[8]>0){
    # print(8)
    optL4_Q = fun_optimQ(stanmodels$mL4_Q, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optL4_Q[[3]]),TRUE,(optL4_Q[[3]]!=0)) | length(optL4_Q)!=9)){
      prior.weights[8] <- 0
      warning('difficulties fitting the Logit model; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      if(data$is_bin == 1) {
        parsL4Q <- parq_extract(optL4_Q, model_name = "L4_Q", pars = c('a', 'b', 'd', 'BMD',
                                                                       paste0('par',1:3)))
      } else {
        parsL4Q <- parq_extract(optL4_Q, model_name = "L4_Q", pars = c('a', 'b', 'd', 'BMD','rho[1]',
                                                                       paste0('par',1:3)),
                                rho = TRUE)
      }

      L4resQ <- quantile(parsL4Q$BMD, pvec, na.rm = T)*data$maxD
      if(data$is_bin == 1){
        L4resQ <- c(L4resQ, apply(parsL4Q[,c('a', 'b', 'd')], 2, median))
      } else {
        L4resQ <- c(L4resQ, apply(parsL4Q[,c('a', 'b', 'd','rho')], 2, median))
      }

      names(L4resQ) <- ifelse(rep(data$is_bin == 1, ifelse(data$is_bin == 1,6,7)),
                              c("BMDL","BMD","BMDU","a","b","d"),
                              c("BMDL","BMD","BMDU","a","b","d","rho")
      )

      if(data$is_bin == 1){
        L4outQ <- outLPQ(parsL4Q, pvec, data$maxD)
      } else {
        L4outQ <- outLPQ(parsL4Q, pvec, data$maxD, rho = TRUE)
      }

      # Covariance between b-d and between BMD-d
      L4covQ = c(cov(parsL4Q[,c("b","d")], use="na.or.complete")["b","d"],
                 cov(parsL4Q[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      L4corrQ = c(cor(parsL4Q[,c("b","d")], use="na.or.complete")["b","d"],
                  cor(parsL4Q[,c("BMD","d")], use="na.or.complete")["BMD","d"])

      DRM_L4_Q = DRM.L4_Q(optL4_Q$par[1:3], data$x, data$q)

      llL4Q <- ifelse(data$is_bin == 1,
                      llfL4_Q(optL4_Q$par[1:3],nvec=data$n,
                              dvec=data$x,
                              yvec=data$y,
                              qval=data$q),
                      llfL42_Q(optL4_Q$par[1:3],
                               nvec=data$n,
                               dvec=data$x,
                               yvec=data$y,
                               qval=data$q,
                               rho = optL4_Q$par[stringr::str_detect(names(optL4_Q$par),'rho') &
                                                   !stringr::str_detect(names(optL4_Q$par),'eta')])
      )

      llQ = c(llQ, llL4Q)
      BMDL = c(BMDL, L4resQ[1])
      BMD = c(BMD, L4resQ[2])
      BMDU = c(BMDU, L4resQ[3])
    }
  }
  if(prior.weights[8] == 0){L4resQ=NULL; llQ = c(llQ,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA);
  BMDU=c(BMDU,NA);
  L4covQ=rep(NA,2); L4corrQ=rep(NA,2); DRM_L4_Q=rep(NA,length(data$x)); parsL4Q <- NA
  L4outQ <- t(data.frame(
    # transformed parameters
    par.at = rep(NA,3),
    par.dt = rep(NA,3),
    par.k = rep(NA,3),
    # parameters on original scale
    par.a = rep(NA,3),
    par.b = rep(NA,3),
    par.c = rep(NA,3),
    par.d = rep(NA,3),
    # natural parameters
    BMD = rep(NA,3),
    min.resp = rep(NA)
  ))
  }


  minll = min(llQ, na.rm=T)

  # the weights
  fun.w.bin <- function(DIH, ll, min.ll, opt, mu, sig, lb, ub, s1, s2, td){

    w1 <- (2*pi)^(2.5)*sqrt(DIH)*exp(ll-min.ll)*
      # d
      truncnorm::dtruncnorm(opt$par[3], b = td, mean = mu[3], sd = sig[3,3])*
      # a
      mc2d::dpert(opt$par[1], min = lb[1], max = ub[1],
                  mode = mu[1], shape = s1)*
      # BMD
      mc2d::dpert(opt$par[2], min = lb[2], max = ub[2],
                  mode = mu[2], shape = s2)

    return(w1)
  }

  fun.w.betabin <- function(DIH, ll, min.ll, opt, mu, sig, lb, ub, s1, s2, td){

    w1 <- (2*pi)^(2.5)*sqrt(DIH)*exp(ll-min.ll)*
      # d
      truncnorm::dtruncnorm(opt$par[3], b = td, mean = mu[3], sd = sig[3,3])*
      # a
      mc2d::dpert(opt$par[1], min = lb[1], max = ub[1],
                  mode = mu[1], shape = s1)*
      # BMD
      mc2d::dpert(opt$par[2], min = lb[2], max = ub[2],
                  mode = mu[2], shape = s2)*
      mc2d::dpert(opt$par[stringr::str_detect(names(opt$par),'rho')],
                  min = 0, max = 1, mode = mu[4], shape = 4)

    return(w1)
  }

  w=c()

  lls <- c(llQ)
  w.msg <- ''

  max.ll = max(lls, na.rm = T)
  if(is.na(lls[which((max.ll-lls[!is.na(lls)]) < 709 & prior.weights>0)][1])){
    lpw <- rep(0, 8)
    lpw[which(lls == max.ll)] <- 1
    w.msg <- 'Laplace weights could not be computed and one model gets all the weight; using another prior for parameter d might help'
    warning('Laplace weights could not be computed and one model gets all the weight; using another prior for parameter d might help')
  }else{

    if(FALSE %in% ((max.ll-lls[!is.na(lls)]) < 709)){
      w.msg <- 'Not all models were used in computation of Laplace weights, some models set to 0; using another prior for parameter d might help'
      warning('Not all models were used in computation of Laplace weights, some models set to 0; using another prior for parameter d might help')
    }

    minll <- min(lls[which((max.ll-lls[!is.na(lls)]) < 709 & prior.weights>0)], na.rm = T)

    if(prior.weights[1]>0){
      DIHE4h=try(det(-pracma::pinv(optE4_Q$hessian)))
      DIHE4=ifelse(DIHE4h<0 | class(DIHE4h)[1] == 'try-error',0,DIHE4h)
      ## Aproximation of marginal (i.e. integrated) likelihood (= 'model evidence')
      if(data$is_bin==1) {
        w = c(w, fun.w.bin(DIHE4, llE4Q, minll, optE4_Q, data$priormu, data$priorSigma,
                           data$priorlb, data$priorub, data$priorgama[1], data$priorgama[2],
                           data$truncd))
      } else {
        w = c(w, fun.w.betabin(DIHE4, llE4Q, minll, optE4_Q, data$priormu, data$priorSigma,
                               data$priorlb, data$priorub, data$priorgama[1], data$priorgama[2],
                               data$truncd))

      }

    }else{w=c(w,0)}

    if(prior.weights[2]>0){
      DIHIE4h=try(det(-pracma::pinv(optIE4_Q$hessian)))
      DIHIE4=ifelse(DIHIE4h<0 | class(DIHIE4h)[1] == 'try-error',0,DIHIE4h)

      if(data$is_bin==1) {
        w = c(w, fun.w.bin(DIHIE4, llIE4Q, minll, optIE4_Q, data$priormu, data$priorSigma,
                           data$priorlb, data$priorub, data$priorgama[1], data$priorgama[2],
                           data$truncd))
      } else {
        w = c(w, fun.w.betabin(DIHIE4, llIE4Q, minll, optIE4_Q, data$priormu, data$priorSigma,
                               data$priorlb, data$priorub, data$priorgama[1], data$priorgama[2],
                               data$truncd))
      }

    }else{w=c(w,0)}

    if(prior.weights[3]>0){
      DIHH4h=try(det(-pracma::pinv(optH4_Q$hessian)))
      DIHH4=ifelse(DIHH4h<0 | class(DIHH4h)[1] == 'try-error',0,DIHH4h)

      if(data$is_bin==1) {
        w = c(w, fun.w.bin(DIHH4, llH4Q, minll, optH4_Q, data$priormu, data$priorSigma,
                           data$priorlb, data$priorub, data$priorgama[1], data$priorgama[2],
                           data$truncd))
      } else {
        w = c(w, fun.w.betabin(DIHH4, llH4Q, minll, optH4_Q, data$priormu, data$priorSigma,
                               data$priorlb, data$priorub, data$priorgama[1], data$priorgama[2],
                               data$truncd))
      }

    }else{w=c(w,0)}

    if(prior.weights[4]>0){
      DIHLN4h=try(det(-pracma::pinv(optLN4_Q$hessian)))
      DIHLN4=ifelse(DIHLN4h<0 | class(DIHLN4h)[1] == 'try-error',0,DIHLN4h)

      if(data$is_bin==1) {
        w = c(w, fun.w.bin(DIHLN4, llLN4Q, minll, optLN4_Q, data$priormu, data$priorSigma,
                           data$priorlb, data$priorub, data$priorgama[1], data$priorgama[2],
                           data$truncd))
      } else {
        w = c(w, fun.w.betabin(DIHLN4, llLN4Q, minll, optLN4_Q, data$priormu, data$priorSigma,
                               data$priorlb, data$priorub, data$priorgama[1], data$priorgama[2],
                               data$truncd))

      }

    }else{w=c(w,0)}

    if(prior.weights[5]>0){
      DIHG4h=try(det(-pracma::pinv(optG4_Q$hessian)))
      DIHG4=ifelse(DIHG4h<0 | class(DIHG4h)[1] == 'try-error',0,DIHG4h)

      if(data$is_bin==1) {
        w = c(w, fun.w.bin(DIHG4, llG4Q, minll, optG4_Q, data$priormu, data$priorSigma,
                           data$priorlb, data$priorub, data$priorgama[1], data$priorgama[2],
                           data$truncd))
      } else {
        w = c(w, fun.w.betabin(DIHG4, llG4Q, minll, optG4_Q, data$priormu, data$priorSigma,
                               data$priorlb, data$priorub, data$priorgama[1], data$priorgama[2],
                               data$truncd))
      }

    }else{w=c(w,0)}

    if(prior.weights[6]>0){
      DIHQE4h=try(det(-pracma::pinv(optQE4_Q$hessian)))
      DIHQE4=ifelse(DIHQE4h<0 | class(DIHQE4h)[1] == 'try-error',0,DIHQE4h)

      if(data$is_bin==1) {
        w = c(w, fun.w.bin(DIHQE4, llQE4Q, minll, optQE4_Q, data$priormuQ, data$priorSigmaQ,
                           data$priorlb, data$priorub, data$priorgama[1], data$priorgama[2],
                           data$truncdQ))
      } else {
        w = c(w, fun.w.betabin(DIHQE4, llQE4Q, minll, optQE4_Q, data$priormuQ, data$priorSigmaQ,
                               data$priorlb, data$priorub, data$priorgama[1], data$priorgama[2],
                               data$truncdQ))
      }

    }else{w=c(w,0)}

    if(prior.weights[7]>0){
      DIHP4h=try(det(-pracma::pinv(optP4_Q$hessian)))
      DIHP4=ifelse(DIHP4h<0 | class(DIHP4h)[1] == 'try-error',0,DIHP4h)

      if(data$is_bin==1) {
        w = c(w, fun.w.bin(DIHP4, llP4Q, minll, optP4_Q, data$priormu, data$priorSigma,
                           data$priorlb, data$priorub, data$priorgama[1], data$priorgama[2],
                           data$truncd))
      } else {
        w = c(w, fun.w.betabin(DIHP4, llP4Q, minll, optP4_Q, data$priormu, data$priorSigma,
                               data$priorlb, data$priorub, data$priorgama[1], data$priorgama[2],
                               data$truncd))
      }

    }else{w=c(w,0)}

    if(prior.weights[8]>0){
      DIHL4h=try(det(-pracma::pinv(optL4_Q$hessian)))
      DIHL4=ifelse(DIHL4h<0 | class(DIHL4h)[1] == 'try-error',0,DIHL4h)

      if(data$is_bin==1) {
        w = c(w, fun.w.bin(DIHL4, llL4Q, minll, optL4_Q, data$priormu, data$priorSigma,
                           data$priorlb, data$priorub, data$priorgama[1], data$priorgama[2],
                           data$truncd))
      } else {
        w = c(w, fun.w.betabin(DIHL4, llL4Q, minll, optL4_Q, data$priormu, data$priorSigma,
                               data$priorlb, data$priorub, data$priorgama[1], data$priorgama[2],
                               data$truncd))
      }

    }else{w=c(w,0)}

    w <- ifelse(w == 'Inf' | is.na(w), 0, w)
    prior.weights = prior.weights/sum(prior.weights==1) # in case this has changed for G4
    lpw=(prior.weights*w)/sum(prior.weights*w)

  }

  # the model average posterior as a mixture
  count=round(lpw*ndraws)
  mabmd=(c(# normal
    if(prior.weights[1]>0) sample(optE4_Q$theta_tilde[,2],count[1],replace=T),
    if(prior.weights[2]>0) sample(optIE4_Q$theta_tilde[,2],count[2],replace=T),
    if(prior.weights[3]>0) sample(optH4_Q$theta_tilde[,2],count[3],replace=T),
    if(prior.weights[4]>0) sample(optLN4_Q$theta_tilde[,2],count[4],replace=T),
    if(prior.weights[5]>0) sample(optG4_Q$theta_tilde[,2],count[5],replace=T),
    if(prior.weights[6]>0) sample(optQE4_Q$theta_tilde[,2],count[6],replace=T),
    if(prior.weights[7]>0) sample(optP4_Q$theta_tilde[,2],count[7],replace=T),
    if(prior.weights[8]>0) sample(optL4_Q$theta_tilde[,2],count[8],replace=T)
  ))

  maci=quantile(mabmd,pvec, na.rm = T)*data$maxD ## original scale
  names(maci)=c("BMDL","BMD","BMDU")

  if(TRUE %in% (mabmd > data$maxD) && data$maxD > 1){
    mabmd = ifelse(mabmd > data$maxD, data$maxD, mabmd)
    p.msg = 'The model averaged posterior distribution has been truncated at max(Dose)^2'
    warnings('The model averaged posterior distribution has been truncated at max(Dose)^2')
  }else{
    p.msg = ''
  }

  BMDq = quantile(mabmd, seq(0,1,0.005), na.rm = T)*data$maxD ## original scale

  BMDL = c(BMDL, maci[1]/data$maxD); BMD = c(BMD, maci[2]/data$maxD);
  BMDU = c(BMDU, maci[3]/data$maxD)

  names(BMDL) <- c("E4_Q","IE4_Q","H4_Q","LN4_Q","G4_Q","QE4_Q","P4_Q","L4_Q","MA")
  model = c("E4_Q","IE4_Q","H4_Q","LN4_Q","G4_Q","QE4_Q","P4_Q","L4_Q","MA")
  model = as.factor(model)
  weight = c(lpw[1], lpw[2], lpw[3], lpw[4], lpw[5], lpw[6], lpw[7], lpw[8], 1)


  names(lpw) = model[1:8]

  ### Model-averaged response per dose level
  dr.MA <- c()
  for(i in 1:length(data$x)){
    xx <- c(DRM_E4_Q[i],DRM_IE4_Q[i],DRM_H4_Q[i],DRM_LN4_Q[i],DRM_G4_Q[i],DRM_QE4_Q[i],
            DRM_P4_Q[i], DRM_L4_Q[i])
    dr.MA[i] = weighted.mean(x = xx,
                             w = lpw,
                             na.rm = T)
  }

  modelnames = c("E4","IE4","H4","LN4","G4","QE4","P4","L4")

  ### Some additional checks

  if(maci[2]/maci[1] > 20){
    warning('BMD/BMDL is larger than 20')
  }
  if(maci[3]/maci[1] > 50){
    warning('BMDU/BMDL is larger than 50')
  }
  if(maci[2] < (data.Q$data$x[2]*data.Q$data$maxD/10)){
    warning('BMD is 10 times lower than the lowest non-zero dose')
  }

  ### best fitting model vs saturated ANOVA model
  best.fit = modelnames[which(weight[1:8] == max(weight[1:8]))][1]
  nrchains = 3; nriterations = 3000; warmup = 1000; delta = 0.8; treedepth = 10
  bfTest <- modelTestQ(best.fit, data.Q, get(paste0('opt', best.fit, '_Q')), type = 'Laplace',
                       seed, ndraws, nrchains, nriterations, warmup, delta, treedepth)
  #warning(bfTest$warn.bf)

  ## Covariances
  covs = t(data.frame(
    E4_Q = E4covQ,
    IE4_Q = IE4covQ,
    H4_Q = H4covQ,
    LN4_Q = LN4covQ,
    G4_Q = G4covQ,
    QE4_Q = QE4covQ,
    P4_Q = P4covQ,
    L4_Q = L4covQ
  ))
  colnames(covs) = c("b-d", "BMD-d")

  corrs = t(data.frame(
    E4_Q = E4corrQ,
    IE4_Q = IE4corrQ,
    H4_Q = H4corrQ,
    LN4_Q = LN4corrQ,
    G4_Q = G4corrQ,
    QE4_Q = QE4corrQ,
    P4_Q = P4corrQ,
    L4_Q = L4corrQ
  ))
  colnames(corrs) = c("b-d", "BMD-d")

  modelnames2 = c("E4_Q","IE4_Q","H4_Q","LN4_Q","G4_Q","QE4_Q","P4_Q","L4_Q")


  ###
  # sd = sqrt(data$s2)
  # N = data$N
  # gmean.a = log(NtoLN(data$m,sd))[1:N]
  # gsd.a = log(NtoLN(data$m,sd))[(N+1):(2*N)]

  ret_results <- list(
    # model parameters
    E4_Q=E4outQ,IE4_Q=IE4outQ,H4_Q=H4outQ,LN4_Q=LN4outQ,G4_Q=G4outQ,QE4_Q=QE4outQ,P4_Q=P4outQ,L4_Q=L4outQ,

    # covariances
    covs = covs, corrs = corrs,
    # weights and MA
    weights=lpw,
    MA=maci,
    llQ=llQ, MA_post = BMDq,
    # dose-response
    MA_dr = dr.MA,
    parsQ = list(parsE4Q, parsIE4Q, parsH4Q, parsLN4Q, parsG4Q,
                 parsQE4Q, parsP4Q, parsL4Q),
    BMDMixture = mabmd*data.Q$data$maxD,
    data = data.frame(
      dose = data.Q$data$x,
      y = data.Q$data$y,
      n = data.Q$data$n),
    max.dose = data.Q$data$maxD,
    models_included_laplace = modelnames2[prior.weights > 0],
    q = data.Q$data$q,
    is_bin = data.Q$data$is_bin,
    is_betabin = data.Q$data$is_betabin,
    bf = bfTest$bayesFactor, gof_check = bfTest$warn.bf,
    # means.SM = bfTest$means.SM, parBestFit = bfTest$par.best,
    # BIC.bestfit = bfTest$BIC.bestfit, BIC.SM = bfTest$BIC.SM,
    w.msg = w.msg, p.msg = p.msg
  )

  attr(ret_results, "class") <- c("BMADRQ", "LP")

  return(ret_results)
}


#' Perform model averaging using Full Laplace method for data with covariate effect
#'
#' @param data the summary data, with columns: dose, response, sd, n, covariate level for continuous; dose, number of adverse events, n, covariate for quantal
#' @param sumstats logical indicating whether summary (T, default) or individual-level (F) data is provided
#' @param sd logical indicating whether standard deviation (T, default) or standard error (F) is provided
#' @param q specified BMR
#' @param prior.d prior distribution for parameter d (on log scale), should be either N11 (default N(1, 1) prior truncated at 5), EPA (N(0.4, sqrt(0.5)) prior) or N05 (for a N(0.5,0.5) prior)
#' @param extended logical indicating whether the dose range should be extended to maxDose^2 (default is TRUE)
#' @param extended.value value for the upper range of BMD prior
#' @param prior.weights a vector specifying which of the 16 (continuous) or 8 (quantal) models should be included (1 = include, 0 = exclude)
#' @param ndraws the number of draws, default 30000
#' @param seed default 123
#' @param pvec vector specifying the three BMD quantiles of interest (default 90% CrI)
#'
#' @description This function performs model averaging using the full Laplace approximation. By default, all 16 models are included for continuous data, and all 8 models for quantal data.
#'              Models can be excluded by setting their respective weight to 0 in \code{prior.weights}. The order of models fitted can be obtained using the \code{get_models()} function.
#'
#'  `full.laplace_MA_Cov` is used for continuous data
#'
#'  `full.laplace_MA_Q_Cov` is used for quantal data
#'
#' @examples
#' summ.data = data.frame(x = data_cont_covariate$Dose, y = data_cont_covariate$Response, s = data_cont_covariate$SD, n = data_cont_covariate$N, cov = data_cont_covariate$Covariate)
#' FLBMD <- full.laplace_MA_Cov(summ.data, sumstats = T, sd = T, q = 0.05, prior.d = 'N11', extended = T, ndraws = 30000, seed = 123, pvec = c(0.05, 0.5, 0.95), prior.weights = rep(1,16))
#'
#' @importFrom truncnorm dtruncnorm
#'
#' @return `MA` Model averaged BMD credible interval per covariate level (in case of an effect)
#' @return `weights` Weights for the best fitting submodels used in model averaging
#' @return `summary` Summary table of results
#' @return `data` Data used in analysis
#' @return `maxDose` Maximum dose level (original scale)
#' @return `BMD_Mixture` 0.5%-percentiles of the model averaged posterior for the BMD
#' @return `parE4_N` Estimated parameters for the (Exponential, Normal) model
#' @return `parIE4_N` Estimated parameters for the (Inverse Exponential, Normal) model
#' @return `parH4_N` Estimated parameters for the (Hill, Normal) model
#' @return `parLN4_N` Estimated parameters for the (Lognormal, Normal) model
#' @return `parG4_N` Estimated parameters for the (Gamma, Normal) model
#' @return `parQE4_N` Estimated parameters for the (Quadratic Exponential, Normal) model
#' @return `parP4_N` Estimated parameters for the (Probit, Normal) model
#' @return `parL4_N` Estimated parameters for the (Logit, Normal) model
#' @return `parE4_LN` Estimated parameters for the (Exponential, Lognormal) model
#' @return `parIE4_LN` Estimated parameters for the (Inverse Exponential, Lognormal) model
#' @return `parH4_LN` Estimated parameters for the (Hill, Lognormal) model
#' @return `parLN4_LN` Estimated parameters for the (Lognormal, Lognormal) model
#' @return `parG4_LN` Estimated parameters for the (Gamma, Lognormal) model
#' @return `parQE4_LN` Estimated parameters for the (Quadratic Exponential, Lognormal) model
#' @return `parP4_LN` Estimated parameters for the (Probit, Lognormal) model
#' @return `parL4_LN` Estimated parameters for the (Logit, Lognormal) model
#' @return `parE4_Q` Estimated parameters for the (Exponential, Quantal) model
#' @return `parIE4_Q` Estimated parameters for the (Inverse Exponential, Quantal) model
#' @return `parH4_Q` Estimated parameters for the (Hill, Quantal) model
#' @return `parLN4_Q` Estimated parameters for the (Lognormal, Quantal) model
#' @return `parG4_Q` Estimated parameters for the (Gamma, Quantal) model
#' @return `parQE4_Q` Estimated parameters for the (Quadratic Exponential, Quantal) model
#' @return `parP4_Q` Estimated parameters for the (Probit, Quantal) model
#' @return `parL4_Q` Estimated parameters for the (Logit, Quantal) model
#' @return `shift` shift in case of negative geometric means
#' @return `q` BMR
#'
#' @export full.laplace_MA_Cov
#'
full.laplace_MA_Cov = function(data, # the summary data
                               sumstats = TRUE,
                               sd = TRUE,
                               q = 0.05,
                               prior.d = 'N11',
                               extended = TRUE, extended.value = 3,
                               prior.weights = rep(1,16),
                               ndraws=30000,seed=123,
                               pvec=c(0.05,0.5,0.95)
){
  #### Set input options (priors, datatype, etc) !!!!

  ##############################################
  ### DATA IN CORRECT FORMAT FOR EACH SUBMODEL

  out.stop = 'ok'

  data_NCOV_all <- PREP_DATA_NCOV(
    data = data,
    sumstats = sumstats,
    sd = sd,
    q = q,
    prior.d = prior.d,
    extended = extended,
    extended.value = extended.value,
    covariate = 'all'
  )

  data_N_noCOV <- PREP_DATA_N(data = data,
                              sumstats = sumstats,
                              sd = sd,
                              q = q,
                              prior.d = prior.d,
                              extended = extended)

  if(1 %in% prior.weights[1:8]){

    data_NCOV_asigma2 <- PREP_DATA_NCOV(
      data = data,
      sumstats = sumstats,
      sd = sd,
      q = q,
      prior.d = prior.d,
      extended = extended,
      extended.value = extended.value,
      covariate = 'a_sigma2'
    )

    data_NCOV_dBMD <- PREP_DATA_NCOV(
      data = data,
      sumstats = sumstats,
      sd = sd,
      q = q,
      prior.d = prior.d,
      extended = extended,
      extended.value = extended.value,
      covariate = 'BMD_d'
    )

  }

  if(1 %in% prior.weights[9:16]){
    data_LNCOV_all <- PREP_DATA_LNCOV(
      data = data,
      sumstats = sumstats,
      sd = sd,
      q = q,
      prior.d = prior.d,
      extended = extended,
      extended.value = extended.value,
      covariate = 'all'
    )

    data_LNCOV_asigma2 <- PREP_DATA_LNCOV(
      data = data,
      sumstats = sumstats,
      sd = sd,
      q = q,
      prior.d = prior.d,
      extended = extended,
      extended.value = extended.value,
      covariate = 'a_sigma2'
    )

    data_LNCOV_dBMD <- PREP_DATA_LNCOV(
      data = data,
      sumstats = sumstats,
      sd = sd,
      q = q,
      prior.d = prior.d,
      extended = extended,
      extended.value = extended.value,
      covariate = 'BMD_d'
    )

    data_LN_noCOV <- PREP_DATA_LN(data = data,
                                  sumstats = sumstats,
                                  sd = sd,
                                  q = q,
                                  prior.d = prior.d,
                                  extended = extended,
                                  extended.value = extended.value)
  }

  #########################################
  ### SELECT BEST SUBMODEL FOR EACH DRM

  best.sub.which <- c()
  best.sub.loglik <- c()
  best.sub.weight <- c()

  nModels = 16

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[1]))
  }
  if(prior.weights[1] > 0){
    # print(1)
    select_E4_N <- fun_cov_selection(model = stanmodels$mE4COV,
                                     model_name = 'E4_N',
                                     model.none = stanmodels$mE4,
                                     loglik = ifelse(data_NCOV_all$data$data_type == 1, llfE4_NI_Cov, llfE4_ND_Cov),
                                     data_asigma2 = data_NCOV_asigma2,
                                     data_dBMD = data_NCOV_dBMD,
                                     data_all = data_NCOV_all,
                                     data_none = data_N_noCOV,
                                     prior.weightsCov = rep(1, 4), # weights for the 4 sub-models
                                     pvec = pvec,
                                     ndraws = ndraws, td = data_NCOV_all$data$truncd, seed = seed)
    if(exists('select_E4_N') & !is.null(select_E4_N)){
      best.sub.which <- c(best.sub.which, select_E4_N$best.submodel.name)
      best.sub.loglik <- c(best.sub.loglik, select_E4_N$Weights$Loglik[select_E4_N$Weights$Covariate == select_E4_N$best.submodel.name])
      best.sub.weight <- c(best.sub.weight, select_E4_N$Weights$Weights[select_E4_N$Weights$Covariate == select_E4_N$best.submodel.name])
      BMD_E4_N <- getBMD(select_E4_N$best.submodel$theta_tilde, pvec, select_E4_N$best.submodel.name, 'E4_N')
    }else{
      best.sub.which <- c(best.sub.which, NA)
      best.sub.loglik <- c(best.sub.loglik, NA)
      BMD_E4_N <- c(NA,NA,NA)
      best.sub.weight <- c(best.sub.weight, NA)
    }
  }else{
    best.sub.which <- c(best.sub.which, NA)
    best.sub.loglik <- c(best.sub.loglik, NA)
    BMD_E4_N <- c(NA,NA,NA)
    best.sub.weight <- c(best.sub.weight, NA)
  }

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[2]))
  }
  if(prior.weights[2] > 0){
    # print(2)
    select_IE4_N <- fun_cov_selection(model = stanmodels$mIE4COV,
                                      model_name = 'IE4_N',
                                      model.none = stanmodels$mIE4,
                                      loglik = ifelse(data_NCOV_all$data$data_type == 1, llfIE4_NI_Cov, llfIE4_ND_Cov),
                                      data_asigma2 = data_NCOV_asigma2,
                                      data_dBMD = data_NCOV_dBMD,
                                      data_all = data_NCOV_all,
                                      data_none = data_N_noCOV,
                                      prior.weightsCov = rep(1, 4), # weights for the 4 sub-models
                                      pvec = pvec,
                                      ndraws = ndraws, td = data_NCOV_all$data$truncd, seed = seed)
    if(exists('select_IE4_N') & !is.null(select_IE4_N)){
      best.sub.which <- c(best.sub.which, select_IE4_N$best.submodel.name)
      best.sub.loglik <- c(best.sub.loglik, select_IE4_N$Weights$Loglik[select_IE4_N$Weights$Covariate == select_IE4_N$best.submodel.name])
      best.sub.weight <- c(best.sub.weight, select_IE4_N$Weights$Weights[select_IE4_N$Weights$Covariate == select_IE4_N$best.submodel.name])
      BMD_IE4_N <- getBMD(select_IE4_N$best.submodel$theta_tilde, pvec, select_IE4_N$best.submodel.name, 'IE4_N')
    }else{
      best.sub.which <- c(best.sub.which, NA)
      best.sub.loglik <- c(best.sub.loglik, NA)
      BMD_IE4_N <- c(NA,NA,NA)
      best.sub.weight <- c(best.sub.weight, NA)
    }
  }else{
    best.sub.which <- c(best.sub.which, NA)
    best.sub.loglik <- c(best.sub.loglik, NA)
    BMD_IE4_N <- c(NA,NA,NA)
    best.sub.weight <- c(best.sub.weight, NA)
  }

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[3]))
  }
  if(prior.weights[3] > 0){
    # print(3)
    select_H4_N <- fun_cov_selection(model = stanmodels$mH4COV,
                                     model_name = 'H4_N',
                                     model.none = stanmodels$mH4,
                                     loglik = ifelse(data_NCOV_all$data$data_type == 1, llfH4_NI_Cov, llfH4_ND_Cov),
                                     data_asigma2 = data_NCOV_asigma2,
                                     data_dBMD = data_NCOV_dBMD,
                                     data_all = data_NCOV_all,
                                     data_none = data_N_noCOV,
                                     prior.weightsCov = rep(1, 4), # weights for the 4 sub-models
                                     pvec = pvec,
                                     ndraws = ndraws, td = data_NCOV_all$data$truncd, seed = seed)
    if(exists('select_H4_N') & !is.null(select_H4_N)){
      best.sub.which <- c(best.sub.which, select_H4_N$best.submodel.name)
      best.sub.loglik <- c(best.sub.loglik, select_H4_N$Weights$Loglik[select_H4_N$Weights$Covariate == select_H4_N$best.submodel.name])
      best.sub.weight <- c(best.sub.weight, select_H4_N$Weights$Weights[select_H4_N$Weights$Covariate == select_H4_N$best.submodel.name])
      BMD_H4_N <- getBMD(select_H4_N$best.submodel$theta_tilde, pvec, select_H4_N$best.submodel.name, 'H4_N')
    }else{
      best.sub.which <- c(best.sub.which, NA)
      best.sub.loglik <- c(best.sub.loglik, NA)
      BMD_H4_N <- c(NA,NA,NA)
      best.sub.weight <- c(best.sub.weight, NA)
    }
  }else{
    best.sub.which <- c(best.sub.which, NA)
    best.sub.loglik <- c(best.sub.loglik, NA)
    BMD_H4_N <- c(NA,NA,NA)
    best.sub.weight <- c(best.sub.weight, NA)
  }

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[4]))
  }
  if(prior.weights[4] > 0){
    # print(4)
    select_LN4_N <- fun_cov_selection(model = stanmodels$mLN4COV,
                                      model_name = 'LN4_N',
                                      model.none = stanmodels$mLN4,
                                      loglik = ifelse(data_NCOV_all$data$data_type == 1, llfLN4_NI_Cov, llfLN4_ND_Cov),
                                      data_asigma2 = data_NCOV_asigma2,
                                      data_dBMD = data_NCOV_dBMD,
                                      data_all = data_NCOV_all,
                                      data_none = data_N_noCOV,
                                      prior.weightsCov = rep(1, 4), # weights for the 4 sub-models
                                      pvec = pvec,
                                      ndraws = ndraws, td = data_NCOV_all$data$truncd, seed = seed)
    if(exists('select_LN4_N') & !is.null(select_LN4_N)){
      best.sub.which <- c(best.sub.which, select_LN4_N$best.submodel.name)
      best.sub.loglik <- c(best.sub.loglik, select_LN4_N$Weights$Loglik[select_LN4_N$Weights$Covariate == select_LN4_N$best.submodel.name])
      best.sub.weight <- c(best.sub.weight, select_LN4_N$Weights$Weights[select_LN4_N$Weights$Covariate == select_LN4_N$best.submodel.name])
      BMD_LN4_N <- getBMD(select_LN4_N$best.submodel$theta_tilde, pvec, select_LN4_N$best.submodel.name, 'LN4_N')
    }else{
      best.sub.which <- c(best.sub.which, NA)
      best.sub.loglik <- c(best.sub.loglik, NA)
      BMD_LN4_N <- c(NA,NA,NA)
      best.sub.weight <- c(best.sub.weight, NA)
    }
  }else{
    best.sub.which <- c(best.sub.which, NA)
    best.sub.loglik <- c(best.sub.loglik, NA)
    BMD_LN4_N <- c(NA,NA,NA)
    best.sub.weight <- c(best.sub.weight, NA)
  }

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[5]))
  }
  if(prior.weights[5] > 0){
    # print(5)
    select_G4_N <- fun_cov_selection(model = stanmodels$mG4COV,
                                     model_name = 'G4_N',
                                     model.none = stanmodels$mG4,
                                     loglik = ifelse(data_NCOV_all$data$data_type == 1, llfG4_NI_Cov, llfG4_ND_Cov),
                                     data_asigma2 = data_NCOV_asigma2,
                                     data_dBMD = data_NCOV_dBMD,
                                     data_all = data_NCOV_all,
                                     data_none = data_N_noCOV,
                                     prior.weightsCov = rep(1, 4), # weights for the 4 sub-models
                                     pvec = pvec,
                                     ndraws = ndraws, td = data_NCOV_all$data$truncd, seed = seed)
    if(exists('select_G4_N') & !is.null(select_G4_N)){
      best.sub.which <- c(best.sub.which, select_G4_N$best.submodel.name)
      best.sub.loglik <- c(best.sub.loglik, select_G4_N$Weights$Loglik[select_G4_N$Weights$Covariate == select_G4_N$best.submodel.name])
      best.sub.weight <- c(best.sub.weight, select_G4_N$Weights$Weights[select_G4_N$Weights$Covariate == select_G4_N$best.submodel.name])
      BMD_G4_N <- getBMD(select_G4_N$best.submodel$theta_tilde, pvec, select_G4_N$best.submodel.name, 'G4_N')
    }else{
      best.sub.which <- c(best.sub.which, NA)
      best.sub.loglik <- c(best.sub.loglik, NA)
      BMD_G4_N <- c(NA,NA,NA)
      best.sub.weight <- c(best.sub.weight, NA)
    }
  }else{
    best.sub.which <- c(best.sub.which, NA)
    best.sub.loglik <- c(best.sub.loglik, NA)
    BMD_G4_N <- c(NA,NA,NA)
    best.sub.weight <- c(best.sub.weight, NA)
  }

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[6]))
  }
  if(prior.weights[6] > 0){
    # print(6)
    select_QE4_N <- fun_cov_selection(model = stanmodels$mQE4COV,
                                      model_name = 'QE4_N',
                                      model.none = stanmodels$mQE4,
                                      loglik = ifelse(data_NCOV_all$data$data_type == 1, llfQE4_NI_Cov, llfQE4_ND_Cov),
                                      data_asigma2 = data_NCOV_asigma2,
                                      data_dBMD = data_NCOV_dBMD,
                                      data_all = data_NCOV_all,
                                      data_none = data_N_noCOV,
                                      prior.weightsCov = rep(1, 4), # weights for the 4 sub-models
                                      pvec = pvec,
                                      ndraws = ndraws, td = data_NCOV_all$data$truncdQ, seed = seed)
    if(exists('select_QE4_N') & !is.null(select_QE4_N)){
      best.sub.which <- c(best.sub.which, select_QE4_N$best.submodel.name)
      best.sub.loglik <- c(best.sub.loglik, select_QE4_N$Weights$Loglik[select_QE4_N$Weights$Covariate == select_QE4_N$best.submodel.name])
      best.sub.weight <- c(best.sub.weight, select_QE4_N$Weights$Weights[select_QE4_N$Weights$Covariate == select_QE4_N$best.submodel.name])
      BMD_QE4_N <- getBMD(select_QE4_N$best.submodel$theta_tilde, pvec, select_QE4_N$best.submodel.name, 'QE4_N')
    }else{
      best.sub.which <- c(best.sub.which, NA)
      best.sub.loglik <- c(best.sub.loglik, NA)
      BMD_QE4_N <- c(NA,NA,NA)
      best.sub.weight <- c(best.sub.weight, NA)
    }
  }else{
    best.sub.which <- c(best.sub.which, NA)
    best.sub.loglik <- c(best.sub.loglik, NA)
    BMD_QE4_N <- c(NA,NA,NA)
    best.sub.weight <- c(best.sub.weight, NA)
  }

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[7]))
  }
  if(prior.weights[7] > 0){
    # print(7)
    select_P4_N <- fun_cov_selection(model = stanmodels$mP4COV,
                                     model_name = 'P4_N',
                                     model.none = stanmodels$mP4,
                                     loglik = ifelse(data_NCOV_all$data$data_type == 1, llfP4_NI_Cov, llfP4_ND_Cov),
                                     data_asigma2 = data_NCOV_asigma2,
                                     data_dBMD = data_NCOV_dBMD,
                                     data_all = data_NCOV_all,
                                     data_none = data_N_noCOV,
                                     prior.weightsCov = rep(1, 4), # weights for the 4 sub-models
                                     pvec = pvec,
                                     ndraws = ndraws, td = data_NCOV_all$data$truncd, seed = seed)
    if(exists('select_P4_N') & !is.null(select_P4_N)){
      best.sub.which <- c(best.sub.which, select_P4_N$best.submodel.name)
      best.sub.loglik <- c(best.sub.loglik, select_P4_N$Weights$Loglik[select_P4_N$Weights$Covariate == select_P4_N$best.submodel.name])
      best.sub.weight <- c(best.sub.weight, select_P4_N$Weights$Weights[select_P4_N$Weights$Covariate == select_P4_N$best.submodel.name])
      BMD_P4_N <- getBMD(select_P4_N$best.submodel$theta_tilde, pvec, select_P4_N$best.submodel.name, 'P4_N')
    }else{
      best.sub.which <- c(best.sub.which, NA)
      best.sub.loglik <- c(best.sub.loglik, NA)
      BMD_P4_N <- c(NA,NA,NA)
      best.sub.weight <- c(best.sub.weight, NA)
    }
  }else{
    best.sub.which <- c(best.sub.which, NA)
    best.sub.loglik <- c(best.sub.loglik, NA)
    BMD_P4_N <- c(NA,NA,NA)
    best.sub.weight <- c(best.sub.weight, NA)
  }

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[8]))
  }
  if(prior.weights[8] > 0){
    # print(8)
    select_L4_N <- fun_cov_selection(model = stanmodels$mL4COV,
                                     model_name = 'L4_N',
                                     model.none = stanmodels$mL4,
                                     loglik = ifelse(data_NCOV_all$data$data_type == 1, llfL4_NI_Cov, llfL4_ND_Cov),
                                     data_asigma2 = data_NCOV_asigma2,
                                     data_dBMD = data_NCOV_dBMD,
                                     data_all = data_NCOV_all,
                                     data_none = data_N_noCOV,
                                     prior.weightsCov = rep(1, 4), # weights for the 4 sub-models
                                     pvec = pvec,
                                     ndraws = ndraws, td = data_NCOV_all$data$truncd, seed = seed)
    if(exists('select_L4_N') & !is.null(select_L4_N)){
      best.sub.which <- c(best.sub.which, select_L4_N$best.submodel.name)
      best.sub.loglik <- c(best.sub.loglik, select_L4_N$Weights$Loglik[select_L4_N$Weights$Covariate == select_L4_N$best.submodel.name])
      best.sub.weight <- c(best.sub.weight, select_L4_N$Weights$Weights[select_L4_N$Weights$Covariate == select_L4_N$best.submodel.name])
      BMD_L4_N <- getBMD(select_L4_N$best.submodel$theta_tilde, pvec, select_L4_N$best.submodel.name, 'L4_N')
    }else{
      best.sub.which <- c(best.sub.which, NA)
      best.sub.loglik <- c(best.sub.loglik, NA)
      BMD_L4_N <- c(NA,NA,NA)
      best.sub.weight <- c(best.sub.weight, NA)
    }
  }else{
    best.sub.which <- c(best.sub.which, NA)
    best.sub.loglik <- c(best.sub.loglik, NA)
    BMD_L4_N <- c(NA,NA,NA)
    best.sub.weight <- c(best.sub.weight, NA)
  }

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[9]))
  }
  if(prior.weights[9] > 0){
    # print(9)
    select_E4_LN <- fun_cov_selection(model = stanmodels$mE4COV,
                                      model_name = 'E4_LN',
                                      model.none = stanmodels$mE4,
                                      loglik = ifelse(data_LNCOV_all$data$data_type == 2, llfE4_LNI_Cov, llfE4_LND_Cov),
                                      data_asigma2 = data_LNCOV_asigma2,
                                      data_dBMD = data_LNCOV_dBMD,
                                      data_all = data_LNCOV_all,
                                      data_none = data_LN_noCOV,
                                      prior.weightsCov = rep(1, 4), # weights for the 4 sub-models
                                      pvec = pvec,
                                      ndraws = ndraws, td = data_LNCOV_all$data$truncd, seed = seed)
    if(exists('select_E4_LN') & !is.null(select_E4_LN)){
      best.sub.which <- c(best.sub.which, select_E4_LN$best.submodel.name)
      best.sub.loglik <- c(best.sub.loglik, select_E4_LN$Weights$Loglik[select_E4_LN$Weights$Covariate == select_E4_LN$best.submodel.name])
      best.sub.weight <- c(best.sub.weight, select_E4_LN$Weights$Weights[select_E4_LN$Weights$Covariate == select_E4_LN$best.submodel.name])
      BMD_E4_LN <- getBMD(select_E4_LN$best.submodel$theta_tilde, pvec, select_E4_LN$best.submodel.name, 'E4_LN')
    }else{
      best.sub.which <- c(best.sub.which, NA)
      best.sub.loglik <- c(best.sub.loglik, NA)
      BMD_E4_LN <- c(NA,NA,NA)
      best.sub.weight <- c(best.sub.weight, NA)
    }
  }else{
    best.sub.which <- c(best.sub.which, NA)
    best.sub.loglik <- c(best.sub.loglik, NA)
    BMD_E4_LN <- c(NA,NA,NA)
    best.sub.weight <- c(best.sub.weight, NA)
  }

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[10]))
  }
  if(prior.weights[10] > 0){
    # print(10)
    select_IE4_LN <- fun_cov_selection(model = stanmodels$mIE4COV,
                                       model_name = 'IE4_LN',
                                       model.none = stanmodels$mIE4,
                                       loglik = ifelse(data_LNCOV_all$data$data_type == 2, llfIE4_LNI_Cov, llfIE4_LND_Cov),
                                       data_asigma2 = data_LNCOV_asigma2,
                                       data_dBMD = data_LNCOV_dBMD,
                                       data_all = data_LNCOV_all,
                                       data_none = data_LN_noCOV,
                                       prior.weightsCov = rep(1, 4), # weights for the 4 sub-models
                                       pvec = pvec,
                                       ndraws = ndraws, td = data_LNCOV_all$data$truncd, seed = seed)
    if(exists('select_IE4_LN') & !is.null(select_IE4_LN)){
      best.sub.which <- c(best.sub.which, select_IE4_LN$best.submodel.name)
      best.sub.loglik <- c(best.sub.loglik, select_IE4_LN$Weights$Loglik[select_IE4_LN$Weights$Covariate == select_IE4_LN$best.submodel.name])
      best.sub.weight <- c(best.sub.weight, select_IE4_LN$Weights$Weights[select_IE4_LN$Weights$Covariate == select_IE4_LN$best.submodel.name])
      BMD_IE4_LN <- getBMD(select_IE4_LN$best.submodel$theta_tilde, pvec, select_IE4_LN$best.submodel.name, 'IE4_LN')
    }else{
      best.sub.which <- c(best.sub.which, NA)
      best.sub.loglik <- c(best.sub.loglik, NA)
      BMD_IE4_LN <- c(NA,NA,NA)
      best.sub.weight <- c(best.sub.weight, NA)
    }
  }else{
    best.sub.which <- c(best.sub.which, NA)
    best.sub.loglik <- c(best.sub.loglik, NA)
    BMD_IE4_LN <- c(NA,NA,NA)
    best.sub.weight <- c(best.sub.weight, NA)
  }

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[11]))
  }
  if(prior.weights[11] > 0){
    # print(11)
    select_H4_LN <- fun_cov_selection(model = stanmodels$mH4COV,
                                      model_name = 'H4_LN',
                                      model.none = stanmodels$mH4,
                                      loglik = ifelse(data_LNCOV_all$data$data_type == 2, llfH4_LNI_Cov, llfH4_LND_Cov),
                                      data_asigma2 = data_LNCOV_asigma2,
                                      data_dBMD = data_LNCOV_dBMD,
                                      data_all = data_LNCOV_all,
                                      data_none = data_LN_noCOV,
                                      prior.weightsCov = rep(1, 4), # weights for the 4 sub-models
                                      pvec = pvec,
                                      ndraws = ndraws, td = data_LNCOV_all$data$truncd, seed = seed)
    if(exists('select_H4_LN') & !is.null(select_H4_LN)){
      best.sub.which <- c(best.sub.which, select_H4_LN$best.submodel.name)
      best.sub.loglik <- c(best.sub.loglik, select_H4_LN$Weights$Loglik[select_H4_LN$Weights$Covariate == select_H4_LN$best.submodel.name])
      best.sub.weight <- c(best.sub.weight, select_H4_LN$Weights$Weights[select_H4_LN$Weights$Covariate == select_H4_LN$best.submodel.name])
      BMD_H4_LN <- getBMD(select_H4_LN$best.submodel$theta_tilde, pvec, select_H4_LN$best.submodel.name, 'H4_LN')
    }else{
      best.sub.which <- c(best.sub.which, NA)
      best.sub.loglik <- c(best.sub.loglik, NA)
      BMD_H4_LN <- c(NA,NA,NA)
      best.sub.weight <- c(best.sub.weight, NA)
    }
  }else{
    best.sub.which <- c(best.sub.which, NA)
    best.sub.loglik <- c(best.sub.loglik, NA)
    BMD_H4_LN <- c(NA,NA,NA)
    best.sub.weight <- c(best.sub.weight, NA)
  }

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[12]))
  }
  if(prior.weights[12] > 0){
    # print(12)
    select_LN4_LN <- fun_cov_selection(model = stanmodels$mLN4COV,
                                       model_name = 'LN4_LN',
                                       model.none = stanmodels$mLN4,
                                       loglik = ifelse(data_LNCOV_all$data$data_type == 2, llfLN4_LNI_Cov, llfLN4_LND_Cov),
                                       data_asigma2 = data_LNCOV_asigma2,
                                       data_dBMD = data_LNCOV_dBMD,
                                       data_all = data_LNCOV_all,
                                       data_none = data_LN_noCOV,
                                       prior.weightsCov = rep(1, 4), # weights for the 4 sub-models
                                       pvec = pvec,
                                       ndraws = ndraws, td = data_LNCOV_all$data$truncd, seed = seed)
    if(exists('select_LN4_LN') & !is.null(select_LN4_LN)){
      best.sub.which <- c(best.sub.which, select_LN4_LN$best.submodel.name)
      best.sub.loglik <- c(best.sub.loglik, select_LN4_LN$Weights$Loglik[select_LN4_LN$Weights$Covariate == select_LN4_LN$best.submodel.name])
      best.sub.weight <- c(best.sub.weight, select_LN4_LN$Weights$Weights[select_LN4_LN$Weights$Covariate == select_LN4_LN$best.submodel.name])
      BMD_LN4_LN <- getBMD(select_LN4_LN$best.submodel$theta_tilde, pvec, select_LN4_LN$best.submodel.name, 'LN4_LN')
    }else{
      best.sub.which <- c(best.sub.which, NA)
      best.sub.loglik <- c(best.sub.loglik, NA)
      BMD_LN4_LN <- c(NA,NA,NA)
      best.sub.weight <- c(best.sub.weight, NA)
    }
  }else{
    best.sub.which <- c(best.sub.which, NA)
    best.sub.loglik <- c(best.sub.loglik, NA)
    BMD_LN4_LN <- c(NA,NA,NA)
    best.sub.weight <- c(best.sub.weight, NA)
  }

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[13]))
  }
  if(prior.weights[13] > 0){
    # print(13)
    select_G4_LN <- fun_cov_selection(model = stanmodels$mG4COV,
                                      model_name = 'G4_LN',
                                      model.none = stanmodels$mG4,
                                      loglik = ifelse(data_LNCOV_all$data$data_type == 2, llfG4_LNI_Cov, llfG4_LND_Cov),
                                      data_asigma2 = data_LNCOV_asigma2,
                                      data_dBMD = data_LNCOV_dBMD,
                                      data_all = data_LNCOV_all,
                                      data_none = data_LN_noCOV,
                                      prior.weightsCov = rep(1, 4), # weights for the 4 sub-models
                                      pvec = pvec,
                                      ndraws = ndraws, td = data_LNCOV_all$data$truncd, seed = seed)
    if(exists('select_G4_LN') & !is.null(select_G4_LN)){
      best.sub.which <- c(best.sub.which, select_G4_LN$best.submodel.name)
      best.sub.loglik <- c(best.sub.loglik, select_G4_LN$Weights$Loglik[select_G4_LN$Weights$Covariate == select_G4_LN$best.submodel.name])
      best.sub.weight <- c(best.sub.weight, select_G4_LN$Weights$Weights[select_G4_LN$Weights$Covariate == select_G4_LN$best.submodel.name])
      BMD_G4_LN <- getBMD(select_G4_LN$best.submodel$theta_tilde, pvec, select_G4_LN$best.submodel.name, 'G4_LN')
    }else{
      best.sub.which <- c(best.sub.which, NA)
      best.sub.loglik <- c(best.sub.loglik, NA)
      BMD_G4_LN <- c(NA,NA,NA)
      best.sub.weight <- c(best.sub.weight, NA)
    }
  }else{
    best.sub.which <- c(best.sub.which, NA)
    best.sub.loglik <- c(best.sub.loglik, NA)
    BMD_G4_LN <- c(NA,NA,NA)
    best.sub.weight <- c(best.sub.weight, NA)
  }

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[14]))
  }
  if(prior.weights[14] > 0){
    # print(14)
    select_QE4_LN <- fun_cov_selection(model = stanmodels$mQE4COV,
                                       model_name = 'QE4_LN',
                                       model.none = stanmodels$mQE4,
                                       loglik = ifelse(data_LNCOV_all$data$data_type == 2, llfQE4_LNI_Cov, llfQE4_LND_Cov),
                                       data_asigma2 = data_LNCOV_asigma2,
                                       data_dBMD = data_LNCOV_dBMD,
                                       data_all = data_LNCOV_all,
                                       data_none = data_LN_noCOV,
                                       prior.weightsCov = rep(1, 4), # weights for the 4 sub-models
                                       pvec = pvec,
                                       ndraws = ndraws, td = data_LNCOV_all$data$truncdQ, seed = seed)
    if(exists('select_QE4_LN') & !is.null(select_QE4_LN)){
      best.sub.which <- c(best.sub.which, select_QE4_LN$best.submodel.name)
      best.sub.loglik <- c(best.sub.loglik, select_QE4_LN$Weights$Loglik[select_QE4_LN$Weights$Covariate == select_QE4_LN$best.submodel.name])
      best.sub.weight <- c(best.sub.weight, select_QE4_LN$Weights$Weights[select_QE4_LN$Weights$Covariate == select_QE4_LN$best.submodel.name])
      BMD_QE4_LN <- getBMD(select_QE4_LN$best.submodel$theta_tilde, pvec, select_QE4_LN$best.submodel.name, 'QE4_LN')
    }else{
      best.sub.which <- c(best.sub.which, NA)
      best.sub.loglik <- c(best.sub.loglik, NA)
      BMD_QE4_LN <- c(NA,NA,NA)
      best.sub.weight <- c(best.sub.weight, NA)
    }
  }else{
    best.sub.which <- c(best.sub.which, NA)
    best.sub.loglik <- c(best.sub.loglik, NA)
    BMD_QE4_LN <- c(NA,NA,NA)
    best.sub.weight <- c(best.sub.weight, NA)
  }

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[15]))
  }
  if(prior.weights[15] > 0){
    # print(15)
    select_P4_LN <- fun_cov_selection(model = stanmodels$mP4COV,
                                      model_name = 'P4_LN',
                                      model.none = stanmodels$mP4,
                                      loglik = ifelse(data_LNCOV_all$data$data_type == 2, llfP4_LNI_Cov, llfP4_LND_Cov),
                                      data_asigma2 = data_LNCOV_asigma2,
                                      data_dBMD = data_LNCOV_dBMD,
                                      data_all = data_LNCOV_all,
                                      data_none = data_LN_noCOV,
                                      prior.weightsCov = rep(1, 4), # weights for the 4 sub-models
                                      pvec = pvec,
                                      ndraws = ndraws, td = data_LNCOV_all$data$truncd, seed = seed)
    if(exists('select_P4_LN') & !is.null(select_P4_LN)){
      best.sub.which <- c(best.sub.which, select_P4_LN$best.submodel.name)
      best.sub.loglik <- c(best.sub.loglik, select_P4_LN$Weights$Loglik[select_P4_LN$Weights$Covariate == select_P4_LN$best.submodel.name])
      best.sub.weight <- c(best.sub.weight, select_P4_LN$Weights$Weights[select_P4_LN$Weights$Covariate == select_P4_LN$best.submodel.name])
      BMD_P4_LN <- getBMD(select_P4_LN$best.submodel$theta_tilde, pvec, select_P4_LN$best.submodel.name, 'P4_LN')
    }else{
      best.sub.which <- c(best.sub.which, NA)
      best.sub.loglik <- c(best.sub.loglik, NA)
      BMD_P4_LN <- c(NA,NA,NA)
      best.sub.weight <- c(best.sub.weight, NA)
    }
  }else{
    best.sub.which <- c(best.sub.which, NA)
    best.sub.loglik <- c(best.sub.loglik, NA)
    BMD_P4_LN <- c(NA,NA,NA)
    best.sub.weight <- c(best.sub.weight, NA)
  }

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[16]))
  }
  if(prior.weights[16] > 0){
    # print(16)
    select_L4_LN <- fun_cov_selection(model = stanmodels$mL4COV,
                                      model_name = 'L4_LN',
                                      model.none = stanmodels$mL4,
                                      loglik = ifelse(data_LNCOV_all$data$data_type == 2, llfL4_LNI_Cov, llfL4_LND_Cov),
                                      data_asigma2 = data_LNCOV_asigma2,
                                      data_dBMD = data_LNCOV_dBMD,
                                      data_all = data_LNCOV_all,
                                      data_none = data_LN_noCOV,
                                      prior.weightsCov = rep(1, 4), # weights for the 4 sub-models
                                      pvec = pvec,
                                      ndraws = ndraws, td = data_LNCOV_all$data$truncd, seed = seed)
    if(exists('select_L4_LN') & !is.null(select_L4_LN)){
      best.sub.which <- c(best.sub.which, select_L4_LN$best.submodel.name)
      best.sub.loglik <- c(best.sub.loglik, select_L4_LN$Weights$Loglik[select_L4_LN$Weights$Covariate == select_L4_LN$best.submodel.name])
      best.sub.weight <- c(best.sub.weight, select_L4_LN$Weights$Weights[select_L4_LN$Weights$Covariate == select_L4_LN$best.submodel.name])
      BMD_L4_LN <- getBMD(select_L4_LN$best.submodel$theta_tilde, pvec, select_L4_LN$best.submodel.name, 'L4_LN')
    }else{
      best.sub.which <- c(best.sub.which, NA)
      best.sub.loglik <- c(best.sub.loglik, NA)
      BMD_L4_LN <- c(NA,NA,NA)
      best.sub.weight <- c(best.sub.weight, NA)
    }
  }else{
    best.sub.which <- c(best.sub.which, NA)
    best.sub.loglik <- c(best.sub.loglik, NA)
    BMD_L4_LN <- c(NA,NA,NA)
    best.sub.weight <- c(best.sub.weight, NA)
  }

  #########################################
  ### Weights based on each selected model

  ## loglikelihoods

  max.ll.best = max(best.sub.loglik, na.rm = T)
  minll.best <- min(best.sub.loglik[which((max.ll.best-best.sub.loglik) < 709)])

  w <- c()

  if(prior.weights[1] > 0 & !is.na(best.sub.which[1])){

    if(best.sub.which[1] == 'a_sigma2'){
      w = c(w, fun.w(select_E4_N$a_sigma2, best.sub.loglik[1], minll.best, data_NCOV_asigma2$data$nlevels,
                     data_NCOV_asigma2, 'a_sigma2'))

    }else if(best.sub.which[1] == 'BMD_d'){
      w = c(w, fun.w(select_E4_N$BMD_d, best.sub.loglik[1], minll.best, data_NCOV_dBMD$data$nlevels,
                     data_NCOV_dBMD, 'BMD_d'))

    }else if(best.sub.which[1] == 'all'){
      w = c(w, fun.w(select_E4_N$all, best.sub.loglik[1], minll.best, data_NCOV_all$data$nlevels,
                     data_NCOV_all, 'all'))

    }else if(best.sub.which[1] == 'none'){
      w = c(w, fun.w(select_E4_N$none, best.sub.loglik[1], minll.best, 1,
                     data_N_noCOV, 'none'))
    }
  }else{
    w = c(w, 0)
  }

  if(prior.weights[2] > 0 & !is.na(best.sub.which[2])){

    if(best.sub.which[2] == 'a_sigma2'){
      w = c(w, fun.w(select_IE4_N$a_sigma2, best.sub.loglik[2], minll.best, data_NCOV_asigma2$data$nlevels,
                     data_NCOV_asigma2, 'a_sigma2'))

    }else if(best.sub.which[2] == 'BMD_d'){
      w = c(w, fun.w(select_IE4_N$BMD_d, best.sub.loglik[2], minll.best, data_NCOV_dBMD$data$nlevels,
                     data_NCOV_dBMD, 'BMD_d'))

    }else if(best.sub.which[2] == 'all'){
      w = c(w, fun.w(select_IE4_N$all, best.sub.loglik[2], minll.best, data_NCOV_all$data$nlevels,
                     data_NCOV_all, 'all'))

    }else if(best.sub.which[2] == 'none'){
      w = c(w, fun.w(select_IE4_N$none, best.sub.loglik[2], minll.best, 1,
                     data_N_noCOV, 'none'))
    }
  }else{
    w = c(w, 0)
  }

  if(prior.weights[3] > 0 & !is.na(best.sub.which[3])){

    if(best.sub.which[3] == 'a_sigma2'){
      w = c(w, fun.w(select_H4_N$a_sigma2, best.sub.loglik[3], minll.best, data_NCOV_asigma2$data$nlevels,
                     data_NCOV_asigma2, 'a_sigma2'))

    }else if(best.sub.which[3] == 'BMD_d'){
      w = c(w, fun.w(select_H4_N$BMD_d, best.sub.loglik[3], minll.best, data_NCOV_dBMD$data$nlevels,
                     data_NCOV_dBMD, 'BMD_d'))

    }else if(best.sub.which[3] == 'all'){
      w = c(w, fun.w(select_H4_N$all, best.sub.loglik[3], minll.best, data_NCOV_all$data$nlevels,
                     data_NCOV_all, 'all'))

    }else if(best.sub.which[3] == 'none'){
      w = c(w, fun.w(select_H4_N$none, best.sub.loglik[3], minll.best, 1,
                     data_N_noCOV, 'none'))
    }
  }else{
    w = c(w, 0)
  }

  if(prior.weights[4] > 0 & !is.na(best.sub.which[4])){

    if(best.sub.which[4] == 'a_sigma2'){
      w = c(w, fun.w(select_LN4_N$a_sigma2, best.sub.loglik[4], minll.best, data_NCOV_asigma2$data$nlevels,
                     data_NCOV_asigma2, 'a_sigma2'))

    }else if(best.sub.which[4] == 'BMD_d'){
      w = c(w, fun.w(select_LN4_N$BMD_d, best.sub.loglik[4], minll.best, data_NCOV_dBMD$data$nlevels,
                     data_NCOV_dBMD, 'BMD_d'))

    }else if(best.sub.which[4] == 'all'){
      w = c(w, fun.w(select_LN4_N$all, best.sub.loglik[4], minll.best, data_NCOV_all$data$nlevels,
                     data_NCOV_all, 'all'))

    }else if(best.sub.which[4] == 'none'){
      w = c(w, fun.w(select_LN4_N$none, best.sub.loglik[4], minll.best, 1,
                     data_N_noCOV, 'none'))
    }
  }else{
    w = c(w, 0)
  }

  if(prior.weights[5] > 0 & !is.na(best.sub.which[5])){

    if(best.sub.which[5] == 'a_sigma2'){
      w = c(w, fun.w(select_G4_N$a_sigma2, best.sub.loglik[5], minll.best, data_NCOV_asigma2$data$nlevels,
                     data_NCOV_asigma2, 'a_sigma2'))

    }else if(best.sub.which[5] == 'BMD_d'){
      w = c(w, fun.w(select_G4_N$BMD_d, best.sub.loglik[5], minll.best, data_NCOV_dBMD$data$nlevels,
                     data_NCOV_dBMD, 'BMD_d'))

    }else if(best.sub.which[5] == 'all'){
      w = c(w, fun.w(select_G4_N$all, best.sub.loglik[5], minll.best, data_NCOV_all$data$nlevels,
                     data_NCOV_all, 'all'))

    }else if(best.sub.which[5] == 'none'){
      w = c(w, fun.w(select_G4_N$none, best.sub.loglik[5], minll.best, 1,
                     data_N_noCOV, 'none'))
    }
  }else{
    w = c(w, 0)
  }

  if(prior.weights[6] > 0 & !is.na(best.sub.which[6])){

    if(best.sub.which[6] == 'a_sigma2'){
      w = c(w, fun.w(select_QE4_N$a_sigma2, best.sub.loglik[6], minll.best, data_NCOV_asigma2$data$nlevels,
                     data_NCOV_asigma2, 'a_sigma2'))

    }else if(best.sub.which[6] == 'BMD_d'){
      w = c(w, fun.w(select_QE4_N$BMD_d, best.sub.loglik[6], minll.best, data_NCOV_dBMD$data$nlevels,
                     data_NCOV_dBMD, 'BMD_d'))

    }else if(best.sub.which[6] == 'all'){
      w = c(w, fun.w(select_QE4_N$all, best.sub.loglik[6], minll.best, data_NCOV_all$data$nlevels,
                     data_NCOV_all, 'all'))

    }else if(best.sub.which[6] == 'none'){
      w = c(w, fun.w(select_QE4_N$none, best.sub.loglik[6], minll.best, 1,
                     data_N_noCOV, 'none'))
    }
  }else{
    w = c(w, 0)
  }

  if(prior.weights[7] > 0 & !is.na(best.sub.which[7])){

    if(best.sub.which[7] == 'a_sigma2'){
      w = c(w, fun.w(select_P4_N$a_sigma2, best.sub.loglik[7], minll.best, data_NCOV_asigma2$data$nlevels,
                     data_NCOV_asigma2, 'a_sigma2'))

    }else if(best.sub.which[7] == 'BMD_d'){
      w = c(w, fun.w(select_P4_N$BMD_d, best.sub.loglik[7], minll.best, data_NCOV_dBMD$data$nlevels,
                     data_NCOV_dBMD, 'BMD_d'))

    }else if(best.sub.which[7] == 'all'){
      w = c(w, fun.w(select_P4_N$all, best.sub.loglik[7], minll.best, data_NCOV_all$data$nlevels,
                     data_NCOV_all, 'all'))

    }else if(best.sub.which[7] == 'none'){
      w = c(w, fun.w(select_P4_N$none, best.sub.loglik[7], minll.best, 1,
                     data_N_noCOV, 'none'))
    }
  }else{
    w = c(w, 0)
  }

  if(prior.weights[8] > 0 & !is.na(best.sub.which[8])){

    if(best.sub.which[8] == 'a_sigma2'){
      w = c(w, fun.w(select_L4_N$a_sigma2, best.sub.loglik[8], minll.best, data_NCOV_asigma2$data$nlevels,
                     data_NCOV_asigma2, 'a_sigma2'))

    }else if(best.sub.which[8] == 'BMD_d'){
      w = c(w, fun.w(select_L4_N$BMD_d, best.sub.loglik[8], minll.best, data_NCOV_dBMD$data$nlevels,
                     data_NCOV_dBMD, 'BMD_d'))

    }else if(best.sub.which[8] == 'all'){
      w = c(w, fun.w(select_L4_N$all, best.sub.loglik[8], minll.best, data_NCOV_all$data$nlevels,
                     data_NCOV_all, 'all'))

    }else if(best.sub.which[8] == 'none'){
      w = c(w, fun.w(select_L4_N$none, best.sub.loglik[8], minll.best, 1,
                     data_N_noCOV, 'none'))
    }
  }else{
    w = c(w, 0)
  }

  if(prior.weights[9] > 0 & !is.na(best.sub.which[9])){

    if(best.sub.which[9] == 'a_sigma2'){
      w = c(w, fun.w(select_E4_LN$a_sigma2, best.sub.loglik[9], minll.best, data_LNCOV_asigma2$data$nlevels,
                     data_LNCOV_asigma2, 'a_sigma2'))

    }else if(best.sub.which[9] == 'BMD_d'){
      w = c(w, fun.w(select_E4_LN$BMD_d, best.sub.loglik[9], minll.best, data_LNCOV_dBMD$data$nlevels,
                     data_LNCOV_dBMD, 'BMD_d'))

    }else if(best.sub.which[9] == 'all'){
      w = c(w, fun.w(select_E4_LN$all, best.sub.loglik[9], minll.best, data_LNCOV_all$data$nlevels,
                     data_LNCOV_all, 'all'))

    }else if(best.sub.which[9] == 'none'){
      w = c(w, fun.w(select_E4_LN$none, best.sub.loglik[9], minll.best, 1,
                     data_LN_noCOV, 'none'))
    }
  }else{
    w = c(w, 0)
  }

  if(prior.weights[10] > 0 & !is.na(best.sub.which[10])){

    if(best.sub.which[10] == 'a_sigma2'){
      w = c(w, fun.w(select_IE4_LN$a_sigma2, best.sub.loglik[10], minll.best, data_LNCOV_asigma2$data$nlevels,
                     data_LNCOV_asigma2, 'a_sigma2'))

    }else if(best.sub.which[10] == 'BMD_d'){
      w = c(w, fun.w(select_IE4_LN$BMD_d, best.sub.loglik[10], minll.best, data_LNCOV_dBMD$data$nlevels,
                     data_LNCOV_dBMD, 'BMD_d'))

    }else if(best.sub.which[10] == 'all'){
      w = c(w, fun.w(select_IE4_LN$all, best.sub.loglik[10], minll.best, data_LNCOV_all$data$nlevels,
                     data_LNCOV_all, 'all'))

    }else if(best.sub.which[10] == 'none'){
      w = c(w, fun.w(select_IE4_LN$none, best.sub.loglik[10], minll.best, 1,
                     data_LN_noCOV, 'none'))
    }
  }else{
    w = c(w, 0)
  }

  if(prior.weights[11] > 0 & !is.na(best.sub.which[11])){

    if(best.sub.which[11] == 'a_sigma2'){
      w = c(w, fun.w(select_H4_LN$a_sigma2, best.sub.loglik[11], minll.best, data_LNCOV_asigma2$data$nlevels,
                     data_LNCOV_asigma2, 'a_sigma2'))

    }else if(best.sub.which[11] == 'BMD_d'){
      w = c(w, fun.w(select_H4_LN$BMD_d, best.sub.loglik[11], minll.best, data_LNCOV_dBMD$data$nlevels,
                     data_LNCOV_dBMD, 'BMD_d'))

    }else if(best.sub.which[11] == 'all'){
      w = c(w, fun.w(select_H4_LN$all, best.sub.loglik[11], minll.best, data_LNCOV_all$data$nlevels,
                     data_LNCOV_all, 'all'))

    }else if(best.sub.which[11] == 'none'){
      w = c(w, fun.w(select_H4_LN$none, best.sub.loglik[11], minll.best, 1,
                     data_LN_noCOV, 'none'))
    }
  }else{
    w = c(w, 0)
  }

  if(prior.weights[12] > 0 & !is.na(best.sub.which[12])){

    if(best.sub.which[12] == 'a_sigma2'){
      w = c(w, fun.w(select_LN4_LN$a_sigma2, best.sub.loglik[12], minll.best, data_LNCOV_asigma2$data$nlevels,
                     data_LNCOV_asigma2, 'a_sigma2'))

    }else if(best.sub.which[12] == 'BMD_d'){
      w = c(w, fun.w(select_LN4_LN$BMD_d, best.sub.loglik[12], minll.best, data_LNCOV_dBMD$data$nlevels,
                     data_LNCOV_dBMD, 'BMD_d'))

    }else if(best.sub.which[12] == 'all'){
      w = c(w, fun.w(select_LN4_LN$all, best.sub.loglik[12], minll.best, data_LNCOV_all$data$nlevels,
                     data_LNCOV_all, 'all'))

    }else if(best.sub.which[12] == 'none'){
      w = c(w, fun.w(select_LN4_LN$none, best.sub.loglik[12], minll.best, 1,
                     data_LN_noCOV, 'none'))
    }
  }else{
    w = c(w, 0)
  }

  if(prior.weights[13] > 0 & !is.na(best.sub.which[13])){

    if(best.sub.which[13] == 'a_sigma2'){
      w = c(w, fun.w(select_G4_LN$a_sigma2, best.sub.loglik[13], minll.best, data_LNCOV_asigma2$data$nlevels,
                     data_LNCOV_asigma2, 'a_sigma2'))

    }else if(best.sub.which[13] == 'BMD_d'){
      w = c(w, fun.w(select_G4_LN$BMD_d, best.sub.loglik[13], minll.best, data_LNCOV_dBMD$data$nlevels,
                     data_LNCOV_dBMD, 'BMD_d'))

    }else if(best.sub.which[13] == 'all'){
      w = c(w, fun.w(select_G4_LN$all, best.sub.loglik[13], minll.best, data_LNCOV_all$data$nlevels,
                     data_LNCOV_all, 'all'))

    }else if(best.sub.which[13] == 'none'){
      w = c(w, fun.w(select_G4_LN$none, best.sub.loglik[13], minll.best, 1,
                     data_LN_noCOV, 'none'))
    }
  }else{
    w = c(w, 0)
  }

  if(prior.weights[14] > 0 & !is.na(best.sub.which[14])){

    if(best.sub.which[14] == 'a_sigma2'){
      w = c(w, fun.w(select_QE4_LN$a_sigma2, best.sub.loglik[14], minll.best, data_LNCOV_asigma2$data$nlevels,
                     data_LNCOV_asigma2, 'a_sigma2'))

    }else if(best.sub.which[14] == 'BMD_d'){
      w = c(w, fun.w(select_QE4_LN$BMD_d, best.sub.loglik[14], minll.best, data_LNCOV_dBMD$data$nlevels,
                     data_LNCOV_dBMD, 'BMD_d'))

    }else if(best.sub.which[14] == 'all'){
      w = c(w, fun.w(select_QE4_LN$all, best.sub.loglik[14], minll.best, data_LNCOV_all$data$nlevels,
                     data_LNCOV_all, 'all'))

    }else if(best.sub.which[14] == 'none'){
      w = c(w, fun.w(select_QE4_LN$none, best.sub.loglik[14], minll.best, 1,
                     data_LN_noCOV, 'none'))
    }
  }else{
    w = c(w, 0)
  }

  if(prior.weights[15] > 0 & !is.na(best.sub.which[15])){

    if(best.sub.which[15] == 'a_sigma2'){
      w = c(w, fun.w(select_P4_LN$a_sigma2, best.sub.loglik[15], minll.best, data_LNCOV_asigma2$data$nlevels,
                     data_LNCOV_asigma2, 'a_sigma2'))

    }else if(best.sub.which[15] == 'BMD_d'){
      w = c(w, fun.w(select_P4_LN$BMD_d, best.sub.loglik[15], minll.best, data_LNCOV_dBMD$data$nlevels,
                     data_LNCOV_dBMD, 'BMD_d'))

    }else if(best.sub.which[15] == 'all'){
      w = c(w, fun.w(select_P4_LN$all, best.sub.loglik[15], minll.best, data_LNCOV_all$data$nlevels,
                     data_LNCOV_all, 'all'))

    }else if(best.sub.which[15] == 'none'){
      w = c(w, fun.w(select_P4_LN$none, best.sub.loglik[15], minll.best, 1,
                     data_LN_noCOV, 'none'))
    }
  }else{
    w = c(w, 0)
  }

  if(prior.weights[16] > 0 & !is.na(best.sub.which[16])){

    if(best.sub.which[16] == 'a_sigma2'){
      w = c(w, fun.w(select_L4_LN$a_sigma2, best.sub.loglik[16], minll.best, data_LNCOV_asigma2$data$nlevels,
                     data_LNCOV_asigma2, 'a_sigma2'))

    }else if(best.sub.which[16] == 'BMD_d'){
      w = c(w, fun.w(select_L4_LN$BMD_d, best.sub.loglik[16], minll.best, data_LNCOV_dBMD$data$nlevels,
                     data_LNCOV_dBMD, 'BMD_d'))

    }else if(best.sub.which[16] == 'all'){
      w = c(w, fun.w(select_L4_LN$all, best.sub.loglik[16], minll.best, data_LNCOV_all$data$nlevels,
                     data_LNCOV_all, 'all'))

    }else if(best.sub.which[16] == 'none'){
      w = c(w, fun.w(select_L4_LN$none, best.sub.loglik[16], minll.best, 1,
                     data_LN_noCOV, 'none'))
    }
  }else{
    w = c(w, 0)
  }

  w <- ifelse(w == 'Inf', 0, w)
  prior.weights = prior.weights/sum(prior.weights==1)
  lpw=(prior.weights*w)/sum(prior.weights*w)

  #####################################
  ### MODEL AVERAGING BMD

  count = round(lpw*ndraws)

  ## If at least one submodel has different BMDs

  # if('all' %in% best.sub.which | 'BMD_d' %in% best.sub.which){

  ## if that model gets count > 0 --> covariate-specific BMD

  if(TRUE %in% (count[best.sub.which=='all' | best.sub.which=='BMD_d'] > 0)){

    max.levels <- data_NCOV_all$data$nlevels

    mabmd <- list()

    for(i in 1:max.levels){

      mabmd[[i]] <- c(

        if(prior.weights[1] > 0 & !is.na(best.sub.which[1])){
          if(best.sub.which[1] %in% c('all', 'BMD_d')){
            sample(select_E4_N$best.submodel$theta_tilde[,grep("par2\\[",colnames(select_E4_N$best.submodel$theta_tilde))][,i],
                   count[1], replace = T)
          }else if(best.sub.which[1] == 'none'){
            sample(select_E4_N$best.submodel$theta_tilde[,'par2'],
                   count[1], replace = T)
          }else{
            sample(select_E4_N$best.submodel$theta_tilde[,'par2[1]'],
                   count[1], replace = T)
          }
        },
        if(prior.weights[2] > 0 & !is.na(best.sub.which[2])){
          if(best.sub.which[2] %in% c('all', 'BMD_d')){
            sample(select_IE4_N$best.submodel$theta_tilde[,grep("par2\\[",colnames(select_IE4_N$best.submodel$theta_tilde))][,i],
                   count[2], replace = T)
          }else if(best.sub.which[2] == 'none'){
            sample(select_IE4_N$best.submodel$theta_tilde[,'par2'],
                   count[2], replace = T)
          }else{
            sample(select_IE4_N$best.submodel$theta_tilde[,'par2[1]'],
                   count[2], replace = T)
          }
        },
        if(prior.weights[3] > 0 & !is.na(best.sub.which[3])){
          if(best.sub.which[3] %in% c('all', 'BMD_d')){
            sample(select_H4_N$best.submodel$theta_tilde[,grep("par2\\[",colnames(select_H4_N$best.submodel$theta_tilde))][,i],
                   count[3], replace = T)
          }else if(best.sub.which[3] == 'none'){
            sample(select_H4_N$best.submodel$theta_tilde[,'par2'],
                   count[3], replace = T)
          }else{
            sample(select_H4_N$best.submodel$theta_tilde[,'par2[1]'],
                   count[3], replace = T)
          }
        },
        if(prior.weights[4] > 0 & !is.na(best.sub.which[4])){
          if(best.sub.which[4] %in% c('all', 'BMD_d')){
            sample(select_LN4_N$best.submodel$theta_tilde[,grep("par2\\[",colnames(select_LN4_N$best.submodel$theta_tilde))][,i],
                   count[4], replace = T)
          }else if(best.sub.which[4] == 'none'){
            sample(select_LN4_N$best.submodel$theta_tilde[,'par2'],
                   count[4], replace = T)
          }else{
            sample(select_LN4_N$best.submodel$theta_tilde[,'par2[1]'],
                   count[4], replace = T)
          }
        },
        if(prior.weights[5] > 0 & !is.na(best.sub.which[5])){
          if(best.sub.which[5] %in% c('all', 'BMD_d')){
            sample(select_G4_N$best.submodel$theta_tilde[,grep("par2\\[",colnames(select_G4_N$best.submodel$theta_tilde))][,i],
                   count[5], replace = T)
          }else if(best.sub.which[5] == 'none'){
            sample(select_G4_N$best.submodel$theta_tilde[,'par2'],
                   count[5], replace = T)
          }else{
            sample(select_G4_N$best.submodel$theta_tilde[,'par2[1]'],
                   count[5], replace = T)
          }
        },
        if(prior.weights[6] > 0 & !is.na(best.sub.which[6])){
          if(best.sub.which[6] %in% c('all', 'BMD_d')){
            sample(select_QE4_N$best.submodel$theta_tilde[,grep("par2\\[",colnames(select_QE4_N$best.submodel$theta_tilde))][,i],
                   count[6], replace = T)
          }else if(best.sub.which[6] == 'none'){
            sample(select_QE4_N$best.submodel$theta_tilde[,'par2'],
                   count[6], replace = T)
          }else{
            sample(select_QE4_N$best.submodel$theta_tilde[,'par2[1]'],
                   count[6], replace = T)
          }
        },
        if(prior.weights[7] > 0 & !is.na(best.sub.which[7])){
          if(best.sub.which[7] %in% c('all', 'BMD_d')){
            sample(select_P4_N$best.submodel$theta_tilde[,grep("par2\\[",colnames(select_P4_N$best.submodel$theta_tilde))][,i],
                   count[7], replace = T)
          }else if(best.sub.which[7] == 'none'){
            sample(select_P4_N$best.submodel$theta_tilde[,'par2'],
                   count[7], replace = T)
          }else{
            sample(select_P4_N$best.submodel$theta_tilde[,'par2[1]'],
                   count[7], replace = T)
          }
        },
        if(prior.weights[8] > 0 & !is.na(best.sub.which[8])){
          if(best.sub.which[8] %in% c('all', 'BMD_d')){
            sample(select_L4_N$best.submodel$theta_tilde[,grep("par2\\[",colnames(select_L4_N$best.submodel$theta_tilde))][,i],
                   count[8], replace = T)
          }else if(best.sub.which[8] == 'none'){
            sample(select_L4_N$best.submodel$theta_tilde[,'par2'],
                   count[8], replace = T)
          }else{
            sample(select_L4_N$best.submodel$theta_tilde[,'par2[1]'],
                   count[8], replace = T)
          }
        },
        if(prior.weights[9] > 0 & !is.na(best.sub.which[9])){
          if(best.sub.which[9] %in% c('all', 'BMD_d')){
            sample(select_E4_LN$best.submodel$theta_tilde[,grep("par2\\[",colnames(select_E4_LN$best.submodel$theta_tilde))][,i],
                   count[9], replace = T)
          }else if(best.sub.which[9] == 'none'){
            sample(select_E4_LN$best.submodel$theta_tilde[,'par2'],
                   count[9], replace = T)
          }else{
            sample(select_E4_LN$best.submodel$theta_tilde[,'par2[1]'],
                   count[9], replace = T)
          }
        },
        if(prior.weights[10] > 0 & !is.na(best.sub.which[10])){
          if(best.sub.which[10] %in% c('all', 'BMD_d')){
            sample(select_IE4_LN$best.submodel$theta_tilde[,grep("par2\\[",colnames(select_IE4_LN$best.submodel$theta_tilde))][,i],
                   count[10], replace = T)
          }else if(best.sub.which[10] == 'none'){
            sample(select_IE4_LN$best.submodel$theta_tilde[,'par2'],
                   count[10], replace = T)
          }else{
            sample(select_IE4_LN$best.submodel$theta_tilde[,'par2[1]'],
                   count[10], replace = T)
          }
        },
        if(prior.weights[11] > 0 & !is.na(best.sub.which[11])){
          if(best.sub.which[11] %in% c('all', 'BMD_d')){
            sample(select_H4_LN$best.submodel$theta_tilde[,grep("par2\\[",colnames(select_H4_LN$best.submodel$theta_tilde))][,i],
                   count[11], replace = T)
          }else if(best.sub.which[11] == 'none'){
            sample(select_H4_LN$best.submodel$theta_tilde[,'par2'],
                   count[11], replace = T)
          }else{
            sample(select_H4_LN$best.submodel$theta_tilde[,'par2[1]'],
                   count[11], replace = T)
          }
        },
        if(prior.weights[12] > 0 & !is.na(best.sub.which[12])){
          if(best.sub.which[12] %in% c('all', 'BMD_d')){
            sample(select_LN4_LN$best.submodel$theta_tilde[,grep("par2\\[",colnames(select_LN4_LN$best.submodel$theta_tilde))][,i],
                   count[12], replace = T)
          }else if(best.sub.which[12] == 'none'){
            sample(select_LN4_LN$best.submodel$theta_tilde[,'par2'],
                   count[12], replace = T)
          }else{
            sample(select_LN4_LN$best.submodel$theta_tilde[,'par2[1]'],
                   count[12], replace = T)
          }
        },
        if(prior.weights[13] > 0 & !is.na(best.sub.which[13])){
          if(best.sub.which[13] %in% c('all', 'BMD_d')){
            sample(select_G4_LN$best.submodel$theta_tilde[,grep("par2\\[",colnames(select_G4_LN$best.submodel$theta_tilde))][,i],
                   count[13], replace = T)
          }else if(best.sub.which[13] == 'none'){
            sample(select_G4_LN$best.submodel$theta_tilde[,'par2'],
                   count[13], replace = T)
          }else{
            sample(select_G4_LN$best.submodel$theta_tilde[,'par2[1]'],
                   count[13], replace = T)
          }
        },
        if(prior.weights[14] > 0 & !is.na(best.sub.which[14])){
          if(best.sub.which[14] %in% c('all', 'BMD_d')){
            sample(select_QE4_LN$best.submodel$theta_tilde[,grep("par2\\[",colnames(select_QE4_LN$best.submodel$theta_tilde))][,i],
                   count[14], replace = T)
          }else if(best.sub.which[14] == 'none'){
            sample(select_QE4_LN$best.submodel$theta_tilde[,'par2'],
                   count[14], replace = T)
          }else{
            sample(select_QE4_LN$best.submodel$theta_tilde[,'par2[1]'],
                   count[14], replace = T)
          }
        },
        if(prior.weights[15] > 0 & !is.na(best.sub.which[15])){
          if(best.sub.which[15] %in% c('all', 'BMD_d')){
            sample(select_P4_LN$best.submodel$theta_tilde[,grep("par2\\[",colnames(select_P4_LN$best.submodel$theta_tilde))][,i],
                   count[15], replace = T)
          }else if(best.sub.which[15] == 'none'){
            sample(select_P4_LN$best.submodel$theta_tilde[,'par2'],
                   count[15], replace = T)
          }else{
            sample(select_P4_LN$best.submodel$theta_tilde[,'par2[1]'],
                   count[15], replace = T)
          }
        },
        if(prior.weights[16] > 0 & !is.na(best.sub.which[16])){
          if(best.sub.which[16] %in% c('all', 'BMD_d')){
            sample(select_L4_LN$best.submodel$theta_tilde[,grep("par2\\[",colnames(select_L4_LN$best.submodel$theta_tilde))][,i],
                   count[16], replace = T)
          }else if(best.sub.which[16] == 'none'){
            sample(select_L4_LN$best.submodel$theta_tilde[,'par2'],
                   count[16], replace = T)
          }else{
            sample(select_L4_LN$best.submodel$theta_tilde[,'par2[1]'],
                   count[16], replace = T)
          }
        }

      )

    }

    maci <- matrix(NA, max.levels, 3)
    for(i in 1:max.levels){
      maci[i, ] <- quantile(mabmd[[i]], pvec, na.rm = T)*data_NCOV_all$data$maxD
    }

    colnames(maci)=c("BMDL","BMD","BMDU")
    rownames(maci)=data_NCOV_all$data$covariate

    BMDq <- list()
    for(i in 1:max.levels){
      BMDq[[i]] = quantile(mabmd[[i]], seq(0,1,0.005), na.rm = T)*data_N_noCOV$data$maxD ## original scale
    }

    model = c("E4_N","IE4_N","H4_N","LN4_N","G4_N","QE4_N","P4_N","L4_N","E4_LN","IE4_LN","H4_LN",
              "LN4_LN","G4_LN","QE4_LN","P4_LN","L4_LN","MA")
    model = as.factor(model)
    weight = c(lpw[1], lpw[2], lpw[3], lpw[4], lpw[5], lpw[6], lpw[7], lpw[8],
               lpw[9], lpw[10], lpw[11], lpw[12], lpw[13], lpw[14], lpw[15], lpw[16], 1)

    names(lpw) = model[1:16]


  }else{
    ## else only one BMD

    mabmd <- c(

      if(prior.weights[1] > 0 & !is.na(best.sub.which[1])){
        if(best.sub.which[1] == 'a_sigma2'){
          sample(select_E4_N$best.submodel$theta_tilde[,'par2[1]'], count[1], replace = T)
        }else if(best.sub.which[1] == 'none'){
          sample(select_E4_N$best.submodel$theta_tilde[,'par2'], count[1], replace = T)
        }
      },
      if(prior.weights[2] > 0 & !is.na(best.sub.which[2])){
        if(best.sub.which[2] == 'a_sigma2'){
          sample(select_IE4_N$best.submodel$theta_tilde[,'par2[1]'], count[2], replace = T)
        }else if(best.sub.which[2] == 'none'){
          sample(select_IE4_N$best.submodel$theta_tilde[,'par2'], count[2], replace = T)
        }
      },
      if(prior.weights[3] > 0 & !is.na(best.sub.which[3])){
        if(best.sub.which[3] == 'a_sigma2'){
          sample(select_H4_N$best.submodel$theta_tilde[,'par2[1]'], count[3], replace = T)
        }else if(best.sub.which[3] == 'none'){
          sample(select_H4_N$best.submodel$theta_tilde[,'par2'], count[3], replace = T)
        }
      },
      if(prior.weights[4] > 0 & !is.na(best.sub.which[4])){
        if(best.sub.which[4] == 'a_sigma2'){
          sample(select_LN4_N$best.submodel$theta_tilde[,'par2[1]'], count[4], replace = T)
        }else if(best.sub.which[4] == 'none'){
          sample(select_LN4_N$best.submodel$theta_tilde[,'par2'], count[4], replace = T)
        }
      },
      if(prior.weights[5] > 0 & !is.na(best.sub.which[5])){
        if(best.sub.which[5] == 'a_sigma2'){
          sample(select_G4_N$best.submodel$theta_tilde[,'par2[1]'], count[5], replace = T)
        }else if(best.sub.which[5] == 'none'){
          sample(select_G4_N$best.submodel$theta_tilde[,'par2'], count[5], replace = T)
        }
      },
      if(prior.weights[6] > 0 & !is.na(best.sub.which[6])){
        if(best.sub.which[6] == 'a_sigma2'){
          sample(select_QE4_N$best.submodel$theta_tilde[,'par2[1]'], count[6], replace = T)
        }else if(best.sub.which[6] == 'none'){
          sample(select_QE4_N$best.submodel$theta_tilde[,'par2'], count[6], replace = T)
        }
      },
      if(prior.weights[7] > 0 & !is.na(best.sub.which[7])){
        if(best.sub.which[7] == 'a_sigma2'){
          sample(select_P4_N$best.submodel$theta_tilde[,'par2[1]'], count[7], replace = T)
        }else if(best.sub.which[7] == 'none'){
          sample(select_P4_N$best.submodel$theta_tilde[,'par2'], count[7], replace = T)
        }
      },
      if(prior.weights[8] > 0 & !is.na(best.sub.which[8])){
        if(best.sub.which[8] == 'a_sigma2'){
          sample(select_L4_N$best.submodel$theta_tilde[,'par2[1]'], count[8], replace = T)
        }else if(best.sub.which[8] == 'none'){
          sample(select_L4_N$best.submodel$theta_tilde[,'par2'], count[8], replace = T)
        }
      },
      if(prior.weights[9] > 0 & !is.na(best.sub.which[9])){
        if(best.sub.which[9] == 'a_sigma2'){
          sample(select_E4_LN$best.submodel$theta_tilde[,'par2[1]'], count[9], replace = T)
        }else if(best.sub.which[9] == 'none'){
          sample(select_E4_LN$best.submodel$theta_tilde[,'par2'], count[9], replace = T)
        }
      },
      if(prior.weights[10] > 0 & !is.na(best.sub.which[10])){
        if(best.sub.which[10] == 'a_sigma2'){
          sample(select_IE4_LN$best.submodel$theta_tilde[,'par2[1]'], count[10], replace = T)
        }else if(best.sub.which[10] == 'none'){
          sample(select_IE4_LN$best.submodel$theta_tilde[,'par2'], count[10], replace = T)
        }
      },
      if(prior.weights[11] > 0 & !is.na(best.sub.which[11])){
        if(best.sub.which[11] == 'a_sigma2'){
          sample(select_H4_LN$best.submodel$theta_tilde[,'par2[1]'], count[11], replace = T)
        }else if(best.sub.which[11] == 'none'){
          sample(select_H4_LN$best.submodel$theta_tilde[,'par2'], count[11], replace = T)
        }
      },
      if(prior.weights[12] > 0 & !is.na(best.sub.which[12])){
        if(best.sub.which[12] == 'a_sigma2'){
          sample(select_LN4_LN$best.submodel$theta_tilde[,'par2[1]'], count[12], replace = T)
        }else if(best.sub.which[12] == 'none'){
          sample(select_LN4_LN$best.submodel$theta_tilde[,'par2'], count[12], replace = T)
        }
      },
      if(prior.weights[13] > 0 & !is.na(best.sub.which[13])){
        if(best.sub.which[13] == 'a_sigma2'){
          sample(select_G4_LN$best.submodel$theta_tilde[,'par2[1]'], count[13], replace = T)
        }else if(best.sub.which[13] == 'none'){
          sample(select_G4_LN$best.submodel$theta_tilde[,'par2'], count[13], replace = T)
        }
      },
      if(prior.weights[14] > 0 & !is.na(best.sub.which[14])){
        if(best.sub.which[14] == 'a_sigma2'){
          sample(select_QE4_LN$best.submodel$theta_tilde[,'par2[1]'], count[14], replace = T)
        }else if(best.sub.which[14] == 'none'){
          sample(select_QE4_LN$best.submodel$theta_tilde[,'par2'], count[14], replace = T)
        }
      },
      if(prior.weights[15] > 0 & !is.na(best.sub.which[15])){
        if(best.sub.which[15] == 'a_sigma2'){
          sample(select_P4_LN$best.submodel$theta_tilde[,'par2[1]'], count[15], replace = T)
        }else if(best.sub.which[15] == 'none'){
          sample(select_P4_LN$best.submodel$theta_tilde[,'par2'], count[15], replace = T)
        }
      },
      if(prior.weights[16] > 0 & !is.na(best.sub.which[16])){
        if(best.sub.which[16] == 'a_sigma2'){
          sample(select_L4_LN$best.submodel$theta_tilde[,'par2[1]'], count[16], replace = T)
        }else if(best.sub.which[16] == 'none'){
          sample(select_L4_LN$best.submodel$theta_tilde[,'par2'], count[16], replace = T)
        }
      }
    )
    maci=quantile(mabmd,pvec)*data_N_noCOV$data$maxD ## original scale
    # names(maci)=c("BMDL","BMD","BMDU")
    # names(maci)=c("BMDL","BMD","BMDU")
    maci=t(as.matrix(maci))
    colnames(maci)=c("BMDL","BMD","BMDU")
    rownames(maci)=""

    BMDq = quantile(mabmd, seq(0,1,0.005), na.rm = T)*data_N_noCOV$data$maxD ## original scale

    # BMDL = c(BMDL, maci[1]); BMD = c(BMD, maci[2]); BMDU = c(BMDU, maci[3])

    # names(BMDL) <- c("E4_N","IE4_N","H4_N","LN4_N","G4_N","QE4_N","P4_N","L4_N","E4_LN","IE4_LN","H4_LN",
    #                  "LN4_LN","G4_LN","QE4_LN","P4_LN","L4_LN","MA")
    model = c("E4_N","IE4_N","H4_N","LN4_N","G4_N","QE4_N","P4_N","L4_N","E4_LN","IE4_LN","H4_LN",
              "LN4_LN","G4_LN","QE4_LN","P4_LN","L4_LN","MA")
    model = as.factor(model)
    weight = c(lpw[1], lpw[2], lpw[3], lpw[4], lpw[5], lpw[6], lpw[7], lpw[8],
               lpw[9], lpw[10], lpw[11], lpw[12], lpw[13], lpw[14], lpw[15], lpw[16], 1)

    names(lpw) = model[1:16]

  }

  if('all' %in% best.sub.which | 'BMD_d' %in% best.sub.which){

    BMDL <- c()
    BMD <- c()
    BMDU <- c()

    for(m in model[1:16]){

      if(!TRUE %in% is.na(get(paste0('BMD_', m)))){
        if(!is.vector(get(paste0('BMD_',m)))){
          # if(dim(get(paste0('BMD_', m)))[1] == data_NCOV_all$data$nlevels){
          for(i in 1:data_NCOV_all$data$nlevels){
            BMDL <- c(BMDL, get(paste0('BMD_',m))[i, 1]*data_N_noCOV$data$maxD)
            BMD <- c(BMD, get(paste0('BMD_',m))[i, 2]*data_N_noCOV$data$maxD)
            BMDU <- c(BMDU, get(paste0('BMD_',m))[i, 3]*data_N_noCOV$data$maxD)
          }
        }else{
          BMDL <- c(BMDL, rep(get(paste0('BMD_',m))[1]*data_N_noCOV$data$maxD, each = data_NCOV_all$data$nlevels))
          BMD <- c(BMD, rep(get(paste0('BMD_',m))[2]*data_N_noCOV$data$maxD, each = data_NCOV_all$data$nlevels))
          BMDU <- c(BMDU, rep(get(paste0('BMD_',m))[3]*data_N_noCOV$data$maxD, each = data_NCOV_all$data$nlevels))
        }
      }else{
        BMDL <- c(BMDL, rep(NA, data_NCOV_all$data$nlevels))
        BMD <- c(BMD, rep(NA, data_NCOV_all$data$nlevels))
        BMDU <- c(BMDU, rep(NA, data_NCOV_all$data$nlevels))
      }
    }

    model.weight <- data.frame(Model = rep(model[1:16], each = data_NCOV_all$data$nlevels),
                               Weight = rep(round(lpw[1:16],4), each = data_NCOV_all$data$nlevels),
                               Submodel = rep(best.sub.which, each = data_NCOV_all$data$nlevels),
                               Submodel.weight = rep(best.sub.weight, each = data_NCOV_all$data$nlevels),
                               Covariate = rep(data_NCOV_all$data$covariate, 16),
                               BMDL = BMDL, BMD = BMD, BMDU = BMDU
                               # LogLik = best.sub.loglik,
    )

  }else{

    BMDL <- c()
    BMD <- c()
    BMDU <- c()

    for(m in model[1:16]){
      if(!is.na(get(paste0('BMD_', m))[1])){
        BMDL <- c(BMDL, get(paste0('BMD_',m))[1]*data_N_noCOV$data$maxD)
        BMD <- c(BMD, get(paste0('BMD_',m))[2]*data_N_noCOV$data$maxD)
        BMDU <- c(BMDU, get(paste0('BMD_',m))[3]*data_N_noCOV$data$maxD)

      }else{
        BMDL <- c(BMDL, NA)
        BMD <- c(BMD, NA)
        BMDU <- c(BMDU, NA)
      }
    }

    model.weight <- data.frame(Model = model[1:16],
                               Weight = round(lpw[1:16],4),
                               Submodel = best.sub.which,
                               Submodel.weight = best.sub.weight,
                               Covariate = "None",
                               BMDL = BMDL, BMD = BMD, BMDU = BMDU
                               # LogLik = best.sub.loglik,
    )

  }

  if(exists('select_E4_N')) { parE4_N = select_E4_N$best.submodel$par} else{parE4_N = NULL}
  if(exists('select_IE4_N')) { parIE4_N = select_IE4_N$best.submodel$par} else{parIE4_N = NULL}
  if(exists('select_H4_N')) { parH4_N = select_H4_N$best.submodel$par} else{parH4_N = NULL}
  if(exists('select_LN4_N')) { parLN4_N = select_LN4_N$best.submodel$par} else{parLN4_N = NULL}
  if(exists('select_G4_N')) { parG4_N = select_G4_N$best.submodel$par} else{parG4_N = NULL}
  if(exists('select_QE4_N')) { parQE4_N = select_QE4_N$best.submodel$par} else{parQE4_N = NULL}
  if(exists('select_P4_N')) { parP4_N = select_P4_N$best.submodel$par} else{parP4_N = NULL}
  if(exists('select_L4_N')) { parL4_N = select_L4_N$best.submodel$par} else{parL4_N = NULL}
  if(exists('select_E4_LN')) { parE4_LN = select_E4_LN$best.submodel$par} else{parE4_LN = NULL}
  if(exists('select_IE4_LN')) { parIE4_LN = select_IE4_LN$best.submodel$par} else{parIE4_LN = NULL}
  if(exists('select_H4_LN')) { parH4_LN = select_H4_LN$best.submodel$par} else{parH4_LN = NULL}
  if(exists('select_LN4_LN')) { parLN4_LN = select_LN4_LN$best.submodel$par} else{parLN4_LN = NULL}
  if(exists('select_G4_LN')) { parG4_LN = select_G4_LN$best.submodel$par} else{parG4_LN = NULL}
  if(exists('select_QE4_LN')) { parQE4_LN = select_QE4_LN$best.submodel$par} else{parQE4_LN = NULL}
  if(exists('select_P4_LN')) { parP4_LN = select_P4_LN$best.submodel$par} else{parP4_LN = NULL}
  if(exists('select_L4_LN')) { parL4_LN = select_L4_LN$best.submodel$par} else{parL4_LN = NULL}

  return(list(MA = maci,
              weights = lpw,
              summary = model.weight,
              # best.submodels = best.sub.which,
              # data = data.frame(x = data_NCOV_all$data$x,
              #                   y = data_NCOV_all$data$m,
              #                   s2 = data_NCOV_all$data$s2,
              #                   n = data_NCOV_all$data$n,
              #                   cov = data_NCOV_all$data$covariate
              # ),
              data = data_NCOV_all$data$org.data,
              maxDose = data_N_noCOV$data$maxD,
              BMD_Mixture = BMDq,
              parE4_N = parE4_N,
              parIE4_N = parIE4_N,
              parH4_N = parH4_N,
              parLN4_N = parLN4_N,
              parG4_N = parG4_N,
              parQE4_N = parQE4_N,
              parP4_N = parP4_N,
              parL4_N = parL4_N,
              parE4_LN = parE4_LN,
              parIE4_LN = parIE4_LN,
              parH4_LN = parH4_LN,
              parLN4_LN = parLN4_LN,
              parG4_LN = parG4_LN,
              parQE4_LN = parQE4_LN,
              parP4_LN = parP4_LN,
              parL4_LN = parL4_LN,
              shift = ifelse((1 %in% prior.weights[9:16]), data_LNCOV_all$data$shift, 0),
              q = data_NCOV_all$data$q
  ))

}


#' @rdname full.laplace_MA_Cov
#' @export
full.laplace_MA_Q_Cov = function(data, # the summary data
                                 sumstats = TRUE,
                                 q = 0.1,
                                 prior.d = 'N11',
                                 extended = TRUE, extended.value = 3,
                                 prior.weights = rep(1,8),
                                 ndraws=30000,seed=123,
                                 pvec=c(0.05,0.5,0.95)
){

  out.stop = 'ok'

  data_all <- PREP_DATA_Q_COV(
    data = data,
    sumstats = sumstats,
    q = q,
    prior.d = prior.d,
    extended = extended,
    extended.value = extended.value,
    covariate = 'all'
  )

  data_bkg <- PREP_DATA_Q_COV(
    data = data,
    sumstats = sumstats,
    q = q,
    prior.d = prior.d,
    extended = extended,
    extended.value = extended.value,
    covariate = 'background'
  )

  data_dBMD <- PREP_DATA_Q_COV(
    data = data,
    sumstats = sumstats,
    q = q,
    prior.d = prior.d,
    extended = extended,
    extended.value = extended.value,
    covariate = 'BMD_d'
  )

  data_noCOV <- PREP_DATA_QA(
    data = data,
    sumstats = sumstats,
    q = q,
    extended = extended,
    extended.value = extended.value,
    prior.d = prior.d
  )

  #########################################
  ### SELECT BEST SUBMODEL FOR EACH DRM

  best.sub.which <- c()
  best.sub.loglik <- c()
  best.sub.weight <- c()

  nModels = 8

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='quantal')[1]))
  }
  if(prior.weights[1] > 0){
    select_E4_Q <- fun_cov_selectionQ(model = stanmodels$mE4_Q_COV,
                                      model_name = 'E4_Q',
                                      model.none = stanmodels$mE4_Q,
                                      loglik = llfE4_Q_Cov,
                                      data_bkg = data_bkg,
                                      data_dBMD = data_dBMD,
                                      data_all = data_all,
                                      data_none = data_noCOV,
                                      prior.weightsCov = rep(1, 4), # weights for the 4 sub-models
                                      pvec = pvec,
                                      ndraws = ndraws, td = data_all$data$truncd, seed = seed)
    if(exists('select_E4_Q') & !is.null(select_E4_Q)){
      best.sub.which <- c(best.sub.which, select_E4_Q$best.submodel.name)
      best.sub.loglik <- c(best.sub.loglik, select_E4_Q$Weights$Loglik[select_E4_Q$Weights$Covariate == select_E4_Q$best.submodel.name])
      best.sub.weight <- c(best.sub.weight, select_E4_Q$Weights$Weights[select_E4_Q$Weights$Covariate == select_E4_Q$best.submodel.name])
      BMD_E4_Q <- getBMD_Q(select_E4_Q$best.submodel$theta_tilde, pvec, select_E4_Q$best.submodel.name, 'E4_Q')
    }else{
      best.sub.which <- c(best.sub.which, NA)
      best.sub.loglik <- c(best.sub.loglik, NA)
      BMD_E4_Q <- c(NA,NA,NA)
      best.sub.weight <- c(best.sub.weight, NA)
    }
  }else{
    best.sub.which <- c(best.sub.which, NA)
    best.sub.loglik <- c(best.sub.loglik, NA)
    BMD_E4_Q <- c(NA,NA,NA)
    best.sub.weight <- c(best.sub.weight, NA)
  }

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='quantal')[2]))
  }
  if(prior.weights[2] > 0){
    select_IE4_Q <- fun_cov_selectionQ(model = stanmodels$mIE4_Q_COV,
                                       model_name = 'IE4_Q',
                                       model.none = stanmodels$mIE4_Q,
                                       loglik = llfIE4_Q_Cov,
                                       data_bkg = data_bkg,
                                       data_dBMD = data_dBMD,
                                       data_all = data_all,
                                       data_none = data_noCOV,
                                       prior.weightsCov = rep(1, 4), # weights for the 4 sub-models
                                       pvec = pvec,
                                       ndraws = ndraws, td = data_all$data$truncd, seed = seed)
    if(exists('select_IE4_Q') & !is.null(select_IE4_Q)){
      best.sub.which <- c(best.sub.which, select_IE4_Q$best.submodel.name)
      best.sub.loglik <- c(best.sub.loglik, select_IE4_Q$Weights$Loglik[select_IE4_Q$Weights$Covariate == select_IE4_Q$best.submodel.name])
      best.sub.weight <- c(best.sub.weight, select_IE4_Q$Weights$Weights[select_IE4_Q$Weights$Covariate == select_IE4_Q$best.submodel.name])
      BMD_IE4_Q <- getBMD_Q(select_IE4_Q$best.submodel$theta_tilde, pvec, select_IE4_Q$best.submodel.name, 'IE4_Q')
    }else{
      best.sub.which <- c(best.sub.which, NA)
      best.sub.loglik <- c(best.sub.loglik, NA)
      BMD_IE4_Q <- c(NA,NA,NA)
      best.sub.weight <- c(best.sub.weight, NA)
    }
  }else{
    best.sub.which <- c(best.sub.which, NA)
    best.sub.loglik <- c(best.sub.loglik, NA)
    BMD_IE4_Q <- c(NA,NA,NA)
    best.sub.weight <- c(best.sub.weight, NA)
  }

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='quantal')[3]))
  }
  if(prior.weights[3] > 0){
    select_H4_Q <- fun_cov_selectionQ(model = stanmodels$mH4_Q_COV,
                                      model_name = 'H4_Q',
                                      model.none = stanmodels$mH4_Q,
                                      loglik = llfH4_Q_Cov,
                                      data_bkg = data_bkg,
                                      data_dBMD = data_dBMD,
                                      data_all = data_all,
                                      data_none = data_noCOV,
                                      prior.weightsCov = rep(1, 4), # weights for the 4 sub-models
                                      pvec = pvec,
                                      ndraws = ndraws, td = data_all$data$truncd, seed = seed)
    if(exists('select_H4_Q') & !is.null(select_H4_Q)){
      best.sub.which <- c(best.sub.which, select_H4_Q$best.submodel.name)
      best.sub.loglik <- c(best.sub.loglik, select_H4_Q$Weights$Loglik[select_H4_Q$Weights$Covariate == select_H4_Q$best.submodel.name])
      best.sub.weight <- c(best.sub.weight, select_H4_Q$Weights$Weights[select_H4_Q$Weights$Covariate == select_H4_Q$best.submodel.name])
      BMD_H4_Q <- getBMD_Q(select_H4_Q$best.submodel$theta_tilde, pvec, select_H4_Q$best.submodel.name, 'H4_Q')
    }else{
      best.sub.which <- c(best.sub.which, NA)
      best.sub.loglik <- c(best.sub.loglik, NA)
      BMD_H4_Q <- c(NA,NA,NA)
      best.sub.weight <- c(best.sub.weight, NA)
    }
  }else{
    best.sub.which <- c(best.sub.which, NA)
    best.sub.loglik <- c(best.sub.loglik, NA)
    BMD_H4_Q <- c(NA,NA,NA)
    best.sub.weight <- c(best.sub.weight, NA)
  }

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='quantal')[4]))
  }
  if(prior.weights[4] > 0){
    select_LN4_Q <- fun_cov_selectionQ(model = stanmodels$mLN4_Q_COV,
                                       model_name = 'LN4_Q',
                                       model.none = stanmodels$mLN4_Q,
                                       loglik = llfLN4_Q_Cov,
                                       data_bkg = data_bkg,
                                       data_dBMD = data_dBMD,
                                       data_all = data_all,
                                       data_none = data_noCOV,
                                       prior.weightsCov = rep(1, 4), # weights for the 4 sub-models
                                       pvec = pvec,
                                       ndraws = ndraws, td = data_all$data$truncd, seed = seed)
    if(exists('select_LN4_Q') & !is.null(select_LN4_Q)){
      best.sub.which <- c(best.sub.which, select_LN4_Q$best.submodel.name)
      best.sub.loglik <- c(best.sub.loglik, select_LN4_Q$Weights$Loglik[select_LN4_Q$Weights$Covariate == select_LN4_Q$best.submodel.name])
      best.sub.weight <- c(best.sub.weight, select_LN4_Q$Weights$Weights[select_LN4_Q$Weights$Covariate == select_LN4_Q$best.submodel.name])
      BMD_LN4_Q <- getBMD_Q(select_LN4_Q$best.submodel$theta_tilde, pvec, select_LN4_Q$best.submodel.name, 'LN4_Q')
    }else{
      best.sub.which <- c(best.sub.which, NA)
      best.sub.loglik <- c(best.sub.loglik, NA)
      BMD_LN4_Q <- c(NA,NA,NA)
      best.sub.weight <- c(best.sub.weight, NA)
    }
  }else{
    best.sub.which <- c(best.sub.which, NA)
    best.sub.loglik <- c(best.sub.loglik, NA)
    BMD_LN4_Q <- c(NA,NA,NA)
    best.sub.weight <- c(best.sub.weight, NA)
  }

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='quantal')[5]))
  }
  if(prior.weights[5] > 0){
    select_G4_Q <- fun_cov_selectionQ(model = stanmodels$mG4_Q_COV,
                                      model_name = 'G4_Q',
                                      model.none = stanmodels$mG4_Q,
                                      loglik = llfG4_Q_Cov,
                                      data_bkg = data_bkg,
                                      data_dBMD = data_dBMD,
                                      data_all = data_all,
                                      data_none = data_noCOV,
                                      prior.weightsCov = rep(1, 4), # weights for the 4 sub-models
                                      pvec = pvec,
                                      ndraws = ndraws, td = data_all$data$truncd, seed = seed)
    if(exists('select_G4_Q') & !is.null(select_G4_Q)){
      best.sub.which <- c(best.sub.which, select_G4_Q$best.submodel.name)
      best.sub.loglik <- c(best.sub.loglik, select_G4_Q$Weights$Loglik[select_G4_Q$Weights$Covariate == select_G4_Q$best.submodel.name])
      best.sub.weight <- c(best.sub.weight, select_G4_Q$Weights$Weights[select_G4_Q$Weights$Covariate == select_G4_Q$best.submodel.name])
      BMD_G4_Q <- getBMD_Q(select_G4_Q$best.submodel$theta_tilde, pvec, select_G4_Q$best.submodel.name, 'G4_Q')
    }else{
      best.sub.which <- c(best.sub.which, NA)
      best.sub.loglik <- c(best.sub.loglik, NA)
      BMD_G4_Q <- c(NA,NA,NA)
      best.sub.weight <- c(best.sub.weight, NA)
    }
  }else{
    best.sub.which <- c(best.sub.which, NA)
    best.sub.loglik <- c(best.sub.loglik, NA)
    BMD_G4_Q <- c(NA,NA,NA)
    best.sub.weight <- c(best.sub.weight, NA)
  }

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='quantal')[6]))
  }
  if(prior.weights[6] > 0){
    select_QE4_Q <- fun_cov_selectionQ(model = stanmodels$mQE4_Q_COV,
                                       model_name = 'QE4_Q',
                                       model.none = stanmodels$mQE4_Q,
                                       loglik = llfQE4_Q_Cov,
                                       data_bkg = data_bkg,
                                       data_dBMD = data_dBMD,
                                       data_all = data_all,
                                       data_none = data_noCOV,
                                       prior.weightsCov = rep(1, 4), # weights for the 4 sub-models
                                       pvec = pvec,
                                       ndraws = ndraws, td = data_all$data$truncdQ, seed = seed)
    if(exists('select_QE4_Q') & !is.null(select_QE4_Q)){
      best.sub.which <- c(best.sub.which, select_QE4_Q$best.submodel.name)
      best.sub.loglik <- c(best.sub.loglik, select_QE4_Q$Weights$Loglik[select_QE4_Q$Weights$Covariate == select_QE4_Q$best.submodel.name])
      best.sub.weight <- c(best.sub.weight, select_QE4_Q$Weights$Weights[select_QE4_Q$Weights$Covariate == select_QE4_Q$best.submodel.name])
      BMD_QE4_Q <- getBMD_Q(select_QE4_Q$best.submodel$theta_tilde, pvec, select_QE4_Q$best.submodel.name, 'QE4_Q')
    }else{
      best.sub.which <- c(best.sub.which, NA)
      best.sub.loglik <- c(best.sub.loglik, NA)
      BMD_QE4_Q <- c(NA,NA,NA)
      best.sub.weight <- c(best.sub.weight, NA)
    }
  }else{
    best.sub.which <- c(best.sub.which, NA)
    best.sub.loglik <- c(best.sub.loglik, NA)
    BMD_QE4_Q <- c(NA,NA,NA)
    best.sub.weight <- c(best.sub.weight, NA)
  }

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='quantal')[7]))
  }
  if(prior.weights[7] > 0){
    select_P4_Q <- fun_cov_selectionQ(model = stanmodels$mP4_Q_COV,
                                      model_name = 'P4_Q',
                                      model.none = stanmodels$mP4_Q,
                                      loglik = llfP4_Q_Cov,
                                      data_bkg = data_bkg,
                                      data_dBMD = data_dBMD,
                                      data_all = data_all,
                                      data_none = data_noCOV,
                                      prior.weightsCov = rep(1, 4), # weights for the 4 sub-models
                                      pvec = pvec,
                                      ndraws = ndraws, td = data_all$data$truncd, seed = seed)
    if(exists('select_P4_Q') & !is.null(select_P4_Q)){
      best.sub.which <- c(best.sub.which, select_P4_Q$best.submodel.name)
      best.sub.loglik <- c(best.sub.loglik, select_P4_Q$Weights$Loglik[select_P4_Q$Weights$Covariate == select_P4_Q$best.submodel.name])
      best.sub.weight <- c(best.sub.weight, select_P4_Q$Weights$Weights[select_P4_Q$Weights$Covariate == select_P4_Q$best.submodel.name])
      BMD_P4_Q <- getBMD_Q(select_P4_Q$best.submodel$theta_tilde, pvec, select_P4_Q$best.submodel.name, 'P4_Q')
    }else{
      best.sub.which <- c(best.sub.which, NA)
      best.sub.loglik <- c(best.sub.loglik, NA)
      BMD_P4_Q <- c(NA,NA,NA)
      best.sub.weight <- c(best.sub.weight, NA)
    }
  }else{
    best.sub.which <- c(best.sub.which, NA)
    best.sub.loglik <- c(best.sub.loglik, NA)
    BMD_P4_Q <- c(NA,NA,NA)
    best.sub.weight <- c(best.sub.weight, NA)
  }

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='quantal')[8]))
  }
  if(prior.weights[8] > 0){
    select_L4_Q <- fun_cov_selectionQ(model = stanmodels$mL4_Q_COV,
                                      model_name = 'L4_Q',
                                      model.none = stanmodels$mL4_Q,
                                      loglik = llfL4_Q_Cov,
                                      data_bkg = data_bkg,
                                      data_dBMD = data_dBMD,
                                      data_all = data_all,
                                      data_none = data_noCOV,
                                      prior.weightsCov = rep(1, 4), # weights for the 4 sub-models
                                      pvec = pvec,
                                      ndraws = ndraws, td = data_all$data$truncd, seed = seed)
    if(exists('select_L4_Q') & !is.null(select_L4_Q)){
      best.sub.which <- c(best.sub.which, select_L4_Q$best.submodel.name)
      best.sub.loglik <- c(best.sub.loglik, select_L4_Q$Weights$Loglik[select_L4_Q$Weights$Covariate == select_L4_Q$best.submodel.name])
      best.sub.weight <- c(best.sub.weight, select_L4_Q$Weights$Weights[select_L4_Q$Weights$Covariate == select_L4_Q$best.submodel.name])
      BMD_L4_Q <- getBMD_Q(select_L4_Q$best.submodel$theta_tilde, pvec, select_L4_Q$best.submodel.name, 'L4_Q')
    }else{
      best.sub.which <- c(best.sub.which, NA)
      best.sub.loglik <- c(best.sub.loglik, NA)
      BMD_L4_Q <- c(NA,NA,NA)
      best.sub.weight <- c(best.sub.weight, NA)
    }
  }else{
    best.sub.which <- c(best.sub.which, NA)
    best.sub.loglik <- c(best.sub.loglik, NA)
    BMD_L4_Q <- c(NA,NA,NA)
    best.sub.weight <- c(best.sub.weight, NA)
  }

  #########################################
  ### Weights based on each selected model

  ## loglikelihoods

  max.ll.best = max(best.sub.loglik, na.rm = T)
  minll.best <- min(best.sub.loglik[which((max.ll.best-best.sub.loglik) < 709)])

  w <- c()

  if(prior.weights[1] > 0 & !is.na(best.sub.which[1])){

    if(best.sub.which[1] == 'background'){
      w = c(w, fun.wQ(select_E4_Q$bkg, best.sub.loglik[1], minll.best, data_bkg$data$nlevels,
                      data_bkg, 'background'))

    }else if(best.sub.which[1] == 'BMD_d'){
      w = c(w, fun.wQ(select_E4_Q$BMD_d, best.sub.loglik[1], minll.best, data_dBMD$data$nlevels,
                      data_dBMD, 'BMD_d'))

    }else if(best.sub.which[1] == 'all'){
      w = c(w, fun.wQ(select_E4_Q$all, best.sub.loglik[1], minll.best, data_all$data$nlevels,
                      data_all, 'all'))

    }else if(best.sub.which[1] == 'none'){
      w = c(w, fun.wQ(select_E4_Q$none, best.sub.loglik[1], minll.best, 1,
                      data_noCOV, 'none'))
    }
  }else{
    w = c(w, 0)
  }

  if(prior.weights[2] > 0 & !is.na(best.sub.which[2])){

    if(best.sub.which[2] == 'background'){
      w = c(w, fun.wQ(select_IE4_Q$bkg, best.sub.loglik[2], minll.best, data_bkg$data$nlevels,
                      data_bkg, 'background'))

    }else if(best.sub.which[2] == 'BMD_d'){
      w = c(w, fun.wQ(select_IE4_Q$BMD_d, best.sub.loglik[2], minll.best, data_dBMD$data$nlevels,
                      data_dBMD, 'BMD_d'))

    }else if(best.sub.which[2] == 'all'){
      w = c(w, fun.wQ(select_IE4_Q$all, best.sub.loglik[2], minll.best, data_all$data$nlevels,
                      data_all, 'all'))

    }else if(best.sub.which[2] == 'none'){
      w = c(w, fun.wQ(select_IE4_Q$none, best.sub.loglik[2], minll.best, 1,
                      data_noCOV, 'none'))
    }
  }else{
    w = c(w, 0)
  }

  if(prior.weights[3] > 0 & !is.na(best.sub.which[3])){

    if(best.sub.which[3] == 'background'){
      w = c(w, fun.wQ(select_H4_Q$bkg, best.sub.loglik[3], minll.best, data_bkg$data$nlevels,
                      data_bkg, 'background'))

    }else if(best.sub.which[3] == 'BMD_d'){
      w = c(w, fun.wQ(select_H4_Q$BMD_d, best.sub.loglik[3], minll.best, data_dBMD$data$nlevels,
                      data_dBMD, 'BMD_d'))

    }else if(best.sub.which[3] == 'all'){
      w = c(w, fun.wQ(select_H4_Q$all, best.sub.loglik[3], minll.best, data_all$data$nlevels,
                      data_all, 'all'))

    }else if(best.sub.which[3] == 'none'){
      w = c(w, fun.wQ(select_H4_Q$none, best.sub.loglik[3], minll.best, 1,
                      data_noCOV, 'none'))
    }
  }else{
    w = c(w, 0)
  }

  if(prior.weights[4] > 0 & !is.na(best.sub.which[4])){

    if(best.sub.which[4] == 'background'){
      w = c(w, fun.wQ(select_LN4_Q$bkg, best.sub.loglik[4], minll.best, data_bkg$data$nlevels,
                      data_bkg, 'background'))

    }else if(best.sub.which[4] == 'BMD_d'){
      w = c(w, fun.wQ(select_LN4_Q$BMD_d, best.sub.loglik[4], minll.best, data_dBMD$data$nlevels,
                      data_dBMD, 'BMD_d'))

    }else if(best.sub.which[4] == 'all'){
      w = c(w, fun.wQ(select_LN4_Q$all, best.sub.loglik[4], minll.best, data_all$data$nlevels,
                      data_all, 'all'))

    }else if(best.sub.which[4] == 'none'){
      w = c(w, fun.wQ(select_LN4_Q$none, best.sub.loglik[4], minll.best, 1,
                      data_noCOV, 'none'))
    }
  }else{
    w = c(w, 0)
  }

  if(prior.weights[5] > 0 & !is.na(best.sub.which[5])){

    if(best.sub.which[5] == 'background'){
      w = c(w, fun.wQ(select_G4_Q$bkg, best.sub.loglik[5], minll.best, data_bkg$data$nlevels,
                      data_bkg, 'background'))

    }else if(best.sub.which[5] == 'BMD_d'){
      w = c(w, fun.wQ(select_G4_Q$BMD_d, best.sub.loglik[5], minll.best, data_dBMD$data$nlevels,
                      data_dBMD, 'BMD_d'))

    }else if(best.sub.which[5] == 'all'){
      w = c(w, fun.wQ(select_G4_Q$all, best.sub.loglik[5], minll.best, data_all$data$nlevels,
                      data_all, 'all'))

    }else if(best.sub.which[5] == 'none'){
      w = c(w, fun.wQ(select_G4_Q$none, best.sub.loglik[5], minll.best, 1,
                      data_noCOV, 'none'))
    }
  }else{
    w = c(w, 0)
  }

  if(prior.weights[6] > 0 & !is.na(best.sub.which[6])){

    if(best.sub.which[6] == 'background'){
      w = c(w, fun.wQ.QE4(select_QE4_Q$bkg, best.sub.loglik[6], minll.best, data_bkg$data$nlevels,
                          data_bkg, 'background'))

    }else if(best.sub.which[6] == 'BMD_d'){
      w = c(w, fun.wQ.QE4(select_QE4_Q$BMD_d, best.sub.loglik[6], minll.best, data_dBMD$data$nlevels,
                          data_dBMD, 'BMD_d'))

    }else if(best.sub.which[6] == 'all'){
      w = c(w, fun.wQ.QE4(select_QE4_Q$all, best.sub.loglik[6], minll.best, data_all$data$nlevels,
                          data_all, 'all'))

    }else if(best.sub.which[6] == 'none'){
      w = c(w, fun.wQ.QE4(select_QE4_Q$none, best.sub.loglik[6], minll.best, 1,
                          data_noCOV, 'none'))
    }
  }else{
    w = c(w, 0)
  }

  if(prior.weights[7] > 0 & !is.na(best.sub.which[7])){

    if(best.sub.which[7] == 'background'){
      w = c(w, fun.wQ(select_P4_Q$bkg, best.sub.loglik[7], minll.best, data_bkg$data$nlevels,
                      data_bkg, 'background'))

    }else if(best.sub.which[7] == 'BMD_d'){
      w = c(w, fun.wQ(select_P4_Q$BMD_d, best.sub.loglik[7], minll.best, data_dBMD$data$nlevels,
                      data_dBMD, 'BMD_d'))

    }else if(best.sub.which[7] == 'all'){
      w = c(w, fun.wQ(select_P4_Q$all, best.sub.loglik[7], minll.best, data_all$data$nlevels,
                      data_all, 'all'))

    }else if(best.sub.which[7] == 'none'){
      w = c(w, fun.wQ(select_P4_Q$none, best.sub.loglik[7], minll.best, 1,
                      data_noCOV, 'none'))
    }
  }else{
    w = c(w, 0)
  }

  if(prior.weights[8] > 0 & !is.na(best.sub.which[8])){

    if(best.sub.which[8] == 'background'){
      w = c(w, fun.wQ(select_L4_Q$bkg, best.sub.loglik[8], minll.best, data_bkg$data$nlevels,
                      data_bkg, 'background'))

    }else if(best.sub.which[8] == 'BMD_d'){
      w = c(w, fun.wQ(select_L4_Q$BMD_d, best.sub.loglik[8], minll.best, data_dBMD$data$nlevels,
                      data_dBMD, 'BMD_d'))

    }else if(best.sub.which[8] == 'all'){
      w = c(w, fun.wQ(select_L4_Q$all, best.sub.loglik[8], minll.best, data_all$data$nlevels,
                      data_all, 'all'))

    }else if(best.sub.which[8] == 'none'){
      w = c(w, fun.wQ(select_L4_Q$none, best.sub.loglik[8], minll.best, 1,
                      data_noCOV, 'none'))
    }
  }else{
    w = c(w, 0)
  }

  w <- ifelse(w == 'Inf', 0, w)
  prior.weights = prior.weights/sum(prior.weights==1)
  lpw=(prior.weights*w)/sum(prior.weights*w)

  #####################################
  ### MODEL AVERAGING BMD

  count = round(lpw*ndraws)

  ## If at least one submodel has different BMDs
  ## and only if that model gets count > 0 --> covariate-specific BMD

  if(TRUE %in% (count[best.sub.which=='all' | best.sub.which=='BMD_d'] > 0)){

    max.levels <- data_all$data$nlevels

    mabmd <- list()

    for(i in 1:max.levels){

      mabmd[[i]] <- c(

        if(prior.weights[1] > 0 & !is.na(best.sub.which[1])){
          if(best.sub.which[1] %in% c('all', 'BMD_d')){
            sample(select_E4_Q$best.submodel$theta_tilde[,grep("par2\\[",colnames(select_E4_Q$best.submodel$theta_tilde))][,i],
                   count[1], replace = T)
          }else if(best.sub.which[1] == 'none'){
            sample(select_E4_Q$best.submodel$theta_tilde[,'par2'],
                   count[1], replace = T)
          }else{
            sample(select_E4_Q$best.submodel$theta_tilde[,'par2[1]'],
                   count[1], replace = T)
          }
        },
        if(prior.weights[2] > 0 & !is.na(best.sub.which[2])){
          if(best.sub.which[2] %in% c('all', 'BMD_d')){
            sample(select_IE4_Q$best.submodel$theta_tilde[,grep("par2\\[",colnames(select_IE4_Q$best.submodel$theta_tilde))][,i],
                   count[2], replace = T)
          }else if(best.sub.which[2] == 'none'){
            sample(select_IE4_Q$best.submodel$theta_tilde[,'par2'],
                   count[2], replace = T)
          }else{
            sample(select_IE4_Q$best.submodel$theta_tilde[,'par2[1]'],
                   count[2], replace = T)
          }
        },
        if(prior.weights[3] > 0 & !is.na(best.sub.which[3])){
          if(best.sub.which[3] %in% c('all', 'BMD_d')){
            sample(select_H4_Q$best.submodel$theta_tilde[,grep("par2\\[",colnames(select_H4_Q$best.submodel$theta_tilde))][,i],
                   count[3], replace = T)
          }else if(best.sub.which[3] == 'none'){
            sample(select_H4_Q$best.submodel$theta_tilde[,'par2'],
                   count[3], replace = T)
          }else{
            sample(select_H4_Q$best.submodel$theta_tilde[,'par2[1]'],
                   count[3], replace = T)
          }
        },
        if(prior.weights[4] > 0 & !is.na(best.sub.which[4])){
          if(best.sub.which[4] %in% c('all', 'BMD_d')){
            sample(select_LN4_Q$best.submodel$theta_tilde[,grep("par2\\[",colnames(select_LN4_Q$best.submodel$theta_tilde))][,i],
                   count[4], replace = T)
          }else if(best.sub.which[4] == 'none'){
            sample(select_LN4_Q$best.submodel$theta_tilde[,'par2'],
                   count[4], replace = T)
          }else{
            sample(select_LN4_Q$best.submodel$theta_tilde[,'par2[1]'],
                   count[4], replace = T)
          }
        },
        if(prior.weights[5] > 0 & !is.na(best.sub.which[5])){
          if(best.sub.which[5] %in% c('all', 'BMD_d')){
            sample(select_G4_Q$best.submodel$theta_tilde[,grep("par2\\[",colnames(select_G4_Q$best.submodel$theta_tilde))][,i],
                   count[5], replace = T)
          }else if(best.sub.which[5] == 'none'){
            sample(select_G4_Q$best.submodel$theta_tilde[,'par2'],
                   count[5], replace = T)
          }else{
            sample(select_G4_Q$best.submodel$theta_tilde[,'par2[1]'],
                   count[5], replace = T)
          }
        },
        if(prior.weights[6] > 0 & !is.na(best.sub.which[6])){
          if(best.sub.which[6] %in% c('all', 'BMD_d')){
            sample(select_QE4_Q$best.submodel$theta_tilde[,grep("par2\\[",colnames(select_QE4_Q$best.submodel$theta_tilde))][,i],
                   count[6], replace = T)
          }else if(best.sub.which[6] == 'none'){
            sample(select_QE4_Q$best.submodel$theta_tilde[,'par2'],
                   count[6], replace = T)
          }else{
            sample(select_QE4_Q$best.submodel$theta_tilde[,'par2[1]'],
                   count[6], replace = T)
          }
        },
        if(prior.weights[7] > 0 & !is.na(best.sub.which[7])){
          if(best.sub.which[7] %in% c('all', 'BMD_d')){
            sample(select_P4_Q$best.submodel$theta_tilde[,grep("par2\\[",colnames(select_P4_Q$best.submodel$theta_tilde))][,i],
                   count[7], replace = T)
          }else if(best.sub.which[7] == 'none'){
            sample(select_P4_Q$best.submodel$theta_tilde[,'par2'],
                   count[7], replace = T)
          }else{
            sample(select_P4_Q$best.submodel$theta_tilde[,'par2[1]'],
                   count[7], replace = T)
          }
        },
        if(prior.weights[8] > 0 & !is.na(best.sub.which[8])){
          if(best.sub.which[8] %in% c('all', 'BMD_d')){
            sample(select_L4_Q$best.submodel$theta_tilde[,grep("par2\\[",colnames(select_L4_Q$best.submodel$theta_tilde))][,i],
                   count[8], replace = T)
          }else if(best.sub.which[8] == 'none'){
            sample(select_L4_Q$best.submodel$theta_tilde[,'par2'],
                   count[8], replace = T)
          }else{
            sample(select_L4_Q$best.submodel$theta_tilde[,'par2[1]'],
                   count[8], replace = T)
          }
        }

      )

    }

    maci <- matrix(NA, max.levels, 3)
    for(i in 1:max.levels){
      maci[i, ] <- quantile(mabmd[[i]], pvec, na.rm = T)*data_all$data$maxD
    }
    colnames(maci)=c("BMDL","BMD","BMDU")
    rownames(maci)=data_all$data$covariate

    BMDq <- list()
    for(i in 1:max.levels){
      BMDq[[i]] = quantile(mabmd[[i]], seq(0,1,0.005), na.rm = T)*data_all$data$maxD ## original scale
    }

    model = c("E4_Q","IE4_Q","H4_Q","LN4_Q","G4_Q","QE4_Q","P4_Q","L4_Q","MA")
    model = as.factor(model)
    weight = c(lpw[1], lpw[2], lpw[3], lpw[4], lpw[5], lpw[6], lpw[7], lpw[8], 1)

    names(lpw) = model[1:8]


  }else{
    ## else only one BMD

    mabmd <- c(

      if(prior.weights[1] > 0 & !is.na(best.sub.which[1])){
        if(best.sub.which[1] == 'background'){
          sample(select_E4_Q$best.submodel$theta_tilde[,'par2[1]'], count[1], replace = T)
        }else if(best.sub.which[1] == 'none'){
          sample(select_E4_Q$best.submodel$theta_tilde[,'par2'], count[1], replace = T)
        }
      },
      if(prior.weights[2] > 0 & !is.na(best.sub.which[2])){
        if(best.sub.which[2] == 'background'){
          sample(select_IE4_Q$best.submodel$theta_tilde[,'par2[1]'], count[2], replace = T)
        }else if(best.sub.which[2] == 'none'){
          sample(select_IE4_Q$best.submodel$theta_tilde[,'par2'], count[2], replace = T)
        }
      },
      if(prior.weights[3] > 0 & !is.na(best.sub.which[3])){
        if(best.sub.which[3] == 'background'){
          sample(select_H4_Q$best.submodel$theta_tilde[,'par2[1]'], count[3], replace = T)
        }else if(best.sub.which[3] == 'none'){
          sample(select_H4_Q$best.submodel$theta_tilde[,'par2'], count[3], replace = T)
        }
      },
      if(prior.weights[4] > 0 & !is.na(best.sub.which[4])){
        if(best.sub.which[4] == 'background'){
          sample(select_LN4_Q$best.submodel$theta_tilde[,'par2[1]'], count[4], replace = T)
        }else if(best.sub.which[4] == 'none'){
          sample(select_LN4_Q$best.submodel$theta_tilde[,'par2'], count[4], replace = T)
        }
      },
      if(prior.weights[5] > 0 & !is.na(best.sub.which[5])){
        if(best.sub.which[5] == 'background'){
          sample(select_G4_Q$best.submodel$theta_tilde[,'par2[1]'], count[5], replace = T)
        }else if(best.sub.which[5] == 'none'){
          sample(select_G4_Q$best.submodel$theta_tilde[,'par2'], count[5], replace = T)
        }
      },
      if(prior.weights[6] > 0 & !is.na(best.sub.which[6])){
        if(best.sub.which[6] == 'background'){
          sample(select_QE4_Q$best.submodel$theta_tilde[,'par2[1]'], count[6], replace = T)
        }else if(best.sub.which[6] == 'none'){
          sample(select_QE4_Q$best.submodel$theta_tilde[,'par2'], count[6], replace = T)
        }
      },
      if(prior.weights[7] > 0 & !is.na(best.sub.which[7])){
        if(best.sub.which[7] == 'background'){
          sample(select_P4_Q$best.submodel$theta_tilde[,'par2[1]'], count[7], replace = T)
        }else if(best.sub.which[7] == 'none'){
          sample(select_P4_Q$best.submodel$theta_tilde[,'par2'], count[7], replace = T)
        }
      },
      if(prior.weights[8] > 0 & !is.na(best.sub.which[8])){
        if(best.sub.which[8] == 'background'){
          sample(select_L4_Q$best.submodel$theta_tilde[,'par2[1]'], count[8], replace = T)
        }else if(best.sub.which[8] == 'none'){
          sample(select_L4_Q$best.submodel$theta_tilde[,'par2'], count[8], replace = T)
        }
      }

    )
    maci=quantile(mabmd,pvec, na.rm = T)*data_all$data$maxD ## original scale
    # names(maci)=c("BMDL","BMD","BMDU")
    maci=t(as.matrix(maci))
    colnames(maci)=c("BMDL","BMD","BMDU")
    rownames(maci)=""

    BMDq = quantile(mabmd, seq(0,1,0.005), na.rm = T)*data_all$data$maxD ## original scale

    # BMDL = c(BMDL, maci[1]); BMD = c(BMD, maci[2]); BMDU = c(BMDU, maci[3])

    # names(BMDL) <- c("E4_N","IE4_N","H4_N","LN4_N","G4_N","QE4_N","P4_N","L4_N","E4_LN","IE4_LN","H4_LN",
    #                  "LN4_LN","G4_LN","QE4_LN","P4_LN","L4_LN","MA")
    model = c("E4_Q","IE4_Q","H4_Q","LN4_Q","G4_Q","QE4_Q","P4_Q","L4_Q","MA")
    model = as.factor(model)
    weight = c(lpw[1], lpw[2], lpw[3], lpw[4], lpw[5], lpw[6], lpw[7], lpw[8], 1)

    names(lpw) = model[1:8]

  }

  if('all' %in% best.sub.which | 'BMD_d' %in% best.sub.which){

    BMDL <- c()
    BMD <- c()
    BMDU <- c()

    for(m in model[1:8]){

      if(!TRUE %in% is.na(get(paste0('BMD_', m)))){

        # if(dim(get(paste0('BMD_', m)))[1] == data_all$data$nlevels){
        if(!is.vector(get(paste0('BMD_',m)))){

          for(i in 1:data_all$data$nlevels){
            BMDL <- c(BMDL, get(paste0('BMD_',m))[i, 1]*data_all$data$maxD)
            BMD <- c(BMD, get(paste0('BMD_',m))[i, 2]*data_all$data$maxD)
            BMDU <- c(BMDU, get(paste0('BMD_',m))[i, 3]*data_all$data$maxD)
          }
        }else{
          BMDL <- c(BMDL, rep(get(paste0('BMD_',m))[1]*data_all$data$maxD, each = data_all$data$nlevels))
          BMD <- c(BMD, rep(get(paste0('BMD_',m))[2]*data_all$data$maxD, each = data_all$data$nlevels))
          BMDU <- c(BMDU, rep(get(paste0('BMD_',m))[3]*data_all$data$maxD, each = data_all$data$nlevels))
        }
      }else{
        BMDL <- c(BMDL, rep(NA, data_all$data$nlevels))
        BMD <- c(BMD, rep(NA, data_all$data$nlevels))
        BMDU <- c(BMDU, rep(NA, data_all$data$nlevels))
      }
    }

    model.weight <- data.frame(Model = rep(model[1:8], each = data_all$data$nlevels),
                               Weight = rep(round(lpw[1:8],4), each = data_all$data$nlevels),
                               Submodel = rep(best.sub.which, each = data_all$data$nlevels),
                               Submodel.weight = rep(best.sub.weight, each = data_all$data$nlevels),
                               Covariate = rep(data_all$data$covariate, 8),
                               BMDL = BMDL, BMD = BMD, BMDU = BMDU
                               # LogLik = best.sub.loglik,
    )

  }else{

    BMDL <- c()
    BMD <- c()
    BMDU <- c()

    for(m in model[1:8]){
      if(!is.na(get(paste0('BMD_', m))[1])){
        BMDL <- c(BMDL, get(paste0('BMD_',m))[1]*data_all$data$maxD)
        BMD <- c(BMD, get(paste0('BMD_',m))[2]*data_all$data$maxD)
        BMDU <- c(BMDU, get(paste0('BMD_',m))[3]*data_all$data$maxD)

      }else{
        BMDL <- c(BMDL, NA)
        BMD <- c(BMD, NA)
        BMDU <- c(BMDU, NA)
      }
    }

    model.weight <- data.frame(Model = model[1:8],
                               Weight = round(lpw[1:8],4),
                               Submodel = best.sub.which,
                               Submodel.weight = best.sub.weight,
                               Covariate = 'None',
                               BMDL = BMDL, BMD = BMD, BMDU = BMDU
                               # LogLik = best.sub.loglik,
    )

  }

  if(exists('select_E4_Q')) { parE4_Q = select_E4_Q$best.submodel$par} else{parE4_Q = NULL}
  if(exists('select_IE4_Q')) { parIE4_Q = select_IE4_Q$best.submodel$par} else{parIE4_Q = NULL}
  if(exists('select_H4_Q')) { parH4_Q = select_H4_Q$best.submodel$par} else{parH4_Q = NULL}
  if(exists('select_LN4_Q')) { parLN4_Q = select_LN4_Q$best.submodel$par} else{parLN4_Q = NULL}
  if(exists('select_G4_Q')) { parG4_Q = select_G4_Q$best.submodel$par} else{parG4_Q = NULL}
  if(exists('select_QE4_Q')) { parQE4_Q = select_QE4_Q$best.submodel$par} else{parQE4_Q = NULL}
  if(exists('select_P4_Q')) { parP4_Q = select_P4_Q$best.submodel$par} else{parP4_Q = NULL}
  if(exists('select_L4_Q')) { parL4_Q = select_L4_Q$best.submodel$par} else{parL4_Q = NULL}


  return(list(MA = maci,
              weights = lpw,
              # best.submodels = best.sub.which,
              # data = data.frame(x = data_all$data$x,
              #                   y = data_all$data$y,
              #                   n = data_all$data$n,
              #                   cov = data_all$data$covariate
              # ),
              data = data_all$data$org.data,
              summary = model.weight,
              maxDose = data_all$data$maxD,
              BMD_Mixture = BMDq,
              parE4_Q = parE4_Q,
              parIE4_Q = parIE4_Q,
              parH4_Q = parH4_Q,
              parLN4_Q = parLN4_Q,
              parG4_Q = parG4_Q,
              parQE4_Q = parQE4_Q,
              parP4_Q = parP4_Q,
              parL4_Q = parL4_Q,
              q = data_all$data$q
  ))

}

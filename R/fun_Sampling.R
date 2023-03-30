
#' Perform model averaging using MCMC methods (Bridge sampling & Partial Laplace)
#'
#' This method assumed data for continuous endpoints.
#'
#' More detailed descriprion
#' @param data.N the input data as returned by function PREP_DATA_N or PREP_DATA_N_C for clustered data
#' @param data.LN the input data as returned by function PREP_DATA_LN or PREP_DATA_LN_C for clustered data
#' @param data.Q the input data as returned by function PREP_DATA_Q
#' @param prior.weights a vector specifying which of the 16 (continuous) or 8 (quantal) models should be included (1 = include, 0 = exclude)
#' @param ndraws the number of draws, default 30000
#' @param nrchains the number of chains to be used in the MCMC, default 3
#' @param nriterations the number of iterations per chain, default 3000
#' @param warmup the number of iterations per chain to be discarded as burnin, default 1000
#' @param delta default 0.8
#' @param treedepth default 10
#' @param seed random seed, default 123
#' @param pvec vector specifying the three BMD quantiles of interest (default 90% CrI)
#'
#' @description Using MCMC, we compute the parameters of each model, perform model averaging using bridge sampling.
#'              We also implemented a Hybrid Laplace method within this function where the model parameters are estimated
#'              with MCMC, but the model weights are computed using Laplace approximation.
#'              By default, all 16 models are included for continuous data, and all 8 models for quantal data.
#'              Models can be excluded by setting their respective weight to 0 in \code{prior.weights}. The order of models fitted can be obtained using the \code{get_models()} function.
#'
#' `sampling_MA` is used for continuous data
#'
#' `sampling_MAc` is used for clustered continuous data
#'
#' `samplingQ_MA` is used for quantal data
#'
#' @importFrom truncnorm dtruncnorm
#'
#' @examples
#'  data_N <- PREP_DATA_N(data = as.data.frame(immunotoxicityData[1:5,]),
#'                        sumstats = TRUE, sd = TRUE, q = 0.1)
#'  data_LN <- PREP_DATA_LN(data = as.data.frame(immunotoxicityData[1:5,]),
#'                          sumstats = TRUE, sd = TRUE, q = 0.1) #'
#'  SBMD <- sampling_MA(data_N, data_LN, prior.weights = c(rep(1,4), rep(0,12)))
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
#' @return `MA_bridge_sampling`  Model averaged BMD credible interval based on Bridge sampling
#' @return `MA_laplace` Model averaged BMD credible interval based on Hybrid Laplace
#' @return `weights_bridge_sampling` Model weights used in the averaging for Bridge sampling
#' @return `weights_laplace` Model weights used in the averaging for Hybrid Laplace
#' @return `convergence` vector indicating whether the models have converged (1) or not (0)
#' @return `divergences` vector containing the proportion of divergent transitions for each model
#' @return `bs_weights_conv` model weights based on bridge sampling including converged models only
#' @return `ls_weights_conv` model weights based on hybrid laplace including converged models only
#' @return `MA_bs_conv` model averaged BMD credible interval based on bridge sampling including converged models only
#' @return `MA_ls_conv` model averaged BMD credible interval based on hybrid laplace including converged models only
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
#' @return `models_included_laplace` vector containing the names of models included in the model averaging using Laplace
#' @return `models_included_bridge` vector containing the names of models included in the model averaging using Bridge sampling
#' @return `q` BMR
#' @return `max.dose` maximum dose level (original scale)
#' @return `dataN` normal summary data used for analysis
#' @return `dataLN` lognormal summary data used for analysis
#' @return `BMDMixture` vector of length \code{ndraws} containing the draws from the model-averaged posterior based on hybrid laplace
#' @return `BMDMixture.conv` vector of length \code{ndraws} containing the draws from the model-averaged posterior based on hybrid laplace, including converged models only
#' @return `BMDMixtureBS` vector of length \code{ndraws} containing the draws from the model-averaged posterior based on bridge sampling
#' @return `BMDMixture.convBS` vector of length \code{ndraws} containing the draws from the model-averaged posterior based on bridge sampling, including converged models only
#' @return `MA_dr_bs` vector containing the model-averaged response at each dose level, for bridge sampling
#' @return `MA_dr_ls` vector containing the model-averaged response at each dose level, for hybrid laplace
#' @return `MA_dr_bs_conv` vector containing the model-averaged response at each dose level, for bridge sampling including converged models only
#' @return `MA_dr_ls_conv` vector containing the model-averaged response at each dose level, for hybrid laplace including converged models only
#' @return `MA_post_bs` vector containing the 0.5%-percentiles of the model-averaged posterior, for bridge sampling
#' @return `MA_post_ls` vector containing the 0.5%-percentiles of the model-averaged posterior, for hybrid laplace
#' @return `MA_post_bs_conv` vector containing the 0.5%-percentiles of the model-averaged posterior, for bridge sampling including converged models only
#' @return `MA_post_ls_conv` vector containing the 0.5%-percentiles of the model-averaged posterior, for hybrid laplace including converged models only
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
#' @export sampling_MA
#'
sampling_MA=function(data.N,data.LN,prior.weights = rep(1,16),
                     ndraws = 30000,nrchains=3,
                     nriterations=3000,warmup=1000,
                     delta=0.8,treedepth=10,seed=123,pvec=c(0.05,0.5,0.95)){


  # if(data.N$increasing == TRUE){

  prior.weights.orig = prior.weights

  ## Normal distribution

  data=data.N$data
  start=data.N$start
  startQ=data.N$startQ

  BMDL=c(); BMD=c(); BMDU=c()
  converged=c()

  nModels = 16

  ## Obtain model parameters via MCMC sampling

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[1]))
  }
  if(prior.weights[1]>0){

    fitstanE4_N = fun_sampling(stanmodels$mE4, data, start,
                               ndraws,nrchains,
                               nriterations,warmup,
                               delta,treedepth,seed,pvec)

    if(is.null(fitstanE4_N)){
      prior.weights[1] <- 0
      warning('difficulties fitting the Exponential (Normal) model using sampling; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsE4N <- par_extract(fitstanE4_N, model_name = "E4_N")
      # parsE4N[,"BMD"] <- parsE4N[,"BMD"]*data$maxD # BMD on original scale

      #diagnostics here
      # posterior_diag(model_stan = fitstanE4_N)
      convergence_stat <- convergence_deci(model_stan = fitstanE4_N, nrchains = nrchains)
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)
      # conv_a_r = 0; if(convergence_stat['a', 'Rhat_Convergence'] == "Convergence") conv_a_r = 1

      # divergent transitions
      div_E4_N <- sum(sapply(rstan::get_sampler_params(fitstanE4_N, inc_warmup = F),
                             function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))

      E4resNI=quantile(as.matrix(fitstanE4_N)[,"par2"]*data$maxD,pvec)
      E4resNI=c(E4resNI,apply(as.matrix(fitstanE4_N),2,median)[c("par1","par2","par3","par4","par5")])
      names(E4resNI)=c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t")

      E4outNI <- outLP(parsE4N, pvec, data$maxD)

      if(data$data_type == 1){
        DRM_E4_N = DRM.E4_NI(E4resNI[4:7], data$x, data$q)
      }else if(data$data_type == 3){
        DRM_E4_N = DRM.E4_ND(E4resNI[4:7], data$x, data$q)
      }


      # Covariance between b-d and between BMD-d
      E4covNI = c(cov(as.matrix(fitstanE4_N)[,c("b","d")], use="complete.obs")["b","d"],
                  cov(as.matrix(fitstanE4_N)[,c("par2","d")], use="complete.obs")["par2","d"])

      E4corrNI = c(cor(as.matrix(fitstanE4_N)[,c("b","d")], use="complete.obs")["b","d"],
                   cor(as.matrix(fitstanE4_N)[,c("par2","d")], use="complete.obs")["par2","d"])

      BMDL=c(BMDL, E4resNI[1]); BMD=c(BMD, E4resNI[2]); BMDU=c(BMDU, E4resNI[3])
      bridgeE4N = bridgesampling::bridge_sampler(fitstanE4_N, silent = T) # compute log marginal likelihood
    }
  }
  if(prior.weights[1] == 0){E4resNI=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeE4N=NA; converged=c(converged, NA); E4covNI=rep(NA,2); E4corrNI=rep(NA,2); DRM_E4_N=rep(NA,length(data$x))
  parsE4N <- NA
  div_E4_N <- NA
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

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[2]))
  }
  if(prior.weights[2]>0){

    fitstanIE4_N = fun_sampling(stanmodels$mIE4, data, start,
                                ndraws,nrchains,
                                nriterations,warmup,
                                delta,treedepth,seed,pvec)

    if(is.null(fitstanIE4_N)){
      prior.weights[2] <- 0
      warning('difficulties fitting the Inverse Exponential (Normal) model using sampling; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsIE4N <- par_extract(fitstanIE4_N, model_name = "IE4_N")
      # parsIE4N[,"BMD"] <- parsIE4N[,"BMD"]*data$maxD # BMD on original scale

      #diagnostics here
      # posterior_diag(model_stan = fitstanIE4_N)
      convergence_stat <- convergence_deci(model_stan = fitstanIE4_N, nrchains = nrchains)
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)

      # divergent transitions
      div_IE4_N <- sum(sapply(rstan::get_sampler_params(fitstanIE4_N, inc_warmup = F),
                              function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))

      IE4resNI=(quantile(as.matrix(fitstanIE4_N)[,"par2"],pvec))*data$maxD
      IE4resNI=c(IE4resNI,apply(as.matrix(fitstanIE4_N),2,median)[c("par1","par2","par3","par4","par5")])
      names(IE4resNI)=c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t")

      IE4outNI <- outLP(parsIE4N, pvec, data$maxD)

      if(data$data_type == 1){
        DRM_IE4_N = DRM.IE4_NI(IE4resNI[4:7], data$x, data$q)
      }else if(data$data_type == 3){
        DRM_IE4_N = DRM.IE4_ND(IE4resNI[4:7], data$x, data$q)
      }

      # Covariance between b-d and between BMD-d
      IE4covNI = c(cov(as.matrix(fitstanIE4_N)[,c("b","d")], use="complete.obs")["b","d"],
                   cov(as.matrix(fitstanIE4_N)[,c("par2","d")], use="complete.obs")["par2","d"])

      IE4corrNI = c(cor(as.matrix(fitstanIE4_N)[,c("b","d")], use="complete.obs")["b","d"],
                    cor(as.matrix(fitstanIE4_N)[,c("par2","d")], use="complete.obs")["par2","d"])

      BMDL=c(BMDL, IE4resNI[1]); BMD=c(BMD, IE4resNI[2]); BMDU=c(BMDU, IE4resNI[3])
      bridgeIE4N = bridgesampling::bridge_sampler(fitstanIE4_N, silent = T)

    }
  }
  if(prior.weights[2] == 0){IE4resNI=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeIE4N=NA; converged=c(converged, NA); IE4covNI=rep(NA,2); IE4corrNI=rep(NA,2); DRM_IE4_N=rep(NA,length(data$x))
  parsIE4N <- NA
  div_IE4_N <- NA

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

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[3]))
  }
  if(prior.weights[3]>0){

    fitstanH4_N = fun_sampling(stanmodels$mH4, data, start,
                               ndraws,nrchains,
                               nriterations,warmup,
                               delta,treedepth,seed,pvec)

    if(is.null(fitstanH4_N)){
      prior.weights[3] <- 0
      warning('difficulties fitting the Hill (Normal) model using sampling; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsH4N <- par_extract(fitstanH4_N, model_name = "H4_N")
      # parsH4N[,"BMD"] <- parsH4N[,"BMD"]*data$maxD # BMD on original scale

      #diagnostics here
      # posterior_diag(model_stan = fitstanH4_N)
      convergence_stat <- convergence_deci(model_stan = fitstanH4_N, nrchains = nrchains)
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)

      # divergent transitions
      div_H4_N <- sum(sapply(rstan::get_sampler_params(fitstanH4_N, inc_warmup = F),
                             function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))

      H4resNI=(quantile(as.matrix(fitstanH4_N)[,"par2"],pvec))*data$maxD
      H4resNI=c(H4resNI,apply(as.matrix(fitstanH4_N),2,median)[c("par1","par2","par3","par4","par5")])
      names(H4resNI)=c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t")

      H4outNI <- outLP(parsH4N, pvec, data$maxD)

      if(data$data_type == 1){
        DRM_H4_N = DRM.H4_NI(H4resNI[4:7], data$x, data$q)
      }else if(data$data_type == 3){
        DRM_H4_N = DRM.H4_ND(H4resNI[4:7], data$x, data$q)
      }
      # Covariance between b-d and between BMD-d
      H4covNI = c(cov(as.matrix(fitstanH4_N)[,c("b","d")], use="complete.obs")["b","d"],
                  cov(as.matrix(fitstanH4_N)[,c("par2","d")], use="complete.obs")["par2","d"])

      H4corrNI = c(cor(as.matrix(fitstanH4_N)[,c("b","d")], use="complete.obs")["b","d"],
                   cor(as.matrix(fitstanH4_N)[,c("par2","d")], use="complete.obs")["par2","d"])

      BMDL=c(BMDL, H4resNI[1]); BMD=c(BMD, H4resNI[2]); BMDU=c(BMDU, H4resNI[3])
      bridgeH4N = bridgesampling::bridge_sampler(fitstanH4_N, silent = T)
    }
  }
  if(prior.weights[3] == 0){H4resNI=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeH4N=NA;
  converged=c(converged, NA); H4covNI=rep(NA,2); H4corrNI=rep(NA,2); DRM_H4_N=rep(NA,length(data$x))
  parsH4N <- NA
  div_H4_N <- NA

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

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[4]))
  }
  if(prior.weights[4]>0){
    fitstanLN4_N = fun_sampling(stanmodels$mLN4, data, start,
                                ndraws,nrchains,
                                nriterations,warmup,
                                delta,treedepth,seed,pvec)

    if(is.null(fitstanLN4_N)){
      prior.weights[4] <- 0
      warning('difficulties fitting the Lognormal (Normal) model using sampling; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsLN4N <- par_extract(fitstanLN4_N, model_name = "LN4_N")
      # parsLN4N[,"BMD"] <- parsLN4N[,"BMD"]*data$maxD # BMD on original scale

      #diagnostics here
      # posterior_diag(model_stan = fitstanLN4_N)
      convergence_stat <- convergence_deci(model_stan = fitstanLN4_N, nrchains = nrchains)
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)

      # divergent transitions
      div_LN4_N <- sum(sapply(rstan::get_sampler_params(fitstanLN4_N, inc_warmup = F),
                              function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))

      LN4resNI=(quantile(as.matrix(fitstanLN4_N)[,"par2"],pvec))*data$maxD
      LN4resNI=c(LN4resNI,apply(as.matrix(fitstanLN4_N),2,median)[c("par1","par2","par3","par4","par5")])
      names(LN4resNI)=c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t")

      LN4outNI <- outLP(parsLN4N, pvec, data$maxD)

      if(data$data_type == 1){
        DRM_LN4_N = DRM.LN4_NI(LN4resNI[4:7], data$x, data$q)
      }else if(data$data_type == 3){
        DRM_LN4_N = DRM.LN4_ND(LN4resNI[4:7], data$x, data$q)
      }

      # Covariance between b-d and between BMD-d
      LN4covNI = c(cov(as.matrix(fitstanLN4_N)[,c("b","d")], use="complete.obs")["b","d"],
                   cov(as.matrix(fitstanLN4_N)[,c("par2","d")], use="complete.obs")["par2","d"])

      LN4corrNI = c(cor(as.matrix(fitstanLN4_N)[,c("b","d")], use="complete.obs")["b","d"],
                    cor(as.matrix(fitstanLN4_N)[,c("par2","d")], use="complete.obs")["par2","d"])

      BMDL=c(BMDL, LN4resNI[1]); BMD=c(BMD, LN4resNI[2]); BMDU=c(BMDU, LN4resNI[3])
      bridgeLN4N = bridgesampling::bridge_sampler(fitstanLN4_N, silent = T)
    }
  }
  if(prior.weights[4] == 0){LN4resNI=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeLN4N=NA; converged=c(converged, NA); LN4covNI=rep(NA,2); LN4corrNI=rep(NA,2); DRM_LN4_N=rep(NA,length(data$x))
  parsLN4N <- NA
  div_LN4_N <- NA

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

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[5]))
  }
  if(prior.weights[5]>0){

    fitstanG4_N = fun_sampling(stanmodels$mG4, data, start,
                               ndraws,nrchains,
                               nriterations,warmup,
                               delta,treedepth,seed,pvec)

    if(is.null(fitstanG4_N)){
      prior.weights[5] <- 0
      warning('difficulties fitting the Gamma (Normal) model using sampling; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsG4N <- par_extract(fitstanG4_N, model_name = "G4_N")
      # parsG4N[,"BMD"] <- parsG4N[,"BMD"]*data$maxD # BMD on original scale

      #diagnostics here
      # posterior_diag(model_stan = fitstanG4_N)
      convergence_stat <- convergence_deci(model_stan = fitstanG4_N, nrchains = nrchains)
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)

      # divergent transitions
      div_G4_N <- sum(sapply(rstan::get_sampler_params(fitstanG4_N, inc_warmup = F),
                             function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))

      G4resNI=(quantile(as.matrix(fitstanG4_N)[,"par2"],pvec))*data$maxD
      G4resNI=c(G4resNI,apply(as.matrix(fitstanG4_N),2,median)[c("par1","par2","par3","par4","par5")])
      names(G4resNI)=c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t")

      G4outNI <- outLP(parsG4N, pvec, data$maxD)

      if(data$data_type == 1){
        DRM_G4_N = DRM.G4_NI(G4resNI[4:7], data$x, data$q)
      }else if(data$data_type == 3){
        DRM_G4_N = DRM.G4_ND(G4resNI[4:7], data$x, data$q)
      }

      # Covariance between b-d and between BMD-d
      G4covNI = c(cov(as.matrix(fitstanG4_N)[,c("b","d")], use="complete.obs")["b","d"],
                  cov(as.matrix(fitstanG4_N)[,c("par2","d")], use="complete.obs")["par2","d"])

      G4corrNI = c(cor(as.matrix(fitstanG4_N)[,c("b","d")], use="complete.obs")["b","d"],
                   cor(as.matrix(fitstanG4_N)[,c("par2","d")], use="complete.obs")["par2","d"])

      BMDL=c(BMDL, G4resNI[1]); BMD=c(BMD, G4resNI[2]); BMDU=c(BMDU, G4resNI[3])
      bridgeG4N = bridgesampling::bridge_sampler(fitstanG4_N, silent = T)
    }
  }
  if(prior.weights[5] == 0){G4resNI=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeG4N=NA; converged=c(converged, NA); G4covNI=rep(NA,2); G4corrNI=rep(NA,2); DRM_G4_N=rep(NA,length(data$x))
  parsG4N <- NA
  div_G4_N <- NA

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

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[6]))
  }
  if(prior.weights[6]>0){
    # print(6)
    fitstanQE4_N = fun_sampling(stanmodels$mQE4, data, startQ,
                                ndraws,nrchains,
                                nriterations,warmup,
                                delta,treedepth,seed,pvec)

    if(is.null(fitstanQE4_N)){
      prior.weights[6] <- 0
      warning('difficulties fitting the Quadratic Exponential (Normal) model using sampling; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsQE4N <- par_extract(fitstanQE4_N, model_name = "QE4_N")
      # parsQE4N[,"BMD"] <- parsQE4N[,"BMD"]*data$maxD # BMD on original scale

      #diagnostics here
      # posterior_diag(model_stan = fitstanQE4_N)
      convergence_stat <- convergence_deci(model_stan = fitstanQE4_N, nrchains = nrchains)
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)

      # divergent transitions
      div_QE4_N <- sum(sapply(rstan::get_sampler_params(fitstanQE4_N, inc_warmup = F),
                              function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))

      QE4resNI=(quantile(as.matrix(fitstanQE4_N)[,"par2"],pvec))*data$maxD
      QE4resNI=c(QE4resNI,apply(as.matrix(fitstanQE4_N),2,median)[c("par1","par2","par3","par4","par5")])
      names(QE4resNI)=c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t")

      QE4outNI <- outLP(parsQE4N, pvec, data$maxD)

      if(data$data_type == 1){
        DRM_QE4_N = DRM.QE4_NI(QE4resNI[4:7], data$x, data$q)
      }else if(data$data_type == 3){
        DRM_QE4_N = DRM.QE4_ND(QE4resNI[4:7], data$x, data$q)
      }

      # Covariance between b-d and between BMD-d
      QE4covNI = c(cov(as.matrix(fitstanQE4_N)[,c("b","d")], use="complete.obs")["b","d"],
                   cov(as.matrix(fitstanQE4_N)[,c("par2","d")], use="complete.obs")["par2","d"])

      QE4corrNI = c(cor(as.matrix(fitstanQE4_N)[,c("b","d")], use="complete.obs")["b","d"],
                    cor(as.matrix(fitstanQE4_N)[,c("par2","d")], use="complete.obs")["par2","d"])

      BMDL=c(BMDL, QE4resNI[1]); BMD=c(BMD, QE4resNI[2]); BMDU=c(BMDU, QE4resNI[3])
      bridgeQE4N = bridgesampling::bridge_sampler(fitstanQE4_N, silent = T)
    }
  }
  if(prior.weights[6] == 0){QE4resNI=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeQE4N=NA; converged=c(converged, NA); QE4covNI=rep(NA,2); QE4corrNI=rep(NA,2); DRM_QE4_N=rep(NA,length(data$x))
  parsQE4N <- NA
  div_QE4_N <- NA

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

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[7]))
  }
  if(prior.weights[7]>0){
    # print(7)
    fitstanP4_N = fun_sampling(stanmodels$mP4, data, start,
                               ndraws,nrchains,
                               nriterations,warmup,
                               delta,treedepth,seed,pvec)

    if(is.null(fitstanP4_N)){
      prior.weights[7] <- 0
      warning('difficulties fitting the Probit (Normal) model using sampling; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsP4N <- par_extract(fitstanP4_N, model_name = "P4_N")
      # parsP4N[,"BMD"] <- parsP4N[,"BMD"]*data$maxD # BMD on original scale

      #diagnostics here
      # posterior_diag(model_stan = fitstanP4_N)
      convergence_stat <- convergence_deci(model_stan = fitstanP4_N, nrchains = nrchains)
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)

      # divergent transitions
      div_P4_N <- sum(sapply(rstan::get_sampler_params(fitstanP4_N, inc_warmup = F),
                             function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))

      P4resNI=(quantile(as.matrix(fitstanP4_N)[,"par2"],pvec))*data$maxD
      P4resNI=c(P4resNI,apply(as.matrix(fitstanP4_N),2,median)[c("par1","par2","par3","par4","par5")])
      names(P4resNI)=c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t")

      P4outNI <- outLP(parsP4N, pvec, data$maxD)

      if(data$data_type == 1){
        DRM_P4_N = DRM.P4_NI(P4resNI[4:7], data$x, data$q)
      }else if(data$data_type == 3){
        DRM_P4_N = DRM.P4_ND(P4resNI[4:7], data$x, data$q)
      }

      # Covariance between b-d and between BMD-d
      P4covNI = c(cov(as.matrix(fitstanP4_N)[,c("b","d")], use="complete.obs")["b","d"],
                  cov(as.matrix(fitstanP4_N)[,c("par2","d")], use="complete.obs")["par2","d"])

      P4corrNI = c(cor(as.matrix(fitstanP4_N)[,c("b","d")], use="complete.obs")["b","d"],
                   cor(as.matrix(fitstanP4_N)[,c("par2","d")], use="complete.obs")["par2","d"])

      BMDL=c(BMDL, P4resNI[1]); BMD=c(BMD, P4resNI[2]); BMDU=c(BMDU, P4resNI[3])
      bridgeP4N = bridgesampling::bridge_sampler(fitstanP4_N, silent = T)
    }
  }
  if(prior.weights[7] == 0){P4resNI=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeP4N=NA; converged=c(converged, NA);
  P4covNI=rep(NA,2); P4corrNI=rep(NA,2);
  DRM_P4_N=rep(NA,length(data$x))
  parsP4N <- NA
  div_P4_N <- NA

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

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[8]))
  }
  if(prior.weights[8]>0){
    # print(8)
    fitstanL4_N = fun_sampling(stanmodels$mL4, data, start,
                               ndraws,nrchains,
                               nriterations,warmup,
                               delta,treedepth,seed,pvec)

    if(is.null(fitstanL4_N)){
      prior.weights[8] <- 0
      warning('difficulties fitting the Logit (Normal) model using sampling; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsL4N <- par_extract(fitstanL4_N, model_name = "L4_N")
      # parsL4N[,"BMD"] <- parsL4N[,"BMD"]*data$maxD # BMD on original scale

      # save(fitstanL4_N, file = "fitL4N.RData")

      #diagnostics here
      # posterior_diag(model_stan = fitstanL4_N)
      convergence_stat <- convergence_deci(model_stan = fitstanL4_N, nrchains = nrchains)
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)

      # divergent transitions
      div_L4_N <- sum(sapply(rstan::get_sampler_params(fitstanL4_N, inc_warmup = F),
                             function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))

      L4resNI=(quantile(as.matrix(fitstanL4_N)[,"par2"],pvec))*data$maxD
      L4resNI=c(L4resNI,apply(as.matrix(fitstanL4_N),2,median)[c("par1","par2","par3","par4","par5")])
      names(L4resNI)=c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t")

      L4outNI <- outLP(parsL4N, pvec, data$maxD)


      if(data$data_type == 1){
        DRM_L4_N = DRM.L4_NI(L4resNI[4:7], data$x, data$q)
      }else if(data$data_type == 3){
        DRM_L4_N = DRM.L4_ND(L4resNI[4:7], data$x, data$q)
      }

      # Covariance between b-d and between BMD-d
      L4covNI = c(cov(as.matrix(fitstanL4_N)[,c("b","d")], use="complete.obs")["b","d"],
                  cov(as.matrix(fitstanL4_N)[,c("par2","d")], use="complete.obs")["par2","d"])

      L4corrNI = c(cor(as.matrix(fitstanL4_N)[,c("b","d")], use="complete.obs")["b","d"],
                   cor(as.matrix(fitstanL4_N)[,c("par2","d")], use="complete.obs")["par2","d"])

      BMDL=c(BMDL, L4resNI[1]); BMD=c(BMD, L4resNI[2]); BMDU=c(BMDU, L4resNI[3])
      bridgeL4N = bridgesampling::bridge_sampler(fitstanL4_N, silent = T)
    }
  }
  if(prior.weights[8] == 0){L4resNI=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeL4N=NA; converged=c(converged, NA); L4covNI=rep(NA,2); L4corrNI=rep(NA,2); DRM_L4_N=rep(NA,length(data$x))
  parsL4N <- NA
  div_L4_N <- NA

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

  ## Lognormal distribution

  data=data.LN$data
  start=data.LN$start
  startQ=data.LN$startQ

  # model specific results
  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[9]))
  }
  if(prior.weights[9]>0){
    # print(9)

    fitstanE4_LN = fun_sampling(stanmodels$mE4, data, start,
                                ndraws,nrchains,
                                nriterations,warmup,
                                delta,treedepth,seed,pvec)

    if(is.null(fitstanE4_LN)){
      prior.weights[9] <- 0
      warning('difficulties fitting the Exponential (Lognormal) model using sampling; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsE4LN <- par_extract(fitstanE4_LN, model_name = "E4_LN")
      # parsE4LN[,"BMD"] <- parsE4LN[,"BMD"]*data$maxD # BMD on original scale

      #diagnostics here
      # posterior_diag(model_stan = fitstanE4_LN)
      convergence_stat <- convergence_deci(model_stan = fitstanE4_LN, nrchains = nrchains)
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)

      # divergent transitions
      div_E4_LN <- sum(sapply(rstan::get_sampler_params(fitstanE4_LN, inc_warmup = F),
                              function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))

      E4resLNI=(quantile(as.matrix(fitstanE4_LN)[,"par2"],pvec))*data$maxD
      E4resLNI=c(E4resLNI,apply(as.matrix(fitstanE4_LN),2,median)[c("par1","par2","par3","par4","par5")])
      names(E4resLNI)=c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t")

      E4outLNI <- outLP(parsE4LN, pvec, data$maxD)

      if(data$data_type == 2){
        DRM_E4_LN = exp(DRM.E4_LNI(E4resLNI[4:7], data$x, data$q, shift=data$shift))
      }else if(data$data_type == 4){
        DRM_E4_LN = exp(DRM.E4_LND(E4resLNI[4:7], data$x, data$q, shift=data$shift))
      }

      # Covariance between b-d and between BMD-d
      E4covLNI = c(cov(as.matrix(fitstanE4_LN)[,c("b","d")], use="complete.obs")["b","d"],
                   cov(as.matrix(fitstanE4_LN)[,c("par2","d")], use="complete.obs")["par2","d"])

      E4corrLNI = c(cor(as.matrix(fitstanE4_LN)[,c("b","d")], use="complete.obs")["b","d"],
                    cor(as.matrix(fitstanE4_LN)[,c("par2","d")], use="complete.obs")["par2","d"])

      BMDL=c(BMDL, E4resLNI[1]); BMD=c(BMD, E4resLNI[2]); BMDU=c(BMDU, E4resLNI[3])
      bridgeE4LN = bridgesampling::bridge_sampler(fitstanE4_LN, silent = T)
    }
  }
  if(prior.weights[9] == 0){E4resLNI=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeE4LN=NA; converged=c(converged, NA); E4covLNI=rep(NA,2); E4corrLNI=rep(NA,2); DRM_E4_LN=rep(NA,length(data$x))
  parsE4LN <- NA
  div_E4_LN <- NA

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

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[10]))
  }
  if(prior.weights[10]>0){
    # print(10)

    fitstanIE4_LN = fun_sampling(stanmodels$mIE4, data, start,
                                 ndraws,nrchains,
                                 nriterations,warmup,
                                 delta,treedepth,seed,pvec)

    if(is.null(fitstanIE4_LN)){
      prior.weights[10] <- 0
      warning('difficulties fitting the Inverse Exponential (Lognormal) model using sampling; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsIE4LN <- par_extract(fitstanIE4_LN, model_name = "IE4_LN")
      # parsIE4LN[,"BMD"] <- parsIE4LN[,"BMD"]*data$maxD # BMD on original scale

      #diagnostics here
      # posterior_diag(model_stan = fitstanIE4_LN)
      convergence_stat <- convergence_deci(model_stan = fitstanIE4_LN, nrchains = nrchains)
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)

      # divergent transitions
      div_IE4_LN <- sum(sapply(rstan::get_sampler_params(fitstanIE4_LN, inc_warmup = F),
                               function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))

      IE4resLNI=(quantile(as.matrix(fitstanIE4_LN)[,"par2"],pvec))*data$maxD
      IE4resLNI=c(IE4resLNI,apply(as.matrix(fitstanIE4_LN),2,median)[c("par1","par2","par3","par4","par5")])
      names(IE4resLNI)=c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t")

      IE4outLNI <- outLP(parsIE4LN, pvec, data$maxD)


      if(data$data_type == 2){
        DRM_IE4_LN = exp(DRM.IE4_LNI(IE4resLNI[4:7], data$x, data$q, shift=data$shift))
      }else if(data$data_type == 4){
        DRM_IE4_LN = exp(DRM.IE4_LND(IE4resLNI[4:7], data$x, data$q, shift=data$shift))
      }

      IE4covLNI = c(cov(as.matrix(fitstanIE4_LN)[,c("b","d")], use="complete.obs")["b","d"],
                    cov(as.matrix(fitstanIE4_LN)[,c("par2","d")], use="complete.obs")["par2","d"])

      IE4corrLNI = c(cor(as.matrix(fitstanIE4_LN)[,c("b","d")], use="complete.obs")["b","d"],
                     cor(as.matrix(fitstanIE4_LN)[,c("par2","d")], use="complete.obs")["par2","d"])

      BMDL=c(BMDL, IE4resLNI[1]); BMD=c(BMD, IE4resLNI[2]); BMDU=c(BMDU, IE4resLNI[3])
      bridgeIE4LN = bridgesampling::bridge_sampler(fitstanIE4_LN, silent = T)
    }
  }
  if(prior.weights[10] == 0){IE4resLNI=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeIE4LN=NA; converged=c(converged, NA); IE4covLNI=rep(NA,2); IE4corrLNI=rep(NA,2); DRM_IE4_LN=rep(NA,length(data$x))
  parsIE4LN <- NA
  div_IE4_LN <- NA

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

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[11]))
  }
  if(prior.weights[11]>0){
    # print(11)

    fitstanH4_LN = fun_sampling(stanmodels$mH4, data, start,
                                ndraws,nrchains,
                                nriterations,warmup,
                                delta,treedepth,seed,pvec)

    if(is.null(fitstanH4_LN)){
      prior.weights[11] <- 0
      warning('difficulties fitting the Hill (Lognormal) model using sampling; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsH4LN <- par_extract(fitstanH4_LN, model_name = "H4_LN")
      # parsH4LN[,"BMD"] <- parsH4LN[,"BMD"]*data$maxD # BMD on original scale

      #diagnostics here
      # posterior_diag(model_stan = fitstanH4_LN)
      convergence_stat <- convergence_deci(model_stan = fitstanH4_LN, nrchains = nrchains)
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)

      # divergent transitions
      div_H4_LN <- sum(sapply(rstan::get_sampler_params(fitstanH4_LN, inc_warmup = F),
                              function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))

      H4resLNI=(quantile(as.matrix(fitstanH4_LN)[,"par2"],pvec))*data$maxD
      H4resLNI=c(H4resLNI,apply(as.matrix(fitstanH4_LN),2,median)[c("par1","par2","par3","par4","par5")])
      names(H4resLNI)=c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t")

      H4outLNI <- outLP(parsH4LN, pvec, data$maxD)

      if(data$data_type == 2){
        DRM_H4_LN = exp(DRM.H4_LNI(H4resLNI[4:7], data$x, data$q, shift=data$shift))
      }else if(data$data_type == 4){
        DRM_H4_LN = exp(DRM.H4_LND(H4resLNI[4:7], data$x, data$q, shift=data$shift))
      }

      H4covLNI = c(cov(as.matrix(fitstanH4_LN)[,c("b","d")], use="complete.obs")["b","d"],
                   cov(as.matrix(fitstanH4_LN)[,c("par2","d")], use="complete.obs")["par2","d"])

      H4corrLNI = c(cor(as.matrix(fitstanH4_LN)[,c("b","d")], use="complete.obs")["b","d"],
                    cor(as.matrix(fitstanH4_LN)[,c("par2","d")], use="complete.obs")["par2","d"])

      BMDL=c(BMDL, H4resLNI[1]); BMD=c(BMD, H4resLNI[2]); BMDU=c(BMDU, H4resLNI[3])
      bridgeH4LN = bridgesampling::bridge_sampler(fitstanH4_LN, silent = T)
    }
  }
  if(prior.weights[11] == 0){H4resLNI=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeH4LN=NA; converged=c(converged, NA); H4covLNI=rep(NA,2); H4corrLNI=rep(NA,2); DRM_H4_LN=rep(NA,length(data$x))
  parsH4LN <- NA
  div_H4_LN <- NA

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

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[12]))
  }
  if(prior.weights[12]>0){
    # print(12)

    fitstanLN4_LN = fun_sampling(stanmodels$mLN4, data, start,
                                 ndraws,nrchains,
                                 nriterations,warmup,
                                 delta,treedepth,seed,pvec)

    if(is.null(fitstanLN4_LN)){
      prior.weights[12] <- 0
      warning('difficulties fitting the Lognormal (Lognormal) model using sampling; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsLN4LN <- par_extract(fitstanLN4_LN, model_name = "LN4_LN")
      # parsLN4LN[,"BMD"] <- parsLN4LN[,"BMD"]*data$maxD # BMD on original scale

      #diagnostics here
      # posterior_diag(model_stan = fitstanLN4_LN)
      convergence_stat <- convergence_deci(model_stan = fitstanLN4_LN, nrchains = nrchains)
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)

      # divergent transitions
      div_LN4_LN <- sum(sapply(rstan::get_sampler_params(fitstanLN4_LN, inc_warmup = F),
                               function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))

      LN4resLNI=(quantile(as.matrix(fitstanLN4_LN)[,"par2"],pvec))*data$maxD
      LN4resLNI=c(LN4resLNI,apply(as.matrix(fitstanLN4_LN),2,median)[c("par1","par2","par3","par4","par5")])
      names(LN4resLNI)=c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t")

      LN4outLNI <- outLP(parsLN4LN, pvec, data$maxD)


      if(data$data_type == 2){
        DRM_LN4_LN = exp(DRM.LN4_LNI(LN4resLNI[4:7], data$x, data$q, shift=data$shift))
      }else if(data$data_type == 4){
        DRM_LN4_LN = exp(DRM.LN4_LND(LN4resLNI[4:7], data$x, data$q, shift=data$shift))
      }

      LN4covLNI = c(cov(as.matrix(fitstanLN4_LN)[,c("b","d")], use="complete.obs")["b","d"],
                    cov(as.matrix(fitstanLN4_LN)[,c("par2","d")], use="complete.obs")["par2","d"])

      LN4corrLNI = c(cor(as.matrix(fitstanLN4_LN)[,c("b","d")], use="complete.obs")["b","d"],
                     cor(as.matrix(fitstanLN4_LN)[,c("par2","d")], use="complete.obs")["par2","d"])

      BMDL=c(BMDL, LN4resLNI[1]); BMD=c(BMD, LN4resLNI[2]); BMDU=c(BMDU, LN4resLNI[3])
      bridgeLN4LN = bridgesampling::bridge_sampler(fitstanLN4_LN, silent = T)
    }
  }
  if(prior.weights[12] == 0){LN4resLNI=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeLN4LN=NA; converged=c(converged, NA); LN4covLNI=rep(NA,2); LN4corrLNI=rep(NA,2); DRM_LN4_LN=rep(NA,length(data$x))
  parsLN4LN <- NA
  div_LN4_LN <- NA

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

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[13]))
  }
  if(prior.weights[13]>0){
    # print(13)

    fitstanG4_LN = fun_sampling(stanmodels$mG4, data, start,
                                ndraws,nrchains,
                                nriterations,warmup,
                                delta,treedepth,seed,pvec)

    if(is.null(fitstanG4_LN)){
      prior.weights[13] <- 0
      warning('difficulties fitting the Gamma (Lognormal) model using sampling; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsG4LN <- par_extract(fitstanG4_LN, model_name = "G4_LN")
      # parsG4LN[,"BMD"] <- parsG4LN[,"BMD"]*data$maxD # BMD on original scale


      #diagnostics here
      # posterior_diag(model_stan = fitstanG4_LN)
      convergence_stat <- convergence_deci(model_stan = fitstanG4_LN, nrchains = nrchains)
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)

      # divergent transitions
      div_G4_LN <- sum(sapply(rstan::get_sampler_params(fitstanG4_LN, inc_warmup = F),
                              function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))

      G4resLNI=(quantile(as.matrix(fitstanG4_LN)[,"par2"],pvec))*data$maxD
      G4resLNI=c(G4resLNI,apply(as.matrix(fitstanG4_LN),2,median)[c("par1","par2","par3","par4","par5")])
      names(G4resLNI)=c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t")

      G4outLNI <- outLP(parsG4LN, pvec, data$maxD)


      if(data$data_type == 2){
        DRM_G4_LN = exp(DRM.G4_LNI(G4resLNI[4:7], data$x, data$q, shift=data$shift))
      }else if(data$data_type == 4){
        DRM_G4_LN = exp(DRM.G4_LND(G4resLNI[4:7], data$x, data$q, shift=data$shift))
      }

      G4covLNI = c(cov(as.matrix(fitstanG4_LN)[,c("b","d")], use="complete.obs")["b","d"],
                   cov(as.matrix(fitstanG4_LN)[,c("par2","d")], use="complete.obs")["par2","d"])

      G4corrLNI = c(cor(as.matrix(fitstanG4_LN)[,c("b","d")], use="complete.obs")["b","d"],
                    cor(as.matrix(fitstanG4_LN)[,c("par2","d")], use="complete.obs")["par2","d"])

      BMDL=c(BMDL, G4resLNI[1]); BMD=c(BMD, G4resLNI[2]); BMDU=c(BMDU, G4resLNI[3])
      bridgeG4LN = bridgesampling::bridge_sampler(fitstanG4_LN, silent = T)
    }
  }
  if(prior.weights[13] == 0){G4resLNI=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeG4LN=NA; converged=c(converged, NA); G4covLNI=rep(NA,2); G4corrLNI=rep(NA,2); DRM_G4_LN=rep(NA,length(data$x))
  parsG4LN <- NA
  div_G4_LN <- NA

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

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[14]))
  }
  if(prior.weights[14]>0){
    # print(14)

    fitstanQE4_LN = fun_sampling(stanmodels$mQE4, data, startQ,
                                 ndraws,nrchains,
                                 nriterations,warmup,
                                 delta,treedepth,seed,pvec)

    if(is.null(fitstanQE4_LN)){
      prior.weights[14] <- 0
      warning('difficulties fitting the Quadratic Exponential (Lognormal) model using sampling; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsQE4LN <- par_extract(fitstanQE4_LN, model_name = "QE4_LN")
      # parsQE4LN[,"BMD"] <- parsQE4LN[,"BMD"]*data$maxD # BMD on original scale

      #diagnostics here
      # posterior_diag(model_stan = fitstanQE4_LN)
      convergence_stat <- convergence_deci(model_stan = fitstanQE4_LN, nrchains = nrchains)
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)

      # divergent transitions
      div_QE4_LN <- sum(sapply(rstan::get_sampler_params(fitstanQE4_LN, inc_warmup = F),
                               function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))

      QE4resLNI=(quantile(as.matrix(fitstanQE4_LN)[,"par2"],pvec))*data$maxD
      QE4resLNI=c(QE4resLNI,apply(as.matrix(fitstanQE4_LN),2,median)[c("par1","par2","par3","par4","par5")])
      names(QE4resLNI)=c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t")

      QE4outLNI <- outLP(parsQE4LN, pvec, data$maxD)

      if(data$data_type == 2){
        DRM_QE4_LN = exp(DRM.QE4_LNI(QE4resLNI[4:7], data$x, data$q, shift=data$shift))
      }else if(data$data_type == 4){
        DRM_QE4_LN = exp(DRM.QE4_LND(QE4resLNI[4:7], data$x, data$q, shift=data$shift))
      }

      QE4covLNI = c(cov(as.matrix(fitstanQE4_LN)[,c("b","d")], use="complete.obs")["b","d"],
                    cov(as.matrix(fitstanQE4_LN)[,c("par2","d")], use="complete.obs")["par2","d"])

      QE4corrLNI = c(cor(as.matrix(fitstanQE4_LN)[,c("b","d")], use="complete.obs")["b","d"],
                     cor(as.matrix(fitstanQE4_LN)[,c("par2","d")], use="complete.obs")["par2","d"])

      BMDL=c(BMDL, QE4resLNI[1]); BMD=c(BMD, QE4resLNI[2]); BMDU=c(BMDU, QE4resLNI[3])
      bridgeQE4LN = bridgesampling::bridge_sampler(fitstanQE4_LN, silent = T)
    }
  }
  if(prior.weights[14] == 0){QE4resLNI=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeQE4LN=NA; converged=c(converged, NA); QE4covLNI=rep(NA,2); QE4corrLNI=rep(NA,2); DRM_QE4_LN=rep(NA,length(data$x))
  parsQE4LN <- NA
  div_QE4_LN <- NA

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

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[15]))
  }
  if(prior.weights[15]>0){
    # print(15)

    fitstanP4_LN = fun_sampling(stanmodels$mP4, data, start,
                                ndraws,nrchains,
                                nriterations,warmup,
                                delta,treedepth,seed,pvec)

    if(is.null(fitstanP4_LN)){
      prior.weights[15] <- 0
      warning('difficulties fitting the Probit (Lognormal) model using sampling; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsP4LN <- par_extract(fitstanP4_LN, model_name = "P4_LN")
      # parsP4LN[,"BMD"] <- parsP4LN[,"BMD"]*data$maxD # BMD on original scale

      #diagnostics here
      # posterior_diag(model_stan = fitstanP4_LN)
      convergence_stat <- convergence_deci(model_stan = fitstanP4_LN, nrchains = nrchains)
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)

      # divergent transitions
      div_P4_LN <- sum(sapply(rstan::get_sampler_params(fitstanP4_LN, inc_warmup = F),
                              function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))

      P4resLNI=(quantile(as.matrix(fitstanP4_LN)[,"par2"],pvec))*data$maxD
      P4resLNI=c(P4resLNI,apply(as.matrix(fitstanP4_LN),2,median)[c("par1","par2","par3","par4","par5")])
      names(P4resLNI)=c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t")

      P4outLNI <- outLP(parsP4LN, pvec, data$maxD)


      if(data$data_type == 2){
        DRM_P4_LN = exp(DRM.P4_LNI(P4resLNI[4:7], data$x, data$q, shift=data$shift))
      }else if(data$data_type == 4){
        DRM_P4_LN = exp(DRM.P4_LND(P4resLNI[4:7], data$x, data$q, shift=data$shift))
      }

      P4covLNI = c(cov(as.matrix(fitstanP4_LN)[,c("b","d")], use="complete.obs")["b","d"],
                   cov(as.matrix(fitstanP4_LN)[,c("par2","d")], use="complete.obs")["par2","d"])

      P4corrLNI = c(cor(as.matrix(fitstanP4_LN)[,c("b","d")], use="complete.obs")["b","d"],
                    cor(as.matrix(fitstanP4_LN)[,c("par2","d")], use="complete.obs")["par2","d"])

      BMDL=c(BMDL, P4resLNI[1]); BMD=c(BMD, P4resLNI[2]); BMDU=c(BMDU, P4resLNI[3])
      bridgeP4LN = bridge_sampler(fitstanP4_LN, silent = T)
    }
  }
  if(prior.weights[15] == 0){P4resLNI=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeP4LN=NA;
  converged=c(converged, NA); P4covLNI=rep(NA,2); P4corrLNI=rep(NA,2); DRM_P4_LN=rep(NA,length(data$x))
  parsP4LN <- NA
  div_P4_LN <- NA

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

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[16]))
  }
  if(prior.weights[16]>0){
    # print(16)

    fitstanL4_LN = fun_sampling(stanmodels$mL4, data, start,
                                ndraws,nrchains,
                                nriterations,warmup,
                                delta,treedepth,seed,pvec)

    if(is.null(fitstanL4_LN)){
      prior.weights[16] <- 0
      warning('difficulties fitting the Logit (Lognormal) model using sampling; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsL4LN <- par_extract(fitstanL4_LN, model_name = "L4_LN")
      # parsL4LN[,"BMD"] <- parsL4LN[,"BMD"]*data$maxD # BMD on original scale


      #diagnostics here
      # posterior_diag(model_stan = fitstanL4_LN)
      convergence_stat <- convergence_deci(model_stan = fitstanL4_LN, nrchains = nrchains)
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)

      # divergent transitions
      div_L4_LN <- sum(sapply(rstan::get_sampler_params(fitstanL4_LN, inc_warmup = F),
                              function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))

      L4resLNI=(quantile(as.matrix(fitstanL4_LN)[,"par2"],pvec))*data$maxD
      L4resLNI=c(L4resLNI,apply(as.matrix(fitstanL4_LN),2,median)[c("par1","par2","par3","par4","par5")])
      names(L4resLNI)=c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t")

      L4outLNI <- outLP(parsL4LN, pvec, data$maxD)


      if(data$data_type == 2){
        DRM_L4_LN = exp(DRM.L4_LNI(L4resLNI[4:7], data$x, data$q, shift=data$shift))
      }else if(data$data_type == 4){
        DRM_L4_LN = exp(DRM.L4_LND(L4resLNI[4:7], data$x, data$q, shift=data$shift))
      }

      L4covLNI = c(cov(as.matrix(fitstanL4_LN)[,c("b","d")], use="complete.obs")["b","d"],
                   cov(as.matrix(fitstanL4_LN)[,c("par2","d")], use="complete.obs")["par2","d"])

      L4corrLNI = c(cor(as.matrix(fitstanL4_LN)[,c("b","d")], use="complete.obs")["b","d"],
                    cor(as.matrix(fitstanL4_LN)[,c("par2","d")], use="complete.obs")["par2","d"])

      BMDL=c(BMDL, L4resLNI[1]); BMD=c(BMD, L4resLNI[2]); BMDU=c(BMDU, L4resLNI[3])
      bridgeL4LN = bridgesampling::bridge_sampler(fitstanL4_LN, silent = T)
    }
  }
  if(prior.weights[16] == 0){L4resLNI=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeL4LN=NA; converged=c(converged, NA); L4covLNI=rep(NA,2); L4corrLNI=rep(NA,2); DRM_L4_LN=rep(NA,length(data$x))
  parsL4LN <- NA
  div_L4_LN <- NA

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

  divergences <- c(div_E4_N, div_IE4_N, div_H4_N, div_LN4_N, div_G4_N, div_QE4_N, div_P4_N, div_L4_N,
                   div_E4_LN, div_IE4_LN, div_H4_LN, div_LN4_LN, div_G4_LN, div_QE4_LN, div_P4_LN, div_L4_LN)


  #-----------------------------------

  ### Weights based on bridge sampling are obtained by computing posterior model probabilities from the marginal likelihoods
  p.weights = prior.weights/sum(prior.weights==1)

  # lpwb = bridgesampling::post_prob(bridgeE4N, bridgeIE4N, bridgeH4N, bridgeLN4N, bridgeG4N, bridgeQE4N, bridgeP4N, bridgeL4N,
  #                  bridgeE4LN, bridgeIE4LN, bridgeH4LN, bridgeLN4LN, bridgeG4LN, bridgeQE4LN, bridgeP4LN, bridgeL4LN,
  #                  prior_prob = p.weights[p.weights>0],
  #                  model_names = c("E4_N","IE4_N","H4_N","LN4_N","G4_N","QE4_N","P4_N","L4_N",
  #                                  "E4_LN","IE4_LN","H4_LN","LN4_LN","G4_LN","QE4_LN","P4_LN","L4_LN")[p.weights>0])


  bridge.mods = c(if(p.weights[1]>0) bridgeE4N$logml,
                  if(p.weights[2]>0) bridgeIE4N$logml,
                  if(p.weights[3]>0) bridgeH4N$logml,
                  if(p.weights[4]>0) bridgeLN4N$logml,
                  if(p.weights[5]>0) bridgeG4N$logml,
                  if(p.weights[6]>0) bridgeQE4N$logml,
                  if(p.weights[7]>0) bridgeP4N$logml,
                  if(p.weights[8]>0) bridgeL4N$logml,
                  if(p.weights[9]>0) bridgeE4LN$logml,
                  if(p.weights[10]>0) bridgeIE4LN$logml,
                  if(p.weights[11]>0) bridgeH4LN$logml,
                  if(p.weights[12]>0) bridgeLN4LN$logml,
                  if(p.weights[13]>0) bridgeG4LN$logml,
                  if(p.weights[14]>0) bridgeQE4LN$logml,
                  if(p.weights[15]>0) bridgeP4LN$logml,
                  if(p.weights[16]>0) bridgeL4LN$logml)

  model_names = c("E4_N","IE4_N","H4_N","LN4_N","G4_N","QE4_N","P4_N","L4_N",
                  "E4_LN","IE4_LN","H4_LN","LN4_LN","G4_LN","QE4_LN","P4_LN","L4_LN")[p.weights>0]

  lpwb = bridgesampling::post_prob(bridge.mods,
                                   prior_prob = p.weights[p.weights>0],
                                   model_names = model_names)

  # # the model average posterior as a mixture
  count=round(lpwb*ndraws)
  names(count) = names(lpwb)
  mabmd1=(c( # normal
    if("E4_N" %in% names(count)) sample(as.matrix(fitstanE4_N)[,2],count[names(count)=="E4_N"],replace=T),
    if("IE4_N" %in% names(count)) sample(as.matrix(fitstanIE4_N)[,2],count[names(count)=="IE4_N"],replace=T),
    if("H4_N" %in% names(count)) sample(as.matrix(fitstanH4_N)[,2],count[names(count)=="H4_N"],replace=T),
    if("LN4_N" %in% names(count)) sample(as.matrix(fitstanLN4_N)[,2],count[names(count)=="LN4_N"],replace=T),
    if("G4_N" %in% names(count)) sample(as.matrix(fitstanG4_N)[,2],count[names(count)=="G4_N"],replace=T),
    if("QE4_N" %in% names(count)) sample(as.matrix(fitstanQE4_N)[,2],count[names(count)=="QE4_N"],replace=T),
    if("P4_N" %in% names(count)) sample(as.matrix(fitstanP4_N)[,2],count[names(count)=="P4_N"],replace=T),
    if("L4_N" %in% names(count)) sample(as.matrix(fitstanL4_N)[,2],count[names(count)=="L4_N"],replace=T),
    # lognormal
    if("E4_LN" %in% names(count)) sample(as.matrix(fitstanE4_LN)[,2],count[names(count)=="E4_LN"],replace=T),
    if("IE4_LN" %in% names(count)) sample(as.matrix(fitstanIE4_LN)[,2],count[names(count)=="IE4_LN"],replace=T),
    if("H4_LN" %in% names(count)) sample(as.matrix(fitstanH4_LN)[,2],count[names(count)=="H4_LN"],replace=T),
    if("LN4_LN" %in% names(count)) sample(as.matrix(fitstanLN4_LN)[,2],count[names(count)=="LN4_LN"],replace=T),
    if("G4_LN" %in% names(count)) sample(as.matrix(fitstanG4_LN)[,2],count[names(count)=="G4_LN"],replace=T),
    if("QE4_LN" %in% names(count)) sample(as.matrix(fitstanQE4_LN)[,2],count[names(count)=="QE4_LN"],replace=T),
    if("P4_LN" %in% names(count)) sample(as.matrix(fitstanP4_LN)[,2],count[names(count)=="P4_LN"],replace=T),
    if("L4_LN" %in% names(count)) sample(as.matrix(fitstanL4_LN)[,2],count[names(count)=="L4_LN"],replace=T)
  ))

  macib=(quantile(mabmd1,pvec))*data$maxD
  names(macib)=c("BMDL","BMD","BMDU") # original scale

  if(TRUE %in% (mabmd1 > data$maxD) && data$maxD > 1){
    mabmd1 = ifelse(mabmd1 > data$maxD, data$maxD, mabmd1)
    p.msg = 'The model averaged posterior distribution has been truncated at max(Dose)^2'
    warnings('The model averaged posterior distribution has been truncated at max(Dose)^2')
  }else{
    p.msg = ''
  }

  BMDq_bs = (quantile(mabmd1, seq(0,1,0.005)))*data$maxD


  mods = c("E4_N","IE4_N","H4_N","LN4_N","G4_N","QE4_N","P4_N","L4_N",
           "E4_LN","IE4_LN","H4_LN","LN4_LN","G4_LN","QE4_LN","P4_LN","L4_LN")
  w.bs = rep(0, 16)
  names(w.bs) = mods
  for(i in mods[mods%in%names(lpwb)]){
    w.bs[i] = lpwb[i]
  }

  ### Model-averaged response per dose level
  dr.MA.bs <- c()
  for(i in 1:length(data$x)){
    dr.MA.bs[i] = weighted.mean(x = c(DRM_E4_N[i],DRM_IE4_N[i],DRM_H4_N[i],DRM_LN4_N[i],DRM_G4_N[i],DRM_QE4_N[i],
                                      DRM_P4_N[i], DRM_L4_N[i] ,
                                      DRM_E4_LN[i], DRM_IE4_LN[i], DRM_H4_LN[i], DRM_LN4_LN[i], DRM_G4_LN[i],
                                      DRM_QE4_LN[i], DRM_P4_LN[i],DRM_L4_LN[i]),
                                w = w.bs,
                                na.rm = T)

  }


  ##### Weights & MA if one of the models is divergent --> this model gets weight 0
  if((0 %in% converged) && (1 %in% converged)){
    prior.weights.new = prior.weights
    div.models = which(converged == 0)
    prior.weights.new[div.models] = 0
    p.weights.new = prior.weights.new/sum(prior.weights.new==1)

    bridge.mods = c(if(p.weights.new[1]>0) bridgeE4N$logml,
                    if(p.weights.new[2]>0) bridgeIE4N$logml,
                    if(p.weights.new[3]>0) bridgeH4N$logml,
                    if(p.weights.new[4]>0) bridgeLN4N$logml,
                    if(p.weights.new[5]>0) bridgeG4N$logml,
                    if(p.weights.new[6]>0) bridgeQE4N$logml,
                    if(p.weights.new[7]>0) bridgeP4N$logml,
                    if(p.weights.new[8]>0) bridgeL4N$logml,
                    if(p.weights.new[9]>0) bridgeE4LN$logml,
                    if(p.weights.new[10]>0) bridgeIE4LN$logml,
                    if(p.weights.new[11]>0) bridgeH4LN$logml,
                    if(p.weights.new[12]>0) bridgeLN4LN$logml,
                    if(p.weights.new[13]>0) bridgeG4LN$logml,
                    if(p.weights.new[14]>0) bridgeQE4LN$logml,
                    if(p.weights.new[15]>0) bridgeP4LN$logml,
                    if(p.weights.new[16]>0) bridgeL4LN$logml)

    model_names = c("E4_N","IE4_N","H4_N","LN4_N","G4_N","QE4_N","P4_N","L4_N",
                    "E4_LN","IE4_LN","H4_LN","LN4_LN","G4_LN","QE4_LN","P4_LN","L4_LN")[p.weights.new>0]

    lpwb.conv = bridgesampling::post_prob(bridge.mods,
                                          prior_prob = p.weights.new[p.weights.new>0],
                                          model_names = model_names)

    # # the model average posterior as a mixture
    count=round(lpwb.conv*ndraws)
    names(count) = names(lpwb.conv)
    mabmd.conv1=(c( # normal
      if("E4_N" %in% names(count)) sample(as.matrix(fitstanE4_N)[,2],count[names(count)=="E4_N"],replace=T),
      if("IE4_N" %in% names(count)) sample(as.matrix(fitstanIE4_N)[,2],count[names(count)=="IE4_N"],replace=T),
      if("H4_N" %in% names(count)) sample(as.matrix(fitstanH4_N)[,2],count[names(count)=="H4_N"],replace=T),
      if("LN4_N" %in% names(count)) sample(as.matrix(fitstanLN4_N)[,2],count[names(count)=="LN4_N"],replace=T),
      if("G4_N" %in% names(count)) sample(as.matrix(fitstanG4_N)[,2],count[names(count)=="G4_N"],replace=T),
      if("QE4_N" %in% names(count)) sample(as.matrix(fitstanQE4_N)[,2],count[names(count)=="QE4_N"],replace=T),
      if("P4_N" %in% names(count)) sample(as.matrix(fitstanP4_N)[,2],count[names(count)=="P4_N"],replace=T),
      if("L4_N" %in% names(count)) sample(as.matrix(fitstanL4_N)[,2],count[names(count)=="L4_N"],replace=T),
      # lognormal
      if("E4_LN" %in% names(count)) sample(as.matrix(fitstanE4_LN)[,2],count[names(count)=="E4_LN"],replace=T),
      if("IE4_LN" %in% names(count)) sample(as.matrix(fitstanIE4_LN)[,2],count[names(count)=="IE4_LN"],replace=T),
      if("H4_LN" %in% names(count)) sample(as.matrix(fitstanH4_LN)[,2],count[names(count)=="H4_LN"],replace=T),
      if("LN4_LN" %in% names(count)) sample(as.matrix(fitstanLN4_LN)[,2],count[names(count)=="LN4_LN"],replace=T),
      if("G4_LN" %in% names(count)) sample(as.matrix(fitstanG4_LN)[,2],count[names(count)=="G4_LN"],replace=T),
      if("QE4_LN" %in% names(count)) sample(as.matrix(fitstanQE4_LN)[,2],count[names(count)=="QE4_LN"],replace=T),
      if("P4_LN" %in% names(count)) sample(as.matrix(fitstanP4_LN)[,2],count[names(count)=="P4_LN"],replace=T),
      if("L4_LN" %in% names(count)) sample(as.matrix(fitstanL4_LN)[,2],count[names(count)=="L4_LN"],replace=T)
    ))

    macib.conv=(quantile(mabmd.conv1,pvec))*data$maxD
    names(macib.conv)=c("BMDL","BMD","BMDU")

    if(TRUE %in% (mabmd.conv1 > data$maxD) && data$maxD > 1){
      mabmd.conv1 = ifelse(mabmd.conv1 > data$maxD, data$maxD, mabmd.conv1)
      p.msg = 'The model averaged posterior distribution has been truncated at max(Dose)^2'
      warnings('The model averaged posterior distribution has been truncated at max(Dose)^2')
    }else{
      p.msg = ''
    }

    BMDq_bs_conv = (quantile(mabmd.conv1, seq(0,1,0.005)))*data$maxD


    mods = c("E4_N","IE4_N","H4_N","LN4_N","G4_N","QE4_N","P4_N","L4_N",
             "E4_LN","IE4_LN","H4_LN","LN4_LN","G4_LN","QE4_LN","P4_LN","L4_LN")
    w.bs.conv = rep(0, 16)
    names(w.bs.conv) = mods
    for(i in mods[mods%in%names(lpwb.conv)]){
      w.bs.conv[i] = lpwb.conv[i]
    }

    ### Model-averaged response per dose level
    dr.MA.bs.conv <- c()
    for(i in 1:length(data$x)){
      dr.MA.bs.conv[i] = weighted.mean(x = c(DRM_E4_N[i],DRM_IE4_N[i],DRM_H4_N[i],DRM_LN4_N[i],DRM_G4_N[i],
                                             DRM_QE4_N[i], DRM_P4_N[i], DRM_L4_N[i] ,
                                             DRM_E4_LN[i], DRM_IE4_LN[i], DRM_H4_LN[i], DRM_LN4_LN[i],
                                             DRM_G4_LN[i],DRM_QE4_LN[i], DRM_P4_LN[i],DRM_L4_LN[i]),
                                       w = w.bs.conv,
                                       na.rm = T)
    }

  }else{ ## If all converged or NONE converged
    w.bs.conv = NULL; macib.conv = NULL; mabmd.conv1 = NA; BMDq_bs_conv = NULL; dr.MA.bs.conv = NULL
  }


  #----------------------------------
  ### weights based on laplace approximation

  prior.weights = prior.weights.orig

  # normal
  data=data.N$data
  start=data.N$start
  startQ=data.N$startQ

  llN = c() # likelihoods

  # getting the posterior modes and the hessian and the model specific posterior distributions
  if(prior.weights[1]>0){

    optE4_NI <- fun_optim(stanmodels$mE4, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optE4_NI[[3]]),TRUE,(optE4_NI[[3]]!=0)) | length(optE4_NI)!=9)){
      prior.weights[1] <- 0
      warning('difficulties fitting the Exponential (Normal) model using Laplace; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      if(data$data_type == 1){
        llE4N=llfE4_NI(optE4_NI$par[c(1,2,9,4,5)],nvec=data$n,
                       dvec=data$x,mvec=data$m,
                       s2vec=data$s2,qval=data$q)
      }else if(data$data_type == 3){
        llE4N=llfE4_ND(optE4_NI$par[c(1,2,9,4,5)],nvec=data$n,
                       dvec=data$x,mvec=data$m,
                       s2vec=data$s2,qval=data$q)
      }


      llN = c(llN, llE4N)
    }
  }
  if(prior.weights[1] == 0){llN = c(llN,NA)}
  #
  if(prior.weights[2]>0){

    optIE4_NI <- fun_optim(stanmodels$mIE4, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optIE4_NI[[3]]),TRUE,(optIE4_NI[[3]]!=0)) | length(optIE4_NI)!=9)){
      prior.weights[2] <- 0
      warning('difficulties fitting the Inverse Exponential (Normal) model using Laplace; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      if(data$data_type == 1){
        llIE4N=llfIE4_NI(optIE4_NI$par[c(1,2,9,4,5)],nvec=data$n,
                         dvec=data$x,mvec=data$m,
                         s2vec=data$s2,qval=data$q)
      }else if(data$data_type == 3){
        llIE4N=llfIE4_ND(optIE4_NI$par[c(1,2,9,4,5)],nvec=data$n,
                         dvec=data$x,mvec=data$m,
                         s2vec=data$s2,qval=data$q)
      }

      llN = c(llN, llIE4N)
    }
  }
  if(prior.weights[2] == 0){llN = c(llN,NA)}
  #
  if(prior.weights[3]>0){

    optH4_NI <- fun_optim(stanmodels$mH4, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optH4_NI[[3]]),TRUE,(optH4_NI[[3]]!=0)) | length(optH4_NI)!=9)){
      prior.weights[3] <- 0
      warning('difficulties fitting the Hill (Normal) model using Laplace; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      if(data$data_type == 1){
        llH4N=llfH4_NI(optH4_NI$par[c(1,2,9,4,5)],nvec=data$n,
                       dvec=data$x,mvec=data$m,
                       s2vec=data$s2,qval=data$q)
      }else if(data$data_type == 3){
        llH4N=llfH4_ND(optH4_NI$par[c(1,2,9,4,5)],nvec=data$n,
                       dvec=data$x,mvec=data$m,
                       s2vec=data$s2,qval=data$q)
      }

      llN = c(llN, llH4N)
    }
  }
  if(prior.weights[3] == 0){llN = c(llN,NA)}
  #
  if(prior.weights[4]>0){

    optLN4_NI <- fun_optim(stanmodels$mLN4, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optLN4_NI[[3]]),TRUE,(optLN4_NI[[3]]!=0)) | length(optLN4_NI)!=9)){
      prior.weights[4] <- 0
      warning('difficulties fitting the Lognormal (Normal) model using Laplace; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      if(data$data_type == 1){
        llLN4N=llfLN4_NI(optLN4_NI$par[c(1,2,9,4,5)],nvec=data$n,
                         dvec=data$x,mvec=data$m,
                         s2vec=data$s2,qval=data$q)
      }else if(data$data_type == 3){
        llLN4N=llfLN4_ND(optLN4_NI$par[c(1,2,9,4,5)],nvec=data$n,
                         dvec=data$x,mvec=data$m,
                         s2vec=data$s2,qval=data$q)
      }


      llN = c(llN, llLN4N)
    }
  }
  if(prior.weights[4] == 0){llN = c(llN,NA)}
  #
  if(prior.weights[5]>0){

    optG4_NI <- fun_optim(stanmodels$mG4, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optG4_NI[[3]]),TRUE,(optG4_NI[[3]]!=0)) | length(optG4_NI)!=9)){
      prior.weights[5] <- 0
      warning('difficulties fitting the Gamma (Normal) model using Laplace; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      if(data$data_type == 1){
        llG4N=llfG4_NI(optG4_NI$par[c(1,2,9,4,5)],nvec=data$n,
                       dvec=data$x,mvec=data$m,
                       s2vec=data$s2,qval=data$q)
      }else if(data$data_type == 3){
        llG4N=llfG4_ND(optG4_NI$par[c(1,2,9,4,5)],nvec=data$n,
                       dvec=data$x,mvec=data$m,
                       s2vec=data$s2,qval=data$q)
      }


      llN = c(llN, llG4N)
    }
  }
  if(prior.weights[5] == 0){llN = c(llN,NA)}
  #
  if(prior.weights[6]>0){

    optQE4_NI <- fun_optim(stanmodels$mQE4, data, startQ, ndraws, 123, pvec)

    if((ifelse(is.na(optQE4_NI[[3]]),TRUE,(optQE4_NI[[3]]!=0)) | length(optQE4_NI)!=9)){
      prior.weights[6] <- 0
      warning('difficulties fitting the Quadratic Exponential (Normal) model using Laplace; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      if(data$data_type == 1){
        llQE4N=llfQE4_NI(optQE4_NI$par[c(1,2,9,4,5)],nvec=data$n,
                         dvec=data$x,mvec=data$m,
                         s2vec=data$s2,qval=data$q)
      }else if(data$data_type == 3){
        llQE4N=llfQE4_ND(optQE4_NI$par[c(1,2,9,4,5)],nvec=data$n,
                         dvec=data$x,mvec=data$m,
                         s2vec=data$s2,qval=data$q)
      }


      llN = c(llN, llQE4N)
    }
  }
  if(prior.weights[6] == 0){llN = c(llN,NA)}
  #
  if(prior.weights[7]>0){

    optP4_NI <- fun_optim(stanmodels$mP4, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optP4_NI[[3]]),TRUE,(optP4_NI[[3]]!=0)) | length(optP4_NI)!=9)){
      prior.weights[7] <- 0
      warning('difficulties fitting the Probit (Normal) model using Laplace; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      if(data$data_type == 1){
        llP4N=llfP4_NI(optP4_NI$par[c(1,2,9,4,5)],nvec=data$n,
                       dvec=data$x,mvec=data$m,
                       s2vec=data$s2,qval=data$q)
      }else if(data$data_type == 3){
        llP4N=llfP4_ND(optP4_NI$par[c(1,2,9,4,5)],nvec=data$n,
                       dvec=data$x,mvec=data$m,
                       s2vec=data$s2,qval=data$q)
      }

      llN = c(llN, llP4N)
    }
  }
  if(prior.weights[7] == 0){llN = c(llN,NA)}
  #
  if(prior.weights[8]>0){

    optL4_NI <- fun_optim(stanmodels$mL4, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optL4_NI[[3]]),TRUE,(optL4_NI[[3]]!=0)) | length(optL4_NI)!=9)){
      prior.weights[8] <- 0
      warning('difficulties fitting the Logit (Normal) model using Laplace; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      if(data$data_type == 1){
        llL4N=llfL4_NI(optL4_NI$par[c(1,2,9,4,5)],nvec=data$n,
                       dvec=data$x,mvec=data$m,
                       s2vec=data$s2,qval=data$q)
      }else if(data$data_type == 3){
        llL4N=llfL4_ND(optL4_NI$par[c(1,2,9,4,5)],nvec=data$n,
                       dvec=data$x,mvec=data$m,
                       s2vec=data$s2,qval=data$q)
      }


      llN = c(llN, llL4N)
    }
  }
  if(prior.weights[8] == 0){llN = c(llN,NA)}


  # lognormal
  data=data.LN$data
  start=data.LN$start
  startQ=data.LN$startQ

  llLN = c()

  if(prior.weights[9]>0){

    optE4_LNI <- fun_optim(stanmodels$mE4, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optE4_LNI[[3]]),TRUE,(optE4_LNI[[3]]!=0)) | length(optE4_LNI)!=9)){
      prior.weights[9] <- 0
      warning('difficulties fitting the Exponential (Lognormal) model using Laplace; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      if(data$data_type == 2){
        llE4LN=llfE4_LNI(optE4_LNI$par[c(1,2,9,4,5)],nvec=data$n,
                         dvec=data$x,mvec=data$m,
                         s2vec=data$s2,qval=data$q, shift=data$shift)
      }else if(data$data_type == 4){
        llE4LN=llfE4_LND(optE4_LNI$par[c(1,2,9,4,5)],nvec=data$n,
                         dvec=data$x,mvec=data$m,
                         s2vec=data$s2,qval=data$q, shift=data$shift)
      }


      llLN = c(llLN, llE4LN)
    }
  }
  if(prior.weights[9] == 0){llLN = c(llLN,NA)}
  #
  if(prior.weights[10]>0){

    optIE4_LNI <- fun_optim(stanmodels$mIE4, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optIE4_LNI[[3]]),TRUE,(optIE4_LNI[[3]]!=0)) | length(optIE4_LNI)!=9)){
      prior.weights[10] <- 0
      warning('difficulties fitting the Inverse Exponential (Lognormal) model using Laplace; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      if(data$data_type == 2){
        llIE4LN=llfIE4_LNI(optIE4_LNI$par[c(1,2,9,4,5)],nvec=data$n,
                           dvec=data$x,mvec=data$m,
                           s2vec=data$s2,qval=data$q, shift=data$shift)
      }else if(data$data_type == 4){
        llIE4LN=llfIE4_LND(optIE4_LNI$par[c(1,2,9,4,5)],nvec=data$n,
                           dvec=data$x,mvec=data$m,
                           s2vec=data$s2,qval=data$q, shift=data$shift)
      }


      llLN = c(llLN, llIE4LN)
    }
  }
  if(prior.weights[10] == 0){llLN = c(llLN,NA)}
  #
  if(prior.weights[11]>0){

    optH4_LNI <- fun_optim(stanmodels$mH4, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optH4_LNI[[3]]),TRUE,(optH4_LNI[[3]]!=0)) | length(optH4_LNI)!=9)){
      prior.weights[11] <- 0
      warning('difficulties fitting the Hill (Lognormal) model using Laplace; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      if(data$data_type == 2){
        llH4LN=llfH4_LNI(optH4_LNI$par[c(1,2,9,4,5)],nvec=data$n,
                         dvec=data$x,mvec=data$m,
                         s2vec=data$s2,qval=data$q, shift=data$shift)
      }else if(data$data_type == 4){
        llH4LN=llfH4_LND(optH4_LNI$par[c(1,2,9,4,5)],nvec=data$n,
                         dvec=data$x,mvec=data$m,
                         s2vec=data$s2,qval=data$q, shift=data$shift)
      }

      llLN = c(llLN, llH4LN)
    }
  }
  if(prior.weights[11] == 0){llLN = c(llLN,NA)}
  #
  if(prior.weights[12]>0){

    optLN4_LNI <- fun_optim(stanmodels$mLN4, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optLN4_LNI[[3]]),TRUE,(optLN4_LNI[[3]]!=0)) | length(optLN4_LNI)!=9)){
      prior.weights[12] <- 0
      warning('difficulties fitting the Lognormal (Lognormal) model using Laplace; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      if(data$data_type == 2){
        llLN4LN=llfLN4_LNI(optLN4_LNI$par[c(1,2,9,4,5)],nvec=data$n,
                           dvec=data$x,mvec=data$m,
                           s2vec=data$s2,qval=data$q, shift=data$shift)
      }else if(data$data_type == 4){
        llLN4LN=llfLN4_LND(optLN4_LNI$par[c(1,2,9,4,5)],nvec=data$n,
                           dvec=data$x,mvec=data$m,
                           s2vec=data$s2,qval=data$q, shift=data$shift)
      }

      llLN = c(llLN, llLN4LN)
    }
  }
  if(prior.weights[12] == 0){llLN = c(llLN,NA)}
  #
  if(prior.weights[13]>0){

    optG4_LNI <- fun_optim(stanmodels$mG4, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optG4_LNI[[3]]),TRUE,(optG4_LNI[[3]]!=0)) | length(optG4_LNI)!=9)){
      prior.weights[13] <- 0
      warning('difficulties fitting the Gamma (Lognormal) model using Laplace; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      if(data$data_type == 2){
        llG4LN=llfG4_LNI(optG4_LNI$par[c(1,2,9,4,5)],nvec=data$n,
                         dvec=data$x,mvec=data$m,
                         s2vec=data$s2,qval=data$q, shift=data$shift)
      }else if(data$data_type == 4){
        llG4LN=llfG4_LND(optG4_LNI$par[c(1,2,9,4,5)],nvec=data$n,
                         dvec=data$x,mvec=data$m,
                         s2vec=data$s2,qval=data$q, shift=data$shift)
      }

      llLN = c(llLN, llG4LN)
    }
  }
  if(prior.weights[13] == 0){llLN = c(llLN,NA)}
  #
  if(prior.weights[14]>0){

    optQE4_LNI <- fun_optim(stanmodels$mQE4, data, startQ, ndraws, 123, pvec)

    if((ifelse(is.na(optQE4_LNI[[3]]),TRUE,(optQE4_LNI[[3]]!=0)) | length(optQE4_LNI)!=9)){
      prior.weights[14] <- 0
      warning('difficulties fitting the Quadratic Exponential (Lognormal) model using Laplace; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      if(data$data_type == 2){
        llQE4LN=llfQE4_LNI(optQE4_LNI$par[c(1,2,9,4,5)],nvec=data$n,
                           dvec=data$x,mvec=data$m,
                           s2vec=data$s2,qval=data$q, shift=data$shift)
      }else if(data$data_type == 4){
        llQE4LN=llfQE4_LND(optQE4_LNI$par[c(1,2,9,4,5)],nvec=data$n,
                           dvec=data$x,mvec=data$m,
                           s2vec=data$s2,qval=data$q, shift=data$shift)
      }

      llLN = c(llLN, llQE4LN)
    }
  }
  if(prior.weights[14] == 0){llLN = c(llLN,NA)}
  #
  if(prior.weights[15]>0){

    optP4_LNI <- fun_optim(stanmodels$mP4, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optP4_LNI[[3]]),TRUE,(optP4_LNI[[3]]!=0)) | length(optP4_LNI)!=9)){
      prior.weights[15] <- 0
      warning('difficulties fitting the Probit (Lognormal) model using Laplace; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      if(data$data_type == 2){
        llP4LN=llfP4_LNI(optP4_LNI$par[c(1,2,9,4,5)],nvec=data$n,
                         dvec=data$x,mvec=data$m,
                         s2vec=data$s2,qval=data$q, shift=data$shift)
      }else if(data$data_type == 4){
        llP4LN=llfP4_LND(optP4_LNI$par[c(1,2,9,4,5)],nvec=data$n,
                         dvec=data$x,mvec=data$m,
                         s2vec=data$s2,qval=data$q, shift=data$shift)
      }

      llLN = c(llLN, llP4LN)
    }
  }
  if(prior.weights[15] == 0){llLN = c(llLN,NA)}
  #
  if(prior.weights[16]>0){

    optL4_LNI <- fun_optim(stanmodels$mL4, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optL4_LNI[[3]]),TRUE,(optL4_LNI[[3]]!=0)) | length(optL4_LNI)!=9)){
      prior.weights[16] <- 0
      warning('difficulties fitting the Logit (Lognormal) model using Laplace; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      if(data$data_type == 2){
        llL4LN=llfL4_LNI(optL4_LNI$par[c(1,2,9,4,5)],nvec=data$n,
                         dvec=data$x,mvec=data$m,
                         s2vec=data$s2,qval=data$q, shift=data$shift)
      }else if(data$data_type == 4){
        llL4LN=llfL4_LND(optL4_LNI$par[c(1,2,9,4,5)],nvec=data$n,
                         dvec=data$x,mvec=data$m,
                         s2vec=data$s2,qval=data$q, shift=data$shift)
      }

      llLN = c(llLN, llL4LN)
    }
  }
  if(prior.weights[16] == 0){llLN = c(llLN,NA)}

  minll = min(llN,llLN,na.rm=T)

  # the weights


  fun.w <- function(DIH, ll, min.ll, opt, mu, sig, lb, ub, s1, s2, s3, td){

    w1 <- (2*pi)^(2.5)*sqrt(DIH)*exp(ll-min.ll)*
      # mvtnorm::dmvnorm(opt$par[c(4,5)],mean=mu[c(4,5)],
      #                  sigma=sig[c(4,5),c(4,5)])*
      # sigma2
      dnorm(opt$par[5], mu[5], sig[5,5])*
      # d
      truncnorm::dtruncnorm(opt$par[4], b = td, mean = mu[4], sd = sig[4,4])*
      #
      mc2d::dpert(opt$par[1], min = lb[1], max = ub[1],
                  mode = mu[1], shape = s1)*
      mc2d::dpert(opt$par[2], min = lb[2], max = ub[2],
                  mode = mu[2], shape = s2)*
      mc2d::dpert(opt$par[9], min = lb[3], max = ub[3],
                  mode = mu[3], shape = s3)
    return(w1)
  }

  w=c()

  # normal

  data=data.N$data
  start=data.N$start
  start=data.N$startQ

  lls <- c(llN, llLN)
  w.msg <- ''

  max.ll = max(lls, na.rm = T)
  if(is.na(lls[which((max.ll-lls[!is.na(lls)]) < 709 & prior.weights>0)][1])){
    lpw <- rep(0, 16)
    lpw[which(lls == max.ll)] <- 1
    w.msg <- 'Laplace weights could not be computed and one model gets all the weight; using another prior for parameter d might help'
    warning('Laplace weights could not be computed and one model gets all the weight; using another prior for parameter d might help')
  }else{

    if(FALSE %in% ((max.ll-lls[!is.na(lls)]) < 709)){
      w.msg <- 'not all models were used in computation of Laplace weights, some models set to 0; using another prior for parameter d might help'
      warning('not all models were used in computation of Laplace weights, some models set to 0; using another prior for parameter d might help')
    }

    minll <- min(lls[which((max.ll-lls[!is.na(lls)]) < 709 & prior.weights>0)], na.rm = T)

    if(prior.weights[1]>0){
      DIHE4h=det(-solve(optE4_NI$hessian))
      DIHE4=ifelse(DIHE4h<0,0,DIHE4h)
      w = c(w, fun.w(DIHE4, llE4N, minll, optE4_NI, data$priormu, data$priorSigma,
                     data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c,
                     data$truncd))

    }else{w=c(w,0)}

    if(prior.weights[2]>0){
      DIHIE4h=det(-solve(optIE4_NI$hessian))
      DIHIE4=ifelse(DIHIE4h<0,0,DIHIE4h)
      w = c(w, fun.w(DIHIE4, llIE4N, minll, optIE4_NI, data$priormu, data$priorSigma,
                     data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c,
                     data$truncd))
    }else{w=c(w,0)}

    if(prior.weights[3]>0){
      DIHH4h=det(-solve(optH4_NI$hessian))
      DIHH4=ifelse(DIHH4h<0,0,DIHH4h)
      w = c(w, fun.w(DIHH4, llH4N, minll, optH4_NI, data$priormu, data$priorSigma,
                     data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c,
                     data$truncd))
    }else{w=c(w,0)}

    if(prior.weights[4]>0){
      DIHLN4h=det(-solve(optLN4_NI$hessian))
      DIHLN4=ifelse(DIHLN4h<0,0,DIHLN4h)
      w = c(w, fun.w(DIHLN4, llLN4N, minll, optLN4_NI, data$priormu, data$priorSigma,
                     data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c,
                     data$truncd))
    }else{w=c(w,0)}

    if(prior.weights[5]>0){
      DIHG4h=det(-solve(optG4_NI$hessian))
      DIHG4=ifelse(DIHG4h<0,0,DIHG4h)
      w = c(w, fun.w(DIHG4, llG4N, minll, optG4_NI, data$priormu, data$priorSigma,
                     data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c,
                     data$truncd))
    }else{w=c(w,0)}

    if(prior.weights[6]>0){
      DIHQE4h=det(-solve(optQE4_NI$hessian))
      DIHQE4=ifelse(DIHQE4h<0,0,DIHQE4h)
      # w = c(w, fun.w(DIHQE4, llQE4N, minll, optQE4_NI, data$priormu, data$priorSigma,
      #                data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c))
      w = c(w, fun.w(DIHQE4, llQE4N, minll, optQE4_NI, data$priormuQ, data$priorSigmaQ,
                     data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c,
                     data$truncdQ))
    }else{w=c(w,0)}

    if(prior.weights[7]>0){
      DIHP4h=det(-solve(optP4_NI$hessian))
      DIHP4=ifelse(DIHP4h<0,0,DIHP4h)
      w = c(w, fun.w(DIHP4, llP4N, minll, optP4_NI, data$priormu, data$priorSigma,
                     data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c,
                     data$truncd))
    }else{w=c(w,0)}

    if(prior.weights[8]>0){
      DIHL4h=det(-solve(optL4_NI$hessian))
      DIHL4=ifelse(DIHL4h<0,0,DIHL4h)
      w = c(w, fun.w(DIHL4, llL4N, minll, optL4_NI, data$priormu, data$priorSigma,
                     data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c,
                     data$truncd))
    }else{w=c(w,0)}


    # lognormal

    data=data.LN$data
    start=data.LN$start
    start=data.LN$startQ

    if(prior.weights[9]>0){
      DIHE4h=det(-solve(optE4_LNI$hessian))
      DIHE4=ifelse(DIHE4h<0,0,DIHE4h)
      w = c(w, fun.w(DIHE4, llE4LN, minll, optE4_LNI, data$priormu, data$priorSigma,
                     data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c,
                     data$truncd))
    }else{w=c(w,0)}

    if(prior.weights[10]>0){
      DIHIE4h=det(-solve(optIE4_LNI$hessian))
      DIHIE4=ifelse(DIHIE4h<0,0,DIHIE4h)
      w = c(w, fun.w(DIHIE4, llIE4LN, minll, optIE4_LNI, data$priormu, data$priorSigma,
                     data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c,
                     data$truncd))
    }else{w=c(w,0)}

    if(prior.weights[11]>0){
      DIHH4h=det(-solve(optH4_LNI$hessian))
      DIHH4=ifelse(DIHH4h<0,0,DIHH4h)
      w = c(w, fun.w(DIHH4, llH4LN, minll, optH4_LNI, data$priormu, data$priorSigma,
                     data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c,
                     data$truncd))
    }else{w=c(w,0)}

    if(prior.weights[12]>0){
      DIHLN4h=det(-solve(optLN4_LNI$hessian))
      DIHLN4=ifelse(DIHLN4h<0,0,DIHLN4h)
      w = c(w, fun.w(DIHLN4, llLN4LN, minll, optLN4_LNI, data$priormu, data$priorSigma,
                     data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c,
                     data$truncd))
    }else{w=c(w,0)}

    if(prior.weights[13]>0){
      DIHG4h=det(-solve(optG4_LNI$hessian))
      DIHG4=ifelse(DIHG4h<0,0,DIHG4h)
      w = c(w, fun.w(DIHG4, llG4LN, minll, optG4_LNI, data$priormu, data$priorSigma,
                     data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c,
                     data$truncd))
    }else{w=c(w,0)}

    if(prior.weights[14]>0){
      DIHQE4h=det(-solve(optQE4_LNI$hessian))
      DIHQE4=ifelse(DIHQE4h<0,0,DIHQE4h)
      # w = c(w, fun.w(DIHQE4, llQE4LN, minll, optQE4_LNI, data$priormu, data$priorSigma,
      #                data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c))
      w = c(w, fun.w(DIHQE4, llQE4LN, minll, optQE4_LNI, data$priormuQ, data$priorSigmaQ,
                     data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c,
                     data$truncdQ))
    }else{w=c(w,0)}

    if(prior.weights[15]>0){
      DIHP4h=det(-solve(optP4_LNI$hessian))
      DIHP4=ifelse(DIHP4h<0,0,DIHP4h)
      w = c(w, fun.w(DIHP4, llP4LN, minll, optP4_LNI, data$priormu, data$priorSigma,
                     data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c,
                     data$truncd))
    }else{w=c(w,0)}

    if(prior.weights[16]>0){
      DIHL4h=det(-solve(optL4_LNI$hessian))
      DIHL4=ifelse(DIHL4h<0,0,DIHL4h)
      w = c(w, fun.w(DIHL4, llL4LN, minll, optL4_LNI, data$priormu, data$priorSigma,
                     data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c,
                     data$truncd))
    }else{w=c(w,0)}
  }

  w <- ifelse(w == 'Inf' | is.na(w), 0, w)
  lpwlp=(prior.weights*w)/sum(prior.weights*w)
  names(lpwlp) = c("E4_N","IE4_N","H4_N","LN4_N","G4_N","QE4_N","P4_N","L4_N","E4_LN","IE4_LN","H4_LN","LN4_LN","G4_LN","QE4_LN","P4_LN","L4_LN")

  # the model average posterior as a mixture
  count=round(lpwlp*ndraws)
  mabmd=(c(# normal
    if(prior.weights[1]>0) sample(as.matrix(fitstanE4_N)[,2],count[1],replace=T),
    if(prior.weights[2]>0) sample(as.matrix(fitstanIE4_N)[,2],count[2],replace=T),
    if(prior.weights[3]>0) sample(as.matrix(fitstanH4_N)[,2],count[3],replace=T),
    if(prior.weights[4]>0) sample(as.matrix(fitstanLN4_N)[,2],count[4],replace=T),
    if(prior.weights[5]>0) sample(as.matrix(fitstanG4_N)[,2],count[5],replace=T),
    if(prior.weights[6]>0) sample(as.matrix(fitstanQE4_N)[,2],count[6],replace=T),
    if(prior.weights[7]>0) sample(as.matrix(fitstanP4_N)[,2],count[7],replace=T),
    if(prior.weights[8]>0) sample(as.matrix(fitstanL4_N)[,2],count[8],replace=T),
    # lognormal
    if(prior.weights[9]>0)  sample(as.matrix(fitstanE4_LN)[,2],count[9],replace=T),
    if(prior.weights[10]>0) sample(as.matrix(fitstanIE4_LN)[,2],count[10],replace=T),
    if(prior.weights[11]>0) sample(as.matrix(fitstanH4_LN)[,2],count[11],replace=T),
    if(prior.weights[12]>0) sample(as.matrix(fitstanLN4_LN)[,2],count[12],replace=T),
    if(prior.weights[13]>0) sample(as.matrix(fitstanG4_LN)[,2],count[13],replace=T),
    if(prior.weights[14]>0) sample(as.matrix(fitstanQE4_LN)[,2],count[14],replace=T),
    if(prior.weights[15]>0) sample(as.matrix(fitstanP4_LN)[,2],count[15],replace=T),
    if(prior.weights[16]>0) sample(as.matrix(fitstanL4_LN)[,2],count[16],replace=T)
  ))

  macilp=(quantile(mabmd,pvec))*data$maxD
  names(macilp)=c("BMDL","BMD","BMDU") # on original scale

  if(TRUE %in% (mabmd > data$maxD) && data$maxD > 1){
    mabmd = ifelse(mabmd > data$maxD, data$maxD, mabmd)
    p.msg = 'The model averaged posterior distribution has been truncated at max(Dose)^2'
    warnings('The model averaged posterior distribution has been truncated at max(Dose)^2')
  }else{
    p.msg = ''
  }

  BMDq_ls = (quantile(mabmd, seq(0,1,0.005)))*data$maxD


  ### Model-averaged response per dose level
  dr.MA.ls <- c()
  for(i in 1:length(data$x)){
    dr.MA.ls[i] = weighted.mean(x = c(DRM_E4_N[i],DRM_IE4_N[i],DRM_H4_N[i],DRM_LN4_N[i],DRM_G4_N[i],DRM_QE4_N[i], DRM_P4_N[i],
                                      DRM_L4_N[i] ,
                                      DRM_E4_LN[i], DRM_IE4_LN[i], DRM_H4_LN[i], DRM_LN4_LN[i], DRM_G4_LN[i],DRM_QE4_LN[i],
                                      DRM_P4_LN[i],DRM_L4_LN[i]),
                                w = lpwlp,
                                na.rm = T)
  }


  ##### Weights & MA if one of the models is divergent --> this model gets weight 0
  if((0 %in% converged) && (1 %in% converged)){

    prior.weights.new = prior.weights
    div.models = which(converged == 0)
    prior.weights.new[div.models] = 0
    p.weights.new = prior.weights.new/sum(prior.weights.new==1)

    lpwlp.conv=(p.weights.new*w)/sum(p.weights.new*w)
    # lpwlp = lpwlp*prior.weights

    # the model average posterior as a mixture
    count=round(lpwlp.conv*ndraws)
    mabmd.conv=(c(# normal
      if(p.weights.new[1]>0) sample(as.matrix(fitstanE4_N)[,2],count[1],replace=T),
      if(p.weights.new[2]>0) sample(as.matrix(fitstanIE4_N)[,2],count[2],replace=T),
      if(p.weights.new[3]>0) sample(as.matrix(fitstanH4_N)[,2],count[3],replace=T),
      if(p.weights.new[4]>0) sample(as.matrix(fitstanLN4_N)[,2],count[4],replace=T),
      if(p.weights.new[5]>0) sample(as.matrix(fitstanG4_N)[,2],count[5],replace=T),
      if(p.weights.new[6]>0) sample(as.matrix(fitstanQE4_N)[,2],count[6],replace=T),
      if(p.weights.new[7]>0) sample(as.matrix(fitstanP4_N)[,2],count[7],replace=T),
      if(p.weights.new[8]>0) sample(as.matrix(fitstanL4_N)[,2],count[8],replace=T),
      # lognormal
      if(p.weights.new[9]>0)  sample(as.matrix(fitstanE4_LN)[,2],count[9],replace=T),
      if(p.weights.new[10]>0) sample(as.matrix(fitstanIE4_LN)[,2],count[10],replace=T),
      if(p.weights.new[11]>0) sample(as.matrix(fitstanH4_LN)[,2],count[11],replace=T),
      if(p.weights.new[12]>0) sample(as.matrix(fitstanLN4_LN)[,2],count[12],replace=T),
      if(p.weights.new[13]>0) sample(as.matrix(fitstanG4_LN)[,2],count[13],replace=T),
      if(p.weights.new[14]>0) sample(as.matrix(fitstanQE4_LN)[,2],count[14],replace=T),
      if(p.weights.new[15]>0) sample(as.matrix(fitstanP4_LN)[,2],count[15],replace=T),
      if(p.weights.new[16]>0) sample(as.matrix(fitstanL4_LN)[,2],count[16],replace=T)
    ))

    macilp.conv=(quantile(mabmd.conv,pvec))*data$maxD
    names(macilp.conv)=c("BMDL","BMD","BMDU")

    if(TRUE %in% (mabmd.conv > data$maxD) && data$maxD > 1){
      mabmd.conv = ifelse(mabmd.conv > data$maxD, data$maxD, mabmd.conv)
      p.msg = 'The model averaged posterior distribution has been truncated at max(Dose)^2'
      warnings('The model averaged posterior distribution has been truncated at max(Dose)^2')
    }else{
      p.msg = ''
    }

    BMDq_ls_conv = (quantile(mabmd.conv, seq(0,1,0.005)))*data$maxD


    ### Model-averaged response per dose level
    dr.MA.ls.conv <- c()
    for(i in 1:length(data$x)){
      dr.MA.ls.conv[i] = weighted.mean(x = c(DRM_E4_N[i],DRM_IE4_N[i],DRM_H4_N[i],DRM_LN4_N[i],DRM_G4_N[i],DRM_QE4_N[i],
                                             DRM_P4_N[i], DRM_L4_N[i] ,
                                             DRM_E4_LN[i], DRM_IE4_LN[i], DRM_H4_LN[i], DRM_LN4_LN[i], DRM_G4_LN[i],DRM_QE4_LN[i],
                                             DRM_P4_LN[i],DRM_L4_LN[i]),
                                       w = lpwlp.conv,
                                       na.rm = T)
    }


  }else{
    lpwlp.conv = NULL; macilp.conv = NULL; BMDq_ls_conv = NULL; dr.MA.ls.conv = NULL; mabmd.conv = NA;
  }

  ## Plot with weights bridge sampling
  BMDL = c(BMDL, macib[1]); BMD = c(BMD, macib[2]); BMDU = c(BMDU, macib[3]) # all on original scale

  names(BMDL) <- c("E4_N","IE4_N","H4_N","LN4_N","G4_N","QE4_N","P4_N","L4_N","E4_LN","IE4_LN","H4_LN","LN4_LN","G4_LN","QE4_LN","P4_LN","L4_LN","MA")
  model = c("E4_N","IE4_N","H4_N","LN4_N","G4_N","QE4_N","P4_N","L4_N","E4_LN","IE4_LN","H4_LN","LN4_LN","G4_LN","QE4_LN","P4_LN","L4_LN","MA")
  model = as.factor(model)

  weight = c(rep(0,16),1)
  names(weight) = c("E4_N","IE4_N","H4_N","LN4_N","G4_N","QE4_N","P4_N","L4_N","E4_LN","IE4_LN","H4_LN","LN4_LN","G4_LN","QE4_LN","P4_LN","L4_LN","MA")
  for(i in names(lpwb)){
    weight[names(weight)==i] = lpwb[names(lpwb)==i]
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

  modelnames <- c("E4_N","IE4_N","H4_N","LN4_N","G4_N","QE4_N","P4_N","L4_N","E4_LN","IE4_LN","H4_LN","LN4_LN","G4_LN","QE4_LN","P4_LN","L4_LN")

  ### Some additional checks


  if(macib[2]/macib[1] > 20){
    warning('BMD/BMDL is larger than 20 for bridge sampling')
  }
  if(macib[3]/macib[1] > 50){
    warning('BMDU/BMDL is larger than 50 for bridge sampling')
  }
  if(macib[2] < (data.N$data$x[2]*data.N$data$maxD/10)){
    warning('BMD is 10 times lower than the lowest non-zero dose for bridge sampling')
  }
  if(macilp[2]/macilp[1] > 20){
    warning('BMD/BMDL is larger than 20 for hybrid Laplace')
  }
  if(macilp[3]/macilp[1] > 50){
    warning('BMDU/BMDL is larger than 50 for hybrid Laplace')
  }
  if(macilp[2] < (data.N$data$x[2]*data.N$data$maxD/10)){
    warning('BMD is 10 times lower than the lowest non-zero dose for hybrid Laplace')
  }

  ### best fitting model vs saturated ANOVA model
  best.fit = modelnames[which(weight[1:16] == max(weight[1:16]))][1]

  bfTest <- modelTest(best.fit, data.N, data.LN, get(paste0('fitstan', best.fit)), type = 'MCMC',
                      seed, ndraws, nrchains, nriterations, warmup, delta, treedepth)
  warning(bfTest$warn.bf)

  ret_results <- list(E4_N=E4outNI,IE4_N=IE4outNI,H4_N=H4outNI,LN4_N=LN4outNI,
                      G4_N=G4outNI,QE4_N=QE4outNI,P4_N=P4outNI,L4_N=L4outNI,
                      E4_LN=E4outLNI,IE4_LN=IE4outLNI,H4_LN=H4outLNI,LN4_LN=LN4outLNI,
                      G4_LN=G4outLNI,QE4_LN=QE4outLNI,P4_LN=P4outLNI,L4_LN=L4outLNI,
                      covs = covs, corrs = corrs,
                      weights_bridge_sampling=w.bs,
                      weights_laplace=lpwlp,
                      MA_bridge_sampling=macib,
                      MA_laplace=macilp,
                      llN=llN, llLN=llLN,
                      convergence=converged, bs_weights_conv=w.bs.conv,
                      ls_weights_conv=lpwlp.conv,
                      MA_bs_conv=macib.conv,
                      MA_ls_conv=macilp.conv,
                      MA_post_bs = BMDq_bs,
                      MA_post_ls = BMDq_ls,
                      MA_post_bs_conv = BMDq_bs_conv,
                      MA_post_ls_conv = BMDq_ls_conv,
                      MA_dr_bs = dr.MA.bs,
                      MA_dr_ls = dr.MA.ls,
                      MA_dr_bs_conv = dr.MA.bs.conv,
                      MA_dr_ls_conv = dr.MA.ls.conv,
                      parsN = list(parsE4N, parsIE4N, parsH4N, parsLN4N, parsG4N,
                                   parsQE4N, parsP4N, parsL4N),
                      parsLN = list(parsE4LN, parsIE4LN, parsH4LN, parsLN4LN, parsG4LN,
                                    parsQE4LN, parsP4LN, parsL4LN),
                      BMDMixture = (mabmd)*data$maxD,
                      BMDMixture.conv = (mabmd.conv)*data$maxD,
                      BMDMixtureBS = (mabmd1)*data$maxD,
                      BMDMixture.convBS = mabmd.conv1*data$maxD,
                      divergences = divergences,
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
                      models_included_bridge = modelnames[p.weights > 0],
                      models_included_laplace = modelnames[lpwlp > 0],
                      bf = bfTest$bayesFactor, gof_check = bfTest$warn.bf,
                      means.SM = bfTest$means.SM, parBestFit = bfTest$par.best,
                      BIC.bestfit = bfTest$BIC.bestfit, BIC.SM = bfTest$BIC.SM,
                      shift = data.LN$data$shift,
                      w.msg = w.msg, p.msg = p.msg
  )

  attr(ret_results, "class") <- c("BMADR", "BS")

  return(ret_results)

}


# Sampling based methods: bridge sampling and laplace approximation
##################################################################################################

#' @rdname sampling_MA
#' @export
sampling_MAc=function(data.N,data.LN,prior.weights = rep(1,16),
                      ndraws = 30000,nrchains=3,
                      nriterations=3000,warmup=1000,
                      delta=0.8,treedepth=10,seed=123,pvec=c(0.05,0.5,0.95)){


  # if(data.N$increasing == TRUE){

  prior.weights.orig = prior.weights

  ## Normal distribution

  data=data.N$data
  start=data.N$start
  startQ=data.N$startQ

  BMDL=c(); BMD=c(); BMDU=c()
  converged=c()

  nModels = 16

  ## Obtain model parameters via MCMC sampling
  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[1]))
  }
  if(prior.weights[1]>0){
    # print(1)

    fitstanE4_N = fun_samplingC(stanmodels$mE4c, data, start,
                                ndraws,nrchains,
                                nriterations,warmup,
                                delta,treedepth,seed,pvec)

    if(is.null(fitstanE4_N)){
      prior.weights[1] <- 0
      warning('difficulties fitting the Exponential (Normal) model using sampling; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsE4N <- par_extractC(fitstanE4_N, model_name = "E4_N")

      #diagnostics here
      convergence_stat <- convergence_deciC(model_stan = fitstanE4_N, nrchains = nrchains)
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)

      # divergent transitions
      div_E4_N <- sum(sapply(rstan::get_sampler_params(fitstanE4_N, inc_warmup = F),
                             function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))

      E4resNI=quantile(as.matrix(fitstanE4_N)[,"par2"]*data$maxD,pvec)
      E4resNI=c(E4resNI,apply(as.matrix(fitstanE4_N),2,median)[c("par1","par2","par3","par4","par5","par6")])
      names(E4resNI)=c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t","rho")

      E4outNI <- outLP(parsE4N, pvec, data$maxD, clustered = T)

      if(data$data_type == 1){
        DRM_E4_N = DRM.E4_NI(E4resNI[4:7], data$x, data$q)
      }else if(data$data_type == 3){
        DRM_E4_N = DRM.E4_ND(E4resNI[4:7], data$x, data$q)
      }


      # Covariance between b-d and between BMD-d
      E4covNI = c(cov(as.matrix(fitstanE4_N)[,c("b","d")], use="complete.obs")["b","d"],
                  cov(as.matrix(fitstanE4_N)[,c("par2","d")], use="complete.obs")["par2","d"])

      E4corrNI = c(cor(as.matrix(fitstanE4_N)[,c("b","d")], use="complete.obs")["b","d"],
                   cor(as.matrix(fitstanE4_N)[,c("par2","d")], use="complete.obs")["par2","d"])

      BMDL=c(BMDL, E4resNI[1]); BMD=c(BMD, E4resNI[2]); BMDU=c(BMDU, E4resNI[3])
      bridgeE4N = bridgesampling::bridge_sampler(fitstanE4_N, silent = T) # compute log marginal likelihood
    }}
  if(prior.weights[1] == 0){E4resNI=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeE4N=NA; converged=c(converged, NA); E4covNI=rep(NA,2); E4corrNI=rep(NA,2); DRM_E4_N=rep(NA,length(data$x))
  parsE4N <- NA
  div_E4_N <- NA
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

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[2]))
  }
  if(prior.weights[2]>0){
    # print(2)

    fitstanIE4_N = fun_samplingC(stanmodels$mIE4c, data, start,
                                 ndraws,nrchains,
                                 nriterations,warmup,
                                 delta,treedepth,seed,pvec)

    if(is.null(fitstanIE4_N)){
      prior.weights[2] <- 0
      warning('difficulties fitting the Inverse Exponential (Normal) model using sampling; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsIE4N <- par_extractC(fitstanIE4_N, model_name = "IE4_N")
      # parsIE4N[,"BMD"] <- parsIE4N[,"BMD"]*data$maxD # BMD on original scale

      #diagnostics here
      # posterior_diag(model_stan = fitstanIE4_N)
      convergence_stat <- convergence_deciC(model_stan = fitstanIE4_N, nrchains = nrchains)
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)

      # divergent transitions
      div_IE4_N <- sum(sapply(rstan::get_sampler_params(fitstanIE4_N, inc_warmup = F),
                              function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))

      IE4resNI=(quantile(as.matrix(fitstanIE4_N)[,"par2"],pvec))*data$maxD
      IE4resNI=c(IE4resNI,apply(as.matrix(fitstanIE4_N),2,median)[c("par1","par2","par3","par4","par5","par6")])
      names(IE4resNI)=c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t","rho")

      IE4outNI <- outLP(parsIE4N, pvec, data$maxD, clustered = T)

      if(data$data_type == 1){
        DRM_IE4_N = DRM.IE4_NI(IE4resNI[4:7], data$x, data$q)
      }else if(data$data_type == 3){
        DRM_IE4_N = DRM.IE4_ND(IE4resNI[4:7], data$x, data$q)
      }

      # Covariance between b-d and between BMD-d
      IE4covNI = c(cov(as.matrix(fitstanIE4_N)[,c("b","d")], use="complete.obs")["b","d"],
                   cov(as.matrix(fitstanIE4_N)[,c("par2","d")], use="complete.obs")["par2","d"])

      IE4corrNI = c(cor(as.matrix(fitstanIE4_N)[,c("b","d")], use="complete.obs")["b","d"],
                    cor(as.matrix(fitstanIE4_N)[,c("par2","d")], use="complete.obs")["par2","d"])

      BMDL=c(BMDL, IE4resNI[1]); BMD=c(BMD, IE4resNI[2]); BMDU=c(BMDU, IE4resNI[3])
      bridgeIE4N = bridgesampling::bridge_sampler(fitstanIE4_N, silent = T)

    }}
  if(prior.weights[2] == 0){IE4resNI=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeIE4N=NA; converged=c(converged, NA); IE4covNI=rep(NA,2); IE4corrNI=rep(NA,2); DRM_IE4_N=rep(NA,length(data$x))
  parsIE4N <- NA
  div_IE4_N <- NA

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

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[3]))
  }
  if(prior.weights[3]>0){
    # print(3)

    fitstanH4_N = fun_samplingC(stanmodels$mH4c, data, start,
                                ndraws,nrchains,
                                nriterations,warmup,
                                delta,treedepth,seed,pvec)

    if(is.null(fitstanH4_N)){
      prior.weights[3] <- 0
      warning('difficulties fitting the Hill (Normal) model using sampling; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsH4N <- par_extractC(fitstanH4_N, model_name = "H4_N")
      # parsH4N[,"BMD"] <- parsH4N[,"BMD"]*data$maxD # BMD on original scale

      #diagnostics here
      # posterior_diag(model_stan = fitstanH4_N)
      convergence_stat <- convergence_deciC(model_stan = fitstanH4_N, nrchains = nrchains)
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)

      # divergent transitions
      div_H4_N <- sum(sapply(rstan::get_sampler_params(fitstanH4_N, inc_warmup = F),
                             function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))

      H4resNI=(quantile(as.matrix(fitstanH4_N)[,"par2"],pvec))*data$maxD
      H4resNI=c(H4resNI,apply(as.matrix(fitstanH4_N),2,median)[c("par1","par2","par3","par4","par5","par6")])
      names(H4resNI)=c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t","rho")

      H4outNI <- outLP(parsH4N, pvec, data$maxD, clustered = T)

      if(data$data_type == 1){
        DRM_H4_N = DRM.H4_NI(H4resNI[4:7], data$x, data$q)
      }else if(data$data_type == 3){
        DRM_H4_N = DRM.H4_ND(H4resNI[4:7], data$x, data$q)
      }
      # Covariance between b-d and between BMD-d
      H4covNI = c(cov(as.matrix(fitstanH4_N)[,c("b","d")], use="complete.obs")["b","d"],
                  cov(as.matrix(fitstanH4_N)[,c("par2","d")], use="complete.obs")["par2","d"])

      H4corrNI = c(cor(as.matrix(fitstanH4_N)[,c("b","d")], use="complete.obs")["b","d"],
                   cor(as.matrix(fitstanH4_N)[,c("par2","d")], use="complete.obs")["par2","d"])

      BMDL=c(BMDL, H4resNI[1]); BMD=c(BMD, H4resNI[2]); BMDU=c(BMDU, H4resNI[3])
      bridgeH4N = bridgesampling::bridge_sampler(fitstanH4_N, silent = T)
    }}
  if(prior.weights[3] == 0){H4resNI=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeH4N=NA;
  converged=c(converged, NA); H4covNI=rep(NA,2); H4corrNI=rep(NA,2); DRM_H4_N=rep(NA,length(data$x))
  parsH4N <- NA
  div_H4_N <- NA

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

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[4]))
  }
  if(prior.weights[4]>0){
    # print(4)
    fitstanLN4_N = fun_samplingC(stanmodels$mLN4c, data, start,
                                 ndraws,nrchains,
                                 nriterations,warmup,
                                 delta,treedepth,seed,pvec)

    if(is.null(fitstanLN4_N)){
      prior.weights[4] <- 0
      warning('difficulties fitting the Lognormal (Normal) model using sampling; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsLN4N <- par_extractC(fitstanLN4_N, model_name = "LN4_N")
      # parsLN4N[,"BMD"] <- parsLN4N[,"BMD"]*data$maxD # BMD on original scale

      #diagnostics here
      # posterior_diag(model_stan = fitstanLN4_N)
      convergence_stat <- convergence_deciC(model_stan = fitstanLN4_N, nrchains = nrchains)
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)

      # divergent transitions
      div_LN4_N <- sum(sapply(rstan::get_sampler_params(fitstanLN4_N, inc_warmup = F),
                              function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))

      LN4resNI=(quantile(as.matrix(fitstanLN4_N)[,"par2"],pvec))*data$maxD
      LN4resNI=c(LN4resNI,apply(as.matrix(fitstanLN4_N),2,median)[c("par1","par2","par3","par4","par5","par6")])
      names(LN4resNI)=c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t","rho")

      LN4outNI <- outLP(parsLN4N, pvec, data$maxD, clustered = T)

      if(data$data_type == 1){
        DRM_LN4_N = DRM.LN4_NI(LN4resNI[4:7], data$x, data$q)
      }else if(data$data_type == 3){
        DRM_LN4_N = DRM.LN4_ND(LN4resNI[4:7], data$x, data$q)
      }

      # Covariance between b-d and between BMD-d
      LN4covNI = c(cov(as.matrix(fitstanLN4_N)[,c("b","d")], use="complete.obs")["b","d"],
                   cov(as.matrix(fitstanLN4_N)[,c("par2","d")], use="complete.obs")["par2","d"])

      LN4corrNI = c(cor(as.matrix(fitstanLN4_N)[,c("b","d")], use="complete.obs")["b","d"],
                    cor(as.matrix(fitstanLN4_N)[,c("par2","d")], use="complete.obs")["par2","d"])

      BMDL=c(BMDL, LN4resNI[1]); BMD=c(BMD, LN4resNI[2]); BMDU=c(BMDU, LN4resNI[3])
      bridgeLN4N = bridgesampling::bridge_sampler(fitstanLN4_N, silent = T)
    }}
  if(prior.weights[4] == 0){LN4resNI=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeLN4N=NA; converged=c(converged, NA); LN4covNI=rep(NA,2); LN4corrNI=rep(NA,2); DRM_LN4_N=rep(NA,length(data$x))
  parsLN4N <- NA
  div_LN4_N <- NA

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

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[5]))
  }
  if(prior.weights[5]>0){
    # print(5)
    fitstanG4_N = fun_samplingC(stanmodels$mG4c, data, start,
                                ndraws,nrchains,
                                nriterations,warmup,
                                delta,treedepth,seed,pvec)

    if(is.null(fitstanG4_N)){
      prior.weights[5] <- 0
      warning('difficulties fitting the Gamma (Normal) model using sampling; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsG4N <- par_extractC(fitstanG4_N, model_name = "G4_N")
      # parsG4N[,"BMD"] <- parsG4N[,"BMD"]*data$maxD # BMD on original scale

      #diagnostics here
      # posterior_diag(model_stan = fitstanG4_N)
      convergence_stat <- convergence_deciC(model_stan = fitstanG4_N, nrchains = nrchains)
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)

      # divergent transitions
      div_G4_N <- sum(sapply(rstan::get_sampler_params(fitstanG4_N, inc_warmup = F),
                             function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))

      G4resNI=(quantile(as.matrix(fitstanG4_N)[,"par2"],pvec))*data$maxD
      G4resNI=c(G4resNI,apply(as.matrix(fitstanG4_N),2,median)[c("par1","par2","par3","par4","par5","par6")])
      names(G4resNI)=c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t","rho")

      G4outNI <- outLP(parsG4N, pvec, data$maxD, clustered = T)

      if(data$data_type == 1){
        DRM_G4_N = DRM.G4_NI(G4resNI[4:7], data$x, data$q)
      }else if(data$data_type == 3){
        DRM_G4_N = DRM.G4_ND(G4resNI[4:7], data$x, data$q)
      }

      # Covariance between b-d and between BMD-d
      G4covNI = c(cov(as.matrix(fitstanG4_N)[,c("b","d")], use="complete.obs")["b","d"],
                  cov(as.matrix(fitstanG4_N)[,c("par2","d")], use="complete.obs")["par2","d"])

      G4corrNI = c(cor(as.matrix(fitstanG4_N)[,c("b","d")], use="complete.obs")["b","d"],
                   cor(as.matrix(fitstanG4_N)[,c("par2","d")], use="complete.obs")["par2","d"])

      BMDL=c(BMDL, G4resNI[1]); BMD=c(BMD, G4resNI[2]); BMDU=c(BMDU, G4resNI[3])
      bridgeG4N = bridgesampling::bridge_sampler(fitstanG4_N, silent = T)
    }}
  if(prior.weights[5] == 0){G4resNI=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeG4N=NA; converged=c(converged, NA); G4covNI=rep(NA,2); G4corrNI=rep(NA,2); DRM_G4_N=rep(NA,length(data$x))
  parsG4N <- NA
  div_G4_N <- NA

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

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[6]))
  }
  if(prior.weights[6]>0){
    # print(6)
    fitstanQE4_N = fun_samplingC(stanmodels$mQE4c, data, startQ,
                                 ndraws,nrchains,
                                 nriterations,warmup,
                                 delta,treedepth,seed,pvec)

    if(is.null(fitstanQE4_N)){
      prior.weights[6] <- 0
      warning('difficulties fitting the Quadratic Exponential (Normal) model using sampling; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsQE4N <- par_extractC(fitstanQE4_N, model_name = "QE4_N")
      # parsQE4N[,"BMD"] <- parsQE4N[,"BMD"]*data$maxD # BMD on original scale

      #diagnostics here
      # posterior_diag(model_stan = fitstanQE4_N)
      convergence_stat <- convergence_deciC(model_stan = fitstanQE4_N, nrchains = nrchains)
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)

      # divergent transitions
      div_QE4_N <- sum(sapply(rstan::get_sampler_params(fitstanQE4_N, inc_warmup = F),
                              function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))

      QE4resNI=(quantile(as.matrix(fitstanQE4_N)[,"par2"],pvec))*data$maxD
      QE4resNI=c(QE4resNI,apply(as.matrix(fitstanQE4_N),2,median)[c("par1","par2","par3","par4","par5","par6")])
      names(QE4resNI)=c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t","rho")

      QE4outNI <- outLP(parsQE4N, pvec, data$maxD, clustered = T)

      if(data$data_type == 1){
        DRM_QE4_N = DRM.QE4_NI(QE4resNI[4:7], data$x, data$q)
      }else if(data$data_type == 3){
        DRM_QE4_N = DRM.QE4_ND(QE4resNI[4:7], data$x, data$q)
      }

      # Covariance between b-d and between BMD-d
      QE4covNI = c(cov(as.matrix(fitstanQE4_N)[,c("b","d")], use="complete.obs")["b","d"],
                   cov(as.matrix(fitstanQE4_N)[,c("par2","d")], use="complete.obs")["par2","d"])

      QE4corrNI = c(cor(as.matrix(fitstanQE4_N)[,c("b","d")], use="complete.obs")["b","d"],
                    cor(as.matrix(fitstanQE4_N)[,c("par2","d")], use="complete.obs")["par2","d"])

      BMDL=c(BMDL, QE4resNI[1]); BMD=c(BMD, QE4resNI[2]); BMDU=c(BMDU, QE4resNI[3])
      bridgeQE4N = bridgesampling::bridge_sampler(fitstanQE4_N, silent = T)
    }}
  if(prior.weights[6] == 0){QE4resNI=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeQE4N=NA; converged=c(converged, NA); QE4covNI=rep(NA,2); QE4corrNI=rep(NA,2); DRM_QE4_N=rep(NA,length(data$x))
  parsQE4N <- NA
  div_QE4_N <- NA

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

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[7]))
  }
  if(prior.weights[7]>0){
    # print(7)
    fitstanP4_N = fun_samplingC(stanmodels$mP4c, data, start,
                                ndraws,nrchains,
                                nriterations,warmup,
                                delta,treedepth,seed,pvec)

    if(is.null(fitstanP4_N)){
      prior.weights[7] <- 0
      warning('difficulties fitting the Probit (Normal) model using sampling; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsP4N <- par_extractC(fitstanP4_N, model_name = "P4_N")
      # parsP4N[,"BMD"] <- parsP4N[,"BMD"]*data$maxD # BMD on original scale

      #diagnostics here
      # posterior_diag(model_stan = fitstanP4_N)
      convergence_stat <- convergence_deciC(model_stan = fitstanP4_N, nrchains = nrchains)
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)

      # divergent transitions
      div_P4_N <- sum(sapply(rstan::get_sampler_params(fitstanP4_N, inc_warmup = F),
                             function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))

      P4resNI=(quantile(as.matrix(fitstanP4_N)[,"par2"],pvec))*data$maxD
      P4resNI=c(P4resNI,apply(as.matrix(fitstanP4_N),2,median)[c("par1","par2","par3","par4","par5","par6")])
      names(P4resNI)=c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t","rho")

      P4outNI <- outLP(parsP4N, pvec, data$maxD, clustered = T)

      if(data$data_type == 1){
        DRM_P4_N = DRM.P4_NI(P4resNI[4:7], data$x, data$q)
      }else if(data$data_type == 3){
        DRM_P4_N = DRM.P4_ND(P4resNI[4:7], data$x, data$q)
      }

      # Covariance between b-d and between BMD-d
      P4covNI = c(cov(as.matrix(fitstanP4_N)[,c("b","d")], use="complete.obs")["b","d"],
                  cov(as.matrix(fitstanP4_N)[,c("par2","d")], use="complete.obs")["par2","d"])

      P4corrNI = c(cor(as.matrix(fitstanP4_N)[,c("b","d")], use="complete.obs")["b","d"],
                   cor(as.matrix(fitstanP4_N)[,c("par2","d")], use="complete.obs")["par2","d"])

      BMDL=c(BMDL, P4resNI[1]); BMD=c(BMD, P4resNI[2]); BMDU=c(BMDU, P4resNI[3])
      bridgeP4N = bridgesampling::bridge_sampler(fitstanP4_N, silent = T)
    }}
  if(prior.weights[7] == 0){P4resNI=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeP4N=NA; converged=c(converged, NA);
  P4covNI=rep(NA,2); P4corrNI=rep(NA,2);
  DRM_P4_N=rep(NA,length(data$x))
  parsP4N <- NA
  div_P4_N <- NA

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

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[8]))
  }
  if(prior.weights[8]>0){
    # print(8)
    fitstanL4_N = fun_samplingC(stanmodels$mL4c, data, start,
                                ndraws,nrchains,
                                nriterations,warmup,
                                delta,treedepth,seed,pvec)

    if(is.null(fitstanL4_N)){
      prior.weights[8] <- 0
      warning('difficulties fitting the Logit (Normal) model using sampling; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsL4N <- par_extractC(fitstanL4_N, model_name = "L4_N")
      # parsL4N[,"BMD"] <- parsL4N[,"BMD"]*data$maxD # BMD on original scale

      # save(fitstanL4_N, file = "fitL4N.RData")

      #diagnostics here
      # posterior_diag(model_stan = fitstanL4_N)
      convergence_stat <- convergence_deciC(model_stan = fitstanL4_N, nrchains = nrchains)
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)

      # divergent transitions
      div_L4_N <- sum(sapply(rstan::get_sampler_params(fitstanL4_N, inc_warmup = F),
                             function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))

      L4resNI=(quantile(as.matrix(fitstanL4_N)[,"par2"],pvec))*data$maxD
      L4resNI=c(L4resNI,apply(as.matrix(fitstanL4_N),2,median)[c("par1","par2","par3","par4","par5","par6")])
      names(L4resNI)=c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t","rho")

      L4outNI <- outLP(parsL4N, pvec, data$maxD, clustered = T)


      if(data$data_type == 1){
        DRM_L4_N = DRM.L4_NI(L4resNI[4:7], data$x, data$q)
      }else if(data$data_type == 3){
        DRM_L4_N = DRM.L4_ND(L4resNI[4:7], data$x, data$q)
      }

      # Covariance between b-d and between BMD-d
      L4covNI = c(cov(as.matrix(fitstanL4_N)[,c("b","d")], use="complete.obs")["b","d"],
                  cov(as.matrix(fitstanL4_N)[,c("par2","d")], use="complete.obs")["par2","d"])

      L4corrNI = c(cor(as.matrix(fitstanL4_N)[,c("b","d")], use="complete.obs")["b","d"],
                   cor(as.matrix(fitstanL4_N)[,c("par2","d")], use="complete.obs")["par2","d"])

      BMDL=c(BMDL, L4resNI[1]); BMD=c(BMD, L4resNI[2]); BMDU=c(BMDU, L4resNI[3])
      bridgeL4N = bridgesampling::bridge_sampler(fitstanL4_N, silent = T)
    }}
  if(prior.weights[8] == 0){L4resNI=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeL4N=NA; converged=c(converged, NA); L4covNI=rep(NA,2); L4corrNI=rep(NA,2); DRM_L4_N=rep(NA,length(data$x))
  parsL4N <- NA
  div_L4_N <- NA

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

  ## Lognormal distribution

  data=data.LN$data
  start=data.LN$start
  startQ=data.LN$startQ

  # model specific results
  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[9]))
  }
  if(prior.weights[9]>0){
    # print(9)

    fitstanE4_LN = fun_samplingC(stanmodels$mE4c, data, start,
                                 ndraws,nrchains,
                                 nriterations,warmup,
                                 delta,treedepth,seed,pvec)

    if(is.null(fitstanE4_LN)){
      prior.weights[9] <- 0
      warning('difficulties fitting the Exponential (Lognormal) model using sampling; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsE4LN <- par_extractC(fitstanE4_LN, model_name = "E4_LN")
      # parsE4LN[,"BMD"] <- parsE4LN[,"BMD"]*data$maxD # BMD on original scale

      #diagnostics here
      # posterior_diag(model_stan = fitstanE4_LN)
      convergence_stat <- convergence_deciC(model_stan = fitstanE4_LN, nrchains = nrchains)
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)

      # divergent transitions
      div_E4_LN <- sum(sapply(rstan::get_sampler_params(fitstanE4_LN, inc_warmup = F),
                              function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))

      E4resLNI=(quantile(as.matrix(fitstanE4_LN)[,"par2"],pvec))*data$maxD
      E4resLNI=c(E4resLNI,apply(as.matrix(fitstanE4_LN),2,median)[c("par1","par2","par3","par4","par5","par6")])
      names(E4resLNI)=c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t","rho")

      E4outLNI <- outLP(parsE4LN, pvec, data$maxD, clustered = T)

      if(data$data_type == 2){
        DRM_E4_LN = exp(DRM.E4_LNI(E4resLNI[4:7], data$x, data$q, shift=data$shift))
      }else if(data$data_type == 4){
        DRM_E4_LN = exp(DRM.E4_LND(E4resLNI[4:7], data$x, data$q, shift=data$shift))
      }

      # Covariance between b-d and between BMD-d
      E4covLNI = c(cov(as.matrix(fitstanE4_LN)[,c("b","d")], use="complete.obs")["b","d"],
                   cov(as.matrix(fitstanE4_LN)[,c("par2","d")], use="complete.obs")["par2","d"])

      E4corrLNI = c(cor(as.matrix(fitstanE4_LN)[,c("b","d")], use="complete.obs")["b","d"],
                    cor(as.matrix(fitstanE4_LN)[,c("par2","d")], use="complete.obs")["par2","d"])

      BMDL=c(BMDL, E4resLNI[1]); BMD=c(BMD, E4resLNI[2]); BMDU=c(BMDU, E4resLNI[3])
      bridgeE4LN = bridgesampling::bridge_sampler(fitstanE4_LN, silent = T)
    }}
  if(prior.weights[9] == 0){E4resLNI=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeE4LN=NA; converged=c(converged, NA); E4covLNI=rep(NA,2); E4corrLNI=rep(NA,2); DRM_E4_LN=rep(NA,length(data$x))
  parsE4LN <- NA
  div_E4_LN <- NA

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

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[10]))
  }
  if(prior.weights[10]>0){
    # print(10)

    fitstanIE4_LN = fun_samplingC(stanmodels$mIE4c, data, start,
                                  ndraws,nrchains,
                                  nriterations,warmup,
                                  delta,treedepth,seed,pvec)

    if(is.null(fitstanIE4_LN)){
      prior.weights[10] <- 0
      warning('difficulties fitting the Inverse Exponential (Lognormal) model using sampling; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsIE4LN <- par_extractC(fitstanIE4_LN, model_name = "IE4_LN")
      # parsIE4LN[,"BMD"] <- parsIE4LN[,"BMD"]*data$maxD # BMD on original scale

      #diagnostics here
      # posterior_diag(model_stan = fitstanIE4_LN)
      convergence_stat <- convergence_deciC(model_stan = fitstanIE4_LN, nrchains = nrchains)
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)

      # divergent transitions
      div_IE4_LN <- sum(sapply(rstan::get_sampler_params(fitstanIE4_LN, inc_warmup = F),
                               function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))

      IE4resLNI=(quantile(as.matrix(fitstanIE4_LN)[,"par2"],pvec))*data$maxD
      IE4resLNI=c(IE4resLNI,apply(as.matrix(fitstanIE4_LN),2,median)[c("par1","par2","par3","par4","par5","par6")])
      names(IE4resLNI)=c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t","rho")

      IE4outLNI <- outLP(parsIE4LN, pvec, data$maxD, clustered = T)


      if(data$data_type == 2){
        DRM_IE4_LN = exp(DRM.IE4_LNI(IE4resLNI[4:7], data$x, data$q, shift=data$shift))
      }else if(data$data_type == 4){
        DRM_IE4_LN = exp(DRM.IE4_LND(IE4resLNI[4:7], data$x, data$q, shift=data$shift))
      }

      IE4covLNI = c(cov(as.matrix(fitstanIE4_LN)[,c("b","d")], use="complete.obs")["b","d"],
                    cov(as.matrix(fitstanIE4_LN)[,c("par2","d")], use="complete.obs")["par2","d"])

      IE4corrLNI = c(cor(as.matrix(fitstanIE4_LN)[,c("b","d")], use="complete.obs")["b","d"],
                     cor(as.matrix(fitstanIE4_LN)[,c("par2","d")], use="complete.obs")["par2","d"])

      BMDL=c(BMDL, IE4resLNI[1]); BMD=c(BMD, IE4resLNI[2]); BMDU=c(BMDU, IE4resLNI[3])
      bridgeIE4LN = bridgesampling::bridge_sampler(fitstanIE4_LN, silent = T)
    }}
  if(prior.weights[10] == 0){IE4resLNI=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeIE4LN=NA; converged=c(converged, NA); IE4covLNI=rep(NA,2); IE4corrLNI=rep(NA,2); DRM_IE4_LN=rep(NA,length(data$x))
  parsIE4LN <- NA
  div_IE4_LN <- NA

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

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[11]))
  }
  if(prior.weights[11]>0){
    # print(11)

    fitstanH4_LN = fun_samplingC(stanmodels$mH4c, data, start,
                                 ndraws,nrchains,
                                 nriterations,warmup,
                                 delta,treedepth,seed,pvec)

    if(is.null(fitstanH4_LN)){
      prior.weights[11] <- 0
      warning('difficulties fitting the Hill (Lognormal) model using sampling; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsH4LN <- par_extractC(fitstanH4_LN, model_name = "H4_LN")
      # parsH4LN[,"BMD"] <- parsH4LN[,"BMD"]*data$maxD # BMD on original scale

      #diagnostics here
      # posterior_diag(model_stan = fitstanH4_LN)
      convergence_stat <- convergence_deciC(model_stan = fitstanH4_LN, nrchains = nrchains)
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)

      # divergent transitions
      div_H4_LN <- sum(sapply(rstan::get_sampler_params(fitstanH4_LN, inc_warmup = F),
                              function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))

      H4resLNI=(quantile(as.matrix(fitstanH4_LN)[,"par2"],pvec))*data$maxD
      H4resLNI=c(H4resLNI,apply(as.matrix(fitstanH4_LN),2,median)[c("par1","par2","par3","par4","par5","par6")])
      names(H4resLNI)=c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t","rho")

      H4outLNI <- outLP(parsH4LN, pvec, data$maxD, clustered = T)

      if(data$data_type == 2){
        DRM_H4_LN = exp(DRM.H4_LNI(H4resLNI[4:7], data$x, data$q, shift=data$shift))
      }else if(data$data_type == 4){
        DRM_H4_LN = exp(DRM.H4_LND(H4resLNI[4:7], data$x, data$q, shift=data$shift))
      }

      H4covLNI = c(cov(as.matrix(fitstanH4_LN)[,c("b","d")], use="complete.obs")["b","d"],
                   cov(as.matrix(fitstanH4_LN)[,c("par2","d")], use="complete.obs")["par2","d"])

      H4corrLNI = c(cor(as.matrix(fitstanH4_LN)[,c("b","d")], use="complete.obs")["b","d"],
                    cor(as.matrix(fitstanH4_LN)[,c("par2","d")], use="complete.obs")["par2","d"])

      BMDL=c(BMDL, H4resLNI[1]); BMD=c(BMD, H4resLNI[2]); BMDU=c(BMDU, H4resLNI[3])
      bridgeH4LN = bridgesampling::bridge_sampler(fitstanH4_LN, silent = T)
    }}
  if(prior.weights[11] == 0){H4resLNI=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeH4LN=NA; converged=c(converged, NA); H4covLNI=rep(NA,2); H4corrLNI=rep(NA,2); DRM_H4_LN=rep(NA,length(data$x))
  parsH4LN <- NA
  div_H4_LN <- NA

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

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[12]))
  }
  if(prior.weights[12]>0){
    # print(12)

    fitstanLN4_LN = fun_samplingC(stanmodels$mLN4c, data, start,
                                  ndraws,nrchains,
                                  nriterations,warmup,
                                  delta,treedepth,seed,pvec)

    if(is.null(fitstanLN4_LN)){
      prior.weights[12] <- 0
      warning('difficulties fitting the Lognormal (Lognormal) model using sampling; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsLN4LN <- par_extractC(fitstanLN4_LN, model_name = "LN4_LN")
      # parsLN4LN[,"BMD"] <- parsLN4LN[,"BMD"]*data$maxD # BMD on original scale

      #diagnostics here
      # posterior_diag(model_stan = fitstanLN4_LN)
      convergence_stat <- convergence_deciC(model_stan = fitstanLN4_LN, nrchains = nrchains)
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)

      # divergent transitions
      div_LN4_LN <- sum(sapply(rstan::get_sampler_params(fitstanLN4_LN, inc_warmup = F),
                               function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))

      LN4resLNI=(quantile(as.matrix(fitstanLN4_LN)[,"par2"],pvec))*data$maxD
      LN4resLNI=c(LN4resLNI,apply(as.matrix(fitstanLN4_LN),2,median)[c("par1","par2","par3","par4","par5","par6")])
      names(LN4resLNI)=c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t","rho")

      LN4outLNI <- outLP(parsLN4LN, pvec, data$maxD, clustered = T)


      if(data$data_type == 2){
        DRM_LN4_LN = exp(DRM.LN4_LNI(LN4resLNI[4:7], data$x, data$q, shift=data$shift))
      }else if(data$data_type == 4){
        DRM_LN4_LN = exp(DRM.LN4_LND(LN4resLNI[4:7], data$x, data$q, shift=data$shift))
      }

      LN4covLNI = c(cov(as.matrix(fitstanLN4_LN)[,c("b","d")], use="complete.obs")["b","d"],
                    cov(as.matrix(fitstanLN4_LN)[,c("par2","d")], use="complete.obs")["par2","d"])

      LN4corrLNI = c(cor(as.matrix(fitstanLN4_LN)[,c("b","d")], use="complete.obs")["b","d"],
                     cor(as.matrix(fitstanLN4_LN)[,c("par2","d")], use="complete.obs")["par2","d"])

      BMDL=c(BMDL, LN4resLNI[1]); BMD=c(BMD, LN4resLNI[2]); BMDU=c(BMDU, LN4resLNI[3])
      bridgeLN4LN = bridgesampling::bridge_sampler(fitstanLN4_LN, silent = T)
    }}
  if(prior.weights[12] == 0){LN4resLNI=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeLN4LN=NA; converged=c(converged, NA); LN4covLNI=rep(NA,2); LN4corrLNI=rep(NA,2); DRM_LN4_LN=rep(NA,length(data$x))
  parsLN4LN <- NA
  div_LN4_LN <- NA

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

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[13]))
  }
  if(prior.weights[13]>0){
    # print(13)

    fitstanG4_LN = fun_samplingC(stanmodels$mG4c, data, start,
                                 ndraws,nrchains,
                                 nriterations,warmup,
                                 delta,treedepth,seed,pvec)

    if(is.null(fitstanG4_LN)){
      prior.weights[13] <- 0
      warning('difficulties fitting the Gamma (Lognormal) model using sampling; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsG4LN <- par_extractC(fitstanG4_LN, model_name = "G4_LN")
      # parsG4LN[,"BMD"] <- parsG4LN[,"BMD"]*data$maxD # BMD on original scale


      #diagnostics here
      # posterior_diag(model_stan = fitstanG4_LN)
      convergence_stat <- convergence_deciC(model_stan = fitstanG4_LN, nrchains = nrchains)
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)

      # divergent transitions
      div_G4_LN <- sum(sapply(rstan::get_sampler_params(fitstanG4_LN, inc_warmup = F),
                              function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))

      G4resLNI=(quantile(as.matrix(fitstanG4_LN)[,"par2"],pvec))*data$maxD
      G4resLNI=c(G4resLNI,apply(as.matrix(fitstanG4_LN),2,median)[c("par1","par2","par3","par4","par5","par6")])
      names(G4resLNI)=c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t","rho")

      G4outLNI <- outLP(parsG4LN, pvec, data$maxD, clustered = T)


      if(data$data_type == 2){
        DRM_G4_LN = exp(DRM.G4_LNI(G4resLNI[4:7], data$x, data$q, shift=data$shift))
      }else if(data$data_type == 4){
        DRM_G4_LN = exp(DRM.G4_LND(G4resLNI[4:7], data$x, data$q, shift=data$shift))
      }

      G4covLNI = c(cov(as.matrix(fitstanG4_LN)[,c("b","d")], use="complete.obs")["b","d"],
                   cov(as.matrix(fitstanG4_LN)[,c("par2","d")], use="complete.obs")["par2","d"])

      G4corrLNI = c(cor(as.matrix(fitstanG4_LN)[,c("b","d")], use="complete.obs")["b","d"],
                    cor(as.matrix(fitstanG4_LN)[,c("par2","d")], use="complete.obs")["par2","d"])

      BMDL=c(BMDL, G4resLNI[1]); BMD=c(BMD, G4resLNI[2]); BMDU=c(BMDU, G4resLNI[3])
      bridgeG4LN = bridgesampling::bridge_sampler(fitstanG4_LN, silent = T)
    }}
  if(prior.weights[13] == 0){G4resLNI=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeG4LN=NA; converged=c(converged, NA); G4covLNI=rep(NA,2); G4corrLNI=rep(NA,2); DRM_G4_LN=rep(NA,length(data$x))
  parsG4LN <- NA
  div_G4_LN <- NA

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

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[14]))
  }
  if(prior.weights[14]>0){
    # print(14)

    fitstanQE4_LN = fun_samplingC(stanmodels$mQE4c, data, startQ,
                                  ndraws,nrchains,
                                  nriterations,warmup,
                                  delta,treedepth,seed,pvec)

    if(is.null(fitstanQE4_LN)){
      prior.weights[14] <- 0
      warning('difficulties fitting the Quadratic Exponential (Lognormal) model using sampling; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsQE4LN <- par_extractC(fitstanQE4_LN, model_name = "QE4_LN")
      # parsQE4LN[,"BMD"] <- parsQE4LN[,"BMD"]*data$maxD # BMD on original scale

      #diagnostics here
      # posterior_diag(model_stan = fitstanQE4_LN)
      convergence_stat <- convergence_deciC(model_stan = fitstanQE4_LN, nrchains = nrchains)
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)

      # divergent transitions
      div_QE4_LN <- sum(sapply(rstan::get_sampler_params(fitstanQE4_LN, inc_warmup = F),
                               function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))

      QE4resLNI=(quantile(as.matrix(fitstanQE4_LN)[,"par2"],pvec))*data$maxD
      QE4resLNI=c(QE4resLNI,apply(as.matrix(fitstanQE4_LN),2,median)[c("par1","par2","par3","par4","par5","par6")])
      names(QE4resLNI)=c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t","rho")

      QE4outLNI <- outLP(parsQE4LN, pvec, data$maxD, clustered = T)

      if(data$data_type == 2){
        DRM_QE4_LN = exp(DRM.QE4_LNI(QE4resLNI[4:7], data$x, data$q, shift=data$shift))
      }else if(data$data_type == 4){
        DRM_QE4_LN = exp(DRM.QE4_LND(QE4resLNI[4:7], data$x, data$q, shift=data$shift))
      }

      QE4covLNI = c(cov(as.matrix(fitstanQE4_LN)[,c("b","d")], use="complete.obs")["b","d"],
                    cov(as.matrix(fitstanQE4_LN)[,c("par2","d")], use="complete.obs")["par2","d"])

      QE4corrLNI = c(cor(as.matrix(fitstanQE4_LN)[,c("b","d")], use="complete.obs")["b","d"],
                     cor(as.matrix(fitstanQE4_LN)[,c("par2","d")], use="complete.obs")["par2","d"])

      BMDL=c(BMDL, QE4resLNI[1]); BMD=c(BMD, QE4resLNI[2]); BMDU=c(BMDU, QE4resLNI[3])
      bridgeQE4LN = bridgesampling::bridge_sampler(fitstanQE4_LN, silent = T)
    }}
  if(prior.weights[14] == 0){QE4resLNI=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeQE4LN=NA; converged=c(converged, NA); QE4covLNI=rep(NA,2); QE4corrLNI=rep(NA,2); DRM_QE4_LN=rep(NA,length(data$x))
  parsQE4LN <- NA
  div_QE4_LN <- NA

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

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[15]))
  }
  if(prior.weights[15]>0){
    # print(15)

    fitstanP4_LN = fun_samplingC(stanmodels$mP4c, data, start,
                                 ndraws,nrchains,
                                 nriterations,warmup,
                                 delta,treedepth,seed,pvec)

    if(is.null(fitstanP4_LN)){
      prior.weights[15] <- 0
      warning('difficulties fitting the Probit (Lognormal) model using sampling; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsP4LN <- par_extractC(fitstanP4_LN, model_name = "P4_LN")
      # parsP4LN[,"BMD"] <- parsP4LN[,"BMD"]*data$maxD # BMD on original scale

      #diagnostics here
      # posterior_diag(model_stan = fitstanP4_LN)
      convergence_stat <- convergence_deciC(model_stan = fitstanP4_LN, nrchains = nrchains)
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)

      # divergent transitions
      div_P4_LN <- sum(sapply(rstan::get_sampler_params(fitstanP4_LN, inc_warmup = F),
                              function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))

      P4resLNI=(quantile(as.matrix(fitstanP4_LN)[,"par2"],pvec))*data$maxD
      P4resLNI=c(P4resLNI,apply(as.matrix(fitstanP4_LN),2,median)[c("par1","par2","par3","par4","par5","par6")])
      names(P4resLNI)=c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t","rho")

      P4outLNI <- outLP(parsP4LN, pvec, data$maxD, clustered = T)


      if(data$data_type == 2){
        DRM_P4_LN = exp(DRM.P4_LNI(P4resLNI[4:7], data$x, data$q, shift=data$shift))
      }else if(data$data_type == 4){
        DRM_P4_LN = exp(DRM.P4_LND(P4resLNI[4:7], data$x, data$q, shift=data$shift))
      }

      P4covLNI = c(cov(as.matrix(fitstanP4_LN)[,c("b","d")], use="complete.obs")["b","d"],
                   cov(as.matrix(fitstanP4_LN)[,c("par2","d")], use="complete.obs")["par2","d"])

      P4corrLNI = c(cor(as.matrix(fitstanP4_LN)[,c("b","d")], use="complete.obs")["b","d"],
                    cor(as.matrix(fitstanP4_LN)[,c("par2","d")], use="complete.obs")["par2","d"])

      BMDL=c(BMDL, P4resLNI[1]); BMD=c(BMD, P4resLNI[2]); BMDU=c(BMDU, P4resLNI[3])
      bridgeP4LN = bridge_sampler(fitstanP4_LN, silent = T)
    }}
  if(prior.weights[15] == 0){P4resLNI=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeP4LN=NA;
  converged=c(converged, NA); P4covLNI=rep(NA,2); P4corrLNI=rep(NA,2); DRM_P4_LN=rep(NA,length(data$x))
  parsP4LN <- NA
  div_P4_LN <- NA

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

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[16]))
  }
  if(prior.weights[16]>0){
    # print(16)

    fitstanL4_LN = fun_samplingC(stanmodels$mL4c, data, start,
                                 ndraws,nrchains,
                                 nriterations,warmup,
                                 delta,treedepth,seed,pvec)

    if(is.null(fitstanL4_LN)){
      prior.weights[16] <- 0
      warning('difficulties fitting the Logit (Lognormal) model using sampling; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsL4LN <- par_extractC(fitstanL4_LN, model_name = "L4_LN")
      # parsL4LN[,"BMD"] <- parsL4LN[,"BMD"]*data$maxD # BMD on original scale


      #diagnostics here
      # posterior_diag(model_stan = fitstanL4_LN)
      convergence_stat <- convergence_deciC(model_stan = fitstanL4_LN, nrchains = nrchains)
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)

      # divergent transitions
      div_L4_LN <- sum(sapply(rstan::get_sampler_params(fitstanL4_LN, inc_warmup = F),
                              function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))

      L4resLNI=(quantile(as.matrix(fitstanL4_LN)[,"par2"],pvec))*data$maxD
      L4resLNI=c(L4resLNI,apply(as.matrix(fitstanL4_LN),2,median)[c("par1","par2","par3","par4","par5","par6")])
      names(L4resLNI)=c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t","rho")

      L4outLNI <- outLP(parsL4LN, pvec, data$maxD, clustered = T)


      if(data$data_type == 2){
        DRM_L4_LN = exp(DRM.L4_LNI(L4resLNI[4:7], data$x, data$q, shift=data$shift))
      }else if(data$data_type == 4){
        DRM_L4_LN = exp(DRM.L4_LND(L4resLNI[4:7], data$x, data$q, shift=data$shift))
      }

      L4covLNI = c(cov(as.matrix(fitstanL4_LN)[,c("b","d")], use="complete.obs")["b","d"],
                   cov(as.matrix(fitstanL4_LN)[,c("par2","d")], use="complete.obs")["par2","d"])

      L4corrLNI = c(cor(as.matrix(fitstanL4_LN)[,c("b","d")], use="complete.obs")["b","d"],
                    cor(as.matrix(fitstanL4_LN)[,c("par2","d")], use="complete.obs")["par2","d"])

      BMDL=c(BMDL, L4resLNI[1]); BMD=c(BMD, L4resLNI[2]); BMDU=c(BMDU, L4resLNI[3])
      bridgeL4LN = bridgesampling::bridge_sampler(fitstanL4_LN, silent = T)
    }}
  if(prior.weights[16] == 0){L4resLNI=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeL4LN=NA; converged=c(converged, NA); L4covLNI=rep(NA,2); L4corrLNI=rep(NA,2); DRM_L4_LN=rep(NA,length(data$x))
  parsL4LN <- NA
  div_L4_LN <- NA

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

  divergences <- c(div_E4_N, div_IE4_N, div_H4_N, div_LN4_N, div_G4_N, div_QE4_N, div_P4_N, div_L4_N,
                   div_E4_LN, div_IE4_LN, div_H4_LN, div_LN4_LN, div_G4_LN, div_QE4_LN, div_P4_LN, div_L4_LN)


  #-----------------------------------

  ### Weights based on bridge sampling are obtained by computing posterior model probabilities from the marginal likelihoods
  p.weights = prior.weights/sum(prior.weights==1)

  # lpwb = bridgesampling::post_prob(bridgeE4N, bridgeIE4N, bridgeH4N, bridgeLN4N, bridgeG4N, bridgeQE4N, bridgeP4N, bridgeL4N,
  #                  bridgeE4LN, bridgeIE4LN, bridgeH4LN, bridgeLN4LN, bridgeG4LN, bridgeQE4LN, bridgeP4LN, bridgeL4LN,
  #                  prior_prob = p.weights[p.weights>0],
  #                  model_names = c("E4_N","IE4_N","H4_N","LN4_N","G4_N","QE4_N","P4_N","L4_N",
  #                                  "E4_LN","IE4_LN","H4_LN","LN4_LN","G4_LN","QE4_LN","P4_LN","L4_LN")[p.weights>0])


  bridge.mods = c(if(p.weights[1]>0) bridgeE4N$logml,
                  if(p.weights[2]>0) bridgeIE4N$logml,
                  if(p.weights[3]>0) bridgeH4N$logml,
                  if(p.weights[4]>0) bridgeLN4N$logml,
                  if(p.weights[5]>0) bridgeG4N$logml,
                  if(p.weights[6]>0) bridgeQE4N$logml,
                  if(p.weights[7]>0) bridgeP4N$logml,
                  if(p.weights[8]>0) bridgeL4N$logml,
                  if(p.weights[9]>0) bridgeE4LN$logml,
                  if(p.weights[10]>0) bridgeIE4LN$logml,
                  if(p.weights[11]>0) bridgeH4LN$logml,
                  if(p.weights[12]>0) bridgeLN4LN$logml,
                  if(p.weights[13]>0) bridgeG4LN$logml,
                  if(p.weights[14]>0) bridgeQE4LN$logml,
                  if(p.weights[15]>0) bridgeP4LN$logml,
                  if(p.weights[16]>0) bridgeL4LN$logml)

  model_names = c("E4_N","IE4_N","H4_N","LN4_N","G4_N","QE4_N","P4_N","L4_N",
                  "E4_LN","IE4_LN","H4_LN","LN4_LN","G4_LN","QE4_LN","P4_LN","L4_LN")[p.weights>0]

  lpwb = bridgesampling::post_prob(bridge.mods,
                                   prior_prob = p.weights[p.weights>0],
                                   model_names = model_names)

  # # the model average posterior as a mixture
  count=round(lpwb*ndraws)
  names(count) = names(lpwb)
  mabmd1=(c( # normal
    if("E4_N" %in% names(count)) sample(as.matrix(fitstanE4_N)[,2],count[names(count)=="E4_N"],replace=T),
    if("IE4_N" %in% names(count)) sample(as.matrix(fitstanIE4_N)[,2],count[names(count)=="IE4_N"],replace=T),
    if("H4_N" %in% names(count)) sample(as.matrix(fitstanH4_N)[,2],count[names(count)=="H4_N"],replace=T),
    if("LN4_N" %in% names(count)) sample(as.matrix(fitstanLN4_N)[,2],count[names(count)=="LN4_N"],replace=T),
    if("G4_N" %in% names(count)) sample(as.matrix(fitstanG4_N)[,2],count[names(count)=="G4_N"],replace=T),
    if("QE4_N" %in% names(count)) sample(as.matrix(fitstanQE4_N)[,2],count[names(count)=="QE4_N"],replace=T),
    if("P4_N" %in% names(count)) sample(as.matrix(fitstanP4_N)[,2],count[names(count)=="P4_N"],replace=T),
    if("L4_N" %in% names(count)) sample(as.matrix(fitstanL4_N)[,2],count[names(count)=="L4_N"],replace=T),
    # lognormal
    if("E4_LN" %in% names(count)) sample(as.matrix(fitstanE4_LN)[,2],count[names(count)=="E4_LN"],replace=T),
    if("IE4_LN" %in% names(count)) sample(as.matrix(fitstanIE4_LN)[,2],count[names(count)=="IE4_LN"],replace=T),
    if("H4_LN" %in% names(count)) sample(as.matrix(fitstanH4_LN)[,2],count[names(count)=="H4_LN"],replace=T),
    if("LN4_LN" %in% names(count)) sample(as.matrix(fitstanLN4_LN)[,2],count[names(count)=="LN4_LN"],replace=T),
    if("G4_LN" %in% names(count)) sample(as.matrix(fitstanG4_LN)[,2],count[names(count)=="G4_LN"],replace=T),
    if("QE4_LN" %in% names(count)) sample(as.matrix(fitstanQE4_LN)[,2],count[names(count)=="QE4_LN"],replace=T),
    if("P4_LN" %in% names(count)) sample(as.matrix(fitstanP4_LN)[,2],count[names(count)=="P4_LN"],replace=T),
    if("L4_LN" %in% names(count)) sample(as.matrix(fitstanL4_LN)[,2],count[names(count)=="L4_LN"],replace=T)
  ))

  macib=(quantile(mabmd1,pvec))*data$maxD
  names(macib)=c("BMDL","BMD","BMDU") # original scale

  if(TRUE %in% (mabmd1 > data$maxD)  && data$maxD > 1){
    mabmd1 = ifelse(mabmd1 > data$maxD, data$maxD, mabmd1)
    p.msg = 'The model averaged posterior distribution has been truncated at max(Dose)^2'
    warnings('The model averaged posterior distribution has been truncated at max(Dose)^2')
  }else{
    p.msg = ''
  }

  BMDq_bs = (quantile(mabmd1, seq(0,1,0.005)))*data$maxD


  mods = c("E4_N","IE4_N","H4_N","LN4_N","G4_N","QE4_N","P4_N","L4_N",
           "E4_LN","IE4_LN","H4_LN","LN4_LN","G4_LN","QE4_LN","P4_LN","L4_LN")
  w.bs = rep(0, 16)
  names(w.bs) = mods
  for(i in mods[mods%in%names(lpwb)]){
    w.bs[i] = lpwb[i]
  }

  ### Model-averaged response per dose level
  dr.MA.bs <- c()
  for(i in 1:length(data$x)){
    dr.MA.bs[i] = weighted.mean(x = c(DRM_E4_N[i],DRM_IE4_N[i],DRM_H4_N[i],DRM_LN4_N[i],DRM_G4_N[i],DRM_QE4_N[i],
                                      DRM_P4_N[i], DRM_L4_N[i] ,
                                      DRM_E4_LN[i], DRM_IE4_LN[i], DRM_H4_LN[i], DRM_LN4_LN[i], DRM_G4_LN[i],
                                      DRM_QE4_LN[i], DRM_P4_LN[i],DRM_L4_LN[i]),
                                w = w.bs,
                                na.rm = T)

  }


  ##### Weights & MA if one of the models is divergent --> this model gets weight 0
  if((0 %in% converged) && (1 %in% converged)){
    prior.weights.new = prior.weights
    div.models = which(converged == 0)
    prior.weights.new[div.models] = 0
    p.weights.new = prior.weights.new/sum(prior.weights.new==1)

    bridge.mods = c(if(p.weights.new[1]>0) bridgeE4N$logml,
                    if(p.weights.new[2]>0) bridgeIE4N$logml,
                    if(p.weights.new[3]>0) bridgeH4N$logml,
                    if(p.weights.new[4]>0) bridgeLN4N$logml,
                    if(p.weights.new[5]>0) bridgeG4N$logml,
                    if(p.weights.new[6]>0) bridgeQE4N$logml,
                    if(p.weights.new[7]>0) bridgeP4N$logml,
                    if(p.weights.new[8]>0) bridgeL4N$logml,
                    if(p.weights.new[9]>0) bridgeE4LN$logml,
                    if(p.weights.new[10]>0) bridgeIE4LN$logml,
                    if(p.weights.new[11]>0) bridgeH4LN$logml,
                    if(p.weights.new[12]>0) bridgeLN4LN$logml,
                    if(p.weights.new[13]>0) bridgeG4LN$logml,
                    if(p.weights.new[14]>0) bridgeQE4LN$logml,
                    if(p.weights.new[15]>0) bridgeP4LN$logml,
                    if(p.weights.new[16]>0) bridgeL4LN$logml)

    model_names = c("E4_N","IE4_N","H4_N","LN4_N","G4_N","QE4_N","P4_N","L4_N",
                    "E4_LN","IE4_LN","H4_LN","LN4_LN","G4_LN","QE4_LN","P4_LN","L4_LN")[p.weights.new>0]

    lpwb.conv = bridgesampling::post_prob(bridge.mods,
                                          prior_prob = p.weights.new[p.weights.new>0],
                                          model_names = model_names)

    # # the model average posterior as a mixture
    count=round(lpwb.conv*ndraws)
    names(count) = names(lpwb.conv)
    mabmd.conv1=(c( # normal
      if("E4_N" %in% names(count)) sample(as.matrix(fitstanE4_N)[,2],count[names(count)=="E4_N"],replace=T),
      if("IE4_N" %in% names(count)) sample(as.matrix(fitstanIE4_N)[,2],count[names(count)=="IE4_N"],replace=T),
      if("H4_N" %in% names(count)) sample(as.matrix(fitstanH4_N)[,2],count[names(count)=="H4_N"],replace=T),
      if("LN4_N" %in% names(count)) sample(as.matrix(fitstanLN4_N)[,2],count[names(count)=="LN4_N"],replace=T),
      if("G4_N" %in% names(count)) sample(as.matrix(fitstanG4_N)[,2],count[names(count)=="G4_N"],replace=T),
      if("QE4_N" %in% names(count)) sample(as.matrix(fitstanQE4_N)[,2],count[names(count)=="QE4_N"],replace=T),
      if("P4_N" %in% names(count)) sample(as.matrix(fitstanP4_N)[,2],count[names(count)=="P4_N"],replace=T),
      if("L4_N" %in% names(count)) sample(as.matrix(fitstanL4_N)[,2],count[names(count)=="L4_N"],replace=T),
      # lognormal
      if("E4_LN" %in% names(count)) sample(as.matrix(fitstanE4_LN)[,2],count[names(count)=="E4_LN"],replace=T),
      if("IE4_LN" %in% names(count)) sample(as.matrix(fitstanIE4_LN)[,2],count[names(count)=="IE4_LN"],replace=T),
      if("H4_LN" %in% names(count)) sample(as.matrix(fitstanH4_LN)[,2],count[names(count)=="H4_LN"],replace=T),
      if("LN4_LN" %in% names(count)) sample(as.matrix(fitstanLN4_LN)[,2],count[names(count)=="LN4_LN"],replace=T),
      if("G4_LN" %in% names(count)) sample(as.matrix(fitstanG4_LN)[,2],count[names(count)=="G4_LN"],replace=T),
      if("QE4_LN" %in% names(count)) sample(as.matrix(fitstanQE4_LN)[,2],count[names(count)=="QE4_LN"],replace=T),
      if("P4_LN" %in% names(count)) sample(as.matrix(fitstanP4_LN)[,2],count[names(count)=="P4_LN"],replace=T),
      if("L4_LN" %in% names(count)) sample(as.matrix(fitstanL4_LN)[,2],count[names(count)=="L4_LN"],replace=T)
    ))

    macib.conv=(quantile(mabmd.conv1,pvec))*data$maxD
    names(macib.conv)=c("BMDL","BMD","BMDU")

    if(TRUE %in% (mabmd.conv1 > data$maxD) && data$maxD > 1){
      mabmd.conv1 = ifelse(mabmd.conv1 > data$maxD, data$maxD, mabmd.conv1)
      p.msg = 'The model averaged posterior distribution has been truncated at max(Dose)^2'
      warnings('The model averaged posterior distribution has been truncated at max(Dose)^2')
    }else{
      p.msg = ''
    }

    BMDq_bs_conv = (quantile(mabmd.conv1, seq(0,1,0.005)))*data$maxD


    mods = c("E4_N","IE4_N","H4_N","LN4_N","G4_N","QE4_N","P4_N","L4_N",
             "E4_LN","IE4_LN","H4_LN","LN4_LN","G4_LN","QE4_LN","P4_LN","L4_LN")
    w.bs.conv = rep(0, 16)
    names(w.bs.conv) = mods
    for(i in mods[mods%in%names(lpwb.conv)]){
      w.bs.conv[i] = lpwb.conv[i]
    }

    ### Model-averaged response per dose level
    dr.MA.bs.conv <- c()
    for(i in 1:length(data$x)){
      dr.MA.bs.conv[i] = weighted.mean(x = c(DRM_E4_N[i],DRM_IE4_N[i],DRM_H4_N[i],DRM_LN4_N[i],DRM_G4_N[i],
                                             DRM_QE4_N[i], DRM_P4_N[i], DRM_L4_N[i] ,
                                             DRM_E4_LN[i], DRM_IE4_LN[i], DRM_H4_LN[i], DRM_LN4_LN[i],
                                             DRM_G4_LN[i],DRM_QE4_LN[i], DRM_P4_LN[i],DRM_L4_LN[i]),
                                       w = w.bs.conv,
                                       na.rm = T)
    }

  }else{
    w.bs.conv = NULL; macib.conv = NULL; mabmd.conv1 = NA; BMDq_bs_conv = NULL; dr.MA.bs.conv = NULL
  }


  #----------------------------------
  ### weights based on laplace approximation

  prior.weights = prior.weights.orig

  # normal
  data=data.N$data
  start=data.N$start
  startQ=data.N$startQ

  llN = c() # likelihoods

  # getting the posterior modes and the hessian and the model specific posterior distributions
  if(prior.weights[1]>0){

    optE4_NI <- fun_optimC(stanmodels$mE4c, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optE4_NI[[3]]),TRUE,(optE4_NI[[3]]!=0)) | length(optE4_NI)!=9)){
      prior.weights[1] <- 0
      warning('difficulties fitting the Exponential (Normal) model using Laplace; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      if(data$data_type == 1){
        llE4N=llfE4_NIc(optE4_NI$par[c(1,2,10,4,5,6)],
                        d=data$x,
                        n=data$n,
                        nij=data$nij,
                        y=data$y,
                        qval=data$q)
      }else if(data$data_type == 3){
        llE4N=llfE4_NDc(optE4_NI$par[c(1,2,10,4,5,6)],
                        d=data$x,
                        n=data$n,
                        nij=data$nij,
                        y=data$y,
                        qval=data$q)
      }


      llN = c(llN, llE4N)
    }}
  if(prior.weights[1] == 0){llN = c(llN,NA)}
  #
  if(prior.weights[2]>0){

    optIE4_NI <- fun_optimC(stanmodels$mIE4c, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optIE4_NI[[3]]),TRUE,(optIE4_NI[[3]]!=0)) | length(optIE4_NI)!=9)){
      prior.weights[2] <- 0
      warning('difficulties fitting the Inverse Exponential (Normal) model using Laplace; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      if(data$data_type == 1){
        llIE4N=llfIE4_NIc(optIE4_NI$par[c(1,2,10,4,5,6)],
                          d=data$x,
                          n=data$n,
                          nij=data$nij,
                          y=data$y,
                          qval=data$q)
      }else if(data$data_type == 3){
        llIE4N=llfIE4_NDc(optIE4_NI$par[c(1,2,10,4,5,6)],
                          d=data$x,
                          n=data$n,
                          nij=data$nij,
                          y=data$y,
                          qval=data$q)
      }

      llN = c(llN, llIE4N)
    }}
  if(prior.weights[2] == 0){llN = c(llN,NA)}
  #
  if(prior.weights[3]>0){

    optH4_NI <- fun_optimC(stanmodels$mH4c, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optH4_NI[[3]]),TRUE,(optH4_NI[[3]]!=0)) | length(optH4_NI)!=9)){
      prior.weights[3] <- 0
      warning('difficulties fitting the Hill (Normal) model using Laplace; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      if(data$data_type == 1){
        llH4N=llfH4_NIc(optH4_NI$par[c(1,2,10,4,5,6)],
                        d=data$x,
                        n=data$n,
                        nij=data$nij,
                        y=data$y,
                        qval=data$q)
      }else if(data$data_type == 3){
        llH4N=llfH4_NDc(optH4_NI$par[c(1,2,10,4,5,6)],
                        d=data$x,
                        n=data$n,
                        nij=data$nij,
                        y=data$y,
                        qval=data$q)
      }

      llN = c(llN, llH4N)
    }}
  if(prior.weights[3] == 0){llN = c(llN,NA)}
  #
  if(prior.weights[4]>0){

    optLN4_NI <- fun_optimC(stanmodels$mLN4c, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optLN4_NI[[3]]),TRUE,(optLN4_NI[[3]]!=0)) | length(optLN4_NI)!=9)){
      prior.weights[4] <- 0
      warning('difficulties fitting the Lognormal (Normal) model using Laplace; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      if(data$data_type == 1){
        llLN4N=llfLN4_NIc(optLN4_NI$par[c(1,2,10,4,5,6)],
                          d=data$x,
                          n=data$n,
                          nij=data$nij,
                          y=data$y,
                          qval=data$q)
      }else if(data$data_type == 3){
        llLN4N=llfLN4_NDc(optLN4_NI$par[c(1,2,10,4,5,6)],
                          d=data$x,
                          n=data$n,
                          nij=data$nij,
                          y=data$y,
                          qval=data$q)
      }


      llN = c(llN, llLN4N)
    }}
  if(prior.weights[4] == 0){llN = c(llN,NA)}
  #
  if(prior.weights[5]>0){

    optG4_NI <- fun_optimC(stanmodels$mG4c, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optG4_NI[[3]]),TRUE,(optG4_NI[[3]]!=0)) | length(optG4_NI)!=9)){
      prior.weights[5] <- 0
      warning('difficulties fitting the Gamma (Normal) model using Laplace; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      if(data$data_type == 1){
        llG4N=llfG4_NIc(optG4_NI$par[c(1,2,10,4,5,6)],
                        d=data$x,
                        n=data$n,
                        nij=data$nij,
                        y=data$y,
                        qval=data$q)
      }else if(data$data_type == 3){
        llG4N=llfG4_NDc(optG4_NI$par[c(1,2,10,4,5,6)],
                        d=data$x,
                        n=data$n,
                        nij=data$nij,
                        y=data$y,
                        qval=data$q)
      }


      llN = c(llN, llG4N)
    }}
  if(prior.weights[5] == 0){llN = c(llN,NA)}
  #
  if(prior.weights[6]>0){

    optQE4_NI <- fun_optimC(stanmodels$mQE4c, data, startQ, ndraws, 123, pvec)

    if((ifelse(is.na(optQE4_NI[[3]]),TRUE,(optQE4_NI[[3]]!=0)) | length(optQE4_NI)!=9)){
      prior.weights[6] <- 0
      warning('difficulties fitting the Quadratic Exponential (Normal) model using Laplace; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      if(data$data_type == 1){
        llQE4N=llfQE4_NIc(optQE4_NI$par[c(1,2,10,4,5,6)],
                          d=data$x,
                          n=data$n,
                          nij=data$nij,
                          y=data$y,
                          qval=data$q)
      }else if(data$data_type == 3){
        llQE4N=llfQE4_NDc(optQE4_NI$par[c(1,2,10,4,5,6)],
                          d=data$x,
                          n=data$n,
                          nij=data$nij,
                          y=data$y,
                          qval=data$q)
      }


      llN = c(llN, llQE4N)
    }}
  if(prior.weights[6] == 0){llN = c(llN,NA)}
  #
  if(prior.weights[7]>0){

    optP4_NI <- fun_optimC(stanmodels$mP4c, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optP4_NI[[3]]),TRUE,(optP4_NI[[3]]!=0)) | length(optP4_NI)!=9)){
      prior.weights[7] <- 0
      warning('difficulties fitting the Probit (Normal) model using Laplace; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      if(data$data_type == 1){
        llP4N=llfP4_NIc(optP4_NI$par[c(1,2,10,4,5,6)],
                        d=data$x,
                        n=data$n,
                        nij=data$nij,
                        y=data$y,
                        qval=data$q)
      }else if(data$data_type == 3){
        llP4N=llfP4_NDc(optP4_NI$par[c(1,2,10,4,5,6)],
                        d=data$x,
                        n=data$n,
                        nij=data$nij,
                        y=data$y,
                        qval=data$q)
      }

      llN = c(llN, llP4N)
    }}
  if(prior.weights[7] == 0){llN = c(llN,NA)}
  #
  if(prior.weights[8]>0){

    optL4_NI <- fun_optimC(stanmodels$mL4c, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optL4_NI[[3]]),TRUE,(optL4_NI[[3]]!=0)) | length(optL4_NI)!=9)){
      prior.weights[8] <- 0
      warning('difficulties fitting the Logit (Normal) model using Laplace; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      if(data$data_type == 1){
        llL4N=llfL4_NIc(optL4_NI$par[c(1,2,10,4,5,6)],
                        d=data$x,
                        n=data$n,
                        nij=data$nij,
                        y=data$y,
                        qval=data$q)
      }else if(data$data_type == 3){
        llL4N=llfL4_NDc(optL4_NI$par[c(1,2,10,4,5,6)],
                        d=data$x,
                        n=data$n,
                        nij=data$nij,
                        y=data$y,
                        qval=data$q)
      }


      llN = c(llN, llL4N)
    }}
  if(prior.weights[8] == 0){llN = c(llN,NA)}


  # lognormal
  data=data.LN$data
  start=data.LN$start
  startQ=data.LN$startQ

  llLN = c()

  if(prior.weights[9]>0){

    optE4_LNI <- fun_optimC(stanmodels$mE4c, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optE4_LNI[[3]]),TRUE,(optE4_LNI[[3]]!=0)) | length(optE4_LNI)!=9)){
      prior.weights[9] <- 0
      warning('difficulties fitting the Exponential (Lognormal) model using Laplace; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      if(data$data_type == 2){
        llE4LN=llfE4_LNIc(optE4_LNI$par[c(1,2,10,4,5,6)],
                          d=data$x,
                          n=data$n,
                          nij=data$nij,
                          y=data$y,
                          qval=data$q,
                          shift=data$shift)
      }else if(data$data_type == 4){
        llE4LN=llfE4_LNDc(optE4_LNI$par[c(1,2,10,4,5,6)],
                          d=data$x,
                          n=data$n,
                          nij=data$nij,
                          y=data$y,
                          qval=data$q,
                          shift=data$shift)
      }


      llLN = c(llLN, llE4LN)
    }}
  if(prior.weights[9] == 0){llLN = c(llLN,NA)}
  #
  if(prior.weights[10]>0){

    optIE4_LNI <- fun_optimC(stanmodels$mIE4c, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optIE4_LNI[[3]]),TRUE,(optIE4_LNI[[3]]!=0)) | length(optIE4_LNI)!=9)){
      prior.weights[10] <- 0
      warning('difficulties fitting the Inverse Exponential (Lognormal) model using Laplace; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      if(data$data_type == 2){
        llIE4LN=llfIE4_LNIc(optIE4_LNI$par[c(1,2,10,4,5,6)],
                            d=data$x,
                            n=data$n,
                            nij=data$nij,
                            y=data$y,
                            qval=data$q,
                            shift=data$shift)
      }else if(data$data_type == 4){
        llIE4LN=llfIE4_LNDc(optIE4_LNI$par[c(1,2,10,4,5,6)],
                            d=data$x,
                            n=data$n,
                            nij=data$nij,
                            y=data$y,
                            qval=data$q,
                            shift=data$shift)
      }


      llLN = c(llLN, llIE4LN)
    }}
  if(prior.weights[10] == 0){llLN = c(llLN,NA)}
  #
  if(prior.weights[11]>0){

    optH4_LNI <- fun_optimC(stanmodels$mH4c, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optH4_LNI[[3]]),TRUE,(optH4_LNI[[3]]!=0)) | length(optH4_LNI)!=9)){
      prior.weights[11] <- 0
      warning('difficulties fitting the Hill (Lognormal) model using Laplace; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      if(data$data_type == 2){
        llH4LN=llfH4_LNIc(optH4_LNI$par[c(1,2,10,4,5,6)],
                          d=data$x,
                          n=data$n,
                          nij=data$nij,
                          y=data$y,
                          qval=data$q,
                          shift=data$shift)
      }else if(data$data_type == 4){
        llH4LN=llfH4_LNDc(optH4_LNI$par[c(1,2,10,4,5,6)],
                          d=data$x,
                          n=data$n,
                          nij=data$nij,
                          y=data$y,
                          qval=data$q,
                          shift=data$shift)
      }

      llLN = c(llLN, llH4LN)
    }}
  if(prior.weights[11] == 0){llLN = c(llLN,NA)}
  #
  if(prior.weights[12]>0){

    optLN4_LNI <- fun_optimC(stanmodels$mLN4c, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optLN4_LNI[[3]]),TRUE,(optLN4_LNI[[3]]!=0)) | length(optLN4_LNI)!=9)){
      prior.weights[12] <- 0
      warning('difficulties fitting the Lognormal (Lognormal) model using Laplace; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      if(data$data_type == 2){
        llLN4LN=llfLN4_LNIc(optLN4_LNI$par[c(1,2,10,4,5,6)],
                            d=data$x,
                            n=data$n,
                            nij=data$nij,
                            y=data$y,
                            qval=data$q,
                            shift=data$shift)
      }else if(data$data_type == 4){
        llLN4LN=llfLN4_LNDc(optLN4_LNI$par[c(1,2,10,4,5,6)],
                            d=data$x,
                            n=data$n,
                            nij=data$nij,
                            y=data$y,
                            qval=data$q,
                            shift=data$shift)
      }

      llLN = c(llLN, llLN4LN)
    }}
  if(prior.weights[12] == 0){llLN = c(llLN,NA)}
  #
  if(prior.weights[13]>0){

    optG4_LNI <- fun_optimC(stanmodels$mG4c, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optG4_LNI[[3]]),TRUE,(optG4_LNI[[3]]!=0)) | length(optG4_LNI)!=9)){
      prior.weights[13] <- 0
      warning('difficulties fitting the Gamma (Lognormal) model using Laplace; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      if(data$data_type == 2){
        llG4LN=llfG4_LNIc(optG4_LNI$par[c(1,2,10,4,5,6)],
                          d=data$x,
                          n=data$n,
                          nij=data$nij,
                          y=data$y,
                          qval=data$q,
                          shift=data$shift)
      }else if(data$data_type == 4){
        llG4LN=llfG4_LNDc(optG4_LNI$par[c(1,2,10,4,5,6)],
                          d=data$x,
                          n=data$n,
                          nij=data$nij,
                          y=data$y,
                          qval=data$q,
                          shift=data$shift)
      }

      llLN = c(llLN, llG4LN)
    }}
  if(prior.weights[13] == 0){llLN = c(llLN,NA)}
  #
  if(prior.weights[14]>0){

    optQE4_LNI <- fun_optimC(stanmodels$mQE4c, data, startQ, ndraws, 123, pvec)

    if((ifelse(is.na(optQE4_LNI[[3]]),TRUE,(optQE4_LNI[[3]]!=0)) | length(optQE4_LNI)!=9)){
      prior.weights[14] <- 0
      warning('difficulties fitting the Quadratic Exponential (Lognormal) model using Laplace; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      if(data$data_type == 2){
        llQE4LN=llfQE4_LNIc(optQE4_LNI$par[c(1,2,10,4,5,6)],
                            d=data$x,
                            n=data$n,
                            nij=data$nij,
                            y=data$y,
                            qval=data$q,
                            shift=data$shift)
      }else if(data$data_type == 4){
        llQE4LN=llfQE4_LNDc(optQE4_LNI$par[c(1,2,10,4,5,6)],
                            d=data$x,
                            n=data$n,
                            nij=data$nij,
                            y=data$y,
                            qval=data$q,
                            shift=data$shift)
      }

      llLN = c(llLN, llQE4LN)
    }}
  if(prior.weights[14] == 0){llLN = c(llLN,NA)}
  #
  if(prior.weights[15]>0){

    optP4_LNI <- fun_optimC(stanmodels$mP4c, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optP4_LNI[[3]]),TRUE,(optP4_LNI[[3]]!=0)) | length(optP4_LNI)!=9)){
      prior.weights[15] <- 0
      warning('difficulties fitting the Probit (Lognormal) model using Laplace; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      if(data$data_type == 2){
        llP4LN=llfP4_LNIc(optP4_LNI$par[c(1,2,10,4,5,6)],
                          d=data$x,
                          n=data$n,
                          nij=data$nij,
                          y=data$y,
                          qval=data$q,
                          shift=data$shift)
      }else if(data$data_type == 4){
        llP4LN=llfP4_LNDc(optP4_LNI$par[c(1,2,10,4,5,6)],
                          d=data$x,
                          n=data$n,
                          nij=data$nij,
                          y=data$y,
                          qval=data$q,
                          shift=data$shift)
      }

      llLN = c(llLN, llP4LN)
    }}
  if(prior.weights[15] == 0){llLN = c(llLN,NA)}
  #
  if(prior.weights[16]>0){

    optL4_LNI <- fun_optimC(stanmodels$mL4c, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optL4_LNI[[3]]),TRUE,(optL4_LNI[[3]]!=0)) | length(optL4_LNI)!=9)){
      prior.weights[16] <- 0
      warning('difficulties fitting the Logit (Lognormal) model using Laplace; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      if(data$data_type == 2){
        llL4LN=llfL4_LNIc(optL4_LNI$par[c(1,2,10,4,5,6)],
                          d=data$x,
                          n=data$n,
                          nij=data$nij,
                          y=data$y,
                          qval=data$q,
                          shift=data$shift)
      }else if(data$data_type == 4){
        llL4LN=llfL4_LNDc(optL4_LNI$par[c(1,2,10,4,5,6)],
                          d=data$x,
                          n=data$n,
                          nij=data$nij,
                          y=data$y,
                          qval=data$q,
                          shift=data$shift)
      }

      llLN = c(llLN, llL4LN)
    }}
  if(prior.weights[16] == 0){llLN = c(llLN,NA)}

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
  if(is.na(lls[which((max.ll-lls[!is.na(lls)]) < 709 & prior.weights>0)][1])){
    lpwlp <- rep(0, 16)
    lpwlp[which(lls == max.ll)] <- 1
    w.msg <- 'Laplace weights could not be computed and one model gets all the weight; using another prior for parameter d might help'
    warning('Laplace weights could not be computed and one model gets all the weight; using another prior for parameter d might help')
  }else{

    if(FALSE %in% ((max.ll-lls[!is.na(lls)]) < 709)){
      w.msg <- 'not all models were used in computation of Laplace weights, some models set to 0; using another prior for parameter d might help'
      warning('not all models were used in computation of Laplace weights, some models set to 0; using another prior for parameter d might help')
    }

    minll <- min(lls[which((max.ll-lls[!is.na(lls)]) < 709 & prior.weights>0)], na.rm = T)

    if(prior.weights[1]>0){
      DIHE4h=det(-solve(optE4_NI$hessian))
      DIHE4=ifelse(DIHE4h<0,0,DIHE4h)
      w = c(w, fun.w(DIHE4, llE4N, minll, optE4_NI, data$priormu, data$priorSigma,
                     data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c,
                     data$truncd))
      # DIHE4h=det(-solve(hessian(func = llE4fN,x=optE4_NI$par[c(1,2,9,3,4)],method.args=list(eps=1e-4, d=0.1, zero.tol=sqrt(.Machine$double.eps/7e-7), r=8, v=2, show.details=FALSE))))
      # DIHE4=ifelse(DIHE4h<0,0,DIHE4h)
      ## Aproximation of marginal (i.e. integrated) likelihood (= 'model evidence')

    }else{w=c(w,0)}

    if(prior.weights[2]>0){
      DIHIE4h=det(-solve(optIE4_NI$hessian))
      DIHIE4=ifelse(DIHIE4h<0,0,DIHIE4h)
      w = c(w, fun.w(DIHIE4, llIE4N, minll, optIE4_NI, data$priormu, data$priorSigma,
                     data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c,
                     data$truncd))
    }else{w=c(w,0)}

    if(prior.weights[3]>0){
      DIHH4h=det(-solve(optH4_NI$hessian))
      DIHH4=ifelse(DIHH4h<0,0,DIHH4h)
      w = c(w, fun.w(DIHH4, llH4N, minll, optH4_NI, data$priormu, data$priorSigma,
                     data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c,
                     data$truncd))
    }else{w=c(w,0)}

    if(prior.weights[4]>0){
      DIHLN4h=det(-solve(optLN4_NI$hessian))
      DIHLN4=ifelse(DIHLN4h<0,0,DIHLN4h)
      w = c(w, fun.w(DIHLN4, llLN4N, minll, optLN4_NI, data$priormu, data$priorSigma,
                     data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c,
                     data$truncd))
    }else{w=c(w,0)}

    if(prior.weights[5]>0){
      DIHG4h=det(-solve(optG4_NI$hessian))
      DIHG4=ifelse(DIHG4h<0,0,DIHG4h)
      w = c(w, fun.w(DIHG4, llG4N, minll, optG4_NI, data$priormu, data$priorSigma,
                     data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c,
                     data$truncd))
    }else{w=c(w,0)}

    if(prior.weights[6]>0){
      DIHQE4h=det(-solve(optQE4_NI$hessian))
      DIHQE4=ifelse(DIHQE4h<0,0,DIHQE4h)
      w = c(w, fun.w(DIHQE4, llQE4N, minll, optQE4_NI, data$priormuQ, data$priorSigmaQ,
                     data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c,
                     data$truncdQ))
    }else{w=c(w,0)}

    if(prior.weights[7]>0){
      DIHP4h=det(-solve(optP4_NI$hessian))
      DIHP4=ifelse(DIHP4h<0,0,DIHP4h)
      w = c(w, fun.w(DIHP4, llP4N, minll, optP4_NI, data$priormu, data$priorSigma,
                     data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c,
                     data$truncd))
    }else{w=c(w,0)}

    if(prior.weights[8]>0){
      DIHL4h=det(-solve(optL4_NI$hessian))
      DIHL4=ifelse(DIHL4h<0,0,DIHL4h)
      w = c(w, fun.w(DIHL4, llL4N, minll, optL4_NI, data$priormu, data$priorSigma,
                     data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c,
                     data$truncd))
    }else{w=c(w,0)}


    # lognormal

    data=data.LN$data
    start=data.LN$start
    startQ=data.LN$startQ

    if(prior.weights[9]>0){
      DIHE4h=det(-solve(optE4_LNI$hessian))
      DIHE4=ifelse(DIHE4h<0,0,DIHE4h)
      w = c(w, fun.w(DIHE4, llE4LN, minll, optE4_LNI, data$priormu, data$priorSigma,
                     data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c,
                     data$truncd))
      # DIHE4h=det(-solve(hessian(func = llE4fLN,x=optE4_LNI$par[c(1,2,9,3,4)],method.args=list(eps=1e-4, d=0.1, zero.tol=sqrt(.Machine$double.eps/7e-7), r=8, v=2, show.details=FALSE))))
      # DIHE4=ifelse(DIHE4h<0,0,DIHE4h)

    }else{w=c(w,0)}

    if(prior.weights[10]>0){
      DIHIE4h=det(-solve(optIE4_LNI$hessian))
      DIHIE4=ifelse(DIHIE4h<0,0,DIHIE4h)
      w = c(w, fun.w(DIHIE4, llIE4LN, minll, optIE4_LNI, data$priormu, data$priorSigma,
                     data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c,
                     data$truncd))
    }else{w=c(w,0)}

    if(prior.weights[11]>0){
      DIHH4h=det(-solve(optH4_LNI$hessian))
      DIHH4=ifelse(DIHH4h<0,0,DIHH4h)
      w = c(w, fun.w(DIHH4, llH4LN, minll, optH4_LNI, data$priormu, data$priorSigma,
                     data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c,
                     data$truncd))
    }else{w=c(w,0)}

    if(prior.weights[12]>0){
      DIHLN4h=det(-solve(optLN4_LNI$hessian))
      DIHLN4=ifelse(DIHLN4h<0,0,DIHLN4h)
      w = c(w, fun.w(DIHLN4, llLN4LN, minll, optLN4_LNI, data$priormu, data$priorSigma,
                     data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c,
                     data$truncd))
    }else{w=c(w,0)}

    if(prior.weights[13]>0){
      DIHG4h=det(-solve(optG4_LNI$hessian))
      DIHG4=ifelse(DIHG4h<0,0,DIHG4h)
      w = c(w, fun.w(DIHG4, llG4LN, minll, optG4_LNI, data$priormu, data$priorSigma,
                     data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c,
                     data$truncd))
    }else{w=c(w,0)}

    if(prior.weights[14]>0){
      DIHQE4h=det(-solve(optQE4_LNI$hessian))
      DIHQE4=ifelse(DIHQE4h<0,0,DIHQE4h)
      w = c(w, fun.w(DIHQE4, llQE4LN, minll, optQE4_LNI, data$priormuQ, data$priorSigmaQ,
                     data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c,
                     data$truncdQ))
    }else{w=c(w,0)}

    if(prior.weights[15]>0){
      DIHP4h=det(-solve(optP4_LNI$hessian))
      DIHP4=ifelse(DIHP4h<0,0,DIHP4h)
      w = c(w, fun.w(DIHP4, llP4LN, minll, optP4_LNI, data$priormu, data$priorSigma,
                     data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c,
                     data$truncd))
    }else{w=c(w,0)}

    if(prior.weights[16]>0){
      DIHL4h=det(-solve(optL4_LNI$hessian))
      DIHL4=ifelse(DIHL4h<0,0,DIHL4h)
      w = c(w, fun.w(DIHL4, llL4LN, minll, optL4_LNI, data$priormu, data$priorSigma,
                     data$priorlb, data$priorub, data$shape.a, data$shape.BMD, data$shape.c,
                     data$truncd))
    }else{w=c(w,0)}

    w <- ifelse(w == 'Inf' | is.na(w), 0, w)
    lpwlp=(prior.weights*w)/sum(prior.weights*w)

  }
  names(lpwlp) = c("E4_N","IE4_N","H4_N","LN4_N","G4_N","QE4_N","P4_N","L4_N","E4_LN","IE4_LN","H4_LN","LN4_LN","G4_LN","QE4_LN","P4_LN","L4_LN")

  # the model average posterior as a mixture
  count=round(lpwlp*ndraws)
  mabmd=(c(# normal
    if(prior.weights[1]>0) sample(as.matrix(fitstanE4_N)[,2],count[1],replace=T),
    if(prior.weights[2]>0) sample(as.matrix(fitstanIE4_N)[,2],count[2],replace=T),
    if(prior.weights[3]>0) sample(as.matrix(fitstanH4_N)[,2],count[3],replace=T),
    if(prior.weights[4]>0) sample(as.matrix(fitstanLN4_N)[,2],count[4],replace=T),
    if(prior.weights[5]>0) sample(as.matrix(fitstanG4_N)[,2],count[5],replace=T),
    if(prior.weights[6]>0) sample(as.matrix(fitstanQE4_N)[,2],count[6],replace=T),
    if(prior.weights[7]>0) sample(as.matrix(fitstanP4_N)[,2],count[7],replace=T),
    if(prior.weights[8]>0) sample(as.matrix(fitstanL4_N)[,2],count[8],replace=T),
    # lognormal
    if(prior.weights[9]>0)  sample(as.matrix(fitstanE4_LN)[,2],count[9],replace=T),
    if(prior.weights[10]>0) sample(as.matrix(fitstanIE4_LN)[,2],count[10],replace=T),
    if(prior.weights[11]>0) sample(as.matrix(fitstanH4_LN)[,2],count[11],replace=T),
    if(prior.weights[12]>0) sample(as.matrix(fitstanLN4_LN)[,2],count[12],replace=T),
    if(prior.weights[13]>0) sample(as.matrix(fitstanG4_LN)[,2],count[13],replace=T),
    if(prior.weights[14]>0) sample(as.matrix(fitstanQE4_LN)[,2],count[14],replace=T),
    if(prior.weights[15]>0) sample(as.matrix(fitstanP4_LN)[,2],count[15],replace=T),
    if(prior.weights[16]>0) sample(as.matrix(fitstanL4_LN)[,2],count[16],replace=T)
  ))

  macilp=(quantile(mabmd,pvec))*data$maxD
  names(macilp)=c("BMDL","BMD","BMDU") # on original scale

  if(TRUE %in% (mabmd > data$maxD) && data$maxD > 1){
    mabmd = ifelse(mabmd > data$maxD, data$maxD, mabmd)
    p.msg = 'The model averaged posterior distribution has been truncated at max(Dose)^2'
    warnings('The model averaged posterior distribution has been truncated at max(Dose)^2')
  }else{
    p.msg = ''
  }

  BMDq_ls = (quantile(mabmd, seq(0,1,0.005)))*data$maxD


  ### Model-averaged response per dose level
  dr.MA.ls <- c()
  for(i in 1:length(data$x)){
    dr.MA.ls[i] = weighted.mean(x = c(DRM_E4_N[i],DRM_IE4_N[i],DRM_H4_N[i],DRM_LN4_N[i],DRM_G4_N[i],DRM_QE4_N[i], DRM_P4_N[i],
                                      DRM_L4_N[i] ,
                                      DRM_E4_LN[i], DRM_IE4_LN[i], DRM_H4_LN[i], DRM_LN4_LN[i], DRM_G4_LN[i],DRM_QE4_LN[i],
                                      DRM_P4_LN[i],DRM_L4_LN[i]),
                                w = lpwlp,
                                na.rm = T)
  }


  ##### Weights & MA if one of the models is divergent --> this model gets weight 0
  if((0 %in% converged) && (1 %in% converged)){

    prior.weights.new = prior.weights
    div.models = which(converged == 0)
    prior.weights.new[div.models] = 0
    p.weights.new = prior.weights.new/sum(prior.weights.new==1)

    lpwlp.conv=(p.weights.new*w)/sum(p.weights.new*w)
    # lpwlp = lpwlp*prior.weights

    # the model average posterior as a mixture
    count=round(lpwlp.conv*ndraws)
    mabmd.conv=(c(# normal
      if(p.weights.new[1]>0) sample(as.matrix(fitstanE4_N)[,2],count[1],replace=T),
      if(p.weights.new[2]>0) sample(as.matrix(fitstanIE4_N)[,2],count[2],replace=T),
      if(p.weights.new[3]>0) sample(as.matrix(fitstanH4_N)[,2],count[3],replace=T),
      if(p.weights.new[4]>0) sample(as.matrix(fitstanLN4_N)[,2],count[4],replace=T),
      if(p.weights.new[5]>0) sample(as.matrix(fitstanG4_N)[,2],count[5],replace=T),
      if(p.weights.new[6]>0) sample(as.matrix(fitstanQE4_N)[,2],count[6],replace=T),
      if(p.weights.new[7]>0) sample(as.matrix(fitstanP4_N)[,2],count[7],replace=T),
      if(p.weights.new[8]>0) sample(as.matrix(fitstanL4_N)[,2],count[8],replace=T),
      # lognormal
      if(p.weights.new[9]>0)  sample(as.matrix(fitstanE4_LN)[,2],count[9],replace=T),
      if(p.weights.new[10]>0) sample(as.matrix(fitstanIE4_LN)[,2],count[10],replace=T),
      if(p.weights.new[11]>0) sample(as.matrix(fitstanH4_LN)[,2],count[11],replace=T),
      if(p.weights.new[12]>0) sample(as.matrix(fitstanLN4_LN)[,2],count[12],replace=T),
      if(p.weights.new[13]>0) sample(as.matrix(fitstanG4_LN)[,2],count[13],replace=T),
      if(p.weights.new[14]>0) sample(as.matrix(fitstanQE4_LN)[,2],count[14],replace=T),
      if(p.weights.new[15]>0) sample(as.matrix(fitstanP4_LN)[,2],count[15],replace=T),
      if(p.weights.new[16]>0) sample(as.matrix(fitstanL4_LN)[,2],count[16],replace=T)
    ))

    macilp.conv=(quantile(mabmd.conv,pvec))*data$maxD
    names(macilp.conv)=c("BMDL","BMD","BMDU")

    if(TRUE %in% (mabmd.conv > data$maxD) && data$maxD > 1){
      mabmd.conv = ifelse(mabmd.conv > data$maxD, data$maxD, mabmd.conv)
      p.msg = 'The model averaged posterior distribution has been truncated at max(Dose)^2'
      warnings('The model averaged posterior distribution has been truncated at max(Dose)^2')
    }else{
      p.msg = ''
    }

    BMDq_ls_conv = (quantile(mabmd.conv, seq(0,1,0.005)))*data$maxD


    ### Model-averaged response per dose level
    dr.MA.ls.conv <- c()
    for(i in 1:length(data$x)){
      dr.MA.ls.conv[i] = weighted.mean(x = c(DRM_E4_N[i],DRM_IE4_N[i],DRM_H4_N[i],DRM_LN4_N[i],DRM_G4_N[i],DRM_QE4_N[i],
                                             DRM_P4_N[i], DRM_L4_N[i] ,
                                             DRM_E4_LN[i], DRM_IE4_LN[i], DRM_H4_LN[i], DRM_LN4_LN[i], DRM_G4_LN[i],DRM_QE4_LN[i],
                                             DRM_P4_LN[i],DRM_L4_LN[i]),
                                       w = lpwlp.conv,
                                       na.rm = T)
    }


  }else{
    lpwlp.conv = NULL; macilp.conv = NULL; BMDq_ls_conv = NULL; dr.MA.ls.conv = NULL; mabmd.conv = NA;
  }

  ## Plot with weights bridge sampling
  BMDL = c(BMDL, macib[1]); BMD = c(BMD, macib[2]); BMDU = c(BMDU, macib[3]) # all on original scale

  names(BMDL) <- c("E4_N","IE4_N","H4_N","LN4_N","G4_N","QE4_N","P4_N","L4_N","E4_LN","IE4_LN","H4_LN","LN4_LN","G4_LN","QE4_LN","P4_LN","L4_LN","MA")
  model = c("E4_N","IE4_N","H4_N","LN4_N","G4_N","QE4_N","P4_N","L4_N","E4_LN","IE4_LN","H4_LN","LN4_LN","G4_LN","QE4_LN","P4_LN","L4_LN","MA")
  model = as.factor(model)

  weight = c(rep(0,16),1)
  names(weight) = c("E4_N","IE4_N","H4_N","LN4_N","G4_N","QE4_N","P4_N","L4_N","E4_LN","IE4_LN","H4_LN","LN4_LN","G4_LN","QE4_LN","P4_LN","L4_LN","MA")
  for(i in names(lpwb)){
    weight[names(weight)==i] = lpwb[names(lpwb)==i]
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

  modelnames <- c("E4_N","IE4_N","H4_N","LN4_N","G4_N","QE4_N","P4_N","L4_N","E4_LN","IE4_LN","H4_LN","LN4_LN","G4_LN","QE4_LN","P4_LN","L4_LN")

  ### Some additional checks


  if(macib[2]/macib[1] > 20){
    warning('BMD/BMDL is larger than 20 for bridge sampling')
  }
  if(macib[3]/macib[1] > 50){
    warning('BMDU/BMDL is larger than 50 for bridge sampling')
  }
  if(macib[2] < (data.N$data$x[2]*data.N$data$maxD/10)){
    warning('BMD is 10 times lower than the lowest non-zero dose for bridge sampling')
  }
  if(macilp[2]/macilp[1] > 20){
    warning('BMD/BMDL is larger than 20 for hybrid Laplace')
  }
  if(macilp[3]/macilp[1] > 50){
    warning('BMDU/BMDL is larger than 50 for hybrid Laplace')
  }
  if(macilp[2] < (data.N$data$x[2]*data.N$data$maxD/10)){
    warning('BMD is 10 times lower than the lowest non-zero dose for hybrid Laplace')
  }

  ### best fitting model vs saturated ANOVA model
  best.fit = modelnames[which(weight[1:16] == max(weight[1:16]))][1]

  bfTest <- modelTestC(best.fit, data.N, data.LN, get(paste0('fitstan', best.fit)), type = 'MCMC',
                       seed, ndraws, nrchains, nriterations, warmup, delta, treedepth)
  warning(bfTest$warn.bf)

  ret_results <- list(E4_N=E4outNI,IE4_N=IE4outNI,H4_N=H4outNI,LN4_N=LN4outNI,
                      G4_N=G4outNI,QE4_N=QE4outNI,P4_N=P4outNI,L4_N=L4outNI,
                      E4_LN=E4outLNI,IE4_LN=IE4outLNI,H4_LN=H4outLNI,LN4_LN=LN4outLNI,
                      G4_LN=G4outLNI,QE4_LN=QE4outLNI,P4_LN=P4outLNI,L4_LN=L4outLNI,
                      covs = covs, corrs = corrs,
                      weights_bridge_sampling=w.bs,
                      weights_laplace=lpwlp,
                      MA_bridge_sampling=macib,
                      MA_laplace=macilp,
                      llN=llN, llLN=llLN,
                      convergence=converged, bs_weights_conv=w.bs.conv,
                      ls_weights_conv=lpwlp.conv,
                      MA_bs_conv=macib.conv,
                      MA_ls_conv=macilp.conv,
                      MA_post_bs = BMDq_bs,
                      MA_post_ls = BMDq_ls,
                      MA_post_bs_conv = BMDq_bs_conv,
                      MA_post_ls_conv = BMDq_ls_conv,
                      MA_dr_bs = dr.MA.bs,
                      MA_dr_ls = dr.MA.ls,
                      MA_dr_bs_conv = dr.MA.bs.conv,
                      MA_dr_ls_conv = dr.MA.ls.conv,
                      parsN = list(parsE4N, parsIE4N, parsH4N, parsLN4N, parsG4N,
                                   parsQE4N, parsP4N, parsL4N),
                      parsLN = list(parsE4LN, parsIE4LN, parsH4LN, parsLN4LN, parsG4LN,
                                    parsQE4LN, parsP4LN, parsL4LN),
                      BMDMixture = (mabmd)*data$maxD,
                      BMDMixture.conv = (mabmd.conv)*data$maxD,
                      BMDMixtureBS = (mabmd1)*data$maxD,
                      BMDMixture.convBS = mabmd.conv1*data$maxD,
                      divergences = divergences,
                      data = data.N$data$data,
                      max.dose = data.N$data$maxD,
                      q = data.N$data$q,
                      # increasing = T,
                      models_included_bridge = modelnames[p.weights > 0],
                      models_included_laplace = modelnames[lpwlp > 0],
                      bf = bfTest$bayesFactor, gof_check = bfTest$warn.bf,
                      # means.SM = bfTest$means.SM, parBestFit = bfTest$par.best,
                      # BIC.bestfit = bfTest$BIC.bestfit, BIC.SM = bfTest$BIC.SM,
                      shift = data.LN$data$shift,
                      w.msg = w.msg, p.msg = p.msg
  )

  attr(ret_results, "class") <- c("BMADR", "BS")

  return(ret_results)

}



#' @rdname sampling_MA
#' @export
samplingQ_MA=function(data.Q,prior.weights = rep(1,8),
                      ndraws = 30000,nrchains=3,
                      nriterations=5000,warmup=1000,
                      delta=0.8,treedepth=10,seed=123,pvec=c(0.025,0.5,0.975)){

  prior.weights.orig = prior.weights

  data=data.Q$data
  start=data.Q$start
  startQ = data.Q$startQ

  BMDL=c(); BMD=c(); BMDU=c()
  converged=c()

  nModels = 8

  ## Obtain model parameters via MCMC sampling

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='quantal')[1]))
  }
  if(prior.weights[1]>0){
    # print(1)
    fitstanE4_Q = fun_samplingQ(stanmodels$mE4_Q, data, start,
                                ndraws,nrchains,
                                nriterations,warmup,
                                delta,treedepth,seed,pvec)

    if(is.null(fitstanE4_Q)){
      prior.weights[1] <- 0
      warning('difficulties fitting the Exponential model using sampling; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      if(data$is_bin == 1) {
        parsE4Q <- parq_extract(fitstanE4_Q, model_name = "E4_Q",
                                pars = c('a', 'b', 'd', 'BMD', 'par1', 'par2', 'par3'))
      } else {
        parsE4Q <- parq_extract(fitstanE4_Q, model_name = "E4_Q",
                                pars = c('a', 'b', 'd', 'BMD', 'rho[1]', 'par1', 'par2', 'par3'),
                                rho = TRUE)
      }
      E4resQ <-  quantile(parsE4Q$BMD, pvec)*data$maxD

      if(data$is_bin == 1){
        E4resQ <- c(E4resQ, apply(parsE4Q[,c('a', 'b', 'd', 'p1', 'p2', 'p3')], 2, median))
      } else {
        E4resQ <- c(E4resQ, apply(parsE4Q[,c('a', 'b', 'd', 'rho', 'p1', 'p2', 'p3')], 2, median))
      }
      names(E4resQ)[1:3] <- c("BMDL","BMD","BMDU")

      if(data$is_bin == 1){
        E4outQ <- outLPQ(parsE4Q, pvec, data$maxD)
      } else {
        E4outQ <- outLPQ(parsE4Q, pvec, data$maxD, rho=TRUE)
      }

      DRM_E4_Q = DRM.E4_Q(E4resQ[names(E4resQ) %in% c('p1', 'p2', 'p3')], data$x, data$q)

      E4covQ = c(cov(as.matrix(fitstanE4_Q)[,c("b","d")], use="complete.obs")["b","d"],
                 cov(as.matrix(fitstanE4_Q)[,c("BMD","d")], use="complete.obs")["BMD","d"])

      E4corrQ = c(cor(as.matrix(fitstanE4_Q)[,c("b","d")], use="complete.obs")["b","d"],
                  cor(as.matrix(fitstanE4_Q)[,c("BMD","d")], use="complete.obs")["BMD","d"])

      BMDL=c(BMDL, E4resQ[1]); BMD=c(BMD, E4resQ[2]); BMDU=c(BMDU, E4resQ[3])
      bridgeE4Q = bridgesampling::bridge_sampler(fitstanE4_Q, silent = T)

      #diagnostics here
      # posterior_diag(model_stan = fitstanE4_Q)
      convergence_stat <- convergence_deci(model_stan = fitstanE4_Q, nrchains = nrchains,
                                           pars = c(paste0('par',1:3), 'a', 'b', 'd', "lp__"))
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)
      # conv_a_r = 0; if(convergence_stat['a', 'Rhat_Convergence'] == "Convergence") conv_a_r = 1

      # divergent transitions
      div_E4_Q <- sum(sapply(rstan::get_sampler_params(fitstanE4_Q, inc_warmup = F),
                             function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))

      # parsE4Q[,"BMD"] <- parsE4Q[,"BMD"]*data$maxD # BMD on original scale
      # Covariance between b-d and between BMD-d
      # compute log marginal likelihood

    }
  }
  if(prior.weights[1] == 0){E4resQ=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA);
  BMDU=c(BMDU,NA); bridgeE4Q=NA; converged=c(converged, NA);
  E4covQ=rep(NA,2); E4corrQ=rep(NA,2); DRM_E4_Q=rep(NA,length(data$x))
  parsE4Q <- NA
  div_E4_Q <- NA
  E4outQ <- t(data.frame(
    # transformed parameters
    par.at = rep(NA,3),
    par.ct = rep(NA,3),
    par.dt = rep(NA,3),
    par.is2t = rep(NA,3),
    par.k = rep(NA,3),
    # parameters on original scale
    par.a = rep(NA,3),
    par.b = rep(NA,3),
    par.c = rep(NA,3),
    par.d = rep(NA,3),
    par.s2 = rep(NA,3),
    # natural parameters
    BMD = rep(NA,3),
    min.resp = rep(NA,3),
    max.resp = rep(NA,3)
  ))
  }

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='quantal')[2]))
  }
  if(prior.weights[2]>0){
    # print(2)
    fitstanIE4_Q = fun_samplingQ(stanmodels$mIE4_Q, data, start,
                                 ndraws,nrchains,
                                 nriterations,warmup,
                                 delta,treedepth,seed,pvec)

    if(is.null(fitstanIE4_Q)){
      prior.weights[2] <- 0
      warning('difficulties fitting the Inverse Exponential model using sampling; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      if(data$is_bin == 1) {
        parsIE4Q <- parq_extract(fitstanIE4_Q, model_name = "IE4_Q",
                                 pars = c('a', 'b', 'd', 'BMD', 'par1', 'par2', 'par3'))
      } else {
        parsIE4Q <- parq_extract(fitstanIE4_Q, model_name = "IE4_Q",
                                 pars = c('a', 'b', 'd', 'BMD', 'rho[1]', 'par1', 'par2', 'par3'),
                                 rho = TRUE)
      }
      IE4resQ <-  quantile(parsIE4Q$BMD, pvec)*data$maxD

      if(data$is_bin == 1){
        IE4resQ <- c(IE4resQ, apply(parsIE4Q[,c('a', 'b', 'd', 'p1', 'p2', 'p3')], 2, median))
      } else {
        IE4resQ <- c(IE4resQ, apply(parsIE4Q[,c('a', 'b', 'd', 'rho', 'p1', 'p2', 'p3')], 2, median))
      }
      names(IE4resQ)[1:3] <- c("BMDL","BMD","BMDU")

      if(data$is_bin == 1){
        IE4outQ <- outLPQ(parsIE4Q, pvec, data$maxD)
      } else {
        IE4outQ <- outLPQ(parsIE4Q, pvec, data$maxD, rho=TRUE)
      }

      DRM_IE4_Q <- DRM.IE4_Q(IE4resQ[names(IE4resQ) %in% c('p1', 'p2', 'p3')],
                             data$x, data$q)

      # Covariance between b-d and between BMD-d
      IE4covQ <- c(cov(parsIE4Q[,c("b","d")], use="complete.obs")["b","d"],
                   cov(parsIE4Q[,c("BMD","d")], use="complete.obs")["BMD","d"])

      IE4corrQ <- c(cor(parsIE4Q[,c("b","d")], use="complete.obs")["b","d"],
                    cor(parsIE4Q[,c("BMD","d")], use="complete.obs")["BMD","d"])
      BMDL=c(BMDL, IE4resQ[1]); BMD=c(BMD, IE4resQ[2]); BMDU=c(BMDU, IE4resQ[3])
      bridgeIE4Q = bridgesampling::bridge_sampler(fitstanIE4_Q, silent = T)

      #diagnostics here
      # posterior_diag(model_stan = fitstanIE4_N)
      convergence_stat <- convergence_deci(model_stan = fitstanIE4_Q, nrchains = nrchains,
                                           pars = c(paste0('par',1:3), 'a', 'b', 'd', "lp__"))
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)

      # divergent transitions
      div_IE4_Q <- sum(sapply(rstan::get_sampler_params(fitstanIE4_Q, inc_warmup = F),
                              function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))

    }
  }
  if(prior.weights[2] == 0){IE4resQ=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeIE4Q=NA;
  converged=c(converged, NA); IE4covQ=rep(NA,2); IE4corrQ=rep(NA,2); DRM_IE4_Q=rep(NA,length(data$x))
  parsIE4Q <- NA
  div_IE4_Q <- NA

  IE4outQ <- t(data.frame(
    # transformed parameters
    par.at = rep(NA,3),
    par.ct = rep(NA,3),
    par.dt = rep(NA,3),
    par.is2t = rep(NA,3),
    par.k = rep(NA,3),
    # parameters on original scale
    par.a = rep(NA,3),
    par.b = rep(NA,3),
    par.c = rep(NA,3),
    par.d = rep(NA,3),
    par.s2 = rep(NA,3),
    # natural parameters
    BMD = rep(NA,3),
    min.resp = rep(NA,3),
    max.resp = rep(NA,3)
  ))
  }

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='quantal')[3]))
  }
  if(prior.weights[3]>0){
    # print(3)
    fitstanH4_Q = fun_samplingQ(stanmodels$mH4_Q, data, start,
                                ndraws,nrchains,
                                nriterations,warmup,
                                delta,treedepth,seed,pvec)

    if(is.null(fitstanH4_Q)){
      prior.weights[3] <- 0
      warning('difficulties fitting the Hill model using sampling; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      if(data$is_bin == 1) {
        parsH4Q <- parq_extract(fitstanH4_Q, model_name = "H4_Q",
                                pars = c('a', 'b', 'd', 'BMD', 'par1', 'par2', 'par3'))
      } else {
        parsH4Q <- parq_extract(fitstanH4_Q, model_name = "H4_Q",
                                pars = c('a', 'b', 'd', 'BMD', 'rho[1]', 'par1', 'par2', 'par3'),
                                rho = TRUE)
      }
      H4resQ <-  quantile(parsH4Q$BMD, pvec)*data$maxD

      if(data$is_bin == 1){
        H4resQ <- c(H4resQ, apply(parsH4Q[,c('a', 'b', 'd', 'p1', 'p2', 'p3')], 2, median))
      } else {
        H4resQ <- c(H4resQ, apply(parsH4Q[,c('a', 'b', 'd', 'rho', 'p1', 'p2', 'p3')], 2, median))
      }
      names(H4resQ)[1:3] <- c("BMDL","BMD","BMDU")

      if(data$is_bin == 1){
        H4outQ <- outLPQ(parsH4Q, pvec, data$maxD)
      } else {
        H4outQ <- outLPQ(parsH4Q, pvec, data$maxD, rho=TRUE)
      }

      DRM_H4_Q = DRM.H4_Q(H4resQ[names(H4resQ) %in% c('p1', 'p2', 'p3')],
                          data$x, data$q)

      # Covariance between b-d and between BMD-d
      H4covQ = c(cov(parsH4Q[,c("b","d")], use="complete.obs")["b","d"],
                 cov(parsH4Q[,c("BMD","d")], use="complete.obs")["BMD","d"])

      H4corrQ = c(cor(parsH4Q[,c("b","d")], use="complete.obs")["b","d"],
                  cor(parsH4Q[,c("BMD","d")], use="complete.obs")["BMD","d"])

      bridgeH4Q = bridge_sampler(fitstanH4_Q, silent = T)
      BMDL = c(BMDL, H4resQ[1]); BMD=c(BMD, H4resQ[2]); BMDU=c(BMDU, H4resQ[3])

      #diagnostics here
      # posterior_diag(model_stan = fitstanH4_N)
      convergence_stat <- convergence_deci(model_stan = fitstanH4_Q, nrchains = nrchains,
                                           pars = c(paste0('par',1:3), 'a', 'b', 'd', "lp__"))
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)

      # divergent transitions
      div_H4_Q <- sum(sapply(rstan::get_sampler_params(fitstanH4_Q, inc_warmup = F),
                             function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))


    }
  }
  if(prior.weights[3] == 0){H4resQ=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeH4Q=NA;
  converged=c(converged, NA); H4covQ=rep(NA,2); H4corrQ=rep(NA,2); DRM_H4_Q=rep(NA,length(data$x))
  parsH4Q <- NA
  div_H4_Q <- NA

  H4outQ <- t(data.frame(
    # transformed parameters
    par.at = rep(NA,3),
    par.ct = rep(NA,3),
    par.dt = rep(NA,3),
    par.is2t = rep(NA,3),
    par.k = rep(NA,3),
    # parameters on original scale
    par.a = rep(NA,3),
    par.b = rep(NA,3),
    par.c = rep(NA,3),
    par.d = rep(NA,3),
    par.s2 = rep(NA,3),
    # natural parameters
    BMD = rep(NA,3),
    min.resp = rep(NA,3),
    max.resp = rep(NA,3)
  ))
  }

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='quantal')[4]))
  }
  if(prior.weights[4]>0){
    # print(4)
    fitstanLN4_Q = fun_samplingQ(stanmodels$mLN4_Q, data, start,
                                 ndraws,nrchains,
                                 nriterations,warmup,
                                 delta,treedepth,seed,pvec)

    if(is.null(fitstanLN4_Q)){
      prior.weights[4] <- 0
      warning('difficulties fitting the Lognormal model using sampling; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      if(data$is_bin == 1) {
        parsLN4Q <- parq_extract(fitstanLN4_Q, model_name = "LN4_Q",
                                 pars = c('a', 'b', 'd', 'BMD', 'par1', 'par2', 'par3'))
      } else {
        parsLN4Q <- parq_extract(fitstanLN4_Q, model_name = "LN4_Q",
                                 pars = c('a', 'b', 'd', 'BMD', 'rho[1]', 'par1', 'par2', 'par3'),
                                 rho = TRUE)
      }
      LN4resQ <-  quantile(parsLN4Q$BMD, pvec)*data$maxD

      if(data$is_bin == 1){
        LN4resQ <- c(LN4resQ, apply(parsLN4Q[,c('a', 'b', 'd', 'p1', 'p2', 'p3')], 2, median))
      } else {
        LN4resQ <- c(LN4resQ, apply(parsLN4Q[,c('a', 'b', 'd', 'rho', 'p1', 'p2', 'p3')], 2, median))
      }
      names(LN4resQ)[1:3] <- c("BMDL","BMD","BMDU")

      if(data$is_bin == 1){
        LN4outQ <- outLPQ(parsLN4Q, pvec, data$maxD)
      } else {
        LN4outQ <- outLPQ(parsLN4Q, pvec, data$maxD, rho=TRUE)
      }

      DRM_LN4_Q = DRM.LN4_Q(LN4resQ[names(LN4resQ) %in% c('p1', 'p2', 'p3')], data$x, data$q)

      BMDL=c(BMDL, LN4resQ[1]); BMD=c(BMD, LN4resQ[2]); BMDU=c(BMDU, LN4resQ[3])
      bridgeLN4Q = bridgesampling::bridge_sampler(fitstanLN4_Q, silent = T)


      # Covariance between b-d and between BMD-d
      LN4covQ = c(cov(parsLN4Q[,c("b","d")], use="complete.obs")["b","d"],
                  cov(parsLN4Q[,c("BMD","d")], use="complete.obs")["BMD","d"])

      LN4corrQ = c(cor(parsLN4Q[,c("b","d")], use="complete.obs")["b","d"],
                   cor(parsLN4Q[,c("BMD","d")], use="complete.obs")["BMD","d"])

      #diagnostics here
      # posterior_diag(model_stan = fitstanLN4_N)
      convergence_stat <- convergence_deci(model_stan = fitstanLN4_Q, nrchains = nrchains,
                                           pars = c(paste0('par',1:3), 'a', 'b', 'd', "lp__"))
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)

      # divergent transitions
      div_LN4_Q <- sum(sapply(rstan::get_sampler_params(fitstanLN4_Q, inc_warmup = F),
                              function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))


    }
  }
  if(prior.weights[4] == 0){LN4resQ=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA);
  bridgeLN4Q=NA; converged=c(converged, NA); LN4covQ=rep(NA,2); LN4corrQ=rep(NA,2);
  DRM_LN4_Q=rep(NA,length(data$x))
  parsLN4Q <- NA
  div_LN4_Q <- NA

  LN4outQ <- t(data.frame(
    # transformed parameters
    par.at = rep(NA,3),
    par.ct = rep(NA,3),
    par.dt = rep(NA,3),
    par.is2t = rep(NA,3),
    par.k = rep(NA,3),
    # parameters on original scale
    par.a = rep(NA,3),
    par.b = rep(NA,3),
    par.c = rep(NA,3),
    par.d = rep(NA,3),
    par.s2 = rep(NA,3),
    # natural parameters
    BMD = rep(NA,3),
    min.resp = rep(NA,3),
    max.resp = rep(NA,3)
  ))
  }

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='quantal')[5]))
  }
  if(prior.weights[5]>0){
    # print(5)

    # data$init_b <- qgamma(data$q, rate=1.0, shape=E4resQ[6])/(E4resQ[2]/data$maxD)

    fitstanG4_Q <- fun_samplingQ(stanmodels$mG4_Q, data, start,
                                 ndraws,nrchains,
                                 nriterations,warmup,
                                 delta,treedepth,seed,pvec)
    if(is.null(fitstanG4_Q)){
      prior.weights[5] <- 0
      warning('difficulties fitting the Gamma model using sampling; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      if(data$is_bin == 1) {
        parsG4Q <- parq_extract(fitstanG4_Q, model_name = "G4_Q",
                                pars = c('a', 'b', 'd', 'BMD', 'par1', 'par2', 'par3'))
      } else {
        parsG4Q <- parq_extract(fitstanG4_Q, model_name = "G4_Q",
                                pars = c('a', 'b', 'd', 'BMD', 'rho[1]', 'par1', 'par2', 'par3'),
                                rho = TRUE)
      }
      G4resQ <-  quantile(parsG4Q$BMD, pvec)*data$maxD

      if(data$is_bin == 1){
        G4resQ <- c(G4resQ, apply(parsG4Q[,c('a', 'b', 'd', 'p1', 'p2', 'p3')], 2, median))
      } else {
        G4resQ <- c(G4resQ, apply(parsG4Q[,c('a', 'b', 'd', 'rho', 'p1', 'p2', 'p3')], 2, median))
      }
      names(G4resQ)[1:3] <- c("BMDL","BMD","BMDU")

      if(data$is_bin == 1){
        G4outQ <- outLPQ(parsG4Q, pvec, data$maxD)
      } else {
        G4outQ <- outLPQ(parsG4Q, pvec, data$maxD, rho=TRUE)
      }

      DRM_G4_Q = DRM.G4_Q(G4resQ[names(G4resQ) %in% c('p1', 'p2', 'p3')], data$x, data$q)

      # Covariance between b-d and between BMD-d
      G4covQ = c(cov(parsG4Q[,c("b","d")], use="complete.obs")["b","d"],
                 cov(parsG4Q[,c("BMD","d")], use="complete.obs")["BMD","d"])

      G4corrQ = c(cor(parsG4Q[,c("b","d")], use="complete.obs")["b","d"],
                  cor(parsG4Q[,c("BMD","d")], use="complete.obs")["BMD","d"])

      BMDL=c(BMDL, G4resQ[1]); BMD=c(BMD, G4resQ[2]); BMDU=c(BMDU, G4resQ[3])
      bridgeG4Q <- bridgesampling::bridge_sampler(fitstanG4_Q, silent = TRUE)

      #diagnostics here
      # posterior_diag(model_stan = fitstanG4_N)
      convergence_stat <- convergence_deci(model_stan = fitstanG4_Q, nrchains = nrchains,
                                           pars = c(paste0('par',1:3), 'a', 'b', 'd', "lp__"))
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)

      # divergent transitions
      div_G4_Q <- sum(sapply(rstan::get_sampler_params(fitstanG4_Q, inc_warmup = F),
                             function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))


    }
  }
  if(prior.weights[5] == 0){G4resNI=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeG4Q=NA;
  converged=c(converged, NA); G4covQ=rep(NA,2); G4corrQ=rep(NA,2); DRM_G4_Q=rep(NA,length(data$x))
  parsG4Q <- NA
  div_G4_Q <- NA

  G4outQ <- t(data.frame(
    # transformed parameters
    par.at = rep(NA,3),
    par.ct = rep(NA,3),
    par.dt = rep(NA,3),
    par.is2t = rep(NA,3),
    par.k = rep(NA,3),
    # parameters on original scale
    par.a = rep(NA,3),
    par.b = rep(NA,3),
    par.c = rep(NA,3),
    par.d = rep(NA,3),
    par.s2 = rep(NA,3),
    # natural parameters
    BMD = rep(NA,3),
    min.resp = rep(NA,3),
    max.resp = rep(NA,3)
  ))
  }

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='quantal')[6]))
  }
  if(prior.weights[6]>0){
    # print(6)
    fitstanQE4_Q = fun_samplingQ(stanmodels$mQE4_Q, data, startQ,
                                 ndraws,nrchains,
                                 nriterations,warmup,
                                 delta,treedepth,seed,pvec)

    if(is.null(fitstanQE4_Q)){
      prior.weights[6] <- 0
      warning('difficulties fitting the Quadratic Exponential model using sampling; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      if(data$is_bin == 1) {
        parsQE4Q <- parq_extract(fitstanQE4_Q, model_name = "QE4_Q",
                                 pars = c('a', 'b', 'd', 'BMD', 'par1', 'par2', 'par3'))
      } else {
        parsQE4Q <- parq_extract(fitstanQE4_Q, model_name = "QE4_Q",
                                 pars = c('a', 'b', 'd', 'BMD', 'rho[1]', 'par1', 'par2', 'par3'))
      }
      QE4resQ <-  quantile(parsQE4Q$BMD, pvec)*data$maxD

      if(data$is_bin == 1){
        QE4resQ <- c(QE4resQ, apply(parsQE4Q[,c('a', 'b', 'd', 'p1', 'p2', 'p3')], 2, median))
      } else {
        QE4resQ <- c(QE4resQ, apply(parsQE4Q[,c('a', 'b', 'd', 'rho', 'p1', 'p2', 'p3')], 2, median))
      }
      names(QE4resQ)[1:3] <- c("BMDL","BMD","BMDU")

      if(data$is_bin == 1){
        QE4outQ <- outLPQ(parsQE4Q, pvec, data$maxD)
      } else {
        QE4outQ <- outLPQ(parsQE4Q, pvec, data$maxD, rho=TRUE)
      }

      DRM_QE4_Q = DRM.QE4_Q(QE4resQ[names(QE4resQ) %in% c('p1', 'p2', 'p3')],
                            data$x, data$q)

      # Covariance between b-d and between BMD-d
      QE4covQ = c(cov(parsQE4Q[,c("b","d")], use="complete.obs")["b","d"],
                  cov(parsQE4Q[,c("BMD","d")], use="complete.obs")["BMD","d"])

      QE4corrQ = c(cor(parsQE4Q[,c("b","d")], use="complete.obs")["b","d"],
                   cor(parsQE4Q[,c("BMD","d")], use="complete.obs")["BMD","d"])

      BMDL=c(BMDL, QE4resQ[1]); BMD=c(BMD, QE4resQ[2]); BMDU=c(BMDU, QE4resQ[3])
      bridgeQE4Q = bridgesampling::bridge_sampler(fitstanQE4_Q, silent = T)

      #diagnostics here
      # posterior_diag(model_stan = fitstanQE4_N)
      convergence_stat <- convergence_deci(model_stan = fitstanQE4_Q, nrchains = nrchains,
                                           pars = c(paste0('par',1:3), 'a', 'b', 'd', "lp__"))
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)

      # divergent transitions
      div_QE4_Q <- sum(sapply(rstan::get_sampler_params(fitstanQE4_Q, inc_warmup = F),
                              function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))

    }
  }
  if(prior.weights[6] == 0){QE4resQ=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeQE4Q=NA;
  converged=c(converged, NA); QE4covQ=rep(NA,2); QE4corrQ=rep(NA,2); DRM_QE4_Q=rep(NA,length(data$x))
  parsQE4Q <- NA
  div_QE4_Q <- NA

  QE4outQ <- t(data.frame(
    # transformed parameters
    par.at = rep(NA,3),
    par.ct = rep(NA,3),
    par.dt = rep(NA,3),
    par.is2t = rep(NA,3),
    par.k = rep(NA,3),
    # parameters on original scale
    par.a = rep(NA,3),
    par.b = rep(NA,3),
    par.c = rep(NA,3),
    par.d = rep(NA,3),
    par.s2 = rep(NA,3),
    # natural parameters
    BMD = rep(NA,3),
    min.resp = rep(NA,3),
    max.resp = rep(NA,3)
  ))
  }

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='quantal')[7]))
  }
  if(prior.weights[7]>0){
    # print(7)
    fitstanP4_Q = fun_samplingQ(stanmodels$mP4_Q, data, start,
                                ndraws,nrchains,
                                nriterations,warmup,
                                delta,treedepth,seed,pvec)

    if(is.null(fitstanP4_Q)){
      prior.weights[7] <- 0
      warning('difficulties fitting the Probit model using sampling; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      if(data$is_bin == 1) {
        parsP4Q <- parq_extract(fitstanP4_Q, model_name = "P4_Q",
                                pars = c('a', 'b', 'd', 'BMD', 'par1', 'par2', 'par3'))
      } else {
        parsP4Q <- parq_extract(fitstanP4_Q, model_name = "P4_Q",
                                pars = c('a', 'b', 'd', 'BMD', 'rho[1]', 'par1', 'par2', 'par3'),
                                rho = TRUE)
      }
      P4resQ <-  quantile(parsP4Q$BMD, pvec)*data$maxD

      if(data$is_bin == 1){
        P4resQ <- c(P4resQ, apply(parsP4Q[,c('a', 'b', 'd', 'p1', 'p2', 'p3')], 2, median))
      } else {
        P4resQ <- c(P4resQ, apply(parsP4Q[,c('a', 'b', 'd', 'rho', 'p1', 'p2', 'p3')], 2, median))
      }
      names(P4resQ)[1:3] <- c("BMDL","BMD","BMDU")

      if(data$is_bin == 1){
        P4outQ <- outLPQ(parsP4Q, pvec, data$maxD)
      } else {
        P4outQ <- outLPQ(parsP4Q, pvec, data$maxD, rho=TRUE)
      }

      DRM_P4_Q <- DRM.P4_Q(E4resQ[names(P4resQ) %in% c('p1', 'p2', 'p3')], data$x, data$q)

      # Covariance between b-d and between BMD-d
      P4covQ = c(cov(parsP4Q[,c("b","d")], use="complete.obs")["b","d"],
                 cov(parsP4Q[,c("BMD","d")], use="complete.obs")["BMD","d"])

      P4corrQ = c(cor(parsP4Q[,c("b","d")], use="complete.obs")["b","d"],
                  cor(parsP4Q[,c("BMD","d")], use="complete.obs")["BMD","d"])

      BMDL=c(BMDL, P4resQ[1]); BMD=c(BMD, P4resQ[2]); BMDU=c(BMDU, P4resQ[3])
      bridgeP4Q = bridgesampling::bridge_sampler(fitstanP4_Q, silent = T)
      #diagnostics here
      # posterior_diag(model_stan = fitstanP4_N)
      convergence_stat <- convergence_deci(model_stan = fitstanP4_Q, nrchains = nrchains,
                                           pars = c(paste0('par',1:3), 'a', 'b', 'd', "lp__"))
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)

      # divergent transitions
      div_P4_Q <- sum(sapply(rstan::get_sampler_params(fitstanP4_Q, inc_warmup = F),
                             function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))

    }
  }
  if(prior.weights[7] == 0){P4resNI=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA);
  bridgeP4Q=NA; converged=c(converged, NA);
  P4covQ=rep(NA,2); P4corrQ=rep(NA,2);
  DRM_P4_Q=rep(NA,length(data$x))
  parsP4Q <- NA
  div_P4_Q <- NA

  P4outQ <- t(data.frame(
    # transformed parameters
    par.at = rep(NA,3),
    par.ct = rep(NA,3),
    par.dt = rep(NA,3),
    par.is2t = rep(NA,3),
    par.k = rep(NA,3),
    # parameters on original scale
    par.a = rep(NA,3),
    par.b = rep(NA,3),
    par.c = rep(NA,3),
    par.d = rep(NA,3),
    par.s2 = rep(NA,3),
    # natural parameters
    BMD = rep(NA,3),
    min.resp = rep(NA,3),
    max.resp = rep(NA,3)
  ))
  }

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='quantal')[8]))
  }
  if(prior.weights[8]>0){
    # print(8)
    fitstanL4_Q = fun_samplingQ(stanmodels$mL4_Q, data, start,
                                ndraws,nrchains,
                                nriterations,warmup,
                                delta,treedepth,seed,pvec)

    if(is.null(fitstanL4_Q)){
      prior.weights[8] <- 0
      warning('difficulties fitting the Logit model using sampling; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      if(data$is_bin == 1) {
        parsL4Q <- parq_extract(fitstanL4_Q, model_name = "L4_Q",
                                pars = c('a', 'b', 'd', 'BMD', 'par1', 'par2', 'par3'))
      } else {
        parsL4Q <- parq_extract(fitstanL4_Q, model_name = "L4_Q",
                                pars = c('a', 'b', 'd', 'BMD', 'rho[1]', 'par1', 'par2', 'par3'),
                                rho = TRUE)
      }
      L4resQ <-  quantile(parsL4Q$BMD, pvec)*data$maxD

      if(data$is_bin == 1){
        L4resQ <- c(L4resQ, apply(parsL4Q[,c('a', 'b', 'd', 'p1', 'p2', 'p3')], 2, median))
      } else {
        L4resQ <- c(L4resQ, apply(parsL4Q[,c('a', 'b', 'd', 'rho', 'p1', 'p2', 'p3')], 2, median))
      }
      names(L4resQ)[1:3] <- c("BMDL","BMD","BMDU")

      if(data$is_bin == 1){
        L4outQ <- outLPQ(parsL4Q, pvec, data$maxD)
      } else {
        L4outQ <- outLPQ(parsL4Q, pvec, data$maxD, rho=TRUE)
      }

      DRM_L4_Q <- DRM.L4_Q(E4resQ[names(P4resQ) %in% c('p1', 'p2', 'p3')],
                           data$x, data$q)

      # Covariance between b-d and between BMD-d
      L4covQ = c(cov(parsL4Q[,c("b","d")], use="complete.obs")["b","d"],
                 cov(parsL4Q[,c("BMD","d")], use="complete.obs")["BMD","d"])

      L4corrQ = c(cor(parsL4Q[,c("b","d")], use="complete.obs")["b","d"],
                  cor(parsL4Q[,c("BMD","d")], use="complete.obs")["BMD","d"])

      BMDL=c(BMDL, L4resQ[1]); BMD=c(BMD, L4resQ[2]); BMDU=c(BMDU, L4resQ[3])
      bridgeL4Q = bridgesampling::bridge_sampler(fitstanL4_Q, silent = T)
      # save(fitstanL4_N, file = "fitL4N.RData")

      #diagnostics here
      # posterior_diag(model_stan = fitstanL4_N)
      convergence_stat <- convergence_deci(model_stan = fitstanL4_Q, nrchains = nrchains,
                                           pars = c(paste0('par',1:3), 'a', 'b', 'd', "lp__"))
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)

      # divergent transitions
      div_L4_Q <- sum(sapply(rstan::get_sampler_params(fitstanL4_Q, inc_warmup = F),
                             function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))

    }
  }
  if(prior.weights[8] == 0){L4resQ=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA);
  bridgeL4Q=NA; converged=c(converged, NA); L4covQ=rep(NA,2); L4corrQ=rep(NA,2);
  DRM_L4_Q=rep(NA,length(data$x))
  parsL4Q <- NA
  div_L4_Q <- NA

  L4outQ <- t(data.frame(
    # transformed parameters
    par.at = rep(NA,3),
    par.ct = rep(NA,3),
    par.dt = rep(NA,3),
    par.is2t = rep(NA,3),
    par.k = rep(NA,3),
    # parameters on original scale
    par.a = rep(NA,3),
    par.b = rep(NA,3),
    par.c = rep(NA,3),
    par.d = rep(NA,3),
    par.s2 = rep(NA,3),
    # natural parameters
    BMD = rep(NA,3),
    min.resp = rep(NA,3),
    max.resp = rep(NA,3)
  ))
  }

  divergences <- c(div_E4_Q, div_IE4_Q, div_H4_Q, div_LN4_Q, div_G4_Q, div_QE4_Q, div_P4_Q, div_L4_Q)


  #-----------------------------------

  ### Weights based on bridge sampling are obtained by computing posterior model probabilities from the marginal likelihoods
  p.weights = prior.weights/sum(prior.weights==1)

  # lpwb = bridgesampling::post_prob(bridgeE4N, bridgeIE4N, bridgeH4N, bridgeLN4N, bridgeG4N, bridgeQE4N, bridgeP4N, bridgeL4N,
  #                  bridgeE4LN, bridgeIE4LN, bridgeH4LN, bridgeLN4LN, bridgeG4LN, bridgeQE4LN, bridgeP4LN, bridgeL4LN,
  #                  prior_prob = p.weights[p.weights>0],
  #                  model_names = c("E4_N","IE4_N","H4_N","LN4_N","G4_N","QE4_N","P4_N","L4_N",
  #                                  "E4_LN","IE4_LN","H4_LN","LN4_LN","G4_LN","QE4_LN","P4_LN","L4_LN")[p.weights>0])


  bridge.mods = c(if(p.weights[1]>0) bridgeE4Q$logml,
                  if(p.weights[2]>0) bridgeIE4Q$logml,
                  if(p.weights[3]>0) bridgeH4Q$logml,
                  if(p.weights[4]>0) bridgeLN4Q$logml,
                  if(p.weights[5]>0) bridgeG4Q$logml,
                  if(p.weights[6]>0) bridgeQE4Q$logml,
                  if(p.weights[7]>0) bridgeP4Q$logml,
                  if(p.weights[8]>0) bridgeL4Q$logml)

  model_names = c("E4_Q","IE4_Q","H4_Q","LN4_Q","G4_Q","QE4_Q","P4_Q","L4_Q")[p.weights>0]

  lpwb = bridgesampling::post_prob(bridge.mods,
                                   prior_prob = p.weights[p.weights>0],
                                   model_names = model_names)

  # # the model average posterior as a mixture
  count=round(lpwb*ndraws)
  names(count) = names(lpwb)
  mabmd1=(c( # normal
    if("E4_Q" %in% names(count)) sample(as.matrix(fitstanE4_Q)[,2],count[names(count)=="E4_Q"],replace=T),
    if("IE4_Q" %in% names(count)) sample(as.matrix(fitstanIE4_Q)[,2],count[names(count)=="IE4_Q"],replace=T),
    if("H4_Q" %in% names(count)) sample(as.matrix(fitstanH4_Q)[,2],count[names(count)=="H4_Q"],replace=T),
    if("LN4_Q" %in% names(count)) sample(as.matrix(fitstanLN4_Q)[,2],count[names(count)=="LN4_Q"],replace=T),
    if("G4_Q" %in% names(count)) sample(as.matrix(fitstanG4_Q)[,2],count[names(count)=="G4_Q"],replace=T),
    if("QE4_Q" %in% names(count)) sample(as.matrix(fitstanQE4_Q)[,2],count[names(count)=="QE4_Q"],replace=T),
    if("P4_Q" %in% names(count)) sample(as.matrix(fitstanP4_Q)[,2],count[names(count)=="P4_Q"],replace=T),
    if("L4_Q" %in% names(count)) sample(as.matrix(fitstanL4_Q)[,2],count[names(count)=="L4_Q"],replace=T)
  ))

  macib=quantile(mabmd1,pvec)*data$maxD
  names(macib)=c("BMDL","BMD","BMDU") # original scale

  if(TRUE %in% (mabmd1 > data$maxD) && data$maxD > 1){
    mabmd1 = ifelse(mabmd1 > data$maxD, data$maxD, mabmd1)
    p.msg = 'The model averaged posterior distribution has been truncated at max(Dose)^2'
    warnings('The model averaged posterior distribution has been truncated at max(Dose)^2')
  }else{
    p.msg = ''
  }

  BMDq_bs = quantile(mabmd1, seq(0,1,0.005))*data$maxD

  mods = c("E4_Q","IE4_Q","H4_Q","LN4_Q","G4_Q","QE4_Q","P4_Q","L4_Q")
  w.bs = rep(0, length(mods))
  names(w.bs) = mods
  for(i in mods[mods%in%names(lpwb)]){
    w.bs[i] = lpwb[i]
  }

  ### Model-averaged response per dose level
  dr.MA.bs <- c()
  for(i in 1:length(data$x)){
    xx <- c(DRM_E4_Q[i],DRM_IE4_Q[i],DRM_H4_Q[i],DRM_LN4_Q[i],DRM_G4_Q[i],DRM_QE4_Q[i],
            DRM_P4_Q[i], DRM_L4_Q[i])
    dr.MA.bs[i] = weighted.mean(x = xx[!is.na(xx)],
                                w = w.bs[!is.na(xx)],
                                na.rm = T)

  }

  ##### Weights & MA if one of the models is divergent --> this model gets weight 0
  if((0 %in% converged) && (1 %in% converged)){
    prior.weights.new = prior.weights
    div.models = which(converged == 0)
    prior.weights.new[div.models] = 0
    p.weights.new = prior.weights.new/sum(prior.weights.new==1)

    bridge.mods = c(if(p.weights.new[1]>0) bridgeE4Q$logml,
                    if(p.weights.new[2]>0) bridgeIE4Q$logml,
                    if(p.weights.new[3]>0) bridgeH4Q$logml,
                    if(p.weights.new[4]>0) bridgeLN4Q$logml,
                    if(p.weights.new[5]>0) bridgeG4Q$logml,
                    if(p.weights.new[6]>0) bridgeQE4Q$logml,
                    if(p.weights.new[7]>0) bridgeP4Q$logml,
                    if(p.weights.new[8]>0) bridgeL4Q$logml)

    model_names = c("E4_Q","IE4_Q","H4_Q","LN4_Q","G4_Q","QE4_Q","P4_Q","L4_Q")[p.weights.new>0]

    lpwb.conv = bridgesampling::post_prob(bridge.mods,
                                          prior_prob = p.weights.new[p.weights.new>0],
                                          model_names = model_names)

    # # the model average posterior as a mixture
    count=round(lpwb.conv*ndraws)
    names(count) = names(lpwb.conv)
    mabmd.conv1=(c( # normal
      if("E4_Q" %in% names(count)) sample(as.matrix(fitstanE4_Q)[,2],count[names(count)=="E4_Q"],replace=T),
      if("IE4_Q" %in% names(count)) sample(as.matrix(fitstanIE4_Q)[,2],count[names(count)=="IE4_Q"],replace=T),
      if("H4_Q" %in% names(count)) sample(as.matrix(fitstanH4_Q)[,2],count[names(count)=="H4_Q"],replace=T),
      if("LN4_Q" %in% names(count)) sample(as.matrix(fitstanLN4_Q)[,2],count[names(count)=="LN4_Q"],replace=T),
      if("G4_Q" %in% names(count)) sample(as.matrix(fitstanG4_Q)[,2],count[names(count)=="G4_Q"],replace=T),
      if("QE4_Q" %in% names(count)) sample(as.matrix(fitstanQE4_Q)[,2],count[names(count)=="QE4_Q"],replace=T),
      if("P4_Q" %in% names(count)) sample(as.matrix(fitstanP4_Q)[,2],count[names(count)=="P4_Q"],replace=T),
      if("L4_Q" %in% names(count)) sample(as.matrix(fitstanL4_Q)[,2],count[names(count)=="L4_Q"],replace=T)
    ))

    macib.conv <- quantile(mabmd.conv1,pvec)*data$maxD
    names(macib.conv) <- c("BMDL","BMD","BMDU")

    if(TRUE %in% (mabmd.conv1 > data$maxD) && data$maxD > 1){
      mabmd.conv1 = ifelse(mabmd.conv1 > data$maxD, data$maxD, mabmd.conv1)
      p.msg = 'The model averaged posterior distribution has been truncated at max(Dose)^2'
      warnings('The model averaged posterior distribution has been truncated at max(Dose)^2')
    }else{
      p.msg = ''
    }

    BMDq_bs_conv <- quantile(mabmd.conv1, seq(0,1,0.005))*data$maxD


    mods = c("E4_Q","IE4_Q","H4_Q","LN4_Q","G4_Q","QE4_Q","P4_Q","L4_Q")
    w.bs.conv = rep(0, length(mods))
    names(w.bs.conv) = mods
    for(i in mods[mods%in%names(lpwb.conv)]){
      w.bs.conv[i] = lpwb.conv[i]
    }

    ### Model-averaged response per dose level
    dr.MA.bs.conv <- c()
    for(i in 1:length(data$x)){
      xx <- c(DRM_E4_Q[i],DRM_IE4_Q[i],DRM_H4_Q[i],DRM_LN4_Q[i],DRM_G4_Q[i],
              DRM_QE4_Q[i], DRM_P4_Q[i], DRM_L4_Q[i])
      dr.MA.bs.conv[i] = weighted.mean(x = xx[!is.na(xx)],
                                       w = w.bs.conv[!is.na(xx)],
                                       na.rm = T)
    }

  }else{
    w.bs.conv = NULL; macib.conv = NULL; BMDq_bs_conv = NULL; dr.MA.bs.conv = NULL; mabmd.conv1 = NA
  }

  #----------------------------------
  ### weights based on laplace approximation

  prior.weights = prior.weights.orig

  llQ = c() # likelihoods

  # getting the posterior modes and the hessian and the model specific posterior distributions
  if(prior.weights[1]>0){

      optE4_Q <- fun_optimQ(stanmodels$mE4_Q, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optE4_Q[[3]]),TRUE,(optE4_Q[[3]]!=0)) | length(optE4_Q)!=9)){
      prior.weights[1] <- 0
      warning('difficulties fitting the Exponential model using Laplace; prior weight was set to 0 and the model is not included in model averaging')
    }else{

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
                               rho = optE4_Q$par[stringr::str_detect(names(optE4_Q$par),'rho')])
      )

      llQ = c(llQ, llE4Q)
    }
  }
  if(prior.weights[1] == 0){llQ = c(llQ,NA)}
  #
  if(prior.weights[2]>0){
    optIE4_Q = fun_optimQ(stanmodels$mIE4_Q, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optIE4_Q[[3]]),TRUE,(optIE4_Q[[3]]!=0)) | length(optIE4_Q)!=9)){
      prior.weights[2] <- 0
      warning('difficulties fitting the Inverse Exponential model using Laplace; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      llIE4Q <- ifelse(data$is_bin == 1,
                       llfIE4_Q(optIE4_Q$par[1:3],
                                nvec=data$n,
                                dvec=data$x,
                                yvec=data$y,
                                qval=data$q),
                       llfIE42_Q(optIE4_Q$par[1:3],
                                 nvec=data$n,
                                 dvec=data$x,
                                 yvec=data$y,
                                 qval=data$q,
                                 rho = optIE4_Q$par[stringr::str_detect(names(optIE4_Q$par),'rho')])
      )

      llQ = c(llQ, llIE4Q)
    }
  }
  if(prior.weights[2] == 0){llQ = c(llQ,NA)}
  #
  if(prior.weights[3]>0){
    optH4_Q = fun_optimQ(stanmodels$mH4_Q, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optH4_Q[[3]]),TRUE,(optH4_Q[[3]]!=0)) | length(optH4_Q)!=9)){
      prior.weights[3] <- 0
      warning('difficulties fitting the Hill model using Laplace; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      llH4Q <- ifelse(data$is_bin == 1,
                      llfH4_Q(optH4_Q$par[1:3],
                              nvec=data$n,
                              dvec=data$x,
                              yvec=data$y,
                              qval=data$q),
                      llfH42_Q(optH4_Q$par[1:3],
                               nvec=data$n,
                               dvec=data$x,
                               yvec=data$y,
                               qval=data$q,
                               rho = optH4_Q$par[stringr::str_detect(names(optH4_Q$par),'rho')])
      )

      llQ = c(llQ, llH4Q)
    }
  }
  if(prior.weights[3] == 0){llQ = c(llQ,NA)}
  #
  if(prior.weights[4]>0){
    optLN4_Q = fun_optimQ(stanmodels$mLN4_Q, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optLN4_Q[[3]]),TRUE,(optLN4_Q[[3]]!=0)) | length(optLN4_Q)!=9)){
      prior.weights[4] <- 0
      warning('difficulties fitting the Lognormal model using Laplace; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      llLN4Q <- ifelse(data$is_bin == 1,
                       llfLN4_Q(optLN4_Q$par[1:3],
                                nvec=data$n,
                                dvec=data$x,
                                yvec=data$y,
                                qval=data$q),
                       llfH42_Q(optLN4_Q$par[1:3],
                                nvec=data$n,
                                dvec=data$x,
                                yvec=data$y,
                                qval=data$q,
                                rho = optLN4_Q$par[stringr::str_detect(names(optLN4_Q$par),'rho')])
      )

      llQ = c(llQ, llLN4Q)
    }
  }
  if(prior.weights[4] == 0){llQ = c(llQ,NA)}
  #
  if(prior.weights[5]>0){

    optG4_Q = fun_optimQ(stanmodels$mG4_Q, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optG4_Q[[3]]),TRUE,(optG4_Q[[3]]!=0)) | length(optG4_Q)!=9)){
      prior.weights[5] <- 0
      warning('difficulties fitting the Gamma model using Laplace; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      llG4Q <- ifelse(data$is_bin == 1,
                      llfG4_Q(optG4_Q$par[1:3],
                              nvec=data$n,
                              dvec=data$x,
                              yvec=data$y,
                              qval=data$q),
                      llfG42_Q(optG4_Q$par[1:3],
                               nvec=data$n,
                               dvec=data$x,
                               yvec=data$y,
                               qval=data$q,
                               rho = optG4_Q$par[stringr::str_detect(names(optG4_Q$par),'rho')])
      )

      llQ = c(llQ, llG4Q)
    }
  }
  if(prior.weights[5] == 0){llQ = c(llQ,NA)}
  #
  if(prior.weights[6]>0){
    optQE4_Q = fun_optimQ(stanmodels$mQE4_Q, data, startQ, ndraws, 123, pvec)

    if((ifelse(is.na(optQE4_Q[[3]]),TRUE,(optQE4_Q[[3]]!=0)) | length(optQE4_Q)!=9)){
      prior.weights[6] <- 0
      warning('difficulties fitting the Quadratic Exponential model using Laplace; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      llQE4Q <- ifelse(data$is_bin == 1,
                       llfQE4_Q(optQE4_Q$par[1:3],
                                nvec=data$n,
                                dvec=data$x,
                                yvec=data$y,
                                qval=data$q),
                       llfQE42_Q(optQE4_Q$par[1:3],
                                 nvec=data$n,
                                 dvec=data$x,
                                 yvec=data$y,
                                 qval=data$q,
                                 rho = optQE4_Q$par[stringr::str_detect(names(optQE4_Q$par),'rho')])
      )

      llQ = c(llQ, llQE4Q)
    }
  }
  if(prior.weights[6] == 0){llQ = c(llQ,NA)}
  #
  if(prior.weights[7]>0){
    optP4_Q = fun_optimQ(stanmodels$mP4_Q, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optP4_Q[[3]]),TRUE,(optP4_Q[[3]]!=0)) | length(optP4_Q)!=9)){
      prior.weights[7] <- 0
      warning('difficulties fitting the Probit model using Laplace; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      llP4Q <- ifelse(data$is_bin == 1,
                      llfP4_Q(optP4_Q$par[1:3],
                              nvec=data$n,
                              dvec=data$x,
                              yvec=data$y,
                              qval=data$q),
                      llfP42_Q(optP4_Q$par[1:3],
                               nvec=data$n,
                               dvec=data$x,
                               yvec=data$y,
                               qval=data$q,
                               rho = optP4_Q$par[stringr::str_detect(names(optP4_Q$par),'rho')])
      )

      llQ = c(llQ, llP4Q)
    }
  }
  if(prior.weights[7] == 0){llQ = c(llQ,NA)}
  #
  if(prior.weights[8]>0){
    optL4_Q = fun_optimQ(stanmodels$mL4_Q, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optL4_Q[[3]]),TRUE,(optL4_Q[[3]]!=0)) | length(optL4_Q)!=9)){
      prior.weights[8] <- 0
      warning('difficulties fitting the Logit model using Laplace; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      llL4Q <- ifelse(data$is_bin == 1,
                      llfL4_Q(optL4_Q$par[1:3],
                              nvec=data$n,
                              dvec=data$x,
                              yvec=data$y,
                              qval=data$q),
                      llfL42_Q(optL4_Q$par[1:3],
                               nvec=data$n,
                               dvec=data$x,
                               yvec=data$y,
                               qval=data$q,
                               rho = optL4_Q$par[stringr::str_detect(names(optL4_Q$par),'rho')])
      )

      llQ = c(llQ, llL4Q)
    }
  }
  if(prior.weights[8] == 0){llQ = c(llQ,NA)}

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
    lpwlp <- rep(0, 8)
    lpwlp[which(lls == max.ll)] <- 1
    w.msg <- 'Laplace weights could not be computed and one model gets all the weight; using another prior for parameter d might help'
    warning('Laplace weights could not be computed and one model gets all the weight; using another prior for parameter d might help')
  }else{

    if(FALSE %in% ((max.ll-lls[!is.na(lls)]) < 709)){
      w.msg <- 'not all models were used in computation of Laplace weights, some models set to 0; using another prior for parameter d might help'
      warning('not all models were used in computation of Laplace weights, some models set to 0; using another prior for parameter d might help')
    }

    minll <- min(lls[which((max.ll-lls[!is.na(lls)]) < 709 & prior.weights>0)], na.rm = T)

  if(prior.weights[1]>0){
    DIHE4h=det(-solve(optE4_Q$hessian))
    DIHE4=ifelse(DIHE4h<0,0,DIHE4h)
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
    DIHIE4h=det(-solve(optIE4_Q$hessian))
    DIHIE4=ifelse(DIHIE4h<0,0,DIHIE4h)

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
    DIHH4h=det(-solve(optH4_Q$hessian))
    DIHH4=ifelse(DIHH4h<0,0,DIHH4h)

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
    DIHLN4h=det(-solve(optLN4_Q$hessian))
    DIHLN4=ifelse(DIHLN4h<0,0,DIHLN4h)

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
    DIHG4h=det(-solve(optG4_Q$hessian))
    DIHG4=ifelse(DIHG4h<0,0,DIHG4h)

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
    DIHQE4h=det(-solve(optQE4_Q$hessian))
    DIHQE4=ifelse(DIHQE4h<0,0,DIHQE4h)

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
    DIHP4h=det(-solve(optP4_Q$hessian))
    DIHP4=ifelse(DIHP4h<0,0,DIHP4h)

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
    DIHL4h=det(-solve(optL4_Q$hessian))
    DIHL4=ifelse(DIHL4h<0,0,DIHL4h)

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
    lpwlp=(prior.weights*w)/sum(prior.weights*w)
  }
  names(lpwlp) = c("E4_Q","IE4_Q","H4_Q","LN4_Q","G4_Q","QE4_Q","P4_Q","L4_Q")

  # the model average posterior as a mixture
  count=round(lpwlp*ndraws)
  mabmd=(c(# normal
    if(prior.weights[1]>0) sample(as.matrix(fitstanE4_Q)[,2],count[1],replace=T),
    if(prior.weights[2]>0) sample(as.matrix(fitstanIE4_Q)[,2],count[2],replace=T),
    if(prior.weights[3]>0) sample(as.matrix(fitstanH4_Q)[,2],count[3],replace=T),
    if(prior.weights[4]>0) sample(as.matrix(fitstanLN4_Q)[,2],count[4],replace=T),
    if(prior.weights[5]>0) sample(as.matrix(fitstanG4_Q)[,2],count[5],replace=T),
    if(prior.weights[6]>0) sample(as.matrix(fitstanQE4_Q)[,2],count[6],replace=T),
    if(prior.weights[7]>0) sample(as.matrix(fitstanP4_Q)[,2],count[7],replace=T),
    if(prior.weights[8]>0) sample(as.matrix(fitstanL4_Q)[,2],count[8],replace=T)

  ))

  macilp=quantile(mabmd,pvec)*data$maxD
  names(macilp)=c("BMDL","BMD","BMDU") # on original scale

  if(TRUE %in% (mabmd > data$maxD) && data$maxD > 1){
    mabmd = ifelse(mabmd > data$maxD, data$maxD, mabmd)
    p.msg = 'The model averaged posterior distribution has been truncated at max(Dose)^2'
    warnings('The model averaged posterior distribution has been truncated at max(Dose)^2')
  }else{
    p.msg = ''
  }

  BMDq_ls = quantile(mabmd, seq(0,1,0.005))*data$maxD


  ### Model-averaged response per dose level
  dr.MA.ls <- c()
  for(i in 1:length(data$x)){
    dr.MA.ls[i] = weighted.mean(x = c(DRM_E4_Q[i],DRM_IE4_Q[i],DRM_H4_Q[i],DRM_LN4_Q[i],DRM_G4_Q[i],
                                      DRM_QE4_Q[i], DRM_P4_Q[i],
                                      DRM_L4_Q[i]),
                                w = lpwlp,
                                na.rm = T)
  }


  ##### Weights & MA if one of the models is divergent --> this model gets weight 0
  if((0 %in% converged) && (1 %in% converged)){

    prior.weights.new = prior.weights
    div.models = which(converged == 0)
    prior.weights.new[div.models] = 0
    p.weights.new = prior.weights.new/sum(prior.weights.new==1)

    lpwlp.conv=(p.weights.new*w)/sum(p.weights.new*w)
    # lpwlp = lpwlp*prior.weights

    # the model average posterior as a mixture
    count=round(lpwlp.conv*ndraws)
    mabmd.conv=(c(# normal
      if(p.weights.new[1]>0) sample(as.matrix(fitstanE4_Q)[,2],count[1],replace=T),
      if(p.weights.new[2]>0) sample(as.matrix(fitstanIE4_Q)[,2],count[2],replace=T),
      if(p.weights.new[3]>0) sample(as.matrix(fitstanH4_Q)[,2],count[3],replace=T),
      if(p.weights.new[4]>0) sample(as.matrix(fitstanLN4_Q)[,2],count[4],replace=T),
      if(p.weights.new[5]>0) sample(as.matrix(fitstanG4_Q)[,2],count[5],replace=T),
      if(p.weights.new[6]>0) sample(as.matrix(fitstanQE4_Q)[,2],count[6],replace=T),
      if(p.weights.new[7]>0) sample(as.matrix(fitstanP4_Q)[,2],count[7],replace=T),
      if(p.weights.new[8]>0) sample(as.matrix(fitstanL4_Q)[,2],count[8],replace=T)
    ))

    macilp.conv = quantile(mabmd.conv,pvec)*data$maxD
    names(macilp.conv)=c("BMDL","BMD","BMDU")

    if(TRUE %in% (mabmd.conv > data$maxD) && data$maxD > 1){
      mabmd.conv = ifelse(mabmd.conv > data$maxD, data$maxD, mabmd.conv)
      p.msg = 'The model averaged posterior distribution has been truncated at max(Dose)^2'
      warnings('The model averaged posterior distribution has been truncated at max(Dose)^2')
    }else{
      p.msg = ''
    }

    BMDq_ls_conv = quantile(mabmd.conv, seq(0,1,0.005))*data$maxD


    ### Model-averaged response per dose level
    dr.MA.ls.conv <- c()
    for(i in 1:length(data$x)){
      dr.MA.ls.conv[i] = weighted.mean(x = c(DRM_E4_Q[i],DRM_IE4_Q[i],DRM_H4_Q[i],DRM_LN4_Q[i],DRM_G4_Q[i],
                                             DRM_QE4_Q[i],
                                             DRM_P4_Q[i], DRM_L4_Q[i]),
                                       w = lpwlp.conv,
                                       na.rm = T)
    }


  }else{
    lpwlp.conv = NULL; macilp.conv = NULL; BMDq_ls_conv = NULL; dr.MA.ls.conv = NULL; mabmd.conv = NA;
  }

  ## Plot with weights bridge sampling
  BMDL = c(BMDL, macib[1]); BMD = c(BMD, macib[2]); BMDU = c(BMDU, macib[3]) # all on original scale

  names(BMDL) <- c("E4_Q","IE4_Q","H4_Q","LN4_Q","G4_Q","QE4_Q","P4_Q","L4_Q", "MA")
  model = c("E4_Q","IE4_Q","H4_Q","LN4_Q","G4_Q","QE4_Q","P4_Q","L4_Q", "MA")
  model = as.factor(model)

  weight = c(rep(0,16),1)
  names(weight) = c("E4_Q","IE4_Q","H4_Q","LN4_Q","G4_Q","QE4_Q","P4_Q","L4_Q","MA")
  for(i in names(lpwb)){
    weight[names(weight)==i] = lpwb[names(lpwb)==i]
  }

  modelnames <- c("E4","IE4","H4","LN4","G4","QE4","P4","L4")

  ### Some additional checks


  if(macib[2]/macib[1] > 20){
    warning('BMD/BMDL is larger than 20 for bridge sampling')
  }
  if(macib[3]/macib[1] > 50){
    warning('BMDU/BMDL is larger than 50 for bridge sampling')
  }
  if(macib[2] < (data.Q$data$x[2]*data.Q$data$maxD/10)){
    warning('BMD is 10 times lower than the lowest non-zero dose for bridge sampling')
  }
  if(macilp[2]/macilp[1] > 20){
    warning('BMD/BMDL is larger than 20 for hybrid Laplace')
  }
  if(macilp[3]/macilp[1] > 50){
    warning('BMDU/BMDL is larger than 50 for hybrid Laplace')
  }
  if(macilp[2] < (data.Q$data$x[2]*data.Q$data$maxD/10)){
    warning('BMD is 10 times lower than the lowest non-zero dose for hybrid Laplace')
  }

  ### best fitting model vs saturated ANOVA model
  best.fit = modelnames[which(weight[1:8] == max(weight[1:8]))]

  bfTest <- modelTestQ(best.fit, data.Q, get(paste0('fitstan', best.fit, '_Q')), type = 'MCMC',
                       seed, ndraws, nrchains, nriterations, warmup, delta, treedepth)
  #print(warning(bfTest$warn.bf))

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

  modelnames2 <- c("E4_Q","IE4_Q","H4_Q","LN4_Q","G4_Q","QE4_Q","P4_Q","L4_Q")

  ret_results <- list(E4_Q=E4outQ,IE4_Q=IE4outQ,H4_Q=H4outQ,LN4_Q=LN4outQ,
                      G4_Q=G4outQ,QE4_Q=QE4outQ,P4_Q=P4outQ,L4_Q=L4outQ,
                      covs = covs, corrs = corrs,
                      weights_bridge_sampling=w.bs,
                      weights_laplace=lpwlp,
                      MA_bridge_sampling=macib,
                      MA_laplace=macilp,
                      llQ=llQ,
                      convergence=converged, bs_weights_conv=w.bs.conv,
                      ls_weights_conv=lpwlp.conv,
                      MA_bs_conv=macib.conv,
                      MA_ls_conv=macilp.conv,
                      MA_post_bs = BMDq_bs,
                      MA_post_ls = BMDq_ls,
                      MA_post_bs_conv = BMDq_bs_conv,
                      MA_post_ls_conv = BMDq_ls_conv,
                      MA_dr_bs = dr.MA.bs,
                      MA_dr_ls = dr.MA.ls,
                      MA_dr_bs_conv = dr.MA.bs.conv,
                      MA_dr_ls_conv = dr.MA.ls.conv,
                      parsQ = list(parsE4Q, parsIE4Q, parsH4Q, parsLN4Q, parsG4Q,
                                   parsQE4Q, parsP4Q, parsL4Q),
                      BMDMixture = mabmd*data$maxD,
                      BMDMixture.conv = mabmd.conv*data$maxD,
                      BMDMixtureBS = (mabmd1)*data$maxD,
                      BMDMixture.convBS = mabmd.conv1*data$maxD,
                      divergences = divergences,
                      data = data.frame(
                        dose = c(data.Q$data$x),
                        y = data.Q$data$y,
                        n = data.Q$data$n),
                      max.dose = data.Q$data$maxD,
                      q = data.Q$data$q,
                      models_included_bridge = modelnames2[p.weights > 0],
                      models_included_laplace = modelnames2[lpwlp > 0],
                      is_bin = data.Q$data$is_bin,
                      is_betabin = data.Q$data$is_betabin,
                      bf = bfTest$bayesFactor, gof_check = bfTest$warn.bf,
                      means.SM = bfTest$means.SM, parBestFit = bfTest$par.best,
                      BIC.bestfit = bfTest$BIC.bestfit, BIC.SM = bfTest$BIC.SM,
                      w.msg = w.msg, p.msg = p.msg
  )

  attr(ret_results, "class") <- c("BMADRQ", "BS")

  return(ret_results)

}


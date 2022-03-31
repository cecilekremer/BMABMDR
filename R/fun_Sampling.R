#' Perform model averaging using MCMC methods (Bridge sampling & Partial Laplace)
#'
#' This method assumed data for continuous endpoints.
#'
#' More detailed descriprion
#' @param data.N the input data as returned by function PREP_DATA_N
#' @param data.LN the input data as returned by function PREP_DATA_LN
#' @param prior.weights a vector specifying which of the 16 models should be included (1 = include, 0 = exclude)
#' @param priordist which prior distribution should be used (currently only "PERT" is implemented)
#' @param prior.BMD logical indicating whether an informative prior for the BMD is used (currently not implemented)
#' @param ndraws the number of draws from the posterior, default 30000
#' @param nrchains the number of chains to be used in the MCMC
#' @param nriterations the number of iterations per chain
#' @param warmup the number of iterations per chain to be discarded as burnin
#' @param delta default 0.8
#' @param treedepth default 10
#' @param seed default 123
#' @param pvec vector specifying the three BMD quantiles of interest
#' @param plot logical indicating whether a simple plot of model fits should be shown
#'
#' @description Using MCMC, we compute the parameters of each model, perform model averaging using bridge sampling.
#'              We also implemented a Laplace-MCMC method within this function where the model parameters are estimated
#'              with MCMC, but the model weights are computed using Laplace approximation.
#'
#' @examples
#'  # we use the first 5 rows because those are observations from subjects belonging to the same group.
#'  data("immunotoxicityData.rda")  #load the immunotoxicity data
#'  data_N <- PREP_DATA_N(data = as.data.frame(immunotoxicityData[1:5,]),
#'                        sumstats = TRUE, sd = TRUE, q = 0.1) #example with default priors
#'  data_LN <- PREP_DATA_LN(data = as.data.frame(immunotoxicityData[1:5,]),
#'                          sumstats = TRUE, sd = TRUE, q = 0.1) #example with default priors
#'  pvec <- c(0.05, 0.5, 0.95)
#'  prior.weights = rep(1, 16)
#'  nrch=3;nriter=3000;wu=1000;dl=0.8;trd=10;sd=123;ndr=30000
#'  SBMD = sampling_MA(data_N,data_LN,prior.weights,
#'                     ndraws=ndr,nrchains=nrch,
#'                     nriterations=nriter,warmup=wu,delta=dl,
#'                    treedepth=trd,seed=sd,pvec=pvec)
#'
#'
#' @return a list containing the following important entries:
#' \enumerate{
#'   \item E4_N parameter estimates from the exponential model
#'   \item IE4_N parameter estimates from the inverse-exponential model
#'   \item H4_N parameter estimates from the Hill model
#'   \item LN4_N parameter estimates from the lognormal model
#'   \item G4_N parameter estimates from the gamma model
#'   \item QE4_N parameter estimates from the quadratic-exponential model
#'   \item P4_N parameter estimates from the probit model
#'   \item L4_N parameter estimates from the logit model
#'   \item E4_LN parameter estimates from the exponential model
#'   \item IE4_LN parameter estimates from the inverse-exponential model
#'   \item H4_LN parameter estimates from the Hill model
#'   \item LN4_LN parameter estimates from the lognormal model
#'   \item G4_LN parameter estimates from the gamma model
#'   \item QE4_LN parameter estimates from the quadratic-exponential model
#'   \item P4_LN parameter estimates from the probit model
#'   \item L4_LN parameter estimates from the logit model
#'   \item MA_laplace Laplace approximation model averaged BMD estimates using all the models
#'   \item MA_bs Bridge sampling model averaged BMD estimates using all the models
#'   \item MA_bs_conv Bridge sampling model averaged BMD estimates using only converged models
#'   \item weights_laplace model weights using Laplace approximation to the posterior
#'   \item weights_bs model weights computed using bridge sampling
#'   \item convergence vector indicating model convergence or not. 1 = converged, 0 otherwise.
#'   \item llN vector of model likelihoods when the distribution is assumed normal
#'    \item llLN vector of model likelihoods when the distribution is assumed lognormal
#'   \item bf Bayes factor comparing the best model against saturated ANOVA model
#' }
#' @export sampling_MA

sampling_MA=function(data.N,data.LN,prior.weights,priordist = "PERT", prior.BMD = FALSE,
                     ndraws = 30000,nrchains=3,
                     nriterations=3000,warmup=1000,
                     delta=0.8,treedepth=10,seed=123,pvec=c(0.05,0.5,0.95),
                     plot = FALSE){

  q = data.N$data.modstan1a$q

  if(priordist == "Normal" & prior.BMD == TRUE){
    stop("Informative prior on BMD can only be used in combination with the PERT prior distribution.")
  }

  if(data.N$increasing == TRUE){

    ## Normal distribution

    data.modstan1a=data.N$data.modstan1a
    data.modstan1bG=data.N$data.modstan1bG
    data.modstan1bQ=data.N$data.modstan1bQ
    data.modstan1aLN=data.N$data.modstan1aLN
    data.modstan2L=data.N$data.modstan2L
    data.modstan2P=data.N$data.modstan2P
    data.modstanE=data.N$data.modstanE
    svF1=data.N$svF1
    svF1LN = data.N$svF1LN
    svF2P=data.N$svF2P
    svF2L=data.N$svF2L

    BMDL=c(); BMD=c(); BMDU=c()
    converged=c()

    ## Obtain model parameters via MCMC sampling
    if(prior.weights[1]>0){
      print(1)

      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          fitstanE4_N = fun_sampling(stanmodels$mE4_NI, data.modstan1a, svF1,
                                     ndraws,nrchains,
                                     nriterations,warmup,
                                     delta,treedepth,seed,pvec)
        }else if(prior.BMD == TRUE){
          fitstanE4_N = fun_sampling(stanmodels$mE4_NI2, data.modstan1a, svF1,
                                     ndraws,nrchains,
                                     nriterations,warmup,
                                     delta,treedepth,seed,pvec)
        }
      }else if(priordist == "Normal"){
        fitstanE4_N = fun_sampling(stanmodels$mE4_NI_n, data.modstan1a, svF1,
                                   ndraws,nrchains,
                                   nriterations,warmup,
                                   delta,treedepth,seed,pvec)
      }

      parsE4N <- par_extract(fitstanE4_N, model_name = "E4_N")
      # parsE4N[,"BMD"] <- parsE4N[,"BMD"]*data.modstan1a$maxD # BMD on original scale

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

      E4resNI=exp(quantile(as.matrix(fitstanE4_N)[,"par[2]"],pvec))*data.modstan1a$maxD
      E4resNI=c(E4resNI,apply(as.matrix(fitstanE4_N),2,median)[c("par[1]","par[2]","par[3]","par[4]","par[5]")])
      names(E4resNI)=c("BMDL","BMD","BMDU","at","k","ct","dt","is2t")

      E4outNI <- outLP(parsE4N, pvec, data.modstan1a$maxD)

      DRM_E4_N = DRM.E4_NI(E4resNI[4:7], data.modstan1a$x, data.modstan1a$q)


      # Covariance between b-d and between BMD-d
      E4covNI = c(cov(as.matrix(fitstanE4_N)[,c("b","d")], use="complete.obs")["b","d"],
                  cov(as.matrix(fitstanE4_N)[,c("BMD","d")], use="complete.obs")["BMD","d"])

      E4corrNI = c(cor(as.matrix(fitstanE4_N)[,c("b","d")], use="complete.obs")["b","d"],
                   cor(as.matrix(fitstanE4_N)[,c("BMD","d")], use="complete.obs")["BMD","d"])

      BMDL=c(BMDL, E4resNI[1]); BMD=c(BMD, E4resNI[2]); BMDU=c(BMDU, E4resNI[3])
      bridgeE4N = bridgesampling::bridge_sampler(fitstanE4_N, silent = T) # compute log marginal likelihood
    }else{E4resNI=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeE4N=NA; converged=c(converged, NA); E4covNI=rep(NA,2); E4corrNI=rep(NA,2); DRM_E4_N=rep(NA,length(data.modstan1a$x))
    parsE4N <- NA
    div_E4_N <- NA
    E4outNI <- t(data.frame(
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

    if(prior.weights[2]>0){
      print(2)
      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          fitstanIE4_N = fun_sampling(stanmodels$mIE4_NI, data.modstan1a, svF1,
                                      ndraws,nrchains,
                                      nriterations,warmup,
                                      delta,treedepth,seed,pvec)
        }else if(prior.BMD == TRUE){
          fitstanIE4_N = fun_sampling(stanmodels$mIE4_NI2, data.modstan1a, svF1,
                                      ndraws,nrchains,
                                      nriterations,warmup,
                                      delta,treedepth,seed,pvec)
        }
      }else if(priordist == "Normal"){
        fitstanIE4_N = fun_sampling(stanmodels$mIE4_NI_n, data.modstan1a, svF1,
                                    ndraws,nrchains,
                                    nriterations,warmup,
                                    delta,treedepth,seed,pvec)
      }
      parsIE4N <- par_extract(fitstanIE4_N, model_name = "IE4_N")
      # parsIE4N[,"BMD"] <- parsIE4N[,"BMD"]*data.modstan1a$maxD # BMD on original scale

      #diagnostics here
      # posterior_diag(model_stan = fitstanIE4_N)
      convergence_stat <- convergence_deci(model_stan = fitstanIE4_N, nrchains = nrchains)
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)

      # divergent transitions
      div_IE4_N <- sum(sapply(rstan::get_sampler_params(fitstanIE4_N, inc_warmup = F),
                              function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))

      IE4resNI=exp(quantile(as.matrix(fitstanIE4_N)[,"par[2]"],pvec))*data.modstan1a$maxD
      IE4resNI=c(IE4resNI,apply(as.matrix(fitstanIE4_N),2,median)[c("par[1]","par[2]","par[3]","par[4]","par[5]")])
      names(IE4resNI)=c("BMDL","BMD","BMDU","at","k","ct","dt","is2t")

      IE4outNI <- outLP(parsIE4N, pvec, data.modstan1a$maxD)

      DRM_IE4_N = DRM.IE4_NI(IE4resNI[4:7], data.modstan1a$x, data.modstan1a$q)


      # Covariance between b-d and between BMD-d
      IE4covNI = c(cov(as.matrix(fitstanIE4_N)[,c("b","d")], use="complete.obs")["b","d"],
                   cov(as.matrix(fitstanIE4_N)[,c("BMD","d")], use="complete.obs")["BMD","d"])

      IE4corrNI = c(cor(as.matrix(fitstanIE4_N)[,c("b","d")], use="complete.obs")["b","d"],
                    cor(as.matrix(fitstanIE4_N)[,c("BMD","d")], use="complete.obs")["BMD","d"])

      BMDL=c(BMDL, IE4resNI[1]); BMD=c(BMD, IE4resNI[2]); BMDU=c(BMDU, IE4resNI[3])
      bridgeIE4N = bridgesampling::bridge_sampler(fitstanIE4_N, silent = T)

    }else{IE4resNI=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeIE4N=NA; converged=c(converged, NA); IE4covNI=rep(NA,2); IE4corrNI=rep(NA,2); DRM_IE4_N=rep(NA,length(data.modstan1a$x))
    parsIE4N <- NA
    div_IE4_N <- NA

    IE4outNI <- t(data.frame(
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

    if(prior.weights[3]>0){
      print(3)
      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          fitstanH4_N = fun_sampling(stanmodels$mH4_NI, data.modstan1a, svF1,
                                     ndraws,nrchains,
                                     nriterations,warmup,
                                     delta,treedepth,seed,pvec)
        }else if(prior.BMD == TRUE){
          fitstanH4_N = fun_sampling(stanmodels$mH4_NI2, data.modstan1a, svF1,
                                     ndraws,nrchains,
                                     nriterations,warmup,
                                     delta,treedepth,seed,pvec)
        }
      }else if(priordist == "Normal"){
        fitstanH4_N = fun_sampling(stanmodels$mH4_NI_n, data.modstan1a, svF1,
                                   ndraws,nrchains,
                                   nriterations,warmup,
                                   delta,treedepth,seed,pvec)
      }
      parsH4N <- par_extract(fitstanH4_N, model_name = "H4_N")
      # parsH4N[,"BMD"] <- parsH4N[,"BMD"]*data.modstan1a$maxD # BMD on original scale

      #diagnostics here
      # posterior_diag(model_stan = fitstanH4_N)
      convergence_stat <- convergence_deci(model_stan = fitstanH4_N, nrchains = nrchains)
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)

      # divergent transitions
      div_H4_N <- sum(sapply(rstan::get_sampler_params(fitstanH4_N, inc_warmup = F),
                             function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))

      H4resNI=exp(quantile(as.matrix(fitstanH4_N)[,"par[2]"],pvec))*data.modstan1a$maxD
      H4resNI=c(H4resNI,apply(as.matrix(fitstanH4_N),2,median)[c("par[1]","par[2]","par[3]","par[4]","par[5]")])
      names(H4resNI)=c("BMDL","BMD","BMDU","at","k","ct","dt","is2t")

      H4outNI <- outLP(parsH4N, pvec, data.modstan1a$maxD)


      DRM_H4_N = DRM.H4_NI(H4resNI[4:7], data.modstan1a$x, data.modstan1a$q)

      # Covariance between b-d and between BMD-d
      H4covNI = c(cov(as.matrix(fitstanH4_N)[,c("b","d")], use="complete.obs")["b","d"],
                  cov(as.matrix(fitstanH4_N)[,c("BMD","d")], use="complete.obs")["BMD","d"])

      H4corrNI = c(cor(as.matrix(fitstanH4_N)[,c("b","d")], use="complete.obs")["b","d"],
                   cor(as.matrix(fitstanH4_N)[,c("BMD","d")], use="complete.obs")["BMD","d"])

      BMDL=c(BMDL, H4resNI[1]); BMD=c(BMD, H4resNI[2]); BMDU=c(BMDU, H4resNI[3])
      bridgeH4N = bridge_sampler(fitstanH4_N, silent = T)
    }else{H4resNI=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeH4N=NA;
    converged=c(converged, NA); H4covNI=rep(NA,2); H4corrNI=rep(NA,2); DRM_H4_N=rep(NA,length(data.modstan1a$x))
    parsH4N <- NA
    div_H4_N <- NA

    H4outNI <- t(data.frame(
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

    if(prior.weights[4]>0){
      print(4)
      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          fitstanLN4_N = fun_sampling(stanmodels$mLN4_NI, data.modstan1aLN, svF1LN,
                                      ndraws,nrchains,
                                      nriterations,warmup,
                                      delta,treedepth,seed,pvec)
        }else if(prior.BMD == TRUE){
          fitstanLN4_N = fun_sampling(stanmodels$mLN4_NI2, data.modstan1aLN, svF1LN,
                                      ndraws,nrchains,
                                      nriterations,warmup,
                                      delta,treedepth,seed,pvec)
        }
      }else if(priordist == "Normal"){
        fitstanLN4_N = fun_sampling(stanmodels$mLN4_NI_n, data.modstan1aLN, svF1LN,
                                    ndraws,nrchains,
                                    nriterations,warmup,
                                    delta,treedepth,seed,pvec)
      }
      parsLN4N <- par_extract(fitstanLN4_N, model_name = "LN4_N")
      # parsLN4N[,"BMD"] <- parsLN4N[,"BMD"]*data.modstan1a$maxD # BMD on original scale

      #diagnostics here
      # posterior_diag(model_stan = fitstanLN4_N)
      convergence_stat <- convergence_deci(model_stan = fitstanLN4_N, nrchains = nrchains)
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)

      # divergent transitions
      div_LN4_N <- sum(sapply(rstan::get_sampler_params(fitstanLN4_N, inc_warmup = F),
                              function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))

      LN4resNI=exp(quantile(as.matrix(fitstanLN4_N)[,"par[2]"],pvec))*data.modstan1a$maxD
      LN4resNI=c(LN4resNI,apply(as.matrix(fitstanLN4_N),2,median)[c("par[1]","par[2]","par[3]","par[4]","par[5]")])
      names(LN4resNI)=c("BMDL","BMD","BMDU","at","k","ct","dt","is2t")

      LN4outNI <- outLP(parsLN4N, pvec, data.modstan1a$maxD)


      DRM_LN4_N = DRM.LN4_NI(LN4resNI[4:7], data.modstan1a$x, data.modstan1a$q)

      # Covariance between b-d and between BMD-d
      LN4covNI = c(cov(as.matrix(fitstanLN4_N)[,c("b","d")], use="complete.obs")["b","d"],
                   cov(as.matrix(fitstanLN4_N)[,c("BMD","d")], use="complete.obs")["BMD","d"])

      LN4corrNI = c(cor(as.matrix(fitstanLN4_N)[,c("b","d")], use="complete.obs")["b","d"],
                    cor(as.matrix(fitstanLN4_N)[,c("BMD","d")], use="complete.obs")["BMD","d"])

      BMDL=c(BMDL, LN4resNI[1]); BMD=c(BMD, LN4resNI[2]); BMDU=c(BMDU, LN4resNI[3])
      bridgeLN4N = bridgesampling::bridge_sampler(fitstanLN4_N, silent = T)
    }else{LN4resNI=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeLN4N=NA; converged=c(converged, NA); LN4covNI=rep(NA,2); LN4corrNI=rep(NA,2); DRM_LN4_N=rep(NA,length(data.modstan1aLN$x))
    parsLN4N <- NA
    div_LN4_N <- NA

    LN4outNI <- t(data.frame(
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

    if(prior.weights[5]>0){
      print(5)
      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          fitstanG4_N = fun_sampling(stanmodels$mG4_NI, data.modstan1bG, svF1,
                                     ndraws,nrchains,
                                     nriterations,warmup,
                                     delta,treedepth,seed,pvec)
        }else if(prior.BMD == TRUE){
          fitstanG4_N = fun_sampling(stanmodels$mG4_NI2, data.modstan1bG, svF1,
                                     ndraws,nrchains,
                                     nriterations,warmup,
                                     delta,treedepth,seed,pvec)
        }
      }else if(priordist == "Normal"){
        fitstanG4_N = fun_sampling(stanmodels$mG4_NI_n, data.modstan1bG, svF1,
                                   ndraws,nrchains,
                                   nriterations,warmup,
                                   delta,treedepth,seed,pvec)
      }

      parsG4N <- par_extract(fitstanG4_N, model_name = "G4_N")
      # parsG4N[,"BMD"] <- parsG4N[,"BMD"]*data.modstan1a$maxD # BMD on original scale

      #diagnostics here
      # posterior_diag(model_stan = fitstanG4_N)
      convergence_stat <- convergence_deci(model_stan = fitstanG4_N, nrchains = nrchains)
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)

      # divergent transitions
      div_G4_N <- sum(sapply(rstan::get_sampler_params(fitstanG4_N, inc_warmup = F),
                             function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))

      G4resNI=exp(quantile(as.matrix(fitstanG4_N)[,"par[2]"],pvec))*data.modstan1a$maxD
      G4resNI=c(G4resNI,apply(as.matrix(fitstanG4_N),2,median)[c("par[1]","par[2]","par[3]","par[4]","par[5]")])
      names(G4resNI)=c("BMDL","BMD","BMDU","at","k","ct","dt","is2t")

      G4outNI <- outLP(parsG4N, pvec, data.modstan1a$maxD)

      DRM_G4_N = DRM.G4_NI(G4resNI[4:7], data.modstan1a$x, data.modstan1a$q)


      # Covariance between b-d and between BMD-d
      G4covNI = c(cov(as.matrix(fitstanG4_N)[,c("b","d")], use="complete.obs")["b","d"],
                  cov(as.matrix(fitstanG4_N)[,c("BMD","d")], use="complete.obs")["BMD","d"])

      G4corrNI = c(cor(as.matrix(fitstanG4_N)[,c("b","d")], use="complete.obs")["b","d"],
                   cor(as.matrix(fitstanG4_N)[,c("BMD","d")], use="complete.obs")["BMD","d"])

      BMDL=c(BMDL, G4resNI[1]); BMD=c(BMD, G4resNI[2]); BMDU=c(BMDU, G4resNI[3])
      bridgeG4N = bridgesampling::bridge_sampler(fitstanG4_N, silent = T)
    }else{G4resNI=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeG4N=NA; converged=c(converged, NA); G4covNI=rep(NA,2); G4corrNI=rep(NA,2); DRM_G4_N=rep(NA,length(data.modstan1bG$x))
    parsG4N <- NA
    div_G4_N <- NA

    G4outNI <- t(data.frame(
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

    if(prior.weights[6]>0){
      print(6)
      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          fitstanQE4_N = fun_sampling(stanmodels$mQE4_NI, data.modstan1bQ, svF1,
                                      ndraws,nrchains,
                                      nriterations,warmup,
                                      delta,treedepth,seed,pvec)
        }else if(prior.BMD == TRUE){
          fitstanQE4_N = fun_sampling(stanmodels$mQE4_NI2, data.modstan1bQ, svF1,
                                      ndraws,nrchains,
                                      nriterations,warmup,
                                      delta,treedepth,seed,pvec)
        }
      }else if(priordist == "Normal"){
        fitstanQE4_N = fun_sampling(stanmodels$mQE4_NI_n, data.modstan1bQ, svF1,
                                    ndraws,nrchains,
                                    nriterations,warmup,
                                    delta,treedepth,seed,pvec)
      }

      parsQE4N <- par_extract(fitstanQE4_N, model_name = "QE4_N")
      # parsQE4N[,"BMD"] <- parsQE4N[,"BMD"]*data.modstan1a$maxD # BMD on original scale

      #diagnostics here
      # posterior_diag(model_stan = fitstanQE4_N)
      convergence_stat <- convergence_deci(model_stan = fitstanQE4_N, nrchains = nrchains)
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)

      # divergent transitions
      div_QE4_N <- sum(sapply(rstan::get_sampler_params(fitstanQE4_N, inc_warmup = F),
                              function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))

      QE4resNI=exp(quantile(as.matrix(fitstanQE4_N)[,"par[2]"],pvec))*data.modstan1a$maxD
      QE4resNI=c(QE4resNI,apply(as.matrix(fitstanQE4_N),2,median)[c("par[1]","par[2]","par[3]","par[4]","par[5]")])
      names(QE4resNI)=c("BMDL","BMD","BMDU","at","k","ct","dt","is2t")

      QE4outNI <- outLP(parsQE4N, pvec, data.modstan1a$maxD)

      DRM_QE4_N = DRM.QE4_NI(QE4resNI[4:7], data.modstan1a$x, data.modstan1a$q)


      # Covariance between b-d and between BMD-d
      QE4covNI = c(cov(as.matrix(fitstanQE4_N)[,c("b","d")], use="complete.obs")["b","d"],
                   cov(as.matrix(fitstanQE4_N)[,c("BMD","d")], use="complete.obs")["BMD","d"])

      QE4corrNI = c(cor(as.matrix(fitstanQE4_N)[,c("b","d")], use="complete.obs")["b","d"],
                    cor(as.matrix(fitstanQE4_N)[,c("BMD","d")], use="complete.obs")["BMD","d"])

      BMDL=c(BMDL, QE4resNI[1]); BMD=c(BMD, QE4resNI[2]); BMDU=c(BMDU, QE4resNI[3])
      bridgeQE4N = bridgesampling::bridge_sampler(fitstanQE4_N, silent = T)
    }else{QE4resNI=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeQE4N=NA; converged=c(converged, NA); QE4covNI=rep(NA,2); QE4corrNI=rep(NA,2); DRM_QE4_N=rep(NA,length(data.modstan1bQ$x))
    parsQE4N <- NA
    div_QE4_N <- NA

    QE4outNI <- t(data.frame(
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

    if(prior.weights[7]>0){
      print(7)
      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          fitstanP4_N = fun_sampling(stanmodels$mP4_NI, data.modstan2P, svF2P,
                                     ndraws,nrchains,
                                     nriterations,warmup,
                                     delta,treedepth,seed,pvec)
        }else if(prior.BMD == TRUE){
          fitstanP4_N = fun_sampling(stanmodels$mP4_NI2, data.modstan2P, svF2P,
                                     ndraws,nrchains,
                                     nriterations,warmup,
                                     delta,treedepth,seed,pvec)
        }
      }else if(priordist == "Normal"){
        fitstanP4_N = fun_sampling(stanmodels$mP4_NI_n, data.modstan2P, svF2P,
                                   ndraws,nrchains,
                                   nriterations,warmup,
                                   delta,treedepth,seed,pvec)
      }
      parsP4N <- par_extract(fitstanP4_N, model_name = "P4_N")
      # parsP4N[,"BMD"] <- parsP4N[,"BMD"]*data.modstan1a$maxD # BMD on original scale

      #diagnostics here
      # posterior_diag(model_stan = fitstanP4_N)
      convergence_stat <- convergence_deci(model_stan = fitstanP4_N, nrchains = nrchains)
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)

      # divergent transitions
      div_P4_N <- sum(sapply(rstan::get_sampler_params(fitstanP4_N, inc_warmup = F),
                             function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))

      P4resNI=exp(quantile(as.matrix(fitstanP4_N)[,"par[2]"],pvec))*data.modstan1a$maxD
      P4resNI=c(P4resNI,apply(as.matrix(fitstanP4_N),2,median)[c("par[1]","par[2]","par[3]","par[4]","par[5]")])
      names(P4resNI)=c("BMDL","BMD","BMDU","at","k","ct","dt","is2t")

      P4outNI <- outLP(parsP4N, pvec, data.modstan1a$maxD)


      DRM_P4_N = DRM.P4_NI(P4resNI[4:7], data.modstan1a$x, data.modstan1a$q)


      # Covariance between b-d and between BMD-d
      P4covNI = c(cov(as.matrix(fitstanP4_N)[,c("b","d")], use="complete.obs")["b","d"],
                  cov(as.matrix(fitstanP4_N)[,c("BMD","d")], use="complete.obs")["BMD","d"])

      P4corrNI = c(cor(as.matrix(fitstanP4_N)[,c("b","d")], use="complete.obs")["b","d"],
                   cor(as.matrix(fitstanP4_N)[,c("BMD","d")], use="complete.obs")["BMD","d"])

      BMDL=c(BMDL, P4resNI[1]); BMD=c(BMD, P4resNI[2]); BMDU=c(BMDU, P4resNI[3])
      bridgeP4N = bridge_sampler(fitstanP4_N, silent = T)
    }else{P4resNI=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeP4N=NA; converged=c(converged, NA);
    P4covNI=rep(NA,2); P4corrNI=rep(NA,2);
    DRM_P4_N=rep(NA,length(data.modstan2P$x))
    parsP4N <- NA
    div_P4_N <- NA

    P4outNI <- t(data.frame(
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

    if(prior.weights[8]>0){
      print(8)
      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          fitstanL4_N = fun_sampling(stanmodels$mL4_NI, data.modstan2L, svF2L,
                                     ndraws,nrchains,
                                     nriterations,warmup,
                                     delta,treedepth,seed,pvec)
        }else if(prior.BMD == TRUE){
          fitstanL4_N = fun_sampling(stanmodels$mL4_NI2, data.modstan2L, svF2L,
                                     ndraws,nrchains,
                                     nriterations,warmup,
                                     delta,treedepth,seed,pvec)
        }
      }else if(priordist == "Normal"){
        fitstanL4_N = fun_sampling(stanmodels$mL4_NI_n, data.modstan2L, svF2L,
                                   ndraws,nrchains,
                                   nriterations,warmup,
                                   delta,treedepth,seed,pvec)
      }
      parsL4N <- par_extract(fitstanL4_N, model_name = "L4_N")
      # parsL4N[,"BMD"] <- parsL4N[,"BMD"]*data.modstan1a$maxD # BMD on original scale

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

      L4resNI=exp(quantile(as.matrix(fitstanL4_N)[,"par[2]"],pvec))*data.modstan1a$maxD
      L4resNI=c(L4resNI,apply(as.matrix(fitstanL4_N),2,median)[c("par[1]","par[2]","par[3]","par[4]","par[5]")])
      names(L4resNI)=c("BMDL","BMD","BMDU","at","k","ct","dt","is2t")

      L4outNI <- outLP(parsL4N, pvec, data.modstan1a$maxD)


      DRM_L4_N = DRM.L4_NI(L4resNI[4:7], data.modstan1a$x, data.modstan1a$q)


      # Covariance between b-d and between BMD-d
      L4covNI = c(cov(as.matrix(fitstanL4_N)[,c("b","d")], use="complete.obs")["b","d"],
                  cov(as.matrix(fitstanL4_N)[,c("BMD","d")], use="complete.obs")["BMD","d"])

      L4corrNI = c(cor(as.matrix(fitstanL4_N)[,c("b","d")], use="complete.obs")["b","d"],
                   cor(as.matrix(fitstanL4_N)[,c("BMD","d")], use="complete.obs")["BMD","d"])

      BMDL=c(BMDL, L4resNI[1]); BMD=c(BMD, L4resNI[2]); BMDU=c(BMDU, L4resNI[3])
      bridgeL4N = bridgesampling::bridge_sampler(fitstanL4_N, silent = T)
    }else{L4resNI=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeL4N=NA; converged=c(converged, NA); L4covNI=rep(NA,2); L4corrNI=rep(NA,2); DRM_L4_N=rep(NA,length(data.modstan2L$x))
    parsL4N <- NA
    div_L4_N <- NA

    L4outNI <- t(data.frame(
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

    ## Lognormal distribution

    data.modstan1a=data.LN$data.modstan1a
    data.modstan1bG=data.LN$data.modstan1bG
    data.modstan1bQ=data.LN$data.modstan1bQ
    data.modstan1aLN=data.LN$data.modstan1aLN
    data.modstan2L=data.LN$data.modstan2L
    data.modstan2P=data.LN$data.modstan2P
    data.modstanE=data.LN$data.modstanE
    svF1=data.LN$svF1
    svF1LN=data.LN$svF1LN
    svF2P=data.LN$svF2P
    svF2L=data.LN$svF2L

    # model specific results
    if(prior.weights[9]>0){
      print(9)
      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          fitstanE4_LN = fun_sampling(stanmodels$mE4_LNI, data.modstan1a, svF1,
                                      ndraws,nrchains,
                                      nriterations,warmup,
                                      delta,treedepth,seed,pvec)
        }else if(prior.BMD == TRUE){
          fitstanE4_LN = fun_sampling(stanmodels$mE4_LNI2, data.modstan1a, svF1,
                                      ndraws,nrchains,
                                      nriterations,warmup,
                                      delta,treedepth,seed,pvec)
        }
      }else if(priordist == "Normal"){
        fitstanE4_LN = fun_sampling(stanmodels$mE4_LNI_n, data.modstan1a, svF1,
                                    ndraws,nrchains,
                                    nriterations,warmup,
                                    delta,treedepth,seed,pvec)
      }

      parsE4LN <- par_extract(fitstanE4_LN, model_name = "E4_LN")
      # parsE4LN[,"BMD"] <- parsE4LN[,"BMD"]*data.modstan1a$maxD # BMD on original scale

      #diagnostics here
      # posterior_diag(model_stan = fitstanE4_LN)
      convergence_stat <- convergence_deci(model_stan = fitstanE4_LN, nrchains = nrchains)
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)

      # divergent transitions
      div_E4_LN <- sum(sapply(rstan::get_sampler_params(fitstanE4_LN, inc_warmup = F),
                              function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))

      E4resLNI=exp(quantile(as.matrix(fitstanE4_LN)[,"par[2]"],pvec))*data.modstan1a$maxD
      E4resLNI=c(E4resLNI,apply(as.matrix(fitstanE4_LN),2,median)[c("par[1]","par[2]","par[3]","par[4]","par[5]")])
      names(E4resLNI)=c("BMDL","BMD","BMDU","at","k","ct","dt","is2t")

      E4outLNI <- outLP(parsE4LN, pvec, data.modstan1a$maxD)

      DRM_E4_LN = exp(DRM.E4_LNI(E4resLNI[4:7], data.modstan1a$x, data.modstan1a$q))

      # Covariance between b-d and between BMD-d
      E4covLNI = c(cov(as.matrix(fitstanE4_LN)[,c("b","d")], use="complete.obs")["b","d"],
                   cov(as.matrix(fitstanE4_LN)[,c("BMD","d")], use="complete.obs")["BMD","d"])

      E4corrLNI = c(cor(as.matrix(fitstanE4_LN)[,c("b","d")], use="complete.obs")["b","d"],
                    cor(as.matrix(fitstanE4_LN)[,c("BMD","d")], use="complete.obs")["BMD","d"])

      BMDL=c(BMDL, E4resLNI[1]); BMD=c(BMD, E4resLNI[2]); BMDU=c(BMDU, E4resLNI[3])
      bridgeE4LN = bridgesampling::bridge_sampler(fitstanE4_LN, silent = T)
    }else{E4resLNI=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeE4LN=NA; converged=c(converged, NA); E4covLNI=rep(NA,2); E4corrLNI=rep(NA,2); DRM_E4_LN=rep(NA,length(data.modstan1a$x))
    parsE4LN <- NA
    div_E4_LN <- NA

    E4outLNI <- t(data.frame(
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

    if(prior.weights[10]>0){
      print(10)
      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          fitstanIE4_LN = fun_sampling(stanmodels$mIE4_LNI, data.modstan1a, svF1,
                                       ndraws,nrchains,
                                       nriterations,warmup,
                                       delta,treedepth,seed,pvec)
        }else if(prior.BMD == TRUE){
          fitstanIE4_LN = fun_sampling(stanmodels$mIE4_LNI2, data.modstan1a, svF1,
                                       ndraws,nrchains,
                                       nriterations,warmup,
                                       delta,treedepth,seed,pvec)
        }
      }else if(priordist == "Normal"){
        fitstanIE4_LN = fun_sampling(stanmodels$mIE4_LNI_n, data.modstan1a, svF1,
                                     ndraws,nrchains,
                                     nriterations,warmup,
                                     delta,treedepth,seed,pvec)
      }

      parsIE4LN <- par_extract(fitstanIE4_LN, model_name = "IE4_LN")
      # parsIE4LN[,"BMD"] <- parsIE4LN[,"BMD"]*data.modstan1a$maxD # BMD on original scale

      #diagnostics here
      # posterior_diag(model_stan = fitstanIE4_LN)
      convergence_stat <- convergence_deci(model_stan = fitstanIE4_LN, nrchains = nrchains)
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)

      # divergent transitions
      div_IE4_LN <- sum(sapply(rstan::get_sampler_params(fitstanIE4_LN, inc_warmup = F),
                               function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))

      IE4resLNI=exp(quantile(as.matrix(fitstanIE4_LN)[,"par[2]"],pvec))*data.modstan1a$maxD
      IE4resLNI=c(IE4resLNI,apply(as.matrix(fitstanIE4_LN),2,median)[c("par[1]","par[2]","par[3]","par[4]","par[5]")])
      names(IE4resLNI)=c("BMDL","BMD","BMDU","at","k","ct","dt","is2t")

      IE4outLNI <- outLP(parsIE4LN, pvec, data.modstan1a$maxD)


      DRM_IE4_LN = exp(DRM.IE4_LNI(IE4resLNI[4:7], data.modstan1a$x, data.modstan1a$q))


      IE4covLNI = c(cov(as.matrix(fitstanIE4_LN)[,c("b","d")], use="complete.obs")["b","d"],
                    cov(as.matrix(fitstanIE4_LN)[,c("BMD","d")], use="complete.obs")["BMD","d"])

      IE4corrLNI = c(cor(as.matrix(fitstanIE4_LN)[,c("b","d")], use="complete.obs")["b","d"],
                     cor(as.matrix(fitstanIE4_LN)[,c("BMD","d")], use="complete.obs")["BMD","d"])

      BMDL=c(BMDL, IE4resLNI[1]); BMD=c(BMD, IE4resLNI[2]); BMDU=c(BMDU, IE4resLNI[3])
      bridgeIE4LN = bridgesampling::bridge_sampler(fitstanIE4_LN, silent = T)
    }else{IE4resLNI=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeIE4LN=NA; converged=c(converged, NA); IE4covLNI=rep(NA,2); IE4corrLNI=rep(NA,2); DRM_IE4_LN=rep(NA,length(data.modstan1a$x))
    parsIE4LN <- NA
    div_IE4_LN <- NA

    IE4outLNI <- t(data.frame(
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

    if(prior.weights[11]>0){
      print(11)
      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          fitstanH4_LN = fun_sampling(stanmodels$mH4_LNI, data.modstan1a, svF1,
                                      ndraws,nrchains,
                                      nriterations,warmup,
                                      delta,treedepth,seed,pvec)
        }else if(prior.BMD == TRUE){
          fitstanH4_LN = fun_sampling(stanmodels$mH4_LNI2, data.modstan1a, svF1,
                                      ndraws,nrchains,
                                      nriterations,warmup,
                                      delta,treedepth,seed,pvec)
        }
      }else if(priordist == "Normal"){
        fitstanH4_LN = fun_sampling(stanmodels$mH4_LNI_n, data.modstan1a, svF1,
                                    ndraws,nrchains,
                                    nriterations,warmup,
                                    delta,treedepth,seed,pvec)
      }

      parsH4LN <- par_extract(fitstanH4_LN, model_name = "H4_LN")
      # parsH4LN[,"BMD"] <- parsH4LN[,"BMD"]*data.modstan1a$maxD # BMD on original scale

      #diagnostics here
      # posterior_diag(model_stan = fitstanH4_LN)
      convergence_stat <- convergence_deci(model_stan = fitstanH4_LN, nrchains = nrchains)
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)

      # divergent transitions
      div_H4_LN <- sum(sapply(rstan::get_sampler_params(fitstanH4_LN, inc_warmup = F),
                              function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))

      H4resLNI=exp(quantile(as.matrix(fitstanH4_LN)[,"par[2]"],pvec))*data.modstan1a$maxD
      H4resLNI=c(H4resLNI,apply(as.matrix(fitstanH4_LN),2,median)[c("par[1]","par[2]","par[3]","par[4]","par[5]")])
      names(H4resLNI)=c("BMDL","BMD","BMDU","at","k","ct","dt","is2t")

      H4outLNI <- outLP(parsH4LN, pvec, data.modstan1a$maxD)

      DRM_H4_LN = exp(DRM.H4_LNI(H4resLNI[4:7], data.modstan1a$x, data.modstan1a$q))

      H4covLNI = c(cov(as.matrix(fitstanH4_LN)[,c("b","d")], use="complete.obs")["b","d"],
                   cov(as.matrix(fitstanH4_LN)[,c("BMD","d")], use="complete.obs")["BMD","d"])

      H4corrLNI = c(cor(as.matrix(fitstanH4_LN)[,c("b","d")], use="complete.obs")["b","d"],
                    cor(as.matrix(fitstanH4_LN)[,c("BMD","d")], use="complete.obs")["BMD","d"])

      BMDL=c(BMDL, H4resLNI[1]); BMD=c(BMD, H4resLNI[2]); BMDU=c(BMDU, H4resLNI[3])
      bridgeH4LN = bridgesampling::bridge_sampler(fitstanH4_LN, silent = T)
    }else{H4resLNI=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeH4LN=NA; converged=c(converged, NA); H4covLNI=rep(NA,2); H4corrLNI=rep(NA,2); DRM_H4_LN=rep(NA,length(data.modstan1a$x))
    parsH4LN <- NA
    div_H4_LN <- NA

    H4outLNI <- t(data.frame(
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

    if(prior.weights[12]>0){
      print(12)
      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          fitstanLN4_LN = fun_sampling(stanmodels$mLN4_LNI, data.modstan1aLN, svF1LN,
                                       ndraws,nrchains,
                                       nriterations,warmup,
                                       delta,treedepth,seed,pvec)
        }else if(prior.BMD == TRUE){
          fitstanLN4_LN = fun_sampling(stanmodels$mLN4_LNI2, data.modstan1aLN, svF1LN,
                                       ndraws,nrchains,
                                       nriterations,warmup,
                                       delta,treedepth,seed,pvec)
        }
      }else if(priordist == "Normal"){
        fitstanLN4_LN = fun_sampling(stanmodels$mLN4_LNI_n, data.modstan1aLN, svF1LN,
                                     ndraws,nrchains,
                                     nriterations,warmup,
                                     delta,treedepth,seed,pvec)
      }

      parsLN4LN <- par_extract(fitstanLN4_LN, model_name = "LN4_LN")
      # parsLN4LN[,"BMD"] <- parsLN4LN[,"BMD"]*data.modstan1a$maxD # BMD on original scale

      #diagnostics here
      # posterior_diag(model_stan = fitstanLN4_LN)
      convergence_stat <- convergence_deci(model_stan = fitstanLN4_LN, nrchains = nrchains)
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)

      # divergent transitions
      div_LN4_LN <- sum(sapply(rstan::get_sampler_params(fitstanLN4_LN, inc_warmup = F),
                               function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))

      LN4resLNI=exp(quantile(as.matrix(fitstanLN4_LN)[,"par[2]"],pvec))*data.modstan1a$maxD
      LN4resLNI=c(LN4resLNI,apply(as.matrix(fitstanLN4_LN),2,median)[c("par[1]","par[2]","par[3]","par[4]","par[5]")])
      names(LN4resLNI)=c("BMDL","BMD","BMDU","at","k","ct","dt","is2t")

      LN4outLNI <- outLP(parsLN4LN, pvec, data.modstan1a$maxD)


      DRM_LN4_LN = exp(DRM.LN4_LNI(LN4resLNI[4:7], data.modstan1a$x, data.modstan1a$q))


      LN4covLNI = c(cov(as.matrix(fitstanLN4_LN)[,c("b","d")], use="complete.obs")["b","d"],
                    cov(as.matrix(fitstanLN4_LN)[,c("BMD","d")], use="complete.obs")["BMD","d"])

      LN4corrLNI = c(cor(as.matrix(fitstanLN4_LN)[,c("b","d")], use="complete.obs")["b","d"],
                     cor(as.matrix(fitstanLN4_LN)[,c("BMD","d")], use="complete.obs")["BMD","d"])

      BMDL=c(BMDL, LN4resLNI[1]); BMD=c(BMD, LN4resLNI[2]); BMDU=c(BMDU, LN4resLNI[3])
      bridgeLN4LN = bridgesampling::bridge_sampler(fitstanLN4_LN, silent = T)
    }else{LN4resLNI=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeLN4LN=NA; converged=c(converged, NA); LN4covLNI=rep(NA,2); LN4corrLNI=rep(NA,2); DRM_LN4_LN=rep(NA,length(data.modstan1aLN$x))
    parsLN4LN <- NA
    div_LN4_LN <- NA

    LN4outLNI <- t(data.frame(
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

    if(prior.weights[13]>0){
      print(13)
      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          fitstanG4_LN = fun_sampling(stanmodels$mG4_LNI, data.modstan1bG, svF1,
                                      ndraws,nrchains,
                                      nriterations,warmup,
                                      delta,treedepth,seed,pvec)
        }else if(prior.BMD == TRUE){
          fitstanG4_LN = fun_sampling(stanmodels$mG4_LNI2, data.modstan1bG, svF1,
                                      ndraws,nrchains,
                                      nriterations,warmup,
                                      delta,treedepth,seed,pvec)
        }
      }else if(priordist == "Normal"){
        fitstanG4_LN = fun_sampling(stanmodels$mG4_LNI_n, data.modstan1bG, svF1,
                                    ndraws,nrchains,
                                    nriterations,warmup,
                                    delta,treedepth,seed,pvec)
      }

      parsG4LN <- par_extract(fitstanG4_LN, model_name = "G4_LN")
      # parsG4LN[,"BMD"] <- parsG4LN[,"BMD"]*data.modstan1a$maxD # BMD on original scale


      #diagnostics here
      # posterior_diag(model_stan = fitstanG4_LN)
      convergence_stat <- convergence_deci(model_stan = fitstanG4_LN, nrchains = nrchains)
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)

      # divergent transitions
      div_G4_LN <- sum(sapply(rstan::get_sampler_params(fitstanG4_LN, inc_warmup = F),
                              function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))

      G4resLNI=exp(quantile(as.matrix(fitstanG4_LN)[,"par[2]"],pvec))*data.modstan1a$maxD
      G4resLNI=c(G4resLNI,apply(as.matrix(fitstanG4_LN),2,median)[c("par[1]","par[2]","par[3]","par[4]","par[5]")])
      names(G4resLNI)=c("BMDL","BMD","BMDU","at","k","ct","dt","is2t")

      G4outLNI <- outLP(parsG4LN, pvec, data.modstan1a$maxD)


      DRM_G4_LN = exp(DRM.G4_LNI(G4resLNI[4:7], data.modstan1a$x, data.modstan1a$q))


      G4covLNI = c(cov(as.matrix(fitstanG4_LN)[,c("b","d")], use="complete.obs")["b","d"],
                   cov(as.matrix(fitstanG4_LN)[,c("BMD","d")], use="complete.obs")["BMD","d"])

      G4corrLNI = c(cor(as.matrix(fitstanG4_LN)[,c("b","d")], use="complete.obs")["b","d"],
                    cor(as.matrix(fitstanG4_LN)[,c("BMD","d")], use="complete.obs")["BMD","d"])

      BMDL=c(BMDL, G4resLNI[1]); BMD=c(BMD, G4resLNI[2]); BMDU=c(BMDU, G4resLNI[3])
      bridgeG4LN = bridgesampling::bridge_sampler(fitstanG4_LN, silent = T)
    }else{G4resLNI=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeG4LN=NA; converged=c(converged, NA); G4covLNI=rep(NA,2); G4corrLNI=rep(NA,2); DRM_G4_LN=rep(NA,length(data.modstan1bG$x))
    parsG4LN <- NA
    div_G4_LN <- NA

    G4outLNI <- t(data.frame(
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

    if(prior.weights[14]>0){
      print(14)
      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          fitstanQE4_LN = fun_sampling(stanmodels$mQE4_LNI, data.modstan1bQ, svF1,
                                       ndraws,nrchains,
                                       nriterations,warmup,
                                       delta,treedepth,seed,pvec)
        }else if(prior.BMD == TRUE){
          fitstanQE4_LN = fun_sampling(stanmodels$mQE4_LNI2, data.modstan1bQ, svF1,
                                       ndraws,nrchains,
                                       nriterations,warmup,
                                       delta,treedepth,seed,pvec)
        }
      }else if(priordist == "Normal"){
        fitstanQE4_LN = fun_sampling(stanmodels$mQE4_LNI_n, data.modstan1bQ, svF1,
                                     ndraws,nrchains,
                                     nriterations,warmup,
                                     delta,treedepth,seed,pvec)
      }

      parsQE4LN <- par_extract(fitstanQE4_LN, model_name = "QE4_LN")
      # parsQE4LN[,"BMD"] <- parsQE4LN[,"BMD"]*data.modstan1a$maxD # BMD on original scale

      #diagnostics here
      # posterior_diag(model_stan = fitstanQE4_LN)
      convergence_stat <- convergence_deci(model_stan = fitstanQE4_LN, nrchains = nrchains)
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)

      # divergent transitions
      div_QE4_LN <- sum(sapply(rstan::get_sampler_params(fitstanQE4_LN, inc_warmup = F),
                               function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))

      QE4resLNI=exp(quantile(as.matrix(fitstanQE4_LN)[,"par[2]"],pvec))*data.modstan1a$maxD
      QE4resLNI=c(QE4resLNI,apply(as.matrix(fitstanQE4_LN),2,median)[c("par[1]","par[2]","par[3]","par[4]","par[5]")])
      names(QE4resLNI)=c("BMDL","BMD","BMDU","at","k","ct","dt","is2t")

      QE4outLNI <- outLP(parsQE4LN, pvec, data.modstan1a$maxD)

      DRM_QE4_LN = exp(DRM.QE4_LNI(QE4resLNI[4:7], data.modstan1a$x, data.modstan1a$q))


      QE4covLNI = c(cov(as.matrix(fitstanQE4_LN)[,c("b","d")], use="complete.obs")["b","d"],
                    cov(as.matrix(fitstanQE4_LN)[,c("BMD","d")], use="complete.obs")["BMD","d"])

      QE4corrLNI = c(cor(as.matrix(fitstanQE4_LN)[,c("b","d")], use="complete.obs")["b","d"],
                     cor(as.matrix(fitstanQE4_LN)[,c("BMD","d")], use="complete.obs")["BMD","d"])

      BMDL=c(BMDL, QE4resLNI[1]); BMD=c(BMD, QE4resLNI[2]); BMDU=c(BMDU, QE4resLNI[3])
      bridgeQE4LN = bridgesampling::bridge_sampler(fitstanQE4_LN, silent = T)
    }else{QE4resLNI=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeQE4LN=NA; converged=c(converged, NA); QE4covLNI=rep(NA,2); QE4corrLNI=rep(NA,2); DRM_QE4_LN=rep(NA,length(data.modstan1bQ$x))
    parsQE4LN <- NA
    div_QE4_LN <- NA

    QE4outLNI <- t(data.frame(
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

    if(prior.weights[15]>0){
      print(15)
      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          fitstanP4_LN = fun_sampling(stanmodels$mP4_LNI, data.modstan2P, svF2P,
                                      ndraws,nrchains,
                                      nriterations,warmup,
                                      delta,treedepth,seed,pvec)
        }else if(prior.BMD == TRUE){
          fitstanP4_LN = fun_sampling(stanmodels$mP4_LNI2, data.modstan2P, svF2P,
                                      ndraws,nrchains,
                                      nriterations,warmup,
                                      delta,treedepth,seed,pvec)
        }
      }else if(priordist == "Normal"){
        fitstanP4_LN = fun_sampling(stanmodels$mP4_LNI_n, data.modstan2P, svF2P,
                                    ndraws,nrchains,
                                    nriterations,warmup,
                                    delta,treedepth,seed,pvec)
      }

      parsP4LN <- par_extract(fitstanP4_LN, model_name = "P4_LN")
      # parsP4LN[,"BMD"] <- parsP4LN[,"BMD"]*data.modstan1a$maxD # BMD on original scale

      #diagnostics here
      # posterior_diag(model_stan = fitstanP4_LN)
      convergence_stat <- convergence_deci(model_stan = fitstanP4_LN, nrchains = nrchains)
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)

      # divergent transitions
      div_P4_LN <- sum(sapply(rstan::get_sampler_params(fitstanP4_LN, inc_warmup = F),
                              function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))

      P4resLNI=exp(quantile(as.matrix(fitstanP4_LN)[,"par[2]"],pvec))*data.modstan1a$maxD
      P4resLNI=c(P4resLNI,apply(as.matrix(fitstanP4_LN),2,median)[c("par[1]","par[2]","par[3]","par[4]","par[5]")])
      names(P4resLNI)=c("BMDL","BMD","BMDU","at","k","ct","dt","is2t")

      P4outLNI <- outLP(parsP4LN, pvec, data.modstan1a$maxD)


      DRM_P4_LN = exp(DRM.P4_LNI(P4resLNI[4:7], data.modstan1a$x, data.modstan1a$q))


      P4covLNI = c(cov(as.matrix(fitstanP4_LN)[,c("b","d")], use="complete.obs")["b","d"],
                   cov(as.matrix(fitstanP4_LN)[,c("BMD","d")], use="complete.obs")["BMD","d"])

      P4corrLNI = c(cor(as.matrix(fitstanP4_LN)[,c("b","d")], use="complete.obs")["b","d"],
                    cor(as.matrix(fitstanP4_LN)[,c("BMD","d")], use="complete.obs")["BMD","d"])

      BMDL=c(BMDL, P4resLNI[1]); BMD=c(BMD, P4resLNI[2]); BMDU=c(BMDU, P4resLNI[3])
      bridgeP4LN = bridge_sampler(fitstanP4_LN, silent = T)
    }else{P4resLNI=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeP4LN=NA;
    converged=c(converged, NA); P4covLNI=rep(NA,2); P4corrLNI=rep(NA,2); DRM_P4_LN=rep(NA,length(data.modstan2P$x))
    parsP4LN <- NA
    div_P4_LN <- NA

    P4outLNI <- t(data.frame(
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

    if(prior.weights[16]>0){
      print(16)
      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          fitstanL4_LN = fun_sampling(stanmodels$mL4_LNI, data.modstan2L, svF2L,
                                      ndraws,nrchains,
                                      nriterations,warmup,
                                      delta,treedepth,seed,pvec)
        }else if(prior.BMD == TRUE){
          fitstanL4_LN = fun_sampling(stanmodels$mL4_LNI2, data.modstan2L, svF2L,
                                      ndraws,nrchains,
                                      nriterations,warmup,
                                      delta,treedepth,seed,pvec)
        }
      }else if(priordist == "Normal"){
        fitstanL4_LN = fun_sampling(stanmodels$mL4_LNI_n, data.modstan2L, svF2L,
                                    ndraws,nrchains,
                                    nriterations,warmup,
                                    delta,treedepth,seed,pvec)
      }
      parsL4LN <- par_extract(fitstanL4_LN, model_name = "L4_LN")
      # parsL4LN[,"BMD"] <- parsL4LN[,"BMD"]*data.modstan1a$maxD # BMD on original scale


      #diagnostics here
      # posterior_diag(model_stan = fitstanL4_LN)
      convergence_stat <- convergence_deci(model_stan = fitstanL4_LN, nrchains = nrchains)
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)

      # divergent transitions
      div_L4_LN <- sum(sapply(rstan::get_sampler_params(fitstanL4_LN, inc_warmup = F),
                              function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))

      L4resLNI=exp(quantile(as.matrix(fitstanL4_LN)[,"par[2]"],pvec))*data.modstan1a$maxD
      L4resLNI=c(L4resLNI,apply(as.matrix(fitstanL4_LN),2,median)[c("par[1]","par[2]","par[3]","par[4]","par[5]")])
      names(L4resLNI)=c("BMDL","BMD","BMDU","at","k","ct","dt","is2t")

      L4outLNI <- outLP(parsL4LN, pvec, data.modstan1a$maxD)


      DRM_L4_LN = exp(DRM.L4_LNI(L4resLNI[4:7], data.modstan1a$x, data.modstan1a$q))


      L4covLNI = c(cov(as.matrix(fitstanL4_LN)[,c("b","d")], use="complete.obs")["b","d"],
                   cov(as.matrix(fitstanL4_LN)[,c("BMD","d")], use="complete.obs")["BMD","d"])

      L4corrLNI = c(cor(as.matrix(fitstanL4_LN)[,c("b","d")], use="complete.obs")["b","d"],
                    cor(as.matrix(fitstanL4_LN)[,c("BMD","d")], use="complete.obs")["BMD","d"])

      BMDL=c(BMDL, L4resLNI[1]); BMD=c(BMD, L4resLNI[2]); BMDU=c(BMDU, L4resLNI[3])
      bridgeL4LN = bridgesampling::bridge_sampler(fitstanL4_LN, silent = T)
    }else{L4resLNI=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeL4LN=NA; converged=c(converged, NA); L4covLNI=rep(NA,2); L4corrLNI=rep(NA,2); DRM_L4_LN=rep(NA,length(data.modstan2L$x))
    parsL4LN <- NA
    div_L4_LN <- NA

    L4outLNI <- t(data.frame(
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
    mabmd=(c( # normal
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
    macib=exp(quantile(mabmd,pvec))*data.modstan1a$maxD
    names(macib)=c("BMDL","BMD","BMDU") # original scale

    BMDq_bs = exp(quantile(mabmd, seq(0,1,0.005)))*data.modstan1a$maxD


    mods = c("E4_N","IE4_N","H4_N","LN4_N","G4_N","QE4_N","P4_N","L4_N",
             "E4_LN","IE4_LN","H4_LN","LN4_LN","G4_LN","QE4_LN","P4_LN","L4_LN")
    w.bs = rep(0, 16)
    names(w.bs) = mods
    for(i in mods[mods%in%names(lpwb)]){
      w.bs[i] = lpwb[i]
    }

    ### Model-averaged response per dose level
    dr.MA.bs <- c()
    for(i in 1:length(data.modstan1a$x)){
      dr.MA.bs[i] = weighted.mean(x = c(DRM_E4_N[i],DRM_IE4_N[i],DRM_H4_N[i],DRM_LN4_N[i],DRM_G4_N[i],DRM_QE4_N[i],
                                        DRM_P4_N[i], DRM_L4_N[i] ,
                                        DRM_E4_LN[i], DRM_IE4_LN[i], DRM_H4_LN[i], DRM_LN4_LN[i], DRM_G4_LN[i],
                                        DRM_QE4_LN[i], DRM_P4_LN[i],DRM_L4_LN[i]),
                                  w = w.bs,
                                  na.rm = T)

    }


    ##### Weights & MA if one of the models is divergent --> this model gets weight 0
    if(0 %in% converged){
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
      mabmd.conv=(c( # normal
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
      macib.conv=exp(quantile(mabmd.conv,pvec))*data.modstan1a$maxD
      names(macib.conv)=c("BMDL","BMD","BMDU")

      BMDq_bs_conv = exp(quantile(mabmd.conv, seq(0,1,0.005)))*data.modstan1a$maxD


      mods = c("E4_N","IE4_N","H4_N","LN4_N","G4_N","QE4_N","P4_N","L4_N",
               "E4_LN","IE4_LN","H4_LN","LN4_LN","G4_LN","QE4_LN","P4_LN","L4_LN")
      w.bs.conv = rep(0, 16)
      names(w.bs.conv) = mods
      for(i in mods[mods%in%names(lpwb.conv)]){
        w.bs.conv[i] = lpwb.conv[i]
      }

      ### Model-averaged response per dose level
      dr.MA.bs.conv <- c()
      for(i in 1:length(data.modstan1a$x)){
        dr.MA.bs.conv[i] = weighted.mean(x = c(DRM_E4_N[i],DRM_IE4_N[i],DRM_H4_N[i],DRM_LN4_N[i],DRM_G4_N[i],
                                               DRM_QE4_N[i], DRM_P4_N[i], DRM_L4_N[i] ,
                                               DRM_E4_LN[i], DRM_IE4_LN[i], DRM_H4_LN[i], DRM_LN4_LN[i],
                                               DRM_G4_LN[i],DRM_QE4_LN[i], DRM_P4_LN[i],DRM_L4_LN[i]),
                                         w = w.bs.conv,
                                         na.rm = T)
      }

    }else{
      w.bs.conv = NULL; macib.conv = NULL; BMDq_bs_conv = NULL; dr.MA.bs.conv = NULL
    }


    #----------------------------------
    ### weights based on laplace approximation

    # normal
    data.modstan1a=data.N$data.modstan1a
    data.modstan1bG=data.N$data.modstan1bG
    data.modstan1bQ=data.N$data.modstan1bQ
    data.modstan1aLN=data.N$data.modstan1aLN
    data.modstan2L=data.N$data.modstan2L
    data.modstan2P=data.N$data.modstan2P
    data.modstanE=data.N$data.modstanE
    svF1=data.N$svF1
    svF1LN=data.N$svF1LN
    svF2P=data.N$svF2P
    svF2L=data.N$svF2L

    llN = c() # likelihoods

    # getting the posterior modes and the hessian and the model specific posterior distributions
    if(prior.weights[1]>0){
      print("pw1")
      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          optE4_NI <- fun_optim(stanmodels$mE4_NI, data.modstan1a, svF1, ndraws, 123, pvec)
        }else if(prior.BMD == TRUE){
          optE4_NI <- fun_optim(stanmodels$mE4_NI2, data.modstan1a, svF1, ndraws, 123, pvec)
        }
      }else if(priordist == "Normal"){
        optE4_NI <- fun_optim(stanmodels$mE4_NI_n, data.modstan1a, svF1, ndraws, 123, pvec)
      }

      llE4N=llfE4_NI(optE4_NI$par[1:5],nvec=data.modstan1a$n,
                     dvec=data.modstan1a$x,mvec=data.modstan1a$m,
                     s2vec=data.modstan1a$s2,qval=data.modstan1a$q)


      llN = c(llN, llE4N)
    }else{llN = c(llN,NA)}
    #
    if(prior.weights[2]>0){
      print("pw2")
      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          optIE4_NI <- fun_optim(stanmodels$mIE4_NI, data.modstan1a, svF1, ndraws, 123, pvec)
        }else if(prior.BMD == TRUE){
          optIE4_NI <- fun_optim(stanmodels$mIE4_NI2, data.modstan1a, svF1, ndraws, 123, pvec)
        }
      }else if(priordist == "Normal"){
        optIE4_NI <- fun_optim(stanmodels$mIE4_NI_n, data.modstan1a, svF1, ndraws, 123, pvec)
      }

      llIE4N=llfIE4_NI(optIE4_NI$par[1:5],nvec=data.modstan1a$n,
                       dvec=data.modstan1a$x,mvec=data.modstan1a$m,
                       s2vec=data.modstan1a$s2,qval=data.modstan1a$q)

      llN = c(llN, llIE4N)
    }else{llN = c(llN,NA)}
    #
    if(prior.weights[3]>0){
      print("pw3")
      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          optH4_NI <- fun_optim(stanmodels$mH4_NI, data.modstan1a, svF1, ndraws, 123, pvec)
        }else if(prior.BMD == TRUE){
          optH4_NI <- fun_optim(stanmodels$mH4_NI2, data.modstan1a, svF1, ndraws, 123, pvec)
        }
      }else if(priordist == "Normal"){
        optH4_NI <- fun_optim(stanmodels$mH4_NI_n, data.modstan1a, svF1, ndraws, 123, pvec)
      }

      llH4N=llfH4_NI(optH4_NI$par[1:5],nvec=data.modstan1a$n,
                     dvec=data.modstan1a$x,mvec=data.modstan1a$m,
                     s2vec=data.modstan1a$s2,qval=data.modstan1a$q)

      llN = c(llN, llH4N)
    }else{llN = c(llN,NA)}
    #
    if(prior.weights[4]>0){
      print("pw4")
      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          optLN4_NI <- fun_optim(stanmodels$mLN4_NI, data.modstan1aLN, svF1LN, ndraws, 123, pvec)
        }else if(prior.BMD == TRUE){
          optLN4_NI <- fun_optim(stanmodels$mLN4_NI2, data.modstan1aLN, svF1LN, ndraws, 123, pvec)
        }
      }else if(priordist == "Normal"){
        optLN4_NI <- fun_optim(stanmodels$mLN4_NI_n, data.modstan1aLN, svF1LN, ndraws, 123, pvec)
      }

      llLN4N=llfLN4_NI(optLN4_NI$par[1:5],nvec=data.modstan1aLN$n,
                       dvec=data.modstan1aLN$x,mvec=data.modstan1aLN$m,
                       s2vec=data.modstan1aLN$s2,qval=data.modstan1aLN$q)


      llN = c(llN, llLN4N)
    }else{llN = c(llN,NA)}
    #
    if(prior.weights[5]>0){
      print("pw5")
      # optG4_NI = fun_optim_G4(modN, stanmodels$mG4_NI, stv = svF1, data.modstan1bG, ndraws, 123, pvec)
      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          optG4_NI <- fun_optim(stanmodels$mG4_NI, data.modstan1bG, svF1, ndraws, 123, pvec)
        }else if(prior.BMD == TRUE){
          optG4_NI <- fun_optim(stanmodels$mG4_NI2, data.modstan1bG, svF1, ndraws, 123, pvec)
        }
      }else if(priordist == "Normal"){
        optG4_NI <- fun_optim(stanmodels$mG4_NI_n, data.modstan1bG, svF1, ndraws, 123, pvec)
      }

      llG4N=llfG4_NI(optG4_NI$par[1:5],nvec=data.modstan1bG$n,
                     dvec=data.modstan1bG$x,mvec=data.modstan1bG$m,
                     s2vec=data.modstan1bG$s2,qval=data.modstan1bG$q)


      llN = c(llN, llG4N)
    }else{llN = c(llN,NA)}
    #
    if(prior.weights[6]>0){
      print("pw6")
      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          optQE4_NI <- fun_optim(stanmodels$mQE4_NI, data.modstan1bQ, svF1, ndraws, 123, pvec)
        }else if(prior.BMD == TRUE){
          optQE4_NI <- fun_optim(stanmodels$mQE4_NI2, data.modstan1bQ, svF1, ndraws, 123, pvec)
        }
      }else if(priordist == "Normal"){
        optQE4_NI <- fun_optim(stanmodels$mQE4_NI_n, data.modstan1bQ, svF1, ndraws, 123, pvec)
      }

      llQE4N=llfQE4_NI(optQE4_NI$par[1:5],nvec=data.modstan1bQ$n,
                       dvec=data.modstan1bQ$x,mvec=data.modstan1bQ$m,
                       s2vec=data.modstan1bQ$s2,qval=data.modstan1bQ$q)


      llN = c(llN, llQE4N)
    }else{llN = c(llN,NA)}
    #
    if(prior.weights[7]>0){
      print("pw7")
      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          optP4_NI <- fun_optim(stanmodels$mP4_NI, data.modstan2P, svF2P, ndraws, 123, pvec)
        }else if(prior.BMD == TRUE){
          optP4_NI <- fun_optim(stanmodels$mP4_NI2, data.modstan2P, svF2P, ndraws, 123, pvec)
        }
      }else if(priordist == "Normal"){
        optP4_NI <- fun_optim(stanmodels$mP4_NI_n, data.modstan2P, svF2P, ndraws, 123, pvec)
      }

      llP4N=llfP4_NI(optP4_NI$par[1:5],nvec=data.modstan2P$n,
                     dvec=data.modstan2P$x,mvec=data.modstan2P$m,
                     s2vec=data.modstan2P$s2,qval=data.modstan2P$q)

      llN = c(llN, llP4N)
    }else{llN = c(llN,NA)}
    #
    if(prior.weights[8]>0){
      print("pw8")
      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          optL4_NI <- fun_optim(stanmodels$mL4_NI, data.modstan2L, svF2L, ndraws, 123, pvec)
        }else if(prior.BMD == TRUE){
          optL4_NI <- fun_optim(stanmodels$mL4_NI2, data.modstan2L, svF2L, ndraws, 123, pvec)
        }
      }else if(priordist == "Normal"){
        optL4_NI <- fun_optim(stanmodels$mL4_NI_n, data.modstan2L, svF2L, ndraws, 123, pvec)
      }

      llL4N=llfL4_NI(optL4_NI$par[1:5],nvec=data.modstan2L$n,
                     dvec=data.modstan2L$x,mvec=data.modstan2L$m,
                     s2vec=data.modstan2L$s2,qval=data.modstan2L$q)


      llN = c(llN, llL4N)
    }else{llN = c(llN,NA)}


    # lognormal
    data.modstan1a=data.LN$data.modstan1a
    data.modstan1bG=data.LN$data.modstan1bG
    data.modstan1bQ=data.LN$data.modstan1bQ
    data.modstan1aLN=data.LN$data.modstan1aLN
    data.modstan2L=data.LN$data.modstan2L
    data.modstan2P=data.LN$data.modstan2P
    data.modstanE=data.LN$data.modstanE
    svF1=data.LN$svF1
    svF1LN=data.LN$svF1LN
    svF2P=data.LN$svF2P
    svF2L=data.LN$svF2L

    llLN = c()

    if(prior.weights[9]>0){
      print("pw9")
      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          optE4_LNI <- fun_optim(stanmodels$mE4_LNI, data.modstan1a, svF1, ndraws, 123, pvec)
        }else if(prior.BMD == TRUE){
          optE4_LNI <- fun_optim(stanmodels$mE4_LNI2, data.modstan1a, svF1, ndraws, 123, pvec)
        }
      }else if(priordist == "Normal"){
        optE4_LNI <- fun_optim(stanmodels$mE4_LNI_n, data.modstan1a, svF1, ndraws, 123, pvec)
      }

      llE4LN=llfE4_LNI(optE4_LNI$par[1:5],nvec=data.modstan1a$n,
                       dvec=data.modstan1a$x,mvec=data.modstan1a$m,
                       s2vec=data.modstan1a$s2,qval=data.modstan1a$q)


      llLN = c(llLN, llE4LN)
    }else{llLN = c(llLN,NA)}
    #
    if(prior.weights[10]>0){
      print("pw10")
      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          optIE4_LNI <- fun_optim(stanmodels$mIE4_LNI, data.modstan1a, svF1, ndraws, 123, pvec)
        }else if(prior.BMD == TRUE){
          optIE4_LNI <- fun_optim(stanmodels$mIE4_LNI2, data.modstan1a, svF1, ndraws, 123, pvec)
        }
      }else if(priordist == "Normal"){
        optIE4_LNI <- fun_optim(stanmodels$mIE4_LNI_n, data.modstan1a, svF1, ndraws, 123, pvec)
      }

      llIE4LN=llfIE4_LNI(optIE4_LNI$par[1:5],nvec=data.modstan1a$n,
                         dvec=data.modstan1a$x,mvec=data.modstan1a$m,
                         s2vec=data.modstan1a$s2,qval=data.modstan1a$q)


      llLN = c(llLN, llIE4LN)
    }else{llLN = c(llLN,NA)}
    #
    if(prior.weights[11]>0){
      print("pw11")
      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          optH4_LNI <- fun_optim(stanmodels$mH4_LNI, data.modstan1a, svF1, ndraws, 123, pvec)
        }else if(prior.BMD == TRUE){
          optH4_LNI <- fun_optim(stanmodels$mH4_LNI2, data.modstan1a, svF1, ndraws, 123, pvec)
        }
      }else if(priordist == "Normal"){
        optH4_LNI <- fun_optim(stanmodels$mH4_LNI_n, data.modstan1a, svF1, ndraws, 123, pvec)
      }

      llH4LN=llfH4_LNI(optH4_LNI$par[1:5],nvec=data.modstan1a$n,
                       dvec=data.modstan1a$x,mvec=data.modstan1a$m,
                       s2vec=data.modstan1a$s2,qval=data.modstan1a$q)

      llLN = c(llLN, llH4LN)
    }else{llLN = c(llLN,NA)}
    #
    if(prior.weights[12]>0){
      print("pw12")
      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          optLN4_LNI <- fun_optim(stanmodels$mLN4_LNI, data.modstan1aLN, svF1LN, ndraws, 123, pvec)
        }else if(prior.BMD == TRUE){
          optLN4_LNI <- fun_optim(stanmodels$mLN4_LNI2, data.modstan1aLN, svF1LN, ndraws, 123, pvec)
        }
      }else if(priordist == "Normal"){
        optLN4_LNI <- fun_optim(stanmodels$mLN4_LNI_n, data.modstan1aLN, svF1LN, ndraws, 123, pvec)
      }

      llLN4LN=llfLN4_LNI(optLN4_LNI$par[1:5],nvec=data.modstan1aLN$n,
                         dvec=data.modstan1aLN$x,mvec=data.modstan1aLN$m,
                         s2vec=data.modstan1aLN$s2,qval=data.modstan1aLN$q)

      llLN = c(llLN, llLN4LN)
    }else{llLN = c(llLN,NA)}
    #
    if(prior.weights[13]>0){

      print("pw13")
      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          optG4_LNI <- fun_optim(stanmodels$mG4_LNI, data.modstan1bG, svF1, ndraws, 123, pvec)
        }else if(prior.BMD == TRUE){
          optG4_LNI <- fun_optim(stanmodels$mG4_LNI2, data.modstan1bG, svF1, ndraws, 123, pvec)
        }
      }else if(priordist == "Normal"){
        optG4_LNI <- fun_optim(stanmodels$mG4_LNI_n, data.modstan1bG, svF1, ndraws, 123, pvec)
      }

      llG4LN=llfG4_LNI(optG4_LNI$par[1:5],nvec=data.modstan1bG$n,
                       dvec=data.modstan1bG$x,mvec=data.modstan1bG$m,
                       s2vec=data.modstan1bG$s2,qval=data.modstan1bG$q)

      llLN = c(llLN, llG4LN)
    }else{llLN = c(llLN,NA)}
    #
    if(prior.weights[14]>0){
      print("pw14")
      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          optQE4_LNI <- fun_optim(stanmodels$mQE4_LNI, data.modstan1bQ, svF1, ndraws, 123, pvec)
        }else if(prior.BMD == TRUE){
          optQE4_LNI <- fun_optim(stanmodels$mQE4_LNI2, data.modstan1bQ, svF1, ndraws, 123, pvec)
        }
      }else if(priordist == "Normal"){
        optQE4_LNI <- fun_optim(stanmodels$mQE4_LNI_n, data.modstan1bQ, svF1, ndraws, 123, pvec)
      }

      llQE4LN=llfQE4_LNI(optQE4_LNI$par[1:5],nvec=data.modstan1bQ$n,
                         dvec=data.modstan1bQ$x,mvec=data.modstan1bQ$m,
                         s2vec=data.modstan1bQ$s2,qval=data.modstan1bQ$q)

      llLN = c(llLN, llQE4LN)
    }else{llLN = c(llLN,NA)}
    #
    if(prior.weights[15]>0){
      print("pw15")
      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          optP4_LNI <- fun_optim(stanmodels$mP4_LNI, data.modstan2P, svF2P, ndraws, 123, pvec)
        }else if(prior.BMD == TRUE){
          optP4_LNI <- fun_optim(stanmodels$mP4_LNI2, data.modstan2P, svF2P, ndraws, 123, pvec)
        }
      }else if(priordist == "Normal"){
        optP4_LNI <- fun_optim(stanmodels$mP4_LNI_n, data.modstan2P, svF2P, ndraws, 123, pvec)
      }

      llP4LN=llfP4_LNI(optP4_LNI$par[1:5],nvec=data.modstan2P$n,
                       dvec=data.modstan2P$x,mvec=data.modstan2P$m,
                       s2vec=data.modstan2P$s2,qval=data.modstan2P$q)

      llLN = c(llLN, llP4LN)
    }else{llLN = c(llLN,NA)}
    #
    if(prior.weights[16]>0){
      print("pw16")
      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          optL4_LNI <- fun_optim(stanmodels$mL4_LNI, data.modstan2L, svF2L, ndraws, 123, pvec)
        }else if(prior.BMD == TRUE){
          optL4_LNI <- fun_optim(stanmodels$mL4_LNI2, data.modstan2L, svF2L, ndraws, 123, pvec)
        }
      }else if(priordist == "Normal"){
        optL4_LNI <- fun_optim(stanmodels$mL4_LNI_n, data.modstan2L, svF2L, ndraws, 123, pvec)
      }

      llL4LN=llfL4_LNI(optL4_LNI$par[1:5],nvec=data.modstan2L$n,
                       dvec=data.modstan2L$x,mvec=data.modstan2L$m,
                       s2vec=data.modstan2L$s2,qval=data.modstan2L$q)

      llLN = c(llLN, llL4LN)
    }else{llLN = c(llLN,NA)}

    minll = min(llN,llLN,na.rm=T)

    # the weights

    fun.w <- function(DIH, ll, min.ll, opt, mu, sig, prior, priorbmd, lb, ub, s1, s2, s3){
      if(prior == "PERT"){
        if(priorbmd == FALSE){
          w1 <- (2*pi)^(2.5)*sqrt(DIH)*exp(ll-min.ll)*
            mvtnorm::dmvnorm(opt$par[c(4,5)],mean=mu[c(4,5)],
                             sigma=sig[c(4,5),c(4,5)])*
            dexp(-opt$par[2],rate=mu[2])*
            mc2d::dpert(opt$par[1], min = lb[1], max = ub[1],
                        mode = mu[1], shape = s1)*
            mc2d::dpert(opt$par[3], min = lb[3], max = ub[3],
                        mode = mu[3], shape = s3)
        }else if(priorbmd == TRUE){
          w1 <- (2*pi)^(2.5)*sqrt(DIH)*exp(ll-min.ll)*
            mvtnorm::dmvnorm(opt$par[c(4,5)],mean=mu[c(4,5)],
                             sigma=sig[c(4,5),c(4,5)])*
            mc2d::dpert(opt$par[1], min = lb[1], max = ub[1],
                        mode = mu[1], shape = s1)*
            mc2d::dpert(opt$par[2], min = lb[2], max = ub[2],
                        mode = mu[2], shape = s2)*
            mc2d::dpert(opt$par[3], min = lb[3], max = ub[3],
                        mode = mu[3], shape = s3)
        }
      }else if(prior == "Normal"){
        w1 <- (2*pi)^(2.5)*sqrt(DIH)*exp(ll-min.ll)*
          mvtnorm::dmvnorm(opt$par[c(1,2,4,5)],mean=mu[c(1,2,4,5)],
                           sigma=sig[c(1,2,4,5),c(1,2,4,5)])*
          dexp(-opt$par[2],rate=mu[2])
      }
      return(w1)
    }

    w=c()

    # normal

    data.modstan1a=data.N$data.modstan1a
    data.modstan1bG=data.N$data.modstan1bG
    data.modstan1bQ=data.N$data.modstan1bQ
    data.modstan1aLN=data.N$data.modstan1aLN
    data.modstan2L=data.N$data.modstan2L
    data.modstan2P=data.N$data.modstan2P
    data.modstanE=data.N$data.modstanE
    svF1=data.N$svF1
    svF1LN = data.N$svF1LN
    svF2P=data.N$svF2P
    svF2L=data.N$svF2L

    if(prior.weights[1]>0){
      DIHE4h=det(-solve(optE4_NI$hessian))
      DIHE4=ifelse(DIHE4h<0,0,DIHE4h)
      ## Aproximation of marginal (i.e. integrated) likelihood (= 'model evidence')
      # if(priordist == "PERT"){
      #   w=c(w,(2*pi)^(2.5)*sqrt(DIHE4)*exp(llE4N-minll)*
      #         mvtnorm::dmvnorm(optE4_NI$par[c(4,5)],mean=data.modstan1a$priormu[c(4,5)],
      #                          sigma=data.modstan1a$priorSigma[c(4,5),c(4,5)])*
      #         dexp(-optE4_NI$par[2],rate=data.modstan1a$priormu[2])*
      #         mc2d::dpert(optE4_NI$par[1], min = data.modstan1a$priorlb[1], max = data.modstan1a$priorub[1],
      #                     mode = data.modstan1a$priormu[1], shape = 4)*
      #         mc2d::dpert(optE4_NI$par[3], min = data.modstan1a$priorlb[3], max = data.modstan1a$priorub[3],
      #                     mode = data.modstan1a$priormu[3], shape = 4)
      #   )
      # }else if(priordist == "Normal"){
      #   w=c(w,(2*pi)^(2.5)*sqrt(DIHE4)*exp(llE4N-minll)*
      #         mvtnorm::dmvnorm(optE4_NI$par[c(1,2,4,5)],mean=data.modstan1a$priormu[c(1,2,4,5)],
      #                          sigma=data.modstan1a$priorSigma[c(1,2,4,5),c(1,2,4,5)])*
      #         dexp(-optE4_NI$par[2],rate=data.modstan1a$priormu[2])
      #   )
      # }
      w = c(w, fun.w(DIHE4, llE4N, minll, optE4_NI, data.modstan1a$priormu, data.modstan1a$priorSigma, priordist, prior.BMD,
                     data.modstan1a$priorlb, data.modstan1a$priorub, data.modstan1a$shape.a, data.modstan1a$shape.BMD, data.modstan1a$shape.c))

    }else{w=c(w,0)}

    if(prior.weights[2]>0){
      DIHIE4h=det(-solve(optIE4_NI$hessian))
      DIHIE4=ifelse(DIHIE4h<0,0,DIHIE4h)
      w = c(w, fun.w(DIHIE4, llIE4N, minll, optIE4_NI, data.modstan1a$priormu, data.modstan1a$priorSigma, priordist, prior.BMD,
                     data.modstan1a$priorlb, data.modstan1a$priorub, data.modstan1a$shape.a, data.modstan1a$shape.BMD, data.modstan1a$shape.c))
    }else{w=c(w,0)}

    if(prior.weights[3]>0){
      DIHH4h=det(-solve(optH4_NI$hessian))
      DIHH4=ifelse(DIHH4h<0,0,DIHH4h)
      w = c(w, fun.w(DIHH4, llH4N, minll, optH4_NI, data.modstan1a$priormu, data.modstan1a$priorSigma, priordist, prior.BMD,
                     data.modstan1a$priorlb, data.modstan1a$priorub, data.modstan1a$shape.a, data.modstan1a$shape.BMD, data.modstan1a$shape.c))
    }else{w=c(w,0)}

    if(prior.weights[4]>0){
      DIHLN4h=det(-solve(optLN4_NI$hessian))
      DIHLN4=ifelse(DIHLN4h<0,0,DIHLN4h)
      w = c(w, fun.w(DIHLN4, llLN4N, minll, optLN4_NI, data.modstan1aLN$priormu, data.modstan1aLN$priorSigma, priordist, prior.BMD,
                     data.modstan1aLN$priorlb, data.modstan1aLN$priorub, data.modstan1a$shape.a, data.modstan1a$shape.BMD, data.modstan1a$shape.c))
    }else{w=c(w,0)}

    if(prior.weights[5]>0){
      DIHG4h=det(-solve(optG4_NI$hessian))
      DIHG4=ifelse(DIHG4h<0,0,DIHG4h)
      w = c(w, fun.w(DIHG4, llG4N, minll, optG4_NI, data.modstan1bG$priormu, data.modstan1bG$priorSigma, priordist, prior.BMD,
                     data.modstan1bG$priorlb, data.modstan1bG$priorub, data.modstan1a$shape.a, data.modstan1a$shape.BMD, data.modstan1a$shape.c))
    }else{w=c(w,0)}

    if(prior.weights[6]>0){
      DIHQE4h=det(-solve(optQE4_NI$hessian))
      DIHQE4=ifelse(DIHQE4h<0,0,DIHQE4h)
      w = c(w, fun.w(DIHQE4, llQE4N, minll, optQE4_NI, data.modstan1bQ$priormu, data.modstan1bQ$priorSigma, priordist, prior.BMD,
                     data.modstan1bQ$priorlb, data.modstan1bQ$priorub, data.modstan1a$shape.a, data.modstan1a$shape.BMD, data.modstan1a$shape.c))
    }else{w=c(w,0)}

    if(prior.weights[7]>0){
      DIHP4h=det(-solve(optP4_NI$hessian))
      DIHP4=ifelse(DIHP4h<0,0,DIHP4h)
      w = c(w, fun.w(DIHP4, llP4N, minll, optP4_NI, data.modstan2P$priormu, data.modstan2P$priorSigma, priordist, prior.BMD,
                     data.modstan2P$priorlb, data.modstan2P$priorub, data.modstan1a$shape.a, data.modstan1a$shape.BMD, data.modstan1a$shape.c))
    }else{w=c(w,0)}

    if(prior.weights[8]>0){
      DIHL4h=det(-solve(optL4_NI$hessian))
      DIHL4=ifelse(DIHL4h<0,0,DIHL4h)
      w = c(w, fun.w(DIHL4, llL4N, minll, optL4_NI, data.modstan2L$priormu, data.modstan2L$priorSigma, priordist, prior.BMD,
                     data.modstan2L$priorlb, data.modstan2L$priorub, data.modstan1a$shape.a, data.modstan1a$shape.BMD, data.modstan1a$shape.c))
    }else{w=c(w,0)}


    # lognormal

    data.modstan1a=data.LN$data.modstan1a
    data.modstan1bG=data.LN$data.modstan1bG
    data.modstan1bQ=data.LN$data.modstan1bQ
    data.modstan1aLN=data.LN$data.modstan1aLN
    data.modstan2L=data.LN$data.modstan2L
    data.modstan2P=data.LN$data.modstan2P
    data.modstanE=data.LN$data.modstanE
    svF1=data.LN$svF1
    svF1LN=data.LN$svF1LN
    svF2P=data.LN$svF2P
    svF2L=data.LN$svF2L

    if(prior.weights[9]>0){
      DIHE4h=det(-solve(optE4_LNI$hessian))
      DIHE4=ifelse(DIHE4h<0,0,DIHE4h)
      w = c(w, fun.w(DIHE4, llE4LN, minll, optE4_LNI, data.modstan1a$priormu, data.modstan1a$priorSigma, priordist, prior.BMD,
                     data.modstan1a$priorlb, data.modstan1a$priorub, data.modstan1a$shape.a, data.modstan1a$shape.BMD, data.modstan1a$shape.c))
    }else{w=c(w,0)}

    if(prior.weights[10]>0){
      DIHIE4h=det(-solve(optIE4_LNI$hessian))
      DIHIE4=ifelse(DIHIE4h<0,0,DIHIE4h)
      w = c(w, fun.w(DIHIE4, llIE4LN, minll, optIE4_LNI, data.modstan1a$priormu, data.modstan1a$priorSigma, priordist, prior.BMD,
                     data.modstan1a$priorlb, data.modstan1a$priorub, data.modstan1a$shape.a, data.modstan1a$shape.BMD, data.modstan1a$shape.c))
    }else{w=c(w,0)}

    if(prior.weights[11]>0){
      DIHH4h=det(-solve(optH4_LNI$hessian))
      DIHH4=ifelse(DIHH4h<0,0,DIHH4h)
      w = c(w, fun.w(DIHH4, llH4LN, minll, optH4_LNI, data.modstan1a$priormu, data.modstan1a$priorSigma, priordist, prior.BMD,
                     data.modstan1a$priorlb, data.modstan1a$priorub, data.modstan1a$shape.a, data.modstan1a$shape.BMD, data.modstan1a$shape.c))
    }else{w=c(w,0)}

    if(prior.weights[12]>0){
      DIHLN4h=det(-solve(optLN4_LNI$hessian))
      DIHLN4=ifelse(DIHLN4h<0,0,DIHLN4h)
      w = c(w, fun.w(DIHLN4, llLN4LN, minll, optLN4_LNI, data.modstan1aLN$priormu, data.modstan1aLN$priorSigma, priordist, prior.BMD,
                     data.modstan1aLN$priorlb, data.modstan1aLN$priorub, data.modstan1a$shape.a, data.modstan1a$shape.BMD, data.modstan1a$shape.c))
    }else{w=c(w,0)}

    if(prior.weights[13]>0){
      DIHG4h=det(-solve(optG4_LNI$hessian))
      DIHG4=ifelse(DIHG4h<0,0,DIHG4h)
      w = c(w, fun.w(DIHG4, llG4LN, minll, optG4_LNI, data.modstan1bG$priormu, data.modstan1bG$priorSigma, priordist, prior.BMD,
                     data.modstan1bG$priorlb, data.modstan1bG$priorub, data.modstan1a$shape.a, data.modstan1a$shape.BMD, data.modstan1a$shape.c))
    }else{w=c(w,0)}

    if(prior.weights[14]>0){
      DIHQE4h=det(-solve(optQE4_LNI$hessian))
      DIHQE4=ifelse(DIHQE4h<0,0,DIHQE4h)
      w = c(w, fun.w(DIHQE4, llQE4LN, minll, optQE4_LNI, data.modstan1bQ$priormu, data.modstan1bQ$priorSigma, priordist, prior.BMD,
                     data.modstan1bQ$priorlb, data.modstan1bQ$priorub, data.modstan1a$shape.a, data.modstan1a$shape.BMD, data.modstan1a$shape.c))
    }else{w=c(w,0)}

    if(prior.weights[15]>0){
      DIHP4h=det(-solve(optP4_LNI$hessian))
      DIHP4=ifelse(DIHP4h<0,0,DIHP4h)
      w = c(w, fun.w(DIHP4, llP4LN, minll, optP4_LNI, data.modstan2P$priormu, data.modstan2P$priorSigma, priordist, prior.BMD,
                     data.modstan2P$priorlb, data.modstan2P$priorub, data.modstan1a$shape.a, data.modstan1a$shape.BMD, data.modstan1a$shape.c))
    }else{w=c(w,0)}

    if(prior.weights[16]>0){
      DIHL4h=det(-solve(optL4_LNI$hessian))
      DIHL4=ifelse(DIHL4h<0,0,DIHL4h)
      w = c(w, fun.w(DIHL4, llL4LN, minll, optL4_LNI, data.modstan2L$priormu, data.modstan2L$priorSigma, priordist, prior.BMD,
                     data.modstan2L$priorlb, data.modstan2L$priorub, data.modstan1a$shape.a, data.modstan1a$shape.BMD, data.modstan1a$shape.c))
    }else{w=c(w,0)}

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
    macilp=exp(quantile(mabmd,pvec))*data.modstan1a$maxD
    names(macilp)=c("BMDL","BMD","BMDU") # on original scale

    BMDq_ls = exp(quantile(mabmd, seq(0,1,0.005)))*data.modstan1a$maxD


    ### Model-averaged response per dose level
    dr.MA.ls <- c()
    for(i in 1:length(data.modstan1a$x)){
      dr.MA.ls[i] = weighted.mean(x = c(DRM_E4_N[i],DRM_IE4_N[i],DRM_H4_N[i],DRM_LN4_N[i],DRM_G4_N[i],DRM_QE4_N[i], DRM_P4_N[i],
                                        DRM_L4_N[i] ,
                                        DRM_E4_LN[i], DRM_IE4_LN[i], DRM_H4_LN[i], DRM_LN4_LN[i], DRM_G4_LN[i],DRM_QE4_LN[i],
                                        DRM_P4_LN[i],DRM_L4_LN[i]),
                                  w = lpwlp,
                                  na.rm = T)
    }


    ##### Weights & MA if one of the models is divergent --> this model gets weight 0
    if(0 %in% converged){

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
      macilp.conv=exp(quantile(mabmd.conv,pvec))*data.modstan1a$maxD
      names(macilp.conv)=c("BMDL","BMD","BMDU")

      BMDq_ls_conv = exp(quantile(mabmd.conv, seq(0,1,0.005)))*data.modstan1a$maxD


      ### Model-averaged response per dose level
      dr.MA.ls.conv <- c()
      for(i in 1:length(data.modstan1a$x)){
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

    if(plot == TRUE){

      data.modstan1a = data.N$data.modstan1a
      data.modstan1bG = data.N$data.modstan1bG
      data.modstanE = data.N$data.modstanE
      data.modstan1bQ = data.N$data.modstan1bQ
      data.modstan1aLN = data.N$data.modstan1aLN
      data.modstan2P = data.N$data.modstan2P
      data.modstan2L = data.N$data.modstan2L
      svF1 = data.N$svF1
      svF1 = data.N$svF1LN
      svF2P = data.N$svF2P
      svF2L = data.N$svF2L

      N = data.N$data.modstan1a$N
      plot2 = function(){
        sd = sqrt(data.modstan1a$s2)
        gmean.a = log(NtoLN(data.modstan1a$m,sd))[1:N]
        gsd.a = log(NtoLN(data.modstan1a$m,sd))[(N+1):(2*N)]
        x = c(log10(data.modstan1a$x[2]/4),log10(data.modstan1a$x)[2:data.modstan1a$N])
        m = log10(exp(gmean.a))
        s = log10(exp(gsd.a))
        bmdl = log10(BMDL/data.modstan1a$maxD)
        x[1] = ifelse(min(bmdl,na.rm=T)<x[1], log10(min(BMDL/data.modstan1a$maxD/2, na.rm=T)), log10(data.modstan1a$x[2]/4))
        g = choose(length(x),2)
        alph = 0.05/(2*g)
        cp = qnorm(1-alph)
        n = data.modstan1a$n
        ci = cp*s/sqrt(n)
        ci2 = 1.96*s/sqrt(n)

        plot(x, m, xlab="log10-Dose", ylab="log10-Response", ylim=c(min(m-2.2*ci),max(m+2.2*ci)),
             xlim=c(ifelse(min(bmdl,na.rm=T)<min(x),min(bmdl,na.rm=T),min(x)),abs(min(x))))#, ylim=c(min(mean.a)-2*sd.a, max(mean.a)+2*sd.a))
        dgr=seq(min(x), abs(min(x)),by=0.01)
        # plot model specific fits
        # normal distribution
        if(prior.weights[1]>0){
          parE4 = E4resNI[4:7]
          DRME4 = log10((DRM.E4_NI(par=parE4,x=10^(dgr[dgr>(x[2]-((x[2]-x[1])/2))]),q=data.modstan1a$q)))
          lines(dgr[dgr>(x[2]-((x[2]-x[1])/2))],DRME4,lwd=1, lty=1, col=1)
          segments(x0=dgr[1], y0=log10(exp(parE4[1])), x1=max(dgr[dgr<=(x[2]-((x[2]-x[1])/2))]), y1=DRME4[1],lty=3,col=1)
        }
        if(prior.weights[2]>0){
          parIE4 = IE4resNI[4:7]
          DRMIE4 = log10((DRM.IE4_NI(par=parIE4,x=10^(dgr[dgr>(x[2]-((x[2]-x[1])/2))]),q=data.modstan1a$q)))
          lines(dgr[dgr>(x[2]-((x[2]-x[1])/2))],DRMIE4,lwd=1, lty=1, col=2)
          segments(x0=dgr[1], y0=log10(exp(parIE4[1])), x1=max(dgr[dgr<=(x[2]-((x[2]-x[1])/2))]), y1=DRMIE4[1],lty=3,col=2)
        }
        if(prior.weights[3]>0){
          parH4 = H4resNI[4:7]
          DRMH4 = log10((DRM.H4_NI(par=parH4,x=10^(dgr[dgr>(x[2]-((x[2]-x[1])/2))]),q=data.modstan1a$q)))
          lines(dgr[dgr>(x[2]-((x[2]-x[1])/2))],DRMH4,lwd=1, lty=1, col=3)
          segments(x0=dgr[1], y0=log10(exp(parH4[1])), x1=max(dgr[dgr<=(x[2]-((x[2]-x[1])/2))]), y1=DRMH4[1],lty=3,col=3)
        }
        if(prior.weights[4]>0){
          parLN4 = LN4resNI[4:7]
          DRMLN4 = log10((DRM.LN4_NI(par=parLN4,x=10^(dgr[dgr>(x[2]-((x[2]-x[1])/2))]),q=data.modstan1a$q)))
          lines(dgr[dgr>(x[2]-((x[2]-x[1])/2))],DRMLN4,lwd=1, lty=1, col=4)
          segments(x0=dgr[1], y0=log10(exp(parLN4[1])), x1=max(dgr[dgr<=(x[2]-((x[2]-x[1])/2))]), y1=DRMLN4[1],lty=3,col=4)
        }
        if(prior.weights[5]>0){
          parG4 = G4resNI[4:7]
          DRMG4 = log10((DRM.G4_NI(par=parG4,x=10^(dgr[dgr>(x[2]-((x[2]-x[1])/2))]),q=data.modstan1a$q)))
          lines(dgr[dgr>(x[2]-((x[2]-x[1])/2))],DRMG4,lwd=1, lty=1, col=5)
          segments(x0=dgr[1], y0=log10(exp(parG4[1])), x1=max(dgr[dgr<=(x[2]-((x[2]-x[1])/2))]), y1=DRMG4[1],lty=3,col=5)
        }
        if(prior.weights[6]>0){
          parQE4 = QE4resNI[4:7]
          DRMQE4 = log10((DRM.QE4_NI(par=parQE4,x=10^(dgr[dgr>(x[2]-((x[2]-x[1])/2))]),q=data.modstan1a$q)))
          lines(dgr[dgr>(x[2]-((x[2]-x[1])/2))],DRMQE4,lwd=1, lty=1, col=6)
          segments(x0=dgr[1], y0=log10(exp(parQE4[1])), x1=max(dgr[dgr<=(x[2]-((x[2]-x[1])/2))]), y1=DRMQE4[1],lty=3,col=6)
        }
        if(prior.weights[7]>0){
          parP4 = P4resNI[4:7]
          DRMP4 = log10((DRM.P4_NI(par=parP4,x=10^(dgr[dgr>(x[2]-((x[2]-x[1])/2))]),q=data.modstan1a$q)))
          lines(dgr[dgr>(x[2]-((x[2]-x[1])/2))],DRMP4,lwd=1, lty=1, col=7)
          segments(x0=dgr[1], y0=log10(exp(parP4[1])*pnorm(qnorm(1/(1+q))-exp(parP4[3]))), x1=max(dgr[dgr<=(x[2]-((x[2]-x[1])/2))]), y1=DRMP4[1],lty=3,col=7)
        }
        if(prior.weights[8]>0){
          parL4 = L4resNI[4:7]
          DRML4 = log10((DRM.L4_NI(par=parL4,x=10^(dgr[dgr>(x[2]-((x[2]-x[1])/2))]),q=data.modstan1a$q)))
          lines(dgr[dgr>(x[2]-((x[2]-x[1])/2))],DRML4,lwd=1, lty=1, col=8)
          segments(x0=dgr[1], y0=log10(exp(parL4[1])*expit(logit(1/(1+q))-exp(parL4[3]))), x1=max(dgr[dgr<=(x[2]-((x[2]-x[1])/2))]), y1=DRML4[1],lty=3,col=8)
        }

        # lognormal distribution
        if(prior.weights[9]>0){
          parE4 = E4resLNI[4:7]
          a=exp(parE4[1])
          DRME4 = log10(exp(DRM.E4_LNI(par=parE4,x=10^(dgr[dgr>(x[2]-((x[2]-x[1])/2))]),q=data.modstan1a$q)))
          if(data.LN$data.modstan1a$shift==T) {DRME4 = log10(exp(DRM.E4_LNI(par=parE4,x=10^(dgr[dgr>(x[2]-((x[2]-x[1])/2))]),q=data.modstan1a$q))) + log10(exp(1.5*min(data.LN$data.modstan1a$m.org)))
          a = exp(parE4[1]) + 1.5*min(data.LN$data.modstan1a$m.org)}
          lines(dgr[dgr>(x[2]-((x[2]-x[1])/2))],DRME4,lwd=1, lty=2, col=1)
          segments(x0=dgr[1], y0=log10(exp(a)), x1=max(dgr[dgr<=(x[2]-((x[2]-x[1])/2))]), y1=min(DRME4),lty=3,col=1)
        }
        if(prior.weights[10]>0){
          parIE4 = IE4resLNI[4:7]
          a=exp(parIE4[1])
          DRMIE4 = log10(exp(DRM.IE4_LNI(par=parIE4,x=10^(dgr[dgr>(x[2]-((x[2]-x[1])/2))]),q=data.modstan1a$q)))
          if(data.LN$data.modstan1a$shift==T) {DRMIE4 = log10(exp(DRM.IE4_LNI(par=parIE4,x=10^(dgr[dgr>(x[2]-((x[2]-x[1])/2))]),q=data.modstan1a$q))) + log10(exp(1.5*min(data.LN$data.modstan1a$m.org)))
          a = exp(parIE4[1]) + 1.5*min(data.LN$data.modstan1a$m.org)}
          lines(dgr[dgr>(x[2]-((x[2]-x[1])/2))],DRMIE4,lwd=1, lty=2, col=2)
          segments(x0=dgr[1], y0=log10(exp(a)), x1=max(dgr[dgr<=(x[2]-((x[2]-x[1])/2))]), y1=min(DRMIE4),lty=3,col=2)
        }
        if(prior.weights[11]>0){
          parH4 = H4resLNI[4:7]
          a=exp(parH4[1])
          DRMH4 = log10(exp(DRM.H4_LNI(par=parH4,x=10^(dgr[dgr>(x[2]-((x[2]-x[1])/2))]),q=data.modstan1a$q)))
          if(data.LN$data.modstan1a$shift==T) {DRMH4 = log10(exp(DRM.H4_LNI(par=parH4,x=10^(dgr[dgr>(x[2]-((x[2]-x[1])/2))]),q=data.modstan1a$q))) + log10(exp(1.5*min(data.LN$data.modstan1a$m.org)))
          a = exp(parH4[1]) + 1.5*min(data.LN$data.modstan1a$m.org)}
          lines(dgr[dgr>(x[2]-((x[2]-x[1])/2))],DRMH4,lwd=1, lty=2, col=3)
          segments(x0=dgr[1], y0=log10(exp(a)), x1=max(dgr[dgr<=(x[2]-((x[2]-x[1])/2))]), y1=min(DRMH4),lty=3,col=3)
        }
        if(prior.weights[12]>0){
          parLN4 = LN4resLNI[4:7]
          a=exp(parLN4[1])
          DRMLN4 = log10(exp(DRM.LN4_LNI(par=parLN4,x=10^(dgr[dgr>(x[2]-((x[2]-x[1])/2))]),q=data.modstan1aLN$q)))
          if(data.LN$data.modstan1aLN$shift==T) {DRMLN4 = log10(exp(DRM.LN4_LNI(par=parLN4,x=10^(dgr[dgr>(x[2]-((x[2]-x[1])/2))]),q=data.modstan1aLN$q))) + log10(exp(1.5*min(data.LN$data.modstan1aLN$m.org)))
          a = exp(parLN4[1]) + 1.5*min(data.LN$data.modstan1aLN$m.org)}
          lines(dgr[dgr>(x[2]-((x[2]-x[1])/2))],DRMLN4,lwd=1, lty=2, col=4)
          segments(x0=dgr[1], y0=log10(exp(a)), x1=max(dgr[dgr<=(x[2]-((x[2]-x[1])/2))]), y1=min(DRMLN4),lty=3,col=4)
        }
        if(prior.weights[13]>0){
          parG4 = G4resLNI[4:7]
          a=exp(parG4[1])
          DRMG4 = log10(exp(DRM.G4_LNI(par=parG4,x=10^(dgr[dgr>(x[2]-((x[2]-x[1])/2))]),q=data.modstan1a$q)))
          if(data.LN$data.modstan1a$shift==T) {DRMG4 = log10(exp(DRM.G4_LNI(par=parG4,x=10^(dgr[dgr>(x[2]-((x[2]-x[1])/2))]),q=data.modstan1a$q))) + log10(exp(1.5*min(data.LN$data.modstan1a$m.org)))
          a = exp(parG4[1]) + 1.5*min(data.LN$data.modstan1a$m.org)}
          lines(dgr[dgr>(x[2]-((x[2]-x[1])/2))],DRMG4,lwd=1, lty=2, col=5)
          segments(x0=dgr[1], y0=log10(exp(a)), x1=max(dgr[dgr<=(x[2]-((x[2]-x[1])/2))]), y1=min(DRMG4),lty=3,col=5)
        }
        if(prior.weights[14]>0){
          parQE4 = QE4resLNI[4:7]
          a=exp(parQE4[1])
          DRMQE4 = log10(exp(DRM.QE4_LNI(par=parQE4,x=10^(dgr[dgr>(x[2]-((x[2]-x[1])/2))]),q=data.modstan1a$q)))
          if(data.LN$data.modstan1a$shift==T) {DRMQE4 = log10(exp(DRM.QE4_LNI(par=parQE4,x=10^(dgr[dgr>(x[2]-((x[2]-x[1])/2))]),q=data.modstan1a$q))) + log10(exp(1.5*min(data.LN$data.modstan1a$m.org)))
          a = exp(parQE4[1]) + 1.5*min(data.LN$data.modstan1a$m.org)}
          lines(dgr[dgr>(x[2]-((x[2]-x[1])/2))],DRMQE4,lwd=1, lty=2, col=6)
          segments(x0=dgr[1], y0=log10(exp(a)), x1=max(dgr[dgr<=(x[2]-((x[2]-x[1])/2))]), y1=min(DRMQE4),lty=3,col=6)
        }
        if(prior.weights[15]>0){
          parP4 = P4resLNI[4:7]
          # a=exp(parP4[1])*pnorm(parP4[3])
          a=exp(parP4[1])*pnorm(qnorm(1 - (log(1+q)/exp(parP4[1])))-exp(parP4[3]))
          DRMP4 = log10(exp(DRM.P4_LNI(par=parP4,x=10^(dgr[dgr>(x[2]-((x[2]-x[1])/2))]),q=data.modstan1a$q)))
          if(data.LN$data.modstan1a$shift==T) {DRMP4 = log10(exp(DRM.P4_LNI(par=parP4,x=10^(dgr[dgr>(x[2]-((x[2]-x[1])/2))]),q=data.modstan1a$q))) + log10(exp(1.5*min(data.LN$data.modstan1a$m.org)))
          a = exp(parP4[1])*pnorm(parP4[3]) + 1.5*min(data.LN$data.modstan1a$m.org)}
          lines(dgr[dgr>(x[2]-((x[2]-x[1])/2))],DRMP4,lwd=1, lty=2, col=7)
          segments(x0=dgr[1], y0=log10(exp(a)), x1=max(dgr[dgr<=(x[2]-((x[2]-x[1])/2))]), y1=min(DRMP4),lty=3,col=7)
        }
        if(prior.weights[16]>0){
          parL4 = L4resLNI[4:7]
          # a=exp(parL4[1])*expit(parL4[3])
          a=exp(parL4[1])*expit(logit(1 - (log(1+q)/exp(parL4[1])))-exp(parL4[3]))
          DRML4 = log10(exp(DRM.L4_LNI(par=parL4,x=10^(dgr[dgr>(x[2]-((x[2]-x[1])/2))]),q=data.modstan1a$q)))
          if(data.LN$data.modstan1a$shift==T) {DRML4 = log10(exp(DRM.L4_LNI(par=parL4,x=10^(dgr[dgr>(x[2]-((x[2]-x[1])/2))]),q=data.modstan1a$q))) + log10(exp(1.5*min(data.LN$data.modstan1a$m.org)))
          a = exp(parL4[1])*expit(parL4[3]) + 1.5*min(data.LN$data.modstan1a$m.org)}
          lines(dgr[dgr>(x[2]-((x[2]-x[1])/2))],DRML4,lwd=1, lty=2, col=8)
          segments(x0=dgr[1], y0=log10(exp(a)), x1=max(dgr[dgr<=(x[2]-((x[2]-x[1])/2))]), y1=min(DRML4),lty=3,col=8)
        }
        leg.all = c(paste0("E4_N (w=",round(weight[1],4),")"),paste0("IE4_N (w=",round(weight[2],4),")"),paste0("H4_N (w=",round(weight[3],4),")"),
                    paste0("LN4_N (w=",round(weight[4],4),")"),paste0("G4_N (w=",round(weight[5],4),")"),paste0("QE4_N (w=",round(weight[6],4),")"),
                    paste0("P4_N (w=",round(weight[7],4),")"),paste0("L4_N (w=",round(weight[8],4),")"),
                    paste0("E4_LN (w=",round(weight[9],4),")"),paste0("IE4_LN (w=",round(weight[10],4),")"),paste0("H4_LN (w=",round(weight[11],4),")"),
                    paste0("LN4_LN (w=",round(weight[12],4),")"),paste0("G4_LN (w=",round(weight[13],4),")"),paste0("QE4_LN (w=",round(weight[14],4),")"),
                    paste0("P4_LN (w=",round(weight[15],4),")"),paste0("L4_LN (w=",round(weight[16],4),")"), "Model averaged BMDL")

        legend(x="bottomright", legend=leg.all[prior.weights>0], col=c(1:8,1:8,1)[prior.weights>0], lwd=c(rep(1,16))[prior.weights>0],
               lty=c(rep(1,8),rep(2,8),NA)[prior.weights>0], pch=c(rep(NA,16),3)[prior.weights>0],
               cex=0.7, title="Model-specific results (fit + BMDL)")
        lines(rbind(x,x,NA),rbind(m-ci,m+ci,NA),lty=3)
        lines(rbind(x,x,NA),rbind(m-ci2,m+ci2,NA))
        points(bmdl[1:16],y=rep(min(m-2*ci),length(bmdl)-1),col=c(1:8,1:8),pch=c(rep(2,8),rep(5,8)))
        points(bmdl[17],y=min(m-1.5*ci),pch=3)
        mtext("Bridge sampling")

      }

      print(plot2())

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
                        BMDMixture = exp(mabmd)*data.modstan1a$maxD,
                        BMDMixture.conv = exp(mabmd.conv)*data.modstan1a$maxD,
                        divergences = divergences,
                        data = data.frame(
                          dose = c(data.N$data.modstan1a$x),
                          sd = sqrt(data.N$data.modstan1a$s2),
                          m = data.N$data.modstan1a$m),
                        max.dose = data.N$data.modstan1a$maxD,
                        q = data.N$data.modstan1a$q,
                        increasing = T,
                        models_included = modelnames[prior.weights > 0]
    )

    attr(ret_results, "class") <- c("BMADR", "BS")

    return(ret_results)

    ### DECREASING

  }else if(data.N$increasing == FALSE){

    ## Normal distribution

    data.modstan1a=data.N$data.modstan1a
    data.modstan1bG=data.N$data.modstan1bG
    data.modstan1bQ=data.N$data.modstan1bQ
    data.modstan1aLN=data.N$data.modstan1aLN
    data.modstan2L=data.N$data.modstan2L
    data.modstan2P=data.N$data.modstan2P
    data.modstanE=data.N$data.modstanE
    svF1=data.N$svF1
    svF1LN = data.N$svF1LN
    svF2P=data.N$svF2P
    svF2L=data.N$svF2L

    BMDL=c(); BMD=c(); BMDU=c()
    converged=c()

    ## Obtain model parameters via MCMC sampling
    if(prior.weights[1]>0){
      print(1)
      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          fitstanE4_N = fun_sampling(stanmodels$mE4_ND, data.modstan1a, svF1,
                                     ndraws,nrchains,
                                     nriterations,warmup,
                                     delta,treedepth,seed,pvec)
        }else if(prior.BMD == TRUE){
          fitstanE4_N = fun_sampling(stanmodels$mE4_ND2, data.modstan1a, svF1,
                                     ndraws,nrchains,
                                     nriterations,warmup,
                                     delta,treedepth,seed,pvec)
        }
      }else if(priordist == "Normal"){
        fitstanE4_N = fun_sampling(stanmodels$mE4_ND_n, data.modstan1a, svF1,
                                   ndraws,nrchains,
                                   nriterations,warmup,
                                   delta,treedepth,seed,pvec)
      }
      parsE4N <- par_extract(fitstanE4_N, model_name = "E4_N")
      # parsE4N[,"BMD"] <- parsE4N[,"BMD"]*data.modstan1a$maxD # BMD on original scale


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

      E4resND=exp(quantile(as.matrix(fitstanE4_N)[,"par[2]"],pvec))*data.modstan1a$maxD
      E4resND=c(E4resND,apply(as.matrix(fitstanE4_N),2,median)[c("par[1]","par[2]","par[3]","par[4]","par[5]")])
      names(E4resND)=c("BMDL","BMD","BMDU","at","k","ct","dt","is2t")

      E4outND <- outLP(parsE4N, pvec, data.modstan1a$maxD)


      DRM_E4_N = DRM.E4_ND(E4resND[4:7], data.modstan1a$x, data.modstan1a$q)


      # Covariance between b-d and between BMD-d
      E4covND = c(cov(as.matrix(fitstanE4_N)[,c("b","d")], use="complete.obs")["b","d"],
                  cov(as.matrix(fitstanE4_N)[,c("BMD","d")], use="complete.obs")["BMD","d"])

      E4corrND = c(cor(as.matrix(fitstanE4_N)[,c("b","d")], use="complete.obs")["b","d"],
                   cor(as.matrix(fitstanE4_N)[,c("BMD","d")], use="complete.obs")["BMD","d"])

      BMDL=c(BMDL, E4resND[1]); BMD=c(BMD, E4resND[2]); BMDU=c(BMDU, E4resND[3])
      bridgeE4N = bridgesampling::bridge_sampler(fitstanE4_N, silent = T) # compute log marginal likelihood
    }else{E4resND=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeE4N=NA; converged=c(converged, NA); E4covND=rep(NA,2); E4corrND=rep(NA,2); DRM_E4_N=rep(NA,length(data.modstan1a$x))
    parsE4N <- NA
    div_E4_N <- NA
    E4outND <- t(data.frame(
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

    if(prior.weights[2]>0){
      print(2)
      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          fitstanIE4_N = fun_sampling(stanmodels$mIE4_ND, data.modstan1a, svF1,
                                     ndraws,nrchains,
                                     nriterations,warmup,
                                     delta,treedepth,seed,pvec)
        }else if(prior.BMD == TRUE){
          fitstanIE4_N = fun_sampling(stanmodels$mIE4_ND2, data.modstan1a, svF1,
                                     ndraws,nrchains,
                                     nriterations,warmup,
                                     delta,treedepth,seed,pvec)
        }
      }else if(priordist == "Normal"){
        fitstanIE4_N = fun_sampling(stanmodels$mIE4_ND_n, data.modstan1a, svF1,
                                    ndraws,nrchains,
                                    nriterations,warmup,
                                    delta,treedepth,seed,pvec)
      }
      parsIE4N <- par_extract(fitstanIE4_N, model_name = "IE4_N")
      # parsIE4N[,"BMD"] <- parsIE4N[,"BMD"]*data.modstan1a$maxD # BMD on original scale

      #diagnostics here
      # posterior_diag(model_stan = fitstanIE4_N)
      convergence_stat <- convergence_deci(model_stan = fitstanIE4_N, nrchains = nrchains)
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)

      # divergent transitions
      div_IE4_N <- sum(sapply(rstan::get_sampler_params(fitstanIE4_N, inc_warmup = F),
                              function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))

      IE4resND=exp(quantile(as.matrix(fitstanIE4_N)[,"par[2]"],pvec))*data.modstan1a$maxD
      IE4resND=c(IE4resND,apply(as.matrix(fitstanIE4_N),2,median)[c("par[1]","par[2]","par[3]","par[4]","par[5]")])
      names(IE4resND)=c("BMDL","BMD","BMDU","at","k","ct","dt","is2t")

      IE4outND <- outLP(parsIE4N, pvec, data.modstan1a$maxD)

      DRM_IE4_N = DRM.IE4_ND(IE4resND[4:7], data.modstan1a$x, data.modstan1a$q)


      # Covariance between b-d and between BMD-d
      IE4covND = c(cov(as.matrix(fitstanIE4_N)[,c("b","d")], use="complete.obs")["b","d"],
                   cov(as.matrix(fitstanIE4_N)[,c("BMD","d")], use="complete.obs")["BMD","d"])

      IE4corrND = c(cor(as.matrix(fitstanIE4_N)[,c("b","d")], use="complete.obs")["b","d"],
                    cor(as.matrix(fitstanIE4_N)[,c("BMD","d")], use="complete.obs")["BMD","d"])

      BMDL=c(BMDL, IE4resND[1]); BMD=c(BMD, IE4resND[2]); BMDU=c(BMDU, IE4resND[3])
      bridgeIE4N = bridgesampling::bridge_sampler(fitstanIE4_N, silent = T)

    }else{IE4resND=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeIE4N=NA; converged=c(converged, NA); IE4covND=rep(NA,2); IE4corrND=rep(NA,2); DRM_IE4_N=rep(NA,length(data.modstan1a$x))
    parsIE4N <- NA
    div_IE4_N <- NA

    IE4outND <- t(data.frame(
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

    if(prior.weights[3]>0){
      print(3)
      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          fitstanH4_N = fun_sampling(stanmodels$mH4_ND, data.modstan1a, svF1,
                                     ndraws,nrchains,
                                     nriterations,warmup,
                                     delta,treedepth,seed,pvec)
        }else if(prior.BMD == TRUE){
          fitstanH4_N = fun_sampling(stanmodels$mH4_ND2, data.modstan1a, svF1,
                                     ndraws,nrchains,
                                     nriterations,warmup,
                                     delta,treedepth,seed,pvec)
        }
      }else if(priordist == "Normal"){
        fitstanH4_N = fun_sampling(stanmodels$mH4_ND_n, data.modstan1a, svF1,
                                   ndraws,nrchains,
                                   nriterations,warmup,
                                   delta,treedepth,seed,pvec)
      }
      parsH4N <- par_extract(fitstanH4_N, model_name = "H4_N")
      # parsH4N[,"BMD"] <- parsH4N[,"BMD"]*data.modstan1a$maxD # BMD on original scale

      #diagnostics here
      # posterior_diag(model_stan = fitstanH4_N)
      convergence_stat <- convergence_deci(model_stan = fitstanH4_N, nrchains = nrchains)
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)

      # divergent transitions
      div_H4_N <- sum(sapply(rstan::get_sampler_params(fitstanH4_N, inc_warmup = F),
                             function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))

      H4resND=exp(quantile(as.matrix(fitstanH4_N)[,"par[2]"],pvec))*data.modstan1a$maxD
      H4resND=c(H4resND,apply(as.matrix(fitstanH4_N),2,median)[c("par[1]","par[2]","par[3]","par[4]","par[5]")])
      names(H4resND)=c("BMDL","BMD","BMDU","at","k","ct","dt","is2t")

      H4outND <- outLP(parsH4N, pvec, data.modstan1a$maxD)


      DRM_H4_N = DRM.H4_ND(H4resND[4:7], data.modstan1a$x, data.modstan1a$q)

      # Covariance between b-d and between BMD-d
      H4covND = c(cov(as.matrix(fitstanH4_N)[,c("b","d")], use="complete.obs")["b","d"],
                  cov(as.matrix(fitstanH4_N)[,c("BMD","d")], use="complete.obs")["BMD","d"])

      H4corrND = c(cor(as.matrix(fitstanH4_N)[,c("b","d")], use="complete.obs")["b","d"],
                   cor(as.matrix(fitstanH4_N)[,c("BMD","d")], use="complete.obs")["BMD","d"])

      BMDL=c(BMDL, H4resND[1]); BMD=c(BMD, H4resND[2]); BMDU=c(BMDU, H4resND[3])
      bridgeH4N = bridgesampling::bridge_sampler(fitstanH4_N, silent = T)
    }else{H4resND=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeH4N=NA; converged=c(converged, NA); H4covND=rep(NA,2); H4corrND=rep(NA,2); DRM_H4_N=rep(NA,length(data.modstan1a$x))
    parsH4N <- NA
    div_H4_N <- NA

    H4outND <- t(data.frame(
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

    if(prior.weights[4]>0){
      print(4)
      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          fitstanLN4_N = fun_sampling(stanmodels$mLN4_ND, data.modstan1aLN, svF1LN,
                                     ndraws,nrchains,
                                     nriterations,warmup,
                                     delta,treedepth,seed,pvec)
        }else if(prior.BMD == TRUE){
          fitstanLN4_N = fun_sampling(stanmodels$mLN4_ND2, data.modstan1aLN, svF1LN,
                                     ndraws,nrchains,
                                     nriterations,warmup,
                                     delta,treedepth,seed,pvec)
        }
      }else if(priordist == "Normal"){
        fitstanLN4_N = fun_sampling(stanmodels$mLN4_ND_n, data.modstan1aLN, svF1LN,
                                    ndraws,nrchains,
                                    nriterations,warmup,
                                    delta,treedepth,seed,pvec)
      }
      parsLN4N <- par_extract(fitstanLN4_N, model_name = "LN4_N")
      # parsLN4N[,"BMD"] <- parsLN4N[,"BMD"]*data.modstan1a$maxD # BMD on original scale

      #diagnostics here
      # posterior_diag(model_stan = fitstanLN4_N)
      convergence_stat <- convergence_deci(model_stan = fitstanLN4_N, nrchains = nrchains)
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)

      # divergent transitions
      div_LN4_N <- sum(sapply(rstan::get_sampler_params(fitstanLN4_N, inc_warmup = F),
                              function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))

      LN4resND=exp(quantile(as.matrix(fitstanLN4_N)[,"par[2]"],pvec))*data.modstan1a$maxD
      LN4resND=c(LN4resND,apply(as.matrix(fitstanLN4_N),2,median)[c("par[1]","par[2]","par[3]","par[4]","par[5]")])
      names(LN4resND)=c("BMDL","BMD","BMDU","at","k","ct","dt","is2t")

      LN4outND <- outLP(parsLN4N, pvec, data.modstan1a$maxD)


      DRM_LN4_N = DRM.LN4_ND(LN4resND[4:7], data.modstan1a$x, data.modstan1a$q)

      # Covariance between b-d and between BMD-d
      LN4covND = c(cov(as.matrix(fitstanLN4_N)[,c("b","d")], use="complete.obs")["b","d"],
                   cov(as.matrix(fitstanLN4_N)[,c("BMD","d")], use="complete.obs")["BMD","d"])

      LN4corrND = c(cor(as.matrix(fitstanLN4_N)[,c("b","d")], use="complete.obs")["b","d"],
                    cor(as.matrix(fitstanLN4_N)[,c("BMD","d")], use="complete.obs")["BMD","d"])

      BMDL=c(BMDL, LN4resND[1]); BMD=c(BMD, LN4resND[2]); BMDU=c(BMDU, LN4resND[3])
      bridgeLN4N = bridgesampling::bridge_sampler(fitstanLN4_N, silent = T)
    }else{LN4resND=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeLN4N=NA; converged=c(converged, NA); LN4covND=rep(NA,2); LN4corrND=rep(NA,2); DRM_LN4_N=rep(NA,length(data.modstan1aLN$x))
    parsLN4N <- NA
    div_LN4_N <- NA

    LN4outND <- t(data.frame(
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

    if(prior.weights[5]>0){
      print(5)
      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          fitstanG4_N = fun_sampling(stanmodels$mG4_ND, data.modstan1bG, svF1,
                                     ndraws,nrchains,
                                     nriterations,warmup,
                                     delta,treedepth,seed,pvec)
        }else if(prior.BMD == TRUE){
          fitstanG4_N = fun_sampling(stanmodels$mG4_ND2, data.modstan1bG, svF1,
                                     ndraws,nrchains,
                                     nriterations,warmup,
                                     delta,treedepth,seed,pvec)
        }
      }else if(priordist == "Normal"){
        fitstanG4_N = fun_sampling(stanmodels$mG4_ND_n, data.modstan1bG, svF1,
                                   ndraws,nrchains,
                                   nriterations,warmup,
                                   delta,treedepth,seed,pvec)
      }

      parsG4N <- par_extract(fitstanG4_N, model_name = "G4_N")
      # parsG4N[,"BMD"] <- parsG4N[,"BMD"]*data.modstan1a$maxD # BMD on original scale

      #diagnostics here
      # posterior_diag(model_stan = fitstanG4_N)
      convergence_stat <- convergence_deci(model_stan = fitstanG4_N, nrchains = nrchains)
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)

      # divergent transitions
      div_G4_N <- sum(sapply(rstan::get_sampler_params(fitstanG4_N, inc_warmup = F),
                             function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))

      G4resND=exp(quantile(as.matrix(fitstanG4_N)[,"par[2]"],pvec))*data.modstan1a$maxD
      G4resND=c(G4resND,apply(as.matrix(fitstanG4_N),2,median)[c("par[1]","par[2]","par[3]","par[4]","par[5]")])
      names(G4resND)=c("BMDL","BMD","BMDU","at","k","ct","dt","is2t")

      G4outND <- outLP(parsG4N, pvec, data.modstan1a$maxD)

      DRM_G4_N = DRM.G4_ND(G4resND[4:7], data.modstan1a$x, data.modstan1a$q)


      # Covariance between b-d and between BMD-d
      G4covND = c(cov(as.matrix(fitstanG4_N)[,c("b","d")], use="complete.obs")["b","d"],
                  cov(as.matrix(fitstanG4_N)[,c("BMD","d")], use="complete.obs")["BMD","d"])

      G4corrND = c(cor(as.matrix(fitstanG4_N)[,c("b","d")], use="complete.obs")["b","d"],
                   cor(as.matrix(fitstanG4_N)[,c("BMD","d")], use="complete.obs")["BMD","d"])

      BMDL=c(BMDL, G4resND[1]); BMD=c(BMD, G4resND[2]); BMDU=c(BMDU, G4resND[3])
      bridgeG4N = bridgesampling::bridge_sampler(fitstanG4_N, silent = T)
    }else{G4resND=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeG4N=NA; converged=c(converged, NA); G4covND=rep(NA,2); G4corrND=rep(NA,2); DRM_G4_N=rep(NA,length(data.modstan1bG$x))
    parsG4N <- NA
    div_G4_N <- NA

    G4outND <- t(data.frame(
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

    if(prior.weights[6]>0){
      print(6)
      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          fitstanQE4_N = fun_sampling(stanmodels$mQE4_ND, data.modstan1bQ, svF1,
                                     ndraws,nrchains,
                                     nriterations,warmup,
                                     delta,treedepth,seed,pvec)
        }else if(prior.BMD == TRUE){
          fitstanQE4_N = fun_sampling(stanmodels$mQE4_ND2, data.modstan1bQ, svF1,
                                     ndraws,nrchains,
                                     nriterations,warmup,
                                     delta,treedepth,seed,pvec)
        }
      }else if(priordist == "Normal"){
        fitstanQE4_N = fun_sampling(stanmodels$mQE4_ND_n, data.modstan1bQ, svF1,
                                    ndraws,nrchains,
                                    nriterations,warmup,
                                    delta,treedepth,seed,pvec)
      }

      parsQE4N <- par_extract(fitstanQE4_N, model_name = "QE4_N")
      # parsQE4N[,"BMD"] <- parsQE4N[,"BMD"]*data.modstan1a$maxD # BMD on original scale

      #diagnostics here
      # posterior_diag(model_stan = fitstanQE4_N)
      convergence_stat <- convergence_deci(model_stan = fitstanQE4_N, nrchains = nrchains)
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)

      # divergent transitions
      div_QE4_N <- sum(sapply(rstan::get_sampler_params(fitstanQE4_N, inc_warmup = F),
                              function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))

      QE4resND=exp(quantile(as.matrix(fitstanQE4_N)[,"par[2]"],pvec))*data.modstan1a$maxD
      QE4resND=c(QE4resND,apply(as.matrix(fitstanQE4_N),2,median)[c("par[1]","par[2]","par[3]","par[4]","par[5]")])
      names(QE4resND)=c("BMDL","BMD","BMDU","at","k","ct","dt","is2t")

      QE4outND <- outLP(parsQE4N, pvec, data.modstan1a$maxD)

      DRM_QE4_N = DRM.QE4_ND(QE4resND[4:7], data.modstan1a$x, data.modstan1a$q)


      # Covariance between b-d and between BMD-d
      QE4covND = c(cov(as.matrix(fitstanQE4_N)[,c("b","d")], use="complete.obs")["b","d"],
                   cov(as.matrix(fitstanQE4_N)[,c("BMD","d")], use="complete.obs")["BMD","d"])

      QE4corrND = c(cor(as.matrix(fitstanQE4_N)[,c("b","d")], use="complete.obs")["b","d"],
                    cor(as.matrix(fitstanQE4_N)[,c("BMD","d")], use="complete.obs")["BMD","d"])

      BMDL=c(BMDL, QE4resND[1]); BMD=c(BMD, QE4resND[2]); BMDU=c(BMDU, QE4resND[3])
      bridgeQE4N = bridgesampling::bridge_sampler(fitstanQE4_N, silent = T)
    }else{QE4resND=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeQE4N=NA; converged=c(converged, NA); QE4covND=rep(NA,2); QE4corrND=rep(NA,2); DRM_QE4_N=rep(NA,length(data.modstan1bQ$x))
    parsQE4N <- NA
    div_QE4_N <- NA

    QE4outND <- t(data.frame(
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

    if(prior.weights[7]>0){
      print(7)
      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          fitstanP4_N = fun_sampling(stanmodels$mP4_ND, data.modstan2P, svF2P,
                                     ndraws,nrchains,
                                     nriterations,warmup,
                                     delta,treedepth,seed,pvec)
        }else if(prior.BMD == TRUE){
          fitstanP4_N = fun_sampling(stanmodels$mP4_ND2, data.modstan2P, svF2P,
                                     ndraws,nrchains,
                                     nriterations,warmup,
                                     delta,treedepth,seed,pvec)
        }
      }else if(priordist == "Normal"){
        fitstanP4_N = fun_sampling(stanmodels$mP4_ND_n, data.modstan2P, svF2P,
                                   ndraws,nrchains,
                                   nriterations,warmup,
                                   delta,treedepth,seed,pvec)
      }
      parsP4N <- par_extract(fitstanP4_N, model_name = "P4_N")
      # parsP4N[,"BMD"] <- parsP4N[,"BMD"]*data.modstan1a$maxD # BMD on original scale

      #diagnostics here
      # posterior_diag(model_stan = fitstanP4_N)
      convergence_stat <- convergence_deci(model_stan = fitstanP4_N, nrchains = nrchains)
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)

      # divergent transitions
      div_P4_N <- sum(sapply(rstan::get_sampler_params(fitstanP4_N, inc_warmup = F),
                             function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))

      P4resND=exp(quantile(as.matrix(fitstanP4_N)[,"par[2]"],pvec))*data.modstan1a$maxD
      P4resND=c(P4resND,apply(as.matrix(fitstanP4_N),2,median)[c("par[1]","par[2]","par[3]","par[4]","par[5]")])
      names(P4resND)=c("BMDL","BMD","BMDU","at","k","ct","dt","is2t")

      P4outND <- outLP(parsP4N, pvec, data.modstan1a$maxD)

      DRM_P4_N = DRM.P4_ND(P4resND[4:7], data.modstan1a$x, data.modstan1a$q)


      # Covariance between b-d and between BMD-d
      P4covND = c(cov(as.matrix(fitstanP4_N)[,c("b","d")], use="complete.obs")["b","d"],
                  cov(as.matrix(fitstanP4_N)[,c("BMD","d")], use="complete.obs")["BMD","d"])

      P4corrND = c(cor(as.matrix(fitstanP4_N)[,c("b","d")], use="complete.obs")["b","d"],
                   cor(as.matrix(fitstanP4_N)[,c("BMD","d")], use="complete.obs")["BMD","d"])

      BMDL=c(BMDL, P4resND[1]); BMD=c(BMD, P4resND[2]); BMDU=c(BMDU, P4resND[3])
      bridgeP4N = bridgesampling::bridge_sampler(fitstanP4_N, silent = T)
    }else{P4resND=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeP4N=NA; converged=c(converged, NA); P4covND=rep(NA,2); P4corrND=rep(NA,2); DRM_P4_N=rep(NA,length(data.modstan2P$x))
    parsP4N <- NA
    div_P4_N <- NA

    P4outND <- t(data.frame(
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

    if(prior.weights[8]>0){
      print(8)
      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          fitstanL4_N = fun_sampling(stanmodels$mL4_ND, data.modstan2L, svF2L,
                                     ndraws,nrchains,
                                     nriterations,warmup,
                                     delta,treedepth,seed,pvec)
        }else if(prior.BMD == TRUE){
          fitstanL4_N = fun_sampling(stanmodels$mL4_ND2, data.modstan2L, svF2L,
                                     ndraws,nrchains,
                                     nriterations,warmup,
                                     delta,treedepth,seed,pvec)
        }
      }else if(priordist == "Normal"){
        fitstanL4_N = fun_sampling(stanmodels$mL4_ND_n, data.modstan2L, svF2L,
                                   ndraws,nrchains,
                                   nriterations,warmup,
                                   delta,treedepth,seed,pvec)
      }
      parsL4N <- par_extract(fitstanL4_N, model_name = "L4_N")
      # parsL4N[,"BMD"] <- parsL4N[,"BMD"]*data.modstan1a$maxD # BMD on original scale

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

      L4resND=exp(quantile(as.matrix(fitstanL4_N)[,"par[2]"],pvec))*data.modstan1a$maxD
      L4resND=c(L4resND,apply(as.matrix(fitstanL4_N),2,median)[c("par[1]","par[2]","par[3]","par[4]","par[5]")])
      names(L4resND)=c("BMDL","BMD","BMDU","at","k","ct","dt","is2t")

      L4outND <- outLP(parsL4N, pvec, data.modstan1a$maxD)


      DRM_L4_N = DRM.L4_ND(L4resND[4:7], data.modstan1a$x, data.modstan1a$q)


      # Covariance between b-d and between BMD-d
      L4covND = c(cov(as.matrix(fitstanL4_N)[,c("b","d")], use="complete.obs")["b","d"],
                  cov(as.matrix(fitstanL4_N)[,c("BMD","d")], use="complete.obs")["BMD","d"])

      L4corrND = c(cor(as.matrix(fitstanL4_N)[,c("b","d")], use="complete.obs")["b","d"],
                   cor(as.matrix(fitstanL4_N)[,c("BMD","d")], use="complete.obs")["BMD","d"])

      BMDL=c(BMDL, L4resND[1]); BMD=c(BMD, L4resND[2]); BMDU=c(BMDU, L4resND[3])
      bridgeL4N = bridgesampling::bridge_sampler(fitstanL4_N, silent = T)
    }else{L4resND=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeL4N=NA; converged=c(converged, NA); L4covND=rep(NA,2); L4corrND=rep(NA,2); DRM_L4_N=rep(NA,length(data.modstan2L$x))
    parsL4N <- NA
    div_L4_N <- NA

    L4outND <- t(data.frame(
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

    ## Lognormal distribution

    data.modstan1a=data.LN$data.modstan1a
    data.modstan1bG=data.LN$data.modstan1bG
    data.modstan1bQ=data.LN$data.modstan1bQ
    data.modstan1aLN=data.LN$data.modstan1aLN
    data.modstan2L=data.LN$data.modstan2L
    data.modstan2P=data.LN$data.modstan2P
    data.modstanE=data.LN$data.modstanE
    svF1=data.LN$svF1
    svF1LN=data.LN$svF1LN
    svF2P=data.LN$svF2P
    svF2L=data.LN$svF2L

    # model specific results
    if(prior.weights[9]>0){
      print(9)
      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          fitstanE4_LN = fun_sampling(stanmodels$mE4_LND, data.modstan1a, svF1,
                                     ndraws,nrchains,
                                     nriterations,warmup,
                                     delta,treedepth,seed,pvec)
        }else if(prior.BMD == TRUE){
          fitstanE4_LN = fun_sampling(stanmodels$mE4_LND2, data.modstan1a, svF1,
                                     ndraws,nrchains,
                                     nriterations,warmup,
                                     delta,treedepth,seed,pvec)
        }
      }else if(priordist == "Normal"){
        fitstanE4_LN = fun_sampling(stanmodels$mE4_LND_n, data.modstan1a, svF1,
                                    ndraws,nrchains,
                                    nriterations,warmup,
                                    delta,treedepth,seed,pvec)
      }

      parsE4LN <- par_extract(fitstanE4_LN, model_name = "E4_LN")
      # parsE4LN[,"BMD"] <- parsE4LN[,"BMD"]*data.modstan1a$maxD # BMD on original scale

      #diagnostics here
      # posterior_diag(model_stan = fitstanE4_LN)
      convergence_stat <- convergence_deci(model_stan = fitstanE4_LN, nrchains = nrchains)
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)

      # divergent transitions
      div_E4_LN <- sum(sapply(rstan::get_sampler_params(fitstanE4_LN, inc_warmup = F),
                              function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))

      E4resLND=exp(quantile(as.matrix(fitstanE4_LN)[,"par[2]"],pvec))*data.modstan1a$maxD
      E4resLND=c(E4resLND,apply(as.matrix(fitstanE4_LN),2,median)[c("par[1]","par[2]","par[3]","par[4]","par[5]")])
      names(E4resLND)=c("BMDL","BMD","BMDU","at","k","ct","dt","is2t")

      E4outLND <- outLP(parsE4LN, pvec, data.modstan1a$maxD)

      DRM_E4_LN = exp(DRM.E4_LND(E4resLND[4:7], data.modstan1a$x, data.modstan1a$q))

      # Covariance between b-d and between BMD-d
      E4covLND = c(cov(as.matrix(fitstanE4_LN)[,c("b","d")], use="complete.obs")["b","d"],
                   cov(as.matrix(fitstanE4_LN)[,c("BMD","d")], use="complete.obs")["BMD","d"])

      E4corrLND = c(cor(as.matrix(fitstanE4_LN)[,c("b","d")], use="complete.obs")["b","d"],
                    cor(as.matrix(fitstanE4_LN)[,c("BMD","d")], use="complete.obs")["BMD","d"])

      BMDL=c(BMDL, E4resLND[1]); BMD=c(BMD, E4resLND[2]); BMDU=c(BMDU, E4resLND[3])
      bridgeE4LN = bridgesampling::bridge_sampler(fitstanE4_LN, silent = T)
    }else{E4resLND=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeE4LN=NA; converged=c(converged, NA); E4covLND=rep(NA,2); E4corrLND=rep(NA,2); DRM_E4_LN=rep(NA,length(data.modstan1a$x))
    parsE4LN <- NA
    div_E4_LN <- NA

    E4outLND <- t(data.frame(
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

    if(prior.weights[10]>0){
      print(10)
      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          fitstanIE4_LN = fun_sampling(stanmodels$mIE4_LND, data.modstan1a, svF1,
                                      ndraws,nrchains,
                                      nriterations,warmup,
                                      delta,treedepth,seed,pvec)
        }else if(prior.BMD == TRUE){
          fitstanIE4_LN = fun_sampling(stanmodels$mIE4_LND2, data.modstan1a, svF1,
                                      ndraws,nrchains,
                                      nriterations,warmup,
                                      delta,treedepth,seed,pvec)
        }
      }else if(priordist == "Normal"){
        fitstanIE4_LN = fun_sampling(stanmodels$mIE4_LND_n, data.modstan1a, svF1,
                                     ndraws,nrchains,
                                     nriterations,warmup,
                                     delta,treedepth,seed,pvec)
      }

      parsIE4LN <- par_extract(fitstanIE4_LN, model_name = "IE4_LN")
      # parsIE4LN[,"BMD"] <- parsIE4LN[,"BMD"]*data.modstan1a$maxD # BMD on original scale

      #diagnostics here
      # posterior_diag(model_stan = fitstanIE4_LN)
      convergence_stat <- convergence_deci(model_stan = fitstanIE4_LN, nrchains = nrchains)
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)

      # divergent transitions
      div_IE4_LN <- sum(sapply(rstan::get_sampler_params(fitstanIE4_LN, inc_warmup = F),
                               function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))

      IE4resLND=exp(quantile(as.matrix(fitstanIE4_LN)[,"par[2]"],pvec))*data.modstan1a$maxD
      IE4resLND=c(IE4resLND,apply(as.matrix(fitstanIE4_LN),2,median)[c("par[1]","par[2]","par[3]","par[4]","par[5]")])
      names(IE4resLND)=c("BMDL","BMD","BMDU","at","k","ct","dt","is2t")

      IE4outLND <- outLP(parsIE4LN, pvec, data.modstan1a$maxD)


      DRM_IE4_LN = exp(DRM.IE4_LND(IE4resLND[4:7], data.modstan1a$x, data.modstan1a$q))


      IE4covLND = c(cov(as.matrix(fitstanIE4_LN)[,c("b","d")], use="complete.obs")["b","d"],
                    cov(as.matrix(fitstanIE4_LN)[,c("BMD","d")], use="complete.obs")["BMD","d"])

      IE4corrLND = c(cor(as.matrix(fitstanIE4_LN)[,c("b","d")], use="complete.obs")["b","d"],
                     cor(as.matrix(fitstanIE4_LN)[,c("BMD","d")], use="complete.obs")["BMD","d"])

      BMDL=c(BMDL, IE4resLND[1]); BMD=c(BMD, IE4resLND[2]); BMDU=c(BMDU, IE4resLND[3])
      bridgeIE4LN = bridgesampling::bridge_sampler(fitstanIE4_LN, silent = T)
    }else{IE4resLND=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeIE4LN=NA; converged=c(converged, NA); IE4covLND=rep(NA,2); IE4corrLND=rep(NA,2); DRM_IE4_LN=rep(NA,length(data.modstan1a$x))
    parsIE4LN <- NA
    div_IE4_LN <- NA

    IE4outLND <- t(data.frame(
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

    if(prior.weights[11]>0){
      print(11)
      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          fitstanH4_LN = fun_sampling(stanmodels$mH4_LND, data.modstan1a, svF1,
                                      ndraws,nrchains,
                                      nriterations,warmup,
                                      delta,treedepth,seed,pvec)
        }else if(prior.BMD == TRUE){
          fitstanH4_LN = fun_sampling(stanmodels$mH4_LND2, data.modstan1a, svF1,
                                      ndraws,nrchains,
                                      nriterations,warmup,
                                      delta,treedepth,seed,pvec)
        }
      }else if(priordist == "Normal"){
        fitstanH4_LN = fun_sampling(stanmodels$mH4_LND_n, data.modstan1a, svF1,
                                    ndraws,nrchains,
                                    nriterations,warmup,
                                    delta,treedepth,seed,pvec)
      }

      parsH4LN <- par_extract(fitstanH4_LN, model_name = "H4_LN")
      # parsH4LN[,"BMD"] <- parsH4LN[,"BMD"]*data.modstan1a$maxD # BMD on original scale

      #diagnostics here
      # posterior_diag(model_stan = fitstanH4_LN)
      convergence_stat <- convergence_deci(model_stan = fitstanH4_LN, nrchains = nrchains)
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)

      # divergent transitions
      div_H4_LN <- sum(sapply(rstan::get_sampler_params(fitstanH4_LN, inc_warmup = F),
                              function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))

      H4resLND=exp(quantile(as.matrix(fitstanH4_LN)[,"par[2]"],pvec))*data.modstan1a$maxD
      H4resLND=c(H4resLND,apply(as.matrix(fitstanH4_LN),2,median)[c("par[1]","par[2]","par[3]","par[4]","par[5]")])
      names(H4resLND)=c("BMDL","BMD","BMDU","at","k","ct","dt","is2t")

      H4outLND <- outLP(parsH4LN, pvec, data.modstan1a$maxD)

      DRM_H4_LN = exp(DRM.H4_LND(H4resLND[4:7], data.modstan1a$x, data.modstan1a$q))

      H4covLND = c(cov(as.matrix(fitstanH4_LN)[,c("b","d")], use="complete.obs")["b","d"],
                   cov(as.matrix(fitstanH4_LN)[,c("BMD","d")], use="complete.obs")["BMD","d"])

      H4corrLND = c(cor(as.matrix(fitstanH4_LN)[,c("b","d")], use="complete.obs")["b","d"],
                    cor(as.matrix(fitstanH4_LN)[,c("BMD","d")], use="complete.obs")["BMD","d"])

      BMDL=c(BMDL, H4resLND[1]); BMD=c(BMD, H4resLND[2]); BMDU=c(BMDU, H4resLND[3])
      bridgeH4LN = bridgesampling::bridge_sampler(fitstanH4_LN, silent = T)
    }else{H4resLND=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeH4LN=NA; converged=c(converged, NA); H4covLND=rep(NA,2); H4corrLND=rep(NA,2); DRM_H4_LN=rep(NA,length(data.modstan1a$x))
    parsH4LN <- NA
    div_H4_LN <- NA

    H4outLND <- t(data.frame(
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

    if(prior.weights[12]>0){
      print(12)
      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          fitstanLN4_LN = fun_sampling(stanmodels$mLN4_LND, data.modstan1aLN, svF1LN,
                                      ndraws,nrchains,
                                      nriterations,warmup,
                                      delta,treedepth,seed,pvec)
        }else if(prior.BMD == TRUE){
          fitstanLN4_LN = fun_sampling(stanmodels$mLN4_LND2, data.modstan1aLN, svF1LN,
                                      ndraws,nrchains,
                                      nriterations,warmup,
                                      delta,treedepth,seed,pvec)
        }
      }else if(priordist == "Normal"){
        fitstanLN4_LN = fun_sampling(stanmodels$mLN4_LND_n, data.modstan1aLN, svF1LN,
                                     ndraws,nrchains,
                                     nriterations,warmup,
                                     delta,treedepth,seed,pvec)
      }

      parsLN4LN <- par_extract(fitstanLN4_LN, model_name = "LN4_LN")
      # parsLN4LN[,"BMD"] <- parsLN4LN[,"BMD"]*data.modstan1a$maxD # BMD on original scale

      #diagnostics here
      # posterior_diag(model_stan = fitstanLN4_LN)
      convergence_stat <- convergence_deci(model_stan = fitstanLN4_LN, nrchains = nrchains)
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)

      # divergent transitions
      div_LN4_LN <- sum(sapply(rstan::get_sampler_params(fitstanLN4_LN, inc_warmup = F),
                               function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))

      LN4resLND=exp(quantile(as.matrix(fitstanLN4_LN)[,"par[2]"],pvec))*data.modstan1a$maxD
      LN4resLND=c(LN4resLND,apply(as.matrix(fitstanLN4_LN),2,median)[c("par[1]","par[2]","par[3]","par[4]","par[5]")])
      names(LN4resLND)=c("BMDL","BMD","BMDU","at","k","ct","dt","is2t")

      LN4outLND <- outLP(parsLN4LN, pvec, data.modstan1a$maxD)


      DRM_LN4_LN = exp(DRM.LN4_LND(LN4resLND[4:7], data.modstan1a$x, data.modstan1a$q))


      LN4covLND = c(cov(as.matrix(fitstanLN4_LN)[,c("b","d")], use="complete.obs")["b","d"],
                    cov(as.matrix(fitstanLN4_LN)[,c("BMD","d")], use="complete.obs")["BMD","d"])

      LN4corrLND = c(cor(as.matrix(fitstanLN4_LN)[,c("b","d")], use="complete.obs")["b","d"],
                     cor(as.matrix(fitstanLN4_LN)[,c("BMD","d")], use="complete.obs")["BMD","d"])

      BMDL=c(BMDL, LN4resLND[1]); BMD=c(BMD, LN4resLND[2]); BMDU=c(BMDU, LN4resLND[3])
      bridgeLN4LN = bridgesampling::bridge_sampler(fitstanLN4_LN, silent = T)
    }else{LN4resLND=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeLN4LN=NA; converged=c(converged, NA); LN4covLND=rep(NA,2); LN4corrLND=rep(NA,2); DRM_LN4_LN=rep(NA,length(data.modstan1aLN$x))
    parsLN4LN <- NA
    div_LN4_LN <- NA

    LN4outLND <- t(data.frame(
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

    if(prior.weights[13]>0){
      print(13)
      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          fitstanG4_LN = fun_sampling(stanmodels$mG4_LND, data.modstan1bG, svF1,
                                      ndraws,nrchains,
                                      nriterations,warmup,
                                      delta,treedepth,seed,pvec)
        }else if(prior.BMD == TRUE){
          fitstanG4_LN = fun_sampling(stanmodels$mG4_LND2, data.modstan1bG, svF1,
                                      ndraws,nrchains,
                                      nriterations,warmup,
                                      delta,treedepth,seed,pvec)
        }
      }else if(priordist == "Normal"){
        fitstanG4_LN = fun_sampling(stanmodels$mG4_LND_n, data.modstan1bG, svF1,
                                    ndraws,nrchains,
                                    nriterations,warmup,
                                    delta,treedepth,seed,pvec)
      }

      parsG4LN <- par_extract(fitstanG4_LN, model_name = "G4_LN")
      # parsG4LN[,"BMD"] <- parsG4LN[,"BMD"]*data.modstan1a$maxD # BMD on original scale


      #diagnostics here
      # posterior_diag(model_stan = fitstanG4_LN)
      convergence_stat <- convergence_deci(model_stan = fitstanG4_LN, nrchains = nrchains)
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)

      # divergent transitions
      div_G4_LN <- sum(sapply(rstan::get_sampler_params(fitstanG4_LN, inc_warmup = F),
                              function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))

      G4resLND=exp(quantile(as.matrix(fitstanG4_LN)[,"par[2]"],pvec))*data.modstan1a$maxD
      G4resLND=c(G4resLND,apply(as.matrix(fitstanG4_LN),2,median)[c("par[1]","par[2]","par[3]","par[4]","par[5]")])
      names(G4resLND)=c("BMDL","BMD","BMDU","at","k","ct","dt","is2t")

      G4outLND <- outLP(parsG4LN, pvec, data.modstan1a$maxD)

      DRM_G4_LN = exp(DRM.G4_LND(G4resLND[4:7], data.modstan1a$x, data.modstan1a$q))


      G4covLND = c(cov(as.matrix(fitstanG4_LN)[,c("b","d")], use="complete.obs")["b","d"],
                   cov(as.matrix(fitstanG4_LN)[,c("BMD","d")], use="complete.obs")["BMD","d"])

      G4corrLND = c(cor(as.matrix(fitstanG4_LN)[,c("b","d")], use="complete.obs")["b","d"],
                    cor(as.matrix(fitstanG4_LN)[,c("BMD","d")], use="complete.obs")["BMD","d"])

      BMDL=c(BMDL, G4resLND[1]); BMD=c(BMD, G4resLND[2]); BMDU=c(BMDU, G4resLND[3])
      bridgeG4LN = bridgesampling::bridge_sampler(fitstanG4_LN, silent = T)
    }else{G4resLND=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeG4LN=NA; converged=c(converged, NA); G4covLND=rep(NA,2); G4corrLND=rep(NA,2); DRM_G4_LN=rep(NA,length(data.modstan1bG$x))
    parsG4LN <- NA
    div_G4_LN <- NA

    G4outLND <- t(data.frame(
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

    if(prior.weights[14]>0){
      print(14)
      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          fitstanQE4_LN = fun_sampling(stanmodels$mQE4_LND, data.modstan1bQ, svF1,
                                      ndraws,nrchains,
                                      nriterations,warmup,
                                      delta,treedepth,seed,pvec)
        }else if(prior.BMD == TRUE){
          fitstanQE4_LN = fun_sampling(stanmodels$mQE4_LND2, data.modstan1bQ, svF1,
                                      ndraws,nrchains,
                                      nriterations,warmup,
                                      delta,treedepth,seed,pvec)
        }
      }else if(priordist == "Normal"){
        fitstanQE4_LN = fun_sampling(stanmodels$mQE4_LND_n, data.modstan1bQ, svF1,
                                     ndraws,nrchains,
                                     nriterations,warmup,
                                     delta,treedepth,seed,pvec)
      }

      parsQE4LN <- par_extract(fitstanQE4_LN, model_name = "QE4_LN")
      # parsQE4LN[,"BMD"] <- parsQE4LN[,"BMD"]*data.modstan1a$maxD # BMD on original scale

      #diagnostics here
      # posterior_diag(model_stan = fitstanQE4_LN)
      convergence_stat <- convergence_deci(model_stan = fitstanQE4_LN, nrchains = nrchains)
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)

      # divergent transitions
      div_QE4_LN <- sum(sapply(rstan::get_sampler_params(fitstanQE4_LN, inc_warmup = F),
                               function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))

      QE4resLND=exp(quantile(as.matrix(fitstanQE4_LN)[,"par[2]"],pvec))*data.modstan1a$maxD
      QE4resLND=c(QE4resLND,apply(as.matrix(fitstanQE4_LN),2,median)[c("par[1]","par[2]","par[3]","par[4]","par[5]")])
      names(QE4resLND)=c("BMDL","BMD","BMDU","at","k","ct","dt","is2t")

      QE4outLND <- outLP(parsQE4LN, pvec, data.modstan1a$maxD)

      DRM_QE4_LN = exp(DRM.QE4_LND(QE4resLND[4:7], data.modstan1a$x, data.modstan1a$q))


      QE4covLND = c(cov(as.matrix(fitstanQE4_LN)[,c("b","d")], use="complete.obs")["b","d"],
                    cov(as.matrix(fitstanQE4_LN)[,c("BMD","d")], use="complete.obs")["BMD","d"])

      QE4corrLND = c(cor(as.matrix(fitstanQE4_LN)[,c("b","d")], use="complete.obs")["b","d"],
                     cor(as.matrix(fitstanQE4_LN)[,c("BMD","d")], use="complete.obs")["BMD","d"])

      BMDL=c(BMDL, QE4resLND[1]); BMD=c(BMD, QE4resLND[2]); BMDU=c(BMDU, QE4resLND[3])
      bridgeQE4LN = bridgesampling::bridge_sampler(fitstanQE4_LN, silent = T)
    }else{QE4resLND=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeQE4LN=NA; converged=c(converged, NA); QE4covLND=rep(NA,2); QE4corrLND=rep(NA,2); DRM_QE4_LN=rep(NA,length(data.modstan1bQ$x))
    parsQE4LN <- NA
    div_QE4_LN <- NA

    QE4outLND <- t(data.frame(
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

    if(prior.weights[15]>0){
      print(15)
      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          fitstanP4_LN = fun_sampling(stanmodels$mP4_LND, data.modstan2P, svF2P,
                                      ndraws,nrchains,
                                      nriterations,warmup,
                                      delta,treedepth,seed,pvec)
        }else if(prior.BMD == TRUE){
          fitstanP4_LN = fun_sampling(stanmodels$mP4_LND2, data.modstan2P, svF2P,
                                      ndraws,nrchains,
                                      nriterations,warmup,
                                      delta,treedepth,seed,pvec)
        }
      }else if(priordist == "Normal"){
        fitstanP4_LN = fun_sampling(stanmodels$mP4_LND_n, data.modstan2P, svF2P,
                                    ndraws,nrchains,
                                    nriterations,warmup,
                                    delta,treedepth,seed,pvec)
      }

      parsP4LN <- par_extract(fitstanP4_LN, model_name = "P4_LN")
      # parsP4LN[,"BMD"] <- parsP4LN[,"BMD"]*data.modstan1a$maxD # BMD on original scale

      #diagnostics here
      # posterior_diag(model_stan = fitstanP4_LN)
      convergence_stat <- convergence_deci(model_stan = fitstanP4_LN, nrchains = nrchains)
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)

      # divergent transitions
      div_P4_LN <- sum(sapply(rstan::get_sampler_params(fitstanP4_LN, inc_warmup = F),
                              function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))

      P4resLND=exp(quantile(as.matrix(fitstanP4_LN)[,"par[2]"],pvec))*data.modstan1a$maxD
      P4resLND=c(P4resLND,apply(as.matrix(fitstanP4_LN),2,median)[c("par[1]","par[2]","par[3]","par[4]","par[5]")])
      names(P4resLND)=c("BMDL","BMD","BMDU","at","k","ct","dt","is2t")

      P4outLND <- outLP(parsP4LN, pvec, data.modstan1a$maxD)

      DRM_P4_LN = exp(DRM.P4_LND(P4resLND[4:7], data.modstan1a$x, data.modstan1a$q))


      P4covLND = c(cov(as.matrix(fitstanP4_LN)[,c("b","d")], use="complete.obs")["b","d"],
                   cov(as.matrix(fitstanP4_LN)[,c("BMD","d")], use="complete.obs")["BMD","d"])

      P4corrLND = c(cor(as.matrix(fitstanP4_LN)[,c("b","d")], use="complete.obs")["b","d"],
                    cor(as.matrix(fitstanP4_LN)[,c("BMD","d")], use="complete.obs")["BMD","d"])

      BMDL=c(BMDL, P4resLND[1]); BMD=c(BMD, P4resLND[2]); BMDU=c(BMDU, P4resLND[3])
      bridgeP4LN = bridgesampling::bridge_sampler(fitstanP4_LN, silent = T)
    }else{P4resLND=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeP4LN=NA; converged=c(converged, NA); P4covLND=rep(NA,2); P4corrLND=rep(NA,2); DRM_P4_LN=rep(NA,length(data.modstan2P$x))
    parsP4LN <- NA
    div_P4_LN <- NA

    P4outLND <- t(data.frame(
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

    if(prior.weights[16]>0){
      print(16)
      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          fitstanL4_LN = fun_sampling(stanmodels$mL4_LND, data.modstan2L, svF2L,
                                      ndraws,nrchains,
                                      nriterations,warmup,
                                      delta,treedepth,seed,pvec)
        }else if(prior.BMD == TRUE){
          fitstanL4_LN = fun_sampling(stanmodels$mL4_LND2, data.modstan2L, svF2L,
                                      ndraws,nrchains,
                                      nriterations,warmup,
                                      delta,treedepth,seed,pvec)
        }
      }else if(priordist == "Normal"){
        fitstanL4_LN = fun_sampling(stanmodels$mL4_LND_n, data.modstan2L, svF2L,
                                    ndraws,nrchains,
                                    nriterations,warmup,
                                    delta,treedepth,seed,pvec)
      }
      parsL4LN <- par_extract(fitstanL4_LN, model_name = "L4_LN")
      # parsL4LN[,"BMD"] <- parsL4LN[,"BMD"]*data.modstan1a$maxD # BMD on original scale


      #diagnostics here
      # posterior_diag(model_stan = fitstanL4_LN)
      convergence_stat <- convergence_deci(model_stan = fitstanL4_LN, nrchains = nrchains)
      conv_lp = 0
      if(convergence_stat$Rhat[convergence_stat$Parameters=="lp__"] < 1.01) conv_lp = 1
      converged = c(converged, conv_lp)

      # divergent transitions
      div_L4_LN <- sum(sapply(rstan::get_sampler_params(fitstanL4_LN, inc_warmup = F),
                              function(x) sum(x[, "divergent__"])))/(nrchains*(nriterations-warmup))

      L4resLND=exp(quantile(as.matrix(fitstanL4_LN)[,"par[2]"],pvec))*data.modstan1a$maxD
      L4resLND=c(L4resLND,apply(as.matrix(fitstanL4_LN),2,median)[c("par[1]","par[2]","par[3]","par[4]","par[5]")])
      names(L4resLND)=c("BMDL","BMD","BMDU","at","k","ct","dt","is2t")

      L4outLND <- outLP(parsL4LN, pvec, data.modstan1a$maxD)

      DRM_L4_LN = exp(DRM.L4_LND(L4resLND[4:7], data.modstan1a$x, data.modstan1a$q))


      L4covLND = c(cov(as.matrix(fitstanL4_LN)[,c("b","d")], use="complete.obs")["b","d"],
                   cov(as.matrix(fitstanL4_LN)[,c("BMD","d")], use="complete.obs")["BMD","d"])

      L4corrLND = c(cor(as.matrix(fitstanL4_LN)[,c("b","d")], use="complete.obs")["b","d"],
                    cor(as.matrix(fitstanL4_LN)[,c("BMD","d")], use="complete.obs")["BMD","d"])

      BMDL=c(BMDL, L4resLND[1]); BMD=c(BMD, L4resLND[2]); BMDU=c(BMDU, L4resLND[3])
      bridgeL4LN = bridgesampling::bridge_sampler(fitstanL4_LN, silent = T)
    }else{L4resLND=NULL; BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); bridgeL4LN=NA; converged=c(converged, NA); L4covLND=rep(NA,2); L4corrLND=rep(NA,2); DRM_L4_LN=rep(NA,length(data.modstan2L$x))
    parsL4LN <- NA
    div_L4_LN <- NA

    L4outLND <- t(data.frame(
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
    mabmd=(c( # normal
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
    macib=exp(quantile(mabmd,pvec))*data.modstan1a$maxD
    names(macib)=c("BMDL","BMD","BMDU") # original scale

    BMDq_bs = exp(quantile(mabmd, seq(0,1,0.005)))*data.modstan1a$maxD


    mods = c("E4_N","IE4_N","H4_N","LN4_N","G4_N","QE4_N","P4_N","L4_N",
             "E4_LN","IE4_LN","H4_LN","LN4_LN","G4_LN","QE4_LN","P4_LN","L4_LN")
    w.bs = rep(0, 16)
    names(w.bs) = mods
    for(i in mods[mods%in%names(lpwb)]){
      w.bs[i] = lpwb[i]
    }

    ### Model-averaged response per dose level
    dr.MA.bs <- c()
    for(i in 1:length(data.modstan1a$x)){
      dr.MA.bs[i] = weighted.mean(x = c(DRM_E4_N[i],DRM_IE4_N[i],DRM_H4_N[i],DRM_LN4_N[i],DRM_G4_N[i],DRM_QE4_N[i], DRM_P4_N[i], DRM_L4_N[i] ,
                                        DRM_E4_LN[i], DRM_IE4_LN[i], DRM_H4_LN[i], DRM_LN4_LN[i], DRM_G4_LN[i],DRM_QE4_LN[i], DRM_P4_LN[i],DRM_L4_LN[i]),
                                  w = w.bs,
                                  na.rm = T)

    }


    ##### Weights & MA if one of the models is divergent --> this model gets weight 0
    if(0 %in% converged){
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
      mabmd.conv=(c( # normal
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
      macib.conv=exp(quantile(mabmd.conv,pvec))*data.modstan1a$maxD
      names(macib.conv)=c("BMDL","BMD","BMDU")

      BMDq_bs_conv = exp(quantile(mabmd.conv, seq(0,1,0.005)))*data.modstan1a$maxD


      mods = c("E4_N","IE4_N","H4_N","LN4_N","G4_N","QE4_N","P4_N","L4_N",
               "E4_LN","IE4_LN","H4_LN","LN4_LN","G4_LN","QE4_LN","P4_LN","L4_LN")
      w.bs.conv = rep(0, 16)
      names(w.bs.conv) = mods
      for(i in mods[mods%in%names(lpwb.conv)]){
        w.bs.conv[i] = lpwb.conv[i]
      }

      ### Model-averaged response per dose level
      dr.MA.bs.conv <- c()
      for(i in 1:length(data.modstan1a$x)){
        dr.MA.bs.conv[i] = weighted.mean(x = c(DRM_E4_N[i],DRM_IE4_N[i],DRM_H4_N[i],DRM_LN4_N[i],DRM_G4_N[i],DRM_QE4_N[i], DRM_P4_N[i], DRM_L4_N[i] ,
                                               DRM_E4_LN[i], DRM_IE4_LN[i], DRM_H4_LN[i], DRM_LN4_LN[i], DRM_G4_LN[i],DRM_QE4_LN[i], DRM_P4_LN[i],DRM_L4_LN[i]),
                                         w = w.bs.conv,
                                         na.rm = T)
      }

    }else{
      w.bs.conv = NULL; macib.conv = NULL; BMDq_bs_conv = NULL; dr.MA.bs.conv = NULL
    }


    #----------------------------------
    ### weights based on laplace approximation

    # normal
    data.modstan1a=data.N$data.modstan1a
    data.modstan1bG=data.N$data.modstan1bG
    data.modstan1bQ=data.N$data.modstan1bQ
    data.modstan1aLN=data.N$data.modstan1aLN
    data.modstan2L=data.N$data.modstan2L
    data.modstan2P=data.N$data.modstan2P
    data.modstanE=data.N$data.modstanE
    svF1=data.N$svF1
    svF1LN=data.N$svF1LN
    svF2P=data.N$svF2P
    svF2L=data.N$svF2L

    llN = c() # likelihoods

    # getting the posterior modes and the hessian and the model specific posterior distributions
    if(prior.weights[1]>0){
      print("pw1")
      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          optE4_ND <- fun_optim(stanmodels$mE4_ND, data.modstan1a, svF1, ndraws, 123, pvec)
        }else if(prior.BMD == TRUE){
          optE4_ND <- fun_optim(stanmodels$mE4_ND2, data.modstan1a, svF1, ndraws, 123, pvec)
        }
        }else if(priordist == "Normal"){
        optE4_ND <- fun_optim(stanmodels$mE4_ND_n, data.modstan1a, svF1, ndraws, 123, pvec)
      }

      llE4N=llfE4_ND(optE4_ND$par[1:5],nvec=data.modstan1a$n,
                     dvec=data.modstan1a$x,mvec=data.modstan1a$m,
                     s2vec=data.modstan1a$s2,qval=data.modstan1a$q)


      llN = c(llN, llE4N)
    }else{llN = c(llN,NA)}
    #
    if(prior.weights[2]>0){
      print("pw2")
      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          optIE4_ND <- fun_optim(stanmodels$mIE4_ND, data.modstan1a, svF1, ndraws, 123, pvec)
        }else if(prior.BMD == TRUE){
          optIE4_ND <- fun_optim(stanmodels$mIE4_ND2, data.modstan1a, svF1, ndraws, 123, pvec)
        }
        }else if(priordist == "Normal"){
        optIE4_ND <- fun_optim(stanmodels$mIE4_ND_n, data.modstan1a, svF1, ndraws, 123, pvec)
      }

      llIE4N=llfIE4_ND(optIE4_ND$par[1:5],nvec=data.modstan1a$n,
                       dvec=data.modstan1a$x,mvec=data.modstan1a$m,
                       s2vec=data.modstan1a$s2,qval=data.modstan1a$q)

      llN = c(llN, llIE4N)
    }else{llN = c(llN,NA)}
    #
    if(prior.weights[3]>0){
      print("pw3")
      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          optH4_ND <- fun_optim(stanmodels$mH4_ND, data.modstan1a, svF1, ndraws, 123, pvec)
        }else if(prior.BMD == TRUE){
          optH4_ND <- fun_optim(stanmodels$mH4_ND2, data.modstan1a, svF1, ndraws, 123, pvec)
        }
        }else if(priordist == "Normal"){
        optH4_ND <- fun_optim(stanmodels$mH4_ND_n, data.modstan1a, svF1, ndraws, 123, pvec)
      }

      llH4N=llfH4_ND(optH4_ND$par[1:5],nvec=data.modstan1a$n,
                     dvec=data.modstan1a$x,mvec=data.modstan1a$m,
                     s2vec=data.modstan1a$s2,qval=data.modstan1a$q)

      llN = c(llN, llH4N)
    }else{llN = c(llN,NA)}
    #
    if(prior.weights[4]>0){
      print("pw4")
      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          optLN4_ND <- fun_optim(stanmodels$mLN4_ND, data.modstan1aLN, svF1LN, ndraws, 123, pvec)
        }else if(prior.BMD == TRUE){
          optLN4_ND <- fun_optim(stanmodels$mLN4_ND2, data.modstan1aLN, svF1LN, ndraws, 123, pvec)
        }
        }else if(priordist == "Normal"){
        optLN4_ND <- fun_optim(stanmodels$mLN4_ND_n, data.modstan1aLN, svF1LN, ndraws, 123, pvec)
      }

      llLN4N=llfLN4_ND(optLN4_ND$par[1:5],nvec=data.modstan1aLN$n,
                       dvec=data.modstan1aLN$x,mvec=data.modstan1aLN$m,
                       s2vec=data.modstan1aLN$s2,qval=data.modstan1aLN$q)


      llN = c(llN, llLN4N)
    }else{llN = c(llN,NA)}
    #
    if(prior.weights[5]>0){
      print("pw5")
      # optG4_ND = fun_optim_G4(modN, stanmodels$mG4_ND, stv = svF1, data.modstan1bG, ndraws, 123, pvec)
      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          optG4_ND <- fun_optim(stanmodels$mG4_ND, data.modstan1bG, svF1, ndraws, 123, pvec)
        }else if(prior.BMD == TRUE){
          optG4_ND <- fun_optim(stanmodels$mG4_ND2, data.modstan1bG, svF1, ndraws, 123, pvec)
        }
        }else if(priordist == "Normal"){
        optG4_ND <- fun_optim(stanmodels$mG4_ND_n, data.modstan1bG, svF1, ndraws, 123, pvec)
      }

      llG4N=llfG4_ND(optG4_ND$par[1:5],nvec=data.modstan1bG$n,
                     dvec=data.modstan1bG$x,mvec=data.modstan1bG$m,
                     s2vec=data.modstan1bG$s2,qval=data.modstan1bG$q)


      llN = c(llN, llG4N)
    }else{llN = c(llN,NA)}
    #
    if(prior.weights[6]>0){
      print("pw6")
      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          optQE4_ND <- fun_optim(stanmodels$mQE4_ND, data.modstan1bQ, svF1, ndraws, 123, pvec)
        }else if(prior.BMD == TRUE){
          optQE4_ND <- fun_optim(stanmodels$mQE4_ND2, data.modstan1bQ, svF1, ndraws, 123, pvec)
        }
        }else if(priordist == "Normal"){
        optQE4_ND <- fun_optim(stanmodels$mQE4_ND_n, data.modstan1bQ, svF1, ndraws, 123, pvec)
      }

      llQE4N=llfQE4_ND(optQE4_ND$par[1:5],nvec=data.modstan1bQ$n,
                       dvec=data.modstan1bQ$x,mvec=data.modstan1bQ$m,
                       s2vec=data.modstan1bQ$s2,qval=data.modstan1bQ$q)


      llN = c(llN, llQE4N)
    }else{llN = c(llN,NA)}
    #
    if(prior.weights[7]>0){
      print("pw7")
      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          optP4_ND <- fun_optim(stanmodels$mP4_ND, data.modstan2P, svF2P, ndraws, 123, pvec)
        }else if(prior.BMD == TRUE){
          optP4_ND <- fun_optim(stanmodels$mP4_ND2, data.modstan2P, svF2P, ndraws, 123, pvec)
        }
        }else if(priordist == "Normal"){
        optP4_ND <- fun_optim(stanmodels$mP4_ND_n, data.modstan2P, svF2P, ndraws, 123, pvec)
      }

      llP4N=llfP4_ND(optP4_ND$par[1:5],nvec=data.modstan2P$n,
                     dvec=data.modstan2P$x,mvec=data.modstan2P$m,
                     s2vec=data.modstan2P$s2,qval=data.modstan2P$q)

      llN = c(llN, llP4N)
    }else{llN = c(llN,NA)}
    #
    if(prior.weights[8]>0){
      print("pw8")
      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          optL4_ND <- fun_optim(stanmodels$mL4_ND, data.modstan2L, svF2L, ndraws, 123, pvec)
        }else if(prior.BMD == TRUE){
          optL4_ND <- fun_optim(stanmodels$mL4_ND2, data.modstan2L, svF2L, ndraws, 123, pvec)
        }
        }else if(priordist == "Normal"){
        optL4_ND <- fun_optim(stanmodels$mL4_ND_n, data.modstan2L, svF2L, ndraws, 123, pvec)
      }

      llL4N=llfL4_ND(optL4_ND$par[1:5],nvec=data.modstan2L$n,
                     dvec=data.modstan2L$x,mvec=data.modstan2L$m,
                     s2vec=data.modstan2L$s2,qval=data.modstan2L$q)


      llN = c(llN, llL4N)
    }else{llN = c(llN,NA)}


    # lognormal
    data.modstan1a=data.LN$data.modstan1a
    data.modstan1bG=data.LN$data.modstan1bG
    data.modstan1bQ=data.LN$data.modstan1bQ
    data.modstan1aLN=data.LN$data.modstan1aLN
    data.modstan2L=data.LN$data.modstan2L
    data.modstan2P=data.LN$data.modstan2P
    data.modstanE=data.LN$data.modstanE
    svF1=data.LN$svF1
    svF1LN=data.LN$svF1LN
    svF2P=data.LN$svF2P
    svF2L=data.LN$svF2L

    llLN = c()

    if(prior.weights[9]>0){
      print("pw9")
      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          optE4_LND <- fun_optim(stanmodels$mE4_LND, data.modstan1a, svF1, ndraws, 123, pvec)
        }else if(prior.BMD == TRUE){
          optE4_LND <- fun_optim(stanmodels$mE4_LND2, data.modstan1a, svF1, ndraws, 123, pvec)
        }
        }else if(priordist == "Normal"){
        optE4_LND <- fun_optim(stanmodels$mE4_LND_n, data.modstan1a, svF1, ndraws, 123, pvec)
      }

      llE4LN=llfE4_LND(optE4_LND$par[1:5],nvec=data.modstan1a$n,
                       dvec=data.modstan1a$x,mvec=data.modstan1a$m,
                       s2vec=data.modstan1a$s2,qval=data.modstan1a$q)


      llLN = c(llLN, llE4LN)
    }else{llLN = c(llLN,NA)}
    #
    if(prior.weights[10]>0){
      print("pw10")
      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          optIE4_LND <- fun_optim(stanmodels$mIE4_LND, data.modstan1a, svF1, ndraws, 123, pvec)
        }else if(prior.BMD == TRUE){
          optIE4_LND <- fun_optim(stanmodels$mIE4_LND2, data.modstan1a, svF1, ndraws, 123, pvec)
        }
        }else if(priordist == "Normal"){
        optIE4_LND <- fun_optim(stanmodels$mIE4_LND_n, data.modstan1a, svF1, ndraws, 123, pvec)
      }

      llIE4LN=llfIE4_LND(optIE4_LND$par[1:5],nvec=data.modstan1a$n,
                         dvec=data.modstan1a$x,mvec=data.modstan1a$m,
                         s2vec=data.modstan1a$s2,qval=data.modstan1a$q)


      llLN = c(llLN, llIE4LN)
    }else{llLN = c(llLN,NA)}
    #
    if(prior.weights[11]>0){
      print("pw11")
      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          optH4_LND <- fun_optim(stanmodels$mH4_LND, data.modstan1a, svF1, ndraws, 123, pvec)
        }else if(prior.BMD == TRUE){
          optH4_LND <- fun_optim(stanmodels$mH4_LND2, data.modstan1a, svF1, ndraws, 123, pvec)
        }
        }else if(priordist == "Normal"){
        optH4_LND <- fun_optim(stanmodels$mH4_LND_n, data.modstan1a, svF1, ndraws, 123, pvec)
      }

      llH4LN=llfH4_LND(optH4_LND$par[1:5],nvec=data.modstan1a$n,
                       dvec=data.modstan1a$x,mvec=data.modstan1a$m,
                       s2vec=data.modstan1a$s2,qval=data.modstan1a$q)

      llLN = c(llLN, llH4LN)
    }else{llLN = c(llLN,NA)}
    #
    if(prior.weights[12]>0){
      print("pw12")
      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          optLN4_LND <- fun_optim(stanmodels$mLN4_LND, data.modstan1aLN, svF1LN, ndraws, 123, pvec)
        }else if(prior.BMD == TRUE){
          optLN4_LND <- fun_optim(stanmodels$mLN4_LND2, data.modstan1aLN, svF1LN, ndraws, 123, pvec)
        }
        }else if(priordist == "Normal"){
        optLN4_LND <- fun_optim(stanmodels$mLN4_LND_n, data.modstan1aLN, svF1LN, ndraws, 123, pvec)
      }

      llLN4LN=llfLN4_LND(optLN4_LND$par[1:5],nvec=data.modstan1aLN$n,
                         dvec=data.modstan1aLN$x,mvec=data.modstan1aLN$m,
                         s2vec=data.modstan1aLN$s2,qval=data.modstan1aLN$q)

      llLN = c(llLN, llLN4LN)
    }else{llLN = c(llLN,NA)}
    #
    if(prior.weights[13]>0){

      print("pw13")
      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          optG4_LND <- fun_optim(stanmodels$mG4_LND, data.modstan1bG, svF1, ndraws, 123, pvec)
        }else if(prior.BMD == TRUE){
          optG4_LND <- fun_optim(stanmodels$mG4_LND2, data.modstan1bG, svF1, ndraws, 123, pvec)
        }
        }else if(priordist == "Normal"){
        optG4_LND <- fun_optim(stanmodels$mG4_LND_n, data.modstan1bG, svF1, ndraws, 123, pvec)
      }

      llG4LN=llfG4_LND(optG4_LND$par[1:5],nvec=data.modstan1bG$n,
                       dvec=data.modstan1bG$x,mvec=data.modstan1bG$m,
                       s2vec=data.modstan1bG$s2,qval=data.modstan1bG$q)

      llLN = c(llLN, llG4LN)
    }else{llLN = c(llLN,NA)}
    #
    if(prior.weights[14]>0){
      print("pw14")
      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          optQE4_LND <- fun_optim(stanmodels$mQE4_LND, data.modstan1bQ, svF1, ndraws, 123, pvec)
        }else if(prior.BMD == TRUE){
          optQE4_LND <- fun_optim(stanmodels$mQE4_LND2, data.modstan1bQ, svF1, ndraws, 123, pvec)
        }
        }else if(priordist == "Normal"){
        optQE4_LND <- fun_optim(stanmodels$mQE4_LND_n, data.modstan1bQ, svF1, ndraws, 123, pvec)
      }

      llQE4LN=llfQE4_LND(optQE4_LND$par[1:5],nvec=data.modstan1bQ$n,
                         dvec=data.modstan1bQ$x,mvec=data.modstan1bQ$m,
                         s2vec=data.modstan1bQ$s2,qval=data.modstan1bQ$q)

      llLN = c(llLN, llQE4LN)
    }else{llLN = c(llLN,NA)}
    #
    if(prior.weights[15]>0){
      print("pw15")
      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          optP4_LND <- fun_optim(stanmodels$mP4_LND, data.modstan2P, svF2P, ndraws, 123, pvec)
        }else if(prior.BMD == TRUE){
          optP4_LND <- fun_optim(stanmodels$mP4_LND2, data.modstan2P, svF2P, ndraws, 123, pvec)
        }
        }else if(priordist == "Normal"){
        optP4_LND <- fun_optim(stanmodels$mP4_LND_n, data.modstan2P, svF2P, ndraws, 123, pvec)
      }

      llP4LN=llfP4_LND(optP4_LND$par[1:5],nvec=data.modstan2P$n,
                       dvec=data.modstan2P$x,mvec=data.modstan2P$m,
                       s2vec=data.modstan2P$s2,qval=data.modstan2P$q)

      llLN = c(llLN, llP4LN)
    }else{llLN = c(llLN,NA)}
    #
    if(prior.weights[16]>0){
      print("pw16")
      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          optL4_LND <- fun_optim(stanmodels$mL4_LND, data.modstan2L, svF2L, ndraws, 123, pvec)
        }else if(prior.BMD == TRUE){
          optL4_LND <- fun_optim(stanmodels$mL4_LND2, data.modstan2L, svF2L, ndraws, 123, pvec)
        }
        }else if(priordist == "Normal"){
        optL4_LND <- fun_optim(stanmodels$mL4_LND_n, data.modstan2L, svF2L, ndraws, 123, pvec)
      }

      llL4LN=llfL4_LND(optL4_LND$par[1:5],nvec=data.modstan2L$n,
                       dvec=data.modstan2L$x,mvec=data.modstan2L$m,
                       s2vec=data.modstan2L$s2,qval=data.modstan2L$q)

      llLN = c(llLN, llL4LN)
    }else{llLN = c(llLN,NA)}

    minll = min(llN,llLN,na.rm=T)

    # the weights

    fun.w <- function(DIH, ll, min.ll, opt, mu, sig, prior, priorbmd, lb, ub, s1, s2, s3){
      if(prior == "PERT"){
        if(priorbmd == FALSE){
          w1 <- (2*pi)^(2.5)*sqrt(DIH)*exp(ll-min.ll)*
            mvtnorm::dmvnorm(opt$par[c(4,5)],mean=mu[c(4,5)],
                             sigma=sig[c(4,5),c(4,5)])*
            dexp(-opt$par[2],rate=mu[2])*
            mc2d::dpert(opt$par[1], min = lb[1], max = ub[1],
                        mode = mu[1], shape = s1)*
            mc2d::dpert(opt$par[3], min = lb[3], max = ub[3],
                        mode = mu[3], shape = s3)
        }else if(priorbmd == TRUE){
          w1 <- (2*pi)^(2.5)*sqrt(DIH)*exp(ll-min.ll)*
            mvtnorm::dmvnorm(opt$par[c(4,5)],mean=mu[c(4,5)],
                             sigma=sig[c(4,5),c(4,5)])*
            mc2d::dpert(opt$par[1], min = lb[1], max = ub[1],
                        mode = mu[1], shape = s1)*
            mc2d::dpert(opt$par[2], min = lb[2], max = ub[2],
                        mode = mu[2], shape = s2)*
            mc2d::dpert(opt$par[3], min = lb[3], max = ub[3],
                        mode = mu[3], shape = s3)
        }
      }else if(prior == "Normal"){
        w1 <- (2*pi)^(2.5)*sqrt(DIH)*exp(ll-min.ll)*
          mvtnorm::dmvnorm(opt$par[c(1,2,4,5)],mean=mu[c(1,2,4,5)],
                           sigma=sig[c(1,2,4,5),c(1,2,4,5)])*
          dexp(-opt$par[2],rate=mu[2])
      }
      return(w1)
    }

    w=c()

    # normal

    data.modstan1a=data.N$data.modstan1a
    data.modstan1bG=data.N$data.modstan1bG
    data.modstan1bQ=data.N$data.modstan1bQ
    data.modstan1aLN=data.N$data.modstan1aLN
    data.modstan2L=data.N$data.modstan2L
    data.modstan2P=data.N$data.modstan2P
    data.modstanE=data.N$data.modstanE
    svF1=data.N$svF1
    svF1LN = data.N$svF1LN
    svF2P=data.N$svF2P
    svF2L=data.N$svF2L

    if(prior.weights[1]>0){
      # determinant of covariance matrix (obtained by examining the curvature of the posterior at the mode)
      DIHE4h = det(-solve(optE4_ND$hessian))
      DIHE4=ifelse(DIHE4h<0,0,DIHE4h)
      ## Aproximation of marginal (i.e. integrated) likelihood (= 'model evidence')
      w = c(w, fun.w(DIHE4, llE4N, minll, optE4_ND, data.modstan1a$priormu, data.modstan1a$priorSigma, priordist, prior.BMD,
                     data.modstan1a$priorlb, data.modstan1a$priorub, data.modstan1a$shape.a, data.modstan1a$shape.BMD, data.modstan1a$shape.c))
    }else{w=c(w,0)}

    if(prior.weights[2]>0){
      DIHIE4h = det(-solve(optIE4_ND$hessian))
      DIHIE4=ifelse(DIHIE4h<0,0,DIHIE4h)
      w = c(w, fun.w(DIHIE4, llIE4N, minll, optIE4_ND, data.modstan1a$priormu, data.modstan1a$priorSigma, priordist, prior.BMD,
                     data.modstan1a$priorlb, data.modstan1a$priorub, data.modstan1a$shape.a, data.modstan1a$shape.BMD, data.modstan1a$shape.c))
    }else{w=c(w,0)}

    if(prior.weights[3]>0){
      DIHH4h = det(-solve(optH4_ND$hessian))
      DIHH4=ifelse(DIHH4h<0,0,DIHH4h)
      w = c(w, fun.w(DIHH4, llH4N, minll, optH4_ND, data.modstan1a$priormu, data.modstan1a$priorSigma, priordist, prior.BMD,
                     data.modstan1a$priorlb, data.modstan1a$priorub, data.modstan1a$shape.a, data.modstan1a$shape.BMD, data.modstan1a$shape.c))
    }else{w=c(w,0)}

    if(prior.weights[4]>0){
      DIHLN4h = det(-solve(optLN4_ND$hessian))
      DIHLN4=ifelse(DIHLN4h<0,0,DIHLN4h)
      w = c(w, fun.w(DIHLN4, llLN4N, minll, optLN4_ND, data.modstan1aLN$priormu, data.modstan1aLN$priorSigma, priordist, prior.BMD,
                     data.modstan1aLN$priorlb, data.modstan1aLN$priorub, data.modstan1a$shape.a, data.modstan1a$shape.BMD, data.modstan1a$shape.c))
    }else{w=c(w,0)}

    if(prior.weights[5]>0){
      DIHG4h=det(-solve(optG4_ND$hessian))
      DIHG4=ifelse(DIHG4h<0,0,DIHG4h)
      w = c(w, fun.w(DIHG4, llG4N, minll, optG4_ND, data.modstan1bG$priormu, data.modstan1bG$priorSigma, priordist, prior.BMD,
                     data.modstan1bG$priorlb, data.modstan1bG$priorub, data.modstan1a$shape.a, data.modstan1a$shape.BMD, data.modstan1a$shape.c))
    }else{w=c(w,0)}

    if(prior.weights[6]>0){
      DIHQE4h = det(-solve(optQE4_ND$hessian))
      DIHQE4=ifelse(DIHQE4h<0,0,DIHQE4h)
      w = c(w, fun.w(DIHQE4, llQE4N, minll, optQE4_ND, data.modstan1bQ$priormu, data.modstan1bQ$priorSigma, priordist, prior.BMD,
                     data.modstan1bQ$priorlb, data.modstan1bQ$priorub, data.modstan1a$shape.a, data.modstan1a$shape.BMD, data.modstan1a$shape.c))
    }else{w=c(w,0)}

    if(prior.weights[7]>0){
      DIHP4h = det(-solve(optP4_ND$hessian))
      DIHP4=ifelse(DIHP4h<0,0,DIHP4h)
      w = c(w, fun.w(DIHP4, llP4N, minll, optP4_ND, data.modstan2P$priormu, data.modstan2P$priorSigma, priordist, prior.BMD,
                     data.modstan2P$priorlb, data.modstan2P$priorub, data.modstan1a$shape.a, data.modstan1a$shape.BMD, data.modstan1a$shape.c))
    }else{w=c(w,0)}

    if(prior.weights[8]>0){
      DIHL4h = det(-solve(optL4_ND$hessian))
      DIHL4=ifelse(DIHL4h<0,0,DIHL4h)
      w = c(w, fun.w(DIHL4, llL4N, minll, optL4_ND, data.modstan2L$priormu, data.modstan2L$priorSigma, priordist, prior.BMD,
                     data.modstan2L$priorlb, data.modstan2L$priorub, data.modstan1a$shape.a, data.modstan1a$shape.BMD, data.modstan1a$shape.c))
    }else{w=c(w,0)}


    # lognormal

    data.modstan1a=data.LN$data.modstan1a
    data.modstan1bG=data.LN$data.modstan1bG
    data.modstan1bQ=data.LN$data.modstan1bQ
    data.modstan1aLN=data.LN$data.modstan1aLN
    data.modstan2L=data.LN$data.modstan2L
    data.modstan2P=data.LN$data.modstan2P
    data.modstanE=data.LN$data.modstanE
    svF1=data.LN$svF1
    svF1LN=data.LN$svF1LN
    svF2P=data.LN$svF2P
    svF2L=data.LN$svF2L

    if(prior.weights[9]>0){
      DIHE4h = det(-solve(optE4_LND$hessian))
      DIHE4=ifelse(DIHE4h<0,0,DIHE4h)
      w = c(w, fun.w(DIHE4, llE4LN, minll, optE4_LND, data.modstan1a$priormu, data.modstan1a$priorSigma, priordist, prior.BMD,
                     data.modstan1a$priorlb, data.modstan1a$priorub, data.modstan1a$shape.a, data.modstan1a$shape.BMD, data.modstan1a$shape.c))
    }else{w=c(w,0)}

    if(prior.weights[10]>0){
      DIHIE4h = det(-solve(optIE4_LND$hessian))
      DIHIE4=ifelse(DIHIE4h<0,0,DIHIE4h)
      w = c(w, fun.w(DIHIE4, llIE4LN, minll, optIE4_LND, data.modstan1a$priormu, data.modstan1a$priorSigma, priordist, prior.BMD,
                     data.modstan1a$priorlb, data.modstan1a$priorub, data.modstan1a$shape.a, data.modstan1a$shape.BMD, data.modstan1a$shape.c))
    }else{w=c(w,0)}

    if(prior.weights[11]>0){
      DIHH4h = det(-solve(optH4_LND$hessian))
      DIHH4=ifelse(DIHH4h<0,0,DIHH4h)
      w = c(w, fun.w(DIHH4, llH4LN, minll, optH4_LND, data.modstan1a$priormu, data.modstan1a$priorSigma, priordist, prior.BMD,
                     data.modstan1a$priorlb, data.modstan1a$priorub, data.modstan1a$shape.a, data.modstan1a$shape.BMD, data.modstan1a$shape.c))
    }else{w=c(w,0)}

    if(prior.weights[12]>0){
      DIHLN4h = det(-solve(optLN4_LND$hessian))
      DIHLN4=ifelse(DIHLN4h<0,0,DIHLN4h)
      w = c(w, fun.w(DIHLN4, llLN4LN, minll, optLN4_LND, data.modstan1aLN$priormu, data.modstan1aLN$priorSigma, priordist, prior.BMD,
                     data.modstan1aLN$priorlb, data.modstan1aLN$priorub, data.modstan1a$shape.a, data.modstan1a$shape.BMD, data.modstan1a$shape.c))
    }else{w=c(w,0)}

    if(prior.weights[13]>0){
      DIHG4h = det(-solve(optG4_LND$hessian))
      DIHG4=ifelse(DIHG4h<0,0,DIHG4h)
      w = c(w, fun.w(DIHG4, llG4LN, minll, optG4_LND, data.modstan1bG$priormu, data.modstan1bG$priorSigma, priordist, prior.BMD,
                     data.modstan1bG$priorlb, data.modstan1bG$priorub, data.modstan1a$shape.a, data.modstan1a$shape.BMD, data.modstan1a$shape.c))
    }else{w=c(w,0)}

    if(prior.weights[14]>0){
      DIHQE4h = det(-solve(optQE4_LND$hessian))
      DIHQE4=ifelse(DIHQE4h<0,0,DIHQE4h)
      w = c(w, fun.w(DIHQE4, llQE4LN, minll, optQE4_LND, data.modstan1bQ$priormu, data.modstan1bQ$priorSigma, priordist, prior.BMD,
                     data.modstan1bQ$priorlb, data.modstan1bQ$priorub, data.modstan1a$shape.a, data.modstan1a$shape.BMD, data.modstan1a$shape.c))
    }else{w=c(w,0)}

    if(prior.weights[15]>0){
      DIHP4h = det(-solve(optP4_LND$hessian))
      DIHP4=ifelse(DIHP4h<0,0,DIHP4h)
      w = c(w, fun.w(DIHP4, llP4LN, minll, optP4_LND, data.modstan2P$priormu, data.modstan2P$priorSigma, priordist, prior.BMD,
                     data.modstan2P$priorlb, data.modstan2P$priorub, data.modstan1a$shape.a, data.modstan1a$shape.BMD, data.modstan1a$shape.c))
    }else{w=c(w,0)}

    if(prior.weights[16]>0){
      DIHL4h = det(-solve(optL4_LND$hessian))
      DIHL4=ifelse(DIHL4h<0,0,DIHL4h)
      w = c(w, fun.w(DIHL4, llL4LN, minll, optL4_LND, data.modstan2L$priormu, data.modstan2L$priorSigma, priordist, prior.BMD,
                     data.modstan2L$priorlb, data.modstan2L$priorub, data.modstan1a$shape.a, data.modstan1a$shape.BMD, data.modstan1a$shape.c))
    }else{w=c(w,0)}

    lpwlp=(prior.weights*w)/sum(prior.weights*w)
    # lpwlp = lpwlp*prior.weights
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
    macilp=exp(quantile(mabmd,pvec))*data.modstan1a$maxD
    names(macilp)=c("BMDL","BMD","BMDU") # on original scale

    BMDq_ls = exp(quantile(mabmd, seq(0,1,0.005)))*data.modstan1a$maxD


    ### Model-averaged response per dose level
    dr.MA.ls <- c()
    for(i in 1:length(data.modstan1a$x)){
      dr.MA.ls[i] = weighted.mean(x = c(DRM_E4_N[i],DRM_IE4_N[i],DRM_H4_N[i],DRM_LN4_N[i],DRM_G4_N[i],DRM_QE4_N[i], DRM_P4_N[i], DRM_L4_N[i] ,
                                        DRM_E4_LN[i], DRM_IE4_LN[i], DRM_H4_LN[i], DRM_LN4_LN[i], DRM_G4_LN[i],DRM_QE4_LN[i], DRM_P4_LN[i],DRM_L4_LN[i]),
                                  w = lpwlp,
                                  na.rm = T)
    }


    ##### Weights & MA if one of the models is divergent --> this model gets weight 0
    if(0 %in% converged){

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
      macilp.conv=exp(quantile(mabmd.conv,pvec))*data.modstan1a$maxD
      names(macilp.conv)=c("BMDL","BMD","BMDU")

      BMDq_ls_conv = exp(quantile(mabmd.conv, seq(0,1,0.005)))*data.modstan1a$maxD


      ### Model-averaged response per dose level
      dr.MA.ls.conv <- c()
      for(i in 1:length(data.modstan1a$x)){
        dr.MA.ls.conv[i] = weighted.mean(x = c(DRM_E4_N[i],DRM_IE4_N[i],DRM_H4_N[i],DRM_LN4_N[i],DRM_G4_N[i],DRM_QE4_N[i], DRM_P4_N[i], DRM_L4_N[i] ,
                                               DRM_E4_LN[i], DRM_IE4_LN[i], DRM_H4_LN[i], DRM_LN4_LN[i], DRM_G4_LN[i],DRM_QE4_LN[i], DRM_P4_LN[i],DRM_L4_LN[i]),
                                         w = lpwlp.conv,
                                         na.rm = T)
      }


    }else{
      lpwlp.conv = NULL; macilp.conv = NULL; BMDq_ls_conv = NULL; dr.MA.ls.conv = NULL; mabmd.conv = NA;
    }

    ## Plot with weights bridge sampling
    BMDL = c(BMDL, macib[1]); BMD = c(BMD, macib[2]); BMDU = c(BMDU, macib[3])

    names(BMDL) <- c("E4_N","IE4_N","H4_N","LN4_N","G4_N","QE4_N","P4_N","L4_N","E4_LN","IE4_LN","H4_LN","LN4_LN","G4_LN","QE4_LN","P4_LN","L4_LN","MA")
    model = c("E4_N","IE4_N","H4_N","LN4_N","G4_N","QE4_N","P4_N","L4_N","E4_LN","IE4_LN","H4_LN","LN4_LN","G4_LN","QE4_LN","P4_LN","L4_LN","MA")
    model = as.factor(model)

    weight = c(rep(0,16),1)
    names(weight) = c("E4_N","IE4_N","H4_N","LN4_N","G4_N","QE4_N","P4_N","L4_N","E4_LN","IE4_LN","H4_LN","LN4_LN","G4_LN","QE4_LN","P4_LN","L4_LN","MA")
    for(i in names(lpwb)){
      weight[names(weight)==i] = lpwb[names(lpwb)==i]
    }

    if(plot == TRUE){

      data.modstan1a = data.N$data.modstan1a
      data.modstan1bG = data.N$data.modstan1bG
      data.modstanE = data.N$data.modstanE
      data.modstan1bQ = data.N$data.modstan1bQ
      data.modstan1aLN = data.N$data.modstan1aLN
      data.modstan2P = data.N$data.modstan2P
      data.modstan2L = data.N$data.modstan2L
      svF1 = data.N$svF1
      svF1LN = data.N$svF1LN
      svF2P = data.N$svF2P
      svF2L = data.N$svF2L

      plot2 = function(){
        sd = sqrt(data.modstan1a$s2)
        gmean.a = log(NtoLN(data.modstan1a$m,sd))[1:N]
        gsd.a = log(NtoLN(data.modstan1a$m,sd))[(N+1):(2*N)]
        x = c(log10(data.modstan1a$x[2]/4),log10(data.modstan1a$x)[2:data.modstan1a$N])
        m = log10(exp(gmean.a))
        s = log10(exp(gsd.a))
        bmdl = log10(BMDL/data.modstan1a$maxD)
        x[1] = ifelse(min(bmdl,na.rm=T)<x[1], log10(min(BMDL/data.modstan1a$maxD/2, na.rm=T)), log10(data.modstan1a$x[2]/4))
        g = choose(length(x),2)
        alph = 0.05/(2*g)
        cp = qnorm(1-alph)
        n = data.modstan1a$n
        ci = cp*s/sqrt(n)
        ci2 = 1.96*s/sqrt(n)

        plot(x, m, xlab="log10-Dose", ylab="log10-Response", ylim=c(min(m-2.2*ci),max(m+2.2*ci)), xlim=c(ifelse(min(bmdl,na.rm=T)<min(x),min(bmdl,na.rm=T),min(x)),abs(min(x))))#, ylim=c(min(mean.a)-2*sd.a, max(mean.a)+2*sd.a))
        dgr=seq(min(x), abs(min(x)),by=0.01)
        # plot model specific fits
        # normal distribution
        if(prior.weights[1]>0){
          parE4 = E4resND[4:7]
          DRME4 = log10((DRM.E4_ND(par=parE4,x=10^(dgr[dgr>(x[2]-((x[2]-x[1])/2))]),q=data.modstan1a$q)))
          lines(dgr[dgr>(x[2]-((x[2]-x[1])/2))],DRME4,lwd=1, lty=1, col=1)
          segments(x0=dgr[1], y0=log10(exp(parE4[1])), x1=max(dgr[dgr<=(x[2]-((x[2]-x[1])/2))]), y1=DRME4[1],lty=3,col=1)
        }
        if(prior.weights[2]>0){
          parIE4 = IE4resND[4:7]
          DRMIE4 = log10((DRM.IE4_ND(par=parIE4,x=10^(dgr[dgr>(x[2]-((x[2]-x[1])/2))]),q=data.modstan1a$q)))
          lines(dgr[dgr>(x[2]-((x[2]-x[1])/2))],DRMIE4,lwd=1, lty=1, col=2)
          segments(x0=dgr[1], y0=log10(exp(parIE4[1])), x1=max(dgr[dgr<=(x[2]-((x[2]-x[1])/2))]), y1=DRMIE4[1],lty=3,col=2)
        }
        if(prior.weights[3]>0){
          parH4 = H4resND[4:7]
          DRMH4 = log10((DRM.H4_ND(par=parH4,x=10^(dgr[dgr>(x[2]-((x[2]-x[1])/2))]),q=data.modstan1a$q)))
          lines(dgr[dgr>(x[2]-((x[2]-x[1])/2))],DRMH4,lwd=1, lty=1, col=3)
          segments(x0=dgr[1], y0=log10(exp(parH4[1])), x1=max(dgr[dgr<=(x[2]-((x[2]-x[1])/2))]), y1=DRMH4[1],lty=3,col=3)
        }
        if(prior.weights[4]>0){
          parLN4 = LN4resND[4:7]
          DRMLN4 = log10((DRM.LN4_ND(par=parLN4,x=10^(dgr[dgr>(x[2]-((x[2]-x[1])/2))]),q=data.modstan1a$q)))
          lines(dgr[dgr>(x[2]-((x[2]-x[1])/2))],DRMLN4,lwd=1, lty=1, col=4)
          segments(x0=dgr[1], y0=log10(exp(parLN4[1])), x1=max(dgr[dgr<=(x[2]-((x[2]-x[1])/2))]), y1=DRMLN4[1],lty=3,col=4)
        }
        if(prior.weights[5]>0){
          parG4 = G4resND[4:7]
          DRMG4 = log10((DRM.G4_ND(par=parG4,x=10^(dgr[dgr>(x[2]-((x[2]-x[1])/2))]),q=data.modstan1a$q)))
          lines(dgr[dgr>(x[2]-((x[2]-x[1])/2))],DRMG4,lwd=1, lty=1, col=5)
          segments(x0=dgr[1], y0=log10(exp(parG4[1])), x1=max(dgr[dgr<=(x[2]-((x[2]-x[1])/2))]), y1=DRMG4[1],lty=3,col=5)
        }
        if(prior.weights[6]>0){
          parQE4 = QE4resND[4:7]
          DRMQE4 = log10((DRM.QE4_ND(par=parQE4,x=10^(dgr[dgr>(x[2]-((x[2]-x[1])/2))]),q=data.modstan1a$q)))
          lines(dgr[dgr>(x[2]-((x[2]-x[1])/2))],DRMQE4,lwd=1, lty=1, col=6)
          segments(x0=dgr[1], y0=log10(exp(parQE4[1])), x1=max(dgr[dgr<=(x[2]-((x[2]-x[1])/2))]), y1=DRMQE4[1],lty=3,col=6)
        }
        if(prior.weights[7]>0){
          parP4 = P4resND[4:7]
          DRMP4 = log10((DRM.P4_ND(par=parP4,x=10^(dgr[dgr>(x[2]-((x[2]-x[1])/2))]),q=data.modstan1a$q)))
          lines(dgr[dgr>(x[2]-((x[2]-x[1])/2))],DRMP4,lwd=1, lty=1, col=7)
          segments(x0=dgr[1], y0=log10(exp(parP4[1])), x1=max(dgr[dgr<=(x[2]-((x[2]-x[1])/2))]), y1=DRMP4[1],lty=3,col=7)
        }
        if(prior.weights[8]>0){
          parL4 = L4resND[4:7]
          DRML4 = log10((DRM.L4_ND(par=parL4,x=10^(dgr[dgr>(x[2]-((x[2]-x[1])/2))]),q=data.modstan1a$q)))
          lines(dgr[dgr>(x[2]-((x[2]-x[1])/2))],DRML4,lwd=1, lty=1, col=8)
          segments(x0=dgr[1], y0=log10(exp(parL4[1])), x1=max(dgr[dgr<=(x[2]-((x[2]-x[1])/2))]), y1=DRML4[1],lty=3,col=8)
        }

        # lognormal distribution
        if(prior.weights[9]>0){
          parE4 = E4resLND[4:7]
          a=exp(parE4[1])
          DRME4 = log10(exp(DRM.E4_LND(par=parE4,x=10^(dgr[dgr>(x[2]-((x[2]-x[1])/2))]),q=data.modstan1a$q)))
          if(data.LN$data.modstan1a$shift==T) {DRME4 = log10(exp(DRM.E4_LND(par=parE4,x=10^(dgr[dgr>(x[2]-((x[2]-x[1])/2))]),q=data.modstan1a$q))) + log10(exp(1.5*min(data.LN$data.modstan1a$m.org)))
          a = exp(parE4[1]) + 1.5*min(data.LN$data.modstan1a$m.org)}
          lines(dgr[dgr>(x[2]-((x[2]-x[1])/2))],DRME4,lwd=1, lty=2, col=1)
          segments(x0=dgr[1], y0=log10(exp(a)), x1=max(dgr[dgr<=(x[2]-((x[2]-x[1])/2))]), y1=DRME4[1],lty=3,col=1)
        }
        if(prior.weights[10]>0){
          parIE4 = IE4resLND[4:7]
          a=exp(parIE4[1])
          DRMIE4 = log10(exp(DRM.IE4_LND(par=parIE4,x=10^(dgr[dgr>(x[2]-((x[2]-x[1])/2))]),q=data.modstan1a$q)))
          if(data.LN$data.modstan1a$shift==T) {DRMIE4 = log10(exp(DRM.IE4_LND(par=parIE4,x=10^(dgr[dgr>(x[2]-((x[2]-x[1])/2))]),q=data.modstan1a$q))) + log10(exp(1.5*min(data.LN$data.modstan1a$m.org)))
          a = exp(parIE4[1]) + 1.5*min(data.LN$data.modstan1a$m.org)}
          lines(dgr[dgr>(x[2]-((x[2]-x[1])/2))],DRMIE4,lwd=1, lty=2, col=2)
          segments(x0=dgr[1], y0=log10(exp(a)), x1=max(dgr[dgr<=(x[2]-((x[2]-x[1])/2))]), y1=DRMIE4[1],lty=3,col=2)
        }
        if(prior.weights[11]>0){
          parH4 = H4resLND[4:7]
          a=exp(parH4[1])
          DRMH4 = log10(exp(DRM.H4_LND(par=parH4,x=10^(dgr[dgr>(x[2]-((x[2]-x[1])/2))]),q=data.modstan1a$q)))
          if(data.LN$data.modstan1a$shift==T) {DRMH4 = log10(exp(DRM.H4_LND(par=parH4,x=10^(dgr[dgr>(x[2]-((x[2]-x[1])/2))]),q=data.modstan1a$q))) + log10(exp(1.5*min(data.LN$data.modstan1a$m.org)))
          a = exp(parH4[1]) + 1.5*min(data.LN$data.modstan1a$m.org)}
          lines(dgr[dgr>(x[2]-((x[2]-x[1])/2))],DRMH4,lwd=1, lty=2, col=3)
          segments(x0=dgr[1], y0=log10(exp(a)), x1=max(dgr[dgr<=(x[2]-((x[2]-x[1])/2))]), y1=DRMH4[1],lty=3,col=3)
        }
        if(prior.weights[12]>0){
          parLN4 = LN4resLND[4:7]
          a=exp(parLN4[1])
          DRMLN4 = log10(exp(DRM.LN4_LND(par=parLN4,x=10^(dgr[dgr>(x[2]-((x[2]-x[1])/2))]),q=data.modstan1a$q)))
          if(data.LN$data.modstan1a$shift==T) {DRMLN4 = log10(exp(DRM.LN4_LND(par=parLN4,x=10^(dgr[dgr>(x[2]-((x[2]-x[1])/2))]),q=data.modstan1a$q))) + log10(exp(1.5*min(data.LN$data.modstan1a$m.org)))
          a = exp(parLN4[1]) + 1.5*min(data.LN$data.modstan1a$m.org)}
          lines(dgr[dgr>(x[2]-((x[2]-x[1])/2))],DRMLN4,lwd=1, lty=2, col=4)
          segments(x0=dgr[1], y0=log10(exp(a)), x1=max(dgr[dgr<=(x[2]-((x[2]-x[1])/2))]), y1=DRMLN4[1],lty=3,col=4)
        }
        if(prior.weights[13]>0){
          parG4 = G4resLND[4:7]
          a=exp(parG4[1])
          DRMG4 = log10(exp(DRM.G4_LND(par=parG4,x=10^(dgr[dgr>(x[2]-((x[2]-x[1])/2))]),q=data.modstan1a$q)))
          if(data.LN$data.modstan1a$shift==T) {DRMG4 = log10(exp(DRM.G4_LND(par=parG4,x=10^(dgr[dgr>(x[2]-((x[2]-x[1])/2))]),q=data.modstan1a$q))) + log10(exp(1.5*min(data.LN$data.modstan1a$m.org)))
          a = exp(parG4[1]) + 1.5*min(data.LN$data.modstan1a$m.org)}
          lines(dgr[dgr>(x[2]-((x[2]-x[1])/2))],DRMG4,lwd=1, lty=2, col=5)
          segments(x0=dgr[1], y0=log10(exp(a)), x1=max(dgr[dgr<=(x[2]-((x[2]-x[1])/2))]), y1=DRMG4[1],lty=3,col=5)
        }
        if(prior.weights[14]>0){
          parQE4 = QE4resLND[4:7]
          a=exp(parQE4[1])
          DRMQE4 = log10(exp(DRM.QE4_LND(par=parQE4,x=10^(dgr[dgr>(x[2]-((x[2]-x[1])/2))]),q=data.modstan1a$q)))
          if(data.LN$data.modstan1a$shift==T) {DRMQE4 = log10(exp(DRM.QE4_LND(par=parQE4,x=10^(dgr[dgr>(x[2]-((x[2]-x[1])/2))]),q=data.modstan1a$q))) + log10(exp(1.5*min(data.LN$data.modstan1a$m.org)))
          a = exp(parQE4[1]) + 1.5*min(data.LN$data.modstan1a$m.org)}
          lines(dgr[dgr>(x[2]-((x[2]-x[1])/2))],DRMQE4,lwd=1, lty=2, col=6)
          segments(x0=dgr[1], y0=log10(exp(a)), x1=max(dgr[dgr<=(x[2]-((x[2]-x[1])/2))]), y1=DRMQE4[1],lty=3,col=6)
        }
        if(prior.weights[15]>0){
          parP4 = P4resLND[4:7]
          a=exp(parP4[1])
          DRMP4 = log10(exp(DRM.P4_LND(par=parP4,x=10^(dgr[dgr>(x[2]-((x[2]-x[1])/2))]),q=data.modstan1a$q)))
          if(data.LN$data.modstan1a$shift==T) {DRMP4 = log10(exp(DRM.P4_LND(par=parP4,x=10^(dgr[dgr>(x[2]-((x[2]-x[1])/2))]),q=data.modstan1a$q))) + log10(exp(1.5*min(data.LN$data.modstan1a$m.org)))
          a = exp(parP4[1]) + 1.5*min(data.LN$data.modstan1a$m.org)}
          lines(dgr[dgr>(x[2]-((x[2]-x[1])/2))],DRMP4,lwd=1, lty=2, col=7)
          segments(x0=dgr[1], y0=log10(exp(a)), x1=max(dgr[dgr<=(x[2]-((x[2]-x[1])/2))]), y1=DRMP4[1],lty=3,col=7)
        }
        if(prior.weights[16]>0){
          parL4 = L4resLND[4:7]
          a=exp(parL4[1])
          DRML4 = log10(exp(DRM.L4_LND(par=parL4,x=10^(dgr[dgr>(x[2]-((x[2]-x[1])/2))]),q=data.modstan1a$q)))
          if(data.LN$data.modstan1a$shift==T) {DRML4 = log10(exp(DRM.L4_LND(par=parL4,x=10^(dgr[dgr>(x[2]-((x[2]-x[1])/2))]),q=data.modstan1a$q))) + log10(exp(1.5*min(data.LN$data.modstan1a$m.org)))
          a = exp(parL4[1]) + 1.5*min(data.LN$data.modstan1a$m.org)}
          lines(dgr[dgr>(x[2]-((x[2]-x[1])/2))],DRML4,lwd=1, lty=2, col=8)
          segments(x0=dgr[1], y0=log10(exp(a)), x1=max(dgr[dgr<=(x[2]-((x[2]-x[1])/2))]), y1=DRML4[1],lty=3,col=8)
        }
        leg.all = c(paste0("E4_N (w=",round(weight[1],4),")"),paste0("IE4_N (w=",round(weight[2],4),")"),paste0("H4_N (w=",round(weight[3],4),")"),
                    paste0("LN4_N (w=",round(weight[4],4),")"),paste0("G4_N (w=",round(weight[5],4),")"),paste0("QE4_N (w=",round(weight[6],4),")"),
                    paste0("P4_N (w=",round(weight[7],4),")"),paste0("L4_N (w=",round(weight[8],4),")"),
                    paste0("E4_LN (w=",round(weight[9],4),")"),paste0("IE4_LN (w=",round(weight[10],4),")"),paste0("H4_LN (w=",round(weight[11],4),")"),
                    paste0("LN4_LN (w=",round(weight[12],4),")"),paste0("G4_LN (w=",round(weight[13],4),")"),paste0("QE4_LN (w=",round(weight[14],4),")"),
                    paste0("P4_LN (w=",round(weight[15],4),")"),paste0("L4_LN (w=",round(weight[16],4),")"), "Model averaged BMDL")

        legend(x="topright", legend=leg.all[prior.weights>0], col=c(1:8,1:8,1)[prior.weights>0], lwd=c(rep(1,16))[prior.weights>0],
               lty=c(rep(1,8),rep(2,8),NA)[prior.weights>0], pch=c(rep(NA,16),3)[prior.weights>0],
               cex=0.7, title="Model-specific results (fit + BMDL)")
        lines(rbind(x,x,NA),rbind(m-ci,m+ci,NA),lty=3)
        lines(rbind(x,x,NA),rbind(m-ci2,m+ci2,NA))
        points(bmdl[1:16],y=rep(min(m-2*ci),length(bmdl)-1),col=c(1:8,1:8),pch=c(rep(2,8),rep(5,8)))
        points(bmdl[17],y=min(m-1.5*ci),pch=3)
        mtext("Bridge sampling")

      }

      print(plot2())

    }

    ## Covariances
    covs = t(data.frame(
      E4_N = E4covND,
      IE4_N = IE4covND,
      H4_N = H4covND,
      LN4_N = LN4covND,
      G4_N = G4covND,
      QE4_N = QE4covND,
      P4_N = P4covND,
      L4_N = L4covND,
      E4_LN = E4covLND,
      IE4_LN = IE4covLND,
      H4_LN = H4covLND,
      LN4_LN = LN4covLND,
      G4_LN = G4covLND,
      QE4_LN = QE4covLND,
      P4_LN = P4covLND,
      L4_LN = L4covLND
    ))
    colnames(covs) = c("b-d", "BMD-d")

    corrs = t(data.frame(
      E4_N = E4corrND,
      IE4_N = IE4corrND,
      H4_N = H4corrND,
      LN4_N = LN4corrND,
      G4_N = G4corrND,
      QE4_N = QE4corrND,
      P4_N = P4corrND,
      L4_N = L4corrND,
      E4_LN = E4corrLND,
      IE4_LN = IE4corrLND,
      H4_LN = H4corrLND,
      LN4_LN = LN4corrLND,
      G4_LN = G4corrLND,
      QE4_LN = QE4corrLND,
      P4_LN = P4corrLND,
      L4_LN = L4corrLND
    ))
    colnames(corrs) = c("b-d", "BMD-d")

    modelnames <- c("E4_N","IE4_N","H4_N","LN4_N","G4_N","QE4_N","P4_N","L4_N","E4_LN","IE4_LN","H4_LN","LN4_LN","G4_LN","QE4_LN","P4_LN","L4_LN")

    ret_results <- list(E4_N=E4outND,IE4_N=IE4outND,H4_N=H4outND,LN4_N=LN4outND,
                        G4_N=G4outND,QE4_N=QE4outND,P4_N=P4outND,L4_N=L4outND,
                        E4_LN=E4outLND,IE4_LN=IE4outLND,H4_LN=H4outLND,LN4_LN=LN4outLND,
                        G4_LN=G4outLND,QE4_LN=QE4outLND,P4_LN=P4outLND,L4_LN=L4outLND,
                        covs = covs, corrs = corrs,
                        weights_bridge_sampling=w.bs,
                        weights_laplace=lpwlp,
                        MA_bridge_sampling=macib,
                        MA_laplace=macilp, llN=llN, llLN=llLN,
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
                        BMDMixture = exp(mabmd)*data.modstan1a$maxD,
                        BMDMixture.conv = exp(mabmd.conv)*data.modstan1a$maxD,
                        divergences = divergences,
                        data = data.frame(
                          dose = c(data.N$data.modstan1a$x),
                          sd = sqrt(data.N$data.modstan1a$s2),
                          m = data.N$data.modstan1a$m),
                        max.dose = data.N$data.modstan1a$maxD,
                        q = data.N$data.modstan1a$q,
                        increasing = F,
                        models_included = modelnames[prior.weights > 0]
    )

    attr(ret_results, "class") <- c("BMADR", "BS")

    return(ret_results)


  }

}

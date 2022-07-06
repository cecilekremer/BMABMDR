
#' function to fit the BMD model using MCMC
#'
#' @param data.Q list containing the data, prior parameters and the starting values
#' @param prior.weights a vector of prior weights for the models. 1 imply a model should be used, 0 otherwise.
#'                      Defaults to rep(1, 8).
#' @param ndraws number of draws from the . Defaults to 30000.
#' @param nrchains number of MCMC chains to run. Defaults to 3.
#' @param nriterations total number of MCMC iterations. Defaults to 5000.
#' @param warmup total number of MCMC iterations to be used as warmup. Defaults to 1000.
#' @param delta adapt value for the MCMC chain. See \code{\link[rstan]{sampling}} for more.
#' @param treedepth adapt value for the MCMC chain. See \code{\link[rstan]{sampling}} for more.
#' @param seed random seed for reproducibility
#' @param pvec probability vector to compute credible interval for the BMD. Defaults to c(0.05,0.5,0.95).
#'
#' @examples
#'
#' @return list containing the following results
#' \enumerate{
#'   \item E4_Q parameter estimates from the exponential model
#'   \item IE4_Q parameter estimates from the inverse-exponential model
#'   \item H4_Q parameter estimates from the Hill model
#'   \item LN4_Q parameter estimates from the lognormal model
#'   \item G4_Q parameter estimates from the gamma model
#'   \item QE4_Q parameter estimates from the quadratic-exponential model
#'   \item P4_Q parameter estimates from the probit model
#'   \item L4_Q parameter estimates from the logit model
#'   \item MA_laplace Laplace approximation model averaged BMD estimates using all the models
#'   \item MA_bs Bridge sampling model averaged BMD estimates using all the models
#'   \item MA_bs_conv Bridge sampling model averaged BMD estimates using only converged models
#'   \item weights_laplace model weights using Laplace approximation to the posterior
#'   \item weights_bs model weights computed using bridge sampling
#'   \item convergence vector indicating model convergence or not. 1 = converged, 0 otherwise.
#'   \item llQ vector of model likelihoods
#'   \item bf Bayes factor comparing the best model against saturated ANOVA model
#' }
#'
#' @export samplingQ_MA
samplingQ_MA=function(data.Q,prior.weights = rep(1,8),
                      ndraws = 30000,nrchains=3,
                      nriterations=5000,warmup=1000,
                      delta=0.8,treedepth=10,seed=123,pvec=c(0.025,0.5,0.975)){

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
      warning('Difficulties fitting the Exponential model; prior weight was set to 0 and the model is not included in model averaging')
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
      warning('Difficulties fitting the Inverse Exponential model; prior weight was set to 0 and the model is not included in model averaging')
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
      warning('Difficulties fitting the Hill model; prior weight was set to 0 and the model is not included in model averaging')
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
      warning('Difficulties fitting the Lognormal model; prior weight was set to 0 and the model is not included in model averaging')
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

    data$init_b <- qgamma(data$q, rate=1.0, shape=E4resQ[6])/(E4resQ[2]/data$maxD)

    fitstanG4_Q <- fun_samplingQ(stanmodels$mG4_Q, data, start,
                                 ndraws,nrchains,
                                 nriterations,warmup,
                                 delta,treedepth,seed,pvec)
    if(is.null(fitstanG4_Q)){
      prior.weights[5] <- 0
      warning('Difficulties fitting the Gamma model; prior weight was set to 0 and the model is not included in model averaging')
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
      warning('Difficulties fitting the Quadratic Exponential model; prior weight was set to 0 and the model is not included in model averaging')
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
      warning('Difficulties fitting the Probit model; prior weight was set to 0 and the model is not included in model averaging')
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
      warning('Difficulties fitting the Logit model; prior weight was set to 0 and the model is not included in model averaging')
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
  mabmd=(c( # normal
    if("E4_Q" %in% names(count)) sample(as.matrix(fitstanE4_Q)[,2],count[names(count)=="E4_Q"],replace=T),
    if("IE4_Q" %in% names(count)) sample(as.matrix(fitstanIE4_Q)[,2],count[names(count)=="IE4_Q"],replace=T),
    if("H4_Q" %in% names(count)) sample(as.matrix(fitstanH4_Q)[,2],count[names(count)=="H4_Q"],replace=T),
    if("LN4_Q" %in% names(count)) sample(as.matrix(fitstanLN4_Q)[,2],count[names(count)=="LN4_Q"],replace=T),
    if("G4_Q" %in% names(count)) sample(as.matrix(fitstanG4_Q)[,2],count[names(count)=="G4_Q"],replace=T),
    if("QE4_Q" %in% names(count)) sample(as.matrix(fitstanQE4_Q)[,2],count[names(count)=="QE4_Q"],replace=T),
    if("P4_Q" %in% names(count)) sample(as.matrix(fitstanP4_Q)[,2],count[names(count)=="P4_Q"],replace=T),
    if("L4_Q" %in% names(count)) sample(as.matrix(fitstanL4_Q)[,2],count[names(count)=="L4_Q"],replace=T)
  ))
  macib=quantile(mabmd,pvec)*data$maxD
  names(macib)=c("BMDL","BMD","BMDU") # original scale

  BMDq_bs = quantile(mabmd, seq(0,1,0.005))*data$maxD

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
  if(0 %in% converged){
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
    mabmd.conv=(c( # normal
      if("E4_Q" %in% names(count)) sample(as.matrix(fitstanE4_Q)[,2],count[names(count)=="E4_Q"],replace=T),
      if("IE4_Q" %in% names(count)) sample(as.matrix(fitstanIE4_Q)[,2],count[names(count)=="IE4_Q"],replace=T),
      if("H4_Q" %in% names(count)) sample(as.matrix(fitstanH4_Q)[,2],count[names(count)=="H4_Q"],replace=T),
      if("LN4_Q" %in% names(count)) sample(as.matrix(fitstanLN4_Q)[,2],count[names(count)=="LN4_Q"],replace=T),
      if("G4_Q" %in% names(count)) sample(as.matrix(fitstanG4_Q)[,2],count[names(count)=="G4_Q"],replace=T),
      if("QE4_Q" %in% names(count)) sample(as.matrix(fitstanQE4_Q)[,2],count[names(count)=="QE4_Q"],replace=T),
      if("P4_Q" %in% names(count)) sample(as.matrix(fitstanP4_Q)[,2],count[names(count)=="P4_Q"],replace=T),
      if("L4_Q" %in% names(count)) sample(as.matrix(fitstanL4_Q)[,2],count[names(count)=="L4_Q"],replace=T)
    ))
    macib.conv <- quantile(mabmd.conv,pvec)*data$maxD
    names(macib.conv) <- c("BMDL","BMD","BMDU")

    BMDq_bs_conv <- quantile(mabmd.conv, seq(0,1,0.005))*data$maxD


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
    w.bs.conv = NULL; macib.conv = NULL; BMDq_bs_conv = NULL; dr.MA.bs.conv = NULL
  }

  #----------------------------------
  ### weights based on laplace approximation

  llQ = c() # likelihoods

  # getting the posterior modes and the hessian and the model specific posterior distributions
  if(prior.weights[1]>0){
    print("pw1")
    optE4_Q <- fun_optimQ(stanmodels$mE4_Q, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optE4_Q[[3]]),TRUE,(optE4_Q[[3]]!=0)) | length(optE4_Q)!=9)){
      prior.weights[1] <- 0
      warning('Difficulties fitting the Exponential model; prior weight was set to 0 and the model is not included in model averaging')
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
    print("pw2")
    optIE4_Q = fun_optimQ(stanmodels$mIE4_Q, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optIE4_Q[[3]]),TRUE,(optIE4_Q[[3]]!=0)) | length(optIE4_Q)!=9)){
      prior.weights[2] <- 0
      warning('Difficulties fitting the Inverse Exponential model; prior weight was set to 0 and the model is not included in model averaging')
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
    print("pw3")
    optH4_Q = fun_optimQ(stanmodels$mH4_Q, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optH4_Q[[3]]),TRUE,(optH4_Q[[3]]!=0)) | length(optH4_Q)!=9)){
      prior.weights[3] <- 0
      warning('Difficulties fitting the Hill model; prior weight was set to 0 and the model is not included in model averaging')
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
    print("pw4")
    optLN4_Q = fun_optimQ(stanmodels$mLN4_Q, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optLN4_Q[[3]]),TRUE,(optLN4_Q[[3]]!=0)) | length(optLN4_Q)!=9)){
      prior.weights[4] <- 0
      warning('Difficulties fitting the Lognormal model; prior weight was set to 0 and the model is not included in model averaging')
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
    print("pw5")

    data$init_b <- qgamma(data$q, rate=1.0, shape=optE4_Q$par[6])/optE4_Q$par[2]

    optG4_Q = fun_optimQ(stanmodels$mG4_Q, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optG4_Q[[3]]),TRUE,(optG4_Q[[3]]!=0)) | length(optG4_Q)!=9)){
      prior.weights[5] <- 0
      warning('Difficulties fitting the Gamma model; prior weight was set to 0 and the model is not included in model averaging')
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
    print("pw6")
    optQE4_Q = fun_optimQ(stanmodels$mQE4_Q, data, startQ, ndraws, 123, pvec)

    if((ifelse(is.na(optQE4_Q[[3]]),TRUE,(optQE4_Q[[3]]!=0)) | length(optQE4_Q)!=9)){
      prior.weights[6] <- 0
      warning('Difficulties fitting the Quadratic Exponential model; prior weight was set to 0 and the model is not included in model averaging')
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
    print("pw7")
    optP4_Q = fun_optimQ(stanmodels$mP4_Q, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optP4_Q[[3]]),TRUE,(optP4_Q[[3]]!=0)) | length(optP4_Q)!=9)){
      prior.weights[7] <- 0
      warning('Difficulties fitting the Probit model; prior weight was set to 0 and the model is not included in model averaging')
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
    print("pw8")
    optL4_Q = fun_optimQ(stanmodels$mL4_Q, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optL4_Q[[3]]),TRUE,(optL4_Q[[3]]!=0)) | length(optL4_Q)!=9)){
      prior.weights[8] <- 0
      warning('Difficulties fitting the Logit model; prior weight was set to 0 and the model is not included in model averaging')
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


  lpwlp=(prior.weights*w)/sum(prior.weights*w)
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
    print(warning('BMD/BMDL is larger than 20 for bridge sampling'))
  }
  if(macib[3]/macib[1] > 50){
    print(warning('BMDU/BMDL is larger than 50 for bridge sampling'))
  }
  if(macib[2] < (data.Q$data$x[2]*data.Q$data$maxD/10)){
    print(warning('BMD is 10 times lower than the lowest non-zero dose for bridge sampling'))
  }
  if(macilp[2]/macilp[1] > 20){
    print(warning('BMD/BMDL is larger than 20 for hybrid Laplace'))
  }
  if(macilp[3]/macilp[1] > 50){
    print(warning('BMDU/BMDL is larger than 50 for hybrid Laplace'))
  }
  if(macilp[2] < (data.Q$data$x[2]*data.Q$data$maxD/10)){
    print(warning('BMD is 10 times lower than the lowest non-zero dose for hybrid Laplace'))
  }

  ### best fitting model vs saturated ANOVA model
  #best.fit = modelnames[which(weight[1:8] == max(weight[1:8]))]

  #bfTest <- modelTestQ(best.fit, data.Q, get(paste0('fitstan', best.fit, '_Q')), type = 'MCMC',
  #                    seed, ndraws, nrchains, nriterations, warmup, delta, treedepth)
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
                      divergences = divergences,
                      data = data.frame(
                        dose = c(data.Q$data$x),
                        y = data.Q$data$y,
                        n = data.Q$data$n),
                      max.dose = data.Q$data$maxD,
                      q = data.Q$data$q,
                      models_included = model[prior.weights > 0],
                      is_bin = data.Q$data$is_bin,
                      is_betabin = data.Q$data$is_betabin#,
                      #bf = bfTest$bayesFactor,
                      #means.SM = bfTest$means.SM, parBestFit = bfTest$par.best,
                      #BIC.bestfit = bfTest$BIC.bestfit, BIC.SM = bfTest$BIC.SM
  )

  attr(ret_results, "class") <- c("BMADRQ", "BS")

  return(ret_results)

}


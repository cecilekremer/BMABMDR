#' Perform model averaging using Full Laplace method
#'
#' This method assumed data for continuous endpoints.
#'
#' More detailed descriprion
#' @param data.N the input data as returned by function PREP_DATA_N
#' @param data.LN the input data as returned by function PREP_DATA_LN
#' @param prior.weights a vector specifying which of the 16 models should be included (1 = include, 0 = exclude)
#' @param ndraws the number of draws, default 30000
#' @param nrchains the number of chains to be used in the MCMC
#' @param nriterations the number of iterations per chain
#' @param warmup the number of iterations per chain to be discarded as burnin
#' @param delta default 0.8
#' @param treedepth default 10
#' @param seed default 123
#' @param pvec vector specifying the three BMD quantiles of interest
#' @param plot logical indicating whether a simple plot of model fits should be shown (defaults to FALSE)
#'
#' @return List with the model-specific results, model weights, the model averaged BMD, BMDL and BMDU. In addition for each model a matrix with covariance/correlation between parameter b/BMD and parameter d is given. Some additional output used in other external functions is also given.
#'
#' @export
#'
full.laplaceQ_MA=function(data.Q, prior.weights,ndraws=30000,nrchains=3,
                          nriterations=3000,warmup=1000,
                          delta=0.8,treedepth=10,seed=123,pvec=c(0.05,0.5,0.95)){

  # prior.weights = prior.weights/sum(prior.weights==1) #this is now done below when calculating weights

  ## Data to use for Normal distribution
  data = data.Q$data
  start = data.Q$start

  llQ = c() # likelihoods
  BMDL = c()
  BMD = c()
  BMDU = c()

  ### Getting the posterior modes and the hessian to obtain the model specific posterior distributions
  # optimizing function: obtain point estimates by maximizing the joint posterior from the Stan model
  ### Additionally save BMD(L/U) estimates and loglikelihood values; set to NA if model has prior weight 0


  if(prior.weights[1]>0){

    print(1)
    optE4_Q <- fun_optimQ(stanmodels$mE4_Q, data, start, ndraws, 123, pvec)

    if(data$is_bin == 1) {
      parsE4Q <- parq_extract(optE4_Q, model_name = "E4_Q", pars = c('a', 'b', 'd', 'BMD',
                                                                     paste0('par',1:3)))
    } else {
      parsE4Q <- parq_extract(optE4_Q, model_name = "E4_Q", pars = c('a', 'b', 'd', 'BMD','rho[1]',
                                                                     paste0('par',1:3)),
                              rho = TRUE)
    }

    E4resQ <- quantile(parsE4Q$BMD, pvec)*data$maxD
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
    E4covQ <- c(cov(parsE4Q[,c("b","d")], use="complete.obs")["b","d"],
                cov(parsE4Q[,c("BMD","d")], use="complete.obs")["BMD","d"])

    E4corrQ <- c(cor(parsE4Q[,c("b","d")], use="complete.obs")["b","d"],
                 cor(parsE4Q[,c("BMD","d")], use="complete.obs")["BMD","d"])

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
  }else{E4resQ=NULL; llQ = c(llQ,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA);
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

  if(prior.weights[2]>0){

    print(2)
    optIE4_Q <- fun_optimQ(stanmodels$mIE4_Q, data, start, ndraws, 123, pvec)

    if(data$is_bin == 1) {
      parsIE4Q <- parq_extract(optIE4_Q, model_name = "IE4_Q", pars = c('a', 'b', 'd', 'BMD',
                                                                        paste0('par',1:3)))
    } else {
      parsIE4Q <- parq_extract(optIE4_Q, model_name = "IE4_Q", pars = c('a', 'b', 'd', 'BMD','rho[1]',
                                                                        paste0('par',1:3)),
                               rho = TRUE)
    }
    IE4resQ <- quantile(parsIE4Q$BMD, pvec)*data$maxD  #exp(quantile(optE4_Q$theta_tilde[,"par[2]"],pvec))

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
    IE4covQ <- c(cov(parsIE4Q[,c("b","d")], use="complete.obs")["b","d"],
                 cov(parsIE4Q[,c("BMD","d")], use="complete.obs")["BMD","d"])

    IE4corrQ <- c(cor(parsIE4Q[,c("b","d")], use="complete.obs")["b","d"],
                  cor(parsIE4Q[,c("BMD","d")], use="complete.obs")["BMD","d"])

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
  }else{
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

  if(prior.weights[3]>0){

    print(3)

    optH4_Q <- fun_optimQ(stanmodels$mH4_Q, data, start, ndraws, 123, pvec)

    if(data$is_bin == 1) {
      parsH4Q <- parq_extract(optH4_Q, model_name = "H4_Q", pars = c('a', 'b', 'd', 'BMD',
                                                                     paste0('par',1:3)))
    } else {
      parsH4Q <- parq_extract(optH4_Q, model_name = "H4_Q", pars = c('a', 'b', 'd', 'BMD','rho[1]',
                                                                     paste0('par',1:3)),
                              rho = TRUE)
    }

    H4resQ <- quantile(parsH4Q$BMD, pvec)*data$maxD
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
    H4covQ = c(cov(parsH4Q[,c("b","d")], use="complete.obs")["b","d"],
               cov(parsH4Q[,c("BMD","d")], use="complete.obs")["BMD","d"])

    H4corrQ = c(cor(parsH4Q[,c("b","d")], use="complete.obs")["b","d"],
                cor(parsH4Q[,c("BMD","d")], use="complete.obs")["BMD","d"])

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
  }else{H4resQ=NULL; llQ = c(llQ,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA);
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

  if(prior.weights[4]>0){

    print(4)
    optLN4_Q <- fun_optimQ(stanmodels$mLN4_Q, data, start, ndraws, 123, pvec)

    if(data$is_bin == 1) {
      parsLN4Q <- parq_extract(optLN4_Q, model_name = "LN4_Q", pars = c('a', 'b', 'd', 'BMD',
                                                                        paste0('par',1:3)))
    } else {
      parsLN4Q <- parq_extract(optLN4_Q, model_name = "LN4_Q", pars = c('a', 'b', 'd', 'BMD','rho[1]',
                                                                        paste0('par',1:3)),
                               rho = TRUE)
    }
    LN4resQ <- quantile(parsLN4Q$BMD, pvec)*data$maxD  #exp(quantile(optE4_Q$theta_tilde[,"par[2]"],pvec))

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
    LN4covQ = c(cov(parsLN4Q[,c("b","d")], use="complete.obs")["b","d"],
                cov(parsLN4Q[,c("BMD","d")], use="complete.obs")["BMD","d"])

    LN4corrQ = c(cor(parsLN4Q[,c("b","d")], use="complete.obs")["b","d"],
                 cor(parsLN4Q[,c("BMD","d")], use="complete.obs")["BMD","d"])

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
  }else{LN4resQ=NULL; llQ = c(llQ,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA);
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

  if(prior.weights[5]>0){
    print(5)
    # optG4_Q = fun_optimQ_G4(modN, stanmodels$mG4_Q, data, stv = start, ndraws, 123, pvec)
    optG4_Q <- fun_optimQ(stanmodels$mG4_Q, data, start, ndraws, 123, pvec)


    if(data$is_bin == 1) {
      parsG4Q <- parq_extract(optG4_Q, model_name = "G4_Q", pars = c('a', 'b', 'd', 'BMD',
                                                                     paste0('par',1:3)))
    } else {
      parsG4Q <- parq_extract(optG4_Q, model_name = "G4_Q", pars = c('a', 'b', 'd', 'BMD','rho[1]',
                                                                     paste0('par',1:3)),
                              rho = TRUE)
    }

    G4resQ <- quantile(parsG4Q$BMD, pvec)*data$maxD
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
    G4covQ = c(cov(parsG4Q[,c("b","d")], use="complete.obs")["b","d"],
               cov(parsG4Q[,c("BMD","d")], use="complete.obs")["BMD","d"])

    G4corrQ = c(cor(parsG4Q[,c("b","d")], use="complete.obs")["b","d"],
                cor(parsG4Q[,c("BMD","d")], use="complete.obs")["BMD","d"])

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
  }else{G4resQ=NULL; llQ = c(llQ,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA);
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

  if(prior.weights[6]>0){
    print(6)
    optQE4_Q = fun_optimQ(stanmodels$mQE4_Q, data, start, ndraws, 123, pvec)

    if(data$is_bin == 1) {
      parsQE4Q <- parq_extract(optQE4_Q, model_name = "QE4_Q", pars = c('a', 'b', 'd', 'BMD',
                                                                        paste0('par',1:3)))
    } else {
      parsQE4Q <- parq_extract(optQE4_Q, model_name = "QE4_Q", pars = c('a', 'b', 'd', 'BMD','rho[1]',
                                                                        paste0('par',1:3)),
                               rho = TRUE)
    }
    QE4resQ <- quantile(parsQE4Q$BMD, pvec)*data$maxD  #exp(quantile(optE4_Q$theta_tilde[,"par[2]"],pvec))

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
    QE4covQ = c(cov(parsQE4Q[,c("b","d")], use="complete.obs")["b","d"],
                cov(parsQE4Q[,c("BMD","d")], use="complete.obs")["BMD","d"])

    QE4corrQ = c(cor(parsQE4Q[,c("b","d")], use="complete.obs")["b","d"],
                 cor(parsQE4Q[,c("BMD","d")], use="complete.obs")["BMD","d"])

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
  }else{QE4resQ=NULL; llQ = c(llQ,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA);
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

  if(prior.weights[7]>0){
    print(7)
    optP4_Q = fun_optimQ(stanmodels$mP4_Q, data, start, ndraws, 123, pvec)

    if(data$is_bin == 1) {
      parsP4Q <- parq_extract(optP4_Q, model_name = "P4_Q", pars = c('a', 'b', 'd', 'BMD',
                                                                     paste0('par',1:3)))
    } else {
      parsP4Q <- parq_extract(optP4_Q, model_name = "P4_Q", pars = c('a', 'b', 'd', 'BMD','rho[1]',
                                                                     paste0('par',1:3)),
                              rho = TRUE)
    }

    P4resQ <- quantile(parsP4Q$BMD, pvec)*data$maxD
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
    P4covQ = c(cov(parsP4Q[,c("b","d")], use="complete.obs")["b","d"],
               cov(parsP4Q[,c("BMD","d")], use="complete.obs")["BMD","d"])

    P4corrQ = c(cor(parsP4Q[,c("b","d")], use="complete.obs")["b","d"],
                cor(parsP4Q[,c("BMD","d")], use="complete.obs")["BMD","d"])

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
  }else{P4resQ=NULL; llQ = c(llQ,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA);
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

  if(prior.weights[8]>0){
    print(8)
    optL4_Q = fun_optimQ(stanmodels$mL4_Q, data, start, ndraws, 123, pvec)

    if(data$is_bin == 1) {
      parsL4Q <- parq_extract(optL4_Q, model_name = "L4_Q", pars = c('a', 'b', 'd', 'BMD',
                                                                     paste0('par',1:3)))
    } else {
      parsL4Q <- parq_extract(optL4_Q, model_name = "L4_Q", pars = c('a', 'b', 'd', 'BMD','rho[1]',
                                                                     paste0('par',1:3)),
                              rho = TRUE)
    }

    L4resQ <- quantile(parsL4Q$BMD, pvec)*data$maxD
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
    L4covQ = c(cov(parsL4Q[,c("b","d")], use="complete.obs")["b","d"],
               cov(parsL4Q[,c("BMD","d")], use="complete.obs")["BMD","d"])

    L4corrQ = c(cor(parsL4Q[,c("b","d")], use="complete.obs")["b","d"],
                cor(parsL4Q[,c("BMD","d")], use="complete.obs")["BMD","d"])

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
  }else{L4resQ=NULL; llQ = c(llQ,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA);
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
  w=c()

  # normal


  if(prior.weights[1]>0){
    DIHE4h=det(-solve(optE4_Q$hessian))
    DIHE4=ifelse(DIHE4h<0,0,DIHE4h)
    ## Aproximation of marginal (i.e. integrated) likelihood (= 'model evidence')
    if(data$is_bin==1) {
      w=c(w,(2*pi)^(2.5)*sqrt(DIHE4)*exp(llE4Q-minll)*
            dnorm(optE4_Q$par[3],mean=data$priormu[3],
                  sd=data$priorSigma[3,3])*
            mc2d::dpert(optE4_Q$par[2], min = data$priorlb[2], max = data$priorub[2],
                        mode = data$priormu[2], shape = data$priorgama[2])*
            mc2d::dpert(optE4_Q$par[1], min = data$priorlb[1], max = data$priorub[1],
                        mode = data$priormu[1], shape = data$priorgama[1])
      )
    } else if(data$is_betabin == 1){
      w=c(w,(2*pi)^(2.5)*sqrt(DIHE4)*exp(llE4Q-minll)*
            dnorm(optE4_Q$par[3],mean=data$priormu[3],
                  sd=data$priorSigma[3,3])*
            dnorm(optE4_Q$par[stringr::str_detect(names(optIE4_Q$par),'etarho')],
                  mean=data$priormu[3],
                  sd=data$priorSigma[3,3])*
            mc2d::dpert(optE4_Q$par[2], min = data$priorlb[2], max = data$priorub[2],
                        mode = data$priormu[2], shape = data$priorgama[2])*
            mc2d::dpert(optE4_Q$par[1], min = data$priorlb[1], max = data$priorub[1],
                        mode = data$priormu[1], shape = data$priorgama[1])
      )
    }

  }else{w=c(w,0)}

  if(prior.weights[2]>0){
    DIHIE4h=det(-solve(optIE4_Q$hessian))
    DIHIE4=ifelse(DIHIE4h<0,0,DIHIE4h)

    if(data$is_bin==1) {
      w=c(w,(2*pi)^(2.5)*sqrt(DIHIE4)*exp(llIE4Q-minll)*
            dnorm(optIE4_Q$par[3],mean=data$priormu[3],
                  sd=data$priorSigma[3,3])*
            mc2d::dpert(optIE4_Q$par[2], min = data$priorlb[2], max = data$priorub[2],
                        mode = data$priormu[2], shape = data$priorgama[2])*
            mc2d::dpert(optIE4_Q$par[1], min = data$priorlb[1], max = data$priorub[1],
                        mode = data$priormu[1], shape = data$priorgama[1])
      )
    } else if(data$is_betabin == 1){
      w=c(w,(2*pi)^(2.5)*sqrt(DIHIE4)*exp(llIE4Q-minll)*
            dnorm(optIE4_Q$par[3],mean=data$priormu[3],
                  sd=data$priorSigma[3,3])*
            dnorm(optIE4_Q$par[stringr::str_detect(names(optIE4_Q$par),'etarho')],
                  mean=data$priormu[3],
                  sd=data$priorSigma[3,3])*
            mc2d::dpert(optIE4_Q$par[2], min = data$priorlb[2], max = data$priorub[2],
                        mode = data$priormu[2], shape = data$priorgama[2])*
            mc2d::dpert(optIE4_Q$par[1], min = data$priorlb[1], max = data$priorub[1],
                        mode = data$priormu[1], shape = data$priorgama[1])
      )
    }

  }else{w=c(w,0)}

  if(prior.weights[3]>0){
    DIHH4h=det(-solve(optH4_Q$hessian))
    DIHH4=ifelse(DIHH4h<0,0,DIHH4h)

    if(data$is_bin==1) {
      w=c(w,(2*pi)^(2.5)*sqrt(DIHH4)*exp(llH4Q-minll)*
            dnorm(optH4_Q$par[3],mean=data$priormu[3],
                  sd=data$priorSigma[3,3])*
            mc2d::dpert(optH4_Q$par[2], min = data$priorlb[2], max = data$priorub[2],
                        mode = data$priormu[2], shape = data$priorgama[2])*
            mc2d::dpert(optH4_Q$par[1], min = data$priorlb[1], max = data$priorub[1],
                        mode = data$priormu[1], shape = data$priorgama[1])
      )
    } else if(data$is_betabin == 1){
      w=c(w,(2*pi)^(2.5)*sqrt(DIHH4)*exp(llH4Q-minll)*
            dnorm(optH4_Q$par[3],mean=data$priormu[3],
                  sd=data$priorSigma[3,3])*
            dnorm(optH4_Q$par[stringr::str_detect(names(optH4_Q$par),'etarho')],
                  mean=data$priormu[3],
                  sd=data$priorSigma[3,3])*
            mc2d::dpert(optH4_Q$par[2], min = data$priorlb[2], max = data$priorub[2],
                        mode = data$priormu[2], shape = data$priorgama[2])*
            mc2d::dpert(optH4_Q$par[1], min = data$priorlb[1], max = data$priorub[1],
                        mode = data$priormu[1], shape = data$priorgama[1])
      )
    }

  }else{w=c(w,0)}

  if(prior.weights[4]>0){
    DIHLN4h=det(-solve(optLN4_Q$hessian))
    DIHLN4=ifelse(DIHLN4h<0,0,DIHLN4h)

    if(data$is_bin==1) {
      w=c(w,(2*pi)^(2.5)*sqrt(DIHLN4)*exp(llLN4Q-minll)*
            dnorm(optLN4_Q$par[3],mean=data$priormu[3],
                  sd=data$priorSigma[3,3])*
            mc2d::dpert(optLN4_Q$par[2], min = data$priorlb[2], max = data$priorub[2],
                        mode = data$priormu[2], shape = data$priorgama[2])*
            mc2d::dpert(optLN4_Q$par[1], min = data$priorlb[1], max = data$priorub[1],
                        mode = data$priormu[1], shape = data$priorgama[1])
      )
    } else if(data$is_betabin == 1){
      w=c(w,(2*pi)^(2.5)*sqrt(DIHLN4)*exp(llLN4Q-minll)*
            dnorm(optLN4_Q$par[3],mean=data$priormu[3],
                  sd=data$priorSigma[3,3])*
            dnorm(optLN4_Q$par[stringr::str_detect(names(optLN4_Q$par),'etarho')],
                  mean=data$priormu[3],
                  sd=data$priorSigma[3,3])*
            mc2d::dpert(optLN4_Q$par[2], min = data$priorlb[2], max = data$priorub[2],
                        mode = data$priormu[2], shape = data$priorgama[2])*
            mc2d::dpert(optLN4_Q$par[1], min = data$priorlb[1], max = data$priorub[1],
                        mode = data$priormu[1], shape = data$priorgama[1])
      )
    }

  }else{w=c(w,0)}

  if(prior.weights[5]>0){
    DIHG4h=det(-solve(optG4_Q$hessian))
    DIHG4=ifelse(DIHG4h<0,0,DIHG4h)

    if(data$is_bin==1) {
      w=c(w,(2*pi)^(2.5)*sqrt(DIHG4)*exp(llG4Q-minll)*
            dnorm(optG4_Q$par[3], mean=data$priormu[3],
                  sd=data$priorSigma[3,3])*
            mc2d::dpert(optG4_Q$par[2], min = data$priorlb[2], max = data$priorub[2],
                        mode = data$priormu[2], shape = data$priorgama[2])*
            mc2d::dpert(optG4_Q$par[1], min = data$priorlb[1], max = data$priorub[1],
                        mode = data$priormu[1], shape = data$priorgama[1])
      )
    } else if(data$is_betabin == 1){
      w=c(w,(2*pi)^(2.5)*sqrt(DIHG4)*exp(llG4Q-minll)*
            dnorm(optG4_Q$par[3], mean=data$priormu[3],
                  sd=data$priorSigma[3,3])*
            dnorm(optG4_Q$par[stringr::str_detect(names(optG4_Q$par),'etarho')],
                  mean=data$priormu[3],
                  sd=data$priorSigma[3,3])*
            mc2d::dpert(optG4_Q$par[2], min = data$priorlb[2], max = data$priorub[2],
                        mode = data$priormu[2], shape = data$priorgama[2])*
            mc2d::dpert(optG4_Q$par[1], min = data$priorlb[1], max = data$priorub[1],
                        mode = data$priormu[1], shape = data$priorgama[1])
      )
    }

  }else{w=c(w,0)}

  if(prior.weights[6]>0){
    DIHQE4h=det(-solve(optQE4_Q$hessian))
    DIHQE4=ifelse(DIHQE4h<0,0,DIHQE4h)

    if(data$is_bin==1) {
      w=c(w,(2*pi)^(2.5)*sqrt(DIHQE4)*exp(llQE4Q-minll)*
            dnorm(optQE4_Q$par[3],mean=data$priormu[3],
                  sd=data$priorSigma[3,3])*
            mc2d::dpert(optQE4_Q$par[2], min = data$priorlb[2], max = data$priorub[2],
                        mode = data$priormu[2], shape = data$priorgama[2])*
            mc2d::dpert(optQE4_Q$par[1], min = data$priorlb[1], max = data$priorub[1],
                        mode = data$priormu[1], shape = data$priorgama[1])
      )
    } else if(data$is_betabin == 1) {
      w=c(w,(2*pi)^(2.5)*sqrt(DIHQE4)*exp(llQE4Q-minll)*
            dnorm(optQE4_Q$par[3],mean=data$priormu[3],
                  sd=data$priorSigma[3,3])*
            dnorm(optQE4_Q$par[stringr::str_detect(names(optQE4_Q$par),'etarho')],
                  mean=data$priormu[3],
                  sd=data$priorSigma[3,3])*
            mc2d::dpert(optQE4_Q$par[2], min = data$priorlb[2], max = data$priorub[2],
                        mode = data$priormu[2], shape = data$priorgama[2])*
            mc2d::dpert(optQE4_Q$par[1], min = data$priorlb[1], max = data$priorub[1],
                        mode = data$priormu[1], shape = data$priorgama[1])
      )
    }

  }else{w=c(w,0)}

  if(prior.weights[7]>0){
    DIHP4h=det(-solve(optP4_Q$hessian))
    DIHP4=ifelse(DIHP4h<0,0,DIHP4h)

    if(data$is_bin==1) {
      w=c(w,(2*pi)^(2.5)*sqrt(DIHP4)*exp(llP4Q-minll)*
            dnorm(optP4_Q$par[3],mean=data$priormu[3],
                  sd=data$priorSigma[3,3])*
            mc2d::dpert(optP4_Q$par[2], min = data$priorlb[2], max = data$priorub[2],
                        mode = data$priormu[2], shape = data$priorgama[2])*
            mc2d::dpert(optP4_Q$par[1], min = data$priorlb[1], max = data$priorub[1],
                        mode = data$priormu[1], shape = data$priorgama[1])
      )
    } else if(data$is_betabin == 1) {
      w=c(w,(2*pi)^(2.5)*sqrt(DIHP4)*exp(llP4Q-minll)*
            dnorm(optP4_Q$par[3],mean=data$priormu[3],
                  sd=data$priorSigma[3,3])*
            dnorm(optP4_Q$par[stringr::str_detect(names(optP4_Q$par),'etarho')],
                  mean=data$priormu[3],
                  sd=data$priorSigma[3,3])*
            mc2d::dpert(optP4_Q$par[2], min = data$priorlb[2], max = data$priorub[2],
                        mode = data$priormu[2], shape = data$priorgama[2])*
            mc2d::dpert(optP4_Q$par[1], min = data$priorlb[1], max = data$priorub[1],
                        mode = data$priormu[1], shape = data$priorgama[1])
      )
    }

  }else{w=c(w,0)}

  if(prior.weights[8]>0){
    DIHL4h=det(-solve(optL4_Q$hessian))
    DIHL4=ifelse(DIHL4h<0,0,DIHL4h)

    if(data$is_bin==1) {
      w=c(w,(2*pi)^(2.5)*sqrt(DIHL4)*exp(llL4Q-minll)*
            dnorm(optL4_Q$par[3],mean=data$priormu[3],
                  sd=data$priorSigma[3,3])*
            mc2d::dpert(optL4_Q$par[2], min = data$priorlb[2], max = data$priorub[2],
                        mode = data$priormu[2], shape = data$priorgama[2])*
            mc2d::dpert(optL4_Q$par[1], min = data$priorlb[1], max = data$priorub[1],
                        mode = data$priormu[1], shape = data$priorgama[1])
      )
    } else if(data$is_betabin == 1) {
      w=c(w,(2*pi)^(2.5)*sqrt(DIHL4)*exp(llL4Q-minll)*
            dnorm(optL4_Q$par[3],mean=data$priormu[3],
                  sd=data$priorSigma[3,3])*
            dnorm(optL4_Q$par[stringr::str_detect(names(optL4_Q$par),'etarho')],
                  mean=data$priormu[3],
                  sd=data$priorSigma[3,3])*
            mc2d::dpert(optL4_Q$par[2], min = data$priorlb[2], max = data$priorub[2],
                        mode = data$priormu[2], shape = data$priorgama[2])*
            mc2d::dpert(optL4_Q$par[1], min = data$priorlb[1], max = data$priorub[1],
                        mode = data$priormu[1], shape = data$priorgama[1])
      )
    }

  }else{w=c(w,0)}


  prior.weights = prior.weights/sum(prior.weights==1) # in case this has changed for G4
  lpw=(prior.weights*w)/sum(prior.weights*w)

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
  maci=quantile(mabmd,pvec)*data$maxD ## original scale
  names(maci)=c("BMDL","BMD","BMDU")

  BMDq = quantile(mabmd, seq(0,1,0.005))*data$maxD ## original scale

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
    print(warning('BMD/BMDL is larger than 20'))
  }
  if(maci[3]/maci[1] > 50){
    print(warning('BMDU/BMDL is larger than 50'))
  }
  if(maci[2] < (data.Q$data$x[2]*data.Q$data$maxD/10)){
    print(warning('BMD is 10 times lower than the lowest non-zero dose'))
  }

  ### best fitting model vs saturated ANOVA model
  best.fit = modelnames[which(weight[1:8] == max(weight[1:8]))]

  bfTest <- modelTestQ(best.fit, data.Q, get(paste0('opt', best.fit, '_Q')), type = 'Laplace',
                      seed, ndraws, nrchains, nriterations, warmup, delta, treedepth)
  print(warning(bfTest$warn.bf))

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
    q = data.Q$data$q,
    is_bin = data.Q$data$is_bin,
    is_betabin = data.Q$data$is_betabin,
    bf = bfTest$bayesFactor,
    means.SM = bfTest$means.SM, parBestFit = bfTest$par.best,
    BIC.bestfit = bfTest$BIC.bestfit, BIC.SM = bfTest$BIC.SM
  )

  attr(ret_results, "class") <- c("BMADRQ", "LP")

  return(ret_results)
}

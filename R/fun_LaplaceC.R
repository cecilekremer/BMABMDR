# Full laplace approximation
##################################################################################################

#' function to perform model-averaging using Laplace approximation
#'
#' @param data.N list containing data values for the normal models
#' @param data.LN list containing data values for the lognormal models
#' @param prior.weights model weights determining if a model is to be included in the model averaging or not.
#'                      1 implies model is included, 0 otherwise. Defaults to rep(1,16).
#' @param ndraws number of draws to be made from the posterior distribution. Defaults to 30000
#' @param seed random seed for reproducibility. Defaults to 123
#' @param pvec probability vector to compute credible interval for the BMD. Defaults to c(0.05,0.5,0.95).
#' @param plot logical variable to determine if the results should be plotted. Defaults to FALSE
#'
#' @examples
#'
#' @return list containing the following results
#' \enumerate{
#'   \item E4_N parameter estimates from the exponential model with normal distribution assumed
#'   \item IE4_N parameter estimates from the inverse-exponential model with normal distribution assumed
#'   \item H4_N parameter estimates from the Hill model with normal distribution assumed
#'   \item LN4_N parameter estimates from the lognormal model with normal distribution assumed
#'   \item G4_N parameter estimates from the gamma model with normal distribution assumed
#'   \item QE4_N parameter estimates from the quadratic-exponential model with normal distribution assumed
#'   \item P4_N parameter estimates from the probit model with normal distribution assumed
#'   \item L4_N parameter estimates from the logit model with normal distribution assumed
#'   \item E4_LN parameter estimates from the exponential model with lognormal distribution assumed
#'   \item IE4_LN parameter estimates from the inverse exponential model with lognormal distribution assumed
#'   \item H4_LN parameter estimates from the Hill model with lognormal distribution assumed
#'   \item LN4_LN parameter estimates from the lognormal model with lognormal distribution assumed
#'   \item G4_LN parameter estimates from the gamma model with lognormal distribution assumed
#'   \item QE4_LN parameter estimates from the quadratic exponential model with lognormal distribution assumed
#'   \item P4_LN parameter estimates from the probit model with lognormal distribution assumed
#'   \item L4_LN parameter estimates from the logit model with lognormal distribution assumed
#'   \item MA_laplace Laplace approximation model averaged BMD estimates using all the models
#'   \item MA_lp_conv Laplace approximation model averaged BMD estimates using only converged models
#'   \item weights_laplace model weights using Laplace approximation to the posterior
#'   \item convergence vector indicating model convergence or not. 1 = converged, 0 otherwise.
#'   \item llN vector of model likelihoods for normal distribution
#'   \item llLN vector of model likelihoods for lognonormal distribution
#'   \item bf Bayes factor comparing the best model against saturated ANOVA model
#' }
#'
#' @export full.laplace_MAc
#'
full.laplace_MAc=function(data.N, data.LN,
                          prior.weights = rep(1,16),
                          ndraws=30000,seed=123,
                          pvec=c(0.05,0.5,0.95)){

  # prior.weights = prior.weights/sum(prior.weights==1) #this is now done below when calculating weights


  data = data.N$data
  start = data.N$start
  startQ = data.N$startQ

  llN = c() # likelihoods
  BMDL = c()
  BMD = c()
  BMDU = c()

  nModels = 16

  if (!is.null(shiny::getDefaultReactiveDomain())){
    shiny::incProgress(amount = 1/(nModels + 1), detail = names(get_models(type='continuous')[1]))
  }
  if(prior.weights[1]>0){

    # print(1)

    optE4_NI <- fun_optimC(stanmodels$mE4c, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optE4_NI[[3]]),TRUE,(optE4_NI[[3]]!=0)) | length(optE4_NI)!=9)){
      prior.weights[1] <- 0
      warning('Difficulties fitting the Exponential (Normal) model; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsE4N <- par_extractC(optE4_NI, model_name = "E4_N")
      E4resNI <- quantile(parsE4N$BMD*data$maxD, pvec)
      E4resNI <- c(E4resNI,apply(parsE4N[,c(paste0("p",1:4), "is2t", "rho")], 2, median))
      names(E4resNI) <- c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t","rho")

      E4outNI <- outLPC(parsE4N, pvec, data$maxD)

      # Covariance between b-d and between BMD-d
      E4covNI <- c(cov(parsE4N[,c("b","d")], use="complete.obs")["b","d"],
                   cov(parsE4N[,c("BMD","d")], use="complete.obs")["BMD","d"])

      E4corrNI <- c(cor(parsE4N[,c("b","d")], use="complete.obs")["b","d"],
                    cor(parsE4N[,c("BMD","d")], use="complete.obs")["BMD","d"])

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
      warning('Difficulties fitting the Inverse Exponential (Normal) model; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsIE4N <- par_extractC(optIE4_NI, model_name = "IE4_N")
      IE4resNI <- quantile(parsIE4N$BMD*data$maxD, pvec)
      IE4resNI <- c(IE4resNI,apply(parsIE4N[,c(paste0("p",1:4), "is2t","rho")], 2, median))
      names(IE4resNI) <- c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t","rho")

      IE4outNI <- outLPC(parsIE4N, pvec, data$maxD)

      # Covariance between b-d and between BMD-d
      IE4covNI <- c(cov(parsIE4N[,c("b","d")], use="complete.obs")["b","d"],
                    cov(parsIE4N[,c("BMD","d")], use="complete.obs")["BMD","d"])

      IE4corrNI <- c(cor(parsIE4N[,c("b","d")], use="complete.obs")["b","d"],
                     cor(parsIE4N[,c("BMD","d")], use="complete.obs")["BMD","d"])

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
      warning('Difficulties fitting the Hill (Normal) model; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsH4N <- par_extractC(optH4_NI, model_name = "H4_N")
      H4resNI <- quantile(parsH4N$BMD*data$maxD, pvec)  #exp(quantile(optE4_NI$theta_tilde[,"par[2]"],pvec))
      H4resNI <- c(H4resNI,apply(parsH4N[,c(paste0("p",1:4), "is2t","rho")], 2, median))
      names(H4resNI) <- c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t","rho")

      H4outNI <- outLPC(parsH4N, pvec, data$maxD)

      # Covariance between b-d and between BMD-d
      H4covNI = c(cov(parsH4N[,c("b","d")], use="complete.obs")["b","d"],
                  cov(parsH4N[,c("BMD","d")], use="complete.obs")["BMD","d"])

      H4corrNI = c(cor(parsH4N[,c("b","d")], use="complete.obs")["b","d"],
                   cor(parsH4N[,c("BMD","d")], use="complete.obs")["BMD","d"])

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
      warning('Difficulties fitting the Lognormal (Normal) model; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsLN4N <- par_extractC(optLN4_NI, model_name = "LN4_N")
      LN4resNI <- quantile(parsLN4N$BMD*data$maxD, pvec)
      LN4resNI <- c(LN4resNI,apply(parsLN4N[,c(paste0("p",1:4), "is2t","rho")], 2, median))
      names(LN4resNI) <- c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t","rho")

      LN4outNI <- outLPC(parsLN4N, pvec, data$maxD)


      # Covariance between b-d and between BMD-d
      LN4covNI = c(cov(parsLN4N[,c("b","d")], use="complete.obs")["b","d"],
                   cov(parsLN4N[,c("BMD","d")], use="complete.obs")["BMD","d"])

      LN4corrNI = c(cor(parsLN4N[,c("b","d")], use="complete.obs")["b","d"],
                    cor(parsLN4N[,c("BMD","d")], use="complete.obs")["BMD","d"])

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

    # if(data$is_increasing == 1){
    #   data$init_b = qgamma(data$q/
    #                          (optE4_NI$par[9]-1), rate=1.0, shape=optE4_NI$par[11])/optE4_NI$par[2]
    # }else if(data$is_decreasing == 1){
    #   data$init_b = qgamma((-data$q)/
    #                          (optE4_NI$par[9]-1), rate=1.0, shape=optE4_NI$par[11])/optE4_NI$par[2]
    # }

    optG4_NI <- fun_optimC(stanmodels$mG4c, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optG4_NI[[3]]),TRUE,(optG4_NI[[3]]!=0)) | length(optG4_NI)!=9)){
      prior.weights[5] <- 0
      warning('Difficulties fitting the Gamma (Normal) model; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsG4N <- par_extractC(optG4_NI, model_name = "G4_N")
      G4resNI <- quantile(parsG4N$BMD*data$maxD, pvec)
      G4resNI <- c(G4resNI,apply(parsG4N[,c(paste0("p",1:4), "is2t", "rho")], 2, median))
      names(G4resNI) <- c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t","rho")

      G4outNI <- outLPC(parsG4N, pvec, data$maxD)

      # Covariance between b-d and between BMD-d
      G4covNI = c(cov(parsG4N[,c("b","d")], use="complete.obs")["b","d"],
                  cov(parsG4N[,c("BMD","d")], use="complete.obs")["BMD","d"])

      G4corrNI = c(cor(parsG4N[,c("b","d")], use="complete.obs")["b","d"],
                   cor(parsG4N[,c("BMD","d")], use="complete.obs")["BMD","d"])

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
      warning('Difficulties fitting the Quadratic Exponential (Normal) model; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsQE4N <- par_extractC(optQE4_NI, model_name = "QE4_N")
      QE4resNI <- quantile(parsQE4N$BMD*data$maxD, pvec)
      QE4resNI <- c(QE4resNI,apply(parsQE4N[,c(paste0("p",1:4), "is2t", "rho")], 2, median))
      names(QE4resNI) <- c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t","rho")

      QE4outNI <- outLPC(parsQE4N, pvec, data$maxD)

      # Covariance between b-d and between BMD-d
      QE4covNI = c(cov(parsQE4N[,c("b","d")], use="complete.obs")["b","d"],
                   cov(parsQE4N[,c("BMD","d")], use="complete.obs")["BMD","d"])

      QE4corrNI = c(cor(parsQE4N[,c("b","d")], use="complete.obs")["b","d"],
                    cor(parsQE4N[,c("BMD","d")], use="complete.obs")["BMD","d"])

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
      warning('Difficulties fitting the Probit (Normal) model; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsP4N <- par_extractC(optP4_NI, model_name = "P4_N")
      P4resNI <- quantile(parsP4N$BMD*data$maxD, pvec)
      P4resNI <- c(P4resNI,apply(parsP4N[,c(paste0("p",1:4), "is2t","rho")], 2, median))
      names(P4resNI) <- c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t","rho")

      P4outNI <- outLPC(parsP4N, pvec, data$maxD)

      # Covariance between b-d and between BMD-d
      P4covNI = c(cov(parsP4N[,c("b","d")], use="complete.obs")["b","d"],
                  cov(parsP4N[,c("BMD","d")], use="complete.obs")["BMD","d"])

      P4corrNI = c(cor(parsP4N[,c("b","d")], use="complete.obs")["b","d"],
                   cor(parsP4N[,c("BMD","d")], use="complete.obs")["BMD","d"])

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
      warning('Difficulties fitting the Logit (Normal) model; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsL4N <- par_extractC(optL4_NI, model_name = "L4_N")
      L4resNI <- quantile(parsL4N$BMD*data$maxD, pvec)
      L4resNI <- c(L4resNI,apply(parsL4N[,c(paste0("p",1:4), "is2t","rho")], 2, median))
      names(L4resNI) <- c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t","rho")

      L4outNI <- outLPC(parsL4N, pvec, data$maxD)


      # Covariance between b-d and between BMD-d
      L4covNI = c(cov(parsL4N[,c("b","d")], use="complete.obs")["b","d"],
                  cov(parsL4N[,c("BMD","d")], use="complete.obs")["BMD","d"])

      L4corrNI = c(cor(parsL4N[,c("b","d")], use="complete.obs")["b","d"],
                   cor(parsL4N[,c("BMD","d")], use="complete.obs")["BMD","d"])

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
      warning('Difficulties fitting the Exponential (Lognormal) model; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsE4LN <- par_extractC(optE4_LNI, model_name = "E4_LN")
      E4resLNI <- quantile(parsE4LN$BMD*data$maxD, pvec)
      E4resLNI <- c(E4resLNI, apply(parsE4LN[,c(paste0("p",1:4), "is2t","rho")], 2, median))
      names(E4resLNI) <- c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t","rho")

      E4outLNI <- outLPC(parsE4LN, pvec, data$maxD)


      # Covariance between b-d and between BMD-d
      E4covLNI = c(cov(parsE4LN[,c("b","d")], use="complete.obs")["b","d"],
                   cov(parsE4LN[,c("BMD","d")], use="complete.obs")["BMD","d"])

      E4corrLNI = c(cor(parsE4LN[,c("b","d")], use="complete.obs")["b","d"],
                    cor(parsE4LN[,c("BMD","d")], use="complete.obs")["BMD","d"])

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
      warning('Difficulties fitting the Inverse Exponential (Lognormal) model; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsIE4LN <- par_extractC(optIE4_LNI, model_name = "IE4_LN")
      IE4resLNI <- quantile(parsIE4LN$BMD*data$maxD, pvec)
      IE4resLNI <- c(IE4resLNI, apply(parsIE4LN[,c(paste0("p",1:4), "is2t","rho")], 2, median))
      names(IE4resLNI) <- c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t","rho")

      IE4outLNI <- outLPC(parsIE4LN, pvec, data$maxD)

      # Covariance between b-d and between BMD-d
      IE4covLNI = c(cov(parsIE4LN[,c("b","d")], use="complete.obs")["b","d"],
                    cov(parsIE4LN[,c("BMD","d")], use="complete.obs")["BMD","d"])

      IE4corrLNI = c(cor(parsIE4LN[,c("b","d")], use="complete.obs")["b","d"],
                     cor(parsIE4LN[,c("BMD","d")], use="complete.obs")["BMD","d"])

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
      warning('Difficulties fitting the Hill (Lognormal) model; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsH4LN <- par_extractC(optH4_LNI, model_name = "H4_LN")
      H4resLNI <- quantile(parsH4LN$BMD*data$maxD, pvec)
      H4resLNI <- c(H4resLNI, apply(parsH4LN[,c(paste0("p",1:4), "is2t","rho")], 2, median))
      names(H4resLNI) <- c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t","rho")

      H4outLNI <- outLPC(parsH4LN, pvec, data$maxD)

      # Covariance between b-d and between BMD-d
      H4covLNI = c(cov(parsH4LN[,c("b","d")], use="complete.obs")["b","d"],
                   cov(parsH4LN[,c("BMD","d")], use="complete.obs")["BMD","d"])

      H4corrLNI = c(cor(parsH4LN[,c("b","d")], use="complete.obs")["b","d"],
                    cor(parsH4LN[,c("BMD","d")], use="complete.obs")["BMD","d"])

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
      warning('Difficulties fitting the Lognormal (Lognormal) model; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsLN4LN <- par_extractC(optLN4_LNI, model_name = "LN4_LN")
      LN4resLNI <- quantile(parsLN4LN$BMD*data$maxD, pvec)
      LN4resLNI <- c(LN4resLNI, apply(parsLN4LN[,c(paste0("p",1:4), "is2t","rho")], 2, median))
      names(LN4resLNI) <- c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t","rho")

      LN4outLNI <- outLPC(parsLN4LN, pvec, data$maxD)

      # Covariance between b-d and between BMD-d
      LN4covLNI = c(cov(parsLN4LN[,c("b","d")], use="complete.obs")["b","d"],
                    cov(parsLN4LN[,c("BMD","d")], use="complete.obs")["BMD","d"])

      LN4corrLNI = c(cor(parsLN4LN[,c("b","d")], use="complete.obs")["b","d"],
                     cor(parsLN4LN[,c("BMD","d")], use="complete.obs")["BMD","d"])

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

    ## this causes problems (not really helping)
    # if(data$is_increasing == 1){
    #   data$init_b = qgamma(log(1+data$q)/
    #                          (optE4_LNI$par[8]*(optE4_LNI$par[9]-1)), rate=1.0, shape=optE4_LNI$par[11])/optE4_LNI$par[2]
    # }else if(data$is_decreasing == 1){
    #   data$init_b = qgamma(log(1-data$q)/
    #                          (optE4_LNI$par[8]*(optE4_LNI$par[9]-1)), rate=1.0, shape=optE4_LNI$par[11])/optE4_LNI$par[2]
    # }

    optG4_LNI <- fun_optimC(stanmodels$mG4c, data, start, ndraws, 123, pvec)

    if((ifelse(is.na(optG4_LNI[[3]]),TRUE,(optG4_LNI[[3]]!=0)) | length(optG4_LNI)!=9)){
      prior.weights[13] <- 0
      warning('Difficulties fitting the Gamma (Lognormal) model; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsG4LN <- par_extractC(optG4_LNI, model_name = "G4_LN")
      G4resLNI <- quantile(parsG4LN$BMD*data$maxD, pvec)
      G4resLNI <- c(G4resLNI, apply(parsG4LN[,c(paste0("p",1:4), "is2t","rho")], 2, median))
      names(G4resLNI) <- c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t","rho")

      G4outLNI <- outLPC(parsG4LN, pvec, data$maxD)

      # Covariance between b-d and between BMD-d
      G4covLNI = c(cov(parsG4LN[,c("b","d")], use="complete.obs")["b","d"],
                   cov(parsG4LN[,c("BMD","d")], use="complete.obs")["BMD","d"])

      G4corrLNI = c(cor(parsG4LN[,c("b","d")], use="complete.obs")["b","d"],
                    cor(parsG4LN[,c("BMD","d")], use="complete.obs")["BMD","d"])

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
      warning('Difficulties fitting the Quadratic Exponential (Lognormal) model; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsQE4LN <- par_extractC(optQE4_LNI, model_name = "QE4_LN")
      QE4resLNI <- quantile(parsQE4LN$BMD*data$maxD, pvec)
      QE4resLNI <- c(QE4resLNI, apply(parsQE4LN[,c(paste0("p",1:4), "is2t","rho")], 2, median))
      names(QE4resLNI) <- c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t","rho")

      QE4outLNI <- outLPC(parsQE4LN, pvec, data$maxD)

      # Covariance between b-d and between BMD-d
      QE4covLNI = c(cov(parsQE4LN[,c("b","d")], use="complete.obs")["b","d"],
                    cov(parsQE4LN[,c("BMD","d")], use="complete.obs")["BMD","d"])

      QE4corrLNI = c(cor(parsQE4LN[,c("b","d")], use="complete.obs")["b","d"],
                     cor(parsQE4LN[,c("BMD","d")], use="complete.obs")["BMD","d"])

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
      warning('Difficulties fitting the Probit (Lognormal) model; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsP4LN <- par_extractC(optP4_LNI, model_name = "P4_LN")
      P4resLNI <- quantile(parsP4LN$BMD*data$maxD, pvec, data$maxD)
      P4resLNI <- c(P4resLNI, apply(parsP4LN[,c(paste0("p",1:4), "is2t","rho")], 2, median))
      names(P4resLNI) <- c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t","rho")

      P4outLNI <- outLPC(parsP4LN, pvec, data$maxD)


      # Covariance between b-d and between BMD-d
      P4covLNI = c(cov(parsP4LN[,c("b","d")], use="complete.obs")["b","d"],
                   cov(parsP4LN[,c("BMD","d")], use="complete.obs")["BMD","d"])

      P4corrLNI = c(cor(parsP4LN[,c("b","d")], use="complete.obs")["b","d"],
                    cor(parsP4LN[,c("BMD","d")], use="complete.obs")["BMD","d"])

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
      warning('Difficulties fitting the Logit (Lognormal) model; prior weight was set to 0 and the model is not included in model averaging')
    }else{

      parsL4LN <- par_extractC(optL4_LNI, model_name = "L4_LN")
      L4resLNI <- quantile(parsL4LN$BMD*data$maxD, pvec, data$maxD)
      L4resLNI <- c(L4resLNI, apply(parsL4LN[,c(paste0("p",1:4), "is2t","rho")], 2, median))
      names(L4resLNI) <- c("BMDL","BMD","BMDU","min_resp","bmd","fold_change","dt","is2t","rho")

      L4outLNI <- outLPC(parsL4LN, pvec, data$maxD)

      # Covariance between b-d and between BMD-d
      L4covLNI = c(cov(parsL4LN[,c("b","d")], use="complete.obs")["b","d"],
                   cov(parsL4LN[,c("BMD","d")], use="complete.obs")["BMD","d"])

      L4corrLNI = c(cor(parsL4LN[,c("b","d")], use="complete.obs")["b","d"],
                    cor(parsL4LN[,c("BMD","d")], use="complete.obs")["BMD","d"])

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
    # w = c(w,
    #       (2*pi)^(2.5)*sqrt(DIHQE4)*exp(llQE4N-minll)*
    #         # sigma2
    #         dnorm(optQE4_NI$par[5],mean=data$priormu[5],
    #                          sd=data$priorSigma[5,5])*
    #         # d
    #         dexp(optQE4_NI$par[4], 1)*
    #         # a
    #         mc2d::dpert(optQE4_NI$par[1], min = data$priorlb[1], max = data$priorub[1],
    #                     mode = data$priormu[1], shape = data$shape.a)*
    #         # BMD
    #         mc2d::dpert(optQE4_NI$par[2], min = data$priorlb[2], max = data$priorub[2],
    #                     mode = data$priormu[2], shape = data$shape.BMD)*
    #         # c
    #         mc2d::dpert(optQE4_NI$par[10], min = data$priorlb[3], max = data$priorub[3],
    #                     mode = data$priormu[3], shape = data$shape.c)*
    #         # rho
    #         mc2d::dpert(optQE4_NI$par[6], min = data$priorlb[6], max = data$priorub[6],
    #                     mode = data$priormu[6], shape = 0.0001)
    # )
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
    # w = c(w,
    #       (2*pi)^(2.5)*sqrt(DIHQE4)*exp(llQE4LN-minll)*
    #         # sigma2
    #         dnorm(optQE4_LNI$par[5],mean=data$priormu[5],
    #                          sd=data$priorSigma[5,5])*
    #         # d
    #         dexp(optQE4_LNI$par[4], 1)*
    #         # a
    #         mc2d::dpert(optQE4_LNI$par[1], min = data$priorlb[1], max = data$priorub[1],
    #                     mode = data$priormu[1], shape = data$shape.a)*
    #         # BMD
    #         mc2d::dpert(optQE4_LNI$par[2], min = data$priorlb[2], max = data$priorub[2],
    #                     mode = data$priormu[2], shape = data$shape.BMD)*
    #         # c
    #         mc2d::dpert(optQE4_LNI$par[10], min = data$priorlb[3], max = data$priorub[3],
    #                     mode = data$priormu[3], shape = data$shape.c)*
    #         # rho
    #         mc2d::dpert(optQE4_LNI$par[6], min = data$priorlb[6], max = data$priorub[6],
    #                     mode = data$priormu[6], shape = 0.0001)
    # )
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

  prior.weights = prior.weights/sum(prior.weights==1)
  lpw=(prior.weights*w)/sum(prior.weights*w)

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
  maci=quantile(mabmd,pvec)*data$maxD ## original scale
  names(maci)=c("BMDL","BMD","BMDU")

  BMDq = quantile(mabmd, seq(0,1,0.005))*data$maxD ## original scale

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
    print(warning('BMD/BMDL is larger than 20'))
  }
  if(maci[3]/maci[1] > 50){
    print(warning('BMDU/BMDL is larger than 50'))
  }
  if(maci[2] < (data.N$data$x[2]*data.N$data$maxD/10)){
    print(warning('BMD is 10 times lower than the lowest non-zero dose'))
  }

  ### best fitting model vs saturated ANOVA model
  best.fit = modelnames[which(weight[1:16] == max(weight[1:16]))][1]

  bfTest <- modelTestC(best.fit, data.N, data.LN, get(paste0('opt', best.fit, 'I')), type = 'Laplace',
                       seed, ndraws, nrchains, nriterations, warmup, delta, treedepth)
  print(warning(bfTest$warn.bf))


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
    models_included = modelnames[prior.weights > 0],
    bf = bfTest$bayesFactor,
    # means.SM = bfTest$means.SM, parBestFit = bfTest$par.best,
    # BIC.bestfit = bfTest$BIC.bestfit, BIC.SM = bfTest$BIC.SM,
    shift = data.LN$data$shift
  )

  attr(ret_results, "class") <- c("BMADR", "LP")

  return(ret_results)

}

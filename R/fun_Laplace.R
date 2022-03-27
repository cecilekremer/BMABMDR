#' Perform model averaging using Full Laplace method
#'
#' This method assumed data for continuous endpoints.
#'
#' More detailed descriprion
#' @param data.N the input data as returned by function PREP_DATA_N
#' @param data.LN the input data as returned by function PREP_DATA_LN
#' @param prior.weights a vector specifying which of the 16 models should be included (1 = include, 0 = exclude)
#' @param priordist which prior distribution should be used (currently only "PERT" is implemented)
#' @param prior.BMD logical indicating whether an informative prior for the BMD is used (currently not implemented)
#' @param ndraws the number of draws, default 30000
#' @param seed default 123
#' @param pvec vector specifying the three BMD quantiles of interest
#' @param plot logical indicating whether a simple plot of model fits should be shown
#'
#' @return List with the model-specific results, model weights, the model averaged BMD, BMDL and BMDU. In addition for each model a matrix with covariance/correlation between parameter b/BMD and parameter d is given. Some additional output used in other external functions is also given.
#'
#' @export
#'
full.laplace_MA=function(data.N, data.LN,
                         prior.weights,
                         priordist = "PERT",
                         prior.BMD = FALSE,
                         ndraws=30000,seed=123,
                         pvec=c(0.05,0.5,0.95),
                         plot=FALSE){

  # prior.weights = prior.weights/sum(prior.weights==1) #this is now done below when calculating weights

  q = data.N$data.modstan1a$q

  if(priordist == "Normal" & prior.BMD == TRUE){
    stop("Informative prior on BMD can only be used in combination with the PERT prior distribution.")
  }

  if(data.N$increasing == TRUE){

    ## Data to use for Normal distribution
    data.modstan1a = data.N$data.modstan1a
    data.modstan1bG = data.N$data.modstan1bG
    data.modstan1aLN = data.N$data.modstan1aLN
    data.modstanE = data.N$data.modstanE
    data.modstan1bQ = data.N$data.modstan1bQ
    data.modstan2P = data.N$data.modstan2P
    data.modstan2L = data.N$data.modstan2L
    svF1 = data.N$svF1
    svF1LN = data.N$svF1LN
    svF2P = data.N$svF2P
    svF2L = data.N$svF2L

    llN = c() # likelihoods
    BMDL = c()
    BMD = c()
    BMDU = c()

    ### Getting the posterior modes and the hessian to obtain the model specific posterior distributions
    # optimizing function: obtain point estimates by maximizing the joint posterior from the Stan model
    ### Additionally save BMD(L/U) estimates and loglikelihood values; set to NA if model has prior weight 0


    if(prior.weights[1]>0){

      print(1)

      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          optE4_NI <- fun_optim(stanmodels$mE4_NI, data.modstan1a, svF1, ndraws, 123, pvec)
        }else if(prior.BMD == TRUE){
          optE4_NI <- fun_optim(stanmodels$mE4_NI2, data.modstan1a, svF1, ndraws, 123, pvec)
        }
      }else if(priordist == "Normal"){
        optE4_NI <- fun_optim(stanmodels$mE4_NI_n, data.modstan1a, svF1, ndraws, 123, pvec)
      }

      parsE4N <- par_extract(optE4_NI, model_name = "E4_N")
      E4resNI <- quantile(parsE4N$BMD, pvec)  #exp(quantile(optE4_NI$theta_tilde[,"par[2]"],pvec))
      E4resNI <- c(E4resNI,apply(parsE4N[,c(paste0("p",1:4), "is2t")], 2, median))
      names(E4resNI) <- c("BMDL","BMD","BMDU","at","k","ct","dt","is2t")

      E4outNI <- outLP(parsE4N, pvec, data.modstan1a$maxD)

      # Covariance between b-d and between BMD-d
      E4covNI <- c(cov(parsE4N[,c("b","d")], use="complete.obs")["b","d"],
                   cov(parsE4N[,c("BMD","d")], use="complete.obs")["BMD","d"])

      E4corrNI <- c(cor(parsE4N[,c("b","d")], use="complete.obs")["b","d"],
                    cor(parsE4N[,c("BMD","d")], use="complete.obs")["BMD","d"])

      DRM_E4_N <- DRM.E4_NI(E4resNI[4:7], data.modstan1a$x, data.modstan1a$q)

      # obtain loglikelihood
      llE4N=llfE4_NI(optE4_NI$par[1:5],
                     nvec=data.modstan1a$n,
                     dvec=data.modstan1a$x,
                     mvec=data.modstan1a$m,
                     s2vec=data.modstan1a$s2,
                     qval=data.modstan1a$q)
      llN = c(llN, llE4N)
      BMDL = c(BMDL, E4resNI[1])
      BMD = c(BMD, E4resNI[2])
      BMDU = c(BMDU, E4resNI[3])
    }else{E4resNI=NULL; llN = c(llN,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA);
    E4covNI=rep(NA,2); E4corrNI=rep(NA,2); DRM_E4_N=rep(NA, length(data.modstan1a$x))
    parsE4N <- NA
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

    #

    if(prior.weights[2]>0){

      print(2)
      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          optIE4_NI <- fun_optim(stanmodels$mIE4_NI, data.modstan1a, svF1, ndraws, 123, pvec)
        }else if(prior.BMD == TRUE){
          optIE4_NI <- fun_optim(stanmodels$mIE4_NI2, data.modstan1a, svF1, ndraws, 123, pvec)
        }
      }else if(priordist == "Normal"){
        optIE4_NI <- fun_optim(stanmodels$mIE4_NI_n, data.modstan1a, svF1, ndraws, 123, pvec)
      }
      parsIE4N <- par_extract(optIE4_NI, model_name = "IE4_N")
      IE4resNI <- quantile(parsIE4N$BMD, pvec)  #exp(quantile(optE4_NI$theta_tilde[,"par[2]"],pvec))
      IE4resNI <- c(IE4resNI,apply(parsIE4N[,c(paste0("p",1:4), "is2t")], 2, median))
      names(IE4resNI) <- c("BMDL","BMD","BMDU","at","k","ct","dt","is2t")

      IE4outNI <- outLP(parsIE4N, pvec, data.modstan1a$maxD)

      # Covariance between b-d and between BMD-d
      IE4covNI <- c(cov(parsIE4N[,c("b","d")], use="complete.obs")["b","d"],
                    cov(parsIE4N[,c("BMD","d")], use="complete.obs")["BMD","d"])

      IE4corrNI <- c(cor(parsIE4N[,c("b","d")], use="complete.obs")["b","d"],
                     cor(parsIE4N[,c("BMD","d")], use="complete.obs")["BMD","d"])

      DRM_IE4_N = DRM.IE4_NI(IE4resNI[4:7], data.modstan1a$x, data.modstan1a$q)

      llIE4N=llfIE4_NI(optIE4_NI$par[1:5],nvec=data.modstan1a$n,
                       dvec=data.modstan1a$x,mvec=data.modstan1a$m,
                       s2vec=data.modstan1a$s2,qval=data.modstan1a$q)

      llN = c(llN, llIE4N)
      BMDL = c(BMDL, IE4resNI[1])
      BMD = c(BMD, IE4resNI[2])
      BMDU = c(BMDU, IE4resNI[3])
    }else{
      IE4resNI=NULL; llN = c(llN,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA);
      IE4covNI=rep(NA,2); IE4corrNI=rep(NA,2); DRM_IE4_N=rep(NA,length(data.modstan1a$x))
      parsIE4N <- NA
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
    #

    if(prior.weights[3]>0){

      print(3)

      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          optH4_NI <- fun_optim(stanmodels$mH4_NI, data.modstan1a, svF1, ndraws, 123, pvec)
        }else if(prior.BMD == TRUE){
          optH4_NI <- fun_optim(stanmodels$mH4_NI2, data.modstan1a, svF1, ndraws, 123, pvec)
        }
      }else if(priordist == "Normal"){
        optH4_NI <- fun_optim(stanmodels$mH4_NI_n, data.modstan1a, svF1, ndraws, 123, pvec)
      }

      parsH4N <- par_extract(optH4_NI, model_name = "H4_N")
      H4resNI <- quantile(parsH4N$BMD, pvec)  #exp(quantile(optE4_NI$theta_tilde[,"par[2]"],pvec))
      H4resNI <- c(H4resNI,apply(parsH4N[,c(paste0("p",1:4), "is2t")], 2, median))
      names(H4resNI) <- c("BMDL","BMD","BMDU","at","k","ct","dt","is2t")

      H4outNI <- outLP(parsH4N, pvec, data.modstan1a$maxD)

      # Covariance between b-d and between BMD-d
      H4covNI = c(cov(parsH4N[,c("b","d")], use="complete.obs")["b","d"],
                  cov(parsH4N[,c("BMD","d")], use="complete.obs")["BMD","d"])

      H4corrNI = c(cor(parsH4N[,c("b","d")], use="complete.obs")["b","d"],
                   cor(parsH4N[,c("BMD","d")], use="complete.obs")["BMD","d"])

      DRM_H4_N = DRM.H4_NI(H4resNI[4:7], data.modstan1a$x, data.modstan1a$q)

      llH4N=llfH4_NI(optH4_NI$par[1:5],nvec=data.modstan1a$n,
                     dvec=data.modstan1a$x,mvec=data.modstan1a$m,
                     s2vec=data.modstan1a$s2,qval=data.modstan1a$q)

      llN = c(llN, llH4N)
      BMDL = c(BMDL, H4resNI[1])
      BMD = c(BMD, H4resNI[2])
      BMDU = c(BMDU, H4resNI[3])
    }else{H4resNI=NULL; llN = c(llN,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA);
    BMDU=c(BMDU,NA); H4covNI=rep(NA,2); H4corrNI=rep(NA,2); DRM_H4_N=rep(NA,length(data.modstan1a$x))
    parsH4N <- NA
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

    #

    if(prior.weights[4]>0){

      print(4)
      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          optLN4_NI <- fun_optim(stanmodels$mLN4_NI, data.modstan1aLN, svF1LN, ndraws, 123, pvec)
        }else if(prior.BMD == TRUE){
          optLN4_NI <- fun_optim(stanmodels$mLN4_NI2, data.modstan1aLN, svF1LN, ndraws, 123, pvec)
        }
      }else if(priordist == "Normal"){
        optLN4_NI <- fun_optim(stanmodels$mLN4_NI_n, data.modstan1aLN, svF1LN, ndraws, 123, pvec)
      }

      parsLN4N <- par_extract(optLN4_NI, model_name = "LN4_N")
      LN4resNI <- quantile(parsLN4N$BMD, pvec)
      LN4resNI <- c(LN4resNI,apply(parsLN4N[,c(paste0("p",1:4), "is2t")], 2, median))
      names(LN4resNI) <- c("BMDL","BMD","BMDU","at","k","ct","dt","is2t")

      LN4outNI <- outLP(parsLN4N, pvec, data.modstan1a$maxD)


      # Covariance between b-d and between BMD-d
      LN4covNI = c(cov(parsLN4N[,c("b","d")], use="complete.obs")["b","d"],
                   cov(parsLN4N[,c("BMD","d")], use="complete.obs")["BMD","d"])

      LN4corrNI = c(cor(parsLN4N[,c("b","d")], use="complete.obs")["b","d"],
                    cor(parsLN4N[,c("BMD","d")], use="complete.obs")["BMD","d"])

      DRM_LN4_N = DRM.LN4_NI(LN4resNI[4:7], data.modstan1a$x, data.modstan1a$q)


      llLN4N=llfLN4_NI(optLN4_NI$par[1:5],nvec=data.modstan1aLN$n,
                       dvec=data.modstan1aLN$x,mvec=data.modstan1aLN$m,
                       s2vec=data.modstan1aLN$s2,qval=data.modstan1aLN$q)

      llN = c(llN, llLN4N)
      BMDL = c(BMDL, LN4resNI[1])
      BMD = c(BMD, LN4resNI[2])
      BMDU = c(BMDU, LN4resNI[3])
    }else{LN4resNI=NULL; llN = c(llN,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA);
    LN4covNI=rep(NA,2); LN4corrNI=rep(NA,2); DRM_LN4_N=rep(NA,length(data.modstan1aLN$x))
    parsLN4N <- NA
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
    #

    if(prior.weights[5]>0){
      print(5)

      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          optG4_NI <- fun_optim(stanmodels$mG4_NI, data.modstan1bG, svF1, ndraws, 123, pvec)
        }else if(prior.BMD == TRUE){
          optG4_NI <- fun_optim(stanmodels$mG4_NI2, data.modstan1bG, svF1, ndraws, 123, pvec)
        }
      }else if(priordist == "Normal"){
        optG4_NI <- fun_optim(stanmodels$mG4_NI_n, data.modstan1bG, svF1, ndraws, 123, pvec)
      }

      parsG4N <- par_extract(optG4_NI, model_name = "G4_N")
      G4resNI <- quantile(parsG4N$BMD, pvec)
      G4resNI <- c(G4resNI,apply(parsG4N[,c(paste0("p",1:4), "is2t")], 2, median))
      names(G4resNI) <- c("BMDL","BMD","BMDU","at","k","ct","dt","is2t")

      G4outNI <- outLP(parsG4N, pvec, data.modstan1a$maxD)

      # Covariance between b-d and between BMD-d
      G4covNI = c(cov(parsG4N[,c("b","d")], use="complete.obs")["b","d"],
                  cov(parsG4N[,c("BMD","d")], use="complete.obs")["BMD","d"])

      G4corrNI = c(cor(parsG4N[,c("b","d")], use="complete.obs")["b","d"],
                   cor(parsG4N[,c("BMD","d")], use="complete.obs")["BMD","d"])

      DRM_G4_N = DRM.G4_NI(G4resNI[4:7], data.modstan1a$x, data.modstan1a$q)


      llG4N=llfG4_NI(optG4_NI$par[1:5],nvec=data.modstan1bG$n,
                     dvec=data.modstan1bG$x,mvec=data.modstan1bG$m,
                     s2vec=data.modstan1bG$s2,qval=data.modstan1bG$q)

      llN = c(llN, llG4N)
      BMDL = c(BMDL, G4resNI[1])
      BMD = c(BMD, G4resNI[2])
      BMDU = c(BMDU, G4resNI[3])
    }else{G4resNI=NULL; llN = c(llN,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA);
    BMDU=c(BMDU,NA); G4covNI=rep(NA,2); G4corrNI=rep(NA,2); DRM_G4_N=rep(NA, length(data.modstan1bG$x))
    parsG4N <-NA
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
    #

    if(prior.weights[6]>0){
      print(6)

      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          optQE4_NI <- fun_optim(stanmodels$mQE4_NI, data.modstan1bQ, svF1, ndraws, 123, pvec)
        }else if(prior.BMD == TRUE){
          optQE4_NI <- fun_optim(stanmodels$mQE4_NI2, data.modstan1bQ, svF1, ndraws, 123, pvec)
        }
      }else if(priordist == "Normal"){
        optQE4_NI <- fun_optim(stanmodels$mQE4_NI_n, data.modstan1bQ, svF1, ndraws, 123, pvec)
      }

      parsQE4N <- par_extract(optQE4_NI, model_name = "QE4_N")
      QE4resNI <- quantile(parsQE4N$BMD, pvec)
      QE4resNI <- c(QE4resNI,apply(parsQE4N[,c(paste0("p",1:4), "is2t")], 2, median))
      names(QE4resNI) <- c("BMDL","BMD","BMDU","at","k","ct","dt","is2t")

      QE4outNI <- outLP(parsQE4N, pvec, data.modstan1a$maxD)

      # Covariance between b-d and between BMD-d
      QE4covNI = c(cov(parsQE4N[,c("b","d")], use="complete.obs")["b","d"],
                   cov(parsQE4N[,c("BMD","d")], use="complete.obs")["BMD","d"])

      QE4corrNI = c(cor(parsQE4N[,c("b","d")], use="complete.obs")["b","d"],
                    cor(parsQE4N[,c("BMD","d")], use="complete.obs")["BMD","d"])

      DRM_QE4_N = DRM.QE4_NI(QE4resNI[4:7], data.modstan1a$x, data.modstan1a$q)


      llQE4N=llfQE4_NI(optQE4_NI$par[1:5],nvec=data.modstan1bQ$n,
                       dvec=data.modstan1bQ$x,mvec=data.modstan1bQ$m,
                       s2vec=data.modstan1bQ$s2,qval=data.modstan1bQ$q)

      llN = c(llN, llQE4N)
      BMDL = c(BMDL, QE4resNI[1])
      BMD = c(BMD, QE4resNI[2])
      BMDU = c(BMDU, QE4resNI[3])
    }else{QE4resNI=NULL; llN = c(llN,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA);
    BMDU=c(BMDU,NA); QE4covNI=rep(NA,2); QE4corrNI=rep(NA,2); DRM_QE4_N=rep(NA,length(data.modstan1bQ$x))
    parsQE4N <- NA
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
    #

    if(prior.weights[7]>0){
      print(7)

      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          optP4_NI <- fun_optim(stanmodels$mP4_NI, data.modstan2P, svF2P, ndraws, 123, pvec)
        }else if(prior.BMD == TRUE){
          optP4_NI <- fun_optim(stanmodels$mP4_NI2, data.modstan2P, svF2P, ndraws, 123, pvec)
        }
      }else if(priordist == "Normal"){
        optP4_NI <- fun_optim(stanmodels$mP4_NI_n, data.modstan2P, svF2P, ndraws, 123, pvec)
      }

      parsP4N <- par_extract(optP4_NI, model_name = "P4_N")
      P4resNI <- quantile(parsP4N$BMD, pvec)
      P4resNI <- c(P4resNI,apply(parsP4N[,c(paste0("p",1:4), "is2t")], 2, median))
      names(P4resNI) <- c("BMDL","BMD","BMDU","at","k","ct","dt","is2t")

      P4outNI <- outLP(parsP4N, pvec, data.modstan1a$maxD)

      # Covariance between b-d and between BMD-d
      P4covNI = c(cov(parsP4N[,c("b","d")], use="complete.obs")["b","d"],
                  cov(parsP4N[,c("BMD","d")], use="complete.obs")["BMD","d"])

      P4corrNI = c(cor(parsP4N[,c("b","d")], use="complete.obs")["b","d"],
                   cor(parsP4N[,c("BMD","d")], use="complete.obs")["BMD","d"])

      DRM_P4_N = DRM.P4_NI(P4resNI[4:7], data.modstan1a$x, data.modstan1a$q)

      llP4N=llfP4_NI(optP4_NI$par[1:5],nvec=data.modstan2P$n,
                     dvec=data.modstan2P$x,mvec=data.modstan2P$m,
                     s2vec=data.modstan2P$s2,qval=data.modstan2P$q)


      llN = c(llN, llP4N)
      BMDL = c(BMDL, P4resNI[1])
      BMD = c(BMD, P4resNI[2])
      BMDU = c(BMDU, P4resNI[3])
    }else{P4resNI=NULL; llN = c(llN,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA);
    P4covNI=rep(NA,2); P4corrNI=rep(NA,2); DRM_P4_N=rep(NA,length(data.modstan2P$x)); parsP4N <- NA
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
    #

    if(prior.weights[8]>0){
      print(8)

      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          optL4_NI <- fun_optim(stanmodels$mL4_NI, data.modstan2L, svF2L, ndraws, 123, pvec)
        }else if(prior.BMD == TRUE){
          optL4_NI <- fun_optim(stanmodels$mL4_NI2, data.modstan2L, svF2L, ndraws, 123, pvec)
        }
      }else if(priordist == "Normal"){
        optL4_NI <- fun_optim(stanmodels$mL4_NI_n, data.modstan2L, svF2L, ndraws, 123, pvec)
      }

      parsL4N <- par_extract(optL4_NI, model_name = "L4_N")
      L4resNI <- quantile(parsL4N$BMD, pvec)
      L4resNI <- c(L4resNI,apply(parsL4N[,c(paste0("p",1:4), "is2t")], 2, median))
      names(L4resNI) <- c("BMDL","BMD","BMDU","at","k","ct","dt","is2t")

      L4outNI <- outLP(parsL4N, pvec, data.modstan1a$maxD)


      # Covariance between b-d and between BMD-d
      L4covNI = c(cov(parsL4N[,c("b","d")], use="complete.obs")["b","d"],
                  cov(parsL4N[,c("BMD","d")], use="complete.obs")["BMD","d"])

      L4corrNI = c(cor(parsL4N[,c("b","d")], use="complete.obs")["b","d"],
                   cor(parsL4N[,c("BMD","d")], use="complete.obs")["BMD","d"])

      DRM_L4_N = DRM.L4_NI(L4resNI[4:7], data.modstan1a$x, data.modstan1a$q)


      llL4N=llfL4_NI(optL4_NI$par[1:5],nvec=data.modstan2L$n,
                     dvec=data.modstan2L$x,mvec=data.modstan2L$m,
                     s2vec=data.modstan2L$s2,qval=data.modstan2L$q)


      llN = c(llN, llL4N)
      BMDL = c(BMDL, L4resNI[1])
      BMD = c(BMD, L4resNI[2])
      BMDU = c(BMDU, L4resNI[3])
    }else{L4resNI=NULL; llN = c(llN,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA);
    L4covNI=rep(NA,2); L4corrNI=rep(NA,2); DRM_L4_N=rep(NA,length(data.modstan2L$x)); parsL4N <- NA
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


    ## Data to use for Lognormal distribution

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
      print(9)

      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          optE4_LNI <- fun_optim(stanmodels$mE4_LNI, data.modstan1a, svF1, ndraws, 123, pvec)
        }else if(prior.BMD == TRUE){
          optE4_LNI <- fun_optim(stanmodels$mE4_LNI2, data.modstan1a, svF1, ndraws, 123, pvec)
        }
      }else if(priordist == "Normal"){
        optE4_LNI <- fun_optim(stanmodels$mE4_LNI_n, data.modstan1a, svF1, ndraws, 123, pvec)
      }

      parsE4LN <- par_extract(optE4_LNI, model_name = "E4_LN")
      E4resLNI <- quantile(parsE4LN$BMD, pvec)
      E4resLNI <- c(E4resLNI, apply(parsE4LN[,c(paste0("p",1:4), "is2t")], 2, median))
      names(E4resLNI) <- c("BMDL","BMD","BMDU","at","k","ct","dt","is2t")

      E4outLNI <- outLP(parsE4LN, pvec, data.modstan1a$maxD)


      # Covariance between b-d and between BMD-d
      E4covLNI = c(cov(parsE4LN[,c("b","d")], use="complete.obs")["b","d"],
                   cov(parsE4LN[,c("BMD","d")], use="complete.obs")["BMD","d"])

      E4corrLNI = c(cor(parsE4LN[,c("b","d")], use="complete.obs")["b","d"],
                    cor(parsE4LN[,c("BMD","d")], use="complete.obs")["BMD","d"])

      DRM_E4_LN = exp(DRM.E4_LNI(E4resLNI[4:7], data.modstan1a$x, data.modstan1a$q))

      llE4LN=llfE4_LNI(optE4_LNI$par[1:5],nvec=data.modstan1a$n,
                       dvec=data.modstan1a$x,mvec=data.modstan1a$m,
                       s2vec=data.modstan1a$s2,qval=data.modstan1a$q)



      llLN = c(llLN, llE4LN)
      BMDL = c(BMDL, E4resLNI[1])
      BMD = c(BMD, E4resLNI[2])
      BMDU = c(BMDU, E4resLNI[3])
    }else{E4resLNI=NULL; llLN = c(llLN,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA);
    E4covLNI=rep(NA,2); E4corrLNI=rep(NA,2); DRM_E4_LN=rep(NA,length(data.modstan1a$x))
    parsE4LN <- NA
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
    #

    if(prior.weights[10]>0){
      print(10)

      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          optIE4_LNI <- fun_optim(stanmodels$mIE4_LNI, data.modstan1a, svF1, ndraws, 123, pvec)
        }else if(prior.BMD == TRUE){
          optIE4_LNI <- fun_optim(stanmodels$mIE4_LNI2, data.modstan1a, svF1, ndraws, 123, pvec)
        }
      }else if(priordist == "Normal"){
        optIE4_LNI <- fun_optim(stanmodels$mIE4_LNI_n, data.modstan1a, svF1, ndraws, 123, pvec)
      }

      parsIE4LN <- par_extract(optIE4_LNI, model_name = "IE4_LN")
      IE4resLNI <- quantile(parsIE4LN$BMD, pvec)
      IE4resLNI <- c(IE4resLNI, apply(parsIE4LN[,c(paste0("p",1:4), "is2t")], 2, median))
      names(IE4resLNI) <- c("BMDL","BMD","BMDU","at","k","ct","dt","is2t")

      IE4outLNI <- outLP(parsIE4LN, pvec, data.modstan1a$maxD)

      # Covariance between b-d and between BMD-d
      IE4covLNI = c(cov(parsIE4LN[,c("b","d")], use="complete.obs")["b","d"],
                    cov(parsIE4LN[,c("BMD","d")], use="complete.obs")["BMD","d"])

      IE4corrLNI = c(cor(parsIE4LN[,c("b","d")], use="complete.obs")["b","d"],
                     cor(parsIE4LN[,c("BMD","d")], use="complete.obs")["BMD","d"])

      DRM_IE4_LN = exp(DRM.IE4_LNI(IE4resLNI[4:7], data.modstan1a$x, data.modstan1a$q))


      llIE4LN=llfIE4_LNI(optIE4_LNI$par[1:5],nvec=data.modstan1a$n,
                         dvec=data.modstan1a$x,mvec=data.modstan1a$m,
                         s2vec=data.modstan1a$s2,qval=data.modstan1a$q)

      llLN = c(llLN, llIE4LN)
      BMDL = c(BMDL, IE4resLNI[1])
      BMD = c(BMD, IE4resLNI[2])
      BMDU = c(BMDU, IE4resLNI[3])
    }else{IE4resLNI=NULL; llLN = c(llLN,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA);
    IE4covLNI=rep(NA,2); IE4corrLNI=rep(NA,2); DRM_IE4_LN=rep(NA,length(data.modstan1a$x))
    parsIE4LN <- NA
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
    #

    if(prior.weights[11]>0){
      print(11)

      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          optH4_LNI <- fun_optim(stanmodels$mH4_LNI, data.modstan1a, svF1, ndraws, 123, pvec)
        }else if(prior.BMD == TRUE){
          optH4_LNI <- fun_optim(stanmodels$mH4_LNI2, data.modstan1a, svF1, ndraws, 123, pvec)
        }
      }else if(priordist == "Normal"){
        optH4_LNI <- fun_optim(stanmodels$mH4_LNI_n, data.modstan1a, svF1, ndraws, 123, pvec)
      }

      parsH4LN <- par_extract(optH4_LNI, model_name = "H4_LN")
      H4resLNI <- quantile(parsH4LN$BMD, pvec)
      H4resLNI <- c(H4resLNI, apply(parsH4LN[,c(paste0("p",1:4), "is2t")], 2, median))
      names(H4resLNI) <- c("BMDL","BMD","BMDU","at","k","ct","dt","is2t")

      H4outLNI <- outLP(parsH4LN, pvec, data.modstan1a$maxD)

      # Covariance between b-d and between BMD-d
      H4covLNI = c(cov(parsH4LN[,c("b","d")], use="complete.obs")["b","d"],
                   cov(parsH4LN[,c("BMD","d")], use="complete.obs")["BMD","d"])

      H4corrLNI = c(cor(parsH4LN[,c("b","d")], use="complete.obs")["b","d"],
                    cor(parsH4LN[,c("BMD","d")], use="complete.obs")["BMD","d"])

      DRM_H4_LN = exp(DRM.H4_LNI(H4resLNI[4:7], data.modstan1a$x, data.modstan1a$q))


      llH4LN=llfH4_LNI(optH4_LNI$par[1:5],nvec=data.modstan1a$n,
                       dvec=data.modstan1a$x,mvec=data.modstan1a$m,
                       s2vec=data.modstan1a$s2,qval=data.modstan1a$q)

      llLN = c(llLN, llH4LN)
      BMDL = c(BMDL, H4resLNI[1])
      BMD = c(BMD, H4resLNI[2])
      BMDU = c(BMDU, H4resLNI[3])
    }else{H4resLNI=NULL; llLN = c(llLN,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA);
    H4covLNI=rep(NA,2); H4corrLNI=rep(NA,2); DRM_H4_LN=rep(NA,length(data.modstan1a$x))
    parsH4LN <- NA
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
    #

    if(prior.weights[12]>0){
      print(12)

      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          optLN4_LNI <- fun_optim(stanmodels$mLN4_LNI, data.modstan1aLN, svF1LN, ndraws, 123, pvec)
        }else if(prior.BMD == TRUE){
          optLN4_LNI <- fun_optim(stanmodels$mLN4_LNI2, data.modstan1aLN, svF1LN, ndraws, 123, pvec)
        }
      }else if(priordist == "Normal"){
        optLN4_LNI <- fun_optim(stanmodels$mLN4_LNI_n, data.modstan1aLN, svF1LN, ndraws, 123, pvec)
      }

      parsLN4LN <- par_extract(optLN4_LNI, model_name = "LN4_LN")
      LN4resLNI <- quantile(parsLN4LN$BMD, pvec)
      LN4resLNI <- c(LN4resLNI, apply(parsLN4LN[,c(paste0("p",1:4), "is2t")], 2, median))
      names(LN4resLNI) <- c("BMDL","BMD","BMDU","at","k","ct","dt","is2t")

      LN4outLNI <- outLP(parsLN4LN, pvec, data.modstan1a$maxD)

      # Covariance between b-d and between BMD-d
      LN4covLNI = c(cov(parsLN4LN[,c("b","d")], use="complete.obs")["b","d"],
                    cov(parsLN4LN[,c("BMD","d")], use="complete.obs")["BMD","d"])

      LN4corrLNI = c(cor(parsLN4LN[,c("b","d")], use="complete.obs")["b","d"],
                     cor(parsLN4LN[,c("BMD","d")], use="complete.obs")["BMD","d"])

      DRM_LN4_LN = exp(DRM.LN4_LNI(LN4resLNI[4:7], data.modstan1a$x, data.modstan1a$q))


      llLN4LN=llfLN4_LNI(optLN4_LNI$par[1:5],nvec=data.modstan1aLN$n,
                         dvec=data.modstan1aLN$x,mvec=data.modstan1aLN$m,
                         s2vec=data.modstan1aLN$s2,qval=data.modstan1aLN$q)


      llLN = c(llLN, llLN4LN)
      BMDL = c(BMDL, LN4resLNI[1])
      BMD = c(BMD, LN4resLNI[2])
      BMDU = c(BMDU, LN4resLNI[3])
    }else{LN4resLNI=NULL; llLN = c(llLN,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA);
    LN4covLNI=rep(NA,2); LN4corrLNI=rep(NA,2); DRM_LN4_LN=rep(NA,length(data.modstan1aLN$x))
    parsLN4LN <- NA
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
    #

    if(prior.weights[13]>0){
      print(13)

      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          optG4_LNI <- fun_optim(stanmodels$mG4_LNI, data.modstan1bG, svF1, ndraws, 123, pvec)
        }else if(prior.BMD == TRUE){
          optG4_LNI <- fun_optim(stanmodels$mG4_LNI2, data.modstan1bG, svF1, ndraws, 123, pvec)
        }
      }else if(priordist == "Normal"){
        optG4_LNI <- fun_optim(stanmodels$mG4_LNI_n, data.modstan1bG, svF1, ndraws, 123, pvec)
      }

      parsG4LN <- par_extract(optG4_LNI, model_name = "G4_LN")
      G4resLNI <- quantile(parsG4LN$BMD, pvec)
      G4resLNI <- c(G4resLNI, apply(parsG4LN[,c(paste0("p",1:4), "is2t")], 2, median))
      names(G4resLNI) <- c("BMDL","BMD","BMDU","at","k","ct","dt","is2t")

      G4outLNI <- outLP(parsG4LN, pvec, data.modstan1a$maxD)

      # Covariance between b-d and between BMD-d
      G4covLNI = c(cov(parsG4LN[,c("b","d")], use="complete.obs")["b","d"],
                   cov(parsG4LN[,c("BMD","d")], use="complete.obs")["BMD","d"])

      G4corrLNI = c(cor(parsG4LN[,c("b","d")], use="complete.obs")["b","d"],
                    cor(parsG4LN[,c("BMD","d")], use="complete.obs")["BMD","d"])

      DRM_G4_LN = exp(DRM.G4_LNI(G4resLNI[4:7], data.modstan1a$x, data.modstan1a$q))


      llG4LN=llfG4_LNI(optG4_LNI$par[1:5],nvec=data.modstan1bG$n,
                       dvec=data.modstan1bG$x,mvec=data.modstan1bG$m,
                       s2vec=data.modstan1bG$s2,qval=data.modstan1bG$q)

      llLN = c(llLN, llG4LN)
      BMDL = c(BMDL, G4resLNI[1])
      BMD = c(BMD, G4resLNI[2])
      BMDU = c(BMDU, G4resLNI[3])
    }else{G4resLNI=NULL; llLN = c(llLN,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA);
    BMDU=c(BMDU,NA); G4covLNI=rep(NA,2); G4corrLNI=rep(NA,2); DRM_G4_LN=rep(NA,length(data.modstan1bG$x))
    parsG4LN <- NA
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
    #

    if(prior.weights[14]>0){
      print(14)

      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          optQE4_LNI <- fun_optim(stanmodels$mQE4_LNI, data.modstan1bQ, svF1, ndraws, 123, pvec)
        }else if(prior.BMD == TRUE){
          optQE4_LNI <- fun_optim(stanmodels$mQE4_LNI2, data.modstan1bQ, svF1, ndraws, 123, pvec)
        }
      }else if(priordist == "Normal"){
        optQE4_LNI <- fun_optim(stanmodels$mQE4_LNI_n, data.modstan1bQ, svF1, ndraws, 123, pvec)
      }

      parsQE4LN <- par_extract(optQE4_LNI, model_name = "QE4_LN")
      QE4resLNI <- quantile(parsQE4LN$BMD, pvec)
      QE4resLNI <- c(QE4resLNI, apply(parsQE4LN[,c(paste0("p",1:4), "is2t")], 2, median))
      names(QE4resLNI) <- c("BMDL","BMD","BMDU","at","k","ct","dt","is2t")

      QE4outLNI <- outLP(parsQE4LN, pvec, data.modstan1a$maxD)

      # Covariance between b-d and between BMD-d
      QE4covLNI = c(cov(parsQE4LN[,c("b","d")], use="complete.obs")["b","d"],
                    cov(parsQE4LN[,c("BMD","d")], use="complete.obs")["BMD","d"])

      QE4corrLNI = c(cor(parsQE4LN[,c("b","d")], use="complete.obs")["b","d"],
                     cor(parsQE4LN[,c("BMD","d")], use="complete.obs")["BMD","d"])

      DRM_QE4_LN = exp(DRM.QE4_LNI(QE4resLNI[4:7], data.modstan1a$x, data.modstan1a$q))


      llQE4LN=llfQE4_LNI(optQE4_LNI$par[1:5],nvec=data.modstan1bQ$n,
                         dvec=data.modstan1bQ$x,mvec=data.modstan1bQ$m,
                         s2vec=data.modstan1bQ$s2,qval=data.modstan1bQ$q)

      llLN = c(llLN, llQE4LN)
      BMDL = c(BMDL, QE4resLNI[1])
      BMD = c(BMD, QE4resLNI[2])
      BMDU = c(BMDU, QE4resLNI[3])
    }else{QE4resLNI=NULL; llLN = c(llLN,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA);
    QE4covLNI=rep(NA,2); QE4corrLNI=rep(NA,2); DRM_QE4_LN=rep(NA,length(data.modstan1bQ$x))
    parsQE4LN <- NA
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
    #

    if(prior.weights[15]>0){
      print(15)

      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          optP4_LNI <- fun_optim(stanmodels$mP4_LNI, data.modstan2P, svF2P, ndraws, 123, pvec)
        }else if(prior.BMD == TRUE){
          optP4_LNI <- fun_optim(stanmodels$mP4_LNI2, data.modstan2P, svF2P, ndraws, 123, pvec)
        }
      }else if(priordist == "Normal"){
        optP4_LNI <- fun_optim(stanmodels$mP4_LNI_n, data.modstan2P, svF2P, ndraws, 123, pvec)
      }

      parsP4LN <- par_extract(optP4_LNI, model_name = "P4_LN")
      P4resLNI <- quantile(parsP4LN$BMD, pvec, data.modstan1a$maxD)
      P4resLNI <- c(P4resLNI, apply(parsP4LN[,c(paste0("p",1:4), "is2t")], 2, median))
      names(P4resLNI) <- c("BMDL","BMD","BMDU","at","k","ct","dt","is2t")

      P4outLNI <- outLP(parsP4LN, pvec, data.modstan1a$maxD)


      # Covariance between b-d and between BMD-d
      P4covLNI = c(cov(parsP4LN[,c("b","d")], use="complete.obs")["b","d"],
                   cov(parsP4LN[,c("BMD","d")], use="complete.obs")["BMD","d"])

      P4corrLNI = c(cor(parsP4LN[,c("b","d")], use="complete.obs")["b","d"],
                    cor(parsP4LN[,c("BMD","d")], use="complete.obs")["BMD","d"])

      DRM_P4_LN = exp(DRM.P4_LNI(P4resLNI[4:7], data.modstan1a$x, data.modstan1a$q))



      llP4LN=llfP4_LNI(optP4_LNI$par[1:5],nvec=data.modstan2P$n,
                       dvec=data.modstan2P$x,mvec=data.modstan2P$m,
                       s2vec=data.modstan2P$s2,qval=data.modstan2P$q)

      llLN = c(llLN, llP4LN)
      BMDL = c(BMDL, P4resLNI[1])
      BMD = c(BMD, P4resLNI[2])
      BMDU = c(BMDU, P4resLNI[3])
    }else{P4resLNI=NULL; llLN = c(llLN,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA);
    P4covLNI=rep(NA,2); P4corrLNI=rep(NA,2); DRM_P4_LN=rep(NA,length(data.modstan2P$x))
    parsP4LN <- NA
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
    #

    if(prior.weights[16]>0){
      print(16)

      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          optL4_LNI <- fun_optim(stanmodels$mL4_LNI, data.modstan2L, svF2L, ndraws, 123, pvec)
        }else if(prior.BMD == TRUE){
          optL4_LNI <- fun_optim(stanmodels$mL4_LNI2, data.modstan2L, svF2L, ndraws, 123, pvec)
        }
      }else if(priordist == "Normal"){
        optL4_LNI <- fun_optim(stanmodels$mL4_LNI_n, data.modstan2L, svF2L, ndraws, 123, pvec)
      }

      parsL4LN <- par_extract(optL4_LNI, model_name = "L4_LN")
      L4resLNI <- quantile(parsL4LN$BMD, pvec, data.modstan1a$maxD)
      L4resLNI <- c(L4resLNI, apply(parsL4LN[,c(paste0("p",1:4), "is2t")], 2, median))
      names(L4resLNI) <- c("BMDL","BMD","BMDU","at","k","ct","dt","is2t")

      L4outLNI <- outLP(parsL4LN, pvec, data.modstan1a$maxD)

      # Covariance between b-d and between BMD-d
      L4covLNI = c(cov(parsL4LN[,c("b","d")], use="complete.obs")["b","d"],
                   cov(parsL4LN[,c("BMD","d")], use="complete.obs")["BMD","d"])

      L4corrLNI = c(cor(parsL4LN[,c("b","d")], use="complete.obs")["b","d"],
                    cor(parsL4LN[,c("BMD","d")], use="complete.obs")["BMD","d"])

      DRM_L4_LN = exp(DRM.L4_LNI(L4resLNI[4:7], data.modstan1a$x, data.modstan1a$q))


      llL4LN=llfL4_LNI(optL4_LNI$par[1:5],nvec=data.modstan2L$n,
                       dvec=data.modstan2L$x,mvec=data.modstan2L$m,
                       s2vec=data.modstan2L$s2,qval=data.modstan2L$q)

      llLN = c(llLN, llL4LN)
      BMDL = c(BMDL, L4resLNI[1])
      BMD = c(BMD, L4resLNI[2])
      BMDU = c(BMDU, L4resLNI[3])
    }else{L4resLNI=NULL; llLN = c(llLN,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA);
    L4covLNI=rep(NA,2); L4corrLNI=rep(NA,2); DRM_L4_LN=rep(NA,length(data.modstan2L$x))
    parsL4LN <- NA
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

    # print(llN); print(llLN); print(w)

    prior.weights = prior.weights/sum(prior.weights==1) # in case this has changed for G4
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
    maci=exp(quantile(mabmd,pvec))*data.modstan1a$maxD ## original scale
    names(maci)=c("BMDL","BMD","BMDU")

    BMDq = exp(quantile(mabmd, seq(0,1,0.005)))*data.modstan1a$maxD ## original scale

    BMDL = c(BMDL, maci[1]/data.modstan1a$maxD); BMD = c(BMD, maci[2]/data.modstan1a$maxD); BMDU = c(BMDU, maci[3]/data.modstan1a$maxD)

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
    for(i in 1:length(data.modstan1a$x)){
      dr.MA[i] = weighted.mean(x = c(DRM_E4_N[i],DRM_IE4_N[i],DRM_H4_N[i],DRM_LN4_N[i],DRM_G4_N[i],DRM_QE4_N[i],
                                     DRM_P4_N[i], DRM_L4_N[i] ,
                                     DRM_E4_LN[i], DRM_IE4_LN[i], DRM_H4_LN[i], DRM_LN4_LN[i], DRM_G4_LN[i],
                                     DRM_QE4_LN[i], DRM_P4_LN[i],DRM_L4_LN[i]),
                               w = lpw,
                               na.rm = T)
    }

    if(plot == TRUE){

      data.modstan1a = data.N$data.modstan1a
      data.modstan1bG = data.N$data.modstan1bG
      data.modstanE = data.N$data.modstanE
      data.modstan1bQ = data.N$data.modstan1bQ
      data.modstan1aLN=data.N$data.modstan1aLN
      data.modstan2P = data.N$data.modstan2P
      data.modstan2L = data.N$data.modstan2L
      svF1 = data.N$svF1
      svF1LN = data.N$svF1LN
      svF2P = data.N$svF2P
      svF2L = data.N$svF2L

      plot2 = function(){
        sd = sqrt(data.modstan1a$s2)
        N = data.modstan1a$N
        gmean.a = log(NtoLN(data.modstan1a$m,sd))[1:N]
        gsd.a = log(NtoLN(data.modstan1a$m,sd))[(N+1):(2*N)]
        # x = c(log10(data.modstan1a$x[2]*data.modstan1a$maxD/4),log10(data.modstan1a$x[2:data.modstan1a$N]*data.modstan1a$maxD))
        x = c(log10(data.modstan1a$x[2]/4),log10(data.modstan1a$x[2:data.modstan1a$N]))
        m = log10(exp(gmean.a))
        s = log10(exp(gsd.a))
        bmdl = log10(BMDL/data.modstan1a$maxD)
        x[1] = ifelse(min(bmdl,na.rm=T)<x[1], log10(min(BMDL/data.modstan1a$maxD/2, na.rm=T)),
                      log10(data.modstan1a$x[2]/4))
        g = choose(length(x),2)
        alph = 0.05/(2*g)
        cp = qnorm(1-alph)
        n = data.modstan1a$n
        ci = cp*s/sqrt(n)
        ci2 = 1.96*s/sqrt(n)

        plot(x, m, xlab="log10-Dose", ylab="log10-Response", ylim=c(min(m-2.2*ci),max(m+2.2*ci)),
             xlim=c(ifelse(min(bmdl,na.rm=T)<min(x),min(bmdl,na.rm=T),min(x)),abs(min(x)))
        )#, ylim=c(min(mean.a)-2*sd.a, max(mean.a)+2*sd.a))
        dgr=seq(min(x), abs(min(x)),by=0.01)
        # dgr = seq(min(x), max(x), by = 0.01)
        # plot model specific fits
        # normal distribution
        if(prior.weights[1]>0){
          parE4 = E4resNI[4:7]
          DRME4 = log10((DRM.E4_NI(par=parE4,x=10^(dgr[dgr>(x[2]-((x[2]-x[1])/2))]),q=data.modstan1a$q)))
          lines(dgr[dgr>(x[2]-((x[2]-x[1])/2))],DRME4,lwd=1, lty=1, col=1)
          segments(x0=dgr[1], y0=log10(exp(parE4[1])), x1=max(dgr[dgr<=(x[2]-((x[2]-x[1])/2))]),
                   y1=DRME4[1],lty=3,col=1)
        }
        if(prior.weights[2]>0){
          parIE4 = IE4resNI[4:7]
          DRMIE4 = log10((DRM.IE4_NI(par=parIE4,x=10^(dgr[dgr>(x[2]-((x[2]-x[1])/2))]),q=data.modstan1a$q)))
          lines(dgr[dgr>(x[2]-((x[2]-x[1])/2))],DRMIE4,lwd=1, lty=1, col=2)
          segments(x0=dgr[1], y0=log10(exp(parIE4[1])), x1=max(dgr[dgr<=(x[2]-((x[2]-x[1])/2))]),
                   y1=DRMIE4[1],lty=3,col=2)
        }
        if(prior.weights[3]>0){
          parH4 = H4resNI[4:7]
          DRMH4 = log10((DRM.H4_NI(par=parH4,x=10^(dgr[dgr>(x[2]-((x[2]-x[1])/2))]),q=data.modstan1a$q)))
          lines(dgr[dgr>(x[2]-((x[2]-x[1])/2))],DRMH4,lwd=1, lty=1, col=3)
          segments(x0=dgr[1], y0=log10(exp(parH4[1])), x1=max(dgr[dgr<=(x[2]-((x[2]-x[1])/2))]),
                   y1=DRMH4[1],lty=3,col=3)
        }
        if(prior.weights[4]>0){
          parLN4 = LN4resNI[4:7]
          DRMLN4 = log10((DRM.LN4_NI(par=parLN4,x=10^(dgr[dgr>(x[2]-((x[2]-x[1])/2))]),q=data.modstan1a$q)))
          lines(dgr[dgr>(x[2]-((x[2]-x[1])/2))],DRMLN4,lwd=1, lty=1, col=4)
          segments(x0=dgr[1], y0=log10(exp(parLN4[1])), x1=max(dgr[dgr<=(x[2]-((x[2]-x[1])/2))]),
                   y1=DRMLN4[1],lty=3,col=4)
        }
        if(prior.weights[5]>0){
          parG4 = G4resNI[4:7]
          DRMG4 = log10((DRM.G4_NI(par=parG4,x=10^(dgr[dgr>(x[2]-((x[2]-x[1])/2))]),q=data.modstan1a$q)))
          lines(dgr[dgr>(x[2]-((x[2]-x[1])/2))],DRMG4,lwd=1, lty=1, col=5)
          segments(x0=dgr[1], y0=log10(exp(parG4[1])), x1=max(dgr[dgr<=(x[2]-((x[2]-x[1])/2))]),
                   y1=DRMG4[1],lty=3,col=5)
        }
        if(prior.weights[6]>0){
          parQE4 = QE4resNI[4:7]
          DRMQE4 = log10((DRM.QE4_NI(par=parQE4,x=10^(dgr[dgr>(x[2]-((x[2]-x[1])/2))]),q=data.modstan1a$q)))
          lines(dgr[dgr>(x[2]-((x[2]-x[1])/2))],DRMQE4,lwd=1, lty=1, col=6)
          segments(x0=dgr[1], y0=log10(exp(parQE4[1])), x1=max(dgr[dgr<=(x[2]-((x[2]-x[1])/2))]),
                   y1=DRMQE4[1],lty=3,col=6)
        }
        if(prior.weights[7]>0){
          parP4 = P4resNI[4:7]
          DRMP4 = log10((DRM.P4_NI(par=parP4,x=10^(dgr[dgr>(x[2]-((x[2]-x[1])/2))]),q=data.modstan1a$q)))
          lines(dgr[dgr>(x[2]-((x[2]-x[1])/2))],DRMP4,lwd=1, lty=1, col=7)
          segments(x0=dgr[1], y0=log10(exp(parP4[1])*pnorm(qnorm(1/(1+q))-exp(parP4[3]))), x1=max(dgr[dgr<=(x[2]-((x[2]-x[1])/2))]), y1=DRMP4[1],lty=3,col=7)
          # segments(x0=dgr[1], y0=log10(exp(parP4[1])*pnorm(parP4[3])), x1=max(dgr[dgr<=(x[2]-((x[2]-x[1])/2))]),
          # y1=DRMP4[1],lty=3,col=7)
        }
        if(prior.weights[8]>0){
          parL4 = L4resNI[4:7]
          DRML4 = log10((DRM.L4_NI(par=parL4,x=10^(dgr[dgr>(x[2]-((x[2]-x[1])/2))]),q=data.modstan1a$q)))
          lines(dgr[dgr>(x[2]-((x[2]-x[1])/2))],DRML4,lwd=1, lty=1, col=8)
          segments(x0=dgr[1], y0=log10(exp(parL4[1])*expit(logit(1/(1+q))-exp(parL4[3]))), x1=max(dgr[dgr<=(x[2]-((x[2]-x[1])/2))]), y1=DRML4[1],lty=3,col=8)
          # segments(x0=dgr[1], y0=log10(exp(parL4[1])*expit(parL4[3])), x1=max(dgr[dgr<=(x[2]-((x[2]-x[1])/2))]),
          # y1=DRML4[1],lty=3,col=8)
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
          if(data.LN$data.modstan1a$shift==T) {DRMIE4 = log10(exp(DRM.IE4_LNI(par=parIE4,
                                                                              x=10^(dgr[dgr>(x[2]-((x[2]-x[1])/2))]),
                                                                              q=data.modstan1a$q))) +
            log10(exp(1.5*min(data.LN$data.modstan1a$m.org)))
          a = exp(parIE4[1]) + 1.5*min(data.LN$data.modstan1a$m.org)}
          lines(dgr[dgr>(x[2]-((x[2]-x[1])/2))],DRMIE4,lwd=1, lty=2, col=2)
          segments(x0=dgr[1], y0=log10(exp(a)), x1=max(dgr[dgr<=(x[2]-((x[2]-x[1])/2))]), y1=min(DRMIE4),
                   lty=3, col=2)
        }
        if(prior.weights[11]>0){
          parH4 = H4resLNI[4:7]
          a=exp(parH4[1])
          DRMH4 = log10(exp(DRM.H4_LNI(par=parH4,x=10^(dgr[dgr>(x[2]-((x[2]-x[1])/2))]),q=data.modstan1a$q)))
          if(data.LN$data.modstan1a$shift==T) {DRMH4 = log10(exp(DRM.H4_LNI(par=parH4,
                                                                            x=10^(dgr[dgr>(x[2]-((x[2]-x[1])/2))]),
                                                                            q=data.modstan1a$q))) +
            log10(exp(1.5*min(data.LN$data.modstan1a$m.org)))
          a = exp(parH4[1]) + 1.5*min(data.LN$data.modstan1a$m.org)}
          lines(dgr[dgr>(x[2]-((x[2]-x[1])/2))],DRMH4,lwd=1, lty=2, col=3)
          segments(x0=dgr[1], y0=log10(exp(a)), x1=max(dgr[dgr<=(x[2]-((x[2]-x[1])/2))]), y1=min(DRMH4),lty=3,col=3)
        }
        if(prior.weights[12]>0){
          parLN4 = LN4resLNI[4:7]
          a=exp(parLN4[1])
          DRMLN4 = log10(exp(DRM.LN4_LNI(par=parLN4,x=10^(dgr[dgr>(x[2]-((x[2]-x[1])/2))]),q=data.modstan1a$q)))
          if(data.LN$data.modstan1aLN$shift==T) {DRMLN4 = log10(exp(DRM.LN4_LNI(par=parLN4,
                                                                                x=10^(dgr[dgr>(x[2]-((x[2]-x[1])/2))]),
                                                                                q=data.modstan1a$q))) +
            log10(exp(1.5*min(data.LN$data.modstan1aLN$m.org)))
          a = exp(parLN4[1]) + 1.5*min(data.LN$data.modstan1aLN$m.org)}
          lines(dgr[dgr>(x[2]-((x[2]-x[1])/2))],DRMLN4,lwd=1, lty=2, col=4)
          segments(x0=dgr[1], y0=log10(exp(a)), x1=max(dgr[dgr<=(x[2]-((x[2]-x[1])/2))]), y1=min(DRMLN4),lty=3,col=4)
        }
        if(prior.weights[13]>0){
          parG4 = G4resLNI[4:7]
          a=exp(parG4[1])
          DRMG4 = log10(exp(DRM.G4_LNI(par=parG4,x=10^(dgr[dgr>(x[2]-((x[2]-x[1])/2))]),q=data.modstan1a$q)))
          if(data.LN$data.modstan1a$shift==T) {DRMG4 = log10(exp(DRM.G4_LNI(par=parG4,
                                                                            x=10^(dgr[dgr>(x[2]-((x[2]-x[1])/2))]),
                                                                            q=data.modstan1a$q))) +
            log10(exp(1.5*min(data.LN$data.modstan1a$m.org)))
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
          a=exp(parP4[1])*pnorm(qnorm(1 - (log(1+q)/exp(parP4[1])))-exp(parP4[3]))
          # a=exp(parP4[1])*pnorm(parP4[3])
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
        leg.all = c(paste0("E4_N (w=",round(weight[1],4),")"),paste0("IE4_N (w=",round(weight[2],4),")"),
                    paste0("H4_N (w=",round(weight[3],4),")"),
                    paste0("LN4_N (w=",round(weight[4],4),")"),paste0("G4_N (w=",round(weight[5],4),")"),
                    paste0("QE4_N (w=",round(weight[6],4),")"),
                    paste0("P4_N (w=",round(weight[7],4),")"),paste0("L4_N (w=",round(weight[8],4),")"),
                    paste0("E4_LN (w=",round(weight[9],4),")"),paste0("IE4_LN (w=",round(weight[10],4),")"),
                    paste0("H4_LN (w=",round(weight[11],4),")"),
                    paste0("LN4_LN (w=",round(weight[12],4),")"),paste0("G4_LN (w=",round(weight[13],4),")"),
                    paste0("QE4_LN (w=",round(weight[14],4),")"),
                    paste0("P4_LN (w=",round(weight[15],4),")"),paste0("L4_LN (w=",round(weight[16],4),")"), "Model averaged BMDL")

        legend(x="bottomright", legend=leg.all[prior.weights>0], col=c(1:8,1:8,1)[prior.weights>0], lwd=c(rep(1,16))[prior.weights>0],
               lty=c(rep(1,8),rep(2,8),NA)[prior.weights>0], pch=c(rep(NA,16),3)[prior.weights>0],
               cex=0.7, title="Model-specific results (fit + BMDL)")
        lines(rbind(x,x,NA),rbind(m-ci,m+ci,NA),lty=3)
        lines(rbind(x,x,NA),rbind(m-ci2,m+ci2,NA))
        points(bmdl[1:16],y=rep(min(m-2*ci),16),col=c(1:8,1:8),pch=c(rep(2,8),rep(5,8)))
        points(bmdl[17],y=min(m-1.5*ci),pch=3)
        mtext("Full Laplace")
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

    modelnames = c("E4_N","IE4_N","H4_N","LN4_N","G4_N","QE4_N","P4_N","L4_N","E4_LN","IE4_LN","H4_LN","LN4_LN","G4_LN","QE4_LN","P4_LN","L4_LN")


    ###
    # sd = sqrt(data.modstan1a$s2)
    # N = data.modstan1a$N
    # gmean.a = log(NtoLN(data.modstan1a$m,sd))[1:N]
    # gsd.a = log(NtoLN(data.modstan1a$m,sd))[(N+1):(2*N)]

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
      BMDMixture = exp(mabmd)*data.modstan1a$maxD,
      data = data.frame(
        dose = c(data.N$data.modstan1a$x),
        sd = sqrt(data.N$data.modstan1a$s2),
        m = data.N$data.modstan1a$m),
      max.dose = data.N$data.modstan1a$maxD,
      q = data.N$data.modstan1a$q,
      increasing = T,
      models_included = modelnames[prior.weights > 0]
    )

    attr(ret_results, "class") <- c("BMADR", "LP")

    return(ret_results)

    ### DECREASING

  }else if(data.N$increasing == FALSE){


    ## Data to use for Normal distribution
    data.modstan1a = data.N$data.modstan1a
    data.modstan1bG = data.N$data.modstan1bG
    data.modstan1aLN = data.N$data.modstan1aLN
    data.modstanE = data.N$data.modstanE
    data.modstan1bQ = data.N$data.modstan1bQ
    data.modstan2P = data.N$data.modstan2P
    data.modstan2L = data.N$data.modstan2L
    svF1 = data.N$svF1
    svF1LN = data.N$svF1LN
    svF2P = data.N$svF2P
    svF2L = data.N$svF2L

    llN = c() # likelihoods
    BMDL = c()
    BMD = c()
    BMDU = c()

    ### Getting the posterior modes and the hessian to obtain the model specific posterior distributions
    # optimizing function: obtain point estimates by maximizing the joint posterior from the Stan model
    ### Additionally save BMD(L/U) estimates and loglikelihood values; set to NA if model has prior weight 0


    if(prior.weights[1]>0){
      print(1)

      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          optE4_ND <- fun_optim(stanmodels$mE4_ND, data.modstan1a, svF1, ndraws, 123, pvec)
        }else if(prior.BMD == TRUE){
          optE4_ND <- fun_optim(stanmodels$mE4_ND2, data.modstan1a, svF1, ndraws, 123, pvec)
        }
      }else if(priordist == "Normal"){
        optE4_ND <- fun_optim(stanmodels$mE4_ND_n, data.modstan1a, svF1, ndraws, 123, pvec)
      }

      parsE4N <- par_extract(optE4_ND, model_name = "E4_N")
      E4resND <- quantile(parsE4N$BMD, pvec)  #exp(quantile(optE4_NI$theta_tilde[,"par[2]"],pvec))
      E4resND <- c(E4resND,apply(parsE4N[,c(paste0("p",1:4), "is2t")], 2, median))
      names(E4resND) <- c("BMDL","BMD","BMDU","at","k","ct","dt","is2t")

      E4outND <- outLP(parsE4N, pvec, data.modstan1a$maxD)

      # Covariance between b-d and between BMD-d
      E4covND = c(cov(optE4_ND$theta_tilde[,c("b","d")], use="complete.obs")["b","d"],
                  cov(optE4_ND$theta_tilde[,c("BMD","d")], use="complete.obs")["BMD","d"])

      E4corrND = c(cor(optE4_ND$theta_tilde[,c("b","d")], use="complete.obs")["b","d"],
                   cor(optE4_ND$theta_tilde[,c("BMD","d")], use="complete.obs")["BMD","d"])

      DRM_E4_N = DRM.E4_ND(E4resND[4:7], data.modstan1a$x, data.modstan1a$q)

      # obtain loglikelihood
      llE4N=llfE4_ND(optE4_ND$par[1:5],
                     nvec=data.modstan1a$n,
                     dvec=data.modstan1a$x,
                     mvec=data.modstan1a$m,
                     s2vec=data.modstan1a$s2,
                     qval=data.modstan1a$q)

      # save LL and BMD(L/U) estimate
      llN = c(llN, llE4N)
      BMDL = c(BMDL, E4resND[1])
      BMD = c(BMD, E4resND[2])
      BMDU = c(BMDU, E4resND[3])
    }else{E4resND=NULL; llN = c(llN,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); E4covND=rep(NA,2); E4corrND=rep(NA,2); DRM_E4_N=rep(NA, length(data.modstan1a$x))
    parsE4N <- NA
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

    #

    if(prior.weights[2]>0){

      print(2)

      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          optIE4_ND <- fun_optim(stanmodels$mIE4_ND, data.modstan1a, svF1, ndraws, 123, pvec)
        }else if(prior.BMD == TRUE){
          optIE4_ND <- fun_optim(stanmodels$mIE4_ND2, data.modstan1a, svF1, ndraws, 123, pvec)
        }
      }else if(priordist == "Normal"){
        optIE4_ND <- fun_optim(stanmodels$mIE4_ND_n, data.modstan1a, svF1, ndraws, 123, pvec)
      }

      parsIE4N <- par_extract(optIE4_ND, model_name = "IE4_N")
      IE4resND <- quantile(parsIE4N$BMD, pvec)  #exp(quantile(optE4_NI$theta_tilde[,"par[2]"],pvec))
      IE4resND <- c(IE4resND,apply(parsIE4N[,c(paste0("p",1:4), "is2t")], 2, median))
      names(IE4resND) <- c("BMDL","BMD","BMDU","at","k","ct","dt","is2t")

      IE4outND <- outLP(parsIE4N, pvec, data.modstan1a$maxD)

      # Covariance between b-d and between BMD-d
      IE4covND = c(cov(optIE4_ND$theta_tilde[,c("b","d")], use="complete.obs")["b","d"],
                   cov(optIE4_ND$theta_tilde[,c("BMD","d")], use="complete.obs")["BMD","d"])

      IE4corrND = c(cor(optIE4_ND$theta_tilde[,c("b","d")], use="complete.obs")["b","d"],
                    cor(optIE4_ND$theta_tilde[,c("BMD","d")], use="complete.obs")["BMD","d"])

      DRM_IE4_N = DRM.IE4_ND(IE4resND[4:7], data.modstan1a$x, data.modstan1a$q)


      llIE4N=llfIE4_ND(optIE4_ND$par[1:5],nvec=data.modstan1a$n,
                       dvec=data.modstan1a$x,mvec=data.modstan1a$m,
                       s2vec=data.modstan1a$s2,qval=data.modstan1a$q)

      llN = c(llN, llIE4N)
      BMDL = c(BMDL, IE4resND[1])
      BMD = c(BMD, IE4resND[2])
      BMDU = c(BMDU, IE4resND[3])
    }else{
      IE4resND=NULL; llN = c(llN,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); IE4covND=rep(NA,2); IE4corrND=rep(NA,2); DRM_IE4_N=rep(NA,length(data.modstan1a$x))
      parsIE4N <- NA
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
    #

    if(prior.weights[3]>0){

      print(3)

      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          optH4_ND <- fun_optim(stanmodels$mH4_ND, data.modstan1a, svF1, ndraws, 123, pvec)
        }else if(prior.BMD == TRUE){
          optH4_ND <- fun_optim(stanmodels$mH4_ND2, data.modstan1a, svF1, ndraws, 123, pvec)
        }
      }else if(priordist == "Normal"){
        optH4_ND <- fun_optim(stanmodels$mH4_ND_n, data.modstan1a, svF1, ndraws, 123, pvec)
      }

      parsH4N <- par_extract(optH4_ND, model_name = "H4_N")
      H4resND <- quantile(parsH4N$BMD, pvec)  #exp(quantile(optH4_NI$theta_tilde[,"par[2]"],pvec))
      H4resND <- c(H4resND,apply(parsH4N[,c(paste0("p",1:4), "is2t")], 2, median))
      names(H4resND) <- c("BMDL","BMD","BMDU","at","k","ct","dt","is2t")

      H4outND <- outLP(parsH4N, pvec, data.modstan1a$maxD)

      # Covariance between b-d and between BMD-d
      H4covND = c(cov(optH4_ND$theta_tilde[,c("b","d")], use="complete.obs")["b","d"],
                  cov(optH4_ND$theta_tilde[,c("BMD","d")], use="complete.obs")["BMD","d"])

      H4corrND = c(cor(optH4_ND$theta_tilde[,c("b","d")], use="complete.obs")["b","d"],
                   cor(optH4_ND$theta_tilde[,c("BMD","d")], use="complete.obs")["BMD","d"])

      DRM_H4_N = DRM.H4_ND(H4resND[4:7], data.modstan1a$x, data.modstan1a$q)


      llH4N=llfH4_ND(optH4_ND$par[1:5],nvec=data.modstan1a$n,
                     dvec=data.modstan1a$x,mvec=data.modstan1a$m,
                     s2vec=data.modstan1a$s2,qval=data.modstan1a$q)

      llN = c(llN, llH4N)
      BMDL = c(BMDL, H4resND[1])
      BMD = c(BMD, H4resND[2])
      BMDU = c(BMDU, H4resND[3])
    }else{H4resND=NULL; llN = c(llN,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); H4covND=rep(NA,2); H4corrND=rep(NA,2); DRM_H4_N=rep(NA,length(data.modstan1a$x))
    parsH4N <- NA
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

    #

    if(prior.weights[4]>0){

      print(4)

      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          optLN4_ND <- fun_optim(stanmodels$mLN4_ND, data.modstan1aLN, svF1LN, ndraws, 123, pvec)
        }else if(prior.BMD == TRUE){
          optLN4_ND <- fun_optim(stanmodels$mLN4_ND2, data.modstan1aLN, svF1LN, ndraws, 123, pvec)
        }
      }else if(priordist == "Normal"){
        optLN4_ND <- fun_optim(stanmodels$mLN4_ND_n, data.modstan1aLN, svF1LN, ndraws, 123, pvec)
      }

      parsLN4N <- par_extract(optLN4_ND, model_name = "LN4_N")
      LN4resND <- quantile(parsLN4N$BMD, pvec)  #exp(quantile(optLN4_NI$theta_tilde[,"par[2]"],pvec))
      LN4resND <- c(LN4resND,apply(parsLN4N[,c(paste0("p",1:4), "is2t")], 2, median))
      names(LN4resND) <- c("BMDL","BMD","BMDU","at","k","ct","dt","is2t")

      LN4outND <- outLP(parsLN4N, pvec, data.modstan1a$maxD)

      # Covariance between b-d and between BMD-d
      LN4covND = c(cov(optLN4_ND$theta_tilde[,c("b","d")], use="complete.obs")["b","d"],
                   cov(optLN4_ND$theta_tilde[,c("BMD","d")], use="complete.obs")["BMD","d"])

      LN4corrND = c(cor(optLN4_ND$theta_tilde[,c("b","d")], use="complete.obs")["b","d"],
                    cor(optLN4_ND$theta_tilde[,c("BMD","d")], use="complete.obs")["BMD","d"])

      DRM_LN4_N = DRM.LN4_ND(LN4resND[4:7], data.modstan1a$x, data.modstan1a$q)


      llLN4N=llfLN4_ND(optLN4_ND$par[1:5],nvec=data.modstan1aLN$n,
                       dvec=data.modstan1aLN$x,mvec=data.modstan1aLN$m,
                       s2vec=data.modstan1aLN$s2,qval=data.modstan1aLN$q)

      llN = c(llN, llLN4N)
      BMDL = c(BMDL, LN4resND[1])
      BMD = c(BMD, LN4resND[2])
      BMDU = c(BMDU, LN4resND[3])
    }else{LN4resND=NULL; llN = c(llN,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); LN4covND=rep(NA,2); LN4corrND=rep(NA,2); DRM_LN4_N=rep(NA,length(data.modstan1aLN$x))
    parsLN4N <- NA
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
    #

    if(prior.weights[5]>0){

      print(5)

      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          optG4_ND <- fun_optim(stanmodels$mG4_ND, data.modstan1bG, svF1, ndraws, 123, pvec)
        }else if(prior.BMD == TRUE){
          optG4_ND <- fun_optim(stanmodels$mG4_ND2, data.modstan1bG, svF1, ndraws, 123, pvec)
        }
      }else if(priordist == "Normal"){
        optG4_ND <- fun_optim(stanmodels$mG4_ND_n, data.modstan1bG, svF1, ndraws, 123, pvec)
      }

      parsG4N <- par_extract(optG4_ND, model_name = "G4_N")
      G4resND <- quantile(parsG4N$BMD, pvec)  #exp(quantile(optE4_NI$theta_tilde[,"par[2]"],pvec))
      G4resND <- c(G4resND,apply(parsG4N[,c(paste0("p",1:4), "is2t")], 2, median))
      names(G4resND) <- c("BMDL","BMD","BMDU","at","k","ct","dt","is2t")

      G4outND <- outLP(parsG4N, pvec, data.modstan1a$maxD)

      # Covariance between b-d and between BMD-d
      G4covND = c(cov(optG4_ND$theta_tilde[,c("b","d")], use="complete.obs")["b","d"],
                  cov(optG4_ND$theta_tilde[,c("BMD","d")], use="complete.obs")["BMD","d"])

      G4corrND = c(cor(optG4_ND$theta_tilde[,c("b","d")], use="complete.obs")["b","d"],
                   cor(optG4_ND$theta_tilde[,c("BMD","d")], use="complete.obs")["BMD","d"])

      DRM_G4_N = DRM.G4_ND(G4resND[4:7], data.modstan1a$x, data.modstan1a$q)


      llG4N=llfG4_ND(optG4_ND$par[1:5],nvec=data.modstan1bG$n,
                     dvec=data.modstan1bG$x,mvec=data.modstan1bG$m,
                     s2vec=data.modstan1bG$s2,qval=data.modstan1bG$q)

      llN = c(llN, llG4N)
      BMDL = c(BMDL, G4resND[1])
      BMD = c(BMD, G4resND[2])
      BMDU = c(BMDU, G4resND[3])
    }else{G4resND=NULL; llN = c(llN,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); G4covND=rep(NA,2); G4corrND=rep(NA,2); DRM_G4_N=rep(NA, length(data.modstan1bG$x))
    parsG4N <- NA
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
    #

    if(prior.weights[6]>0){

      print(6)

      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          optQE4_ND <- fun_optim(stanmodels$mQE4_ND, data.modstan1bQ, svF1, ndraws, 123, pvec)
        }else if(prior.BMD == TRUE){
          optQE4_ND <- fun_optim(stanmodels$mQE4_ND2, data.modstan1bQ, svF1, ndraws, 123, pvec)
        }
      }else if(priordist == "Normal"){
        optQE4_ND <- fun_optim(stanmodels$mQE4_ND_n, data.modstan1bQ, svF1, ndraws, 123, pvec)
      }

      parsQE4N <- par_extract(optQE4_ND, model_name = "QE4_N")
      QE4resND <- quantile(parsQE4N$BMD, pvec)  #exp(quantile(optE4_NI$theta_tilde[,"par[2]"],pvec))
      QE4resND <- c(QE4resND,apply(parsQE4N[,c(paste0("p",1:4), "is2t")], 2, median))
      names(QE4resND) <- c("BMDL","BMD","BMDU","at","k","ct","dt","is2t")

      QE4outND <- outLP(parsQE4N, pvec, data.modstan1a$maxD)

      # Covariance between b-d and between BMD-d
      QE4covND = c(cov(optQE4_ND$theta_tilde[,c("b","d")], use="complete.obs")["b","d"],
                   cov(optQE4_ND$theta_tilde[,c("BMD","d")], use="complete.obs")["BMD","d"])

      QE4corrND = c(cor(optQE4_ND$theta_tilde[,c("b","d")], use="complete.obs")["b","d"],
                    cor(optQE4_ND$theta_tilde[,c("BMD","d")], use="complete.obs")["BMD","d"])

      DRM_QE4_N = DRM.QE4_ND(QE4resND[4:7], data.modstan1a$x, data.modstan1a$q)


      llQE4N=llfQE4_ND(optQE4_ND$par[1:5],nvec=data.modstan1bQ$n,
                       dvec=data.modstan1bQ$x,mvec=data.modstan1bQ$m,
                       s2vec=data.modstan1bQ$s2,qval=data.modstan1bQ$q)

      llN = c(llN, llQE4N)
      BMDL = c(BMDL, QE4resND[1])
      BMD = c(BMD, QE4resND[2])
      BMDU = c(BMDU, QE4resND[3])
    }else{QE4resND=NULL; llN = c(llN,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); QE4covND=rep(NA,2); QE4corrND=rep(NA,2); DRM_QE4_N=rep(NA,length(data.modstan1bQ$x))
    parsQE4N <- NA
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
    #

    if(prior.weights[7]>0){

      print(7)

      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          optP4_ND <- fun_optim(stanmodels$mP4_ND, data.modstan2P, svF2P, ndraws, 123, pvec)
        }else if(prior.BMD == TRUE){
          optP4_ND <- fun_optim(stanmodels$mP4_ND2, data.modstan2P, svF2P, ndraws, 123, pvec)
        }
      }else if(priordist == "Normal"){
        optP4_ND <- fun_optim(stanmodels$mP4_ND_n, data.modstan2P, svF2P, ndraws, 123, pvec)
      }

      parsP4N <- par_extract(optP4_ND, model_name = "P4_N")
      P4resND <- quantile(parsP4N$BMD, pvec)  #exp(quantile(optE4_NI$theta_tilde[,"par[2]"],pvec))
      P4resND <- c(P4resND,apply(parsP4N[,c(paste0("p",1:4), "is2t")], 2, median))
      names(P4resND) <- c("BMDL","BMD","BMDU","at","k","ct","dt","is2t")

      P4outND <- outLP(parsP4N, pvec, data.modstan1a$maxD)

      # Covariance between b-d and between BMD-d
      P4covND = c(cov(optP4_ND$theta_tilde[,c("b","d")], use="complete.obs")["b","d"],
                  cov(optP4_ND$theta_tilde[,c("BMD","d")], use="complete.obs")["BMD","d"])

      P4corrND = c(cor(optP4_ND$theta_tilde[,c("b","d")], use="complete.obs")["b","d"],
                   cor(optP4_ND$theta_tilde[,c("BMD","d")], use="complete.obs")["BMD","d"])

      DRM_P4_N = DRM.P4_ND(P4resND[4:7], data.modstan1a$x, data.modstan1a$q)


      llP4N=llfP4_ND(optP4_ND$par[1:5],nvec=data.modstan2P$n,
                     dvec=data.modstan2P$x,mvec=data.modstan2P$m,
                     s2vec=data.modstan2P$s2,qval=data.modstan2P$q)

      llN = c(llN, llP4N)
      BMDL = c(BMDL, P4resND[1])
      BMD = c(BMD, P4resND[2])
      BMDU = c(BMDU, P4resND[3])
    }else{P4resND=NULL; llN = c(llN,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); P4covND=rep(NA,2); P4corrND=rep(NA,2); DRM_P4_N=rep(NA,length(data.modstan2P$x))
    parsP4N <- NA
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
    #

    if(prior.weights[8]>0){

      print(8)

      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          optL4_ND <- fun_optim(stanmodels$mL4_ND, data.modstan2L, svF2L, ndraws, 123, pvec)
        }else if(prior.BMD == TRUE){
          optL4_ND <- fun_optim(stanmodels$mL4_ND2, data.modstan2L, svF2L, ndraws, 123, pvec)
        }
      }else if(priordist == "Normal"){
        optL4_ND <- fun_optim(stanmodels$mL4_ND_n, data.modstan2L, svF2L, ndraws, 123, pvec)
      }

      parsL4N <- par_extract(optL4_ND, model_name = "L4_N")
      L4resND <- quantile(parsL4N$BMD, pvec)  #exp(quantile(optE4_NI$theta_tilde[,"par[2]"],pvec))
      L4resND <- c(L4resND,apply(parsL4N[,c(paste0("p",1:4), "is2t")], 2, median))
      names(L4resND) <- c("BMDL","BMD","BMDU","at","k","ct","dt","is2t")

      L4outND <- outLP(parsL4N, pvec, data.modstan1a$maxD)

      # Covariance between b-d and between BMD-d
      L4covND = c(cov(optL4_ND$theta_tilde[,c("b","d")], use="complete.obs")["b","d"],
                  cov(optL4_ND$theta_tilde[,c("BMD","d")], use="complete.obs")["BMD","d"])

      L4corrND = c(cor(optL4_ND$theta_tilde[,c("b","d")], use="complete.obs")["b","d"],
                   cor(optL4_ND$theta_tilde[,c("BMD","d")], use="complete.obs")["BMD","d"])

      DRM_L4_N = DRM.L4_ND(L4resND[4:7], data.modstan1a$x, data.modstan1a$q)


      llL4N=llfL4_ND(optL4_ND$par[1:5],nvec=data.modstan2L$n,
                     dvec=data.modstan2L$x,mvec=data.modstan2L$m,
                     s2vec=data.modstan2L$s2,qval=data.modstan2L$q)

      llN = c(llN, llL4N)
      BMDL = c(BMDL, L4resND[1])
      BMD = c(BMD, L4resND[2])
      BMDU = c(BMDU, L4resND[3])
    }else{L4resND=NULL; llN = c(llN,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); L4covND=rep(NA,2); L4corrND=rep(NA,2); DRM_L4_N=rep(NA,length(data.modstan2L$x))
    parsL4N <- NA
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


    ## Data to use for Lognormal distribution

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

      print(9)

      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          optE4_LND <- fun_optim(stanmodels$mE4_LND, data.modstan1a, svF1, ndraws, 123, pvec)
        }else if(prior.BMD == TRUE){
          optE4_LND <- fun_optim(stanmodels$mE4_LND2, data.modstan1a, svF1, ndraws, 123, pvec)
        }
      }else if(priordist == "Normal"){
        optE4_LND <- fun_optim(stanmodels$mE4_LND_n, data.modstan1a, svF1, ndraws, 123, pvec)
      }

      parsE4LN <- par_extract(optE4_LND, model_name = "E4_LN")
      E4resLND <- quantile(parsE4LN$BMD, pvec)  #exp(quantile(optE4_NI$theta_tilde[,"par[2]"],pvec))
      E4resLND <- c(E4resLND,apply(parsE4LN[,c(paste0("p",1:4), "is2t")], 2, median))
      names(E4resLND) <- c("BMDL","BMD","BMDU","at","k","ct","dt","is2t")

      E4outLND <- outLP(parsE4LN, pvec, data.modstan1a$maxD)

      # Covariance between b-d and between BMD-d
      E4covLND = c(cov(optE4_LND$theta_tilde[,c("b","d")], use="complete.obs")["b","d"],
                   cov(optE4_LND$theta_tilde[,c("BMD","d")], use="complete.obs")["BMD","d"])

      E4corrLND = c(cor(optE4_LND$theta_tilde[,c("b","d")], use="complete.obs")["b","d"],
                    cor(optE4_LND$theta_tilde[,c("BMD","d")], use="complete.obs")["BMD","d"])

      DRM_E4_LN = exp(DRM.E4_LND(E4resLND[4:7], data.modstan1a$x, data.modstan1a$q))

      llE4LN=llfE4_LND(optE4_LND$par[1:5],nvec=data.modstan1a$n,
                       dvec=data.modstan1a$x,mvec=data.modstan1a$m,
                       s2vec=data.modstan1a$s2,qval=data.modstan1a$q)

      llLN = c(llLN, llE4LN)
      BMDL = c(BMDL, E4resLND[1])
      BMD = c(BMD, E4resLND[2])
      BMDU = c(BMDU, E4resLND[3])
    }else{E4resLND=NULL; llLN = c(llLN,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); E4covLND=rep(NA,2); E4corrLND=rep(NA,2); DRM_E4_LN=rep(NA,length(data.modstan1a$x))
    parsE4LN <- NA
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
    #

    if(prior.weights[10]>0){

      print(10)

      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          optIE4_LND <- fun_optim(stanmodels$mIE4_LND, data.modstan1a, svF1, ndraws, 123, pvec)
        }else if(prior.BMD == TRUE){
          optIE4_LND <- fun_optim(stanmodels$mIE4_LND2, data.modstan1a, svF1, ndraws, 123, pvec)
        }
        }else if(priordist == "Normal"){
        optIE4_LND <- fun_optim(stanmodels$mIE4_LND_n, data.modstan1a, svF1, ndraws, 123, pvec)
      }

      parsIE4LN <- par_extract(optIE4_LND, model_name = "IE4_LN")
      IE4resLND <- quantile(parsIE4LN$BMD, pvec)  #exp(quantile(optIE4_NI$theta_tilde[,"par[2]"],pvec))
      IE4resLND <- c(IE4resLND,apply(parsIE4LN[,c(paste0("p",1:4), "is2t")], 2, median))
      names(IE4resLND) <- c("BMDL","BMD","BMDU","at","k","ct","dt","is2t")

      IE4outLND <- outLP(parsIE4LN, pvec, data.modstan1a$maxD)

      # Covariance between b-d and between BMD-d
      IE4covLND = c(cov(optIE4_LND$theta_tilde[,c("b","d")], use="complete.obs")["b","d"],
                    cov(optIE4_LND$theta_tilde[,c("BMD","d")], use="complete.obs")["BMD","d"])

      IE4corrLND = c(cor(optIE4_LND$theta_tilde[,c("b","d")], use="complete.obs")["b","d"],
                     cor(optIE4_LND$theta_tilde[,c("BMD","d")], use="complete.obs")["BMD","d"])

      DRM_IE4_LN = exp(DRM.IE4_LND(IE4resLND[4:7], data.modstan1a$x, data.modstan1a$q))


      llIE4LN=llfIE4_LND(optIE4_LND$par[1:5],nvec=data.modstan1a$n,
                         dvec=data.modstan1a$x,mvec=data.modstan1a$m,
                         s2vec=data.modstan1a$s2,qval=data.modstan1a$q)

      llLN = c(llLN, llIE4LN)
      BMDL = c(BMDL, IE4resLND[1])
      BMD = c(BMD, IE4resLND[2])
      BMDU = c(BMDU, IE4resLND[3])
    }else{IE4resLND=NULL; llLN = c(llLN,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); IE4covLND=rep(NA,2); IE4corrLND=rep(NA,2); DRM_IE4_LN=rep(NA,length(data.modstan1a$x))
    parsIE4LN <- NA
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
    #

    if(prior.weights[11]>0){

      print(11)

      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          optH4_LND <- fun_optim(stanmodels$mH4_LND, data.modstan1a, svF1, ndraws, 123, pvec)
        }else if(prior.BMD == TRUE){
          optH4_LND <- fun_optim(stanmodels$mH4_LND2, data.modstan1a, svF1, ndraws, 123, pvec)
        }
        }else if(priordist == "Normal"){
        optH4_LND <- fun_optim(stanmodels$mH4_LND_n, data.modstan1a, svF1, ndraws, 123, pvec)
      }

      parsH4LN <- par_extract(optH4_LND, model_name = "H4_LN")
      H4resLND <- quantile(parsH4LN$BMD, pvec)  #exp(quantile(optH4_NI$theta_tilde[,"par[2]"],pvec))
      H4resLND <- c(H4resLND,apply(parsH4LN[,c(paste0("p",1:4), "is2t")], 2, median))
      names(H4resLND) <- c("BMDL","BMD","BMDU","at","k","ct","dt","is2t")

      H4outLND <- outLP(parsH4LN, pvec, data.modstan1a$maxD)

      # Covariance between b-d and between BMD-d
      H4covLND = c(cov(optH4_LND$theta_tilde[,c("b","d")], use="complete.obs")["b","d"],
                   cov(optH4_LND$theta_tilde[,c("BMD","d")], use="complete.obs")["BMD","d"])

      H4corrLND = c(cor(optH4_LND$theta_tilde[,c("b","d")], use="complete.obs")["b","d"],
                    cor(optH4_LND$theta_tilde[,c("BMD","d")], use="complete.obs")["BMD","d"])

      DRM_H4_LN = exp(DRM.H4_LND(H4resLND[4:7], data.modstan1a$x, data.modstan1a$q))


      llH4LN=llfH4_LND(optH4_LND$par[1:5],nvec=data.modstan1a$n,
                       dvec=data.modstan1a$x,mvec=data.modstan1a$m,
                       s2vec=data.modstan1a$s2,qval=data.modstan1a$q)

      llLN = c(llLN, llH4LN)
      BMDL = c(BMDL, H4resLND[1])
      BMD = c(BMD, H4resLND[2])
      BMDU = c(BMDU, H4resLND[3])
    }else{H4resLND=NULL; llLN = c(llLN,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); H4covLND=rep(NA,2); H4corrLND=rep(NA,2); DRM_H4_LN=rep(NA,length(data.modstan1a$x))
    parsH4LN <- NA
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
    #

    if(prior.weights[12]>0){

      print(12)

      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          optLN4_LND <- fun_optim(stanmodels$mLN4_LND, data.modstan1aLN, svF1LN, ndraws, 123, pvec)
        }else if(prior.BMD == TRUE){
          optLN4_LND <- fun_optim(stanmodels$mLN4_LND2, data.modstan1aLN, svF1LN, ndraws, 123, pvec)
        }
        }else if(priordist == "Normal"){
        optLN4_LND <- fun_optim(stanmodels$mLN4_LND_n, data.modstan1aLN, svF1LN, ndraws, 123, pvec)
      }

      parsLN4LN <- par_extract(optLN4_LND, model_name = "LN4_LN")
      LN4resLND <- quantile(parsLN4LN$BMD, pvec)  #exp(quantile(optLN4_NI$theta_tilde[,"par[2]"],pvec))
      LN4resLND <- c(LN4resLND,apply(parsLN4LN[,c(paste0("p",1:4), "is2t")], 2, median))
      names(LN4resLND) <- c("BMDL","BMD","BMDU","at","k","ct","dt","is2t")

      LN4outLND <- outLP(parsLN4LN, pvec, data.modstan1a$maxD)

      # Covariance between b-d and between BMD-d
      LN4covLND = c(cov(optLN4_LND$theta_tilde[,c("b","d")], use="complete.obs")["b","d"],
                    cov(optLN4_LND$theta_tilde[,c("BMD","d")], use="complete.obs")["BMD","d"])

      LN4corrLND = c(cor(optLN4_LND$theta_tilde[,c("b","d")], use="complete.obs")["b","d"],
                     cor(optLN4_LND$theta_tilde[,c("BMD","d")], use="complete.obs")["BMD","d"])

      DRM_LN4_LN = exp(DRM.LN4_LND(LN4resLND[4:7], data.modstan1a$x, data.modstan1a$q))


      llLN4LN=llfLN4_LND(optLN4_LND$par[1:5],nvec=data.modstan1aLN$n,
                         dvec=data.modstan1aLN$x,mvec=data.modstan1aLN$m,
                         s2vec=data.modstan1aLN$s2,qval=data.modstan1aLN$q)

      llLN = c(llLN, llLN4LN)
      BMDL = c(BMDL, LN4resLND[1])
      BMD = c(BMD, LN4resLND[2])
      BMDU = c(BMDU, LN4resLND[3])
    }else{LN4resLND=NULL; llLN = c(llLN,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); LN4covLND=rep(NA,2); LN4corrLND=rep(NA,2); DRM_LN4_LN=rep(NA,length(data.modstan1aLN$x))
    parsLN4LN <- NA
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
    #

    if(prior.weights[13]>0){

      print(13)

      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          optG4_LND <- fun_optim(stanmodels$mG4_LND, data.modstan1bG, svF1, ndraws, 123, pvec)
        }else if(prior.BMD == TRUE){
          optG4_LND <- fun_optim(stanmodels$mG4_LND2, data.modstan1bG, svF1, ndraws, 123, pvec)
        }
        }else if(priordist == "Normal"){
        optG4_LND <- fun_optim(stanmodels$mG4_LND_n, data.modstan1bG, svF1, ndraws, 123, pvec)
      }

      parsG4LN <- par_extract(optG4_LND, model_name = "G4_LN")
      G4resLND <- quantile(parsG4LN$BMD, pvec)  #exp(quantile(optG4_NI$theta_tilde[,"par[2]"],pvec))
      G4resLND <- c(G4resLND,apply(parsG4LN[,c(paste0("p",1:4), "is2t")], 2, median))
      names(G4resLND) <- c("BMDL","BMD","BMDU","at","k","ct","dt","is2t")

      G4outLND <- outLP(parsG4LN, pvec, data.modstan1a$maxD)

      # Covariance between b-d and between BMD-d
      G4covLND = c(cov(optG4_LND$theta_tilde[,c("b","d")], use="complete.obs")["b","d"],
                   cov(optG4_LND$theta_tilde[,c("BMD","d")], use="complete.obs")["BMD","d"])

      G4corrLND = c(cor(optG4_LND$theta_tilde[,c("b","d")], use="complete.obs")["b","d"],
                    cor(optG4_LND$theta_tilde[,c("BMD","d")], use="complete.obs")["BMD","d"])

      DRM_G4_LN = exp(DRM.G4_LND(G4resLND[4:7], data.modstan1a$x, data.modstan1a$q))


      llG4LN=llfG4_LND(optG4_LND$par[1:5],nvec=data.modstan1bG$n,
                       dvec=data.modstan1bG$x,mvec=data.modstan1bG$m,
                       s2vec=data.modstan1bG$s2,qval=data.modstan1bG$q)

      llLN = c(llLN, llG4LN)
      BMDL = c(BMDL, G4resLND[1])
      BMD = c(BMD, G4resLND[2])
      BMDU = c(BMDU, G4resLND[3])
    }else{G4resLND=NULL; llLN = c(llLN,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); G4covLND=rep(NA,2); G4corrLND=rep(NA,2); DRM_G4_LN=rep(NA,length(data.modstan1bG$x))
    parsG4LN <- NA
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
    #

    if(prior.weights[14]>0){

      print(14)

      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          optQE4_LND <- fun_optim(stanmodels$mQE4_LND, data.modstan1bQ, svF1, ndraws, 123, pvec)
        }else if(prior.BMD == TRUE){
          optQE4_LND <- fun_optim(stanmodels$mQE4_LND2, data.modstan1bQ, svF1, ndraws, 123, pvec)
        }
        }else if(priordist == "Normal"){
        optQE4_LND <- fun_optim(stanmodels$mQE4_LND_n, data.modstan1bQ, svF1, ndraws, 123, pvec)
      }

      parsQE4LN <- par_extract(optQE4_LND, model_name = "QE4_LN")
      QE4resLND <- quantile(parsQE4LN$BMD, pvec)  #exp(quantile(optQE4_NI$theta_tilde[,"par[2]"],pvec))
      QE4resLND <- c(QE4resLND,apply(parsQE4LN[,c(paste0("p",1:4), "is2t")], 2, median))
      names(QE4resLND) <- c("BMDL","BMD","BMDU","at","k","ct","dt","is2t")

      QE4outLND <- outLP(parsQE4LN, pvec, data.modstan1a$maxD)

      # Covariance between b-d and between BMD-d
      QE4covLND = c(cov(optQE4_LND$theta_tilde[,c("b","d")], use="complete.obs")["b","d"],
                    cov(optQE4_LND$theta_tilde[,c("BMD","d")], use="complete.obs")["BMD","d"])

      QE4corrLND = c(cor(optQE4_LND$theta_tilde[,c("b","d")], use="complete.obs")["b","d"],
                     cor(optQE4_LND$theta_tilde[,c("BMD","d")], use="complete.obs")["BMD","d"])

      DRM_QE4_LN = exp(DRM.QE4_LND(QE4resLND[4:7], data.modstan1a$x, data.modstan1a$q))


      llQE4LN=llfQE4_LND(optQE4_LND$par[1:5],nvec=data.modstan1bQ$n,
                         dvec=data.modstan1bQ$x,mvec=data.modstan1bQ$m,
                         s2vec=data.modstan1bQ$s2,qval=data.modstan1bQ$q)

      llLN = c(llLN, llQE4LN)
      BMDL = c(BMDL, QE4resLND[1])
      BMD = c(BMD, QE4resLND[2])
      BMDU = c(BMDU, QE4resLND[3])
    }else{QE4resLND=NULL; llLN = c(llLN,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); QE4covLND=rep(NA,2); QE4corrLND=rep(NA,2); DRM_QE4_LN=rep(NA,length(data.modstan1bQ$x))
    parsQE4LN <- NA
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
    #

    if(prior.weights[15]>0){

      print(15)

      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          optP4_LND <- fun_optim(stanmodels$mP4_LND, data.modstan2P, svF2P, ndraws, 123, pvec)
        }else if(prior.BMD == TRUE){
          optP4_LND <- fun_optim(stanmodels$mP4_LND2, data.modstan2P, svF2P, ndraws, 123, pvec)
        }
        }else if(priordist == "Normal"){
        optP4_LND <- fun_optim(stanmodels$mP4_LND_n, data.modstan2P, svF2P, ndraws, 123, pvec)
      }

      parsP4LN <- par_extract(optP4_LND, model_name = "P4_LN")
      P4resLND <- quantile(parsP4LN$BMD, pvec)  #exp(quantile(optP4_NI$theta_tilde[,"par[2]"],pvec))
      P4resLND <- c(P4resLND,apply(parsP4LN[,c(paste0("p",1:4), "is2t")], 2, median))
      names(P4resLND) <- c("BMDL","BMD","BMDU","at","k","ct","dt","is2t")

      P4outLND <- outLP(parsP4LN, pvec, data.modstan1a$maxD)

      # Covariance between b-d and between BMD-d
      P4covLND = c(cov(optP4_LND$theta_tilde[,c("b","d")], use="complete.obs")["b","d"],
                   cov(optP4_LND$theta_tilde[,c("BMD","d")], use="complete.obs")["BMD","d"])

      P4corrLND = c(cor(optP4_LND$theta_tilde[,c("b","d")], use="complete.obs")["b","d"],
                    cor(optP4_LND$theta_tilde[,c("BMD","d")], use="complete.obs")["BMD","d"])

      DRM_P4_LN = exp(DRM.P4_LND(P4resLND[4:7], data.modstan1a$x, data.modstan1a$q))



      llP4LN=llfP4_LND(optP4_LND$par[1:5],nvec=data.modstan2P$n,
                       dvec=data.modstan2P$x,mvec=data.modstan2P$m,
                       s2vec=data.modstan2P$s2,qval=data.modstan2P$q)

      llLN = c(llLN, llP4LN)
      BMDL = c(BMDL, P4resLND[1])
      BMD = c(BMD, P4resLND[2])
      BMDU = c(BMDU, P4resLND[3])
    }else{P4resLND=NULL; llLN = c(llLN,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); P4covLND=rep(NA,2); P4corrLND=rep(NA,2); DRM_P4_LN=rep(NA,length(data.modstan2P$x))
    parsP4LN <- NA
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
    #

    if(prior.weights[16]>0){

      print(16)

      if(priordist == "PERT"){
        if(prior.BMD == FALSE){
          optL4_LND <- fun_optim(stanmodels$mL4_LND, data.modstan2L, svF2L, ndraws, 123, pvec)
        }else if(prior.BMD == TRUE){
          optL4_LND <- fun_optim(stanmodels$mL4_LND2, data.modstan2L, svF2L, ndraws, 123, pvec)
        }
        }else if(priordist == "Normal"){
        optL4_LND <- fun_optim(stanmodels$mL4_LND_n, data.modstan2L, svF2L, ndraws, 123, pvec)
      }

      parsL4LN <- par_extract(optL4_LND, model_name = "L4_LN")
      L4resLND <- quantile(parsL4LN$BMD, pvec)  #exp(quantile(optL4_NI$theta_tilde[,"par[2]"],pvec))
      L4resLND <- c(L4resLND,apply(parsL4LN[,c(paste0("p",1:4), "is2t")], 2, median))
      names(L4resLND) <- c("BMDL","BMD","BMDU","at","k","ct","dt","is2t")

      L4outLND <- outLP(parsL4LN, pvec, data.modstan1a$maxD)

      # Covariance between b-d and between BMD-d
      L4covLND = c(cov(optL4_LND$theta_tilde[,c("b","d")], use="complete.obs")["b","d"],
                   cov(optL4_LND$theta_tilde[,c("BMD","d")], use="complete.obs")["BMD","d"])

      L4corrLND = c(cor(optL4_LND$theta_tilde[,c("b","d")], use="complete.obs")["b","d"],
                    cor(optL4_LND$theta_tilde[,c("BMD","d")], use="complete.obs")["BMD","d"])

      DRM_L4_LN = exp(DRM.L4_LND(L4resLND[4:7], data.modstan1a$x, data.modstan1a$q))


      llL4LN=llfL4_LND(optL4_LND$par[1:5],nvec=data.modstan2L$n,
                       dvec=data.modstan2L$x,mvec=data.modstan2L$m,
                       s2vec=data.modstan2L$s2,qval=data.modstan2L$q)

      llLN = c(llLN, llL4LN)
      BMDL = c(BMDL, L4resLND[1])
      BMD = c(BMD, L4resLND[2])
      BMDU = c(BMDU, L4resLND[3])
    }else{L4resLND=NULL; llLN = c(llLN,NA); BMDL = c(BMDL,NA); BMD = c(BMD,NA); BMDU=c(BMDU,NA); L4covLND=rep(NA,2); L4corrLND=rep(NA,2); DRM_L4_LN=rep(NA,length(data.modstan2L$x))
    parsL4LN <- NA
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

    # the log-likelihoods functions
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

    # print(llN); print(llLN); print(w)

    prior.weights = prior.weights/sum(prior.weights==1) # in case this has changed for G4
    lpw=(prior.weights*w)/sum(prior.weights*w)

    # the model average posterior as a mixture
    count=round(lpw*ndraws)
    mabmd=(c(# normal
      if(prior.weights[1]>0) sample(optE4_ND$theta_tilde[,2],count[1],replace=T),
      if(prior.weights[2]>0) sample(optIE4_ND$theta_tilde[,2],count[2],replace=T),
      if(prior.weights[3]>0) sample(optH4_ND$theta_tilde[,2],count[3],replace=T),
      if(prior.weights[4]>0) sample(optLN4_ND$theta_tilde[,2],count[4],replace=T),
      if(prior.weights[5]>0) sample(optG4_ND$theta_tilde[,2],count[5],replace=T),
      if(prior.weights[6]>0) sample(optQE4_ND$theta_tilde[,2],count[6],replace=T),
      if(prior.weights[7]>0) sample(optP4_ND$theta_tilde[,2],count[7],replace=T),
      if(prior.weights[8]>0) sample(optL4_ND$theta_tilde[,2],count[8],replace=T),
      # lognormal
      if(prior.weights[9]>0) sample(optE4_LND$theta_tilde[,2],count[9],replace=T),
      if(prior.weights[10]>0) sample(optIE4_LND$theta_tilde[,2],count[10],replace=T),
      if(prior.weights[11]>0) sample(optH4_LND$theta_tilde[,2],count[11],replace=T),
      if(prior.weights[12]>0) sample(optLN4_LND$theta_tilde[,2],count[12],replace=T),
      if(prior.weights[13]>0) sample(optG4_LND$theta_tilde[,2],count[13],replace=T),
      if(prior.weights[14]>0) sample(optQE4_LND$theta_tilde[,2],count[14],replace=T),
      if(prior.weights[15]>0) sample(optP4_LND$theta_tilde[,2],count[15],replace=T),
      if(prior.weights[16]>0) sample(optL4_LND$theta_tilde[,2],count[16],replace=T)
    ))
    maci=exp(quantile(mabmd,pvec))*data.modstan1a$maxD
    names(maci)=c("BMDL","BMD","BMDU")

    BMDq = exp(quantile(mabmd, seq(0,1,0.005)))*data.modstan1a$maxD

    BMDL = c(BMDL, maci[1]/data.modstan1a$maxD); BMD = c(BMD, maci[2]/data.modstan1a$maxD); BMDU = c(BMDU, maci[3]/data.modstan1a$maxD)

    names(BMDL) <- c("E4_N","IE4_N","H4_N","LN4_N","G4_N","QE4_N","P4_N","L4_N","E4_LN","IE4_LN","H4_LN","LN4_LN","G4_LN","QE4_LN","P4_LN","L4_LN","MA")
    model = c("E4_N","IE4_N","H4_N","LN4_N","G4_N","QE4_N","P4_N","L4_N","E4_LN","IE4_LN","H4_LN","LN4_LN","G4_LN","QE4_LN","P4_LN","L4_LN","MA")
    model = as.factor(model)
    weight = c(lpw[1], lpw[2], lpw[3], lpw[4], lpw[5], lpw[6], lpw[7], lpw[8],
               lpw[9], lpw[10], lpw[11], lpw[12], lpw[13], lpw[14], lpw[15], lpw[16], 1)


    names(lpw) = model[1:16]

    ### Model-averaged response per dose level
    dr.MA <- c()
    for(i in 1:length(data.modstan1a$x)){
      dr.MA[i] = weighted.mean(x = c(DRM_E4_N[i],DRM_IE4_N[i],DRM_H4_N[i],DRM_LN4_N[i],DRM_G4_N[i],DRM_QE4_N[i], DRM_P4_N[i], DRM_L4_N[i] ,
                                     DRM_E4_LN[i], DRM_IE4_LN[i], DRM_H4_LN[i], DRM_LN4_LN[i], DRM_G4_LN[i],DRM_QE4_LN[i], DRM_P4_LN[i],DRM_L4_LN[i]),
                               w = lpw,
                               na.rm = T)
    }

    if(plot == TRUE){

      data.modstan1a = data.N$data.modstan1a
      data.modstan1bG = data.N$data.modstan1bG
      data.modstan1aLN = data.N$data.modstan1aLN
      data.modstanE = data.N$data.modstanE
      data.modstan1bQ = data.N$data.modstan1bQ
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
        bmdl = log10(BMDL)
        x[1] = ifelse(min(bmdl,na.rm=T)<x[1], log10(min(BMDL/2, na.rm=T)), log10(data.modstan1a$x[2]/4))
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
          a = exp(parL4[1])*expit(parL4[3]) + 1.5*min(data.LN$data.modstan1a$m.org)}
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
        mtext("Full Laplace")

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

    modelnames = c("E4_N","IE4_N","H4_N","LN4_N","G4_N","QE4_N","P4_N","L4_N","E4_LN","IE4_LN","H4_LN","LN4_LN","G4_LN","QE4_LN","P4_LN","L4_LN")

    ret_results <- list(
      # model parameters
      E4_N=E4outND,IE4_N=IE4outND,H4_N=H4outND,LN4_N=LN4outND,G4_N=G4outND,QE4_N=QE4outND,P4_N=P4outND,L4_N=L4outND,
      E4_LN=E4outLND,IE4_LN=IE4outLND,H4_LN=H4outLND,LN4_LN=LN4outLND,G4_LN=G4outLND,QE4_LN=QE4outLND,
      P4_LN=P4outLND,L4_LN=L4outLND,
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
      BMDMixture = exp(mabmd)*data.modstan1a$maxD,
      data = data.frame(
        dose = c(data.N$data.modstan1a$x),
        sd = sqrt(data.N$data.modstan1a$s2),
        m = data.N$data.modstan1a$m),
      max.dose = data.N$data.modstan1a$maxD,
      q = data.N$data.modstan1a$q,
      increasing = F,
      models_included = modelnames[prior.weights > 0]
    )

    attr(ret_results, "class") <- c("BMADR", "LP")

    return(ret_results)

  }


}

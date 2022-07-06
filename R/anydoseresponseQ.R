#' Function to determine if there is a dose-response effect, based on Bayes Factor (Quantal data version)
#'
#' @param dose.a dose levels
#' @param y.a mean response
#' @param n.a number of observations per dose level
#' @param cluster logical variable indicating if the data is clustered or not.
#'
#' @examples
#'
#' @description This function tests for any dose-response effect using Bayes factor.
#'              It fits a null model and a saturated model and compare these two models using model posterior probabilities.
#' @return list containing estimated Bayes factor and the decision.
#'
#' @export anydoseresponseQ

anydoseresponseQ <- function(dose.a,y.a,n.a, cluster=FALSE){

  nrch=3;nriter=3000;wu=1000;dl=0.8;trd=10;sd=123;delta=0.999;treedepth=15


  if(cluster==TRUE){

    yamax <- tapply(y.a, dose.a, max, na.rm = TRUE)
    dose2.a <- sort(unique(dose.a))
    namax <- tapply(n.a, dose.a, max, na.rm = TRUE)

    yamin <- tapply(y.a, dose.a, min, na.rm = TRUE)
    namin <- tapply(n.a, dose.a, min, na.rm = TRUE)

    yasum <- tapply(y.a, dose.a, mean, na.rm = TRUE)
    nasum <- tapply(n.a, dose.a, mean, na.rm = TRUE)

    ydiff <- abs(diff(yasum/nasum))

    lbs <- ubs <- numeric(length(y.a)-1)

    for(i in 1:length(y.a)){
      if(i == length(y.a)) break
      for(j in 1:length(yasum)){
        if(j == length(yasum)) break
        if(dose.a[i] == dose2.a[j]) {
          tt1 <- prop.test(c(yamin[j+1], yamin[j]),
                           c(namin[j+1], namin[j]))$conf.int
          tt2 <- prop.test(c(yamax[j+1], yamax[j]),
                           c(namax[j+1], namax[j]))$conf.int

        }

        lbs[i] <- tt1[1]
        ubs[i] <- tt2[2]
      }

    }
    #lbs[lbs <= 0] <- .Machine$double.xmin
    #ubs[ubs <= 0] <- max(ubs)

    N <- length(dose.a)
    datf = data.frame(yy = y.a, n.a = n.a, xx = dose.a)
    fpfit2 <- try(gamlss(cbind(yy,n.a-yy)~as.factor(xx), sigma.formula=~1, family=BB, data=datf),
                  silent = TRUE)
    rhohat <- exp(fpfit2$sigma.coefficients)/(exp(fpfit2$sigma.coefficients)+1)
    dim(rhohat) <- 1

    priorH0 = list(
      priormu = c(max(c(y.a[1]/n.a[1], 1/(5*n.a[1]))), rhohat),
      priorlb = ifelse(y.a[1] != 0, max(c(prop.test(y.a[1], n.a[1])$conf.int[1]/2, 1/(10*n.a[1]))),
                       .Machine$double.xmin),
      priorub = min(c(3*prop.test(y.a[1], n.a[1])$conf.int[2]/2, 1 - 1/(10*n.a[1])))
    )

    #ddy <- diff(y.a/n.a)
    #ddy[ddy <= 0] <- max(c(y.a[1]/n.a[1], 1/(5*n.a[1])))
    ddy <- numeric(length(y.a))
    #ddiffy <- diff(yasum/nasum)

    ddy[dose.a== min(dose.a)] <- max(c(yasum[1]/nasum[1], 1/(5*nasum[1])))
    ddy[dose.a == max(dose.a)] <- ydiff[length(ydiff)]

    for(i in 1:length(dose.a)){
      for(j in 1:length(dose2.a)){

        if(dose.a[i] != min(dose.a) & dose.a[i] != max(dose.a) &
           dose.a[i] == dose2.a[j]){
          ddy[i] <- ydiff[j]
        }
      }
    }

    priorSM = list(
      priormu = c(ddy, rhohat),
      priorlb = c(ifelse(yasum[1] != 0, max(c(prop.test(yasum[1], nasum[1])$conf.int[1]/2, 1/(10*nasum[1]))),
                         .Machine$double.xmin), as.numeric(lbs)
      ),
      priorub = c(min(c(3*prop.test(yasum[1], nasum[1])$conf.int[2]/2, 1 - 1/(10*n.a[1]))), as.numeric(ubs))
    )


    svSM = list(par = ddy, # invsigma2
                rho = rhohat
    )

    svH0 = list(par = c(max(c(y.a[1]/n.a[1], 1/(5*n.a[1])))), rho = rhohat)

    data.modstanSM = list(N=N,n=n.a,y=y.a,
                          priormu=priorSM$priormu,
                          priorlb=priorSM$priorlb, priorub=priorSM$priorub,
                          is_bin=0, is_betabin = 1, priorgama = 4, eps = .Machine$double.xmin
    )

    data.modstanH0 = list(N=N,n=n.a,y=y.a, yint=y.a,nint=n.a,
                          priormu=priorH0$priormu,
                          priorlb=priorH0$priorlb, priorub=priorH0$priorub,
                          is_bin=0, is_betabin = 1, priorgama = 4, eps = .Machine$double.xmin
    )

  } else {

    N <- length(dose.a)
    priorH0 = list(
      priormu = max(c(y.a[1]/n.a[1], 1/(5*n.a[1]))),
      priorlb = ifelse(y.a[1] != 0, max(c(prop.test(y.a[1], n.a[1])$conf.int[1]/2, 1/(10*n.a[1]))),
                       .Machine$double.xmin),
      priorub = min(c(3*prop.test(y.a[1], n.a[1])$conf.int[2]/2, 1 - 1/(10*n.a[1])))
    )

    priorSM = list(
      priormu = c(max(c(y.a[1]/n.a[1], 1/(5*n.a[1]))),
                  ydiff, 0.0),
      priorlb = c(ifelse(y.a[1] != 0, max(c(prop.test(y.a[1], n.a[1])$conf.int[1]/2, 1/(10*n.a[1]))),
                         .Machine$double.xmin), lbs),
      priorub = c(min(c(3*prop.test(y.a[1], n.a[1])$conf.int[2]/2, 1 - 1/(10*n.a[1]))),
                  ubs
      )
    )

    svSM = list(par = c(max(c(y.a[1]/n.a[1], 1/(5*n.a[1]))), # background
                        diff(y.a/n.a) # invsigma2
    ))

    svH0 = list(par = c(max(c(y.a[1]/n.a[1], 1/(5*n.a[1])))))

    data.modstanSM = list(N=N,n=n.a,y=y.a, yint=y.a, nint=n.a,
                          priormu = priorSM$priormu,
                          priorlb=priorSM$priorlb, priorub=priorSM$priorub,
                          is_bin=1, is_betabin = 0, priorgama = 4, eps = .Machine$double.xmin
    )+

      data.modstanH0 = list(N=N,n=n.a,y=y.a, yint=y.a, nint=n.a,
                            priormu=c(priorH0$priormu, 0.0),
                            priorlb=priorH0$priorlb, priorub=priorH0$priorub,
                            is_bin=1, is_betabin = 0, priorgama = 4, eps = .Machine$double.xmin
      )
  }





  # Fitting the null model (no effect)
  # svH0 = list(par=log(mean(y.a)))
  # priorH0 = prior_H0(dose.a, y.a, n.a)
  # data.modstanH0 = list(N=N, n=n.a, y=y.a, priormu=priorH0$priormu, priorSigma=priorH0$priorSigma)
  sv=rstan::optimizing(stanmodels$mH0_Q,data = data.modstanH0,init=svH0)$par
  if(data.modstanH0$is_bin == 1){
    initf2 <- function(chain_id = 1) {
      list(par=sv[1] + rnorm(1, sd = 0.01*abs(sv[2])) ,alpha = chain_id)
    }
  } else if(data.modstanH0$is_betabin == 1) {
    initf2 <- function(chain_id = 1) {
      rho = sv[2]; dim(rho)=1
      list(par=sv[1] + rnorm(1, sd = 0.01*abs(sv[2])), rho = rho ,alpha = chain_id)
    }
  }

  init_ll <- lapply(1:nrch, function(id) initf2(chain_id = id))
  fitstanH0=rstan::sampling(stanmodels$mH0_Q,data = data.modstanH0,init=init_ll,iter = nriter,chains = nrch,warmup=wu,seed=sd,
                     control = list(adapt_delta = delta,max_treedepth =treedepth),
                     show_messages = F, refresh = 0)
  while(is.na(dim(fitstanH0)[1])){

    init_ll <- lapply(1:nrch, function(id) initf2(chain_id = id))

    fitstanH0 = rstan::sampling(stanmodels$mH0_Q, data = data.modstanH0, init=init_ll, iter = nriterations,
                                chains = nrch, warmup = warmup, seed = seed,
                                control = list(adapt_delta = delta, max_treedepth =treedepth),
                                show_messages = F, refresh = 0)
  }


  # Fitting the saturated ANOVA model
  # svSM=list(par=log(y.a))
  # priorSM=prior_SM(dose.a, y.a, n.a)
  # data.modstanSM=list(N=N, n=n.a, y=y.a, priormu=priorSM$priormu, priorSigma=priorSM$priorSigma)
  # sampling(stanmodels$mSM_Q,data = data.modstanSM, chains=1, init = list(svSM))
  sv=rstan::optimizing(stanmodels$mSM_Q,data = data.modstanSM,init=svSM)$par

  if(data.modstanH0$is_bin == 1){
    initf2 <- function(chain_id = 1) {
      list(par=sv[1:N] + rnorm(N, sd = 0.01*abs(sv[1:N])) ,alpha = chain_id)
    }
  } else if(data.modstanH0$is_betabin == 1) {
    initf2 <- function(chain_id = 1) {
      rho = par[N+1]; dim(rho)=1
      list(par=sv[1:N] + rnorm(N, sd = 0.01*abs(sv[1:N])), rho = rho ,alpha = chain_id)
    }
  }


  init_ll <- lapply(1:nrch, function(id) initf2(chain_id = id))
  fitstanSM = rstan::sampling(stanmodels$mSM_Q, data = data.modstanSM, init = init_ll, iter = nriter,chains = nrch, warmup=wu, seed=sd,
                       control = list(adapt_delta = delta,max_treedepth =treedepth),
                       show_messages = F, refresh = 0)
  while(is.na(dim(fitstanSM)[1])){

    init_ll <- lapply(1:nrchains, function(id) initf2(chain_id = id))

    fitstanSM = rstan::sampling(stanmodels$mSM_Q, data = data.modstanSM, init=init_ll, iter = nriterations,
                                chains = nrch, warmup = warmup, seed = seed,
                                control = list(adapt_delta = delta, max_treedepth =treedepth),
                                show_messages = F, refresh = 0)
  }

  set.seed(1234)
  bridge_H0 <- bridgesampling::bridge_sampler(fitstanH0, silent=T)
  bridge_SM <- bridgesampling::bridge_sampler(fitstanSM, silent=T)
  bf=bf(bridge_H0,bridge_SM)
  # pb=post_prob(bridge_sampler(fitstanH0, silent = T),bridge_sampler(fitstanSM, silent = T))
  # print(bf)
  # print(pb)

  if (bf$bf<10){
    mess = "there is sufficient evidence that there is a substantial dose-effect"#; therefore models are fitted and the BMDL is calculated"
  } else{
    mess = "attention: there is insufficient evidence that there is a substantial dose-effect"
  }

  return(list(bf = bf, bf.message = mess))

  # if (bf$bf>10) {
  #   plotSM = function(){
  #     bnds=quantile(as.matrix(fitstanSM)[,6],c(0.025,0.975))
  #     for (di in (2:length(dose.a))){
  #       bnds=cbind(bnds,quantile(as.matrix(fitstanSM)[,length(dose.a)+1+di],c(0.025,0.975)))
  #     }
  #     plot(c(log10(0.00001)-1,log10(1)+1,log10(0.00001)-1,log10(1)+1),log10(c(0,1.25*max(bnds),0,1.25*max(bnds))),type="n",xlab="log10(dose.a+0.00001)",ylab="log10(mean)")
  #     points(log10(dose.a+0.00001),log10(y.a))
  #     points(log10(dose.a+0.00001),log10(bnds[1,]),pch = 2, col="blue")
  #     points(log10(dose.a+0.00001),log10(bnds[2,]),pch = 6, col="blue")
  #     for (dj in (1:length(dose.a))){
  #       lines(c(log10(dose.a[dj]+0.00001),log10(dose.a[dj]+0.00001)),c(log10(bnds[1,dj]),log10(bnds[2,dj])),col="blue")
  #     }
  #   }
  #   print(plotSM())
  #   message("Caution: there is insufficient evidence that there is any dose-effect (see plot)")
  #   # therefore no models are fitted and the BMDL is not calculated
  #   # ")
  # }
  # if (bf$bf<10)   return("there is sufficient evidence that there is a substantial dose-effect; therefore models are fitted and the BMDL is calculated")
}

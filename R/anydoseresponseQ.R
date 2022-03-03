#' Function to determine if there is a dose-response effect, based on Bayes Factor (Quantal data version)
#'
#' @param dose.a dose levels
#' @param y.a mean response
#' @param n.a number of observations per dose level
#' @param cluster logical variable indicating if the data is clustered or not.
#'
#' @return list containing estimated Bayes factor and the decision.
#'
#' @export anydoseresponseQ
#'
anydoseresponseQ <- function(dose.a,y.a,n.a, cluster=FALSE){

  nrch=3;nriter=3000;wu=1000;dl=0.8;trd=10;sd=123;delta=0.999;treedepth=15
  N <- length(unique(dose.a))

  ydiff <- abs(diff(y.a/n.a))

  lbs <- ubs <- numeric(length(y.a)-1)
  for(i in 1:length(y.a)){
    if(i == length(y.a)) break
    tt <- prop.test(c(y.a[i+1], y.a[i]),
                    c(n.a[i+1], n.a[i]))$conf.int
    lbs[i] <- tt[1]
    ubs[i] <- tt[2]
  }

  if(cluster==TRUE){

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

    priorSM = list(
      priormu = c(max(c(y.a[1]/n.a[1], 1/(5*n.a[1]))),
                  ydiff, rhohat),
      priorlb = c(ifelse(y.a[1] != 0, max(c(prop.test(y.a[1], n.a[1])$conf.int[1]/2, 1/(10*n.a[1]))),
                         .Machine$double.xmin), lbs),
      priorub = c(min(c(3*prop.test(y.a[1], n.a[1])$conf.int[2]/2, 1 - 1/(10*n.a[1]))),
                  ubs
      )
    )

    svSM = list(par = c(max(c(y.a[1]/n.a[1], 1/(5*n.a[1]))), # background
                        diff(y.a/n.a)
    ), # invsigma2
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
    )

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
  sv=optimizing(stanmodels$mH0_Q,data = data.modstanH0,init=svH0)$par
  if(data.modstanH0$is_bin == 1){
    initf2 <- function(chain_id = 1) {
      list(par=sv[1] + rnorm(1, sd = 0.01*abs(sv[2])) ,alpha = chain_id)
    }
  } else if(data.modstanH0$is_betabin == 1) {
    initf2 <- function(chain_id = 1) {
      rho = par[2]; dim(rho)=1
      list(par=sv[1] + rnorm(1, sd = 0.01*abs(sv[2])), rho = rho ,alpha = chain_id)
    }
  }

  init_ll <- lapply(1:nrch, function(id) initf2(chain_id = id))
  fitstanH0=sampling(stanmodels$mH0_Q,data = data.modstanH0,init=init_ll,iter = nriter,chains = nrch,warmup=wu,seed=sd,
                     control = list(adapt_delta = delta,max_treedepth =treedepth),
                     show_messages = F, refresh = 0)
  while(is.na(dim(fitstanH0)[1])){

    init_ll <- lapply(1:nrch, function(id) initf2(chain_id = id))

    fitstanH0 = rstan::sampling(stanmodels$stanmodels$mH0_Q, data = data.modstanH0, init=init_ll, iter = nriterations,
                                chains = nrch, warmup = warmup, seed = seed,
                                control = list(adapt_delta = delta, max_treedepth =treedepth),
                                show_messages = F, refresh = 0)
  }


  # Fitting the saturated ANOVA model
  # svSM=list(par=log(y.a))
  # priorSM=prior_SM(dose.a, y.a, n.a)
  # data.modstanSM=list(N=N, n=n.a, y=y.a, priormu=priorSM$priormu, priorSigma=priorSM$priorSigma)
  sv=optimizing(stanmodels$mSM_Q,data = data.modstanSM,init=svSM)$par

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
  fitstanSM = sampling(stanmodels$mSM_Q, data = data.modstanSM, init = init_ll, iter = nriter,chains = nrch, warmup=wu, seed=sd,
                       control = list(adapt_delta = delta,max_treedepth =treedepth),
                       show_messages = F, refresh = 0)
  while(is.na(dim(fitstanSM)[1])){

    init_ll <- lapply(1:nrchains, function(id) initf2(chain_id = id))

    fitstanSM = rstan::sampling(stanmodels$stanmodels$mSM_Q, data = data.modstanSM, init=init_ll, iter = nriterations,
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

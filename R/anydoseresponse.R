#' Function to determine if there is a dose-response effect, for the normal distribution
#'
#' @param dose.a ordered dose levels
#' @param mean.a mean response
#' @param sd.a standard deviation
#' @param n.a number of observations per dose level
#'
#' @description This function tests for any dose-response effect using Bayes factor.
#'              It fits a null model and a saturated model and compares these two models using model posterior probabilities obtained via bridge rstan::sampling.
#'
#' @examples
#'
#' # We use the first 5 rows because those are observations from subjects belonging to the same group.
#' data("immunotoxicityData.rda")  #load the immunotoxicity data
#' anydoseresponseN(dose.a = immunotoxicityData$Dose[1:5],
#'                  mean.a = immunotoxicityData$Mean[1:5],
#'                  sd.a = immunotoxicityData$SD[1:5],
#'                  n.a = immunotoxicityData$n[1:5])
#'
#' @return list containing Bayes factor and decision.
#'
#' @export anydoseresponseN
anydoseresponseN=function(dose.a,mean.a,sd.a,n.a){

  if(length(dose.a) != length(unique(dose.a))){
    dose = sort(unique(dose.a))
    N = length(dose)
    mean=rep(NA,N)
    sd=rep(NA,N)
    n=rep(NA,N)
    for (iu in (1:N)){
      mean[iu] = mean(mean.a[dose.a == dose[iu]])
      sd[iu] = mean(sd.a[dose.a == dose[iu]])
      n[iu] = sum(n.a[dose.a == dose[iu]])
    }
    mean.a = mean
    dose.a = dose
    sd.a = sd
    n.a = n
  }

  N = length(dose.a)

  if(mean.a[1] < mean.a[N]){
    data_type = 1
  }else{
    data_type = 3
  }

  nrch=3;nriter=3000;wu=1000;dl=0.8;trd=10;sd=123;delta=0.999;treedepth=15

  priorH0 = list(
    priormu = c(mean.a[1], -2*log(1.5*mean(sd.a))),
    priorSigma = diag(c(1,1)),
    priorlb = 0.001,
    priorub = 2*mean.a[1]
  )

  priorSM = list(
    priormu = c(mean.a[1],
                diff(mean.a),
                -2*log(1.5*mean(sd.a))),
    priorSigma = diag(c(1, rep(1, length(dose.a)-1), 1)),
    priorlb = 0.001,
    priorub = c(2*mean.a[1],
                max(abs(diff(mean.a)))*10
    )
  )

  svSM = list(par = c(mean.a[1], # background
                      diff(mean.a),
                      log(1/mean(sd.a^2))) # invsigma2
  )

  svH0 = list(par = c(mean.a[1], log(1/mean(sd.a^2))))

  data.modstanSM = list(N=N,n=n.a,m=mean.a,s2=sd.a^2,shift=0,
                        priormu=priorSM$priormu, priorSigma=priorSM$priorSigma,
                        priorlb=priorSM$priorlb, priorub=priorSM$priorub,
                        data_type=data_type, priorg = 4
  )

  data.modstanH0 = list(N=N,n=n.a,m=mean.a,s2=sd.a^2,shift=0,
                        priormu=priorH0$priormu, priorSigma=priorH0$priorSigma,
                        priorlb=priorH0$priorlb, priorub=priorH0$priorub,
                        data_type=data_type, priorg = 4
  )

  # Fitting the null model (no effect)
  sv=rstan::optimizing(stanmodels$mH0,data = data.modstanH0,init=svH0)$par
  initf2 <- function(chain_id = 1) {
    list(par=sv[1:2] + rnorm(2, sd = 0.01*abs(sv[2])) ,alpha = chain_id)
  }
  init_ll <- lapply(1:nrch, function(id) initf2(chain_id = id))
  fitstanH0=rstan::sampling(stanmodels$mH0,data = data.modstanH0,init=init_ll,iter = nriter,chains = nrch,warmup=wu,seed=sd,
                            control = list(adapt_delta = delta,max_treedepth =treedepth),
                            show_messages = F, refresh = 0)

  # Fitting the saturated ANOVA model
  sv=rstan::optimizing(stanmodels$mSM,data = data.modstanSM,init=svSM)$par
  initf2 <- function(chain_id = 1) {
    list(par=sv[1:(N+1)] + rnorm(N+1, sd = 0.01*abs(sv[1:(N+1)])) ,alpha = chain_id)
  }
  init_ll <- lapply(1:nrch, function(id) initf2(chain_id = id))
  fitstanSM=rstan::sampling(stanmodels$mSM,data = data.modstanSM,init=init_ll,iter = nriter,chains = nrch,warmup=wu,seed=sd,
                            control = list(adapt_delta = delta,max_treedepth =treedepth),
                            show_messages = F, refresh = 0)

  set.seed(1234)
  bridge_H0 <- bridgesampling::bridge_sampler(fitstanH0, silent=T)
  bridge_SM <- bridgesampling::bridge_sampler(fitstanSM, silent=T)
  bf=bf(bridge_H0,bridge_SM)
  # pb=post_prob(bridge_sampler(fitstanH0, silent=T),bridge_sampler(fitstanSM, silent=T))
  # print(bf)
  # print(pb)

  if (bf$bf<10){
    mess = "there is sufficient evidence that there is a substantial dose-effect"#; therefore models are fitted and the BMDL is calculated"
  } else{
    mess = "attention: there is insufficient evidence that there is a substantial dose-effect"
  }

  return(list(bf = bf, bf.message = mess))

}

#' Function to determine if there is a dose-response effect, for the log-normal distribution
#'
#' @param dose.a ordered dose levels
#' @param mean.a mean response
#' @param sd.a standard deviation
#' @param n.a number of observations per dose level
#'
#' @description This function tests for any dose-response effect using Bayes factor.
#'              It fits a null model and a saturated model and compares these two models using model posterior probabilities obtained via bridge rstan::sampling.
#'
#' @examples
#'
#' # We use the first 5 rows because those are observations from subjects belonging to the same group.
#' data("immunotoxicityData.rda")  #load the immunotoxicity data
#' anydoseresponseLN(dose.a = immunotoxicityData$Dose[1:5],
#'                  mean.a = immunotoxicityData$Mean[1:5],
#'                  sd.a = immunotoxicityData$SD[1:5],
#'                  n.a = immunotoxicityData$n[1:5])
#'
#' @return list containing Bayes factor and decision.
#'
#' @export anydoseresponseLN
anydoseresponseLN=function(dose.a,mean.a,sd.a,n.a){

  if(length(dose.a) != length(unique(dose.a))){
    dose = sort(unique(dose.a))
    N = length(dose)
    mean=rep(NA,N)
    sd=rep(NA,N)
    n=rep(NA,N)
    for (iu in (1:N)){
      mean[iu] = mean(mean.a[dose.a == dose[iu]])
      sd[iu] = mean(sd.a[dose.a == dose[iu]])
      n[iu] = sum(n.a[dose.a == dose[iu]])
    }
    mean.a = mean
    dose.a = dose
    sd.a = sd
    n.a = n
  }

  N = length(dose.a)
  shift=0
  gmean.a2 = log(NtoLN(mean.a,sd.a))[1:N]
  if (min(gmean.a2)<0) {gmean.a = gmean.a2-20*min(gmean.a2); shift=20*min(gmean.a2)}
  if (min(gmean.a2)>=0) gmean.a = gmean.a2
  gsd.a = log(NtoLN(mean.a,sd.a))[(N+1):(2*N)]

  if(gmean.a2[1] < gmean.a2[N]){
    data_type = 2
  }else{
    data_type = 4
  }

  nrch=3;nriter=3000;wu=1000;dl=0.8;trd=10;sd=123;delta=0.999;treedepth=15

  priorH0 = list(
    priormu = c(exp(gmean.a2[1]), -2*log(1.5*mean(sd.a))),
    priorSigma = diag(c(1,1)),
    priorlb = 0.001,
    priorub = 2*exp(gmean.a2[1])
  )

  priorSM = list(
    priormu = c(exp(gmean.a2[1]),
                diff(exp(gmean.a2)),
                -2*log(1.5*mean(sd.a))),
    priorSigma = diag(c(1, rep(1, length(dose.a)-1), 1)),
    priorlb = 0.001,
    priorub = c(2*exp(gmean.a2[1]),
                max(abs(diff(exp(gmean.a2))))*10
    )
  )

  svSM = list(par = c(exp(gmean.a2[1]), # background
                      diff(exp(gmean.a2)),
                      log(1/mean(sd.a^2))) # invsigma2
  )

  svH0 = list(par = c(exp(gmean.a2[1]), log(1/mean(sd.a^2))))

  data.modstanSM = list(N=N,n=n.a,m=gmean.a,s2=gsd.a^2,shift=shift,
                        priormu=priorSM$priormu, priorSigma=priorSM$priorSigma,
                        priorlb=priorSM$priorlb, priorub=priorSM$priorub,
                        data_type=data_type, priorg = 4
  )

  data.modstanH0 = list(N=N,n=n.a,m=gmean.a,s2=gsd.a^2,shift=shift,
                        priormu=priorH0$priormu, priorSigma=priorH0$priorSigma,
                        priorlb=priorH0$priorlb, priorub=priorH0$priorub,
                        data_type=data_type, priorg = 4
  )

  # Fitting the null model (no effect)
  sv=rstan::optimizing(stanmodels$mH0,data = data.modstanH0,init=svH0)$par
  initf2 <- function(chain_id = 1) {
    list(par=sv[1:2] + rnorm(2, sd = 0.01*abs(sv[2])) ,alpha = chain_id)
  }
  init_ll <- lapply(1:nrch, function(id) initf2(chain_id = id))
  fitstanH0=rstan::sampling(stanmodels$mH0,data = data.modstanH0,init=init_ll,iter = nriter,chains = nrch,warmup=wu,seed=sd,
                            control = list(adapt_delta = delta,max_treedepth =treedepth),
                            show_messages = F, refresh = 0)

  # Fitting the saturated ANOVA model
  sv=rstan::optimizing(stanmodels$mSM,data = data.modstanSM,init=svSM)$par
  initf2 <- function(chain_id = 1) {
    list(par=sv[1:(N+1)] + rnorm(N+1, sd = 0.01*abs(sv[1:(N+1)])) ,alpha = chain_id)
  }
  init_ll <- lapply(1:nrch, function(id) initf2(chain_id = id))
  fitstanSM=rstan::sampling(stanmodels$mSM,data = data.modstanSM,init=init_ll,iter = nriter,chains = nrch,warmup=wu,seed=sd,
                            control = list(adapt_delta = delta,max_treedepth =treedepth),
                            show_messages = F, refresh = 0)

  set.seed(1234)
  bridge_H0 <- bridgesampling::bridge_sampler(fitstanH0, silent=T)
  bridge_SM <- bridgesampling::bridge_sampler(fitstanSM, silent=T)
  bf=bf(bridge_H0,bridge_SM)
  # pb=post_prob(bridge_sampler(fitstanH0, silent=T),bridge_sampler(fitstanSM, silent=T))
  # print(bf)
  # print(pb)

  if (bf$bf<10){
    mess = "there is sufficient evidence that there is a substantial dose-effect"#; therefore models are fitted and the BMDL is calculated"
  } else{
    mess = "attention: there is insufficient evidence that there is a substantial dose-effect"
  }

  return(list(bf = bf, bf.message = mess))

}

#' Function to determine if there is a dose-response effect, for the normal distribution
#'
#' @param data Individual data with columns: dose, response, litter
#' @param use.mcmc logical indicating wheter MCMC (TRUE) or Laplace (FALSE, default) should be used
#'
#' @description This function tests for any dose-response effect using Bayes factor.
#'              It fits a null model and a saturated model and compares these two models using model posterior probabilities obtained via bridge sampling.
#'              Currently only based on the normal distribution.
#'
#' @return list containing Bayes factor and decision.
#'
#' @export anydoseresponseC
anydoseresponseC=function(data, use.mcmc = FALSE){

  nrch=3;nriter=300;wu=100;dl=0.8;trd=10;sd=123;delta=0.999;treedepth=15;ndr = 30000

  indiv.data <- data.frame(dose = data[,1],
                           response = data[,2],
                           litter = data[,3])
  indiv.data <- indiv.data %>%
    dplyr::group_by(dose, litter) %>%
    dplyr::arrange(by_group = dose) # order dose groups
  indiv.data <- indiv.data %>%
    dplyr::mutate(cluster = dplyr::cur_group_id(),
                  count = n())
  dose.a = indiv.data$dose
  maxDose = max(dose.a)
  doses = unique(dose.a)
  N = length(unique(dose.a)) # dose groups
  n = c() # number of litters per dose group (vector of size N)
  for(i in 1:N){
    cnt = plyr::count(indiv.data$litter[indiv.data$dose==doses[i]])
    n[i] = length(unique(cnt$x))
  }
  nc = length(unique(indiv.data$cluster)) # number of unique dose x litter combinations (i.e. clusters)
  cid = unique(indiv.data$cluster) # cluster ids
  maxN = max(indiv.data$count) # max number of obs per cluster
  maxNc = max(n) # max number of litters per dose group
  # nij = as.matrix(table(indiv.data$dose, indiv.data$litter))
  nij = matrix(0, nrow = N, ncol = maxNc)
  for(i in 1:N){
    cnt = plyr::count(indiv.data$litter[indiv.data$dose==doses[i]])
    obs = cnt$freq
    if(length(obs) < maxNc){
      obs = c(obs, rep(0, maxNc-length(obs)))
    }
    nij[i, ] = obs
  }
  y = matrix(0, nrow = nc, ncol = maxN)
  for(i in 1:nc){
    obs = indiv.data$response[indiv.data$cluster==cid[i]]
    if(length(obs) < maxN){
      obs = c(obs, rep(0, maxN-length(obs)))
    }
    y[i, ] = obs
  }

  means.all <- indiv.data %>%
    group_by(dose) %>%
    summarise(mresp = mean(response))
  dose.a = unique(indiv.data$dose)
  mean.a = c()
  for(m in 1:length(dose.a)){
    mean.a[m] <- means.all$mresp[means.all$dose == dose.a[m]]
  }

  if(mean.a[1] < mean.a[N]){
    data_type = 1
  }else{
    data_type = 3
  }

  priorH0 = list(
    priormu = c(mean.a[1], -2*log(1.5*sd(y[y!=0]))),
    priorSigma = diag(c(1,1)),
    priorlb = 0.001,0,
    priorub = 2*mean.a[1]
  )

  priorSM = list(
    priormu = c(mean.a[1],
                diff(mean.a),
                -2*log(1.5*sd(y[y!=0])),
                0.5),
    priorSigma = diag(c(1, rep(1, length(dose.a)-1), 1)),
    priorlb = 0.001,
    priorub = c(2*mean.a[1],
                max(abs(diff(mean.a)))*10)
  )

  svSM = list(par = c(mean.a[1],
                      diff(mean.a),
                      log(1/var(y[y!=0])),
                      0.5))

  svH0 = list(par = c(mean.a[1], log(1/var(y[y!=0])), 0.5))



  data.modstanSM = list(N=N, n=n, nc=nc, maxN=maxN, maxNc=maxNc,
                        nij=nij, y=y, q=0, shift=0,
                        priormu=priorSM$priormu, priorSigma=priorSM$priorSigma,
                        priorlb=priorSM$priorlb, priorub=priorSM$priorub,
                        priorg=4, data_type=data_type
  )

  data.modstanH0 = list(N=N, n=n, nc=nc, maxN=maxN, maxNc=maxNc,
                        nij=nij, y=y, q=0, shift=0,
                        priormu=priorH0$priormu, priorSigma=priorH0$priorSigma,
                        priorlb=priorH0$priorlb, priorub=priorH0$priorub,
                        priorg=4, data_type=data_type

  )

  if(use.mcmc == TRUE){

    # Fitting the null model (no effect)
    sv=optimizing(mH0c, data = data.modstanH0, init=svH0)$par
    initf2 <- function(chain_id = 1) {
      list(par=sv[1:3] + rnorm(3, sd = 0.01*abs(sv[1:3])) ,alpha = chain_id)
    }
    init_ll <- lapply(1:nrch, function(id) initf2(chain_id = id))
    fitstanH0=sampling(mH0c,data = data.modstanH0,init=init_ll,iter = nriter,chains = nrch,warmup=wu,seed=sd,
                       control = list(adapt_delta = delta,max_treedepth =treedepth),
                       show_messages = F, refresh = 0)

    # Fitting the saturated ANOVA model
    sv=optimizing(mSMc,data = data.modstanSM,init=svSM)$par
    initf2 <- function(chain_id = 1) {
      list(par=sv[1:(N+2)] + rnorm(N+2, sd = 0.01*abs(sv[1:(N+2)])) ,alpha = chain_id)
    }
    init_ll <- lapply(1:nrch, function(id) initf2(chain_id = id))
    fitstanSM=sampling(mSMc,data = data.modstanSM,init=init_ll,iter = nriter,chains = nrch,warmup=wu,seed=sd,
                       control = list(adapt_delta = delta,max_treedepth =treedepth),
                       show_messages = F, refresh = 0)

    set.seed(1234)
    bridge_H0 <- bridgesampling::bridge_sampler(fitstanH0, silent=T)
    bridge_SM <- bridgesampling::bridge_sampler(fitstanSM, silent=T)
    bf=bf(bridge_H0,bridge_SM)
    # pb=post_prob(bridge_sampler(fitstanH0, silent=T),bridge_sampler(fitstanSM, silent=T))
    # print(bf)
    # print(pb)
    bf = bf$bf


  }else if(use.mcmc == FALSE){

    fitH0 = rstan::optimizing(mH0c, data = data.modstanH0, seed = as.integer(sd), draws = ndr, init = svH0, hessian = T)
    pars.H0 = apply(fitH0$theta_tilde[, c(paste0('a'), paste0('par[2]'),
                                          paste0('rho_cluster'))], 2, median)

    fitSM = rstan::optimizing(mSMc, data = data.modstanSM, seed = as.integer(sd), draws = ndr, init = svSM, hessian = T)
    pars.SM = apply(fitSM$theta_tilde[, c(paste0('a[', 1:N, ']'), paste0('par[', N+1, ']'),
                                          paste0('par[', N+2, ']'))], 2, median)



    llSM = llfSM_Nc(pars.SM,
                    d=doses,
                    n=n,
                    nij=nij,
                    y=y,
                    qval=0)

    llH0 = llfH0_Nc(pars.H0,
                    d=doses,
                    n=n,
                    nij=nij,
                    y=y,
                    qval=0)


    BIC.H0 = - 2 * llH0 + (2 * log(sum(y!=0))) # 2 parameters: a and rho

    BIC.SM = - 2 * llSM + (N+2 * log(sum(y!=0))) # N mean parameters + variance + rho

    bf = exp(-0.5 * (BIC.H0 - BIC.SM))

  }

  if (bf<10){
    mess = "there is sufficient evidence that there is a substantial dose-effect"#; therefore models are fitted and the BMDL is calculated"
  } else{
    mess = "attention: there is insufficient evidence that there is a substantial dose-effect"
  }

  return(list(bf = bf, bf.message = mess))
}

#' Function to determine if there is a dose-response effect, for quantal endpoints
#'
#' @param dose.a ordered dose levels
#' @param y.a number of adverse events
#' @param n.a number of observations per dose level
#' @param cluster logical variable indicating whether the data are clustered or not; if data are clustered, it is assumed that each line of data represents a specific litter
#' @param use.mcmc logical indicating wheter MCMC (TRUE) or Laplace (FALSE, default) should be used for the case of clustered data
#'
#' @description This function tests for any dose-response effect using Bayes factor.
#'              It fits a null model and a saturated model and compares these two models using model posterior probabilities obtained via bridge rstan::sampling.
#'
#' @return list containing estimated Bayes factor and the decision.
#'
#' @export anydoseresponseQ
#'
anydoseresponseQ <- function(dose.a,y.a,n.a, cluster=FALSE, use.mcmc = FALSE){

  ndr=30000;nrch=3;nriter=3000;wu=1000;dl=0.8;trd=10;sd=123;delta=0.999;treedepth=15

  # if(length(dose.a) != length(unique(dose.a))){
  #   dose = sort(unique(dose.a))
  #   N = length(dose)
  #   y=rep(NA,N)
  #   n=rep(NA,N)
  #   for (iu in (1:N)){
  #     y[iu] = sum(y.a[dose.a == dose[iu]])
  #     n[iu] = sum(n.a[dose.a == dose[iu]])
  #   }
  #   y.a = y
  #   dose.a = dose
  #   n.a = n
  # }
  maxDose = max(dose.a)
  N = length(dose.a)
  dose.a = dose.a/maxDose


  if(cluster==TRUE){

    yasum <- tapply(y.a, dose.a, sum, na.rm = TRUE)
    nasum <- tapply(n.a, dose.a, sum, na.rm = TRUE)

    yamean <- yasum/nasum
    ydiff <- diff(yasum/nasum)

    #tt1 <- prop.test(yasum[1], nasum[1])$conf.int
    lbs <- ifelse(yamean[1] != 0, max(c(prop.test(yasum[1], nasum[1])$conf.int[1]/2, 1/(10*nasum[1]))),
                  .Machine$double.xmin)
    ubs <- min(c(3*prop.test(yasum[1], nasum[1])$conf.int[2]/2, 1 - 1/(10*nasum[1])))


    N <- length(dose.a)
    datf = data.frame(yy = y.a, n.a = n.a, xx = dose.a)
    fpfit2 <- try(gamlss(cbind(yy,n.a-yy)~as.factor(xx), sigma.formula=~1, family=BB, data=datf),
                  silent = TRUE)
    rhohat <- exp(fpfit2$sigma.coefficients)/(exp(fpfit2$sigma.coefficients)+1)
    dim(rhohat) <- 1

    priorH0 = list(
      priormu = c(max(c(yasum[1]/nasum[1], 1/(5*nasum[1]))), rhohat),
      priorlb = lbs,
      priorub = ubs
    )

    priorSM = list(
      priormu = c(max(c(yasum[1]/nasum[1], 1/(5*nasum[1]))), rhohat),
      priorlb = lbs,
      priorub = ubs
    )

    ddy <- c(max(c(yasum[1]/nasum[1], 1/(5*nasum[1]))),diff(yamean))
    svSM = list(par = ddy, # invsigma2
                rho = rhohat
    )

    svH0 = list(par = max(c(yasum[1]/nasum[1], 1/(5*nasum[1]))), rho = rhohat)

    data.modstanSM = list(N=N,Ndose=length(unique(dose.a)),n=n.a,y=y.a,
                          priormu=priorSM$priormu,
                          priorlb=priorSM$priorlb, priorub=priorSM$priorub,
                          is_bin=0, is_betabin = 1, priorgama = 4, eps = .Machine$double.xmin
    )

    data.modstanH0 = list(N=N,n=n.a,y=y.a, yint=y.a,nint=n.a,
                          priormu=priorH0$priormu,
                          priorlb=priorH0$priorlb, priorub=priorH0$priorub,
                          is_bin=0, is_betabin = 1, priorgama = 4, eps = .Machine$double.xmin
    )

    if(use.mcmc == TRUE){
      svNM=rstan::optimizing(stanmodels$mH0_Q,data = data.modstanH0,init=svH0)$par

      initf2 <- function(chain_id = 1) {
        rho = svNM[2]; dim(rho)=1
        list(par=svNM[1] + rnorm(1, sd = 0.01*abs(svNM[2])), rho = rho ,alpha = chain_id)
      }

      init_ll <- lapply(1:nrch, function(id) initf2(chain_id = id))
      fitstanH0=rstan::sampling(stanmodels$mH0_Q,data = data.modstanH0,init=init_ll,iter = nriter,chains = nrch,warmup=wu,seed=sd,
                                control = list(adapt_delta = delta,max_treedepth =treedepth),
                                show_messages = F, refresh = 0)
      while(is.na(dim(fitstanH0)[1])){

        init_ll <- lapply(1:nrch, function(id) initf2(chain_id = id))

        fitstanH0 = rstan::sampling(stanmodels$mH0_Q, data = data.modstanH0, init=init_ll, iter = nriterations,
                                    chains = nrch, warmup = warmup, seed = sd,
                                    control = list(adapt_delta = delta, max_treedepth =treedepth),
                                    show_messages = F, refresh = 0)
      }




      # Fitting the saturated ANOVA model
      # svSM=list(par=log(y.a))
      # priorSM=prior_SM(dose.a, y.a, n.a)
      # data.modstanSM=list(N=N, n=n.a, y=y.a, priormu=priorSM$priormu, priorSigma=priorSM$priorSigma)
      #sampling(stanmodels$mSM2_Q,data = data.modstanSM, chains=1, init = list(svSM))
      svH1=rstan::optimizing(stanmodels$mSM_Q,data = data.modstanSM,init=svSM)$par

      initf2 <- function(chain_id = 1) {
        nns <- which(stringr::str_detect(names(svH1),'par'))
        nns_rho <- which(stringr::str_detect(names(svH1),'rho'))

        rho = svH1[nns_rho]; dim(rho)=1
        list(par=svH1[nns] +
               rnorm(length(nns), sd = 0.01*abs(svH1[nns])),
             rho = rho + rnorm(length(nns_rho), sd = 0.01*abs(svH1[nns_rho])), alpha = chain_id)
      }



      init_ll <- lapply(1:nrch, function(id) initf2(chain_id = id))
      fitstanSM = rstan::sampling(stanmodels$mSM_Q, data = data.modstanSM, init = init_ll, iter = nriter,chains = nrch, warmup=wu, seed=sd,
                                  control = list(adapt_delta = delta,max_treedepth =treedepth),
                                  show_messages = F, refresh = 0)
      while(is.na(dim(fitstanSM)[1])){

        init_ll <- lapply(1:nrchains, function(id) initf2(chain_id = id))

        fitstanSM = rstan::sampling(stanmodels$mSM_Q, data = data.modstanSM, init=init_ll, iter = nriterations,
                                    chains = nrch, warmup = warmup, seed = sd,
                                    control = list(adapt_delta = delta, max_treedepth =treedepth),
                                    show_messages = F, refresh = 0)
      }

      set.seed(1234)
      bridge_H0 <- bridgesampling::bridge_sampler(fitstanH0, silent=T)
      bridge_SM <- bridgesampling::bridge_sampler(fitstanSM, silent=T)
      bf=bf(bridge_H0,bridge_SM)
      bf = bf$bf

    }else if(use.mcmc == FALSE){

      fitH0 = rstan::optimizing(stanmodels$mH0_Q, data = data.modstanH0, seed = as.integer(sd), draws = ndr, init = svH0, hessian = T)


      pH0 <- fitH0$par
      pars.H0 = pH0[names(pH0) %in% c(paste0('a[', 1:length(unique(dose.a)), ']'),
                                      paste0('par[', length(unique(dose.a)), ']'),
                                      'rho[1]')]



      fitSM = rstan::optimizing(stanmodels$mSM_Q, data = data.modstanSM, seed = as.integer(sd), draws = ndr, init = svSM, hessian = T)


      pSM <- fitSM$par
      pars.SM = pSM[names(pSM) %in% c(paste0('a[', 1:length(unique(dose.a)), ']'),
                                      paste0('par[', length(unique(dose.a)), ']'),
                                      'rho[1]')]



      llSM = llfSM2_Q(pars.SM[stringr::str_detect(names(pars.SM), 'a\\[')],
                      n.a, dose.a, y.a, 0,
                      pars.SM[stringr::str_detect(names(pars.SM), 'rho')]) ## UPDATE --> rho? !!
      llH0 = llfH02_Q(pars.H0[stringr::str_detect(names(pars.H0), 'a\\[')],
                      n.a, dose.a, y.a, 0,
                      pars.H0[stringr::str_detect(names(pars.H0), 'rho')])


      BIC.H0 = - 2 * llH0 + (2 * log(sum(n.a))) # 2 parameters: a and rho

      BIC.SM = - 2 * llSM + (N+1 * log(sum(n.a)))


      bf = exp(-0.5 * (BIC.H0 - BIC.SM))

    }

  } else {

    N <- length(dose.a)
    priorH0 = list(
      priormu = max(c(y.a[1]/n.a[1], 1/(5*n.a[1]))),
      priorlb = ifelse(y.a[1] != 0, max(c(prop.test(y.a[1], n.a[1])$conf.int[1]/2, 1/(10*n.a[1]))),
                       .Machine$double.xmin),
      priorub = min(c(3*prop.test(y.a[1], n.a[1])$conf.int[2]/2, 1 - 1/(10*n.a[1])))
    )

    priorSM = list(
      priormu = c(max(c(y.a[1]/n.a[1], 1/(5*n.a[1]))), 0.0),
      priorlb = ifelse(y.a[1] != 0, max(c(prop.test(y.a[1], n.a[1])$conf.int[1]/2, 1/(10*n.a[1]))),
                       .Machine$double.xmin),
      priorub = min(c(3*prop.test(y.a[1], n.a[1])$conf.int[2]/2, 1 - 1/(10*n.a[1])))
    )

    svSM = list(par = c(max(c(y.a[1]/n.a[1], 1/(5*n.a[1]))),
                        diff(y.a/n.a)
    ))

    svH0 = list(par = c(max(c(y.a[1]/n.a[1], 1/(5*n.a[1])))))

    data.modstanSM = list(N=N,Ndose=length(unique(dose.a)),n=n.a,y=y.a, yint=y.a, nint=n.a,
                          priormu = priorSM$priormu,
                          priorlb=priorSM$priorlb, priorub=priorSM$priorub,
                          is_bin=1, is_betabin = 0, priorgama = 4, eps = .Machine$double.xmin
    )

    data.modstanH0 = list(N=N,n=n.a,y=y.a, yint=y.a, nint=n.a,
                          priormu=c(priorH0$priormu, 0.0),
                          priorlb=priorH0$priorlb, priorub=priorH0$priorub,
                          is_bin=1, is_betabin = 0, priorgama = 4, eps = .Machine$double.xmin
    )

    svNM=rstan::optimizing(stanmodels$mH0_Q,data = data.modstanH0,init=svH0)$par

    initf2 <- function(chain_id = 1) {
      list(par=svNM[1] + rnorm(1, sd = 0.01*abs(svNM[2])) ,alpha = chain_id)
    }

    init_ll <- lapply(1:nrch, function(id) initf2(chain_id = id))
    fitstanH0=rstan::sampling(stanmodels$mH0_Q,data = data.modstanH0,init=init_ll,iter = nriter,chains = nrch,warmup=wu,seed=sd,
                              control = list(adapt_delta = delta,max_treedepth =treedepth),
                              show_messages = F, refresh = 0)
    while(is.na(dim(fitstanH0)[1])){

      init_ll <- lapply(1:nrch, function(id) initf2(chain_id = id))

      fitstanH0 = rstan::sampling(stanmodels$mH0_Q, data = data.modstanH0, init=init_ll, iter = nriterations,
                                  chains = nrch, warmup = warmup, seed = sd,
                                  control = list(adapt_delta = delta, max_treedepth =treedepth),
                                  show_messages = F, refresh = 0)
    }

    # Fitting the saturated ANOVA model
    # svSM=list(par=log(y.a))
    # priorSM=prior_SM(dose.a, y.a, n.a)
    # data.modstanSM=list(N=N, n=n.a, y=y.a, priormu=priorSM$priormu, priorSigma=priorSM$priorSigma)
    #sampling(stanmodels$mSM2_Q,data = data.modstanSM, chains=1, init = list(svSM))
    svH1=rstan::optimizing(stanmodels$mSM_Q,data = data.modstanSM,init=svSM)$par


    initf2 <- function(chain_id = 1) {
      nns <- which(stringr::str_detect(names(svH1),'par'))
      list(par=svH1[nns] +
             rnorm(length(nns), sd = 0.01*abs(svH1[nns])), alpha = chain_id)
    }



    init_ll <- lapply(1:nrch, function(id) initf2(chain_id = id))
    fitstanSM = rstan::sampling(stanmodels$mSM_Q, data = data.modstanSM, init = init_ll, iter = nriter,chains = nrch, warmup=wu, seed=sd,
                                control = list(adapt_delta = delta,max_treedepth =treedepth),
                                show_messages = F, refresh = 0)
    while(is.na(dim(fitstanSM)[1])){

      init_ll <- lapply(1:nrchains, function(id) initf2(chain_id = id))

      fitstanSM = rstan::sampling(stanmodels$mSM_Q, data = data.modstanSM, init=init_ll, iter = nriterations,
                                  chains = nrch, warmup = warmup, seed = sd,
                                  control = list(adapt_delta = delta, max_treedepth =treedepth),
                                  show_messages = F, refresh = 0)
    }

    set.seed(1234)
    bridge_H0 <- bridgesampling::bridge_sampler(fitstanH0, silent=T)
    bridge_SM <- bridgesampling::bridge_sampler(fitstanSM, silent=T)
    bf=bf(bridge_H0,bridge_SM)
    bf = bf$bf

  }


  # Fitting the null model (no effect)
  # svH0 = list(par=log(mean(y.a)))
  # priorH0 = prior_H0(dose.a, y.a, n.a)
  # data.modstanH0 = list(N=N, n=n.a, y=y.a, priormu=priorH0$priormu, priorSigma=priorH0$priorSigma)


  # pb=post_prob(bridge_sampler(fitstanH0, silent = T),bridge_sampler(fitstanSM, silent = T))
  # print(bf)
  # print(pb)

  if (bf<10){
    mess = "there is sufficient evidence that there is a substantial dose-effect"#; therefore models are fitted and the BMDL is calculated"
  } else{
    mess = "attention: there is insufficient evidence that there is a substantial dose-effect"
  }


  return(list(bf = bf, bf.message = mess))
}

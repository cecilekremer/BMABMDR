#' Function to determine if there is a dose-response effect (Normally distributed data version)
#'
#' @param dose.a dose levels
#' @param mean.a mean response
#' @param sd.a standard deviation
#' @param n.a number of observations per dose level
#'
#' @description This function tests for any dose-response effect using Bayes factor.
#'              It fits a null model and a saturated model and compare these two models using model posterior probabilities.
#'
#' @examples
#'
#' # we use the first 5 rows because those are observations from subjects belonging to the same group.
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

  nrch=3;nriter=3000;wu=1000;dl=0.8;trd=10;sd=123;delta=0.999;treedepth=15
  N = length(dose.a)

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
                        data_type=1, priorg = 4
  )

  data.modstanH0 = list(N=N,n=n.a,m=mean.a,s2=sd.a^2,shift=0,
                        priormu=priorH0$priormu, priorSigma=priorH0$priorSigma,
                        priorlb=priorH0$priorlb, priorub=priorH0$priorub,
                        data_type=1, priorg = 4
  )

  # Fitting the null model (no effect)
  sv=optimizing(stanmodels$mH0,data = data.modstanH0,init=svH0)$par
  initf2 <- function(chain_id = 1) {
    list(par=sv[1:2] + rnorm(2, sd = 0.01*abs(sv[2])) ,alpha = chain_id)
  }
  init_ll <- lapply(1:nrch, function(id) initf2(chain_id = id))
  fitstanH0=sampling(stanmodels$mH0,data = data.modstanH0,init=init_ll,iter = nriter,chains = nrch,warmup=wu,seed=sd,
                     control = list(adapt_delta = delta,max_treedepth =treedepth),
                     show_messages = F, refresh = 0)

  # Fitting the saturated ANOVA model
  sv=optimizing(stanmodels$mSM,data = data.modstanSM,init=svSM)$par
  initf2 <- function(chain_id = 1) {
    list(par=sv[1:(N+1)] + rnorm(N+1, sd = 0.01*abs(sv[1:(N+1)])) ,alpha = chain_id)
  }
  init_ll <- lapply(1:nrch, function(id) initf2(chain_id = id))
  fitstanSM=sampling(stanmodels$mSM,data = data.modstanSM,init=init_ll,iter = nriter,chains = nrch,warmup=wu,seed=sd,
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

  # if (bf$bf>10) {
  #   plotSM = function(){
  #   bnds=quantile(as.matrix(fitstanSM)[,6],c(0.025,0.975))
  #   for (di in (2:length(dose.a))){
  #     bnds=cbind(bnds,quantile(as.matrix(fitstanSM)[,length(dose.a)+1+di],c(0.025,0.975)))
  #   }
  #   plot(c(log10(0.00001)-1,log10(1)+1,log10(0.00001)-1,log10(1)+1),log10(c(0,1.25*max(bnds),0,1.25*max(bnds))),type="n",xlab="log10(dose.a+0.00001)",ylab="log10(mean)")
  #   points(log10(dose.a+0.00001),log10(mean.a))
  #   points(log10(dose.a+0.00001),log10(bnds[1,]),pch = 2, col="blue")
  #   points(log10(dose.a+0.00001),log10(bnds[2,]),pch = 6, col="blue")
  #   for (dj in (1:length(dose.a))){
  #     lines(c(log10(dose.a[dj]+0.00001),log10(dose.a[dj]+0.00001)),c(log10(bnds[1,dj]),log10(bnds[2,dj])),col="blue")
  #   }
  #   }
  #   print(plotSM())
  #   return("there is insufficient evidence that there is any dose-effect; therefore no models are fitted and the BMDL is not calculated")
  # }
}
#' Function to determine if there is a dose-response effect (Lognormal version)
#'
#' @param dose.a ordered dose levels
#' @param mean.a mean response
#' @param sd.a standard deviation
#' @param n.a number of observations per dose level
#'
#' # we use the first 5 rows because those are observations from subjects belonging to the same group.
#' data("immunotoxicityData.rda")  #load the immunotoxicity data
#' anydoseresponseLN(dose.a = immunotoxicityData$Dose[1:5],
#'                  mean.a = immunotoxicityData$Mean[1:5],
#'                  sd.a = immunotoxicityData$SD[1:5],
#'                  n.a = immunotoxicityData$n.a[1:5])
#'
#' @description This function tests for dose response effect using Bayes factor.
#'              It fits a null model and a saturated model and compare these two models using model posterior probabilities.
#' @return list containing Baye's factor and decision.
#'
#' @export anydoseresponseLN
#'
anydoseresponseLN=function(dose.a,mean.a,sd.a,n.a){

  N = length(dose.a)
  shift=0
  gmean.a2 = log(NtoLN(mean.a,sd.a))[1:N]
  if (min(gmean.a2)<0) {gmean.a = gmean.a2-20*min(gmean.a2); shift=20*min(gmean.a2)}
  if (min(gmean.a2)>=0) gmean.a = gmean.a2
  gsd.a = log(NtoLN(mean.a,sd.a))[(N+1):(2*N)]

  nrch=3;nriter=3000;wu=1000;dl=0.8;trd=10;sd=123;delta=0.999;treedepth=15
  N = length(dose.a)

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
                        data_type=2, priorg = 4
  )

  data.modstanH0 = list(N=N,n=n.a,m=gmean.a,s2=gsd.a^2,shift=shift,
                        priormu=priorH0$priormu, priorSigma=priorH0$priorSigma,
                        priorlb=priorH0$priorlb, priorub=priorH0$priorub,
                        data_type=2, priorg = 4
  )

  # Fitting the null model (no effect)
  sv=optimizing(stanmodels$mH0,data = data.modstanH0,init=svH0)$par
  initf2 <- function(chain_id = 1) {
    list(par=sv[1:2] + rnorm(2, sd = 0.01*abs(sv[2])) ,alpha = chain_id)
  }
  init_ll <- lapply(1:nrch, function(id) initf2(chain_id = id))
  fitstanH0=sampling(stanmodels$mH0,data = data.modstanH0,init=init_ll,iter = nriter,chains = nrch,warmup=wu,seed=sd,
                     control = list(adapt_delta = delta,max_treedepth =treedepth),
                     show_messages = F, refresh = 0)

  # Fitting the saturated ANOVA model
  sv=optimizing(stanmodels$mSM,data = data.modstanSM,init=svSM)$par
  initf2 <- function(chain_id = 1) {
    list(par=sv[1:(N+1)] + rnorm(N+1, sd = 0.01*abs(sv[1:(N+1)])) ,alpha = chain_id)
  }
  init_ll <- lapply(1:nrch, function(id) initf2(chain_id = id))
  fitstanSM=sampling(stanmodels$mSM,data = data.modstanSM,init=init_ll,iter = nriter,chains = nrch,warmup=wu,seed=sd,
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

  # if (bf$bf>10) {
  #   plotSM = function(){
  #   bnds=quantile(as.matrix(fitstanSM)[,6],c(0.025,0.975))
  #   for (di in (2:length(dose.a))){
  #     bnds=cbind(bnds,quantile(as.matrix(fitstanSM)[,length(dose.a)+1+di],c(0.025,0.975)))
  #   }
  #   plot(c(log10(0.00001)-1,log10(1)+1,log10(0.00001)-1,log10(1)+1),log10(c(0,1.25*max(bnds),0,1.25*max(bnds))),type="n",xlab="log10(dose.a+0.00001)",ylab="log10(mean)")
  #   points(log10(dose.a+0.00001),log10(mean.a))
  #   points(log10(dose.a+0.00001),log10(bnds[1,]),pch = 2, col="blue")
  #   points(log10(dose.a+0.00001),log10(bnds[2,]),pch = 6, col="blue")
  #   for (dj in (1:length(dose.a))){
  #     lines(c(log10(dose.a[dj]+0.00001),log10(dose.a[dj]+0.00001)),c(log10(bnds[1,dj]),log10(bnds[2,dj])),col="blue")
  #   }
  #   }
  #   print(plotSM())
  #   return("there is insufficient evidence that there is any dose-effect; therefore no models are fitted and the BMDL is not calculated")
  # }
}

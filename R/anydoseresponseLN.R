#' Function to determine if there is a dose-response effect
#'
#' @param dose.a dose levels
#' @param mean.a mean response
#' @param sd.a standard deviation
#' @param n.a number of observations per dose level
#'
#' @return text
#'
#' @export
#'
anydoseresponseLNI=function(dose.a,mean.a,sd.a,n.a){

shift=F
gmean.a2 = log(NtoLN(mean.a,sd.a))[1:N]
if (min(gmean.a2)<0) {gmean.a = gmean.a2-1.5*min(gmean.a2); shift=T}
if (min(gmean.a2)>=0) gmean.a = gmean.a2
gsd.a = log(NtoLN(mean.a,sd.a))[(N+1):(2*N)]

nrch=3;nriter=3000;wu=1000;dl=0.8;trd=10;sd=123;delta=0.999;treedepth=15

prior_H0 = function(dose.a,gmean.a,gsd.a,n.a){
  priormu.m=mean(gmean.a)
  priormu.logis2=0
  priors2.m=1000
  priors2.logis2=1000
  return(list(priormu=c(priormu.m,priormu.logis2),priorSigma=diag(c(priors2.m,priors2.logis2))))
}

prior_SM = function(dose.a,gmean.a,gsd.a,n.a){
  priormu.m=mean(gmean.a)
  priormu.logis2=0
  priors2.m=1000
  priors2.logis2=1000
  return(list(priormu=c(rep(priormu.m,length(dose.a)),priormu.logis2),priorSigma=diag(c(rep(priors2.m,length(dose.a)),priors2.logis2))))
}

# Fitting the null model (no effect)
svH0=list(par=c(mean(gmean.a),log(1/mean(gsd.a^2))))
priorH0=prior_H0(dose.a,gmean.a,gsd.a,n.a)
data.modstanH0=list(N=N,n=n.a,m=gmean.a,s2=gsd.a^2,priormu=priorH0$priormu,priorSigma=priorH0$priorSigma)
sv=optimizing(stanmodels$mH0_NI,data = data.modstanH0,init=svH0)$par
initf2 <- function(chain_id = 1) {
  list(par=sv[1:2],alpha = chain_id)
}
init_ll <- lapply(1:nrch, function(id) initf2(chain_id = id))
fitstanH0=sampling(stanmodels$mH0_NI,data = data.modstanH0,init=init_ll,iter = nriter,chains = nrch,warmup=wu,seed=sd,control = list(adapt_delta = delta,max_treedepth =treedepth))

# Fitting the saturated ANOVA model
svSM=list(par=c(gmean.a,log(1/mean(gsd.a^2))))
priorSM=prior_SM(dose.a,gmean.a,gsd.a,n.a)
data.modstanSM=list(N=N,n=n.a,m=gmean.a,s2=gsd.a^2,priormu=priorSM$priormu,priorSigma=priorSM$priorSigma)
sv=optimizing(stanmodels$mSM_NI,data = data.modstanSM,init=svSM)$par
initf2 <- function(chain_id = 1) {
  list(par=sv[1:(N+1)],alpha = chain_id)
}
init_ll <- lapply(1:nrch, function(id) initf2(chain_id = id))
fitstanSM=sampling(stanmodels$mSM_NI,data = data.modstanSM,init=init_ll,iter = nriter,chains = nrch,warmup=wu,seed=sd,control = list(adapt_delta = delta,max_treedepth =treedepth))

set.seed(1234)
bridge_H0 <- bridge_sampler(fitstanH0)
bridge_SM <- bridge_sampler(fitstanSM)
bf=bf(bridge_H0,bridge_SM)
pb=post_prob(bridge_sampler(fitstanH0),bridge_sampler(fitstanSM))
print(bf)
print(pb)

if (bf$bf>10) {
  plotSM = function(){
  bnds=quantile(as.matrix(fitstanSM)[,6],c(0.025,0.975))
  for (di in (2:length(dose.a))){
    bnds=cbind(bnds,quantile(as.matrix(fitstanSM)[,length(dose.a)+1+di],c(0.025,0.975)))
  }
  plot(c(log10(0.00001)-1,log10(1)+1,log10(0.00001)-1,log10(1)+1),log10(c(0,1.25*max(exp(bnds)),0,1.25*max(exp(bnds)))),type="n",xlab="log10(dose.a+0.00001)",ylab="log10(mean)")
  points(log10(dose.a+0.00001),log10(mean.a))
  points(log10(dose.a+0.00001),log10(exp(bnds[1,])),pch = 2, col="blue")
  points(log10(dose.a+0.00001),log10(exp(bnds[2,])),pch = 6, col="blue")
  for (dj in (1:length(dose.a))){
    lines(c(log10(dose.a[dj]+0.00001),log10(dose.a[dj]+0.00001)),c(log10(exp(bnds[1,dj])),log10(exp(bnds[2,dj]))),col="blue")
  }
  }
  print(plotSM())
  return("there is insufficient evidence that there is any dose-effect; therefore no models are fitted and the BMDL is not calculated")
}
if (bf$bf<10)   return("there is sufficient evidence that there is a substantial dose-effect; therefore models are fitted and the BMDL is calculated")
}
#' @rdname anydoseresponseLNI
#' @export
anydoseresponseLND=function(dose.a,mean.a,sd.a,n.a){

shift=F
gmean.a2 = log(NtoLN(mean.a,sd.a))[1:N]
if (min(gmean.a2)<0) {gmean.a = gmean.a2-1.5*min(gmean.a2); shift=T}
if (min(gmean.a2)>=0) gmean.a = gmean.a2
gsd.a = log(NtoLN(mean.a,sd.a))[(N+1):(2*N)]

nrch=3;nriter=3000;wu=1000;dl=0.8;trd=10;sd=123;delta=0.999;treedepth=15

prior_H0 = function(dose.a,gmean.a,gsd.a,n.a){
  priormu.m=mean(gmean.a)
  priormu.logis2=0
  priors2.m=1000
  priors2.logis2=1000
  return(list(priormu=c(priormu.m,priormu.logis2),priorSigma=diag(c(priors2.m,priors2.logis2))))
}

prior_SM = function(dose.a,gmean.a,gsd.a,n.a){
  priormu.m=mean(gmean.a)
  priormu.logis2=0
  priors2.m=1000
  priors2.logis2=1000
  return(list(priormu=c(rep(priormu.m,length(dose.a)),priormu.logis2),priorSigma=diag(c(rep(priors2.m,length(dose.a)),priors2.logis2))))
}

# Fitting the null model (no effect)
svH0=list(par=c(mean(gmean.a),log(1/mean(gsd.a^2))))
priorH0=prior_H0(dose.a,gmean.a,gsd.a,n.a)
data.modstanH0=list(N=N,n=n.a,m=gmean.a,s2=gsd.a^2,priormu=priorH0$priormu,priorSigma=priorH0$priorSigma)
sv=optimizing(stanmodels$mH0_ND,data = data.modstanH0,init=svH0)$par
initf2 <- function(chain_id = 1) {
  list(par=sv[1:2],alpha = chain_id)
}
init_ll <- lapply(1:nrch, function(id) initf2(chain_id = id))
fitstanH0=sampling(stanmodels$mH0_ND,data = data.modstanH0,init=init_ll,iter = nriter,chains = nrch,warmup=wu,seed=sd,control = list(adapt_delta = delta,max_treedepth =treedepth))

# Fitting the saturated ANOVA model
svSM=list(par=c(gmean.a,log(1/mean(gsd.a^2))))
priorSM=prior_SM(dose.a,gmean.a,gsd.a,n.a)
data.modstanSM=list(N=N,n=n.a,m=gmean.a,s2=gsd.a^2,priormu=priorSM$priormu,priorSigma=priorSM$priorSigma)
sv=optimizing(stanmodels$mSM_ND,data = data.modstanSM,init=svSM)$par
initf2 <- function(chain_id = 1) {
  list(par=sv[1:(N+1)],alpha = chain_id)
}
init_ll <- lapply(1:nrch, function(id) initf2(chain_id = id))
fitstanSM=sampling(stanmodels$mSM_ND,data = data.modstanSM,init=init_ll,iter = nriter,chains = nrch,warmup=wu,seed=sd,control = list(adapt_delta = delta,max_treedepth =treedepth))

set.seed(1234)
bridge_H0 <- bridge_sampler(fitstanH0)
bridge_SM <- bridge_sampler(fitstanSM)
bf=bf(bridge_H0,bridge_SM)
pb=post_prob(bridge_sampler(fitstanH0),bridge_sampler(fitstanSM))
print(bf)
print(pb)

if (bf$bf>10) {
  plotSM = function(){
  bnds=quantile(as.matrix(fitstanSM)[,6],c(0.025,0.975))
  for (di in (2:length(dose.a))){
    bnds=cbind(bnds,quantile(as.matrix(fitstanSM)[,length(dose.a)+1+di],c(0.025,0.975)))
  }
  plot(c(log10(0.00001)-1,log10(1)+1,log10(0.00001)-1,log10(1)+1),log10(c(0,1.25*max(exp(bnds)),0,1.25*max(exp(bnds)))),type="n",xlab="log10(dose.a+0.00001)",ylab="log10(mean)")
  points(log10(dose.a+0.00001),log10(mean.a))
  points(log10(dose.a+0.00001),log10(exp(bnds[1,])),pch = 2, col="blue")
  points(log10(dose.a+0.00001),log10(exp(bnds[2,])),pch = 6, col="blue")
  for (dj in (1:length(dose.a))){
    lines(c(log10(dose.a[dj]+0.00001),log10(dose.a[dj]+0.00001)),c(log10(exp(bnds[1,dj])),log10(exp(bnds[2,dj]))),col="blue")
  }
  }
  print(plotSM())
  return("there is insufficient evidence that there is any dose-effect; therefore no models are fitted and the BMDL is not calculated")
}
if (bf$bf<10)   return("there is sufficient evidence that there is a substantial dose-effect; therefore models are fitted and the BMDL is calculated")
}


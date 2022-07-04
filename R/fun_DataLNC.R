
#' Function to set data in the correct format, and also to generate starting values. Lognormal version
#'
#' This function also generates appropriate start values and uninformative priors for each model
#'
#'
#' @param data dataframe with input data, order of columns should be: dose, response, sd, n.
#' @param sumstats logical. TRUE indicates summary data is provided while FALSE indicates individual-level data.
#' @param geom.stats logical. TRUE if geometric summary data is provided.
#' @param sd logical. TRUE indicates that standard deviation per dose is provided,
#'           FALSE indicates that standard error per dose level is provided
#' @param q the specified BMR
#' @param bkg vector containing informative prior for the background.
#'            It should be specified as minimum, most likely and maximum. Defaults to NULL.
#' @param maxy vector containing informative prior for the maximum response.
#'             It should be specified as minimum, most likely and maximum. Defaults to NULL.
#' @param prior.BMD vector containing informative prior for the BMD.
#'                  It should be specified as minimum, most likely and maximum. Defaults to NULL.
#' @param shape.a shape parameter that determines the flatness of the Pert prior for the background.
#'                 Defaults to 4, which implies a peak at the most likely value.
#' @param shape.c shape parameter that determines the flatness of the Pert prior for c.
#'                Defaults to 4, which implies a peak at the most likely value.
#' @param shape.BMD shape parameter that determines the flatness of the Pert prior for the BMD.
#'                  Defaults to 0.0001, which implies a flat prior.

#'
#' @return List with data and start values in correct format to be directly used within the BMA functions.
#'
#' @export PREP_DATA_LN_C
#'
PREP_DATA_LN_C <- function(data, # a dataframe with input data, order of columns should be: dose, response, litter
                           q, # the BMR
                           bkg = NULL, maxy = NULL, prior.BMD = NULL, # possible expert info on background and max response
                           shape.a = 4, shape.c = 4, shape.BMD = 0.0001, # shape for the PERT distribution
                           prior.d = 'N11',
                           extended = T
){


  ## cluster (ij) = combination of dose (i) and litter (j)
  library(dplyr)
  indiv.data <- data.frame(dose = data[,1],
                           response = data[,2],
                           litter = data[,3])
  indiv.data <- indiv.data %>%
    dplyr::group_by(dose, litter) %>%
    dplyr::arrange(by_group = dose)
  indiv.data <- indiv.data %>%
    dplyr::mutate(cluster = dplyr::cur_group_id(),
                  count = dplyr::n())
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
  yl = matrix(0, nrow = nc, ncol = maxN)
  obs.o = c(); obs.l = c()
  for(i in 1:nc){
    obs = indiv.data$response[indiv.data$cluster==cid[i]]
    if(length(obs) < maxN){
      obs.o = c(obs, rep(0, maxN-length(obs)))
      obs.l = c(log(obs), rep(0, maxN-length(obs)))
    }else{
      obs.o = obs
      obs.l = log(obs)
    }
    y[i, ] = obs.o
    yl[i, ] = obs.l
  }
  if(sum(yl<0)!=0){
    yl2 = yl - 20*min(yl)
    shift = 20*min(yl)
  }else{
    yl2 = yl
    shift = 0
  }

  ## Bartlett test of homoscedasticity
  # on the log scale (constant coefficient of variation)
  # b.test.LN <- bartlett.test(log(indiv.data$response), indiv.data$dose)
  # # b.test.LN <- bartlett(gsd.a, n.a)
  # if(b.test.LN$p.value>=0.05){
  #   test.var = paste0('Distributional assumption of constant variance (on log-scale) are met, Bartlett test p-value is ', round(b.test.LN$p.value, 4))
  #   message(test.var)
  # }else if(b.test.LN$p.value<0.05){
  #   test.var = paste0('Distributional assumption of constant variance (on log-scale) are not met, Bartlett test p-value is ', round(b.test.LN$p.value, 4))
  #   warning(test.var)
  # }

  is_informative_BMD = 0
  is_informative_a = 0
  is_informative_c = 0

  if(!is.null(prior.BMD)) {is_informative_BMD = 1}

  ## Overall mean for test of flatness
  means.all <- indiv.data %>%
    group_by(dose) %>%
    summarise(mresp = mean(response))
  dose.a = unique(indiv.data$dose)
  mean.a = c()
  for(m in 1:length(dose.a)){
    mean.a[m] <- means.all$mresp[means.all$dose == dose.a[m]]
  }


  ######################
  ### PRIORS

  ## Observed background and maximum response

  obs.min = mean(indiv.data$response[indiv.data$dose==0])
  obs.max = mean(indiv.data$response[indiv.data$dose==maxDose])

  ### Family 1

  ## Default range on background and max response

  if(obs.min < obs.max){

    # for a
    min.min = 0.001
    mode.min = obs.min
    max.min = 2*obs.min

    # for c
    min.max = obs.min*(1.01+q)
    max.max = 2*obs.max
    mode.max = obs.max

    if(flat(dose.a, mean.a,inc=T) == F & is.null(maxy)){
      mode.max = 3*obs.max
      min.max  = obs.min*(1.01+q)
      max.max = 2*mode.max

      # warning(
      #   "The data do not contain information on the asymptote, and the default prior for fold change has been set to 3 times the observed maximum.
      #       Please provide prior input on the maximum response by specifying 'maxy', if available."
      # )
    }

    # start value BMD
    datf = data.frame(yy=data$response,xx=(data$dose/max(data$dose))+0.00000000000001)
    fpfit=gamlss::gamlss(yy~fp(xx),family=NO(),data=datf)
    RISK=function(x) (predict(fpfit,newdata=data.frame(xx=c(exp(x))))-predict(fpfit,newdata=data.frame(xx=c(0.00000000000001))))/
      (predict(fpfit,newdata=data.frame(xx=c(0.00000000000001))))-q
    bmd.svh=try(uniroot(RISK, interval=c(0, 1))$root,silent=T)
    # bmd.sv=ifelse((mode(bmd.svh)=="numeric"),exp(bmd.svh),0.05)
    bmd.sv=ifelse((mode(bmd.svh)=="numeric"),bmd.svh,0.5)


    ## Check appropriateness of BMR value
    if(obs.min*(1+q) > obs.max){
      warning('The data do not contain values corresponding to the chosen BMR, lowering the specified value of q may be necessary.')
    }



  }else if(obs.min > obs.max){
    min.min = 0.5*obs.min
    mode.min = obs.min
    max.min = 2*obs.min

    min.max = 0.5*obs.max
    max.max = obs.min*(1-q-0.01)
    # max.max = obs.min
    mode.max = obs.max

    if(flat(dose.a, mean.a,inc=F) == F & is.null(maxy)){
      mode.max = 0.5*obs.max
      min.max  = 0.1*obs.max
      max.max = obs.min*(1-q-0.01)
      # max.max = obs.min

      # warning(
      #   "The data do not contain information on the asymptote, and the default prior for fold change has been set to 3 times the observed maximum.
      #       Please provide prior input on the maximum response by specifying 'maxy', if available."
      # )
    }

    # start value BMD
    datf = data.frame(yy=data$response,xx=(data$dose/max(data$dose))+0.00000000000001)
    fpfit=gamlss::gamlss(yy~fp(xx),family=NO(),data=datf)
    RISK=function(x) (predict(fpfit,newdata=data.frame(xx=c(exp(x))))-predict(fpfit,newdata=data.frame(xx=c(0.00000000000001))))/
      (predict(fpfit,newdata=data.frame(xx=c(0.00000000000001))))+q
    bmd.svh=try(uniroot(RISK, interval=c(0, 1))$root,silent=T)
    bmd.sv=ifelse((mode(bmd.svh)=="numeric"),bmd.svh,0.5)


    ## Check appropriateness of BMR value
    if(obs.min*(1-q) < obs.max){
      warning('The data do not contain values corresponding to the chosen BMR, lowering the specified value of q may be necessary.')
    }

  }

  ## If info on background is given
  if(!is.null(bkg)){

    is_informative_a = 1

    if(!is.na(bkg[2])){
      mode.min = bkg[2]
    }else{
      mode.min = bkg[1] + ((bkg[3]-bkg[1])/2)
    }
    if(!is.na(bkg[1])){
      min.min = bkg[1]
    }
    if(!is.na(bkg[3])){
      max.min = bkg[3]
    }

  }else{
    message("Default prior choices used on background")
  }

  ## If info on max response is given
  if(!is.null(maxy)){

    is_informative_c = 1

    if(!is.na(maxy[2])){
      mode.max = maxy[2]
    }else{
      mode.max = maxy[1] + ((maxy[3]-maxy[1])/2)
    }
    if(!is.na(maxy[1])){
      min.max = maxy[1]
    }
    if(!is.na(maxy[3])){
      max.max = maxy[3]
    }

  }else{
    message("Default prior choices used on fold change")
  }

  ## Default prior BMD
  BMD.min <- .Machine$double.xmin
  if(extended == FALSE){
    BMD.max <- 1
  }else{
    BMD.max <- maxDose
    if(maxDose <= 1){
      BMD.max <- maxDose*1000
    }
  }
  BMD.mode <- 0.5

  ## If info on BMD is given
  if(!is.null(prior.BMD)){

    if(!is.na(prior.BMD[2])){
      BMD.mode = prior.BMD[2]/maxDose
    }else{
      BMD.mode = (prior.BMD[1]/maxDose) + (((prior.BMD[3]/maxDose) - prior.BMD[1]/maxDose)/2)
    }
    if(!is.na(prior.BMD[1])){
      BMD.min = prior.BMD[1]/maxDose
    }
    if(!is.na(prior.BMD[3])){
      BMD.max = prior.BMD[3]/maxDose
    }

  }else {
    message("Default prior choices used on BMD")
  }

  BMD.vec <- c(BMD.min, BMD.mode, BMD.max)

  if(prior.d == 'N11'){
    prvar.d = 1; prmean.d = 1; truncd = 5
  }else if(prior.d == 'EPA'){
    # prvar.d = 0.5^2; prmean.d = 0.4; truncd = 10000
    prvar.d = 0.5; prmean.d = 0.4; truncd = 10000
  }
  prmean.dQE4 = 0; prvar.dQE4 = 1; truncdQ = 10000
  prvar.s=1; prmean.s=-2*log(1.5*sd(yl2[yl2!=0]))

  # Prior on background
  a.vec = c(min.min, mode.min, max.min)

  # Prior on mu(inf)
  c.vec = c(min.max/mode.min, mode.max/mode.min, max.max/mode.min)
  if(c.vec[1] == 0) c.vec[1] = 0.0001
  if(c.vec[2] >= c.vec[3]) c.vec[2] = c.vec[3] - 0.05
  # if(c.vec[2] >= c.vec[3]) c.vec[2] = c.vec[1] + (c.vec[3] - c.vec[1])/2 - 0.01


  priormu1a=c(a.vec[2],BMD.vec[2],c.vec[2],prmean.d,prmean.s,0.5)
  priormu1bQ=c(a.vec[2],BMD.vec[2],c.vec[2],prmean.dQE4,prmean.s,0.5)
  priorSigma1a=diag(c(1,1,1,prvar.d,prvar.s))
  priorSigma1bQ=diag(c(1,1,1,prvar.dQE4,prvar.s))
  priorlb1a = c(a.vec[1],BMD.vec[1],c.vec[1],0,0,0)
  priorub1a = c(a.vec[3],BMD.vec[3],c.vec[3],0,0,1)


  if(!is.null(prior.BMD)) bmd.sv = BMD.vec[2]

  is_increasing = 0; is_decreasing = 0
  if(obs.min < obs.max){ # increasing
    is_increasing = 1
    L = 1+q+0.01
    U = 0
    data_type = 2
    pars3d = numeric()
    pars3i = priormu1a[3] - L
    dim(pars3d)=0
    dim(pars3i)=1
  }else if(obs.min > obs.max){ # decreasing
    is_decreasing = 1
    L = 0.01
    U = 1-q-0.01
    data_type = 4
    pars3d = priormu1a[3] / U
    # pars3d =  L + (U - L) * priormu1a[3];
    pars3i = numeric()
    dim(pars3d)=1
    dim(pars3i)=0
  }


  ## Test for homoscedasticity and normality of residuals, based on saturated ANOVA model

  priorSM = list(
    priormu = c(mean.a[1],
                diff(mean.a),
                -2*log(1.5*sd(yl2[yl2!=0])),
                0.5),
    priorSigma = diag(c(1, rep(1, length(dose.a)-1), 1)),
    priorlb = 0.001,
    priorub = c(2*mean.a[1],
                max(abs(diff(mean.a)))*10)
  )

  svSM = list(par = c(mean.a[1],
                      diff(mean.a),
                      log(1/var(yl2[yl2!=0])),
                      0.5))

  data.modstanSM = list(N=N, n=n, nc=nc, maxN=maxN, maxNc=maxNc,
                        nij=nij, y=yl2, q=q, shift=shift,
                        priormu=priorSM$priormu, priorSigma=priorSM$priorSigma,
                        priorlb=priorSM$priorlb, priorub=priorSM$priorub,
                        priorg=4, data_type=data_type
  )

  optSM = optimizing(stanmodels$mSMc, data = data.modstanSM,
                     seed=as.integer(123), draws = 30000,
                     init = svSM, hessian=TRUE)
  means.SM = apply(optSM$theta_tilde[, paste0('mu[', 1:N, ']')], 2, median)

  means.all$pred = means.SM

  indiv.data$res = NA
  for(i in 1:dim(indiv.data)[1]){
    indiv.data$res[i] <- indiv.data$response[i] - (means.all$pred[which(means.all$dose == indiv.data$dose[i])])
  }

  # homoscedasticity
  b.test.N <- bartlett.test(indiv.data$res, indiv.data$dose)
  # b.test.N <- bartlett(sd.a, n.a)
  if(b.test.N$p.value>=0.05){
    test.var = paste0('Distributional assumption of constant variance for the lognormal distribution is met, Bartlett test p-value is ', round(b.test.N$p.value, 4))
    message(test.var)
  }else if(b.test.N$p.value<0.05){
    test.var = paste0('Distributional assumption of constant variance for the lognormal distribution is not met, Bartlett test p-value is ', round(b.test.N$p.value, 4))
    warning(test.var)
  }

  # normality
  norm.test.N <- shapiro.test(indiv.data$res)
  # b.test.N <- bartlett(sd.a, n.a)
  if(norm.test.N$p.value>=0.05){
    test.var = paste0('Distributional assumption of normality of residuals for the lognormal distribution is met, Shapiro test p-value is ',
                      round(norm.test.N$p.value, 4))
    message(test.var)
  }else if(norm.test.N$p.value<0.05){
    test.var = paste0('Distributional assumption of normality of residuals for the lognormal distribution is not met, Shapiro test p-value is ',
                      round(norm.test.N$p.value, 4))
    warning(test.var)
  }

  ## Data in correct format

  ret.list <- list(data = list(N=N,
                               n=n,
                               x=doses/maxDose,
                               # m=mean.a,
                               nc=nc,
                               maxN=maxN,
                               maxNc=maxNc,
                               nij=nij,
                               y=yl2,
                               shift=shift,
                               data = indiv.data,
                               # s2=gsd.a^2,
                               maxD=maxDose,q=q,priormu=priormu1a,priormuQ = priormu1bQ,
                               shape1 = c(fun.alpha(a = priorlb1a[1], b = priormu1a[1], c = priorub1a[1], g = shape.a),
                                          fun.alpha(a = priorlb1a[2], b = priormu1a[2], c = priorub1a[2], g = shape.BMD),
                                          fun.alpha(a = priorlb1a[3], b = priormu1a[3], c = priorub1a[3], g = shape.c), 0, 0,
                                          fun.alpha(a = priorlb1a[6], b = priormu1a[6], c = priorub1a[6], g = 0.0001)),
                               shape2 = c(fun.beta(a = priorlb1a[1], b = priormu1a[1], c = priorub1a[1], g = shape.a),
                                          fun.beta(a = priorlb1a[2], b = priormu1a[2], c = priorub1a[2], g = shape.BMD),
                                          fun.beta(a = priorlb1a[3], b = priormu1a[3], c = priorub1a[3], g = shape.c), 0, 0,
                                          fun.alpha(a = priorlb1a[6], b = priormu1a[6], c = priorub1a[6], g = 0.0001)),
                               priorlb = priorlb1a, priorub = priorub1a, shape.a = shape.a, shape.c = shape.c, shape.BMD = shape.BMD,
                               priorSigma=priorSigma1a, priorSigmaQ = priorSigma1bQ, truncd = truncd, truncdQ = truncdQ, init_b = 1, data_type = data_type, L = L, U = U, is_increasing = is_increasing,
                               is_decreasing = is_decreasing, is_informative_a = is_informative_a, is_informative_c = is_informative_c,
                               is_informative_BMD = is_informative_BMD),
                   start = list(par1=priormu1a[1],par2=bmd.sv,pars3i=pars3i,pars3d=pars3d,par4=prmean.d,par5=log(1/(var(yl2[yl2!=0]))),
                                par6=0.5),
                   startQ = list(par1=priormu1a[1],par2=bmd.sv,pars3i=pars3i,pars3d=pars3d,par4=prmean.dQE4,par5=log(1/(var(yl2[yl2!=0]))),
                                 par6=0.5),
                   test.var = test.var
  )




  # test for dose-response effect
  # DR.effect = anydoseresponseNI(dose.a,mean.a,sd.a,n.a)

  # data in correct format
  return(ret.list)

}

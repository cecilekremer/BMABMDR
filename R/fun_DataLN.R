
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
#' @export PREP_DATA_LN
#'
PREP_DATA_LN <- function(data, # a dataframe with input data, order of columns should be: dose, response, sd, n
                         sumstats = TRUE, # TRUE if summary data, FALSE if individual data
                         geom.stats = FALSE, # TRUE if geometric summary data is provided
                         sd = TRUE, # TRUE if sd per dose group is given, FALSE is se is given
                         q, # the BMR
                         bkg = NULL, maxy = NULL, prior.BMD = NULL, # possible expert info on background and max response
                         shape.a = 4, shape.c = 4, shape.BMD = 0.0001 # shape for the PERT distribution
){


  if(sumstats == TRUE & geom.stats == FALSE){
    data = data[order(data[, 1]), ]
    dose.a = data[, 1]
    maxDose = max(dose.a)
    mean.a = data[, 2]
    if(sd == TRUE){
      sd.a = data[, 3]
    }else if(sd == FALSE){
      sd.a = data[,3]*sqrt(data[, 4]) # SD = SE * sqrt(n.a)
    }
    n.a = data[, 4]
    N = length(dose.a)
    dose.a = dose.a/maxDose
    # shift if negative means occur
    shift=0
    gmean.a2 = log(NtoLN(mean.a,sd.a))[1:N]
    if (min(gmean.a2)<0) {gmean.a = gmean.a2-20*min(gmean.a2); shift = 20*min(gmean.a2)}
    if (min(gmean.a2)>=0) gmean.a = gmean.a2
    gsd.a = log(NtoLN(mean.a,sd.a))[(N+1):(2*N)]

  }else if(sumstats == TRUE & geom.stats == TRUE){


  }else if(sumstats == FALSE){
    data = data[order(data[, 1]), ]
    doses = data[, 1]
    maxDose = max(doses)
    dose.a=sort(unique(doses))
    N=length(dose.a)
    mean.a=rep(NA,N)
    sd.a=rep(NA,N)
    n.a=rep(NA,N)
    y = data[, 2]
    for (iu in (1:N)){
      mean.a[iu]=exp(mean(log(y[doses==dose.a[iu]])))
      sd.a[iu]=exp(sd(log(y[doses==dose.a[iu]])))
      n.a[iu]=sum(doses==dose.a[iu])
    }
    dose.a = dose.a/maxDose
    # shift if negative means occur
    shift=0
    gmean.a2 = log(mean.a)
    if (min(gmean.a2)<0) {gmean.a = gmean.a2-20*min(gmean.a2); shift = 20*min(gmean.a2)}
    if (min(gmean.a2)>=0) gmean.a = gmean.a2
    gsd.a = log(sd.a)
  }

  ## Bartlett test of homoscedasticity
  # on the log scale (constant coefficient of variation)
  b.test.LN <- bartlett(gsd.a, n.a)
  if(b.test.LN[2]>=0.05){
    test.var = paste0('Distributional assumption of constant variance (on log-scale) are met, Bartlett test p-value is ', round(b.test.LN[2], 4))
    message(test.var)
  }else if(b.test.LN[2]<0.05){
    test.var = paste0('Distributional assumption of constant variance (on log-scale) are not met, Bartlett test p-value is ', round(b.test.LN[2], 4))
    warning(test.var)
  }

  is_informative_BMD = 0
  is_informative_a = 0
  is_informative_c = 0

  if(!is.null(prior.BMD)) {is_informative_BMD = 1}


  ######################
  ### PRIORS

  ## Observed background and maximum response

  obs.min = mean.a[1]
  obs.max = mean.a[N]

  ### Family 1

  ## Default range on background and max response

  if(mean.a[1] < mean.a[N]){
    min.min = 0.001
    mode.min = obs.min
    max.min = 2*obs.min

    min.max = mean.a[1]*(1.01+q)
    max.max = 2*mean.a[N]
    mode.max = obs.max

    if(flat(dose.a, mean.a,inc=T) == F & is.null(maxy)){
      mode.max = 3*mean.a[N]
      min.max  = mean.a[1]*(1.01+q)
      max.max = 2*mode.max

      # warning(
      #   "The data do not contain information on the asymptote, and the default prior for fold change has been set to 3 times the observed maximum.
      #       Please provide prior input on the maximum response by specifying 'maxy', if available."
      # )
    }

    # start value BMD
    datf=data.frame(yy=mean.a,xx=dose.a+0.00000000000001)
    fpfit=gamlss(yy~fp(xx),family=NO(),data=datf)
    RISK=function(x) (predict(fpfit,newdata=data.frame(xx=c(exp(x))))-predict(fpfit,newdata=data.frame(xx=c(0.00000000000001))))/
      (predict(fpfit,newdata=data.frame(xx=c(0.00000000000001))))-q
    bmd.svh=try(uniroot(RISK, interval=c(-5, 0))$root,silent=T)
    bmd.sv=ifelse((mode(bmd.svh)=="numeric"),exp(bmd.svh),0.05)

    ## Check appropriateness of BMR value
    if(mean.a[1]*(1+q) > mean.a[N]){
      warning('The data do not contain values corresponding to the chosen BMR, lowering the specified value of q may be necessary.')
    }

  }else if(mean.a[1] > mean.a[N]){
    min.min = 0.5*obs.min
    mode.min = obs.min
    max.min = 2*obs.min

    mode.max = obs.max
    min.max = 0.5*obs.max
    max.max = obs.min*(1-q-0.01)

    if(flat(dose.a, mean.a,inc=F) == F & is.null(maxy)){
      mode.max = 0.5*obs.max
      min.max  = 0.1*obs.max
      max.max = obs.min*(1-q-0.01)

      # warning(
      #   "The data do not contain information on the asymptote, and the default prior for fold change has been set to 3 times the observed maximum.
      #       Please provide prior input on the maximum response by specifying 'maxy', if available."
      # )
    }

    # start value BMD
    datf=data.frame(yy=mean.a,xx=dose.a+0.00000000000001)
    fpfit=gamlss(yy~fp(xx),family=NO(),data=datf)
    RISK=function(x) (predict(fpfit,newdata=data.frame(xx=c(exp(x))))-predict(fpfit,newdata=data.frame(xx=c(0.00000000000001))))/
      (predict(fpfit,newdata=data.frame(xx=c(0.00000000000001))))+q
    bmd.svh=try(uniroot(RISK, interval=c(-5, 0))$root,silent=T)
    bmd.sv=ifelse((mode(bmd.svh)=="numeric"),exp(bmd.svh),0.05)

    ## Check appropriateness of BMR value
    if(mean.a[1]*(1-q) < mean.a[N]){
      warning('The data do not contain values corresponding to the chosen BMR, lowering the specified value of q may be necessary.')
    }

  }


  ## If info on background is given
  if(!is.null(bkg)){

    is_informative_a = 1

    if(!is.na(bkg[2])){
      mode.min = bkg[2]
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
  BMD.max <- 1
  BMD.mode <- 0.5

  ## If info on BMD is given
  if(!is.null(prior.BMD)){

    if(!is.na(prior.BMD[2])){
      BMD.mode = prior.BMD[2]/maxDose
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

  prvar.d=1; prmean.d=0
  # prvar.d=(exp(sqrt(0.18)))^2; prmean.d=2
  prvar.s=1; prmean.s=-2*log(1.5*mean(gsd.a))

  # Prior on background
  a.vec = c(min.min, mode.min, max.min)

  # Prior on mu(inf)
  c.vec = c(min.max/mode.min, mode.max/mode.min, max.max/mode.min)
  if(c.vec[1] == 0) c.vec[1] = 0.0001
  if(c.vec[2] >= c.vec[3]) c.vec[2] = c.vec[3] - 0.05


  priormu1a=c(a.vec[2],BMD.vec[2],c.vec[2],prmean.d,prmean.s)
  priorSigma1a=diag(c(1,1,1,prvar.d,prvar.s))
  priorlb1a = c(a.vec[1],BMD.vec[1],c.vec[1],0,0)
  priorub1a = c(a.vec[3],BMD.vec[3],c.vec[3],0,0)


  if(!is.null(prior.BMD)) bmd.sv = BMD.vec[2]

  is_increasing = 0; is_decreasing = 0
  if(mean.a[1] < mean.a[N]){ # increasing
    is_increasing = 1
    L = 1+q+0.01
    U = 0
    data_type = 2
    pars3d = numeric()
    pars3i = priormu1a[3] - L
    dim(pars3d)=0
    dim(pars3i)=1
  }else if(mean.a[1] > mean.a[N]){ # decreasing
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

  ## Data in correct format

  ret.list <- list(data = list(N=N,n=n.a,x=dose.a,m=gmean.a,m.org=gmean.a2,shift=shift,s2=gsd.a^2,maxD=maxDose,q=q,priormu=priormu1a,
                               shape1 = c(fun.alpha(a = priorlb1a[1], b = priormu1a[1], c = priorub1a[1], g = shape.a),
                                          fun.alpha(a = priorlb1a[2], b = priormu1a[2], c = priorub1a[2], g = shape.BMD),
                                          fun.alpha(a = priorlb1a[3], b = priormu1a[3], c = priorub1a[3], g = shape.c), 0, 0),
                               shape2 = c(fun.beta(a = priorlb1a[1], b = priormu1a[1], c = priorub1a[1], g = shape.a),
                                          fun.beta(a = priorlb1a[2], b = priormu1a[2], c = priorub1a[2], g = shape.BMD),
                                          fun.beta(a = priorlb1a[3], b = priormu1a[3], c = priorub1a[3], g = shape.c), 0, 0),
                               priorlb = priorlb1a, priorub = priorub1a, shape.a = shape.a, shape.c = shape.c, shape.BMD = shape.BMD,
                               priorSigma=priorSigma1a, init_b = 1, data_type = data_type, L = L, U = U, is_increasing = is_increasing,
                               is_decreasing = is_decreasing, is_informative_a = is_informative_a, is_informative_c = is_informative_c,
                               is_informative_BMD = is_informative_BMD),
                   start = list(par1=priormu1a[1],par2=bmd.sv,pars3i=pars3i,pars3d=pars3d,par4=0,par5=log(1/mean(gsd.a^2))),
                   test.var = test.var
  )


  # data in correct format
  return(ret.list)

}

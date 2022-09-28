
#' Function to set data in the correct input format for model-averaging functions.
#' This function also generates appropriate start values and uninformative priors for each model
#'
#' @param data a dataframe with input data, order of columns should be: dose, response, sd (or se), n, covariate (if given)
#' @param sumstats logical indicating whether summary (T, default) or individual-level (F) data is provided
#' @param geom.stats logicial indicating whether, if summary data are provided, these are geometric (T) or arithmetic (F, default) summary statistics
#' @param sd logical indicating whether standard deviation (T, default) or standard error (F) is provided
#' @param q specified BMR
#' @param bkg vector containing minimum, most likely (optional, can be NA), and maximum value for the background response. Defaults to NULL (non-informative prior)
#' @param maxy vector containing minimum, most likely (optional, can be NA), and maximum value for the response at dose infinity. Defaults to NULL (non-informative prior)
#' @param prior.BMD vector containing minimum, most likely (optional, can be NA), and maximum value for the BMD. Defaults to NULL (non-informative prior)
#' @param shape.a shape parameter for the modified PERT distribution on parameter a, defaults to 4 (peaked at most likely value), a value of 0.0001 implies a uniform distribution
#' @param shape.c shape parameter for the modified PERT distribution on parameter c, defaults to 4 (peaked at most likely value), a value of 0.0001 implies a uniform distribution
#' @param shape.BMD shape parameter for the modified PERT distribution on parameter BMD, defaults to 0.0001 implying a uniform distribution. Can be set to 4 in case of informative prior
#' @param prior.d prior distribution for parameter d, should be either N11 (default), EPA or N05 (for a N(0.5,0.5) prior)
#' @param extended logical indicating whether the dose range should be extended to maxDose^2 (default is TRUE)
#'
#' @description The function a dataset as input and generates the data list and starting values needed by the stan
#'              models to be fitted.
#'
#' @examples
#'
#'  # we use the first 5 rows because those are observations from subjects belonging to the same group.
#'  data("immunotoxicityData.rda")  #load the immunotoxicity data
#'  data_N <- PREP_DATA_N(data = as.data.frame(immunotoxicityData[1:5,]),
#'  sumstats = TRUE, sd = TRUE, q = 0.1) #example with default priors
#'
#'  data_N <- PREP_DATA_N(data = as.data.frame(immunotoxicityData[1:5,]), sumstats = TRUE,
#'                        sd = TRUE, q = 0.1, bkg = c(0.62, 1.34, 2.06)) #example with informative prior on background
#'
#'
#'  data_N <- PREP_DATA_N(data = as.data.frame(immunotoxicityData[1:5,]), sumstats = TRUE,
#'                        sd = TRUE, q = 0.1,
#'                        prior.BMD = c(0.06, 0.25, 1)) #example with informative priors on the BMD
#'
#'
#' @return List with data and start values in correct format to be directly used within the model-averaging functions.
#'
#' @export PREP_DATA_N
#'
PREP_DATA_N <- function(data, # a dataframe with input data, order of columns should be: dose, response, sd, n
                        sumstats = TRUE, # TRUE if summary data, FALSE if individual data
                        geom.stats = FALSE, # TRUE if geometric summary data is provided
                        sd = TRUE, # TRUE if sd per dose group is given, FALSE is se is given
                        q, # the BMR
                        bkg = NULL,
                        maxy = NULL,
                        prior.BMD = NULL, # possible expert info on background and max response
                        shape.a = 4, shape.c = 4, shape.BMD = 0.0001, # shape for the PERT distribution,
                        prior.d = 'N11',
                        extended = TRUE
){

  if(sumstats == TRUE & geom.stats == FALSE){
    data = data[order(data[, 1]), ]
    dose.a = data[, 1]
    maxDose = max(dose.a)

    ## if dose levels not unique
    if(length(dose.a) != length(unique(dose.a))){
      dose.a = sort(unique(dose.a))
      N = length(dose.a)
      mean.a=rep(NA,N)
      sd.a=rep(NA,N)
      n.a=rep(NA,N)
      for (iu in (1:N)){
        mean.a[iu] = mean(data[,2][data[,1] == dose.a[iu]])
        sd.a[iu] = mean(data[,3][data[,1] == dose.a[iu]])
        n.a[iu] = sum(data[,4][data[,1] == dose.a[iu]])
      }
      dose.a = dose.a/maxDose

    }else{
      mean.a = data[, 2]
      if(sd == TRUE){
        sd.a = data[, 3]
      }else if(sd == FALSE){
        sd.a = data[,3]*sqrt(data[, 4]) # SD = SE * sqrt(n.a)
      }
      n.a = data[, 4]
      N = length(dose.a)
      dose.a = dose.a/maxDose
    }
    testNLN <- NA
  }else if(sumstats == TRUE & geom.stats == TRUE){
    data = data[order(data[, 1]), ]
    dose.a = data[, 1]
    maxDose = max(dose.a)

    ## if dose levels not unique
    if(length(dose.a) != length(unique(dose.a))){
      dose.a = sort(unique(dose.a))
      N = length(dose.a)
      gmean.a=rep(NA,N)
      gsd.a=rep(NA,N)
      n.a=rep(NA,N)
      for (iu in (1:N)){
        gmean.a[iu] = mean(data[,2][data[,1] == dose.a[iu]])
        gsd.a[iu] = mean(data[,3][data[,1] == dose.a[iu]])
        n.a[iu] = sum(data[,4][data[,1] == dose.a[iu]])
      }
      dose.a = dose.a/maxDose

    }else{
      gmean.a = data[, 2]
      if(sd == TRUE){
        gsd.a = data[, 3]
      }else if(sd == FALSE){
        gsd.a = data[,3]*sqrt(data[, 4]) # SD = SE * sqrt(n.a)
      }
      n.a = data[, 4]
      N = length(dose.a)
      dose.a = dose.a/maxDose
    }
    testNLN <- NA

    mean.a = LNtoN(gmean.a,gsd.a)[1:N]
    sd.a = LNtoN(gmean.a,gsd.a)[(N+1):(2*N)]

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
      mean.a[iu]=mean(y[doses==dose.a[iu]])
      sd.a[iu]=sd(y[doses==dose.a[iu]])
      n.a[iu]=sum(doses==dose.a[iu])
    }
    dose.a = dose.a/maxDose
    # test normality
    # datind <- data.frame(x = dose.a*maxDose,
    #                      y = y)
    datind <- data.frame(x = doses,
                         y = y)
    testNLN <- NLN_test(datind)
  }

  ## Bartlett test of homoscedasticity
  # on the original scale (constant variance)
  b.test.N <- bartlett(sd.a, n.a)
  if(b.test.N[2]>=0.05){
    test.var = paste0('Distributional assumption of constant variance are met, Bartlett test p-value is ',
                      round(b.test.N[2], 4))
    message(test.var)
  }else if(b.test.N[2]<0.05){
    test.var = paste0('Distributional assumption of constant variance for the normal distribution is not met, Bartlett test p-value is ', round(b.test.N[2], 4))
    warning(test.var)
  }

  is_informative_BMD = 0
  is_informative_a = 0
  is_informative_c = 0

  if(!is.null(prior.BMD)) {is_informative_BMD = 1}


  ######################
  ### PRIORS

  ## Observed background and maximum response

  obs.min = mean.a[dose.a == 0][1]
  obs.max = mean.a[dose.a == max(dose.a)][1]

  ### Family 1

  ## Default range on background and max response

  if(mean.a[1] < mean.a[N]){

    # for a
    min.min = 0.001
    mode.min = obs.min
    max.min = 2*obs.min

    # for c
    min.max = mean.a[1]*(1.01+q)
    max.max = 2*mean.a[N]
    mode.max = obs.max

    if(flat(dose.a, mean.a, n.a, inc=T) == F & is.null(maxy)){
      mode.max = 3*mean.a[N]
      min.max  = mean.a[1]*(1.01+q)
      max.max = 2*mode.max

      warning(
        "The data do not contain information on the asymptote, and the default prior for fold change has been based on 3 times the observed maximum.
            Please provide prior input on the maximum response by specifying 'maxy', if available."
      )
    }


    ## Check appropriateness of BMR value
    if(mean.a[1]*(1+q) > mean.a[N]){
      warning('The data do not contain values corresponding to the chosen BMR, lowering the specified value of q may be necessary.')
    }



  }else if(mean.a[1] > mean.a[N]){
    min.min = 0.5*obs.min
    mode.min = obs.min
    max.min = 2*obs.min

    min.max = 0.5*obs.max
    max.max = obs.min*(1-q-0.01)
    # max.max = obs.min
    mode.max = obs.max

    if(flat(dose.a, mean.a, n.a, inc=F) == F & is.null(maxy)){
      mode.max = 0.5*obs.max
      min.max  = 0.1*obs.max
      max.max = obs.min*(1-q-0.01)
      # max.max = obs.min

      warning(
        "The data do not contain information on the asymptote, and the default prior for fold change has been based on half the observed maximum.
            Please provide prior input on the maximum response by specifying 'maxy', if available."
      )
    }


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
      shape.BMD = 4
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
    prvar.d = 0.5; prmean.d = 0.4; truncd = 10000
  }else if(prior.d == 'N05'){
    prvar.d = 0.25; prmean.d = 0.5; truncd = 10000
  }
  prmean.dQE4 = 0; prvar.dQE4 = 1; truncdQ = 10000
  prvar.s=1; prmean.s=-2*log(1.5*mean(sd.a))

  # Prior on background
  a.vec = c(min.min, mode.min, max.min)

  # Prior on mu(inf)
  c.vec = c(min.max/mode.min, mode.max/mode.min, max.max/mode.min)
  if(c.vec[1] == 0) c.vec[1] = 0.0001
  if(c.vec[2] >= c.vec[3]) c.vec[2] = c.vec[3] - 0.05

  priormu1a=c(a.vec[2],BMD.vec[2],c.vec[2],prmean.d,prmean.s)
  priormu1bQ=c(a.vec[2],BMD.vec[2],c.vec[2],prmean.dQE4,prmean.s)
  priorSigma1a=diag(c(1,1,1,prvar.d,prvar.s))
  priorSigma1bQ=diag(c(1,1,1,prvar.dQE4,prvar.s))
  priorlb1a = c(a.vec[1],BMD.vec[1],c.vec[1],0,0)
  priorub1a = c(a.vec[3],BMD.vec[3],c.vec[3],0,0)

  ## start value BMD

  if(mean.a[1] < mean.a[N]){
    bmr = q
  }else{
    bmr = -q
  }

  datf=data.frame(yy=mean.a,xx=dose.a+0.00000000000001)
  fpfit=gamlss::gamlss(yy~fp(xx),family=NO,data=datf)
  RISK=function(x) (predict(fpfit,newdata=data.frame(xx=c(x)), data = datf)
                    -predict(fpfit,newdata=data.frame(xx=c(0.00000000000001)), data = datf))/
    (predict(fpfit,newdata=data.frame(xx=c(0.00000000000001)), data = datf)) - bmr
  bmd.svh=try(uniroot(RISK, interval=c(0,1))$root,silent=T)
  bmd.sv=ifelse((mode(bmd.svh)=="numeric"),bmd.svh,0.05)

  if(!is.null(prior.BMD)){
    bmd.sv = BMD.vec[2]
  }

  is_increasing = 0; is_decreasing = 0
  if(mean.a[1] < mean.a[N]){ # increasing
    is_increasing = 1
    L = 1+q+0.01
    # L = 1.01
    U = 0
    data_type = 1
    pars3d = numeric()
    pars3i = priormu1a[3] - L
    dim(pars3d)=0
    dim(pars3i)=1
  }else if(mean.a[1] > mean.a[N]){ # decreasing
    is_decreasing = 1
    L = 0.01
    U = 1-q-0.01
    # U = 0.95
    data_type = 3
    pars3d = priormu1a[3] / U
    pars3i = numeric()
    dim(pars3d)=1
    dim(pars3i)=0
  }

  ## Data in correct format

  ret.list <- list(data = list(N=N,n=n.a,x=dose.a,m=mean.a,shift=0,s2=sd.a^2,maxD=maxDose,q=q,priormu=priormu1a,priormuQ=priormu1bQ,
                               shape1 = c(fun.alpha(a = priorlb1a[1], b = priormu1a[1], c = priorub1a[1], g = shape.a),
                                          fun.alpha(a = priorlb1a[2], b = priormu1a[2], c = priorub1a[2], g = shape.BMD),
                                          fun.alpha(a = priorlb1a[3], b = priormu1a[3], c = priorub1a[3], g = shape.c), 0, 0),
                               shape2 = c(fun.beta(a = priorlb1a[1], b = priormu1a[1], c = priorub1a[1], g = shape.a),
                                          fun.beta(a = priorlb1a[2], b = priormu1a[2], c = priorub1a[2], g = shape.BMD),
                                          fun.beta(a = priorlb1a[3], b = priormu1a[3], c = priorub1a[3], g = shape.c), 0, 0),
                               priorlb = priorlb1a, priorub = priorub1a, shape.a = shape.a, shape.c = shape.c, shape.BMD = shape.BMD,
                               priorSigma=priorSigma1a, priorSigmaQ=priorSigma1bQ, init_b = 1, data_type = data_type, L = L, U = U,
                               is_increasing = is_increasing, truncd = truncd, truncdQ = truncdQ,
                               is_decreasing = is_decreasing, is_informative_a = is_informative_a, is_informative_c = is_informative_c,
                               is_informative_BMD = is_informative_BMD),
                   # start values
                   start=list(par1=priormu1a[1],par2=bmd.sv,pars3i=pars3i,pars3d=pars3d,par4=prmean.d,par5=log(1/mean(sd.a^2))),
                   startQ=list(par1=priormu1a[1],par2=bmd.sv,pars3i=pars3i,pars3d=pars3d,par4=prmean.dQE4,par5=log(1/mean(sd.a^2))),
                   test.var = test.var,
                   test.NLN = testNLN
  )


  # data in correct format
  return(ret.list)

}


#' @rdname PREP_DATA_N
#' @export
PREP_DATA_LN <- function(data, # a dataframe with input data, order of columns should be: dose, response, sd, n
                         sumstats = TRUE, # TRUE if summary data, FALSE if individual data
                         geom.stats = FALSE, # TRUE if geometric summary data is provided
                         sd = TRUE, # TRUE if sd per dose group is given, FALSE is se is given
                         q, # the BMR
                         bkg = NULL, maxy = NULL, prior.BMD = NULL, # possible expert info on background and max response
                         shape.a = 4, shape.c = 4, shape.BMD = 0.0001, # shape for the PERT distribution
                         # prmean.d = 1, prmean.dQE4 = 0
                         prior.d = 'N11',
                         extended = TRUE
){


  if(sumstats == TRUE & geom.stats == FALSE){
    data = data[order(data[, 1]), ]
    dose.a = data[, 1]
    maxDose = max(dose.a)
    ## if dose levels not unique
    if(length(dose.a) != length(unique(dose.a))){
      dose.a = sort(unique(dose.a))
      N = length(dose.a)
      mean.a=rep(NA,N)
      sd.a=rep(NA,N)
      n.a=rep(NA,N)
      for (iu in (1:N)){
        mean.a[iu] = mean(data[,2][data[,1] == dose.a[iu]])
        sd.a[iu] = mean(data[,3][data[,1] == dose.a[iu]])
        n.a[iu] = sum(data[,4][data[,1] == dose.a[iu]])
      }
      dose.a = dose.a/maxDose

    }else{

      mean.a = data[, 2]
      if(sd == TRUE){
        sd.a = data[, 3]
      }else if(sd == FALSE){
        sd.a = data[,3]*sqrt(data[, 4]) # SD = SE * sqrt(n.a)
      }
      n.a = data[, 4]
      N = length(dose.a)
      dose.a = dose.a/maxDose
    }
    # shift if negative means occur
    shift=0
    gmean.a2 = log(NtoLN(mean.a,sd.a))[1:N]
    if (min(gmean.a2)<0) {gmean.a = gmean.a2-20*min(gmean.a2); shift = 20*min(gmean.a2)}
    if (min(gmean.a2)>=0) gmean.a = gmean.a2
    gsd.a = log(NtoLN(mean.a,sd.a))[(N+1):(2*N)]

    testNLN <- NA
  }else if(sumstats == TRUE & geom.stats == TRUE){
    data = data[order(data[, 1]), ]
    dose.a = data[, 1]
    maxDose = max(dose.a)
    mean.a = data[, 2]
    if(sd == TRUE){
      sd.a = data[, 3]
    }else if(sd == FALSE){
      sd.a = data[,3]*sqrt(data[, 4]) # SD = SE * sqrt(n.a)
    }
    gsd.a = log(sd.a)
    gmean.a = log(mean.a)
    gmean.a2 = log(mean.a)
    shift = 0
    n.a = data[, 4]
    N = length(dose.a)
    dose.a = dose.a/maxDose
    testNLN <- NA
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
    # datind <- data.frame(x = dose.a*maxDose,
    #                      y = y)
    datind <- data.frame(x = doses,
                         y = y)
    testNLN <- NLN_test(datind)
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

  obs.min = mean.a[dose.a == 0][1]
  obs.max = mean.a[dose.a == max(dose.a)][1]

  ### Family 1

  ## Default range on background and max response

  if(mean.a[1] < mean.a[N]){
    min.min = 0.001
    mode.min = obs.min
    max.min = 2*obs.min

    min.max = mean.a[1]*(1.01+q)
    max.max = 2*mean.a[N]
    mode.max = obs.max

    if(flat(dose.a, mean.a, n.a, inc=T) == F & is.null(maxy)){
      mode.max = 3*mean.a[N]
      min.max  = mean.a[1]*(1.01+q)
      max.max = 2*mode.max

      warning(
        "The data do not contain information on the asymptote, and the default prior for fold change has been based on 3 times the observed maximum.
            Please provide prior input on the maximum response by specifying 'maxy', if available."
      )
    }


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

    if(flat(dose.a, mean.a, n.a, inc=F) == F & is.null(maxy)){
      mode.max = 0.5*obs.max
      min.max  = 0.1*obs.max
      max.max = obs.min*(1-q-0.01)

      warning(
        "The data do not contain information on the asymptote, and the default prior for fold change has been has been based on half the observed maximum.
            Please provide prior input on the maximum response by specifying 'maxy', if available."
      )
    }

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
      shape.BMD = 4
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
    prvar.d = 0.5; prmean.d = 0.4; truncd = 10000
  }else if(prior.d == 'N05'){
    prvar.d = 0.25; prmean.d = 0.5; truncd = 10000
  }
  prmean.dQE4 = 0; prvar.dQE4 = 1; truncdQ = 10000
  prvar.s=1; prmean.s=-2*log(1.5*mean(gsd.a))

  # Prior on background
  a.vec = c(min.min, mode.min, max.min)

  # Prior on mu(inf)
  c.vec = c(min.max/mode.min, mode.max/mode.min, max.max/mode.min)
  if(c.vec[1] == 0) c.vec[1] = 0.0001
  if(c.vec[2] >= c.vec[3]) c.vec[2] = c.vec[3] - 0.05


  priormu1a=c(a.vec[2],BMD.vec[2],c.vec[2],prmean.d,prmean.s)
  priormu1bQ=c(a.vec[2],BMD.vec[2],c.vec[2],prmean.dQE4,prmean.s)
  priorSigma1a=diag(c(1,1,1,prvar.d,prvar.s))
  priorSigma1bQ=diag(c(1,1,1,prvar.dQE4,prvar.s))

  priorlb1a = c(a.vec[1],BMD.vec[1],c.vec[1],0,0)
  priorub1a = c(a.vec[3],BMD.vec[3],c.vec[3],0,0)

  # start value BMD

  if(mean.a[1] < mean.a[N]){
    bmr = q
  }else{
    bmr = -q
  }

  datf=data.frame(yy=mean.a,xx=dose.a+0.00000000000001)
  fpfit=gamlss::gamlss(yy~fp(xx),family=gamlss.dist::NO(),data=datf)
  RISK=function(x) (predict(fpfit,newdata=data.frame(xx=c(x)),data=datf, type = "response")
                    -predict(fpfit,newdata=data.frame(xx=c(0.00000000000001)),data=datf, type = "response"))/
    (predict(fpfit,newdata=data.frame(xx=c(0.00000000000001)),data=datf, type = "response")) - bmr
  bmd.svh=try(uniroot(RISK, interval=c(0,1))$root,silent=T)
  bmd.sv=ifelse((mode(bmd.svh)=="numeric"),bmd.svh,0.05)


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
                               priormuQ = priormu1bQ,
                               shape1 = c(fun.alpha(a = priorlb1a[1], b = priormu1a[1], c = priorub1a[1], g = shape.a),
                                          fun.alpha(a = priorlb1a[2], b = priormu1a[2], c = priorub1a[2], g = shape.BMD),
                                          fun.alpha(a = priorlb1a[3], b = priormu1a[3], c = priorub1a[3], g = shape.c), 0, 0),
                               shape2 = c(fun.beta(a = priorlb1a[1], b = priormu1a[1], c = priorub1a[1], g = shape.a),
                                          fun.beta(a = priorlb1a[2], b = priormu1a[2], c = priorub1a[2], g = shape.BMD),
                                          fun.beta(a = priorlb1a[3], b = priormu1a[3], c = priorub1a[3], g = shape.c), 0, 0),
                               priorlb = priorlb1a, priorub = priorub1a, shape.a = shape.a, shape.c = shape.c, shape.BMD = shape.BMD,
                               priorSigma=priorSigma1a, priorSigmaQ = priorSigma1bQ, truncd = truncd, truncdQ = truncdQ,
                               init_b = 1, data_type = data_type, L = L, U = U, is_increasing = is_increasing,
                               is_decreasing = is_decreasing, is_informative_a = is_informative_a, is_informative_c = is_informative_c,
                               is_informative_BMD = is_informative_BMD),
                   start = list(par1=priormu1a[1],par2=bmd.sv,pars3i=pars3i,pars3d=pars3d,par4=prmean.d,par5=log(1/mean(gsd.a^2))),
                   startQ = list(par1=priormu1a[1],par2=bmd.sv,pars3i=pars3i,pars3d=pars3d,par4=prmean.dQE4,par5=log(1/mean(gsd.a^2))),
                   test.var = test.var,
                   test.NLN = testNLN
  )

  # data in correct format
  return(ret.list)

}




#' Function to set data in the correct input format for model-averaging functions.
#' This function also generates appropriate start values and uninformative priors for each model
#'
#' @param data a dataframe with individual-level input data, order of columns should be: dose, response, litter
#' @param q specified BMR
#' @param bkg vector containing minimum, most likely (optional, can be NA), and maximum value for the background response. Defaults to NULL (non-informative prior)
#' @param maxy vector containing minimum, most likely (optional, can be NA), and maximum value for the response at dose infinity. Defaults to NULL (non-informative prior)
#' @param prior.BMD vector containing minimum, most likely (optional, can be NA), and maximum value for the BMD. Defaults to NULL (non-informative prior)
#' @param shape.a shape parameter for the modified PERT distribution on parameter a, defaults to 4 (peaked at most likely value), a value of 0.0001 implies a uniform distribution
#' @param shape.c shape parameter for the modified PERT distribution on parameter c, defaults to 4 (peaked at most likely value), a value of 0.0001 implies a uniform distribution
#' @param shape.BMD shape parameter for the modified PERT distribution on parameter BMD, defaults to 0.0001 implying a uniform distribution. Can be set to 4 in case of informative prior
#' @param prior.d prior distribution for parameter d, should be either N11 (default), EPA or N05 (for a N(0.5,0.5) prior)
#' @param extended logical indicating whether the dose range should be extended to maxDose^2 (default is TRUE)
#'
#' @description The function takes a dataset as input and generates the data list and starting values needed by the stan
#'              models to be fitted.
#'
#' @return List with data and start values in correct format to be directly used within the model-averaging functions.
#'
#' @export PREP_DATA_N_C
#'
PREP_DATA_N_C <- function(data, # a dataframe with input data, order of columns should be: dose, response, litter
                          q, # the BMR
                          bkg = NULL,
                          maxy = NULL,
                          prior.BMD = NULL, # possible expert info on background and max response
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

  datind <- data.frame(x = indiv.data$dose,
                       y = indiv.data$response)
  # testNLN <- NLN_test(datind)

  ## Observed background and maximum response

  obs.min = mean(indiv.data$response[indiv.data$dose==0])
  obs.max = mean(indiv.data$response[indiv.data$dose==maxDose])

  is_informative_BMD = 0
  is_informative_a = 0
  is_informative_c = 0

  ## Overall mean for test of flatness
  means.all <- indiv.data %>%
    group_by(dose) %>%
    summarise(mresp = mean(response))
  dose.a = unique(indiv.data$dose)
  mean.a = c()
  for(m in 1:length(dose.a)){
    mean.a[m] <- means.all$mresp[means.all$dose == dose.a[m]]
  }

  if(!is.null(prior.BMD)) {is_informative_BMD = 1}

  ######################
  ### PRIORS

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

    if(flatC(dose.a, mean.a,inc=T) == F & is.null(maxy)){
      mode.max = 3*obs.max
      min.max  = obs.min*(1.01+q)
      max.max = 2*mode.max

      warning(
        "The data do not contain information on the asymptote, and the default prior for fold change has been based on 3 times the observed maximum.
            Please provide prior input on the maximum response by specifying 'maxy', if available."
      )
    }

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

    if(flatC(dose.a, mean.a,inc=F) == F & is.null(maxy)){
      mode.max = 0.5*obs.max
      min.max  = 0.1*obs.max
      max.max = obs.min*(1-q-0.01)
      # max.max = obs.min

      warning(
        "The data do not contain information on the asymptote, and the default prior for fold change has been based on half the observed maximum.
            Please provide prior input on the maximum response by specifying 'maxy', if available."
      )
    }

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
      shape.BMD = 4
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
    # prvar.d = 1; prmean.d = 1; truncd = 10000
  }else if(prior.d == 'N05'){
    prvar.d = 0.25; prmean.d = 0.5; truncd = 10000
  }
  # prvar.d=sqrt(0.5); prmean.d = prmean.d
  prmean.dQE4 = 0; prvar.dQE4 = 1; truncdQ = 10000

  prvar.s=1; prmean.s=-2*log(1.5*sd(y[y!=0]))

  # Prior on background
  a.vec = c(min.min, mode.min, max.min)

  # Prior on mu(inf)
  c.vec = c(min.max/mode.min, mode.max/mode.min, max.max/mode.min)
  if(c.vec[1] == 0) c.vec[1] = 0.0001
  if(c.vec[2] >= c.vec[3]) c.vec[2] = c.vec[3] - 0.05

  priormu1a=c(a.vec[2],BMD.vec[2],c.vec[2],prmean.d,prmean.s,0.5)
  priormu1bQ=c(a.vec[2],BMD.vec[2],c.vec[2],prmean.dQE4,prmean.s,0.5)
  priorSigma1a=diag(c(1,1,1,prvar.d,prvar.s))
  priorSigma1bQ=diag(c(1,1,1,prvar.dQE4,prvar.s))
  priorlb1a = c(a.vec[1],BMD.vec[1],c.vec[1],0,0,0)
  priorub1a = c(a.vec[3],BMD.vec[3],c.vec[3],0,0,1)

  # start value BMD
  if(obs.min < obs.max){
    bmr = q
  }else{
    bmr = -q
  }

  # start value BMD
  datf = data.frame(yy=data[,2], xx=data[,1]/max(data[,1])+0.00000000000001)
  fpfit = gamlss(yy~fp(xx),family=NO(),data=datf)
  RISK = function(x) (predict(fpfit,newdata=data.frame(xx=c(x)),data=datf)-predict(fpfit,newdata=data.frame(xx=c(0.00000000000001)),data=datf))/
    (predict(fpfit,newdata=data.frame(xx=c(0.00000000000001)),data=datf)) - bmr
  bmd.svh = try(uniroot(RISK, interval=c(0, 1))$root,silent=T)
  bmd.sv=ifelse((mode(bmd.svh)=="numeric"),bmd.svh,0.05)

  if(!is.null(prior.BMD)) bmd.sv = BMD.vec[2]

  is_increasing = 0; is_decreasing = 0
  if(obs.min < obs.max){ # increasing
    is_increasing = 1
    L = 1+q+0.01
    # L = 1.01
    U = 0
    data_type = 1
    pars3d = numeric()
    pars3i = priormu1a[3] - L
    dim(pars3d)=0
    dim(pars3i)=1
  }else if(obs.min > obs.max){ # decreasing
    is_decreasing = 1
    L = 0.01
    U = 1-q-0.01
    # U = 0.95
    data_type = 3
    pars3d = priormu1a[3] / U
    pars3i = numeric()
    dim(pars3d)=1
    dim(pars3i)=0
  }

  ## Test for homoscedasticity and normality of residuals, based on saturated ANOVA model

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

  data.modstanSM = list(N=N, n=n, nc=nc, maxN=maxN, maxNc=maxNc,
                        nij=nij, y=y, q=q, shift=0,
                        priormu=priorSM$priormu, priorSigma=priorSM$priorSigma,
                        priorlb=priorSM$priorlb, priorub=priorSM$priorub,
                        priorg=4, data_type=data_type
  )

  optSM = rstan::optimizing(stanmodels$mSMc, data = data.modstanSM,
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
    test.var = paste0('Distributional assumption of constant variance for the normal distribution is met, Bartlett test p-value is ', round(b.test.N$p.value, 4))
    message(test.var)
  }else if(b.test.N$p.value<0.05){
    test.var = paste0('Distributional assumption of constant variance for the normal distribution is not met, Bartlett test p-value is ', round(b.test.N$p.value, 4))
    warning(test.var)
  }

  # normality
  norm.test.N <- shapiro.test(indiv.data$res)
  # b.test.N <- bartlett(sd.a, n.a)
  if(norm.test.N$p.value>=0.05){
    test.varN = paste0('Distributional assumption of normality of residuals for the normal distribution is met, Shapiro test p-value is ',
                       round(norm.test.N$p.value, 4))
    message(test.varN)
  }else if(norm.test.N$p.value<0.05){
    test.varN = paste0('Distributional assumption of normality of residuals for the normal distribution is not met, Shapiro test p-value is ',
                       round(norm.test.N$p.value, 4))
    warning(test.varN)
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
                               y=y,
                               shift=0,
                               data = indiv.data,
                               # s2=sd.a^2,
                               maxD=maxDose,q=q,priormu=priormu1a,priormuQ = priormu1bQ,
                               shape1 = c(fun.alpha(a = priorlb1a[1], b = priormu1a[1], c = priorub1a[1], g = shape.a),
                                          fun.alpha(a = priorlb1a[2], b = priormu1a[2], c = priorub1a[2], g = shape.BMD),
                                          fun.alpha(a = priorlb1a[3], b = priormu1a[3], c = priorub1a[3], g = shape.c), 0, 0,
                                          fun.alpha(a = priorlb1a[6], b = priormu1a[6], c = priorub1a[6], g = 0.0001)),
                               shape2 = c(fun.beta(a = priorlb1a[1], b = priormu1a[1], c = priorub1a[1], g = shape.a),
                                          fun.beta(a = priorlb1a[2], b = priormu1a[2], c = priorub1a[2], g = shape.BMD),
                                          fun.beta(a = priorlb1a[3], b = priormu1a[3], c = priorub1a[3], g = shape.c), 0, 0,
                                          fun.beta(a = priorlb1a[6], b = priormu1a[6], c = priorub1a[6], g = 0.0001)),
                               priorlb = priorlb1a, priorub = priorub1a, shape.a = shape.a, shape.c = shape.c, shape.BMD = shape.BMD,
                               priorSigma=priorSigma1a, priorSigmaQ=priorSigma1bQ, init_b = 1, data_type = data_type, L = L, U = U, is_increasing = is_increasing,
                               is_decreasing = is_decreasing, truncd = truncd, truncdQ = truncdQ,
                               is_informative_a = is_informative_a, is_informative_c = is_informative_c,
                               is_informative_BMD = is_informative_BMD),
                   # start values
                   start=list(par1=priormu1a[1],par2=bmd.sv,pars3i=pars3i,pars3d=pars3d,par4=prmean.d,par5=log(1/var(y[y!=0])),
                              par6=0.5),
                   startQ=list(par1=priormu1a[1],par2=bmd.sv,pars3i=pars3i,pars3d=pars3d,par4=prmean.dQE4,par5=log(1/var(y[y!=0])),
                               par6=0.5),
                   # test.var = test.var,
                   # test.NLN = testNLN
                   shapiro.p = norm.test.N$p.value,
                   bartlett.p = b.test.N$p.value,
                   shapiro.msg = test.varN,
                   bartlett.msg = test.var
  )

  # test for dose-response effect
  # DR.effect = anydoseresponseNI(dose.a,mean.a,sd.a,n.a)

  # data in correct format
  return(ret.list)

}


#' @rdname PREP_DATA_N_C
#' @export
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

  datind <- data.frame(x = indiv.data$dose,
                       y = indiv.data$response)
  # testNLN <- NLN_test(datind)

  is_informative_BMD = 0
  is_informative_a = 0
  is_informative_c = 0

  if(!is.null(prior.BMD)) {is_informative_BMD = 1}

  ## Overall mean for test of flatness
  means.all <- indiv.data %>%
    group_by(dose) %>%
    summarise(mresp = mean(log(response)))
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

    if(flatC(dose.a, mean.a,inc=T) == F & is.null(maxy)){
      mode.max = 3*obs.max
      min.max  = obs.min*(1.01+q)
      max.max = 2*mode.max

      warning(
        "The data do not contain information on the asymptote, and the default prior for fold change has been based on 3 times the observed maximum.
            Please provide prior input on the maximum response by specifying 'maxy', if available."
      )
    }


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

    if(flatC(dose.a, mean.a,inc=F) == F & is.null(maxy)){
      mode.max = 0.5*obs.max
      min.max  = 0.1*obs.max
      max.max = obs.min*(1-q-0.01)
      # max.max = obs.min

      warning(
        "The data do not contain information on the asymptote, and the default prior for fold change has been based on half the observed maximum.
            Please provide prior input on the maximum response by specifying 'maxy', if available."
      )
    }


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
      shape.BMD = 4
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
  }else if(prior.d == 'N05'){
    prvar.d = 0.25; prmean.d = 0.5; truncd = 10000
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

  # start value BMD
  if(obs.min < obs.max){
    bmr = q
  }else{
    bmr = -q
  }

  datf = data.frame(yy=data$response,xx=(data$dose/max(data$dose))+0.00000000000001)
  fpfit=gamlss::gamlss(yy~fp(xx),family=NO(),data=datf)
  RISK=function(x) (predict(fpfit,newdata=data.frame(xx=c(x)),data=datf)-predict(fpfit,newdata=data.frame(xx=c(0.00000000000001)),data=datf))/
    (predict(fpfit,newdata=data.frame(xx=c(0.00000000000001)),data=datf)) - bmr
  bmd.svh=try(uniroot(RISK, interval=c(0, 1))$root,silent=T)
  # bmd.sv=ifelse((mode(bmd.svh)=="numeric"),exp(bmd.svh),0.05)
  bmd.sv=ifelse((mode(bmd.svh)=="numeric"),bmd.svh,0.05)

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

  optSM = rstan::optimizing(stanmodels$mSMc, data = data.modstanSM,
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
    test.varN = paste0('Distributional assumption of normality of residuals for the lognormal distribution is met, Shapiro test p-value is ',
                       round(norm.test.N$p.value, 4))
    message(test.varN)
  }else if(norm.test.N$p.value<0.05){
    test.varN = paste0('Distributional assumption of normality of residuals for the lognormal distribution is not met, Shapiro test p-value is ',
                       round(norm.test.N$p.value, 4))
    warning(test.varN)
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
                   # test.var = test.var,
                   # test.NLN = testNLN
                   shapiro.p = norm.test.N$p.value,
                   bartlett.p = b.test.N$p.value,
                   shapiro.msg = test.varN,
                   bartlett.msg = test.var
  )




  # test for dose-response effect
  # DR.effect = anydoseresponseNI(dose.a,mean.a,sd.a,n.a)

  # data in correct format
  return(ret.list)

}



#' Function to set data in the correct format
#'
#' This function also generates appropriate start values and uninformative priors for each model
#' Input should be given as arithmetic mean and standard deviation on the original scale
#'
#' @param data a dataframe with input data, order of columns should be: dose, number of adverse events, n
#' @param sumstats logical indicating whether summary (T, default) or individual-level (F) data is provided. If individual-level data are provided, a litter indicator should be included instead of n (column 3)
#' @param q specified BMR
#' @param bkg vector containing minimum, most likely (optional), and maximum value for the background response. Defaults to NULL (non-informative prior)
#' @param prior.BMD vector containing minimum, most likely (optional), and maximum value for the BMD. Defaults to NULL (non-informative prior)
#' @param shape.a shape parameter for the modified PERT distribution on parameter a, defaults to 4 (peaked at most likely value), a value of 0.0001 implies a uniform distribution
#' @param shape.BMD shape parameter for the modified PERT distribution on parameter BMD, defaults to 0.0001 implying a uniform distribution. Can be set to 4 in case of informative prior
#' @param cluster logical variable to indicate if data is clustered. TRUE = clustered data. Defaults to FALSE
#' @param prior.d prior distribution for parameter d, should be either N11 (default), EPA or N05 (for a N(0.5,0.5) prior)
#' @param extended logical indicating whether the dose range should be extended to maxDose^2 (default is TRUE)
#'
#' @return List with data and start values in correct format to be directly used within the BMA functions.
#'
#' @export PREP_DATA_QA
#'
PREP_DATA_QA <- function(data, # a dataframe with input data, order of columns should be: dose, response, n
                         sumstats = TRUE, # TRUE if summary data, FALSE if individual data
                         q, # the BMR,
                         bkg = NULL, # possible expert info on background response (a),
                         prior.BMD = NULL, # possible expert prior on the BMD,
                         shape.a = 4, #scale parameter for Pert priors for a
                         shape.BMD = 0.0001, #scale parameter for the Pert priors for BMD
                         cluster = FALSE, # indicate if data is clustered
                         prior.d = 'N11',
                         extended = TRUE
){

  if(sumstats == TRUE){
    data = data[order(data[, 1]), ]
    dose.a = data[, 1]
    maxDose = max(dose.a)
    y.a = data[, 2]
    n.a = data[, 3]
    N = length(dose.a)
    dose.a = dose.a/maxDose
  }else if(sumstats == FALSE){ # summarize by dose and cluster ! (NOT TESTED)
    doses = data[, 1]
    maxDose = max(doses)
    # dose.a = sort(unique(doses))
    litter = data[,3]
    # N = length(dose.a)
    y.a = rep(NA,N)
    n.a = rep(NA,N)
    ybin = data[, 2]
    id = 1
    for(iu in unique(litter)){
      y.a[id] = sum(ybin[litter = iu])
      n.a[id] = sum(litter == iu)
      dose.a[id] = doses[litter == iu]
      id = id + 1
    }
  }

  datf = data.frame(yy = y.a, n.a = n.a, xx = dose.a)
  if(cluster == FALSE) {
    fpfit = gamlss(cbind(yy, n.a-yy)~fp(xx),family=BI(mu.link = 'logit'),data=datf)

  } else if(cluster == TRUE) {
    fpfit = gamlss(cbind(yy, n.a-yy)~fp(xx),family=BB,
                   sigma.formula = ~1,data=datf)
    fpfit2 <- try(gamlss(cbind(yy,n.a-yy)~as.factor(xx), sigma.formula=~1, family=BB, data=datf),
                  silent = TRUE)
    rhohat <- exp(fpfit2$sigma.coefficients)/(exp(fpfit2$sigma.coefficients)+1)

  } else stop('provide cluster to be TRUE or FALSE')

  RISK = function(x) (predict(fpfit,newdata=data.frame(xx=c(exp(x))), data=datf, type = "response")-
                        predict(fpfit,newdata=data.frame(xx=c(0.00000000000001)), data=datf, type = "response"))/
    (1 - (predict(fpfit,newdata=data.frame(xx=c(0.00000000000001)), data=datf, type = "response"))) - q
  bmd.svh = try(uniroot(RISK, interval=c(-5,0))$root,silent=TRUE)
  bmd.sv <- ifelse((mode(bmd.svh)=="numeric"),(exp(bmd.svh)+0.5)/2,0.05)

  ## clustered data option
  is_bin <- ifelse(cluster == FALSE, 1, 0)
  is_betabin <- ifelse(cluster == TRUE, 1, 0)

  ######################
  ### PRIORS

  ## Observed background and maximum response

  #obs.min = ifelse(y.a[1]/n.a[1] == 0, 1.0e-03, y.a[1]/n.a[1])

  is_informative_BMD = 0
  is_informative_a = 0


  ## Default range on background
  mindose.a <- which(dose.a == min(dose.a))
  miny.a <- sum(y.a[mindose.a])
  minn.a <- sum(n.a[mindose.a])
  a.min <- ifelse(miny.a != 0, max(c(prop.test(miny.a, minn.a)$conf.int[1]/2, 1/(10*minn.a))),
                  .Machine$double.xmin)
  a.max <- min(c(3*prop.test(miny.a, minn.a)$conf.int[2]/2, 1 - 1/(10*minn.a)))
  a.mode <-  max(c(miny.a/minn.a, 1/(5*minn.a)))

  ## If info on background is given
  if(!is.null(bkg)){

    # is_informative_a = 1

    if(!is.na(bkg[2])){
      a.mode = bkg[2]
    }else{
      a.mode = bkg[1] + ((bkg[3]-bkg[1])/2)
    }
    if(!is.na(bkg[1])){
      a.min = bkg[1]
    }
    if(!is.na(bkg[3])){
      a.max = bkg[3]
    }

    is_informative_a = 1

  }else {
    message("Default prior choices used on background")
  }

  # Prior on a
  a.vec <- c(a.min, a.mode, a.max)

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
      shape.BMD = 4
    }else{
      BMD.mode = (prior.BMD[1]/maxDose) + (((prior.BMD[3]/maxDose) - prior.BMD[1]/maxDose)/2)
    }
    if(!is.na(prior.BMD[1])){
      BMD.min = prior.BMD[1]/maxDose
    }
    if(!is.na(prior.BMD[3])){
      BMD.max = prior.BMD[3]/maxDose
    }

    is_informative_BMD = 1

  }else {
    message("Default prior choices used on BMD")
  }

  BMD.vec <- c(BMD.min, BMD.mode, BMD.max)

  # Default (normal) priors on k, d, Pert on a
  if(prior.d == 'N11'){
    prvar.d = 1; prmean.d = 1; truncd = 5
  }else if(prior.d == 'EPA'){
    prvar.d = 0.5; prmean.d = 0.4; truncd = 10000
  }else if(prior.d == 'N05'){
    prvar.d = 0.25; prmean.d = 0.5; truncd = 10000
  }
  prmean.dQE4 = 0; prvar.dQE4 = 1; truncdQ = 10000

  #prvar.k=1; prmean.k=1
  # family 1a
  priormu1a <- c(a.vec[2], BMD.vec[2], prmean.d, ifelse(is_betabin==1, rhohat, 0))
  priormu1bQ <- c(a.vec[2], BMD.vec[2], prmean.dQE4, ifelse(is_betabin==1, rhohat, 0))
  priorSigma1a <- diag(c(1, 1, prvar.d))
  priorSigma1bQ <- diag(c(1, 1, prvar.dQE4))

  priorlb1a <- c(a.vec[1], BMD.vec[1])
  priorub1a <- c(a.vec[3], BMD.vec[3])

  if(!is.null(prior.BMD)){
    bmd.sv = BMD.vec[2]
  }

  if(is_bin==1) {
    start = list(par1 = priormu1a[1], par2 = bmd.sv, par3 = priormu1a[3])
    startQ = list(par1 = priormu1a[1], par2 = bmd.sv, par3 = priormu1bQ[3])

  } else {
    rho <- rhohat; dim(rho) <- 1
    start = list(par1 = priormu1a[1], par2 = bmd.sv, par3 = priormu1a[3],
                 rho = rho )
    startQ = list(par1 = priormu1a[1], par2 = bmd.sv, par3 = priormu1bQ[3],
                  rho = rho )
  }

  return(list(
    # data and priors
    data = list(N = N, n = n.a, x = dose.a, y = y.a, yint = y.a, nint = n.a, maxD = maxDose,
                q = q, priormu = priormu1a, priormuQ = priormu1bQ, priorSigmaQ=priorSigma1bQ,
                truncd = truncd, truncdQ = truncdQ,
                priorlb = priorlb1a, priorub = priorub1a, priorSigma = priorSigma1a,
                eps = 1.0E-06, priorgama = c(shape.a, shape.BMD),
                init_b = 1, is_informative_a = is_informative_a, is_informative_BMD = is_informative_BMD,
                is_bin = is_bin, is_betabin = is_betabin),
    # start values
    start = start,
    startQ = startQ
  ))

}


#' Function to set data in the correct input format for model-averaging functions.
#' This function also generates appropriate start values and uninformative priors for each model
#'
#' @param data a dataframe with input data, order of columns should be: dose, response, sd (or se), n, covariate (if given)
#' @param sumstats logical indicating whether summary (T, default) or individual-level (F) data is provided
#' @param geom.stats logicial indicating whether, if summary data are provided, these are geometric (T) or arithmetic (F, default) summary statistics
#' @param sd logical indicating whether standard deviation (T, default) or standard error (F) is provided
#' @param q specified BMR
#' @param bkg vector containing minimum, most likely (optional, can be NA), and maximum value for the background response. Defaults to NULL (non-informative prior)
#' @param maxy vector containing minimum, most likely (optional, can be NA), and maximum value for the response at dose infinity. Defaults to NULL (non-informative prior)
#' @param prior.BMD vector containing minimum, most likely (optional, can be NA), and maximum value for the BMD. Defaults to NULL (non-informative prior)
#' @param shape.a shape parameter for the modified PERT distribution on parameter a, defaults to 4 (peaked at most likely value), a value of 0.0001 implies a uniform distribution
#' @param shape.c shape parameter for the modified PERT distribution on parameter c, defaults to 4 (peaked at most likely value), a value of 0.0001 implies a uniform distribution
#' @param shape.BMD shape parameter for the modified PERT distribution on parameter BMD, defaults to 0.0001 implying a uniform distribution. Can be set to 4 in case of informative prior
#' @param prior.d prior distribution for parameter d, should be either N11 (default), EPA or N05 (for a N(0.5,0.5) prior)
#' @param extended logical indicating whether the dose range should be extended to maxDose^2 (default is TRUE)
#' @param covariate on which parameters a covariate effect should be used. Defaults to 'all', other options are 'a_sigma2' and 'BMD_d'
#'
#' @description The function takes a dataset as input and generates the data list and starting values needed by the stan
#'              models to be fitted.

#' @return List with data and start values in correct format to be directly used within the model-averaging functions.
#'
#' @export PREP_DATA_NCOV
#'
PREP_DATA_NCOV <- function(data, # a dataframe with input data, order of columns should be: dose, response, sd, n, covar
                           sumstats = TRUE, # TRUE if summary data, FALSE if individual data
                           geom.stats = FALSE, # TRUE if geometric summary data is provided
                           sd = TRUE, # TRUE if sd per dose group is given, FALSE is se is given
                           q, # the BMR
                           bkg = NULL,
                           maxy = NULL,
                           prior.BMD = NULL, # possible expert info on background and max response
                           shape.a = 4, shape.c = 4, shape.BMD = 0.0001, # shape for the PERT distribution,
                           prior.d = 'N11',
                           extended = TRUE,
                           covariate = 'all' # OPTIONS: 'a_sigma2', 'BMD_d', 'all'; for 'none', use original prep_data
){

  if(sumstats == TRUE & geom.stats == FALSE){
    data = as.data.frame(data[order(data[, 1]), ])
    dose.a = data[, 1]
    maxDose = max(dose.a)
    mean.a = data[, 2]
    if(sd == TRUE){
      sd.a = data[, 3]
    }else if(sd == FALSE){
      sd.a = data[,3]*sqrt(data[, 4]) # SD = SE * sqrt(n.a)
    }
    n.a = data[, 4]
    # N = length(unique(dose.a))
    N = length(dose.a)
    dose.a = dose.a/maxDose
    covar = data[,5]

  }else if(sumstats == TRUE & geom.stats == TRUE){

    data = as.data.frame(data[order(data[, 1]), ])
    dose.a = data[, 1]
    maxDose = max(dose.a)
    gmean.a = data[, 2]
    if(sd == TRUE){
      gsd.a = data[, 3]
    }else if(sd == FALSE){
      gsd.a = data[,3]*sqrt(data[, 4]) # SD = SE * sqrt(n.a)
    }
    n.a = data[, 4]
    # N = length(unique(dose.a))
    N = length(dose.a)
    dose.a = dose.a/maxDose
    covar = data[,5]

    mean.a = LNtoN(gmean.a,gsd.a)[1:N]
    sd.a = LNtoN(gmean.a,gsd.a)[(N+1):(2*N)]

  }else if(sumstats == FALSE){
    # still to be worked on
    # data = data[order(data[, 1]), ]
    # doses = data[, 1]
    # maxDose = max(doses)
    # dose.a=sort(unique(doses))
    # N=length(dose.a)
    # mean.a=rep(NA,N)
    # sd.a=rep(NA,N)
    # n.a=rep(NA,N)
    # y = data[, 2]
    # for (iu in (1:N)){
    #   mean.a[iu]=mean(y[doses==dose.a[iu]])
    #   sd.a[iu]=sd(y[doses==dose.a[iu]])
    #   n.a[iu]=sum(doses==dose.a[iu])
    # }
    # dose.a = dose.a/maxDose
    # covar = data[,5]

    stop('Data should be provided as summary statistics.')
  }

  covar_lvls <- unique(covar)
  nlevels <- length(unique(covar))
  ## Bartlett test of homoscedasticity
  # on the original scale (constant variance)
  #b.test.N <- numeric(nlevels)

  # if(covariate == 'BMD_d' | covariate == 'none'){
  ## get overall data for values of a and sigma
  dose.a2 = sort(unique(dose.a))
  N2 = length(dose.a2)
  mean.a2 = rep(NA, N2)
  sd.a2 = rep(NA, N2)
  n.a2 = rep(NA, N2)
  for(iu in (1:N2)){
    mean.a2[iu] = mean(mean.a[dose.a == dose.a2[iu]])
    sd.a2[iu] = mean(sd.a[dose.a == dose.a2[iu]])
    n.a2[iu] = sum(n.a[dose.a == dose.a2[iu]])
  }
  # }

  if(covariate == 'a_sigma2' | covariate == 'all') {
    test.var <- character(nlevels)
    prmean.s <- par5 <- numeric(nlevels)
    for(i in 1:nlevels){
      b.test.N <- bartlett(sd.a[covar == covar_lvls[i]], n.a[covar == covar_lvls[i]])
      if(b.test.N[2]>=0.05){
        test.var[i] = paste0('Distributional assumption of constant variance are met for group ', covar_lvls[i],
                             ' Bartlett test p-value is ',
                             round(b.test.N[2], 4))
        message(test.var[i])
      }else if(b.test.N[2]<0.05){
        test.var[i] = paste0('Distributional assumption of constant variance for the normal distribution is not met for group ',
                             covar_lvls[i], ' Bartlett test p-value is ', round(b.test.N[2], 4))
        warning(test.var[i])
      }
      prmean.s[i] =-2*log(1.5*mean(sd.a[covar == covar_lvls[i]]))
      par5[i] <- log(1/mean(sd.a[covar == covar_lvls[i]]^2))
    }

  } else {



    b.test.N <- bartlett(sd.a2, n.a2)

    if(b.test.N[2]>=0.05){
      test.var = paste0('Distributional assumption of constant variance are met for ',
                        'Bartlett test p-value is ',
                        round(b.test.N[2], 4))
      message(test.var)
    }else if(b.test.N[2]<0.05){
      test.var = paste0('Distributional assumption of constant variance for the normal distribution is not met for ',
                        ' Bartlett test p-value is ', round(b.test.N[2], 4))
      warning(test.var)
    }

    prmean.s = -2*log(1.5*mean(sd.a2))
    par5 <- log(1/mean(sd.a2^2))
    # dim(par5) <- 1
  }

  nlevels_sigma <- ifelse(covariate == 'a_sigma2' | covariate == 'all', nlevels, 1)
  dim(par5) <- nlevels_sigma

  is_informative_BMD = 0
  is_informative_a = 0
  is_informative_c = 0

  if(!is.null(prior.BMD)) {is_informative_BMD = 1}


  ######################
  ### PRIORS

  ## Observed background and maximum response
  ### Family 1
  if(covariate == 'a_sigma2' | covariate == 'all') {

    obs.min <- obs.max <- min.min <- mode.min <- max.min <- min.max <- max.max <- mode.max <- numeric(nlevels)

    for(i in 1:nlevels) {
      obs.min[i] = mean.a[dose.a == 0 & covar == covar_lvls[i]][1]
      obs.max[i] = mean.a[dose.a == max(dose.a) & covar == covar_lvls[i]][1]

      # obs.min[i] = min(mean.a[covar == covar_lvls[i]])
      # obs.max[i] = max(mean.a[covar == covar_lvls[i]])

      if(obs.min[i] < obs.max[i]) {

        # for a
        min.min = rep(0.001, nlevels)
        mode.min[i] = obs.min[i]
        max.min[i] = 2*obs.min[i]

        # for c
        min.max[i] = obs.min[i]*(1.01+q)
        max.max[i] = 2*obs.max[i]
        mode.max[i] = obs.max[i]

        if(flat(dose.a[covar == covar_lvls[i]], mean.a[covar == covar_lvls[i]],
                n.a[covar == covar_lvls[i]], inc=T) == F & is.null(maxy)){
          mode.max[i] = 3*obs.max[i]
          min.max[i]  = obs.min[i]*(1.01+q)
          max.max[i] = 2*mode.max[i]

        }


        ## Check appropriateness of BMR value
        if(obs.min[i]*(1+q) > obs.max[i]){
          warning(paste0("The data do not contain values corresponding to the chosen BMR for group ", i,
                         " lowering the specified value of q may be necessary."))
        }

      } else if(obs.min[i] > obs.max[i]){

        min.min[i] = 0.5*obs.min[i]
        mode.min[i] = obs.min[i]
        max.min[i] = 2*obs.min[i]

        min.max[i] = 0.5*obs.max[i]
        max.max[i] = obs.min[i]*(1-q-0.01)
        # max.max = obs.min
        mode.max[i] = obs.max[i]

        if(flat(dose.a, mean.a, n.a, inc=F) == F & is.null(maxy)){
          mode.max[i] = 0.5*obs.max[i]
          min.max[i]  = 0.1*obs.max[i]
          max.max[i] = obs.min[i]*(1-q-0.01)
        }

        ## Check appropriateness of BMR value
        if(obs.min[i]*(1-q) < obs.max[i]){
          warning(paste0("The data do not contain values corresponding to the chosen BMR for group ", i,
                         " lowering the specified value of q may be necessary."))
        }

      }

      ## If info on background is given
      if(!is.null(bkg)){

        is_informative_a = 1

        if(!is.na(bkg[2])){
          mode.min[i] = bkg[2]
        }else{
          mode.min[i] = bkg[1] + ((bkg[3]-bkg[1])/2)
        }
        if(!is.na(bkg[1])){
          min.min[i] = bkg[1]
        }
        if(!is.na(bkg[3])){
          max.min[i] = bkg[3]
        }

      }else{
        message("Default prior choices used on background")
      }

      ## If info on max response is given
      if(!is.null(maxy)){

        is_informative_c = 1

        if(!is.na(maxy[2])){
          mode.max[i] = maxy[2]
        }else{
          mode.max[i] = maxy[1] + ((maxy[3]-maxy[1])/2)
        }
        if(!is.na(maxy[1])){
          min.max[i] = maxy[1]
        }
        if(!is.na(maxy[3])){
          max.max[i] = maxy[3]
        }

      }else{
        message("Default prior choices used on fold change")
      }

    }

    # Prior on background
    a.vec = matrix(c(min.min, mode.min, max.min), nrow = nlevels, ncol = 3)
    # Prior on mu(inf)
    c.vec = matrix(c(min.max/mode.min, mode.max/mode.min, max.max/mode.min), nrow = nlevels, ncol = 3)

    shape.a1 <- shape.a2 <- shape.c1 <- shape.c2 <- numeric(nlevels)

    for(i in 1:nlevels){
      if(c.vec[i,1] == 0) c.vec[i,1] = 0.0001
      if(c.vec[i,2] >= c.vec[i,3]) c.vec[i,2] = c.vec[i,3] - 0.05

      shape.a1[i] <- fun.alpha(a = a.vec[i,1], b = a.vec[i,2],
                               c = a.vec[i,3], g = shape.a)
      shape.c1[i] <- fun.alpha(a = c.vec[i,1], b = c.vec[i,2],
                               c = c.vec[i,3], g = shape.c)

      shape.a2[i] <- fun.beta(a = a.vec[i,1], b = a.vec[i,2],
                              c = a.vec[i,3], g = shape.a)
      shape.c2[i] <- fun.beta(a = c.vec[i,1], b = c.vec[i,2],
                              c = c.vec[i,3], g = shape.c)
    }


  } else {

    obs.min = mean.a2[dose.a2 == 0][1]
    obs.max = mean.a2[dose.a2 == max(dose.a)][1]

    if(obs.min < obs.max) {

      # for a
      min.min = 0.001
      mode.min = obs.min
      max.min = 2*obs.min

      # for c
      min.max = obs.min*(1.01+q)
      max.max = 2*obs.max
      mode.max = obs.max

      if(flat(dose.a2, mean.a2,
              n.a2, inc=T) == F & is.null(maxy)){
        mode.max = 3*obs.max
        min.max  = obs.min*(1.01+q)
        max.max = 2*mode.max

      }

      ## Check appropriateness of BMR value
      if(obs.min*(1+q) > obs.max){
        warning('The data do not contain values corresponding to the chosen BMR,
                lowering the specified value of q may be necessary.')
      }

    } else if(obs.min > obs.max){

      min.min = 0.5*obs.min
      mode.min = obs.min
      max.min = 2*obs.min

      min.max = 0.5*obs.max
      max.max = obs.min*(1-q-0.01)
      # max.max = obs.min
      mode.max = obs.max

      if(flat(dose.a2, mean.a2, n.a2, inc=F) == F & is.null(maxy)){
        mode.max = 0.5*obs.max
        min.max = 0.1*obs.max
        max.max = obs.min*(1-q-0.01)
      }


      ## Check appropriateness of BMR value
      if(obs.min*(1-q) < obs.max){
        warning('The data do not contain values corresponding to the chosen BMR,
                lowering the specified value of q may be necessary.')
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


    a.vec = c(min.min, mode.min, max.min)
    # Prior on mu(inf)
    c.vec = c(min.max/mode.min, mode.max/mode.min, max.max/mode.min)

    if(c.vec[1] == 0) c.vec[1] = 0.0001
    if(c.vec[2] >= c.vec[3]) c.vec[2] = c.vec[3] - 0.05

    shape.a1 <- fun.alpha(a = a.vec[1], b = a.vec[2],
                          c = a.vec[3], g = shape.a)
    shape.a2 <- fun.beta(a = a.vec[1], b = a.vec[2],
                         c = a.vec[3], g = shape.a)

    shape.c1 <- fun.alpha(a = c.vec[1], b = c.vec[2],
                          c = c.vec[3], g = shape.c)
    shape.c2 <- fun.beta(a = c.vec[1], b = c.vec[2],
                         c = c.vec[3], g = shape.c)
    # dim(a.vec) <- dim(c.vec) <- 1
  }

  nlevels_a <- nlevels_c <- ifelse(covariate == 'a_sigma2' | covariate == 'all', nlevels, 1)
  # dim(a.vec) <- nlevels_a
  # dim(c.vec) <- nlevels_c


  if(covariate == 'BMD_d' | covariate == 'all') {

    bmd.sv <- numeric(nlevels)

    for(i in 1:nlevels){
      # start value BMD
      datf=data.frame(yy=mean.a[covar == covar_lvls[i]],
                      xx=dose.a[covar == covar_lvls[i]]+0.00000000000001)
      fpfit=gamlss(yy~fp(xx),family=NO(),data=datf)
      RISK=function(x) (predict(fpfit,newdata=data.frame(xx=c((x))),data=datf)-
                          predict(fpfit,newdata=data.frame(xx=c(0.00000000000001)),data=datf))/
        (predict(fpfit,newdata=data.frame(xx=c(0.00000000000001)),data=datf))-q
      bmd.svh=try(uniroot(RISK, interval=c(0,1))$root,silent=T)
      bmd.sv[i]=ifelse((mode(bmd.svh)=="numeric"),(bmd.svh),0.05)
    }

    ## Default prior BMD
    BMD.min <- rep(.Machine$double.xmin, nlevels)
    if(extended == FALSE){
      BMD.max <- rep(1, nlevels)
    }else{
      BMD.max <- rep(maxDose, nlevels)
      if(maxDose <= 1){
        BMD.max <- rep(maxDose*1000, nlevels)
      }
    }

    BMD.mode <- rep(0.5, nlevels)
    ## If info on BMD is given
    if(!is.null(prior.BMD)){

      if(!is.na(prior.BMD[2])){
        BMD.mode = rep(prior.BMD[2]/maxDose, nlevels)
      }else{
        BMD.mode = rep((prior.BMD[1]/maxDose) + (((prior.BMD[3]/maxDose) - prior.BMD[1]/maxDose)/2), nlevels)
      }
      if(!is.na(prior.BMD[1])){
        BMD.min = rep(prior.BMD[1]/maxDose, nlevels)
      }
      if(!is.na(prior.BMD[3])){
        BMD.max = rep(prior.BMD[3]/maxDose, nlevels)
      }

    }else {
      message("Default prior choices used on BMD")
    }

    BMD.vec <- matrix(c(BMD.min, BMD.mode, BMD.max), nrow = nlevels, ncol = 3)
    shape.BMD1 <- shape.BMD2 <- numeric(nlevels)
    for(i in 1:nlevels) {
      shape.BMD1[i] <- fun.alpha(a = BMD.vec[i,1], b = BMD.vec[i,2],
                                 c = BMD.vec[i,3], g = shape.BMD)
      shape.BMD2[i] <- fun.beta(a = BMD.vec[i,1], b = BMD.vec[i,2],
                                c = BMD.vec[i,3], g = shape.BMD)
    }



  } else {


    # start value BMD
    datf=data.frame(yy=mean.a2,
                    xx=dose.a2+0.00000000000001)
    fpfit=gamlss(yy~fp(xx),family=NO(),data=datf)
    RISK=function(x) (predict(fpfit,newdata=data.frame(xx=c((x))),data=datf)-
                        predict(fpfit,newdata=data.frame(xx=c(0.00000000000001)),data=datf))/
      (predict(fpfit,newdata=data.frame(xx=c(0.00000000000001)),data=datf))-q
    bmd.svh=try(uniroot(RISK, interval=c(0,1))$root,silent=T)
    bmd.sv=ifelse((mode(bmd.svh)=="numeric"),(bmd.svh),0.05)
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
    shape.BMD1 <- fun.alpha(a = BMD.vec[1], b = BMD.vec[2],
                            c = BMD.vec[3], g = shape.BMD)
    shape.BMD2 <- fun.beta(a = BMD.vec[1], b = BMD.vec[2],
                           c = BMD.vec[3], g = shape.BMD)
    # dim(bmd.sv) <- 1

  }

  nlevels_BMD <- ifelse(covariate == 'BMD_d' | covariate == 'all', nlevels, 1)
  dim(bmd.sv) <- nlevels_BMD

  if(covariate == 'BMD_d' | covariate == 'all') {

    if(prior.d == 'N11'){
      prvar.d = rep(1, nlevels); prmean.d = rep(1, nlevels); truncd = 5
    }else if(prior.d == 'EPA'){
      # prvar.d = 0.5^2; prmean.d = 0.4; truncd = 10000
      prvar.d = rep(0.5, nlevels); prmean.d = rep(0.4, nlevels); truncd = 10000
      # prvar.d = 1; prmean.d = 1; truncd = 10000
    }
    # prvar.d=sqrt(0.5); prmean.d = prmean.d
    prmean.dQE4 = rep(0, nlevels); prvar.dQE4 = rep(1, nlevels); truncdQ = 10000
    # prvar.d=(exp(sqrt(0.18)))^2; prmean.d=2

  } else {

    if(prior.d == 'N11'){
      prvar.d = 1; prmean.d = 1; truncd = 5
    }else if(prior.d == 'EPA'){
      # prvar.d = 0.5^2; prmean.d = 0.4; truncd = 10000
      prvar.d = 0.5; prmean.d = 0.4; truncd = 10000
      # prvar.d = 1; prmean.d = 1; truncd = 10000
    }else if(prior.d == 'N05'){
      prvar.d = 0.25; prmean.d = 0.5; truncd = 10000
    }
    # prvar.d=sqrt(0.5); prmean.d = prmean.d
    prmean.dQE4 = 0; prvar.dQE4 = 1; truncdQ = 10000
    # prvar.d=(exp(sqrt(0.18)))^2; prmean.d=2
    # dim(prmean.d) <- dim(prmean.dQE4) <- 1
  }

  nlevels_d <- ifelse(covariate == 'BMD_d' | covariate == 'all', nlevels, 1)
  nlevels_b <- ifelse(covariate == 'a_sigma2' | covariate == 'BMD_d' | covariate == 'all', nlevels, 1)
  # dim(prmean.d) <- dim(prmean.dQE4) <- nlevels_d

  prvar.s=1;

  priormu1a <- rbind(ifelse(rep(is.vector(a.vec), nlevels), rep(a.vec[2], nlevels), a.vec[,2]),
                     ifelse(rep(is.vector(BMD.vec), nlevels), rep(BMD.vec[2], nlevels), BMD.vec[,2]),
                     ifelse(rep(is.vector(c.vec), nlevels), rep(c.vec[2], nlevels), c.vec[,2]),
                     ifelse(rep(length(prmean.d) > 1, nlevels), prmean.d, rep(prmean.d, nlevels)),
                     ifelse(rep(length(prmean.s) > 1, nlevels), prmean.s, rep(prmean.s, nlevels))
  )

  priormu1bQ <- rbind(ifelse(rep(is.vector(a.vec), nlevels), rep(a.vec[2], nlevels), a.vec[,2]),
                      ifelse(rep(is.vector(BMD.vec), nlevels), rep(BMD.vec[2], nlevels), BMD.vec[,2]),
                      ifelse(rep(is.vector(c.vec), nlevels), rep(c.vec[2], nlevels), c.vec[,2]),
                      ifelse(rep(length(prmean.dQE4) > 1, nlevels), prmean.dQE4, rep(prmean.dQE4, nlevels)),
                      ifelse(rep(length(prvar.dQE4) > 1, nlevels), prvar.dQE4, rep(prvar.dQE4, nlevels))
  )

  priorSigma1a=diag(c(1,1,1,unique(prvar.d),prvar.s))

  priorSigma1bQ=diag(c(1,1,1,unique(prvar.dQE4),prvar.s))

  priorlb1a = rbind(ifelse(rep(is.vector(a.vec), nlevels), rep(a.vec[1], nlevels),a.vec[,1]),
                    ifelse(rep(is.vector(BMD.vec), nlevels), rep(BMD.vec[1], nlevels),BMD.vec[,1]),
                    ifelse(rep(is.vector(c.vec), nlevels), rep(c.vec[1], nlevels), c.vec[,1]),
                    rep(0, nlevels), rep(0, nlevels))

  priorub1a = rbind(ifelse(rep(is.vector(a.vec), nlevels), rep(a.vec[3], nlevels), a.vec[,3]),
                    ifelse(rep(is.vector(BMD.vec), nlevels), rep(BMD.vec[3], nlevels), BMD.vec[,3]),
                    ifelse(rep(is.vector(c.vec), nlevels), rep(c.vec[3], nlevels), c.vec[,3]),
                    rep(0, nlevels),
                    rep(0, nlevels))
  row.names(priormu1a) <- row.names(priormu1bQ) <- row.names(priorlb1a) <- c('a', 'BMD', 'c', 'd', 's')

  if(!is.null(prior.BMD) & !is.vector(BMD.vec)){
    bmd.sv = BMD.vec[,2]
  }  else if(!is.null(prior.BMD) & is.vector(BMD.vec)) {bmd.sv = BMD.vec[2]}

  is_increasing = 0; is_decreasing = 0
  if((obs.min[1] < obs.max[1]) >= 1){ # increasing
    is_increasing = 1
    L = 1+q+0.01
    # L = 1.01
    U = 0
    data_type = 1
    pars3d = numeric()
    pars3i = max(priormu1a[3,]) - L
    dim(pars3d)=0
    dim(pars3i)=1
  }else if((obs.min[1] > obs.max[1]) >= 1){ # decreasing
    is_decreasing = 1
    L = 0.01
    U = 1-q-0.01
    # U = 0.95
    data_type = 3
    pars3d = max(priormu1a[3,]) / U
    pars3i = numeric()
    dim(pars3d)=1
    dim(pars3i)=0
  } # a third part might be necessary

  trt_ind <- matrix(NA, nrow = N, ncol = nlevels)
  for(i in 1:nlevels) {
    trt_ind[,i] <- as.numeric(covar==covar_lvls[[i]])
  }

  shape1 <- rbind(ifelse(rep(length(shape.a1) > 1, nlevels), shape.a1, rep(shape.a1, nlevels)),
                  ifelse(rep(length(shape.BMD1) > 1, nlevels), shape.BMD1, rep(shape.BMD1, nlevels)),
                  ifelse(rep(length(shape.c1) > 1, nlevels), shape.c1, rep(shape.c1, nlevels)),
                  rep(0, nlevels), rep(0, nlevels))

  shape2 <- rbind(ifelse(rep(length(shape.a2) > 1, nlevels), shape.a2, rep(shape.a2, nlevels)),
                  ifelse(rep(length(shape.BMD2) > 1, nlevels), shape.BMD2, rep(shape.BMD2, nlevels)),
                  ifelse(rep(length(shape.c2) > 1, nlevels), shape.c2, rep(shape.c2, nlevels)),
                  rep(0, nlevels), rep(0, nlevels))

  ## Data in correct format

  sv.a <- if(class(a.vec)[1] == "matrix"){a.vec[,2]}else{a.vec[2]}

  # if(covariate != 'BMD_d'){
  #   N2 = N
  N = N; n = n.a; x = dose.a; m = mean.a; s2 = sd.a^2
  # }else{
  #   N = N2; n = n.a2; x = dose.a2; m = mean.a2; s2 = sd.a2^2
  #   N2 = N*2
  # }

  ret.list <- list(data = list(N=N,
                               # N2=N2,
                               n=n,
                               x=x,
                               m=m,
                               shift=0,
                               s2=s2,
                               maxD=maxDose,q=q,
                               covariate = covar_lvls,
                               priormu=priormu1a,priormuQ=priormu1bQ, trt_ind = trt_ind,
                               nlevels = nlevels,
                               nlevels_a = nlevels_a,
                               nlevels_c = nlevels_c,
                               nlevels_d = nlevels_d,
                               nlevels_BMD = nlevels_BMD,
                               nlevels_sigma = nlevels_sigma,
                               nlevels_b = nlevels_b,
                               shape1 = shape1,
                               shape2 = shape2,
                               priorlb = priorlb1a, priorub = priorub1a, shape.a = shape.a,
                               shape.c = shape.c, shape.BMD = shape.BMD,
                               priorSigma=priorSigma1a, priorSigmaQ=priorSigma1bQ, init_b = 1,
                               data_type = data_type, L = L, U = U,
                               is_increasing = is_increasing, truncd = truncd, truncdQ = truncdQ,
                               is_decreasing = is_decreasing, is_informative_a = is_informative_a,
                               is_informative_c = is_informative_c,
                               is_informative_BMD = is_informative_BMD),
                   # start values
                   start=list(par1=sv.a,
                              par2=bmd.sv, pars3i=pars3i,
                              pars3d=pars3d, par4=prmean.d,
                              par5= par5),
                   startQ=list(par1=sv.a,
                               par2=bmd.sv, pars3i=pars3i, pars3d=pars3d,
                               par4=prmean.dQE4,
                               par5=par5),
                   test.var = test.var
  )
  dim(ret.list$start$par2) <- nlevels_BMD
  dim(ret.list$start$par1) <- nlevels_a
  dim(ret.list$startQ$par2) <- nlevels_BMD
  dim(ret.list$startQ$par1) <- nlevels_a
  dim(ret.list$startQ$par4) <- nlevels_d
  dim(ret.list$start$par4) <- nlevels_d




  # test for dose-response effect
  # DR.effect = anydoseresponseNI(dose.a,mean.a,sd.a,n.a)

  # data in correct format
  return(ret.list)

}

#' @rdname PREP_DATA_NCOV
#' @export
PREP_DATA_LNCOV <- function(data, # a dataframe with input data, order of columns should be: dose, response, sd, n, covar
                            sumstats = TRUE, # TRUE if summary data, FALSE if individual data
                            geom.stats = FALSE, # TRUE if geometric summary data is provided
                            sd = TRUE, # TRUE if sd per dose group is given, FALSE is se is given
                            q, # the BMR
                            bkg = NULL,
                            maxy = NULL,
                            prior.BMD = NULL, # possible expert info on background and max response
                            shape.a = 4, shape.c = 4, shape.BMD = 0.0001, # shape for the PERT distribution,
                            # prmean.d = 1, prmean.dQE4 = 0
                            prior.d = 'N11',
                            extended = TRUE,
                            # covariate = c('a', 'BMD', 'sigma2', 'd')
                            covariate = 'all' # OPTIONS: 'a_sigma2', 'BMD_d', 'all'; for 'none', use original prep_data
){

  if(sumstats == TRUE & geom.stats == FALSE){
    data = as.data.frame(data[order(data[, 1]), ])
    dose.a = data[, 1]
    maxDose = max(dose.a)
    mean.a = data[, 2]
    if(sd == TRUE){
      sd.a = data[, 3]
    }else if(sd == FALSE){
      sd.a = data[,3]*sqrt(data[, 4]) # SD = SE * sqrt(n.a)
    }
    n.a = data[, 4]
    # N = length(unique(dose.a))
    N = length(dose.a)
    dose.a = dose.a/maxDose
    covar = data[,5]
    # shift if negative means occur
    shift=0
    gmean.a2 = log(NtoLN(mean.a,sd.a))[1:N]
    if (min(gmean.a2)<0) {gmean.a = gmean.a2-20*min(gmean.a2); shift = 20*min(gmean.a2)}
    if (min(gmean.a2)>=0) gmean.a = gmean.a2
    gsd.a = log(NtoLN(mean.a,sd.a))[(N+1):(2*N)]

  }else if(sumstats == TRUE & geom.stats == TRUE){

    data = data[order(data[, 1]), ]
    dose.a = data[, 1]
    maxDose = max(dose.a)
    mean.a = data[, 2]
    if(sd == TRUE){
      sd.a = data[, 3]
    }else if(sd == FALSE){
      sd.a = data[,3]*sqrt(data[, 4]) # SD = SE * sqrt(n.a)
    }
    gsd.a = sd.a
    gmean.a = mean.a
    gmean.a2 = mean.a
    shift = 0
    n.a = data[, 4]
    N = length(dose.a)
    covar = data[,5]
    dose.a = dose.a/maxDose

  }else if(sumstats == FALSE){
    # still to be worked on

    stop('Data should be provided as summary statistics.')

  }

  covar_lvls <- unique(covar)
  nlevels <- length(unique(covar))
  ## Bartlett test of homoscedasticity
  # on the original scale (constant variance)
  #b.test.LN <- numeric(nlevels)

  # if(covariate == 'BMD_d' | covariate == 'none'){
  ## get overall data for values of a and sigma
  dose.a2 = sort(unique(dose.a))
  N2 = length(dose.a2)
  mean.a2 = rep(NA, N2)
  sd.a2 = rep(NA, N2)
  n.a2 = rep(NA, N2)
  for(iu in (1:N2)){
    mean.a2[iu] = mean(mean.a[dose.a == dose.a2[iu]])
    sd.a2[iu] = mean(sd.a[dose.a == dose.a2[iu]])
    n.a2[iu] = sum(n.a[dose.a == dose.a2[iu]])
  }
  gmean.a3 = log(NtoLN(mean.a2,sd.a2))[1:N2]
  if (min(gmean.a3)<0) {gmean.a4 = gmean.a3-shift}
  if (min(gmean.a3)>=0) gmean.a4 = gmean.a3
  gsd.a2 = log(NtoLN(mean.a2,sd.a2))[(N2+1):(2*N2)]
  # }

  if(covariate == 'a_sigma2' | covariate == 'all') {
    test.var <- character(nlevels)
    prmean.s <- par5 <- numeric(nlevels)
    for(i in 1:nlevels){
      b.test.LN <- bartlett(gsd.a[covar == covar_lvls[i]], n.a[covar == covar_lvls[i]])
      if(b.test.LN[2]>=0.05){
        test.var[i] = paste0('Distributional assumption of constant coefficient of variation is met for group ', covar_lvls[i],
                             ' Bartlett test p-value is ',
                             round(b.test.LN[2], 4))
        message(test.var[i])
      }else if(b.test.LN[2]<0.05){
        test.var[i] = paste0('Distributional assumption of constant coefficient of variation for the lognormal distribution is not met for group ',
                             covar_lvls[i], ' Bartlett test p-value is ', round(b.test.LN[2], 4))
        warning(test.var[i])
      }
      prmean.s[i] =-2*log(1.5*mean(gsd.a[covar == covar_lvls[i]]))
      par5[i] <- log(1/mean(gsd.a[covar == covar_lvls[i]]^2))
    }

  } else {

    b.test.LN <- bartlett(gsd.a2, n.a2)

    if(b.test.LN[2]>=0.05){
      test.var = paste0('Distributional assumption of constant coefficient of variation are met for ',
                        'Bartlett test p-value is ',
                        round(b.test.LN[2], 4))
      message(test.var)
    }else if(b.test.LN[2]<0.05){
      test.var = paste0('Distributional assumption of constant coefficient variaiton for the lognormal distribution is not met for ',
                        ' Bartlett test p-value is ', round(b.test.LN[2], 4))
      warning(test.var)
    }

    prmean.s = -2*log(1.5*mean(gsd.a2))
    par5 <- log(1/mean(gsd.a2^2))
    # dim(par5) <- 1
  }

  nlevels_sigma <- ifelse(covariate == 'a_sigma2' | covariate == 'all', nlevels, 1)
  dim(par5) <- nlevels_sigma

  is_informative_BMD = 0
  is_informative_a = 0
  is_informative_c = 0

  if(!is.null(prior.BMD)) {is_informative_BMD = 1}


  ######################
  ### PRIORS

  ## Observed background and maximum response
  ### Family 1
  if(covariate == 'a_sigma2' | covariate == 'all') {

    obs.min <- obs.max <- min.min <- mode.min <- max.min <- min.max <- max.max <- mode.max <- numeric(nlevels)
    for(i in 1:nlevels) {
      obs.min[i] = mean.a[dose.a == 0 & covar == covar_lvls[i]][1]
      obs.max[i] = mean.a[dose.a == max(dose.a) & covar == covar_lvls[i]][1]

      if(obs.min[i] < obs.max[i]) {

        # for a
        min.min = rep(0.001, nlevels)
        mode.min[i] = obs.min[i]
        max.min[i] = 2*obs.min[i]

        # for c
        min.max[i] = obs.min[i]*(1.01+q)
        max.max[i] = 2*obs.max[i]
        mode.max[i] = obs.max[i]

        if(flat(dose.a[covar == covar_lvls[i]], mean.a[covar == covar_lvls[i]],
                n.a[covar == covar_lvls[i]], inc=T) == F & is.null(maxy)){
          mode.max[i] = 3*obs.max[i]
          min.max[i]  = obs.min[i]*(1.01+q)
          max.max[i] = 2*mode.max[i]

        }


        ## Check appropriateness of BMR value
        if(obs.min[i]*(1+q) > obs.max[i]){
          warning(paste0("The data do not contain values corresponding to the chosen BMR for group ", i,
                         " lowering the specified value of q may be necessary."))
        }

      } else if(obs.min[i] > obs.max[i]){

        min.min[i] = 0.5*obs.min[i]
        mode.min[i] = obs.min[i]
        max.min[i] = 2*obs.min[i]

        min.max[i] = 0.5*obs.max[i]
        max.max[i] = obs.min[i]*(1-q-0.01)
        # max.max = obs.min
        mode.max[i] = obs.max[i]

        if(flat(dose.a, mean.a, n.a, inc=F) == F & is.null(maxy)){
          mode.max[i] = 0.5*obs.max[i]
          min.max[i]  = 0.1*obs.max[i]
          max.max[i] = obs.min[i]*(1-q-0.01)
        }

        ## Check appropriateness of BMR value
        if(obs.min[i]*(1-q) < obs.max[i]){
          warning(paste0("The data do not contain values corresponding to the chosen BMR for group ", i,
                         " lowering the specified value of q may be necessary."))
        }

      }

      ## If info on background is given
      if(!is.null(bkg)){

        is_informative_a = 1

        if(!is.na(bkg[2])){
          mode.min[i] = bkg[2]
        }else{
          mode.min[i] = bkg[1] + ((bkg[3]-bkg[1])/2)
        }
        if(!is.na(bkg[1])){
          min.min[i] = bkg[1]
        }
        if(!is.na(bkg[3])){
          max.min[i] = bkg[3]
        }

      }else{
        message("Default prior choices used on background")
      }

      ## If info on max response is given
      if(!is.null(maxy)){

        is_informative_c = 1

        if(!is.na(maxy[2])){
          mode.max[i] = maxy[2]
        }else{
          mode.max[i] = maxy[1] + ((maxy[3]-maxy[1])/2)
        }
        if(!is.na(maxy[1])){
          min.max[i] = maxy[1]
        }
        if(!is.na(maxy[3])){
          max.max[i] = maxy[3]
        }

      }else{
        message("Default prior choices used on fold change")
      }

    }

    # Prior on background
    a.vec = matrix(c(min.min, mode.min, max.min), nrow = nlevels, ncol = 3)
    # Prior on mu(inf)
    c.vec = matrix(c(min.max/mode.min, mode.max/mode.min, max.max/mode.min), nrow = nlevels, ncol = 3)

    shape.a1 <- shape.a2 <- shape.c1 <- shape.c2 <- numeric(nlevels)

    for(i in 1:nlevels){
      if(c.vec[i,1] == 0) c.vec[i,1] = 0.0001
      if(c.vec[i,2] >= c.vec[i,3]) c.vec[i,2] = c.vec[i,3] - 0.05

      shape.a1[i] <- fun.alpha(a = a.vec[i,1], b = a.vec[i,2],
                               c = a.vec[i,3], g = shape.a)
      shape.c1[i] <- fun.alpha(a = c.vec[i,1], b = c.vec[i,2],
                               c = c.vec[i,3], g = shape.c)

      shape.a2[i] <- fun.beta(a = a.vec[i,1], b = a.vec[i,2],
                              c = a.vec[i,3], g = shape.a)
      shape.c2[i] <- fun.beta(a = c.vec[i,1], b = c.vec[i,2],
                              c = c.vec[i,3], g = shape.c)
    }


  } else {

    obs.min = mean.a2[dose.a2 == 0][1]
    obs.max = mean.a2[dose.a2 == max(dose.a)][1]

    if(obs.min < obs.max) {

      # for a
      min.min = 0.001
      mode.min = obs.min
      max.min = 2*obs.min

      # for c
      min.max = obs.min*(1.01+q)
      max.max = 2*obs.max
      mode.max = obs.max

      if(flat(dose.a2, mean.a2,
              n.a2, inc=T) == F & is.null(maxy)){
        mode.max = 3*obs.max
        min.max  = obs.min*(1.01+q)
        max.max = 2*mode.max

      }

      ## Check appropriateness of BMR value
      if(obs.min*(1+q) > obs.max){
        warning('The data do not contain values corresponding to the chosen BMR,
                lowering the specified value of q may be necessary.')
      }

    } else if(obs.min > obs.max){

      min.min = 0.5*obs.min
      mode.min = obs.min
      max.min = 2*obs.min

      min.max = 0.5*obs.max
      max.max = obs.min*(1-q-0.01)
      # max.max = obs.min
      mode.max = obs.max

      if(flat(dose.a2, mean.a2, n.a2, inc=F) == F & is.null(maxy)){
        mode.max = 0.5*obs.max
        min.max = 0.1*obs.max
        max.max = obs.min*(1-q-0.01)
      }


      ## Check appropriateness of BMR value
      if(obs.min*(1-q) < obs.max){
        warning('The data do not contain values corresponding to the chosen BMR,
                lowering the specified value of q may be necessary.')
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


    a.vec = c(min.min, mode.min, max.min)
    # Prior on mu(inf)
    c.vec = c(min.max/mode.min, mode.max/mode.min, max.max/mode.min)

    if(c.vec[1] == 0) c.vec[1] = 0.0001
    if(c.vec[2] >= c.vec[3]) c.vec[2] = c.vec[3] - 0.05

    shape.a1 <- fun.alpha(a = a.vec[1], b = a.vec[2],
                          c = a.vec[3], g = shape.a)
    shape.a2 <- fun.beta(a = a.vec[1], b = a.vec[2],
                         c = a.vec[3], g = shape.a)

    shape.c1 <- fun.alpha(a = c.vec[1], b = c.vec[2],
                          c = c.vec[3], g = shape.c)
    shape.c2 <- fun.beta(a = c.vec[1], b = c.vec[2],
                         c = c.vec[3], g = shape.c)
    # dim(a.vec) <- dim(c.vec) <- 1
  }

  nlevels_a <- nlevels_c <- ifelse(covariate == 'a_sigma2' | covariate == 'all', nlevels, 1)
  # dim(a.vec) <- nlevels_a
  # dim(c.vec) <- nlevels_c


  if(covariate == 'BMD_d' | covariate == 'all') {

    bmd.sv <- numeric(nlevels)

    for(i in 1:nlevels){
      # start value BMD
      datf=data.frame(yy=mean.a[covar == covar_lvls[i]],
                      xx=dose.a[covar == covar_lvls[i]]+0.00000000000001)
      fpfit=gamlss(yy~fp(xx),family=NO(),data=datf)
      RISK=function(x) (predict(fpfit,newdata=data.frame(xx=c((x))),data=datf)-
                          predict(fpfit,newdata=data.frame(xx=c(0.00000000000001)),data=datf))/
        (predict(fpfit,newdata=data.frame(xx=c(0.00000000000001)),data=datf))-q
      bmd.svh=try(uniroot(RISK, interval=c(0,1))$root,silent=T)
      bmd.sv[i]=ifelse((mode(bmd.svh)=="numeric"),(bmd.svh),0.05)
    }

    ## Default prior BMD
    BMD.min <- rep(.Machine$double.xmin, nlevels)
    if(extended == FALSE){
      BMD.max <- rep(1, nlevels)
    }else{
      BMD.max <- rep(maxDose, nlevels)
      if(maxDose <= 1){
        BMD.max <- rep(maxDose*1000, nlevels)
      }
    }

    BMD.mode <- rep(0.5, nlevels)
    ## If info on BMD is given
    if(!is.null(prior.BMD)){

      if(!is.na(prior.BMD[2])){
        BMD.mode = rep(prior.BMD[2]/maxDose, nlevels)
      }else{
        BMD.mode = rep((prior.BMD[1]/maxDose) + (((prior.BMD[3]/maxDose) - prior.BMD[1]/maxDose)/2), nlevels)
      }
      if(!is.na(prior.BMD[1])){
        BMD.min = rep(prior.BMD[1]/maxDose, nlevels)
      }
      if(!is.na(prior.BMD[3])){
        BMD.max = rep(prior.BMD[3]/maxDose, nlevels)
      }

    }else {
      message("Default prior choices used on BMD")
    }

    BMD.vec <- matrix(c(BMD.min, BMD.mode, BMD.max), nrow = nlevels, ncol = 3)
    shape.BMD1 <- shape.BMD2 <- numeric(nlevels)
    for(i in 1:nlevels) {
      shape.BMD1[i] <- fun.alpha(a = BMD.vec[i,1], b = BMD.vec[i,2],
                                 c = BMD.vec[i,3], g = shape.BMD)
      shape.BMD2[i] <- fun.beta(a = BMD.vec[i,1], b = BMD.vec[i,2],
                                c = BMD.vec[i,3], g = shape.BMD)
    }



  } else {


    # start value BMD
    datf=data.frame(yy=mean.a2,
                    xx=dose.a2+0.00000000000001)
    fpfit=gamlss(yy~fp(xx),family=NO(),data=datf)
    RISK=function(x) (predict(fpfit,newdata=data.frame(xx=c((x))),data=datf)-
                        predict(fpfit,newdata=data.frame(xx=c(0.00000000000001)),data=datf))/
      (predict(fpfit,newdata=data.frame(xx=c(0.00000000000001)),data=datf))-q
    bmd.svh=try(uniroot(RISK, interval=c(0,1))$root,silent=T)
    bmd.sv=ifelse((mode(bmd.svh)=="numeric"),(bmd.svh),0.05)
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
    shape.BMD1 <- fun.alpha(a = BMD.vec[1], b = BMD.vec[2],
                            c = BMD.vec[3], g = shape.BMD)
    shape.BMD2 <- fun.beta(a = BMD.vec[1], b = BMD.vec[2],
                           c = BMD.vec[3], g = shape.BMD)
    # dim(bmd.sv) <- 1

  }

  nlevels_BMD <- ifelse(covariate == 'BMD_d' | covariate == 'all', nlevels, 1)
  dim(bmd.sv) <- nlevels_BMD

  if(covariate == 'BMD_d' | covariate == 'all') {

    if(prior.d == 'N11'){
      prvar.d = rep(1, nlevels); prmean.d = rep(1, nlevels); truncd = 5
    }else if(prior.d == 'EPA'){
      # prvar.d = 0.5^2; prmean.d = 0.4; truncd = 10000
      prvar.d = rep(0.5, nlevels); prmean.d = rep(0.4, nlevels); truncd = 10000
      # prvar.d = 1; prmean.d = 1; truncd = 10000
    }else if(prior.d == 'N05'){
      prvar.d = rep(0.25, nlevels); prmean.d = rep(0.5, nlevels); truncd = 10000
    }
    # prvar.d=sqrt(0.5); prmean.d = prmean.d
    prmean.dQE4 = rep(0, nlevels); prvar.dQE4 = rep(1, nlevels); truncdQ = 10000
    # prvar.d=(exp(sqrt(0.18)))^2; prmean.d=2

  } else {

    if(prior.d == 'N11'){
      prvar.d = 1; prmean.d = 1; truncd = 5
    }else if(prior.d == 'EPA'){
      # prvar.d = 0.5^2; prmean.d = 0.4; truncd = 10000
      prvar.d = 0.5; prmean.d = 0.4; truncd = 10000
      # prvar.d = 1; prmean.d = 1; truncd = 10000
    }else if(prior.d == 'N05'){
      prvar.d = 0.25; prmean.d = 0.5; truncd = 10000
    }
    # prvar.d=sqrt(0.5); prmean.d = prmean.d
    prmean.dQE4 = 0; prvar.dQE4 = 1; truncdQ = 10000
    # prvar.d=(exp(sqrt(0.18)))^2; prmean.d=2
    # dim(prmean.d) <- dim(prmean.dQE4) <- 1
  }

  nlevels_d <- ifelse(covariate == 'BMD_d' | covariate == 'all', nlevels, 1)
  nlevels_b <- ifelse(covariate == 'a_sigma2' | covariate == 'BMD_d' | covariate == 'all', nlevels, 1)
  # dim(prmean.d) <- dim(prmean.dQE4) <- nlevels_d

  prvar.s=1;

  priormu1a <- rbind(ifelse(rep(is.vector(a.vec), nlevels), rep(a.vec[2], nlevels), a.vec[,2]),
                     ifelse(rep(is.vector(BMD.vec), nlevels), rep(BMD.vec[2], nlevels), BMD.vec[,2]),
                     ifelse(rep(is.vector(c.vec), nlevels), rep(c.vec[2], nlevels), c.vec[,2]),
                     ifelse(rep(length(prmean.d) > 1, nlevels), prmean.d, rep(prmean.d, nlevels)),
                     ifelse(rep(length(prmean.s) > 1, nlevels), prmean.s, rep(prmean.s, nlevels))
  )

  priormu1bQ <- rbind(ifelse(rep(is.vector(a.vec), nlevels), rep(a.vec[2], nlevels), a.vec[,2]),
                      ifelse(rep(is.vector(BMD.vec), nlevels), rep(BMD.vec[2], nlevels), BMD.vec[,2]),
                      ifelse(rep(is.vector(c.vec), nlevels), rep(c.vec[2], nlevels), c.vec[,2]),
                      ifelse(rep(length(prmean.dQE4) > 1, nlevels), prmean.dQE4, rep(prmean.dQE4, nlevels)),
                      ifelse(rep(length(prvar.dQE4) > 1, nlevels), prvar.dQE4, rep(prvar.dQE4, nlevels))
  )

  priorSigma1a=diag(c(1,1,1,unique(prvar.d),prvar.s))

  priorSigma1bQ=diag(c(1,1,1,unique(prvar.dQE4),prvar.s))

  priorlb1a = rbind(ifelse(rep(is.vector(a.vec), nlevels), rep(a.vec[1], nlevels),a.vec[,1]),
                    ifelse(rep(is.vector(BMD.vec), nlevels), rep(BMD.vec[1], nlevels),BMD.vec[,1]),
                    ifelse(rep(is.vector(c.vec), nlevels), rep(c.vec[1], nlevels), c.vec[,1]),
                    rep(0, nlevels), rep(0, nlevels))

  priorub1a = rbind(ifelse(rep(is.vector(a.vec), nlevels), rep(a.vec[3], nlevels), a.vec[,3]),
                    ifelse(rep(is.vector(BMD.vec), nlevels), rep(BMD.vec[3], nlevels), BMD.vec[,3]),
                    ifelse(rep(is.vector(c.vec), nlevels), rep(c.vec[3], nlevels), c.vec[,3]),
                    rep(0, nlevels),
                    rep(0, nlevels))
  row.names(priormu1a) <- row.names(priormu1bQ) <- row.names(priorlb1a) <- c('a', 'BMD', 'c', 'd', 's')

  if(!is.null(prior.BMD) & !is.vector(BMD.vec)){
    bmd.sv = BMD.vec[,2]
  }  else if(!is.null(prior.BMD) & is.vector(BMD.vec)) {bmd.sv = BMD.vec[2]}

  is_increasing = 0; is_decreasing = 0
  if((obs.min[1] < obs.max[1]) >= 1){ # increasing
    is_increasing = 1
    L = 1+q+0.01
    # L = 1.01
    U = 0
    data_type = 2
    pars3d = numeric()
    pars3i = max(priormu1a[3,]) - L
    dim(pars3d)=0
    dim(pars3i)=1
  }else if((obs.min[1] > obs.max[1]) >= 1){ # decreasing
    is_decreasing = 1
    L = 0.01
    U = 1-q-0.01
    # U = 0.95
    data_type = 4
    pars3d = max(priormu1a[3,]) / U
    pars3i = numeric()
    dim(pars3d)=1
    dim(pars3i)=0
  } # a third part might be necessary

  trt_ind <- matrix(NA, nrow = N, ncol = nlevels)
  for(i in 1:nlevels) {
    trt_ind[,i] <- as.numeric(covar==covar_lvls[[i]])
  }

  shape1 <- rbind(ifelse(rep(length(shape.a1) > 1, nlevels), shape.a1, rep(shape.a1, nlevels)),
                  ifelse(rep(length(shape.BMD1) > 1, nlevels), shape.BMD1, rep(shape.BMD1, nlevels)),
                  ifelse(rep(length(shape.c1) > 1, nlevels), shape.c1, rep(shape.c1, nlevels)),
                  rep(0, nlevels), rep(0, nlevels))

  shape2 <- rbind(ifelse(rep(length(shape.a2) > 1, nlevels), shape.a2, rep(shape.a2, nlevels)),
                  ifelse(rep(length(shape.BMD2) > 1, nlevels), shape.BMD2, rep(shape.BMD2, nlevels)),
                  ifelse(rep(length(shape.c2) > 1, nlevels), shape.c2, rep(shape.c2, nlevels)),
                  rep(0, nlevels), rep(0, nlevels))

  ## Data in correct format

  sv.a <- if(class(a.vec)[1] == "matrix"){a.vec[,2]}else{a.vec[2]}

  # if(covariate != 'BMD_d'){
  #   N2 = N
  N = N; n = n.a; x = dose.a; m = gmean.a; m.org = gmean.a2; s2 = gsd.a^2
  # }else{
  #   N = N2; n = n.a2; x = dose.a2; m = mean.a2; s2 = sd.a2^2
  #   N2 = N*2
  # }

  ret.list <- list(data = list(N=N,
                               # N2=N2,
                               n=n,
                               x=x,
                               m=m,
                               shift=0,
                               s2=s2,
                               maxD=maxDose,q=q,
                               covariate = covar_lvls,
                               priormu=priormu1a,priormuQ=priormu1bQ, trt_ind = trt_ind,
                               nlevels = nlevels,
                               nlevels_a = nlevels_a,
                               nlevels_c = nlevels_c,
                               nlevels_d = nlevels_d,
                               nlevels_BMD = nlevels_BMD,
                               nlevels_sigma = nlevels_sigma,
                               nlevels_b = nlevels_b,
                               shape1 = shape1,
                               shape2 = shape2,
                               priorlb = priorlb1a, priorub = priorub1a, shape.a = shape.a,
                               shape.c = shape.c, shape.BMD = shape.BMD,
                               priorSigma=priorSigma1a, priorSigmaQ=priorSigma1bQ, init_b = 1,
                               data_type = data_type, L = L, U = U,
                               is_increasing = is_increasing, truncd = truncd, truncdQ = truncdQ,
                               is_decreasing = is_decreasing, is_informative_a = is_informative_a,
                               is_informative_c = is_informative_c,
                               is_informative_BMD = is_informative_BMD),
                   # start values
                   start=list(par1=sv.a,
                              par2=bmd.sv, pars3i=pars3i,
                              pars3d=pars3d, par4=prmean.d,
                              par5= par5),
                   startQ=list(par1=sv.a,
                               par2=bmd.sv, pars3i=pars3i, pars3d=pars3d,
                               par4=prmean.dQE4,
                               par5=par5),
                   test.var = test.var
  )
  dim(ret.list$start$par2) <- nlevels_BMD
  dim(ret.list$start$par1) <- nlevels_a
  dim(ret.list$startQ$par2) <- nlevels_BMD
  dim(ret.list$startQ$par1) <- nlevels_a
  dim(ret.list$startQ$par4) <- nlevels_d
  dim(ret.list$start$par4) <- nlevels_d




  # test for dose-response effect
  # DR.effect = anydoseresponseNI(dose.a,mean.a,sd.a,n.a)

  # data in correct format
  return(ret.list)

}


#' Function to set data in the correct format
#'
#' This function also generates appropriate start values and uninformative priors for each model
#' Input should be given as arithmetic mean and standard deviation on the original scale
#'
#' @param data a dataframe with input data, order of columns should be: dose, number of adverse events, n
#' @param sumstats logical indicating whether summary (T, default) or individual-level (F) data is provided
#' @param q specified BMR
#' @param bkg vector containing minimum, most likely (optional), and maximum value for the background response. Defaults to NULL (non-informative prior)
#' @param prior.BMD vector containing minimum, most likely (optional), and maximum value for the BMD. Defaults to NULL (non-informative prior)
#' @param shape.a shape parameter for the modified PERT distribution on parameter a, defaults to 4 (peaked at most likely value), a value of 0.0001 implies a uniform distribution
#' @param shape.BMD shape parameter for the modified PERT distribution on parameter BMD, defaults to 0.0001 implying a uniform distribution. Can be set to 4 in case of informative prior
#' @param cluster logical variable to indicate if data is clustered. TRUE = clustered data. Defaults to FALSE
#' @param prior.d prior distribution for parameter d, should be either N11 (default), EPA or N05 (for a N(0.5,0.5) prior)
#' @param extended logical indicating whether the dose range should be extended to maxDose^2 (default is TRUE)
#' @param covariate for which parameter a covariate effect should be included. Defaults to 'all', other options are 'background' or 'BMD_d'
#'
#' @return List with data and start values in correct format to be directly used within the BMA functions.
#'
#' @export PREP_DATA_Q_COV
#'
PREP_DATA_Q_COV <- function(data, # a dataframe with input data, order of columns should be: dose, response (=count?), n
                            sumstats = TRUE, # TRUE if summary data, FALSE if individual data
                            q, # the BMR,
                            bkg = NULL, # possible expert info on background response (a),
                            prior.BMD = NULL, # possible expert prior on the BMD,
                            shape.a = 4, #scale parameter for Pert priors for a
                            shape.BMD = 0.0001, #scale parameter for the Pert priors for BMD
                            cluster = FALSE, # indicate if data is clustered
                            prior.d = 'N11',
                            extended = TRUE,
                            covariate = 'all' # options are 'all', 'BMD_d' or 'background'
){

  if(sumstats == TRUE){
    data = data[order(data[, 1]), ]
    dose.a = data[, 1]
    maxDose = max(dose.a)
    y.a = data[, 2]
    n.a = data[, 3]
    N = length(dose.a)
    dose.a = dose.a/maxDose
    covar = data[, 4]
  } else if(sumstats == FALSE){
    ## STILL TO BE DONE
    doses = data[, 1]
    maxDose = max(doses)
    dose.a = sort(unique(doses))
    N = length(dose.a)
    y.a = rep(NA,N)
    n.a = rep(NA,N)
    ybin = data[, 2]
    for (iu in (1:N)){
      y.a[iu] = sum(ybin[doses == dose.a[iu]])
      n.a[iu] = sum(doses == dose.a[iu])
    }
    dose.a = dose.a/maxDose
  }

  covar_lvls = unique(covar)
  nlevels = length(unique(covar))

  if(covariate == 'BMD_d' | covariate == 'none'){
    # get overall background value
    dose.a2 = sort(unique(dose.a))
    N2 = length(dose.a2)
    y.a2 = rep(NA, N2)
    n.a2 = rep(NA, N2)
    for(iu in (1:N2)){
      y.a2[iu] = sum(y.a[dose.a == dose.a2[iu]])
      n.a2[iu] = sum(n.a[dose.a == dose.a2[iu]])
    }
  }

  ## clustered data option
  is_bin <- ifelse(cluster == FALSE, 1, 0)
  is_betabin <- ifelse(cluster == TRUE, 1, 0)

  ######################
  ### PRIORS

  ## Observed background and maximum response

  #obs.min = ifelse(y.a[1]/n.a[1] == 0, 1.0e-03, y.a[1]/n.a[1])

  is_informative_BMD = 0
  is_informative_a = 0

  if(covariate == 'background' | covariate == 'all'){

    miny.a <- minn.a <- a.mode <- a.min <- a.max <- numeric(nlevels)

    for(i in 1:nlevels){

      ## Default range on background
      miny.a[i] <- sum(y.a[dose.a == 0 & covar == covar_lvls[i]])
      minn.a[i] <- sum(n.a[dose.a == 0 & covar == covar_lvls[i]])

      a.min[i] <- ifelse(miny.a[i] != 0,
                         max(c(prop.test(miny.a[i], minn.a[i])$conf.int[1]/2, 1/(10*minn.a[i]))),
                         .Machine$double.xmin)
      a.max[i] <- min(c(3*prop.test(miny.a[1], minn.a[1])$conf.int[2]/2, 1 - 1/(10*minn.a[i])))
      a.mode[i] <- max(c(miny.a[i]/minn.a[i], 1/(5*minn.a[i])))

    }

    ## If info on background is given (STILL HAS TO BE IMPLEMENTED FOR COVARIATE EFFECT)
    if(!is.null(bkg)){

      # is_informative_a = 1

      if(!is.na(bkg[2])){
        a.mode = bkg[2]
      }
      if(!is.na(bkg[1])){
        a.min = bkg[1]
      }
      if(!is.na(bkg[3])){
        a.max = bkg[3]
      }

      is_informative_a = 1

    }else {
      message("Default prior choices used on background")
    }

    # Prior on a
    a.vec = matrix(c(a.min, a.mode, a.max), nrow = nlevels, ncol = 3)

    shape.a = rep(shape.a, nlevels)

  }else{

    ## Default range on background
    miny.a <- sum(y.a2[dose.a == 0])
    minn.a <- sum(n.a2[dose.a == 0])

    a.min <- ifelse(miny.a != 0,
                    max(c(prop.test(miny.a, minn.a)$conf.int[1]/2, 1/(10*minn.a))),
                    .Machine$double.xmin)
    a.max <- min(c(3*prop.test(miny.a[1], minn.a[1])$conf.int[2]/2, 1 - 1/(10*minn.a)))
    a.mode <- max(c(miny.a/minn.a, 1/(5*minn.a)))


    ## If info on background is given (STILL HAS TO BE IMPLEMENTED FOR COVARIATE EFFECT)
    if(!is.null(bkg)){

      # is_informative_a = 1

      if(!is.na(bkg[2])){
        a.mode = bkg[2]
      }
      if(!is.na(bkg[1])){
        a.min = bkg[1]
      }
      if(!is.na(bkg[3])){
        a.max = bkg[3]
      }

      is_informative_a = 1

    }else {
      message("Default prior choices used on background")
    }

    # Prior on a
    a.vec = c(a.min, a.mode, a.max)
    shape.a = shape.a

  }

  nlevels_a <- ifelse(covariate == 'background' | covariate == 'all', nlevels, 1)

  if(covariate == 'BMD_d' | covariate == 'all'){

    bmd.sv <- numeric(nlevels)

    for(i in 1:nlevels){

      datf = data.frame(yy = y.a[covar == covar_lvls[i]], n.a = n.a[covar == covar_lvls[i]], xx = dose.a[covar == covar_lvls[i]])
      if(cluster == FALSE) {
        fpfit = gamlss(cbind(yy, n.a-yy)~fp(xx),family=BI(mu.link = 'logit'),data=datf)
      }
      # else if(cluster == TRUE) {
      #   fpfit = gamlss(cbind(yy, n.a-yy)~fp(xx),family=BB,
      #                  sigma.formula = ~1,data=datf)
      #   fpfit2 <- try(gamlss(cbind(yy,n.a-yy)~as.factor(xx), sigma.formula=~1, family=BB, data=datf),
      #          silent = TRUE)
      #   rhohat <- exp(fpfit2$sigma.coefficients)/(exp(fpfit2$sigma.coefficients)+1)
      #
      # } else stop('provide cluster to be TRUE or FALSE')

      RISK = function(x) (predict(fpfit,newdata=data.frame(xx=c(exp(x))), data=datf, type = "response")-
                            predict(fpfit,newdata=data.frame(xx=c(0.00000000000001)), data=datf, type = "response"))/
        (1 - (predict(fpfit,newdata=data.frame(xx=c(0.00000000000001)), data=datf, type = "response")))-q
      bmd.svh = try(uniroot(RISK, interval=c(-5, 0))$root,silent=TRUE)
      bmd.sv[i] <- ifelse((mode(bmd.svh)=="numeric"),(exp(bmd.svh)+0.5)/2,0.05) ## why (exp(bmd.svh)+0.5)/2 ?

    }

    ## Default prior BMD
    BMD.min <- rep(.Machine$double.xmin, nlevels)
    if(extended == FALSE){
      BMD.max <- rep(1, nlevels)
    }else{
      BMD.max <- rep(maxDose, nlevels)
      if(maxDose <= 1){
        BMD.max <- rep(maxDose*1000, nlevels)
      }
    }
    BMD.mode <- rep(0.5, nlevels)

    ## If info on BMD is given (STILL HAS TO BE ADAPTED TO COVARIATE EFFECT)
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

      is_informative_BMD = 1

    }else {
      message("Default prior choices used on BMD")
    }

    BMD.vec <- matrix(c(BMD.min, BMD.mode, BMD.max), nrow = nlevels, ncol = 3)
    shape.BMD <- rep(shape.BMD, nlevels)

  }else{

    datf = data.frame(yy = y.a, n.a = n.a, xx = dose.a)
    if(cluster == FALSE) {
      fpfit = gamlss(cbind(yy, n.a-yy)~fp(xx),family=BI(mu.link = 'logit'),data=datf)
    }
    # else if(cluster == TRUE) {
    #   fpfit = gamlss(cbind(yy, n.a-yy)~fp(xx),family=BB,
    #                  sigma.formula = ~1,data=datf)
    #   fpfit2 <- try(gamlss(cbind(yy,n.a-yy)~as.factor(xx), sigma.formula=~1, family=BB, data=datf),
    #          silent = TRUE)
    #   rhohat <- exp(fpfit2$sigma.coefficients)/(exp(fpfit2$sigma.coefficients)+1)
    #
    # } else stop('provide cluster to be TRUE or FALSE')

    RISK = function(x) (predict(fpfit,newdata=data.frame(xx=c(exp(x))), data=datf, type = "response")-
                          predict(fpfit,newdata=data.frame(xx=c(0.00000000000001)), data=datf, type = "response"))/
      (1 - (predict(fpfit,newdata=data.frame(xx=c(0.00000000000001)), data=datf, type = "response")))-q
    bmd.svh = try(uniroot(RISK, interval=c(-5, 0))$root,silent=TRUE)
    bmd.sv <- ifelse((mode(bmd.svh)=="numeric"),(exp(bmd.svh)+0.5)/2,0.05) ## why (exp(bmd.svh)+0.5)/2 ?


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
      }
      if(!is.na(prior.BMD[1])){
        BMD.min = prior.BMD[1]/maxDose
      }
      if(!is.na(prior.BMD[3])){
        BMD.max = prior.BMD[3]/maxDose
      }

      is_informative_BMD = 1

    }else {
      message("Default prior choices used on BMD")
    }

    BMD.vec <- c(BMD.min, BMD.mode, BMD.max)
    shape.BMD <- shape.BMD

  }

  nlevels_BMD <- ifelse(covariate == 'BMD_d' | covariate == 'all', nlevels, 1)
  dim(bmd.sv) <- nlevels_BMD


  ### Parameter d
  if(covariate == 'BMD_d' | covariate == 'all'){

    # Default (normal) priors on k, d, Pert on a
    # prvar.d=1; prmean.d=1;
    # prmean.dQE4=0
    if(prior.d == 'N11'){
      prvar.d = rep(1, nlevels); prmean.d = rep(1, nlevels); truncd = 5
    }else if(prior.d == 'EPA'){
      # prvar.d = 0.5^2; prmean.d = 0.4; truncd = 10000
      prvar.d = rep(0.5, nlevels); prmean.d = rep(0.4, nlevels); truncd = 10000
      # prvar.d = 1; prmean.d = 1; truncd = 10000
    }else if(prior.d == 'N05'){
      prvar.d = rep(0.25, nlevels); prmean.d = rep(0.5, nlevels); truncd = 10000
    }
    # prvar.d=sqrt(0.5); prmean.d = prmean.d
    prmean.dQE4 = rep(0, nlevels); prvar.dQE4 = rep(1, nlevels); truncdQ = 10000


  }else{

    # Default (normal) priors on k, d, Pert on a
    # prvar.d=1; prmean.d=1;
    # prmean.dQE4=0
    if(prior.d == 'N11'){
      prvar.d = 1; prmean.d = 1; truncd = 5
    }else if(prior.d == 'EPA'){
      # prvar.d = 0.5^2; prmean.d = 0.4; truncd = 10000
      prvar.d = 0.5; prmean.d = 0.4; truncd = 10000
      # prvar.d = 1; prmean.d = 1; truncd = 10000
    }else if(prior.d == 'N05'){
      prvar.d = 0.25; prmean.d = 0.5; truncd = 10000
    }
    # prvar.d=sqrt(0.5); prmean.d = prmean.d
    prmean.dQE4 = 0; prvar.dQE4 = 1; truncdQ = 10000

  }

  nlevels_d <- ifelse(covariate == 'BMD_d' | covariate == 'all', nlevels, 1)
  nlevels_b <- ifelse(covariate == 'background' | covariate == 'BMD_d' | covariate == 'all', nlevels, 1)

  priormu1a <- rbind(ifelse(rep(is.vector(a.vec), nlevels), rep(a.vec[2], nlevels), a.vec[ ,2]),
                     ifelse(rep(is.vector(BMD.vec), nlevels), rep(BMD.vec[2], nlevels), BMD.vec[ ,2]),
                     ifelse(rep(length(prmean.d)>1, nlevels), prmean.d, rep(prmean.d, nlevels)),
                     rep(0, nlevels)) # for rho in case of beta_bin

  priormu1bQ <- rbind(ifelse(rep(is.vector(a.vec), nlevels), rep(a.vec[2], nlevels), a.vec[ ,2]),
                      ifelse(rep(is.vector(BMD.vec), nlevels), rep(BMD.vec[2], nlevels), BMD.vec[ ,2]),
                      ifelse(rep(length(prmean.dQE4)>1, nlevels), prmean.dQE4, rep(prmean.dQE4, nlevels)),
                      rep(0, nlevels))

  priorSigma1a <- diag(c(1, 1, prvar.d[1]))
  priorSigma1bQ <- diag(c(1, 1, prvar.dQE4[1]))

  priorlb1a <- rbind(ifelse(rep(is.vector(a.vec), nlevels), rep(a.vec[1], nlevels), a.vec[ ,1]),
                     ifelse(rep(is.vector(BMD.vec), nlevels), rep(BMD.vec[1], nlevels), BMD.vec[ ,1]),
                     rep(0, nlevels),
                     rep(0, nlevels))

  priorub1a <- rbind(ifelse(rep(is.vector(a.vec), nlevels), rep(a.vec[3], nlevels), a.vec[ ,3]),
                     ifelse(rep(is.vector(BMD.vec), nlevels), rep(BMD.vec[3], nlevels), BMD.vec[ ,3]),
                     rep(0, nlevels),
                     rep(0, nlevels))

  priorgama <- rbind(shape.a, shape.BMD)

  row.names(priormu1a) <- row.names(priormu1bQ) <- row.names(priorlb1a) <- c('a', 'BMD', 'd', 'rho')

  if(!is.null(prior.BMD) & !is.vector(BMD.vec)){
    bmd.sv = BMD.vec[,2]
  }  else if(!is.null(prior.BMD) & is.vector(BMD.vec)) {bmd.sv = BMD.vec[2]}

  sv.a <- if(class(a.vec)[1] == "matrix"){a.vec[,2]}else{a.vec[2]}


  if(is_bin==1) {
    start = list(par1 = sv.a, par2 = bmd.sv, par3 = prmean.d)
    startQ = list(par1 = sv.a, par2 = bmd.sv, par3 = prmean.dQE4)

  }


  # else {
  #   rho <- rhohat; dim(rho) <- 1
  #   start = list(par1 = priormu1a[1], par2 = bmd.sv, par3 = priormu1a[3],
  #                rho = rho )
  #   startQ = list(par1 = priormu1a[1], par2 = bmd.sv, par3 = priormu1bQ[3],
  #                rho = rho )
  # }

  trt_ind <- matrix(NA, nrow = N, ncol = nlevels)
  for(i in 1:nlevels) {
    trt_ind[,i] <- as.numeric(covar==covar_lvls[[i]])
  }

  ret.list <- list(
    # data and priors
    data = list(N = N, n = n.a, x = dose.a, y = y.a, yint = y.a, nint = n.a, maxD = maxDose, covariate = covar_lvls,
                q = q, priormu = priormu1a, priormuQ = priormu1bQ, priorSigmaQ=priorSigma1bQ,
                truncd = truncd, truncdQ = truncdQ,
                nlevels = nlevels, nlevels_a = nlevels_a, nlevels_b = nlevels_b, nlevels_BMD = nlevels_BMD, nlevels_d = nlevels_d,
                priorlb = priorlb1a, priorub = priorub1a, priorSigma = priorSigma1a,
                eps = 1.0E-06, priorgama = priorgama, trt_ind = trt_ind,
                init_b = 1, is_informative_a = is_informative_a, is_informative_BMD = is_informative_BMD,
                is_bin = is_bin, is_betabin = is_betabin),
    # start values
    start = start,
    startQ = startQ
  )

  dim(ret.list$start$par2) <- nlevels_BMD
  dim(ret.list$start$par1) <- nlevels_a
  dim(ret.list$startQ$par2) <- nlevels_BMD
  dim(ret.list$startQ$par1) <- nlevels_a
  dim(ret.list$startQ$par3) <- nlevels_d
  dim(ret.list$start$par3) <- nlevels_d

  return(ret.list)

}

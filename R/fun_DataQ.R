
#' Function to set data in the correct format
#'
#' This function also generates appropriate start values and uninformative priors for each model
#' Input should be given as arithmetic mean and standard deviation on the original scale
#'
#' @param data dataframe with input data, order of columns should be: dose, response.
#' @param sumstats number of observations per dose level
#' @param q the specified BMR
#' @param bkg vector containing informative prior for the background.
#'            It should be specified as minimum, most likely and maximum. Defaults to NULL.
#' @param prior.BMD vector containing informative prior for the BMD.
#'                  It should be specified as minimum, most likely and maximum. Defaults to NULL.
#' @param shape.a shape parameter that determines the flatness of the Pert prior for the background.
#'                 Defaults to 4, which implies a peak at the most likely value.
#' @param shape.BMD shape parameter that determines the flatness of the Pert prior for the BMD.
#'                  Defaults to 0.0001, which implies a flat prior.
#' @param cluster logical variable to indicate if data is clustered. TRUE = clustered data. Defaults to FALSE
#'
#' @examples
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
                         # prmean.d = 1, prmean.dQE4 = 0
                         prior.d = c('N11', 'EPA'),
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
  } else if(sumstats == FALSE){
    doses = data[, 1]
    maxDose = max(doses)
    dose.a = sort(unique(doses))
    N = length(dose.a)
    y.a = rep(NA,N)
    n.a = rep(NA,N)
    ybin = data[, 2]
    for (iu in (1:N)){
      y.a[iu] = mean(ybin[doses == dose.a[iu]])
      n.a[iu] = sum(doses == dose.a[iu])
    }
    dose.a = dose.a/maxDose
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
    (1 - (predict(fpfit,newdata=data.frame(xx=c(0.00000000000001)), data=datf, type = "response")))-q
  bmd.svh = try(uniroot(RISK, interval=c(-5, 0))$root,silent=TRUE)
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
  # prvar.d=1; prmean.d=1;
  # prmean.dQE4=0
  if(prior.d == 'N11'){
    prvar.d = 1; prmean.d = 1; truncd = 5
  }else if(prior.d == 'EPA'){
    # prvar.d = 0.5^2; prmean.d = 0.4; truncd = 10000
    prvar.d = 0.5; prmean.d = 0.4; truncd = 10000
    # prvar.d = 1; prmean.d = 1; truncd = 10000
  }
  # prvar.d=sqrt(0.5); prmean.d = prmean.d
  prmean.dQE4 = 0; prvar.dQE4 = 1; truncdQ = 10000

  #prvar.k=1; prmean.k=1
  # family 1a
  priormu1a <- c(a.vec[2], BMD.vec[2], prmean.d, ifelse(is_betabin==1, rhohat, 0))
  priormu1bQ <- c(a.vec[2], BMD.vec[2], prmean.dQE4, ifelse(is_betabin==1, rhohat, 0))
  priorSigma1a <- diag(c(1, 1, prvar.d))
  priorSigma1bQ <- diag(c(1, 1, prvar.dQE4))

  priorlb1a <- c(a.vec[1], BMD.vec[1])
  priorub1a <- c(a.vec[3], BMD.vec[3])

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

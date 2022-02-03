#' Function to set data in the correct format
#'
#' This function also generates appropriate start values and uninformative priors for each model
#'
#' @param data a dataframe with input data, order of columns should be: dose, response, sd (or se), n
#' @param sumstats logical indicating whether summary (T, default) or individual-level (F) data is provided
#' @param q specified BMR
#' @param shape.a shape parameter for the modified PERT distribution on parameter a, defaults to 4, a value of 0.0001 implies a uniform distribution
#' @param shape.BMD shape parameter for the modified PERT distribution on parameter BMD, defaults to 0.0001 (implies a uniform distribution)
#' @param cluster logical indicating if data is clustered (defaults to FALSE)
#' @param bkg vector containing minimum, most likely, and maximum value for the background response
#' @param BMD_p vector containing minimum, most likely, and maximum value for the BMD
#'
#' @return List with data and start values in correct format to be directly used within the BMA functions.
#'
#' @export
#'
PREP_DATA_QA <- function(data, # a dataframe with input data, order of columns should be: dose, response, n
                         sumstats = TRUE, # TRUE if summary data, FALSE if individual data
                         q, # the BMR,
                         shape.a = 4, #scale parameter for Pert priors for a
                         shape.BMD = 0.0001, #scale parameter for the Pert priors for BMD
                         cluster = FALSE, # indicate if data is clustered
                         bkg = NULL, # possible expert info on background response (a),
                         BMD_p = NULL # possible expert prior on the BMD,
){

  if(sumstats == TRUE){
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
  a.min <- ifelse(y.a[1] != 0, max(c(prop.test(y.a[1], n.a[1])$conf.int[1]/2, 1/(10*n.a[1]))),
                  .Machine$double.xmin)
  a.max <- min(c(3*prop.test(y.a[1], n.a[1])$conf.int[2]/2, 1 - 1/(10*n.a[1])))
  a.mode <-  max(c(y.a[1]/n.a[1], 1/(5*n.a[1])))

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
  BMD.max <- 1
  BMD.mode <- 0.5

  ## If info on BMD is given
  if(!is.null(BMD_p)){

    if(!is.na(BMD_p[2])){
      BMD.mode = BMD_p[2]/maxDose
    }
    if(!is.na(BMD_p[1])){
      BMD.min = BMD_p[1]/maxDose
    }
    if(!is.na(BMD_p[3])){
      BMD.max = BMD_p[3]/maxDose
    }

    is_informative_BMD = 1

  }else {
    message("Default prior choices used on BMD")
  }

  BMD.vec <- c(BMD.min, BMD.mode, BMD.max)

  # Default (normal) priors on k, d, Pert on a
  prvar.d=1; prmean.d=0;

  #prvar.k=1; prmean.k=1
  # family 1a
  priormu1a <- c(a.vec[2], BMD.vec[2], prmean.d)
  priorSigma1a <- diag(c(1, 1, prvar.d))
  priorlb1a <- c(a.vec[1], BMD.vec[1])
  priorub1a <- c(a.vec[3], BMD.vec[3])

  if(is_bin==1) {
    start = list(par1 = priormu1a[1], par2 = bmd.sv, par3 = priormu1a[3])
  } else if(is_betabin == 1) {
    etarho <- 0.6; dim(etarho) <- 1
    start = list(par1 = priormu1a[1], par2 = bmd.sv, par3 = priormu1a[3],
                 etarho = etarho )
  }


  return(list(
    # data and priors
    data = list(N = N, n = n.a, x = dose.a, y = y.a, maxD = maxDose, q = q, priormu = priormu1a,
                priorlb = priorlb1a, priorub = priorub1a, priorSigma = priorSigma1a,
                eps = .Machine$double.xmin, priorgama = c(shape.a, shape.BMD),
                init_b = 1, is_informative_a = is_informative_a, is_informative_BMD = is_informative_BMD,
                is_bin = is_bin, is_betabin = is_betabin),
    # start values
    start = start
  ))

}

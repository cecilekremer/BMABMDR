#' Function to set data in the correct format
#'
#' This function also generates appropriate start values and uninformative priors for each model
#'
#' @param data a dataframe with input data, order of columns should be: dose, response, sd (or se), n
#' @param sumstats logical indicating whether summary (T, default) or individual-level (F) data is provided
#' @param geom.stats logicial indicating whether, if summary data are provided, these are geometric (T) or arithmetic (F, default) summary statistics
#' @param sd logical indicating whether standard deviation (T, default) or standard error (F) is provided
#' @param q specified BMR
#' @param prior which prior distribution will be used, defaults to "PERT" (and currently this is the only one implemented)
#' @param bkg vector containing minimum, most likely, and maximum value for the background response
#' @param maxy vector containing minimum, most likely, and maximum value for the response at dose infinity
#' @param prior.BMD vector containing minimum, most likely, and maximum value for the BMD
#' @param shape.a shape parameter for the modified PERT distribution on parameter a, defaults to 4, a value of 0.0001 implies a uniform distribution
#' @param shape.c shape parameter for the modified PERT distribution on parameter c, defaults to 4, a value of 0.0001 implies a uniform distribution
#' @param shape.BMD shape parameter for the modified PERT distribution on parameter BMD, defaults to 4, a value of 0.0001 implies a uniform distribution
#'
#' @return List with data and start values in correct format to be directly used within the BMA functions.
#'
#' @export
#'
PREP_DATA_LN <- function(data, # a dataframe with input data, order of columns should be: dose, response, sd, n
                         sumstats = TRUE, # TRUE if summary data, FALSE if individual data
                         geom.stats = FALSE, # TRUE if geom summ stats are provided
                         sd = TRUE, # TRUE if sd per dose group is given, FALSE is se is given
                         q, # the BMR
                         prior = "PERT",
                         bkg = NULL, maxy = NULL, prior.BMD = NULL, # possible expert info on background and max response
                         shape.a = 4, shape.c = 4, shape.BMD = 4 # shape for the PERT distribution
){

  if(sumstats == TRUE & geom.stats == FALSE){

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
    shift=F
    gmean.a2 = log(NtoLN(mean.a,sd.a))[1:N]
    if (min(gmean.a2)<0) {gmean.a = gmean.a2-1.5*min(gmean.a2); shift=T}
    if (min(gmean.a2)>=0) gmean.a = gmean.a2
    gsd.a = log(NtoLN(mean.a,sd.a))[(N+1):(2*N)]

  }else if(sumstats == TRUE & geom.stats == TRUE){


  }else if(sumstats == FALSE){

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
    shift=F
    gmean.a2 = log(mean.a)
    if (min(gmean.a2)<0) {gmean.a = gmean.a2-1.5*min(gmean.a2); shift=T}
    if (min(gmean.a2)>=0) gmean.a = gmean.a2
    gsd.a = log(sd.a)
  }

  # if(gmean.a[N] < gmean.a[1]) stop("Dose-response curve is not increasing!")

  #####---------------######
  ## For INCREASING

  if(gmean.a[N] > gmean.a[1]){

    # start values
    datf=data.frame(yy=exp(gmean.a),xx=dose.a+0.00000000000001) # exp(mi) = geometric mean
    fpfit=gamlss(yy~fp(xx),family=LOGNO(),data=datf)
    RISK=function(x) (predict(fpfit,newdata=data.frame(xx=c(exp(x))))-predict(fpfit,newdata=data.frame(xx=c(0.00000000000001))))/
      (predict(fpfit,newdata=data.frame(xx=c(0.00000000000001))))-q
    bmd.svh=try(uniroot(RISK, interval=c(-5, 0))$root,silent=T)
    bmd.sv=ifelse((mode(bmd.svh)=="numeric"),bmd.svh,log(0.05))

    ######################
    ### PRIORS

    if(!is.null(bkg)) bkg.ln = log(bkg)
    if(!is.null(maxy)) maxy.ln = log(maxy)

    if(prior == "PERT"){

      ### Family 1

      obs.min = gmean.a[1]
      obs.max = gmean.a[N]

      # min.min = log(0.5*exp(obs.min))
      min.min = 0.001
      mode.min = obs.min
      max.min = log(2*exp(obs.min))

      min.max = gmean.a[1]
      max.max = log(2*exp(gmean.a[N]))
      mode.max = obs.max

      # if(flat2(dose.a, mean.a)$decision == F  & is.null(maxy)){
      if(flat(dose.a, mean.a,inc=T) == F & is.null(maxy)){
        mode.max = log(3*exp(gmean.a[N]))
        # min.max = log(0.5*exp(mode.max))
        # min.max = gmean.a[1]
        min.max  = log(0.5*exp(gmean.a[N]))
        max.max = log(6*exp(gmean.a[N]))

        # warning(
        #   "The data do not contain information on the asymptote, and the default prior for fold change has been set to 3 times the observed maximum.
        #       Please provide prior input on the maximum response by specifying 'maxy', if available."
        # )
      }

      ## If info on background is given
      if(exists("bkg.ln")){
        if(!is.null(bkg.ln)){

          if(!is.na(bkg.ln[2])){
            mode.min = bkg.ln[2]
          }
          if(!is.na(bkg.ln[1])){
            min.min = bkg.ln[1]
          }
          if(!is.na(bkg.ln[3])){
            max.min = bkg.ln[3]
          }

        }
      }

      ## If info on max response is given
      if(exists("maxy.ln")){
        if(!is.null(maxy.ln)){

          if(!is.na(maxy.ln[2])){
            mode.max = maxy.ln[2]
          }
          if(!is.na(maxy.ln[1])){
            min.max = maxy.ln[1]
          }
          if(!is.na(maxy.ln[3])){
            max.max = maxy.ln[3]
          }

        }
      }

      if(!is.null(prior.BMD)){
        mode.BMD = prior.BMD[2]/maxDose
        min.BMD = prior.BMD[1]/maxDose
        if(min.BMD == 0) min.BMD = 0.0001
        max.BMD = prior.BMD[3]/maxDose
      }

      # Default priors on k, d, sigma
      # prvar.k=1; prmean.k=1
      if(!is.null(prior.BMD)){
        BMD.vec = log(c(min.BMD, mode.BMD, max.BMD))
      }else{
        BMD.vec = c(0.01,1,2) # to avoid NA values in shape
      }
      prvar.d=1; prmean.d=0
      # prvar.d=(exp(sqrt(0.18)))^2; prmean.d=2
      prvar.s=1; prmean.s=-2*log(1.5*mean(sd.a))


      # Default prior on a  a.vec = c(min.min, mode.min, max.min)

      a.vec = c(min.min, mode.min, max.min)

      a.vec = log(a.vec)

      # Default prior on c

      c.vec = c(min.max/mode.min, mode.max/mode.min, max.max/mode.min)
      if(c.vec[1] < 1.001) c.vec[1] = 1.001

      c.vec = log(c.vec - 1)

      # family 1a
      priormu1a=c(a.vec[2],BMD.vec[2],c.vec[2],prmean.d,prmean.s)
      priorSigma1a=diag(c(1,1,1,prvar.d,prvar.s))
      priorlb1a = c(a.vec[1],BMD.vec[1],c.vec[1],0,0)
      priorub1a = c(a.vec[3],BMD.vec[3],c.vec[3],0,0)

      # family 1b Gamma
      priormu1bG=c(a.vec[2],BMD.vec[2],c.vec[2],prmean.d,prmean.s)
      priorSigma1bG=diag(c(1,1,1,prvar.d,prvar.s))
      priorlb1bG = c(a.vec[1],BMD.vec[1],c.vec[1],0,0)
      priorub1bG = c(a.vec[3],BMD.vec[3],c.vec[3],0,0)

      # family 1b QE
      priormu1bQ=c(a.vec[2],BMD.vec[2],c.vec[2],prmean.d,prmean.s)
      priorSigma1bQ=diag(c(1,1,1,prvar.d,prvar.s))
      priorlb1bQ = c(a.vec[1],BMD.vec[1],c.vec[1],0,0)
      priorub1bQ = c(a.vec[3],BMD.vec[3],c.vec[3],0,0)

      # family 1a LN
      c.vec = c(min.max/mode.min, mode.max/mode.min, max.max/mode.min)
      if(c.vec[1] < (1+(log(1+q)/gmean.a[1]))) c.vec[1] = 1.001 + (log(1+q)/gmean.a[1])
      c.vec = log(c.vec-(1+(log(1+q)/gmean.a[1])))


      # if(c.vec[1] <= log(log(1+q)/gmean.a[1]) ) c.vec[1] = log(log(1+q)/gmean.a[1]) + 0.0001
      priormu1aLN=c(a.vec[2],BMD.vec[2],c.vec[2],prmean.d,prmean.s)
      priorSigma1aLN=diag(c(1,1,1,prvar.d,prvar.s))
      priorlb1aLN = c(a.vec[1],BMD.vec[1],c.vec[1],0,0)
      priorub1aLN = c(a.vec[3],BMD.vec[3],c.vec[3],0,0)

      # }

      ### Family 2

      # Default priors on k, d, sigma
      # prvar.k=1; prmean.k=1
      prvar.d=1; prmean.d=0
      # prvar.d=(exp(sqrt(0.18)))^2; prmean.d=2
      prvar.s=1; prmean.s=-2*log(1.5*mean(sd.a))

      ## Default range on background and max response

      # min.min = log(0.5*exp(obs.min))
      min.min = 0.001
      mode.min = obs.min
      max.min = log(2*exp(obs.min))

      min.max = gmean.a[1]
      mode.max = obs.max
      max.max = log(2*exp(gmean.a[N]))

      # if(flat2(dose.a, mean.a)$decision == F  & is.null(maxy)) {
      if(flat(dose.a, mean.a,inc=T) == F & is.null(maxy)){
        # min.max = gmean.a[N]
        min.max  = log(0.5*exp(gmean.a[N]))
        mode.max = log(3*exp(gmean.a[N]))
        # max.max = 10*obs.max
        max.max = log(6*exp(gmean.a[N]))

      }

      if(exists("bkg.ln")){
        if(!is.null(bkg.ln)){

          if(!is.na(bkg.ln[2])){
            mode.min = bkg.ln[2]
          }
          if(!is.na(bkg.ln[1])){
            min.min = bkg.ln[1]
          }
          if(!is.na(bkg.ln[3])){
            max.min = bkg.ln[3]
          }

        }
      }

      ## If info on max response is given
      if(exists("maxy.ln")){
        if(!is.null(maxy.ln)){

          if(!is.na(maxy.ln[2])){
            mode.max = maxy.ln[2]
          }
          if(!is.na(maxy.ln[1])){
            min.max = maxy.ln[1]
          }
          if(!is.na(maxy.ln[3])){
            max.max = maxy.ln[3]
          }

        }
      }

      # Prior on a

      a.vec = c(min.max, mode.max, max.max)
      a.vec = log(a.vec)

      # Prior on c

      if(max.min/mode.max >= 1){
        c.vec.P = c(qnorm(min.min/mode.max), qnorm(mode.min/mode.max),qnorm(max.min/(1.001*max.min)))
      }else{
        c.vec.P = c(qnorm(min.min/mode.max), qnorm(mode.min/mode.max), qnorm(max.min/mode.max))
      }
      if(c.vec.P[3] >= qnorm(1 - (log(1+q)/gmean.a[1]))) c.vec.P[3] = qnorm(1 - (log(1+q)/gmean.a[1])) - 0.0001
      c.vec.P = log(qnorm(1-(log(1+q)/gmean.a[1])) - c.vec.P)
      if(c.vec.P[1] > c.vec.P[3]){
        cc.vec.P <- c()
        cc.vec.P[2] = c.vec.P[2]
        cc.vec.P[1] = c.vec.P[3]; cc.vec.P[3] = c.vec.P[1]
        c.vec.P = cc.vec.P
      }

      if(max.min/mode.max >= 1){
        c.vec.L = c(logit(min.min/mode.max), logit(mode.min/mode.max),logit(max.min/(1.001*max.min)))
      }else{
        c.vec.L = c(logit(min.min/mode.max), logit(mode.min/mode.max), logit(max.min/mode.max))
      }
      if(c.vec.L[3] >= logit(1 - (log(1+q)/gmean.a[1]))) c.vec.L[3] = logit(1 - (log(1+q)/gmean.a[1])) - 0.0001
      c.vec.L = log(logit(1-(log(1+q)/gmean.a[1])) - c.vec.L)
      if(c.vec.L[1] > c.vec.L[3]){
        cc.vec.L <- c()
        cc.vec.L[2] = c.vec.L[2]
        cc.vec.L[1] = c.vec.L[3]; cc.vec.L[3] = c.vec.L[1]
        c.vec.L = cc.vec.L
      }


      # family 2P
      priormu2P=c(a.vec[2],BMD.vec[2],c.vec.P[2],prmean.d,prmean.s)
      priorSigma2P=diag(c(1,1,1,prvar.d,prvar.s))
      priorlb2P = c(a.vec[1],BMD.vec[1],c.vec.P[1],0,0)
      priorub2P = c(a.vec[3],BMD.vec[3],c.vec.P[3],0,0)

      # family 2L
      priormu2L=c(a.vec[2],BMD.vec[2],c.vec.L[2],prmean.d,prmean.s)
      priorSigma2L=diag(c(1,1,1,prvar.d,prvar.s))
      priorlb2L = c(a.vec[1],BMD.vec[1],c.vec.L[1],0,0)
      priorub2L = c(a.vec[3],BMD.vec[3],c.vec.L[3],0,0)

      if(!is.null(prior.BMD)) bmd.sv = BMD.vec[2]

      ret.list <- list(
        # data and priors
        data.modstan1a=list(N=N,n=n.a,x=dose.a,m=gmean.a,m.org=gmean.a2,shift=shift,s2=gsd.a^2,maxD=maxDose,q=q,priormu=priormu1a,
                            shape1 = c(fun.alpha(a = priorlb1a[1], b = priormu1a[1], c = priorub1a[1], g = shape.a),
                                       fun.alpha(a = priorlb1a[2], b = priormu1a[2], c = priorub1a[2], g = shape.BMD),
                                       fun.alpha(a = priorlb1a[3], b = priormu1a[3], c = priorub1a[3], g = shape.c), 0, 0),
                            shape2 = c(fun.beta(a = priorlb1a[1], b = priormu1a[1], c = priorub1a[1], g = shape.a),
                                       fun.beta(a = priorlb1a[2], b = priormu1a[2], c = priorub1a[2], g = shape.BMD),
                                       fun.beta(a = priorlb1a[3], b = priormu1a[3], c = priorub1a[3], g = shape.c), 0, 0),
                            priorlb = priorlb1a, priorub = priorub1a, shape.a = shape.a, shape.c = shape.c, shape.BMD = shape.BMD,
                            priorSigma=priorSigma1a),

        data.modstan1bG=list(N=N,n=n.a,x=dose.a,m=gmean.a,m.org=gmean.a2,shift=shift,s2=gsd.a^2,maxD=maxDose,q=q,priormu=priormu1bG,
                             shape1 = c(fun.alpha(a = priorlb1bG[1], b = priormu1bG[1], c = priorub1bG[1], g = shape.a),
                                        fun.alpha(a = priorlb1bG[2], b = priormu1bG[2], c = priorub1bG[2], g = shape.BMD),
                                        fun.alpha(a = priorlb1bG[3], b = priormu1bG[3], c = priorub1bG[3], g = shape.c), 0, 0),
                             shape2 = c(fun.beta(a = priorlb1bG[1], b = priormu1bG[1], c = priorub1bG[1], g = shape.a),
                                        fun.beta(a = priorlb1bG[2], b = priormu1bG[2], c = priorub1bG[2], g = shape.BMD),
                                        fun.beta(a = priorlb1bG[3], b = priormu1bG[3], c = priorub1bG[3], g = shape.c), 0, 0),
                             priorlb = priorlb1bG, priorub = priorub1bG, shape.a = shape.a, shape.c = shape.c, shape.BMD = shape.BMD,
                             priorSigma=priorSigma1bG,
                             init_b=1),

        data.modstan1bQ=list(N=N,n=n.a,x=dose.a,m=gmean.a,m.org=gmean.a2,shift=shift,s2=gsd.a^2,maxD=maxDose,q=q,priormu=priormu1bQ,
                             shape1 = c(fun.alpha(a = priorlb1bQ[1], b = priormu1bQ[1], c = priorub1bQ[1], g = shape.a),
                                        fun.alpha(a = priorlb1bQ[2], b = priormu1bQ[2], c = priorub1bQ[2], g = shape.BMD),
                                        fun.alpha(a = priorlb1bQ[3], b = priormu1bQ[3], c = priorub1bQ[3], g = shape.c), 0, 0),
                             shape2 = c(fun.beta(a = priorlb1bQ[1], b = priormu1bQ[1], c = priorub1bQ[1], g = shape.a),
                                        fun.beta(a = priorlb1bQ[2], b = priormu1bQ[2], c = priorub1bQ[2], g = shape.BMD),
                                        fun.beta(a = priorlb1bQ[3], b = priormu1bQ[3], c = priorub1bQ[3], g = shape.c), 0, 0),
                             priorlb = priorlb1bQ, priorub = priorub1bQ, shape.a = shape.a, shape.c = shape.c, shape.BMD = shape.BMD,
                             priorSigma=priorSigma1bQ),

        data.modstan1aLN=list(N=N,n=n.a,x=dose.a,m=gmean.a,m.org=gmean.a2,shift=shift,s2=gsd.a^2,maxD=maxDose,q=q,priormu=priormu1aLN,
                              shape1 = c(fun.alpha(a = priorlb1aLN[1], b = priormu1aLN[1], c = priorub1aLN[1], g = shape.a),
                                         fun.alpha(a = priorlb1aLN[2], b = priormu1aLN[2], c = priorub1aLN[2], g = shape.BMD),
                                         fun.alpha(a = priorlb1aLN[3], b = priormu1aLN[3], c = priorub1aLN[3], g = shape.c), 0, 0),
                              shape2 = c(fun.beta(a = priorlb1aLN[1], b = priormu1aLN[1], c = priorub1aLN[1], g = shape.a),
                                         fun.beta(a = priorlb1aLN[2], b = priormu1aLN[2], c = priorub1aLN[2], g = shape.BMD),
                                         fun.beta(a = priorlb1aLN[3], b = priormu1aLN[3], c = priorub1aLN[3], g = shape.c), 0, 0),
                              priorlb = priorlb1aLN, priorub = priorub1aLN, shape.a = shape.a, shape.c = shape.c, shape.BMD = shape.BMD,
                              priorSigma=priorSigma1aLN),

        data.modstan2P=list(N=N,n=n.a,x=dose.a,m=gmean.a,m.org=gmean.a2,shift=shift,s2=gsd.a^2,maxD=maxDose,q=q,priormu=priormu2P,
                            shape1 = c(fun.alpha(a = priorlb2P[1], b = priormu2P[1], c = priorub2P[1], g = shape.a),
                                       fun.alpha(a = priorlb2P[2], b = priormu2P[2], c = priorub2P[2], g = shape.BMD),
                                       fun.alpha(a = priorlb2P[3], b = priormu2P[3], c = priorub2P[3], g = shape.c), 0, 0),
                            shape2 = c(fun.beta(a = priorlb2P[1], b = priormu2P[1], c = priorub2P[1], g = shape.a),
                                       fun.beta(a = priorlb2P[2], b = priormu2P[2], c = priorub2P[2], g = shape.BMD),
                                       fun.beta(a = priorlb2P[3], b = priormu2P[3], c = priorub2P[3], g = shape.c), 0, 0),
                            priorlb = priorlb2P, priorub = priorub2P, shape.a = shape.a, shape.c = shape.c, shape.BMD = shape.BMD,
                            priorSigma=priorSigma2P),

        data.modstan2L=list(N=N,n=n.a,x=dose.a,m=gmean.a,m.org=gmean.a2,shift=shift,s2=gsd.a^2,maxD=maxDose,q=q,priormu=priormu2L,
                            shape1 = c(fun.alpha(a = priorlb2L[1], b = priormu2L[1], c = priorub2L[1], g = shape.a),
                                       fun.alpha(a = priorlb2L[2], b = priormu2L[2], c = priorub2L[2], g = shape.BMD),
                                       fun.alpha(a = priorlb2L[3], b = priormu2L[3], c = priorub2L[3], g = shape.c), 0, 0),
                            shape2 = c(fun.beta(a = priorlb2L[1], b = priormu2L[1], c = priorub2L[1], g = shape.a),
                                       fun.beta(a = priorlb2L[2], b = priormu2L[2], c = priorub2L[2], g = shape.BMD),
                                       fun.beta(a = priorlb2L[3], b = priormu2L[3], c = priorub2L[3], g = shape.c), 0, 0),
                            priorlb = priorlb2L, priorub = priorub2L, shape.a = shape.a, shape.c = shape.c, shape.BMD = shape.BMD,
                            priorSigma=priorSigma2L),
        # start values
        svF1=list(par=c(priormu1a[1],bmd.sv, priormu1a[3], 0,log(1/mean(gsd.a^2)))),
        svF1LN=list(par=c(priormu1aLN[1],bmd.sv, priormu1aLN[3],0,log(1/mean(gsd.a^2)))),
        svF2P=list(par=c(priormu2P[1], bmd.sv, priormu2P[3], 0, log(1/mean(gsd.a^2)))), # probit
        svF2L=list(par=c(priormu2L[1], bmd.sv, priormu2L[3], 0, log(1/mean(gsd.a^2)))),
        # ,
        # DR effect
        # print(DR.effect)
        increasing = TRUE
      )

    }else if(prior == "Normal"){

      # priors

      # family 1a
      prvar.uninf.a=1; prmean.uninf.a=log(mean(gmean.a))
      prvar.uninf.k=1; prmean.uninf.k=1
      prvar.uninf.c=1; prmean.uninf.c=log(round(gmean.a[N]/gmean.a[1],1)-1)
      if(flat(dose.a,mean.a,inc=T)==F) prmean.uninf.c=log(3*round(gmean.a[N]/gmean.a[1],1)-1)
      prvar.uninf.d=1; prmean.uninf.d=0
      prvar.uninf.s=1; prmean.uninf.s=-2*log(1.5*mean(gsd.a))
      priormu1a=c(prmean.uninf.a,prmean.uninf.k,prmean.uninf.c,prmean.uninf.d,prmean.uninf.s)
      priorSigma1a=diag(c(prvar.uninf.a,prvar.uninf.k,prvar.uninf.c,prvar.uninf.d,prvar.uninf.s))

      # family 1a LN
      prvar.uninf.c=1; prmean.uninf.c=log(round(gmean.a[N]/gmean.a[1],1)-(1+(log(1+q)/mean(gmean.a))))
      if(flat(dose.a,mean.a,inc=T)==F) prmean.uninf.c=log(3*round(gmean.a[N]/gmean.a[1],1)-(1+(log(1+q)/mean(gmean.a))))
      priormu1aLN=c(prmean.uninf.a,prmean.uninf.k,prmean.uninf.c,prmean.uninf.d,prmean.uninf.s)
      priorSigma1aLN=diag(c(prvar.uninf.a,prvar.uninf.k,prvar.uninf.c,prvar.uninf.d,prvar.uninf.s))

      # family 1b Gamma
      prvar.uninf.a=1; prmean.uninf.a=log(mean(gmean.a))
      prvar.uninf.k=1; prmean.uninf.k=1
      prvar.uninf.c=1; prmean.uninf.c=log(1.5*round(gmean.a[N]/gmean.a[1],1)-1)
      if(flat(dose.a,mean.a,inc=T)==F) prmean.uninf.c=log(3*round(gmean.a[N]/gmean.a[1],1)-1)
      prvar.uninf.d=1; prmean.uninf.d=0
      prvar.uninf.s=1; prmean.uninf.s=-2*log(1.5*mean(gsd.a))
      priormu1bG=c(prmean.uninf.a,prmean.uninf.k,prmean.uninf.c,prmean.uninf.d,prmean.uninf.s)
      priorSigma1bG=diag(c(prvar.uninf.a,prvar.uninf.k,prvar.uninf.c,prvar.uninf.d,prvar.uninf.s))
      # Prior for E4 to use as start values for G4
      # prvar.uninf.d=0.001; prmean.uninf.d=0
      # priormuE=c(prmean.uninf.a,prmean.uninf.k,prmean.uninf.c,prmean.uninf.d,prmean.uninf.s)
      # priorSigmaE=diag(c(prvar.uninf.a,prvar.uninf.k,prvar.uninf.c,prvar.uninf.d,prvar.uninf.s))

      # family 1b QE
      prvar.uninf.a=1; prmean.uninf.a=log(mean(gmean.a))
      prvar.uninf.k=1; prmean.uninf.k=1
      prvar.uninf.c=1; prmean.uninf.c=log(1.5*round(gmean.a[N]/gmean.a[1],1)-1)
      if(flat(dose.a,mean.a,inc=T)==F) prmean.uninf.c=log(3*round(gmean.a[N]/gmean.a[1],1)-1)
      prvar.uninf.d=1; prmean.uninf.d=0
      prvar.uninf.s=1; prmean.uninf.s=-2*log(1.5*mean(gsd.a))
      priormu1bQ=c(prmean.uninf.a,prmean.uninf.k,prmean.uninf.c,prmean.uninf.d,prmean.uninf.s)
      priorSigma1bQ=diag(c(prvar.uninf.a,prvar.uninf.k,prvar.uninf.c,prvar.uninf.d,prvar.uninf.s))

      # family 2P
      prvar.uninf.a=0.1; prmean.uninf.a=log(gmean.a[N])
      if(flat(dose.a,mean.a,inc=T)==F) prmean.uninf.a=log(3*gmean.a[N])
      prvar.uninf.k=1; prmean.uninf.k=1
      # prvar.uninf.c=1; prmean.uninf.c=0.5*round(qnorm(rnd(gmean.a[1]/gmean.a[N])),1)
      prvar.uninf.c = 1; prmean.uninf.c = log(qnorm(1-(log(1+q)/gmean.a[N])) - qnorm(gmean.a[1]/gmean.a[N]))
      prvar.uninf.d=1; prmean.uninf.d=0
      prvar.uninf.s=1; prmean.uninf.s=-2*log(1.5*mean(gsd.a))
      priormu2P=c(prmean.uninf.a,prmean.uninf.k,prmean.uninf.c,prmean.uninf.d,prmean.uninf.s)
      priorSigma2P=diag(c(prvar.uninf.a,prvar.uninf.k,prvar.uninf.c,prvar.uninf.d,prvar.uninf.s))

      # family 2L
      prvar.uninf.c = 1; prmean.uninf.c = log(logit(1-(log(1+q)/gmean.a[N])) - logit(gmean.a[1]/gmean.a[N]))
      priormu2L=c(prmean.uninf.a,prmean.uninf.k,prmean.uninf.c,prmean.uninf.d,prmean.uninf.s)
      priorSigma2L=diag(c(prvar.uninf.a,prvar.uninf.k,prvar.uninf.c,prvar.uninf.d,prvar.uninf.s))

      ## Data in correct format
      ret.list <- list(
        # data and priors
        data.modstan1a=list(N=N,n=n.a,x=dose.a,m=gmean.a,m.org=gmean.a2,shift=shift,s2=gsd.a^2,maxD=maxDose,q=q,priormu=priormu1a,priorSigma=priorSigma1a),
        data.modstan1bG=list(N=N,n=n.a,x=dose.a,m=gmean.a,m.org=gmean.a2,shift=shift,s2=gsd.a^2,maxD=maxDose,q=q,priormu=priormu1bG,priorSigma=priorSigma1bG,init_b=1),
        data.modstan1bQ=list(N=N,n=n.a,x=dose.a,m=gmean.a,m.org=gmean.a2,shift=shift,s2=gsd.a^2,maxD=maxDose,q=q,priormu=priormu1bQ,priorSigma=priorSigma1bQ),
        data.modstan1aLN=list(N=N,n=n.a,x=dose.a,m=gmean.a,m.org=gmean.a2,shift=shift,s2=gsd.a^2,maxD=maxDose,q=q,priormu=priormu1aLN,priorSigma=priorSigma1aLN),
        # data.modstanE=list(N=N,n=n.a,x=dose.a,m=gmean.a,m.org=gmean.a2,shift=shift,s2=gsd.a^2,maxD=maxDose,q=q,priormu=priormuE,priorSigma=priorSigmaE),
        data.modstan2P=list(N=N,n=n.a,x=dose.a,m=gmean.a,m.org=gmean.a2,shift=shift,s2=gsd.a^2,maxD=maxDose,q=q,priormu=priormu2P,priorSigma=priorSigma2P),
        data.modstan2L=list(N=N,n=n.a,x=dose.a,m=gmean.a,m.org=gmean.a2,shift=shift,s2=gsd.a^2,maxD=maxDose,q=q,priormu=priormu2L,priorSigma=priorSigma2L),
        # start values
        svF1=list(par=c(log(gmean.a[1]),bmd.sv, log(gmean.a[N]/gmean.a[1]-1),0,log(1/mean(gsd.a^2)))),
        svF1LN=list(par=c(log(gmean.a[1]),bmd.sv, log(gmean.a[N]/gmean.a[1]- -(1+(log(1+q)/mean(gmean.a)))),0,log(1/mean(gsd.a^2)))),
        # svF2P=list(par=c(log(gmean.a[N]), bmd.sv, qnorm(gmean.a[1]/gmean.a[N]), 0, log(1/mean(gsd.a^2)))), # probit
        # svF2L=list(par=c(log(gmean.a[N]), bmd.sv, logit(gmean.a[1]/gmean.a[N]), 0, log(1/mean(gsd.a^2))))
        svF2P=list(par=c(log(gmean.a[N]), bmd.sv, log(qnorm(1-(log(1+q)/gmean.a[N])) - qnorm(gmean.a[1]/gmean.a[N])), 0, log(1/mean(gsd.a^2)))), # probit
        svF2L=list(par=c(log(gmean.a[N]), bmd.sv, log(logit(1-(log(1+q)/gmean.a[N])) - logit(gmean.a[1]/gmean.a[N])), 0, log(1/mean(gsd.a^2))))
        ,
        # DR effect
        # print(DR.effect)
        increasing = TRUE
      )


    }

    # test for dose-response effect
    # DR.effect = anydoseresponseLNI(dose.a,mean.a,sd.a,n.a)
    # DR.effect = anydoseresponseLNI(dose.a,mean.a,gsd.a,n.a)

    # data in correct format
    return(ret.list)

    ##########------------############
    ## For DECREASING

  }else if(gmean.a[N] < gmean.a[1]){

    # start value BMD
    datf=data.frame(yy=gmean.a,xx=dose.a+0.00000000000001)
    fpfit=gamlss(yy~fp(xx),family=NO(),data=datf)
    RISK=function(x) (predict(fpfit,newdata=data.frame(xx=c(exp(x))))-predict(fpfit,newdata=data.frame(xx=c(0.00000000000001))))/
      (predict(fpfit,newdata=data.frame(xx=c(0.00000000000001))))+q
    bmd.svh=try(uniroot(RISK, interval=c(-5, 0))$root,silent=T)
    bmd.sv=ifelse((mode(bmd.svh)=="numeric"),bmd.svh,log(0.05))

    ######################
    ### PRIORS

    if(!is.null(bkg)) bkg.ln = log(bkg)
    if(!is.null(maxy)) maxy.ln = log(maxy)

    ## Observed background and maximum response

    obs.min = gmean.a[1]
    obs.max = gmean.a[N]

    if(prior == "PERT"){

      ### Family 1

      ## Default range on background and max response

      min.min = log(0.5*exp(obs.min))
      mode.min = obs.min
      max.min = log(2*exp(obs.min))

      min.max = log(0.5*exp(obs.max))
      max.max = obs.min
      mode.max = obs.max

      if(flat(dose.a, mean.a,inc=F) == F & is.null(maxy)){
        mode.max = log(0.5*exp(obs.max))
        min.max  = 0.001
        max.max = obs.min

        # warning(
        #   "The data do not contain information on the asymptote, and the default prior for fold change has been set to 3 times the observed maximum.
        #       Please provide prior input on the maximum response by specifying 'maxy', if available."
        # )
      }


      ## If info on background is given
      if(exists("bkg.ln")){
        if(!is.null(bkg.ln)){

          if(!is.na(bkg.ln[2])){
            mode.min = bkg.ln[2]
          }
          if(!is.na(bkg.ln[1])){
            min.min = bkg.ln[1]
          }
          if(!is.na(bkg.ln[3])){
            max.min = bkg.ln[3]
          }

        }
      }

      ## If info on max response is given
      if(exists("maxy.ln")){
        if(!is.null(maxy.ln)){

          if(!is.na(maxy.ln[2])){
            mode.max = maxy.ln[2]
          }
          if(!is.na(maxy.ln[1])){
            min.max = maxy.ln[1]
          }
          if(!is.na(maxy.ln[3])){
            max.max = maxy.ln[3]
          }

        }
      }

      if(!is.null(prior.BMD)){
        mode.BMD = prior.BMD[2]/maxDose
        min.BMD = prior.BMD[1]/maxDose
        if(min.BMD == 0) min.BMD = 0.0001
        max.BMD = prior.BMD[3]/maxDose
      }

      # Default (normal) priors on k, d, sigma
      # prvar.k=1; prmean.k=1
      if(!is.null(prior.BMD)){
        BMD.vec = log(c(min.BMD, mode.BMD, max.BMD))
      }else{
        BMD.vec = c(0.01,1,2) # to avoid NA values in shape
      }
      prvar.d=1; prmean.d=0
      prvar.s=1; prmean.s=-2*log(1.5*mean(sd.a))


      # Prior on a
      a.vec = c(min.min, mode.min, max.min)
      a.vec = log(a.vec)

      # Prior on c
      c.vec = c(min.max/mode.min, mode.max/mode.min, max.max/mode.min)
      if(c.vec[3] >= 1) c.vec[3] = 0.999
      c.vec = log(c.vec/(1-c.vec))

      # family 1a
      priormu1a=c(a.vec[2],BMD.vec[2],c.vec[2],prmean.d,prmean.s)
      priorSigma1a=diag(c(1,1,1,prvar.d,prvar.s))
      priorlb1a = c(a.vec[1],BMD.vec[1],c.vec[1],0,0)
      priorub1a = c(a.vec[3],BMD.vec[3],c.vec[3],0,0)

      # family 1b Gamma
      priormu1bG=c(a.vec[2],BMD.vec[2],c.vec[2],prmean.d,prmean.s)
      priorSigma1bG=diag(c(1,1,1,prvar.d,prvar.s))
      priorlb1bG = c(a.vec[1],BMD.vec[1],c.vec[1],0,0)
      priorub1bG = c(a.vec[3],BMD.vec[3],c.vec[3],0,0)

      # family 1b QE
      priormu1bQ=c(a.vec[2],BMD.vec[2],c.vec[2],prmean.d,prmean.s)
      priorSigma1bQ=diag(c(1,1,1,prvar.d,prvar.s))
      priorlb1bQ = c(a.vec[1],BMD.vec[1],c.vec[1],0,0)
      priorub1bQ = c(a.vec[3],BMD.vec[3],c.vec[3],0,0)

      # family 1a LN
      c.vec = c(min.max/mode.min, mode.max/mode.min, max.max/mode.min)
      c.vec = c.vec/(1+(log(1-q)/mean.a[1]))
      if(c.vec[3] >= 1) c.vec[3] = 0.999
      c.vec = log(c.vec/(1-c.vec))
      priormu1aLN=c(a.vec[2],BMD.vec[2],c.vec[2],prmean.d,prmean.s)
      priorSigma1aLN=diag(c(1,1,1,prvar.d,prvar.s))
      priorlb1aLN = c(a.vec[1],BMD.vec[1],c.vec[1],0,0)
      priorub1aLN = c(a.vec[3],BMD.vec[3],c.vec[3],0,0)

      ### Family 2

      # Default priors on k, d, sigma
      # prvar.k=1; prmean.k=1
      prvar.d=1; prmean.d=0
      prvar.s=1; prmean.s=-2*log(1.5*mean(sd.a))

      ## Default range on background and max response

      min.min = log(0.5*exp(obs.min))
      mode.min = obs.min
      max.min = log(2*exp(obs.min))

      min.max = log(0.5*exp(obs.max))
      max.max = obs.min
      mode.max = obs.max

      if(flat(dose.a, mean.a,inc=F) == F & is.null(maxy)){
        mode.max = log(0.5*exp(obs.max))
        min.max  = 0.001
        max.max = obs.min

        # warning(
        #   "The data do not contain information on the asymptote, and the default prior for fold change has been set to 3 times the observed maximum.
        #       Please provide prior input on the maximum response by specifying 'maxy', if available."
        # )
      }

      ## If info on background is given
      if(exists("bkg.ln")){

        if(!is.null(bkg.ln)){

          if(!is.na(bkg.ln[2])){
            mode.min = bkg.ln[2]
          }
          if(!is.na(bkg.ln[1])){
            min.min = bkg.ln[1]
          }
          if(!is.na(bkg.ln[3])){
            max.min = bkg.ln[3]
          }

        }
      }

      ## If info on max response is given
      if(exists("maxy.ln")){

        if(!is.null(maxy.ln)){

          if(!is.na(maxy.ln[2])){
            mode.max = maxy.ln[2]
          }
          if(!is.na(maxy.ln[1])){
            min.max = maxy.ln[1]
          }
          if(!is.na(maxy.ln[3])){
            max.max = maxy.ln[3]
          }

        }
      }

      # Prior on a
      a.vec = c(min.min, mode.min, max.min)
      a.vec = log(a.vec)

      # Prior on c
      if(max.max/mode.min >= 1){
        c.vec.P = c(qnorm(min.max/mode.min), qnorm(mode.max/mode.min), qnorm(max.max/(1.001*max.max)))
      }else{
        c.vec.P = c(qnorm(min.max/mode.min), qnorm(mode.max/mode.min), qnorm(max.max/mode.min))
      }
      if(c.vec.P[3] >= qnorm(1 + (log(1-q)/gmean.a[1]))) c.vec.P[3] = qnorm(1 + (log(1-q)/gmean.a[1])) - 0.0001
      c.vec.P = log(qnorm(1+(log(1-q)/gmean.a[1])) - c.vec.P)
      if(c.vec.P[1] > c.vec.P[3]){
        cc.vec.P <- c()
        cc.vec.P[2] = c.vec.P[2]
        cc.vec.P[1] = c.vec.P[3]; cc.vec.P[3] = c.vec.P[1]
        c.vec.P = cc.vec.P
      }

      if(max.max/mode.min >= 1){
        c.vec.L = c(logit(min.max/mode.min), logit(mode.max/mode.min), logit(max.max/(1.001*max.max)))
      }else{
        c.vec.L = c(logit(min.max/mode.min), logit(mode.max/mode.min), logit(max.max/mode.min))
      }
      if(c.vec.L[3] >= logit(1 + (log(1-q)/gmean.a[1]))) c.vec.L[3] = logit(1 + (log(1-q)/gmean.a[1])) - 0.0001
      c.vec.L = log(logit(1+(log(1-q)/gmean.a[1])) - c.vec.L)
      if(c.vec.L[1] > c.vec.L[3]){
        cc.vec.L <- c()
        cc.vec.L[2] = c.vec.L[2]
        cc.vec.L[1] = c.vec.L[3]; cc.vec.L[3] = c.vec.L[1]
        c.vec.L = cc.vec.L
      }

      # family 2P
      priormu2P=c(a.vec[2],BMD.vec[2],c.vec.P[2],prmean.d,prmean.s)
      priorSigma2P=diag(c(1,1,1,prvar.d,prvar.s))
      priorlb2P = c(a.vec[1],BMD.vec[1],c.vec.P[1],0,0)
      priorub2P = c(a.vec[3],BMD.vec[3],c.vec.P[3],0,0)

      # family 2L
      priormu2L=c(a.vec[2],BMD.vec[2],c.vec.L[2],prmean.d,prmean.s)
      priorSigma2L=diag(c(1,1,1,prvar.d,prvar.s))
      priorlb2L = c(a.vec[1],BMD.vec[1],c.vec.L[1],0,0)
      priorub2L = c(a.vec[3],BMD.vec[3],c.vec.L[3],0,0)

      if(!is.null(prior.BMD)) bmd.sv = BMD.vec[2]

      ret.list <- list(
        # data and priors
        data.modstan1a=list(N=N,n=n.a,x=dose.a,m=gmean.a,m.org=gmean.a2,shift=shift,s2=gsd.a^2,maxD=maxDose,q=q,priormu=priormu1a,
                            shape1 = c(fun.alpha(a = priorlb1a[1], b = priormu1a[1], c = priorub1a[1], g = shape.a),
                                       fun.alpha(a = priorlb1a[2], b = priormu1a[2], c = priorub1a[2], g = shape.BMD),
                                       fun.alpha(a = priorlb1a[3], b = priormu1a[3], c = priorub1a[3], g = shape.c), 0, 0),
                            shape2 = c(fun.beta(a = priorlb1a[1], b = priormu1a[1], c = priorub1a[1], g = shape.a),
                                       fun.beta(a = priorlb1a[2], b = priormu1a[2], c = priorub1a[2], g = shape.BMD),
                                       fun.beta(a = priorlb1a[3], b = priormu1a[3], c = priorub1a[3], g = shape.c), 0, 0),
                            priorlb = priorlb1a, priorub = priorub1a, shape.a = shape.a, shape.c = shape.c, shape.BMD = shape.BMD,
                            priorSigma=priorSigma1a),

        data.modstan1bG=list(N=N,n=n.a,x=dose.a,m=gmean.a,m.org=gmean.a2,shift=shift,s2=gsd.a^2,maxD=maxDose,q=q,priormu=priormu1bG,
                             shape1 = c(fun.alpha(a = priorlb1bG[1], b = priormu1bG[1], c = priorub1bG[1], g = shape.a),
                                        fun.alpha(a = priorlb1bG[2], b = priormu1bG[2], c = priorub1bG[2], g = shape.BMD),
                                        fun.alpha(a = priorlb1bG[3], b = priormu1bG[3], c = priorub1bG[3], g = shape.c), 0, 0),
                             shape2 = c(fun.beta(a = priorlb1bG[1], b = priormu1bG[1], c = priorub1bG[1], g = shape.a),
                                        fun.beta(a = priorlb1bG[2], b = priormu1bG[2], c = priorub1bG[2], g = shape.BMD),
                                        fun.beta(a = priorlb1bG[3], b = priormu1bG[3], c = priorub1bG[3], g = shape.c), 0, 0),
                             priorlb = priorlb1bG, priorub = priorub1bG, shape.a = shape.a, shape.c = shape.c, shape.BMD = shape.BMD,
                             priorSigma=priorSigma1bG,
                             init_b=1),

        data.modstan1bQ=list(N=N,n=n.a,x=dose.a,m=gmean.a,m.org=gmean.a2,shift=shift,s2=gsd.a^2,maxD=maxDose,q=q,priormu=priormu1bQ,
                             shape1 = c(fun.alpha(a = priorlb1bQ[1], b = priormu1bQ[1], c = priorub1bQ[1], g = shape.a),
                                        fun.alpha(a = priorlb1bQ[2], b = priormu1bQ[2], c = priorub1bQ[2], g = shape.BMD),
                                        fun.alpha(a = priorlb1bQ[3], b = priormu1bQ[3], c = priorub1bQ[3], g = shape.c), 0, 0),
                             shape2 = c(fun.beta(a = priorlb1bQ[1], b = priormu1bQ[1], c = priorub1bQ[1], g = shape.a),
                                        fun.beta(a = priorlb1bQ[2], b = priormu1bQ[2], c = priorub1bQ[2], g = shape.BMD),
                                        fun.beta(a = priorlb1bQ[3], b = priormu1bQ[3], c = priorub1bQ[3], g = shape.c), 0, 0),
                             priorlb = priorlb1bQ, priorub = priorub1bQ, shape.a = shape.a, shape.c = shape.c, shape.BMD = shape.BMD,
                             priorSigma=priorSigma1bQ),

        data.modstan1aLN=list(N=N,n=n.a,x=dose.a,m=gmean.a,m.org=gmean.a2,shift=shift,s2=gsd.a^2,maxD=maxDose,q=q,priormu=priormu1aLN,
                              shape1 = c(fun.alpha(a = priorlb1aLN[1], b = priormu1aLN[1], c = priorub1aLN[1], g = shape.a),
                                         fun.alpha(a = priorlb1aLN[2], b = priormu1aLN[2], c = priorub1aLN[2], g = shape.BMD),
                                         fun.alpha(a = priorlb1aLN[3], b = priormu1aLN[3], c = priorub1aLN[3], g = shape.c), 0, 0),
                              shape2 = c(fun.beta(a = priorlb1aLN[1], b = priormu1aLN[1], c = priorub1aLN[1], g = shape.a),
                                         fun.beta(a = priorlb1aLN[2], b = priormu1aLN[2], c = priorub1aLN[2], g = shape.BMD),
                                         fun.beta(a = priorlb1aLN[3], b = priormu1aLN[3], c = priorub1aLN[3], g = shape.c), 0, 0),
                              priorlb = priorlb1aLN, priorub = priorub1aLN, shape.a = shape.a, shape.c = shape.c, shape.BMD = shape.BMD,
                              priorSigma=priorSigma1aLN),

        data.modstan2P=list(N=N,n=n.a,x=dose.a,m=gmean.a,m.org=gmean.a2,shift=shift,s2=gsd.a^2,maxD=maxDose,q=q,priormu=priormu2P,
                            shape1 = c(fun.alpha(a = priorlb2P[1], b = priormu2P[1], c = priorub2P[1], g = shape.a),
                                       fun.alpha(a = priorlb2P[2], b = priormu2P[2], c = priorub2P[2], g = shape.BMD),
                                       fun.alpha(a = priorlb2P[3], b = priormu2P[3], c = priorub2P[3], g = shape.c), 0, 0),
                            shape2 = c(fun.beta(a = priorlb2P[1], b = priormu2P[1], c = priorub2P[1], g = shape.a),
                                       fun.beta(a = priorlb2P[2], b = priormu2P[2], c = priorub2P[2], g = shape.BMD),
                                       fun.beta(a = priorlb2P[3], b = priormu2P[3], c = priorub2P[3], g = shape.c), 0, 0),
                            priorlb = priorlb2P, priorub = priorub2P, shape.a = shape.a, shape.c = shape.c, shape.BMD = shape.BMD,
                            priorSigma=priorSigma2P),

        data.modstan2L=list(N=N,n=n.a,x=dose.a,m=gmean.a,m.org=gmean.a2,shift=shift,s2=gsd.a^2,maxD=maxDose,q=q,priormu=priormu2L,
                            shape1 = c(fun.alpha(a = priorlb2L[1], b = priormu2L[1], c = priorub2L[1], g = shape.a),
                                       fun.alpha(a = priorlb2L[2], b = priormu2L[2], c = priorub2L[2], g = shape.BMD),
                                       fun.alpha(a = priorlb2L[3], b = priormu2L[3], c = priorub2L[3], g = shape.c), 0, 0),
                            shape2 = c(fun.beta(a = priorlb2L[1], b = priormu2L[1], c = priorub2L[1], g = shape.a),
                                       fun.beta(a = priorlb2L[2], b = priormu2L[2], c = priorub2L[2], g = shape.BMD),
                                       fun.beta(a = priorlb2L[3], b = priormu2L[3], c = priorub2L[3], g = shape.c), 0, 0),
                            priorlb = priorlb2L, priorub = priorub2L, shape.a = shape.a, shape.c = shape.c, shape.BMD = shape.BMD,
                            priorSigma=priorSigma2L),
        # start values
        svF1=list(par=c(priormu1a[1],bmd.sv, priormu1a[3], 1,log(1/mean(gsd.a^2)))),
        svF1LN=list(par=c(priormu1aLN[1],bmd.sv, priormu1aLN[3],1,log(1/mean(gsd.a^2)))),
        svF2P=list(par=c(priormu2P[1], bmd.sv, priormu2P[3], 0, log(1/mean(gsd.a^2)))), # probit
        svF2L=list(par=c(priormu2L[1], bmd.sv, priormu2L[3], 0, log(1/mean(gsd.a^2)))),
        # ,
        # DR effect
        # print(DR.effect)
        increasing = FALSE
      )


    }else if(prior == "Normal"){

      # priors

      # family 1a
      prvar.uninf.a=1; prmean.uninf.a=log(mean(gmean.a))
      prvar.uninf.k=1; prmean.uninf.k=1
      prvar.uninf.c=1; prmean.uninf.c=log(round(gmean.a[N]/(gmean.a[1]-gmean.a[N]),2))
      if(flat(dose.a,mean.a,inc=F)==F) prmean.uninf.c=log((1/3)*round(gmean.a[N]/(gmean.a[1]-gmean.a[N]),2))
      prvar.uninf.d=1; prmean.uninf.d=0
      prvar.uninf.s=1; prmean.uninf.s=-2*log(1.5*mean(gsd.a))
      priormu1a=c(prmean.uninf.a,prmean.uninf.k,prmean.uninf.c,prmean.uninf.d,prmean.uninf.s)
      priorSigma1a=diag(c(prvar.uninf.a,prvar.uninf.k,prvar.uninf.c,prvar.uninf.d,prvar.uninf.s))

      # family 1a LN
      prvar.uninf.c=1; prmean.uninf.c=logit(round(gmean.a[N]/(gmean.a[1]),2)/(1+(log(1-q)/mean(gmean.a))))/(sqrt(3)/pi)
      if(flat(dose.a,mean.a,inc=F)==F) prmean.uninf.c=logit((1/3)*round(gmean.a[N]/(gmean.a[1]),2)/(1+(log(1-q)/mean(gmean.a))))/(sqrt(3)/pi)
      priormu1aLN=c(prmean.uninf.a,prmean.uninf.k,prmean.uninf.c,prmean.uninf.d,prmean.uninf.s)
      priorSigma1aLN=diag(c(prvar.uninf.a,prvar.uninf.k,prvar.uninf.c,prvar.uninf.d,prvar.uninf.s))

      # family 1b Gamma
      prvar.uninf.a=1; prmean.uninf.a=log(mean(gmean.a))
      # prvar.uninf.c=1; prmean.uninf.c=logit(round(gmean.a[N]/(gmean.a[1]),2)/(1+(log(1-q)/mean(gmean.a))))/(sqrt(3)/pi)
      # if(flat(dose.a,mean.a,inc=F)==F) prmean.uninf.c=logit((1/3)*round(gmean.a[N]/(gmean.a[1]),2)/(1+(log(1-q)/mean(gmean.a))))/(sqrt(3)/pi)
      prvar.uninf.c=1; prmean.uninf.c=log(1.5*round(gmean.a[N]/(gmean.a[1]-gmean.a[N]),2))
      if(flat(dose.a,mean.a,inc=F)==F) prmean.uninf.c=log(0.5*round(gmean.a[N]/(gmean.a[1]-gmean.a[N]),2))
      prvar.uninf.d=1; prmean.uninf.d=0
      prvar.uninf.s=1; prmean.uninf.s=-2*log(1.5*mean(gsd.a))
      priormu1bG=c(prmean.uninf.a,prmean.uninf.k,prmean.uninf.c,prmean.uninf.d,prmean.uninf.s)
      priorSigma1bG=diag(c(prvar.uninf.a,prvar.uninf.k,prvar.uninf.c,prvar.uninf.d,prvar.uninf.s))
      # Prior for E4 to use as start values for G4
      # prvar.uninf.d=0.001; prmean.uninf.d=0
      # priormuE=c(prmean.uninf.a,prmean.uninf.k,prmean.uninf.c,prmean.uninf.d,prmean.uninf.s)
      # priorSigmaE=diag(c(prvar.uninf.a,prvar.uninf.k,prvar.uninf.c,prvar.uninf.d,prvar.uninf.s))

      # family 1b QE
      prvar.uninf.a=1; prmean.uninf.a=log(mean(gmean.a))
      # prvar.uninf.c=1; prmean.uninf.c=log(round(gmean.a[N]/(gmean.a[1]-gmean.a[N]),2))
      # if(flat(dose.a,mean.a)==F) prmean.uninf.c=logit(1/3*round(gmean.a[N]/gmean.a[1],2))
      # if(flat(dose.a,mean.a,inc=F)==F) prmean.uninf.c=log((1/3)*round(gmean.a[N]/(gmean.a[1]-gmean.a[N]),2))
      prvar.uninf.s=1; prmean.uninf.s=-2*log(1.5*mean(gsd.a))
      priormu1bQ=c(prmean.uninf.a,prmean.uninf.k,prmean.uninf.c,prmean.uninf.d,prmean.uninf.s)
      priorSigma1bQ=diag(c(prvar.uninf.a,prvar.uninf.k,prvar.uninf.c,prvar.uninf.d,prvar.uninf.s))

      # family 2P
      prvar.uninf.a=0.1; prmean.uninf.a=log(gmean.a[1])
      prvar.uninf.c=1; prmean.uninf.c=log((qnorm(1+(log(1-q)/gmean.a[1])) - round(qnorm(rnd(gmean.a[N]/gmean.a[1])),2)))
      if(flat(dose.a,mean.a,inc=F)==F) prmean.uninf.c=log(0.5*(qnorm(1+(log(1-q)/gmean.a[1])) - round(qnorm(rnd(gmean.a[N]/gmean.a[1])),2)))
      prvar.uninf.d=1; prmean.uninf.d=0
      prvar.uninf.s=1; prmean.uninf.s=-2*log(1.5*mean(gsd.a))
      priormu2P=c(prmean.uninf.a,prmean.uninf.k,prmean.uninf.c,prmean.uninf.d,prmean.uninf.s)
      priorSigma2P=diag(c(prvar.uninf.a,prvar.uninf.k,prvar.uninf.c,prvar.uninf.d,prvar.uninf.s))

      # family 2L
      prvar.uninf.c=1; prmean.uninf.c=log((logit(1+(log(1-q)/gmean.a[1])) - round(logit(rnd(gmean.a[N]/gmean.a[1])),2)))
      if(flat(dose.a,mean.a,inc=F)==F) prmean.uninf.c=log(0.5*(logit(1+(log(1-q)/gmean.a[1])) - round(logit(rnd(gmean.a[N]/gmean.a[1])),2)))
      priormu2L=c(prmean.uninf.a,prmean.uninf.k,prmean.uninf.c,prmean.uninf.d,prmean.uninf.s)
      priorSigma2L=diag(c(prvar.uninf.a,prvar.uninf.k,prvar.uninf.c,prvar.uninf.d,prvar.uninf.s))

      ## data in correct format

      ret.list <- list(
        # data and priors
        data.modstan1a=list(N=N,n=n.a,x=dose.a,m=gmean.a,m.org=gmean.a2,shift=shift,s2=gsd.a^2,maxD=maxDose,q=q,priormu=priormu1a,priorSigma=priorSigma1a),
        data.modstan1bG=list(N=N,n=n.a,x=dose.a,m=gmean.a,m.org=gmean.a2,shift=shift,s2=gsd.a^2,maxD=maxDose,q=q,priormu=priormu1bG,priorSigma=priorSigma1bG,init_b=1),
        data.modstan1bQ=list(N=N,n=n.a,x=dose.a,m=gmean.a,m.org=gmean.a2,shift=shift,s2=gsd.a^2,maxD=maxDose,q=q,priormu=priormu1bQ,priorSigma=priorSigma1bQ),
        data.modstan1aLN=list(N=N,n=n.a,x=dose.a,m=gmean.a,m.org=gmean.a2,shift=shift,s2=gsd.a^2,maxD=maxDose,q=q,priormu=priormu1aLN,priorSigma=priorSigma1aLN),
        # data.modstanE=list(N=N,n=n.a,x=dose.a,m=gmean.a,m.org=gmean.a2,shift=shift,s2=gsd.a^2,maxD=maxDose,q=q,priormu=priormuE,priorSigma=priorSigmaE),
        data.modstan2P=list(N=N,n=n.a,x=dose.a,m=gmean.a,m.org=gmean.a2,shift=shift,s2=gsd.a^2,maxD=maxDose,q=q,priormu=priormu2P,priorSigma=priorSigma2P),
        data.modstan2L=list(N=N,n=n.a,x=dose.a,m=gmean.a,m.org=gmean.a2,shift=shift,s2=gsd.a^2,maxD=maxDose,q=q,priormu=priormu2L,priorSigma=priorSigma2L),
        # start values
        svF1=list(par=c(log(gmean.a[1]),bmd.sv, log(gmean.a[N]/(gmean.a[1]-gmean.a[N])),0,log(1/mean(gsd.a^2)))),
        svF1LN=list(par=c(log(gmean.a[1]),bmd.sv, logit((gmean.a[N]/(gmean.a[1]))/(1+(log(1-q)/mean(gmean.a))))/(sqrt(3)/pi),0,log(1/mean(gsd.a^2)))),
        svF2P=list(par=c(log(gmean.a[1]), bmd.sv, log(qnorm(1+(log(1-q)/gmean.a[1])) - qnorm(min(mean.a)/max(mean.a))), 0, log(1/mean(gsd.a^2)))), # probit
        svF2L=list(par=c(log(gmean.a[1]), bmd.sv, log(logit(1+(log(1-q)/gmean.a[1])) - logit(min(mean.a)/max(mean.a))), 0, log(1/mean(gsd.a^2)))),
        # DR effect
        # print(DR.effect)
        increasing = FALSE
      )

    }



    # test for dose-response effect
    # DR.effect = anydoseresponseNI(dose.a,mean.a,sd.a,n.a)


    # data in correct format
    return(ret.list)


  }

}

#' Function for internal use
#'
#' @param x value
#'
#' @return .
#'
#' @export expit
#'
expit=function(x) 1/(1+exp(-pi*x/sqrt(3)))

#' Function for internal use
#'
#' @param x value
#'
#' @return .
#'
#' @export logit
#'
logit=function(x) sqrt(3)/pi*log(x/(1-x))

#' Function for internal use
#'
#' @param p value
#'
#' @return .
#'
#' @export rnd
#'
rnd=function(p) (p<0.5)*round(p,1)+(p>=0.5)*floor(p*10)/10

#' Function for internal use
#'
#' @param dose value
#' @param mean value
#' @param n value
#' @param inc logical variable to indicate if the dose-resonse curve is increasing or decreasing
#' @return .logical value indicating if the dose-response curve is flat or not
#'
#' @export flat
#'
flat = function(dose,mean,n,inc){ # To determine if DR curve flattens or not
  flat=F
  if(inc==TRUE){
    if(length(unique(dose)) == length(dose)){
      dat = data.frame(dose,mean)
      datf=data.frame(yy=mean,xx=dose+0.0001)
      fpfit=gamlss::gamlss(yy~fp(xx),family=NO,data=datf)
      # maxdiff=max(abs(diff(predict(fpfit),lag=1,differences=1))/diff(log(dose+0.0001),lag=1,differences=1))
      # lastdiff=(abs(diff(predict(fpfit),lag=1,differences=1))/diff(log(dose+0.0001),lag=1,differences=1))[length(dose)-1]
      maxdiff=max((diff(predict(fpfit),lag=1,differences=1))/(diff(log(dose+0.0001),lag=1,differences=1)))
      lastdiff=((diff(predict(fpfit),lag=1,differences=1))/(diff(log(dose+0.0001),lag=1,differences=1)))[length(dose)-1]
      if (lastdiff/maxdiff<(0.5)) flat=T # flat if last incremental change smaller than 50% of the maximal change
    }
    else{
      mean.i = c()
      dose.i = c()
      j = 1
      for(i in unique(dose)){
        dose.i[j] = i
        mean.i[j] = weighted.mean(x = mean[dose == i], w = n[dose == i])
        j = j+1
      }
      dat = data.frame(dose.i,mean.i)
      datf=data.frame(yy=mean.i,xx=dose.i+0.0001)
      fpfit=gamlss::gamlss(yy~fp(xx),family=NO,data=datf)
      maxdiff=max((diff(predict(fpfit),lag=1,differences=1))/(diff(log(dose.i+0.0001),lag=1,differences=1)))
      lastdiff=((diff(predict(fpfit),lag=1,differences=1))/(diff(log(dose.i+0.0001),lag=1,differences=1)))[length(dose.i)-1]
      if (lastdiff/maxdiff<(0.5)) flat=T # flat if last incremental change smaller than 50% of the maximal change
    }

    return(flat)

  }else if(inc==FALSE){
    if(length(unique(dose)) == length(dose)){
      dat = data.frame(dose,mean)
      datf=data.frame(yy=mean,xx=dose+0.0001)
      fpfit=gamlss::gamlss(yy~fp(xx),family=NO,data=datf)
      # maxdiff=max(abs(diff(predict(fpfit),lag=1,differences=1))/diff(log(dose+0.0001),lag=1,differences=1))
      # lastdiff=(abs(diff(predict(fpfit),lag=1,differences=1))/diff(log(dose+0.0001),lag=1,differences=1))[length(dose)-1]
      maxdiff=max(abs(diff(predict(fpfit),lag=1,differences=1))/(diff(log(dose+0.0001),lag=1,differences=1)))
      lastdiff=(abs(diff(predict(fpfit),lag=1,differences=1))/(diff(log(dose+0.0001),lag=1,differences=1)))[length(dose)-1]
      if (lastdiff/maxdiff<(0.5)) flat=T # flat if last incremental change smaller than 50% of the maximal change
    }else{
      mean.i = c()
      dose.i = c()
      j = 1
      for(i in unique(dose)){
        dose.i[j] = i
        mean.i[j] = weighted.mean(x = mean[dose == i], w = n[dose == i])
        j = j+1
      }
      dat = data.frame(dose.i,mean.i)
      datf=data.frame(yy=mean.i,xx=dose.i+0.0001)
      fpfit=gamlss::gamlss(yy~fp(xx),family=NO,data=datf)
      maxdiff=max((diff(predict(fpfit),lag=1,differences=1))/(diff(log(dose.i+0.0001),lag=1,differences=1)))
      lastdiff=((diff(predict(fpfit),lag=1,differences=1))/(diff(log(dose.i+0.0001),lag=1,differences=1)))[length(dose.i)-1]
      if (lastdiff/maxdiff<(0.5)) flat=T # flat if last incremental change smaller than 50% of the maximal change
    }

    return(flat)
  }

}
#' Functions to get parameters for PERT prior
#'
#' @param a minimum
#' @param b most likely value
#' @param c maximum
#' @param g shape
#'
#' @return .
#'
#' @export fun.alpha
#'
fun.alpha = function(a,b,c,g){
  # if((b-a) == (c-b)){
  #   c = c + 0.0001
  # }
  # m = (a + g*b + c)/(g + 2)
  # s = ((m - a)*(2*b - a - c))/((b - m)*(c - a))
  s = 1 + (g*((b-a)/(c-a)))
  if(s<0) stop("Specified values lead to negative shape parameter for the PERT distribution;
               try increasing the prior range")
  return(s)
  # return((shape*b + c - 5*a)/(c-a))
}
#' Functions to get parameters for PERT prior
#'
#' @param a minimum
#' @param b most likely value
#' @param c maximum
#' @param g shape
#'
#' @return .
#'
#' @export fun.beta
#'
fun.beta = function(a,b,c,g){
  # if((b-a) == (c-b)){
  #   c = c + 0.0001
  # }
  # m = (a + g*b + c)/(g + 2)
  # alpha = fun.alpha(a,b,c,g)
  # s = (alpha * (c - m))/(m - a)
  s = 1 + (g*((c-b)/(c-a)))
  if(s<0) stop("Specified values lead to negative shape parameter for the PERT distribution; try increasing the prior range")
  return(s)
  # return((5*c - a - shape*b)/(c-a))
}
#' Functions to get parameters for PERT prior
#'
#' @param x value
#' @param lb lower bound
#' @param ub upper bound
#' @param s1 shape
#' @param s2 shape
#'
#' @return .
#'
#' @export fun.alpha
#'
pert_dist = function(x,lb,ub,s1,s2){
  # return(
  # dbeta((x-lb)/(ub-lb), shape1 = s1, shape2 = s2, log = T) - log((ub - lb))
  x1 = (s1-1) * log((x - lb))
  x2 = (s2-1) * log((ub - x))
  x3 = (s1+s2-1) * log((ub - lb))
  x4 = lbeta(s1, s2)
  return( x1 + x2 - x3 - x4)

  # )
}

#' Bartlett function for testing constant variance and co-efficient of variation; H0 equal variances
#'
#' @param sd standard deviation per dose group
#' @param n number of observations per dose group
#'
#' @return vector with test statistic and pvalue for the Bartlett's test
#' @export bartlett
#'
bartlett <- function(sd,n){
  N<-sum(n);k<-length(n)

  sp2<-sum((n-1)*sd^2)/(N-k)
  B<-(
    (N-k)*log(sp2)-sum((n-1)*log(sd^2))
  )/(
    1+(sum(1/(n-1))-1/(N-k))/(3*(k-1))
  )
  c(B,pchisq(B,df=(k-1),lower.tail = F))
}

#' Obtain a list of all available models
#'
#' @param type one of 'continuous' or 'quantal'
#'
#' @return .
#' @export get_models
#'
get_models <- function(type = c('continuous', 'quantal')){
  if(type == 'continuous'){
    mods <- c("E4_N","IE4_N","H4_N","LN4_N","G4_N","QE4_N","P4_N","L4_N","E4_LN","IE4_LN","H4_LN","LN4_LN","G4_LN","QE4_LN","P4_LN","L4_LN")
    names(mods) <- c('Exponential Normal',
                     'Inverse Exponential Normal',
                     'Hill Normal',
                     'Lognormal Normal',
                     'Gamma Normal',
                     'Quadratic Exponential Normal',
                     'Probit Normal',
                     'Logit Normal',
                     'Exponential Lognormal',
                     'Inverse Exponential Lognormal',
                     'Hill Lognormal',
                     'Lognormal Lognormal',
                     'Gamma Lognormal',
                     'Quadratic Exponential Lognormal',
                     'Probit Lognormal',
                     'Logit Lognormal')
  }else if(type == 'quantal'){
    mods <- c("E4_Q","IE4_Q","H4_Q","LN4_Q","G4_Q","QE4_Q","P4_Q","L4_Q")
    names(mods) <- c('Exponential',
                     'Inverse Exponential',
                     'Hill',
                     'Lognormal',
                     'Gamma',
                     'Quadratic Exponential',
                     'Probit',
                     'Logit')
  }
  return(mods)
}



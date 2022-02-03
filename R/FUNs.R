#' Function for internal use
#'
#' @param x value
#'
#' @return .
#'
expit=function(x) 1/(1+exp(-pi*x/sqrt(3)))

#' Function for internal use
#'
#' @param x value
#'
#' @return .
#'
logit=function(x) sqrt(3)/pi*log(x/(1-x))

#' Function for internal use
#'
#' @param p value
#'
#' @return .
#'
rnd=function(p) (p<0.5)*round(p,1)+(p>=0.5)*floor(p*10)/10

#' Function for internal use
#'
#' @param dose value
#' @param mean value
#' @param inc logical
#'
#' @return .
#'
flat = function(dose,mean,inc){ # To determine if DR curve flattens or not
  flat=F
  if(inc==TRUE){
    dat = data.frame(dose,mean)
    datf=data.frame(yy=mean,xx=dose+0.0001)
    fpfit=gamlss::gamlss(yy~fp(xx),family=NO,data=datf)
    # maxdiff=max(abs(diff(predict(fpfit),lag=1,differences=1))/diff(log(dose+0.0001),lag=1,differences=1))
    # lastdiff=(abs(diff(predict(fpfit),lag=1,differences=1))/diff(log(dose+0.0001),lag=1,differences=1))[length(dose)-1]
    maxdiff=max((diff(predict(fpfit),lag=1,differences=1))/diff(log(dose+0.0001),lag=1,differences=1))
    lastdiff=((diff(predict(fpfit),lag=1,differences=1))/diff(log(dose+0.0001),lag=1,differences=1))[length(dose)-1]
    if (lastdiff/maxdiff<(0.5)) flat=T # flat if last incremental change smaller than 50% of the maximal change
    return(flat)
  }else if(inc==FALSE){
    dat = data.frame(dose,mean)
    datf=data.frame(yy=mean,xx=dose+0.0001)
    fpfit=gamlss::gamlss(yy~fp(xx),family=NO,data=datf)
    # maxdiff=max(abs(diff(predict(fpfit),lag=1,differences=1))/diff(log(dose+0.0001),lag=1,differences=1))
    # lastdiff=(abs(diff(predict(fpfit),lag=1,differences=1))/diff(log(dose+0.0001),lag=1,differences=1))[length(dose)-1]
    maxdiff=max(abs(diff(predict(fpfit),lag=1,differences=1))/diff(log(dose+0.0001),lag=1,differences=1))
    lastdiff=(abs(diff(predict(fpfit),lag=1,differences=1))/diff(log(dose+0.0001),lag=1,differences=1))[length(dose)-1]
    if (lastdiff/maxdiff<(0.5)) flat=T # flat if last incremental change smaller than 50% of the maximal change
    return(flat)
  }

}


#' Functions to get parameters for PERT prior
#'
#' @param a minimum
#' @param b most likely value
#' @param c maximum
#' @param g shape of modified PERT
#'
#' @return .
#'
fun.alpha = function(a,b,c,g){
  if((b-a) == (c-b)){
    c = c + 0.0001
  }
  m = (a + g*b + c)/(g + 2)
  s = ((m - a)*(2*b - a - c))/((b - m)*(c - a))
  if(s<0) stop("Specified values lead to negative shape parameter for the PERT distribution; try increasing the prior range")
  return(s)
  # return((shape*b + c - 5*a)/(c-a))
}
#' @rdname fun.alpha
fun.beta = function(a,b,c,g){
  if((b-a) == (c-b)){
    c = c + 0.0001
  }
  m = (a + g*b + c)/(g + 2)
  alpha = fun.alpha(a,b,c,g)
  s = (alpha * (c - m))/(m - a)
  if(s<0) stop("Specified values lead to negative shape parameter for the PERT distribution; try increasing the prior range")
  return(s)
  # return((5*c - a - shape*b)/(c-a))
}
#' @rdname fun.alpha
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
#' @return .
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

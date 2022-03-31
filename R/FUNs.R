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
  return(((m - a)*(2*b - a - c))/((b - m)*(c - a)))
  # return((shape*b + c - 5*a)/(c-a))
}
#' @rdname fun.alpha
fun.beta = function(a,b,c,g){
  if((b-a) == (c-b)){
    c = c + 0.0001
  }
  m = (a + g*b + c)/(g + 2)
  alpha = fun.alpha(a,b,c,g)
  return(
    (alpha * (c - m))/(m - a)
  )
  # return((5*c - a - shape*b)/(c-a))
}
#' @rdname fun.alpha
pert_dist = function(x,lb,ub,s1,s2){
  return(
    dbeta((x-lb)/(ub-lb), shape1 = s1, shape2 = s2, log = T) - log((ub - lb))
  )
}


#' Function to check if the mean for the Dose-Response data is approaching or at its asymptote
#' for internal use
#'
#' @param dose value
#' @param response value
#' @param h step size for the derivative
#' @param B number of bootstraps
#'
#' @return .
#'
#'
flat2 <- function(dose, response, h = 1/1000, B = 200) {

  md_data <- data.frame(dose = (dose/max(dose)) + 1.0e-5,
                        response = response)

  pseudo_dose <- seq(0, 1, by = h) + 1.0e-5

  dd <- matrix(NA, nrow = length(pseudo_dose)-1, ncol = B)

  for(i in 1:B) {
    ss <- sample(1:nrow(md_data), replace = TRUE)

    p_data <- md_data[ss,]

    fpm <- try(gamlss::gamlss(response ~ fp(dose),
                              data = p_data, family = NO,
                              control = gamlss::gamlss.control(trace = FALSE)),
               silent = TRUE)


    if(class(fpm)[1] != 'try-error') {

      ffp_pred2 <- try(predict(fpm, what = 'mu',
                               newdata = data.frame(dose = pseudo_dose),
                               data = p_data), silent = TRUE)
      #compute the diff
      dd[,i] <- diff(ffp_pred2, lag = 1, differences = 1)/h

    } else dd[,i] <- NA

  }

  #compute the quantiles of h
  qdiff <- apply(dd, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
  #check if the quantile range contains 0
  test <- apply(qdiff, 2, function(x) {
    return(x[1] <= 0 & x[2] >= 0)
  })
  flat <- ifelse(sum(test)==ncol(qdiff), T, F)
  return(list(last3_derivative = dd, qts =  qdiff, decision = flat))

}



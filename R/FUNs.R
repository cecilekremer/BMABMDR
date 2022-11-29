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

#' Function to convert arithmetic to geometric mean/standard deviation
#'
#'
#'
#' @param am Arithmetic means per odered dose level, on original scale
#' @param asd Arithmetic standard deviations per odered dose level, on original scale
#'
#' @description Converst arithmetic mean (sd) on original scale to geometric mean (sd) to be used in log-normal distribution
#'
#' @examples
#'  data("immunotoxicityData.rda")
#'  NtoLN(am = immunotoxicityData$Mean, asd = immunotoxicityData$SD)
#'
#' @return Vector containing the geometric means and standard deviations per ordered dose level.
#'
#' @export NtoLN
#'
NtoLN=function(am,asd){
  gm=am/sqrt(asd^2/am^2+1)
  gsd=exp(sqrt(log(asd^2/am^2+1)))
  return(c(gm,gsd))
}

#' Function to convert geometric to arithmetic mean/standard deviation
#'
#'
#' @param gm Geometric means per odered dose level, on original scale
#' @param gsd Geometric standard deviations per odered dose level, on original scale
#'
#' @description Converst geometric mean (sd) on original scale to arithmetic mean (sd) to be used in normal distribution
#'
#' @examples
#'  data("immunotoxicityData.rda")
#'  LNtoN(gm = immunotoxicityData$Mean, gsd = immunotoxicityData$SD)
#'
#' @return Vector containing the arithmetic means and standard deviations per ordered dose level.
#'
#' @export LNtoN
#'
LNtoN=function(gm,gsd){
  am=exp(log(gm)+log(gsd)^2/2)
  asd=sqrt((exp(log(gsd)^2)-1)*exp(2*log(gm)+log(gsd)^2))
  return(c(am,asd))
}


#' Function for internal use
#'
#' @param dose value
#' @param mean value
#' @param n value
#' @param inc logical variable to indicate if the dose-resonse curve is increasing or decreasing
#'
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
      fpfit=try(gamlss::gamlss(yy~fp(xx),family=NO,data=datf), silent = T)
      if(class(fpfit)[1] == 'try-error'){
        flat = F
      }else{
        maxdiff=max((diff(predict(fpfit),lag=1,differences=1))/(diff(log(dose+0.0001),lag=1,differences=1)))
        lastdiff=((diff(predict(fpfit),lag=1,differences=1))/(diff(log(dose+0.0001),lag=1,differences=1)))[length(dose)-1]
        if (lastdiff/maxdiff<(0.5)) flat=T # flat if last incremental change smaller than 50% of the maximal change
      }
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
      fpfit=try(gamlss::gamlss(yy~fp(xx),family=NO,data=datf), silent = T)
      if(class(fpfit)[1] == 'try-error'){
        flat = F
      }else{
        maxdiff=max((diff(predict(fpfit),lag=1,differences=1))/(diff(log(dose.i+0.0001),lag=1,differences=1)))
        lastdiff=((diff(predict(fpfit),lag=1,differences=1))/(diff(log(dose.i+0.0001),lag=1,differences=1)))[length(dose.i)-1]
        if (lastdiff/maxdiff<(0.5)) flat=T # flat if last incremental change smaller than 50% of the maximal change
      }
    }

    return(flat)

  }else if(inc==FALSE){
    if(length(unique(dose)) == length(dose)){
      dat = data.frame(dose,mean)
      datf=data.frame(yy=mean,xx=dose+0.0001)
      fpfit=try(gamlss::gamlss(yy~fp(xx),family=NO,data=datf), silent = T)
      if(class(fpfit)[1] == 'try-error'){
        flat = F
      }else{
        maxdiff=max(abs(diff(predict(fpfit),lag=1,differences=1))/(diff(log(dose+0.0001),lag=1,differences=1)))
        lastdiff=(abs(diff(predict(fpfit),lag=1,differences=1))/(diff(log(dose+0.0001),lag=1,differences=1)))[length(dose)-1]
        if (lastdiff/maxdiff<(0.5)) flat=T # flat if last incremental change smaller than 50% of the maximal change
      }
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
      fpfit=try(gamlss::gamlss(yy~fp(xx),family=NO,data=datf), silent = T)
      if(class(fpfit)[1] == 'try-error'){
        flat = F
      }else{
        maxdiff=max((diff(predict(fpfit),lag=1,differences=1))/(diff(log(dose.i+0.0001),lag=1,differences=1)))
        lastdiff=((diff(predict(fpfit),lag=1,differences=1))/(diff(log(dose.i+0.0001),lag=1,differences=1)))[length(dose.i)-1]
        if (lastdiff/maxdiff<(0.5)) flat=T # flat if last incremental change smaller than 50% of the maximal change
      }
    }

    return(flat)
  }

}

#' Function for internal use
#'
#' @param dose value
#' @param mean value
#' @param inc logical variable to indicate if the dose-resonse curve is increasing or decreasing
#' @return .logical value indicating if the dose-response curve is flat or not
#'
#' @export flatC
#'
flatC = function(dose.a,mean.a,inc){ # To determine if DR curve flattens or not
  flat=F
  if(inc==TRUE){
    # dat = data.frame(dose,mean)
    # datf=data.frame(yy=mean.a,xx=dose.a+0.0001)
    # fpfit=gamlss::gamlss(yy~fp(xx),family=NO,data=datf)
    fpfit=try(gamlss::gamlss(mean.a~fp(dose.a+0.0001),family=NO), silent = T)
    if(class(fpfit)[1] == 'try-error'){
      flat = F
    }else{
      maxdiff=max((diff(predict(fpfit),lag=1,differences=1))/diff(log(dose.a+0.0001),lag=1,differences=1))
      lastdiff=((diff(predict(fpfit),lag=1,differences=1))/diff(log(dose.a+0.0001),lag=1,differences=1))[length(dose.a)-1]
      if (lastdiff/maxdiff<(0.5)) flat=T # flat if last incremental change smaller than 50% of the maximal change
    }
    return(flat)
  }else if(inc==FALSE){
    # dat = data.frame(dose,mean)
    # datf=data.frame(yy=mean.a,xx=dose.a+0.0001)
    # fpfit=gamlss::gamlss(yy~fp(xx),family=NO,data=datf)
    fpfit=try(gamlss::gamlss(mean.a~fp(dose.a+0.0001),family=NO), silent = T)
    if(class(fpfit)[1] == 'try-error'){
      flat = F
    }else{
      maxdiff=max(abs(diff(predict(fpfit),lag=1,differences=1))/diff(log(dose.a+0.0001),lag=1,differences=1))
      lastdiff=(abs(diff(predict(fpfit),lag=1,differences=1))/diff(log(dose.a+0.0001),lag=1,differences=1))[length(dose.a)-1]
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
  if(s<0) stop("Specified values lead to negative shape parameter for the PERT distribution; try increasing the prior range")
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
#' @param lb value
#' @param ub value
#' @param s1 value
#' @param s2 value
#'
#' @return .
#'
#' @export pert_dist
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

#' Function to check if the mean for the Dose-Response data is approaching or at its asymptote
#' for internal use
#'
#' @param dose value
#' @param response value
#' @param h step size for the derivative
#' @param B number of bootstraps
#'
#' @return .
#' @export flat2
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

#' Bartlett function for testing constant variance and co-efficient of variation; H0 equal variances
#'
#' @param sd standard deviation per dose group
#' @param n number of observations per dose group
#'
#' @return pvalue for the Bartlett's test
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
#' @return list of all available models
#'
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

#' Performs Shapiro-Wilk's test for normality on individual dose-response data
#'
#' @param datind a dataframe containing dose (must be named x) and response named y
#'
#'
#' @return a list of p-values and messages confirming the normality or lognormality of the data
#'
#' @export NLN_test
#'
NLN_test <- function(datind) {

  doseunique <- sort(unique(datind$x))
  normtestres <- c()
  # for each dose separately
  for (dd in doseunique){
    seldd <- (datind$x==dd)
    if(sum(seldd)>2){
      normtestres <- c(normtestres,shapiro.test(datind$y[seldd])$p.value)
    }else{
      normtestres <- c(normtestres, NA)
    }
  }
  # global test of anova residuals, more data but pooling and possibly missing local deviations
  normtestres <- c(normtestres, shapiro.test(lm(datind$y~as.factor(datind$x))$residuals)$p.value)
  #normtestres<-data.frame(normtestres)
  test05 <- (normtestres<0.05)
  test10 <- (normtestres<0.10)
  normtestres <- cbind(normtestres,test05,test10)
  colnames(normtestres) <- c("p-value Shapiro normality test","reject at 5%","reject at 10%")
  rnames<-c()
  for (dd in doseunique) rnames <- c(rnames,as.character(dd))
  rnames <- c(rnames,"global")
  rownames(normtestres) <- rnames


  # log-normality tests

  doseunique <- sort(unique(datind$x))
  lognormtestres <- c()
  # for each dose separately
  for (dd in doseunique){
    seldd <- (datind$x==dd)
    if(sum(seldd)>2){
      lognormtestres <- c(lognormtestres, shapiro.test(log(datind$y[seldd]))$p.value)
    }else{
      lognormtestres <- c(lognormtestres, NA)
    }
  }
  # global test of anova residuals, more data but pooling and possibly missing local deviations
  lognormtestres <- c(lognormtestres,shapiro.test(lm(log(datind$y)~as.factor(datind$x))$residuals)$p.value)
  #lognormtestres<-data.frame(lognormtestres)
  test05 <- (lognormtestres<0.05)
  test10 <- (lognormtestres<0.10)
  lognormtestres <- cbind(lognormtestres,test05,test10)
  colnames(lognormtestres) <- c("p-value Shapiro log-normality test","reject at 5%","reject at 10%")
  rnames <- c()
  for (dd in doseunique) rnames<-c(rnames,as.character(dd))
  rnames <- c(rnames,"global")
  rownames(lognormtestres) <- rnames

  #normtestres
  if (sum(normtestres[,2], na.rm = T)>0){
    warning("there is evidence against normality at level 5%")
    msg_5N <- "there is evidence against normality at level 5%"
  } else {
    warning("there is no evidence against normality at level 5%")
    msg_5N <- "there is no evidence against normality at level 5%"
  }
  if (sum(normtestres[,3], na.rm = T)>0) {
    warning("there is evidence against normality at level 10%")
    msg_10N <- "there is evidence against normality at level 10%"
  } else {
    warning("there is no evidence against normality at level 10%")
    msg_10N <- "there is no evidence against normality at level 10%"
  }

  #lognormtestres
  if (sum(lognormtestres[,2], na.rm = T)>0) {
    warning("there is evidence against log-normality at level 5%")
    msg_5LN <- "there is evidence against log-normality at level 5%"
  } else {
    warning("there is no evidence against log-normality at level 5%")
    msg_5LN <- "there is no evidence against log-normality at level 5%"
  }
  if (sum(lognormtestres[,3], na.rm = T)>0) {
    warning("there is evidence against log-normality at level 10%")
    msg_10LN <- "there is evidence against log-normality at level 10%"
  } else {
    warning("there is no evidence against log-normality at level 10%")
    msg_10LN <- "there is no evidence against log-normality at level 10%"
  }

  return(list(
    msg_5N = msg_5N,
    msg_10N = msg_10N,
    msg_5LN = msg_5LN,
    msg_10LN = msg_10LN,
    normal_test = normtestres,
    lognormal_test = lognormtestres

  ))
}





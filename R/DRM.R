#' The dose response models
#'
#' These are the dose response models used internally in the BMA functions.
#'
#' @param par parameters in the order a, c, d, b
#' @param x unique ordered dose levels
#' @param q specified BMR
#'
#' @return Vector containing the expected response at each dose level
#'
#' @example
#' data("immunotoxicityData.rda")  #load the immunotoxicity data
#' DRM.E4_NI(par = c(1.06, 0.015, 0.2, 1), x = immunotoxicityData$Dose[1:5], q = 0.1)
#'
#' @description the mean/median dose-response model per dose-level.
#'
#' @export
#'
DRM.E4_NI=function(par,x,q){
  a=exp(par[1])
  c=exp(par[3])+1 # c > 1, so take log to be >0 and then add 1 ?
  # a=par[1]
  # c=par[3]
  d=exp(par[4])
  b=-exp(-par[2]*d)*log(1-q/(c-1)) # par[2] is k (BMD)
  a*(1+(c-1)*(1-exp(-b*x^d)))
}
#' @rdname DRM.E4_NI
#' @export
DRM.IE4_NI=function(par,x,q){
  a=exp(par[1])
  c=exp(par[3])+1
  d=exp(par[4])
  b=-exp(par[2]*d)*log(q/(c-1))
  a*(1+(c-1)*exp(-b*x^(-d)))
}
#' @rdname DRM.E4_NI
#' @export
DRM.H4_NI=function(par,x,q){
  a=exp(par[1])
  c=exp(par[3])+1
  d=exp(par[4])
  b=exp(par[2]*d)*(((c-1)/q)-1)
  a*(1+(c-1)*(1-(b/(b+x^d))))
}
#' @rdname DRM.E4_NI
#' @export
DRM.LN4_NI=function(par,x,q){
  a=exp(par[1])
  c=exp(par[3])+1+q
  # c=exp(par[3])+1
  d=exp(par[4])
  b=exp(qnorm(q/(c-1))-(par[2]*d))
  a*(1+(c-1)*pnorm(log(b)+(d*log(x))))
}
#' @rdname DRM.E4_NI
#' @export
DRM.G4_NI=function(par,x,q){
  a=exp(par[1])
  c=exp(par[3])+1
  d=exp(par[4])
  b=qgamma(q/(c-1), rate=1.0, shape=d)/exp(par[2])
  a*(1+(c-1)*pgamma(x,shape=d,rate=b))
}
#' @rdname DRM.E4_NI
#' @export
DRM.QE4_NI=function(par,x,q){
  a=exp(par[1])
  c=exp(par[3])+1
  d=exp(par[4])
  b=(-log(1-(q/(c-1)))) / (exp(par[2])+((exp(par[2])*(exp(par[2])-1))/exp(d)))
  a*(1+(c-1)*(1-exp((-b*x)-(b/exp(d)*x*(x-1)))))
}
#' @rdname DRM.E4_NI
#' @export
DRM.P4_NI=function(par,x,q){
  a=exp(par[1])
  # c=par[3]
  c=qnorm(1/(1+q)) - exp(par[3])
  d=exp(par[4])
  b=exp(-par[2]*d)*(qnorm(pnorm(c)*(1+q))-c)
  a*pnorm(c+(b*x^d))
}
#' @rdname DRM.E4_NI
#' @export
DRM.L4_NI=function(par,x,q){
  a=exp(par[1])
  # c=par[3]
  c=logit(1/(1+q)) - exp(par[3])
  d=exp(par[4])
  b=exp(-par[2]*d)*(logit(expit(c)*(1+q))-c)
  a*expit(c+(b*x^d))
}
#' @rdname DRM.E4_NI
#' @export
DRM.E4_LNI=function(par,x,q){
  a=exp(par[1])
  c=exp(par[3])+1
  d=exp(par[4])
  b=-exp(-par[2]*d)*log(1-(log(1+q)/(a*(c-1))))
  a*(1+(c-1)*(1-exp(-b*x^d)))
}
#' @rdname DRM.E4_NI
#' @export
DRM.IE4_LNI=function(par,x,q){
  a=exp(par[1]) # no log transformation in this case
  c=exp(par[3])+1
  d=exp(par[4])
  b=-exp(par[2]*d)*log((log(1+q))/(a*(c-1)))
  a*(1+(c-1)*exp(-b*x^(-d)))
}
#' @rdname DRM.E4_NI
#' @export
DRM.H4_LNI=function(par,x,q){
  a=exp(par[1])
  c=exp(par[3])+1
  d=exp(par[4])
  b=exp(par[2]*d)*(((a*(c-1))/log(1+q))-1)
  a*(1+(c-1)*(1-(b/(b+x^d))))
}
#' @rdname DRM.E4_NI
#' @export
DRM.LN4_LNI=function(par,x,q){
  a=exp(par[1])
  c=exp(par[3])+(1+(log(1+q)/a))
  # c=exp(par[3])+1
  d=exp(par[4])
  b=exp(qnorm((log(1+q))/(a*(c-1)))-(par[2]*d))
  a*(1+(c-1)*pnorm(log(b)+(d*log(x))))
}

#' @rdname DRM.E4_NI
#' @export
DRM.G4_LNI=function(par,x,q){
  a=exp(par[1])
  c=exp(par[3])+1
  d=exp(par[4])
  b=qgamma(log(1+q)/(a*(c-1)), rate=1.0, shape=d)/exp(par[2])
  a*(1+(c-1)*pgamma(x,shape=d,rate=b))
}
#' @rdname DRM.E4_NI
#' @export
DRM.QE4_LNI=function(par,x,q){
  a=exp(par[1])
  c=exp(par[3])+1
  d=exp(par[4])
  b=(-log(1-(log(1+q)/(a*(c-1))))) / (exp(par[2])+((exp(par[2])*(exp(par[2])-1))/exp(d)))
  a*(1+(c-1)*(1-exp((-b*x)-(b/exp(d)*x*(x-1)))))
}
#' @rdname DRM.E4_NI
#' @export
DRM.P4_LNI=function(par,x,q){
  a=exp(par[1]) # no log transformation in this case
  c=qnorm(1 - (log(1+q)/a)) - exp(par[3])
  # c=par[3]
  d=exp(par[4])
  b=exp(-par[2]*d)*(qnorm(pnorm(c)+log(1+q)/a)-c)
  a*pnorm(c+(b*(x^d)))
}
#' @rdname DRM.E4_NI
#' @export
DRM.L4_LNI=function(par,x,q){
  a=exp(par[1]) # no log transformation in this case
  c=logit(1 - (log(1+q)/a)) - exp(par[3])
  # c=par[3]
  d=exp(par[4])
  b=exp(-par[2]*d)*(logit(expit(c)+log(1+q)/a)-c)
  a*expit(c+(b*(x^d)))
}
#' @rdname DRM.E4_NI
#' @export
DRM.E4_ND=function(par,x,q){
  a=exp(par[1])
  c=1/(1+exp(-par[3]))
  d=exp(par[4])
  b=-exp(-par[2]*d)*log(1+q/(c-1)) # par[2] is k (BMD)
  a*(1+(c-1)*(1-exp(-b*x^d)))
}
#' @rdname DRM.E4_NI
#' @export
DRM.IE4_ND=function(par,x,q){
  a=exp(par[1])
  c=1/(1+exp(-par[3]))
  d=exp(par[4])
  b=-exp(par[2]*d)*log(-q/(c-1))
  a*(1+(c-1)*exp(-b*x^(-d)))
}
#' @rdname DRM.E4_NI
#' @export
DRM.H4_ND=function(par,x,q){
  a=exp(par[1])
  c=1/(1+exp(-par[3]))
  d=exp(par[4])
  b=exp(par[2]*d)*(((c-1)/-q)-1)
  a*(1+(c-1)*(1-(b/(b+x^d))))
}
#' @rdname DRM.E4_NI
#' @export
DRM.LN4_ND=function(par,x,q){
  a=exp(par[1])
  # c=(1-q)*(1/(1+exp(-par[3])))
  c=1/(1+exp(-par[3]))
  d=exp(par[4])
  b=exp(qnorm(-q/(c-1))-(par[2]*d))
  a*(1+(c-1)*pnorm(log(b)+(d*log(x))))
}
#' @rdname DRM.E4_NI
#' @export
DRM.G4_ND=function(par,x,q){
  a=exp(par[1])
  # c=(1-q)*(1/(1+exp(-par[3])))
  c=1/(1+exp(-par[3]))
  d=exp(par[4])
  b=qgamma((-q)/(c-1), rate=1.0, shape=d)/exp(par[2])
  # fct = function(x) pgamma(exp(par[2]),shape=d,rate=x)*gamma(d) - (((-q)*gamma(d))/(c-1))
  # bh=try(uniroot(f=fct, interval=c(0,1000))$root,silent=T)
  # b=ifelse((is.numeric(bh) & bh>0),bh,0)
  a*(1+(c-1)*pgamma(x,shape=d,rate=b))
}
#' @rdname DRM.E4_NI
#' @export
DRM.QE4_ND=function(par,x,q){
  a=exp(par[1])
  c=1/(1+exp(-par[3]))
  d=exp(par[4])
  b=(-log(1-(-q/(c-1)))) / (exp(par[2])+((exp(par[2])*(exp(par[2])-1))/exp(d)))
  # a*(1+(c-1)*(1-exp((-b*x)-(b/exp(d)*x*(x-1)))))
  a*(1+(c-1)*(1-exp(
    (-b*x) - ((b*x^2)/exp(d)) + ((b*x)/exp(d))
  )
  ))
}
#' @rdname DRM.E4_NI
#' @export
DRM.P4_ND=function(par,x,q){
  a=exp(par[1])
  # c=par[3]
  c=qnorm(1-q) - exp(par[3])
  d=exp(par[4])
  b=exp(-par[2]*d)*(qnorm(pnorm(c)+q)-c)
  a*(1+pnorm(c))-(a*pnorm(c+(b*x^d)))
}
#' @rdname DRM.E4_NI
#' @export
DRM.L4_ND=function(par,x,q){
  a=exp(par[1])
  # c=par[3]
  c=logit(1-q) - exp(par[3])
  d=exp(par[4])
  b=exp(-par[2]*d)*(logit(expit(c)+q)-c)
  (a*(1+expit(c)))-(a*expit(c+(b*x^d)))
}
#' @rdname DRM.E4_NI
#' @export
DRM.E4_LND=function(par,x,q){
  a=exp(par[1])
  c=1/(1+exp(-par[3]))
  d=exp(par[4])
  b=-exp(-par[2]*d)*log(1-(log(1-q)/(a*(c-1))))
  a*(1+(c-1)*(1-exp(-b*x^d))) # DR model on log scale
}
#' @rdname DRM.E4_NI
#' @export
DRM.IE4_LND=function(par,x,q){
  a=exp(par[1]) # no log transformation in this case
  c=1/(1+exp(-par[3]))
  d=exp(par[4])
  b=-exp(par[2]*d)*log((log(1-q))/(a*(c-1)))
  a*(1+(c-1)*exp(-b*x^(-d)))
}
#' @rdname DRM.E4_NI
#' @export
DRM.H4_LND=function(par,x,q){
  a=exp(par[1]) # no log transformation in this case
  c=1/(1+exp(-par[3]))
  d=exp(par[4])
  b=exp(par[2]*d)*(((a*(c-1))/log(1-q))-1)
  a*(1+(c-1)*(1-(b/(b+x^d))))
}
#' @rdname DRM.E4_NI
#' @export
DRM.LN4_LND=function(par,x,q){
  a=exp(par[1])
  c=(1+(log(1-q)/a))*(1/(1+exp(-par[3])))
  # c=1/(1+exp(-par[3]))
  d=exp(par[4])
  b=exp(qnorm((log(1-q))/(a*(c-1)))-(par[2]*d))
  a*(1+(c-1)*pnorm(log(b)+(d*log(x))))
}
#' @rdname DRM.E4_NI
#' @export
DRM.G4_LND=function(par,x,q){
  a=exp(par[1])
  c=1/(1+exp(-par[3]))
  # c=(1+(log(1-q)/a))*(1/(1+exp(-par[3])))
  d=exp(par[4])
  b=qgamma(log(1-q)/(a*(c-1)), rate=1.0, shape=d)/exp(par[2])
  # fct = function(x) pgamma(exp(par[2]),shape=d,rate=x)*gamma(d) - ((log(1-q)*gamma(d))/(a*(c-1)))
  # bh=try(uniroot(f=fct, interval=c(0,1000))$root,silent=T)
  # b=ifelse((is.numeric(bh) & bh>0),bh,0)
  a*(1+(c-1)*pgamma(x,shape=d,rate=b))
}
#' @rdname DRM.E4_NI
#' @export
DRM.QE4_LND=function(par,x,q){
  a=exp(par[1])
  c=1/(1+exp(-par[3]))
  d=exp(par[4])
  b=(-log(1-(log(1-q)/(a*(c-1))))) / (exp(par[2])+((exp(par[2])*(exp(par[2])-1))/exp(d)))
  a*(1+(c-1)*(1-exp((-b*x)-(b/exp(d)*x*(x-1)))))
}
#' @rdname DRM.E4_NI
#' @export
DRM.P4_LND=function(par,x,q){
  a=exp(par[1]) # no log transformation in this case
  # c=par[3]
  c=qnorm(1+(log(1-q)/a)) - exp(par[3])
  d=exp(par[4])
  b=exp(-par[2]*d)*(qnorm(pnorm(c)-(log(1-q)/a))-c)
  (a*(1+pnorm(c)))-(a*pnorm(c+(b*(x^d))))
}
#' @rdname DRM.E4_NI
#' @export
DRM.L4_LND=function(par,x,q){
  a=exp(par[1]) # no log transformation in this case
  # c=par[3]
  c=logit(1+(log(1-q)/a)) - exp(par[3])
  d=exp(par[4])
  b=exp(-par[2]*d)*(logit(expit(c)-log(1-q)/a)-c)
  (a*(1+expit(c)))-(a*expit(c+(b*(x^d))))
}





#' These are the dose response models used internally in the BMA functions.
#'
#' @param par parameters in the order a, c, d, b
#' @param x unique ordered dose levels
#' @param q specified BMR
#'
#' @return Vector containing the expected response at each dose level
#'
#' @export
#'
DRM.E4_Q=function(par,x,q){
  a = par[1]
  d = exp(par[3])
  k = log(par[2])
  b = -exp(-k*d)*log(1-q) # k is k (BMD)
  a + (1 - a)*(1 - exp(-b*x^d)) - .Machine$double.xmin
}

#' @rdname DRM.E4_Q
#' @export
#'
DRM.IE4_Q=function(par,x,q){
  a = par[1]
  d = exp(par[3])
  k = log(par[2])
  b = -exp(k*d)*log(q)
  a + (1 - a)*exp(-b*x^-d) - .Machine$double.xmin
}
#' @rdname DRM.E4_Q
#' @export
DRM.H4_Q=function(par,x,q){
  a = par[1]
  d = exp(par[3])
  k = log(par[2])
  b = exp(k*d)*((1/q)-1)
  a + (1 - a)*(1 - (b / (b + x^d))) - .Machine$double.xmin
}
#' @rdname DRM.E4_Q
#' @export
DRM.LN4_Q=function(par,x,q){
  a = par[1]
  d = exp(par[3])
  k = log(par[2])
  b = exp(qnorm(q)-(k*d))
  a + (1 - a)*pnorm(log(b) + d*log(x)) - .Machine$double.xmin
}
#' @rdname DRM.E4_Q
#' @export
DRM.G4_Q=function(par,x,q){
  a = par[1]
  d = exp(par[3])+0.000000001
  k = log(par[2])
  b = qgamma(q, rate=1.0, shape=d)/exp(k)
  a + (1 - a)*pgamma(x, shape = d, rate = b) - .Machine$double.xmin
}
#' @rdname DRM.E4_Q
#' @export
DRM.QE4_Q=function(par,x,q){
  a = par[1]
  d = exp(par[3])
  k = log(par[2])
  b = (-log(1-q)) / (exp(k)+((exp(k)*(exp(k)-1))/exp(d)))
  a + (1 - a)*(1 - exp(-b*x - ( (b/exp(d)) * x * (x - 1)))) - .Machine$double.xmin
}
#' @rdname DRM.E4_Q
#' @export
DRM.P4_Q=function(par,x,q){
  a = par[1]
  d = exp(par[3])
  k = log(par[2])
  b = (qnorm(q*(1-a) + a) - qnorm(a))/exp(k*d);
  pnorm(qnorm(a) + b*x^d) - .Machine$double.xmin

}
#' @rdname DRM.E4_Q
#' @export
DRM.L4_Q=function(par,x,q){
  a = par[1]
  d = exp(par[3])
  k = log(par[2])
  b = exp(-k*d)*(logit(q*(1-a)+a) - logit(a))
  expit(logit(a) + b*x^d) - .Machine$double.xmin

}

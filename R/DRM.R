#' The dose response models
#'
#' These are the dose response models used internally in the BMA functions.
#'
#' @param par parameters in the order min response, BMD, fold change (for continuous data only), d
#' @param x unique rescaled ordered dose levels
#' @param q specified BMR
#' @param shift shift in the response in case of negative geometric means
#'
#' @return Vector containing the expected response at each (rescaled) dose level
#'
#' @examples
#' DRM.E4_NI(par = c(1.06, 0.3, 5.6, 1), x = immunotoxicityData$Dose[1:5]/max(immunotoxicityData$Dose[1:5]), q = 0.1)
#'
#' @description the mean/median response per dose-level
#'
#' @export
#'
DRM.E4_NI=function(par,x,q,shift=0){
  a = par[1]
  mu_inf = par[1]*par[3];
  c = mu_inf/par[1];
  d = exp(par[4])
  k = log(par[2])
  b=-exp(-k*d)*log(1-q/(c-1))
  a*(1+(c-1)*(1-exp(-b*x^d)))
}
#' @rdname DRM.E4_NI
#' @export
DRM.IE4_NI=function(par,x,q,shift=0){
  a = par[1]
  mu_inf = par[1]*par[3];
  c = mu_inf/par[1];
  d = exp(par[4])
  k = log(par[2])
  b=-exp(k*d)*log(q/(c-1))
  a*(1+(c-1)*exp(-b*x^(-d)))
}
#' @rdname DRM.E4_NI
#' @export
DRM.H4_NI=function(par,x,q,shift=0){
  a = par[1]
  mu_inf = par[1]*par[3];
  c = mu_inf/par[1];
  d = exp(par[4])
  k = log(par[2])
  b=exp(k*d)*(((c-1)/q)-1)
  a*(1+(c-1)*(1-(b/(b+x^d))))
}
#' @rdname DRM.E4_NI
#' @export
DRM.LN4_NI=function(par,x,q,shift=0){
  a = par[1]
  mu_inf = par[1]*par[3];
  c = mu_inf/par[1];
  d = exp(par[4])
  k = log(par[2])
  b=exp(qnorm(q/(c-1))-(k*d))
  a*(1+(c-1)*pnorm(log(b)+(d*log(x))))
}
#' @rdname DRM.E4_NI
#' @export
DRM.G4_NI=function(par,x,q,shift=0){
  a = par[1]
  mu_inf = par[1]*par[3];
  c = mu_inf/par[1];
  d = exp(par[4])
  k = log(par[2])
  b=qgamma(q/(c-1), rate=1.0, shape=d)/exp(k)
  a*(1+(c-1)*pgamma(x,shape=d,rate=b))
}
#' @rdname DRM.E4_NI
#' @export
DRM.QE4_NI=function(par,x,q,shift=0){
  a = par[1]
  mu_inf = par[1]*par[3];
  c = mu_inf/par[1];
  d = exp(par[4])
  k = log(par[2])
  b=(-log(1-(q/(c-1)))) / (exp(k)+(exp(k)*((exp(k)-1))/exp(d)))
  # fct = function(x) log(1-(q/(c-1))) + x*exp(k) + d*exp(2*k)
  #b=max(try(uniroot(f=fct, c(-10000,1000))$root),0)
  # bh=try(uniroot(f=fct, interval=c(-1000,1000))$root,silent=T)
  # b=ifelse((is.numeric(bh) & bh>0),bh,0)
  a*(1+(c-1)*(1-exp((-b*x)-(b/exp(d)*x*(x-1)))))
  # a*(1+(c-1)*(1-exp(-b*x-d*(x^2))))
}
#' @rdname DRM.E4_NI
#' @export
DRM.P4_NI=function(par,x,q,shift=0){
  mu_0 = par[1]
  mu_inf = par[1]*par[3]
  a = mu_inf
  c = qnorm(mu_0/mu_inf);
  d = exp(par[4])
  k = log(par[2])
  b=exp(-k*d)*(qnorm(pnorm(c)*(1+q))-c)
  a*pnorm(c+(b*x^d))
}
#' @rdname DRM.E4_NI
#' @export
DRM.L4_NI=function(par,x,q,shift=0){
  mu_0 = par[1]
  mu_inf = par[1]*par[3]
  a = mu_inf
  c = logit(mu_0/mu_inf);
  d = exp(par[4])
  k = log(par[2])
  b=exp(-k*d)*(logit(expit(c)*(1+q))-c)
  a*expit(c+(b*x^d))
}
#' @rdname DRM.E4_NI
#' @export
DRM.E4_LNI=function(par,x,q,shift){
  a = log(par[1])-shift
  mu_inf = par[1]*par[3];
  c = (log(mu_inf)-shift)/a;
  d = exp(par[4])
  k = log(par[2])
  b=-exp(-k*d)*log(1-(log(1+q)/(a*(c-1))))
  a*(1+(c-1)*(1-exp(-b*x^d))) + shift
}
#' @rdname DRM.E4_NI
#' @export
DRM.IE4_LNI=function(par,x,q,shift){
  a = log(par[1])-shift
  mu_inf = par[1]*par[3];
  c = (log(mu_inf)-shift)/a;
  d = exp(par[4])
  k = log(par[2])
  b=-exp(k*d)*log((log(1+q))/(a*(c-1)))
  a*(1+(c-1)*exp(-b*x^(-d))) + shift
}
#' @rdname DRM.E4_NI
#' @export
DRM.H4_LNI=function(par,x,q,shift){
  a = log(par[1])-shift
  mu_inf = par[1]*par[3];
  c = (log(mu_inf)-shift)/a;
  d = exp(par[4])
  k = log(par[2])
  b=exp(k*d)*(((a*(c-1))/log(1+q))-1)
  a*(1+(c-1)*(1-(b/(b+x^d)))) + shift
}
#' @rdname DRM.E4_NI
#' @export
DRM.LN4_LNI=function(par,x,q,shift){
  a = log(par[1])-shift
  mu_inf = par[1]*par[3];
  c = (log(mu_inf)-shift)/a;
  d = exp(par[4])
  k = log(par[2])
  b=exp(qnorm((log(1+q))/(a*(c-1)))-(k*d))
  a*(1+(c-1)*pnorm(log(b)+(d*log(x)))) + shift
}

#' @rdname DRM.E4_NI
#' @export
DRM.G4_LNI=function(par,x,q,shift){
  a = log(par[1])-shift
  mu_inf = par[1]*par[3];
  c = (log(mu_inf)-shift)/a;
  d = exp(par[4])
  k = log(par[2])
  b=qgamma(log(1+q)/(a*(c-1)), rate=1.0, shape=d)/exp(k)
  a*(1+(c-1)*pgamma(x,shape=d,rate=b)) + shift
}
#' @rdname DRM.E4_NI
#' @export
DRM.QE4_LNI=function(par,x,q,shift){
  a = log(par[1])-shift
  mu_inf = par[1]*par[3];
  c = (log(mu_inf)-shift)/a;
  d = exp(par[4])
  k = log(par[2])
  b=(-log(1-(log(1+q)/(a*(c-1))))) / (exp(k)+((exp(k)*(exp(k)-1))/exp(d)))
  a*(1+(c-1)*(1-exp((-b*x)-(b/exp(d)*x*(x-1))))) + shift
}
#' @rdname DRM.E4_NI
#' @export
DRM.P4_LNI=function(par,x,q,shift){
  mu_0 = par[1]
  mu_inf = par[1]*par[3]
  a = log(mu_inf)-shift
  c = qnorm((log(mu_0)-shift)/(log(mu_inf)-shift));
  d = exp(par[4])
  k = log(par[2])
  b=exp(-k*d)*(qnorm(pnorm(c)+log(1+q)/a)-c)
  a*pnorm(c+(b*(x^d))) + shift
}
#' @rdname DRM.E4_NI
#' @export
DRM.L4_LNI=function(par,x,q,shift){
  mu_0 = par[1]
  mu_inf = par[1]*par[3]
  a = log(mu_inf)-shift
  c = logit((log(mu_0)-shift)/(log(mu_inf)-shift));
  d = exp(par[4])
  k = log(par[2])
  b=exp(-k*d)*(logit(expit(c)+log(1+q)/a)-c)
  a*expit(c+(b*(x^d))) + shift
}
#' @rdname DRM.E4_NI
#' @export
DRM.E4_ND=function(par,x,q,shift=0){
  a = par[1]
  mu_inf = par[1]*par[3];
  c = mu_inf/par[1];
  d = exp(par[4])
  k = log(par[2])
  b=-exp(-k*d)*log(1+q/(c-1)) # par[2] is k (BMD)
  a*(1+(c-1)*(1-exp(-b*x^d)))
}
#' @rdname DRM.E4_NI
#' @export
DRM.IE4_ND=function(par,x,q,shift=0){
  a = par[1]
  mu_inf = par[1]*par[3];
  c = mu_inf/par[1];
  d = exp(par[4])
  k = log(par[2])
  b=-exp(k*d)*log(-q/(c-1))
  a*(1+(c-1)*exp(-b*x^(-d)))
}
#' @rdname DRM.E4_NI
#' @export
DRM.H4_ND=function(par,x,q,shift=0){
  a = par[1]
  mu_inf = par[1]*par[3];
  c = mu_inf/par[1];
  d = exp(par[4])
  k = log(par[2])
  b=exp(k*d)*(((c-1)/-q)-1)
  a*(1+(c-1)*(1-(b/(b+x^d))))
}
#' @rdname DRM.E4_NI
#' @export
DRM.LN4_ND=function(par,x,q,shift=0){
  a = par[1]
  mu_inf = par[1]*par[3];
  c = mu_inf/par[1];
  d = exp(par[4])
  k = log(par[2])
  b=exp(qnorm(-q/(c-1))-(k*d))
  a*(1+(c-1)*pnorm(log(b)+(d*log(x))))
}
#' @rdname DRM.E4_NI
#' @export
DRM.G4_ND=function(par,x,q,shift=0){
  a = par[1]
  mu_inf = par[1]*par[3];
  c = mu_inf/par[1];
  d = exp(par[4])
  k = log(par[2])
  b=qgamma((-q)/(c-1), rate=1.0, shape=d)/exp(k)
  a*(1+(c-1)*pgamma(x,shape=d,rate=b))
}
#' @rdname DRM.E4_NI
#' @export
DRM.QE4_ND=function(par,x,q,shift=0){
  a = par[1]
  mu_inf = par[1]*par[3];
  c = mu_inf/par[1];
  d = exp(par[4])
  k = log(par[2])
  b=(-log(1-(-q/(c-1)))) / (exp(k)+((exp(k)*(exp(k)-1))/exp(d)))
  # a*(1+(c-1)*(1-exp((-b*x)-(b/exp(d)*x*(x-1)))))
  a*(1+(c-1)*(1-exp(
    (-b*x) - ((b*x^2)/exp(d)) + ((b*x)/exp(d))
  )
  ))
}
#' @rdname DRM.E4_NI
#' @export
DRM.P4_ND=function(par,x,q,shift=0){
  mu_0 = par[1]
  mu_inf = par[1]*par[3]
  a = mu_0
  c = qnorm(mu_inf/mu_0);
  d = exp(par[4])
  k = log(par[2])
  b=exp(-k*d)*(qnorm(pnorm(c)+q)-c)
  a*(1+pnorm(c))-(a*pnorm(c+(b*x^d)))
}
#' @rdname DRM.E4_NI
#' @export
DRM.L4_ND=function(par,x,q,shift=0){
  mu_0 = par[1]
  mu_inf = par[1]*par[3]
  a = mu_0
  c = logit(mu_inf/mu_0);
  d = exp(par[4])
  k = log(par[2])
  b=exp(-k*d)*(logit(expit(c)+q)-c)
  (a*(1+expit(c)))-(a*expit(c+(b*x^d)))
}
#' @rdname DRM.E4_NI
#' @export
DRM.E4_LND=function(par,x,q,shift){
  a = log(par[1])-shift
  mu_inf = par[1]*par[3];
  c = (log(mu_inf)-shift)/a;
  d = exp(par[4])
  k = log(par[2])
  b=-exp(-k*d)*log(1-(log(1-q)/(a*(c-1))))
  a*(1+(c-1)*(1-exp(-b*x^d))) + shift
}
#' @rdname DRM.E4_NI
#' @export
DRM.IE4_LND=function(par,x,q,shift){
  a = log(par[1])-shift
  mu_inf = par[1]*par[3];
  c = (log(mu_inf)-shift)/a;
  d = exp(par[4])
  k = log(par[2])
  b=-exp(k*d)*log((log(1-q))/(a*(c-1)))
  a*(1+(c-1)*exp(-b*x^(-d))) + shift
}
#' @rdname DRM.E4_NI
#' @export
DRM.H4_LND=function(par,x,q,shift){
  a = log(par[1])-shift
  mu_inf = par[1]*par[3];
  c = (log(mu_inf)-shift)/a;
  d = exp(par[4])
  k = log(par[2])
  b=exp(k*d)*(((a*(c-1))/log(1-q))-1)
  a*(1+(c-1)*(1-(b/(b+x^d)))) + shift
}
#' @rdname DRM.E4_NI
#' @export
DRM.LN4_LND=function(par,x,q,shift){
  a = log(par[1])-shift
  mu_inf = par[1]*par[3];
  c = (log(mu_inf)-shift)/a;
  d = exp(par[4])
  k = log(par[2])
  b=exp(qnorm((log(1-q))/(a*(c-1)))-(k*d))
  a*(1+(c-1)*pnorm(log(b)+(d*log(x)))) + shift
}
#' @rdname DRM.E4_NI
#' @export
DRM.G4_LND=function(par,x,q,shift){
  a = log(par[1]) - shift
  mu_inf = par[1]*par[3];
  c = (log(mu_inf)-shift)/a;
  d = exp(par[4])
  k = log(par[2])
  b=qgamma(log(1-q)/(a*(c-1)), rate=1.0, shape=d)/exp(k)
  a*(1+(c-1)*pgamma(x,shape=d,rate=b)) + shift
}
#' @rdname DRM.E4_NI
#' @export
DRM.QE4_LND=function(par,x,q,shift){
  a = log(par[1])-shift
  mu_inf = par[1]*par[3];
  c = (log(mu_inf)-shift)/a;
  d = exp(par[4])
  k = log(par[2])
  b=(-log(1-(log(1-q)/(a*(c-1))))) / (exp(k)+((exp(k)*(exp(k)-1))/exp(d)))
  a*(1+(c-1)*(1-exp((-b*x)-(b/exp(d)*x*(x-1))))) + shift
}
#' @rdname DRM.E4_NI
#' @export
DRM.P4_LND=function(par,x,q,shift){
  mu_0 = par[1]
  mu_inf = par[1]*par[3]
  a = log(mu_0)-shift
  c = qnorm((log(mu_inf)-shift)/(log(mu_0)-shift));
  d = exp(par[4])
  k = log(par[2])
  b=exp(-k*d)*(qnorm(pnorm(c)-(log(1-q)/a))-c)
  (a*(1+pnorm(c)))-(a*pnorm(c+(b*(x^d)))) + shift
}
#' @rdname DRM.E4_NI
#' @export
DRM.L4_LND=function(par,x,q,shift){
  mu_0 = par[1]
  mu_inf = par[1]*par[3]
  a = log(mu_0) - shift
  c = logit((log(mu_inf)-shift)/(log(mu_0)-shift));
  d = exp(par[4])
  k = log(par[2])
  b=exp(-k*d)*(logit(expit(c)-log(1-q)/a)-c)
  (a*(1+expit(c)))-(a*expit(c+(b*(x^d)))) + shift
}


#' @rdname DRM.E4_NI
#' @export
DRM.E4_Q=function(par,x,q){
  a = par[1]
  d = exp(par[3])
  k = log(par[2])
  b = -exp(-k*d)*log(1-q) # k is k (BMD)
  a + (1 - a)*(1 - exp(-b*x^d)) - .Machine$double.xmin
}

#' @rdname DRM.E4_NI
#' @export
#'
DRM.IE4_Q=function(par,x,q){
  a = par[1]
  d = exp(par[3])
  k = log(par[2])
  b = -exp(k*d)*log(q)
  a + (1 - a)*exp(-b*x^-d) - .Machine$double.xmin
}
#' @rdname DRM.E4_NI
#' @export
DRM.H4_Q=function(par,x,q){
  a = par[1]
  d = exp(par[3])
  k = log(par[2])
  b = exp(k*d)*((1/q)-1)
  a + (1 - a)*(1 - (b / (b + x^d))) - .Machine$double.xmin
}
#' @rdname DRM.E4_NI
#' @export
DRM.LN4_Q=function(par,x,q){
  a = par[1]
  d = exp(par[3])
  k = log(par[2])
  b = exp(qnorm(q)-(k*d))
  a + (1 - a)*pnorm(log(b) + d*log(x)) - .Machine$double.xmin
}
#' @rdname DRM.E4_NI
#' @export
DRM.G4_Q=function(par,x,q){
  a = par[1]
  d = exp(par[3])+0.000000001
  k = log(par[2])
  b = qgamma(q, rate=1.0, shape=d)/exp(k)
  a + (1 - a)*pgamma(x, shape = d, rate = b) - .Machine$double.xmin
}
#' @rdname DRM.E4_NI
#' @export
DRM.QE4_Q=function(par,x,q){
  a = par[1]
  d = exp(par[3])
  k = log(par[2])
  b = (-log(1-q)) / (exp(k)+((exp(k)*(exp(k)-1))/exp(d)))
  a + (1 - a)*(1 - exp((-b*x) - ( (b/exp(d)) * x * (x - 1)))) - .Machine$double.xmin
}

#' @rdname DRM.E4_NI
#' @export
DRM.P4_Q=function(par,x,q){
  a = par[1]
  d = exp(par[3])
  k = log(par[2])
  b = (qnorm(q*(1-a) + a) - qnorm(a))/exp(k*d);
  pnorm(qnorm(a) + b*x^d) - .Machine$double.xmin

}
#' @rdname DRM.E4_NI
#' @export
DRM.L4_Q=function(par,x,q){
  a = par[1]
  d = exp(par[3])
  k = log(par[2])
  b = exp(-k*d)*(logit(q*(1-a)+a) - logit(a))
  expit(logit(a) + b*x^d) - .Machine$double.xmin
}


#' The dose response models
#'
#' These are the dose response models used internally in the BMA functions.
#'
#' @param par parameters in the order background, BMD, fold change, log(d), invsigma2
#' @param x unique ordered dose levels
#' @param q specified BMR
#'
#' @return Vector containing the expected response at each dose level
#'
#' @export
#'
DRM.E4_NI=function(par,x,q){
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
DRM.IE4_NI=function(par,x,q){
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
DRM.H4_NI=function(par,x,q){
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
DRM.LN4_NI=function(par,x,q){
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
DRM.G4_NI=function(par,x,q){
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
DRM.QE4_NI=function(par,x,q){
  a = par[1]
  mu_inf = par[1]*par[3];
  c = mu_inf/par[1];
  d = exp(par[4])
  k = log(par[2])
  b=(-log(1-(q/(c-1)))) / (exp(k)+((exp(k)*(exp(k)-1))/exp(d)))
  a*(1+(c-1)*(1-exp((-b*x)-(b/exp(d)*x*(x-1)))))
}
#' @rdname DRM.E4_NI
#' @export
DRM.P4_NI=function(par,x,q){
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
DRM.L4_NI=function(par,x,q){
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
DRM.E4_LNI=function(par,x,q){
  a = log(par[1])
  mu_inf = par[1]*par[3];
  c = log(mu_inf)/a;
  d = exp(par[4])
  k = log(par[2])
  b=-exp(-k*d)*log(1-(log(1+q)/(a*(c-1))))
  a*(1+(c-1)*(1-exp(-b*x^d)))
}
#' @rdname DRM.E4_NI
#' @export
DRM.IE4_LNI=function(par,x,q){
  a = log(par[1])
  mu_inf = par[1]*par[3];
  c = log(mu_inf)/a;
  d = exp(par[4])
  k = log(par[2])
  b=-exp(k*d)*log((log(1+q))/(a*(c-1)))
  a*(1+(c-1)*exp(-b*x^(-d)))
}
#' @rdname DRM.E4_NI
#' @export
DRM.H4_LNI=function(par,x,q){
  a = log(par[1])
  mu_inf = par[1]*par[3];
  c = log(mu_inf)/a;
  d = exp(par[4])
  k = log(par[2])
  b=exp(k*d)*(((a*(c-1))/log(1+q))-1)
  a*(1+(c-1)*(1-(b/(b+x^d))))
}
#' @rdname DRM.E4_NI
#' @export
DRM.LN4_LNI=function(par,x,q){
  a = log(par[1])
  mu_inf = par[1]*par[3];
  c = log(mu_inf)/a;
  d = exp(par[4])
  k = log(par[2])
  b=exp(qnorm((log(1+q))/(a*(c-1)))-(k*d))
  a*(1+(c-1)*pnorm(log(b)+(d*log(x))))
}

#' @rdname DRM.E4_NI
#' @export
DRM.G4_LNI=function(par,x,q){
  a = log(par[1])
  mu_inf = par[1]*par[3];
  c = log(mu_inf)/a;
  d = exp(par[4])
  k = log(par[2])
  b=qgamma(log(1+q)/(a*(c-1)), rate=1.0, shape=d)/exp(k)
  a*(1+(c-1)*pgamma(x,shape=d,rate=b))
}
#' @rdname DRM.E4_NI
#' @export
DRM.QE4_LNI=function(par,x,q){
  a = log(par[1])
  mu_inf = par[1]*par[3];
  c = log(mu_inf)/a;
  d = exp(par[4])
  k = log(par[2])
  b=(-log(1-(log(1+q)/(a*(c-1))))) / (exp(k)+((exp(k)*(exp(k)-1))/exp(d)))
  a*(1+(c-1)*(1-exp((-b*x)-(b/exp(d)*x*(x-1)))))
}
#' @rdname DRM.E4_NI
#' @export
DRM.P4_LNI=function(par,x,q){
  mu_0 = par[1]
  mu_inf = par[1]*par[3]
  a = log(mu_inf)
  c = qnorm(log(mu_0)/log(mu_inf));
  d = exp(par[4])
  k = log(par[2])
  b=exp(-k*d)*(qnorm(pnorm(c)+log(1+q)/a)-c)
  a*pnorm(c+(b*(x^d)))
}
#' @rdname DRM.E4_NI
#' @export
DRM.L4_LNI=function(par,x,q){
  mu_0 = par[1]
  mu_inf = par[1]*par[3]
  a = log(mu_inf)
  c = logit(log(mu_0)/log(mu_inf));
  d = exp(par[4])
  k = log(par[2])
  b=exp(-k*d)*(logit(expit(c)+log(1+q)/a)-c)
  a*expit(c+(b*(x^d)))
}
#' @rdname DRM.E4_NI
#' @export
DRM.E4_ND=function(par,x,q){
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
DRM.IE4_ND=function(par,x,q){
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
DRM.H4_ND=function(par,x,q){
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
DRM.LN4_ND=function(par,x,q){
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
DRM.G4_ND=function(par,x,q){
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
DRM.QE4_ND=function(par,x,q){
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
DRM.P4_ND=function(par,x,q){
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
DRM.L4_ND=function(par,x,q){
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
DRM.E4_LND=function(par,x,q){
  a = log(par[1])
  mu_inf = par[1]*par[3];
  c = log(mu_inf)/a;
  d = exp(par[4])
  k = log(par[2])
  b=-exp(-k*d)*log(1-(log(1-q)/(a*(c-1))))
  a*(1+(c-1)*(1-exp(-b*x^d))) # DR model on log scale
}
#' @rdname DRM.E4_NI
#' @export
DRM.IE4_LND=function(par,x,q){
  a = log(par[1])
  mu_inf = par[1]*par[3];
  c = log(mu_inf)/a;
  d = exp(par[4])
  k = log(par[2])
  b=-exp(k*d)*log((log(1-q))/(a*(c-1)))
  a*(1+(c-1)*exp(-b*x^(-d)))
}
#' @rdname DRM.E4_NI
#' @export
DRM.H4_LND=function(par,x,q){
  a = log(par[1])
  mu_inf = par[1]*par[3];
  c = log(mu_inf)/a;
  d = exp(par[4])
  k = log(par[2])
  b=exp(k*d)*(((a*(c-1))/log(1-q))-1)
  a*(1+(c-1)*(1-(b/(b+x^d))))
}
#' @rdname DRM.E4_NI
#' @export
DRM.LN4_LND=function(par,x,q){
  a = log(par[1])
  mu_inf = par[1]*par[3];
  c = log(mu_inf)/a;
  d = exp(par[4])
  k = log(par[2])
  b=exp(qnorm((log(1-q))/(a*(c-1)))-(k*d))
  a*(1+(c-1)*pnorm(log(b)+(d*log(x))))
}
#' @rdname DRM.E4_NI
#' @export
DRM.G4_LND=function(par,x,q){
  a = log(par[1])
  mu_inf = par[1]*par[3];
  c = log(mu_inf)/a;
  d = exp(par[4])
  k = log(par[2])
  b=qgamma(log(1-q)/(a*(c-1)), rate=1.0, shape=d)/exp(k)
  a*(1+(c-1)*pgamma(x,shape=d,rate=b))
}
#' @rdname DRM.E4_NI
#' @export
DRM.QE4_LND=function(par,x,q){
  a = log(par[1])
  mu_inf = par[1]*par[3];
  c = log(mu_inf)/a;
  d = exp(par[4])
  k = log(par[2])
  b=(-log(1-(log(1-q)/(a*(c-1))))) / (exp(k)+((exp(k)*(exp(k)-1))/exp(d)))
  a*(1+(c-1)*(1-exp((-b*x)-(b/exp(d)*x*(x-1)))))
}
#' @rdname DRM.E4_NI
#' @export
DRM.P4_LND=function(par,x,q){
  mu_0 = par[1]
  mu_inf = par[1]*par[3]
  a = log(mu_0)
  c = qnorm(log(mu_inf)/log(mu_0));
  d = exp(par[4])
  k = log(par[2])
  b=exp(-k*d)*(qnorm(pnorm(c)-(log(1-q)/a))-c)
  (a*(1+pnorm(c)))-(a*pnorm(c+(b*(x^d))))
}
#' @rdname DRM.E4_NI
#' @export
DRM.L4_LND=function(par,x,q){
  mu_0 = par[1]
  mu_inf = par[1]*par[3]
  a = log(mu_0)
  c = logit(log(mu_inf)/log(mu_0));
  d = exp(par[4])
  k = log(par[2])
  b=exp(-k*d)*(logit(expit(c)-log(1-q)/a)-c)
  (a*(1+expit(c)))-(a*expit(c+(b*(x^d))))
}

#' The dose response models
#'
#' These are the dose response models used internally in the BMA functions.
#'
#' @param par parameters in the order background, BMD, log(d)
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



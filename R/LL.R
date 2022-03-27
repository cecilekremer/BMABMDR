#' Loglikelihood functions for the different dose response models
#'
#' @param x parameter values
#' @param nvec number of observations
#' @param dvec dose levels
#' @param mvec response
#' @param s2vec variance
#' @param qval BMR
#'
#' @return .
#'
llfE4_NI=function(x,nvec,dvec,mvec,s2vec,qval){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[5]-0.5*(nvec-1)*s2vec*exp(x[5])-0.5*nvec*((mvec-DRM.E4_NI(x[1:4],dvec,qval))^2)*exp(x[5]))
}
#' @rdname llfE4_NI
llfIE4_NI=function(x,nvec,dvec,mvec,s2vec,qval){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[5]-0.5*(nvec-1)*s2vec*exp(x[5])-0.5*nvec*((mvec-DRM.IE4_NI(x[1:4],dvec,qval))^2)*exp(x[5]))
}
#' @rdname llfE4_NI
llfH4_NI=function(x,nvec,dvec,mvec,s2vec,qval){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[5]-0.5*(nvec-1)*s2vec*exp(x[5])-0.5*nvec*((mvec-DRM.H4_NI(x[1:4],dvec,qval))^2)*exp(x[5]))
}
#' @rdname llfE4_NI
llfLN4_NI=function(x,nvec,dvec,mvec,s2vec,qval){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[5]-0.5*(nvec-1)*s2vec*exp(x[5])-0.5*nvec*((mvec-DRM.LN4_NI(x[1:4],dvec,qval))^2)*exp(x[5]))
}
#' @rdname llfE4_NI
llfG4_NI=function(x,nvec,dvec,mvec,s2vec,qval){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[5]-0.5*(nvec-1)*s2vec*exp(x[5])-0.5*nvec*((mvec-DRM.G4_NI(x[1:4],dvec,qval))^2)*exp(x[5]))
}
#' @rdname llfE4_NI
llfQE4_NI=function(x,nvec,dvec,mvec,s2vec,qval){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[5]-0.5*(nvec-1)*s2vec*exp(x[5])-0.5*nvec*((mvec-DRM.QE4_NI(x[1:4],dvec,qval))^2)*exp(x[5]))
}
#' @rdname llfE4_NI
llfP4_NI=function(x,nvec,dvec,mvec,s2vec,qval){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[5]-0.5*(nvec-1)*s2vec*exp(x[5])-0.5*nvec*((mvec-DRM.P4_NI(x[1:4],dvec,qval))^2)*exp(x[5]))
}
#' @rdname llfE4_NI
llfL4_NI=function(x,nvec,dvec,mvec,s2vec,qval){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[5]-0.5*(nvec-1)*s2vec*exp(x[5])-0.5*nvec*((mvec-DRM.L4_NI(x[1:4],dvec,qval))^2)*exp(x[5]))
}
#' @rdname llfE4_NI
llfE4_LNI=function(x,nvec,dvec,mvec,s2vec,qval){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[5]-0.5*(nvec-1)*s2vec*exp(x[5])-0.5*nvec*((mvec-DRM.E4_LNI(x[1:4],dvec,qval))^2)*exp(x[5])) - sum(mvec*nvec)
}
#' @rdname llfE4_NI
llfIE4_LNI=function(x,nvec,dvec,mvec,s2vec,qval){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[5]-0.5*(nvec-1)*s2vec*exp(x[5])-0.5*nvec*((mvec-DRM.IE4_LNI(x[1:4],dvec,qval))^2)*exp(x[5])) - sum(mvec*nvec)
}
#' @rdname llfE4_NI
llfH4_LNI=function(x,nvec,dvec,mvec,s2vec,qval){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[5]-0.5*(nvec-1)*s2vec*exp(x[5])-0.5*nvec*((mvec-DRM.H4_LNI(x[1:4],dvec,qval))^2)*exp(x[5])) - sum(mvec*nvec)
}
#' @rdname llfE4_NI
llfLN4_LNI=function(x,nvec,dvec,mvec,s2vec,qval){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[5]-0.5*(nvec-1)*s2vec*exp(x[5])-0.5*nvec*((mvec-DRM.LN4_LNI(x[1:4],dvec,qval))^2)*exp(x[5])) - sum(mvec*nvec)
}
#' @rdname llfE4_NI
llfG4_LNI=function(x,nvec,dvec,mvec,s2vec,qval){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[5]-0.5*(nvec-1)*s2vec*exp(x[5])-0.5*nvec*((mvec-DRM.G4_LNI(x[1:4],dvec,qval))^2)*exp(x[5])) - sum(mvec*nvec)
}
#' @rdname llfE4_NI
llfQE4_LNI=function(x,nvec,dvec,mvec,s2vec,qval){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[5]-0.5*(nvec-1)*s2vec*exp(x[5])-0.5*nvec*((mvec-DRM.QE4_LNI(x[1:4],dvec,qval))^2)*exp(x[5])) - sum(mvec*nvec)
}
#' @rdname llfE4_NI
llfP4_LNI=function(x,nvec,dvec,mvec,s2vec,qval){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[5]-0.5*(nvec-1)*s2vec*exp(x[5])-0.5*nvec*((mvec-DRM.P4_LNI(x[1:4],dvec,qval))^2)*exp(x[5])) - sum(mvec*nvec)
}
#' @rdname llfE4_NI
llfL4_LNI=function(x,nvec,dvec,mvec,s2vec,qval){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[5]-0.5*(nvec-1)*s2vec*exp(x[5])-0.5*nvec*((mvec-DRM.L4_LNI(x[1:4],dvec,qval))^2)*exp(x[5])) - sum(mvec*nvec)
}
#' @rdname llfE4_NI
llfE4_ND=function(x,nvec,dvec,mvec,s2vec,qval){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[5]-0.5*(nvec-1)*s2vec*exp(x[5])-0.5*nvec*((mvec-DRM.E4_ND(x[1:4],dvec,qval))^2)*exp(x[5]))
}
#' @rdname llfE4_NI
llfIE4_ND=function(x,nvec,dvec,mvec,s2vec,qval){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[5]-0.5*(nvec-1)*s2vec*exp(x[5])-0.5*nvec*((mvec-DRM.IE4_ND(x[1:4],dvec,qval))^2)*exp(x[5]))
}
#' @rdname llfE4_NI
llfH4_ND=function(x,nvec,dvec,mvec,s2vec,qval){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[5]-0.5*(nvec-1)*s2vec*exp(x[5])-0.5*nvec*((mvec-DRM.H4_ND(x[1:4],dvec,qval))^2)*exp(x[5]))
}
#' @rdname llfE4_NI
llfLN4_ND=function(x,nvec,dvec,mvec,s2vec,qval){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[5]-0.5*(nvec-1)*s2vec*exp(x[5])-0.5*nvec*((mvec-DRM.LN4_ND(x[1:4],dvec,qval))^2)*exp(x[5]))
}
#' @rdname llfE4_NI
llfG4_ND=function(x,nvec,dvec,mvec,s2vec,qval){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[5]-0.5*(nvec-1)*s2vec*exp(x[5])-0.5*nvec*((mvec-DRM.G4_ND(x[1:4],dvec,qval))^2)*exp(x[5]))
}
#' @rdname llfE4_NI
llfQE4_ND=function(x,nvec,dvec,mvec,s2vec,qval){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[5]-0.5*(nvec-1)*s2vec*exp(x[5])-0.5*nvec*((mvec-DRM.QE4_ND(x[1:4],dvec,qval))^2)*exp(x[5]))
}
#' @rdname llfE4_NI
llfP4_ND=function(x,nvec,dvec,mvec,s2vec,qval){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[5]-0.5*(nvec-1)*s2vec*exp(x[5])-0.5*nvec*((mvec-DRM.P4_ND(x[1:4],dvec,qval))^2)*exp(x[5]))
}
#' @rdname llfE4_NI
llfL4_ND=function(x,nvec,dvec,mvec,s2vec,qval){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[5]-0.5*(nvec-1)*s2vec*exp(x[5])-0.5*nvec*((mvec-DRM.L4_ND(x[1:4],dvec,qval))^2)*exp(x[5]))
}
#' @rdname llfE4_NI
llfE4_LND=function(x,nvec,dvec,mvec,s2vec,qval){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[5]-0.5*(nvec-1)*s2vec*exp(x[5])-0.5*nvec*((mvec-DRM.E4_LND(x[1:4],dvec,qval))^2)*exp(x[5])) - sum(mvec*nvec)
}
#' @rdname llfE4_NI
llfIE4_LND=function(x,nvec,dvec,mvec,s2vec,qval){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[5]-0.5*(nvec-1)*s2vec*exp(x[5])-0.5*nvec*((mvec-DRM.IE4_LND(x[1:4],dvec,qval))^2)*exp(x[5])) - sum(mvec*nvec)
}
#' @rdname llfE4_NI
llfH4_LND=function(x,nvec,dvec,mvec,s2vec,qval){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[5]-0.5*(nvec-1)*s2vec*exp(x[5])-0.5*nvec*((mvec-DRM.H4_LND(x[1:4],dvec,qval))^2)*exp(x[5])) - sum(mvec*nvec)
}
#' @rdname llfE4_NI
llfLN4_LND=function(x,nvec,dvec,mvec,s2vec,qval){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[5]-0.5*(nvec-1)*s2vec*exp(x[5])-0.5*nvec*((mvec-DRM.LN4_LND(x[1:4],dvec,qval))^2)*exp(x[5])) - sum(mvec*nvec)
}
#' @rdname llfE4_NI
llfG4_LND=function(x,nvec,dvec,mvec,s2vec,qval){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[5]-0.5*(nvec-1)*s2vec*exp(x[5])-0.5*nvec*((mvec-DRM.G4_LND(x[1:4],dvec,qval))^2)*exp(x[5])) - sum(mvec*nvec)
}
#' @rdname llfE4_NI
llfQE4_LND=function(x,nvec,dvec,mvec,s2vec,qval){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[5]-0.5*(nvec-1)*s2vec*exp(x[5])-0.5*nvec*((mvec-DRM.QE4_LND(x[1:4],dvec,qval))^2)*exp(x[5])) - sum(mvec*nvec)
}
#' @rdname llfE4_NI
llfP4_LND=function(x,nvec,dvec,mvec,s2vec,qval){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[5]-0.5*(nvec-1)*s2vec*exp(x[5])-0.5*nvec*((mvec-DRM.P4_LND(x[1:4],dvec,qval))^2)*exp(x[5])) - sum(mvec*nvec)
}
#' @rdname llfE4_NI
llfL4_LND=function(x,nvec,dvec,mvec,s2vec,qval){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[5]-0.5*(nvec-1)*s2vec*exp(x[5])-0.5*nvec*((mvec-DRM.L4_LND(x[1:4],dvec,qval))^2)*exp(x[5])) - sum(mvec*nvec)
}





#' Loglikelihood functions for the different dose response models
#'
#' @param x parameter values
#' @param nvec number of observations
#' @param dvec dose levels
#' @param yvec response
#' @param qval BMR
#'
#' @return loglikelihood value at x.
#'

### Binomial
llfE4_Q=function(x,nvec,dvec,yvec,qval){
  sum(lchoose(nvec, yvec) + yvec*log(DRM.E4_Q(x[1:3], dvec, qval) + .Machine$double.xmin) +
        (nvec - yvec)*log(1 - DRM.E4_Q(x, dvec, qval) + .Machine$double.xmin))
}
#' @rdname llfE4_Q
llfIE4_Q=function(x,nvec,dvec,yvec,qval){
  sum(lchoose(nvec, yvec) + yvec*log(DRM.IE4_Q(x[1:3], dvec, qval) + .Machine$double.xmin) +
        (nvec - yvec)*log(1 - DRM.IE4_Q(x, dvec, qval) + .Machine$double.xmin))
}
#' @rdname llfE4_Q
llfH4_Q=function(x,nvec,dvec,yvec,qval){
  sum(lchoose(nvec, yvec) + yvec*log(DRM.H4_Q(x[1:3], dvec, qval) + .Machine$double.xmin) +
        (nvec - yvec)*log(1 - DRM.H4_Q(x, dvec, qval) + .Machine$double.xmin))
}
#' @rdname llfE4_Q
llfLN4_Q=function(x,nvec,dvec,yvec,qval){
  sum(lchoose(nvec, yvec) + yvec*log(DRM.LN4_Q(x[1:3], dvec, qval) + .Machine$double.xmin) +
        (nvec - yvec)*log(1 - DRM.LN4_Q(x, dvec, qval) + .Machine$double.xmin))
}
#' @rdname llfE4_Q
llfG4_Q=function(x,nvec,dvec,yvec,qval){
  sum(lchoose(nvec, yvec) + yvec*log(DRM.G4_Q(x[1:3], dvec, qval) + .Machine$double.xmin) +
        (nvec - yvec)*log(1 - DRM.G4_Q(x, dvec, qval) + .Machine$double.xmin))
}
#' @rdname llfE4_Q
llfQE4_Q=function(x,nvec,dvec,yvec,qval){
  sum(lchoose(nvec, yvec) + yvec*log(DRM.QE4_Q(x[1:3], dvec, qval) + .Machine$double.xmin) +
        (nvec - yvec)*log(1 - DRM.QE4_Q(x, dvec, qval) + .Machine$double.xmin))
}
#' @rdname llfE4_Q
llfP4_Q=function(x,nvec,dvec,yvec,qval){
  sum(lchoose(nvec, yvec) + yvec*log(DRM.P4_Q(x[1:3], dvec, qval) + .Machine$double.xmin) +
        (nvec - yvec)*log(1 - DRM.P4_Q(x, dvec, qval)))
}
#' @rdname llfE4_Q
llfL4_Q=function(x,nvec,dvec,yvec,qval){
  sum(lchoose(nvec, yvec) + yvec*log(DRM.L4_Q(x[1:3], dvec, qval) + .Machine$double.xmin) +
        (nvec - yvec)*log(1 - DRM.L4_Q(x, dvec, qval) + .Machine$double.xmin))
}
#' @rdname llfE4_Q
llfSM_Q=function(x,nvec,dvec,yvec,qval){
  sum(lchoose(nvec, yvec) + yvec*log(x[1:length(dvec)] + .Machine$double.xmin) +
        (nvec - yvec)*log(1 - x[1:length(dvec)] + .Machine$double.xmin))
}


### Beta-Binomial

#' Loglikelihood functions for the different dose response models (Beta-Binomial)
#'
#' @param x parameter values
#' @param nvec number of observations
#' @param dvec dose levels
#' @param yvec response
#' @param qval BMR
#' @param rho intra-cluster correlation parameter
#'
#' @return loglikelihood value at x.

llfE42_Q=function(x,nvec,dvec,yvec,qval,rho){
  m = DRM.E4_Q(x[1:3], dvec, qval)
  abet = m*((1/rho)-1)
  bbet = (1.0 - m)*((1/rho)-1)

  sum(lchoose(nvec, yvec) + lgamma(abet+yvec+.Machine$double.xmin) + lgamma(bbet+nvec-yvec+.Machine$double.xmin) -
        lgamma(abet+bbet+nvec+.Machine$double.xmin) - lgamma(abet+.Machine$double.xmin) -
        lgamma(bbet+.Machine$double.xmin) +
        lgamma(abet+bbet+.Machine$double.xmin))
}

#' @rdname llfE42_Q
llfIE42_Q=function(x,nvec,dvec,yvec,qval,rho){

  m = DRM.IE4_Q(x[1:3], dvec, qval)
  abet = m*((1/rho)-1)
  bbet = (1.0 - m)*((1/rho)-1)

  sum(lchoose(nvec, yvec) + lgamma(abet+yvec+.Machine$double.xmin) + lgamma(bbet+nvec-yvec+.Machine$double.xmin) -
        lgamma(abet+bbet+nvec+.Machine$double.xmin) - lgamma(abet+.Machine$double.xmin) -
        lgamma(bbet+.Machine$double.xmin) +
        lgamma(abet+bbet+.Machine$double.xmin))
}

#' @rdname llfE42_Q
llfH42_Q=function(x,nvec,dvec,yvec,qval,rho){
  m = DRM.H4_Q(x[1:3], dvec, qval)
  abet = m*((1/rho)-1)
  bbet = (1.0 - m)*((1/rho)-1)

  sum(lchoose(nvec, yvec) + lgamma(abet+yvec+.Machine$double.xmin) + lgamma(bbet+nvec-yvec+.Machine$double.xmin) -
        lgamma(abet+bbet+nvec+.Machine$double.xmin) - lgamma(abet+.Machine$double.xmin) -
        lgamma(bbet+.Machine$double.xmin) +
        lgamma(abet+bbet+.Machine$double.xmin))
}

#' @rdname llfE42_Q
llfLN42_Q=function(x,nvec,dvec,yvec,qval,rho){
  m = DRM.LN4_Q(x[1:3], dvec, qval)
  abet = m*((1/rho)-1)
  bbet = (1.0 - m)*((1/rho)-1)

  sum(lchoose(nvec, yvec) + lgamma(abet+yvec+.Machine$double.xmin) + lgamma(bbet+nvec-yvec+.Machine$double.xmin) -
        lgamma(abet+bbet+nvec+.Machine$double.xmin) - lgamma(abet+.Machine$double.xmin) -
        lgamma(bbet+.Machine$double.xmin) +
        lgamma(abet+bbet+.Machine$double.xmin))
}

#' @rdname llfE42_Q
llfG42_Q=function(x,nvec,dvec,yvec,qval,rho){
  m = DRM.G4_Q(x[1:3], dvec, qval)
  abet = m*((1/rho)-1)
  bbet = (1.0 - m)*((1/rho)-1)

  sum(lchoose(nvec, yvec) + lgamma(abet+yvec+.Machine$double.xmin) + lgamma(bbet+nvec-yvec+.Machine$double.xmin) -
        lgamma(abet+bbet+nvec+.Machine$double.xmin) - lgamma(abet+.Machine$double.xmin) -
        lgamma(bbet+.Machine$double.xmin) +
        lgamma(abet+bbet+.Machine$double.xmin))
}


#' @rdname llfE42_Q
llfQE42_Q=function(x,nvec,dvec,yvec,qval,rho){
  m = DRM.QE4_Q(x[1:3], dvec, qval)
  abet = m*((1/rho)-1)
  bbet = (1.0 - m)*((1/rho)-1)

  sum(lchoose(nvec, yvec) + lgamma(abet+yvec+.Machine$double.xmin) + lgamma(bbet+nvec-yvec+.Machine$double.xmin) -
        lgamma(abet+bbet+nvec+.Machine$double.xmin) - lgamma(abet+.Machine$double.xmin) -
        lgamma(bbet+.Machine$double.xmin) +
        lgamma(abet+bbet+.Machine$double.xmin))
}

llfP42_Q=function(x,nvec,dvec,yvec,qval,rho){
  m = DRM.P4_Q(x[1:3], dvec, qval)
  abet = m*((1/rho)-1)
  bbet = (1.0 - m)*((1/rho)-1)

  sum(lchoose(nvec, yvec) + lgamma(abet+yvec+.Machine$double.xmin) + lgamma(bbet+nvec-yvec+.Machine$double.xmin) -
        lgamma(abet+bbet+nvec+.Machine$double.xmin) - lgamma(abet+.Machine$double.xmin) -
        lgamma(bbet+.Machine$double.xmin) +
        lgamma(abet+bbet+.Machine$double.xmin))
}


#' @rdname llfE42_Q
llfL42_Q=function(x,nvec,dvec,yvec,qval,rho){
  m = DRM.IE4_Q(x[1:3], dvec, qval)
  abet = m*((1/rho)-1)
  bbet = (1.0 - m)*((1/rho)-1)

  sum(lchoose(nvec, yvec) + lgamma(abet+yvec+.Machine$double.xmin) + lgamma(bbet+nvec-yvec+.Machine$double.xmin) -
        lgamma(abet+bbet+nvec+.Machine$double.xmin) - lgamma(abet+.Machine$double.xmin) -
        lgamma(bbet+.Machine$double.xmin) +
        lgamma(abet+bbet+.Machine$double.xmin))
}

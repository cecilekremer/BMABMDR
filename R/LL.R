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

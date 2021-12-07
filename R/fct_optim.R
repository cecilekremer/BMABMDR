#' Function for internal use
#'
#' @param mod stan model
#' @param data input data
#' @param stv start values
#' @param ndraws ndraws
#' @param seed seed
#' @param pvec pvec
#'
#' @return .
#'
fun_optim = function(mod, data, stv,
                        ndraws = ndraws,seed=seed,pvec){

  opt=try(rstan::optimizing(mod,data = data,
                          seed=as.integer(seed),draws=ndraws,
                          init=stv, hessian=TRUE),silent=T)
  if(length(opt)!=9) {

    # svh = try(list(par=as.vector(rstan::optimizing(mod,data = data)$par[1:5])),silent=T)
    opt=try(rstan::optimizing(mod,data = data,seed=as.integer(seed),draws=ndraws,hessian=TRUE),silent=T)
  }
  while(ifelse(is.na(opt[3]),TRUE,(opt[3]!=0)) | length(opt)!=9){
    # print('stuck here')
    # svF2P$par[3]=svF2P$par[3]+runif(1,-0.1,0.1)
    svh=stv
    svh$par=svh$par+runif(5,-0.1,0.1)
    opt=try(rstan::optimizing(mod,data = data,
                            seed=as.integer(seed),draws=ndraws,
                            init=svh,hessian=TRUE),silent=T)
  }

  attr(opt, "class") <- "stanfitOptim"
  return(opt)

}

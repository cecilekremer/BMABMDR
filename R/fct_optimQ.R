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
fun_optimQ = function(mod, data, stv,
                      ndraws = ndraws,seed=seed,pvec){

  opt=try(rstan::optimizing(mod,data = data,
                            seed=as.integer(seed),draws=ndraws,
                            init=stv, hessian=TRUE),silent=T)
  if(length(opt)!=9) {

    # svh = try(list(par=as.vector(rstan::optimizing(mod,data = data)$par[1:5])),silent=T)
    opt=try(rstan::optimizing(mod,data = data,seed=as.integer(seed),draws=ndraws,hessian=TRUE),silent=T)
  }

  while(ifelse(is.na(opt[3]),TRUE,(opt[3]!=0)) | length(opt)!=9){
    svh=stv
    par = unname(unlist(svh))

    if(data$is_informative_a == 1){
      par[1] = runif(1, data$priorlb[1], data$priorub[1])
    }else{
      par[1] = par[1] + rnorm(1, sd = 0.01*abs(par[1]))
    }
    if(data$is_informative_BMD == 1){
      par[2] = runif(1, data$priorlb[2], data$priorub[2])
    }else{
      par[2] = par[2] + rnorm(1, sd = 0.01*abs(par[2]))
    }
    par[3] = par[3] + rnorm(1, sd = 0.01*abs(par[3]))

    opt=try(rstan::optimizing(mod,data = data,
                              seed=as.integer(seed),draws=ndraws,
                              init=svh,hessian=TRUE),silent=T)
  }

  attr(opt, "class") <- "stanfitOptim"
  return(opt)

}

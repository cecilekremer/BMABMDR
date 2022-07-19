#'  function to perform Laplace approximation for a given dose-response model.
#'
#' @param mod stan model
#' @param data list containing data values to be passed to the stan model file
#' @param stv list of starting values
#' @param ndraws number of draws to be made from the posterior distribution
#' @param seed random seed for reproducibility
#' @param pvec probability vector to compute credible interval for the BMD
#'
#' @examples
#'
#' @return a stan results object
#'
#' @export fun_optim
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
    svh=stv
    par = unname(unlist(svh))

    if(data$is_informative_a == 1){
      par[1] = runif(1, data$priorlb[1], data$priorub[1])
    }else{
      par[1] = par[1] + rnorm(1, sd = 0.01*abs(par[1]))
    }
    if(data$is_informative_c == 1){
      if(data$is_increasing == 1){
        par[3] = runif(1, data$priorlb[3] - data$L, data$priorub[3] - data$L)
      }else{
        par[3] = runif(1, data$priorlb[3] / data$U, data$priorub[3] / data$U)
      }
    }else{
      par[3] = par[3] + rnorm(1, sd = 0.01*abs(par[3]))
    }

    ## for decreasing, par[3] < 1
    if(data$is_decreasing == 1 & par[3] > 1){
      par[3] = 0.9999
    }

    if(data$is_informative_BMD == 1){
      par[2] = runif(1, data$priorlb[2], data$priorub[2])
    }else{
      par[2] = par[2] + rnorm(1, sd = 0.01*abs(par[2]))
    }
    par[4:5] = par[4:5] + rnorm(2, sd = 0.01*abs(par[4:5]))

    ## BMD between 0 and 1
    if(par[2] > 1) par[2] = 0.9

    if(data$data_type==1|data$data_type==2){
      pars3d = numeric()
      pars3i = par[3]
      dim(pars3d)=0
      dim(pars3i)=1
      svh = list(par1=par[1],par2=par[2],pars3i=pars3i,pars3d=pars3d,par4=par[4],par5=par[5])
    }else if(data$data_type==3|data$data_type==4){
      pars3d = par[3]
      pars3i = numeric()
      dim(pars3d)=1
      dim(pars3i)=0
      svh = list(par1=par[1],par2=par[2],pars3i=pars3i,pars3d=pars3d,par4=par[4],par5=par[5])
    }
    opt=try(rstan::optimizing(mod,data = data,
                              seed=as.integer(seed),draws=ndraws,
                              init=svh,hessian=TRUE),silent=T)
  }

  attr(opt, "class") <- "stanfitOptim"
  return(opt)

}

#'  Function to perform MCMC sampling for a given dose-response model.
#'
#' @param mod stan model
#' @param data list containing data values to be passed to the stan model file
#' @param stv list of starting values
#' @param ndraws number of draws to be made from the posterior distribution
#' @param nrchains number of MCMC chains
#' @param nriterations number of MCMC iterations
#' @param warmup  number of MCMC iterations for warmup
#' @param delta adapt_delta value for the HMC in stan. See \code{\link[rstan]{sampling}} for more.
#' @param treedepth tree_depth value for the HMC in stan. See \code{\link[rstan]{sampling}} for more
#' @param seed random seed for reproducibility
#' @param pvec probability vector to compute credible interval for the BMD
#'
#' @examples
#'
#' @return a stan results object
#'
#' @export fun_sampling
#'
fun_sampling = function(mod, data, stv,
                        ndraws = ndraws,nrchains=nrchains,
                        nriterations=nriterations,warmup=warmup,
                        delta=delta,treedepth=treedepth,seed=seed,pvec){

  if(exists("opt")) rm(opt)
  opt = try(rstan::optimizing(mod,data = data,init=stv), silent = T)
  if(ifelse(is.na(opt[3]),TRUE,(opt[3]!=0))){
    opt = try(rstan::optimizing(mod, data = data), silent = T)
  }

  n.attempts <- 1
  while(ifelse(is.na(opt[3]),TRUE,(opt[3]!=0)) & n.attempts < 100){
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
    if(data$is_informative_BMD == 1){
      par[2] = runif(1, data$priorlb[2], data$priorub[2])
    }else{
      par[2] = par[2] + rnorm(1, sd = 0.01*abs(par[2]))
    }
    par[4:5] = par[4:5] + rnorm(2, sd = 0.01*abs(par[4:5]))

    ## for decreasing, par[3] < 1
    if(data$is_decreasing == 1 & par[3] > 1){
      par[3] = 0.9999
    }
    ## for increasing, par[3] > 0

    ## BMD between 0 and 1
    # if(par[2] > 1) par[2] = 0.9


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
                              seed=as.integer(seed),init=svh),silent=T)

    n.attempts <- n.attempts + 1

  }
  sv = opt$par

  initf2 <- function(chain_id = 1) {
    par = sv[1:5]
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
    # if(par[2] > 1) par[2] = 0.9

    if(data$data_type==1|data$data_type==2){
      pars3d = numeric()
      pars3i = par[3]
      dim(pars3d)=0
      dim(pars3i)=1
      list(par1=par[1],par2=par[2],pars3i=pars3i,pars3d=pars3d,par4=par[4],par5=par[5],
           alpha = chain_id)
    }else if(data$data_type==3|data$data_type==4){
      pars3d = par[3]
      pars3i = numeric()
      dim(pars3d)=1
      dim(pars3i)=0
      list(par1=par[1],par2=par[2],pars3i=pars3i,pars3d=pars3d,par4=par[4],par5=par[5],
           alpha = chain_id)
    }
  }
  init_ll <- lapply(1:nrchains, function(id) initf2(chain_id = id))
  fitstan <- rstan::sampling(mod, data = data,
                             init=init_ll,iter = nriterations,
                             chains = nrchains,warmup=warmup,
                             seed=123,
                             control = list(adapt_delta = delta,
                                            max_treedepth =treedepth)
                             # )
                             ,
                             show_messages = F, refresh = 0, verbose = F)

  if(is.na(dim(fitstan)[1])){

    fitstan <- NULL

  }


  return(fitstan)

}

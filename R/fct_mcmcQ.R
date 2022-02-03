#' Function for internal use
#'
#' @param mod stan model
#' @param data input data
#' @param stv start values
#' @param ndraws ndraws
#' @param nrchains nrchains
#' @param nriterations nriterations
#' @param warmup burnin
#' @param delta delta
#' @param treedepth treedepth
#' @param seed seed
#' @param pvec pvec
#'
#' @return .
#'
fun_samplingQ = function(mod, data, stv,
                         ndraws = ndraws,nrchains=nrchains,
                         nriterations=nriterations,warmup=warmup,
                         delta=delta,treedepth=treedepth,seed=seed,pvec){

  if(exists("opt")) rm(opt)
  opt = try(rstan::optimizing(mod,data = data,init=stv), silent = T)
  if(ifelse(is.na(opt[3]),TRUE,(opt[3]!=0))){
    opt = try(rstan::optimizing(mod, data = data), silent = T)
  }
  while(ifelse(is.na(opt[3]),TRUE,(opt[3]!=0))){
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
                              seed=as.integer(seed),init=svh),silent=T)
  }

  if(data$is_bin == 1) {

    sv <- opt$par[names(opt$par) %in% c('par1', 'par2', 'par3')]
    initf2 <- function(chain_id = 1) {
      par = sv[1:3]
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

      list(par1=par[1],par2=par[2],par3=par[3],
           alpha = chain_id)
    }

  } else if(data$is_betabin == 1) {

    sv <- opt$par[names(opt$par) %in% c('par1', 'par2', 'par3', 'etarho[1]')]
    initf2 <- function(chain_id = 1) {
      par = sv[1:4]
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
      par[4] = par[4] + rnorm(1, sd = 0.01*abs(par[4]))
      etarho = par[4]; dim(etarho)=1
      list(par1=par[1],par2=par[2],par3=par[3], etarho = etarho,
           alpha = chain_id)
    }
  }



  init_ll <- lapply(1:nrchains, function(id) initf2(chain_id = id))
  fitstan <- rstan::sampling(mod, data = data,
                             init=init_ll,iter = nriterations,
                             chains = nrchains,warmup=warmup,
                             seed=123,
                             control = list(adapt_delta = delta,
                                            max_treedepth =treedepth),
                             show_messages = F, refresh = 0)

  while(is.na(dim(fitstan)[1])){

    init_ll <- lapply(1:nrchains, function(id) initf2(chain_id = id))

    fitstan=try(rstan::sampling(mod,data = data,
                                init=init_ll,iter = nriterations,
                                chains = nrchains,warmup=warmup,seed=123,
                                control = list(adapt_delta = delta,
                                               max_treedepth =treedepth),
                                show_messages = F, refresh = 0),silent=T)
  }


  return(fitstan)

}

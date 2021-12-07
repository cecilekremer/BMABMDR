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
fun_sampling = function(mod, data, stv,
                        ndraws = ndraws,nrchains=nrchains,
                        nriterations=nriterations,warmup=warmup,
                        delta=delta,treedepth=treedepth,seed=seed,pvec){

  if(exists("opt")) rm(opt)
  opt = try(rstan::optimizing(mod,data = data,init=stv), silent = T)
  if(ifelse(is.na(opt[3]),TRUE,(opt[3]!=0))){
    opt = try(rstan::optimizing(mod, data = data), silent = T)
  }
  while(ifelse(is.na(opt[3]),TRUE,(opt[3]!=0))){
    svF1h = stv
    svF1h$par=stv$par+runif(5,-0.1,0.1)
    opt=try(rstan::optimizing(mod,data = data,
                       seed=as.integer(seed),init=svF1h),silent=T)
  }
  sv = opt$par

  initf2 <- function(chain_id = 1) {
    list(par=sv[1:5]+rnorm(5,sd=0.01*abs(sv[1:5])),alpha = chain_id)
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
                      show_messages = F, refresh = 0)

  while(is.na(dim(fitstan)[1])){

    initf2 <- function(chain_id = 1) {
      list(par=sv[1:5]+rnorm(5,sd=0.01*abs(sv[1:5])),alpha = chain_id)
    }
    init_ll <- lapply(1:nrchains, function(id) initf2(chain_id = id))

    # save(init_ll, file = "init.RData")
    fitstan=try(rstan::sampling(mod,data = data,
                         init=init_ll,iter = nriterations,
                         chains = nrchains,warmup=warmup,seed=123,
                         control = list(adapt_delta = delta,
                                        max_treedepth =treedepth),
                         show_messages = F, refresh = 0),silent=T)
  }


  return(fitstan)

}

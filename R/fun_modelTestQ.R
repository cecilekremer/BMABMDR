#' Test whether the best fit model fits equally well as the saturated ANOVA model
#'
#'
#' More detailed descriprion
#' @param best.fit name of the best fitting model
#' @param data.N the input data as returned by function PREP_DATA_N
#' @param data.LN the input data as returned by function PREP_DATA_LN
#' @param stanBest fit of the best fitting model as returned by rstan's optimizing or sampling routine
#' @param type either 'MCMC' or 'Laplace'
#' @param seed default 123
#' @param ndraws the number of draws from the posterior, default 30000
#' @param nrchains the number of chains to be used in the MCMC
#' @param nriterations the number of iterations per chain
#' @param warmup the number of iterations per chain to be discarded as burnin
#' @param delta default 0.8
#' @param treedepth default 10
#'
#' @return Bayes factor for the best fitting model compared to the saturated ANOVA model. A Bayes factor < 10 indicates the best fitting model fits equally well as the ANOVA model.
#'
#' @export
#'
modelTestQ <- function(best.fit, data.Q, stanBest, type, seed, ndraws, nrchains, nriterations, warmup, delta, treedepth){

  N = data.Q$data$N

  priorSM = list(
    priormu = c(max(c(data.Q$data$y[1]/data.Q$data$n[1], 1/(5*data.Q$data$n[1]))),
                diff(data.Q$data$y/data.Q$data$n)),
    priorlb = ifelse(data.Q$data$y[1] != 0, max(c(prop.test(data.Q$data$y[1], data.Q$data$n[1])$conf.int[1]/2, 1/(10*data.Q$data$n[1]))),
                     .Machine$double.xmin),
    priorub = c(min(c(3*prop.test(data.Q$data$y[1], data.Q$data$n[1])$conf.int[2]/2, 1 - 1/(10*data.Q$data$n[1]))),
                max(abs(diff(data.Q$data$y/data.Q$data$n)))
    )
  )

  svSM = list(par = c(max(c(data.Q$data$y[1]/data.Q$data$n[1], 1/(5*data.Q$data$n[1]))), # background
                      diff(data.Q$data$y/data.Q$data$n)) # invsigma2
  )

  data.modstanSM = list(N=N,n=data.Q$data$n,y=data.Q$data$y,
                        priormu=priorSM$priormu,
                        priorlb=priorSM$priorlb, priorub=priorSM$priorub,
                        is_bin=1, is_betabin = 0, priorgama = 4, eps = .Machine$double.xmin
  )

  if(type == 'MCMC'){

    sv=rstan::optimizing(stanmodels$mSM_Q,data = data.modstanSM,init=svSM)$par

    initf2 <- function(chain_id = 1) {
      list(par=sv[1:(data.Q$data$N)] + rnorm(data.Q$data$N, sd = 0.01*abs(sv[1:(data.Q$data$N)])) ,alpha = chain_id)
    }
    init_ll <- lapply(1:nrchains, function(id) initf2(chain_id = id))
    fitstanSM = rstan::sampling(stanmodels$mSM_Q, data = data.modstanSM, init=init_ll, iter = nriterations,
                                chains = nrchains, warmup = warmup, seed = seed,
                                control = list(adapt_delta = delta, max_treedepth =treedepth),
                                show_messages = F, refresh = 0)

    while(is.na(dim(fitstanSM)[1])){

      init_ll <- lapply(1:nrchains, function(id) initf2(chain_id = id))

      fitstanSM = rstan::sampling(stanmodels$mSM_Q, data = data.modstanSM, init=init_ll, iter = nriterations,
                                  chains = nrchains, warmup = warmup, seed = seed,
                                  control = list(adapt_delta = delta, max_treedepth =treedepth),
                                  show_messages = F, refresh = 0)
    }

    pars.bestfit = apply(as.matrix(stanBest),2,median)[c("par1","par2","par3")]

    parsSM = as.matrix(fitstanSM)
    means.SM = apply(parsSM[, paste0('a[', 1:data.Q$data$N, ']')], 2, median)
    pars.SM = apply(parsSM[, c(paste0('a[', 1:data.Q$data$N, ']'), paste0('par[', data.Q$data$N, ']'))], 2, median)

  }else if(type == 'Laplace'){

    all.pars.bestfit = parq_extract(stanBest, model_name = paste0(best.fit,'_Q'), pars = c('a', 'b', 'd', 'BMD',
                                                                                           paste0('par',1:3)))
    pars.bestfit = apply(all.pars.bestfit[,c(paste0("p",1:3))], 2, median)

    optSM = optimizing(stanmodels$mSM_Q, data = data.modstanSM,
                       seed=as.integer(seed), draws = ndraws,
                       init = svSM, hessian=TRUE)
    pars.SM = apply(optSM$theta_tilde[, c(paste0('a[', 1:data.Q$data$N, ']'), paste0('par[', data.Q$data$N, ']'))], 2, median)
    means.SM = apply(optSM$theta_tilde[, paste0('a[', 1:data.Q$data$N, ']')], 2, median)

  }

  if(data.Q$data$is_bin == 1){
    llfun = paste0('llf',best.fit,'_Q')
  }else if(data.Q$data$is_betabin == 1){
    llfun = paste0('llf',best.fit,'2_Q')
  }
  llBestfitf = get(llfun)
  llBestfit = llBestfitf(x = pars.bestfit, data.Q$data$n, data.Q$data$x, data.Q$data$y, data.Q$data$q)

  if(data.Q$data$is_bin == 1){
    llSM = llfSM_Q(pars.SM, data.Q$data$n, data.Q$data$x, data.Q$data$y)
  }else if(data.Q$data$is_betabin == 1){
    llSM = llfSM2_Q(pars.SM, data.Q$data$n, data.Q$data$x, data.Q$data$y, rho) ## UPDATE --> rho? !!
  }

  BIC.bestfit = - 2 * llBestfit + (3 * log(sum(data.Q$data$n)))
  BIC.SM = - 2 * llSM + ((data.Q$data$N) * log(sum(data.Q$data$n)))

  bf = exp(-0.5 * (BIC.bestfit - BIC.SM))

  if(bf < 1/10){
    warn.bf = 'None of the models provide an adequate fit do the data.'
  }else if(bf > 1/10 & bf < 10){
    warn.bf = 'Best fitting model fits well.'
  }else if(bf > 10){
    warn.bf = 'attention: bayes factor is larger than 10 in favor of the best fitting model'
  }

  return(list(bayesFactor = bf,
              means.SM = means.SM,
              par.best = pars.bestfit,
              BIC.bestfit = BIC.bestfit,
              BIC.SM = BIC.SM,
              warn.bf = warn.bf)
  )

}

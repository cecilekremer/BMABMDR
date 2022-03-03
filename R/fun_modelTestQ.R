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
#' @export modelTestQ
#'
modelTestQ <- function(best.fit, data.Q, stanBest, type, seed, ndraws, nrchains, nriterations,
                       warmup, delta, treedepth){

  N = data.Q$data$N

  ydiff <- abs(diff(data.Q$data$y/data.Q$data$n))

  tts <- sapply(1:(length(data.Q$data$y)-1), function(i){
    tt <- prop.test(c(data.Q$data$y[i+1], data.Q$data$y[i]),
                    c(data.Q$data$n[i+1], data.Q$data$n[i]))$conf.int
    return(c(lbs = tt[1], ubs = tt[2]))
  })
  lbs <- as.numeric(tts[1,])
  ubs <- as.numeric(tts[2,])

  pls <- which(lbs==0 & ubs==0)
  lbs[pls] <- min(lbs, na.rm = TRUE)
  ubs[pls] <- max(ubs, na.rm = TRUE)

  if(data.Q$data$is_betabin == 1){

    datf = data.frame(yy = data.Q$data$y, n.a = data.Q$data$n, xx = data.Q$data$x)
    fpfit2 <- try(gamlss(cbind(yy,n.a-yy)~as.factor(xx), sigma.formula=~1, family=BB, data=datf),
                  silent = TRUE)
    rhohat <- exp(fpfit2$sigma.coefficients)/(exp(fpfit2$sigma.coefficients)+1)
    dim(rhohat) <- 1

    priorSM = list(
      priormu = c(max(c(data.Q$data$y[1]/data.Q$data$n[1], 1/(5*data.Q$data$n[1]))),
                  diff(data.Q$data$y/data.Q$data$n), rhohat),
      priorlb = c(ifelse(data.Q$data$y[1] != 0, max(c(prop.test(data.Q$data$y[1],
                                                                data.Q$data$n[1])$conf.int[1]/2,
                                                      1/(10*data.Q$data$n[1]))),
                         .Machine$double.xmin), lbs),
      priorub = c(min(c(3*prop.test(data.Q$data$y[1], data.Q$data$n[1])$conf.int[2]/2, 1 - 1/(10*data.Q$data$n[1]))),
                  ubs
      )
    )

    svSM = list(par = c(max(c(data.Q$data$y[1]/data.Q$data$n[1], 1/(5*data.Q$data$n[1]))), # background
                        diff(data.Q$data$y/data.Q$data$n)), rho = rhohat # invsigma2
    )

    data.modstanSM = list(N=N,n=data.Q$data$n,y=data.Q$data$y,yint = data.Q$data$y,nint = data.Q$data$n,
                          priormu=priorSM$priormu,
                          priorlb=priorSM$priorlb, priorub=priorSM$priorub,
                          is_bin=0, is_betabin = 1, priorgama = 4, eps = .Machine$double.xmin
    )
  } else if(data.Q$data$is_bin == 1){
    priorSM = list(
      priormu = c(max(c(data.Q$data$y[1]/data.Q$data$n[1], 1/(5*data.Q$data$n[1]))),
                  diff(data.Q$data$y/data.Q$data$n), 0.0),
      priorlb = c(ifelse(data.Q$data$y[1] != 0, max(c(prop.test(data.Q$data$y[1],
                                                                data.Q$data$n[1])$conf.int[1]/2,
                                                      1/(10*data.Q$data$n[1]))),
                         .Machine$double.xmin), lbs),
      priorub = c(min(c(3*prop.test(data.Q$data$y[1], data.Q$data$n[1])$conf.int[2]/2, 1 - 1/(10*data.Q$data$n[1]))),
                  ubs
      )
    )

    svSM = list(par = c(max(c(data.Q$data$y[1]/data.Q$data$n[1], 1/(5*data.Q$data$n[1]))), # background
                        diff(data.Q$data$y/data.Q$data$n)) # invsigma2
    )

    data.modstanSM = list(N=N,n=data.Q$data$n,y=data.Q$data$y,yint = data.Q$data$y,nint = data.Q$data$n,
                          priormu=priorSM$priormu,
                          priorlb=priorSM$priorlb, priorub=priorSM$priorub,
                          is_bin=1, is_betabin = 0, priorgama = 4, eps = .Machine$double.xmin
    )
  } else stop("data must be either clustered or independent")

  #max(abs(diff(data.Q$data$y/data.Q$data$n))),


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
                                  control = list(adapt_delta = delta, max_treedepth = treedepth),
                                  show_messages = F, refresh = 0)
    }

    if(data.Q$data$is_bin == 1){
      pars.bestfit = apply(as.matrix(stanBest),2,median)[c("par1","par2","par3")]
    } else if(data.Q$data$is_betabin == 1){
      pars.bestfit = apply(as.matrix(stanBest),2,median)[c("par1","par2","par3", "rho[1]")]
    } else stop("data must be either clustered or independent")

    parsSM = as.matrix(fitstanSM)
    means.SM = apply(parsSM[, c(paste0('a[', 1:data.Q$data$N, ']'), "rho[1]")], 2, median)
    pars.SM = apply(parsSM[, c(paste0('a[', 1:data.Q$data$N, ']'),
                               paste0('par[', data.Q$data$N, ']'), "rho[1]")], 2, median)

  }else if(type == 'Laplace'){

    if(data.Q$data$is_bin == 1){
      all.pars.bestfit = parq_extract(stanBest, model_name = paste0(best.fit,'_Q'), pars = c('a', 'b', 'd', 'BMD',
                                                                                             paste0('par',1:3)))
      pars.bestfit = apply(all.pars.bestfit[,c(paste0("p",1:3))], 2, median)

      optSM = optimizing(stanmodels$mSM_Q, data = data.modstanSM,
                         seed=as.integer(seed), draws = ndraws,
                         init = svSM, hessian=TRUE)
      pars.SM = apply(optSM$theta_tilde[, c(paste0('a[', 1:data.Q$data$N, ']'),
                                            paste0('par[', data.Q$data$N, ']'))], 2, median)
      means.SM = apply(optSM$theta_tilde[, paste0('a[', 1:data.Q$data$N, ']')], 2, median)
    } else if(data.Q$data$is_betabin == 1) {
      all.pars.bestfit = parq_extract(stanBest, model_name = paste0(best.fit,'_Q'),
                                      pars = c('a', 'b', 'd', 'rho[1]','BMD', paste0('par',1:3)),
                                      rho = TRUE)
      pars.bestfit = apply(all.pars.bestfit[,c(paste0("p",1:3),"rho")], 2, median)

      optSM = optimizing(stanmodels$mSM_Q, data = data.modstanSM,
                         seed=as.integer(seed), draws = ndraws,
                         init = svSM, hessian=TRUE)
      pars.SM = apply(optSM$theta_tilde[, c(paste0('a[', 1:data.Q$data$N, ']'),
                                            paste0('par[', data.Q$data$N, ']'),'rho[1]')], 2, median)
      means.SM = apply(optSM$theta_tilde[, c(paste0('a[', 1:data.Q$data$N, ']'),'rho[1]')], 2, median)
    }

  }

  if(data.Q$data$is_bin == 1){
    llfun = paste0('llf',best.fit,'_Q')
  }else if(data.Q$data$is_betabin == 1){
    llfun = paste0('llf',best.fit,'22_Q')
  }
  llBestfitf = get(llfun)
  llBestfit = llBestfitf(x = pars.bestfit, data.Q$data$n, data.Q$data$x, data.Q$data$y, data.Q$data$q)

  if(data.Q$data$is_bin == 1){
    llSM = llfSM_Q(pars.SM, data.Q$data$n, data.Q$data$x, data.Q$data$y)
  }else if(data.Q$data$is_betabin == 1){
    llSM = llfSM2_Q(pars.SM, data.Q$data$n, data.Q$data$x, data.Q$data$y,
                    pars.SM[stringr::str_detect(names(pars.SM), 'rho')]) ## UPDATE --> rho? !!
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

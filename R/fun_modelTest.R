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
modelTest <- function(best.fit, data.N, data.LN, stanBest, type, seed, ndraws, nrchains, nriterations, warmup, delta, treedepth){

  if(grepl('_N', best.fit)){

    svSM=list(par = c(data.N$data$m[1], # background
                      diff(data.N$data$m), # increments
                      log(1/mean(data.N$data$s2))) # invsigma2
    )

    priorSM = list(
      priormu = c(data.N$data$m[1],
                  diff(data.N$data$m), # increments
                  -2*log(1.5*mean(sqrt(data.N$data$s2)))),
      priorSigma = diag(c(1, rep(1, data.N$data$N-1), 1)),
      priorlb = data.N$data$priorlb[1],
      priorub = c(data.N$data$priorub[1],
                  max(abs(diff(data.N$data$m)))*10
      )
    )

    data.modstanSM=list(N=data.N$data$N, n=data.N$data$n, m=data.N$data$m, s2=data.N$data$s2,
                        priormu=priorSM$priormu, priorSigma=priorSM$priorSigma,
                        priorlb=priorSM$priorlb, priorub=priorSM$priorub,
                        data_type=data.N$data$data_type, priorg = data.N$data$shape.a
    )

    if(type == 'MCMC'){

      sv=rstan::optimizing(stanmodels$mSM,data = data.modstanSM,init=svSM)$par

      initf2 <- function(chain_id = 1) {
        list(par=sv[1:(data.N$data$N+1)] + rnorm(data.N$data$N+1, sd = 0.01*abs(sv[1:(data.N$data$N+1)])) ,alpha = chain_id)
      }
      init_ll <- lapply(1:nrchains, function(id) initf2(chain_id = id))
      fitstanSM = rstan::sampling(stanmodels$mSM, data = data.modstanSM, init=init_ll, iter = nriterations,
                                  chains = nrchains, warmup = warmup, seed = seed,
                                  control = list(adapt_delta = delta, max_treedepth =treedepth),
                                  show_messages = F, refresh = 0)

      pars.bestfit = apply(as.matrix(stanBest),2,median)[c("par1","par2","par3","par4","par5")]

      parsSM = as.matrix(fitstanSM)
      means.SM = apply(parsSM[, paste0('mu[', 1:data.N$data$N, ']')], 2, median)
      pars.SM = apply(parsSM[, c(paste0('a[', 1:data.N$data$N, ']'), paste0('par[', data.N$data$N+1, ']'))], 2, median)

    }else if(type == 'Laplace'){

      all.pars.bestfit = par_extract(stanBest, model_name = best.fit)
      pars.bestfit = apply(all.pars.bestfit[,c(paste0("p",1:4), "is2t")], 2, median)

      optSM = optimizing(stanmodels$mSM, data = data.modstanSM,
                         seed=as.integer(seed), draws = ndraws,
                         init = svSM, hessian=TRUE)
      pars.SM = apply(optSM$theta_tilde[, c(paste0('a[', 1:data.N$data$N, ']'), paste0('par[', data.N$data$N+1, ']'))], 2, median)
      means.SM = apply(optSM$theta_tilde[, paste0('mu[', 1:data.N$data$N, ']')], 2, median)

    }

    if(data.N$data$is_increasing == 1){
      llfun = paste0('llf',best.fit,'I')
    }else if(data.N$data$is_decreasing == 1){
      llfun = paste0('llf',best.fit,'D')
    }
    llBestfitf = get(llfun)
    llBestfit = llBestfitf(x = pars.bestfit, data.N$data$n, data.N$data$x, data.N$data$m, data.N$data$s2, data.N$data$q)

    llSM = llfSM_N(pars.SM, data.N$data$n, data.N$data$x, data.N$data$m, data.N$data$s2)

    BIC.bestfit = - 2 * llBestfit + (5 * log(sum(data.N$data$n)))
    BIC.SM = - 2 * llSM + ((data.N$data$N + 1) * log(sum(data.N$data$n)))

    bf = exp(-0.5 * (BIC.bestfit - BIC.SM))

    if(bf < 1/10){
      warn.bf = 'None of the models provide an adequate fit do the data.'
    }else if(bf > 1/10 & bf < 10){
      warn.bf = 'Best fitting model fits well.'
    }else if(bf > 10){
      warn.bf = 'attention: bayes factor is larger than 10 in favor of the best fitting model'
    }


  #############--LOGNORMAL--##################################

  }else if(grepl('_LN', best.fit)){

    svSM=list(par = c(exp(data.LN$data$m[1]), # background
                      diff(exp(data.LN$data$m)),
                      log(1/mean(data.LN$data$s2))) # invsigma2
    )

    priorSM = list(
      priormu = c(exp(data.LN$data$m[1]),
                  diff(exp(data.LN$data$m)),
                  -2*log(1.5*mean(sqrt(data.LN$data$s2)))),
      priorSigma = diag(c(1, rep(1, data.LN$data$N-1), 1)),
      priorlb = data.LN$data$priorlb[1],
      priorub = c(data.LN$data$priorub[1],
                  max(abs(diff(exp(data.LN$data$m))))*10
      )
    )

    data.modstanSM=list(N=data.LN$data$N, n=data.LN$data$n, m=data.LN$data$m, s2=data.LN$data$s2,
                        priormu=priorSM$priormu, priorSigma=priorSM$priorSigma,
                        priorlb=priorSM$priorlb, priorub=priorSM$priorub,
                        data_type=data.LN$data$data_type, priorg = data.LN$data$shape.a
    )

    if(type == 'MCMC'){

      sv=rstan::optimizing(stanmodels$mSM,data = data.modstanSM,init=svSM)$par

      initf2 <- function(chain_id = 1) {
        list(par=sv[1:(data.LN$data$N+1)] + rnorm(data.LN$data$N+1, sd = 0.01*abs(sv[1:(data.LN$data$N+1)])) ,alpha = chain_id)
      }
      init_ll <- lapply(1:nrchains, function(id) initf2(chain_id = id))
      fitstanSM = rstan::sampling(stanmodels$mSM, data = data.modstanSM, init=init_ll, iter = nriterations,
                                  chains = nrchains, warmup = warmup, seed = seed,
                                  control = list(adapt_delta = delta, max_treedepth =treedepth),
                                  show_messages = F, refresh = 0)

      pars.bestfit = apply(as.matrix(stanBest),2,median)[c("par1","par2","par3","par4","par5")]

      parsSM = as.matrix(fitstanSM)
      means.SM = apply(parsSM[, paste0('mu[', 1:data.LN$data$N, ']')], 2, median)
      pars.SM = apply(parsSM[, c(paste0('a[', 1:data.LN$data$N, ']'), paste0('par[', data.LN$data$N+1, ']'))], 2, median)


    }else if(type == 'Laplace'){

      all.pars.bestfit = par_extract(stanBest, model_name = best.fit)
      pars.bestfit = apply(all.pars.bestfit[,c(paste0("p",1:4), "is2t")], 2, median)

      optSM = optimizing(stanmodels$mSM, data = data.modstanSM,
                         seed=as.integer(seed), draws = ndraws,
                         init = svSM, hessian=TRUE)
      pars.SM = apply(optSM$theta_tilde[, c(paste0('a[', 1:data.LN$data$N, ']'), paste0('par[', data.LN$data$N+1, ']'))], 2, median)
      means.SM = apply(optSM$theta_tilde[, paste0('mu[', 1:data.LN$data$N, ']')], 2, median)

    }

    if(data.LN$data$is_increasing == 1){
      llfun = paste0('llf',best.fit,'I')
    }else if(data.LN$data$is_decreasing == 1){
      llfun = paste0('llf',best.fit,'D')
    }
    llBestfitf = get(llfun)
    llBestfit = llBestfitf(x = pars.bestfit, data.LN$data$n, data.LN$data$x, data.LN$data$m, data.LN$data$s2, data.LN$data$q)

    llSM = llfSM_LN(x = pars.SM, data.LN$data$n, data.LN$data$x, data.LN$data$m, data.LN$data$s2)

    BIC.bestfit = - 2 * llBestfit + (5 * log(sum(data.LN$data$n)))
    BIC.SM = - 2 * llSM + ((data.LN$data$N + 1) * log(sum(data.LN$data$n)))

    bf = exp(-0.5 * (BIC.bestfit - BIC.SM))

    if(bf < 1/10){
      warn.bf = 'None of the models provide an adequate fit do the data.'
    }else if(bf > 1/10 & bf < 10){
      warn.bf = 'Best fitting model fits well.'
    }else if(bf > 10){
      warn.bf = 'attention: bayes factor is larger than 10 in favor of the best fitting model'
    }

  }

  return(list(bayesFactor = bf,
              means.SM = means.SM,
              par.best = pars.bestfit,
              BIC.bestfit = BIC.bestfit,
              BIC.SM = BIC.SM,
              warn.bf = warn.bf)
  )

}

#'  function to perform MCMC sampling for all the dose-response model.
#'
#' @param best.fit best fitting model
#' @param data.N list containing data values for the normal models
#' @param data.LN list containing data values for the lognormal models
#' @param stanBest stan object for the best fitting model
#' @param type estimation type. Laplace approximation or McMC
#' @param ndraws number of draws to be made from the posterior distribution. Defaults to 30000
#' @param nrchains number of MCMC chains. Defaults to 3
#' @param nriterations number of MCMC iterations.Defaults to 3000
#' @param warmup  number of MCMC iterations for warmup. Defaults to 1000
#' @param delta adapt_delta value for the HMC in stan. See \code{\link[rstan]{sampling}} for more.
#'              Defaults to 0.8.
#' @param treedepth tree_depth value for the HMC in stan. See \code{\link[rstan]{sampling}} for more.
#'                  Defaults to 10.
#' @param seed random seed for reproducibility. Defaults to 123
#'
#' @examples{
#'
#'
#' }
#'
#'
#' @export modelTestC
#'
modelTestC <- function(best.fit, data.N, data.LN, stanBest, type, seed,
                       ndraws, nrchains, nriterations, warmup, delta, treedepth){

  if(grepl('_N', best.fit)){

    means.all <- data.N$data$data %>%
      group_by(dose) %>%
      summarise(mresp = mean(response))
    dose.a = unique(data.N$data$data$dose)
    mean.a = c()
    for(m in 1:length(dose.a)){
      mean.a[m] <- means.all$mresp[means.all$dose == dose.a[m]]
    }

    svSM = list(par = c(mean.a[1],
                        diff(mean.a),
                        log(1/var(data.N$data$y[data.N$data$y!=0])),
                        0.5))

    priorSM = list(
      priormu = c(mean.a[1],
                  diff(mean.a),
                  -2*log(1.5*sd(data.N$data$y[data.N$data$y!=0])),
                  0.5),
      priorSigma = diag(c(1, rep(1, length(dose.a)-1), 1)),
      priorlb = 0.001,
      priorub = c(2*mean.a[1],
                  max(abs(diff(mean.a)))*10)
    )

    data.modstanSM = list(N=data.N$data$N, n=data.N$data$n, nc=data.N$data$nc, maxN=data.N$data$maxN, maxNc=data.N$data$maxNc,
                          nij=data.N$data$nij, y=data.N$data$y, q=data.N$data$q, shift=data.N$data$shift,
                          priormu=priorSM$priormu, priorSigma=priorSM$priorSigma,
                          priorlb=priorSM$priorlb, priorub=priorSM$priorub,
                          priorg=4, data_type=data.N$data$data_type
    )

    if(type == 'MCMC'){

      sv=rstan::optimizing(stanmodels$mSMc,data = data.modstanSM,init=svSM)$par

      initf2 <- function(chain_id = 1) {
        list(par=sv[1:(data.N$data$N+2)] + rnorm(data.N$data$N+2, sd = 0.01*abs(sv[1:(data.N$data$N+2)])) ,alpha = chain_id)
      }
      init_ll <- lapply(1:nrchains, function(id) initf2(chain_id = id))
      fitstanSM = rstan::sampling(stanmodels$mSMc, data = data.modstanSM, init=init_ll, iter = nriterations,
                                  chains = nrchains, warmup = warmup, seed = seed,
                                  control = list(adapt_delta = delta, max_treedepth =treedepth),
                                  show_messages = F, refresh = 0)

      pars.bestfit = apply(as.matrix(stanBest),2,median)[c("par1","par2","par3","par4","par5","par6")]

      parsSM = as.matrix(fitstanSM)
      means.SM = apply(parsSM[, paste0('mu[', 1:data.N$data$N, ']')], 2, median)
      pars.SM = apply(parsSM[, c(paste0('a[', 1:data.N$data$N, ']'), paste0('par[', data.N$data$N+1, ']'),
                                 paste0('par[', data.N$data$N+2, ']'))], 2, median)

    }else if(type == 'Laplace'){

      all.pars.bestfit = par_extractC(stanBest, model_name = best.fit)
      pars.bestfit = apply(all.pars.bestfit[,c(paste0("p",1:4), "is2t", "rho")], 2, median)

      optSM = optimizing(stanmodels$mSMc, data = data.modstanSM,
                         seed=as.integer(seed), draws = ndraws,
                         init = svSM, hessian=TRUE)
      pars.SM = apply(optSM$theta_tilde[, c(paste0('a[', 1:data.N$data$N, ']'), paste0('par[', data.N$data$N+1, ']'),
                                            paste0('par[', data.N$data$N+2, ']'))], 2, median)
      means.SM = apply(optSM$theta_tilde[, paste0('mu[', 1:data.N$data$N, ']')], 2, median)

    }

    if(data.N$data$is_increasing == 1){
      llfun = paste0('llf',best.fit,'Ic')
    }else if(data.N$data$is_decreasing == 1){
      llfun = paste0('llf',best.fit,'Dc')
    }
    llBestfitf = get(llfun)
    llBestfit = llBestfitf(x = pars.bestfit,
                           d=data.N$data$x,
                           n=data.N$data$n,
                           nij=data.N$data$nij,
                           y=data.N$data$y,
                           qval=data.N$data$q)

    llSM = llfSM_Nc(pars.SM,
                    d=data.N$data$x,
                    n=data.N$data$n,
                    nij=data.N$data$nij,
                    y=data.N$data$y,
                    qval=data.N$data$q)

    BIC.bestfit = - 2 * llBestfit + (6 * log(sum(data.N$data$y!=0)))
    BIC.SM = - 2 * llSM + ((data.N$data$N + 2) * log(sum(data.N$data$y!=0)))

    bf = exp(-0.5 * (BIC.bestfit - BIC.SM))

    if(bf < 1/10){
      warn.bf = 'None of the models provide an adequate fit do the data.'
    }else if(bf > 1/10){
      warn.bf = 'Best fitting model fits sufficiently well.'
    }


    #############--LOGNORMAL--##################################

  }else if(grepl('_LN', best.fit)){

    means.all <- data.LN$data$data %>%
      group_by(dose) %>%
      summarise(mresp = mean(response))
    dose.a = unique(data.LN$data$data$dose)
    mean.a = c()
    for(m in 1:length(dose.a)){
      mean.a[m] <- means.all$mresp[means.all$dose == dose.a[m]]
    }

    svSM = list(par = c(mean.a[1],
                        diff(mean.a),
                        log(1/var(data.LN$data$y[data.LN$data$y!=0])),
                        0.5))

    priorSM = list(
      priormu = c(mean.a[1],
                  diff(mean.a),
                  -2*log(1.5*sd(data.LN$data$y[data.LN$data$y!=0])),
                  0.5),
      priorSigma = diag(c(1, rep(1, length(dose.a)-1), 1)),
      priorlb = 0.001,
      priorub = c(2*mean.a[1],
                  max(abs(diff(mean.a)))*10)
    )

    data.modstanSM = list(N=data.LN$data$N, n=data.LN$data$n, nc=data.LN$data$nc, maxN=data.LN$data$maxN, maxNc=data.LN$data$maxNc,
                          nij=data.LN$data$nij, y=data.LN$data$y, q=data.LN$data$q, shift=data.LN$data$shift,
                          priormu=priorSM$priormu, priorSigma=priorSM$priorSigma,
                          priorlb=priorSM$priorlb, priorub=priorSM$priorub,
                          priorg=4, data_type=data.LN$data$data_type
    )

    if(type == 'MCMC'){

      sv=rstan::optimizing(stanmodels$mSMc,data = data.modstanSM,init=svSM)$par

      initf2 <- function(chain_id = 1) {
        list(par=sv[1:(data.LN$data$N+2)] + rnorm(data.LN$data$N+2, sd = 0.01*abs(sv[1:(data.LN$data$N+2)])) ,alpha = chain_id)
      }
      init_ll <- lapply(1:nrchains, function(id) initf2(chain_id = id))
      fitstanSM = rstan::sampling(stanmodels$mSMc, data = data.modstanSM, init=init_ll, iter = nriterations,
                                  chains = nrchains, warmup = warmup, seed = seed,
                                  control = list(adapt_delta = delta, max_treedepth =treedepth),
                                  show_messages = F, refresh = 0)

      pars.bestfit = apply(as.matrix(stanBest),2,median)[c("par1","par2","par3","par4","par5","par6")]

      parsSM = as.matrix(fitstanSM)
      means.SM = apply(parsSM[, paste0('mu[', 1:data.LN$data$N, ']')], 2, median)
      pars.SM = apply(parsSM[, c(paste0('a[', 1:data.LN$data$N, ']'), paste0('par[', data.LN$data$N+1, ']'),
                                 paste0('par[', data.LN$data$N+2, ']'))], 2, median)


    }else if(type == 'Laplace'){

      all.pars.bestfit = par_extractC(stanBest, model_name = best.fit)
      pars.bestfit = apply(all.pars.bestfit[,c(paste0("p",1:4), "is2t", "rho")], 2, median)

      optSM = optimizing(stanmodels$mSMc, data = data.modstanSM,
                         seed=as.integer(seed), draws = ndraws,
                         init = svSM, hessian=TRUE)
      pars.SM = apply(optSM$theta_tilde[, c(paste0('a[', 1:data.LN$data$N, ']'), paste0('par[', data.LN$data$N+1, ']'),
                                            paste0('par[', data.LN$data$N+2, ']'))], 2, median)
      means.SM = apply(optSM$theta_tilde[, paste0('mu[', 1:data.LN$data$N, ']')], 2, median)

    }

    if(data.LN$data$is_increasing == 1){
      llfun = paste0('llf',best.fit,'Ic')
    }else if(data.LN$data$is_decreasing == 1){
      llfun = paste0('llf',best.fit,'Dc')
    }
    llBestfitf = get(llfun)
    llBestfit = llBestfitf(x = pars.bestfit,
                           d=data.LN$data$x,
                           n=data.LN$data$n,
                           nij=data.LN$data$nij,
                           y=data.LN$data$y,
                           qval=data.LN$data$q,
                           shift=data.LN$data$shift)

    llSM = llfSM_LNc(pars.SM,
                     d=data.LN$data$x,
                     n=data.LN$data$n,
                     nij=data.LN$data$nij,
                     y=data.LN$data$y,
                     qval=data.LN$data$q,
                     shift=data.LN$data$shift)

    BIC.bestfit = - 2 * llBestfit + (6 * log(sum(data.LN$data$y!=0)))
    BIC.SM = - 2 * llSM + ((data.LN$data$N + 2) * log(sum(data.LN$data$y!=0)))

    bf = exp(-0.5 * (BIC.bestfit - BIC.SM))

    if(bf < 1/10){
      warn.bf = 'None of the models provide an adequate fit do the data.'
    }else if(bf > 1/10){
      warn.bf = 'Best fitting model fits sufficiently well.'
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

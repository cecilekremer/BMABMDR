
#'  Function used internally to test whether the best fitting model fits the data well, compared to a saturated (non-monotone) model
#'
#' @param best.fit model name
#' @param data.N object as given by PREP_DATA_N
#' @param data.LN object as given by PREP_DATA_LN
#' @param data.Q object as given by PREP_DATA_Q
#' @param stanBest stan object for the model given in \code{best.fit}
#' @param type estimation type. Laplace or MCMC
#' @param seed random seed for reproducibility
#' @param ndraws number of draws to be made from the posterior distribution. Defaults to 30000
#' @param nrchains number of MCMC chains. Defaults to 3
#' @param nriterations number of MCMC iterations.Defaults to 3000
#' @param warmup  number of MCMC iterations for warmup. Defaults to 1000
#' @param delta adapt_delta value for the HMC in stan. See \code{\link[rstan]{sampling}} for more.
#'              Defaults to 0.8.
#' @param treedepth tree_depth value for the HMC in stan. See \code{\link[rstan]{sampling}} for more.
#'                  Defaults to 10.
#'
#' @description The function compares a dose-response model to the saturated model using an approximation to BIC.
#'
#' @return Bayes factor for the best fitting model compared to the saturated ANOVA model.
#'
#'
#' @export modelTest
#'
modelTest <- function(best.fit, data.N, data.LN, stanBest, type, seed,
                      ndraws, nrchains, nriterations, warmup, delta, treedepth){

  if(grepl('_N', best.fit)){

    svSM=list(par = c(data.N$data$priormu[1], # background
                      diff(data.N$data$m), # increments
                      log(1/mean(data.N$data$s2))) # invsigma2
    )

    priorSM = list(
      priormu = c(data.N$data$priormu[1],
                  diff(data.N$data$m), # increments
                  -2*log(1.5*mean(sqrt(data.N$data$s2)))),
      priorSigma = diag(c(1, rep(1, data.N$data$N-1), 1)),
      priorlb = data.N$data$priorlb[1],
      priorub = c(data.N$data$priorub[1],
                  max(abs(diff(data.N$data$m)))*10
      )
    )

    data.modstanSM=list(N=data.N$data$N, n=data.N$data$n, m=data.N$data$m, s2=data.N$data$s2, shift=data.N$data$shift,
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
      means.SM = apply(parsSM[, paste0('mu[', 1:data.N$data$N, ']')], 2, median, na.rm = T)
      pars.SM = apply(parsSM[, c(paste0('a[', 1:data.N$data$N, ']'), paste0('par[', data.N$data$N+1, ']'))], 2, median, na.rm = T)

    }else if(type == 'Laplace'){

      # all.pars.bestfit = par_extract(stanBest, model_name = best.fit)
      # pars.bestfit = apply(all.pars.bestfit[,c(paste0("p",1:4), "is2t")], 2, median)
      pars.bestfit = stanBest$par[c(1,2,9,4,5)]

      optSM = optimizing(stanmodels$mSM, data = data.modstanSM,
                         seed=as.integer(seed), draws = ndraws,
                         init = svSM, hessian=TRUE)
      pars.SM = apply(as.data.frame(optSM$theta_tilde)[, c(paste0('a[', 1:data.N$data$N, ']'), paste0('par[', data.N$data$N+1, ']'))], 2, median, na.rm = T)
      means.SM = apply(as.data.frame(optSM$theta_tilde)[, paste0('mu[', 1:data.N$data$N, ']')], 2, median, na.rm = T)

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
      warn.bf = paste0('None of the models provide an adequate fit do the data (Bayes factor is ', formatC(1/bf, digits=2, format='e'), ').')
      # warn.bf = paste0('None of the models provide an adequate fit do the data (Bayes factor is ', formatC(bf, digits=2, format='e'), ').')
    }else if(bf >= 1/10){
      warn.bf = paste0('Best fitting model fits sufficiently well (Bayes factor is ', formatC(1/bf, digits=2, format='e'), ').')
      # warn.bf = paste0('Best fitting model fits sufficiently well (Bayes factor is ', formatC(bf, digits=2, format='e'), ').')
    }


    #############--LOGNORMAL--##################################

  }else if(grepl('_LN', best.fit)){

    svSM=list(par = c(data.LN$data$priormu[1], # background
                      # exp(data.LN$data$m.org[1]),
                      diff(exp(data.LN$data$m.org)),
                      log(1/mean(data.LN$data$s2))) # invsigma2
    )

    priorSM = list(
      priormu = c(data.LN$data$priormu[1],
                  diff(exp(data.LN$data$m.org)),
                  -2*log(1.5*mean(sqrt(data.LN$data$s2)))),
      # priormu = c(exp(data.LN$data$m.org[1]),
      #             diff(exp(data.LN$data$m.org)),
      #             -2*log(1.5*mean(sqrt(data.LN$data$s2)))),
      priorSigma = diag(c(1, rep(1, data.LN$data$N-1), 1)),
      priorlb = data.LN$data$priorlb[1],
      priorub = c(data.LN$data$priorub[1],
                  max(abs(diff(exp(data.LN$data$m.org))))*3
      )
    )

    data.modstanSM=list(N=data.LN$data$N, n=data.LN$data$n, m=data.LN$data$m, s2=data.LN$data$s2, shift=data.LN$data$shift,
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
      fitstanSM = try(rstan::sampling(stanmodels$mSM, data = data.modstanSM, init=init_ll, iter = nriterations,
                                      chains = nrchains, warmup = warmup, seed = seed,
                                      control = list(adapt_delta = delta, max_treedepth =treedepth),
                                      show_messages = F, refresh = 0))
      # n.at <- 1
      if(class(fitstanSM) == 'try-error'){# && n.at < 11){
        # priorgam <- data.LN$data$shape.a + runif(1, -0.01, 0.01)
        # if(priorgam < 0) priorgam <- 0.0001
        priorgam <- 4
        data.modstanSM=list(N=data.LN$data$N, n=data.LN$data$n, m=data.LN$data$m, s2=data.LN$data$s2, shift=data.LN$data$shift,
                            priormu=priorSM$priormu, priorSigma=priorSM$priorSigma,
                            priorlb=priorSM$priorlb, priorub=priorSM$priorub,
                            data_type=data.LN$data$data_type, priorg = priorgam
        )
        sv=rstan::optimizing(stanmodels$mSM,data = data.modstanSM,init=svSM)$par
        initf2 <- function(chain_id = 1) {
          list(par=sv[1:(data.LN$data$N+1)] + rnorm(data.LN$data$N+1, sd = 0.01*abs(sv[1:(data.LN$data$N+1)])) ,alpha = chain_id)
        }
        init_ll <- lapply(1:nrchains, function(id) initf2(chain_id = id))
        fitstanSM = try(rstan::sampling(stanmodels$mSM, data = data.modstanSM, init=init_ll, iter = nriterations,
                                        chains = nrchains, warmup = warmup, seed = seed,
                                        control = list(adapt_delta = delta, max_treedepth =treedepth),
                                        show_messages = F, refresh = 0))
        # n.at <- n.at + 1
      }

      pars.bestfit = apply(as.matrix(stanBest),2,median)[c("par1","par2","par3","par4","par5")]

      parsSM = as.matrix(fitstanSM)
      means.SM = apply(parsSM[, paste0('mu[', 1:data.LN$data$N, ']')], 2, median, na.rm = T)
      pars.SM = apply(parsSM[, c(paste0('a[', 1:data.LN$data$N, ']'), paste0('par[', data.LN$data$N+1, ']'))], 2, median, na.rm = T)


    }else if(type == 'Laplace'){

      # all.pars.bestfit = par_extract(stanBest, model_name = best.fit)
      # pars.bestfit = apply(all.pars.bestfit[,c(paste0("p",1:4), "is2t")], 2, median, na.rm = T)
      pars.bestfit = stanBest$par[c(1,2,9,4,5)]

      optSM = optimizing(stanmodels$mSM, data = data.modstanSM,
                         seed=as.integer(seed), draws = ndraws,
                         init = svSM, hessian=TRUE)
      pars.SM = apply(as.data.frame(optSM$theta_tilde)[, c(paste0('a[', 1:data.LN$data$N, ']'), paste0('par[', data.LN$data$N+1, ']'))], 2, median, na.rm = T)
      means.SM = apply(as.data.frame(optSM$theta_tilde)[, paste0('mu[', 1:data.LN$data$N, ']')], 2, median, na.rm = T)

    }

    if(data.LN$data$is_increasing == 1){
      llfun = paste0('llf',best.fit,'I')
    }else if(data.LN$data$is_decreasing == 1){
      llfun = paste0('llf',best.fit,'D')
    }
    llBestfitf = get(llfun)
    llBestfit = llBestfitf(x = pars.bestfit, data.LN$data$n, data.LN$data$x, data.LN$data$m, data.LN$data$s2, data.LN$data$q, data.LN$data$shift)

    llSM = llfSM_LN(x = pars.SM, data.LN$data$n, data.LN$data$x, data.LN$data$m, data.LN$data$s2, data.LN$data$shift)

    BIC.bestfit = - 2 * llBestfit + (5 * log(sum(data.LN$data$n)))
    BIC.SM = - 2 * llSM + ((data.LN$data$N + 1) * log(sum(data.LN$data$n)))

    # print(means.SM);  print(pars.SM); print(pars.bestfit); print(llSM); print(llBestfit); print(best.fit)

    bf = exp(-0.5 * (BIC.bestfit - BIC.SM)) # bf in factor of SM if bf < 1/10

    if(bf < 1/10){
      warn.bf = paste0('None of the models provide an adequate fit do the data (Bayes factor is ', formatC(1/bf, digits=2, format='e'), ').')
      # warn.bf = paste0('None of the models provide an adequate fit do the data (Bayes factor is ', formatC(bf, digits=2, format='e'), ').')
    }else if(bf >= 1/10){
      warn.bf = paste0('Best fitting model fits sufficiently well (Bayes factor is ', formatC(1/bf, digits=2, format='e'), ').')
      # warn.bf = paste0('Best fitting model fits sufficiently well (Bayes factor is ', formatC(bf, digits=2, format='e'), ').')
    }

  }

  return(list(bayesFactor = 1/bf,
              # means.SM = means.SM,
              # par.best = pars.bestfit,
              # BIC.bestfit = BIC.bestfit,
              # BIC.SM = BIC.SM,
              pars.SM = pars.SM,
              # llSM = llSM,
              warn.bf = warn.bf)
  )

}

#' @rdname modelTest
#' @export
modelTestC <- function(best.fit, data.N, data.LN, stanBest, type, seed,
                       ndraws, nrchains, nriterations, warmup, delta, treedepth){

  if(grepl('_N', best.fit)){

    means.all <- data.N$data$data %>%
      dplyr::group_by(dose) %>%
      dplyr::summarise(mresp = mean(response))
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
      means.SM = apply(parsSM[, paste0('mu[', 1:data.N$data$N, ']')], 2, median, na.rm = T)
      pars.SM = apply(parsSM[, c(paste0('a[', 1:data.N$data$N, ']'), paste0('par[', data.N$data$N+1, ']'),
                                 paste0('par[', data.N$data$N+2, ']'))], 2, median, na.rm = T)

    }else if(type == 'Laplace'){

      # all.pars.bestfit = par_extractC(stanBest, model_name = best.fit)
      # pars.bestfit = apply(all.pars.bestfit[,c(paste0("p",1:4), "is2t", "rho")], 2, median)
      pars.bestfit = stanBest$par[c(1,2,10,4,5,6)]

      optSM = optimizing(stanmodels$mSMc, data = data.modstanSM,
                         seed=as.integer(seed), draws = ndraws,
                         init = svSM, hessian=TRUE)
      pars.SM = apply(as.data.frame(optSM$theta_tilde)[, c(paste0('a[', 1:data.N$data$N, ']'), paste0('par[', data.N$data$N+1, ']'),
                                                           paste0('par[', data.N$data$N+2, ']'))], 2, median, na.rm = T)
      means.SM = apply(as.data.frame(optSM$theta_tilde)[, paste0('mu[', 1:data.N$data$N, ']')], 2, median, na.rm = T)

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
      warn.bf = paste0('None of the models provide an adequate fit do the data (Bayes factor is ', formatC(1/bf, digits=2, format='e'), ').')
      # warn.bf = paste0('None of the models provide an adequate fit do the data (Bayes factor is ', formatC(bf, digits=2, format='e'), ').')
    }else if(bf >= 1/10){
      warn.bf = paste0('Best fitting model fits sufficiently well (Bayes factor is ', formatC(1/bf, digits=2, format='e'), ').')
      # warn.bf = paste0('Best fitting model fits sufficiently well (Bayes factor is ', formatC(bf, digits=2, format='e'), ').')
    }


    #############--LOGNORMAL--##################################

  }else if(grepl('_LN', best.fit)){

    means.all <- data.LN$data$data %>%
      dplyr::group_by(dose) %>%
      dplyr::summarise(mresp = mean(response))
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
      means.SM = apply(parsSM[, paste0('mu[', 1:data.LN$data$N, ']')], 2, median, na.rm = T)
      pars.SM = apply(parsSM[, c(paste0('a[', 1:data.LN$data$N, ']'), paste0('par[', data.LN$data$N+1, ']'),
                                 paste0('par[', data.LN$data$N+2, ']'))], 2, median, na.rm = T)


    }else if(type == 'Laplace'){

      # all.pars.bestfit = par_extractC(stanBest, model_name = best.fit)
      # pars.bestfit = apply(all.pars.bestfit[,c(paste0("p",1:4), "is2t", "rho")], 2, median, na.rm = T)
      pars.bestfit = stanBest$par[c(1,2,10,4,5,6)]

      optSM = optimizing(stanmodels$mSMc, data = data.modstanSM,
                         seed=as.integer(seed), draws = ndraws,
                         init = svSM, hessian=TRUE)
      pars.SM = apply(as.data.frame(optSM$theta_tilde)[, c(paste0('a[', 1:data.LN$data$N, ']'), paste0('par[', data.LN$data$N+1, ']'),
                                                           paste0('par[', data.LN$data$N+2, ']'))], 2, median, na.rm = T)
      means.SM = apply(as.data.frame(optSM$theta_tilde)[, paste0('mu[', 1:data.LN$data$N, ']')], 2, median, na.rm = T)

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
      warn.bf = paste0('None of the models provide an adequate fit do the data (Bayes factor is ', formatC(1/bf, digits=2, format='e'), ').')
      # warn.bf = paste0('None of the models provide an adequate fit do the data (Bayes factor is ', formatC(bf, digits=2, format='e'), ').')
    }else if(bf >= 1/10){
      warn.bf = paste0('Best fitting model fits sufficiently well (Bayes factor is ', formatC(1/bf, digits=2, format='e'), ').')
      # warn.bf = paste0('Best fitting model fits sufficiently well (Bayes factor is ', formatC(bf, digits=2, format='e'), ').')
    }

  }

  return(list(bayesFactor = 1/bf,
              # means.SM = means.SM,
              # par.best = pars.bestfit,
              # BIC.bestfit = BIC.bestfit,
              # BIC.SM = BIC.SM,
              llSM = llSM,
              warn.bf = warn.bf)
  )

}

#' @rdname modelTest
#' @export
modelTestQ <- function(best.fit, data.Q, stanBest, type, seed, ndraws, nrchains, nriterations,
                       warmup, delta, treedepth){

  N = data.Q$data$N


  if(data.Q$data$is_betabin == 1){

    data.modstanSM = list(N = data.Q$data$N, Y = data.Q$data$y, trials = data.Q$data$n, K = 2,
                           X = cbind(rep(1,data.Q$data$N), data.Q$data$x), Kc = 1,
                           N_1 = length(data.Q$data$x), M_1 = 1, J_1 = 1:length(data.Q$data$x), Z_1_1 = data.Q$data$x)

  } else if(data.Q$data$is_bin == 1){

    dose.a = data.Q$data$x
    y.a = data.Q$data$y
    n.a = data.Q$data$n

    if(length(data.Q$data$x) != length(unique(data.Q$data$x))){
      dose = sort(unique(data.Q$data$x))
      N = length(dose)
      y=rep(NA,N)
      n=rep(NA,N)
      for (iu in (1:N)){
        y[iu] = sum(data.Q$data$y[data.Q$data$x == dose[iu]])
        n[iu] = sum(data.Q$data$n[data.Q$data$x == dose[iu]])
      }
      y.a = y
      dose.a = dose
      n.a = n
    }

    N <- length(unique(data.Q$data$x))
    Ndose <- length(unique(dose.a))

    data.modstanSM = list(N = N, Y = y.a, trials = n.a, K = 2, X = cbind(rep(1, N), dose.a), Kc = 1)

  } else stop("data must be either clustered or independent")

  if(type == 'MCMC'){

    if(data.Q$data$is_bin == 1){

      fitstanSM <- rstan::sampling(stanmodels$mSM_Q, data = data.modstanSM, iter = nriterations,
                                   chains = nrchains, warmup = warmup, seed = seed,
                                   control = list(adapt_delta = delta, max_treedepth =treedepth),
                                   show_messages = F, refresh = 0)






    }else if(data.Q$data$is_betabin == 1){

      fitstanSM = rstan::sampling(stanmodels$mSM_Qc, data = data.modstanSM, iter = nriterations,
                                  chains = nrchains, warmup = warmup, seed = seed,
                                  control = list(adapt_delta = delta, max_treedepth =treedepth),
                                  show_messages = F, refresh = 0)

    }

  bridge_best <- bridgesampling::bridge_sampler(stanBest, silent=T)
  bridge_SM <- bridgesampling::bridge_sampler(fitstanSM, silent=T)
  # BF_brms_bridge = bridgesampling::bf(bridge_H0,bridge_SM)
  BF_brms_bridge = bridgesampling::bf(bridge_SM, bridge_best) # BF in favor of SM
  bf = BF_brms_bridge$bf

  } else if(type == 'Laplace'){

    if(data.Q$data$is_bin == 1){

      optSM <- rstan::optimizing(stanmodels$mSM_Q, data = data.modstanSM, hessian = T, draws = ndraws)
      llSM <- llfSM_Q(b = optSM$par[1], Intercept = optSM$par[4], # non-centered intercept
                      Y = data.modstanSM$Y, trials = data.modstanSM$trials, Xc = as.matrix(data.Q$data$x))
      llfun = paste0('llf',best.fit,'_Q')
      llBestfitf = get(llfun)
      llBestfit <- llBestfitf(x = stanBest$par[1:3], nvec = data.Q$data$n, dvec = data.Q$data$x,
                              yvec = data.Q$data$y, qval = data.Q$data$q)

      BIC.bestfit = - 2 * llBestfit + (3 * log(sum(data.Q$data$n))) # parms: a, b, d
      BIC.SM = - 2 * llSM + (2 * log(sum(data.Q$data$n)))


    }else if(data.Q$data$is_betabin == 1){

      optSM <- rstan::optimizing(stanmodels$mSM_Qc, data = data.modstanSM, hessian = T, draws = ndraws)
      llSM <- llfSM2_Q(b = optSM$par[1], Intercept = optSM$par[length(optSM$par)],
                        r_1_1 = optSM$par[(data.Q$data$N + 4):(2*data.Q$data$N + 3)],
                        Y = data.Q$data$y,
                        trials = data.Q$data$n, Xc = as.matrix(data.Q$data$x),
                        J_1 = 1:length(data.Q$data$x), Z_1_1 = data.Q$data$x, N = data.Q$data$N)

      llfun = paste0('llf',best.fit,'2_Q')
      llBestfitf = get(llfun)
      llBestfit <- llBestfitf(x = stanBest$par[1:3], nvec = data.Q$data$n, dvec = data.Q$data$x, yvec = data.Q$data$y,
                            qval = data.Q$data$q, rho = stanBest$par[stringr::str_detect(names(stanBest$par),'rho') &
                                                                           !stringr::str_detect(names(stanBest$par),'eta')])

      BIC.bestfit = - 2 * llBestfit + (4 * log(sum(data.Q$data$n)))
      BIC.SM = - 2 * llSM + ((2+data.Q$data$N) * log(sum(data.Q$data$n)))

    }

    bf = exp(0.5 * (BIC.bestfit - BIC.SM))

  }

  if(bf > 10){
    warn.bf = paste0('None of the models provide an adequate fit do the data (Bayes factor in favor of saturated model is ', formatC(bf, digits = 2, format = 'e'), ').')
  }else{
    warn.bf = paste0('Best fitting model fits sufficiently well (Bayes factor in favor of saturated model is ', formatC(bf, digits = 2, format = 'e'), ').')
  }

  return(list(bayesFactor = bf,
              # means.SM = means.SM,
              # par.best = pars.bestfit,
              # BIC.bestfit = BIC.bestfit,
              # BIC.SM = BIC.SM,
              llSM = llSM,
              warn.bf = warn.bf)
  )

}


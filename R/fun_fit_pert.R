
#' Function to get the mode
#'
#' @param v vector of observations
#'
#' @return numeric value
#'
#' @export
#'
getmode <- function(v){
  dx <- density(v)
  dx$x[dx$y == max(dx$y)]
}

#' Function for loglikelihood of PERT distribution
#'
#' @param m mode
#' @param s shape
#' @param lb minimum
#' @param ub maximum
#'
#' @return loglikelihood
#'
#' @export
#'
fun_loglik <- function(m, s, lb, ub){
  return(
    -sum(dpert(y, min = lb, mode = m, max = ub, shape = s,
               log = TRUE))
  )
}

#' Function to get 1000 quantiles
#'
#' @param x vector of observations
#'
#' @return vector of quantiles
#'
#' @export
fun.quantile.dist <- function(x){
  qnr <- 1000
  q.y <- quantile(x, seq(0, 1, 1/qnr))
  qdat_l <- q.y[1:qnr]
  qdat_r <- q.y[2:(qnr+1)]
  qdat_m <- (qdat_l + qdat_r)/2
  return(qdat_m)
}


#' Function to fit PERT distribution
#'
#' This function can be used to estimate the shape parameter of a PERT distribution by fitting to a posterior
#'
#' @param mod.obj an object of class BMADR
#' @param param which parameter's posterior the PERT should be fitted to (BMD, bkg, maxy)
#' @param method estimation method to use (MLE, JKL)
#'
#' @return Vector containing the parameters of a PERT distribution
#'
#' @export
#'
fit_pert <- function(mod.obj, param = 'BMD', method = 'MLE'){
  ## Estimate shape of PERT for each parameter

  require(bbmle)
  require(mc2d)
  require(kldest)

  if(class(mod.obj)[2] == "LP"){
    if(param == "BMD"){
      q_post <- fun.quantile.dist(mod.obj$BMDMixture)
    }else if(param == "bkg"){
      q_post <- fun.quantile.dist(mod.obj$bkg_post)
    }else if(param == "maxy"){
      q_post <- fun.quantile.dist(mod.obj$maxy_post)
    }else{
      stop(" param must be one of BMD, bkg, or maxy ")
    }
  }else if(class(mod.obj)[2] == "BS"){
    if(param == "BMD"){
      q_post <- fun.quantile.dist(mod.obj$MA_post_full_bs)
    }else if(param == "bkg"){
      q_post <- fun.quantile.dist(mod.obj$bkg_post_bs)
    }else if(param == "maxy"){
      q_post <- fun.quantile.dist(mod.obj$maxy_post_bs)
    }else{
      stop(" param must be one of BMD, bkg, or maxy ")
    }
  }else{
    stop("mod.obj must be an object of class BMADR")
  }

  y <- q_post

  if(method == "MLE"){
    # estimate shape
    ml_fit <- try(bbmle::mle2(minuslogl = fun_loglik,
                              data = list(y = y),
                              fixed = list(lb = min(y)-0.000000001, ub = max(y)+0.000000001, m = getmode(y)),
                              lower = list(s = 0.000000001),
                              upper = list(s = 20),
                              start = list(s = 0.1), method = 'L-BFGS-B'), silent = T)
    if(class(ml_fit)[1] == 'try-error'){
      # BMD.pars <- c(dataN$data$priorlb[2], dataN$data$priormu[2], dataN$data$priorub[2], 0.0001)
      stop("Could not fit a PERT to this posterior")
    }else{
      pars <- c(coef(ml_fit)[3], # min
                coef(ml_fit)[1], # mode
                coef(ml_fit)[4], # max
                coef(ml_fit)[2] # shape
      )
    }
    pert.l <- pars[1]
    pert.u <- pars[3]
    pert.m <- pars[2]
    pert.s <- pars[4]
    pert.mu <- (pert.l + pert.s*pert.m + pert.u) / (pert.s + 2)
    pert.var <- ((pert.mu - pert.l)*(pert.u - pert.mu)) / (pert.s + 3)
    pert.var.discounted <- pert.var * 2
    pert.s.discounted <- (((pert.mu - pert.l)*(pert.u - pert.mu))/pert.var.discounted) - 3
    if(pert.s.discounted < 0.0001) { pert.s.discounted = 0.0001 }

  }else if(method == "JKL"){

    true.dist <- as.numeric(q_post)  ## the true distribution
    mintd <- min(true.dist)-0.000000001
    maxtd <- max(true.dist)+0.000000001
    modetd <- getmode(true.dist)
    fun_KLdistance <- function(s){
      ensemble_jeffreys_kl <- mean(replicate(10, kld_est_nn(rpert(1000, min=mintd, mode=modetd, max=maxtd, shape=s), true.dist, k=5)+kld_est_nn(true.dist, rpert(1000, min=mintd, mode=modetd, max=maxtd, shape=s), k=5)))
      return(ensemble_jeffreys_kl)  # X true distribution, # Y approximate distribution
    }
    KL_fit <- try(optim(par=c(1),fun_KLdistance,method="Brent",lower = 0.000000001,
                        upper = 50), silent = T)
    shape.approx <- KL_fit$par
    pars <- c(mintd,modetd,maxtd,shape.approx)
    pert.l <- pars[1]
    pert.u <- pars[3]
    pert.m <- pars[2]
    pert.s <- pars[4]
    pert.mu <- (pert.l + pert.s*pert.m + pert.u) / (pert.s + 2)
    pert.var <- ((pert.mu - pert.l)*(pert.u - pert.mu)) / (pert.s + 3)
    pert.var.discounted <- pert.var * 2
    pert.s.discounted <- (((pert.mu - pert.l)*(pert.u - pert.mu))/pert.var.discounted) - 3
    if(pert.s.discounted < 0.0001) { pert.s.discounted = 0.0001 }

  }else{
    stop("method must be one of MLE, JKL")
  }

  ## Save output
  pars.est <- as.data.frame(t(c(pert.l, pert.u, pert.m, pert.s, pert.mu, pert.var, pert.var.discounted, pert.s.discounted)))
  names(pars.est) <- c('min','max','mode','shape','mean','variance','variance.discounted','shape.discounted')

  return(pars.est)

}

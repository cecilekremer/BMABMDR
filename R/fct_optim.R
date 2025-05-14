
#' Function to perform Laplace approximation for a given dose-response model (used internally).
#'
#' @param mod stan model
#' @param data list containing the data to be passed to the stan model
#' @param stv list of starting values
#' @param ndraws number of draws to be made from the posterior distribution
#' @param seed random seed for reproducibility
#' @param pvec probability vector to compute credible interval for the BMD
#'
#' @examples
#' \dontrun{
#' fun_optim(stanmodels$mE4, data, start, 30000, 123, c(0.05,0.5,0.95))
#' }
#'
#' @return a stan results object, which is a list containing point estimates, value of the log-posterior, value of the return code from the optimizer, Hessian matrix, matrix of parameter draws
#'
#' @export fun_optim
#'
fun_optim = function(mod, data, stv,
                     ndraws,seed,pvec){

  opt=try(rstan::optimizing(mod,data = data,
                            seed=as.integer(seed),draws=ndraws,
                            init=stv, hessian=TRUE),silent=T)
  if(length(opt)!=9) {

    # svh = try(list(par=as.vector(rstan::optimizing(mod,data = data)$par[1:5])),silent=T)
    opt=try(rstan::optimizing(mod,data = data,seed=as.integer(seed),draws=ndraws,hessian=TRUE),silent=T)
  }
  if(class(opt) == 'try-error'){
    opt = c(NA, NA, NA)
  }
  n.attempts <- 1
  while((ifelse(is.na(opt[3]),TRUE,(opt[3]!=0)) | length(opt)!=9) & n.attempts < 100){
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
                              seed=as.integer(seed),draws=ndraws,
                              init=svh,hessian=TRUE),silent=T)

    n.attempts <- n.attempts + 1
  }
  if(class(opt) == 'try-error'){
    opt = c(NA, NA, NA)
  }
  # else{
  #   message('All errors related to \'chol.default()\' can be ignored')
  # }

  attr(opt, "class") <- "stanfitOptim"
  return(opt)

}
#' @rdname fun_optim
#' @export
fun_optimC = function(mod, data, stv,
                      ndraws,seed,pvec){

  opt=try(rstan::optimizing(mod,data = data,
                            seed=as.integer(seed),draws=ndraws,
                            init=stv, hessian=TRUE),silent=T)
  if(length(opt)!=9) {

    # svh = try(list(par=as.vector(rstan::optimizing(mod,data = data)$par[1:5])),silent=T)
    opt=try(rstan::optimizing(mod,data = data,seed=as.integer(seed),draws=ndraws,hessian=TRUE),silent=T)
  }
  if(class(opt) == 'try-error'){
    opt = c(NA, NA, NA)
  }
  n.attempts <- 1
  while((ifelse(is.na(opt[3]),TRUE,(opt[3]!=0)) | length(opt)!=9) & n.attempts < 100){
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
    par[4:6] = par[4:6] + rnorm(3, sd = 0.01*abs(par[4:5]))

    ## BMD between 0 and 1
    # if(par[2] > 1) par[2] = 0.9

    if(data$data_type==1|data$data_type==2){
      pars3d = numeric()
      pars3i = par[3]
      dim(pars3d)=0
      dim(pars3i)=1
      svh = list(par1=par[1],par2=par[2],pars3i=pars3i,pars3d=pars3d,par4=par[4],par5=par[5],par6=par[6])
    }else if(data$data_type==3|data$data_type==4){
      pars3d = par[3]
      pars3i = numeric()
      dim(pars3d)=1
      dim(pars3i)=0
      svh = list(par1=par[1],par2=par[2],pars3i=pars3i,pars3d=pars3d,par4=par[4],par5=par[5],par6=par[6])
    }
    opt=try(rstan::optimizing(mod,data = data,
                              seed=as.integer(seed),draws=ndraws,
                              init=svh,hessian=TRUE),silent=T)
    n.attempts <- n.attempts + 1
  }
  if(class(opt) == 'try-error'){
    opt = c(NA, NA, NA)
  }

  attr(opt, "class") <- "stanfitOptim"
  return(opt)

}
#' @rdname fun_optim
#' @export
fun_optimQ = function(mod, data, stv,
                      ndraws, seed, pvec){

  opt=try(rstan::optimizing(mod,data = data,
                            seed=as.integer(seed),draws=ndraws,
                            init=stv, hessian=TRUE),silent=T)
  # if(is.na(opt$return_code) | opt$return_code != 0) {
  if(length(opt)!=9){
    # svh = try(list(par=as.vector(rstan::optimizing(mod,data = data)$par[1:5])),silent=T)
    opt=try(rstan::optimizing(mod,data = data,seed=as.integer(seed),draws=ndraws,hessian=TRUE),silent=T)
  }
  if(class(opt) == 'try-error'){
    opt = c(NA, NA, NA)
  }
  n.attempts <- 1
  if(class(opt) == 'try-error') {

    while(TRUE){ # set limit at 100 because takes longer if there are difficulties
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

      svh = list(par1=par[1], par2=par[2], par3=par[3])

      #opt=try(rstan::optimizing(mod,data = data,
      #                          seed=as.integer(seed),draws=ndraws,
      #                          init=svh,hessian=TRUE),silent=T)
      # seed <- runif(1, 1, 10000)
      opt=try(rstan::optimizing(mod,data = data, init = svh,
                                seed=as.integer(seed), draws=ndraws, hessian=TRUE),silent=T)
      n.attempts <- n.attempts + 1

      if(class(opt) == 'try-error') {
        next
      } else if((!is.na(opt[[3]]) | opt[[3]]==0 & length(opt)==9) | n.attempts == 100){
        break
      } else next
      # print(n.attempts)
    }

  } else {

    while((ifelse(is.na(opt[[3]]),TRUE,(opt[[3]]!=0)) | length(opt)!=9) & n.attempts < 100){ # set limit at 100 because takes longer if there are difficulties
      svh = stv
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

      svh = list(par1=par[1], par2=par[2], par3=par[3])

      #opt=try(rstan::optimizing(mod,data = data,
      #                          seed=as.integer(seed),draws=ndraws,
      #                          init=svh,hessian=TRUE),silent=T)
      # seed <- runif(1, 1, 10000)
      opt=try(rstan::optimizing(mod,data = data, init = svh,
                                seed=as.integer(seed), draws=ndraws, hessian=TRUE),silent=T)
      n.attempts <- n.attempts + 1
      # print(n.attempts)
    }

  }
  if(class(opt) == 'try-error'){
    opt = c(NA, NA, NA)
  }

  attr(opt, "class") <- "stanfitOptim"
  return(opt)

}
#' @rdname fun_optim
#' @export
fun_optimCov = function(mod, data, stv,
                        ndraws,seed,pvec){

  opt=try(rstan::optimizing(mod,data = data,
                            seed=as.integer(seed),draws=ndraws,
                            init=stv, hessian=TRUE),silent=T)
  if(length(opt)!=9) {

    # svh = try(list(par=as.vector(rstan::optimizing(mod,data = data)$par[1:5])),silent=T)
    opt=try(rstan::optimizing(mod,data = data,seed=as.integer(seed),draws=ndraws,hessian=TRUE),silent=T)
  }

  if(class(opt)!='try-error'){
    if(opt$return_code != 0){
      par = opt$theta_tilde
      a = par[1:data$nlevels_a]
      bmd = par[(data$nlevels_a+1):(data$nlevels_a+data$nlevels_BMD)]
      c = par[(data$nlevels_a+data$nlevels_BMD+1+data$nlevels_d+data$nlevels_a+data$nlevels+data$nlevels_a+1):
                (data$nlevels_a+data$nlevels_BMD+1+data$nlevels_d+data$nlevels_a+data$nlevels+data$nlevels_a+1)]
      d = par[(data$nlevels_a+data$nlevels_BMD+2):(data$nlevels_a+data$nlevels_BMD+1+data$nlevels_d)]
      s = par[(data$nlevels_a+data$nlevels_BMD+1+data$nlevels_d+1):(data$nlevels_a+data$nlevels_BMD+1+data$nlevels_d+data$nlevels_a)]

      if(data$data_type==1|data$data_type==2){
        pars3d = numeric()
        pars3i = c
        dim(pars3d)=0
        dim(pars3i)=1
        svh = list(par1=a,par2=bmd,pars3i=pars3i,pars3d=pars3d,par4=d,par5=s)
      }else if(data$data_type==3|data$data_type==4){
        pars3d = par[3]
        pars3i = numeric()
        dim(pars3d)=1
        dim(pars3i)=0
        svh = list(par1=a,par2=bmd,pars3i=pars3i,pars3d=pars3d,par4=d,par5=s)
      }

      dim(svh$par2) <- data$nlevels_BMD
      dim(svh$par1) <- data$nlevels_a
      dim(svh$par4) <- data$nlevels_d
      dim(svh$par5) <- data$nlevels_sigma

      opt=try(rstan::optimizing(mod,data = data,
                                seed=as.integer(seed),draws=ndraws,
                                init=svh,hessian=TRUE),silent=T)
    }
  }
  if(class(opt) == 'try-error'){
    opt = c(NA, NA, NA)
  }
  n.attempts <- 1
  while((ifelse(is.na(opt[3]),TRUE,(opt[3]!=0)) | length(opt)!=9) & n.attempts < 10){
    svh=stv
    par = unname(unlist(svh))

    a = par[1:data$nlevels_a]
    bmd = par[(data$nlevels_a+1):(data$nlevels_a+data$nlevels_BMD)]
    c = par[(data$nlevels_a+data$nlevels_BMD+1):(data$nlevels_a+data$nlevels_BMD+1)]
    d = par[(data$nlevels_a+data$nlevels_BMD+2):(data$nlevels_a+data$nlevels_BMD+1+data$nlevels_d)]
    s = par[(data$nlevels_a+data$nlevels_BMD+1+data$nlevels_d+1):(length(par))]



    if(data$is_informative_a == 1){
      a = runif(length(a), data$priorlb[1], data$priorub[1])
    }else{
      a = a + rnorm(length(a), sd = 0.01*abs(a))
    }

    if(data$is_informative_c == 1){
      if(data$is_increasing == 1){
        c = runif(1, data$priorlb[3] - data$L, data$priorub[3] - data$L)
      }else{
        c = runif(1, data$priorlb[3] / data$U, data$priorub[3] / data$U)
      }
    }else{
      c = c + rnorm(1, sd = 0.01*abs(c))
    }

    ## for decreasing, par[3] < 1
    if(data$is_decreasing == 1 & c > 1){
      c = 0.9999
    }

    if(data$is_informative_BMD == 1){
      bmd = runif(length(bmd), data$priorlb[2], data$priorub[2])
    }else{
      bmd = bmd + rnorm(length(bmd), sd = 0.01*abs(bmd))
    }

    d = d + rnorm(length(d), sd = 0.01*abs(d))
    s = s + rnorm(length(s), sd = 0.01*abs(s))

    if(data$data_type==1|data$data_type==2){
      pars3d = numeric()
      pars3i = c
      dim(pars3d)=0
      dim(pars3i)=1
      svh = list(par1=a,par2=bmd,pars3i=pars3i,pars3d=pars3d,par4=d,par5=s)
    }else if(data$data_type==3|data$data_type==4){
      pars3d = par[3]
      pars3i = numeric()
      dim(pars3d)=1
      dim(pars3i)=0
      svh = list(par1=a,par2=bmd,pars3i=pars3i,pars3d=pars3d,par4=d,par5=s)
    }

    dim(svh$par2) <- data$nlevels_BMD
    dim(svh$par1) <- data$nlevels_a
    dim(svh$par4) <- data$nlevels_d
    dim(svh$par5) <- data$nlevels_sigma

    opt=try(rstan::optimizing(mod,data = data,
                              seed=as.integer(seed),draws=ndraws,
                              init=svh,hessian=TRUE),silent=T)

    n.attempts <- n.attempts + 1
  }
  if(class(opt) == 'try-error'){
    opt = c(NA, NA, NA)
  }

  attr(opt, "class") <- "stanfitOptim"
  return(opt)

}
#' @rdname fun_optim
#' @export
fun_optimQCov = function(mod, data, stv,
                         ndraws,seed,pvec){

  opt=try(rstan::optimizing(mod,data = data,
                            seed=as.integer(seed),draws=ndraws,
                            init=stv, hessian=TRUE),silent=T)
  if(length(opt)!=9) {
    opt=try(rstan::optimizing(mod,data = data,seed=as.integer(seed),draws=ndraws,hessian=TRUE),silent=T)
  }

  if(class(opt)!='try-error'){
    if(opt$return_code != 0){

      par = opt$par

      a = par[1:data$nlevels_a]
      bmd = par[(data$nlevels_a+1):(data$nlevels_a+data$nlevels_BMD)]
      d = par[(data$nlevels_a+data$nlevels_BMD+1):(data$nlevels_a+data$nlevels_BMD+data$nlevels_d)]

      svh = list(par1=a,par2=bmd,par3=d)

      dim(svh$par2) <- data$nlevels_BMD
      dim(svh$par1) <- data$nlevels_a
      dim(svh$par3) <- data$nlevels_d

      opt=try(rstan::optimizing(mod,data = data,
                                seed=as.integer(seed),draws=ndraws,
                                init=svh,hessian=TRUE),silent=T)
    }
  }
  if(class(opt) == 'try-error'){
    opt = c(NA, NA, NA)
  }
  n.attempts <- 1
  while((ifelse(is.na(opt[3]),TRUE,(opt[3]!=0)) | length(opt)!=9) & n.attempts < 10){

    svh=stv
    par = unname(unlist(svh))

    a = par[1:data$nlevels_a]
    bmd = par[(data$nlevels_a+1):(data$nlevels_a+data$nlevels_BMD)]
    d = par[(data$nlevels_a+data$nlevels_BMD+1):(data$nlevels_a+data$nlevels_BMD+data$nlevels_d)]

    if(data$is_informative_a == 1){
      a = runif(length(a), data$priorlb[1], data$priorub[1])
    }else{
      a = a + rnorm(length(a), sd = 0.01*abs(a))
    }

    if(data$is_informative_BMD == 1){
      bmd = runif(length(bmd), data$priorlb[2], data$priorub[2])
    }else{
      bmd = bmd + rnorm(length(bmd), sd = 0.01*abs(bmd))
    }

    d = d + rnorm(length(d), sd = 0.01*abs(d))

    svh = list(par1=a,par2=bmd,par3=d)


    dim(svh$par2) <- data$nlevels_BMD
    dim(svh$par1) <- data$nlevels_a
    dim(svh$par3) <- data$nlevels_d

    opt=try(rstan::optimizing(mod,data = data,
                              seed=as.integer(seed),draws=ndraws,
                              init=svh,hessian=TRUE),silent=T)

    n.attempts <- n.attempts + 1
  }
  if(class(opt) == 'try-error'){
    opt = c(NA, NA, NA)
  }

  attr(opt, "class") <- "stanfitOptim"
  return(opt)

}

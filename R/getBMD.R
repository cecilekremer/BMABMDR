
#' Function for internal use
#'
#' @param pars value
#' @param pvec value
#' @param covar value
#' @param model_name value
#'
#' @return .
#' 
#' @export getBMD
#'
getBMD <- function(pars, pvec, covar, model_name){
  
  covar = match.arg(covar, c('a_sigma2', 'BMD_d', 'all', 'none'))
  nmpar <- colnames(pars)
  
  if(covar != 'none'){
    bmd <- pars[,grep("par2\\[", nmpar)]
  } else {
    bmd <- pars[,nmpar == "par2"]
  }
  
  if(covar == 'BMD_d' | covar == 'all'){
    bmd.mat <- t(apply(bmd, 2, quantile, pvec))
    colnames(bmd.mat) <- c('BMDL','BMD','BMDU')
    rownames(bmd.mat) <- rep(model_name, dim(bmd.mat)[1])    
  }else{
    bmd.mat <- quantile(bmd, pvec)
    names(bmd.mat) <- c('BMDL','BMD','BMDU')
    # rownames(bmd.mat) <- model_name    
  }

  
  return(bmd.mat)
}

#' @rdname getBMD
#' @export
getBMD_Q <- function(pars, pvec, covar, model_name){
  
  covar = match.arg(covar, c('background', 'BMD_d', 'all', 'none'))
  nmpar <- colnames(pars)
  
  if(covar != 'none'){
    bmd <- pars[,grep("par2\\[", nmpar)]
  } else {
    bmd <- pars[,nmpar == "par2"]
  }
  
  if(covar == 'BMD_d' | covar == 'all'){
    bmd.mat <- t(apply(bmd, 2, quantile, pvec))
    colnames(bmd.mat) <- c('BMDL','BMD','BMDU')
    rownames(bmd.mat) <- rep(model_name, dim(bmd.mat)[1])    
  }else{
    bmd.mat <- quantile(bmd, pvec)
    names(bmd.mat) <- c('BMDL','BMD','BMDU')
    # rownames(bmd.mat) <- model_name    
  }
  return(bmd.mat)
}
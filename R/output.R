
#' output function for estimated model parameter from each model
#'
#' @param par_obj list containing parameter draws from either Laplace or MCMC fit
#' @param pvec probability vector to compute credible interval for the BMD.
#' @param max.dose maximum dose tested
#'
#' @return dataframe containing model estimates
#'
#'
#' @export outLP
#'
outLP <- function(par_obj, pvec, max.dose){
  ret <- t(data.frame(
    # model parameters
    # transformed
    par.dt = quantile(par_obj$p4, c(0.025,0.5,0.975), na.rm=T),
    par.is2t = quantile(par_obj$is2t, c(0.025,0.5,0.975), na.rm=T),
    par.k = quantile(par_obj$k, pvec, na.rm=T),
    # original scale
    par.a = quantile(par_obj$a, c(0.025,0.5,0.975), na.rm=T),
    par.b = quantile(par_obj$b, c(0.025,0.5,0.975), na.rm=T),
    par.c = quantile(par_obj$c, c(0.025,0.5,0.975), na.rm=T),
    par.d = quantile(par_obj$d, c(0.025,0.5,0.975), na.rm=T),
    par.s2 = quantile(par_obj$s2, c(0.025,0.5,0.975), na.rm=T),
    # natural parameters
    BMD = quantile(par_obj$BMD*max.dose, pvec, na.rm=T),
    min.resp = quantile(par_obj$min_response, c(0.025,0.5,0.975), na.rm=T),
    max.resp = quantile(par_obj$max_response, c(0.025,0.5,0.975), na.rm=T),
    fold.change = quantile(par_obj$p3, c(0.025, 0.5, 0.975), na.rm=T)
  ))

  return(ret)
}


#### function to extract the sampled values for the model parameters and the BMD
#' function to extract the sampled values for the model parameters and the BMD
#'
#' @param mod_obj model object either of class stanfit or stanfitOptim
#' @param pars parameters to be extracted.
#' @param model_name name of the model to be extracted
#'
#' @return dataframe containing model estimates
#'
#' @export par_extract
#'
par_extract <- function(mod_obj, pars = c(letters[1:4],"k",
                                          paste0("par",1:5,""),
                                          "invsigma2","mu_inf","mu_0"), model_name) {

  if(is.stanfit(mod_obj)){

    pars_samples <- rstan::extract(mod_obj, pars)
    pars_samplesd <- data.frame(ModelName = model_name,
                                a = pars_samples$a,
                                b = pars_samples$b,
                                c = pars_samples$c,
                                d = pars_samples$d,
                                k = pars_samples$k,
                                BMD = pars_samples$`par2`,
                                p1 = pars_samples$`par1`, # mu(0)
                                p2 = pars_samples$`par2`, # BMD
                                p3 = pars_samples$`par3`, # fold change
                                p4 = pars_samples$`par4`, # log(d)
                                is2t = pars_samples$`par5`,
                                s2 = 1/pars_samples$invsigma2,
                                min_response = pars_samples$`mu_0`,
                                max_response = pars_samples$`mu_inf`
    )

  } else if(is.stanfitOptim(mod_obj)) {

    pars_samples <- mod_obj$theta_tilde[,pars]
    pars_samplesd <- data.frame(ModelName = model_name,
                                a = pars_samples[,"a"],
                                b = pars_samples[,"b"],
                                c = pars_samples[,"c"],
                                d = pars_samples[,"d"],
                                k = pars_samples[,"k"],
                                BMD = pars_samples[,"par2"],
                                p1 = pars_samples[,"par1"],
                                p2 = pars_samples[,"par2"],
                                p3 = pars_samples[,"par3"],
                                p4 = pars_samples[,"par4"],
                                is2t = pars_samples[,"par5"],
                                s2 = 1/pars_samples[,"invsigma2"],
                                min_response = pars_samples[,"mu_0"],
                                max_response = pars_samples[,"mu_inf"]
    )

    return(pars_samplesd)

  } else stop("please check if input is of class stanfit or stanfitOptim")



}

# function to turn the vector of weights into dataframe
#' function to extract the sampled values for the model parameters and the BMD
#'
#' @param mod_obj model object either of class stanfit or stanfitOptim
#' @param type type of weight to be extracted. It can be either of BS (bridge sampling), LP (Laplace) and
#'             both to extract bridge sampling and laplace approximation.
#' @param model_name name of the model to be extracted
#'
#' @return dataframe containing model estimates
#'
#' @export weights_extract
#'
weights_extract <- function(mod.obj, type = c("BS", "LP", "both")) {
  type <- match.arg(type)
  if(is.BMADR2(mod.obj)[3] == 2 & type == 'BS'){

    return(data.frame(Model = names(mod.obj$weights_bridge_sampling),
                      BS_Weights = mod.obj$weights_bridge_sampling))

  } else if(is.BMADR2(mod.obj)[3] == 2 & type == 'LP') {

    return(data.frame(Model = names(mod.obj$weights_bridge_sampling),
                      LP_Weights = mod.obj$weights_laplace))

  } else if(is.BMADR2(mod.obj)[3] == 2 & type == "both"){

    return(data.frame(Model = names(mod.obj$weights_bridge_sampling),
                      BS_Weights = mod.obj$weights_bridge_sampling,
                      LP_Weights = mod.obj$weights_laplace))

  } else if(is.BMADR2(mod.obj)[2] == 2){

    return(data.frame(Model = names(mod.obj$weights),
                      LP_Weights = mod.obj$weights))

  }else stop('please check mod.obj or provide type')

}

# function to extract mixture
#' function to extract BMD mixture
#'
#' @param mod_obj BMDBMA model object
#' @param conv logical to indicate if mixture should be only based on converged
#'             models or all the models regardless of the convergence status
#' @return dataframe containing Mixture per model
#' @export BMDmixture_extract
#'
BMDmixture_extract <- function(mod.obj, conv = FALSE){

  if(is.BMADR2(mod.obj)[3] == 2 & conv == TRUE) {
    BMDMixture <- data.frame(Model = "Model Averaged",
                             BMDMixture = mod.obj$BMDMixture.conv
    )
  } else if(is.BMADR2(mod.obj)[3] == 2 & conv == FALSE) {
    BMDMixture <- data.frame(Model = "Model Averaged",
                             BMDMixture = mod.obj$BMDMixture
    )
  } else if(is.BMADR2(mod.obj)[2] == 2) {

    BMDMixture <- data.frame(Model = "Model Averaged",
                             BMDMixture = mod.obj$BMDMixture
    )

  } else stop("please check mod.obj or supply conv")

  return(BMDMixture)

}


#' function to extract model averaged BMD
#'
#' @param mod_obj BMDBMA model object
#' @param conv logical to indicate if mixture should be only based on converged
#'             models or all the models regardless of the convergence status
#' @return dataframe containing BMDL, BMD and BMDU
#' @export BMDMA_extract
#'
BMDMA_extract <- function(mod.obj, conv = FALSE) {

  if(is.BMADR2(mod.obj)[3] == 2 & conv == TRUE) {
    BMDMA <- data.frame(Model = "Model Averaged",
                        Type = c("BS", "LP"),
                        BMDL = c(mod.obj$MA_bs_conv[1], mod.obj$MA_ls_conv[1]),
                        BMD = c(mod.obj$MA_bs_conv[2], mod.obj$MA_ls_conv[2]),
                        BMDU = c(mod.obj$MA_bs_conv[3], mod.obj$MA_ls_conv[3])
    )
  } else if(is.BMADR2(mod.obj)[3] == 2 & conv == FALSE) {
    BMDMA <- data.frame(Model = "Model Averaged",
                        Type = c("BS", "LP"),
                        BMDL = c(mod.obj$MA_bridge_sampling[1], mod.obj$MA_laplace[1]),
                        BMD = c(mod.obj$MA_bridge_sampling[2], mod.obj$MA_laplace[2]),
                        BMDU = c(mod.obj$MA_bridge_sampling[3], mod.obj$MA_laplace[3])
    )
  } else if(is.BMADR2(mod.obj)[2] == 2) {

    BMDMA <- data.frame(Model = "Model Averaged",
                        Type = "LP",
                        BMDL = mod.obj$MA[1],
                        BMD = mod.obj$MA[2],
                        BMDU = mod.obj$MA[3]
    )

  } else stop("check mod.obj or supply conv")

  return(BMDMA)

}


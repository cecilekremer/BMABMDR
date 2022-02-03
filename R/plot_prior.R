#' Function to plot prior vs posterior distribution
#'
#' @param mod.obj the model object returned by full laplace or sampling method
#' @param data the data for a specific model used by the above method
#' @param model_name the model for which plots should be generated
#' @param parms if TRUE, model parameters (background, fold change, BMD, d) are shown; if FALSE, background, maximum response, BMD and d are shown
#'
#' @return .
#'
#' @export
#'
plot_prior <- function(mod.obj, data, model_name, parms = T){ # pars = T for parameters, F for background & fold change

  if(model_name %in% mod.obj$models_included){

    names(mod.obj$parsN) = c("E4_N","IE4_N","H4_N","LN4_N","G4_N","QE4_N","P4_N","L4_N")
    names(mod.obj$parsLN) = c("E4_LN","IE4_LN","H4_LN","LN4_LN","G4_LN","QE4_LN","P4_LN","L4_LN")

    pars = c(mod.obj$parsN, mod.obj$parsLN)

    df.mod <- (as.data.frame(pars[[model_name]]))[, c("ModelName","a","p3","d","BMD","min_response","max_response")]
    q <- mod.obj$q

    obs.bkg <- mod.obj$data$m[1]

    obs.fold <- mod.obj$data$m[length(mod.obj$data$m)] / mod.obj$data$m[1]

    obs.max <- mod.obj$data$m[length(mod.obj$data$m)]

    bkg.prior <- rpert(dim(df.mod)[1], min = (data$priorlb[1]),
                       mode = (data$priormu[1]), max = (data$priorub[1]),
                       shape = data$shape.a)

    BMD.prior <- rpert(dim(df.mod)[1], min = (data$priorlb[2]),
                       mode = (data$priormu[2]), max = (data$priorub[2]),
                       shape = data$shape.BMD)

    fold.prior <- rpert(dim(df.mod)[1], min = (data$priorlb[3]),
                        mode = (data$priormu[3]), max = (data$priorub[3]),
                        shape = data$shape.c)

    max.prior <- fold.prior*bkg.prior


    bkg.post <- df.mod$min_response[df.mod$ModelName == model_name]

    BMD.post <- df.mod$BMD[df.mod$ModelName == model_name]

    max.post <- df.mod$max_response[df.mod$ModelName == model_name]

    fold.post <- df.mod$p3[df.mod$ModelName == model_name]


    d.post <- df.mod$d[df.mod$ModelName == model_name]

    d.prior <- rnorm(dim(df.mod)[1], mean = data$priormu[4], sd = sqrt(data$priorSigma[4,4]))

    d.prior <- exp(d.prior)


    df.par.bkg <- data.frame(dist = c(rep("posterior",length(bkg.post)), rep("prior", length(bkg.prior))),
                             value = c(bkg.post, bkg.prior))

    plot.bkg <- ggplot(data = df.par.bkg,
                       aes(x = value, group = as.factor(dist), fill = as.factor(dist))) +
      geom_density(alpha = .3, color = NA) +
      geom_vline(xintercept = obs.bkg, linetype = "dashed", alpha = 0.5) +
      xlab("Background response") + ylab("Density") +
      xlim(quantile(bkg.prior, 0.01), quantile(bkg.prior, 0.99)) +
      scale_fill_manual(name = "", breaks = c("posterior","prior"), labels = c("Posterior","Prior"), values = c("coral","lightblue")) +
      theme_minimal()

    df.par.max <- data.frame(dist = c(rep("posterior",length(max.post)), rep("prior", length(max.prior))),
                             value = c(max.post, max.prior))

    plot.max <- ggplot(data = df.par.max,
                       aes(x = value, group = as.factor(dist), fill = as.factor(dist))) +
      geom_density(alpha = .3, color = NA) +
      geom_vline(xintercept = obs.max, linetype = "dashed", alpha = 0.5) +
      xlab("Maximum response") + ylab("Density") +
      xlim(min(c(quantile(max.post, 0),quantile(max.prior, 0.01))), min(c(quantile(max.post,0.99),quantile(max.prior, 0.99)))) +
      scale_fill_manual(name = "", breaks = c("posterior","prior"), labels = c("Posterior","Prior"), values = c("coral","lightblue")) +
      theme_minimal()

    df.par.fold <- data.frame(dist = c(rep("posterior",length(fold.post)), rep("prior", length(fold.prior))),
                              value = c(fold.post, fold.prior))

    plot.fold <- ggplot(data = df.par.fold,
                        aes(x = value, group = as.factor(dist), fill = as.factor(dist))) +
      geom_density(alpha = .3, color = NA) +
      geom_vline(xintercept = obs.fold, linetype = "dashed", alpha = 0.5) +
      xlab("Fold change") + ylab("Density") +
      xlim(quantile(fold.prior, 0.01), quantile(fold.prior, 0.99)) +
      scale_fill_manual(name = "", breaks = c("posterior","prior"), labels = c("Posterior","Prior"), values = c("coral","lightblue")) +
      theme_minimal()

    df.par.d <- data.frame(dist = c(rep("posterior",length(d.post)), rep("prior", length(d.prior))),
                           value = c(d.post, d.prior))

    plot.d <- ggplot(data = df.par.d,
                     aes(x = value, group = as.factor(dist), fill = as.factor(dist))) +
      geom_density(alpha = .3, color = NA) +
      xlab("Parameter d") + ylab("Density") +
      xlim(quantile(d.prior, 0.01), quantile(d.prior, 0.99)) +
      scale_fill_manual(name = "", breaks = c("posterior","prior"), labels = c("Posterior","Prior"), values = c("coral","lightblue")) +
      theme_minimal()

    df.par.BMD <- data.frame(dist = c(rep("posterior",length(BMD.post)), rep("prior", length(BMD.prior))),
                             value = c(BMD.post*mod.obj$max.dose, BMD.prior*mod.obj$max.dose))

    plot.BMD <- ggplot(data = df.par.BMD,
                       aes(x = value, group = as.factor(dist), fill = as.factor(dist))) +
      # geom_density(aes(x=value, y=..scaled.., group = as.factor(dist), fill = as.factor(dist)), alpha = .3, color = NA) +
      geom_density(alpha = .3, color = NA) +
      xlab("BMD") + ylab("Density") +
      xlim(0,mod.obj$max.dose) +
      scale_fill_manual(name = "", breaks = c("posterior","prior"), labels = c("Posterior","Prior"), values = c("coral","lightblue")) +
      theme_minimal()

    if(parms == T){

      pts <- ggpubr::ggarrange(plot.bkg, plot.fold,
                               plot.d, plot.BMD,
                               nrow = 2, ncol = 2,
                               common.legend = T, legend = "bottom")

    }else if(parms == F){


      pts <- ggpubr::ggarrange(plot.bkg, plot.max,
                               plot.d, plot.BMD,
                               nrow = 2, ncol = 2,
                               common.legend = T, legend = "bottom")

    }


    return(ggpubr::annotate_figure(pts, top = text_grob(paste0("Model ", model_name),
                                                        color = "red", face = "bold", size = 14)))
  }else{
    stop("The requested model was not included in the original analysis.")
  }



}

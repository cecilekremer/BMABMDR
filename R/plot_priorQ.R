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
plot_priorQ <- function(mod.obj, data, model_name){ # pars = T for parameters, F for background & fold change

  names(mod.obj$parsQ) <- c("E4_Q","IE4_Q","H4_Q","LN4_Q","G4_Q","QE4_Q","P4_Q","L4_Q")

  pars <- mod.obj$parsQ
  parsQ_class <- sapply(mod.obj$parsQ, class)
  pos_miss <- which(parsQ_class != "data.frame")
  if(length(pos_miss) != 0) {
    pars <- pars[-pos_miss]
  }

  df.mod <- (as.data.frame(pars[[model_name]]))[, c("ModelName","a","d","BMD")]
  q <- mod.obj$q

  if(model_name %in% c("E4_Q","IE4_Q","H4_Q","LN4_Q","G4_Q","QE4_Q", "P4_Q","L4_Q")){

    obs.a <- mod.obj$data$y[1]/mod.obj$data$n[1]

    a.prior <- rpert(dim(df.mod)[1], min = data$priorlb[1],
                     mode = data$priormu[1], max = data$priorub[1],
                     shape = data$priorgama[1])

    a.post <- df.mod$a[df.mod$ModelName == model_name]
    d.post <- df.mod$d[df.mod$ModelName == model_name]
    d.prior <- rnorm(dim(df.mod)[1], mean = data$priormu[3], sd = sqrt(data$priorSigma[3,3]))
    d.prior <- exp(d.prior)

    BMD.post <- df.mod$BMD[df.mod$ModelName == model_name]
    BMD.prior <- rpert(dim(df.mod)[1], min = data$priorlb[2], mode = data$priormu[2],
                       max = data$priorub[2], shape = data$priorgama[2])

    dens.a.post <- density(a.post, from = min(a.post), to = max(a.post))
    dens.a.prior <- density(a.prior, from = min(a.prior), to = max(a.prior))
    df.par.a <- data.frame(dist = rep(c("posterior", "prior"), each = length(dens.a.post$x)),
                           x = c(dens.a.post$x, dens.a.prior$x),
                           value = c(dens.a.post$y, dens.a.prior$y))

    plot.a <- ggplot(data = df.par.a,
                     aes(x = x, y = value, fill = dist)) +
      geom_area(alpha = .3) +
      xlab("Background response") + ylab("Density") +
      geom_vline(xintercept = obs.a, linetype = "dashed", alpha = 0.5, size = 1) +
      scale_fill_manual(name = "", breaks = c("posterior","prior"),
                        labels = c("Posterior","Prior"), values = c("coral","lightblue"),
                        aesthetics = 'fill') +
      # coord_cartesian(xlim = c(quantile(c(a.prior, a.post), 0.01), quantile(c(a.prior, a.post), 0.99)),
      # ylim = NULL)+
      xlim(c(quantile(c(a.prior, a.post), 0.01), quantile(c(a.prior, a.post), 0.99))) +
      theme_minimal()

    dens.d.post <- density(d.post, from = min(d.post), to = max(d.post))
    dens.d.prior <- density(d.prior, from = min(d.prior), to = max(d.prior))
    df.par.d <- data.frame(dist = rep(c("posterior", "prior"), each = length(dens.d.post$x)),
                           x = c(dens.d.post$x, dens.d.prior$x),
                           value = c(dens.d.post$y, dens.d.prior$y))
    plot.d <- ggplot(data = df.par.d,
                     aes(x = x, y = value, fill = dist)) +
      geom_area(alpha = .3) +
      xlab("Parameter d") + ylab("Density") +
      scale_fill_manual(name = "", breaks = c("posterior","prior"), labels = c("Posterior","Prior"),
                        values = c("coral","lightblue")) +
      # coord_cartesian(xlim = c(quantile(c(d.prior, d.post), 0.01), quantile(d.prior, 0.99))) +
      xlim(quantile(d.prior, 0.01), quantile(d.prior, 0.99)) +
      theme_minimal()


    dens.BMD.post <- density(BMD.post*mod.obj$max.dose,
                             from = min(BMD.post*mod.obj$max.dose),
                             to = max(BMD.post*mod.obj$max.dose))
    dens.BMD.prior <- density(BMD.prior*mod.obj$max.dose,
                              from = min(BMD.prior*mod.obj$max.dose),
                              to = max(BMD.prior*mod.obj$max.dose))

    df.par.BMD <- data.frame(dist = rep(c("posterior", "prior"), each = length(dens.BMD.post$x)),
                             x = c(dens.BMD.post$x, dens.BMD.prior$x),
                             value = c(dens.BMD.post$y, dens.BMD.prior$y))

    plot.BMD <- ggplot(data = df.par.BMD,
                       aes(x = x, y = value, fill = dist)) +
      geom_area(alpha = .3) +
      xlab("Parameter BMD") + ylab("Density") +
      scale_fill_manual(name = "", breaks = c("posterior","prior"), labels = c("Posterior","Prior"),
                        values = c("coral","lightblue")) +
      # coord_cartesian(xlim = c(quantile(c(BMD.prior*mod.obj$max.dose, BMD.post*mod.obj$max.dose), 0.01),
      # quantile(BMD.prior*mod.obj$max.dose, 0.99))) +
      xlim(min(BMD.prior*mod.obj$max.dose), max(BMD.prior*mod.obj$max.dose)) +
      theme_minimal()

  } else stop('model_name not correct. Please check model_name')

  if(mod.obj$is_bin == 1){

    pts <- ggpubr::ggarrange(plot.a,
                             plot.d, plot.BMD,
                             nrow = 2, ncol = 2,
                             common.legend = TRUE, legend = "bottom")

    return(ggpubr::annotate_figure(pts, top = text_grob(paste0("Model ", model_name),
                                                        color = "red", face = "bold", size = 14)))
  } else if(mod.obj$isbetabin == 1){

    rho.post <- df.mod$rho[df.mod$ModelName == model_name]
    dens.rho.post <- density(rho.post*mod.obj$max.dose,
                             from = min(rho.post*mod.obj$max.dose),
                             to = max(rho.post*mod.obj$max.dose))
    # dens.rho.prior <- density(rho.prior*mod.obj$max.dose,
    #                           from = min(rho.prior*mod.obj$max.dose),
    #                           to = max(rho.prior*mod.obj$max.dose))
    df.par.rho <- data.frame(dist = "posterior",
                             x = dens.rho.post$x,
                             value = dens.rho.post$y)
    # df.par.rho <- data.frame(dist = rep(c("posterior", "prior"), each = length(dens.rho.post$x)),
    #                          x = c(dens.rho.post$x, dens.rho.prior$x),
    #                          value = c(dens.rho.post$y, dens.rho.prior$y))

    pts <- ggpubr::ggarrange(plot.a,
                             plot.d, plot.BMD, plot.rho,
                             nrow = 2, ncol = 2,
                             common.legend = TRUE, legend = "bottom")

    return(ggpubr::annotate_figure(pts, top = text_grob(paste0("Model ", model_name),
                                                        color = "red", face = "bold", size = 14)))
  } else stop("The requested model was not included in the original analysis.")

}

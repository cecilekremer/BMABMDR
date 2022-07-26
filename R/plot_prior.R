#' plot function for the prior and posterior distributio of parameter estimates
#'
#' @param mod.obj BMDBMA model object
#' @param data data list from the \code{\link{PREP_DATA_N}} or \code{\link{PREP_DATA_LN}}
#' @param model_name name of the model whose parameters are to be plotted
#' It can be either of c("E4_N","IE4_N","H4_N","LN4_N","G4_N",
#'                                    "QE4_N","P4_N","L4_N",
#'                                     "E4_LN","IE4_LN","H4_LN","LN4_LN",
#'                                     "G4_LN","QE4_LN","P4_LN","L4_LN")
#' @param parms logical indicating if foldchange or maximum response should be plotted. Defaults to TRUE
#' @param clustered logical indicating whether the data are clustered (T) or not (F, default)
#'              
#' @return object of class ggplot
#' 
#' @export plot_prior
#' 
plot_prior <- function(mod.obj, data, model_name, 
                       parms = TRUE, clustered = FALSE){ # pars = T for parameters, F for background & fold change
  
  if(model_name %in% mod.obj$models_included){
    
    names(mod.obj$parsN) = c("E4_N","IE4_N","H4_N","LN4_N","G4_N","QE4_N","P4_N","L4_N")
    names(mod.obj$parsLN) = c("E4_LN","IE4_LN","H4_LN","LN4_LN","G4_LN","QE4_LN","P4_LN","L4_LN")
    
    pars = c(mod.obj$parsN, mod.obj$parsLN)
    

    if(clustered == T){
      df.mod <- (as.data.frame(pars[[model_name]]))[, c("ModelName","a","p3","d","BMD","min_response","max_response","rho")]
      q <- mod.obj$q
      
      obs.bkg <- mean(data$data$response[data$data$dose == data$data$dose[1]])
      
      obs.fold <- mean(data$data$response[data$data$dose == data$data$dose[dim(data$data)[1]]]) / mean(data$data$response[data$data$dose == data$data$dose[1]])
      
      obs.max <- mean(data$data$response[data$data$dose == data$data$dose[dim(data$data)[1]]])
    }else{
      df.mod <- (as.data.frame(pars[[model_name]]))[, c("ModelName","a","p3","d","BMD","min_response","max_response")]
      q <- mod.obj$q
      
      obs.bkg <- mod.obj$data$m[1]
      
      obs.fold <- mod.obj$data$m[length(mod.obj$data$m)] / mod.obj$data$m[1]
      
      obs.max <- mod.obj$data$m[length(mod.obj$data$m)]      
    }
    
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
    
    if(clustered == T){
      rho.post <- df.mod$rho[df.mod$ModelName == model_name]
      
      rho.prior <- rpert(dim(df.mod)[1], min = (data$priorlb[6]),
                         mode = (data$priormu[6]), max = (data$priorub[6]),
                         shape = 0.0001)
    }
    
    
    d.post <- df.mod$d[df.mod$ModelName == model_name]
    
    d.prior <- rnorm(dim(df.mod)[1], mean = data$priormu[4], sd = sqrt(data$priorSigma[4,4]))
    
    if(grepl('QE4_', model_name)){
      d.prior <- rnorm(dim(df.mod)[1], mean = data$priormuQ[4], sd = sqrt(data$priorSigma[4,4]))
    }
    
    d.prior <- exp(d.prior)
    
    
    df.par.bkg <- data.frame(dist = c(rep("posterior",length(bkg.post)), rep("prior", length(bkg.prior))),
                             value = c(bkg.post, bkg.prior))
    
    plot.bkg <- ggplot(data = df.par.bkg,
                       aes(x = value, group = as.factor(dist), fill = as.factor(dist))) +
      geom_density(alpha = .3, color = NA) +
      geom_vline(xintercept = obs.bkg, linetype = "dashed", alpha = 0.5) +
      xlab("Background response") + ylab("Density") +
      #xlim(quantile(bkg.prior, 0.01), quantile(bkg.prior, 0.99)) +
      coord_cartesian(xlim = c(quantile(bkg.prior, 0.01), quantile(bkg.prior, 0.99))) +
      scale_fill_manual(name = "", breaks = c("posterior","prior"), 
                        labels = c("Posterior","Prior"), values = c("coral","lightblue")) +
      theme_minimal()
    
    df.par.max <- data.frame(dist = c(rep("posterior",length(max.post)), rep("prior", length(max.prior))),
                             value = c(max.post, max.prior))
    
    plot.max <- ggplot(data = df.par.max,
                       aes(x = value, group = as.factor(dist), fill = as.factor(dist))) +
      geom_density(alpha = .3, color = NA) +
      geom_vline(xintercept = obs.max, linetype = "dashed", alpha = 0.5) +
      xlab("Maximum response") + ylab("Density") +
      #xlim(min(c(quantile(max.post, 0),quantile(max.prior, 0.01))), min(c(quantile(max.post,0.99),
      #                                                                    quantile(max.prior, 0.99)))) +
      coord_cartesian(xlim = c(min(c(quantile(max.post, 0),quantile(max.prior, 0.01))), 
                               min(c(quantile(max.post,0.99),
                                   quantile(max.prior, 0.99))))) +
      scale_fill_manual(name = "", breaks = c("posterior","prior"), labels = c("Posterior","Prior"), 
                        values = c("coral","lightblue")) +
      theme_minimal()
    
    df.par.fold <- data.frame(dist = c(rep("posterior",length(fold.post)), rep("prior", length(fold.prior))),
                              value = c(fold.post, fold.prior))
    
    plot.fold <- ggplot(data = df.par.fold,
                        aes(x = value, group = as.factor(dist), fill = as.factor(dist))) +
      geom_density(alpha = .3, color = NA) +
      geom_vline(xintercept = obs.fold, linetype = "dashed", alpha = 0.5) +
      xlab("Fold change") + ylab("Density") +
      #xlim(quantile(fold.prior, 0.01), quantile(fold.prior, 0.99)) +
      coord_cartesian(xlim = c(quantile(fold.prior, 0.01), quantile(fold.prior, 0.99))) +
      scale_fill_manual(name = "", breaks = c("posterior","prior"), labels = c("Posterior","Prior"), 
                        values = c("coral","lightblue")) +
      theme_minimal()
    
    df.par.d <- data.frame(dist = c(rep("posterior",length(d.post)), rep("prior", length(d.prior))),
                           value = c(d.post, d.prior))
    
    plot.d <- ggplot(data = df.par.d,
                     aes(x = value, group = as.factor(dist), fill = as.factor(dist))) +
      geom_density(alpha = .3, color = NA) +
      xlab("Parameter d") + ylab("Density") +
      #xlim(quantile(d.prior, 0.01), quantile(d.prior, 0.99)) +
      coord_cartesian(xlim = c(quantile(d.prior, 0.01), quantile(d.prior, 0.99))) +
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
      # coord_cartesian(xlim = c(0,mod.obj$max.dose*2)) +
      scale_fill_manual(name = "", breaks = c("posterior","prior"), labels = c("Posterior","Prior"), 
                        values = c("coral","lightblue")) +
      theme_minimal()
    
    if(clustered == T){
      df.par.rho <- data.frame(dist = c(rep("posterior",length(rho.post)), rep("prior", length(rho.prior))),
                               value = c(rho.post, rho.prior))
      
      plot.rho <- ggplot(data = df.par.rho,
                       aes(x = value, group = as.factor(dist), fill = as.factor(dist))) +
        geom_density(alpha = .3, color = NA) +
        xlab("Correlation rho") + ylab("Density") +
        #xlim(quantile(d.prior, 0.01), quantile(d.prior, 0.99)) +
        coord_cartesian(xlim = c(quantile(rho.prior, 0.01), quantile(rho.prior, 0.99))) +
        scale_fill_manual(name = "", breaks = c("posterior","prior"), labels = c("Posterior","Prior"), values = c("coral","lightblue")) +
        theme_minimal()
    }
    
    if(parms == T){
      
      pts <- ggpubr::ggarrange(plot.bkg, plot.fold,  
                               plot.d, plot.BMD,
                               nrow = 2, ncol = 2,
                               common.legend = T, legend = "bottom")  
      
      if(clustered == T){
        pts <- ggpubr::ggarrange(plot.bkg, plot.fold,  
                                 plot.d, plot.BMD,
                                 plot.rho,
                                 nrow = 3, ncol = 2,
                                 common.legend = T, legend = "bottom")  
      }
      
    }else if(parms == F){
      
      
      pts <- ggpubr::ggarrange(plot.bkg, plot.max,  
                               plot.d, plot.BMD,
                               nrow = 2, ncol = 2,
                               common.legend = T, legend = "bottom")  
      
      if(clustered == T){
        pts <- ggpubr::ggarrange(plot.bkg, plot.max,  
                                 plot.d, plot.BMD,
                                 plot.rho,
                                 nrow = 3, ncol = 2,
                                 common.legend = T, legend = "bottom")  
      }
      
    }
    
    
    return(ggpubr::annotate_figure(pts, top = text_grob(names(get_models('continuous')[which(get_models('continuous') == model_name)]), 
                                                        color = "red", face = "bold", size = 14)))
  }else{
    stop("The requested model was not included in the original analysis.")
  }
  
  
  
}


#' @rdname plot_prior
#' @export
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
    
    a.prior <- mc2d::rpert(dim(df.mod)[1], min = data$priorlb[1],
                           mode = data$priormu[1], max = data$priorub[1],
                           shape = data$priorgama[1])
    
    a.post <- df.mod$a[df.mod$ModelName == model_name]
    d.post <- df.mod$d[df.mod$ModelName == model_name]
    d.prior <- rnorm(dim(df.mod)[1], mean = data$priormu[3], sd = sqrt(data$priorSigma[3,3]))
    d.prior <- exp(d.prior)
    
    BMD.post <- df.mod$BMD[df.mod$ModelName == model_name]
    BMD.prior <- mc2d::rpert(dim(df.mod)[1], min = data$priorlb[2], mode = data$priormu[2], 
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
      coord_cartesian(xlim = c(quantile(c(a.prior, a.post), 0.01), quantile(c(a.prior, a.post), 0.99)),
                      ylim = NULL)+
      #xlim(c(quantile(c(a.prior, a.post), 0.01), quantile(c(a.prior, a.post), 0.99))) + 
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
      coord_cartesian(xlim = c(quantile(c(d.prior, d.post), 0.01), quantile(d.prior, 0.99))) +
      #xlim(quantile(d.prior, 0.01), quantile(d.prior, 0.99)) +
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
      coord_cartesian(xlim = c(quantile(c(BMD.prior*mod.obj$max.dose, BMD.post*mod.obj$max.dose), 0.01),
                               min(quantile(c(BMD.prior*mod.obj$max.dose, BMD.post*mod.obj$max.dose), 0.99), 
                                   2*mod.obj$max.dose)
      )
      ) +
      #xlim(min(BMD.prior*mod.obj$max.dose), max(BMD.prior*mod.obj$max.dose)) + 
      theme_minimal()
    
  } else stop('model_name not correct. Please check model_name')
  
  if(mod.obj$is_bin == 1){
    
    pts <- ggpubr::ggarrange(plot.a,  
                             plot.d, plot.BMD,
                             nrow = 2, ncol = 2,
                             common.legend = TRUE, legend = "bottom")
    
    return(ggpubr::annotate_figure(pts, top = text_grob(paste0("Model ", model_name), 
                                                        color = "red", face = "bold", size = 14)))
  } else if(mod.obj$is_bin == 0){
    
    rho.prior <- mc2d::rpert(dim(df.mod)[1], min = 0,
                             mode = data$priormu[4], max = 1,
                             shape = 4)
    dens.rho.prior <- density(rho.prior, from = min(rho.prior), to = max(rho.prior))
    
    
    df.mod <- (as.data.frame(pars[[model_name]]))[, c("ModelName","a","d","BMD", "rho")]
    rho.post <- df.mod$rho[df.mod$ModelName == model_name]
    dens.rho.post <- density(rho.post, 
                             from = min(rho.post), 
                             to = max(rho.post))
    
    df.par.rho <- data.frame(dist = rep(c("posterior", "prior"), each = length(dens.rho.post$x)),
                             x = c(dens.rho.post$x, dens.rho.prior$x),
                             value = c(dens.rho.post$y, dens.rho.prior$y))
    
    plot.rho <- ggplot(data = df.par.rho,
                       aes(x = x, y = value, fill = dist)) +
      geom_area(alpha = .3) +
      xlab("Parameter rho") + ylab("Density") +
      scale_fill_manual(name = "", breaks = c("posterior","prior"), labels = c("Posterior","Prior"), 
                        values = c("coral","lightblue")) +
      coord_cartesian(xlim = c(quantile(c(rho.prior, rho.post), 0.01),
                               quantile(rho.prior, 0.99))) +
      #xlim(min(BMD.prior*mod.obj$max.dose), max(BMD.prior*mod.obj$max.dose)) + 
      theme_minimal()
    
    pts <- ggpubr::ggarrange(plot.a,  
                             plot.d, plot.BMD, plot.rho,
                             nrow = 2, ncol = 2,
                             common.legend = TRUE, legend = "bottom")
    
    return(ggpubr::annotate_figure(pts, top = text_grob(names(get_models('quantal')[
      which(get_models('quantal') == model_name)]), color = "red", face = "bold", size = 14)))
  } else stop("The requested model was not included in the original analysis.")
  
}

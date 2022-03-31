#' Function to plot prior vs posterior distribution
#'
#' @param mod.obj the model object returned by full laplace or sampling method
#' @param data the data for a specific model used by the above method
#' @param model_name the model for which plots should be generated
#' @param priordist the prior distribution used
#' @param priorbmd logical indicating if informative prior for BMD was used
#'
#' @return .
#'
#' @export
#'
plot_prior <- function(mod.obj, data, model_name, parms, priordist, priorbmd){ # pars = T for parameters, F for background & fold change

  if(model_name %in% mod.obj$models_included){

    names(mod.obj$parsN) = c("E4_N","IE4_N","H4_N","LN4_N","G4_N","QE4_N","P4_N","L4_N")
    names(mod.obj$parsLN) = c("E4_LN","IE4_LN","H4_LN","LN4_LN","G4_LN","QE4_LN","P4_LN","L4_LN")

    pars = c(mod.obj$parsN, mod.obj$parsLN)

    df.mod <- (as.data.frame(pars[[model_name]]))[, c("ModelName","a","c","d","BMD","min_response","max_response")]
    q <- mod.obj$q

    if(mod.obj$increasing == TRUE){

      if(model_name %in% c("E4_N","IE4_N","H4_N","LN4_N","G4_N","QE4_N","E4_LN","IE4_LN","H4_LN","LN4_LN","G4_LN","QE4_LN")){

        obs.a <- mod.obj$data$m[1]
        if(grepl("_LN", model_name)){
          obs.a <- log(obs.a)
        }

        obs.c <- mod.obj$data$m[length(mod.obj$data$m)] / mod.obj$data$m[1]
        if(grepl("_LN", model_name)){
          obs.c <- log(mod.obj$data$m[length(mod.obj$data$m)]) / log(mod.obj$data$m[1])
        }

        obs.bkg <- obs.a

        obs.fold <- obs.c

        if(priordist == "PERT"){
          a.prior <- rpert(dim(df.mod)[1], min = (data$priorlb[1]),
                           mode = (data$priormu[1]), max = (data$priorub[1]),
                           shape = 4)

          c.prior <- rpert(dim(df.mod)[1], min = (data$priorlb[3]),
                           mode = (data$priormu[3]), max = (data$priorub[3]),
                           shape = 4)
        }else if(priordist == "Normal"){
          a.prior <- rnorm(dim(df.mod)[1], mean = data$priormu[1], sd = sqrt(data$priorSigma[1,1]))

          c.prior <- rnorm(dim(df.mod)[1], mean = data$priormu[3], sd = sqrt(data$priorSigma[3,3]))

        }

        a.post <- df.mod$a[df.mod$ModelName == model_name]

        c.post <- df.mod$c[df.mod$ModelName == model_name]

        if(model_name == "LN4_N"){

          c.prior <- exp(c.prior)+1+q

        }else if(model_name == "LN4_LN"){

          c.prior <- exp(c.prior)+(1+(log(1+q)/exp(data$priormu[1])))

        }else{

          c.prior <- exp(c.prior) + 1

        }

        a.prior <- exp(a.prior)

        d.post <- df.mod$d[df.mod$ModelName == model_name]

        d.prior <- rnorm(dim(df.mod)[1], mean = data$priormu[4], sd = sqrt(data$priorSigma[4,4]))

        d.prior <- exp(d.prior)

        BMD.post <- df.mod$BMD[df.mod$ModelName == model_name]

        if(priorbmd == FALSE){
          BMD.prior <- rexp(dim(df.mod)[1], 1)

          BMD.prior <- exp(-BMD.prior)
        }else if(priorbmd == TRUE){
          BMD.prior <- rpert(dim(df.mod)[1], min = data$priorlb[2], mode = data$priormu[2], max = data$priorub[2], shape = data$shape.BMD)
          BMD.prior <- exp(BMD.prior)
        }


        ### Other

        bkg.post <- df.mod$min_response[df.mod$ModelName == model_name]

        fold.post <- df.mod$max_response[df.mod$ModelName == model_name] / df.mod$min_response[df.mod$ModelName == model_name]

        bkg.prior <- a.prior

        fold.prior <- c.prior

      }else if(model_name %in% c("P4_N","L4_N","P4_LN","L4_LN")){

        if(grepl("P4", model_name)){
          funF2 <- qnorm
          funF2inv <- pnorm
        }else if(grepl("L4", model_name)){
          funF2 <- logit
          funF2inv <- expit
        }

        obs.a <- mod.obj$data$m[length(mod.obj$data$m)]
        if(grepl("_LN", model_name)){
          obs.a <- log(obs.a)
        }

        if(grepl("_N", model_name)){
          obs.c <- funF2(mod.obj$data$m[1] / mod.obj$data$m[length(mod.obj$data$m)])
        }else if(grepl("_LN", model_name)){
          obs.c <- funF2(log(mod.obj$data$m[1]) / log(mod.obj$data$m[length(mod.obj$data$m)]))
        }

        obs.bkg <- mod.obj$data$m[1]
        obs.fold <- mod.obj$data$m[length(mod.obj$data$m)] / mod.obj$data$m[1]
        if(grepl("_LN", model_name)){
          obs.fold <- log(mod.obj$data$m[length(mod.obj$data$m)]) / log(mod.obj$data$m[1])
          obs.bkg <- log(mod.obj$data$m[1])
        }

        if(priordist == "PERT"){
          a.prior <- rpert(dim(df.mod)[1], min = (data$priorlb[1]),
                           mode = (data$priormu[1]), max = (data$priorub[1]),
                           shape = 4)

          c.prior <- rpert(dim(df.mod)[1], min = (data$priorlb[3]),
                           mode = (data$priormu[3]), max = (data$priorub[3]),
                           shape = 4)
        }else if(priordist == "Normal"){

          a.prior <- rnorm(dim(df.mod)[1], mean = data$priormu[1], sd = sqrt(data$priorSigma[1,1]))

          c.prior <- rnorm(dim(df.mod)[1], mean = data$priormu[3], sd = sqrt(data$priorSigma[3,3]))

        }

        a.post <- df.mod$a[df.mod$ModelName == model_name]

        c.post <- df.mod$c[df.mod$ModelName == model_name]

        a.prior <- exp(a.prior)

        if(grepl("_N", model_name)){
          c.prior <- funF2(1/(1+q)) - exp(c.prior)
        }else if(grepl("_LN", model_name)){
          c.prior <- funF2(1 - (log(1+q)/exp(data$priormu[1]))) - exp(c.prior)
        }


        d.post <- df.mod$d[df.mod$ModelName == model_name]

        d.prior <- rnorm(dim(df.mod)[1], mean = data$priormu[4], sd = sqrt(data$priorSigma[4,4]))

        d.prior <- exp(d.prior)

        BMD.post <- df.mod$BMD[df.mod$ModelName == model_name]

        if(priorbmd == FALSE){
          BMD.prior <- rexp(dim(df.mod)[1], 1)
          BMD.prior <- exp(-BMD.prior)
        }else if(priorbmd == TRUE){
          BMD.prior <- rpert(dim(df.mod)[1], min = data$priorlb[2], mode = data$priormu[2], max = data$priorub[2], shape = data$shape.BMD)
          BMD.prior <- exp(BMD.prior)
        }

        ### Other

        bkg.post <- df.mod$min_response[df.mod$ModelName == model_name]

        fold.post <- df.mod$max_response[df.mod$ModelName == model_name] / df.mod$min_response[df.mod$ModelName == model_name]

        bkg.prior <- a.prior*funF2inv(c.prior)

        fold.prior <- a.prior/bkg.prior

      }

    }else if(mod.obj$increasing == FALSE){


      if(model_name %in% c("E4_N","IE4_N","H4_N","LN4_N","G4_N","QE4_N","E4_LN","IE4_LN","H4_LN","LN4_LN","G4_LN","QE4_LN")){

        obs.a <- mod.obj$data$m[1]
        if(grepl("_LN", model_name)){
          obs.a <- log(obs.a)
        }

        obs.c <- mod.obj$data$m[length(mod.obj$data$m)] / mod.obj$data$m[1]
        if(grepl("_LN", model_name)){
          obs.c <- log(mod.obj$data$m[length(mod.obj$data$m)]) / log(mod.obj$data$m[1])
        }

        obs.bkg <- obs.a

        obs.fold <- obs.c


        if(priordist == "PERT"){
          a.prior <- rpert(dim(df.mod)[1], min = (data$priorlb[1]),
                           mode = (data$priormu[1]), max = (data$priorub[1]),
                           shape = 4)

          c.prior <- rpert(dim(df.mod)[1], min = (data$priorlb[3]),
                           mode = (data$priormu[3]), max = (data$priorub[3]),
                           shape = 4)
        }else if(priordist == "Normal"){

          a.prior <- rnorm(dim(df.mod)[1], mean = data$priormu[1], sd = sqrt(data$priorSigma[1,1]))

          c.prior <- rnorm(dim(df.mod)[1], mean = data$priormu[3], sd = sqrt(data$priorSigma[3,3]))

        }

        a.post <- df.mod$a[df.mod$ModelName == model_name]

        c.post <- df.mod$c[df.mod$ModelName == model_name]

        if(model_name == "LN4_N"){

          c.prior <- (1-q)*(1/(1+exp(-c.prior)))

        }else if(model_name == "LN4_LN"){

          c.prior <- (1+(log(1-q)/exp(data$priormu[1])))*(1/(1+exp(-c.prior)))

        }else{

          c.prior <- 1/(1+exp(-c.prior))

        }

        a.prior <- exp(a.prior)

        d.post <- df.mod$d[df.mod$ModelName == model_name]

        d.prior <- rnorm(dim(df.mod)[1], mean = data$priormu[4], sd = sqrt(data$priorSigma[4,4]))

        d.prior <- exp(d.prior)

        BMD.post <- df.mod$BMD[df.mod$ModelName == model_name]

        if(priorbmd == FALSE){
          BMD.prior <- rexp(dim(df.mod)[1], 1)
          BMD.prior <- exp(-BMD.prior)
        }else if(priorbmd == TRUE){
          BMD.prior <- rpert(dim(df.mod)[1], min = data$priorlb[2], mode = data$priormu[2], max = data$priorub[2], shape = data$shape.BMD)
          BMD.prior <- exp(BMD.prior)
        }

        ### Other

        bkg.post <- df.mod$min_response[df.mod$ModelName == model_name]

        fold.post <- df.mod$max_response[df.mod$ModelName == model_name] / df.mod$min_response[df.mod$ModelName == model_name]

        bkg.prior <- a.prior

        fold.prior <- c.prior

      }else if(model_name %in% c("P4_N","L4_N","P4_LN","L4_LN")){

        if(grepl("P4", model_name)){
          funF2 <- qnorm
          funF2inv <- pnorm
        }else if(grepl("L4", model_name)){
          funF2 <- logit
          funF2inv <- expit
        }

        obs.a <- mod.obj$data$m[1]
        if(grepl("_LN", model_name)){
          obs.a <- log(obs.a)
        }

        if(grepl("_N", model_name)){
          obs.c <- funF2(mod.obj$data$m[length(mod.obj$data$m)] / mod.obj$data$m[1])
        }else if(grepl("_LN", model_name)){
          obs.c <- funF2(log(mod.obj$data$m[length(mod.obj$data$m)]) / log(mod.obj$data$m[1]))
        }

        obs.bkg <- obs.a

        obs.fold <- mod.obj$data$m[length(mod.obj$data$m)] / mod.obj$data$m[1]
        if(grepl("_LN", model_name)){
          obs.fold <- log(mod.obj$data$m[length(mod.obj$data$m)]) / log(mod.obj$data$m[1])
        }

        if(priordist == "PERT"){
          a.prior <- rpert(dim(df.mod)[1], min = (data$priorlb[1]),
                           mode = (data$priormu[1]), max = (data$priorub[1]),
                           shape = 4)

          c.prior <- rpert(dim(df.mod)[1], min = (data$priorlb[3]),
                           mode = (data$priormu[3]), max = (data$priorub[3]),
                           shape = 4)

        }else if(priordist == "Normal"){

          a.prior <- rnorm(dim(df.mod)[1], mean = data$priormu[1], sd = sqrt(data$priorSigma[1,1]))

          c.prior <- rnorm(dim(df.mod)[1], mean = data$priormu[3], sd = sqrt(data$priorSigma[3,3]))

        }

        a.post <- df.mod$a[df.mod$ModelName == model_name]

        c.post <- df.mod$c[df.mod$ModelName == model_name]

        a.prior <- exp(a.prior)

        if(grepl("_N", model_name)){
          c.prior <- funF2(1-q) - exp(c.prior)
        }else if(grepl("_LN", model_name)){
          c.prior <- funF2(1 + (log(1-q)/exp(data$priormu[1]))) - exp(c.prior)
        }


        d.post <- df.mod$d[df.mod$ModelName == model_name]

        d.prior <- rnorm(dim(df.mod)[1], mean = data$priormu[4], sd = sqrt(data$priorSigma[4,4]))

        d.prior <- exp(d.prior)

        BMD.post <- df.mod$BMD[df.mod$ModelName == model_name]

        if(priorbmd == FALSE){
          BMD.prior <- rexp(dim(df.mod)[1], 1)
          BMD.prior <- exp(-BMD.prior)
        }else if(priorbmd == TRUE){
          BMD.prior <- rpert(dim(df.mod)[1], min = data$priorlb[2], mode = data$priormu[2], max = data$priorub[2], shape = data$shape.BMD)
          BMD.prior <- exp(BMD.prior)
        }

        ### Other

        bkg.post <- df.mod$min_response[df.mod$ModelName == model_name]

        fold.post <- df.mod$max_response[df.mod$ModelName == model_name] / df.mod$min_response[df.mod$ModelName == model_name]

        bkg.prior <- a.prior

        fold.prior <- (a.prior*funF2inv(c.prior))/a.prior

      }


    }


    df.par.a <- data.frame(dist = c(rep("posterior",length(a.post)), rep("prior", length(a.prior))),
                           value = c(a.post, a.prior))

    plot.a <- ggplot(data = df.par.a,
                     aes(x = value, group = as.factor(dist), fill = as.factor(dist))) +
      geom_density(alpha = .3, color = NA) +
      geom_vline(xintercept = obs.a, linetype = "dashed", alpha = 0.5) +
      xlab("Parameter a") + ylab("Density") +
      xlim(quantile(a.prior, 0.01), quantile(a.prior, 0.99)) +
      scale_fill_manual(name = "", breaks = c("posterior","prior"), labels = c("Posterior","Prior"), values = c("coral","lightblue")) +
      theme_minimal()

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

    df.par.c <- data.frame(dist = c(rep("posterior",length(c.post)), rep("prior", length(c.prior))),
                           value = c(c.post, c.prior))

    plot.c <- ggplot(data = df.par.c,
                     aes(x = value, group = as.factor(dist), fill = as.factor(dist))) +
      geom_density(alpha = .3, color = NA) +
      geom_vline(xintercept = obs.c, linetype = "dashed", alpha = 0.5) +
      xlab("Parameter c") + ylab("Density") +
      xlim(quantile(c.prior, 0.01), quantile(c.prior, 0.99)) +
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
                             value = c(BMD.post, BMD.prior))

    # plot.BMD <- ggplot(data = df.par.BMD,
    #                    aes(x = value, group = as.factor(dist), fill = as.factor(dist))) +
    #   # geom_density(aes(x=value, y=..scaled.., group = as.factor(dist), fill = as.factor(dist)), alpha = .3, color = NA) +
    #   geom_density(alpha = .3, color = NA) +
    #   xlab("BMD") + ylab("Density") +
    #   xlim(0,1) +
    #   scale_fill_manual(name = "", breaks = c("posterior","prior"), labels = c("Posterior","Prior"), values = c("coral","lightblue")) +
    #   theme_minimal()

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

      pts <- ggpubr::ggarrange(plot.a, plot.c,
                               plot.d, plot.BMD,
                               nrow = 2, ncol = 2,
                               common.legend = T, legend = "bottom")

    }else if(parms == F){

      # pts <- ggpubr::ggarrange(plot.bkg, NULL, plot.fold, NULL, plot.BMD, NULL, ncol = 3, nrow = 2,
      #                          widths = c(2,0,2,1,2,1),common.legend = T)

      pts <- ggpubr::ggarrange(plot.bkg, plot.fold,
                               plot.BMD,
                               nrow = 2, ncol = 2,
                               common.legend = T, legend = "bottom")

    }


    return(ggpubr::annotate_figure(pts, top = ggpubr::text_grob(paste0("Model ", model_name),
                                                        color = "red", face = "bold", size = 14)))
  }else{
    stop("The requested model was not included in the original analysis.")
  }



}

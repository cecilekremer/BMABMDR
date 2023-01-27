
#' Function for internal use
#'
#' @param DR_df dataframe containing individual-level dose response data
#' @param DRM_df dataframe for dose-response model
#' @param x_DR_df name of the x varaibel in the individual-level dose response data.
#'                Typically dose.
#' @param y_DR_df name of the y variable in the individual-level dose response data.
#' @param x_DRM_df name of the x variable in the dose response curve data. Typically dose
#' @param y_DRM_df name of the y variable in the individual-level dose response data.
#' @param col_DRM color scheme
#'
#' @return .
#'
#' @export DR_scatter
#'
DR_scatter <- function(DR_df = NULL, DRM_df = NULL,
                       x_DR_df, y_DR_df,
                       x_DRM_df, y_DRM_df,
                       col_DRM = brewer.pal(3, "Set1")[1]) {

  if(!is.null(DR_df)) {

    base_plt <- ggplot(data = DR_df, aes_string(x = x_DR_df,
                                                y = y_DR_df)) +
      geom_point(size = 4, shape = 19,
                 alpha = 0.4) +
      labs(x = 'Dose', y = 'Response') +
      theme_minimal()
  }

  mx_x <- max(DR_df[, x_DR_df])

  if(!is.null(DR_df) & is.null(DRM_df)) {

    p1 <-  base_plt +
      #scale_x_continuous(breaks = seq(0, mx_x, leng)) +
      theme(strip.text = element_text(size = 15, face = "bold"),
            axis.text = element_text(size = 10, face = "bold"),
            axis.title = element_text(size = 15, face = "bold"),
            legend.text = element_text(size = 10, face = "bold"),
            legend.title = element_text(size = 15, face = "bold"),
            panel.spacing = unit(5, "lines"),
            legend.position = "top",
            legend.direction = "horizontal",
            title = element_text(size = 15, face = "bold"))

  } else if(is.null(DR_df) & !is.null(DRM_df)) {

    mx_x <- max(DRM_df[, x_DRM_df])

    p1 <- ggplot(data = DRM_df, aes_string(x = x_DRM_df,
                                           y = y_DRM_df,
                                           group = 1)) +
      geom_line(size = 3.5, color = col_DRM) +
      labs(x = 'Dose', y = 'Response') +
      theme_minimal() +
      #scale_x_continuous(breaks = seq(0, mx_x, by = 0.1)) +
      theme(strip.text = element_text(size = 15, face = "bold"),
            axis.text = element_text(size = 10, face = "bold"),
            axis.title = element_text(size = 20, face = "bold"),
            legend.text = element_text(size = 10, face = "bold"),
            legend.title = element_text(size = 20, face = "bold"),
            panel.spacing = unit(5, "lines"),
            legend.position = "top",
            legend.direction = "horizontal",
            title = element_text(size = 20, face = "bold"))

  } else if(!is.null(DR_df) & !is.null(DRM_df)) {

    p1 <- base_plt +
      geom_line(data = DRM_df,
                aes_string(x = x_DRM_df,
                           y = y_DRM_df,
                           group = 1),
                color = col_DRM, size = 3.5) +
      labs(x = 'Dose', y = 'Response') +
      theme_minimal() +
      #scale_x_continuous(breaks = seq(0, mx_x, by = 0.1)) +
      theme(strip.text = element_text(size = 15, face = "bold"),
            axis.text = element_text(size = 10, face = "bold"),
            axis.title = element_text(size = 20, face = "bold"),
            legend.text = element_text(size = 10, face = "bold"),
            legend.title = element_text(size = 20, face = "bold"),
            panel.spacing = unit(5, "lines"),
            legend.position = "top",
            legend.direction = "horizontal",
            plot.title = element_text(hjust = 0.5),
            title = element_text(size = 20, face = "bold"))

  } else return(stop("Check the inputs"))

  return(p1)

}

#' Function for internal use
#'
#' @param BMD_DF dataframe containing BMD, BMDL and BMDU
#' @param BMD.true true value for BMD. Only necessary if such exist
#' @param x value
#' @param BMD BMD
#' @param BMDL BMDL
#' @param BMDU BMDU
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @param title title
#'
#' @return .
#'
#' @export plot_BMD
#'
plot_BMD <- function(BMD_DF, BMD.true = NULL,
                     x, BMD, BMDL, BMDU,
                     xlab = "Dataset",
                     ylab = expression(log[10](BMD)),
                     title = "BMD Estimate") {

  p1 <- ggplot(data = BMD_DF, aes_string(x = x, y = BMD, group = 1)) +
    geom_errorbar(aes_string(ymin = BMDL,
                             ymax = BMDU),
                  position = position_dodge(0.4),
                  width = 0.1, size = 1.3) +
    geom_point(size = 6,
               position = position_dodge(0.8),
               color = RColorBrewer::brewer.pal(4, "Set1")[2]) +
    labs(x = xlab, y = ylab, title = title) +
    theme_minimal() +
    theme(strip.text = element_text(size = 15, face = "bold"),
          axis.text = element_text(size = 10, face = "bold"),
          axis.title = element_text(size = 15, face = "bold"),
          legend.text = element_text(size = 10, face = "bold"),
          legend.title = element_text(size = 15, face = "bold"),
          panel.spacing = unit(5, "lines"),
          legend.position = "top",
          legend.direction = "horizontal",
          plot.title = element_text(hjust = 0.5),
          title = element_text(size = 15, face = "bold"))

  if(!is.null(BMD.true)) {

    p1 <- p1 + geom_hline(yintercept = BMD.true,
                          linetype = "dashed", size = 3)  #+
    #scale_color_manual(values = colls)
  }

  return(p1)
}

#' Function for internal use
#'
#' @param BMD_DF dataframe containing BMD, BMDL and BMDU
#' @param BMD.true true value for BMD. Only necessary if such exist
#' @param to_plot which value to plot
#' @param xlab x-axis label
#' @param title title
#'
#' @return .
#'
#' @export density_BMD
#'
density_BMD <- function(BMD_DF, BMD.true = NULL,
                        to_plot = c("BMD",
                                    "BMDL",
                                    "BMDU"),
                        xlab = expression(log[10](BMD)),
                        title = "MLE") {
  # meds <- median(BMD_DF[,  to_plot])
  p1 <- ggplot(data = BMD_DF, aes_string(x = to_plot)) +
    geom_density(fill = "#377EB8") +
    theme_minimal() +
    labs(x = xlab, y = "Density",
         title = title) +
    # geom_vline(xintercept = 1,
    #            linetype = "dotted", size = 3,
    #            color = "#4DAF4A") +
    theme(strip.text = element_text(size = 15, face = "bold"),
          axis.text = element_text(size = 10, face = "bold"),
          axis.title = element_text(size = 15, face = "bold"),
          legend.text = element_text(size = 10, face = "bold"),
          legend.title = element_text(size = 15, face = "bold"),
          panel.spacing = unit(5, "lines"),
          legend.position = "top",
          legend.direction = "horizontal",
          plot.title = element_text(hjust = 0.5),
          title = element_text(size = 15, face = "bold"))

  if(!is.null(BMD.true)) {
    p1 <- p1 + geom_vline(xintercept = BMD.true,
                          linetype = "dashed", size = 3)

  }#else{
  #  p1 <- p1 #+ geom_vline(xintercept = 1,
  #          #                linetype = "dotted", size = 3,
  #          #                color = "#4DAF4A")
  #}

  return(p1)

}


#' Function for internal use
#'
#' @param Model which model
#' @param type increasing or decreasing
#'
#' @return .
#'
#' @export DRM
#'
DRM <- function(Model, type = c("increasing", "decreasing", "quantal")) {

  type <- match.arg(type)
  if(type == "increasing") {

    DRM <- switch (Model,
                   E4_N = DRM.E4_NI,
                   E4_LN = DRM.E4_LNI,
                   IE4_N = DRM.IE4_NI,
                   IE4_LN = DRM.IE4_LNI,
                   H4_N = DRM.H4_NI,
                   H4_LN = DRM.H4_LNI,
                   LN4_N = DRM.LN4_NI,
                   LN4_LN = DRM.LN4_LNI,
                   G4_N = DRM.G4_NI,
                   G4_LN = DRM.G4_LNI,
                   QE4_N = DRM.QE4_NI,
                   QE4_LN = DRM.QE4_LNI,
                   P4_N = DRM.P4_NI,
                   P4_LN = DRM.P4_LNI,
                   L4_N = DRM.L4_NI,
                   L4_LN = DRM.L4_LNI
    )

  } else if(type == "decreasing") {

    DRM <- switch (Model,
                   E4_N = DRM.E4_ND,
                   E4_LN = DRM.E4_LND,
                   IE4_N = DRM.IE4_ND,
                   IE4_LN = DRM.IE4_LND,
                   H4_N = DRM.H4_ND,
                   H4_LN = DRM.H4_LND,
                   LN4_N = DRM.LN4_ND,
                   LN4_LN = DRM.LN4_LND,
                   G4_N = DRM.G4_ND,
                   G4_LN = DRM.G4_LND,
                   QE4_N = DRM.QE4_ND,
                   QE4_LN = DRM.QE4_LND,
                   P4_N = DRM.P4_ND,
                   P4_LN = DRM.P4_LND,
                   L4_N = DRM.L4_ND,
                   L4_LN = DRM.L4_LND
    )

  } else if(type == 'quantal') {

    DRM <- switch (Model,
                   E4_Q = DRM.E4_Q,
                   IE4_Q = DRM.IE4_Q,
                   H4_Q = DRM.H4_Q,
                   LN4_Q = DRM.LN4_Q,
                   G4_Q = DRM.G4_Q,
                   QE4_Q = DRM.QE4_Q,
                   P4_Q = DRM.P4_Q,
                   L4_Q = DRM.L4_Q
    )


  } else stop("Provide type to be either 'increasing', 'decreasing' or 'quantal'")
  return(DRM)
}


#' Function to plot the results for continuous endpoints
#'
#' @param mod.obj BMDBMA model object
#' @param type dose-response type. It can either be "increasing" or "decreasing"
#' @param clustered logical indicating whether clustered data is used
#' @param weight_type type of model weight to be plotted. It can either be "BS" or "LP"
#' @param include_data logical argument to indicate if data should be included in the plots. Defaults to TRUE
#' @param all logical argument to indicate if all plots should be displayed
#' @param title title of the plot
#' @param log logical whether the fit of Normal models should be on log10-scale (log = T) or original scale (log = F)
#'
#' @return .
#'
#' @export plot.BMADR
#'
plot.BMADR <- function(mod.obj,
                       type = c("increasing", "decreasing"),
                       clustered = FALSE,
                       weight_type = c("BS", "LP"),
                       include_data = TRUE,
                       all = TRUE, title, log = FALSE
) {
  type <- match.arg(type)
  weight_type <- match.arg(weight_type)
  q <- mod.obj$q
  # if(mod.obj$increasing == T){
  #   type = "increasing"
  # }else{type = "decreasing"}

  refactor <- function(x, type){

    if(type == "quantal") {

      x <- dplyr::mutate(x, Model = factor(x$Model,
                                           levels = c("Model Averaged", "E4_Q", "IE4_Q", "H4_Q",
                                                      "LN4_Q", "G4_Q",  "QE4_Q", "P4_Q", "L4_Q"),
                                           labels = c("Model Averaged", "Exp", "InvExp", "Hill",
                                                      "LogNormal", "Gamma",  "QuadExp", "Probit",
                                                      "Logistic")
      ))

    } else {

      x <- dplyr::mutate(x, Model = factor(x$Model,
                                           levels = c("Model Averaged", "E4_N", "IE4_N", "H4_N",
                                                      "LN4_N", "G4_N",  "QE4_N", "P4_N", "L4_N",
                                                      "E4_LN", "IE4_LN", "H4_LN", "LN4_LN",
                                                      "G4_LN", "QE4_LN", "P4_LN", "L4_LN"),
                                           labels = c("Model Averaged", "Exp(N)", "InvExp(N)", "Hill(N)",
                                                      "LogNormal(N)", "Gamma(N)",  "QuadExp(N)", "Probit(N)",
                                                      "Logistic(N)", "Exp(LN)", "InvExp(LN)", "Hill4(LN)",
                                                      "LogNormal(LN)", "Gamma(LN)", "QuadExp(LN)",
                                                      "Probit(LN)", "Logistic(LN)")
      ))
    }

    return(x)

  }

  BMDW <- BMDWeights(mod.obj)
  BMDBMA <- BMDMA_extract(mod.obj, conv = FALSE)

  if(clustered == TRUE){
    mod.obj$dataN <- mod.obj$data
    mod.obj$dataLN <- mod.obj$data
  }

  mod.obj$dataN <- mod.obj$dataN[order(mod.obj$dataN$dose), ] #order the data by dose
  mod.obj$dataLN <- mod.obj$dataLN[order(mod.obj$dataLN$dose), ] #order the data by dose

  dose <- sort(unique(mod.obj$dataN$dose)/max(mod.obj$dataN$dose))
  if(min(dose) == 0){
    ddd <- c(min(dose[dose > 0])/4, dose[2:length(dose)])
    lg10d <- c(log10(min(dose[dose > 0])/4), log10(dose[2:length(dose)]))
    bmdl <- log10(BMDW$BMDL/mod.obj$max.dose)
    bmdlo <- BMDW$BMDL/mod.obj$max.dose
    lg10d[1] <- ifelse((min(bmdl, na.rm=T) < lg10d[1] & min(bmdl, na.rm=T)!='-Inf' ),
                       log10(min(BMDW$BMDL/mod.obj$max.dose/2, na.rm=T)),
                       log10(min(dose[dose > 0])/4))
    ddd[1] <- ifelse(min(bmdlo, na.rm=T) < ddd[1],
                     min(BMDW$BMDL/mod.obj$max.dose/2, na.rm=T),
                     min(dose[dose > 0])/4)
  }else{
    ddd <- c(min(dose)/4, dose)
    lg10d <- c(log10(min(dose)/4), log10(dose))
    bmdl <- log10(BMDW$BMDL/mod.obj$max.dose)
    bmdlo <- BMDW$BMDL/mod.obj$max.dose
    lg10d[1] <- ifelse((min(bmdl, na.rm=T) < lg10d[1] & min(bmdl, na.rm=T)!='-Inf' ),
                       log10(min(BMDW$BMDL/mod.obj$max.dose/2, na.rm=T)),
                       log10(min(dose)/4))
    ddd[1] <- ifelse(min(bmdlo, na.rm=T) < ddd[1],
                     min(BMDW$BMDL/mod.obj$max.dose/2, na.rm=T),
                     min(dose)/4)
  }


  if(min(dose) == 0){
    if(clustered == F){
      orig.y <- log(NtoLN(mod.obj$dataN$m,mod.obj$dataN$sd))[1:length(mod.obj$dataN$dose)]
      orig.s <- log(NtoLN(mod.obj$dataN$m,mod.obj$dataN$sd))[(length(mod.obj$dataN$dose)+1):(2*length(mod.obj$dataN$dose))]
      orig_ptdataN <- data.frame(dose = (mod.obj$dataN$dose),
                                 dose2 = rep(ddd, times = table(mod.obj$dataN$dose)),
                                 lg10d = rep(lg10d, times = table(mod.obj$dataN$dose)),
                                 m = (mod.obj$dataN$m),
                                 sd = (mod.obj$dataN$sd),
                                 log10m = log10(exp(orig.y)),
                                 log10s = log10(exp(orig.s))
      )
      orig_ptdataLN <- data.frame(dose = mod.obj$dataLN$dose,
                                  dose2 = rep(ddd, times = table(mod.obj$dataLN$dose)),
                                  lg10d = rep(lg10d, times = table(mod.obj$dataLN$dose)),
                                  m = (mod.obj$dataLN$m),
                                  sd = (mod.obj$dataLN$sd),
                                  log10m = log10(exp(mod.obj$dataLN$m)),
                                  log10s = log10(exp(mod.obj$dataLN$sd))
      )

    }else if(clustered == T){
      orig_ptdata <- data.frame(dose = (mod.obj$dataN$dose),
                                dose2 = rep(ddd, times = table(mod.obj$dataN$dose)),
                                lg10d = rep(lg10d, times = table(mod.obj$dataN$dose)),
                                y = mod.obj$dataN$response,
                                yl = log(mod.obj$dataN$response),
                                log10y = log10(mod.obj$dataN$response))
    }

    if(clustered == T){

      ## overall mean
      means.all <- mod.obj$dataN %>%
        dplyr::group_by(dose) %>%
        dplyr::summarise(mresp = mean(response),
                         mlresp = mean(log(response)))
      means.litter <- mod.obj$dataN %>%
        dplyr::group_by(dose, litter) %>%
        dplyr::summarise(mresp = mean(response),
                         mlresp = mean(log(response)))

      means.all$dose = means.all$dose/max(means.all$dose)
      means.litter$dose = means.litter$dose/max(means.litter$dose)

    }

    mod.obj$dataN$lg10d <- rep(lg10d, times = table(mod.obj$dataN$dose))
    mod.obj$dataN$dose2 <- rep(ddd, times = table(mod.obj$dataN$dose))
    mod.obj$dataLN$lg10d <- rep(lg10d, times = table(mod.obj$dataLN$dose))
    mod.obj$dataLN$dose2 <- rep(ddd, times = table(mod.obj$dataLN$dose))

    if(clustered == F){
      mod.obj$dataN$lg10m <- log10(exp(orig.y))
      mod.obj$dataN$lg10s <- log10(exp(orig.s))
      mod.obj$dataLN$lg10m <- log10(exp(mod.obj$dataLN$m))
      mod.obj$dataLN$lg10s <- log10(exp(mod.obj$dataLN$s))
    }else if(clustered == T){
      # means.all$lg10m <- log10(means.all$mresp)
      # means.litter$lg10m <- log10(means.litter$mresp)
      means.all$dose[means.all$dose==0] = ddd[1]
      means.litter$dose[means.litter$dose==0] = ddd[1]
    }

  }else{

    if(clustered == F){
      orig.y <- log(NtoLN(mod.obj$dataN$m,mod.obj$dataN$sd))[1:length(mod.obj$dataN$dose)]
      orig.s <- log(NtoLN(mod.obj$dataN$m,mod.obj$dataN$sd))[(length(mod.obj$dataN$dose)+1):(2*length(mod.obj$dataN$dose))]
      orig_ptdataN <- data.frame(dose = (mod.obj$dataN$dose),
                                 dose2 = rep(ddd[2:(length(unique(mod.obj$dataN$dose))+1)], times = table(mod.obj$dataN$dose)),
                                 lg10d = rep(lg10d[2:(length(unique(mod.obj$dataN$dose))+1)], times = table(mod.obj$dataN$dose)),
                                 m = (mod.obj$dataN$m),
                                 sd = (mod.obj$dataN$sd),
                                 log10m = log10(exp(orig.y)),
                                 log10s = log10(exp(orig.s))
      )
      orig_ptdataLN <- data.frame(dose = mod.obj$dataLN$dose,
                                  dose2 = rep(ddd[2:(length(unique(mod.obj$dataLN$dose))+1)], times = table(mod.obj$dataLN$dose)),
                                  lg10d = rep(lg10d[2:(length(unique(mod.obj$dataLN$dose))+1)], times = table(mod.obj$dataLN$dose)),
                                  m = (mod.obj$dataLN$m),
                                  sd = (mod.obj$dataLN$sd),
                                  log10m = log10(exp(mod.obj$dataLN$m)),
                                  log10s = log10(exp(mod.obj$dataLN$sd))
      )

    }else if(clustered == T){
      orig_ptdata <- data.frame(dose = (mod.obj$dataN$dose),
                                dose2 = rep(ddd[2:(length(unique(mod.obj$dataN$dose))+1)], times = table(mod.obj$dataN$dose)),
                                lg10d = rep(lg10d[2:(length(unique(mod.obj$dataN$dose))+1)], times = table(mod.obj$dataN$dose)),
                                y = mod.obj$dataN$response,
                                yl = log(mod.obj$dataN$response),
                                log10y = log10(mod.obj$dataN$response))
    }

    if(clustered == T){

      ## overall mean
      means.all <- mod.obj$dataN %>%
        dplyr::group_by(dose) %>%
        dplyr::summarise(mresp = mean(response),
                         mlresp = mean(log(response)))
      means.litter <- mod.obj$dataN %>%
        dplyr::group_by(dose, litter) %>%
        dplyr::summarise(mresp = mean(response),
                         mlresp = mean(log(response)))

      means.all$dose = means.all$dose/max(means.all$dose)
      means.litter$dose = means.litter$dose/max(means.litter$dose)

    }

    mod.obj$dataN$lg10d <- rep(lg10d[2:(length(unique(mod.obj$dataN$dose))+1)], times = table(mod.obj$dataN$dose))
    mod.obj$dataN$dose2 <- rep(ddd[2:(length(unique(mod.obj$dataN$dose))+1)], times = table(mod.obj$dataN$dose))
    mod.obj$dataLN$lg10d <- rep(lg10d[2:(length(unique(mod.obj$dataLN$dose))+1)], times = table(mod.obj$dataLN$dose))
    mod.obj$dataLN$dose2 <- rep(ddd[2:(length(unique(mod.obj$dataLN$dose))+1)], times = table(mod.obj$dataLN$dose))

    if(clustered == F){
      mod.obj$dataN$lg10m <- log10(exp(orig.y))
      mod.obj$dataN$lg10s <- log10(exp(orig.s))
      mod.obj$dataLN$lg10m <- log10(exp(mod.obj$dataLN$m))
      mod.obj$dataLN$lg10s <- log10(exp(mod.obj$dataLN$s))
    }
    # else if(clustered == T){
    #   # means.all$lg10m <- log10(means.all$mresp)
    #   # means.litter$lg10m <- log10(means.litter$mresp)
    #   means.all$dose[means.all$dose==0] = ddd[1]
    #   means.litter$dose[means.litter$dose==0] = ddd[1]
    # }

  }



  dgr <- seq(min(lg10d), abs(min(lg10d)), by=0.01)
  dose2 <- 10^dgr[(dgr > (lg10d[2]-((lg10d[2]-lg10d[1])/2)))]
  if(clustered == F){
    preds <- predict.BMADR(mod.obj, dose = dose2,
                           type = type, what = "predicted",
                           model_averaged = TRUE,
                           weight_type = weight_type)
  }else if(clustered == T){
    preds <- predict.BMADR(mod.obj, dose = dose2,
                           type = type, what = "predicted",
                           model_averaged = TRUE,
                           clustered = TRUE,
                           weight_type = weight_type)
  }


  preds$predicted$dgr <- dgr[(dgr > (lg10d[2]-((lg10d[2]-lg10d[1])/2)))]
  preds$model_averaged$dgr <- dgr[(dgr > (lg10d[2]-((lg10d[2]-lg10d[1])/2)))]

  preds$predicted$lg10pred <- log10(preds$predicted$predicted)

  dgrprime <- dgr[(dgr <= (lg10d[2]-((lg10d[2]-lg10d[1])/2)))]
  preds_min <- tidyr::nest(dplyr::group_by(dplyr::filter(preds$predicted, Dose == min(Dose)), Model))
  if(clustered == F){
    as_per_model <- dplyr::select(par_med(mod.obj), Model, min_response)
  }else if(clustered == T){
    as_per_model <- dplyr::select(par_med(mod.obj, clustered = TRUE), Model, min_response)
  }
  preds_min$data <- lapply(preds_min$data, function(x) data.frame(dgrprime = dgrprime,
                                                                  Dose = 10^dgrprime,
                                                                  #Dose = x$Dose,
                                                                  predicted = x$predicted,
                                                                  lg10predicted = x$lg10pred))
  preds_min <- tidyr::unnest(preds_min, cols = c('data'))
  preds_min <- merge(preds_min, as_per_model, by = 'Model', sort = FALSE)
  preds_min$min_response[stringr::str_detect(preds_min$Model, "_LN")] <- (
    preds_min$min_response[stringr::str_detect(preds_min$Model, "_LN")])

  preds_min2 <- data.frame(lg10d = dgrprime,
                           Dose = 10^dgrprime,
                           log10MA = rep(log10(preds$model_averaged$model_averaged[
                             preds$model_averaged$Dose==min(preds$model_averaged$Dose)]),
                             length(dgrprime)),
                           MA = rep(preds$model_averaged$model_averaged[
                             preds$model_averaged$Dose==min(preds$model_averaged$Dose)],
                             length(dgrprime))
  )

  #preds$predicted$Dose <- preds$predicted$Dose*mod.obj$max.dose
  #preds$model_averaged$Dose <- preds$model_averaged$Dose*mod.obj$max.dose

  if(clustered == F){
    respBMD <- predict.BMADR(mod.obj, type = type,
                             what = "resp_at_BMD",
                             model_averaged = TRUE,
                             weight_type = weight_type)
  }else if(clustered == T){
    respBMD <- predict.BMADR(mod.obj, type = type,
                             what = "resp_at_BMD",
                             model_averaged = TRUE,
                             clustered = TRUE,
                             weight_type = weight_type)
  }

  respBMD$resp_at_BMD$lg10bmd <- log10(respBMD$resp_at_BMD$BMD)
  #respBMD$resp_at_BMD[,"BMD"] = respBMD$resp_at_BMD[,"BMD"]*mod.obj$max.dose
  #print(respBMD)
  #respBMD <- refactor(respBMD)

  BMDMixture <- BMDmixture_extract(mod.obj, weight_type, conv=FALSE) # BMD values
  BMDMixture$BMDMixture2 <- log10(BMDMixture$BMDMixture/mod.obj$max.dose)

  gghst2 <- hist(BMDMixture$BMDMixture, breaks = sqrt(nrow(BMDMixture)), plot = FALSE) #hist on original scale

  lais <- log10(gghst2$breaks) # log10 of the breakpoints
  dlais <- abs(diff(lais)) # width of the interval on log scale
  glais <- (gghst2$counts/dlais) # divide the original counts by the width on the log scale
  glais2 <- glais/sum(glais) * nrow(BMDMixture) # normalise the new frequencie

  BMDMixture2 <- data.frame(Model = unique(BMDMixture$Model),#rep(unique(BMDMixture$Model), length(gghst2$counts)),
                            # Dose = gghst2$mids, #midpoints
                            Dose = 10**lais[1:length(gghst2$counts)], #midpoints
                            y = gghst2$counts, #frequencies
                            y2 = glais2
  )

  #BMDMixtureD <- density(BMDMixture$BMDMixture2, # density of BMD mixture on log10 scale
  #                       from = min(BMDMixture$BMDMixture2),
  #                       to = max(BMDMixture$BMDMixture2))

  #BMDMixtureD2 <- density(BMDMixture$BMDMixture, # density of BMD mixture on original dose scale
  #                        from = min(BMDMixture$BMDMixture),
  #                        to = max(BMDMixture$BMDMixture))

  #BMDMixture2 <- data.frame(Model = unique(BMDMixture$Model), # y = Density
  #                          Dose10 = BMDMixtureD$x,
  #                          Dose102 = 10^BMDMixtureD$x,
  #                          y10 = BMDMixtureD$y,
  #                          Dose = BMDMixtureD2$x,
  #                          y = BMDMixtureD2$y)
  #BMDMixture <- refactor(BMDMixture)

  #BMDs and Weights
  respBMDBMDW <- merge(respBMD$resp_at_BMD, BMDW, by = c("Model", "BMD"), sort = FALSE)
  respBMDBMDW <- refactor(respBMDBMDW, type)
  respBMDBMDW$Distribution <- ifelse(stringr::str_detect(respBMDBMDW$Model, "(LN)"), "LN", "N")
  pd <- position_dodge(0.5)
  # dist_fills <- c("N" = "#EFC000FF", "LN" = "#0073C2FF")
  # dist_names <- c("N" = "Normal", "LN" = "LogNormal")
  dist_fills <- c("#EFC000FF", "#0073C2FF")
  names(dist_fills) <- c("N","LN")
  dist_names <- c("Normal","LogNormal")
  names(dist_names) <- c("N","LN")

  ## If BMDU > maxD^2 set to maxD^2
  # if(mod.obj$max.dose > 1){
  #   respBMDBMDW$BMDU <- ifelse(respBMDBMDW$BMDU > mod.obj$max.dose^2, mod.obj$max.dose^2, respBMDBMDW$BMDU)
  # }
  if(mod.obj$max.dose > 1){
    respBMDBMDW$BMDU <- ifelse(((respBMDBMDW$BMDU > mod.obj$max.dose^2) &
                                  ((respBMDBMDW$BMDU/respBMDBMDW$BMD)/(respBMDBMDW$BMD/respBMDBMDW$BMDL)) > 10) |
                                 respBMDBMDW$BMDU > 10*mod.obj$max.dose^2,
                               # mod.obj$max.dose^2,
                               2*(respBMDBMDW$BMD - respBMDBMDW$BMDL),
                               respBMDBMDW$BMDU)
  }

  if(clustered == F){
    mod.obj$dataN$response = mod.obj$dataN$m
    mod.obj$dataLN$response = mod.obj$dataLN$m
  }


  #BMD plots
  pBMDs <- ggplot(data = respBMDBMDW, aes(x = BMD, y = Model, group = Model)) +
    geom_errorbarh(aes(xmin = BMDL, xmax = BMDU, y = Model),
                   linetype = "solid", show.legend = FALSE,
                   size = 3, height = 1.2, color = brewer.pal(8, "Set2")[8]) +
    geom_point(aes(x = BMD, y = Model, fill = Distribution),
               size = 7, color = 1, shape = 21) +
    theme_minimal() +
    labs(x = expression(BMD), y = "", fill = "Distribution") +
    scale_fill_manual(values = dist_fills,
                      labels = dist_names) +
    labs(x="BMD on original scale") +
    #labs(x = "BMD", y = "Models", title = "BMDs per Model") +
    # scale_x_continuous(breaks = seq(min(log10(respBMDBMDW$BMDL)),
    #                                 max(log10(respBMDBMDW$BMDU)), by = 0.08)
    #                    ) +
    theme(strip.text = element_text(size = 15, face = "bold"),
          axis.text = element_text(size = 10, face = "bold"),
          axis.title = element_text(size = 15, face = "bold"),
          legend.text = element_text(size = 10, face = "bold"),
          legend.title = element_text(size = 15, face = "bold"),
          panel.spacing = unit(5, "lines"),
          legend.position = "top",
          legend.direction = "horizontal",
          title = element_text(size = 15, face = "bold")) +
    geom_errorbarh(data = BMDBMA[BMDBMA$Type == weight_type,],
                   aes(xmin = BMDL, xmax = BMDU, y = Model),
                   linetype = "solid", show.legend = FALSE,
                   size = 3, height = 1.2, color = brewer.pal(9, "Set1")[3]) +
    geom_point(data = BMDBMA[BMDBMA$Type == weight_type,],
               aes(x = BMD, y = Model),
               size = 7, color = 1, shape = 21,
               fill = brewer.pal(9, "Set1")[1]) +
    #scale_x_log10() +
    scale_y_discrete(limits = rev(c("Model Averaged", "Exp(N)", "InvExp(N)", "Hill(N)",
                                    "LogNormal(N)", "Gamma(N)",  "QuadExp(N)", "Probit(N)",
                                    "Logistic(N)", "Exp(LN)", "InvExp(LN)", "Hill4(LN)",
                                    "LogNormal(LN)", "Gamma(LN)", "QuadExp(LN)",
                                    "Probit(LN)", "Logistic(LN)")),
                     labels = rev(c("Model Averaged", "Exp(N)", "InvExp(N)", "Hill(N)",
                                    "LogNormal(N)", "Gamma(N)",  "QuadExp(N)", "Probit(N)",
                                    "Logistic(N)", "Exp(LN)", "InvExp(LN)", "Hill4(LN)",
                                    "LogNormal(LN)", "Gamma(LN)", "QuadExp(LN)",
                                    "Probit(LN)", "Logistic(LN)")))
  #Weights plot
  needed_data <- tidyr::separate(BMDW,
                                 col = "Model", sep = "_",
                                 into = c("Model", "Distribution"))

  if(weight_type == "LP") {

    pWeights <- ggpubr::ggdotchart(
      needed_data, x = "Model", y = "LP_Weights",
      group = "Distribution", color = "Distribution",
      add = "segment", position = position_dodge(0.3),
      sorting = "descending", size = 3, dot.size = 7) +
      labs(x = "", y = "Weight", color = "Distribution",
           title = "") +
      theme_minimal() +
      theme(strip.text = element_text(size = 15, face = "bold"),
            axis.text = element_text(size = 10, face = "bold"),
            axis.title = element_text(size = 15, face = "bold"),
            legend.text = element_text(size = 10, face = "bold"),
            legend.title = element_text(size = 15, face = "bold"),
            panel.spacing = unit(5, "lines"),
            legend.position = "top",
            legend.direction = "horizontal",
            title = element_text(size = 15, face = "bold")) +
      scale_x_discrete(
        limits = c("E4", "IE4", "H4", "LN4", "G4", "QE4", "P4", "L4"),
        labels = c("Exp", "InvExp", "Hill", "LogNormal", "Gamma",
                   "QuadExp", "Probit", "Logistic")
      ) +
      scale_color_manual(values = dist_fills,
                         labels = dist_names)
    needed_data2 <- dplyr::arrange(needed_data, desc(LP_Weights))
    needed_data2$CWeights <- cumsum(needed_data2$LP_Weights)
    needed_data2 <- dplyr::rename(needed_data2, Weights = LP_Weights)
  } else {

    pWeights <- ggpubr::ggdotchart(
      needed_data, x = "Model", y = "BS_Weights",
      group = "Distribution", color = "Distribution",
      add = "segment", position = position_dodge(0.3),
      sorting = "descending", size = 3, dot.size = 7) +
      labs(x = "", y = "Weight", color = "Distribution",
           title = "") +
      theme_minimal() +
      theme(strip.text = element_text(size = 15, face = "bold"),
            axis.text = element_text(size = 10, face = "bold"),
            axis.title = element_text(size = 15, face = "bold"),
            legend.text = element_text(size = 10, face = "bold"),
            legend.title = element_text(size = 15, face = "bold"),
            panel.spacing = unit(5, "lines"),
            legend.position = "top",
            legend.direction = "horizontal",
            title = element_text(size = 15, face = "bold")) +
      scale_x_discrete(
        limits = c("E4", "IE4", "H4", "LN4", "G4", "QE4", "P4", "L4"),
        labels = c("Exp", "InvExp", "Hill", "LogNormal", "Gamma",
                   "QuadExp", "Probit", "Logistic")
      ) +
      scale_color_manual(values = dist_fills,
                         labels = dist_names)

    needed_data2 <- dplyr::arrange(needed_data, desc(BS_Weights))
    needed_data2$CWeights <- cumsum(needed_data2$BS_Weights)
    needed_data2 <- dplyr::rename(needed_data2, Weights = BS_Weights)

  }

  ##### Mixture Distribution
  # needed_data3 <- needed_data2[which(needed_data2$CWeights <= 0.85), ]
  # needed_data3$mod_dist <- paste(needed_data3$Model, needed_data3$Distribution, sep = '_')
  #
  # names(mod.obj$parsN) = c("E4_N","IE4_N","H4_N","LN4_N","G4_N","QE4_N","P4_N","L4_N")
  # names(mod.obj$parsLN) = c("E4_LN","IE4_LN","H4_LN","LN4_LN","G4_LN","QE4_LN","P4_LN","L4_LN")
  #
  # pars = c(mod.obj$parsN, mod.obj$parsLN)
  #
  # ddens <- lapply(needed_data3$mod_dist, function(i){
  #   xx <- pars[[i]][,c("ModelName","BMD")]
  #   yy <- density(xx$BMD)
  #   data.frame(x = yy$x, y = yy$y,
  #              yw = yy$y*needed_data3$Weights[needed_data3$mod_dist==i],
  #              ModelName = i)
  # })
  # df.mod <- do.call(rbind.data.frame, ddens)




  #Prediction plot
  mods_fills <- c("Model-averaged BMD" = "coral")
  #print(BMDBMA)
  if(is.BMADR2(mod.obj)[3]==2) {

    BMDBMA$Response <- c(sum(respBMDBMDW$resp_at_BMD*respBMDBMDW$BS_Weights),
                         sum(respBMDBMDW$resp_at_BMD*respBMDBMDW$LP_Weights))

  } else if(is.BMADR2(mod.obj)[2]==2) {

    BMDBMA$Response <- sum(respBMDBMDW$resp_at_BMD*respBMDBMDW$LP_Weights)

  } else {

    pts <- ggpubr::ggarrange(pBMDs, pWeights, nrow = 1, ncol = 2)
    return(pts)
    stop("cannot compute model averaged response at BMD. Please check the inputs.")
  }

  #BMDBMA$Response <- c(sum(respBMDBMDW$resp_at_BMD*respBMDBMDW$BS_Weights),
  #sum(respBMDBMDW$resp_at_BMD*respBMDBMDW$LP_Weights))

  preds2 <- tidyr::separate(preds$predicted,
                            col = "Model", sep = "_",
                            into = c("Model", "Distribution"))
  preds_min <- tidyr::separate(preds_min,
                               col = "Model", sep = "_",
                               into = c("Model", "Distribution"))
  lty <- c("LN" = "twodash",
           "N" = "solid")
  md_cls <- RColorBrewer::brewer.pal(8, "Dark2")
  names(md_cls) <-  c("E4", "IE4", "H4",
                      "LN4", "G4", "QE4",
                      "P4", "L4")


  if(min(mod.obj$dataN$dose) == 0){
    plot.labs = dose*mod.obj$max.dose
  }else{
    plot.labs = c(0, dose*mod.obj$max.dose)
  }

  if(clustered == F){

    ymin = 10^min(mod.obj$dataN$lg10m - 2*mod.obj$dataN$lg10s, na.rm=T)
    ymax = 10^max(mod.obj$dataN$lg10m + 2*mod.obj$dataN$lg10s, na.rm=T)

    # plot for Normal distribution
    pplotN <- ggplot(data = preds2[preds2$Distribution=="N",],
                     aes(x = Dose*mod.obj$max.dose, y = predicted, group = Model,
                         color = Model)) +
      geom_line(alpha = 0.6,
                size = 1,
                show.legend = TRUE, linetype = 1) +
      labs(color = "Model", x = expression(dose),
           y = expression(response), title = "Normal distribution",
           caption = "data and vertical bars based on arithmetic sample means and standard deviations") +

      geom_segment(data = preds_min[preds_min$Distribution=="N",],
                   mapping = aes(x = Dose[1]*mod.obj$max.dose, y = min_response,
                                 xend = max(Dose*mod.obj$max.dose),#max(dgr[(dgr <= (lg10d[2]-((lg10d[2]-lg10d[1])/2)))]),
                                 yend = predicted,
                                 group = interaction(Model, Distribution),
                                 color = Model),
                   linetype = "dotted", alpha = 0.6,
                   size = 0.8, inherit.aes = FALSE, show.legend = FALSE) +
      #ylim(ymin, ymax) +
      geom_errorbar(data = orig_ptdataN,
                    mapping = aes(x = dose2*mod.obj$max.dose, ymin = m-sd,
                                  ymax = m+sd),
                    width = NA, position = pd, size = 1,
                    inherit.aes = FALSE) +
      geom_point(data = orig_ptdataN, mapping = aes(x = dose2*mod.obj$max.dose, y = m) ,
                 size = 2, color = 1, shape = 21,
                 fill = brewer.pal(9, "Set1")[2],
                 inherit.aes = FALSE) +
      # geom_jitter(data = orig_ptdata, mapping = aes(x = dose2*mod.obj$max.dose, y = m) ,
      #              size = 2, color = 1, shape = 21,
      #              fill = brewer.pal(9, "Set1")[2], position = position_jitter(h = 0.01, width = 0),
      #              inherit.aes = FALSE) +
      geom_errorbarh(data = dplyr::filter(BMDBMA, Type == weight_type),
                     aes(xmin = BMDL, xmax = BMDU,
                         y = Response,
                         group = Model),
                     linetype = "solid", show.legend = FALSE,
                     size = 2, height = 0.01*log10(mean(mod.obj$dataN$m)), inherit.aes = FALSE,
                     color = brewer.pal(9, "Set1")[3]) +
      geom_point(data = dplyr::filter(BMDBMA, Type == weight_type),
                 aes(x = BMD, y = Response, group = Model),
                 size = 5, shape = 19,
                 color = brewer.pal(9, "Set1")[1],
                 show.legend = FALSE,
                 inherit.aes = FALSE) +
      coord_cartesian(xlim = c(min(preds_min$Dose*mod.obj$max.dose),
                               2*mod.obj$max.dose) ) +
      # scale_x_continuous(trans = 'log10', labels = scales::comma,
      #                    breaks = orig_ptdata$dose2[2:length(orig_ptdata$dose2)]*mod.obj$max.dose) +
      scale_x_continuous(trans = 'log10', labels = plot.labs,
                         breaks = ddd*mod.obj$max.dose) +
      # scale_linetype_manual(values = lty,
      # labels = c('LogNormal', 'Normal')) +
      scale_color_manual(values = md_cls,
                         labels = c("Exp", "InvExp", "Hill", "LogNormal", "Gamma",
                                    "QuadExp", "Probit", "Logistic")) +
      theme_minimal() +
      theme(strip.text = element_text(size = 15, face = "bold"),
            axis.text = element_text(size = 10, face = "bold"),
            axis.title = element_text(size = 15, face = "bold"),
            legend.text = element_text(size = 8),
            legend.title = element_text(size = 12, face = "bold"),
            panel.spacing = unit(5, "lines"),
            legend.position = "top",
            legend.direction = "horizontal",
            title = element_text(size = 15, face = "bold"))

    if(log == F){
      pplotN <- pplotN + scale_y_continuous(trans = 'identity', labels = scales::comma)
    }else{
      pplotN <- pplotN + scale_y_continuous(trans = 'log10', labels = scales::comma)
    }

    ymin = 10^min(mod.obj$dataLN$lg10m - 2*mod.obj$dataLN$lg10s, na.rm=T)
    ymax = 10^max(mod.obj$dataLN$lg10m + 2*mod.obj$dataLN$lg10s, na.rm=T)

    # plot for LogNormal distribution
    pplotLN <- ggplot(data = preds2[preds2$Distribution=="LN",],
                      aes(x = Dose*mod.obj$max.dose, y = predicted, group = Model,
                          color = Model)) +
      geom_line(alpha = 0.6,
                size = 1,
                show.legend = TRUE, linetype = 2) +
      labs(color = "Model", title = "LogNormal distribution", x = expression(dose),
           y = expression(response),
           caption = "data and vertical bars based on geometric sample means and standard deviations") +
      geom_segment(data = preds_min[preds_min$Distribution=="LN",],
                   mapping = aes(x = Dose[1]*mod.obj$max.dose, y = min_response,
                                 xend = max(Dose*mod.obj$max.dose), #max(dgr[(dgr <= (lg10d[2]-((lg10d[2]-lg10d[1])/2)))]),
                                 yend = predicted,
                                 group = interaction(Model, Distribution),
                                 color = Model),
                   linetype = "dotted", alpha = 0.6,
                   size = 0.8, inherit.aes = FALSE, show.legend = FALSE) +
      #ylim(ymin, ymax) +
      geom_errorbar(data = orig_ptdataLN,
                    mapping = aes(x = dose2*mod.obj$max.dose, ymin = 10^(log10m-log10s),
                                  ymax = 10^(log10m+log10s)),
                    width = NA, position = pd, size = 1,
                    inherit.aes = FALSE) +
      geom_point(data = orig_ptdataLN, mapping = aes(x = dose2*mod.obj$max.dose, y = 10^log10m) ,
                 size = 2, color = 1, shape = 21,
                 fill = brewer.pal(9, "Set1")[2],
                 inherit.aes = FALSE) +
      geom_errorbarh(data = dplyr::filter(BMDBMA, Type == weight_type),
                     aes(xmin = BMDL, xmax = BMDU,
                         y = Response,
                         group = Model),
                     linetype = "solid", show.legend = FALSE,
                     size = 2, height = 0.01*log10(mean(mod.obj$dataN$m)), inherit.aes = FALSE,
                     color = brewer.pal(9, "Set1")[3]) +
      geom_point(data = dplyr::filter(BMDBMA, Type == weight_type),
                 aes(x = BMD, y = Response, group = Model),
                 size = 5, shape = 19,
                 color = brewer.pal(9, "Set1")[1],
                 show.legend = FALSE,
                 inherit.aes = FALSE) +
      # scale_linetype_manual(values = lty,
      # labels = c('LogNormal')) +
      coord_cartesian(xlim = c(min(preds_min$Dose*mod.obj$max.dose),
                               2*mod.obj$max.dose) ) +
      # scale_x_continuous(trans = 'log10', labels = scales::comma,
      #                    breaks = orig_ptdata$dose2[2:length(orig_ptdata$dose2)]*mod.obj$max.dose) +
      scale_x_continuous(trans = 'log10', labels = plot.labs,
                         breaks = ddd*mod.obj$max.dose) +
      scale_y_continuous(trans = 'log10', labels = scales::comma) +
      scale_colour_manual(values = md_cls,
                         labels = c("Exp", "InvExp", "Hill", "LogNormal", "Gamma",
                                    "QuadExp", "Probit", "Logistic")) +
      theme_minimal() +
      theme(strip.text = element_text(size = 15, face = "bold"),
            axis.text = element_text(size = 10, face = "bold"),
            axis.title = element_text(size = 15, face = "bold"),
            legend.text = element_text(size = 8),
            legend.title = element_text(size = 12, face = "bold"),
            panel.spacing = unit(5, "lines"),
            legend.position = "top",
            legend.direction = "horizontal",
            title = element_text(size = 15, face = "bold"))

    ## Plot for both distributions
    pplot <- ggplot(data = preds2,
                    aes(x = Dose*mod.obj$max.dose, y = predicted, group = interaction(Model, Distribution),
                        color = Model, linetype = Distribution)) +
      geom_line(alpha = 0.6,
                size = 1,
                show.legend = TRUE) +
      labs(color = "Model", linetype = "Distribution", x = expression(dose),
           y = expression(response), title = "") +

      geom_segment(data = preds_min, mapping = aes(x = Dose[1]*mod.obj$max.dose, y = min_response,
                                                   xend = max(Dose*mod.obj$max.dose), #max(dgr[(dgr <= (lg10d[2]-((lg10d[2]-lg10d[1])/2)))]),
                                                   yend = predicted,
                                                   group = interaction(Model, Distribution),
                                                   color = Model),
                   linetype = "dotted", alpha = 0.6,
                   size = 0.8, inherit.aes = FALSE, show.legend = FALSE) +

      # geom_line(data = preds_min, mapping = aes(x = dgrprime, y = log10(min_response),
      #                                           group = interaction(Model, Distribution),
      #                                           color = Model),
      #           linetype = "dotted", alpha = 0.6,
      #           size = 0.8, inherit.aes = FALSE, show.legend = FALSE) +
      #ylim(ymin, ymax) +
      geom_errorbarh(data = dplyr::filter(BMDBMA, Type == weight_type),
                     aes(xmin = BMDL, xmax = BMDU,
                         y = Response,
                         group = Model),
                     linetype = "solid", show.legend = FALSE,
                     size = 2, height = 0.01*log10(mean(mod.obj$dataN$m)), inherit.aes = FALSE,
                     color = brewer.pal(9, "Set1")[3]) +
      geom_point(data = dplyr::filter(BMDBMA, Type == weight_type),
                 aes(x = BMD, y = Response, group = Model),
                 size = 5, shape = 19,
                 color = brewer.pal(9, "Set1")[1],
                 show.legend = FALSE,
                 inherit.aes = FALSE) +
      scale_linetype_manual(values = lty,
                            labels = c('LogNormal', 'Normal')) +
      scale_color_manual(values = md_cls,
                         labels = c("Exp", "InvExp", "Hill", "LogNormal", "Gamma",
                                    "QuadExp", "Probit", "Logistic")) +
      theme_minimal() +
      coord_cartesian(xlim = c(min(preds_min$Dose*mod.obj$max.dose),
                               2*mod.obj$max.dose) ) +
      # scale_x_continuous(trans = 'log10', labels = scales::comma,
      #                    breaks = orig_ptdata$dose2[2:length(orig_ptdata$dose2)]*mod.obj$max.dose) +
      scale_x_continuous(trans = 'log10', labels = plot.labs,
                         breaks = ddd*mod.obj$max.dose) +
      scale_y_continuous(trans = 'log10', labels = scales::comma) +
      theme(strip.text = element_text(size = 15, face = "bold"),
            axis.text = element_text(size = 10, face = "bold"),
            axis.title = element_text(size = 15, face = "bold"),
            legend.text = element_text(size = 8),
            legend.title = element_text(size = 12, face = "bold"),
            panel.spacing = unit(5, "lines"),
            legend.position = "top",
            legend.direction = "horizontal",
            title = element_text(size = 15, face = "bold"))

  }else if(clustered == T){

    # plot for Normal distribution
    pplotN <- ggplot(data = preds2[preds2$Distribution=="N",],
                     aes(x = Dose*mod.obj$max.dose, y = predicted, group = Model,
                         color = Model)) +
      geom_line(alpha = 0.6,
                size = 1,
                show.legend = TRUE, linetype = 1) +
      labs(color = "Model", x = expression(dose),
           y = expression(response), title = "Normal distribution",
           caption = "diamonds represent the arithmetic sample mean") +
      geom_segment(data = preds_min[preds_min$Distribution=="N",],
                   mapping = aes(x = Dose[1]*mod.obj$max.dose, y = min_response,
                                 xend = max(Dose*mod.obj$max.dose),#max(dgr[(dgr <= (lg10d[2]-((lg10d[2]-lg10d[1])/2)))]),
                                 yend = predicted,
                                 group = interaction(Model, Distribution),
                                 color = Model),
                   linetype = "dotted", alpha = 0.6,
                   size = 0.8, inherit.aes = FALSE, show.legend = FALSE) +
      scale_x_continuous(trans = 'log10', labels = plot.labs,
                         breaks = ddd*mod.obj$max.dose) +

      geom_jitter(data = orig_ptdata, mapping = aes(x = dose2*mod.obj$max.dose, y = y),
                  size = 1, color = 3, shape = 20,
                  # fill = brewer.pal(9, "Set1")[2],
                  position = position_jitter(h = 0.01, width = 0),
                  inherit.aes = FALSE) +

      geom_point(data = means.litter, mapping = aes(x = dose*mod.obj$max.dose, y = mresp),
                 size = 2, color = 1, shape = 21,
                 fill = 1,
                 inherit.aes = FALSE) +
      geom_point(data = means.all, mapping = aes(x = dose*mod.obj$max.dose, y = mresp),
                 size = 3, color = 2, shape = 23,
                 fill = 2,
                 inherit.aes = FALSE) +
      geom_errorbarh(data = dplyr::filter(BMDBMA, Type == weight_type),
                     aes(xmin = BMDL, xmax = BMDU,
                         y = Response,
                         group = Model),
                     linetype = "solid", show.legend = FALSE,
                     size = 2, height = 0.01*log10(mean(mod.obj$dataN$response)), inherit.aes = FALSE,
                     color = brewer.pal(9, "Set1")[3]) +
      geom_point(data = dplyr::filter(BMDBMA, Type == weight_type),
                 aes(x = BMD, y = Response, group = Model),
                 size = 5, shape = 19,
                 color = brewer.pal(9, "Set1")[1],
                 show.legend = FALSE,
                 inherit.aes = FALSE) +
      coord_cartesian(xlim = c(min(preds_min$Dose*mod.obj$max.dose),
                               2*mod.obj$max.dose) ) +
      scale_color_manual(values = md_cls,
                         labels = c("Exp", "InvExp", "Hill", "LogNormal", "Gamma",
                                    "QuadExp", "Probit", "Logistic")) +
      theme_minimal() +
      theme(strip.text = element_text(size = 15, face = "bold"),
            axis.text = element_text(size = 10, face = "bold"),
            axis.title = element_text(size = 15, face = "bold"),
            legend.text = element_text(size = 8),
            legend.title = element_text(size = 12, face = "bold"),
            panel.spacing = unit(5, "lines"),
            legend.position = "top",
            legend.direction = "horizontal",
            title = element_text(size = 15, face = "bold"))

    if(log == F){
      pplotN <- pplotN + scale_y_continuous(trans = 'identity', labels = scales::comma)
    }else{
      pplotN <- pplotN + scale_y_continuous(trans = 'log10', labels = scales::comma)
    }

    # plot for LogNormal distribution
    pplotLN <- ggplot(data = preds2[preds2$Distribution=="LN",],
                      aes(x = Dose*mod.obj$max.dose, y = predicted, group = Model,
                          color = Model)) +
      geom_line(alpha = 0.6,
                size = 1,
                show.legend = TRUE, linetype = 2) +
      labs(color = "Model", title = "LogNormal distribution", x = expression(dose),
           y = expression(response),
           caption = "diamonds represent the geometric sample mean") +

      geom_segment(data = preds_min[preds_min$Distribution=="LN",],
                   mapping = aes(x = Dose[1]*mod.obj$max.dose, y = min_response,
                                 xend = max(Dose*mod.obj$max.dose), #max(dgr[(dgr <= (lg10d[2]-((lg10d[2]-lg10d[1])/2)))]),
                                 yend = predicted,
                                 group = interaction(Model, Distribution),
                                 color = Model),
                   linetype = "dotted", alpha = 0.6,
                   size = 0.8, inherit.aes = FALSE, show.legend = FALSE) +

      geom_jitter(data = orig_ptdata, mapping = aes(x = dose2*mod.obj$max.dose, y = 10^log10(exp(yl))),
                  size = 1, color = 3, shape = 20,
                  # fill = brewer.pal(9, "Set1")[2],
                  position = position_jitter(h = 0.01, width = 0),
                  inherit.aes = FALSE) +

      geom_point(data = means.litter, mapping = aes(x = dose*mod.obj$max.dose, y = 10^log10(exp(mlresp))),
                 size = 2, color = 1, shape = 21,
                 fill = 1,
                 inherit.aes = FALSE) +
      geom_point(data = means.all, mapping = aes(x = dose*mod.obj$max.dose, y = 10^log10(exp(mlresp))),
                 size = 3, color = 2, shape = 23,
                 fill = 2,
                 inherit.aes = FALSE) +

      # geom_point(data = orig_ptdata, mapping = aes(x = dose2*mod.obj$max.dose, y = 10^log10m) ,
      #            size = 2, color = 1, shape = 21,
      #            fill = brewer.pal(9, "Set1")[2],
      #            inherit.aes = FALSE) +
      geom_errorbarh(data = dplyr::filter(BMDBMA, Type == weight_type),
                     aes(xmin = BMDL, xmax = BMDU,
                         y = Response,
                         group = Model),
                     linetype = "solid", show.legend = FALSE,
                     size = 2, height = 0.01*log10(mean(mod.obj$dataN$response)), inherit.aes = FALSE,
                     color = brewer.pal(9, "Set1")[3]) +
      geom_point(data = dplyr::filter(BMDBMA, Type == weight_type),
                 aes(x = BMD, y = Response, group = Model),
                 size = 5, shape = 19,
                 color = brewer.pal(9, "Set1")[1],
                 show.legend = FALSE,
                 inherit.aes = FALSE) +
      # scale_linetype_manual(values = lty,
      # labels = c('LogNormal')) +
      coord_cartesian(xlim = c(min(preds_min$Dose*mod.obj$max.dose),
                               2*mod.obj$max.dose) ) +
      # scale_x_continuous(trans = 'log10', labels = scales::comma,
      #                    breaks = orig_ptdata$dose2[2:length(orig_ptdata$dose2)]*mod.obj$max.dose) +
      scale_x_continuous(trans = 'log10', labels = plot.labs,
                         breaks = ddd*mod.obj$max.dose) +
      scale_y_continuous(trans = 'log10', labels = scales::comma) +
      scale_color_manual(values = md_cls,
                         labels = c("Exp", "InvExp", "Hill", "LogNormal", "Gamma",
                                    "QuadExp", "Probit", "Logistic")) +
      theme_minimal() +
      theme(strip.text = element_text(size = 15, face = "bold"),
            axis.text = element_text(size = 10, face = "bold"),
            axis.title = element_text(size = 15, face = "bold"),
            legend.text = element_text(size = 8),
            legend.title = element_text(size = 12, face = "bold"),
            panel.spacing = unit(5, "lines"),
            legend.position = "top",
            legend.direction = "horizontal",
            title = element_text(size = 15, face = "bold"))

    ## Plot for both distributions
    pplot <- ggplot(data = preds2,
                    aes(x = Dose*mod.obj$max.dose, y = predicted, group = interaction(Model, Distribution),
                        color = Model, linetype = Distribution)) +
      geom_line(alpha = 0.6,
                size = 1,
                show.legend = TRUE) +
      labs(color = "Model", linetype = "Distribution", x = expression(dose),
           y = expression(response), title = "") +

      geom_segment(data = preds_min, mapping = aes(x = Dose[1]*mod.obj$max.dose, y = min_response,
                                                   xend = max(Dose*mod.obj$max.dose), #max(dgr[(dgr <= (lg10d[2]-((lg10d[2]-lg10d[1])/2)))]),
                                                   yend = predicted,
                                                   group = interaction(Model, Distribution),
                                                   color = Model),
                   linetype = "dotted", alpha = 0.6,
                   size = 0.8, inherit.aes = FALSE, show.legend = FALSE) +

      # geom_line(data = preds_min, mapping = aes(x = dgrprime, y = log10(min_response),
      #                                           group = interaction(Model, Distribution),
      #                                           color = Model),
      #           linetype = "dotted", alpha = 0.6,
      #           size = 0.8, inherit.aes = FALSE, show.legend = FALSE) +
      #ylim(ymin, ymax) +
      geom_errorbarh(data = dplyr::filter(BMDBMA, Type == weight_type),
                     aes(xmin = BMDL, xmax = BMDU,
                         y = Response,
                         group = Model),
                     linetype = "solid", show.legend = FALSE,
                     size = 2, height = 0.01*log10(mean(mod.obj$dataN$response)), inherit.aes = FALSE,
                     color = brewer.pal(9, "Set1")[3]) +
      geom_point(data = dplyr::filter(BMDBMA, Type == weight_type),
                 aes(x = BMD, y = Response, group = Model),
                 size = 5, shape = 19,
                 color = brewer.pal(9, "Set1")[1],
                 show.legend = FALSE,
                 inherit.aes = FALSE) +
      scale_linetype_manual(values = lty,
                            labels = c('LogNormal', 'Normal')) +
      scale_color_manual(values = md_cls,
                         labels = c("Exp", "InvExp", "Hill", "LogNormal", "Gamma",
                                    "QuadExp", "Probit", "Logistic")) +
      theme_minimal() +
      coord_cartesian(xlim = c(min(preds_min$Dose*mod.obj$max.dose),
                               2*mod.obj$max.dose) ) +
      # scale_x_continuous(trans = 'log10', labels = scales::comma,
      #                    breaks = orig_ptdata$dose2[2:length(orig_ptdata$dose2)]*mod.obj$max.dose) +
      scale_x_continuous(trans = 'log10', labels = plot.labs,
                         breaks = ddd*mod.obj$max.dose) +
      scale_y_continuous(trans = 'log10', labels = scales::comma) +
      theme(strip.text = element_text(size = 15, face = "bold"),
            axis.text = element_text(size = 10, face = "bold"),
            axis.title = element_text(size = 15, face = "bold"),
            legend.text = element_text(size = 8),
            legend.title = element_text(size = 12, face = "bold"),
            panel.spacing = unit(5, "lines"),
            legend.position = "top",
            legend.direction = "horizontal",
            title = element_text(size = 15, face = "bold"))

  }

  cmax <- max(preds$model_averaged$model_averaged)
  cmin <- min(preds$model_averaged$model_averaged)
  # BMDMixture$yres <- (BMDMixture$y/max(BMDMixture$y)) * (cmax * 1.2) # Density rescaled

  ylim.prim <- c(cmin, cmax)
  # ylim.sec <- c(min(BMDMixture2$y2), max(BMDMixture2$y2))
  ylim.sec <- c(min(BMDMixture2$y), max(BMDMixture2$y))
  b <- diff(ylim.prim)/diff(ylim.sec)
  a <- ylim.prim[1] - b*ylim.sec[1]

  # BMDMixture2$yres = a + BMDMixture2$y2*b
  BMDMixture2$yres = a + BMDMixture2$y*b

  dplot <- ggplot(data = preds$model_averaged, aes(x = Dose*mod.obj$max.dose, y = model_averaged,
                                                   group = 1)) +
    geom_line(show.legend = FALSE, linetype = "dashed", size = 3) +
    #geom_col(data = BMDMixture2, aes(x = Dose2, y = yres, fill = Model), alpha = 0.6,
    #         inherit.aes = FALSE, color = NA) +
    #geom_area(data = BMDMixture2,
    #           aes(x = Dose, y = yres, group = Model, fill = Model),
    #           color = NA, inherit.aes = FALSE, alpha = 0.5) +
    #geom_histogram(data = BMDMixture,
    #               aes(x = BMDMixture, fill = Model), #y = ..density..,
    #               color = NA, alpha = 0.5, inherit.aes = FALSE, bins = sqrt(nrow(BMDMixture))) +
    geom_ribbon(data = BMDMixture2, aes(x = Dose, y = yres, ymin = min(yres), ymax = yres, fill = Model),
                color=NA, alpha = 0.5) +
    # geom_line(data = preds_min, mapping = aes(x = dgrprime, y = lg10predicted,
    #                                           group = interaction(Model, Distribution)),
    #           linetype = "dashed", alpha = 0.6,
    #           size = 0.8, inherit.aes = FALSE, show.legend = FALSE, color = 1) +

    # geom_line(data = preds_min2, mapping = aes(x = Dose, y = log10MA),
    #           linetype = "dotted", alpha = 1,
    #           size = 3, inherit.aes = FALSE, show.legend = FALSE) +

    geom_segment(data = preds_min2, mapping = aes(x = Dose[1]*mod.obj$max.dose, y = MA,
                                                  xend = max(Dose*mod.obj$max.dose),
                                                  #max(dgr[(dgr <= (lg10d[2]-((lg10d[2]-lg10d[1])/2)))]),
                                                  yend = preds$model_averaged$model_averaged[1],
                                                  # group = interaction(Model, Distribution),
                                                  # color = Model
    ),
    linetype = "dotted",
    size = 3, inherit.aes = FALSE, show.legend = FALSE) +
    geom_errorbarh(data = dplyr::filter(BMDBMA, Type == weight_type),
                   aes(xmin = BMDL, xmax = BMDU,
                       y = respBMD$model_averaged$model_averaged_response,#min(preds$model_averaged$model_averaged),
                       group = Model),
                   linetype = "solid", show.legend = FALSE,
                   size = 2, height = 0.01*mean(mod.obj$dataN$response), inherit.aes = FALSE,
                   color = brewer.pal(9, "Set1")[3]) +
    geom_point(data = dplyr::filter(BMDBMA, Type == weight_type),
               aes(x = BMD,
                   y = respBMD$model_averaged$model_averaged_response,#min(preds$model_averaged$model_averaged),
                   group = Model),
               size = 5, shape = 19,
               color = brewer.pal(9, "Set1")[1],
               show.legend = FALSE,
               inherit.aes = FALSE)  +
    scale_y_continuous(expression(response),
                       sec.axis = sec_axis(#~.*1.2,
                         # ~ ((. /max(BMDMixture$y)) *(cmax*1.2)),
                         ~ (. - a)/b,
                         name = "Rescaled Density", #trans = 'log10',
                         labels = scales::comma)
    ) +
    labs(x = expression(dose)) +
    theme_minimal() +
    coord_cartesian(xlim = c(min(preds_min$Dose*mod.obj$max.dose),
                             2*mod.obj$max.dose),
                    ylim =  c(min(BMDMixture2$yres)*0.9, max(BMDMixture2$yres)*1.1)) +
    # scale_x_continuous(trans = 'log10', labels = scales::comma,
    #                    breaks = orig_ptdata$dose2[2:length(orig_ptdata$dose2)]*mod.obj$max.dose) +
    scale_x_continuous(trans = 'log10', labels = plot.labs,
                       breaks = ddd*mod.obj$max.dose) +
    theme(strip.text = element_text(size = 15, face = "bold"),
          axis.text = element_text(size = 10, face = "bold"),
          axis.title = element_text(size = 15, face = "bold"),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 15, face = "bold"),
          panel.spacing = unit(5, "lines"),
          legend.position = "top",
          legend.direction = "horizontal",
          # title = element_text(size = 15, face = "bold")
    ) +
    scale_fill_manual(values = unname(mods_fills), name = NULL)#, labels = c(""))


  if(include_data == TRUE) {

    # total weight for Normal and Lognormal dist

    if(is.BMADR2(mod.obj)[3]==2 & weight_type =="BS"){
      weightsDist = c(sum(mod.obj$weights_bridge_sampling[grepl("_N", names(mod.obj$weights_bridge_sampling))]),
                      sum(mod.obj$weights_bridge_sampling[grepl("_LN", names(mod.obj$weights_bridge_sampling))]))
    }else if(is.BMADR2(mod.obj)[3]==2 & weight_type =="LP"){
      weightsDist = c(sum(mod.obj$weights_laplace[grepl("_N", names(mod.obj$weights_laplace))]),
                      sum(mod.obj$weights_laplace[grepl("_LN", names(mod.obj$weights_laplace))]))
    }else{
      weightsDist = c(sum(mod.obj$weights[grepl("_N", names(mod.obj$weights))]),
                      sum(mod.obj$weights[grepl("_LN", names(mod.obj$weights))]))
    }

    if(clustered == F){
      # decide which data to show based on highest weights
      if(weightsDist[1] > weightsDist[2]){
        data.plot <- data.frame(dose = orig_ptdataN$dose,
                                dose2 = orig_ptdataN$dose2,
                                lg10d = orig_ptdataN$lg10d,
                                m = orig_ptdataN$m,
                                min.y = orig_ptdataN$m - orig_ptdataN$sd,
                                max.y = orig_ptdataN$m + orig_ptdataN$sd
        )
        w.data = "arithmetic"
      }else if(weightsDist[2] > weightsDist[1]){
        data.plot <- data.frame(dose = orig_ptdataLN$dose,
                                dose2 = orig_ptdataLN$dose2,
                                lg10d = orig_ptdataLN$lg10d,
                                m = 10^orig_ptdataLN$log10m,
                                min.y = 10^(orig_ptdataLN$log10m - orig_ptdataLN$log10s),
                                max.y = 10^(orig_ptdataLN$log10m + orig_ptdataLN$log10s)
        )
        w.data = "geometric"
      }
    }else if(clustered == T){
      # decide which data to show based on highest weights
      if(weightsDist[1] > weightsDist[2]){
        data.plot <- data.frame(dose = orig_ptdata$dose,
                                dose2 = orig_ptdata$dose2,
                                lg10d = orig_ptdata$lg10d,
                                y = orig_ptdata$y
        )
        data.litter <- data.frame(dose = means.litter$dose,
                                  resp = means.litter$mresp)
        data.all <- data.frame(dose = means.all$dose,
                               resp = means.all$mresp)
        w.data = "arithmetic"
      }else if(weightsDist[2] > weightsDist[1]){
        data.plot <- data.frame(dose = orig_ptdata$dose,
                                dose2 = orig_ptdata$dose2,
                                lg10d = orig_ptdata$lg10d,
                                y = 10^log10(exp(orig_ptdata$yl))

        )
        data.litter <- data.frame(dose = means.litter$dose,
                                  resp = 10^log10(exp(means.litter$mlresp)))
        data.all <- data.frame(dose = means.all$dose,
                               resp = 10^log10(exp(means.all$mlresp)))
        w.data = "geometric"
      }
    }

    if(clustered == F){
      dplot2 <- dplot +
        geom_errorbar(data = data.plot,
                      mapping = aes(x = dose2*mod.obj$max.dose, ymin = min.y,
                                    ymax = max.y),
                      width = NA, position = pd, size = 1,
                      inherit.aes = FALSE) +
        geom_point(data = data.plot, mapping = aes(x = dose2*mod.obj$max.dose, y = m) ,
                   size = 2, color = 1, shape = 21,
                   fill = brewer.pal(9, "Set1")[2],
                   inherit.aes = FALSE) +
        labs(caption = paste0("data and vertical bars based on ", w.data, " sample means and standard deviations"))
    }else if(clustered == T){
      dplot2 <- dplot +
        geom_jitter(data = orig_ptdata, mapping = aes(x = dose2*mod.obj$max.dose, y = y),
                    size = 1, color = 3, shape = 20,
                    # fill = brewer.pal(9, "Set1")[2],
                    position = position_jitter(h = 0.01, width = 0),
                    inherit.aes = FALSE) +
        geom_point(data = data.litter, mapping = aes(x = dose*mod.obj$max.dose, y = resp),
                   size = 2, color = 1, shape = 21,
                   fill = 1,
                   inherit.aes = FALSE) +
        geom_point(data = data.all, mapping = aes(x = dose*mod.obj$max.dose, y = resp),
                   size = 3, color = 2, shape = 23,
                   fill = 2,
                   inherit.aes = FALSE) +
        labs(caption = paste0("diamonds represent the ", w.data, " sample mean"))
    }

    if((TRUE %in% grepl('_LN', mod.obj$models_included)) && (TRUE %in% grepl('_N', mod.obj$models_included))){
      pts1 <- ggpubr::ggarrange(pBMDs, pWeights, pplotN, pplotLN, pplot,
                                dplot2,
                                nrow = 3, ncol = 2)
      pts2 <- ggpubr::annotate_figure(pts1, top = text_grob(title,
                                                            color = "red", face = "bold", size = 14))
    }else if(!(TRUE %in% grepl('_LN', mod.obj$models_included)) && (TRUE %in% grepl('_N', mod.obj$models_included))){
      pts1 <- ggpubr::ggarrange(pBMDs, pWeights, pplotN, pplot,
                                dplot2,
                                nrow = 3, ncol = 2)
      pts2 <- ggpubr::annotate_figure(pts1, top = text_grob(title,
                                                            color = "red", face = "bold", size = 14))
    }else if((TRUE %in% grepl('_LN', mod.obj$models_included)) && !(TRUE %in% grepl('_N', mod.obj$models_included))){
      pts1 <- ggpubr::ggarrange(pBMDs, pWeights, pplotLN, pplot,
                                dplot2,
                                nrow = 3, ncol = 2)
      pts2 <- ggpubr::annotate_figure(pts1, top = text_grob(title,
                                                            color = "red", face = "bold", size = 14))
    }


  } else {

    if((TRUE %in% grepl('_LN', mod.obj$models_included)) && (TRUE %in% grepl('_N', mod.obj$models_included))){
      pts1 <- ggpubr::ggarrange(pBMDs, pWeights, pplotN, pplotLN, pplot,
                                dplot2,
                                nrow = 3, ncol = 2)
      pts2 <- ggpubr::annotate_figure(pts1, top = text_grob(title,
                                                            color = "red", face = "bold", size = 14))
    }else if(!(TRUE %in% grepl('_LN', mod.obj$models_included)) && (TRUE %in% grepl('_N', mod.obj$models_included))){
      pts1 <- ggpubr::ggarrange(pBMDs, pWeights, pplotN, pplot,
                                dplot2,
                                nrow = 3, ncol = 2)
      pts2 <- ggpubr::annotate_figure(pts1, top = text_grob(title,
                                                            color = "red", face = "bold", size = 14))
    }else if((TRUE %in% grepl('_LN', mod.obj$models_included)) && !(TRUE %in% grepl('_N', mod.obj$models_included))){
      pts1 <- ggpubr::ggarrange(pBMDs, pWeights, pplotLN, pplot,
                                dplot2,
                                nrow = 3, ncol = 2)
      pts2 <- ggpubr::annotate_figure(pts1, top = text_grob(title,
                                                            color = "red", face = "bold", size = 14))
    }

  }


  if(all == TRUE){
    return(pts2)
  }else{
    if(include_data == TRUE){
      if((TRUE %in% grepl('_LN', mod.obj$models_included)) && (TRUE %in% grepl('_N', mod.obj$models_included))){
        return(list(BMDs = pBMDs, weights = pWeights, model_fit_N = pplotN,
                    model_fit_LN = pplotLN, model_fit = pplot, MA_fit = dplot2))
      }else if(!(TRUE %in% grepl('_LN', mod.obj$models_included)) && (TRUE %in% grepl('_N', mod.obj$models_included))){
        return(list(BMDs = pBMDs, weights = pWeights, model_fit_N = pplotN,
                    model_fit = pplot, MA_fit = dplot2))
      }else if((TRUE %in% grepl('_LN', mod.obj$models_included)) && !(TRUE %in% grepl('_N', mod.obj$models_included))){
        return(list(BMDs = pBMDs, weights = pWeights,
                    model_fit_LN = pplotLN, model_fit = pplot, MA_fit = dplot2))
      }

    }else{
      if((TRUE %in% grepl('_LN', mod.obj$models_included)) && (TRUE %in% grepl('_N', mod.obj$models_included))){
        return(list(BMDs = pBMDs, weights = pWeights, model_fit_N = pplotN,
                    model_fit_LN = pplotLN, model_fit = pplot, MA_fit = dplot2))
      }else if(!(TRUE %in% grepl('_LN', mod.obj$models_included)) && (TRUE %in% grepl('_N', mod.obj$models_included))){
        return(list(BMDs = pBMDs, weights = pWeights, model_fit_N = pplotN,
                    model_fit = pplot, MA_fit = dplot2))
      }else if((TRUE %in% grepl('_LN', mod.obj$models_included)) && !(TRUE %in% grepl('_N', mod.obj$models_included))){
        return(list(BMDs = pBMDs, weights = pWeights,
                    model_fit_LN = pplotLN, model_fit = pplot, MA_fit = dplot2))
      }
    }
  }

}



#' plot function for the BMD estimates, model weights, model predictions and model average predictions (Quantal dat)
#'
#' @param mod.obj BMDBMA model object
#' @param weight_type type of model weight to be plotted. It can either be "BS" or "LP"
#' @param include_data logical argument to indicate if data should be included in the plots. Defaults to TRUE
#' @param all logical argument to indicate if all plots should be displayed
#' @param title title of the plot
#'
#' @return object of class ggplot
#' @export plot.BMADRQ
#'
plot.BMADRQ <- function(mod.obj,
                        weight_type = c("BS", "LP"),
                        include_data = TRUE,
                        all = TRUE, title
) {
  type <- 'quantal'
  weight_type <- match.arg(weight_type)
  q <- mod.obj$q

  refactor <- function(x, type){

    x <- dplyr::mutate(x, Model = factor(x$Model,
                                         levels = c("Model Averaged", "E4_Q", "IE4_Q", "H4_Q",
                                                    "LN4_Q", "G4_Q",  "QE4_Q", "P4_Q", "L4_Q"),
                                         labels = c("Model Averaged", "Exp", "InvExp", "Hill",
                                                    "LogNormal", "Gamma",  "QuadExp", "Probit",
                                                    "Logistic")
    ))

    return(x)

  }

  BMDW <- BMDWeights(mod.obj, type = type)
  BMDBMA <- BMDMAQ_extract(mod.obj, conv = FALSE)
  mod.obj$data <- mod.obj$data[order(mod.obj$data$dose), ]

  dose <- sort(unique(mod.obj$data$dose)/max(mod.obj$data$dose))
  if(min(dose) == 0){
    ddd <- c(min(dose[dose > 0])/4, dose[2:length(dose)])
    lg10d <- c(log10(min(dose[dose > 0])/4), log10(dose[2:length(dose)]))
    bmdl <- log10(BMDW$BMDL/mod.obj$max.dose)
    bmdlo <- BMDW$BMDL/mod.obj$max.dose

    lg10d[1] <- ifelse((min(bmdl, na.rm=T) < lg10d[1] & min(bmdl, na.rm=T)!='-Inf' ),
                       log10(min(BMDW$BMDL/mod.obj$max.dose/2, na.rm=T)),
                       log10(min(dose[dose > 0])/4))
    ddd[1] <- ifelse(min(bmdlo, na.rm=T) < ddd[1],
                     min(BMDW$BMDL/mod.obj$max.dose/2, na.rm=T),
                     min(dose[dose > 0])/4)
  }else{
    ddd <- c(min(dose)/4, dose)
    lg10d <- c(log10(min(dose)/4), log10(dose))
    bmdl <- log10(BMDW$BMDL/mod.obj$max.dose)
    bmdlo <- BMDW$BMDL/mod.obj$max.dose
    lg10d[1] <- ifelse((min(bmdl, na.rm=T) < lg10d[1] & min(bmdl, na.rm=T)!='-Inf' ),
                       log10(min(BMDW$BMDL/mod.obj$max.dose/2, na.rm=T)),
                       log10(min(dose)/4))
    ddd[1] <- ifelse(min(bmdlo, na.rm=T) < ddd[1],
                     min(BMDW$BMDL/mod.obj$max.dose/2, na.rm=T),
                     min(dose)/4)
  }


  orig.p <- (mod.obj$data$y/mod.obj$data$n)
  orig.ps <- log(orig.p + .Machine$double.xmin)
  orig.s <- log(sqrt(orig.p*(1-orig.p)/mod.obj$data$n))

  mps <- dplyr::ungroup(dplyr::summarise(dplyr::group_by(mod.obj$data, dose), mp = sum(y)/sum(n),
                                         n = sum(n)) )

  #dplyr::summarise(dplyr::group_by(mod.obj$data, dose), n = sum(n))
  if(min(dose) == 0){
    orig_ptdata <- merge(data.frame(dose = mod.obj$data$dose,
                                    dose2 = rep(ddd, times = table(mod.obj$data$dose)),
                                    lg10d = rep(lg10d, times = table(mod.obj$data$dose))
    ), mps, by = "dose")
    mod.obj$data$lg10d <- rep(lg10d, times = table(mod.obj$data$dose))
    mod.obj$data$dose2 <- rep(ddd, times = table(mod.obj$data$dose))
  }else{
    orig_ptdata <- merge(data.frame(dose = mod.obj$data$dose,
                                    dose2 = rep(ddd[2:(length(unique(mod.obj$data$dose))+1)], times = table(mod.obj$data$dose)),
                                    lg10d = rep(lg10d[2:(length(unique(mod.obj$data$dose))+1)], times = table(mod.obj$data$dose))
    ), mps, by = "dose")
    mod.obj$data$lg10d <- rep(lg10d[2:(length(unique(mod.obj$data$dose))+1)], times = table(mod.obj$data$dose))
    mod.obj$data$dose2 <- rep(ddd[2:(length(unique(mod.obj$data$dose))+1)], times = table(mod.obj$data$dose))
  }


  dgr <- seq(min(lg10d), abs(min(lg10d)), by=0.01)
  #dgr2 <- seq(min(ddd), abs(max(ddd)), by=0.01)
  dose2 <- 10^dgr[(dgr > (lg10d[2]-((lg10d[2]-lg10d[1])/2)))]
  preds <- predict.BMADRQ(mod.obj, dose = dose2,
                          what = "predicted",
                          model_averaged = TRUE,
                          weight_type = weight_type)

  preds$predicted$dgr <- dgr[(dgr > (lg10d[2]-((lg10d[2]-lg10d[1])/2)))]
  preds$model_averaged$dgr <- dgr[(dgr > (lg10d[2]-((lg10d[2]-lg10d[1])/2)))]

  preds$predicted$lg10pred <- log10(preds$predicted$predicted)

  dgrprime <- dgr[(dgr <= (lg10d[2]-((lg10d[2]-lg10d[1])/2)))]
  preds_min <- tidyr::nest(dplyr::group_by(dplyr::filter(preds$predicted, Dose == min(Dose)), Model))
  as_per_model <- dplyr::select(par_med(mod.obj, type = 'quantal'), Model, a)
  preds_min$data <- lapply(preds_min$data, function(x) data.frame(dgrprime = dgrprime,
                                                                  Dose = 10^dgrprime,
                                                                  predicted = x$predicted,
                                                                  lg10predicted = x$lg10pred))
  preds_min <- tidyr::unnest(preds_min, cols = c('data'))
  preds_min <- merge(preds_min, as_per_model, by = 'Model', sort = FALSE)

  preds_min2 <- data.frame(Dose = 10^dgrprime,
                           lg10d = dgrprime,
                           log10MA = rep(log10(preds$model_averaged$model_averaged[
                             preds$model_averaged$Dose==min(preds$model_averaged$Dose)]),
                             length(dgrprime)),
                           MA = rep(preds$model_averaged$model_averaged[
                             preds$model_averaged$Dose==min(preds$model_averaged$Dose)],
                             length(dgrprime)))

  #preds$predicted$Dose <- preds$predicted$Dose*mod.obj$max.dose
  #preds$model_averaged$Dose <- preds$model_averaged$Dose*mod.obj$max.dose

  respBMD <- predict.BMADRQ(mod.obj,
                            what = "resp_at_BMD",
                            model_averaged = TRUE,
                            weight_type = weight_type)
  respBMD$resp_at_BMD$lg10bmd <- log10(respBMD$resp_at_BMD$BMD)
  #respBMD$resp_at_BMD[,"BMD"] = respBMD$resp_at_BMD[,"BMD"]*mod.obj$max.dose
  #print(respBMD)
  #respBMD <- refactor(respBMD)

  BMDMixture <- BMDQmixture_extract(mod.obj, weight_type, conv=FALSE) # BMD values
  BMDMixture$BMDMixture2 <- log10(BMDMixture$BMDMixture/mod.obj$max.dose)

  gghst2 <- hist(BMDMixture$BMDMixture, breaks = sqrt(nrow(BMDMixture)), plot = FALSE) #hist on original scale

  lais <- log10(gghst2$breaks) # log10 of the breaks
  dlais <- abs(diff(lais)) # width of the interval on log scale
  glais <- (gghst2$counts/dlais) # divide the original counts by the width on the log scale
  glais2 <- glais/sum(glais) * nrow(BMDMixture) # normalise the new frequencies

  BMDMixture2 <- data.frame(Model = unique(BMDMixture$Model),#rep(unique(BMDMixture$Model), length(gghst2$counts)),
                            # Dose = gghst2$mids, #midpoints
                            Dose = 10**lais[1:length(gghst2$counts)], #midpoints
                            y = gghst2$counts, #frequencies
                            y2 = glais2
  )
  # BMDMixtureD <- density(BMDMixture$BMDMixture2, # density of BMD mixture on log10 scale
  #                        from = min(BMDMixture$BMDMixture2),
  #                        to = max(BMDMixture$BMDMixture2))
  #
  # BMDMixtureD2 <- density(BMDMixture$BMDMixture, # density of BMD mixture on 0-1 scale
  #                         from = min(BMDMixture$BMDMixture),
  #                         to = max(BMDMixture$BMDMixture))
  #
  #
  # BMDMixture2 <- data.frame(Model = unique(BMDMixture$Model), # y = Density
  #                           Dose10 = BMDMixtureD$x,
  #                           Dose102 = 10^BMDMixtureD$x,
  #                           y10 = BMDMixtureD$y,
  #                           Dose = BMDMixtureD2$x,
  #                           y = BMDMixtureD2$y
  # )
  #BMDMixture <- refactor(BMDMixture)
  #print(BMDMixture)
  #BMDs and Weights
  respBMDBMDW <- merge(respBMD$resp_at_BMD, BMDW, by = c("Model", "BMD"), sort = FALSE)
  respBMDBMDW <- refactor(respBMDBMDW, type)
  pd <- position_dodge(0.5)
  # dist_fills <- c("N" = "#EFC000FF", "LN" = "#0073C2FF")
  # dist_names <- c("N" = "Normal", "LN" = "LogNormal")
  #dist_fills <- c("#EFC000FF", "#0073C2FF")
  #names(dist_fills) <- c("N","LN")
  #dist_names <- c("Normal","LogNormal")
  #names(dist_names) <- c("N","LN")
  if(mod.obj$max.dose > 1){
    respBMDBMDW$BMDU <- ifelse(((respBMDBMDW$BMDU > mod.obj$max.dose^2) &
                                  ((respBMDBMDW$BMDU/respBMDBMDW$BMD)/(respBMDBMDW$BMD/respBMDBMDW$BMDL)) > 10) |
                                 respBMDBMDW$BMDU > 10*mod.obj$max.dose^2,
                               # mod.obj$max.dose^2,
                               2*(respBMDBMDW$BMD - respBMDBMDW$BMDL),
                               respBMDBMDW$BMDU)
  }
  #BMD plots
  pBMDs <- ggplot(data = respBMDBMDW, aes(x = BMD, y = Model, group = Model)) +
    geom_errorbarh(aes(xmin = BMDL, xmax = BMDU, y = Model),
                   linetype = "solid", show.legend = FALSE,
                   size = 3, height = 1.2, color = brewer.pal(8, "Set2")[8]) +
    geom_point(aes(x = BMD, y = Model), fill = "#0073C2FF",
               size = 7, color = 1, shape = 21) +
    theme_minimal() +
    labs(x = expression(BMD), y = "") +
    labs(x="BMD on original scale") +
    theme(strip.text = element_text(size = 15, face = "bold"),
          axis.text = element_text(size = 10, face = "bold"),
          axis.title = element_text(size = 15, face = "bold"),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 15, face = "bold"),
          panel.spacing = unit(5, "lines"),
          legend.position = "top",
          legend.direction = "horizontal",
          title = element_text(size = 15, face = "bold")) +
    geom_errorbarh(data = BMDBMA[BMDBMA$Type == weight_type,],
                   aes(xmin = BMDL, xmax = BMDU, y = Model),
                   linetype = "solid", show.legend = FALSE,
                   size = 3, height = 1.2, color = brewer.pal(9, "Set1")[3]) +
    geom_point(data = BMDBMA[BMDBMA$Type == weight_type,],
               aes(x = BMD, y = Model),
               size = 7, color = 1, shape = 21,
               fill = brewer.pal(9, "Set1")[1]) +
    #scale_x_log10() +
    scale_y_discrete(limits = rev(c("Model Averaged", "Exp", "InvExp", "Hill",
                                    "LogNormal", "Gamma",  "QuadExp", "Probit",
                                    "Logistic")),
                     labels = rev(c("Model Averaged", "Exp(Q)", "InvExp(Q)", "Hill(Q)",
                                    "LogNormal(Q)", "Gamma(Q)",  "QuadExp(Q)", "Probit(Q)",
                                    "Logistic(Q)")))
  #Weights plot
  needed_data <- tidyr::separate(BMDW,
                                 col = "Model", sep = "_",
                                 into = c("Model", "Distribution"))

  if(weight_type == "LP") {

    pWeights <- ggpubr::ggdotchart(
      needed_data, x = "Model", y = "LP_Weights",
      group = "Distribution", color = "#0073C2FF",
      add = "segment", position = position_dodge(0.3),
      sorting = "descending", size = 3, dot.size = 7) +
      labs(x = "", y = "Weight",
           title = "") +
      theme_minimal() +
      theme(strip.text = element_text(size = 15, face = "bold"),
            axis.text = element_text(size = 10, face = "bold"),
            axis.title = element_text(size = 15, face = "bold"),
            legend.text = element_text(size = 10),
            legend.title = element_text(size = 15, face = "bold"),
            panel.spacing = unit(5, "lines"),
            legend.position = "top",
            legend.direction = "horizontal",
            title = element_text(size = 15, face = "bold")) +
      scale_x_discrete(
        limits = c("E4", "IE4", "H4", "LN4", "G4", "QE4", "P4", "L4"),
        labels = c("Exp", "InvExp", "Hill", "LogNormal", "Gamma",
                   "QuadExp", "Probit", "Logistic")
      )

  } else {

    pWeights <- ggpubr::ggdotchart(
      needed_data, x = "Model", y = "BS_Weights",
      group = "Distribution", color = "#0073C2FF",
      add = "segment", position = position_dodge(0.3),
      sorting = "descending", size = 3, dot.size = 7) +
      labs(x = "", y = "Weight", title = "") +
      theme_minimal() +
      theme(strip.text = element_text(size = 15, face = "bold"),
            axis.text = element_text(size = 10, face = "bold"),
            axis.title = element_text(size = 15, face = "bold"),
            legend.text = element_text(size = 10),
            legend.title = element_text(size = 15, face = "bold"),
            panel.spacing = unit(5, "lines"),
            legend.position = "top",
            legend.direction = "horizontal",
            title = element_text(size = 15, face = "bold")) +
      scale_x_discrete(
        limits = c("E4", "IE4", "H4", "LN4", "G4", "QE4", "P4", "L4"),
        labels = c("Exp", "InvExp", "Hill", "LogNormal", "Gamma",
                   "QuadExp", "Probit", "Logistic")
      )

  }

  ## Mixture density plots

  #Prediction plot
  mods_fills <- c("Model-averaged BMD" = "coral")
  #print(BMDBMA)
  if(is.BMADRQ2(mod.obj)[3]==2) {

    BMDBMA$Response <- c(sum(respBMDBMDW$resp_at_BMD*respBMDBMDW$BS_Weights),
                         sum(respBMDBMDW$resp_at_BMD*respBMDBMDW$LP_Weights))

  } else if(is.BMADRQ2(mod.obj)[2]==2) {

    BMDBMA$Response <- sum(respBMDBMDW$resp_at_BMD*respBMDBMDW$LP_Weights)

  } else {

    pts <- ggpubr::ggarrange(pBMDs, pWeights, nrow = 1, ncol = 2)
    return(pts)
    stop("cannot compute model averaged response at BMD. Please check the inputs.")
  }

  #BMDBMA$Response <- c(sum(respBMDBMDW$resp_at_BMD*respBMDBMDW$BS_Weights),
  #sum(respBMDBMDW$resp_at_BMD*respBMDBMDW$LP_Weights))

  preds2 <- tidyr::separate(preds$predicted,
                            col = "Model", sep = "_",
                            into = c("Model", "Distribution"))
  preds_min <- tidyr::separate(preds_min,
                               col = "Model", sep = "_",
                               into = c("Model", "Distribution"))
  lty <- c("Q" = "solid")
  md_cls <- RColorBrewer::brewer.pal(8, "Dark2")
  names(md_cls) <-  c("E4", "IE4", "H4",
                      "LN4", "G4", "QE4",
                      "P4", "L4")

  #ymin = min(10^mod.obj$data$lg10m - 2*exp(orig.s), na.rm=T)
  #ymax = max(10^mod.obj$data$lg10m + 2*exp(orig.s), na.rm=T)

  if(min(mod.obj$data$dose) == 0){
    plot.labs = dose*mod.obj$max.dose
  }else{
    plot.labs = c(0, dose*mod.obj$max.dose)
  }

  pplot <- ggplot(data = preds2,
                  aes(x = Dose*mod.obj$max.dose, y = predicted, group = Model,
                      color = Model)) +
    geom_line(#alpha = 0.6,
      size = 1,
      show.legend = TRUE) +
    labs(color = "Model",  x = expression(dose),
         y = expression(p(y==1)), title = "") +

    geom_segment(data = preds_min, mapping = aes(x = Dose[1]*mod.obj$max.dose, y = a,
                                                 xend = max(Dose*mod.obj$max.dose),
                                                 yend = predicted,
                                                 group = interaction(Model, Distribution),
                                                 color = Model),
                 linetype = "dotted", alpha = 0.6,
                 size = 0.8, inherit.aes = FALSE, show.legend = FALSE) +
    geom_errorbarh(data = dplyr::filter(BMDBMA, Type == weight_type),
                   aes(xmin = BMDL, xmax = BMDU,
                       y = Response,
                       group = Model),
                   linetype = "solid", show.legend = FALSE,
                   size = 2, height = 0.03, inherit.aes = FALSE,
                   color = brewer.pal(9, "Set1")[3]) +
    geom_point(data = dplyr::filter(BMDBMA, Type == weight_type),
               aes(x = BMD, y = Response, group = Model),
               size = 5, shape = 19,
               color = brewer.pal(9, "Set1")[1],
               show.legend = FALSE,
               inherit.aes = FALSE) +
    coord_cartesian(xlim = c(min(preds_min$Dose*mod.obj$max.dose),
                             2*mod.obj$max.dose)) +
    scale_color_manual(values = md_cls,
                       labels = c("Exp", "InvExp", "Hill", "LogNormal", "Gamma",
                                  "QuadExp", "Probit", "Logistic")) +
     # scale_x_continuous(trans = 'log10', labels = scales::comma,
     #                   breaks = orig_ptdata$dose2[2:length(orig_ptdata$dose2)]*mod.obj$max.dose) +
    scale_x_continuous(trans = 'log10', labels = plot.labs,
                       breaks = ddd*mod.obj$max.dose) +
    theme_minimal() +
    theme(strip.text = element_text(size = 15, face = "bold"),
          axis.text = element_text(size = 10, face = "bold"),
          axis.title = element_text(size = 15, face = "bold"),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 12, face = "bold"),
          panel.spacing = unit(5, "lines"),
          legend.position = "top",
          legend.direction = "horizontal",
          title = element_text(size = 15, face = "bold"))

  cmax <- max(preds$model_averaged$model_averaged)
  cmin <- min(preds$model_averaged$model_averaged)
  # BMDMixture$yres <- (BMDMixture$y/max(BMDMixture$y)) * (cmax * 1.2) # Density rescaled

  ylim.prim <- c(cmin, cmax)
  # ylim.sec <- c(min(BMDMixture2$y2), max(BMDMixture2$y2))
  ylim.sec <- c(min(BMDMixture2$y), max(BMDMixture2$y))
  b <- diff(ylim.prim)/diff(ylim.sec)
  a <- ylim.prim[1] - b*ylim.sec[1]

  # BMDMixture2$yres = a + BMDMixture2$y2*b
  BMDMixture2$yres = a + BMDMixture2$y*b

  dplot <- ggplot(data = preds$model_averaged, aes(x = Dose*mod.obj$max.dose, y = model_averaged,
                                                   group = 1)) +
    geom_line(show.legend = FALSE, linetype = "dashed", size = 3) +
    #geom_col(data = BMDMixture2, aes(x = Dose, y = yres, fill = Model), alpha = 0.6,
    #         inherit.aes = FALSE, color = NA) +
    #geom_histogram(data = BMDMixture,
    #               aes(x = BMDMixture, y = ..density.., fill = Model),
    #               color = NA, alpha = 0.5, inherit.aes = FALSE, bins = sqrt(nrow(BMDMixture))) +
    geom_ribbon(data = BMDMixture2, aes(x = Dose,
                                        y = yres,
                                        ymin = min(yres), ymax = yres, fill = Model),
                color=NA, alpha = 0.5) +
    #geom_area(data = BMDMixture, aes(x = Dose,
    #                                 y = yres, fill = Model),
    #            color=NA, alpha = 0.5) +
    geom_segment(data = preds_min2, mapping = aes(x = Dose[1]*mod.obj$max.dose, y = MA,
                                                  xend = max(Dose*mod.obj$max.dose),
                                                  yend = preds$model_averaged$model_averaged[1]
                                                  # group = interaction(Model, Distribution),
                                                  # color = Model
    ),
    linetype = "dotted",
    size = 3, inherit.aes = FALSE, show.legend = FALSE) +
    geom_errorbarh(data = dplyr::filter(BMDBMA, Type == weight_type),
                   aes(xmin = BMDL, xmax = BMDU,
                       y = respBMD$model_averaged$model_averaged_response,#min(preds$model_averaged$model_averaged),
                       group = Model),
                   linetype = "solid", show.legend = FALSE,
                   size = 2, height = 0.01*mean(mod.obj$data$y), inherit.aes = FALSE,
                   color = brewer.pal(9, "Set1")[3]) +
    geom_point(data = dplyr::filter(BMDBMA, Type == weight_type),
               aes(x = BMD,
                   y = respBMD$model_averaged$model_averaged_response,#min(preds$model_averaged$model_averaged),
                   group = Model),
               size = 5, shape = 19,
               color = brewer.pal(9, "Set1")[1],
               show.legend = FALSE,
               inherit.aes = FALSE)  +
    scale_y_continuous(expression(p(y==1)),
                       sec.axis = sec_axis(#~.*1.2,
                         # ~ ((. /max(BMDMixture$y)) *(cmax*1.2)),
                         ~ (. - a)/b,
                         name = "Rescaled Density",
                         labels = scales::comma)
    ) +
    labs(x = expression(dose)) +
    theme_minimal() +
    coord_cartesian(xlim = c(min(preds_min$Dose*mod.obj$max.dose),
                             2*mod.obj$max.dose),
                    ylim = c(min(BMDMixture2$yres), max(BMDMixture2$yres)) ) +
    # scale_x_continuous(trans = 'log10', labels = scales::comma,
    #                    breaks = orig_ptdata$dose2[2:length(orig_ptdata$dose2)]*mod.obj$max.dose) +
    # scale_x_continuous(trans = 'log10', labels = dose*mod.obj$max.dose,
    #                    breaks = ddd*mod.obj$max.dose) +
    scale_x_continuous(trans = 'log10', labels = plot.labs,
                       breaks = ddd*mod.obj$max.dose) +
    theme(strip.text = element_text(size = 15, face = "bold"),
          axis.text = element_text(size = 10, face = "bold"),
          axis.title = element_text(size = 15, face = "bold"),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 15, face = "bold"),
          panel.spacing = unit(5, "lines"),
          legend.position = "top",
          legend.direction = "horizontal"#,
          # title = element_text(size = 15, face = "bold")
    ) +
    scale_fill_manual(values = unname(mods_fills), name = NULL)#, labels = c(""))

  if(include_data == TRUE) {


    orig_ptdata$ymin <- orig_ptdata$mp - sqrt(orig_ptdata$mp*(1-orig_ptdata$mp))
    orig_ptdata$ymax <- orig_ptdata$mp + sqrt(orig_ptdata$mp*(1-orig_ptdata$mp))


    pplot2 <- pplot +
      geom_point(data = orig_ptdata, mapping = aes(x = (dose2*mod.obj$max.dose)+.Machine$double.xmin, y = mp) ,
                 size = 4, color = 1, shape = 21,
                 fill = brewer.pal(9, "Set1")[2],
                 inherit.aes = FALSE)

    dplot2 <- dplot +
      geom_point(data = orig_ptdata, mapping = aes(x = (dose2*mod.obj$max.dose)+.Machine$double.xmin, y = mp) ,
                 size = 6, color = 1, shape = 21,
                 fill = brewer.pal(9, "Set1")[2],
                 inherit.aes = FALSE)

    if(mod.obj$is_bin == 0){

      pplot3 <- pplot2 +
        geom_jitter(data = mod.obj$data, mapping = aes(x = dose2*mod.obj$max.dose, y = y/n) ,
                    size = 2.5, color = 1, shape = 23,
                    fill = brewer.pal(9, "Set1")[3], width = 0.1,
                    inherit.aes = FALSE)

      dplot3 <- dplot2 +
        geom_jitter(data = mod.obj$data, mapping = aes(x = dose2*mod.obj$max.dose, y = y/n) ,
                    size = 2.5, color = 1, shape = 23,
                    fill = brewer.pal(9, "Set1")[3], width = 0.1,
                    inherit.aes = FALSE)

      pts1 <- ggpubr::ggarrange(pBMDs, pWeights, pplot2, pplot3,
                                dplot2, dplot3,
                                nrow = 2, ncol = 2)
      pts2 <- ggpubr::annotate_figure(pts1, top = text_grob(title,
                                                            color = "red", face = "bold", size = 14))

    } else {
      pts1 <- ggpubr::ggarrange(pBMDs, pWeights, pplot2,
                                dplot2,
                                nrow = 2, ncol = 2)
      pts2 <- ggpubr::annotate_figure(pts1, top = text_grob(title,
                                                            color = "red", face = "bold", size = 14))
    }




  } else {

    pts1 <- ggpubr::ggarrange(pBMDs, pWeights, pplot, dplot, nrow = 2, ncol = 2)
    pts2 <- ggpubr::annotate_figure(pts1, top = text_grob(title,
                                                          color = "red", face = "bold", size = 14))

  }

  if(all == TRUE){
    return(pts2)
  }else{
    if(include_data == TRUE & mod.obj$is_bin == 0){
      return(list(BMDs = pBMDs, weights = pWeights, model_fit2 = pplot2, model_fit = pplot3,
                  MA_fit2 = dplot2, MA_fit = dplot3))
    } else if(include_data == TRUE & mod.obj$is_bin != 0){
      return(list(BMDs = pBMDs, weights = pWeights, model_fit = pplot2, MA_fit = dplot2))
    } else {
      return(list(BMDs = pBMDs, weights = pWeights, model_fit = pplot, MA_fit = dplot))
    }
  }

}


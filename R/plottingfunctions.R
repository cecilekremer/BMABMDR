#' Function for internal use
#'
#' @param DR_df dataframe containing individual-level dose response data
#' @param DRM_df dataframe for dose-response model
#' @param x_DR_df name of the x varaibel in the individual-level dose response data.
#'                Typically dose.
#' @param y_DR_df name of the y variable in the individual-level dose response data.
#' @param x_DRM_df name of the x variable in the dose response curve data. Typically dose
#' @param y_DRM_df name of the y variable in the individual-level dose response data.
#'
#' @return .
#'
DR_scatter <- function(DR_df = NULL, DRM_df = NULL,
                       x_DR_df, y_DR_df,
                       x_DRM_df, y_DRM_df,
                       col_DRM = RColorBrewer::brewer.pal(3, "Set1")[1]) {

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
#' @param x doses
#' @param BMD column name for BMD
#' @param BMDL column name for BMDL
#' @param BMDU column name for BMDU.
#'
#' @return .
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
#' @param to_plot which quantity to plot
#' @param xlab x-axis label
#' @param title title
#'
#' @return .
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
DRM <- function(Model, type = c("increasing", "decreasing")) {

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

  } else stop("Provide type to be either 'increasing' or 'decreasing'")
  return(DRM)
}

#' Function for internal use
#'
#' @param mod.obj fitted models
#' @param weight_type BS or LP
#' @param include_data include data points in plot
#' @param all all 4 plots in one plot
#' @param title title for plot if all = TRUE
#'
#' @return .
#'
plot.BMADR <- function(mod.obj,
                       # type = c("increasing", "decreasing"),
                       weight_type = c("BS", "LP"),
                       include_data = TRUE,
                       all = TRUE, title
) {
  # type <- match.arg(type)
  weight_type <- match.arg(weight_type)
  q <- mod.obj$q
  if(mod.obj$increasing == T){
    type = "increasing"
  }else{type = "decreasing"}

  refactor <- function(x){

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

    return(x)

  }

  BMDW <- BMDWeights(mod.obj)
  BMDBMA <- BMDMA_extract(mod.obj, conv = FALSE)

  dose <- mod.obj$data$dose
  lg10d <- c(log10(dose[2]/4), log10(mod.obj$data$dose[2:nrow(mod.obj$data)]))
  bmdl <- log10(BMDW$BMDL/mod.obj$max.dose)
  lg10d[1] <- ifelse(min(bmdl, na.rm=T) < lg10d[1],
                     log10(min(BMDW$BMDL/mod.obj$max.dose/2, na.rm=T)),
                     log10(dose[2]/4))
  orig.y <- log(NtoLN(mod.obj$data$m,mod.obj$data$sd))[1:length(mod.obj$data$dose)]
  orig.s <- log(NtoLN(mod.obj$data$m,mod.obj$data$sd))[(length(mod.obj$data$dose)+1):(2*length(mod.obj$data$dose))]
  orig_ptdata <- data.frame(dose = mod.obj$data$dose,
                            lg10d = lg10d, m = mod.obj$data$m,
                            sd = mod.obj$data$sd,
                            log10m = log10(exp(orig.y)),
                            log10s = log10(exp(orig.s)))
  mod.obj$data$lg10d <- lg10d
  mod.obj$data$lg10m <- log10(exp(orig.y))
  mod.obj$data$lg10s <- log10(exp(orig.s))
  dgr <- seq(min(lg10d), abs(min(lg10d)), by=0.01)
  dose2 <- 10^dgr[(dgr > (lg10d[2]-((lg10d[2]-lg10d[1])/2)))]
  preds <- predict.BMADR(mod.obj, dose = dose2,
                         type = type, what = "predicted",
                         model_averaged = TRUE,
                         weight_type = weight_type)

  preds$predicted$dgr <- dgr[(dgr > (lg10d[2]-((lg10d[2]-lg10d[1])/2)))]
  preds$model_averaged$dgr <- dgr[(dgr > (lg10d[2]-((lg10d[2]-lg10d[1])/2)))]

  preds$predicted$lg10pred = log10(preds$predicted$predicted)

  dgrprime <- dgr[(dgr <= (lg10d[2]-((lg10d[2]-lg10d[1])/2)))]
  preds_min <- tidyr::nest(dplyr::group_by(dplyr::filter(preds$predicted, Dose == min(Dose)), Model))
  as_per_model <- dplyr::select(par_med(mod.obj), Model, min_response)
  preds_min$data <- lapply(preds_min$data, function(x) data.frame(dgrprime = dgrprime,
                                                                  predicted = x$predicted,
                                                                  lg10predicted = x$lg10pred))
  preds_min <- tidyr::unnest(preds_min, cols = c('data'))
  preds_min <- merge(preds_min, as_per_model, by = 'Model', sort = FALSE)
  preds_min$min_response[stringr::str_detect(preds_min$Model, "_LN")] <- exp(
    preds_min$min_response[stringr::str_detect(preds_min$Model, "_LN")])

  preds_min2 <- data.frame(Dose = dgrprime,
                           log10MA = rep(log10(preds$model_averaged$model_averaged[
                             preds$model_averaged$Dose==min(preds$model_averaged$Dose)]),
                             length(dgrprime)))

  #preds$predicted$Dose <- preds$predicted$Dose*mod.obj$max.dose
  #preds$model_averaged$Dose <- preds$model_averaged$Dose*mod.obj$max.dose

  respBMD <- predict(mod.obj, type = type,
                     what = "resp_at_BMD",
                     model_averaged = TRUE,
                     weight_type = weight_type)
  respBMD$resp_at_BMD$lg10bmd <- log10(respBMD$resp_at_BMD$BMD)
  #respBMD$resp_at_BMD[,"BMD"] = respBMD$resp_at_BMD[,"BMD"]*mod.obj$max.dose
  #print(respBMD)
  #respBMD <- refactor(respBMD)

  BMDMixture <- BMDmixture_extract(mod.obj, conv=FALSE) # BMD values
  BMDMixture$BMDMixture2 <- log10(BMDMixture$BMDMixture/mod.obj$max.dose)

  BMDMixtureD <- density(BMDMixture$BMDMixture2, # density of BMD mixture
                         from = min(BMDMixture$BMDMixture2),
                         to = max(BMDMixture$BMDMixture2))
  BMDMixture <- data.frame(Model = unique(BMDMixture$Model), # y = Density
                           Dose = BMDMixtureD$x,
                           y = BMDMixtureD$y)
  #BMDMixture <- refactor(BMDMixture)

  #BMDs and Weights
  respBMDBMDW <- merge(respBMD$resp_at_BMD, BMDW, by = c("Model", "BMD"), sort = FALSE)
  respBMDBMDW <- refactor(respBMDBMDW)
  respBMDBMDW$Distribution <- ifelse(stringr::str_detect(respBMDBMDW$Model, "(LN)"), "LN", "N")
  pd <- position_dodge(0.5)

  #BMD plots
  pBMDs <- ggplot(data = respBMDBMDW, aes(x = BMD, y = Model, group = Model)) +
    geom_errorbarh(aes(xmin = BMDL, xmax = BMDU, y = Model),
                   linetype = "solid", show.legend = FALSE,
                   size = 3, height = 1.2, color = RColorBrewer::brewer.pal(8, "Set2")[8]) +
    geom_point(aes(x = BMD, y = Model, fill = Distribution),
               size = 7, color = 1, shape = 21) +
    theme_minimal() +
    labs(x = expression(BMD), y = "", fill = "Distribution") +
    scale_fill_manual(values = c("#0073C2FF", "#EFC000FF"),
                      labels = c('LogNormal', 'Normal')) +
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
                   size = 3, height = 1.2, color = RColorBrewer::brewer.pal(9, "Set1")[3]) +
    geom_point(data = BMDBMA[BMDBMA$Type == weight_type,],
               aes(x = BMD, y = Model),
               size = 7, color = 1, shape = 21,
               fill = RColorBrewer::brewer.pal(9, "Set1")[1]) +
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
      scale_color_manual(values = c("#0073C2FF", "#EFC000FF"),
                         labels = c('LogNormal', 'Normal'))

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
      scale_color_manual(values = c("#0073C2FF", "#EFC000FF"),
                         labels = c('LogNormal', 'Normal'))

  }


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
  lty <- c("LN" = "dashed",
           "N" = "solid")
  md_cls <- RColorBrewer::brewer.pal(8, "Dark2")
  names(md_cls) <-  c("E4", "IE4", "H4",
                      "LN4", "G4", "QE4",
                      "P4", "L4")

  ymin = min(mod.obj$data$lg10m - 2*mod.obj$data$lg10s, na.rm=T)
  ymax = max(mod.obj$data$lg10m + 2*mod.obj$data$lg10s, na.rm=T)

  pplot <- ggplot(data = preds2,
                  aes(x = dgr, y = lg10pred, group = interaction(Model, Distribution),
                      color = Model, linetype = Distribution)) +
    geom_line(alpha = 0.6,
              size = 1,
              show.legend = TRUE) +
    labs(color = "Model", linetype = "Distribution", x = expression(log[10](dose)),
         y = expression(log[10](response)), title = "") +

    geom_segment(data = preds_min, mapping = aes(x = dgr[1], y = log10(min_response),
                                                 xend = max(dgr[(dgr <= (lg10d[2]-((lg10d[2]-lg10d[1])/2)))]),
                                                 yend = lg10predicted,
                                                 group = interaction(Model, Distribution),
                                                 color = Model),
                 linetype = "dotted", alpha = 0.6,
                 size = 0.8, inherit.aes = FALSE, show.legend = FALSE) +

    # geom_line(data = preds_min, mapping = aes(x = dgrprime, y = log10(min_response),
    #                                           group = interaction(Model, Distribution),
    #                                           color = Model),
    #           linetype = "dotted", alpha = 0.6,
    #           size = 0.8, inherit.aes = FALSE, show.legend = FALSE) +
    ylim(ymin, ymax) +
    # geom_errorbar(data = orig_ptdata,
    #               mapping = aes(x = lg10d, ymin = log10m - log10s,
    #                             ymax = log10m + log10s),
    #               width = NA, position = pd, size = 1,
    #               inherit.aes = FALSE) +
    # geom_point(data = orig_ptdata, mapping = aes(x = lg10d, y = log10m) ,
    #            size = 2, color = 1, shape = 21,
    #            fill = RColorBrewer::brewer.pal(9, "Set1")[2],
    #            inherit.aes = FALSE) +
    geom_errorbarh(data = dplyr::filter(BMDBMA, Type == weight_type),
                   aes(xmin = log10(BMDL/mod.obj$max.dose), xmax = log10(BMDU/mod.obj$max.dose),
                       y = log10(Response),
                       group = Model),
                   linetype = "solid", show.legend = FALSE,
                   size = 2, height = 0.03, inherit.aes = FALSE,
                   color = RColorBrewer::brewer.pal(9, "Set1")[3]) +
    geom_point(data = dplyr::filter(BMDBMA, Type == weight_type),
               aes(x = log10(BMD/mod.obj$max.dose), y = log10(Response), group = Model),
               size = 5, shape = 19,
               color = RColorBrewer::brewer.pal(9, "Set1")[1],
               show.legend = FALSE,
               inherit.aes = FALSE) +
    scale_linetype_manual(values = lty,
                          labels = c('LogNormal', 'Normal')) +
    scale_color_manual(values = md_cls,
                       labels = c("Exp", "InvExp", "Hill", "LogNormal", "Gamma",
                                  "QuadExp", "Probit", "Logistic")) +
    theme_minimal() +
    theme(strip.text = element_text(size = 15, face = "bold"),
          axis.text = element_text(size = 10, face = "bold"),
          axis.title = element_text(size = 15, face = "bold"),
          legend.text = element_text(size = 8, face = "bold"),
          legend.title = element_text(size = 12, face = "bold"),
          panel.spacing = unit(5, "lines"),
          legend.position = "top",
          legend.direction = "horizontal",
          title = element_text(size = 15, face = "bold"))

  cmax <- log10(max(preds$model_averaged$model_averaged))
  cmin <- log10(min(preds$model_averaged$model_averaged))
  # BMDMixture$yres <- (BMDMixture$y/max(BMDMixture$y)) * (cmax * 1.2) # Density rescaled

  ylim.prim <- c(cmin, cmax)
  ylim.sec <- c(min(BMDMixture$y), max(BMDMixture$y))
  b <- diff(ylim.prim)/diff(ylim.sec)
  a <- ylim.prim[1] - b*ylim.sec[1]

  BMDMixture$yres = a + BMDMixture$y*b

  dplot <- ggplot(data = preds$model_averaged, aes(x = dgr, y = log10(model_averaged),
                                                   group = 1)) +
    geom_line(show.legend = FALSE, linetype = "dashed", size = 3) +
    # geom_area(data = BMDMixture,
    #           aes(x = Dose, y = yres, group = Model, fill = Model),
    #           color = NA, inherit.aes = FALSE, alpha = 0.5) +
    geom_ribbon(data = BMDMixture, aes(x = Dose, y = yres, ymin = min(yres), ymax = yres, fill = Model),
                color=NA, alpha = 0.5) +
    # geom_line(data = preds_min, mapping = aes(x = dgrprime, y = lg10predicted,
    #                                           group = interaction(Model, Distribution)),
    #           linetype = "dashed", alpha = 0.6,
    #           size = 0.8, inherit.aes = FALSE, show.legend = FALSE, color = 1) +

    # geom_line(data = preds_min2, mapping = aes(x = Dose, y = log10MA),
    #           linetype = "dotted", alpha = 1,
    #           size = 3, inherit.aes = FALSE, show.legend = FALSE) +

    geom_segment(data = preds_min2, mapping = aes(x = dgr[1], y = log10MA,
                                                  xend = max(dgr[(dgr <= (lg10d[2]-((lg10d[2]-lg10d[1])/2)))]),
                                                  yend = log10(preds$model_averaged$model_averaged[1]),
                                                  # group = interaction(Model, Distribution),
                                                  # color = Model
    ),
    linetype = "dotted",
    size = 3, inherit.aes = FALSE, show.legend = FALSE) +

    geom_errorbarh(data = dplyr::filter(BMDBMA, Type == weight_type),
                   aes(xmin = log10(BMDL/mod.obj$max.dose), xmax = log10(BMDU/mod.obj$max.dose),
                       y = log10(min(preds$model_averaged$model_averaged)),
                       group = Model),
                   linetype = "solid", show.legend = FALSE,
                   size = 2, height = 0.01, inherit.aes = FALSE,
                   color = RColorBrewer::brewer.pal(9, "Set1")[3]) +
    geom_point(data = dplyr::filter(BMDBMA, Type == weight_type),
               aes(x = log10(BMD/mod.obj$max.dose),
                   y = log10(min(preds$model_averaged$model_averaged)), group = Model),
               size = 5, shape = 19,
               color = RColorBrewer::brewer.pal(9, "Set1")[1],
               show.legend = FALSE,
               inherit.aes = FALSE)  +
    scale_y_continuous(expression(log[10](response)),
                       sec.axis = sec_axis(#~.*1.2,
                         # ~ ((. /max(BMDMixture$y)) *(cmax*1.2)),
                         ~ (. - a)/b,
                         name = "Rescaled Density")
    ) +
    labs(x = expression(log[10](dose))) +
    theme_minimal() +
    theme(strip.text = element_text(size = 15, face = "bold"),
          axis.text = element_text(size = 10, face = "bold"),
          axis.title = element_text(size = 15, face = "bold"),
          legend.text = element_text(size = 7, face = "bold"),
          legend.title = element_text(size = 15, face = "bold"),
          panel.spacing = unit(5, "lines"),
          legend.position = "top",
          legend.direction = "horizontal",
          # title = element_text(size = 15, face = "bold")
    ) +
    scale_fill_manual(values = unname(mods_fills), name = NULL)#, labels = c(""))

  if(include_data == TRUE) {

    pplot2 <- pplot +
      geom_errorbar(data = orig_ptdata,
                    mapping = aes(x = lg10d, ymin = log10m - log10s,
                                  ymax = log10m + log10s),
                    width = NA, position = pd, size = 1,
                    inherit.aes = FALSE) +
      geom_point(data = orig_ptdata, mapping = aes(x = lg10d, y = log10m) ,
                 size = 2, color = 1, shape = 21,
                 fill = RColorBrewer::brewer.pal(9, "Set1")[2],
                 inherit.aes = FALSE)

    dplot2 <- dplot +
      geom_errorbar(data = orig_ptdata,
                    mapping = aes(x = lg10d, ymin = log10m - log10s,
                                  ymax = log10m + log10s),
                    width = NA, position = pd, size = 1,
                    inherit.aes = FALSE) +
      geom_point(data = orig_ptdata, mapping = aes(x = lg10d, y = log10m) ,
                 size = 2, color = 1, shape = 21,
                 fill = RColorBrewer::brewer.pal(9, "Set1")[2],
                 inherit.aes = FALSE)

    pts1 <- ggpubr::ggarrange(pBMDs, pWeights, pplot2,
                             dplot2,
                             nrow = 2, ncol = 2)
    pts2 <- ggpubr::annotate_figure(pts1, top = ggpubr::text_grob(title,
                                                           color = "red", face = "bold", size = 14))

  } else {

    pts1 <- ggpubr::ggarrange(pBMDs, pWeights, pplot, dplot, nrow = 2, ncol = 2)
    pts2 <- ggpubr::annotate_figure(pts1, top = ggpubr::text_grob(title,
                                                          color = "red", face = "bold", size = 14))

  }

  if(all == TRUE){
    return(pts2)
  }else{
    if(include_data == TRUE){
      return(list(BMDs = pBMDs, weights = pWeights, model_fit = pplot2, MA_fit = dplot2))
    }else{
      return(list(BMDs = pBMDs, weights = pWeights, model_fit = pplot, MA_fit = dplot))
    }
  }

}


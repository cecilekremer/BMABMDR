#' Function for internal use
#'
#' @param model_stan stan model
#' @param pars parameters
#'
#' @return .
#'
posterior_diag <- function(model_stan,
                           pars = c('par[2]', letters[c(1, 3, 4)],
                                    "invsigma2", "lp__")) {

  color_scheme_set("darkgray")
  if(class(model_stan)[1] == "CmdStanMCMC") {


    posterior_draws <- as_draws_matrix(model_stan$draws())
    lp_draws <- lapply(pars, function(i) {
      as_draws_matrix(model_stan$draws(i))
    })

    rhats <- unlist(lapply(lp_draws, Rhat))

  } else {

    posterior_draws <- as_draws_matrix(model_stan)
    rhats <- rhat(model_stan)
    neffs <- neff_ratio(model_stan, pars)

  }

  lp_draws <- log_posterior(model_stan)
  neffs <- neff_ratio(model_stan, pars)
  nts_model <- nuts_params(model_stan)

  ## Convergence traceplots
  tplt <- mcmc_trace(posterior_draws, pars = pars) +
    labs(x = 'Iteration', y = 'Draws') +
    theme(axis.text = element_text(size = 10, face = "bold"))
  # print(tplt)

  ## MCMC intervals
  # int_plt1 <- mcmc_intervals(posterior_draws, pars = pars,
  #                   prob = 0.60,
  #                   prob_outer = 0.95) +
  #   theme(axis.text = element_text(size = 10, face = "bold"))
  # print(int_plt1)

  ## MCMC Density
  int_plt2 <- mcmc_areas(posterior_draws, pars = pars,
                         prob = 0.60,
                         prob_outer = 0.95) +
    theme(axis.text = element_text(size = 10, face = "bold"))
  # print(int_plt2)

  ## Univarite plot per chain
  # vio_plt <- mcmc_violin(posterior_draws, pars = pars,
  #             probs = c(0.1, 0.5, 0.9)) +
  #   labs(x = "Chains", y = "Draws") +
  #   theme(axis.text = element_text(size = 10, face = "bold"))+
  #   yaxis_text(hjust = 0)
  # print(vio_plt)


  ## Scatter plot
  plt_pairs <- mcmc_pairs(posterior_draws, np = nts_model,
                          pars = pars,
                          off_diag_args = list(size = 0.75)) #+
  #labs()
  #theme(axis.text = element_text(size = 10, face = "bold"))
  # print(plt_pairs)


  ## Divergence Transitions
  #plt_nuts <- mcmc_nuts_divergence(nts_model, lp = lp_draws) #+
  #theme(axis.text = element_text(size = 10, face = "bold"))
  #print(plt_nuts)

  ########## Traduitional Convergence
  color_scheme_set("brightblue")
  ## Autocorrelation
  autcoor <- mcmc_acf_bar(posterior_draws, pars = pars,
                          lags = 10) +
    theme(axis.text = element_text(size = 10, face = "bold")) +
    ggtitle('Autocorrelation')
  #print(autcoor)

  ## Rhat
  color_scheme_set("brightblue")
  rhat_plot <- mcmc_rhat(rhats) +
    yaxis_text(hjust = 0) +
    geom_vline(yintercept = 1.01, linetype = 'dashed') +
    theme(axis.text = element_text(size = 10, face = "bold")) +
    ggtitle(expression(hat(R)))
  #print(rhat_plot)

  ## N efficiency
  color_scheme_set("brightblue")
  neff_plot <- mcmc_neff(neffs) +
    annotate("text", x = 0.5, y = mean(1:length(pars)),
             parse = TRUE, size = 4,
             label = " frac(N[eff], N) < 0.1 * 'is Problematic'") +
    yaxis_text(hjust = 0) +
    theme(axis.text = element_text(size = 10, face = "bold")) +
    ggtitle(expression(N[eff]))

  #print(neff_plot)
  classic_plts <- ggarrange(autcoor, rhat_plot, neff_plot,
                            ncol = 2, nrow = 2)

  # print(classic_plts)

  #return(NULL)

}

#' Function for internal use
#'
#' @param model_stan stan model
#' @param pars parameters
#'
#' @return .
#'
convergence_deci <- function(model_stan,
                             nrchains = 3,
                             pars = c('par[2]', letters[c(1, 3, 4)],

                                      "invsigma2", "lp__")) {


  #posterior_draws <- as.array(model_stan)
  #lp_draws <- log_posterior(model_stan)
  if("CmdStanMCMC" %in% class(model_stan)) {

    lp_draws <- lapply(pars, function(i) {
      as_draws_matrix(model_stan$draws(i))
    })


    #rhats <- rhats[names(rhats) %in% pars]

  } else {

    lp_draws <- lapply(pars, function(i) {
      extract_variable_matrix(model_stan, i)
    })

    #rhats <- bayesplot::rhat(model_stan)
    #rhats <- rhats[names(rhats) %in% pars]
    #rhats <- rhat(model_stan)
    #rhats <- rhats[names(rhats) %in% pars]


  }

  rhats <- unlist(lapply(lp_draws, Rhat))
  ess_b <- unlist(lapply(lp_draws, ess_bulk))
  ess_t <- unlist(lapply(lp_draws, ess_tail))

  #neffs <- neff_ratio(model_stan)
  #neffs <- neffs[names(neffs) %in% pars]



  #Rhat > 1.01 or 1.05 is problematic
  test_rhat <- (rhats < 1.01)
  test_rhat[is.na(test_rhat)] <- FALSE
  report_rhat <- ifelse(test_rhat == TRUE, 1, 0)

  if(sum(test_rhat) < length(pars)) {

    problem_rhat <- pars[which(test_rhat == FALSE)]
    err_msg_rhat <- paste0("Rhat > 1.01 for the following parameters; ",
                           problem_rhat, "\n")
  } else err_msg_rhat <- ""

  #bulk and tail ESS should be above nchain*100

  test_bulk <- ess_b >= (100 * nrchains)
  test_bulk[is.na(test_bulk)] <- FALSE
  report_bulk <- ifelse(test_bulk == TRUE, 1, 0)

  test_tail <- ess_t >= (100 * nrchains)
  test_tail[is.na(test_tail)] <- FALSE
  report_tail <- ifelse(test_tail == TRUE, 1, 0)

  if(sum(test_bulk) < length(pars) & sum(test_tail) < length(pars)) {

    problem_bulk <- pars[which(test_bulk == FALSE)]
    err_msg_bulk <- paste0("Sampling not efficeintly done for the bulk of the posterior distributions of ",
                           problem_bulk, ".",
                           " Measures of central tendency cannot be trusted.\n")

    problem_tail <- pars[which(test_tail == FALSE)]
    err_msg_tail <- paste0("Sampling not efficeintly done at the tail of the posterior distributions of ",
                           problem_tail, ".",
                           " Quantile estimates cannot be trusted.\n")

    # stop(paste(err_msg_bulk, err_msg_tail, "Increasing the number
    #      of iterations might help.", sep = " "))

  } else if(sum(test_bulk) == length(pars) & sum(test_tail) < length(pars)) {

    problem_tail <- pars[which(test_tail == FALSE)]
    err_msg_tail <- paste0("Sampling not efficeintly done at the tail of the posterior distributions of ",
                           problem_tail, ".",
                           " Quantile estimates cannot be trusted.\n")
    err_msg_bulk <- ""

    # stop(paste(err_msg_tail, "Increasing the number
    #      of iterations might help.", sep = " "))

  } else if(sum(test_bulk) < length(pars) & sum(test_tail) == length(pars)) {

    problem_bulk <- pars[which(test_bulk == FALSE)]
    err_msg_bulk <- paste0("Sampling not efficeintly done for the bulk of the posterior distributions of ",
                           problem_bulk, ".",
                           " Measures of central tendency cannot be trusted.\n")
    err_msg_tail <- ""

    # stop(paste(err_msg_bulk, "Increasing the number
    #      of iterations might help.", sep = " "))

  } else {
    err_msg_bulk <- ""
    err_msg_tail <- ""
  }


  if(nzchar(err_msg_rhat[1]) & nzchar(err_msg_bulk[1]) & nzchar(err_msg_tail[1])) {

    warning(paste0(err_msg_rhat, err_msg_bulk, err_msg_tail,
                   "Increasing the number of iterations might help.\n"))
    msg_label <- "Severe Convergence Issues.\n"

  } else if(nzchar(err_msg_rhat[1]) & nzchar(err_msg_bulk[1]) & !nzchar(err_msg_tail[1])) {

    warning(paste0(err_msg_rhat, err_msg_bulk,
                   "Increasing the number of iterations might help.\n"))
    msg_label <- "Problems with Rhat and Bulk of the posterior distribution.\n"

  } else if(nzchar(err_msg_rhat[1]) & !nzchar(err_msg_bulk[1]) & nzchar(err_msg_tail[1])) {

    warning(paste0(err_msg_rhat, err_msg_tail,
                   "Increasing the number of iterations might help.\n"))
    msg_label <- 'Problems with Rhat and Tail of the posterior distribution.\n'

  } else if(!nzchar(err_msg_rhat[1]) & nzchar(err_msg_bulk[1]) & nzchar(err_msg_tail[1])) {

    warning(paste0(err_msg_bulk, err_msg_tail,
                   "Please check if Rhat < 1.01. Also, check the trace and autocorrelation plot to see if measures of central tendency and quantiles can be trusted.\n"))
    msg_label <- 'Problems with Bulk and Tail of the posterior distribution.\n'

  } else if(nzchar(err_msg_rhat[1]) & !nzchar(err_msg_bulk[1]) & !nzchar(err_msg_tail[1])) {

    warning(paste0(err_msg_rhat,
                   "Please check the trace and autocorrelation plot to see if measures of central tendency and quantiles can be trusted.\n"))
    msg_label <- 'Problems with Rhat.\n'

  } else if(!nzchar(err_msg_rhat[1]) & nzchar(err_msg_bulk[1]) & !nzchar(err_msg_tail[1])) {

    warning(paste0(err_msg_bulk,
                   "Problems with Bulk of the posterior distribution. Please check if Rhat < 1.01.  Also check the trace and autocorrelation plot to see if measures of central tendency can be trusted.\n"))
    #msg_label <- 'Problems with Bulk of the posterior distribution.\n'

  } else if(!nzchar(err_msg_rhat[1]) & !nzchar(err_msg_bulk[1]) & nzchar(err_msg_tail[1])) {

    warning(paste0(err_msg_tail,
                   "Problems with Tail of the posterior distribution. Please check if Rhat < 1.01. Also check the trace and autocorrelation plot to see if quantiles can be trusted.\n"))
    #msg_label <- 'Problems with Tail of the posterior distribution.\n'

  } else {
    msg_label <- 'convergence achieved.\n'
    message('convergence achieved.\n')
  }


  err_table <- data.frame(Parameters = pars,
                          Rhat = rhats,
                          Bulk_ESS = ess_b,
                          Tail_ESS = ess_t,
                          Rhat_Convergence = report_rhat,
                          Bulk_Convergence = report_bulk,
                          Tail_Convergence = report_tail
  )


  return(err_table)


}







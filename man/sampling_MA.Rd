% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fun_Sampling.R
\name{sampling_MA}
\alias{sampling_MA}
\alias{sampling_MAc}
\alias{samplingQ_MA}
\title{Perform model averaging using MCMC methods (Bridge sampling & Partial Laplace)}
\usage{
sampling_MA(
  data.N,
  data.LN,
  prior.weights = rep(1, 16),
  ndraws = 30000,
  nrchains = 3,
  nriterations = 3000,
  warmup = 1000,
  delta = 0.8,
  treedepth = 10,
  seed = 123,
  pvec = c(0.05, 0.5, 0.95)
)

sampling_MAc(
  data.N,
  data.LN,
  prior.weights = rep(1, 16),
  ndraws = 30000,
  nrchains = 3,
  nriterations = 3000,
  warmup = 1000,
  delta = 0.8,
  treedepth = 10,
  seed = 123,
  pvec = c(0.05, 0.5, 0.95)
)

samplingQ_MA(
  data.Q,
  prior.weights = rep(1, 8),
  ndraws = 30000,
  nrchains = 3,
  nriterations = 5000,
  warmup = 1000,
  delta = 0.8,
  treedepth = 10,
  seed = 123,
  pvec = c(0.025, 0.5, 0.975)
)
}
\arguments{
\item{data.N}{the input data as returned by function PREP_DATA_N or PREP_DATA_N_C for clustered data}

\item{data.LN}{the input data as returned by function PREP_DATA_LN or PREP_DATA_LN_C for clustered data}

\item{prior.weights}{a vector specifying which of the 16 (continuous) or 8 (quantal) models should be included (1 = include, 0 = exclude)}

\item{ndraws}{the number of draws, default 30000}

\item{nrchains}{the number of chains to be used in the MCMC, default 3}

\item{nriterations}{the number of iterations per chain, default 3000}

\item{warmup}{the number of iterations per chain to be discarded as burnin, default 1000}

\item{delta}{default 0.8}

\item{treedepth}{default 10}

\item{seed}{random seed, default 123}

\item{pvec}{vector specifying the three BMD quantiles of interest (default 90\% CrI)}

\item{data.Q}{the input data as returned by function PREP_DATA_Q}
}
\value{
The function returns a BMDBMA model object, which is a list containing the following objects:

\code{E4_N} parameter estimates from the (Exponential, Normal) model

\code{IE4_N} parameter estimates from the (Inverse Exponential, Normal) model

\code{H4_N} parameter estimates from the (Hill, Normal) model

\code{LN4_N} parameter estimates from the (Log-normal, Normal)

\code{G4_N} parameter estimates from the (Gamma, Normal) model

\code{QE4_N} parameter estimates from the (Quadratic Exponential, Normal) model

\code{P4_N} parameter estimates from the (Probit, Normal) model

\code{L4_N} parameter estimates from the (Logit, Normal) model

\code{E4_LN} parameter estimates from the (Exponential, Log-normal) model

\code{IE4_LN} parameter estimates from the (Inverse Exponential, Log-normal) model

\code{H4_LN} parameter estimates from the (Hill, Log-normal) model

\code{LN4_LN} parameter estimates from the (Log-normal, Log-normal)

\code{G4_LN} parameter estimates from the (Gamma, Log-normal) model

\code{QE4_LN} parameter estimates from the (Quadratic Exponential, Log-normal) model

\code{P4_LN} parameter estimates from the (Probit, Log-normal) model

\code{L4_LN} parameter estimates from the (Logit, Log-normal) model

\code{E4_Q} parameter estimates from the (Exponential, Quantal) model

\code{IE4_Q} parameter estimates from the (Inverse Exponential, Quantal) model

\code{H4_Q} parameter estimates from the (Hill, Quantal) model

\code{LN4_Q} parameter estimates from the (Lognormal, Quantal) model

\code{G4_Q} parameter estimates from the (Gamma, Quantal) model

\code{QE4_Q} parameter estimates from the (Quadratic Exponential, Quantal) model

\code{P4_Q} parameter estimates from the (Probit, Quantal) model

\code{L4_Q} parameter estimates from the (Logit, Quantal) model

\code{MA_bridge_sampling}  Model averaged BMD credible interval based on Bridge sampling

\code{MA_laplace} Model averaged BMD credible interval based on Hybrid Laplace

\code{weights_bridge_sampling} Model weights used in the averaging for Bridge sampling

\code{weights_laplace} Model weights used in the averaging for Hybrid Laplace

\code{convergence} vector indicating whether the models have converged (1) or not (0)

\code{divergences} vector containing the proportion of divergent transitions for each model

\code{bs_weights_conv} model weights based on bridge sampling including converged models only

\code{ls_weights_conv} model weights based on hybrid laplace including converged models only

\code{MA_bs_conv} model averaged BMD credible interval based on bridge sampling including converged models only

\code{MA_ls_conv} model averaged BMD credible interval based on hybrid laplace including converged models only

\code{bf} Bayes factor comparing the best model against saturated ANOVA model

\code{covs} matrix with covariances between parameters b-d and BMD-d

\code{corrs} matrix with correlation between parameters b-d and BMD-d

\code{p.msg} warning message if model averaged posterior has been truncated

\code{w.msg} warning message if Laplace weights could not be computed and one model gets all the weight

\code{shift} shift value for lognormal data

\code{BIC.SM} BIC value of saturated model, used to test for GOF

\code{BIC.bestfit} BIC value of best fitting model, used to test for GOF

\code{means.SM} mean response per dose level estimated from the saturated model, used to test for GOF

\code{gof_check} GOF message

\code{models_included_laplace} vector containing the names of models included in the model averaging using Laplace

\code{models_included_bridge} vector containing the names of models included in the model averaging using Bridge sampling

\code{q} BMR

\code{max.dose} maximum dose level (original scale)

\code{dataN} normal summary data used for analysis

\code{dataLN} lognormal summary data used for analysis

\code{BMDMixture} vector of length \code{ndraws} containing the draws from the model-averaged posterior based on hybrid laplace

\code{BMDMixture.conv} vector of length \code{ndraws} containing the draws from the model-averaged posterior based on hybrid laplace, including converged models only

\code{BMDMixtureBS} vector of length \code{ndraws} containing the draws from the model-averaged posterior based on bridge sampling

\code{BMDMixture.convBS} vector of length \code{ndraws} containing the draws from the model-averaged posterior based on bridge sampling, including converged models only

\code{MA_dr_bs} vector containing the model-averaged response at each dose level, for bridge sampling

\code{MA_dr_ls} vector containing the model-averaged response at each dose level, for hybrid laplace

\code{MA_dr_bs_conv} vector containing the model-averaged response at each dose level, for bridge sampling including converged models only

\code{MA_dr_ls_conv} vector containing the model-averaged response at each dose level, for hybrid laplace including converged models only

\code{MA_post_bs} vector containing the 0.5\%-percentiles of the model-averaged posterior, for bridge sampling

\code{MA_post_ls} vector containing the 0.5\%-percentiles of the model-averaged posterior, for hybrid laplace

\code{MA_post_bs_conv} vector containing the 0.5\%-percentiles of the model-averaged posterior, for bridge sampling including converged models only

\code{MA_post_ls_conv} vector containing the 0.5\%-percentiles of the model-averaged posterior, for hybrid laplace including converged models only

\code{llN} vector containing the loglikelihood values of the Normal models

\code{llLN} vector containing the loglikelihood values of the Lognormal models

\code{parsN} list containing the fitted Normal models

\code{parsLN} list containing the fitted Lognormal models

\code{is_bin} logical indicating whether binomial (no litter) model was used

\code{is_betabin} logical indicating whether betabinomial (with litter) model was used

\code{data} quantal summary data used for analysis

\code{parsQ} list containing the fitted Quantal models

\code{llQ} vector containing the loglikelihood values of the Quantal models
}
\description{
Using MCMC, we compute the parameters of each model, perform model averaging using bridge sampling.
We also implemented a Hybrid Laplace method within this function where the model parameters are estimated
with MCMC, but the model weights are computed using Laplace approximation.
By default, all 16 models are included for continuous data, and all 8 models for quantal data.
Models can be excluded by setting their respective weight to 0 in \code{prior.weights}. The order of models fitted can be obtained using the \code{get_models()} function.

\code{sampling_MA} is used for continuous data

\code{sampling_MAc} is used for clustered continuous data

\code{samplingQ_MA} is used for quantal data
}
\details{
This method assumed data for continuous endpoints.

More detailed descriprion
}
\examples{
 data_N <- PREP_DATA_N(data = as.data.frame(immunotoxicityData[1:5,]),
                       sumstats = TRUE, sd = TRUE, q = 0.1)
 data_LN <- PREP_DATA_LN(data = as.data.frame(immunotoxicityData[1:5,]),
                         sumstats = TRUE, sd = TRUE, q = 0.1) #'
 SBMD <- sampling_MA(data_N, data_LN, prior.weights = c(rep(1,4), rep(0,12)))

}

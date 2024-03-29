% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diagnostics.R
\name{convergence_deciC}
\alias{convergence_deciC}
\title{Function for internal use; decides if model has converged based on Rhat, bulk and tail ESS}
\usage{
convergence_deciC(
  model_stan,
  nrchains = 3,
  pars = c("par1", "par2", "par3", "d", "invsigma2", "rho_cluster", "lp__")
)
}
\arguments{
\item{model_stan}{stan model}

\item{nrchains}{default 3; number of MCMC chains}

\item{pars}{parameters}
}
\value{
a table with convergence diagnostics
}
\description{
Function for internal use; decides if model has converged based on Rhat, bulk and tail ESS
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fct_optim.R
\name{fun_optim}
\alias{fun_optim}
\alias{fun_optimC}
\alias{fun_optimQ}
\alias{fun_optimCov}
\alias{fun_optimQCov}
\title{Function to perform Laplace approximation for a given dose-response model (used internally).}
\usage{
fun_optim(mod, data, stv, ndraws, seed, pvec)

fun_optimC(mod, data, stv, ndraws, seed, pvec)

fun_optimQ(mod, data, stv, ndraws, seed, pvec)

fun_optimCov(mod, data, stv, ndraws, seed, pvec)

fun_optimQCov(mod, data, stv, ndraws, seed, pvec)
}
\arguments{
\item{mod}{stan model}

\item{data}{list containing the data to be passed to the stan model}

\item{stv}{list of starting values}

\item{ndraws}{number of draws to be made from the posterior distribution}

\item{seed}{random seed for reproducibility}

\item{pvec}{probability vector to compute credible interval for the BMD}
}
\value{
a stan results object, which is a list containing point estimates, value of the log-posterior, value of the return code from the optimizer, Hessian matrix, matrix of parameter draws
}
\description{
Function to perform Laplace approximation for a given dose-response model (used internally).
}
\examples{
\dontrun{
fun_optim(stanmodels$mE4, data, start, 30000, 123, c(0.05,0.5,0.95))
}

}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fun_summData.R
\name{summarize.indiv.data}
\alias{summarize.indiv.data}
\title{Function to summarize individual data}
\usage{
summarize.indiv.data(
  data,
  type = c("continuous", "quantal"),
  cluster = FALSE,
  covar = FALSE
)
}
\arguments{
\item{data}{data in order of: dose, response, covariate OR dose, response, litter}

\item{type}{either 'continuous' or 'quantal'}

\item{cluster}{logical indicating whether the data include a litter effect (defaults to FALSE)}

\item{covar}{logical indicating whether the data include a covariate (defaults to FALSE)}
}
\value{
dataframe of summary data
}
\description{
Function to summarize individual data
}
\examples{
summarize.indiv.data(das5$data[,c(1,11)], type = 'continuous', cluster = FALSE, covar = FALSE)
summarize.indiv.data(das5$data[,c(1,11,2)], type = 'continuous', cluster = FALSE, covar = TRUE)
summarize.indiv.data(das5$data[,c(1,11,10)], type = 'continuous', cluster = TRUE, covar = FALSE)
\dontrun{
# This gives error because covariates not implemented for clustered data
summarize.indiv.data(das5$data[,c(1,11,10,2)], type = 'continuous', cluster = TRUE, covar = TRUE)
}

}

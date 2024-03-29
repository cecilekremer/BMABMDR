% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FUNs.R
\name{flat}
\alias{flat}
\title{Function for internal use to determine if DR curve flattens out for continuous data}
\usage{
flat(dose, mean, n, inc)
}
\arguments{
\item{dose}{ordered dose levels}

\item{mean}{mean response per ordered dose level}

\item{n}{value}

\item{inc}{logical variable to indicate if the dose-resonse curve is increasing or decreasing}
}
\value{
logical value indicating if the dose-response curve is flat or not
}
\description{
Function for internal use to determine if DR curve flattens out for continuous data
}
\examples{
flat(immunotoxicityData$Dose[1:5], immunotoxicityData$Mean[1:5], immunotoxicityData$n[1:5], inc = TRUE)

}

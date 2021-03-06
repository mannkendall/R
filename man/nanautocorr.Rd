% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nanautocorr.R
\name{nanautocorr}
\alias{nanautocorr}
\title{Autocorrelation function}
\usage{
nanautocorr(data, nlags, R)
}
\arguments{
\item{data}{vector of time series}

\item{nlags}{maximum number of lags in the computation of ACF}

\item{R}{number of lags until the model is  supposed to have a significant AC coefficient (optional)}
}
\value{
out vector of the ACF
}
\description{
`nanautocorr` computes autocorrelation coefficient dealing with missing values
}
\details{
Calculates the `nlag` autocorrelation coefficient for a data vector containing missing values
couples of data including `NA` are excluded from the computation.
Here the ACF is calculated using the Pearson\'s correlation coefficient for each lag.
}
\author{
Fabio Oriani, fabio.oriani@unine.ch, University of Neuchatel (CH) original author; Alessandro Bigi (abigi@unimore.it), University of Modena and Reggio Emilia (IT) for the R version
}

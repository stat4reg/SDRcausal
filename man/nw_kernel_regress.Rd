% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nw_kernel_regress.R
\name{nw_kernel_regress}
\alias{nw_kernel_regress}
\title{The Nadaraya-Watson kernel estimator}
\usage{
nw_kernel_regress(
  y,
  x,
  bandwidth = 1,
  kernel = "EPAN",
  gauss_cutoff = 0.001,
  verbose = FALSE
)
}
\arguments{
\item{y}{Y in E(Y|X)}

\item{x}{X in E(Y|X)}

\item{bandwidth}{Kernel bandwidth}

\item{kernel}{Indicates which kernel function to be used}

\item{gauss_cutoff}{Cutoff value for Gaussian kernel}

\item{verbose}{Specifies if the program should print output while
running.}
}
\value{
Value of kernel regression
}
\description{
Gives the expected value of Y given X = x by kernel regression
             according to the Nadaraya-Watson kernel estimator to get
             E(Y|X). Note that y and x may be vectors or matrices, as long
             as dim(x)[1] == dim(y)[1].
}
\examples{
# Using example data from package SDRcausal
library(SDRcausal)

# Import example data
x <- SDRcausal::covariates
y <- SDRcausal::outcomes

# Extimating y given x, E(y | x)
k <- nw_kernel_regress(y, x, bandwidth = 1)

}
\references{
Ghosh, T., Ma, Y., & De Luna, X. (2020). Sufficient dimension
reduction for feasible and robust estimation of average causal effect.
Statistica Sinica, accepted.
}

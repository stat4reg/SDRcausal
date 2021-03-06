% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inf.ate.R
\name{inf.ate}
\alias{inf.ate}
\title{Performs Estimations of Average Treatment Effect and Infrences}
\usage{
inf.ate(
  x,
  y,
  treated,
  beta_guess1,
  beta_guess0,
  imp.solver = "optim",
  imp.kernel = "EPAN",
  imp.explicit_bandwidth = FALSE,
  imp.recalc_bandwidth = TRUE,
  bwc_dim_red1 = 1,
  bwc_impute1 = 1.25,
  bwc_dim_red0 = 1,
  bwc_impute0 = 1.25,
  imp.gauss_cutoff = 0.001,
  imp.penalty = 10,
  imp.n_before_pen = 5,
  imp.to_extrapolate = TRUE,
  imp.to_truncate = TRUE,
  imp.extrapolation_basis = 5,
  alpha_initial,
  ipw.solver = "optim",
  ipw.kernel = "EPAN",
  ipw.explicit_bandwidth = FALSE,
  ipw.recalc_bandwidth = TRUE,
  bwc_dim_red = 1,
  bwc_prop_score = 10,
  ipw.gauss_cutoff = 0.001,
  ipw.penalty = 10,
  ipw.n_before_pen = 1,
  n_threads = 1,
  verbose = TRUE,
  imp.solver.options = NA,
  ipw.solver.options = NA
)
}
\arguments{
\item{x}{Covariate matrix}

\item{y}{Response vector}

\item{treated}{A binary vector indicating treatment status}

\item{beta_guess1}{Initial guess for \eqn{\beta_1}}

\item{beta_guess0}{Initial guess for \eqn{\beta_0}}

\item{imp.solver}{Specifies which solver is to be used. Current options are \code{optim} and \code{cobyla} (from \code{nloptr} package). The diffault value is \code{"optim"}.}

\item{imp.kernel}{Specifies which kernel function is to be used, current options are: \code{"EPAN"}, \code{"QUARTIC"}, and \code{"GAUSSIAN"}. The default value is \code{"EPAN"}.}

\item{imp.explicit_bandwidth}{Specifies if \code{bandwidth_scale} will be used as the bandwidth or if it will be calculated as \code{bandwidth_scale} * sd(\eqn{\beta^T x}) * \eqn{n^{(1/5)}}. The default value is \code{FALSE}.}

\item{imp.recalc_bandwidth}{Specifies whether the bandwidth should be recalculated after the first stage (the estimations of dimension reduction step). If \code{explicit_bandwidth} is \code{TRUE}, \code{recalc_bandwidth} is not used, but if \code{explicit_bandwidth} is \code{FALSE}, then if \code{recalc_bandwidth} is \code{TRUE}, bandwidths are recalculated at the beginning of the second step based on \code{bwc_impute0} and \code{bwc_impute1}. If \code{recalc_bandwidth} is \code{FALSE}, the first step bandwidths are used. The default value is \code{FALSE}.}

\item{bwc_dim_red1}{Scaling of calculated bandwidth, or if \code{explicit_bandwidth = TRUE} used as the bandwidth. It is used in the dimension reduction step for \eqn{\hat{m}_1(\beta_1^T x)}. The default value is \code{1}.}

\item{bwc_impute1}{Scaling of calculated bandwidth, or if \code{explicit_bandwidth = TRUE} used as the bandwidth. It is used in the imputation step for \eqn{\hat{m}_1(\beta_1^T x)}. The default  value is \code{1.25}.}

\item{bwc_dim_red0}{Scaling of calculated bandwidth, or if \code{explicit_bandwidth = TRUE} used as the bandwidth. It is used in the dimension reduction step for \eqn{\hat{m}_0(\beta_0^T x)}. The default value is \code{1}.}

\item{bwc_impute0}{Scaling of calculated bandwidth, or if \code{explicit_bandwidth = TRUE} used as the bandwidth. It is used in the imputation step for \eqn{\hat{m}_0(\beta_0^T x)}. The default value is \code{1.25}.}

\item{imp.gauss_cutoff}{The cutoff value for Gaussian kernel. The default value is \code{1e-3}.}

\item{imp.penalty}{Penalty for the optimizer if local linear regression fails. Added to the function value in solver as \code{penalty}^(n - \code{n_before_pen}), where n is the number of times local linear regression fails. The default value is \code{10}.}

\item{imp.n_before_pen}{The number of acceptable local linear regression failures during the estimation of \eqn{\beta_0} and \eqn{\beta_1} phase. The default value is \code{5}.}

\item{imp.to_extrapolate}{Specifies whether to extrapolate or not. Since in \eqn{\hat{m}_0(\beta_0^T x)} and \eqn{\hat{m}_1(\beta_1^T x)} estimates in terms of \eqn{\beta_0} and \eqn{\beta_1}, local linear regression at the boundaries of \eqn{\beta_0^Tx} and \eqn{\beta_1^Tx} can be very volatile, it is recommended to use extrapolation on those points instead of local linear regression. The default value is \code{TRUE}.}

\item{imp.to_truncate}{Specifies whether to truncate \eqn{\hat{m}_0(\beta_0^T x)} and \eqn{\hat{m}_1(\beta_1^T x)} or not. After estimating \eqn{\hat{m}_0(\beta_0^T x)} and \eqn{\hat{m}_1(\beta_1^T x)}, if they are outside the range of observed outputs, they are replaced with the minimum and maximum observed outputs. The default value is \code{TRUE}.}

\item{imp.extrapolation_basis}{The number of data points to base extrapolation on. Extrapolation at border points can be done based on a different number of neighborhood points. \code{extrapolation_basis} is how many neighborhood points are used. The default value is \code{5}.}

\item{alpha_initial}{Initial guess for \eqn{\alpha}}

\item{ipw.solver}{Specifies which solver is to be used. Current options are \code{optim} and \code{cobyla} (from \code{nloptr} package). The diffault value is \code{"optim"}.}

\item{ipw.kernel}{Specifies which kernel function is to be used, current options are: \code{"EPAN"}, \code{"QUARTIC"}, and \code{"GAUSSIAN"}. The default value is \code{"EPAN"}.}

\item{ipw.explicit_bandwidth}{Specifies if \code{bandwidth_scale} will be used as the bandwidth or if it will be calculated as \code{bandwidth_scale} * sd(\eqn{\alpha^T x}) * \eqn{n^{(1/5)}}. The default value is \code{FALSE}.}

\item{ipw.recalc_bandwidth}{Specifies whether the bandwidth should be recalculated after the estimations of \eqn{\alpha}. If \code{explicit_bandwidth} is \code{TRUE}, \code{recalc_bandwidth} is not used, but if \code{explicit_bandwidth} is \code{FALSE}, then if \code{recalc_bandwidth} is \code{TRUE}, bandwidths are recalculated at the beginning of the second step based on \code{bwc_prop_score}. If \code{recalc_bandwidth} is \code{FALSE}, the first step bandwidths are used. The default value is \code{FALSE}.}

\item{bwc_dim_red}{Scaling of calculated bandwidth, or if \code{explicit_bandwidth = TRUE} used as the bandwidth. It is used in the dimension reduction step for \eqn{\alpha^T x}. The default value is \code{1}.}

\item{bwc_prop_score}{Scaling of calculated bandwidth, or if \code{explicit_bandwidth = TRUE} used as the bandwidth. It is used for the estimation of the propensity score. The default value is \code{10}.}

\item{ipw.gauss_cutoff}{The cutoff value for Gaussian kernel. The default value is \code{1e-3}.}

\item{ipw.penalty}{Penalty for the optimizer if a probability is outside (0, 1) during the estimation of \eqn{\alpha} phase. Added to the function value in solver as \code{penalty}^(n - \code{n_before_pen}), where n is the number of probabilities outside (0, 1). The default value is \code{10}.}

\item{ipw.n_before_pen}{The number of probabilities outside the range (0, 1) to accept during the estimation of \eqn{\alpha} phase. The default value is \code{1}.}

\item{n_threads}{Sets the number of threads for parallel computing. Set to 1 serial. If \code{n_threads} exceeds the maximum number of threads, sets \code{n_threads} to max_threads - 1. To use max_threads, set to \code{n_threads} to max_threads of system. The default value is \code{1}.}

\item{verbose}{Specifies if the program should print output while running. The default value is \code{TRUE}.}

\item{imp.solver.options}{Additional parameters passed to \code{optim} or \code{cobyla} for \code{imp.ate}.}

\item{ipw.solver.options}{Additional parameters passed to \code{optim} or \code{cobyla} for \code{ipw.ate}.}
}
\value{
A list containing the average treatment effect of the
        combination of observed and imputed values (ate), the average
        treatment effect based on the imputed values only (ate2), the
        imputed values for treated (m1) and untreated treated (m0), the and
        the output from optim (op).
}
\description{
Semiparametric estimation of the average treatment effect based
             on the all methods described in Ghosh, Ma, & De Luna
             (2020) and all infrences.
}
\examples{
# Using example data from package SDRcausal
library(SDRcausal)

# Import example data
x <- SDRcausal::covariates
y <- SDRcausal::outcomes
trt <- SDRcausal::treated
b1 <- SDRcausal::beta1_guess
b0 <- SDRcausal::beta0_guess
a <- SDRcausal::alpha_guess

# Perform semiparametric imputation
inf.ate <- SDRcausal::inf.ate(x, y, trt, b1, b0, alpha_initial = a)

}
\references{
Ghosh, T., Ma, Y., & De Luna, X. (2020). Sufficient dimension
reduction for feasible and robust estimation of average causal effect.
Statistica Sinica, accepted.
}
\seealso{
[stats::optim]
}

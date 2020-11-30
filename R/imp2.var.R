#' @title Estimates IMP2 variance
#'
#' @description Variance of IMP2 as in Ghosh, Ma, & De Luna (2020).
#'
#' @param x                   Covariate matrix
#' @param y                   Response vector
#' @param treated             A binary vetor indicating treatment
#' @param imp                 imp output object from \code{imp.ate}
#' @param ipw                 ipw output object from \code{ipw.ate}
#' @param bandwidth_scale1    Scaling of the calculated bandwidth, or in case of \code{explicit_bandwidth = TRUE}, the actual bandwidth for the estimation of \eqn{E(.|\beta_1^T X)}. The default value is \code{imp$bw1}. If this default value is used, one should use the default value \code{TRUE} for \code{explicit_bandwidth}.
#' @param bandwidth_scale0    Scaling of the calculated bandwidth, or in case of \code{explicit_bandwidth = TRUE}, the actual bandwidth for the estimation of \eqn{E(.|\beta_0^T X)}. The default value is \code{imp$bw0}. If this default value is used, one should use the default value \code{TRUE} for \code{explicit_bandwidth}.
#' @param kernel              Specifies which kernel function to be used. The default is \code{"EPAN"}.
#' @param explicit_bandwidth  Specifies if bandwidth_scale will be used as the bandwidth or if it will be calculated as \code{bandwidth_scale} * sd(\eqn{\beta_t^T x}) * \eqn{n^{(1/5)}}. The default value is \code{TRUE}.
#' @param gauss_cutoff        The cutoff value for Gaussian kernel. The default value is \code{1e-3}.
#'
#' @return Variance of IMP2
#'
#' @importFrom stats sd
#'
#' @export
#'
#' @references Ghosh, T., Ma, Y., & De Luna, X. (2020). Sufficient dimension
#' reduction for feasible and robust estimation of average causal effect.
#' Statistica Sinica, accepted.
#'
#' @examples
#' # Using example data from package SDRcausal
#' library(SDRcausal)
#'
#' # Import example data
#' x <- SDRcausal::covariates
#' y <- SDRcausal::outcomes
#' trt <- SDRcausal::treated
#' b1 <- SDRcausal::beta1_guess
#' b0 <- SDRcausal::beta0_guess
#' alp <- SDRcausal::alpha_guess
#'
#' # Perform semiparametric imputation
#' imp <- SDRcausal::imp.ate(x, y, trt, b1, b0,
#'            explicit_bandwidth = TRUE, bwc_dim_red1 = 1, bwc_impute1 = 1,
#'            bwc_dim_red0 = 1, bwc_impute0 = 1)
#'
#' # Perform semiparametric inverse probability weighting
#' ipw <- SDRcausal::ipw.ate(x, y, trt, alp, bwc_dim_red = 10,
#'            bwc_prop_score = 18)
#'
#' # Calculate the variance of the IMP2 estimator.
#' var <- SDRcausal::imp2.var(x, y, trt, imp, ipw,
#'            bandwidth_scale1 = imp$bw1, bandwidth_scale0 = imp$bw0)
#'
imp2.var <- function(x,
                          y,
                          treated,
                          imp,
                          ipw,
                          bandwidth_scale1 = imp$bw1,
                          bandwidth_scale0 = imp$bw0,
                          kernel = "EPAN",
                          explicit_bandwidth = TRUE,
                          gauss_cutoff = 1e-3)
{


  # Number of observations
  n <- as.integer(dim(x)[1])
  n_ones <- rep(1, times = n)
  p <- as.integer(dim(x)[2])


  # Boolean treatement vector
  tbl <- as.logical(treated)


  # Imputation input
  beta1 <- imp$beta1_hat
  m1 <- imp$m1$m
  dm1 <- imp$m1$dm
  xb1 <- x %*% beta1
  d <- as.integer(dim(dm1)[2])
  # Lower p-d x matrix
  x_lower <- x[,(d+1):p]


  beta0 <- imp$beta0_hat
  m0 <- imp$m0$m
  dm0 <- imp$m0$dm
  xb0 <- x %*% beta0

# Deriving parameters from input
  # Checking if explicit bandwidth
  if (explicit_bandwidth) {
    # Setting explicit bandwidths
    bw1 <- bandwidth_scale1
    bw0 <- bandwidth_scale0
  } else {
    # Calculating bandwidths
    sd_xb1 <- sd(xb1[as.logical(treated)])
    bw1 <- bandwidth_scale1 * sd_xb1 * sum(treated)**(-1/5)

    sd_xb0 <- sd(xb0[as.logical(!treated)])
    bw0 <- bandwidth_scale0 * sd_xb0 * sum(!treated)**(-1/5)

    # Setting explicit_bandwidth to TRUE so bw1 and bw0 is used
  }


  # IPW input
  pr <- ipw$pr

  # Calculating naive estimators
  #e1 <- rep(sum(y[tbl]) / sum(treated), n)
  #e0 <- rep(sum(y[!tbl]) / (n - sum(treated)), n)
  e1 = mean(m1)
  e0 = mean(m0)


  # Term 1
  term1 <- m1 - m0 - (e1 - e0)

  # Term 2
  k <- nw_kernel_regress(n_ones / pr, x %*% beta1, bandwidth = bw1)
  term2 <- k * treated * (y - m1)

  # Term 3
  k <- nw_kernel_regress(n_ones / (n_ones - pr), x %*% beta0, bandwidth = bw0)
  term3 <- k * (n_ones - treated) * (y - m0)

  # Term 4
  b1 <- b10_fun(x = x,
                treated = treated,
                dm = dm1,
                beta = beta1,
                kernel = "EPAN",
                bandwidth = bw1,
                gauss_cutoff = gauss_cutoff)

  xc1 <- x_lower - nw_kernel_regress(x_lower, x %*% beta1)




  part1 <- ( t(rep(1, d)) %x% x_lower * (dm1 %x% t(rep(1, p - d)) )/ n)
  part2 <- part1 %*% b1
  part3 <- sweep(part2, MARGIN = 1, treated * (y - m1), "*")
  part4 <- (t(rep(1, d)) %x% xc1 ) * (dm1 %x% t(rep(1, p - d)))

  # The matrix multiplication and subsequent summation over the columns
  # represents the scalar product (1d case) of each observation (i)
  term4 <- colSums(part3 %*% t(part4))

  # Term 5
  b0 <- b10_fun(x = x,
                treated = (n_ones - treated),
                dm = dm0,
                beta = beta0,
                kernel = kernel,
                bandwidth = bw0,
                gauss_cutoff = gauss_cutoff)

  xc0 <- x_lower - nw_kernel_regress(x_lower, x %*% beta0)

  part1 <- ( t(rep(1, d)) %x% x_lower * (dm0 %x% t(rep(1, p - d)) )/ n)
  part2 <- part1 %*% b0
  part3 <- sweep(part2, MARGIN = 1, (n_ones - treated) * (y - m0), "*")
  part4 <- (t(rep(1, d)) %x% xc0 ) * (dm0 %x% t(rep(1, p - d)))

  # The elementwise multiplication and subsequent summation over the columns
  # represents the scalar product (1d case) of each observation (i)
  term5 <- colSums(part3 %*% t(part4))


  output <- mean((term1 + term2 - term3 - term4 + term5)**2)/n

  return(output)
}

#' @title Estimates IPW variance
#'
#' @description Variance of the IPW as in Ghosh, Ma, & De Luna (2020).
#'
#' @param x                  Covariate matrix
#' @param y                  Response vector
#' @param treated            Binary vetor indicating treatment
#' @param imp                imp_output object from imp.ate()
#' @param ipw                ipw_output object from ipw.ate()
#' @param treated            Binary vetor indicating treatment
#' @param bandwidth_scale    Scaling of the calculated bandwidth, or in case of
#' @param kernel             Specifies which kernel function to be used
#' @param explicit_bandwidth Specifies if bandwidth_scale will be used as the
#'                           bandwidth or if it will be calculated as bw =
#'                           bandwidth_scale * sd(x * beta) * n^(1/3).
#' @param gauss_cutoff       Cutoff value for Gaussian kernel
#' @param num_deriv_h        Step size of numerical derivative.
#' @param verbose            Specifies if the program should print output
#'                           while running.
#'
#' @return The variance of IPW
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
#' # Calculate the variance of the Augmented IPW (AIPW)
#' var <- SDRcausal::ipw.var(x, y, trt, imp, ipw,
#'            bandwidth_scale = ipw$bw_pr)
#'
ipw.var <- function(x,
                         y,
                         treated,
                         imp,
                         ipw,
                         bandwidth_scale,
                         kernel = "EPAN",
                         explicit_bandwidth = TRUE,
                         gauss_cutoff = 1e-3,
                         num_deriv_h = 1e-3,
                         verbose = FALSE)
{
  # Deriving parameters from input
  # Number of observations
  n <- as.integer(dim(x)[1])
  n_ones <- rep(1, times = n)
  p <- as.integer(dim(x)[2])

  alpha_hat <- ipw$alpha_hat
  # Extracting the lower p-d matrix of x
  d <- as.integer(dim(alpha_hat)[2])
  x_lower <- x[,(d+1):p]

  # Boolean treatement vector
  tbl <- as.logical(treated)

  # Imputation input
  m1 <- imp$m1$m
  m0 <- imp$m0$m

  # IPW input

  pr <- ipw$pr
  d_pr <- ipw$d_pr
  eta <- log(pr / (1 - pr))
  xa <- x %*% alpha_hat
  d_eta <- d_pr / (pr * (1 - pr)) # Makes sense
  bw_pr <- ipw$bw_pr
  # Calculating bandwidth
  if (explicit_bandwidth) {
    # Setting explicit bandwidths
    bw <- bandwidth_scale
  } else {
    # Calculating bandwidths
    sd_xa <- sd(xa)
    bw <- bandwidth_scale * sd_xa * n**(-1/5)
  }

  e1 = mean(treated * y / pr)
  e0 = mean((n_ones - treated) * y / (n_ones - pr))
  # Calculating the terms that make out the variance
  # Term 1
  term1 <-(treated * y / pr -
           e1 -
           (n_ones - treated) * y / (n_ones - pr) +
           e0)


  # Term 2

  k <- nw_kernel_regress(m1, xa, bandwidth = bw)
  term2 <- (n_ones - (treated / pr)) * k

  # Term 3
  k <- nw_kernel_regress(m0, xa, bandwidth = bw)
  term3 <- ((treated - pr) / (n_ones - pr)) * k

  # Term 4
  b <- b_fun(x,
             treated,
             alpha_hat,
             num_deriv_h,
             kernel,
             bw,
			 bw_pr,
             verbose)

  xc <- x_lower - nw_kernel_regress(x_lower, xa, bandwidth = bw)

  part1 <- colSums(   (t(rep(1, d)) %x%  x_lower) * ( ((m1 * (n_ones - pr ) + m0 * pr)* d_eta)    %x%   t(rep(1, p-d)) ) / n)


  part2 <- part1 %*% b
  part3 <- (t(rep(1, d)) %x% xc) * ( ((treated - pr) * d_eta) %x%   t(rep(1, p - d)) )

  term4 <- rowSums(sweep(part3, MARGIN = 2, part2, "*"))

  # Adding the terms and calculating the variance
  var <- mean((term1 + term2 - term3 + term4)**2)/n

  return(var)
}

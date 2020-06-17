#' @title Estimates Augmented Inverse Probability variance
#'
#' @description Variance of the Augmented IPW as in Ghosh, Ma, & De Luna
#'              (2020).
#'
#' @param x                   Covariate matrix
#' @param y                   Response vector
#' @param treated             Binary vetor indicating treatment
#' @param imp                 imp_output object from semipar_imputation()
#' @param ipw                 ipw_output object from semipar_ipw()
#' @param treated             Binary vetor indicating treatment
#' @param bandwidth_scale1    Scaling of the calculated bandwidth, m1
#' @param bandwidth_scale0    Scaling of the calculated bandwidth, m0
#' @param bandwidth_scale_pr  Scaling of the calculated bandwidth, pr
#' @param kernel              Specifies which kernel function to be used
#' @param explicit_bandwidth  Specifies if bandwidth_scale will be used as the
#'                            bandwidth or if it will be calculated as bw =
#'                            bandwidth_scale * sd(x * beta) * n^(1/5).
#' @param gauss_cutoff        Cutoff value for Gaussian kernel
#' @param num_deriv_h         Step size of numerical derivative.
#' @param verbose             Specifies if the program should print output
#'                            while running.
#'
#' @return The variance of Augmented IPW
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
#' imp <- SDRcausal::semipar_imputation(x, y, trt, b1, b0,
#'            explicit_bandwidth = TRUE, bwc_dim_red1 = 1, bwc_impute1 = 1,
#'            bwc_dim_red0 = 1, bwc_impute0 = 1)
#'
#' # Perform semiparametric inverse probability weighting
#' ipw <- SDRcausal::semipar_ipw(x, y, trt, alp, bwc_dim_red = 10,
#'            bwc_prop_score = 18)
#'
#' # Calculate the variance of the Augmented IPW (AIPW)
#' var <- SDRcausal::aipw_variance(x, y, trt, imp, ipw,
#'            bandwidth_scale1 = imp$bw1, bandwidth_scale0 = imp$bw0,
#'            bandwidth_scale_pr = ipw$bw_pr)
#'
aipw_variance <- function(x,
                          y,
                          treated,
                          imp,
                          ipw,
                          bandwidth_scale1,
                          bandwidth_scale0,
                          bandwidth_scale_pr,
                          kernel = "EPAN",
                          explicit_bandwidth = TRUE,
                          gauss_cutoff = 1e-3,
                          num_deriv_h = 1e-8,
                          verbose = FALSE)
{
  # Deriving parameters from input
  # Number of observations
  n <- as.integer(dim(x)[1])
  n_ones <- rep(1, times = n)
  p <- as.integer(dim(x)[2])
  d <- as.integer(dim(imp$beta0_hat)[2])


  # Boolean treatement vector
  tbl <- as.logical(treated)

  # Lower p-d x matrix
  x_lower <- x[,(d+1):p]

  # Imputation input
  beta1 <- imp$beta1_hat
  m1 <- imp$m1$m
  dm1 <- imp$m1$dm

  beta0 <- imp$beta0_hat
  m0 <- imp$m0$m
  dm0 <- imp$m0$dm

  # CMS
  xb1 <- x %*% beta1
  xb0 <- x %*% beta0

  # IPW input
  alpha_hat <- ipw$alpha_hat
  xa <- x %*% alpha_hat
  pr <- ipw$pr
  d_pr <- ipw$d_pr
  eta <- log(pr / (1 - pr))
  d_eta <- d_pr / (pr * (1 - pr) ) # Makes sense
  # Calculating expected value of combined observed and imputed outcomes (e0,
  # e1)
  #e1 <- rep(sum(y[tbl]) / sum(treated), n)
  #e0 <- rep(sum(y[!tbl]) / (n - sum(treated)), n)


  # Calculating bandwidth
  if (explicit_bandwidth) {
    # Setting explicit bandwidths
    bw1 <- bandwidth_scale1
    bw0 <- bandwidth_scale0
    bw_pr <- bandwidth_scale_pr
  } else {
    # Calculating bandwidths
    sd_xb1 <- sd(xb1[as.logical(treated)])
    bw1 <- bandwidth_scale1 * sd_xb1 * sum(treated)**(-1/5)

    sd_xb0 <- sd(xb0[as.logical(!treated)])
    bw0 <- bandwidth_scale0 * sd_xb0 * sum(!treated)**(-1/5)

    sd_xa <- sd(xa)
    bw_pr <- bandwidth_scale_pr * sd_xa * n**(-1/5)
  }

  # Calculating parameters B B1,0 C1,0 D1,0
  b1 <- b10_fun(x = x,
                treated = treated,
                dm = dm1,
                beta = beta1,
                kernel = kernel,
                bandwidth = bw1,
                gauss_cutoff = gauss_cutoff)

  c1 <- matrix(nrow = 1, ncol = d*(p-d))
  for (i in 1:d){
    c1[1,((i-1)*(p-d) +1):(i*(p-d))] <- colSums(sweep(x_lower * dm1[,i], MARGIN = 1, (n_ones - treated * (n_ones / pr)), "*")) / n

  }

  # Calculating C1B1 for term2
  c1b1 <- c1 %*% b1

  b0 <- b10_fun(x= x,
                treated = (n_ones - treated),
                dm = dm0,
                beta = beta0,
                kernel = kernel,
                bandwidth = bw0,
                gauss_cutoff = gauss_cutoff)


  c0 <- matrix(nrow = 1, ncol = d*(p-d))
  for (i in 1:d){
      c0[1,((i-1)*(p-d) +1):(i*(p-d))] <- colSums(sweep(x_lower * dm0[,i], MARGIN = 1, (n_ones - (n_ones - treated)* (n_ones / (n_ones - pr))), "*")) / n
  }

  # Calculating C0B0 for term3
  c0b0 <- c0 %*% b0

  # Calculating b
  b <- b_fun(x,
             treated,
             alpha_hat,
             num_deriv_h,
             kernel,
             bw_pr,
			       ipw$bw_pr,
             verbose)

  # Calculating Ds
  d1 <- matrix(nrow = 1, ncol = d*(p-d))
  for (i in 1:d){
      d1[1,((i-1)*(p-d) +1):(i*(p-d))] <- colSums(sweep(x_lower, MARGIN = 1, (y - m1) * treated * ( (n_ones/pr) - n_ones) * d_eta[,i], "*")) / n
  }

  # Calculating D1B for term4
  d1b <- d1 %*% b

  d0 <- matrix(nrow = 1, ncol = d*(p-d))
  for (i in 1:d){
      d0[1,((i-1)*(p-d) +1):(i*(p-d))] <- colSums(sweep(x_lower, MARGIN = 1, (y - m0) * (n_ones - treated) * (pr/(n_ones - pr))* d_eta[,i], "*")) / n
  }

  # Calculating D1B for term4
  d0b <- d0 %*% b

  # Calculating expected value of combined observed and imputed outcomes (e0,
  # e1)
  e1 =  mean((y - m1) * treated * (n_ones / pr) + m1)
  e0 =  mean((y - m0) * (n_ones - treated) * (- pr / (n_ones - pr)) + m0)

  # Calculating terms to be summed up
  # Term 1
  term1 <- (y - m1) * treated * (n_ones / pr) + m1 - e1

  # Term 2
  k <- nw_kernel_regress(x_lower, xb1, bandwidth = bw1)
  xc <- x_lower - k

  #term2 <- c1b1 %*% t(sweep(xc, MARGIN = 1, treated * (y - m1) * dm1, "*"))
  term2 <- treated * (y - m1) *t( ((dm1 %x% t(rep(1,(p-d)))) * ( t(rep(1,d)) %x% xc)) %*% t(c1b1) )
  # Term 3
  k <- nw_kernel_regress(x_lower, xa, bandwidth = bw_pr)
  xc <- x_lower - k

  #term3 <- d1b %*% t(sweep(xc, MARGIN = 1, (treated - pr) * d_eta))
  term3 <- (treated - pr) * t( ((d_eta %x% t(rep(1,(p-d)))) * ( t(rep(1,d)) %x% xc)) %*% t(d1b) )

  # Term 4
  term4 <- (y - m0) * (n_ones - treated) * (pr / (n_ones - pr)) - m0 + e0

  # Term 5
  k <- nw_kernel_regress(x_lower, xb0, bandwidth = bw0)
  xc <- x_lower - k

  term5 <- (n_ones - treated) * (y - m0) * t( ((dm0 %x% t(rep(1,(p-d)))) * ( t(rep(1,d)) %x% xc))  %*% t(c0b0) )
  #term5 <- c0b0 %*% t(sweep(xc, MARGIN = 1, (n_ones - treated) * (y - m0) * dm0, "*"))

  # Term 6
  k <- nw_kernel_regress(x_lower, xa, bandwidth = bw_pr)
  xc <- x_lower - k

  #term6 <- d0b %*% t(sweep(xc, MARGIN = 1, (treated - pr) * d_eta))
  term6 <- (treated - pr) * t( ((d_eta %x% t(rep(1,(p-d)))) * ( t(rep(1,d)) %x% xc)) %*% t(d0b)  )

  var <- mean((term1 - term2 + term3 - term4 + term5 + term6)**2)/n

  return(var)
}

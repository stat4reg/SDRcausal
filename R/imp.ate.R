#' @title Estimates Average Treatment Effect (ATE) by imputation (IMP)
#'
#' @description Semiparametric estimation of the average treatment effect based
#'              on the imputation method described in Ghosh, Ma, & De Luna
#'              (2020).
#'
#' @param x                   Covariate matrix
#' @param y                   Response vector
#' @param treated1            A binary vector indicating treatment.
#' @param beta_guess1         Initial guess for \eqn{\beta_1}
#' @param beta_guess0         Initial guess for \eqn{\beta_0}
#' @param solver              Specifies which solver is to be used. Current options are \code{optim} and \code{cobyla} (from \code{nloptr} package). The diffault value is \code{"optim"}.
#' @param kernel              Specifies which kernel function is to be used, current options are: \code{"EPAN"}, \code{"QUARTIC"}, and \code{"GAUSSIAN"}. The default value is \code{"EPAN"}.
#' @param explicit_bandwidth  Specifies if \code{bandwidth_scale} will be used as the bandwidth or if it will be calculated as \code{bandwidth_scale} * sd(\eqn{\beta^T x}) * \eqn{n^{(1/5)}}. The default value is \code{FALSE}.
#' @param recalc_bandwidth    Specifies whether the bandwidth should be recalculated after the first stage (the estimations of dimension reduction step). If \code{explicit_bandwidth} is \code{TRUE}, \code{recalc_bandwidth} is not used, but if \code{explicit_bandwidth} is \code{FALSE}, then if \code{recalc_bandwidth} is \code{TRUE}, bandwidths are recalculated at the beginning of the second step based on \code{bwc_impute0} and \code{bwc_impute1}. If \code{recalc_bandwidth} is \code{FALSE}, the first step bandwidths are used. The default value is \code{FALSE}. 
#' @param bwc_dim_red1        Scaling of calculated bandwidth, or if \code{explicit_bandwidth = TRUE} used as the bandwidth. It is used in the dimension reduction step for \eqn{\hat{m}_1(\beta_1^T x)}. The default value is \code{1}.
#' @param bwc_dim_red0        Scaling of calculated bandwidth, or if \code{explicit_bandwidth = TRUE} used as the bandwidth. It is used in the dimension reduction step for \eqn{\hat{m}_0(\beta_0^T x)}. The default value is \code{1}.
#' @param bwc_impute1         Scaling of calculated bandwidth, or if \code{explicit_bandwidth = TRUE} used as the bandwidth. It is used in the imputation step for \eqn{\hat{m}_1(\beta_1^T x)}. The default  value is \code{1.25}.
#' @param bwc_impute0         Scaling of calculated bandwidth, or if \code{explicit_bandwidth = TRUE} used as the bandwidth. It is used in the imputation step for \eqn{\hat{m}_0(\beta_0^T x)}. The default value is \code{1.25}.
#' @param gauss_cutoff        The cutoff value for Gaussian kernel. The default value is \code{1e-3}.
#' @param penalty             Penalty for the optimizer if local linear regression fails. Added to the function value in solver as \code{penalty}^(n - \code{n_before_pen}), where n is the number of times local linear regression fails. The default value is \code{10}.
#' @param n_before_pen        The number of acceptable local linear regression failures during the estimation of \eqn{\beta_0} and \eqn{\beta_1} phase. The default value is \code{5}.
#' @param to_extrapolate      Specifies whether to extrapolate or not. Since in \eqn{\hat{m}_0(\beta_0^T x)} and \eqn{\hat{m}_1(\beta_1^T x)} estimates in terms of \eqn{\beta_0} and \eqn{\beta_1}, local linear regression at the boundaries of \eqn{\beta_0^Tx} and \eqn{\beta_1^Tx} can be very volatile, it is recommended to use extrapolation on those points instead of local linear regression. The default value is \code{TRUE}. 
#' @param extrapolation_basis The number of data points to base extrapolation on. Extrapolation at border points can be done based on a different number of neighborhood points. \code{extrapolation_basis} is how many neighborhood points are used. The default value is \code{5}.
#' @param to_truncate         Specifies whether to truncate \eqn{\hat{m}_0(\beta_0^T x)} and \eqn{\hat{m}_1(\beta_1^T x)} or not. After estimating \eqn{\hat{m}_0(\beta_0^T x)} and \eqn{\hat{m}_1(\beta_1^T x)}, if they are outside the range of observed outputs, they are replaced with the minimum and maximum observed outputs. The default value is \code{TRUE}. 
#' @param n_threads           Sets the number of threads for parallel computing. Set to 1 serial. If \code{n_threads} exceeds the maximum number of threads, sets \code{n_threads} to max_threads - 1. To use max_threads, set to \code{n_threads} to max_threads of system. The default value is \code{1}.
#' @param verbose             Specifies if the program should print output while running. The default value is \code{TRUE}.
#' @param ...                 Additional parameters passed to \code{optim} or \code{cobyla}.
#' 
#'
#' @return A list containing the average treatment effect of the
#'         combination of observed and imputed values (ate), the average
#'         treatment effect based on the imputed values only (ate2), the
#'         imputed values for treated (m1) and untreated treated (m0), the and
#'         the output from optim (op).
#'
#' @seealso [stats::optim]
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
#'
#' # Perform semiparametric imputation
#' imp <- SDRcausal::imp.ate(x, y, trt, b1, b0,
#'            explicit_bandwidth = TRUE, bwc_dim_red1 = 1, bwc_impute1 = 1,
#'            bwc_dim_red0 = 1, bwc_impute0 = 1)
#'
imp.ate <- function(x,
                               y,
                               treated1,
                               beta_guess1,
                               beta_guess0,
                               solver = "optim",
                               kernel = "EPAN",
                               explicit_bandwidth = FALSE,
                               recalc_bandwidth = TRUE,
                               bwc_dim_red1 = 1,
                               bwc_impute1 = 1.25,
                               bwc_dim_red0 = 1,
                               bwc_impute0 = 1.25,
                               gauss_cutoff = 1e-3,
                               penalty = 10,
                               n_before_pen = 5,
                               to_extrapolate = TRUE,
                               to_truncate = TRUE,
                               extrapolation_basis = 5,
                               n_threads = 1,
                               verbose = TRUE,
                               ...)
{
  # Printing input
  if (verbose) {
    cat("\nPerforming semi parametric imputation\n")
    cat("Input:\n")
    cat("Verbose\t\t", verbose, "\n")
    cat("Kernel\t\t", kernel, "\n")
    cat("Explicit bw\t", explicit_bandwidth, "\n")
    cat("Recalculate bw\t", recalc_bandwidth, "\n")
  }

  # Converting input input
  class(x) <- "numeric"
  class(y) <- "numeric"
  class(treated1) <- "integer"
  class(beta_guess1) <- "numeric"
  class(beta_guess0) <- "numeric"
  class(bwc_dim_red1) <- "numeric"
  class(bwc_impute1) <- "numeric"
  class(bwc_dim_red0) <- "numeric"
  class(bwc_impute0) <- "numeric"
  class(gauss_cutoff) <- "numeric"
  class(penalty) <- "numeric"
  class(n_before_pen) <- "integer"
  class(extrapolation_basis) <- "integer"
  class(n_threads) <- "integer"

  # Deriving dimensions
  n <- as.integer(dim(x)[1])
  p <- as.integer(dim(x)[2])
  # Setting d as dimensions of CMS and testing length
  if (is.matrix(beta_guess1) && is.matrix(beta_guess0)) {
    d <- as.integer(dim(beta_guess1)[2])
    stopifnot(dim(beta_guess1)[1] == p)
    stopifnot(dim(beta_guess0)[1] == p)
    stopifnot(dim(beta_guess0)[2] == d)
  } else {
    d <- as.integer(1)
    stopifnot(length(beta_guess1) == p)
    stopifnot(length(beta_guess0) == p)
  }

  # Creating inverted treatment vector for beta0 and m0
  treated0 <- as.integer(rep(1, times = n) - treated1)

  # Checking so arguments are valid
  kernel <- match.arg(kernel, c("EPAN", "QUARTIC", "GAUSSIAN"))
  if ((bwc_dim_red1 <= 0) ||
      (bwc_impute1 <= 0) ||
      (bwc_dim_red0 <= 0) ||
      (bwc_impute0 <= 0)) {
    stop("Invalid bandwidth scaling. bwcs must be > 0.")
  }
  # Testing data types of input
  stopifnot(is.numeric(x),
            is.numeric(y),
            length(y) == n,
            is.integer(treated1),
            is.integer(treated0),
            is.numeric(beta_guess1),
            is.numeric(beta_guess0),
            is.numeric(bwc_dim_red1),
            is.numeric(bwc_impute1),
            is.numeric(bwc_dim_red0),
            is.numeric(bwc_impute0),
            is.numeric(gauss_cutoff),
            is.numeric(penalty),
            is.integer(n_before_pen),
            is.integer(extrapolation_basis),
            is.integer(n_threads))

  if (verbose) {
    cat("\nReducing dimension for treated (beta 1)\n")
  }
  cms1 <- cms.semi(x = x,
                      y = y,
                      treated = treated1,
                      beta_initial = beta_guess1,
                      solver = solver,
                      kernel = kernel,
                      explicit_bandwidth = explicit_bandwidth,
                      bandwidth_scale = bwc_dim_red1,
                      gauss_cutoff = gauss_cutoff,
                      penalty = penalty,
                      n_before_pen = n_before_pen,
                      n_threads = n_threads,
                      verbose = verbose,
                      ...)

  if (verbose) {
    cat("\nReducing dimension for untreated (beta 0)\n")
  }
  cms0 <- cms.semi(x = x,
                      y = y,
                      treated = treated0,
                      beta_initial = beta_guess0,
                      solver = solver,
                      kernel = kernel,
                      explicit_bandwidth = explicit_bandwidth,
                      bandwidth_scale = bwc_dim_red0,
                      gauss_cutoff = gauss_cutoff,
                      penalty = penalty,
                      n_before_pen = n_before_pen,
                      n_threads = n_threads,
                      verbose = verbose,
                      ...)

  beta1_hat <- cms1$fb
  beta0_hat <- cms0$fb

  # Sets explicit_bandwidth and bandwidth_scale so it is recalculated if specified
  if (!recalc_bandwidth && !explicit_bandwidth) {
    explicit_bandwidth = TRUE
    bwc_impute1 = cms1$bw
    bwc_impute0 = cms0$bw
  }

  if (verbose) {
    cat("\nImputing treated values\n")
  }
  m1 <- imp.val(x = x,
               y = y,
               treated = treated1,
               beta_hat = beta1_hat,
               kernel = kernel,
               explicit_bandwidth = explicit_bandwidth,
               bandwidth_scale = bwc_impute1,
               gauss_cutoff = gauss_cutoff,
               to_extrapolate = to_extrapolate,
               to_truncate = to_truncate,
               extrapolation_basis = extrapolation_basis,
               verbose = verbose)

  if (verbose) {
    cat("\nImputing untreated values\n")
  }
  m0 <- imp.val(x = x,
               y = y,
               treated = treated0,
               beta_hat = beta0_hat,
               kernel = kernel,
               explicit_bandwidth = explicit_bandwidth,
               bandwidth_scale = bwc_impute0,
               gauss_cutoff = gauss_cutoff,
               to_extrapolate = to_extrapolate,
               to_truncate = to_truncate,
               extrapolation_basis = extrapolation_basis,
               verbose = verbose)

  # Calculating ATE
  t_logical <- as.logical(treated1)
  y1 <- m1$m
  y1[t_logical] <- y[t_logical]
  y0 <- m0$m
  y0[!t_logical] <- y[!t_logical]

  ate <- mean(y1 - y0)
  ate2 <- mean(m1$m - m0$m)

  # Creating output list
  output <- list(ate = ate,
                 ate2 = ate2,
                 m1 = m1,
                 m0 = m0,
                 beta1_hat = beta1_hat,
                 beta0_hat = beta0_hat,
                 bw1 = m1$bw,
                 bw0 = m0$bw)
  class(output) <- 'imp'
  return(output)
}

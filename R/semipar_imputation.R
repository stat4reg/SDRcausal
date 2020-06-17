#' @title Estimates Average Treatment Effect (ATE) by imputation (IMP)
#'
#' @description Semiparametric estimation of the average treatment effect based
#'              on the imputation method described in Ghosh, Ma, & De Luna
#'              (2020).
#'
#' @param x                   Covariate matrix
#' @param y                   Response vector
#' @param treated1            Binary vector indicating treatment.
#' @param beta_guess1         Initial guess of beta for m1
#' @param beta_guess0         Initial guess of beta for m0
#' @param solver              Specifies which solver to be used. Current options
#'                            optim and cobyla (from nloptr package).
#' @param kernel              Specifies which kernel function to be used,
#'                            current options are: "EPAN", "QUARTIC", and
#'                            "GAUSSIAN".
#' @param explicit_bandwidth  Specifies if bandwidth_scale will be used as the
#'                            bandwidth or if it will be calculated as bw =
#'                            bandwidth_scale * sd(x * beta) * n^(1/3).
#' @param recalc_bandwidth    Specifies wheter the bandwidth should be
#'                            recalculated after the estimation of alpha
#'                            (ipw_dim_red).
#' @param bwc_dim_red1        Scaling of calculated bandwidth, or if
#'                            explicit_bandwidth = TRUE used as the banddwidth.
#'                            For dimension reduction (imp_dim_red).
#' @param bwc_dim_red0        See bwc_dim_red1
#' @param bwc_impute1         Scaling of calculated bandwidth, or if
#'                            explicit_bandwidth = TRUE used as the banddwidth.
#'                            Recalculated if explicit_bandwidth = FALSE and
#'                            recalc_bandwidth = TRUE. For imputation.
#' @param bwc_impute0         See bwc_impute1
#' @param gauss_cutoff        Cutoff value for Gaussian kernel
#' @param penalty             Penalty for the optimizer if local linear
#'                            regression fails. Added to the function value in
#'                            solver as: penalty^(n - n_before_pen), where n is
#'                            the number of llr fails.
#' @param n_before_pen        Number of probabilities outside the range (0, 1)
#'                            to accept during dimension reduction.
#' @param to_extrapolate      Specifies wheter to extrapolate or not
#' @param to_truncate         Specifies wheter to extrapolate or not
#' @param extrapolation_basis Number of data point to base extrapolation on.
#' @param n_threads           Sets number of threads for parallel run. Set to 0
#'                            serial. If n_threads exceeds maximum number of
#'                            threads, sets n_threads to max_threads - 1. To
#'                            use max_threads, set to n_threads to max_threads
#'                            of system.
#' @param verbose             Specifies if the program should print output
#'                            while running.
#' @param ...                Additional parameters passed to optim.
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
#' imp <- SDRcausal::semipar_imputation(x, y, trt, b1, b0,
#'            explicit_bandwidth = TRUE, bwc_dim_red1 = 1, bwc_impute1 = 1,
#'            bwc_dim_red0 = 1, bwc_impute0 = 1)
#'
semipar_imputation <- function(x,
                               y,
                               treated1,
                               beta_guess1,
                               beta_guess0,
                               solver = "optim",
                               kernel = "EPAN",
                               explicit_bandwidth = FALSE,
                               recalc_bandwidth = FALSE,
                               bwc_dim_red1 = 1,
                               bwc_impute1 = 1,
                               bwc_dim_red0 = 1,
                               bwc_impute0 = 1,
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
  cms1 <- imp_dim_red(x = x,
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
  cms0 <- imp_dim_red(x = x,
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
  m1 <- impute(x = x,
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
  m0 <- impute(x = x,
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

  return(output)
}

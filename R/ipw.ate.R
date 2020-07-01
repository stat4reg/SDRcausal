#' @title Estimates average treatment effect through IPW
#'
#' @description Semiparametric estimation of the average treatment effect based
#'              on the IPW method described in Ghosh, Ma, & De Luna (2020).
#'
#' @param x                  Covariate matrix
#' @param y                  Response vector
#' @param treated            Binary vector indicating treatment.
#' @param alpha_initial      Initial guess of beta for m1
#' @param kernel              Specifies which kernel function to be used,
#'                            current options are: "EPAN", "QUARTIC", and
#'                            "GAUSSIAN".
#' @param explicit_bandwidth Specifies if bandwidth_scale will be used as the
#'                           bandwidth or if it will be calculated as bw =
#'                           bandwidth_scale * sd(x * beta) * n^(1/3).
#' @param recalc_bandwidth   Specifies wheter the bandwidth should be
#'                           recalculated after the estimation of alpha
#'                           (cms.ps.semi)
#' @param bwc_dim_red        Scaling of calculated bandwidth, or if
#'                           explicit_bandwidth = TRUE used as the banddwidth.
#'                           For dimension reduction (cms.ps.semi).
#' @param bwc_prop_score     Scaling of calculated bandwidth, or if
#'                           explicit_bandwidth = TRUE used as the banddwidth.
#'                           Recalculated if explicit_bandwidth = FALSE and
#'                           recalc_bandwidth = TRUE. For propensity score.
#' @param gauss_cutoff       cutoff value for Gaussian kernel
#' @param penalty            Penalty for the optimizer if a probability is
#'                           outside (0, 1) during dimension reduction. Added
#'                           to the function value in solver as: penalty^(n -
#'                           n_before_pen), where n is the number of
#'                           probabilities outside (0, 1).
#' @param n_before_pen       Number of probabilities outside the range (0, 1)
#'                           to accept during dimension reduction.
#' @param n_threads          Sets number of threads for parallel run. Set to 0
#'                           serial. If n_threads exceeds maximum number of
#'                           threads, sets n_threads to max_threads - 1. To
#'                           use max_threads, set to n_threads to max_threads
#'                           of system.
#' @param verbose            Specifies if the program should print output while
#'                           running.
#' @param ...                Additional parameters passed to optim.
#'
#' @return A list containing the average treatment effect (ate), the propensity
#'         score (pr), the final alpha (fa), and the output from optim (op).
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
#' alp <- SDRcausal::alpha_guess
#'
#' # Perform semiparametric inverse probability weighting
#' ipw <- SDRcausal::ipw.ate(x, y, trt, alp, bwc_dim_red = 8,
#'            bwc_prop_score = 8)
#'
ipw.ate <- function(x,
                        y,
                        treated,
                        alpha_initial,
                        kernel = "EPAN",
                        explicit_bandwidth = FALSE,
                        recalc_bandwidth = TRUE,
                        bwc_dim_red = 1,
                        bwc_prop_score = 10,
                        gauss_cutoff = 1e-3,
                        penalty = 10,
                        n_before_pen = 1,
                        n_threads = 1,
                        verbose = TRUE,
                        ...)
{
  # Printing input
  if (verbose) {
    cat("\nPerforming augmented inverse probability weighting\n")
    cat("Input:\n")
    cat("Verbose\t\t", verbose, "\n")
    cat("Kernel\t\t", kernel, "\n")
    cat("Explicit bw\t", explicit_bandwidth, "\n")
    cat("Recalculate bw\t", recalc_bandwidth, "\n")
  }

  # Converting input
  class(x) <- "numeric"
  class(y) <- "numeric"
  class(treated) <- "integer"
  class(alpha_initial) <- "numeric"
  class(bwc_dim_red) <- "numeric"
  class(bwc_prop_score) <- "numeric"
  class(penalty) <- "numeric"
  class(n_before_pen) <- "integer"
  class(n_threads) <- "integer"

  # Deriving dimensions
  n <- as.integer(dim(x)[1])
  p <- as.integer(dim(x)[2])
  # Setting d as dimensions of CMS and testing length
  if (is.matrix(alpha_initial)) {
    d <- as.integer(dim(alpha_initial)[2])
    stopifnot(dim(alpha_initial)[1] == p)
  } else {
    d <- as.integer(1)
    stopifnot(length(alpha_initial) == p)
  }

  # Checking so arguments are valid
  kernel <- match.arg(kernel, c("EPAN", "QUARTIC", "GAUSSIAN"))
  if ((bwc_dim_red <= 0) ||
      (bwc_prop_score <= 0)) {
    stop("Invalid bandwidth scaling. bwcs must be > 0.")
  }

  # Testing data types of input
  stopifnot(is.numeric(x),
            is.numeric(y),
            length(y) == n,
            is.integer(treated),
            is.numeric(alpha_initial),
            is.numeric(bwc_dim_red),
            is.numeric(bwc_prop_score),
            is.numeric(penalty),
            is.integer(n_before_pen),
            is.integer(n_threads),
            is.numeric(gauss_cutoff))

  if (verbose) {
    cat("\nReducing dimension (alpha)\n")
  }
  cms <- cms.ps.semi(x,
                     treated,
                     alpha_initial = alpha_initial,
                     kernel = "EPAN",
                     explicit_bandwidth = explicit_bandwidth,
                     bandwidth_scale = bwc_dim_red,
                     gauss_cutoff = gauss_cutoff,
                     penalty = penalty,
                     n_before_pen = n_before_pen,
                     n_threads = n_threads,
                     verbose = verbose,
                     ...)

  if (verbose) {
    cat("\nEstimating propensity score (pr)\n")
  }

  if (!recalc_bandwidth && !explicit_bandwidth) {
    explicit_bandwidth = TRUE
    bwc_prop_score = cms$bw
  }

  prop_score <- ps.semi(x,
                                 treated,
                                 cms$fa,
                                 kernel = "EPAN",
                                 explicit_bandwidth = explicit_bandwidth,
                                 bandwidth_scale = bwc_prop_score,
                                 verbose = verbose)

  pr <- prop_score$pr
  d_pr <- prop_score$d_pr

  tau_est <- rep(0, length = n)
  tbl <- as.logical(treated)
  tau_est[tbl] <- y[tbl] / pr[tbl]
  tau_est[!tbl] <- -1 * y[!tbl] / (1 - pr[!tbl])

  ate <- mean(tau_est[!is.infinite(tau_est)])

  output <- list(ate = ate,
                 pr = pr,
                 d_pr = d_pr,
                 alpha_hat = cms$fa,
                 optim_op = cms$op,
                 bw_dr = cms$bw,
                 bw_pr = prop_score$bw)
  class(output) <- 'ipw'
  return(output)
}

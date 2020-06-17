#' @title Estimates propensity score
#'
#' @description Semiparametric estimation of the propensity score as in Ghosh,
#'              Ma, & De Luna (2020). To be used with SDRcausal::ipw_dim_red().
#'
#' @param x                  Covariate matrix
#' @param treated            Binary vetor indicating treatment
#' @param alpha_hat          Locally efficient CMS
#' @param kernel             Kernel specification
#' @param explicit_bandwidth Specifies if bandwidth_scale will be used as the
#'                           bandwidth or if it will be calculated as bw =
#'                           bandwidth_scale * sd(x * beta) * n^(1/3).
#' @param bandwidth_scale    Scaling of calculated bandwidth, or if
#'                           explicit_bandwidth = TRUE used as the banddwidth.
#' @param verbose            Specifies if the program should print output while
#'                           running.
#'
#' @return A list containing the estimated propensity scores values and their
#'         derivatives, and the bandwidth used.
#'
#' @useDynLib SDRcausal C_propensity_score
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
#' # Perform semiparametric dimension reduction
#' cms <- SDRcausal::ipw_dim_red(x, trt1, alp,
#'            explicit_bandwidth = TRUE, bandwidth_scale = 8)
#'
#' # Estimate propensity score
#' pr_score <- SDRcausal::propensity_score(x, trt, cms$fa,
#'                bandwidth_scale = 8)
#'
propensity_score <- function(x,
                             treated,
                             alpha_hat,
                             kernel = "EPAN",
                             explicit_bandwidth = FALSE,
                             bandwidth_scale = 1,
                             verbose = FALSE)
{
  # Converting input input
  class(x) <- "numeric"
  class(treated) <- "integer"
  class(alpha_hat) <- "numeric"
  class(bandwidth_scale) <- "numeric"

  # Deriving dimensions
  n <- as.integer(dim(x)[1])
  p <- as.integer(dim(x)[2])
  # Setting d as dimensions of CMS and testing length
  if (is.matrix(alpha_hat)) {
    d <- as.integer(dim(alpha_hat)[2])
    stopifnot(dim(alpha_hat)[1] == p)
  } else {
    d <- as.integer(1)
    stopifnot(length(alpha_hat) == p)
  }

  # Testing so kernel function is valid
  valid_kernels <- list("EPAN", "QUARTIC", "GAUSSIAN")
  stopifnot(kernel %in% valid_kernels)
  kernel <- switch(kernel,
                        EPAN = 1,
                        QUARTIC = 2,
                        GAUSSIAN = 3)
  class(kernel) <- "integer"

  # Testing data types of input
  stopifnot(is.numeric(x),
            is.integer(treated),
            is.numeric(alpha_hat),
            is.numeric(bandwidth_scale))

  # Calculating CMS projection of x and transposing since C code is row major
  xa <- t(x %*% alpha_hat)

  # Creating vectors to recieve C output
  pr <- as.numeric(rep(0, length = n))
  d_pr <- as.numeric(rep(0, length = n * d))
  #d_pr <- as.numeric(matrix(rep(0, length = n * d), nrow = d, ncol = n))
  # Calculating bandwidth
  if (explicit_bandwidth) {
    bw <- bandwidth_scale
  } else {
    sd_xa <- sd(x %*% alpha_hat)
    bw <- bandwidth_scale * sd_xa * dim(x)[1]**(-1/5)
  }

  if (verbose) {
    cat("Using:\n")
    cat("Bandwidth\t", bw, "\n")
    t1 <- Sys.time()
  }

  .Call("C_propensity_score",
        n,
        d,
        xa,
        treated,
        kernel,
        bw,
        pr,
        d_pr)

  output <- list(alpha_hat = alpha_hat,
                 pr = pr,
                 d_pr = matrix(d_pr, ncol = d,byrow = T),
                 bw = bw)

  if (verbose) {
    cat("Finished\n")
    t2 <- Sys.time()
    print(t2 - t1)
  }

  return(output)
}

#' @title Estimates imputed values based on CMS
#'
#' @description Performs semiparametric imputation based on the CMS calculated
#'              by imp_dim_red, as in Ghosh, Ma, & De Luna (2020).
#'
#' @param x                   Covariate matrix
#' @param y                   Response vector
#' @param treated             Binary vetor indicating treatment
#' @param beta_hat            Locally efficient CMS
#' @param kernel              Specifies which kernel function to be used
#' @param explicit_bandwidth  Specifies if bandwidth_scale will be used as the bandwidth or
#'                            if it will be calculated as bw = bandwidth_scale
#'                            * sd(x * beta) * n^(1/3)
#' @param bandwidth_scale     Kernel bandwidth
#' @param gauss_cutoff        Cutoff value for Gaussian kernel
#' @param to_extrapolate      Specifies wheter to extrapolate or not
#' @param to_truncate         Specifies wheter to extrapolate or not
#' @param extrapolation_basis Number of data point to base extrapolation on.
#' @param verbose             Specifies if the program should print output while
#'                            running
#'
#' @return A list containing the reduced space xb, the imputed values and their
#'         derivatives.
#'
#' @useDynLib SDRcausal C_impute
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
#' # Using example data from package SDRcausal
#' library(SDRcausal)
#'
#' # Import example data
#' x <- SDRcausal::covariates
#' y <- SDRcausal::outcomes
#' trt1 <- SDRcausal::treated
#' n <- as.integer(dim(x)[1])
#' trt0 <- as.integer(rep(1, times = n) - trt1)
#' b1 <- SDRcausal::beta1_guess
#' b0 <- SDRcausal::beta0_guess
#'
#' # Perform semiparametric dimension reduction for treated
#' cms1 <- SDRcausal::imp_dim_red(x, y, trt1, b1,
#'            explicit_bandwidth = 1, bandwidth_scale = 1)
#'
#' # Perform semiparametric dimension reduction for untreated
#' cms0 <- SDRcausal::imp_dim_red(x, y, trt0, b0,
#'            explicit_bandwidth = 1, bandwidth_scale = 1)
#'
#' # Perform semiparametric imputation for treated
#' m1 <- SDRcausal::impute(x, y, trt1, cms1$fb,
#'          explicit_bandwidth = 1, bandwidth_scale = cms1$bw)
#'
#' # Perform semiparametric imputation for untreated
#' m0 <- SDRcausal::impute(x, y, trt0, cms0$fb,
#'          explicit_bandwidth = 1, bandwidth_scale = cms0$bw)
#'
impute <- function(x,
                   y,
                   treated,
                   beta_hat,
                   kernel = "EPAN",
                   explicit_bandwidth = FALSE,
                   bandwidth_scale = 1,
                   gauss_cutoff = 1e-3,
                   to_extrapolate = TRUE,
                   to_truncate = TRUE,
                   extrapolation_basis = as.integer(5),
                   verbose = FALSE)
{
  # Converting input input
  class(x) <- "numeric"
  class(y) <- "numeric"
  class(treated) <- "integer"
  class(beta_hat) <- "numeric"
  class(bandwidth_scale) <- "numeric"
  class(gauss_cutoff) <- "numeric"
  class(extrapolation_basis) <- "integer"

  # Deriving dimensions
  n <- as.integer(dim(x)[1])
  p <- as.integer(dim(x)[2])
  # Setting d as dimensions of CMS and testing length
  if (is.matrix(beta_hat)) {
    d <- as.integer(dim(beta_hat)[2])
    stopifnot(dim(beta_hat)[1] == p)
  } else {
    d <- as.integer(1)
    stopifnot(length(beta_hat) == p)
  }

  # Testing so kernel function is valid
  kernel <- match.arg(kernel, c("EPAN", "QUARTIC", "GAUSSIAN"))
  kernel <- switch(kernel,
                   EPAN = 1,
                   QUARTIC = 2,
                   GAUSSIAN = 3)
  class(kernel) <- "integer"

  # Testing data types of input
  stopifnot(is.numeric(x),
            is.numeric(y),
            is.integer(treated),
            is.numeric(beta_hat),
            is.numeric(bandwidth_scale),
            is.numeric(gauss_cutoff),
            is.integer(extrapolation_basis),
            length(y) == n)

  # Transpose since C code is in row major
  xb <- t(x %*% beta_hat)

  # Calculating bandwidth
  if (explicit_bandwidth) {
    bw <- bandwidth_scale
  } else {
    sd_xb <- sd(xb[as.logical(treated)])
    bw <- bandwidth_scale * sd_xb * sum(treated)**(-1/5)
  }

  # Creating vectors to recieve C output
  m <- as.numeric(rep(0, length = n))
  dm <- as.numeric(rep(0, length = n*d))

  if (verbose) {
    cat("Using:\n")
    cat("Bandwidth\t", bw, "\n")
    t1 <- Sys.time()
  }

  .Call("C_impute",
        n,
        d,
        xb,
        y,
        treated,
        kernel,
        bw,
        bw,
        bw,
        bw,
        gauss_cutoff,
        as.integer(to_extrapolate),
        as.integer(to_truncate),
        extrapolation_basis,
        m,
        dm)

  output <- list(m = m,
                 dm = matrix(dm, ncol = d, byrow = T),
                 bw = bw)

  if (verbose) {
    cat("Finished\n")
    t2 <- Sys.time()
    print(t2 - t1)
  }

  return(output)
}

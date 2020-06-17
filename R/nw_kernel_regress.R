#' @title The Nadaraya-Watson kernel estimator
#'
#' @description Gives the expected value of Y given X = x by kernel regression
#'              according to the Nadaraya-Watson kernel estimator to get
#'              E(Y|X). Note that y and x may be vectors or matrices, as long
#'              as dim(x)[1] == dim(y)[1].
#'
#' @param y            Y in E(Y|X)
#' @param x            X in E(Y|X)
#' @param kernel       Indicates which kernel function to be used
#' @param bandwidth    Kernel bandwidth
#' @param gauss_cutoff Cutoff value for Gaussian kernel
#' @param verbose      Specifies if the program should print output while
#'                     running.
#'
#' @return Value of kernel regression
#'
#' @useDynLib SDRcausal C_nw_kernel_regress
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
#'
#' # Extimating y given x, E(y | x)
#' k <- nw_kernel_regress(x, y, bandwidth = 1)
#'
nw_kernel_regress <- function(y,
                              x,
                              bandwidth = 1,
                              kernel = "EPAN",
                              gauss_cutoff = 1e-3,
                              verbose = FALSE)
{
  # Converting input input
  class(x) <- "numeric"
  class(y) <- "numeric"
  class(bandwidth) <- "numeric"
  class(gauss_cutoff) <- "numeric"

  # Deriving dimensions
  n <- as.integer(dim(x)[1])
  p <- as.integer(dim(x)[2])
  if (is.matrix(y)) {
    m <- as.integer(dim(y)[2])
  } else if (is.vector(y)) {
    m <- as.integer(1)
  }


  # Testing so kernel function is valid
  kernel <- match.arg(kernel, c("EPAN", "QUARTIC", "GAUSSIAN"))
  kernel <- switch(kernel,
                   EPAN = 1,
                   QUARTIC = 2,
                   GAUSSIAN = 3)
  class(kernel) <- "integer"

  stopifnot(is.numeric(x),
            is.numeric(y),
            is.numeric(bandwidth),
            is.numeric(gauss_cutoff),
            dim(y)[1] == n)

  if (verbose) {
    cat("Performing kernel regression\n")
    cat("Using:\n")
    cat("Bandwidth\t", bandwidth, "\n")
  }

  k <- as.numeric(matrix(rep(0, times = n*m), nrow = n, ncol = m))
  .Call("C_nw_kernel_regress",
        n,
        p,
        m,
        t(x),
        y,
        kernel,
        bandwidth,
        gauss_cutoff,
        k)


  return(k)
}

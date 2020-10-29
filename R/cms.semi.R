#' @title Estimates the Central Mean Space (CMS)
#'
#' @description Semiparametric estimation of the Central Mean Space (CMS) as in
#'              Ghosh, Ma, & De Luna (2020). To be used with
#'              SDRcausal::imp.val().
#'
#' @param x                  Covariate matrix
#' @param y                  Response vector
#' @param treated            Binary vetor indicating treatment
#' @param beta_initial       Initial guess of CMS
#' @param solver             Specifies which solver to be used. Current options
#'                           optim and cobyla (from nloptr package).
#' @param kernel             Specifies which kernel function to be used,
#'                           current options are: "EPAN", "QUARTIC", and
#'                           "GAUSSIAN".
#' @param explicit_bandwidth Specifies if bandwidth_scale will be used as the
#'                           bandwidth or if it will be calculated as
#'                           bw = bandwidth_scale * sd(x * beta) * n^(1/5)
#' @param bandwidth_scale    Scaling of the bandwidth or the actual bandwidth
#'                           if explicit bandwidth.
#' @param gauss_cutoff       cutoff value for Gaussian kernel
#' @param penalty            Penalty for the optimizer if local linear
#'                           regression fails. Added to the function value in
#'                           solver as: penalty^(n - n_before_pen), where n is
#'                           the number of llr fails.
#' @param n_before_pen       Number of probabilities outside the range (0, 1)
#'                           to accept during dimension reduction.
#' @param root_tol           Tolerance which makes the program warn if optim
#'                           stops at at a value higher than root_tol.
#' @param n_threads          Sets number of threads for parallel run. Set to 0
#'                           serial. If n_threads exceeds maximum number of
#'                           threads, sets n_threads to max_threads - 1. To use
#'                           max_threads, set to n_threads to max_threads of
#'                           system.
#' @param verbose            Specifies if the program should print output while
#'                           running.
#' @param ...                Additional parameters passed to optim.
#'
#' @return A list containing the final beta, the bandwidth used, a warning
#'         if optim does not converge or converges to a value that is larger
#'         than root_tol, and the output of optim.
#'
#' @importFrom stats optim
#' @seealso [stats::optim]
#'
#' @useDynLib SDRcausal C_imp_dim_red
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
#' trt1 <- SDRcausal::treated
#' trt0 <- rep(1, length(trt1)) - trt1
#' b1 <- SDRcausal::beta1_guess
#' b0 <- SDRcausal::beta0_guess
#'
#' # Perform semiparametric dimension reduction for treated
#' cms1 <- SDRcausal::cms.semi(x, y, trt1, b1,
#'            explicit_bandwidth = TRUE, bandwidth_scale = 1)
#'
#' # Perform semiparametric dimension reduction for untreated
#' cms0 <- SDRcausal::cms.semi(x, y, trt0, b0,
#'            explicit_bandwidth = TRUE, bandwidth_scale = 1)
#'
cms.semi <- function(x,
                        y,
                        treated,
                        beta_initial,
                        solver = "optim",
                        kernel = "EPAN",
                        explicit_bandwidth = FALSE,
                        bandwidth_scale = 1,
                        gauss_cutoff = 1e-3,
                        penalty = 10,
                        n_before_pen = 1,
                        root_tol = 1e-3,
                        n_threads = 1,
                        verbose = FALSE,
                        ...)
{
  # Converting input input
  class(x) <- "numeric"
  class(y) <- "numeric"
  class(treated) <- "integer"
  class(beta_initial) <- "numeric"
  class(bandwidth_scale) <- "numeric"
  class(gauss_cutoff) <- "numeric"
  class(penalty) <- "numeric"
  class(n_before_pen) <- "integer"
  class(n_threads) <- "integer"

  # Deriving dimensions
  n <- as.integer(dim(x)[1])
  p <- as.integer(dim(x)[2])
  # Setting d as dimensions of CMS and testing length
  if (is.matrix(beta_initial)) {
    d <- as.integer(dim(beta_initial)[2])
    stopifnot(dim(beta_initial)[1] == p)
  } else {
    d <- as.integer(1)
    stopifnot(length(beta_initial) == p)
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
            is.numeric(beta_initial),
            is.numeric(bandwidth_scale),
            is.numeric(gauss_cutoff),
            is.numeric(penalty),
            is.integer(n_before_pen),
            is.integer(n_threads),
            length(y) == n)

  # Calculating bandwidth
  if (explicit_bandwidth) {
    bw <- bandwidth_scale
  } else {
    xb <- x %*% beta_initial
    sd_xb<- sd(xb[as.logical(treated)])
    bw <- bandwidth_scale * sd_xb * sum(treated)**(-1/5)
  }

  # Parameter list for optimization
    params <- list(n = n,
                   p = p,
                   d = d,
                   x = t(x),
                   y = y,
                   t = treated,
                   kspec = kernel,
                   bw = bw,
                   b0 = beta_initial[1],
                   n_thr = n_threads,
                   pen = penalty,
                   n_pen = n_before_pen,
                   gco = gauss_cutoff)

    # Initial guess of lower p-d beta
    if (d == 1) {
      beta_guess = beta_initial[2:length(beta_initial)]
    } else {
      beta_guess = as.vector(beta_initial[(d+1):p,])
      params$b0 = beta_initial[1:d,]
    }


    # Function to be minimized
    dr_fun <- function(beta_guess, params)
    {
      # Concatenating beta with identity matrix of dimension d
      beta_guess <- matrix(beta_guess, ncol = params$d)
      beta_guess <- rbind(params$b0, beta_guess)

      # Multi d case
      # beta_guess <- rbind(diag(params$d), beta_guess)

      # 1d implementation (multi d would need transpose)
      xb <- x %*% beta_guess

      fval <- 0
      .Call("C_imp_dim_red",
            params$n,
            params$p,
            params$d,
            params$x,
            xb,
            params$y,
            params$t,
            params$kspec,
            params$bw,
            params$bw,
            params$bw,
            params$bw,
            params$bw,
            params$gco,
            params$pen,
            params$n_pen,
            params$n_thr,
            fval)
      if (verbose) {
        cat("\rfval =", fval)
      }

      return(fval)
    }


  if (solver == "optim") {
    if (verbose) {
      cat("Using:\n")
      cat("N threads\t", n_threads, "\n")
      if ("method" %in% names(list(...))) {
        cat("Method\t\t", list(...)$method, "\n")
      }
      else {
        cat("Method\t\t", "Nelder-Mead", "\n")
      }
      cat("Initial beta\t")
      cat(paste0("[", paste0(trimws(format(round(beta_initial, digits = 2),
                                           nsmall = 2)), collapse = " "), "]\n"))
      cat("Bandwidth\t", bw, "\n")
      t1 <- Sys.time()
    }

    # Performing minimization
    optim_output <- optim(beta_guess,
                          dr_fun,
                          params = params,
                          ...)

    # Checking if optim converges and warns in case not
    if (optim_output$convergence != 0) {
      cat("\nWarning: Optim did not converge\n")
      cat("Cause:")
      cat(optim_output$convergence)
      cat("\n\n")
      warning = 1
    } else if (optim_output$value > root_tol) {
      cat("\nWarning: Final function value high\n")
      cat("Function value:")
      cat(optim_output$value)
      cat("\n")
      warning = 2
    } else {
      warning = 0
    }

    # 1d implementation
    # Concatenating beta with identity matrix of dimension d
    beta_final <- matrix(optim_output$par, ncol = params$d)
    beta_final <- rbind(params$b0, beta_final)

    # Multi d case
    # beta_final <- rbind(diag(params$d), as.matrix(beta_final))

    output <- list(fb = beta_final,
                   bw = bw,
                   op = optim_output,
                   warning = warning)

    if (verbose) {
      cat("\nFinished\n")
      t2 <- Sys.time()
      print(t2 - t1)
      cat("\n")
      cat("Final beta\t")
      cat(paste0("[", paste0(trimws(format(round(beta_final, digits = 2), nsmall
                                           = 2)), collapse = " "), "]\n"))
    }

    return(output)
  }
  else if (solver == "cobyla") {

    stopifnot("nloptr" %in% (.packages()))

    if (verbose) {
      cat("Using:\n")
      cat("N threads\t", n_threads, "\n")
      cat("Method\t\t", "COBYLA", "\n")
      cat("Initial beta\t")
      cat(paste0("[", paste0(trimws(format(round(beta_initial, digits = 2),
                                           nsmall = 2)), collapse = " "), "]\n"))
      cat("Bandwidth\t", bw, "\n")
      t1 <- Sys.time()
    }

    # Performing minimization
    cobyla_output <- nloptr::cobyla(x0 = beta_guess,
                            fn = dr_fun,
                            params = params,
                            ...)

    # Checking if cobyla converges and warns in case not
    if (cobyla_output$convergence < 0) {
      cat("\n\nWarning: COBYLA did not converge\n")
      cat("Cause:")
      cat(cobyla_output$convergence)
      cat("\n")
      warning = 1
    } else if (cobyla_output$value > root_tol) {
      cat("\nWarning: Final function value high\n")
      cat("Function value:")
      cat(cobyla_output$value)
      cat("\n\n")
      warning = 2
    } else {
      warning = 0
    }

    # 1d implementation
    # Concatenating beta with identity matrix of dimension d
    beta_final <- matrix(cobyla_output$par, ncol = params$d)
    beta_final <- rbind(params$b0, beta_final)

    # Multi d case
    # beta_final <- rbind(diag(params$d), as.matrix(beta_final))

    output <- list(fb = beta_final,
                   bw = bw,
                   op = cobyla_output,
                   warning = warning)

    if (verbose) {
      cat("\nFinished\n")
      t2 <- Sys.time()
      print(t2 - t1)
      cat("\n")
      cat("Final beta\t")
      cat(paste0("[", paste0(trimws(format(round(beta_final, digits = 2), nsmall
                                           = 2)), collapse = " "), "]\n"))
    }

    return(output)


  }
  else {
    cat("Unknown solver", solver, "\n\n")
    stop()
  }
}

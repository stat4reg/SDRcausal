#' @title Estimates the Central Mean Space (CMS)
#'
#' @description Semiparametric estimation of the Central Mean Space (CMS) as in
#'              Ghosh, Ma, & De Luna (2020). To be used with SDRcausal::propensity_score().
#'
#' @param x                  Covariate matrix
#' @param treated            Binary vetor indicating treatment
#' @param alpha_initial      Initial guess of CMS
#' @param solver             Specifies which solver to be used. Current options
#'                           optim and cobyla (from nloptr package).
#' @param kernel             Specifies which kernel function to be used,
#'                           current options are: "EPAN", "QUARTIC", and
#'                           "GAUSSIAN".
#' @param explicit_bandwidth Specifies if bandwidth_scale will be used
#'                           explicitly as the bandwidth.
#' @param bandwidth_scale    Scaling of the calculated bandwidth, or in case of
#'                           explicit_bandwidth = TRUE the bandwidth.
#' @param gauss_cutoff       cutoff value for Gaussian kernel
#' @param penalty            Penalty for the optimizer if a probability is
#'                           outside (0, 1). Added to the function value in
#'                           optim as: penalty^(n), where n is the number of
#'                           probabilities outside (0, 1).
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
#' @param ...                Additional parameters passed to solver.
#'
#' @return A list containing the final alpha, bandwwidth used, and the output
#'         of optim
#'
#' @importFrom stats optim
#' @seealso [stats::optim]
#'
#' @useDynLib SDRcausal C_ipw_dim_red
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
#' # Perform semiparametric dimension reduction for treated
#' cms <- SDRcausal::ipw_dim_red(x, trt1, alp,
#'            explicit_bandwidth = TRUE, bandwidth_scale = 8)
#'
ipw_dim_red <- function(x,
                        treated,
                        alpha_initial,
                        solver = "optim",
                        kernel = "EPAN",
                        explicit_bandwidth = FALSE,
                        bandwidth_scale = 1,
                        gauss_cutoff = 1e-3,
                        penalty = 10,
                        n_before_pen = 5,
                        root_tol = 1e-3,
                        n_threads = 1,
                        verbose = FALSE,
                        ...)
{
  # Converting input input
  class(x) <- "numeric"
  class(treated) <- "integer"
  class(alpha_initial) <- "numeric"
  class(bandwidth_scale) <- "numeric"
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

  # Testing so kernel function is valid
  kernel <- match.arg(kernel, c("EPAN", "QUARTIC", "GAUSSIAN"))
  kernel <- switch(kernel,
                   EPAN = 1,
                   QUARTIC = 2,
                   GAUSSIAN = 3)
  class(kernel) <- "integer"

  # Testing data types of input
  stopifnot(is.numeric(x),
            is.integer(treated),
            is.numeric(alpha_initial),
            is.numeric(bandwidth_scale),
            is.numeric(gauss_cutoff),
            is.numeric(penalty),
            is.integer(n_before_pen),
            is.integer(n_threads))

  # Calculating bandwidth
  if (explicit_bandwidth) {
    bw <- bandwidth_scale
  } else {
    sd_xa <- sd(x %*% alpha_initial)
    bw <- bandwidth_scale * sd_xa * dim(x)[1]**(-1/5)
  }

  # Parameter list for optimization
  params <- list(n = n,
                 p = p,
                 d = d,
                 # Transpose since C code is row major
                 x = t(x),
                 t = treated,
                 kernel = kernel,
                 bw = bw,
                 gco = gauss_cutoff,
                 a0 = alpha_initial[1],
                 n_thr = n_threads,
                 pen = penalty,
                 n_pen = n_before_pen,
                 n_thr = n_threads)


  # Initial guess of lower p-d alpha
  if (d == 1) {
    alpha_guess = alpha_initial[2:length(alpha_initial)]
  } else {
    alpha_guess = as.vector(alpha_initial[(d+1):p,])
    params$a0 = alpha_initial[1:d,]
  }

  # Function to be minimized
  dr_fun <- function(alpha_guess, params)
  {
    # 1d case
    alpha_guess <- matrix(alpha_guess, ncol = params$d)
    alpha_guess <- rbind(params$a0, alpha_guess)

    # Multi d case
    # alpha_guess <- rbind(diag(params$d), alpha_guess)

    # 1d implementation (multi d would need transpose)
    xa <- x %*% alpha_guess
    fval <- as.numeric(0)

    .Call("C_ipw_dim_red",
          params$n,
          params$p,
          params$d,
          params$x,
          xa,
          params$t,
          params$kernel,
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

  # The code below is divided into two if statements, one for optim and one for
  # cobyla.
  # Optim solver
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
      cat("Bandwidth\t", bw, "\n")
      cat("Initial alpha\t")
      cat(paste0("[", paste0(trimws(format(round(alpha_initial, digits = 2),
                                           nsmall = 2)), collapse = " "), "]\n"))
      t1 <- Sys.time()
    }

    # Performing minimization
    optim_output <- optim(alpha_guess,
                          dr_fun,
                          params = params,
                          ...)

    # Checking if optim converges and warns in case not
    if (optim_output$convergence != 0) {
      cat("\n\nWarning: Optim did not converge\n")
      cat("Cause:")
      cat(optim_output$convergence)
      cat("\n")
      warning = 1
    } else if (optim_output$value > root_tol) {
      cat("\nWarning: Final function value high\n")
      cat("Function value:")
      cat(optim_output$value)
      cat("\n\n")
      warning = 2
    } else {
      warning = 0
    }


    # 1d case
    alpha_final <- matrix(optim_output$par, ncol = params$d)
    alpha_final <- rbind(params$a0, alpha_final)

    # Multi d case
    # alpha_final <- matrix(optim_output$par, ncol = params$d)
    # alpha_final <- rbind(diag(params$d), as.matrix(alpha_final))

    output <- list(fa = alpha_final,
                   bw = bw,
                   op = optim_output,
                   warning = warning)

    if (verbose) {
      cat("\nFinished\n")
      t2 <- Sys.time()
      print(t2 - t1)
      cat("\n")
      cat("Final alpha\t")
      cat(paste0("[", paste0(trimws(format(round(alpha_final, digits = 2), nsmall
                                           = 2)), collapse = " "), "]\n"))
    }

    return(output)
  }
  # Cobyla solver
  else if (solver == "cobyla") {

    stopifnot("nloptr" %in% (.packages()))

    if (verbose) {
      cat("Using:\n")
      cat("N threads\t", n_threads, "\n")
      cat("Method\t\t", "COBYLA", "\n")
      cat("Bandwidth\t", bw, "\n")
      cat("Initial alpha\t")
      cat(paste0("[", paste0(trimws(format(round(alpha_initial, digits = 2),
                                           nsmall = 2)), collapse = " "), "]\n"))
      t1 <- Sys.time()
    }

    # Performing minimization
    cobyla_output <- cobyla(x0 = alpha_guess,
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


    # 1d case
    alpha_final <- matrix(cobyla_output$par, ncol = params$d)
    alpha_final <- rbind(params$a0, alpha_final)

    # Multi d case
    # alpha_final <- matrix(cobyla_output$par, ncol = params$d)
    # alpha_final <- rbind(diag(params$d), as.matrix(alpha_final))

    output <- list(fa = alpha_final,
                   bw = bw,
                   op = cobyla_output,
                   warning = warning)

    if (verbose) {
      cat("\nFinished\n")
      t2 <- Sys.time()
      print(t2 - t1)
      cat("\n")
      cat("Final alpha\t")
      cat(paste0("[", paste0(trimws(format(round(alpha_final, digits = 2), nsmall
                                           = 2)), collapse = " "), "]\n"))
    }

    return(output)

  }
  else {
    cat("Unknown solver", solver, "\n\n")
    stop()
  }
}

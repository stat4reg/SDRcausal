#' @title Performs Estimations of Average Treatment Effect and Infrences 
#'
#' @description Semiparametric estimation of the average treatment effect based
#'              on the all methods described in Ghosh, Ma, & De Luna
#'              (2020) and all infrences.
#'
#' @param x                   Covariate matrix
#' @param y                   Response vector
#' @param treated            Binary vector indicating treatment.
#' @param beta_guess1         Initial guess of beta for m1
#' @param beta_guess0         Initial guess of beta for m0
#' @param imp.solver              Specifies which solver to be used. Current options
#'                            optim and cobyla (from nloptr package).
#' @param imp.kernel              Specifies which kernel function to be used,
#'                            current options are: "EPAN", "QUARTIC", and
#'                            "GAUSSIAN".
#' @param imp.explicit_bandwidth  Specifies if bandwidth_scale will be used as the
#'                            bandwidth or if it will be calculated as bw =
#'                            bandwidth_scale * sd(x * beta) * n^(1/3).
#' @param imp.recalc_bandwidth    Specifies wheter the bandwidth should be
#'                            recalculated after the estimation of alpha
#'                            (cms.ps.semi).
#' @param bwc_dim_red1        Scaling of calculated bandwidth, or if
#'                            explicit_bandwidth = TRUE used as the banddwidth.
#'                            For dimension reduction (cms.semi).
#' @param bwc_dim_red0        See bwc_dim_red1
#' @param bwc_impute1         Scaling of calculated bandwidth, or if
#'                            explicit_bandwidth = TRUE used as the banddwidth.
#'                            Recalculated if explicit_bandwidth = FALSE and
#'                            recalc_bandwidth = TRUE. For imputation.
#' @param bwc_impute0         See bwc_impute1
#' @param imp.gauss_cutoff        Cutoff value for Gaussian kernel
#' @param imp.penalty             Penalty for the optimizer if local linear
#'                            regression fails. Added to the function value in
#'                            solver as: penalty^(n - n_before_pen), where n is
#'                            the number of llr fails.
#' @param imp.n_before_pen        Number of probabilities outside the range (0, 1)
#'                            to accept during dimension reduction.
#' @param imp.to_extrapolate      Specifies wheter to extrapolate or not
#' @param imp.to_truncate         Specifies wheter to extrapolate or not
#' @param imp.extrapolation_basis Number of data point to base extrapolation on.
#' @param n_threads           Sets number of threads for parallel run. Set to 0
#'                            serial. If n_threads exceeds maximum number of
#'                            threads, sets n_threads to max_threads - 1. To
#'                            use max_threads, set to n_threads to max_threads
#'                            of system.
#' @param verbose             Specifies if the program should print output
#'                            while running.
#' 
#' @param alpha_initial      Initial guess of beta for m1
#' @param ipw.kernel              Specifies which kernel function to be used,
#'                            current options are: "EPAN", "QUARTIC", and
#'                            "GAUSSIAN".
#' @param ipw.explicit_bandwidth Specifies if bandwidth_scale will be used as the
#'                           bandwidth or if it will be calculated as bw =
#'                           bandwidth_scale * sd(x * beta) * n^(1/3).
#' @param ipw.recalc_bandwidth   Specifies wheter the bandwidth should be
#'                           recalculated after the estimation of alpha
#'                           (cms.ps.semi)
#' @param bwc_dim_red        Scaling of calculated bandwidth, or if
#'                           explicit_bandwidth = TRUE used as the banddwidth.
#'                           For dimension reduction (cms.ps.semi).
#' @param bwc_prop_score     Scaling of calculated bandwidth, or if
#'                           explicit_bandwidth = TRUE used as the banddwidth.
#'                           Recalculated if explicit_bandwidth = FALSE and
#'                           recalc_bandwidth = TRUE. For propensity score.
#' @param ipw.gauss_cutoff       cutoff value for Gaussian kernel
#' @param ipw.penalty            Penalty for the optimizer if a probability is
#'                           outside (0, 1) during dimension reduction. Added
#'                           to the function value in solver as: penalty^(n -
#'                           n_before_pen), where n is the number of
#'                           probabilities outside (0, 1).
#' @param ipw.n_before_pen       Number of probabilities outside the range (0, 1)
#'                           to accept during dimension reduction.
#' @param ipw.solver         Specifies which solver to be used. Current options
#'                            optim and cobyla (from nloptr package).
#' @param imp.solver.options Additional parameters passed to optim for imp.
#' @param ipw.solver.options Additional parameters passed to optim for ipw.
#'
#'
#'
#'
#'
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
#' a <- SDRcausal::alpha_guess
#'
#' # Perform semiparametric imputation
#' inf.ate <- SDRcausal::inf.ate(x, y, trt, b1, b0, alpha_initial = a)
#'
inf.ate <- function(x,
                    y,
                    treated,
                    beta_guess1,
                    beta_guess0,
                    imp.solver = "optim",
                    imp.kernel = "EPAN",
                    imp.explicit_bandwidth = FALSE,
                    imp.recalc_bandwidth = TRUE,
                    bwc_dim_red1 = 1,
                    bwc_impute1 = 1.25,
                    bwc_dim_red0 = 1,
                    bwc_impute0 = 1.25,
                    imp.gauss_cutoff = 1e-3,
                    imp.penalty = 10,
                    imp.n_before_pen = 5,
                    imp.to_extrapolate = TRUE,
                    imp.to_truncate = TRUE,
                    imp.extrapolation_basis = 5,
                    alpha_initial,
                    ipw.solver = "optim",
                    ipw.kernel = "EPAN",
                    ipw.explicit_bandwidth = FALSE,
                    ipw.recalc_bandwidth = TRUE,
                    bwc_dim_red = 1,
                    bwc_prop_score = 10,
                    ipw.gauss_cutoff = 1e-3,
                    ipw.penalty = 10,
                    ipw.n_before_pen = 1,
                    n_threads = 1,
                    verbose = TRUE,
                    imp.solver.options = NA,
                    ipw.solver.options = NA)
{
  
  imp.list = list("x"= x,
          'y' = y,
          'treated1' = treated,
          'beta_guess1' = beta_guess1,
          'beta_guess0' = beta_guess0,
          'solver' = imp.solver,
          'kernel' = imp.kernel,
          'explicit_bandwidth' = imp.explicit_bandwidth,
          'recalc_bandwidth' = imp.recalc_bandwidth,
          'bwc_dim_red1' = bwc_dim_red1,
          'bwc_impute1' = bwc_impute1,
          'bwc_dim_red0' = bwc_dim_red0,
          'bwc_impute0' = bwc_impute0,
          'gauss_cutoff' = imp.gauss_cutoff,
          'penalty' = imp.penalty,
          'n_before_pen' = imp.n_before_pen,
          'to_extrapolate' = imp.to_extrapolate,
          'to_truncate' = imp.to_truncate,
          'extrapolation_basis' = imp.extrapolation_basis,
          'n_threads' = n_threads,
          'verbose' = verbose)
  if (!is.na(imp.solver.options)){
    imp.list = append(imp.list, imp.solver.options )
  }
  imp = do.call(imp.ate , args = imp.list)
  
  
  ipw.list =list('x' = x, 'y' = y,'treated' = treated,
                      'alpha_initial' = alpha_initial,
                      'solver' = ipw.solver,
                      'kernel' = ipw.kernel,
                      'explicit_bandwidth' = ipw.explicit_bandwidth,
                      'recalc_bandwidth' = ipw.recalc_bandwidth,
                      'bwc_dim_red' = bwc_dim_red,
                      'bwc_prop_score' = bwc_prop_score,
                      'gauss_cutoff' = ipw.gauss_cutoff,
                      'penalty' = ipw.penalty,
                      'n_before_pen' = ipw.n_before_pen,
                      'n_threads' = n_threads,
                      'verbose' = verbose)
  if (!is.na(ipw.solver.options)){
    ipw.list = append(ipw.list, ipw.solver.options )
  }
  ipw = do.call(ipw.ate , args = ipw.list)
  
  aipw = aipw.ate(y,treated,
                       imp,
                       ipw)
  aipw2 = aipw2.ate(y,treated,
                        imp,
                        ipw)
  #### var
  imp_var =  imp.var(x,y,treated,imp,
                                 ipw,
                                 kernel = imp.kernel,
                                 explicit_bandwidth = imp.explicit_bandwidth,
                                 gauss_cutoff = imp.gauss_cutoff)
    
  imp2_var = imp2.var(x,
                      y,
                      treated,
                      imp,
                      ipw,
                      kernel = imp.kernel,
                      explicit_bandwidth = imp.explicit_bandwidth,
                      gauss_cutoff = imp.gauss_cutoff)
  ipw_var = ipw.var(x,y,treated,
                                imp,
                                ipw,
                                kernel = ipw.kernel,
                                explicit_bandwidth = ipw.explicit_bandwidth,
                                gauss_cutoff = ipw.gauss_cutoff,
                                verbose = verbose)
  aipw_var = aipw.var(x,y,treated,
                                  imp,
                                  ipw,
                                  kernel = ipw.kernel,
                                  explicit_bandwidth = ipw.explicit_bandwidth,
                                  gauss_cutoff = ipw.gauss_cutoff,
                                  verbose = verbose)
  
  return(list('imp'=imp, 'ipw'=ipw, 'aipw'=aipw, 'aipw2'=aipw2, 'imp_var'=imp_var, 'imp2_var'=imp2_var, 'ipw_var'=ipw_var, 'aipw_var'=aipw_var))
    
}

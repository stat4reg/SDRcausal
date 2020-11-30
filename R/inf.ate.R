#' @title Performs Estimations of Average Treatment Effect and Infrences 
#'
#' @description Semiparametric estimation of the average treatment effect based
#'              on the all methods described in Ghosh, Ma, & De Luna
#'              (2020) and all infrences.
#'
#' @param x                   Covariate matrix
#' @param y                   Response vector
#' @param treated             A binary vector indicating treatment status
#' @param beta_guess1         Initial guess for \eqn{\beta_1}
#' @param beta_guess0         Initial guess for \eqn{\beta_0}
#' @param imp.solver          Specifies which solver is to be used. Current options are \code{optim} and \code{cobyla} (from \code{nloptr} package). The diffault value is \code{"optim"}.
#' @param imp.kernel          Specifies which kernel function is to be used, current options are: \code{"EPAN"}, \code{"QUARTIC"}, and \code{"GAUSSIAN"}. The default value is \code{"EPAN"}.
#' @param imp.explicit_bandwidth  Specifies if \code{bandwidth_scale} will be used as the bandwidth or if it will be calculated as \code{bandwidth_scale} * sd(\eqn{\beta^T x}) * \eqn{n^{(1/5)}}. The default value is \code{FALSE}.
#' @param imp.recalc_bandwidth    Specifies whether the bandwidth should be recalculated after the first stage (the estimations of dimension reduction step). If \code{explicit_bandwidth} is \code{TRUE}, \code{recalc_bandwidth} is not used, but if \code{explicit_bandwidth} is \code{FALSE}, then if \code{recalc_bandwidth} is \code{TRUE}, bandwidths are recalculated at the beginning of the second step based on \code{bwc_impute0} and \code{bwc_impute1}. If \code{recalc_bandwidth} is \code{FALSE}, the first step bandwidths are used. The default value is \code{FALSE}. 
#' @param bwc_dim_red1        Scaling of calculated bandwidth, or if \code{explicit_bandwidth = TRUE} used as the bandwidth. It is used in the dimension reduction step for \eqn{\hat{m}_1(\beta_1^T x)}. The default value is \code{1}.
#' @param bwc_dim_red0        Scaling of calculated bandwidth, or if \code{explicit_bandwidth = TRUE} used as the bandwidth. It is used in the dimension reduction step for \eqn{\hat{m}_0(\beta_0^T x)}. The default value is \code{1}.
#' @param bwc_impute1         Scaling of calculated bandwidth, or if \code{explicit_bandwidth = TRUE} used as the bandwidth. It is used in the imputation step for \eqn{\hat{m}_1(\beta_1^T x)}. The default  value is \code{1.25}.
#' @param bwc_impute0         Scaling of calculated bandwidth, or if \code{explicit_bandwidth = TRUE} used as the bandwidth. It is used in the imputation step for \eqn{\hat{m}_0(\beta_0^T x)}. The default value is \code{1.25}.
#' @param imp.gauss_cutoff    The cutoff value for Gaussian kernel. The default value is \code{1e-3}.
#' @param imp.penalty         Penalty for the optimizer if local linear regression fails. Added to the function value in solver as \code{penalty}^(n - \code{n_before_pen}), where n is the number of times local linear regression fails. The default value is \code{10}.
#' @param imp.n_before_pen    The number of acceptable local linear regression failures during the estimation of \eqn{\beta_0} and \eqn{\beta_1} phase. The default value is \code{5}.
#' @param imp.to_extrapolate  Specifies whether to extrapolate or not. Since in \eqn{\hat{m}_0(\beta_0^T x)} and \eqn{\hat{m}_1(\beta_1^T x)} estimates in terms of \eqn{\beta_0} and \eqn{\beta_1}, local linear regression at the boundaries of \eqn{\beta_0^Tx} and \eqn{\beta_1^Tx} can be very volatile, it is recommended to use extrapolation on those points instead of local linear regression. The default value is \code{TRUE}. 
#' @param imp.extrapolation_basis The number of data points to base extrapolation on. Extrapolation at border points can be done based on a different number of neighborhood points. \code{extrapolation_basis} is how many neighborhood points are used. The default value is \code{5}.
#' @param imp.to_truncate     Specifies whether to truncate \eqn{\hat{m}_0(\beta_0^T x)} and \eqn{\hat{m}_1(\beta_1^T x)} or not. After estimating \eqn{\hat{m}_0(\beta_0^T x)} and \eqn{\hat{m}_1(\beta_1^T x)}, if they are outside the range of observed outputs, they are replaced with the minimum and maximum observed outputs. The default value is \code{TRUE}. 
#' @param n_threads           Sets the number of threads for parallel computing. Set to 1 serial. If \code{n_threads} exceeds the maximum number of threads, sets \code{n_threads} to max_threads - 1. To use max_threads, set to \code{n_threads} to max_threads of system. The default value is \code{1}.
#' @param verbose             Specifies if the program should print output while running. The default value is \code{TRUE}.
#' 
#' @param alpha_initial      Initial guess for \eqn{\alpha}
#' @param ipw.kernel         Specifies which kernel function is to be used, current options are: \code{"EPAN"}, \code{"QUARTIC"}, and \code{"GAUSSIAN"}. The default value is \code{"EPAN"}.
#' @param ipw.explicit_bandwidth Specifies if \code{bandwidth_scale} will be used as the bandwidth or if it will be calculated as \code{bandwidth_scale} * sd(\eqn{\alpha^T x}) * \eqn{n^{(1/5)}}. The default value is \code{FALSE}.
#' @param ipw.recalc_bandwidth   Specifies whether the bandwidth should be recalculated after the estimations of \eqn{\alpha}. If \code{explicit_bandwidth} is \code{TRUE}, \code{recalc_bandwidth} is not used, but if \code{explicit_bandwidth} is \code{FALSE}, then if \code{recalc_bandwidth} is \code{TRUE}, bandwidths are recalculated at the beginning of the second step based on \code{bwc_prop_score}. If \code{recalc_bandwidth} is \code{FALSE}, the first step bandwidths are used. The default value is \code{FALSE}. 
#' @param bwc_dim_red        Scaling of calculated bandwidth, or if \code{explicit_bandwidth = TRUE} used as the bandwidth. It is used in the dimension reduction step for \eqn{\alpha^T x}. The default value is \code{1}.
#' @param bwc_prop_score     Scaling of calculated bandwidth, or if \code{explicit_bandwidth = TRUE} used as the bandwidth. It is used for the estimation of the propensity score. The default value is \code{10}.
#' @param ipw.gauss_cutoff   The cutoff value for Gaussian kernel. The default value is \code{1e-3}.
#' @param ipw.penalty        Penalty for the optimizer if a probability is outside (0, 1) during the estimation of \eqn{\alpha} phase. Added to the function value in solver as \code{penalty}^(n - \code{n_before_pen}), where n is the number of probabilities outside (0, 1). The default value is \code{10}.
#' @param ipw.n_before_pen   The number of probabilities outside the range (0, 1) to accept during the estimation of \eqn{\alpha} phase. The default value is \code{1}.
#' @param ipw.solver         Specifies which solver is to be used. Current options are \code{optim} and \code{cobyla} (from \code{nloptr} package). The diffault value is \code{"optim"}.
#' @param imp.solver.options Additional parameters passed to \code{optim} or \code{cobyla} for \code{imp.ate}.
#' @param ipw.solver.options Additional parameters passed to \code{optim} or \code{cobyla} for \code{ipw.ate}.
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

#' @title Combines IPW and IMP estimators to form the augmented IPW, AIPW
#'
#' @description Augmented IPW (AIPW) as in Ghosh, Ma, & De Luna (2020).
#'
#' @param y       Observed response
#' @param treated Binary vetor indicating treatment
#' @param imp     imp output object from \code{imp.ate()}
#' @param ipw     ipw output object from \code{ipw.ate()}
#'
#' @return Average treatment effect (ATE) for the augmented IPW (AIPW)
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
#' alp <- SDRcausal::alpha_guess
#'
#' # Perform semiparametric imputation
#' imp <- SDRcausal::imp.ate(x, y, trt, b1, b0,
#'            explicit_bandwidth = TRUE, bwc_dim_red1 = 1, bwc_impute1 = 1,
#'            bwc_dim_red0 = 1, bwc_impute0 = 1)
#'
#' # Perform semiparametric inverse probability weighting
#' ipw <- SDRcausal::ipw.ate(x, y, trt, alp, bwc_dim_red = 8,
#'            bwc_prop_score = 8)
#'
#' # Calculate the Augmented IPW (AIPW)
#' aipw <- SDRcausal::aipw.ate(y, trt, imp, ipw)
#'
aipw.ate <- function(y,
                        treated,
                        imp,
                        ipw)
{
  m1 <- imp$m1$m
  m0 <- imp$m0$m
  pr <- ipw$pr

  n_ones <- rep(1, length(y))

  y1 <- (treated * y / pr) + (n_ones - (treated / pr)) * m1

  y0 <- ((n_ones - treated) * (y / (n_ones - pr))  + (n_ones - ((n_ones- treated) / (n_ones - pr))) * m0)

  ATE <- mean(y1) - mean(y0)

  return(ATE)
}

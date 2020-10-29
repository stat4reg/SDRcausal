#' @title Improved Augmented IPW (AIPW2)
#'
#' @description Combines IPW and IMP estimators to form the improved augmented
#'              IPW, AIPW2 as in Ghosh, Ma, & De Luna (2020).
#'
#' @param y       Observed response
#' @param treated Binary vetor indicating treatment
#' @param imp     imp output object from \code{imp.ate()}
#' @param ipw     ipw output object from \code{ipw.ate()}
#'
#' @return Average treatment effect (ATE) for the improved augmented IPW
#'         (AIPW2)
#'
#' @export
#'
#' @importFrom stats cov
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
#' ipw <- SDRcausal::ipw.ate(x, y, trt, alp, bwc_dim_red = 10,
#'            bwc_prop_score = 18)
#'
#' # Calculate the Improved Augmented IPW (AIPW2)
#' aipw2 <- SDRcausal::aipw2.ate(y, trt, imp, ipw)
#'
aipw2.ate <- function(y,
                      treated,
                      imp,
                      ipw)
{
  m1 <- imp$m1$m
  m0 <- imp$m0$m
  pr <- ipw$pr
  ones <- rep(1, length(y))

  g1a <- cov(m1 * treated / pr, m1 * (ones - treated / pr))
  g1b <- cov(treated * y / pr, (ones - treated / pr) * m1)
  gamma_1 <- rep(g1b / g1a, length(y))

  y1 <- treated * y / pr + gamma_1 * (ones - treated / pr) * m1

  g0a <- cov(m0 * (ones - treated) / (ones - pr), m0 *
             (ones - (ones - treated) / (ones - pr)))
  g0b <- cov((ones - treated) * y / (ones - pr), m0 *
             (ones - (ones - treated) / (ones - pr)))
  gamma_0 <- rep(g0b / g0a, length(y))

  y0 <- (ones - treated) * y / (ones - pr) + gamma_0 *
    (ones - (ones - treated) / (ones - pr)) * m0

  ATE <- mean(y1) - mean(y0)

  return(ATE)
}

#' @title Plots imputation output
#'
#' @description Plot function for visualisation of imputation output from
#'              semipar_imputation. Note: The function requires ggplot2.
#'
#' @param x                   Covariate matrix
#' @param y                   Response vector
#' @param treated             Binary vetor indicating treatment
#' @param imp                 imp_output object from semipar_imputation()
#'
#' @return A list of ggplot plots of observed and imputed values (pl_imp),
#'         imputed treated values vs CMS (pl_m1), and imputed untreated values
#'         vs CMS (pl_m0).
#'
#' @export
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
#' imp <- SDRcausal::semipar_imputation(x, y, trt, b1, b0,
#'            explicit_bandwidth = TRUE, bwc_dim_red1 = 1, bwc_impute1 = 1,
#'            bwc_dim_red0 = 1, bwc_impute0 = 1)
#'
#' # Plotting
#' plots <- plot_imp(x, y, trt, imp)
#'
plot_imp <- function(x,
                     y,
                     treated,
                     imp)
{
  stopifnot("ggplot2" %in% (.packages()))

  # Number of observations
  n <- length(y)
  n_seq <- seq(1, n)

  # CMS projections of covariate matrix
  xb1 <- x %*% imp$beta1_hat
  xb0 <- x %*% imp$beta0_hat

  # Treated vector as boolean
  tbl <- as.logical(treated)

  # Creating outcome vectors based on observed and imputet values
  y1 <- imp$m1$m
  y1[tbl] <- y[tbl]
  y0 <- imp$m0$m
  y0[!tbl] <- y[!tbl]

  # Creating vector to label observed or imputed.
  tpl <- vector(, length = n)
  tpl[tbl] <- "treated"
  tpl[!tbl] <- "untreated"
  # Creating vector to label observed or imputed.
  ipl1 <- vector(, length = n)
  ipl1[tbl] <- "observed"
  ipl1[!tbl] <- "imputed"
  ipl0 <- vector(, length = n)
  ipl0[tbl] <- "imputed"
  ipl0[!tbl] <- "observed"

  # Plot of outcome vs observation
  p1 <- ggplot(,aes(x = n_seq)) +
    geom_point(aes(y = y1, color = "treated", shape = ipl1)) +
    geom_point(aes(y = y0, color = "untreated", shape = ipl0)) +
    scale_shape_manual(values=c(17, 16)) +
    guides(shape = guide_legend(title = "Obs/Imp"),
           color = guide_legend(title = "Treatment")) +
    xlab("Observation") +
    ylab("Response") +
    ggtitle("Imputed and observed values VS index")

  # Plot of imputed treated outcomes vs CMS projection
  p2 <- ggplot(, aes(x = xb1)) +
    geom_point(aes(y = imp$m1$m, shape = "imputed", color = tpl), alpha = 0.4) +
    geom_point(aes(y = y, shape = "observed", color = tpl), alpha = 0.4) +
    geom_line(aes(y = imp$m1$m, color = tpl), alpha = 0.6) +
    scale_shape_manual("Source", values = c(2, 0, 1, 3)) +
    guides(shape = guide_legend(title = "Obs/Imp"),
           color = guide_legend(title = "Treatment"),
           linetype = guide_legend(title = " ")) +
    ggtitle("Imputated and observed outcomes VS subspace beta[t]^T x", ) +
    xlab("CMS") +
    ylab("Outcome")

  # Plot of imputed untreated outcomes vs CMS projection
  p3 <- ggplot(, aes(x = xb0)) +
    geom_point(aes(y = imp$m0$m, shape = "imputed", color = tpl), alpha = 0.4) +
    geom_point(aes(y = y, shape = "observed", color = tpl), alpha = 0.4) +
    geom_line(aes(y = imp$m0$m, color = tpl), alpha = 0.6) +
    scale_shape_manual("Source", values = c(2, 0, 1, 3)) +
    guides(shape = guide_legend(title = "Obs/Imp"),
           color = guide_legend(title = "Treatment"),
           linetype = guide_legend(title = " ")) +
    ggtitle("Imputation of treated outcomes") +
    xlab("CMS") +
    ylab("Outcome")

  return(list(pl_imp = p1, pl_m1 = p2, pl_m0 = p3))

}

#' @title Plots IPW output
#'
#' @description Plot function for visualisation of IPW output from semipar_ipw.
#'              Note: The function requires ggplot2.
#'
#' @param treated             Binary vetor indicating treatment
#' @param ipw                 ipw_output object from semipar_ipw()
#'
#' @return ggplot plot of the propensity score vs CMS.
#'
#' @export
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
#' imp <- SDRcausal::semipar_imputation(x, y, trt, b1, b0,
#'            explicit_bandwidth = TRUE, bwc_dim_red1 = 1, bwc_impute1 = 1,
#'            bwc_dim_red0 = 1, bwc_impute0 = 1)
#'
#' # Plotting
#' plots <- plot_imp(x, y, trt, imp)
#'
plot_ipw <- function(treated,
                     ipw)
{
  stopifnot("ggplot2" %in% (.packages()))

  # CMS projections of covariate matrix
  xa <- x %*% ipw$alpha_hat

  pl <- ggplot(, aes(x = xa)) +
    geom_point(aes(y = ipw$pr, color = "propensity score"), alpha = 0.4,
               size = 3) +
    geom_line(aes(y = ipw$pr), alpha = 0.6) +
    geom_point(aes(y = treated, color = "treatment"), color = "black",
               alpha = 0.03, size = 4) +
    xlab("CMS") +
    ylab("Propensity score")

  return(pl)
}

#' @title Example data
#'
#' @description Data generated as in paper, study 1. Using the betas in betas
#'              data. Use beta1/0 for imputation as the initial guess of the
#'              central mean space (CMS) and alpha as the initial guess of the
#'              CMS for IPW.
#'
#' @name example_data
#' @docType data
#' @keywords data
#'
#' @format Data used in examples of the SDRcausal package
#' \describe{
#'   \item{covariates}{covariate matrix}
#'   \item{outcomes}{observed outcome vector}
#'   \item{treated}{binary treatment vector}
#'   \item{beta1_guess}{Starting guess for CMS for treated}
#'   \item{beta0_guess}{Starting guess for CMS for untreated}
#'   \item{alpha_guess}{Starting guess for CMS for propensity score}
#' }
NULL

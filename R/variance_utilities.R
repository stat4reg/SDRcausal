#' @title Calculates B1/0
#'
#' @description Calculates Eq 2.8 or 2.10 in Ghosh, Ma, & De Luna (2020).
#'
#' @param x_lower      Projection of covariate matrix on CMS
#' @param treated      Binary vector indicating treatment.
#' @param dm           Derivative of imputed values
#' @param beta         CMS
#' @param kernel       Specifies which kernel function to be used
#' @param bandwidth    Specifies if bandwidth_scale will be used as the
#' @param gauss_cutoff Cutoff value for Gaussian kernel
#'
#' @return B1/0 matrix
#'
b10_fun <- function(x,
                    treated,
                    dm,
                    beta,
                    kernel,
                    bandwidth,
                    gauss_cutoff)
{
  # Derive dimensions
  n <- as.integer(dim(x)[1])
  p <- as.integer(dim(x)[2])
  d <- as.integer(dim(dm)[2])
  p0 <- p - d

  # Extracts lower p-d elements of X and performs kernel regression to get
  # the expectation of the same with respect to beta x X
  x_lower <- x[,(d+1):p]
  x_lower_est <- nw_kernel_regress(x_lower,
                                   x %*% beta,
                                   bandwidth = bandwidth,
                                   kernel = kernel,
                                   gauss_cutoff = gauss_cutoff)
  xc <- x_lower - x_lower_est

  part1 <- - ((dm ) %x% t(rep(1,p0)) ) * (t(rep(1,d)) %x% xc)

  part2 <- matrix(rep(0, p0*p0*d*d), nrow = p0*d, ncol = p0*d)
  for (i in seq(1, n)) {
    if (treated[i]) {
      part2 <- part2 + treated[i] * (dm[i,] %x% x_lower[i,]) %o% part1[i,]
    }
  }

  part2 <- part2 / n

  b_out <- solve(part2 )

  return(b_out)
}

#' @title Calculates B1/0
#'
#' @description Calculates Eq 2.8 or 2.10 in Ghosh, Ma, & De Luna (2020).
#'
#' @param x         Projection of covariate matrix on CMS
#' @param treated   Treated
#' @param alpha_hat Derivative of imputed values
#' @param h         CMS
#' @param kernel    Specifies which kernel function to be used
#' @param bandwidth Kernel bandwidth
#' @param verbose   Specifies if the program should print output while running.
#'
#' @return B1/0 matrix
#'
b_fun <- function(x,
                  treated,
                  alpha_hat,
                  h,
                  kernel,
                  bandwidth,
				  bandwidth_pr,
                  verbose)
{
  # Deriving dimensions
  n <- as.integer(dim(x)[1])
  p <- as.integer(dim(x)[2])
  d <- as.integer(dim(alpha_hat)[2])
  p0 <- p - d

  x_lower <- x[,(d+1):p]

  xa <- x %*% alpha_hat
  xc <- x_lower - nw_kernel_regress(x_lower, xa, bandwidth = bandwidth)

  # Extracting the p-d lower elements of alpha for numerical integration
  if (d==1){
      alpha_low <- alpha_hat[(d+1):p]
      alpha_high <- alpha_hat[1:d]
  }else{
      alpha_low <- as.vector(alpha_hat[(d+1):p,])
      alpha_high <- alpha_hat[1:d,]
  }


  # Creating matrices where each row represents a small step in alpha[i] to the
  # left and right respectively, where i is the row of the matrix.
  h_mat <- h * diag(p0 *d)
  alpha_left <- (matrix(rep(alpha_low, times = p0*d), nrow = p0*d, ncol = p0*d)
                 - h_mat)
  alpha_right <- (matrix(rep(alpha_low, times = p0*d), nrow = p0*d, ncol = p0*d)
                  + h_mat)

  # Creating matrices to store the left and right results
  b_left <- matrix(, nrow = p0*d, ncol = p0*d)
  b_right <- matrix(, nrow = p0*d, ncol = p0*d)

  # Loop over all the rows
  for (i in seq(1, p0*d)) {
    # Calculating propensity score based on beta_left
    gamma = h/2
    xa_left <- if (d==1){x %*% c(alpha_hat[1], alpha_left[,i])}else{x %*% rbind(alpha_high, matrix(alpha_left[,i], ncol= d))}


    d_eta_left = apply(xa_left,2 ,function(xa_left_v){

      xa_left_v = matrix(xa_left_v)
      xa_left_gamma = xa_left_v *(1 + gamma)
      pr_left <- nw_kernel_regress(as.numeric(treated), xa_left_v, bandwidth = bandwidth)
      pr_left_gamma <- nw_kernel_regress(as.numeric(treated), xa_left_gamma, bandwidth = bandwidth)
      pr_left[pr_left >.99] <- 0.99
      pr_left_gamma[pr_left_gamma >.99] <- 0.99
      pr_left[pr_left <.01] <- 0.01
      pr_left_gamma[pr_left_gamma <.01] <- 0.01
      d_pr_left = (pr_left_gamma - pr_left)/(xa_left_v * gamma)
      return(d_pr_left/(pr_left * (1 - pr_left)))
    })
    pr_left <- nw_kernel_regress(as.numeric(treated), xa_left, bandwidth = bandwidth)
    xc_left <- x_lower - nw_kernel_regress(x_lower, xa_left, bandwidth = bandwidth)
    #b_left[i,] <- colSums(sweep(xc, MARGIN = 1, (treated - pr_left) * d_eta,"*")) / n
    b_left[i,] <- colSums( (t(rep(1, d)) %x% xc_left) * ( ((treated - pr_left) * d_eta_left) %x% t(rep(1, p-d)) ) ) / n



    # Calculating propensity score based on beta_right
    gamma = h/2

    xa_right <- if (d==1){x %*% c(alpha_hat[1], alpha_right[,i])}else{x %*% rbind(alpha_high, matrix(alpha_right[,i], ncol= d))}
    d_eta_right = apply(xa_right,2 ,function(xa_right){

    xa_right = matrix(xa_right)
    xa_right_gamma = xa_right *(1 + gamma)
    pr_right <- nw_kernel_regress(as.numeric(treated), xa_right, bandwidth = bandwidth)
    pr_right_gamma <- nw_kernel_regress(as.numeric(treated), xa_right_gamma, bandwidth = bandwidth)
    pr_right[pr_right >.99] <- 0.99
    pr_right_gamma[pr_right_gamma >.99] <- 0.99
    pr_right[pr_right <.01] <- 0.01
    pr_right_gamma[pr_right_gamma <.01] <- 0.01
    d_pr_right = (pr_right_gamma - pr_right)/(xa_right * gamma)
    return(d_pr_right/(pr_right * (1 - pr_right)))
    })

    pr_right <- nw_kernel_regress(as.numeric(treated), xa_right, bandwidth = bandwidth)
    xc_right <- x_lower - nw_kernel_regress(x_lower, xa_right, bandwidth = bandwidth)
    #b_right[i,] <- colSums(sweep(xc, MARGIN = 1, (treated - pr_right) * d_eta,"*")) / n
    #b_right[i,] <- colSums(sweep(xc_right, MARGIN = 1, (treated - pr_right) * d_eta_right,"*")) / n
    b_right[i,] <- colSums( (t(rep(1, d)) %x% xc_right) * ( ((treated - pr_right) * d_eta_right) %x% t(rep(1, p-d)) ) ) / n
  }

  # Performing the numerical derivative, creating a jacobian matrix, and
  # takes the inverse of the matrix.
  b_out <- solve(t(b_right - b_left) / (2 * h))

  return(b_out)
}

//     ------------------------------------------------------------------------
//
//     This file is part of SDRcausal.
//
//     SDRcausal is free software: you can redistribute it and/or modify it
//     under the terms of the GNU General Public License as published by the
//     Free Software Foundation, either version 3 of the License, or (at your
//     option) any later version.
//
//     SDRcausal is distributed in the hope that it will be useful, but WITHOUT
//     ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//     FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//     for more details.
//
//     You should have received a copy of the GNU General Public License along
//     with SDRcausal.  If not, see <https://www.gnu.org/licenses/>.
//
//     ------------------------------------------------------------------------
//
#ifndef C_IPW_DIM_RED_H
#define C_IPW_DIM_RED_H

/**
 * @brief Wrapper function for ipw_dim_red which is used to estimate alpha.
 *
 * @param n_in         Number of observations
 * @param p_in         Number of covariates
 * @param d_in         Structural dimension
 * @param x            Covariate matrix
 * @param x_alpha      CMS projection of x
 * @param treated      Binary treatment vector
 * @param kernel_spec  Indicates which kernel function to be used
 * @param h1           Kernel bandwidth
 * @param h2           Kernel bandwidth
 * @param gauss_cutoff Cutoff value for Gaussian kernel
 * @param penalty      Penalty for probabilities outside of (0, 1)
 * @param n_pen_in     Number of points to accept before penalising
 * @param n_threads_in Number of threads for parallel computing
 * @param fval         Norm squared of output from cseff
 *
 * @return R_NilValue
 */
SEXP
C_ipw_dim_red (SEXP       n_in,
               SEXP       p_in,
               SEXP       d_in,
               const SEXP x,
               const SEXP x_alpha,
               const SEXP treated,
               SEXP       kernel_spec,
               SEXP       h1,
               SEXP       h2,
               SEXP       gauss_cutoff,
               SEXP       penalty,
               SEXP       n_pen_in,
               SEXP       n_threads_in,
               SEXP       fval);

#endif

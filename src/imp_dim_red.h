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
#ifndef IMP_DIM_RED_H
#define IMP_DIM_RED_H

/**
 * @brief Calculates Eq. 2.4 from paper.
 *
 *
 * @param n            Number of observations
 * @param p            Number of covariates
 * @param d            Structural dimension
 * @param x            Covariate matrix
 * @param x_beta       CMS projection of x
 * @param y            Response vector
 * @param treated      Binary treatment vector
 * @param kernel_spec  Indicates which kernel function to be used
 * @param h_0          Kernel bandwidth
 * @param h11          Kernel bandwidth
 * @param h12          Kernel bandwidth
 * @param h13          Kernel bandwidth
 * @param h14          Kernel bandwidth
 * @param gauss_cutoff Cutoff value for Gaussian kernel
 * @param n_threads    Number of threads for parallel computing
 * @param n_llr_fail        Number of local linear regression fails
 * @param y_out        Equation 2.4 before summation
 */
void
imp_dim_red (int            n,
             int            p,
             int            d,
             const double   x[n*p],
             const double   x_beta[n*d],
             const double   y[n],
             const int      treated[n],
             int            kernel_spec,
             double         h0,
             double         h11,
             double         h12,
             double         h13,
             double         h14,
             double         gauss_cutoff,
             int            n_threads,
             int           *n_llr_fail,
             double       (*y_out)[d*(p-d)]);

#endif

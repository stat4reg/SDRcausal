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
#ifndef IPW_DIM_RED_H
#define IPW_DIM_RED_H

/**
 * @brief
 *
 * @param n            Number of observations
 * @param d            Structural dimension
 * @param x            Matrix of observations
 * @param x_alpha      CMS projection of x
 * @param treated      Binary treatment vector
 * @param kernel_spec  Indicates which kernel function to be used
 * @param h1           Kernel bandwidth
 * @param h2           Kernel bandwidth
 * @param gauss_cutoff Cutoff value for Gaussian kernel
 * @param n_pr_nan     Number of probabilities that is outside of (0.01, 0.99)
 * @param col_sum      Summation of the columns of the estimating equation
 */
void
ipw_dim_red (int            n,
             int            p,
             int            d,
             const double   x[n*p],
             const double   x_alpha[n*d],
             const int      treated[n],
             int            kernel_spec,
             double         h1,
             double         h2,
             double         gauss_cutoff,
             int           *n_pr_nan,
             double       (*col_sum)[(p-d)*d]);

#endif

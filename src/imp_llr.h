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
#ifndef IMP_LLR_H
#define IMP_LLR_H

/**
 * @brief Performs local linear regression by kernel estimation.
 *
 * @param n            Number of observations
 * @param d            Structural dimension
 * @param x            CMS projection of covariate matrix
 * @param x0           Current x_beta value
 * @param y            Response vector
 * @param treated      Binary treatment vector
 * @param x0           Reference observation, one row of x_beta
 * @param kernel_spec  Indicates which kernel function to be used
 * @param h11          Kernel bandwidth
 * @param h12          Kernel bandwidth
 * @param h13          Kernel bandwidth
 * @param h14          Kernel bandwidth
 * @param gauss_cutoff Cutoff value for Gaussian kernel
 * @param m            Imputed values
 * @param dm           Derivative of imputed values
 */
void
imp_llr (int            n,
         int            d,
         const double   x[n*d],
         const double   x0[d],
         const double   y[n],
         const int      treated[n],
         int            kernel_spec,
         double         h11,
         double         h12,
         double         h13,
         double         h14,
         double         gauss_cutoff,
         double        *m,
         double       (*dm)[d]);

#endif

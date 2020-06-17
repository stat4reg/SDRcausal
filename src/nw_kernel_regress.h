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
#ifndef NW_KERNEL_REGRESS_H
#define NW_KERNEL_REGRESS_H

/**
 * @brief Estimates E(Y | X) using the Nadaraya-Watson kernel estimator.
 *
 * @param n            number of rows of X and Y
 * @param m            number of columns of Y
 * @param u            distances at to evaluate kernel function
 * @param y            Y in E(Y | X)
 * @param kernel_spec  Indicates which kernel function to be used
 * @param bandwidth    kernel bandwidth
 * @param gauss_cutoff cutoff value for Gaussian kernel
 * @param k            the estimated values
 */
void
nw_kernel_regress (int            n,
                   int            m,
                   const double   u[n],
                   const double   y[n*m],
                   int            kernel_spec,
                   double         bandwidth,
                   double         gauss_cutoff,
                   double       (*k)[m]);

/**
 * @brief Estimates E(Y | X) using the Nadaraya-Watson kernel estimator,
 * accounting for treatment.
 *
 * @param n            number of rows of X and Y
 * @param m            number of columns of Y
 * @param u            distances at to evaluate kernel function
 * @param y            Y in E(Y | X)
 * @param treated      binary treatment vector
 * @param kernel_spec  Indicates which kernel function to be used
 * @param bandwidth    kernel bandwidth
 * @param gauss_cutoff cutoff value for Gaussian kernel
 * @param k            the estimated values
 */
void
nw_kernel_regress_treatment (int            n,
                             int            m,
                             const double   u[n],
                             const double   y[n*m],
                             const int      treated[n],
                             int            kernel_spec,
                             double         bandwidth,
                             double         gauss_cutoff,
                             double       (*k)[m]);
#endif

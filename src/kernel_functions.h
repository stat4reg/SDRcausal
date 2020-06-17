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
#ifndef KERNEL_FUNCTIONS_H
#define KERNEL_FUNCTIONS_H

/**
 * @brief Calculates the Epanechnikov kernel for u with bandwidth h.
 *
 * @param u input to kernel function
 * @param h bandwidth
 *
 * @return Kernel function value
 */
double
epanechnikov_kernel (double u,
                     double h);

/**
 * @brief Calculates the biweight quartic kernel for u with bandwidth h.
 *
 * @param u input to kernel function
 * @param h bandwidth
 *
 * @return Kernel function value
 */
double
quartic_biweight_kernel (double u,
                         double h);

/**
 * @brief Calculates the Gaussian kernel for u.
 *
 * @param u input to kernel function
 * @param h bandwidth
 * @param c cutoff
 *
 * @return Kernel function value
 */
double
gaussian_kernel (double u,
                 double h,
                 double c);

/**
 * @brief Chooses kernel function and evaluates it
 *
 * @param u           input to kernel function
 * @param kernel_spec indicates which kernel function to be used
 * @param h           kernel bandwidth
 * @param c           gaussian cutoff
 */
double
kernel_evaluator (double u,
                  int    kernel_spec,
                  double h,
                  double c);

#endif

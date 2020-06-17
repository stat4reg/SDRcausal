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
#ifndef C_IMPUTE_H
#define C_IMPUTE_H

/**
 * @brief Wrapper function for impute which computes the imputed values.
 *
 * @param n_in                Number of observations
 * @param d_in                Structural dimension
 * @param x_beta              CMS projection of x
 * @param y                   Response vector
 * @param treated             Binary treatment vector
 * @param kernel_spec         Indicates which kernel function to be used
 * @param h11                 Kernel bandwidth
 * @param h12                 Kernel bandwidth
 * @param h13                 Kernel bandwidth
 * @param h14                 Kernel bandwidth
 * @param gauss_cutoff        Cutoff value for Gaussian kernel
 * @param extrapolation_basis Number of point to base extrapolation on.
 * @param to_extrapolate      Specifies wheter to extrapolate or not
 * @param to_truncate         Specifies wheter to truncate or not
 * @param m                   Imputed values
 * @param dm                  Derivatives of imputed values
 *
 * @return R_NilValue
 */
SEXP
C_impute (SEXP       n_in,
          SEXP       d_in,
          const SEXP x_beta,
          const SEXP y,
          const SEXP treated,
          SEXP       kernel_spec,
          SEXP       h11,
          SEXP       h12,
          SEXP       h13,
          SEXP       h14,
          SEXP       extrapolation_basis,
          SEXP       gauss_cutoff,
          SEXP       to_extrapolate,
          SEXP       to_truncate,
          SEXP       m,
          SEXP       dm);

#endif

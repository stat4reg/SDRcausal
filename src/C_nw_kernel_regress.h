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
#ifndef C_NW_KERNEL_REGRESS_H
#define C_NW_KERNEL_REGRESS_H
/**
 * @brief Wrapper function for nw_kernel_regress.
 *
 * @param n_in         Number of rows of X and Y
 * @param d_in         Structural dimension
 * @param m_in         Number of columns of Y
 * @param u            Distances at to evaluate kernel function
 * @param y            Y in E(Y | X)
 * @param kernel_spec  Indicates which kernel function to be used
 * @param bandwidth    Kernel bandwidth
 * @param gauss_cutoff Cutoff value for Gaussian kernel
 * @param k            The estimated values
 * @param fval         Pointer to store output value
 *
 * @return R_NilValue
 */
SEXP
C_nw_kernel_regress (SEXP       n_in,
                     SEXP       d_in,
                     SEXP       m_in,
                     const SEXP u,
                     const SEXP y,
                     SEXP       kernel_spec,
                     SEXP       bandwidth,
                     SEXP       gauss_cutoff,
                     SEXP       k);

#endif

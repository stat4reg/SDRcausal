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
#ifndef C_PROPENSITY_SCORE_H
#define C_PROPENSITY_SCORE_H

/**
 * @brief Wrapper function for propensity_score.
 *
 * @param n_in         Number of observations
 * @param d_in         Structural dimension
 * @param x_beta_in    CMS projection of x
 * @param treated      Binary treatment vector
 * @param kernel_spec  Indicates which kernel function to be used
 * @param h1           Kernel bandwidth
 * @param pr_out       Probabilities
 * @param d_pr_out     Derivative of probabilities
 *
 * @return R_NilValue
 */
SEXP
C_propensity_score (SEXP       n_in,
                    SEXP       d_in,
                    const SEXP x_beta_in,
                    const SEXP treated,
                    SEXP       kernel_spec,
                    SEXP       h1,
                    SEXP       pr_out,
                    SEXP       d_pr_out);

#endif

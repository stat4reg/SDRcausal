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
#ifndef PROPENSITY_SCORE_H
#define PROPENSITY_SCORE_H

/**
 * @brief Performs local linear regression by kernel estimation.
 *
 * @param n           Number of observations
 * @param d           Structural dimension
 * @param x_beta      CMS projection of x
 * @param x0          Reference observation, one row of x_beta
 * @param treated     Binary treatment vector
 * @param kernel_spec Indicates which kernel function to be used
 * @param h1          Kernel bandwidth
 * @param pr          Propensity score estimation
 * @param d_pr        Derivative of propensity score estimation
 */
void
propensity_score (int            n,
                  int            d,
                  const double   x_beta[n*d],
                  const double   x0[d],
                  const int      treated[n],
                  int            kernel_spec,
                  double         h1,
                  double        *pr,
                  double       (*d_pr)[d]);

#endif

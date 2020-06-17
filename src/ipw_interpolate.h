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
#ifndef IPW_INTERPOLATE_H
#define IPW_INTERPOLATE_H

/**
 * @brief Interpolates values of pr.
 *
 * @param n       Number of observations
 * @param d       Structural dimension
 * @param x_alpha CMS projection of covariate matrix
 * @param pr      Propensity score
 */
void
ipw_interpolate (int            n,
                 int            d,
                 const double   x_alpha[n*d],
                 double       (*pr)[n],
				 double       (*d_pr)[n*d]);

#endif

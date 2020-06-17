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
#ifndef IMP_EXTRAPOLATE_H
#define IMP_EXTRAPOLATE_H

/**
 * @brief Extrapolates the edges of imputed values (m).
 *
 * @param n                   Number of observations
 * @param d                   Structural dimension
 * @param x_beta              CMS projection of covariate matrix
 * @param treated             Binary treatment vector
 * @param extrapolation_basis Number of points to base extrapolation on.
 * @param m                   Imputed values
 * @param dm                  Derivatives of imputed values
 * @param n_extrapolated      Number of extrapolated elements
 */
void
imp_extrapolate (int            n,
                 int            d,
                 const double   x_beta[n*d],
                 const int      treated[n],
                 int            extrapolation_basis,
                 double       (*m)[n],
                 double       (*dm)[n*d],
                 int           *n_extrapolated);

#endif

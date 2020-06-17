//     ------------------------------------------------------------------------
//
//     This file is part of SemiparEstimators.
//
//     SemiparEstimators is free software: you can redistribute it and/or
//     modify it under the terms of the GNU General Public License as published
//     by the Free Software Foundation, either version 3 of the License, or (at
//     your option) any later version.
//
//     SemiparEstimators is distributed in the hope that it will be useful, but
//     WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//     General Public License for more details.
//
//     You should have received a copy of the GNU General Public License
//     along with SemiparEstimators.  If not, see
//     <https://www.gnu.org/licenses/>.
//
//     ------------------------------------------------------------------------
//
#ifndef IMP_TRUNCATE_H
#define IMP_TRUNCATE_H

/**
 * @brief Truncating values of m with respect to y.
 *
 * @param n           number of observations
 * @param y           response vector
 * @param m           imputed values
 * @param n_truncated number of extrapolated elements
 */
void
imp_truncate (int            n,
              const double   y[n],
              double       (*m)[n],
              int           *n_truncated);

#endif

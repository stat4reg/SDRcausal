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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "matrix_utilities.h"
#include "imp_local_linear_regression.h"
#include "nw_kernel_regress.h"
#include "imp_dim_red.h"

void
imp_dim_red (int            n,
             int            p,
             int            d,
             const double   x[n*p],
             const double   x_beta[n*d],
             const double   y[n],
             const int      treated[n],
             int            kernel_spec,
             double         h0,
             double         h11,
             double         h12,
             double         h13,
             double         h14,
             double         gauss_cutoff,
             double       (*y_out)[d*(p-d)])
{
  // Creating lower n x (p - d) submatrix of x
  int p0 = p - d;
  double x_lower[n*p0];
  for (int i=0; i<n; i++)
    for (int j=d; j<p; j++)
      x_lower[i*p0 + (j-d)] = x[i*p + j];

  for (int i=0; i<d*p0; i++)
    (*y_out)[i] = 0;

  // Calculating eq 2.4 in paper excluding the summation over i
  for (int i=0; i<n; i++)
    {
      if (treated[i])
        {
          double x0[d];
          for (int j=0; j<d; j++)
            x0[j] = x_beta[i*d + j];

          double m, dm[d];
          imp_local_linear_regression ( n,
                                        d,
                                        x_beta,
                                        x0,
                                        y,
                                        treated,
                                        kernel_spec,
                                        h11,
                                        h12,
                                        h13,
                                        h14,
                                        gauss_cutoff,
                                       &m,
                                       &dm);

          double x_distance[n];
          for (int j=0; j<n; j++)
            {
              double delta_x = 0;
              for (int k=0; k<d; k++)
                delta_x += pow (x_beta[j*d + k] - x0[k], 2);
              x_distance[j] = sqrt (delta_x);
            }

          double x_lower_estimate[p0];
          nw_kernel_regress ( n,
                              p0,
                              x_distance,
                              x_lower,
                              kernel_spec,
                              h0,
                              gauss_cutoff,
                             &x_lower_estimate);

          if (d == 1)
            {
              double tmp[p0];
              for (int j=0; j<p0; j++)
                tmp[j] = dm[0] * (x_lower[i*p0 + j] - x_lower_estimate[j]);

              // Row of eq 2.4 before summation
              for (int j=0; j<p0; j++)
                (*y_out)[j] += (y[i] - m) * tmp[j];
            }
          else
            {
              double tmp1[p0], tmp2[d*p0];
              for (int j=0; j<p0; j++)
                tmp1[j] = x_lower[i*p0 + j] - x_lower_estimate[j];

              matrix_multiplication ( d,
                                      1,
                                      p0,
                                      dm,
                                      tmp1,
                                     &tmp2);

              for (int j=0; j<d*p0; j++)
                (*y_out)[j] += (y[i] - m) * tmp2[j];

            }
        }
    }

  return;
}

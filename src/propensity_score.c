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
#include <stdio.h>
#include <math.h>
#include "kernel_functions.h"
#include "matrix_utilities.h"
#include "propensity_score.h"

void
propensity_score (int            n,
                  int            d,
                  const double   x_beta[n*d],
                  const double   x0[d],
                  const int      treated[n],
                  int            kernel_spec,
                  double         h1,
                  double        *pr,
                  double       (*d_pr)[d])
{
  int d1 = d + 1, d12 = (d + 1) * (d + 1);

  // Looping to compute the summations in Eq 2.11, where u1 is a vector
  // u1 = [1, dx] and uu is a matrix Am = sum([1, dx; dx, dx^2]), where ;
  // signifies a row break and the rows are the rows in the system of equations
  // in eq 2.11.
  double u1[n*d1], w[n], Am[d12];
  for (int i=0; i<d12; i++)
    Am[i] = 0;
  for (int i=0; i<n; i++)
    {
      double u_tmp[d1], u_norm2 = 0;
      u1[i] = 1;
      u_tmp[0] = 1;
      for (int j=0; j<d; j++)
        {
          double dx = x_beta[i*d + j] - x0[j];
          u1[i + (j+1)*n] = dx;
          u_tmp[j + 1] = dx;
          u_norm2 += pow (dx, 2);
        }

      double uu_tmp[d12];
      vector_outer_product ( d1,
                             u_tmp,
                             u_tmp,
                            &uu_tmp);

      // Euclidean distance for kernel evaluation
      double u = sqrt (u_norm2);
      w[i] = kernel_evaluator (u,
                               kernel_spec,
                               h1,
                               1e-3);

      for (int j=0; j<d12; j++)
        Am[j] += w[i] * uu_tmp[j];
    }

  double Am_inv[d12];
  if (d == 1)
    matrix_2d_inversion (Am, &Am_inv);
  else if (d == 2)
    matrix_3d_inversion (Am, &Am_inv);

  double bv[d1];
  for (int i=0; i<d1; i++)
    bv[i] = 0;
  for (int i=0; i<n; i++)
    if (treated[i])
      for (int j=0; j<d1; j++)
        bv[j] += w[i] * u1[i + j*n];

  double f[d1];
  matrix_multiplication ( d1,
                          d1,
                          1,
                          Am_inv,
						  bv,
                         &f);

  // Checks if probability is okay, otherwise sets to nan
  if ((f[0] >= 0.99) ||
      (f[0] <= 0.01))
    {
      *pr = NAN;
    }
  else
    {
      *pr = f[0];
    }

  for (int i=0; i<d; i++)
    (*d_pr)[i] = f[i + 1];

  return;
}

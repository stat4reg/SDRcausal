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
#include <string.h>
#include "omp_test.h"
#include "nw_kernel_regress.h"
#include "matrix_utilities.h"
#include "propensity_score.h"
#include "ipw_dim_red.h"

void
ipw_dim_red (int            n,
             int            p,
             int            d,
             const double   x[n*p],
             const double   x_alpha[n*d],
             const int      treated[n],
             int            kernel_spec,
             double         h1,
             double         h2,
             double         gauss_cutoff,
             int            n_threads,
             int           *n_pr_nan,
             double       (*col_sum)[(p-d)*d])
{
  // Creating lower n x (p - d) submatrix of x
  int p0 = p - d;
  double x_lower[n*p0];
  for (int i=0; i<n; i++)
    for (int j=0; j<p0; j++)
      x_lower[i*p0 + j] = x[i*p + (j+d)];

  int n_nan_tmp = 0;
  double col_sum_tmp[d*p0];
  for (int i=0; i<d*p0; i++)
    col_sum_tmp[i] = 0;
  #if defined(_OPENMP)
    #pragma omp parallel for schedule (static) num_threads (n_threads)
  #endif
  for (int i=0; i<n; i++)
    {
      double x0[d];
      for (int j=0; j<d; j++)
        x0[j] = x_alpha[i*d + j];

      double pr, d_pr[d], d_eta_est[d];
      propensity_score ( n,
                         d,
                         x_alpha,
                         x0,
                         treated,
                         kernel_spec,
                         h1,
                        &pr,
                        &d_pr);

      // Checks if probability is outside of range. In that case the row is
      // ignored but the optimization will be punished in the parent script.
      if (isnan (pr))
        {
          #if defined(_OPENMP)
            #pragma omp atomic
          #endif
          n_nan_tmp++;
          continue;
        }
      else
        {
          for (int j=0; j<d; j++)
            d_eta_est[j] = d_pr[j] / (pr * (1 - pr));
        }

      double x_distance[n*d];
      for (int j=0; j<n; j++)
        for (int k=0; k<d; k++)
          x_distance[j*d + k] = x_alpha[j*d + k] - x0[k];

      // Estimates E(X_i | X * alpha)
      double e_hat[p0];
      nw_kernel_regress ( n,
                          p0,
                          x_distance,
                          x_lower,
                          kernel_spec,
                          h2,
                          gauss_cutoff,
                         &e_hat);

      // Calculating the n:th row of the estimating equation
      for (int j=0; j<p0; j++)
        {
          // 1d implementation
          #if defined(_OPENMP)
            #pragma omp atomic
          #endif
          col_sum_tmp[j] += (x_lower[i*p0 + j] - e_hat[j])
            * ( (double) treated[i] - pr) * d_eta_est[0];
        }
    }

  *n_pr_nan = n_nan_tmp;
  memmove (col_sum, col_sum_tmp, d*p0*sizeof (double));

  return;
}

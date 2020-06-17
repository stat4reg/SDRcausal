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
#include <math.h>
#include <R.h>
#include <Rdefines.h>
#include "matrix_utilities.h"
#include "ipw_interpolate.h"
#include "propensity_score.h"
#include "C_propensity_score.h"

SEXP
C_propensity_score (SEXP       n_in,
                    SEXP       d_in,
                    const SEXP x_beta_in,
                    const SEXP treated,
                    SEXP       kernel_spec,
                    SEXP       h1,
                    SEXP       pr_out,
                    SEXP       d_pr_out)
{
  // Storing in local variables for more convenient code, escaping the use of
  // REAL and INTEGER
  int n = *INTEGER (n_in), d = *INTEGER (d_in);

  // Creates pointer to x_beta_in to escape the use of REAL in every iteration.
  double *x_beta = REAL (x_beta_in);

  double pr[n], d_pr[n*d];
  for (int i=0; i<n; i++)
    {
      double x0[d];
      for (int j=0; j<d; j++)
        x0[j] = x_beta[i*d + j];

      double d_pr_tmp[d];
      propensity_score ( n,
                         d,
                         x_beta,
                         x0,
                         INTEGER (treated),
                        *INTEGER (kernel_spec),
                        *REAL (h1),
                        &pr[i],
                        &d_pr_tmp);

      for (int j=0; j<d; j++)
        d_pr[i*d + j] = d_pr_tmp[j];
    }

  ipw_interpolate ( n,
                    d,
                    x_beta,
                   &pr , &d_pr);

  memmove (REAL (pr_out), pr, n*sizeof (double));
  memmove (REAL (d_pr_out), d_pr, n*d*sizeof (double));

  return R_NilValue;
}

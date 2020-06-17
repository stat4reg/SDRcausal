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
#include "nw_kernel_regress.h"
#include "C_nw_kernel_regress.h"

SEXP
C_nw_kernel_regress (SEXP       n_in,
                     SEXP       p_in,
                     SEXP       m_in,
                     const SEXP u,
                     const SEXP y,
                     SEXP       kernel_spec,
                     SEXP       bandwidth,
                     SEXP       gauss_cutoff,
                     SEXP       k)
{
  // Storing in local variables for more convenient code, escaping the use of
  // INTEGER
  int n = *INTEGER (n_in), m = *INTEGER (m_in), p = *INTEGER (p_in);
  double *x_ptr = REAL (u);

  double k_out[n*m];
  for (int i=0; i<n; i++)
    {
      double x0[p];
      for (int j=0; j<p; j++)
        x0[j] = x_ptr[i*p + j];

      double u[n];
      for (int j=0; j<n; j++)
        {
          double delta_x = 0;
          for (int k=0; k<p; k++)
            delta_x += pow (x_ptr[j*p + k] - x0[k], 2);
          u[j] = sqrt (delta_x);
        }

      double k_tmp[m];
      nw_kernel_regress ( n,
                          m,
                          u,
                          REAL (y),
                         *INTEGER (kernel_spec),
                         *REAL (bandwidth),
                         *REAL (gauss_cutoff),
                         &k_tmp);

      for (int j=0; j<m; j++)
        k_out[i*m + j] = k_tmp[j];
    }

  memmove (REAL (k), k_out, n*m*sizeof (double));

  return R_NilValue;
}

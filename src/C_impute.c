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
#include <string.h>
#include <math.h>
#include <R.h>
#include <Rdefines.h>
#include "matrix_utilities.h"
#include "impute.h"
#include "imp_extrapolate.h"
#include "imp_truncate.h"
#include "C_impute.h"

SEXP
C_impute (SEXP       n_in,
          SEXP       d_in,
          const SEXP x_beta,
          const SEXP y,
          const SEXP treated,
          SEXP       kernel_spec,
          SEXP       h11,
          SEXP       h12,
          SEXP       h13,
          SEXP       h14,
          SEXP       gauss_cutoff,
          SEXP       to_extrapolate,
          SEXP       to_truncate,
          SEXP       extrapolation_basis,
          SEXP       m,
          SEXP       dm)
{
  // Storing in local variables for more convenient code, escaping the use of
  // INTEGER
  int n = *INTEGER (n_in), d = *INTEGER (d_in);

  // Creating temporary variables to be used inside C
  double m_tmp[n], dm_tmp[n*d];

  // Performing the local linear regression to get the imputed values
  // int n_extrapolated = 0, n_truncated = 0;
  impute ( n,
           d,
           REAL (x_beta),
           REAL (y),
           INTEGER (treated),
          *INTEGER (kernel_spec),
          *REAL (h11),
          *REAL (h12),
          *REAL (h13),
          *REAL (h14),
          *REAL (gauss_cutoff),
          &m_tmp,
          &dm_tmp);

  // Extrapolating the edges of the imputed values, if needed.
  if (*INTEGER (to_extrapolate))
    {
      int n_extrapolated = 0;
      imp_extrapolate ( n,
                       d,
                       REAL (x_beta),
                       INTEGER (treated),
                       *INTEGER (extrapolation_basis),
                       &m_tmp,
                       &dm_tmp,
                       &n_extrapolated);
    }

  // Truncating the imputed values to within the range of y.
  if (*INTEGER (to_truncate))
    {
      int n_truncated = 0;
      imp_truncate ( n,
                    REAL (y),
                    &m_tmp,
                    &n_truncated);
    }

  memmove (REAL (m), m_tmp, n*sizeof (double));
  memmove (REAL (dm), dm_tmp, n*d*sizeof (double));

  return R_NilValue;
}

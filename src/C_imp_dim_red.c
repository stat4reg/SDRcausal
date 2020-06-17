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
#include "omp_test.h"
#include "matrix_utilities.h"
#include "imp_dim_red.h"
#include "C_imp_dim_red.h"

SEXP
C_imp_dim_red (SEXP       n_in,
               SEXP       p_in,
               SEXP       d_in,
               const SEXP x,
               const SEXP x_beta,
               const SEXP y,
               const SEXP treated,
               SEXP       kernel_spec,
               SEXP       h0,
               SEXP       h11,
               SEXP       h12,
               SEXP       h13,
               SEXP       h14,
               SEXP       gauss_cutoff,
               SEXP       penalty,
               SEXP       n_pen_in,
               SEXP       n_threads_in,
               SEXP       fval)
{
  // Storing in local variables for more convenient code, escaping the use of
  // REAL and INTEGER
  int n = *INTEGER (n_in), p = *INTEGER (p_in), d = *INTEGER (d_in);
  int n_pen = *INTEGER (n_pen_in);

  // Handles number of threads
  int max_threads = omp_get_max_threads(), n_threads = *INTEGER (n_threads_in);
  if (n_threads > 1)
    {
      OMPMSG(1);
    }
  if (n_threads <= 0)
    n_threads = max_threads;

  if (n_threads > max_threads)
  {
      Rprintf("\nWarning: n_threads (%d) exceeds maximum number of threads"
              "available (%d). \nUsing maximum number of threads - 1 (%d)."
              "\nTo use maximum number of threads, set n_threads to %d.\n\n",
              n_threads, max_threads, max_threads, max_threads);

      n_threads = max_threads;
  }

  int p0 = p - d, n_llr_fail = 0;
  double y_out[d*p0];
  imp_dim_red ( n,
                p,
                d,
                REAL (x),
                REAL (x_beta),
                REAL (y),
                INTEGER (treated),
               *INTEGER (kernel_spec),
               *REAL (h0),
               *REAL (h11),
               *REAL (h12),
               *REAL (h13),
               *REAL (h14),
               *REAL (gauss_cutoff),
                n_threads,
               &n_llr_fail,
               &y_out);

  // Norm of y_out
  double sum = 0;
  for (int i=0; i<p0; i++)
    sum += pow(y_out[i], 2);

  // Penalising for local linear approximation fails
  if (n_llr_fail > n_pen)
    sum = sum + pow(*REAL (penalty), n_llr_fail - n_pen);
  *REAL (fval) = sum;

  return R_NilValue;
}

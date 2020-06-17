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
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "matrix_utilities.h"
#include "imp_llr.h"
#include "impute.h"

void
impute (int              n,
        int              d,
        const double     x_beta[n*d],
        const double     y[n],
        const int        treated[n],
        int              kernel_spec,
        double           h11,
        double           h12,
        double           h13,
        double           h14,
        double           gauss_cutoff,
        double         (*m)[n],
        double         (*dm)[n*d])
{
  // Calculating the imputed values
  for (int i=0; i<n; i++)
    {
      double x0[d], dm_tmp[d];

      for (int j=0; j<d; j++)
        x0[j] = (x_beta)[i*d + j];

      imp_llr ( n,
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
               &(*m)[i],
               &dm_tmp);

      for (int j=0; j<d; j++)
        {
          if (isnan (dm_tmp[j]))
            (*dm)[i*d + j] = 0;
          else
            (*dm)[i*d + j] = dm_tmp[j];
        }

    }

  return;
}

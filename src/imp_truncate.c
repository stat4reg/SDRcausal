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
#include "imp_truncate.h"

void
imp_truncate (int            n,
              const double   y[n],
              double       (*m)[n],
              int           *n_truncated)
{
  // Identifying y_min and y_max of y
  double y_min = INFINITY, y_max = -INFINITY;
  for (int i=0; i<n; i++)
    {
      if (y[i] < y_min)
        y_min = y[i];
      if (y[i] > y_max)
        y_max = y[i];
    }

  // Truncating values
  for (int i=0; i<n; i++)
    {
      if ((*m)[i] < y_min)
        {
          // Truncating value
          (*m)[i] = y_min;

          (*n_truncated)++;
        }
      if ((*m)[i] > y_max)
        {
          // Truncating value
          (*m)[i] = y_max;

          (*n_truncated)++;
        }
    }

  return;
}

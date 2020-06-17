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
#include "index_element.h"
#include "imp_extrapolate.h"

void
imp_extrapolate (int            n,
                 int            d,
                 const double   x_beta[n*d],
                 const int      treated[n],
                 int            extrapolation_basis,
                 double       (*m)[n],
                 double       (*dm)[n*d],
                 int           *n_extrapolated)
{
  // Creating vector of index elements to keep track of original indices
  IndexElement xb_sorted[n*d];
  for (int i=0; i<n*d; i++)
    {
      xb_sorted[i].index = i;
      xb_sorted[i].value = x_beta[i];
    }

  // Sorting x_beta
  qsort (xb_sorted,
         (size_t) n*d,
         sizeof (IndexElement),
         compare_IndexElement);

  // Creating temporary vectors of m, dm, and treated that matches sorted x_beta
  int treated_tmp[n];
  double m_tmp[n];
  for (int i=0; i<n; i++)
    {
      m_tmp[i] = (*m)[xb_sorted[i].index];
      treated_tmp[i] = treated[xb_sorted[i].index];
    }

  // Finding values to extrapolate and counts number of values to be
  // extrapolated
  int idx_left = 0, idx_right = n - 1, n_to_extrapolate = 0;
  while ((idx_left < n-1) &&
         (!treated_tmp[idx_left]))
    {
      idx_left++;
      n_to_extrapolate++;
    }
  while ((idx_right > 0) &&
         (!treated_tmp[idx_right]))
    {
      idx_right--;
      n_to_extrapolate++;
    }

  // Exits if no values to extrapolate
  if (n_to_extrapolate == 0)
    {
      *n_extrapolated = 0;
      return;
    }

    // Extrapolating left side
    {
      int idx_step = idx_left + 1, n_based = 0;
      double x, y, dydx = 0;
      x = xb_sorted[idx_left].value;
      y = m_tmp[idx_left];

      // Basing extrapolation on treated values only
      while (n_based < extrapolation_basis)
        {
          if (!treated_tmp[idx_step])
            {
              // Adding current points slope
              double dx, dy;
              dx = xb_sorted[idx_step].value - x;
              dy = m_tmp[idx_step] - y;

              if (dx == 0)
                dydx += 0;
              else
                dydx += dy / dx;

              n_based++;
            }
          idx_step++;
        }

      // Averaging the derivative
      dydx = dydx / extrapolation_basis;

      for (int i=0; i<idx_left; i++)
        {
          // Linear extrapolation based on dm
          double y_new = y + dydx * (xb_sorted[i].value - x);
          (*m)[xb_sorted[i].index] = y_new;

          // Setting derivative
          (*dm)[xb_sorted[i].index] = dydx;

        }
    }

    // Extrapolating right side
    {
      int idx_step = idx_right - 1, n_based = 0;
      double x, y, dydx = 0;
      x = xb_sorted[idx_right].value;
      y = m_tmp[idx_right];

      // Basing extrapolation on treated values only
      while (n_based < extrapolation_basis)
        {
          if (!treated_tmp[idx_step])
            {
              // Adding current point slope
              double dx, dy;
              dx = xb_sorted[idx_step].value - x;
              dy = m_tmp[idx_step] - y;

              if (dx == 0)
                dydx += 0;
              else
                dydx += dy / dx;

              n_based++;
            }
          idx_step--;
        }

      // Averaging the derivative
      dydx = dydx / (double) extrapolation_basis;

      for (int i=n-1; i>idx_right; i--)
        {
          // Linear extrapolation
          double y_new = y + dydx * (xb_sorted[i].value - x);
          (*m)[xb_sorted[i].index] = y_new;

          // Setting derivative
          (*dm)[xb_sorted[i].index] = dydx;

        }
    }

  *n_extrapolated = n_to_extrapolate;

  return;
}

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
#include "ipw_interpolate.h"

void
ipw_interpolate (int            n,
                 int            d,
                 const double   x_alpha[n*d],
                 double       (*pr)[n],
				 double       (*d_pr)[n*d])
{
  int extrapolation_basis = 3;

  // Creating vector of index elements to keep track of original indices
  IndexElement xa_sorted[n*d];
  for (int i=0; i<n*d; i++)
    {
      xa_sorted[i].index = i;
      xa_sorted[i].value = x_alpha[i];
    }

  // Sorting x_alpha
  qsort (xa_sorted,
         (size_t) n*d,
         sizeof (IndexElement),
         compare_IndexElement);

  // Creating temporary vector of pr that matches sorted x_alpha
  double pr_sorted[n];
  for (int i=0; i<n; i++)
    pr_sorted[i] = (*pr)[xa_sorted[i].index];
	
  double d_pr_sorted[n*d];
  for (int i=0; i<n; i++)
  {
    for (int k=0; k<d; k++)
		d_pr_sorted[i+k] = (*d_pr)[xa_sorted[i].index + k];
  }
  // Extrapolating
  int idx_left = 0, idx_right = n - 1, n_to_extrapolate = 0;
  while ((idx_left < n-1) &&
         (isnan (pr_sorted[idx_left])))
    {
      idx_left++;
      n_to_extrapolate++;
    }
  while ((idx_right > 0) &&
         (isnan (pr_sorted[idx_right])))
    {
      idx_right--;
      n_to_extrapolate++;
    }

    // Extrapolating left side
    {
      int idx_step = idx_left + 1, n_based = 0;
      double x, y, dydx = 0;
      x = xa_sorted[idx_left].value;
      y = pr_sorted[idx_left];
      // Creating a linear fit based on the extrapolation_basis leftmost points
      while (n_based < extrapolation_basis)
        {
          // Uses only points with values within (0, 1)
          if (isnan (pr_sorted[idx_step]))
            {
              // Skips point if not in (0, 1)
              idx_step++;
            }
          else
            {
              // Adding current points slope
              double dx, dy;
              dx = xa_sorted[idx_step].value - x;
              dy = pr_sorted[idx_step] - y;

              if (dx == 0)
                dydx += 0;
              else
                dydx += dy / dx;

              n_based++;
              idx_step++;
            }
        }
      dydx = dydx / (double) extrapolation_basis;

      for (int i=0; i<idx_left; i++)
        {
          // Calculating extrapolated value
          double y_new = y + dydx * (xa_sorted[i].value - x);
		  
          for (int k=0; k<d; k++)
		  	d_pr_sorted[i+k] = dydx;  
          // Checks that extrapolated value is within range
          if (y_new <= 0.01 ||
              y_new >= 0.99)
            {
              // If extrapolation gives a value outside (0, 1) the value of the
              // outmost ok point gets assigned to y_new
              y_new = y;
			  for (int k=0; k<d; k++)
		  		d_pr_sorted[i+k] = 0;  
            }

          pr_sorted[i] = y_new;
        }
    }

    // Extrapolating right side
    {
      int idx_step = idx_right - 1;
      int n_based = 0;
      double x, y, dydx = 0;
      x = xa_sorted[idx_right].value;
      y = pr_sorted[idx_right];
      // Creating a linear fit based on the extrapolation_basis leftmost points

      while (n_based < extrapolation_basis)
        {
          // Uses only points with values within (0, 1)
          if (isnan (pr_sorted[idx_step]))
            {
              // Skips point if not in (0, 1)
              idx_step--;
            }
          else
            {
              // Adding current point slope
              double dx, dy;
              dx = xa_sorted[idx_step].value - x;
              dy = pr_sorted[idx_step] - y;

              if (dx == 0)
                dydx += 0;
              else
                dydx += dy / dx;

              n_based++;
              idx_step--;
            }
        }
      dydx = dydx / (double) extrapolation_basis;

      for (int i=n-1; i>idx_right; i--)
        {
          // Calculating extrapolated value
          double y_new = y + dydx * (xa_sorted[i].value - x);

          for (int k=0; k<d; k++)
		  	d_pr_sorted[i+k] = dydx; 
          // Checks that extrapolated value is within range
          if (y_new <= 0.01 ||
              y_new >= 0.99)
            {
              // If extrapolation gives a value outside (0, 1) the value of the
              // outmost ok point gets assigned to y_new
              y_new = y;
			  
			  for (int k=0; k<d; k++)
		  		d_pr_sorted[i+k] = 0; 
            }

          pr_sorted[i] = y_new;
        }
    }

  // Interpolating
  int n_interpolated = 0;
  for (int i=1; i<n-1; i++)
    {
      if (isnan (pr_sorted[i]))
        {
          double sum = 0;
		  double dy = 0;
		  double dx = 0;
          sum += pr_sorted[i - 1];
		  dx = xa_sorted[i - 1].value;
          dy = pr_sorted[i - 1];
		  
          int right_done = 0, idr = 1;
          while (!right_done)
            {
              // Checks so next value is not also outside the boundary. In that
              // case moves on to the next value.
              if (isnan (pr_sorted[i+idr]))
                {
                  idr++;
                }
              else
                {
                  sum += pr_sorted[i+idr];
                  dx = dx - xa_sorted[i+idr].value;
                  dy = dy - pr_sorted[i+idr];
				  right_done = 1;
                }
            }
          pr_sorted[i] = sum / 2;
		  for (int k=0; k<d; k++)
		  	d_pr_sorted[i+k] = dy/dx; 
          n_interpolated++;
        }
    }

  for (int i=0; i<n; i++){
    (*pr)[xa_sorted[i].index] = pr_sorted[i];
    for (int k=0; k<d; k++)
			(*d_pr)[xa_sorted[i+k].index] = d_pr_sorted[i+k];
 }
  return;
}

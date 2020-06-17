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
#include <stdio.h>
#include "kernel_functions.h"

double
epanechnikov_kernel (double u,
                     double h)
{
  u = u / h;
  double k = 0;
  if (fabs(u) < 1)
    k = (double) 0.75 * ( (double) 1.0 - pow(u, 2)) / h;

  return k;
}

double
quartic_biweight_kernel (double u,
                         double h)
{
  u = u / h;
  double k = 0;
  if (fabs(u) < 1)
    k = (double) 0.9375 * pow((double) 1 - pow(u, 2), 2) / h;

  return k;
}

double
gaussian_kernel (double u,
                 double h,
                 double c)
{
  double k = 0, b = 2.506628274631000; // \approx sqrt((double) 2 * M_PI)
  u = u / h;

  k = exp((double) -0.5 * pow(u, 2)) / b;

  if (k < c)
    return (double) 0;
  else
    return k;
}

double
kernel_evaluator (double u,
                  int    kernel_spec,
                  double h,
                  double c)
{
  double k = 0;
  switch (kernel_spec)
    {
    case 1:
      k = epanechnikov_kernel (u, h);
      break;
    case 2:
      k = quartic_biweight_kernel (u, h);
      break;
    case 3:
      k = gaussian_kernel (u, h, c);
      break;
    }

  return k;
}

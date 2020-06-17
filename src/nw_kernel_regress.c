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
#include "kernel_functions.h"
#include "nw_kernel_regress.h"

void
nw_kernel_regress (int            n,
                   int            m,
                   const double   u[n],
                   const double   y[n*m],
                   int            kernel_spec,
                   double         bandwidth,
                   double         gauss_cutoff,
                   double       (*k)[m])
{

  double denominator = 0, nominator[m];
  for (int i=0; i<m; i++)
    nominator[i] = 0;

  // Calculating the estimation
  for (int i=0; i<n; i++)
    {
      double s = 0;
      s = kernel_evaluator (u[i],
                            kernel_spec,
                            bandwidth,
                            gauss_cutoff);

      denominator += s;
      for (int j=0; j<m; j++)
        nominator[j] += y[i*m + j] * s;
    }

  // Checking if division by 0, sets value to 0 in that case.
  if (denominator == 0)
    for (int i=0; i<m; i++)
      (*k)[i] = 0;
  else
    for (int i=0; i<m; i++)
      (*k)[i] = nominator[i] / denominator;

  return;
}

void
nw_kernel_regress_treatment (int            n,
                             int            m,
                             const double   u[n],
                             const double   y[n*m],
                             const int      treated[n],
                             int            kernel_spec,
                             double         bandwidth,
                             double         gauss_cutoff,
                             double       (*k)[m])
{

  double denominator = 0, nominator[m];
  for (int i=0; i<m; i++)
    nominator[i] = 0;

  // Calculating the estimation
  for (int i=0; i<n; i++)
    {
      double s = 0;
      // Using only treated values.
      if (treated[i])
        s = kernel_evaluator (u[i],
                              kernel_spec,
                              bandwidth,
                              gauss_cutoff);

      denominator += s;
      for (int j=0; j<m; j++)
        nominator[j] += y[i*m + j] * s;
    }

  // Checking if division by 0, sets value to 0 in that case.
  if (denominator == 0)
    for (int i=0; i<m; i++)
      (*k)[i] = 0;
  else
    for (int i=0; i<m; i++)
      (*k)[i] = nominator[i] / denominator;

  return;
}

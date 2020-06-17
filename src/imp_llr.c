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
#include <R.h>
#include "matrix_utilities.h"
#include "nw_kernel_regress.h"
#include "imp_llr.h"

void
imp_llr (int            n,
         int            d,
         const double   x[n*d],
         const double   x0[d],
         const double   y[n],
         const int      treated[n],
         int            kernel_spec,
         double         h11,
         double         h12,
         double         h13,
         double         h14,
         double         gauss_cutoff,
         double        *m,
         double       (*dm)[d])
{
  double u[n], y11[n], y12[n*d], y13[n*d], y14[n*d*d];

  for (int i=0; i<n; i++)
    {
      y11[i] = y[i];
      double u_tmp[d], u_tmp_sum = 0;
      for (int j=0; j<d; j++)
        {
          u_tmp_sum += pow (x[i*d + j] - x0[j], 2);

          u_tmp[j] = x[i*d + j] - x0[j];
          y12[i*d + j] = y[i]*u_tmp[j];
          y13[i*d + j] = u_tmp[j];
        }
      // Euclidean distance for kernel evaluation
      u[i] = sqrt (u_tmp_sum);

      // This needs to be checked
      for (int j=0; j<d; j++)
        for (int k=0; k<d; k++)
          y14[i*d*d + j*d + k] = u_tmp[j] * u_tmp[k];
    }

  double a11[1], a12[d], a13[d], a14[d*d];
  nw_kernel_regress_treatment ( n,
                                1,
                                u,
                                y11,
                                treated,
                                kernel_spec,
                                h11,
                                gauss_cutoff,
                               &a11);

  nw_kernel_regress_treatment ( n,
                                d,
                                u,
                                y12,
                                treated,
                                kernel_spec,
                                h12,
                                gauss_cutoff,
                               &a12);

  nw_kernel_regress_treatment ( n,
                                d,
                                u,
                                y13,
                                treated,
                                kernel_spec,
                                h13,
                                gauss_cutoff,
                               &a13);

  nw_kernel_regress_treatment ( n,
                                d*d,
                                u,
                                y14,
                                treated,
                                kernel_spec,
                                h14,
                                gauss_cutoff,
                               &a14);

  // Calculating m and dm according to the following formulas. The
  // corresponding steps in the 2, 3 d code are shown in the parentheses
  // below the expression as step 1 = s1 and so on.
  // m  = A_11 - A_13^T * (A_14 - A_13 * A_13^T)^(-1) * (A_12 - A_11 * A_13)
  //                              (     s1    )
  //                      (          s2        )
  //                      (          s3             )
  //             (             s4                   )
  //                                                            (    s5    )
  //                                                    (        s6        )
  //             (                           s7                            )
  // dm = (A_14 - A_13 * A_13^T)^(-1) * (A_12 - A_11 * A_13)
  //      (           s3            )
  //                                    (        s6        )
  //
  if (d == 1)
    {
      double m_a, m_b, dm_a, dm_b;
      m_a = a13[0] * (a12[0] - a11[0] * a13[0]);
      m_b = a14[0] - a13[0] * a13[0];

      // Checks if denominator is 0 and in that case uses the value of a11[0]
      // which is the same as the funciton value for treated and 0 for not
      // treated.
      if (m_b == 0)
        *m = a11[0] - m_a / (m_b+ 1e-04);
      else
        *m = a11[0] - m_a / (m_b+ 1e-04);

      dm_a = a12[0] - a11[0] * a13[0];
      dm_b = a14[0] - a13[0] * a13[0];
      // Checks if denominator is 0 and in that case sets the derivative to 0
      // Should set to large number instead of zero???
      if (dm_b == 0)
        (*dm)[0] = dm_a / (dm_b + 1e-04);
      else
        (*dm)[0] = dm_a / (dm_b + 1e-04);
    }
  else
    {
      double s1[d*d], s2[d*d], s3[d*d], s4[d], s5[d], s6[d], s7;

      vector_outer_product ( d,
                             a13,
                             a13,
                            &s1);

      matrix_subtraction ( d,
                           d,
                           a14,
                           s1,
                          &s2);

      if (d == 2)
        {
          int inversion_fail = matrix_2d_inversion (s2, &s3);
          if (inversion_fail)
            {
              *m = a11[0];
              for (int i=0; i<d; i++)
                (*dm)[i] = 0;

              return;
            }
        }
      else
        {
          int inversion_fail = matrix_3d_inversion (s2, &s3);
          if (inversion_fail)
            {
              *m = a11[0];
              for (int i=0; i<d; i++)
                (*dm)[i] = 0;

              return;
            }
        }

      matrix_multiplication ( 1,
                              d,
                              d,
                              a13,
                              s3,
                             &s4);

      matrix_scaling ( 1,
                       d,
                       a11[0],
                       a13,
                      &s5);

      matrix_subtraction ( 1,
                           d,
                           a12,
                           s5,
                          &s6);

      vector_scalar_product ( d,
                              s4,
                              s6,
                             &s7);

      // Assigning m value
      *m = a11[0] - s7;

      // Assigning dm value
      matrix_multiplication(d,
                            d,
                            1,
                            s3,
                            s6,
                            dm);
    }

  return;
}

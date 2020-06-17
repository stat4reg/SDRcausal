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

void
vector_scalar_product(int          n,
                      const double  vec_a[n],
                      const double  vec_b[n],
                      double       *out)
{
  *out = 0;
  for (int i=0; i<n; i++)
    *out += vec_a[i] * vec_b[i];

  return;
}


void
vector_outer_product(int          n,
                     const double vec_a[n],
                     const double vec_b[n],
                     double       (*mat_out)[n*n])
{
  for (int i=0; i<n; i++)
    for (int j=0; j<n; j++)
      (*mat_out)[i*n + j] = vec_a[i] * vec_b[j];

  return;
}

void
vector_3d_cross_product(const double vec_a[3],
                        const double vec_b[3],
                        double       (*vec_out)[3])
{
  (*vec_out)[0] = vec_a[1] * vec_b[2] - vec_a[2] * vec_b[1];
  (*vec_out)[1] = vec_a[2] * vec_b[0] - vec_a[0] * vec_b[2];
  (*vec_out)[2] = vec_a[0] * vec_b[1] - vec_a[1] * vec_b[0];

  return;
}

void
matrix_scaling(int          n,
               int          m,
               double       f,
               const double mat[n*m],
               double       (*mat_out)[n*m])
{
  for (int i=0; i<n*m; i++)
    (*mat_out)[i] = f * mat[i];

  return;
}

void
matrix_addition(int          n,
                int          m,
                const double mat_a[n*m],
                const double mat_b[n*m],
                double       (*mat_out)[n*m])
{
  for (int i=0; i<n*m; i++)
    (*mat_out)[i] = mat_a[i] + mat_b[i];

  return;
}

void
matrix_subtraction(int          n,
                   int          m,
                   const double mat_a[n*m],
                   const double mat_b[n*m],
                   double       (*mat_out)[n*m])
{
  for (int i=0; i<n*m; i++)
    (*mat_out)[i] = mat_a[i] - mat_b[i];

  return;
}

void
matrix_elementwise_sum(int           n,
                       int           m,
                       const double  mat[n*m],
                       double       *sum)
{
  *sum = 0;
  for (int i=0; i<n*m; i++)
    *sum += mat[i];

  return;
}

void
matrix_multiplication(int          n,
                      int          m,
                      int          p,
                      const double mat_a[n*m],
                      const double mat_b[m*p],
                      double       (*mat_out)[n*p])
{
  // Zeroing all values of output matrix
  for (int i=0; i<n*p; i++)
    (*mat_out)[i] = 0;

  for (int i=0; i<n; i++)
    for (int j=0; j<p; j++)
      for (int k=0; k<m; k++)
        (*mat_out)[i*p + j] += mat_a[i*m + k] * mat_b[j + k*p];

  return;
}

int
matrix_2d_inversion(const double mat[4],
                    double       (*mat_out)[4])
{
  // Computing determinant
  double det = mat[0] * mat[3] - mat[1] * mat[2];

  if (det == 0)
    {
      for (int i=0; i<4; i++)
        (*mat_out)[i] = 0;
      return 1;
    }

  (*mat_out)[0] = (1 / det) * mat[3];
  (*mat_out)[1] = (1 / det) * (- mat[1]);
  (*mat_out)[2] = (1 / det) * (- mat[2]);
  (*mat_out)[3] = (1 / det) * mat[0];

  return 0;
}

int
matrix_3d_inversion(const double mat[9],
                    double       (*mat_out)[9])
{

  double a[3], b[3], c[3];
  for (int i=0; i<3; i++)
    {
      a[i] = mat[i*3];
      b[i] = mat[i*3 + 1];
      c[i] = mat[i*3 + 2];
    }

  double bc[3], ca[3], ab[3];
  vector_3d_cross_product(b, c, &bc);
  vector_3d_cross_product(c, a, &ca);
  vector_3d_cross_product(a, b, &ab);

  // Computing determinant by Sarrus rule
  double det = mat[0] * mat[4] * mat[8]
             + mat[1] * mat[5] * mat[6]
             + mat[2] * mat[3] * mat[7]
             - mat[6] * mat[4] * mat[2]
             - mat[7] * mat[5] * mat[0]
             - mat[8] * mat[3] * mat[1];

  // double det = 0;
  // vector_scalar_product(3, a, bc, &det);

  if (det == 0)
    {
      for (int i=0; i<9; i++)
        (*mat_out)[i] = 0;
      return 1;
    }


  for (int i=0; i<3; i++)
    {
      (*mat_out)[i] = (1 / det) * bc[i];
      (*mat_out)[3 + i] = (1 / det) * ca[i];
      (*mat_out)[6 + i] = (1 / det) * ab[i];
    }

  return 0;
}

void
matrix_row_to_col_major(int          n,
                        int          m,
                        const double mat[n*m],
                        double       (*mat_out)[n*m])
{
  for (int i=0; i<m; i++)
    for (int j=0; j<n; j++)
      (*mat_out)[i*n + j] = mat[i + j*m];

  return;
}

void
matrix_col_to_row_major(int          n,
                        int          m,
                        const double mat[n*m],
                        double       (*mat_out)[n*m])
{
  for (int i=0; i<n; i++)
    for (int j=0; j<m; j++)
      (*mat_out)[i*m + j] = mat[i + j*n];

  return;
}

// Print functions for debugging. Commented to not trigger warnings when
// compiling.
// void
// print_array (int           n,
//              int           m,
//              const double *mat,
//              const char   *fmt)
// {
//   for (int i=0; i<n; i++)
//     {
//     printf("\t");
//     for (int j=0; j<m; j++)
//       printf(fmt, mat[i*m + j]);
//
//     printf("\n");
//     }
//
//   return;
// }
//
// void
// print_array_col_major (int           n,
//                        int           m,
//                        const double *mat,
//                        const char   *fmt)
// {
//   for (int i=0; i<n; i++)
//     {
//     printf("\t");
//     for (int j=0; j<m; j++)
//       printf(fmt, mat[i + j*n]);
//
//     printf("\n");
//     }
//
//   return;
// }
//
// void
// print_array_int (int         n,
//                  int         m,
//                  const int  *mat,
//                  const char *fmt)
// {
//   for (int i=0; i<n; i++)
//     {
//     printf("\t");
//     for (int j=0; j<m; j++)
//       printf(fmt, mat[i*m + j]);
//
//     printf("\n");
//     }
//
//   return;
// }

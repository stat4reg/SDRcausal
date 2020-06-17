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
#ifndef MATRIX_UTILITIES_H
#define MATRIX_UTILITIES_H

/**
 * @brief Calculates the scalar product of two vectors.
 *
 * @param n    length of vectors
 * @param vec_a first vector
 * @param vec_b second vector
 * @param out  resulting value
 */
void
vector_scalar_product(int           n,
                      const double  vec_a[n],
                      const double  vec_b[n],
                      double       *out);

/**
 * @brief Calculates the outer product of two vectors.
 *
 * @param n       length of vectors
 * @param vec_a    first vector
 * @param vec_b    second vector
 * @param mat_out resulting matrix
 */
void
vector_outer_product(int          n,
                     const double vec_a[n],
                     const double vec_b[n],
                     double       (*mat_out)[n*n]);

/**
 * @brief Calculates the vector (cross) product of two vectors.
 *
 * @param vec_a    first vector
 * @param vec_b    second vector
 * @param vec_out resulting vector
 */
void
vector_3d_cross_product(const double vec_a[3],
                        const double vec_b[3],
                        double       (*vec_out)[3]);

/**
 * @brief Performs scales a matrix by a factor
 *
 * @param n       first dimension of matrices
 * @param m       second dimension of matrices
 * @param mat     matrix
 * @param f       factor
 * @param mat_out resulting matrix
 */
void
matrix_scaling(int          n,
               int          m,
               double       f,
               const double mat[n*m],
               double       (*mat_out)[n*m]);

/**
 * @brief Performs elementwise addition of two matrices.
 *
 * @param n       first dimension of matrices
 * @param m       second dimension of matrices
 * @param mat_a   first matrix
 * @param mat_b   second matrix
 * @param mat_out resulting matrix
 */
void
matrix_addition(int          n,
                int          m,
                const double mat_a[n*m],
                const double mat_b[n*m],
                double       (*mat_out)[n*m]);

/**
 * @brief Performs elementwise subtraction of two matrices.
 *
 * @param n       first dimension of matrices
 * @param m       second dimension of matrices
 * @param mat_a   first matrix
 * @param mat_b   second matrix
 * @param mat_out resulting matrix
 */
void
matrix_subtraction(int          n,
                   int          m,
                   const double mat_a[n*m],
                   const double mat_b[n*m],
                   double       (*mat_out)[n*m]);

/**
 * @brief Performs elementwise sum of a 2d matrix
 *
 * @param n       first dimension of matrix
 * @param m       second dimension of matrix
 * @param mat     matrix
 * @param sum     resulting sum
 */
void
matrix_elementwise_sum(int           n,
                       int           m,
                       const double  mat[n*m],
                       double       *sum);

/**
 * @brief Performs matrix multiplication. Sets output matrix to 0 before
 * multiplication.
 *
 * @param n       number of rows of mat_a
 * @param m       number of cols of mat_a, number of rows of mat_b
 * @param p       number of cols of mat_b
 * @param mat_a   n x m matrix
 * @param mat_b   m x p matrix
 * @param mat_out resulting n x p matrix
 */
void
matrix_multiplication(int          n,
                      int          m,
                      int          p,
                      const double mat_a[n*m],
                      const double mat_b[m*n],
                      double       (*mat_out)[n*n]);

/**
 * @brief Computes the inverse of a 2d matrix.
 *
 * Computes the inverse of a 2d matrix by straiht forward calculations.
 *
 * If the determinant is 0 the function returns without performing any
 * additional calculations.
 *
 * @param mat     matrix to be inversed
 * @param mat_out resulting matrix
 *
 * @return 0 on success, 1 if inversion fails
 */
int
matrix_2d_inversion(const double mat[4],
                    double       (*mat_out)[4]);

/**
 * @brief Computes the inverse of a 3d matrix.
 *
 * Computes the inverse of a 3d matrix by dividing the matrix in three column
 * vectors which are then used in a triple product to get the determinant and a
 * series of cross products and transposes to generate the inverted matrix.
 *
 * If the determinant is 0 the function returns without performing any
 * additional calculations.
 *
 * @param mat     matrix to be inversed
 * @param mat_out resulting matrix
 *
 * @return 0 on success, 1 if inversion fails
 */
int
matrix_3d_inversion(const double mat[9],
                    double       (*mat_out)[9]);

/**
 * @brief Converts the storage of a flat matrix from row major to col major
 *
 * @param n       first dimension of matrices
 * @param m       second dimension of matrices
 * @param mat     matrix
 * @param mat_out resulting matrix
 */
void
matrix_row_to_col_major(int          n,
                        int          m,
                        const double mat[n*m],
                        double       (*mat_out)[n*m]);

/**
 * @brief Converts the storage of a flat matrix from col major to row major
 *
 * @param n       first dimension of matrices
 * @param m       second dimension of matrices
 * @param mat     matrix
 * @param mat_out resulting matrix
 */
void
matrix_col_to_row_major(int          n,
                        int          m,
                        const double mat[n*m],
                        double       (*mat_out)[n*m]);

/**
 * @brief Prints an array of doubles in a given format
 *
 * @param n   number of rows of mat
 * @param m   number of cols of mat
 * @param mat matrix to be printed
 * @param fmt printing format
 */
void
print_array (int           m,
             int           n,
             const double *mat,
             const char   *fmt);

/**
 * @brief Prints an array of doubles in a given format, column major
 *
 * @param n   number of rows of mat
 * @param m   number of cols of mat
 * @param mat matrix to be printed
 * @param fmt printing format
 */
void
print_array_col_major (int           m,
                       int           n,
                       const double *mat,
                       const char   *fmt);

/**
 * @brief Prints an array of integers in a given format
 *
 * @param n   number of rows of mat
 * @param m   number of cols of mat
 * @param mat matrix to be printed
 * @param fmt printing format
 */
void
print_array_int (int         m,
                 int         n,
                 const int  *mat,
                 const char *fmt);
#endif

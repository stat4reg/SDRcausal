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
#ifndef IMPUTE_H
#define IMPUTE_H

/**
 * @breif Imputes values based on treated outcomes.
 *
 *
 * @param n            Number of observations
 * @param d            Structural dimension
 * @param x_beta       CMS projection of covariate matrix
 * @param y            Response vector
 * @param treated      Binary treatment vector
 * @param beta_hat     Locally efficient estimator of beta
 * @param kernel_spec  Indicates which kernel function to be used
 * @param h11          Kernel bandwidth
 * @param h12          Kernel bandwidth
 * @param h13          Kernel bandwidth
 * @param h14          Kernel bandwidth
 * @param gauss_cutoff Cutoff value for Gaussian kernel
 * @param m            Imputed values
 * @param dm           Derivatives of imputed values
 */
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
        double         (*dm)[n*d]);

#endif

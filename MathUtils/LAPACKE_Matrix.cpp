/* 
 * LAPACKE_Matrix.cpp, part of the RotationShield program
 * Copyright (c) 2007-2014 Michael P. Mendenhall
 *
 * This code uses the LAPACKE C interface to LAPACK;
 * see http://www.netlib.org/lapack/lapacke.html
 * and the GSL interface to CBLAS, https://www.gnu.org/software/gsl/
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

#ifdef WITH_LAPACKE
#include "LAPACKE_Matrix.hh"

template<>
auto Mat_Ops_Real<float>::f_gemm = &cblas_sgemm;
template<>
auto Mat_Ops_Real<double>::f_gemm = &cblas_dgemm;

template<>
auto Mat_Ops_Complex<lapack_complex_float>::f_gemm = &cblas_cgemm;
template<>
auto Mat_Ops_Complex<lapack_complex_double>::f_gemm = &cblas_zgemm;

template<>
auto LAPACKE_Matrix_SVD<double,double>::f_gebrd = &LAPACKE_dgebrd;
template<>
auto LAPACKE_Matrix_SVD<double,double>::f_bdsqr = &LAPACKE_dbdsqr;
template<>
auto LAPACKE_Matrix_SVD<double,double>::f_orgbr = &LAPACKE_dorgbr;
template<>
Mat_Ops<double>* LAPACKE_Matrix_SVD<double,double>::myOps = new Mat_Ops_Real<double>();

template<>
auto LAPACKE_Matrix_SVD<double, lapack_complex_double >::f_gebrd = &LAPACKE_zgebrd;
template<>
auto LAPACKE_Matrix_SVD<double, lapack_complex_double >::f_bdsqr = &LAPACKE_zbdsqr;
template<>
auto LAPACKE_Matrix_SVD<double, lapack_complex_double >::f_orgbr = &LAPACKE_zungbr;
template<>
Mat_Ops<lapack_complex_double>* LAPACKE_Matrix_SVD<double, lapack_complex_double >::myOps = new Mat_Ops_Complex<lapack_complex_double>();
#endif

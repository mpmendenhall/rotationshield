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
void (*Mat_Ops_Real<float>::f_gemm)(const enum CBLAS_ORDER, const enum CBLAS_TRANSPOSE, const enum CBLAS_TRANSPOSE,
                                                  const int, const int, const int, const float, const float*, const int lda,
                                                  const float*, const int, const float, float*, const int) = &cblas_sgemm;
template<>
void (*Mat_Ops_Real<double>::f_gemm)(const enum CBLAS_ORDER, const enum CBLAS_TRANSPOSE, const enum CBLAS_TRANSPOSE,
                                     const int, const int, const int, const double, const double*, const int lda,
                                     const double*, const int, const double, double*, const int) = &cblas_dgemm;
                                                                                                              
template<>
void (*Mat_Ops_Complex<lapack_complex_float>::f_gemm)(const enum CBLAS_ORDER, const enum CBLAS_TRANSPOSE, const enum CBLAS_TRANSPOSE,
                                                      const int, const int, const int, const void*, const void*, const int lda,
                                                      const void*, const int, const void*, void*, const int) = &cblas_cgemm;
template<>
void (*Mat_Ops_Complex<lapack_complex_double>::f_gemm)(const enum CBLAS_ORDER, const enum CBLAS_TRANSPOSE, const enum CBLAS_TRANSPOSE,
                                                       const int, const int, const int, const void*, const void*, const int lda,
                                                       const void*, const int, const void*, void*, const int) = &cblas_zgemm;

                                                
template<>
lapack_int (*LAPACKE_Matrix_SVD<double,double>::f_gebrd)(int, lapack_int, lapack_int, double*, lapack_int, double*, double*, double*, double*) = &LAPACKE_dgebrd;
template<>
lapack_int (*LAPACKE_Matrix_SVD<double,double>::f_bdsqr)(int, char, lapack_int, lapack_int, lapack_int, lapack_int, double*, double*, double*, lapack_int, double*, lapack_int, double*, lapack_int) = &LAPACKE_dbdsqr;
template<>
lapack_int (*LAPACKE_Matrix_SVD<double,double>::f_orgbr)(int, char, lapack_int, lapack_int, lapack_int, double*, lapack_int, const double*) = &LAPACKE_dorgbr;
template<>
Mat_Ops<double>* LAPACKE_Matrix_SVD<double,double>::myOps = new Mat_Ops_Real<double>();

template<>
lapack_int (*LAPACKE_Matrix_SVD<double,lapack_complex_double>::f_gebrd)(int, lapack_int, lapack_int, lapack_complex_double*, lapack_int, double*, double*, lapack_complex_double*, lapack_complex_double*) = &LAPACKE_zgebrd;
template<>
lapack_int (*LAPACKE_Matrix_SVD<double,lapack_complex_double>::f_bdsqr)(int, char, lapack_int, lapack_int, lapack_int, lapack_int, double*, double*, lapack_complex_double*, lapack_int, lapack_complex_double*, lapack_int, lapack_complex_double*, lapack_int) = &LAPACKE_zbdsqr;
template<>
lapack_int (*LAPACKE_Matrix_SVD<double,lapack_complex_double>::f_orgbr)(int, char, lapack_int, lapack_int, lapack_int, lapack_complex_double*, lapack_int, const lapack_complex_double*) = &LAPACKE_zungbr;
template<>
Mat_Ops<lapack_complex_double>* LAPACKE_Matrix_SVD<double, lapack_complex_double >::myOps = new Mat_Ops_Complex<lapack_complex_double>();
#endif

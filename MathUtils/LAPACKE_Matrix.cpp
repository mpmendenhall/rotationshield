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

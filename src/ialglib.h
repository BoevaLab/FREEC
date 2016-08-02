/********************************************************************
optimized ALGLIB subroutines.
********************************************************************/
#ifndef IALGLIB_H
#define IALGLIB_H

#include "ap.h"
#include <string>
#define ALGLIB_INTERCEPTS_ABLAS

namespace ialglib
{
void mv_32(const double *a, const double *x, double *y, int stride, double alpha, double beta);
void mv(int m, int n, const double *a, const double *x, double *y, int stride, double alpha, double beta);
void mv_generic(int m, int n, const double *a, const double *x, double *y, int stride, double alpha, double beta);
void mv_complex(int m, int n, const double *a, const double *x, ap::complex *cy, double *dy, int stride, ap::complex alpha, ap::complex beta);
void mv_complex_generic(int m, int n, const double *a, const double *x, ap::complex *cy, double *dy, int stride, ap::complex alpha, ap::complex beta);

void vzero(int n, double *p, int stride);
void vzero_complex(int n, ap::complex *p, int stride);
void vcopy(int n, const double *a, int stridea, double *b, int strideb);
void vcopy_complex(int n, const ap::complex *a, int stridea, double *b, int strideb, char *conj);
void vcopy_complex(int n, const double *a, int stridea, double *b, int strideb, char *conj);
void mcopyblock(int m, int n, const double *a, int op, int stride, double *b);
void mcopyunblock(int m, int n, const double *a, int op, double *b, int stride);
void mcopyblock_complex(int m, int n, const ap::complex *a, int op, int stride, double *b);
void mcopyunblock_complex(int m, int n, const double *a, int op, ap::complex* b, int stride);

bool _i_rmatrixgemmf(int m,
     int n,
     int k,
     double alpha,
     const ap::real_2d_array& a,
     int ia,
     int ja,
     int optypea,
     const ap::real_2d_array& b,
     int ib,
     int jb,
     int optypeb,
     double beta,
     ap::real_2d_array& c,
     int ic,
     int jc);
bool _i_cmatrixgemmf(int m,
     int n,
     int k,
     ap::complex alpha,
     const ap::complex_2d_array& a,
     int ia,
     int ja,
     int optypea,
     const ap::complex_2d_array& b,
     int ib,
     int jb,
     int optypeb,
     ap::complex beta,
     ap::complex_2d_array& c,
     int ic,
     int jc);
bool _i_cmatrixrighttrsmf(int m,
     int n,
     const ap::complex_2d_array& a,
     int i1,
     int j1,
     bool isupper,
     bool isunit,
     int optype,
     ap::complex_2d_array& x,
     int i2,
     int j2);
bool _i_rmatrixrighttrsmf(int m,
     int n,
     const ap::real_2d_array& a,
     int i1,
     int j1,
     bool isupper,
     bool isunit,
     int optype,
     ap::real_2d_array& x,
     int i2,
     int j2);
bool _i_cmatrixlefttrsmf(int m,
     int n,
     const ap::complex_2d_array& a,
     int i1,
     int j1,
     bool isupper,
     bool isunit,
     int optype,
     ap::complex_2d_array& x,
     int i2,
     int j2);
bool _i_rmatrixlefttrsmf(int m,
     int n,
     const ap::real_2d_array& a,
     int i1,
     int j1,
     bool isupper,
     bool isunit,
     int optype,
     ap::real_2d_array& x,
     int i2,
     int j2);
bool _i_cmatrixsyrkf(int n,
     int k,
     double alpha,
     const ap::complex_2d_array& a,
     int ia,
     int ja,
     int optypea,
     double beta,
     ap::complex_2d_array& c,
     int ic,
     int jc,
     bool isupper);
bool _i_rmatrixsyrkf(int n,
     int k,
     double alpha,
     const ap::real_2d_array& a,
     int ia,
     int ja,
     int optypea,
     double beta,
     ap::real_2d_array& c,
     int ic,
     int jc,
     bool isupper);
bool _i_cmatrixrank1f(int m,
     int n,
     ap::complex_2d_array& a,
     int ia,
     int ja,
     ap::complex_1d_array& u,
     int uoffs,
     ap::complex_1d_array& v,
     int voffs);
bool _i_rmatrixrank1f(int m,
     int n,
     ap::real_2d_array& a,
     int ia,
     int ja,
     ap::real_1d_array& u,
     int uoffs,
     ap::real_1d_array& v,
     int voffs);

}
#endif

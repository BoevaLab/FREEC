/*************************************************************************
Copyright (c) 2009, Sergey Bochkanov (ALGLIB project).

>>> SOURCE LICENSE >>>
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation (www.fsf.org); either version 2 of the 
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License is available at
http://www.fsf.org/licensing/licenses

>>> END OF LICENSE >>>
*************************************************************************/

//#include <stdafx.h>
#include "ablasf.h"

/*************************************************************************
Fast kernel

  -- ALGLIB routine --
     19.01.2010
     Bochkanov Sergey
*************************************************************************/
bool cmatrixrank1f(int m,
     int n,
     ap::complex_2d_array& a,
     int ia,
     int ja,
     ap::complex_1d_array& u,
     int iu,
     ap::complex_1d_array& v,
     int iv)
{
#ifndef ALGLIB_INTERCEPTS_ABLAS
    bool result;

    result = false;
    return result;
#else
    return ialglib::_i_cmatrixrank1f(m, n, a, ia, ja, u, iu, v, iv);
#endif
}


/*************************************************************************
Fast kernel

  -- ALGLIB routine --
     19.01.2010
     Bochkanov Sergey
*************************************************************************/
bool rmatrixrank1f(int m,
     int n,
     ap::real_2d_array& a,
     int ia,
     int ja,
     ap::real_1d_array& u,
     int iu,
     ap::real_1d_array& v,
     int iv)
{
#ifndef ALGLIB_INTERCEPTS_ABLAS
    bool result;

    result = false;
    return result;
#else
    return ialglib::_i_rmatrixrank1f(m, n, a, ia, ja, u, iu, v, iv);
#endif
}


/*************************************************************************
Fast kernel

  -- ALGLIB routine --
     19.01.2010
     Bochkanov Sergey
*************************************************************************/
bool cmatrixmvf(int m,
     int n,
     ap::complex_2d_array& a,
     int ia,
     int ja,
     int opa,
     ap::complex_1d_array& x,
     int ix,
     ap::complex_1d_array& y,
     int iy)
{
    bool result;

    result = false;
    return result;
}


/*************************************************************************
Fast kernel

  -- ALGLIB routine --
     19.01.2010
     Bochkanov Sergey
*************************************************************************/
bool rmatrixmvf(int m,
     int n,
     ap::real_2d_array& a,
     int ia,
     int ja,
     int opa,
     ap::real_1d_array& x,
     int ix,
     ap::real_1d_array& y,
     int iy)
{
    bool result;

    result = false;
    return result;
}


/*************************************************************************
Fast kernel

  -- ALGLIB routine --
     19.01.2010
     Bochkanov Sergey
*************************************************************************/
bool cmatrixrighttrsmf(int m,
     int n,
     const ap::complex_2d_array& a,
     int i1,
     int j1,
     bool isupper,
     bool isunit,
     int optype,
     ap::complex_2d_array& x,
     int i2,
     int j2)
{
#ifndef ALGLIB_INTERCEPTS_ABLAS
    bool result;

    result = false;
    return result;
#else
    return ialglib::_i_cmatrixrighttrsmf(m, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
#endif
}


/*************************************************************************
Fast kernel

  -- ALGLIB routine --
     19.01.2010
     Bochkanov Sergey
*************************************************************************/
bool cmatrixlefttrsmf(int m,
     int n,
     const ap::complex_2d_array& a,
     int i1,
     int j1,
     bool isupper,
     bool isunit,
     int optype,
     ap::complex_2d_array& x,
     int i2,
     int j2)
{
#ifndef ALGLIB_INTERCEPTS_ABLAS
    bool result;

    result = false;
    return result;
#else
    return ialglib::_i_cmatrixlefttrsmf(m, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
#endif
}


/*************************************************************************
Fast kernel

  -- ALGLIB routine --
     19.01.2010
     Bochkanov Sergey
*************************************************************************/
bool rmatrixrighttrsmf(int m,
     int n,
     const ap::real_2d_array& a,
     int i1,
     int j1,
     bool isupper,
     bool isunit,
     int optype,
     ap::real_2d_array& x,
     int i2,
     int j2)
{
#ifndef ALGLIB_INTERCEPTS_ABLAS
    bool result;

    result = false;
    return result;
#else
    return ialglib::_i_rmatrixrighttrsmf(m, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
#endif
}


/*************************************************************************
Fast kernel

  -- ALGLIB routine --
     19.01.2010
     Bochkanov Sergey
*************************************************************************/
bool rmatrixlefttrsmf(int m,
     int n,
     const ap::real_2d_array& a,
     int i1,
     int j1,
     bool isupper,
     bool isunit,
     int optype,
     ap::real_2d_array& x,
     int i2,
     int j2)
{
#ifndef ALGLIB_INTERCEPTS_ABLAS
    bool result;

    result = false;
    return result;
#else
    return ialglib::_i_rmatrixlefttrsmf(m, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
#endif
}


/*************************************************************************
Fast kernel

  -- ALGLIB routine --
     19.01.2010
     Bochkanov Sergey
*************************************************************************/
bool cmatrixsyrkf(int n,
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
     bool isupper)
{
#ifndef ALGLIB_INTERCEPTS_ABLAS
    bool result;

    result = false;
    return result;
#else
    return ialglib::_i_cmatrixsyrkf(n, k, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper);
#endif
}


/*************************************************************************
Fast kernel

  -- ALGLIB routine --
     19.01.2010
     Bochkanov Sergey
*************************************************************************/
bool rmatrixsyrkf(int n,
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
     bool isupper)
{
#ifndef ALGLIB_INTERCEPTS_ABLAS
    bool result;

    result = false;
    return result;
#else
    return ialglib::_i_rmatrixsyrkf(n, k, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper);
#endif
}


/*************************************************************************
Fast kernel

  -- ALGLIB routine --
     19.01.2010
     Bochkanov Sergey
*************************************************************************/
bool rmatrixgemmf(int m,
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
     int jc)
{
#ifndef ALGLIB_INTERCEPTS_ABLAS
    bool result;

    result = false;
    return result;
#else
    return ialglib::_i_rmatrixgemmf(m, n, k, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc);
#endif
}


/*************************************************************************
Fast kernel

  -- ALGLIB routine --
     19.01.2010
     Bochkanov Sergey
*************************************************************************/
bool cmatrixgemmf(int m,
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
     int jc)
{
#ifndef ALGLIB_INTERCEPTS_ABLAS
    bool result;

    result = false;
    return result;
#else
    return ialglib::_i_cmatrixgemmf(m, n, k, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc);
#endif
}





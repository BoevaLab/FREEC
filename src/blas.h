/*************************************************************************
Copyright (c) 2005-2007, Sergey Bochkanov (ALGLIB project).

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

#ifndef _blas_h
#define _blas_h

#include "ap.h"
#include "ialglib.h"

double vectornorm2(const ap::real_1d_array& x, int i1, int i2);


int vectoridxabsmax(const ap::real_1d_array& x, int i1, int i2);


int columnidxabsmax(const ap::real_2d_array& x, int i1, int i2, int j);


int rowidxabsmax(const ap::real_2d_array& x, int j1, int j2, int i);


double upperhessenberg1norm(const ap::real_2d_array& a,
     int i1,
     int i2,
     int j1,
     int j2,
     ap::real_1d_array& work);


void copymatrix(const ap::real_2d_array& a,
     int is1,
     int is2,
     int js1,
     int js2,
     ap::real_2d_array& b,
     int id1,
     int id2,
     int jd1,
     int jd2);


void inplacetranspose(ap::real_2d_array& a,
     int i1,
     int i2,
     int j1,
     int j2,
     ap::real_1d_array& work);


void copyandtranspose(const ap::real_2d_array& a,
     int is1,
     int is2,
     int js1,
     int js2,
     ap::real_2d_array& b,
     int id1,
     int id2,
     int jd1,
     int jd2);


void matrixvectormultiply(const ap::real_2d_array& a,
     int i1,
     int i2,
     int j1,
     int j2,
     bool trans,
     const ap::real_1d_array& x,
     int ix1,
     int ix2,
     double alpha,
     ap::real_1d_array& y,
     int iy1,
     int iy2,
     double beta);


double pythag2(double x, double y);


void matrixmatrixmultiply(const ap::real_2d_array& a,
     int ai1,
     int ai2,
     int aj1,
     int aj2,
     bool transa,
     const ap::real_2d_array& b,
     int bi1,
     int bi2,
     int bj1,
     int bj2,
     bool transb,
     double alpha,
     ap::real_2d_array& c,
     int ci1,
     int ci2,
     int cj1,
     int cj2,
     double beta,
     ap::real_1d_array& work);


#endif


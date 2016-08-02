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

#ifndef _svd_h
#define _svd_h

#include "ap.h"
#include "ialglib.h"

#include "hblas.h"
#include "reflections.h"
#include "creflections.h"
#include "sblas.h"
#include "ablasf.h"
#include "ablas.h"
#include "ortfac.h"
#include "blas.h"
#include "rotations.h"
#include "bdsvd.h"


/*************************************************************************
Singular value decomposition of a rectangular matrix.

The algorithm calculates the singular value decomposition of a matrix of
size MxN: A = U * S * V^T

The algorithm finds the singular values and, optionally, matrices U and V^T.
The algorithm can find both first min(M,N) columns of matrix U and rows of
matrix V^T (singular vectors), and matrices U and V^T wholly (of sizes MxM
and NxN respectively).

Take into account that the subroutine does not return matrix V but V^T.

Input parameters:
    A           -   matrix to be decomposed.
                    Array whose indexes range within [0..M-1, 0..N-1].
    M           -   number of rows in matrix A.
    N           -   number of columns in matrix A.
    UNeeded     -   0, 1 or 2. See the description of the parameter U.
    VTNeeded    -   0, 1 or 2. See the description of the parameter VT.
    AdditionalMemory -
                    If the parameter:
                     * equals 0, the algorithm doesn’t use additional
                       memory (lower requirements, lower performance).
                     * equals 1, the algorithm uses additional
                       memory of size min(M,N)*min(M,N) of real numbers.
                       It often speeds up the algorithm.
                     * equals 2, the algorithm uses additional
                       memory of size M*min(M,N) of real numbers.
                       It allows to get a maximum performance.
                    The recommended value of the parameter is 2.

Output parameters:
    W           -   contains singular values in descending order.
    U           -   if UNeeded=0, U isn't changed, the left singular vectors
                    are not calculated.
                    if Uneeded=1, U contains left singular vectors (first
                    min(M,N) columns of matrix U). Array whose indexes range
                    within [0..M-1, 0..Min(M,N)-1].
                    if UNeeded=2, U contains matrix U wholly. Array whose
                    indexes range within [0..M-1, 0..M-1].
    VT          -   if VTNeeded=0, VT isn’t changed, the right singular vectors
                    are not calculated.
                    if VTNeeded=1, VT contains right singular vectors (first
                    min(M,N) rows of matrix V^T). Array whose indexes range
                    within [0..min(M,N)-1, 0..N-1].
                    if VTNeeded=2, VT contains matrix V^T wholly. Array whose
                    indexes range within [0..N-1, 0..N-1].

  -- ALGLIB --
     Copyright 2005 by Bochkanov Sergey
*************************************************************************/
bool rmatrixsvd(ap::real_2d_array a,
     int m,
     int n,
     int uneeded,
     int vtneeded,
     int additionalmemory,
     ap::real_1d_array& w,
     ap::real_2d_array& u,
     ap::real_2d_array& vt);


#endif


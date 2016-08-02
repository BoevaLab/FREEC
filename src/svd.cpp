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

//#include <stdafx.h>
#include "svd.h"

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
     ap::real_2d_array& vt)
{
    bool result;
    ap::real_1d_array tauq;
    ap::real_1d_array taup;
    ap::real_1d_array tau;
    ap::real_1d_array e;
    ap::real_1d_array work;
    ap::real_2d_array t2;
    bool isupper;
    int minmn;
    int ncu;
    int nrvt;
    int nru;
    int ncvt;
    int i;
    int j;

    result = true;
    if( m==0||n==0 )
    {
        return result;
    }
    ap::ap_error::make_assertion(uneeded>=0&&uneeded<=2, "SVDDecomposition: wrong parameters!");
    ap::ap_error::make_assertion(vtneeded>=0&&vtneeded<=2, "SVDDecomposition: wrong parameters!");
    ap::ap_error::make_assertion(additionalmemory>=0&&additionalmemory<=2, "SVDDecomposition: wrong parameters!");
    
    //
    // initialize
    //
    minmn = ap::minint(m, n);
    w.setbounds(1, minmn);
    ncu = 0;
    nru = 0;
    if( uneeded==1 )
    {
        nru = m;
        ncu = minmn;
        u.setbounds(0, nru-1, 0, ncu-1);
    }
    if( uneeded==2 )
    {
        nru = m;
        ncu = m;
        u.setbounds(0, nru-1, 0, ncu-1);
    }
    nrvt = 0;
    ncvt = 0;
    if( vtneeded==1 )
    {
        nrvt = minmn;
        ncvt = n;
        vt.setbounds(0, nrvt-1, 0, ncvt-1);
    }
    if( vtneeded==2 )
    {
        nrvt = n;
        ncvt = n;
        vt.setbounds(0, nrvt-1, 0, ncvt-1);
    }
    
    //
    // M much larger than N
    // Use bidiagonal reduction with QR-decomposition
    //
    if( ap::fp_greater(m,1.6*n) )
    {
        if( uneeded==0 )
        {
            
            //
            // No left singular vectors to be computed
            //
            rmatrixqr(a, m, n, tau);
            for(i = 0; i <= n-1; i++)
            {
                for(j = 0; j <= i-1; j++)
                {
                    a(i,j) = 0;
                }
            }
            rmatrixbd(a, n, n, tauq, taup);
            rmatrixbdunpackpt(a, n, n, taup, nrvt, vt);
            rmatrixbdunpackdiagonals(a, n, n, isupper, w, e);
            result = rmatrixbdsvd(w, e, n, isupper, false, u, 0, a, 0, vt, ncvt);
            return result;
        }
        else
        {
            
            //
            // Left singular vectors (may be full matrix U) to be computed
            //
            rmatrixqr(a, m, n, tau);
            rmatrixqrunpackq(a, m, n, tau, ncu, u);
            for(i = 0; i <= n-1; i++)
            {
                for(j = 0; j <= i-1; j++)
                {
                    a(i,j) = 0;
                }
            }
            rmatrixbd(a, n, n, tauq, taup);
            rmatrixbdunpackpt(a, n, n, taup, nrvt, vt);
            rmatrixbdunpackdiagonals(a, n, n, isupper, w, e);
            if( additionalmemory<1 )
            {
                
                //
                // No additional memory can be used
                //
                rmatrixbdmultiplybyq(a, n, n, tauq, u, m, n, true, false);
                result = rmatrixbdsvd(w, e, n, isupper, false, u, m, a, 0, vt, ncvt);
            }
            else
            {
                
                //
                // Large U. Transforming intermediate matrix T2
                //
                work.setbounds(1, ap::maxint(m, n));
                rmatrixbdunpackq(a, n, n, tauq, n, t2);
                copymatrix(u, 0, m-1, 0, n-1, a, 0, m-1, 0, n-1);
                inplacetranspose(t2, 0, n-1, 0, n-1, work);
                result = rmatrixbdsvd(w, e, n, isupper, false, u, 0, t2, n, vt, ncvt);
                matrixmatrixmultiply(a, 0, m-1, 0, n-1, false, t2, 0, n-1, 0, n-1, true, 1.0, u, 0, m-1, 0, n-1, 0.0, work);
            }
            return result;
        }
    }
    
    //
    // N much larger than M
    // Use bidiagonal reduction with LQ-decomposition
    //
    if( ap::fp_greater(n,1.6*m) )
    {
        if( vtneeded==0 )
        {
            
            //
            // No right singular vectors to be computed
            //
            rmatrixlq(a, m, n, tau);
            for(i = 0; i <= m-1; i++)
            {
                for(j = i+1; j <= m-1; j++)
                {
                    a(i,j) = 0;
                }
            }
            rmatrixbd(a, m, m, tauq, taup);
            rmatrixbdunpackq(a, m, m, tauq, ncu, u);
            rmatrixbdunpackdiagonals(a, m, m, isupper, w, e);
            work.setbounds(1, m);
            inplacetranspose(u, 0, nru-1, 0, ncu-1, work);
            result = rmatrixbdsvd(w, e, m, isupper, false, a, 0, u, nru, vt, 0);
            inplacetranspose(u, 0, nru-1, 0, ncu-1, work);
            return result;
        }
        else
        {
            
            //
            // Right singular vectors (may be full matrix VT) to be computed
            //
            rmatrixlq(a, m, n, tau);
            rmatrixlqunpackq(a, m, n, tau, nrvt, vt);
            for(i = 0; i <= m-1; i++)
            {
                for(j = i+1; j <= m-1; j++)
                {
                    a(i,j) = 0;
                }
            }
            rmatrixbd(a, m, m, tauq, taup);
            rmatrixbdunpackq(a, m, m, tauq, ncu, u);
            rmatrixbdunpackdiagonals(a, m, m, isupper, w, e);
            work.setbounds(1, ap::maxint(m, n));
            inplacetranspose(u, 0, nru-1, 0, ncu-1, work);
            if( additionalmemory<1 )
            {
                
                //
                // No additional memory available
                //
                rmatrixbdmultiplybyp(a, m, m, taup, vt, m, n, false, true);
                result = rmatrixbdsvd(w, e, m, isupper, false, a, 0, u, nru, vt, n);
            }
            else
            {
                
                //
                // Large VT. Transforming intermediate matrix T2
                //
                rmatrixbdunpackpt(a, m, m, taup, m, t2);
                result = rmatrixbdsvd(w, e, m, isupper, false, a, 0, u, nru, t2, m);
                copymatrix(vt, 0, m-1, 0, n-1, a, 0, m-1, 0, n-1);
                matrixmatrixmultiply(t2, 0, m-1, 0, m-1, false, a, 0, m-1, 0, n-1, false, 1.0, vt, 0, m-1, 0, n-1, 0.0, work);
            }
            inplacetranspose(u, 0, nru-1, 0, ncu-1, work);
            return result;
        }
    }
    
    //
    // M<=N
    // We can use inplace transposition of U to get rid of columnwise operations
    //
    if( m<=n )
    {
        rmatrixbd(a, m, n, tauq, taup);
        rmatrixbdunpackq(a, m, n, tauq, ncu, u);
        rmatrixbdunpackpt(a, m, n, taup, nrvt, vt);
        rmatrixbdunpackdiagonals(a, m, n, isupper, w, e);
        work.setbounds(1, m);
        inplacetranspose(u, 0, nru-1, 0, ncu-1, work);
        result = rmatrixbdsvd(w, e, minmn, isupper, false, a, 0, u, nru, vt, ncvt);
        inplacetranspose(u, 0, nru-1, 0, ncu-1, work);
        return result;
    }
    
    //
    // Simple bidiagonal reduction
    //
    rmatrixbd(a, m, n, tauq, taup);
    rmatrixbdunpackq(a, m, n, tauq, ncu, u);
    rmatrixbdunpackpt(a, m, n, taup, nrvt, vt);
    rmatrixbdunpackdiagonals(a, m, n, isupper, w, e);
    if( additionalmemory<2||uneeded==0 )
    {
        
        //
        // We cant use additional memory or there is no need in such operations
        //
        result = rmatrixbdsvd(w, e, minmn, isupper, false, u, nru, a, 0, vt, ncvt);
    }
    else
    {
        
        //
        // We can use additional memory
        //
        t2.setbounds(0, minmn-1, 0, m-1);
        copyandtranspose(u, 0, m-1, 0, minmn-1, t2, 0, minmn-1, 0, m-1);
        result = rmatrixbdsvd(w, e, minmn, isupper, false, u, 0, t2, m, vt, ncvt);
        copyandtranspose(t2, 0, minmn-1, 0, m-1, u, 0, m-1, 0, minmn-1);
    }
    return result;
}





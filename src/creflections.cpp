/*************************************************************************
Copyright (c) 1992-2007 The University of Tennessee.  All rights reserved.

Contributors:
    * Sergey Bochkanov (ALGLIB project). Translation from FORTRAN to
      pseudocode.

See subroutines comments for additional copyrights.

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
#include "creflections.h"

/*************************************************************************
Generation of an elementary complex reflection transformation

The subroutine generates elementary complex reflection H of  order  N,  so
that, for a given X, the following equality holds true:

     ( X(1) )   ( Beta )
H' * (  ..  ) = (  0   ),   H'*H = I,   Beta is a real number
     ( X(n) )   (  0   )

where

              ( V(1) )
H = 1 - Tau * (  ..  ) * ( conj(V(1)), ..., conj(V(n)) )
              ( V(n) )

where the first component of vector V equals 1.

Input parameters:
    X   -   vector. Array with elements [1..N].
    N   -   reflection order.

Output parameters:
    X   -   components from 2 to N are replaced by vector V.
            The first component is replaced with parameter Beta.
    Tau -   scalar value Tau.

This subroutine is the modification of CLARFG subroutines  from the LAPACK
library. It has similar functionality except for the fact that it  doesn’t
handle errors when intermediate results cause an overflow.

  -- LAPACK auxiliary routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     September 30, 1994
*************************************************************************/
void complexgeneratereflection(ap::complex_1d_array& x,
     int n,
     ap::complex& tau)
{
    int j;
    ap::complex alpha;
    double alphi;
    double alphr;
    double beta;
    double xnorm;
    double mx;
    ap::complex t;
    double s;
    ap::complex v;

    if( n<=0 )
    {
        tau = 0;
        return;
    }
    
    //
    // Scale if needed (to avoid overflow/underflow during intermediate
    // calculations).
    //
    mx = 0;
    for(j = 1; j <= n; j++)
    {
        mx = ap::maxreal(ap::abscomplex(x(j)), mx);
    }
    s = 1;
    if( ap::fp_neq(mx,0) )
    {
        if( ap::fp_less(mx,1) )
        {
            s = sqrt(ap::minrealnumber);
            v = 1/s;
            ap::vmul(&x(1), 1, ap::vlen(1,n), v);
        }
        else
        {
            s = sqrt(ap::maxrealnumber);
            v = 1/s;
            ap::vmul(&x(1), 1, ap::vlen(1,n), v);
        }
    }
    
    //
    // calculate
    //
    alpha = x(1);
    mx = 0;
    for(j = 2; j <= n; j++)
    {
        mx = ap::maxreal(ap::abscomplex(x(j)), mx);
    }
    xnorm = 0;
    if( ap::fp_neq(mx,0) )
    {
        for(j = 2; j <= n; j++)
        {
            t = x(j)/mx;
            xnorm = xnorm+(t*ap::conj(t)).x;
        }
        xnorm = sqrt(xnorm)*mx;
    }
    alphr = alpha.x;
    alphi = alpha.y;
    if( ap::fp_eq(xnorm,0)&&ap::fp_eq(alphi,0) )
    {
        tau = 0;
        x(1) = x(1)*s;
        return;
    }
    mx = ap::maxreal(fabs(alphr), fabs(alphi));
    mx = ap::maxreal(mx, fabs(xnorm));
    beta = -mx*sqrt(ap::sqr(alphr/mx)+ap::sqr(alphi/mx)+ap::sqr(xnorm/mx));
    if( ap::fp_less(alphr,0) )
    {
        beta = -beta;
    }
    tau.x = (beta-alphr)/beta;
    tau.y = -alphi/beta;
    alpha = 1/(alpha-beta);
    if( n>1 )
    {
        ap::vmul(&x(2), 1, ap::vlen(2,n), alpha);
    }
    alpha = beta;
    x(1) = alpha;
    
    //
    // Scale back
    //
    x(1) = x(1)*s;
}


/*************************************************************************
Application of an elementary reflection to a rectangular matrix of size MxN

The  algorithm  pre-multiplies  the  matrix  by  an  elementary reflection
transformation  which  is  given  by  column  V  and  scalar  Tau (see the
description of the GenerateReflection). Not the whole matrix  but  only  a
part of it is transformed (rows from M1 to M2, columns from N1 to N2). Only
the elements of this submatrix are changed.

Note: the matrix is multiplied by H, not by H'.   If  it  is  required  to
multiply the matrix by H', it is necessary to pass Conj(Tau) instead of Tau.

Input parameters:
    C       -   matrix to be transformed.
    Tau     -   scalar defining transformation.
    V       -   column defining transformation.
                Array whose index ranges within [1..M2-M1+1]
    M1, M2  -   range of rows to be transformed.
    N1, N2  -   range of columns to be transformed.
    WORK    -   working array whose index goes from N1 to N2.

Output parameters:
    C       -   the result of multiplying the input matrix C by the
                transformation matrix which is given by Tau and V.
                If N1>N2 or M1>M2, C is not modified.

  -- LAPACK auxiliary routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     September 30, 1994
*************************************************************************/
void complexapplyreflectionfromtheleft(ap::complex_2d_array& c,
     ap::complex tau,
     const ap::complex_1d_array& v,
     int m1,
     int m2,
     int n1,
     int n2,
     ap::complex_1d_array& work)
{
    ap::complex t;
    int i;
    int vm;

    if( tau==0||n1>n2||m1>m2 )
    {
        return;
    }
    
    //
    // w := C^T * conj(v)
    //
    vm = m2-m1+1;
    for(i = n1; i <= n2; i++)
    {
        work(i) = 0;
    }
    for(i = m1; i <= m2; i++)
    {
        t = ap::conj(v(i+1-m1));
        ap::vadd(&work(n1), 1, &c(i, n1), 1, "N", ap::vlen(n1,n2), t);
    }
    
    //
    // C := C - tau * v * w^T
    //
    for(i = m1; i <= m2; i++)
    {
        t = v(i-m1+1)*tau;
        ap::vsub(&c(i, n1), 1, &work(n1), 1, "N", ap::vlen(n1,n2), t);
    }
}


/*************************************************************************
Application of an elementary reflection to a rectangular matrix of size MxN

The  algorithm  post-multiplies  the  matrix  by  an elementary reflection
transformation  which  is  given  by  column  V  and  scalar  Tau (see the
description  of  the  GenerateReflection). Not the whole matrix but only a
part  of  it  is  transformed (rows from M1 to M2, columns from N1 to N2).
Only the elements of this submatrix are changed.

Input parameters:
    C       -   matrix to be transformed.
    Tau     -   scalar defining transformation.
    V       -   column defining transformation.
                Array whose index ranges within [1..N2-N1+1]
    M1, M2  -   range of rows to be transformed.
    N1, N2  -   range of columns to be transformed.
    WORK    -   working array whose index goes from M1 to M2.

Output parameters:
    C       -   the result of multiplying the input matrix C by the
                transformation matrix which is given by Tau and V.
                If N1>N2 or M1>M2, C is not modified.

  -- LAPACK auxiliary routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     September 30, 1994
*************************************************************************/
void complexapplyreflectionfromtheright(ap::complex_2d_array& c,
     ap::complex tau,
     ap::complex_1d_array& v,
     int m1,
     int m2,
     int n1,
     int n2,
     ap::complex_1d_array& work)
{
    ap::complex t;
    int i;
    int vm;

    if( tau==0||n1>n2||m1>m2 )
    {
        return;
    }
    
    //
    // w := C * v
    //
    vm = n2-n1+1;
    for(i = m1; i <= m2; i++)
    {
        t = ap::vdotproduct(&c(i, n1), 1, "N", &v(1), 1, "N", ap::vlen(n1,n2));
        work(i) = t;
    }
    
    //
    // C := C - w * conj(v^T)
    //
    ap::vmove(&v(1), 1, &v(1), 1, "Conj", ap::vlen(1,vm));
    for(i = m1; i <= m2; i++)
    {
        t = work(i)*tau;
        ap::vsub(&c(i, n1), 1, &v(1), 1, "N", ap::vlen(n1,n2), t);
    }
    ap::vmove(&v(1), 1, &v(1), 1, "Conj", ap::vlen(1,vm));
}





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
#include "sblas.h"

void symmetricmatrixvectormultiply(const ap::real_2d_array& a,
     bool isupper,
     int i1,
     int i2,
     const ap::real_1d_array& x,
     double alpha,
     ap::real_1d_array& y)
{
    int i;
    int ba1;
    int ba2;
    int by1;
    int by2;
    int bx1;
    int bx2;
    int n;
    double v;

    n = i2-i1+1;
    if( n<=0 )
    {
        return;
    }
    
    //
    // Let A = L + D + U, where
    //  L is strictly lower triangular (main diagonal is zero)
    //  D is diagonal
    //  U is strictly upper triangular (main diagonal is zero)
    //
    // A*x = L*x + D*x + U*x
    //
    // Calculate D*x first
    //
    for(i = i1; i <= i2; i++)
    {
        y(i-i1+1) = a(i,i)*x(i-i1+1);
    }
    
    //
    // Add L*x + U*x
    //
    if( isupper )
    {
        for(i = i1; i <= i2-1; i++)
        {
            
            //
            // Add L*x to the result
            //
            v = x(i-i1+1);
            by1 = i-i1+2;
            by2 = n;
            ba1 = i+1;
            ba2 = i2;
            ap::vadd(&y(by1), 1, &a(i, ba1), 1, ap::vlen(by1,by2), v);
            
            //
            // Add U*x to the result
            //
            bx1 = i-i1+2;
            bx2 = n;
            ba1 = i+1;
            ba2 = i2;
            v = ap::vdotproduct(&x(bx1), 1, &a(i, ba1), 1, ap::vlen(bx1,bx2));
            y(i-i1+1) = y(i-i1+1)+v;
        }
    }
    else
    {
        for(i = i1+1; i <= i2; i++)
        {
            
            //
            // Add L*x to the result
            //
            bx1 = 1;
            bx2 = i-i1;
            ba1 = i1;
            ba2 = i-1;
            v = ap::vdotproduct(&x(bx1), 1, &a(i, ba1), 1, ap::vlen(bx1,bx2));
            y(i-i1+1) = y(i-i1+1)+v;
            
            //
            // Add U*x to the result
            //
            v = x(i-i1+1);
            by1 = 1;
            by2 = i-i1;
            ba1 = i1;
            ba2 = i-1;
            ap::vadd(&y(by1), 1, &a(i, ba1), 1, ap::vlen(by1,by2), v);
        }
    }
    ap::vmul(&y(1), 1, ap::vlen(1,n), alpha);
}


void symmetricrank2update(ap::real_2d_array& a,
     bool isupper,
     int i1,
     int i2,
     const ap::real_1d_array& x,
     const ap::real_1d_array& y,
     ap::real_1d_array& t,
     double alpha)
{
    int i;
    int tp1;
    int tp2;
    double v;

    if( isupper )
    {
        for(i = i1; i <= i2; i++)
        {
            tp1 = i+1-i1;
            tp2 = i2-i1+1;
            v = x(i+1-i1);
            ap::vmove(&t(tp1), 1, &y(tp1), 1, ap::vlen(tp1,tp2), v);
            v = y(i+1-i1);
            ap::vadd(&t(tp1), 1, &x(tp1), 1, ap::vlen(tp1,tp2), v);
            ap::vmul(&t(tp1), 1, ap::vlen(tp1,tp2), alpha);
            ap::vadd(&a(i, i), 1, &t(tp1), 1, ap::vlen(i,i2));
        }
    }
    else
    {
        for(i = i1; i <= i2; i++)
        {
            tp1 = 1;
            tp2 = i+1-i1;
            v = x(i+1-i1);
            ap::vmove(&t(tp1), 1, &y(tp1), 1, ap::vlen(tp1,tp2), v);
            v = y(i+1-i1);
            ap::vadd(&t(tp1), 1, &x(tp1), 1, ap::vlen(tp1,tp2), v);
            ap::vmul(&t(tp1), 1, ap::vlen(tp1,tp2), alpha);
            ap::vadd(&a(i, i1), 1, &t(tp1), 1, ap::vlen(i1,i));
        }
    }
}





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
#include "rotations.h"

/*************************************************************************
Application of a sequence of  elementary rotations to a matrix

The algorithm pre-multiplies the matrix by a sequence of rotation
transformations which is given by arrays C and S. Depending on the value
of the IsForward parameter either 1 and 2, 3 and 4 and so on (if IsForward=true)
rows are rotated, or the rows N and N-1, N-2 and N-3 and so on, are rotated.

Not the whole matrix but only a part of it is transformed (rows from M1 to
M2, columns from N1 to N2). Only the elements of this submatrix are changed.

Input parameters:
    IsForward   -   the sequence of the rotation application.
    M1,M2       -   the range of rows to be transformed.
    N1, N2      -   the range of columns to be transformed.
    C,S         -   transformation coefficients.
                    Array whose index ranges within [1..M2-M1].
    A           -   processed matrix.
    WORK        -   working array whose index ranges within [N1..N2].

Output parameters:
    A           -   transformed matrix.

Utility subroutine.
*************************************************************************/
void applyrotationsfromtheleft(bool isforward,
     int m1,
     int m2,
     int n1,
     int n2,
     const ap::real_1d_array& c,
     const ap::real_1d_array& s,
     ap::real_2d_array& a,
     ap::real_1d_array& work)
{
    int j;
    int jp1;
    double ctemp;
    double stemp;
    double temp;

    if( m1>m2||n1>n2 )
    {
        return;
    }
    
    //
    // Form  P * A
    //
    if( isforward )
    {
        if( n1!=n2 )
        {
            
            //
            // Common case: N1<>N2
            //
            for(j = m1; j <= m2-1; j++)
            {
                ctemp = c(j-m1+1);
                stemp = s(j-m1+1);
                if( ap::fp_neq(ctemp,1)||ap::fp_neq(stemp,0) )
                {
                    jp1 = j+1;
                    ap::vmove(&work(n1), 1, &a(jp1, n1), 1, ap::vlen(n1,n2), ctemp);
                    ap::vsub(&work(n1), 1, &a(j, n1), 1, ap::vlen(n1,n2), stemp);
                    ap::vmul(&a(j, n1), 1, ap::vlen(n1,n2), ctemp);
                    ap::vadd(&a(j, n1), 1, &a(jp1, n1), 1, ap::vlen(n1,n2), stemp);
                    ap::vmove(&a(jp1, n1), 1, &work(n1), 1, ap::vlen(n1,n2));
                }
            }
        }
        else
        {
            
            //
            // Special case: N1=N2
            //
            for(j = m1; j <= m2-1; j++)
            {
                ctemp = c(j-m1+1);
                stemp = s(j-m1+1);
                if( ap::fp_neq(ctemp,1)||ap::fp_neq(stemp,0) )
                {
                    temp = a(j+1,n1);
                    a(j+1,n1) = ctemp*temp-stemp*a(j,n1);
                    a(j,n1) = stemp*temp+ctemp*a(j,n1);
                }
            }
        }
    }
    else
    {
        if( n1!=n2 )
        {
            
            //
            // Common case: N1<>N2
            //
            for(j = m2-1; j >= m1; j--)
            {
                ctemp = c(j-m1+1);
                stemp = s(j-m1+1);
                if( ap::fp_neq(ctemp,1)||ap::fp_neq(stemp,0) )
                {
                    jp1 = j+1;
                    ap::vmove(&work(n1), 1, &a(jp1, n1), 1, ap::vlen(n1,n2), ctemp);
                    ap::vsub(&work(n1), 1, &a(j, n1), 1, ap::vlen(n1,n2), stemp);
                    ap::vmul(&a(j, n1), 1, ap::vlen(n1,n2), ctemp);
                    ap::vadd(&a(j, n1), 1, &a(jp1, n1), 1, ap::vlen(n1,n2), stemp);
                    ap::vmove(&a(jp1, n1), 1, &work(n1), 1, ap::vlen(n1,n2));
                }
            }
        }
        else
        {
            
            //
            // Special case: N1=N2
            //
            for(j = m2-1; j >= m1; j--)
            {
                ctemp = c(j-m1+1);
                stemp = s(j-m1+1);
                if( ap::fp_neq(ctemp,1)||ap::fp_neq(stemp,0) )
                {
                    temp = a(j+1,n1);
                    a(j+1,n1) = ctemp*temp-stemp*a(j,n1);
                    a(j,n1) = stemp*temp+ctemp*a(j,n1);
                }
            }
        }
    }
}


/*************************************************************************
Application of a sequence of  elementary rotations to a matrix

The algorithm post-multiplies the matrix by a sequence of rotation
transformations which is given by arrays C and S. Depending on the value
of the IsForward parameter either 1 and 2, 3 and 4 and so on (if IsForward=true)
rows are rotated, or the rows N and N-1, N-2 and N-3 and so on are rotated.

Not the whole matrix but only a part of it is transformed (rows from M1
to M2, columns from N1 to N2). Only the elements of this submatrix are changed.

Input parameters:
    IsForward   -   the sequence of the rotation application.
    M1,M2       -   the range of rows to be transformed.
    N1, N2      -   the range of columns to be transformed.
    C,S         -   transformation coefficients.
                    Array whose index ranges within [1..N2-N1].
    A           -   processed matrix.
    WORK        -   working array whose index ranges within [M1..M2].

Output parameters:
    A           -   transformed matrix.

Utility subroutine.
*************************************************************************/
void applyrotationsfromtheright(bool isforward,
     int m1,
     int m2,
     int n1,
     int n2,
     const ap::real_1d_array& c,
     const ap::real_1d_array& s,
     ap::real_2d_array& a,
     ap::real_1d_array& work)
{
    int j;
    int jp1;
    double ctemp;
    double stemp;
    double temp;

    
    //
    // Form A * P'
    //
    if( isforward )
    {
        if( m1!=m2 )
        {
            
            //
            // Common case: M1<>M2
            //
            for(j = n1; j <= n2-1; j++)
            {
                ctemp = c(j-n1+1);
                stemp = s(j-n1+1);
                if( ap::fp_neq(ctemp,1)||ap::fp_neq(stemp,0) )
                {
                    jp1 = j+1;
                    ap::vmove(&work(m1), 1, &a(m1, jp1), a.getstride(), ap::vlen(m1,m2), ctemp);
                    ap::vsub(&work(m1), 1, &a(m1, j), a.getstride(), ap::vlen(m1,m2), stemp);
                    ap::vmul(&a(m1, j), a.getstride(), ap::vlen(m1,m2), ctemp);
                    ap::vadd(&a(m1, j), a.getstride(), &a(m1, jp1), a.getstride(), ap::vlen(m1,m2), stemp);
                    ap::vmove(&a(m1, jp1), a.getstride(), &work(m1), 1, ap::vlen(m1,m2));
                }
            }
        }
        else
        {
            
            //
            // Special case: M1=M2
            //
            for(j = n1; j <= n2-1; j++)
            {
                ctemp = c(j-n1+1);
                stemp = s(j-n1+1);
                if( ap::fp_neq(ctemp,1)||ap::fp_neq(stemp,0) )
                {
                    temp = a(m1,j+1);
                    a(m1,j+1) = ctemp*temp-stemp*a(m1,j);
                    a(m1,j) = stemp*temp+ctemp*a(m1,j);
                }
            }
        }
    }
    else
    {
        if( m1!=m2 )
        {
            
            //
            // Common case: M1<>M2
            //
            for(j = n2-1; j >= n1; j--)
            {
                ctemp = c(j-n1+1);
                stemp = s(j-n1+1);
                if( ap::fp_neq(ctemp,1)||ap::fp_neq(stemp,0) )
                {
                    jp1 = j+1;
                    ap::vmove(&work(m1), 1, &a(m1, jp1), a.getstride(), ap::vlen(m1,m2), ctemp);
                    ap::vsub(&work(m1), 1, &a(m1, j), a.getstride(), ap::vlen(m1,m2), stemp);
                    ap::vmul(&a(m1, j), a.getstride(), ap::vlen(m1,m2), ctemp);
                    ap::vadd(&a(m1, j), a.getstride(), &a(m1, jp1), a.getstride(), ap::vlen(m1,m2), stemp);
                    ap::vmove(&a(m1, jp1), a.getstride(), &work(m1), 1, ap::vlen(m1,m2));
                }
            }
        }
        else
        {
            
            //
            // Special case: M1=M2
            //
            for(j = n2-1; j >= n1; j--)
            {
                ctemp = c(j-n1+1);
                stemp = s(j-n1+1);
                if( ap::fp_neq(ctemp,1)||ap::fp_neq(stemp,0) )
                {
                    temp = a(m1,j+1);
                    a(m1,j+1) = ctemp*temp-stemp*a(m1,j);
                    a(m1,j) = stemp*temp+ctemp*a(m1,j);
                }
            }
        }
    }
}


/*************************************************************************
The subroutine generates the elementary rotation, so that:

[  CS  SN  ]  .  [ F ]  =  [ R ]
[ -SN  CS  ]     [ G ]     [ 0 ]

CS**2 + SN**2 = 1
*************************************************************************/
void generaterotation(double f, double g, double& cs, double& sn, double& r)
{
    double f1;
    double g1;

    if( ap::fp_eq(g,0) )
    {
        cs = 1;
        sn = 0;
        r = f;
    }
    else
    {
        if( ap::fp_eq(f,0) )
        {
            cs = 0;
            sn = 1;
            r = g;
        }
        else
        {
            f1 = f;
            g1 = g;
            if( ap::fp_greater(fabs(f1),fabs(g1)) )
            {
                r = fabs(f1)*sqrt(1+ap::sqr(g1/f1));
            }
            else
            {
                r = fabs(g1)*sqrt(1+ap::sqr(f1/g1));
            }
            cs = f1/r;
            sn = g1/r;
            if( ap::fp_greater(fabs(f),fabs(g))&&ap::fp_less(cs,0) )
            {
                cs = -cs;
                sn = -sn;
                r = -r;
            }
        }
    }
}





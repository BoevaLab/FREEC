/*************************************************************************
Copyright (c) 2009-2010, Sergey Bochkanov (ALGLIB project).

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
#include "ablas.h"

static void ablasinternalsplitlength(int n, int nb, int& n1, int& n2);
static void cmatrixrighttrsm2(int m,
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
static void cmatrixlefttrsm2(int m,
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
static void rmatrixrighttrsm2(int m,
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
static void rmatrixlefttrsm2(int m,
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
static void cmatrixsyrk2(int n,
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
static void rmatrixsyrk2(int n,
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
static void cmatrixgemmk(int m,
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
static void rmatrixgemmk(int m,
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

/*************************************************************************
Splits matrix length in two parts, left part should match ABLAS block size

INPUT PARAMETERS
    A   -   real matrix, is passed to ensure that we didn't split
            complex matrix using real splitting subroutine.
            matrix itself is not changed.
    N   -   length, N>0

OUTPUT PARAMETERS
    N1  -   length
    N2  -   length

N1+N2=N, N1>=N2, N2 may be zero

  -- ALGLIB routine --
     15.12.2009
     Bochkanov Sergey
*************************************************************************/
void ablassplitlength(const ap::real_2d_array& a, int n, int& n1, int& n2)
{

    if( n>ablasblocksize(a) )
    {
        ablasinternalsplitlength(n, ablasblocksize(a), n1, n2);
    }
    else
    {
        ablasinternalsplitlength(n, ablasmicroblocksize(), n1, n2);
    }
}


/*************************************************************************
Complex ABLASSplitLength

  -- ALGLIB routine --
     15.12.2009
     Bochkanov Sergey
*************************************************************************/
void ablascomplexsplitlength(const ap::complex_2d_array& a,
     int n,
     int& n1,
     int& n2)
{

    if( n>ablascomplexblocksize(a) )
    {
        ablasinternalsplitlength(n, ablascomplexblocksize(a), n1, n2);
    }
    else
    {
        ablasinternalsplitlength(n, ablasmicroblocksize(), n1, n2);
    }
}


/*************************************************************************
Returns block size - subdivision size where  cache-oblivious  soubroutines
switch to the optimized kernel.

INPUT PARAMETERS
    A   -   real matrix, is passed to ensure that we didn't split
            complex matrix using real splitting subroutine.
            matrix itself is not changed.

  -- ALGLIB routine --
     15.12.2009
     Bochkanov Sergey
*************************************************************************/
int ablasblocksize(const ap::real_2d_array& a)
{
    int result;

    result = 32;
    return result;
}


/*************************************************************************
Block size for complex subroutines.

  -- ALGLIB routine --
     15.12.2009
     Bochkanov Sergey
*************************************************************************/
int ablascomplexblocksize(const ap::complex_2d_array& a)
{
    int result;

    result = 24;
    return result;
}


/*************************************************************************
Microblock size

  -- ALGLIB routine --
     15.12.2009
     Bochkanov Sergey
*************************************************************************/
int ablasmicroblocksize()
{
    int result;

    result = 8;
    return result;
}


/*************************************************************************
Cache-oblivous complex "copy-and-transpose"

Input parameters:
    M   -   number of rows
    N   -   number of columns
    A   -   source matrix, MxN submatrix is copied and transposed
    IA  -   submatrix offset (row index)
    JA  -   submatrix offset (column index)
    A   -   destination matrix
    IB  -   submatrix offset (row index)
    JB  -   submatrix offset (column index)
*************************************************************************/
void cmatrixtranspose(int m,
     int n,
     const ap::complex_2d_array& a,
     int ia,
     int ja,
     ap::complex_2d_array& b,
     int ib,
     int jb)
{
    int i;
    int s1;
    int s2;

    if( m<=2*ablascomplexblocksize(a)&&n<=2*ablascomplexblocksize(a) )
    {
        
        //
        // base case
        //
        for(i = 0; i <= m-1; i++)
        {
            ap::vmove(&b(ib, jb+i), b.getstride(), &a(ia+i, ja), 1, "N", ap::vlen(ib,ib+n-1));
        }
    }
    else
    {
        
        //
        // Cache-oblivious recursion
        //
        if( m>n )
        {
            ablascomplexsplitlength(a, m, s1, s2);
            cmatrixtranspose(s1, n, a, ia, ja, b, ib, jb);
            cmatrixtranspose(s2, n, a, ia+s1, ja, b, ib, jb+s1);
        }
        else
        {
            ablascomplexsplitlength(a, n, s1, s2);
            cmatrixtranspose(m, s1, a, ia, ja, b, ib, jb);
            cmatrixtranspose(m, s2, a, ia, ja+s1, b, ib+s1, jb);
        }
    }
}


/*************************************************************************
Cache-oblivous real "copy-and-transpose"

Input parameters:
    M   -   number of rows
    N   -   number of columns
    A   -   source matrix, MxN submatrix is copied and transposed
    IA  -   submatrix offset (row index)
    JA  -   submatrix offset (column index)
    A   -   destination matrix
    IB  -   submatrix offset (row index)
    JB  -   submatrix offset (column index)
*************************************************************************/
void rmatrixtranspose(int m,
     int n,
     const ap::real_2d_array& a,
     int ia,
     int ja,
     ap::real_2d_array& b,
     int ib,
     int jb)
{
    int i;
    int s1;
    int s2;

    if( m<=2*ablasblocksize(a)&&n<=2*ablasblocksize(a) )
    {
        
        //
        // base case
        //
        for(i = 0; i <= m-1; i++)
        {
            ap::vmove(&b(ib, jb+i), b.getstride(), &a(ia+i, ja), 1, ap::vlen(ib,ib+n-1));
        }
    }
    else
    {
        
        //
        // Cache-oblivious recursion
        //
        if( m>n )
        {
            ablassplitlength(a, m, s1, s2);
            rmatrixtranspose(s1, n, a, ia, ja, b, ib, jb);
            rmatrixtranspose(s2, n, a, ia+s1, ja, b, ib, jb+s1);
        }
        else
        {
            ablassplitlength(a, n, s1, s2);
            rmatrixtranspose(m, s1, a, ia, ja, b, ib, jb);
            rmatrixtranspose(m, s2, a, ia, ja+s1, b, ib+s1, jb);
        }
    }
}


/*************************************************************************
Copy

Input parameters:
    M   -   number of rows
    N   -   number of columns
    A   -   source matrix, MxN submatrix is copied and transposed
    IA  -   submatrix offset (row index)
    JA  -   submatrix offset (column index)
    B   -   destination matrix
    IB  -   submatrix offset (row index)
    JB  -   submatrix offset (column index)
*************************************************************************/
void cmatrixcopy(int m,
     int n,
     const ap::complex_2d_array& a,
     int ia,
     int ja,
     ap::complex_2d_array& b,
     int ib,
     int jb)
{
    int i;

    for(i = 0; i <= m-1; i++)
    {
        ap::vmove(&b(ib+i, jb), 1, &a(ia+i, ja), 1, "N", ap::vlen(jb,jb+n-1));
    }
}


/*************************************************************************
Copy

Input parameters:
    M   -   number of rows
    N   -   number of columns
    A   -   source matrix, MxN submatrix is copied and transposed
    IA  -   submatrix offset (row index)
    JA  -   submatrix offset (column index)
    B   -   destination matrix
    IB  -   submatrix offset (row index)
    JB  -   submatrix offset (column index)
*************************************************************************/
void rmatrixcopy(int m,
     int n,
     const ap::real_2d_array& a,
     int ia,
     int ja,
     ap::real_2d_array& b,
     int ib,
     int jb)
{
    int i;

    for(i = 0; i <= m-1; i++)
    {
        ap::vmove(&b(ib+i, jb), 1, &a(ia+i, ja), 1, ap::vlen(jb,jb+n-1));
    }
}


/*************************************************************************
Rank-1 correction: A := A + u*v'

INPUT PARAMETERS:
    M   -   number of rows
    N   -   number of columns
    A   -   target matrix, MxN submatrix is updated
    IA  -   submatrix offset (row index)
    JA  -   submatrix offset (column index)
    U   -   vector #1
    IU  -   subvector offset
    V   -   vector #2
    IV  -   subvector offset
*************************************************************************/
void cmatrixrank1(int m,
     int n,
     ap::complex_2d_array& a,
     int ia,
     int ja,
     ap::complex_1d_array& u,
     int iu,
     ap::complex_1d_array& v,
     int iv)
{
    int i;
    ap::complex s;

    if( m==0||n==0 )
    {
        return;
    }
    if( cmatrixrank1f(m, n, a, ia, ja, u, iu, v, iv) )
    {
        return;
    }
    for(i = 0; i <= m-1; i++)
    {
        s = u(iu+i);
        ap::vadd(&a(ia+i, ja), 1, &v(iv), 1, "N", ap::vlen(ja,ja+n-1), s);
    }
}


/*************************************************************************
Rank-1 correction: A := A + u*v'

INPUT PARAMETERS:
    M   -   number of rows
    N   -   number of columns
    A   -   target matrix, MxN submatrix is updated
    IA  -   submatrix offset (row index)
    JA  -   submatrix offset (column index)
    U   -   vector #1
    IU  -   subvector offset
    V   -   vector #2
    IV  -   subvector offset
*************************************************************************/
void rmatrixrank1(int m,
     int n,
     ap::real_2d_array& a,
     int ia,
     int ja,
     ap::real_1d_array& u,
     int iu,
     ap::real_1d_array& v,
     int iv)
{
    int i;
    double s;

    if( m==0||n==0 )
    {
        return;
    }
    if( rmatrixrank1f(m, n, a, ia, ja, u, iu, v, iv) )
    {
        return;
    }
    for(i = 0; i <= m-1; i++)
    {
        s = u(iu+i);
        ap::vadd(&a(ia+i, ja), 1, &v(iv), 1, ap::vlen(ja,ja+n-1), s);
    }
}


/*************************************************************************
Matrix-vector product: y := op(A)*x

INPUT PARAMETERS:
    M   -   number of rows of op(A)
            M>=0
    N   -   number of columns of op(A)
            N>=0
    A   -   target matrix
    IA  -   submatrix offset (row index)
    JA  -   submatrix offset (column index)
    OpA -   operation type:
            * OpA=0     =>  op(A) = A
            * OpA=1     =>  op(A) = A^T
            * OpA=2     =>  op(A) = A^H
    X   -   input vector
    IX  -   subvector offset
    IY  -   subvector offset

OUTPUT PARAMETERS:
    Y   -   vector which stores result

if M=0, then subroutine does nothing.
if N=0, Y is filled by zeros.


  -- ALGLIB routine --

     28.01.2010
     Bochkanov Sergey
*************************************************************************/
void cmatrixmv(int m,
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
    int i;
    ap::complex v;

    if( m==0 )
    {
        return;
    }
    if( n==0 )
    {
        for(i = 0; i <= m-1; i++)
        {
            y(iy+i) = 0;
        }
        return;
    }
    if( cmatrixmvf(m, n, a, ia, ja, opa, x, ix, y, iy) )
    {
        return;
    }
    if( opa==0 )
    {
        
        //
        // y = A*x
        //
        for(i = 0; i <= m-1; i++)
        {
            v = ap::vdotproduct(&a(ia+i, ja), 1, "N", &x(ix), 1, "N", ap::vlen(ja,ja+n-1));
            y(iy+i) = v;
        }
        return;
    }
    if( opa==1 )
    {
        
        //
        // y = A^T*x
        //
        for(i = 0; i <= m-1; i++)
        {
            y(iy+i) = 0;
        }
        for(i = 0; i <= n-1; i++)
        {
            v = x(ix+i);
            ap::vadd(&y(iy), 1, &a(ia+i, ja), 1, "N", ap::vlen(iy,iy+m-1), v);
        }
        return;
    }
    if( opa==2 )
    {
        
        //
        // y = A^H*x
        //
        for(i = 0; i <= m-1; i++)
        {
            y(iy+i) = 0;
        }
        for(i = 0; i <= n-1; i++)
        {
            v = x(ix+i);
            ap::vadd(&y(iy), 1, &a(ia+i, ja), 1, "Conj", ap::vlen(iy,iy+m-1), v);
        }
        return;
    }
}


/*************************************************************************
Matrix-vector product: y := op(A)*x

INPUT PARAMETERS:
    M   -   number of rows of op(A)
    N   -   number of columns of op(A)
    A   -   target matrix
    IA  -   submatrix offset (row index)
    JA  -   submatrix offset (column index)
    OpA -   operation type:
            * OpA=0     =>  op(A) = A
            * OpA=1     =>  op(A) = A^T
    X   -   input vector
    IX  -   subvector offset
    IY  -   subvector offset

OUTPUT PARAMETERS:
    Y   -   vector which stores result

if M=0, then subroutine does nothing.
if N=0, Y is filled by zeros.


  -- ALGLIB routine --

     28.01.2010
     Bochkanov Sergey
*************************************************************************/
void rmatrixmv(int m,
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
    int i;
    double v;

    if( m==0 )
    {
        return;
    }
    if( n==0 )
    {
        for(i = 0; i <= m-1; i++)
        {
            y(iy+i) = 0;
        }
        return;
    }
    if( rmatrixmvf(m, n, a, ia, ja, opa, x, ix, y, iy) )
    {
        return;
    }
    if( opa==0 )
    {
        
        //
        // y = A*x
        //
        for(i = 0; i <= m-1; i++)
        {
            v = ap::vdotproduct(&a(ia+i, ja), 1, &x(ix), 1, ap::vlen(ja,ja+n-1));
            y(iy+i) = v;
        }
        return;
    }
    if( opa==1 )
    {
        
        //
        // y = A^T*x
        //
        for(i = 0; i <= m-1; i++)
        {
            y(iy+i) = 0;
        }
        for(i = 0; i <= n-1; i++)
        {
            v = x(ix+i);
            ap::vadd(&y(iy), 1, &a(ia+i, ja), 1, ap::vlen(iy,iy+m-1), v);
        }
        return;
    }
}


/*************************************************************************
This subroutine calculates X*op(A^-1) where:
* X is MxN general matrix
* A is NxN upper/lower triangular/unitriangular matrix
* "op" may be identity transformation, transposition, conjugate transposition

Multiplication result replaces X.
Cache-oblivious algorithm is used.

INPUT PARAMETERS
    N   -   matrix size, N>=0
    M   -   matrix size, N>=0
    A       -   matrix, actial matrix is stored in A[I1:I1+N-1,J1:J1+N-1]
    I1      -   submatrix offset
    J1      -   submatrix offset
    IsUpper -   whether matrix is upper triangular
    IsUnit  -   whether matrix is unitriangular
    OpType  -   transformation type:
                * 0 - no transformation
                * 1 - transposition
                * 2 - conjugate transposition
    C   -   matrix, actial matrix is stored in C[I2:I2+M-1,J2:J2+N-1]
    I2  -   submatrix offset
    J2  -   submatrix offset

  -- ALGLIB routine --
     15.12.2009
     Bochkanov Sergey
*************************************************************************/
void cmatrixrighttrsm(int m,
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
    int s1;
    int s2;
    int bs;

    bs = ablascomplexblocksize(a);
    if( m<=bs&&n<=bs )
    {
        cmatrixrighttrsm2(m, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
        return;
    }
    if( m>=n )
    {
        
        //
        // Split X: X*A = (X1 X2)^T*A
        //
        ablascomplexsplitlength(a, m, s1, s2);
        cmatrixrighttrsm(s1, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
        cmatrixrighttrsm(s2, n, a, i1, j1, isupper, isunit, optype, x, i2+s1, j2);
    }
    else
    {
        
        //
        // Split A:
        //               (A1  A12)
        // X*op(A) = X*op(       )
        //               (     A2)
        //
        // Different variants depending on
        // IsUpper/OpType combinations
        //
        ablascomplexsplitlength(a, n, s1, s2);
        if( isupper&&optype==0 )
        {
            
            //
            //                  (A1  A12)-1
            // X*A^-1 = (X1 X2)*(       )
            //                  (     A2)
            //
            cmatrixrighttrsm(m, s1, a, i1, j1, isupper, isunit, optype, x, i2, j2);
            cmatrixgemm(m, s2, s1, -1.0, x, i2, j2, 0, a, i1, j1+s1, 0, 1.0, x, i2, j2+s1);
            cmatrixrighttrsm(m, s2, a, i1+s1, j1+s1, isupper, isunit, optype, x, i2, j2+s1);
            return;
        }
        if( isupper&&optype!=0 )
        {
            
            //
            //                  (A1'     )-1
            // X*A^-1 = (X1 X2)*(        )
            //                  (A12' A2')
            //
            cmatrixrighttrsm(m, s2, a, i1+s1, j1+s1, isupper, isunit, optype, x, i2, j2+s1);
            cmatrixgemm(m, s1, s2, -1.0, x, i2, j2+s1, 0, a, i1, j1+s1, optype, 1.0, x, i2, j2);
            cmatrixrighttrsm(m, s1, a, i1, j1, isupper, isunit, optype, x, i2, j2);
            return;
        }
        if( !isupper&&optype==0 )
        {
            
            //
            //                  (A1     )-1
            // X*A^-1 = (X1 X2)*(       )
            //                  (A21  A2)
            //
            cmatrixrighttrsm(m, s2, a, i1+s1, j1+s1, isupper, isunit, optype, x, i2, j2+s1);
            cmatrixgemm(m, s1, s2, -1.0, x, i2, j2+s1, 0, a, i1+s1, j1, 0, 1.0, x, i2, j2);
            cmatrixrighttrsm(m, s1, a, i1, j1, isupper, isunit, optype, x, i2, j2);
            return;
        }
        if( !isupper&&optype!=0 )
        {
            
            //
            //                  (A1' A21')-1
            // X*A^-1 = (X1 X2)*(        )
            //                  (     A2')
            //
            cmatrixrighttrsm(m, s1, a, i1, j1, isupper, isunit, optype, x, i2, j2);
            cmatrixgemm(m, s2, s1, -1.0, x, i2, j2, 0, a, i1+s1, j1, optype, 1.0, x, i2, j2+s1);
            cmatrixrighttrsm(m, s2, a, i1+s1, j1+s1, isupper, isunit, optype, x, i2, j2+s1);
            return;
        }
    }
}


/*************************************************************************
This subroutine calculates op(A^-1)*X where:
* X is MxN general matrix
* A is MxM upper/lower triangular/unitriangular matrix
* "op" may be identity transformation, transposition, conjugate transposition

Multiplication result replaces X.
Cache-oblivious algorithm is used.

INPUT PARAMETERS
    N   -   matrix size, N>=0
    M   -   matrix size, N>=0
    A       -   matrix, actial matrix is stored in A[I1:I1+M-1,J1:J1+M-1]
    I1      -   submatrix offset
    J1      -   submatrix offset
    IsUpper -   whether matrix is upper triangular
    IsUnit  -   whether matrix is unitriangular
    OpType  -   transformation type:
                * 0 - no transformation
                * 1 - transposition
                * 2 - conjugate transposition
    C   -   matrix, actial matrix is stored in C[I2:I2+M-1,J2:J2+N-1]
    I2  -   submatrix offset
    J2  -   submatrix offset

  -- ALGLIB routine --
     15.12.2009
     Bochkanov Sergey
*************************************************************************/
void cmatrixlefttrsm(int m,
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
    int s1;
    int s2;
    int bs;

    bs = ablascomplexblocksize(a);
    if( m<=bs&&n<=bs )
    {
        cmatrixlefttrsm2(m, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
        return;
    }
    if( n>=m )
    {
        
        //
        // Split X: op(A)^-1*X = op(A)^-1*(X1 X2)
        //
        ablascomplexsplitlength(x, n, s1, s2);
        cmatrixlefttrsm(m, s1, a, i1, j1, isupper, isunit, optype, x, i2, j2);
        cmatrixlefttrsm(m, s2, a, i1, j1, isupper, isunit, optype, x, i2, j2+s1);
    }
    else
    {
        
        //
        // Split A
        //
        ablascomplexsplitlength(a, m, s1, s2);
        if( isupper&&optype==0 )
        {
            
            //
            //           (A1  A12)-1  ( X1 )
            // A^-1*X* = (       )   *(    )
            //           (     A2)    ( X2 )
            //
            cmatrixlefttrsm(s2, n, a, i1+s1, j1+s1, isupper, isunit, optype, x, i2+s1, j2);
            cmatrixgemm(s1, n, s2, -1.0, a, i1, j1+s1, 0, x, i2+s1, j2, 0, 1.0, x, i2, j2);
            cmatrixlefttrsm(s1, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
            return;
        }
        if( isupper&&optype!=0 )
        {
            
            //
            //          (A1'     )-1 ( X1 )
            // A^-1*X = (        )  *(    )
            //          (A12' A2')   ( X2 )
            //
            cmatrixlefttrsm(s1, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
            cmatrixgemm(s2, n, s1, -1.0, a, i1, j1+s1, optype, x, i2, j2, 0, 1.0, x, i2+s1, j2);
            cmatrixlefttrsm(s2, n, a, i1+s1, j1+s1, isupper, isunit, optype, x, i2+s1, j2);
            return;
        }
        if( !isupper&&optype==0 )
        {
            
            //
            //          (A1     )-1 ( X1 )
            // A^-1*X = (       )  *(    )
            //          (A21  A2)   ( X2 )
            //
            cmatrixlefttrsm(s1, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
            cmatrixgemm(s2, n, s1, -1.0, a, i1+s1, j1, 0, x, i2, j2, 0, 1.0, x, i2+s1, j2);
            cmatrixlefttrsm(s2, n, a, i1+s1, j1+s1, isupper, isunit, optype, x, i2+s1, j2);
            return;
        }
        if( !isupper&&optype!=0 )
        {
            
            //
            //          (A1' A21')-1 ( X1 )
            // A^-1*X = (        )  *(    )
            //          (     A2')   ( X2 )
            //
            cmatrixlefttrsm(s2, n, a, i1+s1, j1+s1, isupper, isunit, optype, x, i2+s1, j2);
            cmatrixgemm(s1, n, s2, -1.0, a, i1+s1, j1, optype, x, i2+s1, j2, 0, 1.0, x, i2, j2);
            cmatrixlefttrsm(s1, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
            return;
        }
    }
}


/*************************************************************************
Same as CMatrixRightTRSM, but for real matrices

OpType may be only 0 or 1.

  -- ALGLIB routine --
     15.12.2009
     Bochkanov Sergey
*************************************************************************/
void rmatrixrighttrsm(int m,
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
    int s1;
    int s2;
    int bs;

    bs = ablasblocksize(a);
    if( m<=bs&&n<=bs )
    {
        rmatrixrighttrsm2(m, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
        return;
    }
    if( m>=n )
    {
        
        //
        // Split X: X*A = (X1 X2)^T*A
        //
        ablassplitlength(a, m, s1, s2);
        rmatrixrighttrsm(s1, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
        rmatrixrighttrsm(s2, n, a, i1, j1, isupper, isunit, optype, x, i2+s1, j2);
    }
    else
    {
        
        //
        // Split A:
        //               (A1  A12)
        // X*op(A) = X*op(       )
        //               (     A2)
        //
        // Different variants depending on
        // IsUpper/OpType combinations
        //
        ablassplitlength(a, n, s1, s2);
        if( isupper&&optype==0 )
        {
            
            //
            //                  (A1  A12)-1
            // X*A^-1 = (X1 X2)*(       )
            //                  (     A2)
            //
            rmatrixrighttrsm(m, s1, a, i1, j1, isupper, isunit, optype, x, i2, j2);
            rmatrixgemm(m, s2, s1, -1.0, x, i2, j2, 0, a, i1, j1+s1, 0, 1.0, x, i2, j2+s1);
            rmatrixrighttrsm(m, s2, a, i1+s1, j1+s1, isupper, isunit, optype, x, i2, j2+s1);
            return;
        }
        if( isupper&&optype!=0 )
        {
            
            //
            //                  (A1'     )-1
            // X*A^-1 = (X1 X2)*(        )
            //                  (A12' A2')
            //
            rmatrixrighttrsm(m, s2, a, i1+s1, j1+s1, isupper, isunit, optype, x, i2, j2+s1);
            rmatrixgemm(m, s1, s2, -1.0, x, i2, j2+s1, 0, a, i1, j1+s1, optype, 1.0, x, i2, j2);
            rmatrixrighttrsm(m, s1, a, i1, j1, isupper, isunit, optype, x, i2, j2);
            return;
        }
        if( !isupper&&optype==0 )
        {
            
            //
            //                  (A1     )-1
            // X*A^-1 = (X1 X2)*(       )
            //                  (A21  A2)
            //
            rmatrixrighttrsm(m, s2, a, i1+s1, j1+s1, isupper, isunit, optype, x, i2, j2+s1);
            rmatrixgemm(m, s1, s2, -1.0, x, i2, j2+s1, 0, a, i1+s1, j1, 0, 1.0, x, i2, j2);
            rmatrixrighttrsm(m, s1, a, i1, j1, isupper, isunit, optype, x, i2, j2);
            return;
        }
        if( !isupper&&optype!=0 )
        {
            
            //
            //                  (A1' A21')-1
            // X*A^-1 = (X1 X2)*(        )
            //                  (     A2')
            //
            rmatrixrighttrsm(m, s1, a, i1, j1, isupper, isunit, optype, x, i2, j2);
            rmatrixgemm(m, s2, s1, -1.0, x, i2, j2, 0, a, i1+s1, j1, optype, 1.0, x, i2, j2+s1);
            rmatrixrighttrsm(m, s2, a, i1+s1, j1+s1, isupper, isunit, optype, x, i2, j2+s1);
            return;
        }
    }
}


/*************************************************************************
Same as CMatrixLeftTRSM, but for real matrices

OpType may be only 0 or 1.

  -- ALGLIB routine --
     15.12.2009
     Bochkanov Sergey
*************************************************************************/
void rmatrixlefttrsm(int m,
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
    int s1;
    int s2;
    int bs;

    bs = ablasblocksize(a);
    if( m<=bs&&n<=bs )
    {
        rmatrixlefttrsm2(m, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
        return;
    }
    if( n>=m )
    {
        
        //
        // Split X: op(A)^-1*X = op(A)^-1*(X1 X2)
        //
        ablassplitlength(x, n, s1, s2);
        rmatrixlefttrsm(m, s1, a, i1, j1, isupper, isunit, optype, x, i2, j2);
        rmatrixlefttrsm(m, s2, a, i1, j1, isupper, isunit, optype, x, i2, j2+s1);
    }
    else
    {
        
        //
        // Split A
        //
        ablassplitlength(a, m, s1, s2);
        if( isupper&&optype==0 )
        {
            
            //
            //           (A1  A12)-1  ( X1 )
            // A^-1*X* = (       )   *(    )
            //           (     A2)    ( X2 )
            //
            rmatrixlefttrsm(s2, n, a, i1+s1, j1+s1, isupper, isunit, optype, x, i2+s1, j2);
            rmatrixgemm(s1, n, s2, -1.0, a, i1, j1+s1, 0, x, i2+s1, j2, 0, 1.0, x, i2, j2);
            rmatrixlefttrsm(s1, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
            return;
        }
        if( isupper&&optype!=0 )
        {
            
            //
            //          (A1'     )-1 ( X1 )
            // A^-1*X = (        )  *(    )
            //          (A12' A2')   ( X2 )
            //
            rmatrixlefttrsm(s1, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
            rmatrixgemm(s2, n, s1, -1.0, a, i1, j1+s1, optype, x, i2, j2, 0, 1.0, x, i2+s1, j2);
            rmatrixlefttrsm(s2, n, a, i1+s1, j1+s1, isupper, isunit, optype, x, i2+s1, j2);
            return;
        }
        if( !isupper&&optype==0 )
        {
            
            //
            //          (A1     )-1 ( X1 )
            // A^-1*X = (       )  *(    )
            //          (A21  A2)   ( X2 )
            //
            rmatrixlefttrsm(s1, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
            rmatrixgemm(s2, n, s1, -1.0, a, i1+s1, j1, 0, x, i2, j2, 0, 1.0, x, i2+s1, j2);
            rmatrixlefttrsm(s2, n, a, i1+s1, j1+s1, isupper, isunit, optype, x, i2+s1, j2);
            return;
        }
        if( !isupper&&optype!=0 )
        {
            
            //
            //          (A1' A21')-1 ( X1 )
            // A^-1*X = (        )  *(    )
            //          (     A2')   ( X2 )
            //
            rmatrixlefttrsm(s2, n, a, i1+s1, j1+s1, isupper, isunit, optype, x, i2+s1, j2);
            rmatrixgemm(s1, n, s2, -1.0, a, i1+s1, j1, optype, x, i2+s1, j2, 0, 1.0, x, i2, j2);
            rmatrixlefttrsm(s1, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
            return;
        }
    }
}


/*************************************************************************
This subroutine calculates  C=alpha*A*A^H+beta*C  or  C=alpha*A^H*A+beta*C
where:
* C is NxN Hermitian matrix given by its upper/lower triangle
* A is NxK matrix when A*A^H is calculated, KxN matrix otherwise

Additional info:
* cache-oblivious algorithm is used.
* multiplication result replaces C. If Beta=0, C elements are not used in
  calculations (not multiplied by zero - just not referenced)
* if Alpha=0, A is not used (not multiplied by zero - just not referenced)
* if both Beta and Alpha are zero, C is filled by zeros.

INPUT PARAMETERS
    N       -   matrix size, N>=0
    K       -   matrix size, K>=0
    Alpha   -   coefficient
    A       -   matrix
    IA      -   submatrix offset
    JA      -   submatrix offset
    OpTypeA -   multiplication type:
                * 0 - A*A^H is calculated
                * 2 - A^H*A is calculated
    Beta    -   coefficient
    C       -   matrix
    IC      -   submatrix offset
    JC      -   submatrix offset
    IsUpper -   whether C is upper triangular or lower triangular

  -- ALGLIB routine --
     16.12.2009
     Bochkanov Sergey
*************************************************************************/
void cmatrixsyrk(int n,
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
    int s1;
    int s2;
    int bs;

    bs = ablascomplexblocksize(a);
    if( n<=bs&&k<=bs )
    {
        cmatrixsyrk2(n, k, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper);
        return;
    }
    if( k>=n )
    {
        
        //
        // Split K
        //
        ablascomplexsplitlength(a, k, s1, s2);
        if( optypea==0 )
        {
            cmatrixsyrk(n, s1, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper);
            cmatrixsyrk(n, s2, alpha, a, ia, ja+s1, optypea, 1.0, c, ic, jc, isupper);
        }
        else
        {
            cmatrixsyrk(n, s1, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper);
            cmatrixsyrk(n, s2, alpha, a, ia+s1, ja, optypea, 1.0, c, ic, jc, isupper);
        }
    }
    else
    {
        
        //
        // Split N
        //
        ablascomplexsplitlength(a, n, s1, s2);
        if( optypea==0&&isupper )
        {
            cmatrixsyrk(s1, k, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper);
            cmatrixgemm(s1, s2, k, alpha, a, ia, ja, 0, a, ia+s1, ja, 2, beta, c, ic, jc+s1);
            cmatrixsyrk(s2, k, alpha, a, ia+s1, ja, optypea, beta, c, ic+s1, jc+s1, isupper);
            return;
        }
        if( optypea==0&&!isupper )
        {
            cmatrixsyrk(s1, k, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper);
            cmatrixgemm(s2, s1, k, alpha, a, ia+s1, ja, 0, a, ia, ja, 2, beta, c, ic+s1, jc);
            cmatrixsyrk(s2, k, alpha, a, ia+s1, ja, optypea, beta, c, ic+s1, jc+s1, isupper);
            return;
        }
        if( optypea!=0&&isupper )
        {
            cmatrixsyrk(s1, k, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper);
            cmatrixgemm(s1, s2, k, alpha, a, ia, ja, 2, a, ia, ja+s1, 0, beta, c, ic, jc+s1);
            cmatrixsyrk(s2, k, alpha, a, ia, ja+s1, optypea, beta, c, ic+s1, jc+s1, isupper);
            return;
        }
        if( optypea!=0&&!isupper )
        {
            cmatrixsyrk(s1, k, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper);
            cmatrixgemm(s2, s1, k, alpha, a, ia, ja+s1, 2, a, ia, ja, 0, beta, c, ic+s1, jc);
            cmatrixsyrk(s2, k, alpha, a, ia, ja+s1, optypea, beta, c, ic+s1, jc+s1, isupper);
            return;
        }
    }
}


/*************************************************************************
Same as CMatrixSYRK, but for real matrices

OpType may be only 0 or 1.

  -- ALGLIB routine --
     16.12.2009
     Bochkanov Sergey
*************************************************************************/
void rmatrixsyrk(int n,
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
    int s1;
    int s2;
    int bs;

    bs = ablasblocksize(a);
    if( n<=bs&&k<=bs )
    {
        rmatrixsyrk2(n, k, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper);
        return;
    }
    if( k>=n )
    {
        
        //
        // Split K
        //
        ablassplitlength(a, k, s1, s2);
        if( optypea==0 )
        {
            rmatrixsyrk(n, s1, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper);
            rmatrixsyrk(n, s2, alpha, a, ia, ja+s1, optypea, 1.0, c, ic, jc, isupper);
        }
        else
        {
            rmatrixsyrk(n, s1, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper);
            rmatrixsyrk(n, s2, alpha, a, ia+s1, ja, optypea, 1.0, c, ic, jc, isupper);
        }
    }
    else
    {
        
        //
        // Split N
        //
        ablassplitlength(a, n, s1, s2);
        if( optypea==0&&isupper )
        {
            rmatrixsyrk(s1, k, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper);
            rmatrixgemm(s1, s2, k, alpha, a, ia, ja, 0, a, ia+s1, ja, 1, beta, c, ic, jc+s1);
            rmatrixsyrk(s2, k, alpha, a, ia+s1, ja, optypea, beta, c, ic+s1, jc+s1, isupper);
            return;
        }
        if( optypea==0&&!isupper )
        {
            rmatrixsyrk(s1, k, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper);
            rmatrixgemm(s2, s1, k, alpha, a, ia+s1, ja, 0, a, ia, ja, 1, beta, c, ic+s1, jc);
            rmatrixsyrk(s2, k, alpha, a, ia+s1, ja, optypea, beta, c, ic+s1, jc+s1, isupper);
            return;
        }
        if( optypea!=0&&isupper )
        {
            rmatrixsyrk(s1, k, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper);
            rmatrixgemm(s1, s2, k, alpha, a, ia, ja, 1, a, ia, ja+s1, 0, beta, c, ic, jc+s1);
            rmatrixsyrk(s2, k, alpha, a, ia, ja+s1, optypea, beta, c, ic+s1, jc+s1, isupper);
            return;
        }
        if( optypea!=0&&!isupper )
        {
            rmatrixsyrk(s1, k, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper);
            rmatrixgemm(s2, s1, k, alpha, a, ia, ja+s1, 1, a, ia, ja, 0, beta, c, ic+s1, jc);
            rmatrixsyrk(s2, k, alpha, a, ia, ja+s1, optypea, beta, c, ic+s1, jc+s1, isupper);
            return;
        }
    }
}


/*************************************************************************
This subroutine calculates C = alpha*op1(A)*op2(B) +beta*C where:
* C is MxN general matrix
* op1(A) is MxK matrix
* op2(B) is KxN matrix
* "op" may be identity transformation, transposition, conjugate transposition

Additional info:
* cache-oblivious algorithm is used.
* multiplication result replaces C. If Beta=0, C elements are not used in
  calculations (not multiplied by zero - just not referenced)
* if Alpha=0, A is not used (not multiplied by zero - just not referenced)
* if both Beta and Alpha are zero, C is filled by zeros.

INPUT PARAMETERS
    N       -   matrix size, N>0
    M       -   matrix size, N>0
    K       -   matrix size, K>0
    Alpha   -   coefficient
    A       -   matrix
    IA      -   submatrix offset
    JA      -   submatrix offset
    OpTypeA -   transformation type:
                * 0 - no transformation
                * 1 - transposition
                * 2 - conjugate transposition
    B       -   matrix
    IB      -   submatrix offset
    JB      -   submatrix offset
    OpTypeB -   transformation type:
                * 0 - no transformation
                * 1 - transposition
                * 2 - conjugate transposition
    Beta    -   coefficient
    C       -   matrix
    IC      -   submatrix offset
    JC      -   submatrix offset

  -- ALGLIB routine --
     16.12.2009
     Bochkanov Sergey
*************************************************************************/
void cmatrixgemm(int m,
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
    int s1;
    int s2;
    int bs;

    bs = ablascomplexblocksize(a);
    if( m<=bs&&n<=bs&&k<=bs )
    {
        cmatrixgemmk(m, n, k, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc);
        return;
    }
    if( m>=n&&m>=k )
    {
        
        //
        // A*B = (A1 A2)^T*B
        //
        ablascomplexsplitlength(a, m, s1, s2);
        cmatrixgemm(s1, n, k, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc);
        if( optypea==0 )
        {
            cmatrixgemm(s2, n, k, alpha, a, ia+s1, ja, optypea, b, ib, jb, optypeb, beta, c, ic+s1, jc);
        }
        else
        {
            cmatrixgemm(s2, n, k, alpha, a, ia, ja+s1, optypea, b, ib, jb, optypeb, beta, c, ic+s1, jc);
        }
        return;
    }
    if( n>=m&&n>=k )
    {
        
        //
        // A*B = A*(B1 B2)
        //
        ablascomplexsplitlength(a, n, s1, s2);
        if( optypeb==0 )
        {
            cmatrixgemm(m, s1, k, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc);
            cmatrixgemm(m, s2, k, alpha, a, ia, ja, optypea, b, ib, jb+s1, optypeb, beta, c, ic, jc+s1);
        }
        else
        {
            cmatrixgemm(m, s1, k, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc);
            cmatrixgemm(m, s2, k, alpha, a, ia, ja, optypea, b, ib+s1, jb, optypeb, beta, c, ic, jc+s1);
        }
        return;
    }
    if( k>=m&&k>=n )
    {
        
        //
        // A*B = (A1 A2)*(B1 B2)^T
        //
        ablascomplexsplitlength(a, k, s1, s2);
        if( optypea==0&&optypeb==0 )
        {
            cmatrixgemm(m, n, s1, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc);
            cmatrixgemm(m, n, s2, alpha, a, ia, ja+s1, optypea, b, ib+s1, jb, optypeb, 1.0, c, ic, jc);
        }
        if( optypea==0&&optypeb!=0 )
        {
            cmatrixgemm(m, n, s1, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc);
            cmatrixgemm(m, n, s2, alpha, a, ia, ja+s1, optypea, b, ib, jb+s1, optypeb, 1.0, c, ic, jc);
        }
        if( optypea!=0&&optypeb==0 )
        {
            cmatrixgemm(m, n, s1, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc);
            cmatrixgemm(m, n, s2, alpha, a, ia+s1, ja, optypea, b, ib+s1, jb, optypeb, 1.0, c, ic, jc);
        }
        if( optypea!=0&&optypeb!=0 )
        {
            cmatrixgemm(m, n, s1, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc);
            cmatrixgemm(m, n, s2, alpha, a, ia+s1, ja, optypea, b, ib, jb+s1, optypeb, 1.0, c, ic, jc);
        }
        return;
    }
}


/*************************************************************************
Same as CMatrixGEMM, but for real numbers.
OpType may be only 0 or 1.

  -- ALGLIB routine --
     16.12.2009
     Bochkanov Sergey
*************************************************************************/
void rmatrixgemm(int m,
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
    int s1;
    int s2;
    int bs;

    bs = ablasblocksize(a);
    if( m<=bs&&n<=bs&&k<=bs )
    {
        rmatrixgemmk(m, n, k, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc);
        return;
    }
    if( m>=n&&m>=k )
    {
        
        //
        // A*B = (A1 A2)^T*B
        //
        ablassplitlength(a, m, s1, s2);
        if( optypea==0 )
        {
            rmatrixgemm(s1, n, k, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc);
            rmatrixgemm(s2, n, k, alpha, a, ia+s1, ja, optypea, b, ib, jb, optypeb, beta, c, ic+s1, jc);
        }
        else
        {
            rmatrixgemm(s1, n, k, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc);
            rmatrixgemm(s2, n, k, alpha, a, ia, ja+s1, optypea, b, ib, jb, optypeb, beta, c, ic+s1, jc);
        }
        return;
    }
    if( n>=m&&n>=k )
    {
        
        //
        // A*B = A*(B1 B2)
        //
        ablassplitlength(a, n, s1, s2);
        if( optypeb==0 )
        {
            rmatrixgemm(m, s1, k, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc);
            rmatrixgemm(m, s2, k, alpha, a, ia, ja, optypea, b, ib, jb+s1, optypeb, beta, c, ic, jc+s1);
        }
        else
        {
            rmatrixgemm(m, s1, k, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc);
            rmatrixgemm(m, s2, k, alpha, a, ia, ja, optypea, b, ib+s1, jb, optypeb, beta, c, ic, jc+s1);
        }
        return;
    }
    if( k>=m&&k>=n )
    {
        
        //
        // A*B = (A1 A2)*(B1 B2)^T
        //
        ablassplitlength(a, k, s1, s2);
        if( optypea==0&&optypeb==0 )
        {
            rmatrixgemm(m, n, s1, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc);
            rmatrixgemm(m, n, s2, alpha, a, ia, ja+s1, optypea, b, ib+s1, jb, optypeb, 1.0, c, ic, jc);
        }
        if( optypea==0&&optypeb!=0 )
        {
            rmatrixgemm(m, n, s1, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc);
            rmatrixgemm(m, n, s2, alpha, a, ia, ja+s1, optypea, b, ib, jb+s1, optypeb, 1.0, c, ic, jc);
        }
        if( optypea!=0&&optypeb==0 )
        {
            rmatrixgemm(m, n, s1, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc);
            rmatrixgemm(m, n, s2, alpha, a, ia+s1, ja, optypea, b, ib+s1, jb, optypeb, 1.0, c, ic, jc);
        }
        if( optypea!=0&&optypeb!=0 )
        {
            rmatrixgemm(m, n, s1, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc);
            rmatrixgemm(m, n, s2, alpha, a, ia+s1, ja, optypea, b, ib, jb+s1, optypeb, 1.0, c, ic, jc);
        }
        return;
    }
}


/*************************************************************************
Complex ABLASSplitLength

  -- ALGLIB routine --
     15.12.2009
     Bochkanov Sergey
*************************************************************************/
static void ablasinternalsplitlength(int n, int nb, int& n1, int& n2)
{
    int r;

    if( n<=nb )
    {
        
        //
        // Block size, no further splitting
        //
        n1 = n;
        n2 = 0;
    }
    else
    {
        
        //
        // Greater than block size
        //
        if( n%nb!=0 )
        {
            
            //
            // Split remainder
            //
            n2 = n%nb;
            n1 = n-n2;
        }
        else
        {
            
            //
            // Split on block boundaries
            //
            n2 = n/2;
            n1 = n-n2;
            if( n1%nb==0 )
            {
                return;
            }
            r = nb-n1%nb;
            n1 = n1+r;
            n2 = n2-r;
        }
    }
}


/*************************************************************************
Level 2 variant of CMatrixRightTRSM
*************************************************************************/
static void cmatrixrighttrsm2(int m,
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
    int i;
    int j;
    ap::complex vc;
    ap::complex vd;

    
    //
    // Special case
    //
    if( n*m==0 )
    {
        return;
    }
    
    //
    // Try to call fast TRSM
    //
    if( cmatrixrighttrsmf(m, n, a, i1, j1, isupper, isunit, optype, x, i2, j2) )
    {
        return;
    }
    
    //
    // General case
    //
    if( isupper )
    {
        
        //
        // Upper triangular matrix
        //
        if( optype==0 )
        {
            
            //
            // X*A^(-1)
            //
            for(i = 0; i <= m-1; i++)
            {
                for(j = 0; j <= n-1; j++)
                {
                    if( isunit )
                    {
                        vd = 1;
                    }
                    else
                    {
                        vd = a(i1+j,j1+j);
                    }
                    x(i2+i,j2+j) = x(i2+i,j2+j)/vd;
                    if( j<n-1 )
                    {
                        vc = x(i2+i,j2+j);
                        ap::vsub(&x(i2+i, j2+j+1), 1, &a(i1+j, j1+j+1), 1, "N", ap::vlen(j2+j+1,j2+n-1), vc);
                    }
                }
            }
            return;
        }
        if( optype==1 )
        {
            
            //
            // X*A^(-T)
            //
            for(i = 0; i <= m-1; i++)
            {
                for(j = n-1; j >= 0; j--)
                {
                    vc = 0;
                    vd = 1;
                    if( j<n-1 )
                    {
                        vc = ap::vdotproduct(&x(i2+i, j2+j+1), 1, "N", &a(i1+j, j1+j+1), 1, "N", ap::vlen(j2+j+1,j2+n-1));
                    }
                    if( !isunit )
                    {
                        vd = a(i1+j,j1+j);
                    }
                    x(i2+i,j2+j) = (x(i2+i,j2+j)-vc)/vd;
                }
            }
            return;
        }
        if( optype==2 )
        {
            
            //
            // X*A^(-H)
            //
            for(i = 0; i <= m-1; i++)
            {
                for(j = n-1; j >= 0; j--)
                {
                    vc = 0;
                    vd = 1;
                    if( j<n-1 )
                    {
                        vc = ap::vdotproduct(&x(i2+i, j2+j+1), 1, "N", &a(i1+j, j1+j+1), 1, "Conj", ap::vlen(j2+j+1,j2+n-1));
                    }
                    if( !isunit )
                    {
                        vd = ap::conj(a(i1+j,j1+j));
                    }
                    x(i2+i,j2+j) = (x(i2+i,j2+j)-vc)/vd;
                }
            }
            return;
        }
    }
    else
    {
        
        //
        // Lower triangular matrix
        //
        if( optype==0 )
        {
            
            //
            // X*A^(-1)
            //
            for(i = 0; i <= m-1; i++)
            {
                for(j = n-1; j >= 0; j--)
                {
                    if( isunit )
                    {
                        vd = 1;
                    }
                    else
                    {
                        vd = a(i1+j,j1+j);
                    }
                    x(i2+i,j2+j) = x(i2+i,j2+j)/vd;
                    if( j>0 )
                    {
                        vc = x(i2+i,j2+j);
                        ap::vsub(&x(i2+i, j2), 1, &a(i1+j, j1), 1, "N", ap::vlen(j2,j2+j-1), vc);
                    }
                }
            }
            return;
        }
        if( optype==1 )
        {
            
            //
            // X*A^(-T)
            //
            for(i = 0; i <= m-1; i++)
            {
                for(j = 0; j <= n-1; j++)
                {
                    vc = 0;
                    vd = 1;
                    if( j>0 )
                    {
                        vc = ap::vdotproduct(&x(i2+i, j2), 1, "N", &a(i1+j, j1), 1, "N", ap::vlen(j2,j2+j-1));
                    }
                    if( !isunit )
                    {
                        vd = a(i1+j,j1+j);
                    }
                    x(i2+i,j2+j) = (x(i2+i,j2+j)-vc)/vd;
                }
            }
            return;
        }
        if( optype==2 )
        {
            
            //
            // X*A^(-H)
            //
            for(i = 0; i <= m-1; i++)
            {
                for(j = 0; j <= n-1; j++)
                {
                    vc = 0;
                    vd = 1;
                    if( j>0 )
                    {
                        vc = ap::vdotproduct(&x(i2+i, j2), 1, "N", &a(i1+j, j1), 1, "Conj", ap::vlen(j2,j2+j-1));
                    }
                    if( !isunit )
                    {
                        vd = ap::conj(a(i1+j,j1+j));
                    }
                    x(i2+i,j2+j) = (x(i2+i,j2+j)-vc)/vd;
                }
            }
            return;
        }
    }
}


/*************************************************************************
Level-2 subroutine
*************************************************************************/
static void cmatrixlefttrsm2(int m,
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
    int i;
    int j;
    ap::complex vc;
    ap::complex vd;

    
    //
    // Special case
    //
    if( n*m==0 )
    {
        return;
    }
    
    //
    // Try to call fast TRSM
    //
    if( cmatrixlefttrsmf(m, n, a, i1, j1, isupper, isunit, optype, x, i2, j2) )
    {
        return;
    }
    
    //
    // General case
    //
    if( isupper )
    {
        
        //
        // Upper triangular matrix
        //
        if( optype==0 )
        {
            
            //
            // A^(-1)*X
            //
            for(i = m-1; i >= 0; i--)
            {
                for(j = i+1; j <= m-1; j++)
                {
                    vc = a(i1+i,j1+j);
                    ap::vsub(&x(i2+i, j2), 1, &x(i2+j, j2), 1, "N", ap::vlen(j2,j2+n-1), vc);
                }
                if( !isunit )
                {
                    vd = 1/a(i1+i,j1+i);
                    ap::vmul(&x(i2+i, j2), 1, ap::vlen(j2,j2+n-1), vd);
                }
            }
            return;
        }
        if( optype==1 )
        {
            
            //
            // A^(-T)*X
            //
            for(i = 0; i <= m-1; i++)
            {
                if( isunit )
                {
                    vd = 1;
                }
                else
                {
                    vd = 1/a(i1+i,j1+i);
                }
                ap::vmul(&x(i2+i, j2), 1, ap::vlen(j2,j2+n-1), vd);
                for(j = i+1; j <= m-1; j++)
                {
                    vc = a(i1+i,j1+j);
                    ap::vsub(&x(i2+j, j2), 1, &x(i2+i, j2), 1, "N", ap::vlen(j2,j2+n-1), vc);
                }
            }
            return;
        }
        if( optype==2 )
        {
            
            //
            // A^(-H)*X
            //
            for(i = 0; i <= m-1; i++)
            {
                if( isunit )
                {
                    vd = 1;
                }
                else
                {
                    vd = 1/ap::conj(a(i1+i,j1+i));
                }
                ap::vmul(&x(i2+i, j2), 1, ap::vlen(j2,j2+n-1), vd);
                for(j = i+1; j <= m-1; j++)
                {
                    vc = ap::conj(a(i1+i,j1+j));
                    ap::vsub(&x(i2+j, j2), 1, &x(i2+i, j2), 1, "N", ap::vlen(j2,j2+n-1), vc);
                }
            }
            return;
        }
    }
    else
    {
        
        //
        // Lower triangular matrix
        //
        if( optype==0 )
        {
            
            //
            // A^(-1)*X
            //
            for(i = 0; i <= m-1; i++)
            {
                for(j = 0; j <= i-1; j++)
                {
                    vc = a(i1+i,j1+j);
                    ap::vsub(&x(i2+i, j2), 1, &x(i2+j, j2), 1, "N", ap::vlen(j2,j2+n-1), vc);
                }
                if( isunit )
                {
                    vd = 1;
                }
                else
                {
                    vd = 1/a(i1+j,j1+j);
                }
                ap::vmul(&x(i2+i, j2), 1, ap::vlen(j2,j2+n-1), vd);
            }
            return;
        }
        if( optype==1 )
        {
            
            //
            // A^(-T)*X
            //
            for(i = m-1; i >= 0; i--)
            {
                if( isunit )
                {
                    vd = 1;
                }
                else
                {
                    vd = 1/a(i1+i,j1+i);
                }
                ap::vmul(&x(i2+i, j2), 1, ap::vlen(j2,j2+n-1), vd);
                for(j = i-1; j >= 0; j--)
                {
                    vc = a(i1+i,j1+j);
                    ap::vsub(&x(i2+j, j2), 1, &x(i2+i, j2), 1, "N", ap::vlen(j2,j2+n-1), vc);
                }
            }
            return;
        }
        if( optype==2 )
        {
            
            //
            // A^(-H)*X
            //
            for(i = m-1; i >= 0; i--)
            {
                if( isunit )
                {
                    vd = 1;
                }
                else
                {
                    vd = 1/ap::conj(a(i1+i,j1+i));
                }
                ap::vmul(&x(i2+i, j2), 1, ap::vlen(j2,j2+n-1), vd);
                for(j = i-1; j >= 0; j--)
                {
                    vc = ap::conj(a(i1+i,j1+j));
                    ap::vsub(&x(i2+j, j2), 1, &x(i2+i, j2), 1, "N", ap::vlen(j2,j2+n-1), vc);
                }
            }
            return;
        }
    }
}


/*************************************************************************
Level 2 subroutine

  -- ALGLIB routine --
     15.12.2009
     Bochkanov Sergey
*************************************************************************/
static void rmatrixrighttrsm2(int m,
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
    int i;
    int j;
    double vr;
    double vd;

    
    //
    // Special case
    //
    if( n*m==0 )
    {
        return;
    }
    
    //
    // Try to use "fast" code
    //
    if( rmatrixrighttrsmf(m, n, a, i1, j1, isupper, isunit, optype, x, i2, j2) )
    {
        return;
    }
    
    //
    // General case
    //
    if( isupper )
    {
        
        //
        // Upper triangular matrix
        //
        if( optype==0 )
        {
            
            //
            // X*A^(-1)
            //
            for(i = 0; i <= m-1; i++)
            {
                for(j = 0; j <= n-1; j++)
                {
                    if( isunit )
                    {
                        vd = 1;
                    }
                    else
                    {
                        vd = a(i1+j,j1+j);
                    }
                    x(i2+i,j2+j) = x(i2+i,j2+j)/vd;
                    if( j<n-1 )
                    {
                        vr = x(i2+i,j2+j);
                        ap::vsub(&x(i2+i, j2+j+1), 1, &a(i1+j, j1+j+1), 1, ap::vlen(j2+j+1,j2+n-1), vr);
                    }
                }
            }
            return;
        }
        if( optype==1 )
        {
            
            //
            // X*A^(-T)
            //
            for(i = 0; i <= m-1; i++)
            {
                for(j = n-1; j >= 0; j--)
                {
                    vr = 0;
                    vd = 1;
                    if( j<n-1 )
                    {
                        vr = ap::vdotproduct(&x(i2+i, j2+j+1), 1, &a(i1+j, j1+j+1), 1, ap::vlen(j2+j+1,j2+n-1));
                    }
                    if( !isunit )
                    {
                        vd = a(i1+j,j1+j);
                    }
                    x(i2+i,j2+j) = (x(i2+i,j2+j)-vr)/vd;
                }
            }
            return;
        }
    }
    else
    {
        
        //
        // Lower triangular matrix
        //
        if( optype==0 )
        {
            
            //
            // X*A^(-1)
            //
            for(i = 0; i <= m-1; i++)
            {
                for(j = n-1; j >= 0; j--)
                {
                    if( isunit )
                    {
                        vd = 1;
                    }
                    else
                    {
                        vd = a(i1+j,j1+j);
                    }
                    x(i2+i,j2+j) = x(i2+i,j2+j)/vd;
                    if( j>0 )
                    {
                        vr = x(i2+i,j2+j);
                        ap::vsub(&x(i2+i, j2), 1, &a(i1+j, j1), 1, ap::vlen(j2,j2+j-1), vr);
                    }
                }
            }
            return;
        }
        if( optype==1 )
        {
            
            //
            // X*A^(-T)
            //
            for(i = 0; i <= m-1; i++)
            {
                for(j = 0; j <= n-1; j++)
                {
                    vr = 0;
                    vd = 1;
                    if( j>0 )
                    {
                        vr = ap::vdotproduct(&x(i2+i, j2), 1, &a(i1+j, j1), 1, ap::vlen(j2,j2+j-1));
                    }
                    if( !isunit )
                    {
                        vd = a(i1+j,j1+j);
                    }
                    x(i2+i,j2+j) = (x(i2+i,j2+j)-vr)/vd;
                }
            }
            return;
        }
    }
}


/*************************************************************************
Level 2 subroutine
*************************************************************************/
static void rmatrixlefttrsm2(int m,
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
    int i;
    int j;
    double vr;
    double vd;

    
    //
    // Special case
    //
    if( n*m==0 )
    {
        return;
    }
    
    //
    // Try fast code
    //
    if( rmatrixlefttrsmf(m, n, a, i1, j1, isupper, isunit, optype, x, i2, j2) )
    {
        return;
    }
    
    //
    // General case
    //
    if( isupper )
    {
        
        //
        // Upper triangular matrix
        //
        if( optype==0 )
        {
            
            //
            // A^(-1)*X
            //
            for(i = m-1; i >= 0; i--)
            {
                for(j = i+1; j <= m-1; j++)
                {
                    vr = a(i1+i,j1+j);
                    ap::vsub(&x(i2+i, j2), 1, &x(i2+j, j2), 1, ap::vlen(j2,j2+n-1), vr);
                }
                if( !isunit )
                {
                    vd = 1/a(i1+i,j1+i);
                    ap::vmul(&x(i2+i, j2), 1, ap::vlen(j2,j2+n-1), vd);
                }
            }
            return;
        }
        if( optype==1 )
        {
            
            //
            // A^(-T)*X
            //
            for(i = 0; i <= m-1; i++)
            {
                if( isunit )
                {
                    vd = 1;
                }
                else
                {
                    vd = 1/a(i1+i,j1+i);
                }
                ap::vmul(&x(i2+i, j2), 1, ap::vlen(j2,j2+n-1), vd);
                for(j = i+1; j <= m-1; j++)
                {
                    vr = a(i1+i,j1+j);
                    ap::vsub(&x(i2+j, j2), 1, &x(i2+i, j2), 1, ap::vlen(j2,j2+n-1), vr);
                }
            }
            return;
        }
    }
    else
    {
        
        //
        // Lower triangular matrix
        //
        if( optype==0 )
        {
            
            //
            // A^(-1)*X
            //
            for(i = 0; i <= m-1; i++)
            {
                for(j = 0; j <= i-1; j++)
                {
                    vr = a(i1+i,j1+j);
                    ap::vsub(&x(i2+i, j2), 1, &x(i2+j, j2), 1, ap::vlen(j2,j2+n-1), vr);
                }
                if( isunit )
                {
                    vd = 1;
                }
                else
                {
                    vd = 1/a(i1+j,j1+j);
                }
                ap::vmul(&x(i2+i, j2), 1, ap::vlen(j2,j2+n-1), vd);
            }
            return;
        }
        if( optype==1 )
        {
            
            //
            // A^(-T)*X
            //
            for(i = m-1; i >= 0; i--)
            {
                if( isunit )
                {
                    vd = 1;
                }
                else
                {
                    vd = 1/a(i1+i,j1+i);
                }
                ap::vmul(&x(i2+i, j2), 1, ap::vlen(j2,j2+n-1), vd);
                for(j = i-1; j >= 0; j--)
                {
                    vr = a(i1+i,j1+j);
                    ap::vsub(&x(i2+j, j2), 1, &x(i2+i, j2), 1, ap::vlen(j2,j2+n-1), vr);
                }
            }
            return;
        }
    }
}


/*************************************************************************
Level 2 subroutine
*************************************************************************/
static void cmatrixsyrk2(int n,
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
    int i;
    int j;
    int j1;
    int j2;
    ap::complex v;

    
    //
    // Fast exit (nothing to be done)
    //
    if( (ap::fp_eq(alpha,0)||k==0)&&ap::fp_eq(beta,1) )
    {
        return;
    }
    
    //
    // Try to call fast SYRK
    //
    if( cmatrixsyrkf(n, k, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper) )
    {
        return;
    }
    
    //
    // SYRK
    //
    if( optypea==0 )
    {
        
        //
        // C=alpha*A*A^H+beta*C
        //
        for(i = 0; i <= n-1; i++)
        {
            if( isupper )
            {
                j1 = i;
                j2 = n-1;
            }
            else
            {
                j1 = 0;
                j2 = i;
            }
            for(j = j1; j <= j2; j++)
            {
                if( ap::fp_neq(alpha,0)&&k>0 )
                {
                    v = ap::vdotproduct(&a(ia+i, ja), 1, "N", &a(ia+j, ja), 1, "Conj", ap::vlen(ja,ja+k-1));
                }
                else
                {
                    v = 0;
                }
                if( ap::fp_eq(beta,0) )
                {
                    c(ic+i,jc+j) = alpha*v;
                }
                else
                {
                    c(ic+i,jc+j) = beta*c(ic+i,jc+j)+alpha*v;
                }
            }
        }
        return;
    }
    else
    {
        
        //
        // C=alpha*A^H*A+beta*C
        //
        for(i = 0; i <= n-1; i++)
        {
            if( isupper )
            {
                j1 = i;
                j2 = n-1;
            }
            else
            {
                j1 = 0;
                j2 = i;
            }
            if( ap::fp_eq(beta,0) )
            {
                for(j = j1; j <= j2; j++)
                {
                    c(ic+i,jc+j) = 0;
                }
            }
            else
            {
                ap::vmul(&c(ic+i, jc+j1), 1, ap::vlen(jc+j1,jc+j2), beta);
            }
        }
        for(i = 0; i <= k-1; i++)
        {
            for(j = 0; j <= n-1; j++)
            {
                if( isupper )
                {
                    j1 = j;
                    j2 = n-1;
                }
                else
                {
                    j1 = 0;
                    j2 = j;
                }
                v = alpha*ap::conj(a(ia+i,ja+j));
                ap::vadd(&c(ic+j, jc+j1), 1, &a(ia+i, ja+j1), 1, "N", ap::vlen(jc+j1,jc+j2), v);
            }
        }
        return;
    }
}


/*************************************************************************
Level 2 subrotuine
*************************************************************************/
static void rmatrixsyrk2(int n,
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
    int i;
    int j;
    int j1;
    int j2;
    double v;

    
    //
    // Fast exit (nothing to be done)
    //
    if( (ap::fp_eq(alpha,0)||k==0)&&ap::fp_eq(beta,1) )
    {
        return;
    }
    
    //
    // Try to call fast SYRK
    //
    if( rmatrixsyrkf(n, k, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper) )
    {
        return;
    }
    
    //
    // SYRK
    //
    if( optypea==0 )
    {
        
        //
        // C=alpha*A*A^H+beta*C
        //
        for(i = 0; i <= n-1; i++)
        {
            if( isupper )
            {
                j1 = i;
                j2 = n-1;
            }
            else
            {
                j1 = 0;
                j2 = i;
            }
            for(j = j1; j <= j2; j++)
            {
                if( ap::fp_neq(alpha,0)&&k>0 )
                {
                    v = ap::vdotproduct(&a(ia+i, ja), 1, &a(ia+j, ja), 1, ap::vlen(ja,ja+k-1));
                }
                else
                {
                    v = 0;
                }
                if( ap::fp_eq(beta,0) )
                {
                    c(ic+i,jc+j) = alpha*v;
                }
                else
                {
                    c(ic+i,jc+j) = beta*c(ic+i,jc+j)+alpha*v;
                }
            }
        }
        return;
    }
    else
    {
        
        //
        // C=alpha*A^H*A+beta*C
        //
        for(i = 0; i <= n-1; i++)
        {
            if( isupper )
            {
                j1 = i;
                j2 = n-1;
            }
            else
            {
                j1 = 0;
                j2 = i;
            }
            if( ap::fp_eq(beta,0) )
            {
                for(j = j1; j <= j2; j++)
                {
                    c(ic+i,jc+j) = 0;
                }
            }
            else
            {
                ap::vmul(&c(ic+i, jc+j1), 1, ap::vlen(jc+j1,jc+j2), beta);
            }
        }
        for(i = 0; i <= k-1; i++)
        {
            for(j = 0; j <= n-1; j++)
            {
                if( isupper )
                {
                    j1 = j;
                    j2 = n-1;
                }
                else
                {
                    j1 = 0;
                    j2 = j;
                }
                v = alpha*a(ia+i,ja+j);
                ap::vadd(&c(ic+j, jc+j1), 1, &a(ia+i, ja+j1), 1, ap::vlen(jc+j1,jc+j2), v);
            }
        }
        return;
    }
}


/*************************************************************************
GEMM kernel

  -- ALGLIB routine --
     16.12.2009
     Bochkanov Sergey
*************************************************************************/
static void cmatrixgemmk(int m,
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
    int i;
    int j;
    ap::complex v;

    
    //
    // Special case
    //
    if( m*n==0 )
    {
        return;
    }
    
    //
    // Try optimized code
    //
    if( cmatrixgemmf(m, n, k, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc) )
    {
        return;
    }
    
    //
    // Another special case
    //
    if( k==0 )
    {
        if( beta!=0 )
        {
            for(i = 0; i <= m-1; i++)
            {
                for(j = 0; j <= n-1; j++)
                {
                    c(ic+i,jc+j) = beta*c(ic+i,jc+j);
                }
            }
        }
        else
        {
            for(i = 0; i <= m-1; i++)
            {
                for(j = 0; j <= n-1; j++)
                {
                    c(ic+i,jc+j) = 0;
                }
            }
        }
        return;
    }
    
    //
    // General case
    //
    if( optypea==0&&optypeb!=0 )
    {
        
        //
        // A*B'
        //
        for(i = 0; i <= m-1; i++)
        {
            for(j = 0; j <= n-1; j++)
            {
                if( k==0||alpha==0 )
                {
                    v = 0;
                }
                else
                {
                    if( optypeb==1 )
                    {
                        v = ap::vdotproduct(&a(ia+i, ja), 1, "N", &b(ib+j, jb), 1, "N", ap::vlen(ja,ja+k-1));
                    }
                    else
                    {
                        v = ap::vdotproduct(&a(ia+i, ja), 1, "N", &b(ib+j, jb), 1, "Conj", ap::vlen(ja,ja+k-1));
                    }
                }
                if( beta==0 )
                {
                    c(ic+i,jc+j) = alpha*v;
                }
                else
                {
                    c(ic+i,jc+j) = beta*c(ic+i,jc+j)+alpha*v;
                }
            }
        }
        return;
    }
    if( optypea==0&&optypeb==0 )
    {
        
        //
        // A*B
        //
        for(i = 0; i <= m-1; i++)
        {
            if( beta!=0 )
            {
                ap::vmul(&c(ic+i, jc), 1, ap::vlen(jc,jc+n-1), beta);
            }
            else
            {
                for(j = 0; j <= n-1; j++)
                {
                    c(ic+i,jc+j) = 0;
                }
            }
            if( alpha!=0 )
            {
                for(j = 0; j <= k-1; j++)
                {
                    v = alpha*a(ia+i,ja+j);
                    ap::vadd(&c(ic+i, jc), 1, &b(ib+j, jb), 1, "N", ap::vlen(jc,jc+n-1), v);
                }
            }
        }
        return;
    }
    if( optypea!=0&&optypeb!=0 )
    {
        
        //
        // A'*B'
        //
        for(i = 0; i <= m-1; i++)
        {
            for(j = 0; j <= n-1; j++)
            {
                if( alpha==0 )
                {
                    v = 0;
                }
                else
                {
                    if( optypea==1 )
                    {
                        if( optypeb==1 )
                        {
                            v = ap::vdotproduct(&a(ia, ja+i), a.getstride(), "N", &b(ib+j, jb), 1, "N", ap::vlen(ia,ia+k-1));
                        }
                        else
                        {
                            v = ap::vdotproduct(&a(ia, ja+i), a.getstride(), "N", &b(ib+j, jb), 1, "Conj", ap::vlen(ia,ia+k-1));
                        }
                    }
                    else
                    {
                        if( optypeb==1 )
                        {
                            v = ap::vdotproduct(&a(ia, ja+i), a.getstride(), "Conj", &b(ib+j, jb), 1, "N", ap::vlen(ia,ia+k-1));
                        }
                        else
                        {
                            v = ap::vdotproduct(&a(ia, ja+i), a.getstride(), "Conj", &b(ib+j, jb), 1, "Conj", ap::vlen(ia,ia+k-1));
                        }
                    }
                }
                if( beta==0 )
                {
                    c(ic+i,jc+j) = alpha*v;
                }
                else
                {
                    c(ic+i,jc+j) = beta*c(ic+i,jc+j)+alpha*v;
                }
            }
        }
        return;
    }
    if( optypea!=0&&optypeb==0 )
    {
        
        //
        // A'*B
        //
        if( beta==0 )
        {
            for(i = 0; i <= m-1; i++)
            {
                for(j = 0; j <= n-1; j++)
                {
                    c(ic+i,jc+j) = 0;
                }
            }
        }
        else
        {
            for(i = 0; i <= m-1; i++)
            {
                ap::vmul(&c(ic+i, jc), 1, ap::vlen(jc,jc+n-1), beta);
            }
        }
        if( alpha!=0 )
        {
            for(j = 0; j <= k-1; j++)
            {
                for(i = 0; i <= m-1; i++)
                {
                    if( optypea==1 )
                    {
                        v = alpha*a(ia+j,ja+i);
                    }
                    else
                    {
                        v = alpha*ap::conj(a(ia+j,ja+i));
                    }
                    ap::vadd(&c(ic+i, jc), 1, &b(ib+j, jb), 1, "N", ap::vlen(jc,jc+n-1), v);
                }
            }
        }
        return;
    }
}


/*************************************************************************
GEMM kernel

  -- ALGLIB routine --
     16.12.2009
     Bochkanov Sergey
*************************************************************************/
static void rmatrixgemmk(int m,
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
    int i;
    int j;
    double v;

    
    //
    // if matrix size is zero
    //
    if( m*n==0 )
    {
        return;
    }
    
    //
    // Try optimized code
    //
    if( rmatrixgemmf(m, n, k, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc) )
    {
        return;
    }
    
    //
    // if K=0, then C=Beta*C
    //
    if( k==0 )
    {
        if( ap::fp_neq(beta,1) )
        {
            if( ap::fp_neq(beta,0) )
            {
                for(i = 0; i <= m-1; i++)
                {
                    for(j = 0; j <= n-1; j++)
                    {
                        c(ic+i,jc+j) = beta*c(ic+i,jc+j);
                    }
                }
            }
            else
            {
                for(i = 0; i <= m-1; i++)
                {
                    for(j = 0; j <= n-1; j++)
                    {
                        c(ic+i,jc+j) = 0;
                    }
                }
            }
        }
        return;
    }
    
    //
    // General case
    //
    if( optypea==0&&optypeb!=0 )
    {
        
        //
        // A*B'
        //
        for(i = 0; i <= m-1; i++)
        {
            for(j = 0; j <= n-1; j++)
            {
                if( k==0||ap::fp_eq(alpha,0) )
                {
                    v = 0;
                }
                else
                {
                    v = ap::vdotproduct(&a(ia+i, ja), 1, &b(ib+j, jb), 1, ap::vlen(ja,ja+k-1));
                }
                if( ap::fp_eq(beta,0) )
                {
                    c(ic+i,jc+j) = alpha*v;
                }
                else
                {
                    c(ic+i,jc+j) = beta*c(ic+i,jc+j)+alpha*v;
                }
            }
        }
        return;
    }
    if( optypea==0&&optypeb==0 )
    {
        
        //
        // A*B
        //
        for(i = 0; i <= m-1; i++)
        {
            if( ap::fp_neq(beta,0) )
            {
                ap::vmul(&c(ic+i, jc), 1, ap::vlen(jc,jc+n-1), beta);
            }
            else
            {
                for(j = 0; j <= n-1; j++)
                {
                    c(ic+i,jc+j) = 0;
                }
            }
            if( ap::fp_neq(alpha,0) )
            {
                for(j = 0; j <= k-1; j++)
                {
                    v = alpha*a(ia+i,ja+j);
                    ap::vadd(&c(ic+i, jc), 1, &b(ib+j, jb), 1, ap::vlen(jc,jc+n-1), v);
                }
            }
        }
        return;
    }
    if( optypea!=0&&optypeb!=0 )
    {
        
        //
        // A'*B'
        //
        for(i = 0; i <= m-1; i++)
        {
            for(j = 0; j <= n-1; j++)
            {
                if( ap::fp_eq(alpha,0) )
                {
                    v = 0;
                }
                else
                {
                    v = ap::vdotproduct(&a(ia, ja+i), a.getstride(), &b(ib+j, jb), 1, ap::vlen(ia,ia+k-1));
                }
                if( ap::fp_eq(beta,0) )
                {
                    c(ic+i,jc+j) = alpha*v;
                }
                else
                {
                    c(ic+i,jc+j) = beta*c(ic+i,jc+j)+alpha*v;
                }
            }
        }
        return;
    }
    if( optypea!=0&&optypeb==0 )
    {
        
        //
        // A'*B
        //
        if( ap::fp_eq(beta,0) )
        {
            for(i = 0; i <= m-1; i++)
            {
                for(j = 0; j <= n-1; j++)
                {
                    c(ic+i,jc+j) = 0;
                }
            }
        }
        else
        {
            for(i = 0; i <= m-1; i++)
            {
                ap::vmul(&c(ic+i, jc), 1, ap::vlen(jc,jc+n-1), beta);
            }
        }
        if( ap::fp_neq(alpha,0) )
        {
            for(j = 0; j <= k-1; j++)
            {
                for(i = 0; i <= m-1; i++)
                {
                    v = alpha*a(ia+j,ja+i);
                    ap::vadd(&c(ic+i, jc), 1, &b(ib+j, jb), 1, ap::vlen(jc,jc+n-1), v);
                }
            }
        }
        return;
    }
}





/*************************************************************************
Copyright (c) 2005-2010 Sergey Bochkanov.

Additional copyrights:
    1992-2007 The University of Tennessee (as indicated in subroutines
    comments).

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
#include "ortfac.h"

static void rmatrixqrbasecase(ap::real_2d_array& a,
     int m,
     int n,
     ap::real_1d_array& work,
     ap::real_1d_array& t,
     ap::real_1d_array& tau);
static void rmatrixlqbasecase(ap::real_2d_array& a,
     int m,
     int n,
     ap::real_1d_array& work,
     ap::real_1d_array& t,
     ap::real_1d_array& tau);
static void cmatrixqrbasecase(ap::complex_2d_array& a,
     int m,
     int n,
     ap::complex_1d_array& work,
     ap::complex_1d_array& t,
     ap::complex_1d_array& tau);
static void cmatrixlqbasecase(ap::complex_2d_array& a,
     int m,
     int n,
     ap::complex_1d_array& work,
     ap::complex_1d_array& t,
     ap::complex_1d_array& tau);
static void rmatrixblockreflector(ap::real_2d_array& a,
     ap::real_1d_array& tau,
     bool columnwisea,
     int lengtha,
     int blocksize,
     ap::real_2d_array& t,
     ap::real_1d_array& work);
static void cmatrixblockreflector(ap::complex_2d_array& a,
     ap::complex_1d_array& tau,
     bool columnwisea,
     int lengtha,
     int blocksize,
     ap::complex_2d_array& t,
     ap::complex_1d_array& work);

/*************************************************************************
QR decomposition of a rectangular matrix of size MxN

Input parameters:
    A   -   matrix A whose indexes range within [0..M-1, 0..N-1].
    M   -   number of rows in matrix A.
    N   -   number of columns in matrix A.

Output parameters:
    A   -   matrices Q and R in compact form (see below).
    Tau -   array of scalar factors which are used to form
            matrix Q. Array whose index ranges within [0.. Min(M-1,N-1)].

Matrix A is represented as A = QR, where Q is an orthogonal matrix of size
MxM, R - upper triangular (or upper trapezoid) matrix of size M x N.

The elements of matrix R are located on and above the main diagonal of
matrix A. The elements which are located in Tau array and below the main
diagonal of matrix A are used to form matrix Q as follows:

Matrix Q is represented as a product of elementary reflections

Q = H(0)*H(2)*...*H(k-1),

where k = min(m,n), and each H(i) is in the form

H(i) = 1 - tau * v * (v^T)

where tau is a scalar stored in Tau[I]; v - real vector,
so that v(0:i-1) = 0, v(i) = 1, v(i+1:m-1) stored in A(i+1:m-1,i).

  -- ALGLIB routine --
     17.02.2010
     Bochkanov Sergey
*************************************************************************/
void rmatrixqr(ap::real_2d_array& a, int m, int n, ap::real_1d_array& tau)
{
    ap::real_1d_array work;
    ap::real_1d_array t;
    ap::real_1d_array taubuf;
    int minmn;
    ap::real_2d_array tmpa;
    ap::real_2d_array tmpt;
    ap::real_2d_array tmpr;
    int blockstart;
    int blocksize;
    int rowscount;
    int i;
//    int j;
//    int k;
//    double v;

    if( m<=0||n<=0 )
    {
        return;
    }
    minmn = ap::minint(m, n);
    work.setlength(ap::maxint(m, n)+1);
    t.setlength(ap::maxint(m, n)+1);
    tau.setlength(minmn);
    taubuf.setlength(minmn);
    tmpa.setlength(m, ablasblocksize(a));
    tmpt.setlength(ablasblocksize(a), 2*ablasblocksize(a));
    tmpr.setlength(2*ablasblocksize(a), n);
    
    //
    // Blocked code
    //
    blockstart = 0;
    while(blockstart!=minmn)
    {
        
        //
        // Determine block size
        //
        blocksize = minmn-blockstart;
        if( blocksize>ablasblocksize(a) )
        {
            blocksize = ablasblocksize(a);
        }
        rowscount = m-blockstart;
        
        //
        // QR decomposition of submatrix.
        // Matrix is copied to temporary storage to solve
        // some TLB issues arising from non-contiguous memory
        // access pattern.
        //
        rmatrixcopy(rowscount, blocksize, a, blockstart, blockstart, tmpa, 0, 0);
        rmatrixqrbasecase(tmpa, rowscount, blocksize, work, t, taubuf);
        rmatrixcopy(rowscount, blocksize, tmpa, 0, 0, a, blockstart, blockstart);
        ap::vmove(&tau(blockstart), 1, &taubuf(0), 1, ap::vlen(blockstart,blockstart+blocksize-1));
        
        //
        // Update the rest, choose between:
        // a) Level 2 algorithm (when the rest of the matrix is small enough)
        // b) blocked algorithm, see algorithm 5 from  'A storage efficient WY
        //    representation for products of Householder transformations',
        //    by R. Schreiber and C. Van Loan.
        //
        if( blockstart+blocksize<=n-1 )
        {
            if( n-blockstart-blocksize>=2*ablasblocksize(a)||rowscount>=4*ablasblocksize(a) )
            {
                
                //
                // Prepare block reflector
                //
                rmatrixblockreflector(tmpa, taubuf, true, rowscount, blocksize, tmpt, work);
                
                //
                // Multiply the rest of A by Q'.
                //
                // Q  = E + Y*T*Y'  = E + TmpA*TmpT*TmpA'
                // Q' = E + Y*T'*Y' = E + TmpA*TmpT'*TmpA'
                //
                rmatrixgemm(blocksize, n-blockstart-blocksize, rowscount, 1.0, tmpa, 0, 0, 1, a, blockstart, blockstart+blocksize, 0, 0.0, tmpr, 0, 0);
                rmatrixgemm(blocksize, n-blockstart-blocksize, blocksize, 1.0, tmpt, 0, 0, 1, tmpr, 0, 0, 0, 0.0, tmpr, blocksize, 0);
                rmatrixgemm(rowscount, n-blockstart-blocksize, blocksize, 1.0, tmpa, 0, 0, 0, tmpr, blocksize, 0, 0, 1.0, a, blockstart, blockstart+blocksize);
            }
            else
            {
                
                //
                // Level 2 algorithm
                //
                for(i = 0; i <= blocksize-1; i++)
                {
                    ap::vmove(&t(1), 1, &tmpa(i, i), tmpa.getstride(), ap::vlen(1,rowscount-i));
                    t(1) = 1;
                    applyreflectionfromtheleft(a, taubuf(i), t, blockstart+i, m-1, blockstart+blocksize, n-1, work);
                }
            }
        }
        
        //
        // Advance
        //
        blockstart = blockstart+blocksize;
    }
}


/*************************************************************************
LQ decomposition of a rectangular matrix of size MxN

Input parameters:
    A   -   matrix A whose indexes range within [0..M-1, 0..N-1].
    M   -   number of rows in matrix A.
    N   -   number of columns in matrix A.

Output parameters:
    A   -   matrices L and Q in compact form (see below)
    Tau -   array of scalar factors which are used to form
            matrix Q. Array whose index ranges within [0..Min(M,N)-1].

Matrix A is represented as A = LQ, where Q is an orthogonal matrix of size
MxM, L - lower triangular (or lower trapezoid) matrix of size M x N.

The elements of matrix L are located on and below  the  main  diagonal  of
matrix A. The elements which are located in Tau array and above  the  main
diagonal of matrix A are used to form matrix Q as follows:

Matrix Q is represented as a product of elementary reflections

Q = H(k-1)*H(k-2)*...*H(1)*H(0),

where k = min(m,n), and each H(i) is of the form

H(i) = 1 - tau * v * (v^T)

where tau is a scalar stored in Tau[I]; v - real vector, so that v(0:i-1)=0,
v(i) = 1, v(i+1:n-1) stored in A(i,i+1:n-1).

  -- ALGLIB routine --
     17.02.2010
     Bochkanov Sergey
*************************************************************************/
void rmatrixlq(ap::real_2d_array& a, int m, int n, ap::real_1d_array& tau)
{
    ap::real_1d_array work;
    ap::real_1d_array t;
    ap::real_1d_array taubuf;
    int minmn;
    ap::real_2d_array tmpa;
    ap::real_2d_array tmpt;
    ap::real_2d_array tmpr;
    int blockstart;
    int blocksize;
    int columnscount;
    int i;
//    int j;
//    int k;
//    double v;

    if( m<=0||n<=0 )
    {
        return;
    }
    minmn = ap::minint(m, n);
    work.setlength(ap::maxint(m, n)+1);
    t.setlength(ap::maxint(m, n)+1);
    tau.setlength(minmn);
    taubuf.setlength(minmn);
    tmpa.setlength(ablasblocksize(a), n);
    tmpt.setlength(ablasblocksize(a), 2*ablasblocksize(a));
    tmpr.setlength(m, 2*ablasblocksize(a));
    
    //
    // Blocked code
    //
    blockstart = 0;
    while(blockstart!=minmn)
    {
        
        //
        // Determine block size
        //
        blocksize = minmn-blockstart;
        if( blocksize>ablasblocksize(a) )
        {
            blocksize = ablasblocksize(a);
        }
        columnscount = n-blockstart;
        
        //
        // LQ decomposition of submatrix.
        // Matrix is copied to temporary storage to solve
        // some TLB issues arising from non-contiguous memory
        // access pattern.
        //
        rmatrixcopy(blocksize, columnscount, a, blockstart, blockstart, tmpa, 0, 0);
        rmatrixlqbasecase(tmpa, blocksize, columnscount, work, t, taubuf);
        rmatrixcopy(blocksize, columnscount, tmpa, 0, 0, a, blockstart, blockstart);
        ap::vmove(&tau(blockstart), 1, &taubuf(0), 1, ap::vlen(blockstart,blockstart+blocksize-1));
        
        //
        // Update the rest, choose between:
        // a) Level 2 algorithm (when the rest of the matrix is small enough)
        // b) blocked algorithm, see algorithm 5 from  'A storage efficient WY
        //    representation for products of Householder transformations',
        //    by R. Schreiber and C. Van Loan.
        //
        if( blockstart+blocksize<=m-1 )
        {
            if( m-blockstart-blocksize>=2*ablasblocksize(a) )
            {
                
                //
                // Prepare block reflector
                //
                rmatrixblockreflector(tmpa, taubuf, false, columnscount, blocksize, tmpt, work);
                
                //
                // Multiply the rest of A by Q.
                //
                // Q  = E + Y*T*Y'  = E + TmpA'*TmpT*TmpA
                //
                rmatrixgemm(m-blockstart-blocksize, blocksize, columnscount, 1.0, a, blockstart+blocksize, blockstart, 0, tmpa, 0, 0, 1, 0.0, tmpr, 0, 0);
                rmatrixgemm(m-blockstart-blocksize, blocksize, blocksize, 1.0, tmpr, 0, 0, 0, tmpt, 0, 0, 0, 0.0, tmpr, 0, blocksize);
                rmatrixgemm(m-blockstart-blocksize, columnscount, blocksize, 1.0, tmpr, 0, blocksize, 0, tmpa, 0, 0, 0, 1.0, a, blockstart+blocksize, blockstart);
            }
            else
            {
                
                //
                // Level 2 algorithm
                //
                for(i = 0; i <= blocksize-1; i++)
                {
                    ap::vmove(&t(1), 1, &tmpa(i, i), 1, ap::vlen(1,columnscount-i));
                    t(1) = 1;
                    applyreflectionfromtheright(a, taubuf(i), t, blockstart+blocksize, m-1, blockstart+i, n-1, work);
                }
            }
        }
        
        //
        // Advance
        //
        blockstart = blockstart+blocksize;
    }
}


/*************************************************************************
QR decomposition of a rectangular complex matrix of size MxN

Input parameters:
    A   -   matrix A whose indexes range within [0..M-1, 0..N-1]
    M   -   number of rows in matrix A.
    N   -   number of columns in matrix A.

Output parameters:
    A   -   matrices Q and R in compact form
    Tau -   array of scalar factors which are used to form matrix Q. Array
            whose indexes range within [0.. Min(M,N)-1]

Matrix A is represented as A = QR, where Q is an orthogonal matrix of size
MxM, R - upper triangular (or upper trapezoid) matrix of size MxN.

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     September 30, 1994
*************************************************************************/
void cmatrixqr(ap::complex_2d_array& a,
     int m,
     int n,
     ap::complex_1d_array& tau)
{
    ap::complex_1d_array work;
    ap::complex_1d_array t;
    ap::complex_1d_array taubuf;
    int minmn;
    ap::complex_2d_array tmpa;
    ap::complex_2d_array tmpt;
    ap::complex_2d_array tmpr;
    int blockstart;
    int blocksize;
    int rowscount;
    int i;
//    int j;
//    int k;
    ap::complex v;

    if( m<=0||n<=0 )
    {
        return;
    }
    minmn = ap::minint(m, n);
    work.setlength(ap::maxint(m, n)+1);
    t.setlength(ap::maxint(m, n)+1);
    tau.setlength(minmn);
    taubuf.setlength(minmn);
    tmpa.setlength(m, ablascomplexblocksize(a));
    tmpt.setlength(ablascomplexblocksize(a), ablascomplexblocksize(a));
    tmpr.setlength(2*ablascomplexblocksize(a), n);
    
    //
    // Blocked code
    //
    blockstart = 0;
    while(blockstart!=minmn)
    {
        
        //
        // Determine block size
        //
        blocksize = minmn-blockstart;
        if( blocksize>ablascomplexblocksize(a) )
        {
            blocksize = ablascomplexblocksize(a);
        }
        rowscount = m-blockstart;
        
        //
        // QR decomposition of submatrix.
        // Matrix is copied to temporary storage to solve
        // some TLB issues arising from non-contiguous memory
        // access pattern.
        //
        cmatrixcopy(rowscount, blocksize, a, blockstart, blockstart, tmpa, 0, 0);
        cmatrixqrbasecase(tmpa, rowscount, blocksize, work, t, taubuf);
        cmatrixcopy(rowscount, blocksize, tmpa, 0, 0, a, blockstart, blockstart);
        ap::vmove(&tau(blockstart), 1, &taubuf(0), 1, "N", ap::vlen(blockstart,blockstart+blocksize-1));
        
        //
        // Update the rest, choose between:
        // a) Level 2 algorithm (when the rest of the matrix is small enough)
        // b) blocked algorithm, see algorithm 5 from  'A storage efficient WY
        //    representation for products of Householder transformations',
        //    by R. Schreiber and C. Van Loan.
        //
        if( blockstart+blocksize<=n-1 )
        {
            if( n-blockstart-blocksize>=2*ablascomplexblocksize(a) )
            {
                
                //
                // Prepare block reflector
                //
                cmatrixblockreflector(tmpa, taubuf, true, rowscount, blocksize, tmpt, work);
                
                //
                // Multiply the rest of A by Q'.
                //
                // Q  = E + Y*T*Y'  = E + TmpA*TmpT*TmpA'
                // Q' = E + Y*T'*Y' = E + TmpA*TmpT'*TmpA'
                //
                cmatrixgemm(blocksize, n-blockstart-blocksize, rowscount, 1.0, tmpa, 0, 0, 2, a, blockstart, blockstart+blocksize, 0, 0.0, tmpr, 0, 0);
                cmatrixgemm(blocksize, n-blockstart-blocksize, blocksize, 1.0, tmpt, 0, 0, 2, tmpr, 0, 0, 0, 0.0, tmpr, blocksize, 0);
                cmatrixgemm(rowscount, n-blockstart-blocksize, blocksize, 1.0, tmpa, 0, 0, 0, tmpr, blocksize, 0, 0, 1.0, a, blockstart, blockstart+blocksize);
            }
            else
            {
                
                //
                // Level 2 algorithm
                //
                for(i = 0; i <= blocksize-1; i++)
                {
                    ap::vmove(&t(1), 1, &tmpa(i, i), tmpa.getstride(), "N", ap::vlen(1,rowscount-i));
                    t(1) = 1;
                    complexapplyreflectionfromtheleft(a, ap::conj(taubuf(i)), t, blockstart+i, m-1, blockstart+blocksize, n-1, work);
                }
            }
        }
        
        //
        // Advance
        //
        blockstart = blockstart+blocksize;
    }
}


/*************************************************************************
LQ decomposition of a rectangular complex matrix of size MxN

Input parameters:
    A   -   matrix A whose indexes range within [0..M-1, 0..N-1]
    M   -   number of rows in matrix A.
    N   -   number of columns in matrix A.

Output parameters:
    A   -   matrices Q and L in compact form
    Tau -   array of scalar factors which are used to form matrix Q. Array
            whose indexes range within [0.. Min(M,N)-1]

Matrix A is represented as A = LQ, where Q is an orthogonal matrix of size
MxM, L - lower triangular (or lower trapezoid) matrix of size MxN.

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     September 30, 1994
*************************************************************************/
void cmatrixlq(ap::complex_2d_array& a,
     int m,
     int n,
     ap::complex_1d_array& tau)
{
    ap::complex_1d_array work;
    ap::complex_1d_array t;
    ap::complex_1d_array taubuf;
    int minmn;
    ap::complex_2d_array tmpa;
    ap::complex_2d_array tmpt;
    ap::complex_2d_array tmpr;
    int blockstart;
    int blocksize;
    int columnscount;
    int i;
//    int j;
//    int k;
    ap::complex v;

    if( m<=0||n<=0 )
    {
        return;
    }
    minmn = ap::minint(m, n);
    work.setlength(ap::maxint(m, n)+1);
    t.setlength(ap::maxint(m, n)+1);
    tau.setlength(minmn);
    taubuf.setlength(minmn);
    tmpa.setlength(ablascomplexblocksize(a), n);
    tmpt.setlength(ablascomplexblocksize(a), ablascomplexblocksize(a));
    tmpr.setlength(m, 2*ablascomplexblocksize(a));
    
    //
    // Blocked code
    //
    blockstart = 0;
    while(blockstart!=minmn)
    {
        
        //
        // Determine block size
        //
        blocksize = minmn-blockstart;
        if( blocksize>ablascomplexblocksize(a) )
        {
            blocksize = ablascomplexblocksize(a);
        }
        columnscount = n-blockstart;
        
        //
        // LQ decomposition of submatrix.
        // Matrix is copied to temporary storage to solve
        // some TLB issues arising from non-contiguous memory
        // access pattern.
        //
        cmatrixcopy(blocksize, columnscount, a, blockstart, blockstart, tmpa, 0, 0);
        cmatrixlqbasecase(tmpa, blocksize, columnscount, work, t, taubuf);
        cmatrixcopy(blocksize, columnscount, tmpa, 0, 0, a, blockstart, blockstart);
        ap::vmove(&tau(blockstart), 1, &taubuf(0), 1, "N", ap::vlen(blockstart,blockstart+blocksize-1));
        
        //
        // Update the rest, choose between:
        // a) Level 2 algorithm (when the rest of the matrix is small enough)
        // b) blocked algorithm, see algorithm 5 from  'A storage efficient WY
        //    representation for products of Householder transformations',
        //    by R. Schreiber and C. Van Loan.
        //
        if( blockstart+blocksize<=m-1 )
        {
            if( m-blockstart-blocksize>=2*ablascomplexblocksize(a) )
            {
                
                //
                // Prepare block reflector
                //
                cmatrixblockreflector(tmpa, taubuf, false, columnscount, blocksize, tmpt, work);
                
                //
                // Multiply the rest of A by Q.
                //
                // Q  = E + Y*T*Y'  = E + TmpA'*TmpT*TmpA
                //
                cmatrixgemm(m-blockstart-blocksize, blocksize, columnscount, 1.0, a, blockstart+blocksize, blockstart, 0, tmpa, 0, 0, 2, 0.0, tmpr, 0, 0);
                cmatrixgemm(m-blockstart-blocksize, blocksize, blocksize, 1.0, tmpr, 0, 0, 0, tmpt, 0, 0, 0, 0.0, tmpr, 0, blocksize);
                cmatrixgemm(m-blockstart-blocksize, columnscount, blocksize, 1.0, tmpr, 0, blocksize, 0, tmpa, 0, 0, 0, 1.0, a, blockstart+blocksize, blockstart);
            }
            else
            {
                
                //
                // Level 2 algorithm
                //
                for(i = 0; i <= blocksize-1; i++)
                {
                    ap::vmove(&t(1), 1, &tmpa(i, i), 1, "Conj", ap::vlen(1,columnscount-i));
                    t(1) = 1;
                    complexapplyreflectionfromtheright(a, taubuf(i), t, blockstart+blocksize, m-1, blockstart+i, n-1, work);
                }
            }
        }
        
        //
        // Advance
        //
        blockstart = blockstart+blocksize;
    }
}


/*************************************************************************
Partial unpacking of matrix Q from the QR decomposition of a matrix A

Input parameters:
    A       -   matrices Q and R in compact form.
                Output of RMatrixQR subroutine.
    M       -   number of rows in given matrix A. M>=0.
    N       -   number of columns in given matrix A. N>=0.
    Tau     -   scalar factors which are used to form Q.
                Output of the RMatrixQR subroutine.
    QColumns -  required number of columns of matrix Q. M>=QColumns>=0.

Output parameters:
    Q       -   first QColumns columns of matrix Q.
                Array whose indexes range within [0..M-1, 0..QColumns-1].
                If QColumns=0, the array remains unchanged.

  -- ALGLIB routine --
     17.02.2010
     Bochkanov Sergey
*************************************************************************/
void rmatrixqrunpackq(const ap::real_2d_array& a,
     int m,
     int n,
     const ap::real_1d_array& tau,
     int qcolumns,
     ap::real_2d_array& q)
{
    ap::real_1d_array work;
    ap::real_1d_array t;
    ap::real_1d_array taubuf;
    int minmn;
    int refcnt;
    ap::real_2d_array tmpa;
    ap::real_2d_array tmpt;
    ap::real_2d_array tmpr;
    int blockstart;
    int blocksize;
    int rowscount;
    int i;
    int j;
//    int k;
//    double v;

    ap::ap_error::make_assertion(qcolumns<=m, "UnpackQFromQR: QColumns>M!");
    if( m<=0||n<=0||qcolumns<=0 )
    {
        return;
    }
    
    //
    // init
    //
    minmn = ap::minint(m, n);
    refcnt = ap::minint(minmn, qcolumns);
    q.setlength(m, qcolumns);
    for(i = 0; i <= m-1; i++)
    {
        for(j = 0; j <= qcolumns-1; j++)
        {
            if( i==j )
            {
                q(i,j) = 1;
            }
            else
            {
                q(i,j) = 0;
            }
        }
    }
    work.setlength(ap::maxint(m, qcolumns)+1);
    t.setlength(ap::maxint(m, qcolumns)+1);
    taubuf.setlength(minmn);
    tmpa.setlength(m, ablasblocksize(a));
    tmpt.setlength(ablasblocksize(a), 2*ablasblocksize(a));
    tmpr.setlength(2*ablasblocksize(a), qcolumns);
    
    //
    // Blocked code
    //
    blockstart = ablasblocksize(a)*(refcnt/ablasblocksize(a));
    blocksize = refcnt-blockstart;
    while(blockstart>=0)
    {
        rowscount = m-blockstart;
        
        //
        // Copy current block
        //
        rmatrixcopy(rowscount, blocksize, a, blockstart, blockstart, tmpa, 0, 0);
        ap::vmove(&taubuf(0), 1, &tau(blockstart), 1, ap::vlen(0,blocksize-1));
        
        //
        // Update, choose between:
        // a) Level 2 algorithm (when the rest of the matrix is small enough)
        // b) blocked algorithm, see algorithm 5 from  'A storage efficient WY
        //    representation for products of Householder transformations',
        //    by R. Schreiber and C. Van Loan.
        //
        if( qcolumns>=2*ablasblocksize(a) )
        {
            
            //
            // Prepare block reflector
            //
            rmatrixblockreflector(tmpa, taubuf, true, rowscount, blocksize, tmpt, work);
            
            //
            // Multiply matrix by Q.
            //
            // Q  = E + Y*T*Y'  = E + TmpA*TmpT*TmpA'
            //
            rmatrixgemm(blocksize, qcolumns, rowscount, 1.0, tmpa, 0, 0, 1, q, blockstart, 0, 0, 0.0, tmpr, 0, 0);
            rmatrixgemm(blocksize, qcolumns, blocksize, 1.0, tmpt, 0, 0, 0, tmpr, 0, 0, 0, 0.0, tmpr, blocksize, 0);
            rmatrixgemm(rowscount, qcolumns, blocksize, 1.0, tmpa, 0, 0, 0, tmpr, blocksize, 0, 0, 1.0, q, blockstart, 0);
        }
        else
        {
            
            //
            // Level 2 algorithm
            //
            for(i = blocksize-1; i >= 0; i--)
            {
                ap::vmove(&t(1), 1, &tmpa(i, i), tmpa.getstride(), ap::vlen(1,rowscount-i));
                t(1) = 1;
                applyreflectionfromtheleft(q, taubuf(i), t, blockstart+i, m-1, 0, qcolumns-1, work);
            }
        }
        
        //
        // Advance
        //
        blockstart = blockstart-ablasblocksize(a);
        blocksize = ablasblocksize(a);
    }
}


/*************************************************************************
Unpacking of matrix R from the QR decomposition of a matrix A

Input parameters:
    A       -   matrices Q and R in compact form.
                Output of RMatrixQR subroutine.
    M       -   number of rows in given matrix A. M>=0.
    N       -   number of columns in given matrix A. N>=0.

Output parameters:
    R       -   matrix R, array[0..M-1, 0..N-1].

  -- ALGLIB routine --
     17.02.2010
     Bochkanov Sergey
*************************************************************************/
void rmatrixqrunpackr(const ap::real_2d_array& a,
     int m,
     int n,
     ap::real_2d_array& r)
{
    int i;
    int k;

    if( m<=0||n<=0 )
    {
        return;
    }
    k = ap::minint(m, n);
    r.setlength(m, n);
    for(i = 0; i <= n-1; i++)
    {
        r(0,i) = 0;
    }
    for(i = 1; i <= m-1; i++)
    {
        ap::vmove(&r(i, 0), 1, &r(0, 0), 1, ap::vlen(0,n-1));
    }
    for(i = 0; i <= k-1; i++)
    {
        ap::vmove(&r(i, i), 1, &a(i, i), 1, ap::vlen(i,n-1));
    }
}


/*************************************************************************
Partial unpacking of matrix Q from the LQ decomposition of a matrix A

Input parameters:
    A       -   matrices L and Q in compact form.
                Output of RMatrixLQ subroutine.
    M       -   number of rows in given matrix A. M>=0.
    N       -   number of columns in given matrix A. N>=0.
    Tau     -   scalar factors which are used to form Q.
                Output of the RMatrixLQ subroutine.
    QRows   -   required number of rows in matrix Q. N>=QRows>=0.

Output parameters:
    Q       -   first QRows rows of matrix Q. Array whose indexes range
                within [0..QRows-1, 0..N-1]. If QRows=0, the array remains
                unchanged.

  -- ALGLIB routine --
     17.02.2010
     Bochkanov Sergey
*************************************************************************/
void rmatrixlqunpackq(const ap::real_2d_array& a,
     int m,
     int n,
     const ap::real_1d_array& tau,
     int qrows,
     ap::real_2d_array& q)
{
    ap::real_1d_array work;
    ap::real_1d_array t;
    ap::real_1d_array taubuf;
    int minmn;
    int refcnt;
    ap::real_2d_array tmpa;
    ap::real_2d_array tmpt;
    ap::real_2d_array tmpr;
    int blockstart;
    int blocksize;
    int columnscount;
    int i;
    int j;
//    int k;
//    double v;

    ap::ap_error::make_assertion(qrows<=n, "RMatrixLQUnpackQ: QRows>N!");
    if( m<=0||n<=0||qrows<=0 )
    {
        return;
    }
    
    //
    // init
    //
    minmn = ap::minint(m, n);
    refcnt = ap::minint(minmn, qrows);
    work.setlength(ap::maxint(m, n)+1);
    t.setlength(ap::maxint(m, n)+1);
    taubuf.setlength(minmn);
    tmpa.setlength(ablasblocksize(a), n);
    tmpt.setlength(ablasblocksize(a), 2*ablasblocksize(a));
    tmpr.setlength(qrows, 2*ablasblocksize(a));
    q.setlength(qrows, n);
    for(i = 0; i <= qrows-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            if( i==j )
            {
                q(i,j) = 1;
            }
            else
            {
                q(i,j) = 0;
            }
        }
    }
    
    //
    // Blocked code
    //
    blockstart = ablasblocksize(a)*(refcnt/ablasblocksize(a));
    blocksize = refcnt-blockstart;
    while(blockstart>=0)
    {
        columnscount = n-blockstart;
        
        //
        // Copy submatrix
        //
        rmatrixcopy(blocksize, columnscount, a, blockstart, blockstart, tmpa, 0, 0);
        ap::vmove(&taubuf(0), 1, &tau(blockstart), 1, ap::vlen(0,blocksize-1));
        
        //
        // Update matrix, choose between:
        // a) Level 2 algorithm (when the rest of the matrix is small enough)
        // b) blocked algorithm, see algorithm 5 from  'A storage efficient WY
        //    representation for products of Householder transformations',
        //    by R. Schreiber and C. Van Loan.
        //
        if( qrows>=2*ablasblocksize(a) )
        {
            
            //
            // Prepare block reflector
            //
            rmatrixblockreflector(tmpa, taubuf, false, columnscount, blocksize, tmpt, work);
            
            //
            // Multiply the rest of A by Q'.
            //
            // Q'  = E + Y*T'*Y'  = E + TmpA'*TmpT'*TmpA
            //
            rmatrixgemm(qrows, blocksize, columnscount, 1.0, q, 0, blockstart, 0, tmpa, 0, 0, 1, 0.0, tmpr, 0, 0);
            rmatrixgemm(qrows, blocksize, blocksize, 1.0, tmpr, 0, 0, 0, tmpt, 0, 0, 1, 0.0, tmpr, 0, blocksize);
            rmatrixgemm(qrows, columnscount, blocksize, 1.0, tmpr, 0, blocksize, 0, tmpa, 0, 0, 0, 1.0, q, 0, blockstart);
        }
        else
        {
            
            //
            // Level 2 algorithm
            //
            for(i = blocksize-1; i >= 0; i--)
            {
                ap::vmove(&t(1), 1, &tmpa(i, i), 1, ap::vlen(1,columnscount-i));
                t(1) = 1;
                applyreflectionfromtheright(q, taubuf(i), t, 0, qrows-1, blockstart+i, n-1, work);
            }
        }
        
        //
        // Advance
        //
        blockstart = blockstart-ablasblocksize(a);
        blocksize = ablasblocksize(a);
    }
}


/*************************************************************************
Unpacking of matrix L from the LQ decomposition of a matrix A

Input parameters:
    A       -   matrices Q and L in compact form.
                Output of RMatrixLQ subroutine.
    M       -   number of rows in given matrix A. M>=0.
    N       -   number of columns in given matrix A. N>=0.

Output parameters:
    L       -   matrix L, array[0..M-1, 0..N-1].

  -- ALGLIB routine --
     17.02.2010
     Bochkanov Sergey
*************************************************************************/
void rmatrixlqunpackl(const ap::real_2d_array& a,
     int m,
     int n,
     ap::real_2d_array& l)
{
    int i;
    int k;

    if( m<=0||n<=0 )
    {
        return;
    }
    l.setlength(m, n);
    for(i = 0; i <= n-1; i++)
    {
        l(0,i) = 0;
    }
    for(i = 1; i <= m-1; i++)
    {
        ap::vmove(&l(i, 0), 1, &l(0, 0), 1, ap::vlen(0,n-1));
    }
    for(i = 0; i <= m-1; i++)
    {
        k = ap::minint(i, n-1);
        ap::vmove(&l(i, 0), 1, &a(i, 0), 1, ap::vlen(0,k));
    }
}


/*************************************************************************
Partial unpacking of matrix Q from QR decomposition of a complex matrix A.

Input parameters:
    A           -   matrices Q and R in compact form.
                    Output of CMatrixQR subroutine .
    M           -   number of rows in matrix A. M>=0.
    N           -   number of columns in matrix A. N>=0.
    Tau         -   scalar factors which are used to form Q.
                    Output of CMatrixQR subroutine .
    QColumns    -   required number of columns in matrix Q. M>=QColumns>=0.

Output parameters:
    Q           -   first QColumns columns of matrix Q.
                    Array whose index ranges within [0..M-1, 0..QColumns-1].
                    If QColumns=0, array isn't changed.

  -- ALGLIB routine --
     17.02.2010
     Bochkanov Sergey
*************************************************************************/
void cmatrixqrunpackq(const ap::complex_2d_array& a,
     int m,
     int n,
     const ap::complex_1d_array& tau,
     int qcolumns,
     ap::complex_2d_array& q)
{
    ap::complex_1d_array work;
    ap::complex_1d_array t;
    ap::complex_1d_array taubuf;
    int minmn;
    int refcnt;
    ap::complex_2d_array tmpa;
    ap::complex_2d_array tmpt;
    ap::complex_2d_array tmpr;
    int blockstart;
    int blocksize;
    int rowscount;
    int i;
    int j;
//    int k;
    ap::complex v;

    ap::ap_error::make_assertion(qcolumns<=m, "UnpackQFromQR: QColumns>M!");
    if( m<=0||n<=0 )
    {
        return;
    }
    
    //
    // init
    //
    minmn = ap::minint(m, n);
    refcnt = ap::minint(minmn, qcolumns);
    work.setlength(ap::maxint(m, n)+1);
    t.setlength(ap::maxint(m, n)+1);
    taubuf.setlength(minmn);
    tmpa.setlength(m, ablascomplexblocksize(a));
    tmpt.setlength(ablascomplexblocksize(a), ablascomplexblocksize(a));
    tmpr.setlength(2*ablascomplexblocksize(a), qcolumns);
    q.setlength(m, qcolumns);
    for(i = 0; i <= m-1; i++)
    {
        for(j = 0; j <= qcolumns-1; j++)
        {
            if( i==j )
            {
                q(i,j) = 1;
            }
            else
            {
                q(i,j) = 0;
            }
        }
    }
    
    //
    // Blocked code
    //
    blockstart = ablascomplexblocksize(a)*(refcnt/ablascomplexblocksize(a));
    blocksize = refcnt-blockstart;
    while(blockstart>=0)
    {
        rowscount = m-blockstart;
        
        //
        // QR decomposition of submatrix.
        // Matrix is copied to temporary storage to solve
        // some TLB issues arising from non-contiguous memory
        // access pattern.
        //
        cmatrixcopy(rowscount, blocksize, a, blockstart, blockstart, tmpa, 0, 0);
        ap::vmove(&taubuf(0), 1, &tau(blockstart), 1, "N", ap::vlen(0,blocksize-1));
        
        //
        // Update matrix, choose between:
        // a) Level 2 algorithm (when the rest of the matrix is small enough)
        // b) blocked algorithm, see algorithm 5 from  'A storage efficient WY
        //    representation for products of Householder transformations',
        //    by R. Schreiber and C. Van Loan.
        //
        if( qcolumns>=2*ablascomplexblocksize(a) )
        {
            
            //
            // Prepare block reflector
            //
            cmatrixblockreflector(tmpa, taubuf, true, rowscount, blocksize, tmpt, work);
            
            //
            // Multiply the rest of A by Q.
            //
            // Q  = E + Y*T*Y'  = E + TmpA*TmpT*TmpA'
            //
            cmatrixgemm(blocksize, qcolumns, rowscount, 1.0, tmpa, 0, 0, 2, q, blockstart, 0, 0, 0.0, tmpr, 0, 0);
            cmatrixgemm(blocksize, qcolumns, blocksize, 1.0, tmpt, 0, 0, 0, tmpr, 0, 0, 0, 0.0, tmpr, blocksize, 0);
            cmatrixgemm(rowscount, qcolumns, blocksize, 1.0, tmpa, 0, 0, 0, tmpr, blocksize, 0, 0, 1.0, q, blockstart, 0);
        }
        else
        {
            
            //
            // Level 2 algorithm
            //
            for(i = blocksize-1; i >= 0; i--)
            {
                ap::vmove(&t(1), 1, &tmpa(i, i), tmpa.getstride(), "N", ap::vlen(1,rowscount-i));
                t(1) = 1;
                complexapplyreflectionfromtheleft(q, taubuf(i), t, blockstart+i, m-1, 0, qcolumns-1, work);
            }
        }
        
        //
        // Advance
        //
        blockstart = blockstart-ablascomplexblocksize(a);
        blocksize = ablascomplexblocksize(a);
    }
}


/*************************************************************************
Unpacking of matrix R from the QR decomposition of a matrix A

Input parameters:
    A       -   matrices Q and R in compact form.
                Output of CMatrixQR subroutine.
    M       -   number of rows in given matrix A. M>=0.
    N       -   number of columns in given matrix A. N>=0.

Output parameters:
    R       -   matrix R, array[0..M-1, 0..N-1].

  -- ALGLIB routine --
     17.02.2010
     Bochkanov Sergey
*************************************************************************/
void cmatrixqrunpackr(const ap::complex_2d_array& a,
     int m,
     int n,
     ap::complex_2d_array& r)
{
    int i;
    int k;

    if( m<=0||n<=0 )
    {
        return;
    }
    k = ap::minint(m, n);
    r.setlength(m, n);
    for(i = 0; i <= n-1; i++)
    {
        r(0,i) = 0;
    }
    for(i = 1; i <= m-1; i++)
    {
        ap::vmove(&r(i, 0), 1, &r(0, 0), 1, "N", ap::vlen(0,n-1));
    }
    for(i = 0; i <= k-1; i++)
    {
        ap::vmove(&r(i, i), 1, &a(i, i), 1, "N", ap::vlen(i,n-1));
    }
}


/*************************************************************************
Partial unpacking of matrix Q from LQ decomposition of a complex matrix A.

Input parameters:
    A           -   matrices Q and R in compact form.
                    Output of CMatrixLQ subroutine .
    M           -   number of rows in matrix A. M>=0.
    N           -   number of columns in matrix A. N>=0.
    Tau         -   scalar factors which are used to form Q.
                    Output of CMatrixLQ subroutine .
    QRows       -   required number of rows in matrix Q. N>=QColumns>=0.

Output parameters:
    Q           -   first QRows rows of matrix Q.
                    Array whose index ranges within [0..QRows-1, 0..N-1].
                    If QRows=0, array isn't changed.

  -- ALGLIB routine --
     17.02.2010
     Bochkanov Sergey
*************************************************************************/
void cmatrixlqunpackq(const ap::complex_2d_array& a,
     int m,
     int n,
     const ap::complex_1d_array& tau,
     int qrows,
     ap::complex_2d_array& q)
{
    ap::complex_1d_array work;
    ap::complex_1d_array t;
    ap::complex_1d_array taubuf;
    int minmn;
    int refcnt;
    ap::complex_2d_array tmpa;
    ap::complex_2d_array tmpt;
    ap::complex_2d_array tmpr;
    int blockstart;
    int blocksize;
    int columnscount;
    int i;
    int j;
//    int k;
    ap::complex v;

    if( m<=0||n<=0 )
    {
        return;
    }
    
    //
    // Init
    //
    minmn = ap::minint(m, n);
    refcnt = ap::minint(minmn, qrows);
    work.setlength(ap::maxint(m, n)+1);
    t.setlength(ap::maxint(m, n)+1);
    taubuf.setlength(minmn);
    tmpa.setlength(ablascomplexblocksize(a), n);
    tmpt.setlength(ablascomplexblocksize(a), ablascomplexblocksize(a));
    tmpr.setlength(qrows, 2*ablascomplexblocksize(a));
    q.setlength(qrows, n);
    for(i = 0; i <= qrows-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            if( i==j )
            {
                q(i,j) = 1;
            }
            else
            {
                q(i,j) = 0;
            }
        }
    }
    
    //
    // Blocked code
    //
    blockstart = ablascomplexblocksize(a)*(refcnt/ablascomplexblocksize(a));
    blocksize = refcnt-blockstart;
    while(blockstart>=0)
    {
        columnscount = n-blockstart;
        
        //
        // LQ decomposition of submatrix.
        // Matrix is copied to temporary storage to solve
        // some TLB issues arising from non-contiguous memory
        // access pattern.
        //
        cmatrixcopy(blocksize, columnscount, a, blockstart, blockstart, tmpa, 0, 0);
        ap::vmove(&taubuf(0), 1, &tau(blockstart), 1, "N", ap::vlen(0,blocksize-1));
        
        //
        // Update matrix, choose between:
        // a) Level 2 algorithm (when the rest of the matrix is small enough)
        // b) blocked algorithm, see algorithm 5 from  'A storage efficient WY
        //    representation for products of Householder transformations',
        //    by R. Schreiber and C. Van Loan.
        //
        if( qrows>=2*ablascomplexblocksize(a) )
        {
            
            //
            // Prepare block reflector
            //
            cmatrixblockreflector(tmpa, taubuf, false, columnscount, blocksize, tmpt, work);
            
            //
            // Multiply the rest of A by Q'.
            //
            // Q'  = E + Y*T'*Y'  = E + TmpA'*TmpT'*TmpA
            //
            cmatrixgemm(qrows, blocksize, columnscount, 1.0, q, 0, blockstart, 0, tmpa, 0, 0, 2, 0.0, tmpr, 0, 0);
            cmatrixgemm(qrows, blocksize, blocksize, 1.0, tmpr, 0, 0, 0, tmpt, 0, 0, 2, 0.0, tmpr, 0, blocksize);
            cmatrixgemm(qrows, columnscount, blocksize, 1.0, tmpr, 0, blocksize, 0, tmpa, 0, 0, 0, 1.0, q, 0, blockstart);
        }
        else
        {
            
            //
            // Level 2 algorithm
            //
            for(i = blocksize-1; i >= 0; i--)
            {
                ap::vmove(&t(1), 1, &tmpa(i, i), 1, "Conj", ap::vlen(1,columnscount-i));
                t(1) = 1;
                complexapplyreflectionfromtheright(q, ap::conj(taubuf(i)), t, 0, qrows-1, blockstart+i, n-1, work);
            }
        }
        
        //
        // Advance
        //
        blockstart = blockstart-ablascomplexblocksize(a);
        blocksize = ablascomplexblocksize(a);
    }
}


/*************************************************************************
Unpacking of matrix L from the LQ decomposition of a matrix A

Input parameters:
    A       -   matrices Q and L in compact form.
                Output of CMatrixLQ subroutine.
    M       -   number of rows in given matrix A. M>=0.
    N       -   number of columns in given matrix A. N>=0.

Output parameters:
    L       -   matrix L, array[0..M-1, 0..N-1].

  -- ALGLIB routine --
     17.02.2010
     Bochkanov Sergey
*************************************************************************/
void cmatrixlqunpackl(const ap::complex_2d_array& a,
     int m,
     int n,
     ap::complex_2d_array& l)
{
    int i;
    int k;

    if( m<=0||n<=0 )
    {
        return;
    }
    l.setlength(m, n);
    for(i = 0; i <= n-1; i++)
    {
        l(0,i) = 0;
    }
    for(i = 1; i <= m-1; i++)
    {
        ap::vmove(&l(i, 0), 1, &l(0, 0), 1, "N", ap::vlen(0,n-1));
    }
    for(i = 0; i <= m-1; i++)
    {
        k = ap::minint(i, n-1);
        ap::vmove(&l(i, 0), 1, &a(i, 0), 1, "N", ap::vlen(0,k));
    }
}


/*************************************************************************
Reduction of a rectangular matrix to  bidiagonal form

The algorithm reduces the rectangular matrix A to  bidiagonal form by
orthogonal transformations P and Q: A = Q*B*P.

Input parameters:
    A       -   source matrix. array[0..M-1, 0..N-1]
    M       -   number of rows in matrix A.
    N       -   number of columns in matrix A.

Output parameters:
    A       -   matrices Q, B, P in compact form (see below).
    TauQ    -   scalar factors which are used to form matrix Q.
    TauP    -   scalar factors which are used to form matrix P.

The main diagonal and one of the  secondary  diagonals  of  matrix  A  are
replaced with bidiagonal  matrix  B.  Other  elements  contain  elementary
reflections which form MxM matrix Q and NxN matrix P, respectively.

If M>=N, B is the upper  bidiagonal  MxN  matrix  and  is  stored  in  the
corresponding  elements  of  matrix  A.  Matrix  Q  is  represented  as  a
product   of   elementary   reflections   Q = H(0)*H(1)*...*H(n-1),  where
H(i) = 1-tau*v*v'. Here tau is a scalar which is stored  in  TauQ[i],  and
vector v has the following  structure:  v(0:i-1)=0, v(i)=1, v(i+1:m-1)  is
stored   in   elements   A(i+1:m-1,i).   Matrix   P  is  as  follows:  P =
G(0)*G(1)*...*G(n-2), where G(i) = 1 - tau*u*u'. Tau is stored in TauP[i],
u(0:i)=0, u(i+1)=1, u(i+2:n-1) is stored in elements A(i,i+2:n-1).

If M<N, B is the  lower  bidiagonal  MxN  matrix  and  is  stored  in  the
corresponding   elements  of  matrix  A.  Q = H(0)*H(1)*...*H(m-2),  where
H(i) = 1 - tau*v*v', tau is stored in TauQ, v(0:i)=0, v(i+1)=1, v(i+2:m-1)
is    stored    in   elements   A(i+2:m-1,i).    P = G(0)*G(1)*...*G(m-1),
G(i) = 1-tau*u*u', tau is stored in  TauP,  u(0:i-1)=0, u(i)=1, u(i+1:n-1)
is stored in A(i,i+1:n-1).

EXAMPLE:

m=6, n=5 (m > n):               m=5, n=6 (m < n):

(  d   e   u1  u1  u1 )         (  d   u1  u1  u1  u1  u1 )
(  v1  d   e   u2  u2 )         (  e   d   u2  u2  u2  u2 )
(  v1  v2  d   e   u3 )         (  v1  e   d   u3  u3  u3 )
(  v1  v2  v3  d   e  )         (  v1  v2  e   d   u4  u4 )
(  v1  v2  v3  v4  d  )         (  v1  v2  v3  e   d   u5 )
(  v1  v2  v3  v4  v5 )

Here vi and ui are vectors which form H(i) and G(i), and d and e -
are the diagonal and off-diagonal elements of matrix B.

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     September 30, 1994.
     Sergey Bochkanov, ALGLIB project, translation from FORTRAN to
     pseudocode, 2007-2010.
*************************************************************************/
void rmatrixbd(ap::real_2d_array& a,
     int m,
     int n,
     ap::real_1d_array& tauq,
     ap::real_1d_array& taup)
{
    ap::real_1d_array work;
    ap::real_1d_array t;
    int minmn;
    int maxmn;
    int i;
    double ltau;

    
    //
    // Prepare
    //
    if( n<=0||m<=0 )
    {
        return;
    }
    minmn = ap::minint(m, n);
    maxmn = ap::maxint(m, n);
    work.setlength(maxmn+1);
    t.setlength(maxmn+1);
    if( m>=n )
    {
        tauq.setlength(n);
        taup.setlength(n);
    }
    else
    {
        tauq.setlength(m);
        taup.setlength(m);
    }
    if( m>=n )
    {
        
        //
        // Reduce to upper bidiagonal form
        //
        for(i = 0; i <= n-1; i++)
        {
            
            //
            // Generate elementary reflector H(i) to annihilate A(i+1:m-1,i)
            //
            ap::vmove(&t(1), 1, &a(i, i), a.getstride(), ap::vlen(1,m-i));
            generatereflection(t, m-i, ltau);
            tauq(i) = ltau;
            ap::vmove(&a(i, i), a.getstride(), &t(1), 1, ap::vlen(i,m-1));
            t(1) = 1;
            
            //
            // Apply H(i) to A(i:m-1,i+1:n-1) from the left
            //
            applyreflectionfromtheleft(a, ltau, t, i, m-1, i+1, n-1, work);
            if( i<n-1 )
            {
                
                //
                // Generate elementary reflector G(i) to annihilate
                // A(i,i+2:n-1)
                //
                ap::vmove(&t(1), 1, &a(i, i+1), 1, ap::vlen(1,n-i-1));
                generatereflection(t, n-1-i, ltau);
                taup(i) = ltau;
                ap::vmove(&a(i, i+1), 1, &t(1), 1, ap::vlen(i+1,n-1));
                t(1) = 1;
                
                //
                // Apply G(i) to A(i+1:m-1,i+1:n-1) from the right
                //
                applyreflectionfromtheright(a, ltau, t, i+1, m-1, i+1, n-1, work);
            }
            else
            {
                taup(i) = 0;
            }
        }
    }
    else
    {
        
        //
        // Reduce to lower bidiagonal form
        //
        for(i = 0; i <= m-1; i++)
        {
            
            //
            // Generate elementary reflector G(i) to annihilate A(i,i+1:n-1)
            //
            ap::vmove(&t(1), 1, &a(i, i), 1, ap::vlen(1,n-i));
            generatereflection(t, n-i, ltau);
            taup(i) = ltau;
            ap::vmove(&a(i, i), 1, &t(1), 1, ap::vlen(i,n-1));
            t(1) = 1;
            
            //
            // Apply G(i) to A(i+1:m-1,i:n-1) from the right
            //
            applyreflectionfromtheright(a, ltau, t, i+1, m-1, i, n-1, work);
            if( i<m-1 )
            {
                
                //
                // Generate elementary reflector H(i) to annihilate
                // A(i+2:m-1,i)
                //
                ap::vmove(&t(1), 1, &a(i+1, i), a.getstride(), ap::vlen(1,m-1-i));
                generatereflection(t, m-1-i, ltau);
                tauq(i) = ltau;
                ap::vmove(&a(i+1, i), a.getstride(), &t(1), 1, ap::vlen(i+1,m-1));
                t(1) = 1;
                
                //
                // Apply H(i) to A(i+1:m-1,i+1:n-1) from the left
                //
                applyreflectionfromtheleft(a, ltau, t, i+1, m-1, i+1, n-1, work);
            }
            else
            {
                tauq(i) = 0;
            }
        }
    }
}


/*************************************************************************
Unpacking matrix Q which reduces a matrix to bidiagonal form.

Input parameters:
    QP          -   matrices Q and P in compact form.
                    Output of ToBidiagonal subroutine.
    M           -   number of rows in matrix A.
    N           -   number of columns in matrix A.
    TAUQ        -   scalar factors which are used to form Q.
                    Output of ToBidiagonal subroutine.
    QColumns    -   required number of columns in matrix Q.
                    M>=QColumns>=0.

Output parameters:
    Q           -   first QColumns columns of matrix Q.
                    Array[0..M-1, 0..QColumns-1]
                    If QColumns=0, the array is not modified.

  -- ALGLIB --
     2005-2010
     Bochkanov Sergey
*************************************************************************/
void rmatrixbdunpackq(const ap::real_2d_array& qp,
     int m,
     int n,
     const ap::real_1d_array& tauq,
     int qcolumns,
     ap::real_2d_array& q)
{
    int i;
    int j;

    ap::ap_error::make_assertion(qcolumns<=m, "RMatrixBDUnpackQ: QColumns>M!");
    ap::ap_error::make_assertion(qcolumns>=0, "RMatrixBDUnpackQ: QColumns<0!");
    if( m==0||n==0||qcolumns==0 )
    {
        return;
    }
    
    //
    // prepare Q
    //
    q.setlength(m, qcolumns);
    for(i = 0; i <= m-1; i++)
    {
        for(j = 0; j <= qcolumns-1; j++)
        {
            if( i==j )
            {
                q(i,j) = 1;
            }
            else
            {
                q(i,j) = 0;
            }
        }
    }
    
    //
    // Calculate
    //
    rmatrixbdmultiplybyq(qp, m, n, tauq, q, m, qcolumns, false, false);
}


/*************************************************************************
Multiplication by matrix Q which reduces matrix A to  bidiagonal form.

The algorithm allows pre- or post-multiply by Q or Q'.

Input parameters:
    QP          -   matrices Q and P in compact form.
                    Output of ToBidiagonal subroutine.
    M           -   number of rows in matrix A.
    N           -   number of columns in matrix A.
    TAUQ        -   scalar factors which are used to form Q.
                    Output of ToBidiagonal subroutine.
    Z           -   multiplied matrix.
                    array[0..ZRows-1,0..ZColumns-1]
    ZRows       -   number of rows in matrix Z. If FromTheRight=False,
                    ZRows=M, otherwise ZRows can be arbitrary.
    ZColumns    -   number of columns in matrix Z. If FromTheRight=True,
                    ZColumns=M, otherwise ZColumns can be arbitrary.
    FromTheRight -  pre- or post-multiply.
    DoTranspose -   multiply by Q or Q'.

Output parameters:
    Z           -   product of Z and Q.
                    Array[0..ZRows-1,0..ZColumns-1]
                    If ZRows=0 or ZColumns=0, the array is not modified.

  -- ALGLIB --
     2005-2010
     Bochkanov Sergey
*************************************************************************/
void rmatrixbdmultiplybyq(const ap::real_2d_array& qp,
     int m,
     int n,
     const ap::real_1d_array& tauq,
     ap::real_2d_array& z,
     int zrows,
     int zcolumns,
     bool fromtheright,
     bool dotranspose)
{
    int i;
    int i1;
    int i2;
    int istep;
    ap::real_1d_array v;
    ap::real_1d_array work;
    int mx;

    if( m<=0||n<=0||zrows<=0||zcolumns<=0 )
    {
        return;
    }
    ap::ap_error::make_assertion((fromtheright&&zcolumns==m)||((!fromtheright)&&(zrows==m)), "RMatrixBDMultiplyByQ: incorrect Z size!");
    
    //
    // init
    //
    mx = ap::maxint(m, n);
    mx = ap::maxint(mx, zrows);
    mx = ap::maxint(mx, zcolumns);
    v.setlength(mx+1);
    work.setlength(mx+1);
    if( m>=n )
    {
        
        //
        // setup
        //
        if( fromtheright )
        {
            i1 = 0;
            i2 = n-1;
            istep = +1;
        }
        else
        {
            i1 = n-1;
            i2 = 0;
            istep = -1;
        }
        if( dotranspose )
        {
            i = i1;
            i1 = i2;
            i2 = i;
            istep = -istep;
        }
        
        //
        // Process
        //
        i = i1;
        do
        {
            ap::vmove(&v(1), 1, &qp(i, i), qp.getstride(), ap::vlen(1,m-i));
            v(1) = 1;
            if( fromtheright )
            {
                applyreflectionfromtheright(z, tauq(i), v, 0, zrows-1, i, m-1, work);
            }
            else
            {
                applyreflectionfromtheleft(z, tauq(i), v, i, m-1, 0, zcolumns-1, work);
            }
            i = i+istep;
        }
        while(i!=i2+istep);
    }
    else
    {
        
        //
        // setup
        //
        if( fromtheright )
        {
            i1 = 0;
            i2 = m-2;
            istep = +1;
        }
        else
        {
            i1 = m-2;
            i2 = 0;
            istep = -1;
        }
        if( dotranspose )
        {
            i = i1;
            i1 = i2;
            i2 = i;
            istep = -istep;
        }
        
        //
        // Process
        //
        if( m-1>0 )
        {
            i = i1;
            do
            {
                ap::vmove(&v(1), 1, &qp(i+1, i), qp.getstride(), ap::vlen(1,m-i-1));
                v(1) = 1;
                if( fromtheright )
                {
                    applyreflectionfromtheright(z, tauq(i), v, 0, zrows-1, i+1, m-1, work);
                }
                else
                {
                    applyreflectionfromtheleft(z, tauq(i), v, i+1, m-1, 0, zcolumns-1, work);
                }
                i = i+istep;
            }
            while(i!=i2+istep);
        }
    }
}


/*************************************************************************
Unpacking matrix P which reduces matrix A to bidiagonal form.
The subroutine returns transposed matrix P.

Input parameters:
    QP      -   matrices Q and P in compact form.
                Output of ToBidiagonal subroutine.
    M       -   number of rows in matrix A.
    N       -   number of columns in matrix A.
    TAUP    -   scalar factors which are used to form P.
                Output of ToBidiagonal subroutine.
    PTRows  -   required number of rows of matrix P^T. N >= PTRows >= 0.

Output parameters:
    PT      -   first PTRows columns of matrix P^T
                Array[0..PTRows-1, 0..N-1]
                If PTRows=0, the array is not modified.

  -- ALGLIB --
     2005-2010
     Bochkanov Sergey
*************************************************************************/
void rmatrixbdunpackpt(const ap::real_2d_array& qp,
     int m,
     int n,
     const ap::real_1d_array& taup,
     int ptrows,
     ap::real_2d_array& pt)
{
    int i;
    int j;

    ap::ap_error::make_assertion(ptrows<=n, "RMatrixBDUnpackPT: PTRows>N!");
    ap::ap_error::make_assertion(ptrows>=0, "RMatrixBDUnpackPT: PTRows<0!");
    if( m==0||n==0||ptrows==0 )
    {
        return;
    }
    
    //
    // prepare PT
    //
    pt.setlength(ptrows, n);
    for(i = 0; i <= ptrows-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            if( i==j )
            {
                pt(i,j) = 1;
            }
            else
            {
                pt(i,j) = 0;
            }
        }
    }
    
    //
    // Calculate
    //
    rmatrixbdmultiplybyp(qp, m, n, taup, pt, ptrows, n, true, true);
}


/*************************************************************************
Multiplication by matrix P which reduces matrix A to  bidiagonal form.

The algorithm allows pre- or post-multiply by P or P'.

Input parameters:
    QP          -   matrices Q and P in compact form.
                    Output of RMatrixBD subroutine.
    M           -   number of rows in matrix A.
    N           -   number of columns in matrix A.
    TAUP        -   scalar factors which are used to form P.
                    Output of RMatrixBD subroutine.
    Z           -   multiplied matrix.
                    Array whose indexes range within [0..ZRows-1,0..ZColumns-1].
    ZRows       -   number of rows in matrix Z. If FromTheRight=False,
                    ZRows=N, otherwise ZRows can be arbitrary.
    ZColumns    -   number of columns in matrix Z. If FromTheRight=True,
                    ZColumns=N, otherwise ZColumns can be arbitrary.
    FromTheRight -  pre- or post-multiply.
    DoTranspose -   multiply by P or P'.

Output parameters:
    Z - product of Z and P.
                Array whose indexes range within [0..ZRows-1,0..ZColumns-1].
                If ZRows=0 or ZColumns=0, the array is not modified.

  -- ALGLIB --
     2005-2010
     Bochkanov Sergey
*************************************************************************/
void rmatrixbdmultiplybyp(const ap::real_2d_array& qp,
     int m,
     int n,
     const ap::real_1d_array& taup,
     ap::real_2d_array& z,
     int zrows,
     int zcolumns,
     bool fromtheright,
     bool dotranspose)
{
    int i;
    ap::real_1d_array v;
    ap::real_1d_array work;
    int mx;
    int i1;
    int i2;
    int istep;

    if( m<=0||n<=0||zrows<=0||zcolumns<=0 )
    {
        return;
    }
    ap::ap_error::make_assertion((fromtheright&&(zcolumns==n))||((!fromtheright)&&(zrows==n)), "RMatrixBDMultiplyByP: incorrect Z size!");
    
    //
    // init
    //
    mx = ap::maxint(m, n);
    mx = ap::maxint(mx, zrows);
    mx = ap::maxint(mx, zcolumns);
    v.setlength(mx+1);
    work.setlength(mx+1);
    if( m>=n )
    {
        
        //
        // setup
        //
        if( fromtheright )
        {
            i1 = n-2;
            i2 = 0;
            istep = -1;
        }
        else
        {
            i1 = 0;
            i2 = n-2;
            istep = +1;
        }
        if( !dotranspose )
        {
            i = i1;
            i1 = i2;
            i2 = i;
            istep = -istep;
        }
        
        //
        // Process
        //
        if( n-1>0 )
        {
            i = i1;
            do
            {
                ap::vmove(&v(1), 1, &qp(i, i+1), 1, ap::vlen(1,n-1-i));
                v(1) = 1;
                if( fromtheright )
                {
                    applyreflectionfromtheright(z, taup(i), v, 0, zrows-1, i+1, n-1, work);
                }
                else
                {
                    applyreflectionfromtheleft(z, taup(i), v, i+1, n-1, 0, zcolumns-1, work);
                }
                i = i+istep;
            }
            while(i!=i2+istep);
        }
    }
    else
    {
        
        //
        // setup
        //
        if( fromtheright )
        {
            i1 = m-1;
            i2 = 0;
            istep = -1;
        }
        else
        {
            i1 = 0;
            i2 = m-1;
            istep = +1;
        }
        if( !dotranspose )
        {
            i = i1;
            i1 = i2;
            i2 = i;
            istep = -istep;
        }
        
        //
        // Process
        //
        i = i1;
        do
        {
            ap::vmove(&v(1), 1, &qp(i, i), 1, ap::vlen(1,n-i));
            v(1) = 1;
            if( fromtheright )
            {
                applyreflectionfromtheright(z, taup(i), v, 0, zrows-1, i, n-1, work);
            }
            else
            {
                applyreflectionfromtheleft(z, taup(i), v, i, n-1, 0, zcolumns-1, work);
            }
            i = i+istep;
        }
        while(i!=i2+istep);
    }
}


/*************************************************************************
Unpacking of the main and secondary diagonals of bidiagonal decomposition
of matrix A.

Input parameters:
    B   -   output of RMatrixBD subroutine.
    M   -   number of rows in matrix B.
    N   -   number of columns in matrix B.

Output parameters:
    IsUpper -   True, if the matrix is upper bidiagonal.
                otherwise IsUpper is False.
    D       -   the main diagonal.
                Array whose index ranges within [0..Min(M,N)-1].
    E       -   the secondary diagonal (upper or lower, depending on
                the value of IsUpper).
                Array index ranges within [0..Min(M,N)-1], the last
                element is not used.

  -- ALGLIB --
     2005-2010
     Bochkanov Sergey
*************************************************************************/
void rmatrixbdunpackdiagonals(const ap::real_2d_array& b,
     int m,
     int n,
     bool& isupper,
     ap::real_1d_array& d,
     ap::real_1d_array& e)
{
    int i;

    isupper = m>=n;
    if( m<=0||n<=0 )
    {
        return;
    }
    if( isupper )
    {
        d.setlength(n);
        e.setlength(n);
        for(i = 0; i <= n-2; i++)
        {
            d(i) = b(i,i);
            e(i) = b(i,i+1);
        }
        d(n-1) = b(n-1,n-1);
    }
    else
    {
        d.setlength(m);
        e.setlength(m);
        for(i = 0; i <= m-2; i++)
        {
            d(i) = b(i,i);
            e(i) = b(i+1,i);
        }
        d(m-1) = b(m-1,m-1);
    }
}


/*************************************************************************
Reduction of a square matrix to  upper Hessenberg form: Q'*A*Q = H,
where Q is an orthogonal matrix, H - Hessenberg matrix.

Input parameters:
    A       -   matrix A with elements [0..N-1, 0..N-1]
    N       -   size of matrix A.

Output parameters:
    A       -   matrices Q and P in  compact form (see below).
    Tau     -   array of scalar factors which are used to form matrix Q.
                Array whose index ranges within [0..N-2]

Matrix H is located on the main diagonal, on the lower secondary  diagonal
and above the main diagonal of matrix A. The elements which are used to
form matrix Q are situated in array Tau and below the lower secondary
diagonal of matrix A as follows:

Matrix Q is represented as a product of elementary reflections

Q = H(0)*H(2)*...*H(n-2),

where each H(i) is given by

H(i) = 1 - tau * v * (v^T)

where tau is a scalar stored in Tau[I]; v - is a real vector,
so that v(0:i) = 0, v(i+1) = 1, v(i+2:n-1) stored in A(i+2:n-1,i).

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     October 31, 1992
*************************************************************************/
void rmatrixhessenberg(ap::real_2d_array& a, int n, ap::real_1d_array& tau)
{
    int i;
    double v;
    ap::real_1d_array t;
    ap::real_1d_array work;

    ap::ap_error::make_assertion(n>=0, "RMatrixHessenberg: incorrect N!");
    
    //
    // Quick return if possible
    //
    if( n<=1 )
    {
        return;
    }
    tau.setbounds(0, n-2);
    t.setbounds(1, n);
    work.setbounds(0, n-1);
    for(i = 0; i <= n-2; i++)
    {
        
        //
        // Compute elementary reflector H(i) to annihilate A(i+2:ihi,i)
        //
        ap::vmove(&t(1), 1, &a(i+1, i), a.getstride(), ap::vlen(1,n-i-1));
        generatereflection(t, n-i-1, v);
        ap::vmove(&a(i+1, i), a.getstride(), &t(1), 1, ap::vlen(i+1,n-1));
        tau(i) = v;
        t(1) = 1;
        
        //
        // Apply H(i) to A(1:ihi,i+1:ihi) from the right
        //
        applyreflectionfromtheright(a, v, t, 0, n-1, i+1, n-1, work);
        
        //
        // Apply H(i) to A(i+1:ihi,i+1:n) from the left
        //
        applyreflectionfromtheleft(a, v, t, i+1, n-1, i+1, n-1, work);
    }
}


/*************************************************************************
Unpacking matrix Q which reduces matrix A to upper Hessenberg form

Input parameters:
    A   -   output of RMatrixHessenberg subroutine.
    N   -   size of matrix A.
    Tau -   scalar factors which are used to form Q.
            Output of RMatrixHessenberg subroutine.

Output parameters:
    Q   -   matrix Q.
            Array whose indexes range within [0..N-1, 0..N-1].

  -- ALGLIB --
     2005-2010
     Bochkanov Sergey
*************************************************************************/
void rmatrixhessenbergunpackq(const ap::real_2d_array& a,
     int n,
     const ap::real_1d_array& tau,
     ap::real_2d_array& q)
{
    int i;
    int j;
    ap::real_1d_array v;
    ap::real_1d_array work;

    if( n==0 )
    {
        return;
    }
    
    //
    // init
    //
    q.setbounds(0, n-1, 0, n-1);
    v.setbounds(0, n-1);
    work.setbounds(0, n-1);
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            if( i==j )
            {
                q(i,j) = 1;
            }
            else
            {
                q(i,j) = 0;
            }
        }
    }
    
    //
    // unpack Q
    //
    for(i = 0; i <= n-2; i++)
    {
        
        //
        // Apply H(i)
        //
        ap::vmove(&v(1), 1, &a(i+1, i), a.getstride(), ap::vlen(1,n-i-1));
        v(1) = 1;
        applyreflectionfromtheright(q, tau(i), v, 0, n-1, i+1, n-1, work);
    }
}


/*************************************************************************
Unpacking matrix H (the result of matrix A reduction to upper Hessenberg form)

Input parameters:
    A   -   output of RMatrixHessenberg subroutine.
    N   -   size of matrix A.

Output parameters:
    H   -   matrix H. Array whose indexes range within [0..N-1, 0..N-1].

  -- ALGLIB --
     2005-2010
     Bochkanov Sergey
*************************************************************************/
void rmatrixhessenbergunpackh(const ap::real_2d_array& a,
     int n,
     ap::real_2d_array& h)
{
    int i;
    int j;
    ap::real_1d_array v;
    ap::real_1d_array work;

    if( n==0 )
    {
        return;
    }
    h.setbounds(0, n-1, 0, n-1);
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= i-2; j++)
        {
            h(i,j) = 0;
        }
        j = ap::maxint(0, i-1);
        ap::vmove(&h(i, j), 1, &a(i, j), 1, ap::vlen(j,n-1));
    }
}


/*************************************************************************
Reduction of a symmetric matrix which is given by its higher or lower
triangular part to a tridiagonal matrix using orthogonal similarity
transformation: Q'*A*Q=T.

Input parameters:
    A       -   matrix to be transformed
                array with elements [0..N-1, 0..N-1].
    N       -   size of matrix A.
    IsUpper -   storage format. If IsUpper = True, then matrix A is given
                by its upper triangle, and the lower triangle is not used
                and not modified by the algorithm, and vice versa
                if IsUpper = False.

Output parameters:
    A       -   matrices T and Q in  compact form (see lower)
    Tau     -   array of factors which are forming matrices H(i)
                array with elements [0..N-2].
    D       -   main diagonal of symmetric matrix T.
                array with elements [0..N-1].
    E       -   secondary diagonal of symmetric matrix T.
                array with elements [0..N-2].


  If IsUpper=True, the matrix Q is represented as a product of elementary
  reflectors

     Q = H(n-2) . . . H(2) H(0).

  Each H(i) has the form

     H(i) = I - tau * v * v'

  where tau is a real scalar, and v is a real vector with
  v(i+1:n-1) = 0, v(i) = 1, v(0:i-1) is stored on exit in
  A(0:i-1,i+1), and tau in TAU(i).

  If IsUpper=False, the matrix Q is represented as a product of elementary
  reflectors

     Q = H(0) H(2) . . . H(n-2).

  Each H(i) has the form

     H(i) = I - tau * v * v'

  where tau is a real scalar, and v is a real vector with
  v(0:i) = 0, v(i+1) = 1, v(i+2:n-1) is stored on exit in A(i+2:n-1,i),
  and tau in TAU(i).

  The contents of A on exit are illustrated by the following examples
  with n = 5:

  if UPLO = 'U':                       if UPLO = 'L':

    (  d   e   v1  v2  v3 )              (  d                  )
    (      d   e   v2  v3 )              (  e   d              )
    (          d   e   v3 )              (  v0  e   d          )
    (              d   e  )              (  v0  v1  e   d      )
    (                  d  )              (  v0  v1  v2  e   d  )

  where d and e denote diagonal and off-diagonal elements of T, and vi
  denotes an element of the vector defining H(i).

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     October 31, 1992
*************************************************************************/
void smatrixtd(ap::real_2d_array& a,
     int n,
     bool isupper,
     ap::real_1d_array& tau,
     ap::real_1d_array& d,
     ap::real_1d_array& e)
{
    int i;
    double alpha;
    double taui;
    double v;
    ap::real_1d_array t;
    ap::real_1d_array t2;
    ap::real_1d_array t3;

    if( n<=0 )
    {
        return;
    }
    t.setbounds(1, n);
    t2.setbounds(1, n);
    t3.setbounds(1, n);
    if( n>1 )
    {
        tau.setbounds(0, n-2);
    }
    d.setbounds(0, n-1);
    if( n>1 )
    {
        e.setbounds(0, n-2);
    }
    if( isupper )
    {
        
        //
        // Reduce the upper triangle of A
        //
        for(i = n-2; i >= 0; i--)
        {
            
            //
            // Generate elementary reflector H() = E - tau * v * v'
            //
            if( i>=1 )
            {
                ap::vmove(&t(2), 1, &a(0, i+1), a.getstride(), ap::vlen(2,i+1));
            }
            t(1) = a(i,i+1);
            generatereflection(t, i+1, taui);
            if( i>=1 )
            {
                ap::vmove(&a(0, i+1), a.getstride(), &t(2), 1, ap::vlen(0,i-1));
            }
            a(i,i+1) = t(1);
            e(i) = a(i,i+1);
            if( ap::fp_neq(taui,0) )
            {
                
                //
                // Apply H from both sides to A
                //
                a(i,i+1) = 1;
                
                //
                // Compute  x := tau * A * v  storing x in TAU
                //
                ap::vmove(&t(1), 1, &a(0, i+1), a.getstride(), ap::vlen(1,i+1));
                symmetricmatrixvectormultiply(a, isupper, 0, i, t, taui, t3);
                ap::vmove(&tau(0), 1, &t3(1), 1, ap::vlen(0,i));
                
                //
                // Compute  w := x - 1/2 * tau * (x'*v) * v
                //
                v = ap::vdotproduct(&tau(0), 1, &a(0, i+1), a.getstride(), ap::vlen(0,i));
                alpha = -0.5*taui*v;
                ap::vadd(&tau(0), 1, &a(0, i+1), a.getstride(), ap::vlen(0,i), alpha);
                
                //
                // Apply the transformation as a rank-2 update:
                //    A := A - v * w' - w * v'
                //
                ap::vmove(&t(1), 1, &a(0, i+1), a.getstride(), ap::vlen(1,i+1));
                ap::vmove(&t3(1), 1, &tau(0), 1, ap::vlen(1,i+1));
                symmetricrank2update(a, isupper, 0, i, t, t3, t2, double(-1));
                a(i,i+1) = e(i);
            }
            d(i+1) = a(i+1,i+1);
            tau(i) = taui;
        }
        d(0) = a(0,0);
    }
    else
    {
        
        //
        // Reduce the lower triangle of A
        //
        for(i = 0; i <= n-2; i++)
        {
            
            //
            // Generate elementary reflector H = E - tau * v * v'
            //
            ap::vmove(&t(1), 1, &a(i+1, i), a.getstride(), ap::vlen(1,n-i-1));
            generatereflection(t, n-i-1, taui);
            ap::vmove(&a(i+1, i), a.getstride(), &t(1), 1, ap::vlen(i+1,n-1));
            e(i) = a(i+1,i);
            if( ap::fp_neq(taui,0) )
            {
                
                //
                // Apply H from both sides to A
                //
                a(i+1,i) = 1;
                
                //
                // Compute  x := tau * A * v  storing y in TAU
                //
                ap::vmove(&t(1), 1, &a(i+1, i), a.getstride(), ap::vlen(1,n-i-1));
                symmetricmatrixvectormultiply(a, isupper, i+1, n-1, t, taui, t2);
                ap::vmove(&tau(i), 1, &t2(1), 1, ap::vlen(i,n-2));
                
                //
                // Compute  w := x - 1/2 * tau * (x'*v) * v
                //
                v = ap::vdotproduct(&tau(i), 1, &a(i+1, i), a.getstride(), ap::vlen(i,n-2));
                alpha = -0.5*taui*v;
                ap::vadd(&tau(i), 1, &a(i+1, i), a.getstride(), ap::vlen(i,n-2), alpha);
                
                //
                // Apply the transformation as a rank-2 update:
                //     A := A - v * w' - w * v'
                //
                //
                ap::vmove(&t(1), 1, &a(i+1, i), a.getstride(), ap::vlen(1,n-i-1));
                ap::vmove(&t2(1), 1, &tau(i), 1, ap::vlen(1,n-i-1));
                symmetricrank2update(a, isupper, i+1, n-1, t, t2, t3, double(-1));
                a(i+1,i) = e(i);
            }
            d(i) = a(i,i);
            tau(i) = taui;
        }
        d(n-1) = a(n-1,n-1);
    }
}


/*************************************************************************
Unpacking matrix Q which reduces symmetric matrix to a tridiagonal
form.

Input parameters:
    A       -   the result of a SMatrixTD subroutine
    N       -   size of matrix A.
    IsUpper -   storage format (a parameter of SMatrixTD subroutine)
    Tau     -   the result of a SMatrixTD subroutine

Output parameters:
    Q       -   transformation matrix.
                array with elements [0..N-1, 0..N-1].

  -- ALGLIB --
     Copyright 2005-2010 by Bochkanov Sergey
*************************************************************************/
void smatrixtdunpackq(const ap::real_2d_array& a,
     const int& n,
     const bool& isupper,
     const ap::real_1d_array& tau,
     ap::real_2d_array& q)
{
    int i;
    int j;
    ap::real_1d_array v;
    ap::real_1d_array work;

    if( n==0 )
    {
        return;
    }
    
    //
    // init
    //
    q.setbounds(0, n-1, 0, n-1);
    v.setbounds(1, n);
    work.setbounds(0, n-1);
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            if( i==j )
            {
                q(i,j) = 1;
            }
            else
            {
                q(i,j) = 0;
            }
        }
    }
    
    //
    // unpack Q
    //
    if( isupper )
    {
        for(i = 0; i <= n-2; i++)
        {
            
            //
            // Apply H(i)
            //
            ap::vmove(&v(1), 1, &a(0, i+1), a.getstride(), ap::vlen(1,i+1));
            v(i+1) = 1;
            applyreflectionfromtheleft(q, tau(i), v, 0, i, 0, n-1, work);
        }
    }
    else
    {
        for(i = n-2; i >= 0; i--)
        {
            
            //
            // Apply H(i)
            //
            ap::vmove(&v(1), 1, &a(i+1, i), a.getstride(), ap::vlen(1,n-i-1));
            v(1) = 1;
            applyreflectionfromtheleft(q, tau(i), v, i+1, n-1, 0, n-1, work);
        }
    }
}


/*************************************************************************
Reduction of a Hermitian matrix which is given  by  its  higher  or  lower
triangular part to a real  tridiagonal  matrix  using  unitary  similarity
transformation: Q'*A*Q = T.

Input parameters:
    A       -   matrix to be transformed
                array with elements [0..N-1, 0..N-1].
    N       -   size of matrix A.
    IsUpper -   storage format. If IsUpper = True, then matrix A is  given
                by its upper triangle, and the lower triangle is not  used
                and not modified by the algorithm, and vice versa
                if IsUpper = False.

Output parameters:
    A       -   matrices T and Q in  compact form (see lower)
    Tau     -   array of factors which are forming matrices H(i)
                array with elements [0..N-2].
    D       -   main diagonal of real symmetric matrix T.
                array with elements [0..N-1].
    E       -   secondary diagonal of real symmetric matrix T.
                array with elements [0..N-2].


  If IsUpper=True, the matrix Q is represented as a product of elementary
  reflectors

     Q = H(n-2) . . . H(2) H(0).

  Each H(i) has the form

     H(i) = I - tau * v * v'

  where tau is a complex scalar, and v is a complex vector with
  v(i+1:n-1) = 0, v(i) = 1, v(0:i-1) is stored on exit in
  A(0:i-1,i+1), and tau in TAU(i).

  If IsUpper=False, the matrix Q is represented as a product of elementary
  reflectors

     Q = H(0) H(2) . . . H(n-2).

  Each H(i) has the form

     H(i) = I - tau * v * v'

  where tau is a complex scalar, and v is a complex vector with
  v(0:i) = 0, v(i+1) = 1, v(i+2:n-1) is stored on exit in A(i+2:n-1,i),
  and tau in TAU(i).

  The contents of A on exit are illustrated by the following examples
  with n = 5:

  if UPLO = 'U':                       if UPLO = 'L':

    (  d   e   v1  v2  v3 )              (  d                  )
    (      d   e   v2  v3 )              (  e   d              )
    (          d   e   v3 )              (  v0  e   d          )
    (              d   e  )              (  v0  v1  e   d      )
    (                  d  )              (  v0  v1  v2  e   d  )

where d and e denote diagonal and off-diagonal elements of T, and vi
denotes an element of the vector defining H(i).

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     October 31, 1992
*************************************************************************/
void hmatrixtd(ap::complex_2d_array& a,
     int n,
     bool isupper,
     ap::complex_1d_array& tau,
     ap::real_1d_array& d,
     ap::real_1d_array& e)
{
    int i;
    ap::complex alpha;
    ap::complex taui;
    ap::complex v;
    ap::complex_1d_array t;
    ap::complex_1d_array t2;
    ap::complex_1d_array t3;

    if( n<=0 )
    {
        return;
    }
    for(i = 0; i <= n-1; i++)
    {
        ap::ap_error::make_assertion(ap::fp_eq(a(i,i).y,0), "");
    }
    if( n>1 )
    {
        tau.setbounds(0, n-2);
        e.setbounds(0, n-2);
    }
    d.setbounds(0, n-1);
    t.setbounds(0, n-1);
    t2.setbounds(0, n-1);
    t3.setbounds(0, n-1);
    if( isupper )
    {
        
        //
        // Reduce the upper triangle of A
        //
        a(n-1,n-1) = a(n-1,n-1).x;
        for(i = n-2; i >= 0; i--)
        {
            
            //
            // Generate elementary reflector H = I+1 - tau * v * v'
            //
            alpha = a(i,i+1);
            t(1) = alpha;
            if( i>=1 )
            {
                ap::vmove(&t(2), 1, &a(0, i+1), a.getstride(), "N", ap::vlen(2,i+1));
            }
            complexgeneratereflection(t, i+1, taui);
            if( i>=1 )
            {
                ap::vmove(&a(0, i+1), a.getstride(), &t(2), 1, "N", ap::vlen(0,i-1));
            }
            alpha = t(1);
            e(i) = alpha.x;
            if( taui!=0 )
            {
                
                //
                // Apply H(I+1) from both sides to A
                //
                a(i,i+1) = 1;
                
                //
                // Compute  x := tau * A * v  storing x in TAU
                //
                ap::vmove(&t(1), 1, &a(0, i+1), a.getstride(), "N", ap::vlen(1,i+1));
                hermitianmatrixvectormultiply(a, isupper, 0, i, t, taui, t2);
                ap::vmove(&tau(0), 1, &t2(1), 1, "N", ap::vlen(0,i));
                
                //
                // Compute  w := x - 1/2 * tau * (x'*v) * v
                //
                v = ap::vdotproduct(&tau(0), 1, "Conj", &a(0, i+1), a.getstride(), "N", ap::vlen(0,i));
                alpha = -0.5*taui*v;
                ap::vadd(&tau(0), 1, &a(0, i+1), a.getstride(), "N", ap::vlen(0,i), alpha);
                
                //
                // Apply the transformation as a rank-2 update:
                //    A := A - v * w' - w * v'
                //
                ap::vmove(&t(1), 1, &a(0, i+1), a.getstride(), "N", ap::vlen(1,i+1));
                ap::vmove(&t3(1), 1, &tau(0), 1, "N", ap::vlen(1,i+1));
                hermitianrank2update(a, isupper, 0, i, t, t3, t2, -1);
            }
            else
            {
                a(i,i) = a(i,i).x;
            }
            a(i,i+1) = e(i);
            d(i+1) = a(i+1,i+1).x;
            tau(i) = taui;
        }
        d(0) = a(0,0).x;
    }
    else
    {
        
        //
        // Reduce the lower triangle of A
        //
        a(0,0) = a(0,0).x;
        for(i = 0; i <= n-2; i++)
        {
            
            //
            // Generate elementary reflector H = I - tau * v * v'
            //
            ap::vmove(&t(1), 1, &a(i+1, i), a.getstride(), "N", ap::vlen(1,n-i-1));
            complexgeneratereflection(t, n-i-1, taui);
            ap::vmove(&a(i+1, i), a.getstride(), &t(1), 1, "N", ap::vlen(i+1,n-1));
            e(i) = a(i+1,i).x;
            if( taui!=0 )
            {
                
                //
                // Apply H(i) from both sides to A(i+1:n,i+1:n)
                //
                a(i+1,i) = 1;
                
                //
                // Compute  x := tau * A * v  storing y in TAU
                //
                ap::vmove(&t(1), 1, &a(i+1, i), a.getstride(), "N", ap::vlen(1,n-i-1));
                hermitianmatrixvectormultiply(a, isupper, i+1, n-1, t, taui, t2);
                ap::vmove(&tau(i), 1, &t2(1), 1, "N", ap::vlen(i,n-2));
                
                //
                // Compute  w := x - 1/2 * tau * (x'*v) * v
                //
                v = ap::vdotproduct(&tau(i), 1, "Conj", &a(i+1, i), a.getstride(), "N", ap::vlen(i,n-2));
                alpha = -0.5*taui*v;
                ap::vadd(&tau(i), 1, &a(i+1, i), a.getstride(), "N", ap::vlen(i,n-2), alpha);
                
                //
                // Apply the transformation as a rank-2 update:
                // A := A - v * w' - w * v'
                //
                ap::vmove(&t(1), 1, &a(i+1, i), a.getstride(), "N", ap::vlen(1,n-i-1));
                ap::vmove(&t2(1), 1, &tau(i), 1, "N", ap::vlen(1,n-i-1));
                hermitianrank2update(a, isupper, i+1, n-1, t, t2, t3, -1);
            }
            else
            {
                a(i+1,i+1) = a(i+1,i+1).x;
            }
            a(i+1,i) = e(i);
            d(i) = a(i,i).x;
            tau(i) = taui;
        }
        d(n-1) = a(n-1,n-1).x;
    }
}


/*************************************************************************
Unpacking matrix Q which reduces a Hermitian matrix to a real  tridiagonal
form.

Input parameters:
    A       -   the result of a HMatrixTD subroutine
    N       -   size of matrix A.
    IsUpper -   storage format (a parameter of HMatrixTD subroutine)
    Tau     -   the result of a HMatrixTD subroutine

Output parameters:
    Q       -   transformation matrix.
                array with elements [0..N-1, 0..N-1].

  -- ALGLIB --
     Copyright 2005-2010 by Bochkanov Sergey
*************************************************************************/
void hmatrixtdunpackq(const ap::complex_2d_array& a,
     const int& n,
     const bool& isupper,
     const ap::complex_1d_array& tau,
     ap::complex_2d_array& q)
{
    int i;
    int j;
    ap::complex_1d_array v;
    ap::complex_1d_array work;

    if( n==0 )
    {
        return;
    }
    
    //
    // init
    //
    q.setbounds(0, n-1, 0, n-1);
    v.setbounds(1, n);
    work.setbounds(0, n-1);
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            if( i==j )
            {
                q(i,j) = 1;
            }
            else
            {
                q(i,j) = 0;
            }
        }
    }
    
    //
    // unpack Q
    //
    if( isupper )
    {
        for(i = 0; i <= n-2; i++)
        {
            
            //
            // Apply H(i)
            //
            ap::vmove(&v(1), 1, &a(0, i+1), a.getstride(), "N", ap::vlen(1,i+1));
            v(i+1) = 1;
            complexapplyreflectionfromtheleft(q, tau(i), v, 0, i, 0, n-1, work);
        }
    }
    else
    {
        for(i = n-2; i >= 0; i--)
        {
            
            //
            // Apply H(i)
            //
            ap::vmove(&v(1), 1, &a(i+1, i), a.getstride(), "N", ap::vlen(1,n-i-1));
            v(1) = 1;
            complexapplyreflectionfromtheleft(q, tau(i), v, i+1, n-1, 0, n-1, work);
        }
    }
}


/*************************************************************************
Base case for real QR

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     September 30, 1994.
     Sergey Bochkanov, ALGLIB project, translation from FORTRAN to
     pseudocode, 2007-2010.
*************************************************************************/
static void rmatrixqrbasecase(ap::real_2d_array& a,
     int m,
     int n,
     ap::real_1d_array& work,
     ap::real_1d_array& t,
     ap::real_1d_array& tau)
{
    int i;
    int k;
    int minmn;
    double tmp;

    minmn = ap::minint(m, n);
    
    //
    // Test the input arguments
    //
    k = minmn;
    for(i = 0; i <= k-1; i++)
    {
        
        //
        // Generate elementary reflector H(i) to annihilate A(i+1:m,i)
        //
        ap::vmove(&t(1), 1, &a(i, i), a.getstride(), ap::vlen(1,m-i));
        generatereflection(t, m-i, tmp);
        tau(i) = tmp;
        ap::vmove(&a(i, i), a.getstride(), &t(1), 1, ap::vlen(i,m-1));
        t(1) = 1;
        if( i<n )
        {
            
            //
            // Apply H(i) to A(i:m-1,i+1:n-1) from the left
            //
            applyreflectionfromtheleft(a, tau(i), t, i, m-1, i+1, n-1, work);
        }
    }
}


/*************************************************************************
Base case for real LQ

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     September 30, 1994.
     Sergey Bochkanov, ALGLIB project, translation from FORTRAN to
     pseudocode, 2007-2010.
*************************************************************************/
static void rmatrixlqbasecase(ap::real_2d_array& a,
     int m,
     int n,
     ap::real_1d_array& work,
     ap::real_1d_array& t,
     ap::real_1d_array& tau)
{
    int i;
    int k;
    int minmn;
    double tmp;

    minmn = ap::minint(m, n);
    k = ap::minint(m, n);
    for(i = 0; i <= k-1; i++)
    {
        
        //
        // Generate elementary reflector H(i) to annihilate A(i,i+1:n-1)
        //
        ap::vmove(&t(1), 1, &a(i, i), 1, ap::vlen(1,n-i));
        generatereflection(t, n-i, tmp);
        tau(i) = tmp;
        ap::vmove(&a(i, i), 1, &t(1), 1, ap::vlen(i,n-1));
        t(1) = 1;
        if( i<n )
        {
            
            //
            // Apply H(i) to A(i+1:m,i:n) from the right
            //
            applyreflectionfromtheright(a, tau(i), t, i+1, m-1, i, n-1, work);
        }
    }
}


/*************************************************************************
Base case for complex QR

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     September 30, 1994.
     Sergey Bochkanov, ALGLIB project, translation from FORTRAN to
     pseudocode, 2007-2010.
*************************************************************************/
static void cmatrixqrbasecase(ap::complex_2d_array& a,
     int m,
     int n,
     ap::complex_1d_array& work,
     ap::complex_1d_array& t,
     ap::complex_1d_array& tau)
{
    int i;
    int k;
    int mmi;
    int minmn;
    ap::complex tmp;

    minmn = ap::minint(m, n);
    if( minmn<=0 )
    {
        return;
    }
    
    //
    // Test the input arguments
    //
    k = ap::minint(m, n);
    for(i = 0; i <= k-1; i++)
    {
        
        //
        // Generate elementary reflector H(i) to annihilate A(i+1:m,i)
        //
        mmi = m-i;
        ap::vmove(&t(1), 1, &a(i, i), a.getstride(), "N", ap::vlen(1,mmi));
        complexgeneratereflection(t, mmi, tmp);
        tau(i) = tmp;
        ap::vmove(&a(i, i), a.getstride(), &t(1), 1, "N", ap::vlen(i,m-1));
        t(1) = 1;
        if( i<n-1 )
        {
            
            //
            // Apply H'(i) to A(i:m,i+1:n) from the left
            //
            complexapplyreflectionfromtheleft(a, ap::conj(tau(i)), t, i, m-1, i+1, n-1, work);
        }
    }
}


/*************************************************************************
Base case for complex LQ

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     September 30, 1994.
     Sergey Bochkanov, ALGLIB project, translation from FORTRAN to
     pseudocode, 2007-2010.
*************************************************************************/
static void cmatrixlqbasecase(ap::complex_2d_array& a,
     int m,
     int n,
     ap::complex_1d_array& work,
     ap::complex_1d_array& t,
     ap::complex_1d_array& tau)
{
    int i;
    int minmn;
    ap::complex tmp;

    minmn = ap::minint(m, n);
    if( minmn<=0 )
    {
        return;
    }
    
    //
    // Test the input arguments
    //
    for(i = 0; i <= minmn-1; i++)
    {
        
        //
        // Generate elementary reflector H(i)
        //
        // NOTE: ComplexGenerateReflection() generates left reflector,
        // i.e. H which reduces x by applyiong from the left, but we
        // need RIGHT reflector. So we replace H=E-tau*v*v' by H^H,
        // which changes v to conj(v).
        //
        ap::vmove(&t(1), 1, &a(i, i), 1, "Conj", ap::vlen(1,n-i));
        complexgeneratereflection(t, n-i, tmp);
        tau(i) = tmp;
        ap::vmove(&a(i, i), 1, &t(1), 1, "Conj", ap::vlen(i,n-1));
        t(1) = 1;
        if( i<m-1 )
        {
            
            //
            // Apply H'(i)
            //
            complexapplyreflectionfromtheright(a, tau(i), t, i+1, m-1, i, n-1, work);
        }
    }
}


/*************************************************************************
Generate block reflector:
* fill unused parts of reflectors matrix by zeros
* fill diagonal of reflectors matrix by ones
* generate triangular factor T

PARAMETERS:
    A           -   either LengthA*BlockSize (if ColumnwiseA) or
                    BlockSize*LengthA (if not ColumnwiseA) matrix of
                    elementary reflectors.
                    Modified on exit.
    Tau         -   scalar factors
    ColumnwiseA -   reflectors are stored in rows or in columns
    LengthA     -   length of largest reflector
    BlockSize   -   number of reflectors
    T           -   array[BlockSize,2*BlockSize]. Left BlockSize*BlockSize
                    submatrix stores triangular factor on exit.
    WORK        -   array[BlockSize]
    
  -- ALGLIB routine --
     17.02.2010
     Bochkanov Sergey
*************************************************************************/
static void rmatrixblockreflector(ap::real_2d_array& a,
     ap::real_1d_array& tau,
     bool columnwisea,
     int lengtha,
     int blocksize,
     ap::real_2d_array& t,
     ap::real_1d_array& work)
{
    int i;
    int j;
    int k;
    double v;

    
    //
    // fill beginning of new column with zeros,
    // load 1.0 in the first non-zero element
    //
    for(k = 0; k <= blocksize-1; k++)
    {
        if( columnwisea )
        {
            for(i = 0; i <= k-1; i++)
            {
                a(i,k) = 0;
            }
        }
        else
        {
            for(i = 0; i <= k-1; i++)
            {
                a(k,i) = 0;
            }
        }
        a(k,k) = 1;
    }
    
    //
    // Calculate Gram matrix of A
    //
    for(i = 0; i <= blocksize-1; i++)
    {
        for(j = 0; j <= blocksize-1; j++)
        {
            t(i,blocksize+j) = 0;
        }
    }
    for(k = 0; k <= lengtha-1; k++)
    {
        for(j = 1; j <= blocksize-1; j++)
        {
            if( columnwisea )
            {
                v = a(k,j);
                if( ap::fp_neq(v,0) )
                {
                    ap::vadd(&t(j, blocksize), 1, &a(k, 0), 1, ap::vlen(blocksize,blocksize+j-1), v);
                }
            }
            else
            {
                v = a(j,k);
                if( ap::fp_neq(v,0) )
                {
                    ap::vadd(&t(j, blocksize), 1, &a(0, k), a.getstride(), ap::vlen(blocksize,blocksize+j-1), v);
                }
            }
        }
    }
    
    //
    // Prepare Y (stored in TmpA) and T (stored in TmpT)
    //
    for(k = 0; k <= blocksize-1; k++)
    {
        
        //
        // fill non-zero part of T, use pre-calculated Gram matrix
        //
        ap::vmove(&work(0), 1, &t(k, blocksize), 1, ap::vlen(0,k-1));
        for(i = 0; i <= k-1; i++)
        {
            v = ap::vdotproduct(&t(i, i), 1, &work(i), 1, ap::vlen(i,k-1));
            t(i,k) = -tau(k)*v;
        }
        t(k,k) = -tau(k);
        
        //
        // Rest of T is filled by zeros
        //
        for(i = k+1; i <= blocksize-1; i++)
        {
            t(i,k) = 0;
        }
    }
}


/*************************************************************************
Generate block reflector (complex):
* fill unused parts of reflectors matrix by zeros
* fill diagonal of reflectors matrix by ones
* generate triangular factor T


  -- ALGLIB routine --
     17.02.2010
     Bochkanov Sergey
*************************************************************************/
static void cmatrixblockreflector(ap::complex_2d_array& a,
     ap::complex_1d_array& tau,
     bool columnwisea,
     int lengtha,
     int blocksize,
     ap::complex_2d_array& t,
     ap::complex_1d_array& work)
{
    int i;
    int k;
    ap::complex v;

    
    //
    // Prepare Y (stored in TmpA) and T (stored in TmpT)
    //
    for(k = 0; k <= blocksize-1; k++)
    {
        
        //
        // fill beginning of new column with zeros,
        // load 1.0 in the first non-zero element
        //
        if( columnwisea )
        {
            for(i = 0; i <= k-1; i++)
            {
                a(i,k) = 0;
            }
        }
        else
        {
            for(i = 0; i <= k-1; i++)
            {
                a(k,i) = 0;
            }
        }
        a(k,k) = 1;
        
        //
        // fill non-zero part of T,
        //
        for(i = 0; i <= k-1; i++)
        {
            if( columnwisea )
            {
                v = ap::vdotproduct(&a(k, i), a.getstride(), "Conj", &a(k, k), a.getstride(), "N", ap::vlen(k,lengtha-1));
            }
            else
            {
                v = ap::vdotproduct(&a(i, k), 1, "N", &a(k, k), 1, "Conj", ap::vlen(k,lengtha-1));
            }
            work(i) = v;
        }
        for(i = 0; i <= k-1; i++)
        {
            v = ap::vdotproduct(&t(i, i), 1, "N", &work(i), 1, "N", ap::vlen(i,k-1));
            t(i,k) = -tau(k)*v;
        }
        t(k,k) = -tau(k);
        
        //
        // Rest of T is filled by zeros
        //
        for(i = k+1; i <= blocksize-1; i++)
        {
            t(i,k) = 0;
        }
    }
}





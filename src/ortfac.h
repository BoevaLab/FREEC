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

#ifndef _ortfac_h
#define _ortfac_h

#include "ap.h"
#include "ialglib.h"

#include "hblas.h"
#include "reflections.h"
#include "creflections.h"
#include "sblas.h"
#include "ablasf.h"
#include "ablas.h"


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
void rmatrixqr(ap::real_2d_array& a, int m, int n, ap::real_1d_array& tau);


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
void rmatrixlq(ap::real_2d_array& a, int m, int n, ap::real_1d_array& tau);


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
     ap::complex_1d_array& tau);


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
     ap::complex_1d_array& tau);


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
     ap::real_2d_array& q);


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
     ap::real_2d_array& r);


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
     ap::real_2d_array& q);


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
     ap::real_2d_array& l);


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
     ap::complex_2d_array& q);


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
     ap::complex_2d_array& r);


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
     ap::complex_2d_array& q);


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
     ap::complex_2d_array& l);


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
     ap::real_1d_array& taup);


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
     ap::real_2d_array& q);


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
     bool dotranspose);


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
     ap::real_2d_array& pt);


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
     bool dotranspose);


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
     ap::real_1d_array& e);


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
void rmatrixhessenberg(ap::real_2d_array& a, int n, ap::real_1d_array& tau);


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
     ap::real_2d_array& q);


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
     ap::real_2d_array& h);


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
     ap::real_1d_array& e);


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
     ap::real_2d_array& q);


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
     ap::real_1d_array& e);


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
     ap::complex_2d_array& q);


#endif


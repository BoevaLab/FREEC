/*************************************************************************
Copyright (c) 2007-2008, Sergey Bochkanov (ALGLIB project).

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

#ifndef _linreg_h
#define _linreg_h

#include "ap.h"
#include "ialglib.h"

#include "descriptivestatistics.h"
#include "gammafunc.h"
#include "normaldistr.h"
#include "igammaf.h"
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
#include "svd.h"


struct linearmodel
{
    ap::real_1d_array w;
};


/*************************************************************************
LRReport structure contains additional information about linear model:
* C             -   covariation matrix,  array[0..NVars,0..NVars].
                    C[i,j] = Cov(A[i],A[j])
* RMSError      -   root mean square error on a training set
* AvgError      -   average error on a training set
* AvgRelError   -   average relative error on a training set (excluding
                    observations with zero function value).
* CVRMSError    -   leave-one-out cross-validation estimate of
                    generalization error. Calculated using fast algorithm
                    with O(NVars*NPoints) complexity.
* CVAvgError    -   cross-validation estimate of average error
* CVAvgRelError -   cross-validation estimate of average relative error

All other fields of the structure are intended for internal use and should
not be used outside ALGLIB.
*************************************************************************/
struct lrreport
{
    ap::real_2d_array c;
    double rmserror;
    double avgerror;
    double avgrelerror;
    double cvrmserror;
    double cvavgerror;
    double cvavgrelerror;
    int ncvdefects;
    ap::integer_1d_array cvdefects;
};




/*************************************************************************
Linear regression

Subroutine builds model:

    Y = A(0)*X[0] + ... + A(N-1)*X[N-1] + A(N)

and model found in ALGLIB format, covariation matrix, training set  errors
(rms,  average,  average  relative)   and  leave-one-out  cross-validation
estimate of the generalization error. CV  estimate calculated  using  fast
algorithm with O(NPoints*NVars) complexity.

When  covariation  matrix  is  calculated  standard deviations of function
values are assumed to be equal to RMS error on the training set.

INPUT PARAMETERS:
    XY          -   training set, array [0..NPoints-1,0..NVars]:
                    * NVars columns - independent variables
                    * last column - dependent variable
    NPoints     -   training set size, NPoints>NVars+1
    NVars       -   number of independent variables

OUTPUT PARAMETERS:
    Info        -   return code:
                    * -255, in case of unknown internal error
                    * -4, if internal SVD subroutine haven't converged
                    * -1, if incorrect parameters was passed (NPoints<NVars+2, NVars<1).
                    *  1, if subroutine successfully finished
    LM          -   linear model in the ALGLIB format. Use subroutines of
                    this unit to work with the model.
    AR          -   additional results


  -- ALGLIB --
     Copyright 02.08.2008 by Bochkanov Sergey
*************************************************************************/
void lrbuild(const ap::real_2d_array& xy,
     int npoints,
     int nvars,
     int& info,
     linearmodel& lm,
     lrreport& ar);


/*************************************************************************
Linear regression

Variant of LRBuild which uses vector of standatd deviations (errors in
function values).

INPUT PARAMETERS:
    XY          -   training set, array [0..NPoints-1,0..NVars]:
                    * NVars columns - independent variables
                    * last column - dependent variable
    S           -   standard deviations (errors in function values)
                    array[0..NPoints-1], S[i]>0.
    NPoints     -   training set size, NPoints>NVars+1
    NVars       -   number of independent variables

OUTPUT PARAMETERS:
    Info        -   return code:
                    * -255, in case of unknown internal error
                    * -4, if internal SVD subroutine haven't converged
                    * -1, if incorrect parameters was passed (NPoints<NVars+2, NVars<1).
                    * -2, if S[I]<=0
                    *  1, if subroutine successfully finished
    LM          -   linear model in the ALGLIB format. Use subroutines of
                    this unit to work with the model.
    AR          -   additional results


  -- ALGLIB --
     Copyright 02.08.2008 by Bochkanov Sergey
*************************************************************************/
void lrbuilds(const ap::real_2d_array& xy,
     const ap::real_1d_array& s,
     int npoints,
     int nvars,
     int& info,
     linearmodel& lm,
     lrreport& ar);


/*************************************************************************
Like LRBuildS, but builds model

    Y = A(0)*X[0] + ... + A(N-1)*X[N-1]

i.e. with zero constant term.

  -- ALGLIB --
     Copyright 30.10.2008 by Bochkanov Sergey
*************************************************************************/
void lrbuildzs(const ap::real_2d_array& xy,
     const ap::real_1d_array& s,
     int npoints,
     int nvars,
     int& info,
     linearmodel& lm,
     lrreport& ar);


/*************************************************************************
Like LRBuild but builds model

    Y = A(0)*X[0] + ... + A(N-1)*X[N-1]

i.e. with zero constant term.

  -- ALGLIB --
     Copyright 30.10.2008 by Bochkanov Sergey
*************************************************************************/
void lrbuildz(const ap::real_2d_array& xy,
     int npoints,
     int nvars,
     int& info,
     linearmodel& lm,
     lrreport& ar);


/*************************************************************************
Unpacks coefficients of linear model.

INPUT PARAMETERS:
    LM          -   linear model in ALGLIB format

OUTPUT PARAMETERS:
    V           -   coefficients, array[0..NVars]
    NVars       -   number of independent variables (one less than number
                    of coefficients)

  -- ALGLIB --
     Copyright 30.08.2008 by Bochkanov Sergey
*************************************************************************/
void lrunpack(const linearmodel& lm, ap::real_1d_array& v, int& nvars);


/*************************************************************************
"Packs" coefficients and creates linear model in ALGLIB format (LRUnpack
reversed).

INPUT PARAMETERS:
    V           -   coefficients, array[0..NVars]
    NVars       -   number of independent variables

OUTPUT PAREMETERS:
    LM          -   linear model.

  -- ALGLIB --
     Copyright 30.08.2008 by Bochkanov Sergey
*************************************************************************/
void lrpack(const ap::real_1d_array& v, int nvars, linearmodel& lm);


/*************************************************************************
Procesing

INPUT PARAMETERS:
    LM      -   linear model
    X       -   input vector,  array[0..NVars-1].

Result:
    value of linear model regression estimate

  -- ALGLIB --
     Copyright 03.09.2008 by Bochkanov Sergey
*************************************************************************/
double lrprocess(const linearmodel& lm, const ap::real_1d_array& x);


/*************************************************************************
RMS error on the test set

INPUT PARAMETERS:
    LM      -   linear model
    XY      -   test set
    NPoints -   test set size

RESULT:
    root mean square error.

  -- ALGLIB --
     Copyright 30.08.2008 by Bochkanov Sergey
*************************************************************************/
double lrrmserror(const linearmodel& lm,
     const ap::real_2d_array& xy,
     int npoints);


/*************************************************************************
Average error on the test set

INPUT PARAMETERS:
    LM      -   linear model
    XY      -   test set
    NPoints -   test set size

RESULT:
    average error.

  -- ALGLIB --
     Copyright 30.08.2008 by Bochkanov Sergey
*************************************************************************/
double lravgerror(const linearmodel& lm,
     const ap::real_2d_array& xy,
     int npoints);


/*************************************************************************
RMS error on the test set

INPUT PARAMETERS:
    LM      -   linear model
    XY      -   test set
    NPoints -   test set size

RESULT:
    average relative error.

  -- ALGLIB --
     Copyright 30.08.2008 by Bochkanov Sergey
*************************************************************************/
double lravgrelerror(const linearmodel& lm,
     const ap::real_2d_array& xy,
     int npoints);


/*************************************************************************
Copying of LinearModel strucure

INPUT PARAMETERS:
    LM1 -   original

OUTPUT PARAMETERS:
    LM2 -   copy

  -- ALGLIB --
     Copyright 15.03.2009 by Bochkanov Sergey
*************************************************************************/
void lrcopy(const linearmodel& lm1, linearmodel& lm2);


/*************************************************************************
Serialization of LinearModel strucure

INPUT PARAMETERS:
    LM      -   original

OUTPUT PARAMETERS:
    RA      -   array of real numbers which stores model,
                array[0..RLen-1]
    RLen    -   RA lenght

  -- ALGLIB --
     Copyright 15.03.2009 by Bochkanov Sergey
*************************************************************************/
void lrserialize(const linearmodel& lm, ap::real_1d_array& ra, int& rlen);


/*************************************************************************
Unserialization of DecisionForest strucure

INPUT PARAMETERS:
    RA      -   real array which stores decision forest

OUTPUT PARAMETERS:
    LM      -   unserialized structure

  -- ALGLIB --
     Copyright 15.03.2009 by Bochkanov Sergey
*************************************************************************/
void lrunserialize(const ap::real_1d_array& ra, linearmodel& lm);


void lrlines(const ap::real_2d_array& xy,
     const ap::real_1d_array& s,
     int n,
     int& info,
     double& a,
     double& b,
     double& vara,
     double& varb,
     double& covab,
     double& corrab,
     double& p);


void lrline(const ap::real_2d_array& xy,
     int n,
     int& info,
     double& a,
     double& b);


#endif


/*************************************************************************
Copyright (c) 2007, Sergey Bochkanov (ALGLIB project).

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

#ifndef _descriptivestatistics_h
#define _descriptivestatistics_h

#include "ap.h"
#include "ialglib.h"

/*************************************************************************
Calculation of the distribution moments: mean, variance, slewness, kurtosis.

Input parameters:
    X       -   sample. Array with whose indexes range within [0..N-1]
    N       -   sample size.
    
Output parameters:
    Mean    -   mean.
    Variance-   variance.
    Skewness-   skewness (if variance<>0; zero otherwise).
    Kurtosis-   kurtosis (if variance<>0; zero otherwise).

  -- ALGLIB --
     Copyright 06.09.2006 by Bochkanov Sergey
*************************************************************************/
void calculatemoments(const ap::real_1d_array& x,
     int n,
     double& mean,
     double& variance,
     double& skewness,
     double& kurtosis);


/*************************************************************************
ADev

Input parameters:
    X   -   sample (array indexes: [0..N-1])
    N   -   sample size
    
Output parameters:
    ADev-   ADev

  -- ALGLIB --
     Copyright 06.09.2006 by Bochkanov Sergey
*************************************************************************/
void calculateadev(const ap::real_1d_array& x, int n, double& adev);


/*************************************************************************
Median calculation.

Input parameters:
    X   -   sample (array indexes: [0..N-1])
    N   -   sample size

Output parameters:
    Median

  -- ALGLIB --
     Copyright 06.09.2006 by Bochkanov Sergey
*************************************************************************/
void calculatemedian(ap::real_1d_array x, int n, double& median);


/*************************************************************************
Percentile calculation.

Input parameters:
    X   -   sample (array indexes: [0..N-1])
    N   -   sample size, N>1
    P   -   percentile (0<=P<=1)

Output parameters:
    V   -   percentile

  -- ALGLIB --
     Copyright 01.03.2008 by Bochkanov Sergey
*************************************************************************/
void calculatepercentile(ap::real_1d_array x, int n, double p, double& v);


#endif


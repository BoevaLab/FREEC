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

////#include <stdafx.h>
#include "descriptivestatistics.h"

static void internalstatheapsort(ap::real_1d_array& arr, int n);

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
     double& kurtosis)
{
    int i;
    double v;
    double v1;
    double v2;
    double stddev;

    mean = 0;
    variance = 0;
    skewness = 0;
    kurtosis = 0;
    stddev = 0;
    if( n<=0 )
    {
        return;
    }
    
    //
    // Mean
    //
    for(i = 0; i <= n-1; i++)
    {
        mean = mean+x(i);
    }
    mean = mean/n;
    
    //
    // Variance (using corrected two-pass algorithm)
    //
    if( n!=1 )
    {
        v1 = 0;
        for(i = 0; i <= n-1; i++)
        {
            v1 = v1+ap::sqr(x(i)-mean);
        }
        v2 = 0;
        for(i = 0; i <= n-1; i++)
        {
            v2 = v2+(x(i)-mean);
        }
        v2 = ap::sqr(v2)/n;
        variance = (v1-v2)/(n-1);
        if( ap::fp_less(variance,0) )
        {
            variance = 0;
        }
        stddev = sqrt(variance);
    }
    
    //
    // Skewness and kurtosis
    //
    if( ap::fp_neq(stddev,0) )
    {
        for(i = 0; i <= n-1; i++)
        {
            v = (x(i)-mean)/stddev;
            v2 = ap::sqr(v);
            skewness = skewness+v2*v;
            kurtosis = kurtosis+ap::sqr(v2);
        }
        skewness = skewness/n;
        kurtosis = kurtosis/n-3;
    }
}


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
void calculateadev(const ap::real_1d_array& x, int n, double& adev)
{
    int i;
    double mean;

    mean = 0;
    adev = 0;
    if( n<=0 )
    {
        return;
    }
    
    //
    // Mean
    //
    for(i = 0; i <= n-1; i++)
    {
        mean = mean+x(i);
    }
    mean = mean/n;
    
    //
    // ADev
    //
    for(i = 0; i <= n-1; i++)
    {
        adev = adev+fabs(x(i)-mean);
    }
    adev = adev/n;
}


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
void calculatemedian(ap::real_1d_array x, int n, double& median)
{
    int i;
    int ir;
    int j;
    int l;
    int midp;
    int k;
    double a;
    double tval;

    
    //
    // Some degenerate cases
    //
    median = 0;
    if( n<=0 )
    {
        return;
    }
    if( n==1 )
    {
        median = x(0);
        return;
    }
    if( n==2 )
    {
        median = 0.5*(x(0)+x(1));
        return;
    }
    
    //
    // Common case, N>=3.
    // Choose X[(N-1)/2]
    //
    l = 0;
    ir = n-1;
    k = (n-1)/2;
    while(true)
    {
        if( ir<=l+1 )
        {
            
            //
            // 1 or 2 elements in partition
            //
            if( ir==l+1&&ap::fp_less(x(ir),x(l)) )
            {
                tval = x(l);
                x(l) = x(ir);
                x(ir) = tval;
            }
            break;
        }
        else
        {
            midp = (l+ir)/2;
            tval = x(midp);
            x(midp) = x(l+1);
            x(l+1) = tval;
            if( ap::fp_greater(x(l),x(ir)) )
            {
                tval = x(l);
                x(l) = x(ir);
                x(ir) = tval;
            }
            if( ap::fp_greater(x(l+1),x(ir)) )
            {
                tval = x(l+1);
                x(l+1) = x(ir);
                x(ir) = tval;
            }
            if( ap::fp_greater(x(l),x(l+1)) )
            {
                tval = x(l);
                x(l) = x(l+1);
                x(l+1) = tval;
            }
            i = l+1;
            j = ir;
            a = x(l+1);
            while(true)
            {
                do
                {
                    i = i+1;
                }
                while(ap::fp_less(x(i),a));
                do
                {
                    j = j-1;
                }
                while(ap::fp_greater(x(j),a));
                if( j<i )
                {
                    break;
                }
                tval = x(i);
                x(i) = x(j);
                x(j) = tval;
            }
            x(l+1) = x(j);
            x(j) = a;
            if( j>=k )
            {
                ir = j-1;
            }
            if( j<=k )
            {
                l = i;
            }
        }
    }
    
    //
    // If N is odd, return result
    //
    if( n%2==1 )
    {
        median = x(k);
        return;
    }
    a = x(n-1);
    for(i = k+1; i <= n-1; i++)
    {
        if( ap::fp_less(x(i),a) )
        {
            a = x(i);
        }
    }
    median = 0.5*(x(k)+a);
}


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
void calculatepercentile(ap::real_1d_array x, int n, double p, double& v)
{
    int i1;
    double t;

    ap::ap_error::make_assertion(n>1, "CalculatePercentile: N<=1!");
    ap::ap_error::make_assertion(ap::fp_greater_eq(p,0)&&ap::fp_less_eq(p,1), "CalculatePercentile: incorrect P!");
    internalstatheapsort(x, n);
    if( ap::fp_eq(p,0) )
    {
        v = x(0);
        return;
    }
    if( ap::fp_eq(p,1) )
    {
        v = x(n-1);
        return;
    }
    t = p*(n-1);
    i1 = ap::ifloor(t);
    t = t-ap::ifloor(t);
    v = x(i1)*(1-t)+x(i1+1)*t;
}


static void internalstatheapsort(ap::real_1d_array& arr, int n)
{
    int i;
    int k;
    int t;
    double tmp;

    if( n==1 )
    {
        return;
    }
    i = 2;
    do
    {
        t = i;
        while(t!=1)
        {
            k = t/2;
            if( ap::fp_greater_eq(arr(k-1),arr(t-1)) )
            {
                t = 1;
            }
            else
            {
                tmp = arr(k-1);
                arr(k-1) = arr(t-1);
                arr(t-1) = tmp;
                t = k;
            }
        }
        i = i+1;
    }
    while(i<=n);
    i = n-1;
    do
    {
        tmp = arr(i);
        arr(i) = arr(0);
        arr(0) = tmp;
        t = 1;
        while(t!=0)
        {
            k = 2*t;
            if( k>i )
            {
                t = 0;
            }
            else
            {
                if( k<i )
                {
                    if( ap::fp_greater(arr(k),arr(k-1)) )
                    {
                        k = k+1;
                    }
                }
                if( ap::fp_greater_eq(arr(t-1),arr(k-1)) )
                {
                    t = 0;
                }
                else
                {
                    tmp = arr(k-1);
                    arr(k-1) = arr(t-1);
                    arr(t-1) = tmp;
                    t = k;
                }
            }
        }
        i = i-1;
    }
    while(i>=1);
}





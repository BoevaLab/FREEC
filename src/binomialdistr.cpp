/*************************************************************************
Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1995, 2000 by Stephen L. Moshier

Contributors:
    * Sergey Bochkanov (ALGLIB project). Translation from C to
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
#include "binomialdistr.h"

/*************************************************************************
Binomial distribution

Returns the sum of the terms 0 through k of the Binomial
probability density:

  k
  --  ( n )   j      n-j
  >   (   )  p  (1-p)
  --  ( j )
 j=0

The terms are not summed directly; instead the incomplete
beta integral is employed, according to the formula

y = bdtr( k, n, p ) = incbet( n-k, k+1, 1-p ).

The arguments must be positive, with p ranging from 0 to 1.

ACCURACY:

Tested at random points (a,b,p), with p between 0 and 1.

              a,b                     Relative error:
arithmetic  domain     # trials      peak         rms
 For p between 0.001 and 1:
   IEEE     0,100       100000      4.3e-15     2.6e-16

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1995, 2000 by Stephen L. Moshier
*************************************************************************/
double binomialdistribution(int k, int n, double p)
{
    double result;
    double dk;
    double dn;

    ap::ap_error::make_assertion(ap::fp_greater_eq(p,0)&&ap::fp_less_eq(p,1), "Domain error in BinomialDistribution");
    ap::ap_error::make_assertion(k>=-1&&k<=n, "Domain error in BinomialDistribution");
    if( k==-1 )
    {
        result = 0;
        return result;
    }
    if( k==n )
    {
        result = 1;
        return result;
    }
    dn = n-k;
    if( k==0 )
    {
        dk = pow(1.0-p, dn);
    }
    else
    {
        dk = k+1;
        dk = incompletebeta(dn, dk, 1.0-p);
    }
    result = dk;
    return result;
}


/*************************************************************************
Complemented binomial distribution

Returns the sum of the terms k+1 through n of the Binomial
probability density:

  n
  --  ( n )   j      n-j
  >   (   )  p  (1-p)
  --  ( j )
 j=k+1

The terms are not summed directly; instead the incomplete
beta integral is employed, according to the formula

y = bdtrc( k, n, p ) = incbet( k+1, n-k, p ).

The arguments must be positive, with p ranging from 0 to 1.

ACCURACY:

Tested at random points (a,b,p).

              a,b                     Relative error:
arithmetic  domain     # trials      peak         rms
 For p between 0.001 and 1:
   IEEE     0,100       100000      6.7e-15     8.2e-16
 For p between 0 and .001:
   IEEE     0,100       100000      1.5e-13     2.7e-15

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1995, 2000 by Stephen L. Moshier
*************************************************************************/
double binomialcdistribution(int k, int n, double p)
{
    double result;
    double dk;
    double dn;
    //std::cout << k << " of "<< n<<"\n";
    ap::ap_error::make_assertion(ap::fp_greater_eq(p,0)&&ap::fp_less_eq(p,1), "Domain error in BinomialDistributionC");
    ap::ap_error::make_assertion(k>=-1&&k<=n, "Domain error in BinomialDistributionC");
    if( k==-1 )
    {
        result = 1;
        return result;
    }
    if( k==n )
    {
        result = 0;
        return result;
    }
    dn = n-k;
    if( k==0 )
    {
        if( ap::fp_less(p,0.01) )
        {
            dk = -expm1(dn*log1p(-p));
        }
        else
        {
            dk = 1.0-pow(1.0-p, dn);
        }
    }
    else
    {
        dk = k+1;
        dk = incompletebeta(dk, dn, p);
    }
    result = dk;
    return result;
}


/*************************************************************************
Inverse binomial distribution

Finds the event probability p such that the sum of the
terms 0 through k of the Binomial probability density
is equal to the given cumulative probability y.

This is accomplished using the inverse beta integral
function and the relation

1 - p = incbi( n-k, k+1, y ).

ACCURACY:

Tested at random points (a,b,p).

              a,b                     Relative error:
arithmetic  domain     # trials      peak         rms
 For p between 0.001 and 1:
   IEEE     0,100       100000      2.3e-14     6.4e-16
   IEEE     0,10000     100000      6.6e-12     1.2e-13
 For p between 10^-6 and 0.001:
   IEEE     0,100       100000      2.0e-12     1.3e-14
   IEEE     0,10000     100000      1.5e-12     3.2e-14

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1995, 2000 by Stephen L. Moshier
*************************************************************************/
double invbinomialdistribution(int k, int n, double y)
{
    double result;
    double dk;
    double dn;
    double p;

    ap::ap_error::make_assertion(k>=0&&k<n, "Domain error in InvBinomialDistribution");
    dn = n-k;
    if( k==0 )
    {
        if( ap::fp_greater(y,0.8) )
        {
            p = -expm1(log1p(y-1.0)/dn);
        }
        else
        {
            p = 1.0-pow(y, 1.0/dn);
        }
    }
    else
    {
        dk = k+1;
        p = incompletebeta(dn, dk, 0.5);
        if( ap::fp_greater(p,0.5) )
        {
            p = invincompletebeta(dk, dn, 1.0-y);
        }
        else
        {
            p = 1.0-invincompletebeta(dn, dk, y);
        }
    }
    result = p;
    return result;
}





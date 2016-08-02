#ifndef HEADER_362C1E424D7DDEE3
#define HEADER_362C1E424D7DDEE3

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

#ifndef _binomialdistr_h
#define _binomialdistr_h

#include "ap.h"
#include "ialglib.h"

#include "gammafunc.h"
#include "normaldistr.h"
#include "ibetaf.h"
#include <iostream>
#include <cmath>
#include "myFunc.h"


//#include "nearunityunit.h"


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
double binomialdistribution(int k, int n, double p);


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
double binomialcdistribution(int k, int n, double p);


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
double invbinomialdistribution(int k, int n, double y);


#endif


#endif // header guard

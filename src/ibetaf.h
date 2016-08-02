/*************************************************************************
Cephes Math Library Release 2.8:  June, 2000
Copyright by Stephen L. Moshier

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

#ifndef _ibetaf_h
#define _ibetaf_h

#include "ap.h"
#include "ialglib.h"

#include "gammafunc.h"
#include "normaldistr.h"


/*************************************************************************
Incomplete beta integral

Returns incomplete beta integral of the arguments, evaluated
from zero to x.  The function is defined as

                 x
    -            -
   | (a+b)      | |  a-1     b-1
 -----------    |   t   (1-t)   dt.
  -     -     | |
 | (a) | (b)   -
                0

The domain of definition is 0 <= x <= 1.  In this
implementation a and b are restricted to positive values.
The integral from x to 1 may be obtained by the symmetry
relation

   1 - incbet( a, b, x )  =  incbet( b, a, 1-x ).

The integral is evaluated by a continued fraction expansion
or, when b*x is small, by a power series.

ACCURACY:

Tested at uniformly distributed random points (a,b,x) with a and b
in "domain" and x between 0 and 1.
                                       Relative error
arithmetic   domain     # trials      peak         rms
   IEEE      0,5         10000       6.9e-15     4.5e-16
   IEEE      0,85       250000       2.2e-13     1.7e-14
   IEEE      0,1000      30000       5.3e-12     6.3e-13
   IEEE      0,10000    250000       9.3e-11     7.1e-12
   IEEE      0,100000    10000       8.7e-10     4.8e-11
Outputs smaller than the IEEE gradual underflow threshold
were excluded from these statistics.

Cephes Math Library, Release 2.8:  June, 2000
Copyright 1984, 1995, 2000 by Stephen L. Moshier
*************************************************************************/
double incompletebeta(double a, double b, double x);


/*************************************************************************
Inverse of imcomplete beta integral

Given y, the function finds x such that

 incbet( a, b, x ) = y .

The routine performs interval halving or Newton iterations to find the
root of incbet(a,b,x) - y = 0.


ACCURACY:

                     Relative error:
               x     a,b
arithmetic   domain  domain  # trials    peak       rms
   IEEE      0,1    .5,10000   50000    5.8e-12   1.3e-13
   IEEE      0,1   .25,100    100000    1.8e-13   3.9e-15
   IEEE      0,1     0,5       50000    1.1e-12   5.5e-15
With a and b constrained to half-integer or integer values:
   IEEE      0,1    .5,10000   50000    5.8e-12   1.1e-13
   IEEE      0,1    .5,100    100000    1.7e-14   7.9e-16
With a = .5, b constrained to half-integer or integer values:
   IEEE      0,1    .5,10000   10000    8.3e-11   1.0e-11

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1996, 2000 by Stephen L. Moshier
*************************************************************************/
double invincompletebeta(double a, double b, double y);


#endif


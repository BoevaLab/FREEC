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

#ifndef _igammaf_h
#define _igammaf_h

#include "ap.h"
#include "ialglib.h"

#include "gammafunc.h"
#include "normaldistr.h"


/*************************************************************************
Incomplete gamma integral

The function is defined by

                          x
                           -
                  1       | |  -t  a-1
 igam(a,x)  =   -----     |   e   t   dt.
                 -      | |
                | (a)    -
                          0


In this implementation both arguments must be positive.
The integral is evaluated by either a power series or
continued fraction expansion, depending on the relative
values of a and x.

ACCURACY:

                     Relative error:
arithmetic   domain     # trials      peak         rms
   IEEE      0,30       200000       3.6e-14     2.9e-15
   IEEE      0,100      300000       9.9e-14     1.5e-14

Cephes Math Library Release 2.8:  June, 2000
Copyright 1985, 1987, 2000 by Stephen L. Moshier
*************************************************************************/
double incompletegamma(double a, double x);


/*************************************************************************
Complemented incomplete gamma integral

The function is defined by


 igamc(a,x)   =   1 - igam(a,x)

                           inf.
                             -
                    1       | |  -t  a-1
              =   -----     |   e   t   dt.
                   -      | |
                  | (a)    -
                            x


In this implementation both arguments must be positive.
The integral is evaluated by either a power series or
continued fraction expansion, depending on the relative
values of a and x.

ACCURACY:

Tested at random a, x.
               a         x                      Relative error:
arithmetic   domain   domain     # trials      peak         rms
   IEEE     0.5,100   0,100      200000       1.9e-14     1.7e-15
   IEEE     0.01,0.5  0,100      200000       1.4e-13     1.6e-15

Cephes Math Library Release 2.8:  June, 2000
Copyright 1985, 1987, 2000 by Stephen L. Moshier
*************************************************************************/
double incompletegammac(double a, double x);


/*************************************************************************
Inverse of complemented imcomplete gamma integral

Given p, the function finds x such that

 igamc( a, x ) = p.

Starting with the approximate value

        3
 x = a t

 where

 t = 1 - d - ndtri(p) sqrt(d)

and

 d = 1/9a,

the routine performs up to 10 Newton iterations to find the
root of igamc(a,x) - p = 0.

ACCURACY:

Tested at random a, p in the intervals indicated.

               a        p                      Relative error:
arithmetic   domain   domain     # trials      peak         rms
   IEEE     0.5,100   0,0.5       100000       1.0e-14     1.7e-15
   IEEE     0.01,0.5  0,0.5       100000       9.0e-14     3.4e-15
   IEEE    0.5,10000  0,0.5        20000       2.3e-13     3.8e-14

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1995, 2000 by Stephen L. Moshier
*************************************************************************/
double invincompletegammac(double a, double y0);


#endif


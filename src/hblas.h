/*************************************************************************
Copyright (c) 1992-2007 The University of Tennessee.  All rights reserved.

Contributors:
    * Sergey Bochkanov (ALGLIB project). Translation from FORTRAN to
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

#ifndef _hblas_h
#define _hblas_h

#include "ap.h"
#include "ialglib.h"

void hermitianmatrixvectormultiply(const ap::complex_2d_array& a,
     bool isupper,
     int i1,
     int i2,
     const ap::complex_1d_array& x,
     ap::complex alpha,
     ap::complex_1d_array& y);


void hermitianrank2update(ap::complex_2d_array& a,
     bool isupper,
     int i1,
     int i2,
     const ap::complex_1d_array& x,
     const ap::complex_1d_array& y,
     ap::complex_1d_array& t,
     ap::complex alpha);


#endif


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

#ifndef _sblas_h
#define _sblas_h

#include "ap.h"
#include "ialglib.h"

void symmetricmatrixvectormultiply(const ap::real_2d_array& a,
     bool isupper,
     int i1,
     int i2,
     const ap::real_1d_array& x,
     double alpha,
     ap::real_1d_array& y);


void symmetricrank2update(ap::real_2d_array& a,
     bool isupper,
     int i1,
     int i2,
     const ap::real_1d_array& x,
     const ap::real_1d_array& y,
     ap::real_1d_array& t,
     double alpha);


#endif


/*************************************************************************
Copyright (c) 2010-2011, Valentina BOEVA.

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


#include "KernelVector.h"


KernelVector::KernelVector(int r)  {
	radius_ = r;
	std::vector<float> V(r);
	sum_ = 0;
	for (int i = 0; i<r; i++) {
		V[i] = EpanechnikovKernel(float(i)/r);
		//std::cout << i << " " << V[i]<< "\n";
		sum_ += 2*V[i];
	}
	sum_ -= EpanechnikovKernel(0);
	std::cout << sum_ << "\n";
	vector_ = V;
	
}

float KernelVector::getK(int i , int j) {
	int l = abs(i-j);
	if (l<radius_)
		return vector_[l];
	return 0;
}



float KernelVector::getSum() {
	return sum_;
}

int KernelVector::getRadius() {
	return radius_;
}

KernelVector::~KernelVector(void) {
	vector_.clear();
	radius_ = 0;
}

float EpanechnikovKernel (float u) {
	if (-1<=u && u<=1) return 3*(1-u*u)/4 ;
	return 0;
}

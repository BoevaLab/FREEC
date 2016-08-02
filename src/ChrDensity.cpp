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


#include "ChrDensity.h"

using namespace std ;

ChrDensity::ChrDensity(int length)
{
	length_ = length;	
	coverage_ = vector <int>(length,0); 	
	
	//coverage_.resize(length);
}

int ChrDensity::getLength() {
	return this->length_;
}

int ChrDensity::getCoverageAtI(int i) {
	if (i>=length_)
		return 0;

	if (i<0)
		return -1;
	return coverage_[i];
}
float ChrDensity::getDensityAtI(int i) {
	if (i>=length_)
		return 0;

	if (i<0)
		return -1;
	return density_[i];
}

void ChrDensity::coveragePlusOneAtI(int i) {
	if (i>=0 && i<length_)
		coverage_[i]++; 
}

void ChrDensity::coveragePlusOneAtI(int start, int end) {
	if (start>=0 && end < length_)
		for (int j = start; j<end; j++)
			coverage_[j]++; 
}

void ChrDensity::smooth(KernelVector kernelVector) {	
	
	vector <float> third (length_,0); 
	density_ = third;
	//KernelVector kernelVector(r); 
	int r = kernelVector.getRadius();

	//left point 
	for (int i = 0; i<r-1; i++) {
		float sumKY = 0;
		float sumK = 0;
		for (int j = 0; j<i+r; j++) {
			sumKY+=kernelVector.getK(j,i)*coverage_[j];		
			sumK +=kernelVector.getK(j,i);
		}
		density_[i]=sumKY/sumK;
	}

	//right point
	for (int i = length_-r+1; i<length_; i++) {
		float sumKY = 0;
		float sumK = 0;
		for (int j = i-r+1; j<length_; j++) {
			sumKY+=kernelVector.getK(j,i)*coverage_[j];		
			sumK +=kernelVector.getK(j,i);
		}		
		density_[i]=sumKY/sumK;
	}

	//all other points
	for (int i = r-1; i<length_-r+1; i++) {
		float sumKY = 0;
		float sumK = kernelVector.getSum();
		for (int j = i-r+1; j<i+r; j++) 
			sumKY+=kernelVector.getK(j,i)*coverage_[j];	
		density_[i]=sumKY/sumK;

		if (i % 1000000 == 0)
			cout << i << "\t" << density_[i] << "\t" << coverage_[i] << "\n";
	}
}

void ChrDensity::calculateLogRatio(ChrDensity controlCD, std::string type) {
	typeOfLogRatio_ = type;
	if (length_!=controlCD.getLength())
		return; //throw error TODO
	std::vector<float> second(length_);
	logRatio = second;
	if (type.compare("norm")==0) {	
		for (int i = 0; i<length_; i++) {
			int y = controlCD.getCoverageAtI(i);
			if (y!=0)
				logRatio[i] = float(coverage_[i])/y;
			else 
				logRatio[i] = -1;
		
			if (i % 50000 == 0) {
			//	std::cout << i << "\t" << logRatio[i] << "\t" << coverage_[i] << "\t" << y << "\n";
				std::cout << i << "\t" << exp(-(float)y)*pow((float)y,coverage_[i])/factorial(coverage_[i]) << "\t" << -log(exp((float)y)*pow((float)y,coverage_[i])/factorial(coverage_[i])) << "\t" << coverage_[i] << "\t" << y << "\n";
			}			
		}
	}
	else {	
		for (int i = 0; i<length_; i++) {
			float y = controlCD.getDensityAtI(i);
			if (y!=0)
				logRatio[i] = float(density_[i])/y;
			else 
				logRatio[i] = -1;
		
			//if (i % 50000 == 0)
			//	std::cout << i << "\t" << logRatio[i] << "\t" << density_[i] << "\t" << y << "\n";
			
		}
	}

}

ChrDensity::~ChrDensity(void)
{
	coverage_.clear();
	density_.clear();
	length_ = 0;
}


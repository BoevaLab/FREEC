#pragma once
#ifndef _CHR_DENSITY_H
#define _CHR_DENSITY_H

#include <iostream>
#include <vector>

#include "myFunc.h"
#include "KernelVector.h"
//#include "poissondistr.h"

class ChrDensity
{
public:
	ChrDensity(int length);
	~ChrDensity(void);
	int getLength();
	int getCoverageAtI(int i);
	float getDensityAtI(int i);
	void coveragePlusOneAtI(int i) ;
	void coveragePlusOneAtI(int start, int end);
	void smooth(KernelVector kernelVector) ;
	void calculateLogRatio(ChrDensity controlCD, std::string type);
private:
	int length_;
	std::vector<int> coverage_;
	std::vector<float> density_;
	std::vector<float> logRatio;
	std::string typeOfLogRatio_;
};
#endif 
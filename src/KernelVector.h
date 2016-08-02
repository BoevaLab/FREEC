#pragma once
#ifndef __KERNEL_VECTOR_H__
#define __KERNEL_VECTOR_H__

#include <vector>
#include <stdlib.h>
#include <iostream>

class KernelVector
{
public:
	KernelVector(int radius);
	~KernelVector(void);
	float getK(int i , int j);
	float getSum();
	int getRadius();
private:
	std::vector<float> vector_;
	int radius_;
	float sum_;
};

float EpanechnikovKernel (float u);

#endif 
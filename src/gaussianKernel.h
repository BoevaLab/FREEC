//
// Created by gajaques on 6/1/22.
//

#ifndef SRC_GAUSSIANKERNEL_H
#define SRC_GAUSSIANKERNEL_H

#include <vector>

class gaussianKernel
{
public:
    gaussianKernel();
    ~gaussianKernel();
};

std::vector<float> computeGaussianKernel(const std::vector<float>& x, float x_value, float cov);
float applyGaussianKernelToValues(const std::vector<float>& kernel, const std::vector<float>& interval_x_values);
std::vector<float> applyGaussianKernel(const std::vector<float>& x, float cov, int size = 4);



#endif //SRC_GAUSSIANKERNEL_H

//
// Created by gajaques on 6/1/22.
//

#include "gaussianKernel.h"
#include <iostream>
#include <cmath>

gaussianKernel::gaussianKernel() {}
gaussianKernel::~gaussianKernel() {}

std::vector<float> computeGaussianKernel(const std::vector<float>& x, float x_value, float cov)
{
    int n = x.size();
    std::vector<float> kernel(n);
    float sum_kernel = 0;
    for (int i = 0; i < n; i++) {
        kernel[i] = exp(-pow(x[i] - x_value, 2)/(2*pow(cov,2)));
        sum_kernel += kernel[i];
    }
    for (int i = 0; i < n; i++) {
        kernel[i] = kernel[i]/sum_kernel;
    }
    return kernel;
}

float applyGaussianKernelToValues(const std::vector<float>& kernel, const std::vector<float>& interval_x_values)
{
    if (kernel.size() != interval_x_values.size()) {
        std::cerr << "Error: the kernel needs to have the same size as the vector passed by argument" << std::endl;
    }
    float x_filtered = 0;
    for (unsigned int i = 0; i < interval_x_values.size(); i++) {
        x_filtered += kernel[i] * interval_x_values[i];
    }
    return x_filtered;
}

std::vector<float> applyGaussianKernel(const std::vector<float>& x, float cov, int size)
{
    int n = x.size();
    // Round size to the nearest odd
    int kernel_width = 1;
    if (floor(size*cov) != 0) {
        kernel_width = floor(size * cov);
    } else {
        kernel_width = ceil(size * cov);
    }

    int hkernel_width = floor(kernel_width/2);
    std::vector<float> interval(kernel_width);
    std::vector<float> kernel(kernel_width);
    std::vector<float> x_filtered(n);
    for (int i = 0; i < n; i++) {
        if ((i >= hkernel_width) & ((n-i) > hkernel_width)) {
            interval = {x.begin() + i - 2, x.begin() + i + 2};
        } else if ((i < hkernel_width) & ((n-i) > hkernel_width)) {
            interval = {x.begin(), x.begin() + i + 2};
        } else if ((i >= hkernel_width) & ((n-i) <= hkernel_width)) {
            interval = {x.begin() + i - 2, x.end()};
        } else {
            interval = x;
        }
        kernel = computeGaussianKernel(interval, x[i], cov);
        x_filtered[i] = applyGaussianKernelToValues(kernel, interval);
    }
    return x_filtered;
}
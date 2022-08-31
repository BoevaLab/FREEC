//
// Created by Garance Jaques on 31.05.22.
//

#ifndef SRC_GENOMEGMM_H
#define SRC_GENOMEGMM_H

#include <string>
#include <vector>
#include <Eigen/Dense>
#include "ap.h"
#include "MixtureComponent.h"


class GenomeGMM {
public:
    GenomeGMM();
    ~GenomeGMM();
    
    // Methods for shifting expected MAF
    float estimateSigmaAB(const std::vector<float>& ratio, const std::vector<float>& maf, int ploidy, float purity);
    float estimateSigmaAABB(const std::vector<float>& ratio, const std::vector<float>& maf, int ploidy, float purity);
    float calculateSigma(int ploidy, float purity, int i, float sigma2 = -1, float sigma4 = -1);
    Eigen::MatrixXd subclonalClusters(float s, const std::vector<float> &ratio, const std::vector<float> &maf, float alpha, float tau);
    void correctMAF(std::vector<float> sigma);
    void shiftMAF(const std::vector<float>& ratio, const std::vector<float>& maf_smoothed, int ploidy, float purity);
    // For a given pair of ploidy and purity, the means of ratio and MAF are fixed
    void calculateExpectedRatioAndMAF(int ploidy, float purity);
    // Print the cluster centres
    void printExpectedRatioAndMAF(std::string const& outFile, int ploidy, float purity);
    
    // Helper function that transforms std::vector into Eigen::Vector so that I can use matrices
    // in the GMM fitting part without having to transform them every time
    Eigen::VectorXd transformVectorToMatrix(const std::vector<double>& v);
    
    // Methods for model fitting
    void initialiseGMMParameters(int n_init = 1, int max_iter = 600, double tol = 5.0e-3, double reg_cov = 1.0e-7, double noise = 1, bool b = false, double lbr = 0.0, double ubr = 0.01, double lbm = 0.0, double ubm = 0.015);
    Eigen::Matrix2d initialiseCovariance();
    double initialiseWeight();
    std::vector<double> calculateMarginalResponsabilities(int n);
    double computeLoss(int n);
    double fitPredict(int tau, float alpha, const std::vector<float>& ratio, const std::vector<float>& maf);
    
    
    // Get methods
    std::vector<MixtureComponent> getMixtures();
    double getFinalL();
    int getNumberOfIterations();
    
    // Set methods
    void setNinit(int n);
    void setMaxIter(int m);
    void setTolerance(double t);
    void setRegCovar(double c);
    void setNoise(double n);
    void setBound(bool b);
    void setBoundsRatio(double l, double u);
    void setBoundsMAF(double l, double u);

private:
    int length_expectation_;
    std::vector<float> expected_ratio_;
    std::vector<float> expected_maf_;
    std::vector<float> shifted_expected_maf_;
    bool shifted_MAF_;
    int K_;
    double finalL_;
    int number_of_iterations_;
    
    // GMM parameters
    int n_init_;
    int max_iter_;
    double tolerance_;
    double reg_covar_;
    double noise_;
    bool bound_;
    double lower_bound_ratio_;
    double upper_bound_ratio_;
    double lower_bound_maf_;
    double upper_bound_maf_;
    
    // GMM
    std::vector<MixtureComponent> mixture_;
};


#endif //SRC_GENOMEGMM_H

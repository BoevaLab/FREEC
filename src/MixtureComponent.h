//
// Created by Garance Jaques on 17.06.22.
//

#ifndef SRC_MIXTURECOMPONENT_H
#define SRC_MIXTURECOMPONENT_H

#include <string>
#include <vector>
#include <Eigen/Dense>
#include "myFunc.h"

class MixtureComponent {
public:
    MixtureComponent();
    MixtureComponent(const Eigen::VectorXd& ratio, const Eigen::VectorXd& maf, const Eigen::Vector2d& mu, float w, const Eigen::Matrix2d& cov, int N);
    ~MixtureComponent();
    
    // EM methods
    Eigen::MatrixXd computeShift(const Eigen::VectorXd& ratio, const Eigen::VectorXd& maf);
    void calculateSampleLikelihoods(int n);
    void calculateComponentResponsabilities(const std::vector<double>& marginal_responsabilities);
    double calculateTotalResponsability();
    void updateWeight(double total_resp, int n);
    void updateCovariance(bool bound, double total_resp, int n, double lbr, double ubr, double lbm, double ubm);
    void updateCovarianceAndWeight(bool bound, double lbr, double ubr, double lbm, double ubm);
    
    // Get methods
    double getMixtureComponentWeight();
    Eigen::Matrix2d getMixtureComponentCov();
    //std::vector<float> getMixtureComponentLikelihood();
    double getMixtureComponentLikelihoodAt(int i);
    
    // Set methods
    void setMixtureComponentCovariance(const Eigen::Matrix2d& cov);

    Eigen::Vector2d getMixtureComponentMu();

private:
	Eigen::Vector2d mu_;
    Eigen::MatrixXd shift_;
	Eigen::Matrix2d covariance_;
    double weight_;
    std::vector<double> likelihoods_;
    std::vector<double> responsabilities_;

};


#endif //SRC_MIXTURECOMPONENT_H

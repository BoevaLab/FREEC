//
// Created by Garance Jaques on 17.06.22.
//

#include "MixtureComponent.h"
#include <iostream>

MixtureComponent::MixtureComponent() {}
MixtureComponent::MixtureComponent(const Eigen::VectorXd& ratio, const Eigen::VectorXd& maf, const Eigen::Vector2d& mu, float w, const Eigen::Matrix2d& cov, int N) {
	mu_ = mu;
	weight_ = w;
	covariance_ = cov;
	shift_.resize(ratio.size(), 2);
	shift_ << computeShift(ratio, maf);
	responsabilities_.resize(N);
}
MixtureComponent::~MixtureComponent() {}


// EM methods
Eigen::MatrixXd MixtureComponent::computeShift(const Eigen::VectorXd& ratio, const Eigen::VectorXd& maf)
{
	int n = ratio.size();
	Eigen::MatrixXd shift(n, 2);
	double r, m = 0.0;
	if (n != maf.size()) {
		std::cerr << "Error: ratio and MAF vectors must have the same size" << std::endl;
	} else {
		for (int i = 0; i < n; i++) {
			r = ratio(i) - mu_(0);
			m = maf(i) - mu_(1);
			shift.row(i) << r, m;
		}
	}
    return shift;
}

void MixtureComponent::calculateSampleLikelihoods(int n) {
	int D = 2;
	likelihoods_.resize(n);
    double det = covariance_(0,0)*covariance_(1,1);
    double denom = std::pow(std::sqrt(2*PI_), -D/2.0) * std::pow(det, -0.5);
    double inv11 = covariance_(1,1)/det;
    double inv22 = covariance_(0,0)/det;
    Eigen::Matrix2d inverse;
    inverse << inv11, 0.0, 0.0, inv22;
	for (int i = 0; i < n; i++) {
        Eigen::RowVector2d shift;
        shift << shift_.row(i);
		double exp_term = -0.5 * shift * inverse * shift.transpose();
		likelihoods_[i] = log(denom) + exp_term;
	}
}

void MixtureComponent::calculateComponentResponsabilities(const std::vector<double>& marginal_responsabilities)
{
	int n = marginal_responsabilities.size();
	if (marginal_responsabilities.size() != n) {
        std::cerr << "Error: the length of marginal responsabilities vector needs to be the same as the length of the likelihoods vector" << std::endl;
    }
    for (int i = 0; i < n; i++) {
        double joint = log(weight_) + likelihoods_[i];
	responsabilities_[i] = exp(joint - marginal_responsabilities[i]);
	}
}

double MixtureComponent::calculateTotalResponsability()
{
	double total_responsability = 0;
	int n = responsabilities_.size();
	for (int i = 0; i < n; i++) {
		total_responsability += responsabilities_[i];
	}
	//double log_total_responsability = log(total_responsability);
	return total_responsability;
}

void MixtureComponent::updateWeight(double total_resp, int n)
{
	weight_ = (double)total_resp/(double)n;
}

void MixtureComponent::updateCovariance(bool bound, double total_resp, int n, double lbr, double ubr, double lbm, double ubm)
{
	Eigen::Matrix2d new_cov;
    new_cov << 0, 0, 0, 0;
	for (int i = 0; i < n; i++) {
		Eigen::RowVector2d shift;
		shift << shift_.row(i);
		new_cov += responsabilities_[i] * shift.transpose() * shift;
	}
	new_cov = new_cov/total_resp;
	covariance_(0,0) = new_cov(0,0);
	covariance_(1,1) = new_cov(1,1);
	if (bound) {
		if (new_cov(0,0) > ubr) {
                covariance_(0,0) = ubr;
            } else if (new_cov(0,0) < lbr) {
        	covariance_(0,0) = lbr;
            }
            if (new_cov(1,1) > ubm) {
                covariance_(1,1) = ubm;
            } else if (new_cov(1,1) < lbm) {
                covariance_(1,1) = lbm;
            }
	}
}

void MixtureComponent::updateCovarianceAndWeight(bool bound, double lbr, double ubr, double lbm, double ubm)
{
	int n = responsabilities_.size();
	double total_resp = calculateTotalResponsability();
	updateWeight(total_resp, n);
	updateCovariance(bound, total_resp, n, lbr, ubr, lbm, ubm);
}
    

// Get methods
double MixtureComponent::getMixtureComponentWeight()
{
	return weight_;
}
Eigen::Matrix2d MixtureComponent::getMixtureComponentCov()
{
	return covariance_;
}
Eigen::Vector2d MixtureComponent::getMixtureComponentMu()
{
	return mu_;
}

double MixtureComponent::getMixtureComponentLikelihoodAt(int i)
{
	return likelihoods_[i];
}
    

// Set methods
void MixtureComponent::setMixtureComponentCovariance(const Eigen::Matrix2d& cov)
{
	covariance_ = cov;
}

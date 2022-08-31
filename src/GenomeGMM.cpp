//
// Created by Garance Jaques on 31.05.22.
//

#include <iostream>
#include <cmath>
#include <algorithm>
#include "GenomeGMM.h"
#include "myFunc.h"
#include "cumNormalDistr.h"

#include "ConfigFile.h"

GenomeGMM::GenomeGMM() {
    shifted_MAF_ = false;
    finalL_ = 0.0;
    number_of_iterations_ = 0;
}
GenomeGMM::~GenomeGMM() {}

float GenomeGMM::estimateSigmaAB(const std::vector<float>& ratio, const std::vector<float>& maf, int ploidy, float purity)
{
    std::cout << "..Estimate sigma of cluster AB" << std::endl;
    float estimated_sigma = -1;
    float ratio1 = (purity + 2*(1 - purity))/(purity*ploidy + 2*(1 - purity));
    float ratio2 = 2/(purity*ploidy + 2*(1 - purity));
    float ratio3 = (purity*3 + 2*(1 - purity))/(purity*ploidy + 2*(1 - purity));
    if (ratio.size() != maf.size()) {
        std::cerr << "Error: ratio and MAF vectors should be of the same size." << std::endl;
    }
    std::vector<float> maf_diploid;
    for (unsigned int i = 0; i < ratio.size(); i++) {
        if((ratio[i] >= ((ratio1 + ratio2)/2)) & (ratio[i] <= ((ratio2 + ratio3)/2))) {
	    float maf_val = maf[i];
            maf_diploid.push_back(maf_val);
        }
    }
    if (!maf_diploid.empty()) {
	    float max_maf = *std::max_element(maf_diploid.begin(), maf_diploid.end());
	    std::cout << "Max element of MAF diploid 2: " << max_maf << std::endl;
        int n = ceil((max_maf+0.001-0.500)/0.001);
        std::vector<float> bins(n);
        float x = 0.499;
        std::generate(bins.begin(), bins.end(), [&]{return x+=0.001;});
	    std::vector<float> counts = makeHistogramOfCounts(maf_diploid, bins);
	    std::map<int, float> peaks = findPeaks(counts, 2, 50, 20);
        if (peaks.size() != 0) {
		    // Find the peak that is closest to 0.5
            float maf_empirical_mean = max_maf + 0.001;
	        int peak_pos = 0;
	        float peak_freq = max_maf + 0.001;
	        float diff = abs(maf_empirical_mean - 0.500);
            for (std::map<int, float>::iterator it = peaks.begin(); it != peaks.end(); ++it) {
                peak_pos = it->first;
		        peak_freq = bins[peak_pos];
		        if (abs(peak_freq - 0.5) < diff) {
			        diff = abs(peak_freq - 0.500);
                    maf_empirical_mean = peak_freq;
                }
            }
	    std::cout << "MAF empirical mean 2: " << maf_empirical_mean << std::endl;
            if (maf_empirical_mean <= 0.6) {
		        estimated_sigma = (maf_empirical_mean - 0.5)*std::pow(2*PI_, 0.5)*0.5;
            }
        }
    }
    std::cout << "Estimated sigma: " << estimated_sigma << std::endl;
    return estimated_sigma;
}

float GenomeGMM::estimateSigmaAABB(const std::vector<float>& ratio, const std::vector<float>& maf, int ploidy, float purity)
{
    std::cout << "..Estimate sigma using cluster AABB" << std::endl;
    float estimated_sigma = -1;
    float ratio3 = (purity*3 + 2*(1 - purity))/(purity*ploidy + 2*(1 - purity));
    float ratio4 = (2 + 2*purity)/(purity*ploidy + 2*(1 - purity));
    float ratio5 = (purity*5 + 2*(1 - purity))/(purity*ploidy + 2*(1 - purity));
    if (ratio.size() != maf.size()) {
        std::cerr << "Error: ratio and MAF vectors should be of the same size." << std::endl;
    }
    std::vector<float> maf_diploid;
    for (unsigned int i = 0; i < ratio.size(); i++) {
        if((ratio[i] >= ((ratio3 + ratio4)/2)) & (ratio[i] <= ((ratio4 + ratio5)/2))) {
	        float maf_val = maf[i];
            maf_diploid.push_back(maf_val);
        }
    }
    if (!maf_diploid.empty()) {
        float max_maf = *std::max_element(maf_diploid.begin(), maf_diploid.end());
	std::cout << "Max element MAF diploid 4: " << max_maf << std::endl;
        int n = ceil((max_maf+0.001-0.500)/0.001);
        std::vector<float> bins(n);
        float x = 0.500;
        std::generate(bins.begin(), bins.end(), [&]{return x+=0.001;});
        std::vector<float> counts = makeHistogramOfCounts(maf_diploid, bins);
        std::map<int, float> peaks = findPeaks(counts, 2, 50, 20);
        if (peaks.size() != 0) {
            // Find the peak that is closest to 0.5
            float maf_empirical_mean = max_maf + 0.001;
            int peak_pos = 0;
            float peak_freq = maf_empirical_mean;
	        float diff = abs(maf_empirical_mean - 0.500);
            for (std::map<int, float>::iterator it = peaks.begin(); it != peaks.end(); ++it) {
                peak_pos = it->first;
                peak_freq = bins[peak_pos];
                if (abs(peak_freq - 0.5) < diff) {
                    diff = abs(peak_freq - 0.500);
                    maf_empirical_mean = peak_freq;
                }
            }
	    std::cout << "MAF empirical mean 4: " << maf_empirical_mean << std::endl;
            if (maf_empirical_mean <= 0.6) {
                estimated_sigma = (maf_empirical_mean - 0.5)*pow(2*PI_, 0.5)*0.5;
            }
        }
    }
    std::cout << "Estimated sigma: " << estimated_sigma << std::endl;
    return estimated_sigma;
}

float GenomeGMM::calculateSigma(int ploidy, float purity, int i, float sigma2, float sigma4)
{
    float threshold = 1.0/3.0;
	if (sigma2 != -1) {
        float sigma = sigma2*pow(1/expected_ratio_[i], 0.5) * pow(expected_maf_[i]*(1 - expected_maf_[i]), 0.5)/0.5;
        if (sigma < threshold) {
            return sigma;
        } else {
            std::cerr << "Warning: the estimated sigma is very high and is thus set to 1/3" << std::endl;
            sigma = threshold;
            return sigma;
        }
    } else if (sigma4 != -1) {
        float sigma = sigma4*pow(1/expected_ratio_[i], 0.5) * pow(expected_maf_[i]*(1-expected_maf_[i]), 0.5)/0.5;
        if (sigma < threshold) {
            return sigma;
        } else {
            std::cerr << "Warning: the estimated sigma is very high and is thus set to 1/3" << std::endl;
            sigma = threshold;
            return sigma;
        }
    } else {
        std::cerr << "Error: either sigma2 or sigma4 should be passed as non-default arguments!" << std::endl;
    }
}

void GenomeGMM::correctMAF(std::vector<float> sigma)
{
    std::cout << "..Correct MAF expected values" << std::endl;
    shifted_expected_maf_.resize(length_expectation_);
    float exp_BAF1, exp_BAF2 = 0;
    float pdf1, pdf2, cdf1, cdf2, term1, term2 = 0;
    for (int i = 0; i < length_expectation_; i++) {
        exp_BAF1 = expected_maf_[i];
        exp_BAF2 = 1 - exp_BAF1;
        pdf1 = NormalDistributionDensity(0.5, exp_BAF1, sigma[i]);
        pdf2 = NormalDistributionDensity(0.5, exp_BAF2, sigma[i]);
        cdf1 = calculateNormalCDF(0.5, exp_BAF1, sigma[i]);
        cdf2 = calculateNormalCDF(0.5, exp_BAF2, sigma[i]);
        term1 = exp_BAF1 * (1 - cdf1) + pow(sigma[i], 2) * pdf1;
        term2 = exp_BAF2 * (1 - cdf2) + pow(sigma[i], 2) * pdf2;
        shifted_expected_maf_[i] = term1 + term2;
	std::cout << "Shifted MAF: " << shifted_expected_maf_[i] << " and corresponding sigma: " << sigma[i] << std::endl;
    }
}

void GenomeGMM::shiftMAF(const std::vector<float>& ratio, const std::vector<float>& maf_smoothed, int ploidy, float purity)
{
    std::cout << "..Shift MAF centers" << std::endl;
    std::vector<float> estimated_sigma(length_expectation_);

    float sigma2 = estimateSigmaAB(ratio, maf_smoothed, ploidy, purity);
    if ((sigma2 != -1) & (ploidy < 4)) {
        for (int i = 0; i < length_expectation_; i++) {
            estimated_sigma[i] = calculateSigma(ploidy, purity, i, sigma2, -1);
        }
        correctMAF(estimated_sigma);
        shifted_MAF_ = true;
    } else {
        float sigma4 = estimateSigmaAABB(ratio, maf_smoothed, ploidy, purity);
        if (sigma4 != -1) {
            for (int i = 0; i < length_expectation_; i++) {
                estimated_sigma[i] = calculateSigma(ploidy, purity, i, -1, sigma4);
            }
            correctMAF(estimated_sigma);
            shifted_MAF_ = true;
        } else {
            std::cerr << "Warning! Cluster centers were not corrected for shifting" << std::endl;
	    shifted_expected_maf_ = expected_maf_;
        }
    }
}

void GenomeGMM::calculateExpectedRatioAndMAF(int ploidy, float purity)
{
    std::cout << "..Calculate expected ratio and expected MAF" << std::endl;
    int max_copy_number = 0;
    switch(ploidy) {
        case 1:
            max_copy_number = 5;
            break;
        case 2:
            max_copy_number = 7;
            break;
        case 3:
            max_copy_number = 10;
            break;
        case 4:
            max_copy_number = 12;
            break;
        default:
            max_copy_number = 15;
    }
    std::vector<int> Q;
    std::vector<int> Qh;
    for(int i = 0; i <= max_copy_number; i++) {
        for(int j = ceil((float)i/(float)2); j <= i; j++) {
            Q.push_back(i);
            Qh.push_back(j);
        }
    }
    int n = Q.size();
    length_expectation_ = n;
    expected_ratio_.reserve(n);
    expected_maf_.reserve(n);
    float normalFraction = 2*(1 - purity);
    for(int i = 0; i < n; i++) {
        expected_ratio_.push_back((Q[i]*purity + normalFraction)/(purity*ploidy + normalFraction + 1e-8));
        expected_maf_.push_back((purity*Qh[i] + 1 - purity)/(purity*Q[i] + normalFraction + 1e-8));
    }
    K_ = expected_ratio_.size();
}


Eigen::VectorXd GenomeGMM::transformVectorToMatrix(const std::vector<double>& v)
{
	const double* cpv = v.data();
	Eigen::VectorXd vxd = Eigen::VectorXd::Map(cpv, v.size());
	return vxd;
}

void GenomeGMM::initialiseGMMParameters(int n_init, int max_iter, double tol, double reg_cov, double noise, bool b, double lbr, double ubr, double lbm, double ubm)
{
	setNinit(n_init);
	setMaxIter(max_iter);
	setTolerance(tol);
	setRegCovar(reg_cov);
	setNoise(noise);
	setBound(b);
	setBoundsRatio(lbr, ubr);
	setBoundsMAF(lbm, ubm);
}

Eigen::Matrix2d GenomeGMM::initialiseCovariance() {
	Eigen::Matrix2d cov;
    cov << 0.005, 0.0, 0.0, 0.01;
	return cov;
}

double GenomeGMM::initialiseWeight() {
	return (double)1/(double)K_;
}

std::vector<double> GenomeGMM::calculateMarginalResponsabilities(int n)
{
	std::vector<double> marginal_responsabilities;
    marginal_responsabilities.resize(n);
	for (int i = 0; i < n; i++) {
		double marginal_responsability = 0;
		for (int k = 0; k < K_; k++) {
			marginal_responsability += exp(log(mixture_[k].getMixtureComponentWeight()) + mixture_[k].getMixtureComponentLikelihoodAt(i));
		}
                marginal_responsabilities[i] = log(marginal_responsability);
	}
	return marginal_responsabilities;
}

double GenomeGMM::computeLoss(int n)
{
    double likelihood = 0.0;
    for (int i = 0; i < n; i++) {
        double marg_likelihood = 0.0;
        for (int k = 0; k < K_; k++) {
	    double component_weight = mixture_[k].getMixtureComponentWeight();
	    double l = exp(log(component_weight) + mixture_[k].getMixtureComponentLikelihoodAt(i));
            marg_likelihood += l;
        }
	double log_marg = log(marg_likelihood);
	likelihood += log_marg;
    }
    return likelihood;
}

double GenomeGMM::fitPredict(int tau, float alpha, const std::vector<float>& ratio, const std::vector<float>& maf)
{
    // TMP
    #ifdef CLUSTER
        ConfigFile cf("/cluster/home/kzaitse/custom_configs/freec_configs/meta-config.txt");
    #else
        ConfigFile cf("/local/home/kzaitse/lab/custom_configs/freec_configs/meta-config.txt");
    #endif // !CLUSTER
	int maxIter = (int)cf.Value("gmm", "maxIter", 1);
    initialiseGMMParameters(1, maxIter, 2000, 1.0e-7, 1, false, 0.0, 0.01, 0.0, 0.015);
	
    if (noise_ == 1) {
		expected_ratio_.push_back(1.0);
		expected_maf_.push_back(0.95);
		shifted_expected_maf_.push_back(0.95);
		K_ = expected_ratio_.size();
	}
	mixture_.resize(K_);
	int N = ratio.size();
	
	// Transform std::vectors in Eigen::vectors for matrix operations
    std::vector<double> d_ratio(ratio.begin(), ratio.end());
    std::vector<double> d_maf(maf.begin(), maf.end());
	Eigen::VectorXd e_ratio = transformVectorToMatrix(d_ratio);
    Eigen::VectorXd e_maf = transformVectorToMatrix(d_maf);
	
	// Initialise all mixture components
	Eigen::Matrix2d cov = initialiseCovariance();
	double w = initialiseWeight();
	for (int k = 0; k < K_; k++) {
		std::vector<double> mu = {expected_ratio_[k], shifted_expected_maf_[k]};
		Eigen::VectorXd ex_mu;
		ex_mu = transformVectorToMatrix(mu);
        	Eigen::Vector2d e_mu;
        	if (ex_mu.size() == 2) {
			e_mu << ex_mu;
        	} else {
            		std::cerr << "Error: mu needs to be of size 2!" << std::endl;
        	}
		MixtureComponent m(e_ratio, e_maf, e_mu, w, cov, N);
        	mixture_[k] = m;
	}
	
	// Initialise variables used in the EM algorithm
	std::vector<double> marginal_responsabilities(N);
    double currentL = 0.0;
    double previousL = 0.0;
	
	// EM algorithm
	int iter = 0;
	double diff = std::numeric_limits<float>::infinity();
	std::cout << "Begin of EM algorithm" << std::endl;
	while((iter < max_iter_) & (abs(diff) > tolerance_)) {
		std::cout << tau << "/" << alpha << ":\tIteration " << iter << std::endl;
		// E-step
		std::cout << "..Calculate sample likelihood" << std::endl;
		for (int k = 0; k < K_; ++k) {
			mixture_[k].calculateSampleLikelihoods(N);
		}
		marginal_responsabilities = calculateMarginalResponsabilities(N);
		for (int k = 0; k < K_; k++) {
			mixture_[k].calculateComponentResponsabilities(marginal_responsabilities);
			mixture_[k].updateCovarianceAndWeight(bound_, lower_bound_ratio_, upper_bound_ratio_, lower_bound_maf_, upper_bound_maf_);
			mixture_[k].calculateSampleLikelihoods(N);
		}
        previousL = currentL;
        currentL = computeLoss(N);
		std::cout << tau << "/" << alpha << ":\tCurrent likelihood: " << currentL << std::endl;
        diff = (currentL - previousL);

        if (diff > 1e8){
            std::cout << "Clipping difference between current and previous likelihood to 1e8" << diff << std::endl;
            diff = 1e8;
        } else if (diff < -1e8){
            std::cout << "Clipping difference between current and previous likelihood to -1e8" << diff << std::endl;
            diff = -1e8;
        } else{
		    std::cout << tau << "/" << alpha << ":\tDifference between current and previous likelihood: " << diff << std::endl;
        }

		iter += 1;
		if (iter >= max_iter_) {
			std::cout << "The algorithm stopped because it reached the maximum number of iterations " << max_iter_ << std::endl;
		} else if (abs(diff) <= tolerance_) {
			std::cout << "The algorithm stopped because the difference between the current and the previous likelihood is below the given tolerance " << tolerance_ << std::endl;
		}
	}
	finalL_ = currentL;
	number_of_iterations_ = iter;
    return finalL_;
}


void GenomeGMM::printExpectedRatioAndMAF(std::string const& outFile, int ploidy, float purity) {
    std::string myName = outFile + "_cluster-centers_ploidy" + std::to_string(ploidy) + "_purity" + std::to_string(purity) + ".txt";
    std::ofstream file (myName.c_str());

    file << "Cluster centers coordinates for ploidy " << ploidy << " and purity " << purity;
    file << "\n" << "\n" << "\n";
    file << "ExpectedRatio\tExpectedMAF";
    if (shifted_MAF_) {
        file << "\tShiftedExpectedMAF";
    }
    file << "\n";

    for (int i = 0; i < expected_ratio_.size(); i++) {
        file << expected_ratio_[i] << "\t" << expected_maf_[i];
        if (shifted_MAF_) {
            file << "\t" << shifted_expected_maf_[i];
        }
        file << "\n";
    }
    file.close();
    std::cout << "Wrote CLUSTER CENTERS to file:\t" << myName << std::endl;
}

double GenomeGMM::getFinalL()
{
	return finalL_;
}

int GenomeGMM::getNumberOfIterations()
{
	return number_of_iterations_;
}

void GenomeGMM::setNinit(int n) 
{
	n_init_ = n;
}

void GenomeGMM::setMaxIter(int m)
{
	max_iter_ = m;
}

void GenomeGMM::setTolerance(double t)
{
	tolerance_ = t;
}
void GenomeGMM::setRegCovar(double c)
{
	reg_covar_ = c;
}

void GenomeGMM::setNoise(double n)
{
	noise_ = n;
}

void GenomeGMM::setBound(bool b)
{
	bound_ = b;
}

void GenomeGMM::setBoundsRatio(double l, double u)
{
	lower_bound_ratio_ = l;
	upper_bound_ratio_ = u;
}

void GenomeGMM::setBoundsMAF(double l, double u)
{
	lower_bound_maf_ = l;
	upper_bound_maf_ = u;
}

std::vector<MixtureComponent> GenomeGMM::getMixtures() {
    return mixture_;
}

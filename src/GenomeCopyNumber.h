#ifndef HEADER_4D7F42FDD0838B93
#define HEADER_4D7F42FDD0838B93

#pragma once
#ifndef _GENOME_CPN_H
#define _GENOME_CPN_H

#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <map>
#include <fstream>
#include <bitset>
#include "SVfinder.h"
#include <utility>

#include "myFunc.h"
#include "ChrCopyNumber.h"
#include "SVfinder.h"
#include "chisquaredistr.h" //to calculate chisquare distribution
#include "EntryCNV.h"
#include "SNPinGenome.h"
#include "GenomeGMM.h"
#include "bayes_opt.h"
#include <Eigen/Dense>


class GenomeCopyNumber
{
public:
	GenomeCopyNumber(void);
	~GenomeCopyNumber(void);

	void readCopyNumber(std::string const& mateFileName ,std::string const& inputFormat, std::string const& matesOrientation, std::string const& chrLenFileName, int windowSize, int step, std::string targetBed = "" );
	void readCopyNumber(std::string const& mateFileName ,std::string const& inputFormat, std::string const& matesOrientation, std::string const& chrLenFileName, float coefficientOfVariation );
	void readCopyNumber(std::string const& inFile);
	int readCGprofile(std::string const& inFile);
	void readGemMappabilityFile(std::string const& inFile);
    int processRead(InputFormat inputFormat, MateOrientation matesOrientation, const char* line_buffer,  int& prevInd, std::string targetBed = "", std::string mateFileName = "");
	int processReadWithBowtie(std::string const& inputFormat, std::string const& matesOrientation,std::string const line,std::string const line2);
    int focusOnCapture (std::string const& captureFile);
    float removeLargeExons(float iqrToKeep);
	void initCopyNumber(std::string const& chrLenFileName, int windowSize , int step, std::string targetBed);
	void finishCopyNumber(long normalCount);
    void addBAFinfo(SNPinGenome & snpingenome);
    void removeLowReadCountWindows(GenomeCopyNumber & controlCopyNumber, int RCThresh);
    void removeLowReadCountWindowsFromControl (int RCThresh);

    int fillInRatio();
	int calculateRatio( GenomeCopyNumber & controlCopyNumber, int degree, bool intercept) ;
    void calculateRatioUsingCG( GenomeCopyNumber & controlCopyNumber) ;
    void calculateRatioUsingCG_Regression( GenomeCopyNumber & controlCopyNumber) ;
	float calculateNormalizationConstant(GenomeCopyNumber & controlCopyNumber);
	void calculateBreakpoints(double breakPointThreshold, int breakPointType);
	void calculateBAFBreakpoints(double breakPointThreshold, int breakPointType);
	void calculateSDAndMed(int ploidy,std::map<float,float> &sds,std::map<float,float> &meds);
	void calculateSDs(int ploidy, std::map <float,float> &sds) ; //only SDs
	float calculateVarianceForNormalCopy(int ploidy); //OLD
	void calculatePloidy(int minCNAlength);
	double calculateXiSum(int ploidy); //calculates sum_{i}{(med_i-supposedValue_i) / SQRT(var_i)}   var_i = pi/2/n*sd_i*2 which should be distributed as Xi_sqaure if the null hypo is correct
	void calculateCopyNumberProbs_and_genomeLength(int breakPointType) ;
	void calculateCopyNumberProbs_and_exomeLength(int breakPointType) ;
	void deleteFlanks(int telo_centromeric_flanks);
	void recalcFlanks(int telo_centromeric_flanks, int minNumberOfWindows);
	int calculateRatioUsingCG (bool intercept, float minExpectedGC, float maxExpectedGC) ; //will try different degrees; returns 1 if #interation < max
	int calculateRatioUsingCG (int degree, bool intercept, float minExpectedGC, float maxExpectedGC) ;

    void recalculateRatioUsingCG (int degree, bool intercept, float minExpectedGC, float maxExpectedGC) ;
	void recalculateRatio (float contamination);
	void calculateCopyNumberMedians (int minCNAlength, bool noisyData, bool CompleteGenomicsData);
	double calculateMedianRatioAround (float interval, float around);
    double calculateMedianAround (float interval, float around);
    void calculateSomaticCNVs (std::vector <EntryCNV> controlCNVs, int controlPloidy);
    int calculateMedianReadCountPerWindow();
    int calculateSDReadCountPerWindow(int mean);
    void calculateReadCountQuartiles(int& lower, int& upper);

	double calculateMedianAround (GenomeCopyNumber & controlCopyNumber, float interval, float around);
	float evaluateContamination();
    float evaluateContaminationwithLR();

	void printRatio(std::string const& outFile, bool ifBedGraphOutPut, bool printNA, bool has_smoothed_profile = false, bool has_MAF = false);
	void printRatio(std::string const& chr, std::string const& outFile, bool printNA, bool has_smoothed_profile, bool has_MAF);
	void printRatio(std::string const& chr, std::ofstream & file, bool printNA, bool has_smoothed_profile, bool has_MAF);
    void printRatioBedGraph(std::string const& chr, std::ofstream & file, std::string const& typeCNA);
	void printPloidy(std::string const& outFile) ;
	void printPloidy(std::ofstream & file);
	void printCopyNumber(std::string const& outFile);
	void printCopyNumber(std::string const& chr, std::ofstream & file) ;
	void printCopyNumber(std::string const& chr, std::string const& outFile);
	void printCGprofile(std::string const& outFile);
	void printCGprofile(std::string const& chr, std::ofstream & file);
	void printCNVs (std::string const& outFile);
	void printBAF(std::string const& outFile, SNPinGenome& snpingenome, std::string const& matefile = "");
    void printBAF(std::string const& chr, std::ofstream & file, SNPatChr& snpAtChrom, std::string myName, std::string const& matefile);
    void printInfo(std::ofstream & file) ;

    void shiftNeutalRatioTo1();


	int findIndex (std::string const& chr);
	void fillCGprofile(std::string const& chrFolder);

	ChrCopyNumber getChrCopyNumber(std::string const& chr);
	float getMedianCopyNumber();
	long getTotalNumberOfPairs();
	float getMedianRatio();
	int getWindowSize(void);
	long getNormalNumberOfPairs();
	std::vector <EntryCNV> getCNVs ();
    int getPloidy();
	ChrCopyNumber& getChrCopyNumberAt(int index) {return chrCopyNumber_[index];}
	int getStep() const {return step_;}
    bool getWESanalysis();
    double getGenomeRefSize();
    int getNumberOfChromosomes();
    bool ifHasRatio();


	void setPloidy(int ploidy);
	void setStep(int step);
    void setBAFtrue();
    void setNormalContamination(float normalContamination) ;
    void setAllNormal ();
    void setSamtools(std::string const& pathToSamtools, std::string const& SambambaThreads_);
    void setSambamba(std::string const& pathToSambamba, std::string const& SambambaThreads_);
    void setReference(std::string const& pathToReference);
    bool ifHasBAF();
    void setSex(std::string sex);
    void setSeekSubclones(bool seekSubclones);
    void setDataSubsamplingRateInPuritySearch(int dataSubsamplingRateInPuritySearch);
    void setDataSubsamplingRateInPloidyEvaluation(int dataSubsamplingRateInPloidyEvaluation);
    void setDoRescaleRatio(bool doRescaleRatio);

    int findWinNumber(int position, std::string myName, std::string const& matefile);
    void setWESanalysis(bool WESgiven);
    void setmakingPileup(bool makingPileup_given);
    void setIfLogged(bool);

    double Percentage_GenomeExplained(int &);
    long double calculateRSS(int ploidy);
    bool isMappUsed();
    int getDataSubsamplingRateInPuritySearch();
    int getDataSubsamplingRateInPloidyEvaluation();

    // Added by Gara
    float findMaximumPeak(std::vector<float> total_ratio);
    void rescaleRatio();
    void preprocessing(SNPinGenome& snpingenome);
    // Print preprocessing outputs
    void printRatioAfterPreprocessing(std::string const& outFile);
    // Get ratio and MAF vectors that span on the whole genome, instead of a single chromosome
    float getFilteredRatio(int index, int i);
    float getFilteredMAF(int index, int i);
    // For a given pair of ploidy tau and purity alpha, the means of ratio and MAF are fixed
    void fixExpectedRatioAndMAF(GenomeGMM& gmm, const std::vector<float>& ratio, const std::vector<float>& sm_maf, int ploidy, float purity);
    // Fit the GMM
    double processGMM(std::string const& outFile, int tau, float alpha, bool take_subsample_flag, int take_each_nth_element_for_gmm_fit, bool usePenalty = true);

    void printGMMFinalLikelihoods(std::string const& outFile, const std::vector<std::pair<int, double>>& ploidy_purity_pairs, const std::vector<double> likelihoods, const std::vector<int> number_of_iterations);

private:

    std::vector<ChrCopyNumber> chrCopyNumber_; //should stay private !!!
	std::map<std::string, int> chromosomesInd_; //should stay private!!!
    bool WESanalysis;
    bool makingPileup;
    bool SeekingSubc_;
    bool isMappUsed_;
    bool isRatioLogged_;

	void fillMyHash(std::string const& mateFileName , std::string const& inputFormat, std::string const& matesOrientation, int windowSize, int step, std::string targetBed = "");
    int  dataSubsamplingRateInPuritySearch;
    int  dataSubsamplingRateInPloidyEvaluation;
    bool doRescaleRatio;
	int windowSize_;
	int step_;
	long totalNumberOfPairs_;
	long normalNumberOfPairs_;
	double refGenomeSize_;
	int ploidy_;
	double ploidy_pvalue_;
	double estimationOfGenomeSize_;
	double estimationOfExomeSize_;
	std::map<int, double> copyNumberProbs_;
	int telo_centromeric_flanks_;
	std::vector <EntryCNV> CNVs_;
	bool hasBAF_;
	bool ifUsedControl_;
	float normalContamination_;
	std::string sex_;
	std::string pathToSamtools_;
	std::string pathToSambamba_;
	std::string pathToReference_;
	std::string SambambaThreads_;

    float eval_objective(GenomeGMM sampleGMM, int tau, float alpha, std::vector<float>& genome_ratio_smoothed,
                         std::vector<float>& genome_maf_smoothed);

    float
    compute_penalty(const std::vector<double> &weights, const std::vector<Eigen::Matrix2d> &covs, Eigen::MatrixXd mu,
                    int N,
                    float threshold_empty = 5e-3,
                    float threshold_noise = 0.1,
                    int C = 7500);
};
#endif

//
// Multi Thread Support
//

struct GenomeCopyNumberCalculateBreakpointArgWrapper : public ThreadArg {
  GenomeCopyNumber& genomeCopyNumber;
  double breakPointThreshold;
  int breakPointType;

  GenomeCopyNumberCalculateBreakpointArgWrapper(GenomeCopyNumber& genomeCopyNumber, double breakPointThreshold, int breakPointType) : genomeCopyNumber(genomeCopyNumber), breakPointThreshold(breakPointThreshold), breakPointType(breakPointType) { }
};

extern void* GenomeCopyNumber_calculateBreakpoint_wrapper(void *arg);

struct GenomeCopyNumberReadMateFileArgWrapper : public ThreadArg {
  SNPinGenome& snpInGenome;
  std::string mateFile;
  std::string inputFormat;
  int minimalTotalLetterCountPerPosition;
  int minimalQualityPerPosition;
  GenomeCopyNumber* p_genomeCopyNumber;
  std::string chrLenFileName;
  int windowSize;
  int step;
  std::string targetBed;

  GenomeCopyNumberReadMateFileArgWrapper(SNPinGenome& snpInGenome, std::string const& mateFile, const std::string& inputFormat, int minimalTotalLetterCountPerPosition, int minimalQualityPerPosition) : snpInGenome(snpInGenome), mateFile(mateFile), inputFormat(inputFormat), minimalTotalLetterCountPerPosition(minimalTotalLetterCountPerPosition), minimalQualityPerPosition(minimalQualityPerPosition), p_genomeCopyNumber(NULL), windowSize(0), step(0) { }

  GenomeCopyNumberReadMateFileArgWrapper(SNPinGenome& snpInGenome, std::string const& mateFile, std::string const& inputFormat, int minimalTotalLetterCountPerPosition, int minimalQualityPerPosition, GenomeCopyNumber& genomeCopyNumber, std::string const& chrLenFileName, int windowSize, int step, std::string targetBed) :  snpInGenome(snpInGenome), mateFile(mateFile), inputFormat(inputFormat), minimalTotalLetterCountPerPosition(minimalTotalLetterCountPerPosition), minimalQualityPerPosition(minimalQualityPerPosition), p_genomeCopyNumber(&genomeCopyNumber), chrLenFileName(chrLenFileName), windowSize(windowSize), step(step), targetBed(targetBed) { }
};

extern void* GenomeCopyNumber_readMateFile_wrapper(void *arg);

#endif // header guard



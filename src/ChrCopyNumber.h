#ifndef HEADER_6AE83A503DDCA686
#define HEADER_6AE83A503DDCA686

#pragma once
#ifndef _CHR_CPN_H
#define _CHR_CPN_H

#include <vector>
#include <map>
#include <string>
#include <stdlib.h>
#include <iostream>
#include "myFunc.h"
#include "SVfinder.h"
#include "SNPatChr.h"
#include "ThreadPool.h"

class ChrCopyNumber
{
public:
	ChrCopyNumber(void);
	ChrCopyNumber(std::string const& chrName);
	ChrCopyNumber(int windowSize, int chrLength, std::string const& chrName);
	ChrCopyNumber(int windowSize, int chrLength, std::string const& chrName, int step, std::string targetBed = "");
	~ChrCopyNumber(void);
	void mappedPlusOneAtI(int i, int step, int l = -1);

	void addBAFinfo(SNPinGenome & snpingenome,int indexSNP);

	void fillInRatio(bool islog);
	void calculateRatio(ChrCopyNumber control, float normalizationConst) ;
	void recalculateRatio (float constant);
	void recalculateLogRatio (float constant) ;
	void recalculateRatioWithContam(float contamination, float normGenytype, bool isLogged);
    void recalculateRatio(ChrCopyNumber control);
	void calculateRatio(ChrCopyNumber control, double a0, double a1);
	void calculateRatio(ChrCopyNumber control, const double * a, const int degree);
	void calculateRatioLog(ChrCopyNumber control, const double * a, const int degree);
	int calculateBreakpoints(double breakPointThreshold, int firstChrLength, int breakPointType);
    int calculateBAFBreakpoints(double breakPointThreshold, int firstChrLength, int breakPointType);
	double calculateXiSum(int ploidy, std::map <float,float> &sds, std::map <float,float> &meds);
	double calculateXiSum(int ploidy, std::map <float,float> &sds);
	void calculateCopyNumberMedian(); //create median profiles using 'bpfinal_' and store them in medianProfile_, info about medians themselves is stored in medianValues_ and about SD in sd_, lengths of fragments in bpLengths_
	void calculateCopyNumberMedian(int ploidy, int minCNAlength, bool noisyData, bool CompleteGenomicsData, bool isLogged); //create median profiles as calculateCopyNumberMedian(), but merges close regions (roundByPloidy(median))
    void recalcFlanksForIndeces (int i, int j);
    void recalcFlanks(int telo_centromeric_flanks, int minNumberOfWindows); //merge short notNA-segments around NA-segments

	//void calculateRatio(double a,double b, double c, double d); //depreciated
	void calculateRatio(double *a, int degree);
    void recalculateRatio(double *a, int degree);

    void removeLowReadCountWindows(ChrCopyNumber control, const int RCThresh);
    void removeLowReadCountWindows(const int RCThresh) ;

	void deleteFlanks(int telo_centromeric_flanks);
	void deleteFragment(int i) ;
    int removeLargeExons(float threshold);

    std::string getGeneNameAtBin(int i);
	float		getValueAt(int i);
	int			getCoordinateAtBin(int i);
	int			getEndAtBin(int i);
	//void printLog2Ratio(std::ofstream const& file) ;
	float		getRatioAtBin(int i);
	std::vector <float>	getRatio();
	int			getLength();
    int         getChrLength();
	std::string getChromosome();
	int getCoveredPart(int breakPointStart, int breakPointEnd);//for exome-seq: get length of the genome covered by the targeted region (from breakPointStart to breakPointEnd)
	std::vector <float> getValues(); //get readCount_
	float		getMedianProfileAtI (int i) ;
	std::vector <float> getMedianValues () ;
	float getMedianValuesAt (int i) ;
	std::vector <float> getSDs ();//unUsed
	std::vector <int> getFragmentLengths_notNA ();
	std::vector <int> getFragmentLengths ();
	int getFragmentLengthsAt(int i);
    int getFragmentLengths_notNA_At (int i);
	int			getNumberOfFragments();
	int			getNumberOfGoodFragments();
	double		getXiSum(int ploidy, float minSD); //OLD unUsed
	std::vector <int> getBreakPoints();
	float getCGprofileAt(int i);
	float getMappabilityProfileAt(int i);
	float getNotNprofileAt(int i);
	float getSmoothedProfileAtI(int i);
	float getBAFat (int i);
	float getBAFProfileAt (int i);
	std::string getBAFsymbolAt (int i);
	float getEstimatedBAFuncertaintyAtI(int i) ;
	std::string getBAFsymbPerFrg (int i) ;
    float getEstimatedBAFuncertaintyAtBin(int i);
    float getFittedBAFProfileAt (int i);


	int getEndsSize();
	int	getMappabilityLength();

	bool isMedianCalculated();
	bool isSmoothed();
    bool ifHasRatio();


	void fillCGprofile(std::string const& chrFolder);

	float nextNoNAMedian(int i1, int ploidy);
	int nextNoNALength(int i1, int ploidy);
	int nextNoNAIndex(int i1, int ploidy, int min_fragment);

	void addToReadCount(float);
	void addToCoordinates(int);
	void addToEnds(int i);
	void addToCGcontent (float valueToAdd);
	void addToNonNpercent (float valueToAdd);
    void addToMappabilityProfile(float valueToAdd);
    void addToGenes_name(std::string i);

    void clearCGcontent () ;
    void clearNonNpercent () ;
    void clearMappabilityProfile ();

	void setWindowSize(int);
	void setVectorLength(int);
	void setChrLength(int);
	void setIsSmoothed(bool value);
	void setNotNprofileAt(int i, float value);
	void setMappabilityProfileAt(int i, float value);
	void setBAFat(int i, float value);
	void setStep(int step);
	void setPloidy(int ploidy);
	void setNormalContamination(float normalContamination);
    void setRCountToZeroForNNNN();
    void setAllNormal();
    void setValueAt(int i, float val) ;


    void createMappabilityProfile();
    void createBAF(float);
    void createBAFvalues();
    void checkOrCreateNotNprofileWithZeros();

	void pushSmoothedProfile(float value);

	float getLevelAt(int unsigned i, int ploidy);

	int getExons_Countchr();
	void setCN_subc(int i, int CN_subc);
	int getCN_subc(int i);
	void setPopulation_subc(int i, float pop_subc);
	float getPopulation_subc(int i);
	void setCN_subcLength(int len);
	void setpop_subcLength(int len);
    void setLookingForSubclones(bool);
    float getSmoothedForInterval(int start , int end);

private:
   // std::vector <std::string> coordinatesTmp_;
//	std::vector <std::string> endsTmp_;
//	std::vector <std::string> chr_namestmp;
 //   std::vector <std::string> chr_names;
	std::vector <std::string> genes_names;
  //  int exons_Count;
	int exons_Countchr_;
	std::vector <int> copy_number_subc_;
	std::vector <float> population_subc_;

	bool isLookingForSubclones_;


    int ploidy_;
    float normalContamination_;
	int chrLength_;
	int length_;
	int windowSize_;
	int step_;
	bool isMedianCalculated_;
	bool isSmoothed_; //medianProfileHasBeenSmoothed
	std::string chromosome_;
	std::vector <int> coordinates_;
	std::vector <int> ends_;
	std::vector <float> readCount_;
	std::vector <float> ratio_;
	std::vector <int> bpfinal_;
	std::vector <int> fragmentNotNA_lengths_;
	std::vector <int> fragment_lengths_;
	std::vector <float> medianValues_; //medianValues for each segment
	std::vector <float> medianProfile_; //medianValues for each position
	std::vector <float> sd_;
	std::vector <float> GCprofile_;//CG-content in a window
	std::vector <float> notNprofile_;//percentage of not 'N' in a window
	std::vector <float> mappabilityProfile_;//percentage of mappable positions in a window
	std::vector <float> smoothedProfile_;
    std::vector <float> BAF_; //median value of abs(BAF-0.5) in a window (only for heterozygous SNPs)
    std::vector <float> medianBAFProfile_; //median values of BAF_ for each segment)
    std::vector <std::string> BAFvalues_; // values of abs(BAF-0.5) in a window
    std::vector <float> estimatedBAFProfile_; //estimation of BAF for each segment (value per window)
    std::vector <float> fittedBAFProfile_;
    std::vector <std::string> medianBAFSymbol_; //estimation of BAF for each segment: AA,AAB;AB,AABB etc
    std::vector <float> estimatedBAFuncertainty_;//uncertainty of estimation of BAF for each segment (1./(LL_best-LL_secondBest)))
    std::vector <std::string> BAFsymbPerFrag_; //estimation of BAF for each segment: AA,AAB;AB,AABB etc - one value per fragment
    std::vector <float> estBAFuncertaintyPerFrag_;//uncertainty of estimation of BAF for each segment (1./(LL_best-LL_secondBest))) - one value per fragment


};
#endif

//
// Multi Thread Support
//

struct ChrCopyNumberCalculateBreakpointArgWrapper : public ThreadArg {
  ChrCopyNumber& chrCopyNumber;
  double breakPointThreshold;
  int breakPointType;

  ChrCopyNumberCalculateBreakpointArgWrapper(ChrCopyNumber& chrCopyNumber, double breakPointThreshold, int breakPointType) : chrCopyNumber(chrCopyNumber), breakPointThreshold(breakPointThreshold), breakPointType(breakPointType) { }
};

extern void* ChrCopyNumber_calculateBreakpoint_wrapper(void *arg);
extern void* ChrCopyNumber_calculateBAFBreakpoint_wrapper(void *arg);

#endif // header guard

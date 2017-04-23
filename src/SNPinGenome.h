#ifndef SNPINGENOME_H
#define SNPINGENOME_H

#include <string>
#include <vector>
#include <iostream>

#include "SNPatChr.h"
#include "binomialdistr.h"
#include "ThreadPool.h"

#define ERROR_PER_POS 0.01

class GenomeCopyNumber;

class SNPinGenome
{
    public:
	    SNPinGenome();
  	    void setSNPChr(std::vector<SNPatChr>* SNP_atChr_);
		void readSNPs(std::string const& inFile);
		void perform(std::string const& mateFile, const std::string& inputFormat, int minimalTotalLetterCountPerPosition, int minimalQualityPerPosition, bool noisyData, bool CompleteGenomicsData, GenomeCopyNumber& genomeCopyNumber, double breakPointThreshold, int breakPointType, int minCNAlength, const char* what);
        virtual ~SNPinGenome();
        void assignValues(std::string const& inFile, std::string inputFormat,int minimalTotalLetterCountPerPosition, int minimalQualityPerPosition, GenomeCopyNumber* p_genomeCopyNumber = NULL);

        void setBinAt(int indexSNP,int SNPcount,int valueToSet);
        void setWESanalysis(bool WESanalysis);
        void setCopyNumberFromPileup(bool CopyNumberFromPileup);

        int findIndex (std::string chromosome) const;
        const SNPatChr& SNP_atChr(int) const;
        SNPatChr& SNP_atChr(int);
		void readMateFile(std::string const& mateFile, const std::string& inputFormat, int minimalTotalLetterCountPerPosition, int minimalQualityPerPosition);
		void readMateFile(std::string const& mateFile, std::string const& inputFormat, int minimalTotalLetterCountPerPosition, int minimalQualityPerPosition, GenomeCopyNumber& genomeCopyNumber, std::string const& chrLenFileName, int windowSize, int step,  std::string targetBed = "");
        std::vector <SNPatChr>* getSNPChr() {return SNP_atChr_;}

    protected:
    private:
        std::vector <SNPatChr>* SNP_atChr_;
        float addInfoFromAPileUp (int totalLetterCount, int minimalTotalLetterCountPerPosition,char whatToLook,
                                  int index,int &positionCount, int &sNPpositionToProceed,const char * pileup,int minimalQualityPerPosition,const char * quality); //returns BAF in case on heterozygous SNP (NA otherwise)
		void readPileUP(FILE* stream, int minimalTotalLetterCountPerPosition, int minimalQualityPerPosition, GenomeCopyNumber* p_genomeCopyNumber);
        long processPileUPLine(int & positionCount, char* line, std::string & oldChr, int & sNPpositionToProceed,int minimalTotalLetterCountPerPosition, int & index, int minimalQualityPerPosition, GenomeCopyNumber* p_genomeCopyNumber);
		int processSNPLine(bool isVCF, char * line, std::string & myChr, int & index,int &previousPos) ;
		bool pileup_read;
        bool WESanalysis_;
        bool CopyNumberFromPileup_;
};


//
// Multi Thread Support
//

struct SNPinGenomePerformArgWrapper : public ThreadArg {
  SNPinGenome& snpingenome;
  std::string mateFile;
  std::string inputFormat;
  int minimalTotalLetterCountPerPosition;
  int minimalQualityPerPosition;
  bool noisyData;
  bool CompleteGenomicsData;
  GenomeCopyNumber& genomeCopyNumber;
  double breakPointThreshold;
  int breakPointType;
  int minCNAlength;
  const char* what;

  SNPinGenomePerformArgWrapper(SNPinGenome& snpingenome, std::string const& mateFile, const std::string& inputFormat, int minimalTotalLetterCountPerPosition, int minimalQualityPerPosition, bool noisyData, bool CompleteGenomicsData, GenomeCopyNumber& genomeCopyNumber, double breakPointThreshold, int breakPointType, int minCNAlength, const char* what) : snpingenome(snpingenome), mateFile(mateFile), inputFormat(inputFormat), minimalTotalLetterCountPerPosition(minimalTotalLetterCountPerPosition), minimalQualityPerPosition(minimalQualityPerPosition), noisyData(noisyData),CompleteGenomicsData(CompleteGenomicsData),genomeCopyNumber(genomeCopyNumber), breakPointThreshold(breakPointThreshold), breakPointType(breakPointType), minCNAlength(minCNAlength), what(what) { }
};

extern void* SNPinGenome_perform_wrapper(void *arg);

#endif // SNPINGENOME_H

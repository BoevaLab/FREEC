/*************************************************************************
Copyright (c) 2010-2011, Valentina BOEVA.

>>> SOURCE LICENSE >>>
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation (www.fsf.org); either version 2 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License is available at
http://www.fsf.org/licensing/licenses

>>> END OF LICENSE >>>
*************************************************************************/

#include "SNPinGenome.h"
#include "GenomeCopyNumber.h"

using namespace std;

SNPinGenome::SNPinGenome() : SNP_atChr_(NULL)
{
  pileup_read = false;
}

void SNPinGenome::setSNPChr(std::vector<SNPatChr>* SNP_atChr_)
{
  this->SNP_atChr_ = new std::vector<SNPatChr>(*SNP_atChr_);
}

int SNPinGenome::processSNPLine(bool isVCF, char * line, string & myChr, int & index,int &previousPos) {

//chr1	11273	C/G	rs72481019

//#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
//1	13372	.	G	C	608.91

    #define MAX_COLS 32

    char* strs[MAX_COLS];
    unsigned int strs_cnt = split(line, '\t', strs);
    if (strs_cnt>=3) {
        string chr = strs[0];
        processChrName(chr);
        if (chr.compare(myChr)!=0) {
            index++;
            SNP_atChr_->push_back(SNPatChr(chr));
            myChr = chr;
        }
        int position = atoi(strs[1]);
        if (position!=previousPos) {
            previousPos = position;
            if (int(SNP_atChr_->size())>index) {
                if (!isVCF) {
                    if (strlen(strs[4])==1) {
                        (*SNP_atChr_)[index].push_SNP(SNPposition(position,strs[2],strs[3],strs[4]));
                        return 1;
                    }
                    return 0;
                } else if (strlen(strs[3])==1 && strlen(strs[4])<10) {
                        (*SNP_atChr_)[index].push_SNP(SNPposition(position,strs[4])); //if VCF
                        return 1;
                } else {return 0;}
            } else {
                cerr << "something is wrong with reading SNP positions"<< endl;
                return 0;
            }

        } else { return 0;}
    }
    return 0;
}


void SNPinGenome::readSNPs(std::string const& inFile)
{
    cout << "..Starting reading "<< inFile << " to get SNP positions" << std::endl;

	std::ifstream fileSNP (inFile.c_str());

    if (!fileSNP.is_open()) {
	    cerr << "Error: unable to open "+inFile+"\n" ;
	    exit(-1);
	}
    SNP_atChr_ = new std::vector <SNPatChr>();

    int count = 0;
    string line;

    int index = 0;
    string myChr = "1";
    SNP_atChr_->push_back(SNPatChr("1"));

    int previousPos = NA;

#ifdef PROFILE_TRACE
	time_t t0 = time(NULL);
#endif

	/*if (makingpileup != true)*/
	{
        //check whether the file is compressed:
        bool ifGZ = 0;
        if (inFile.substr(inFile.size()-3,3).compare(".gz")==0){ifGZ=1;}
        //check whether the file is in VCF format:
        bool ifVCF = 0;
        std::size_t found = inFile.find(".vcf");
        if (found!=std::string::npos) {ifVCF=1;}

        if (ifGZ) {
            fileSNP.close();
            FILE *stream;
            char buffer[MAX_BUFFER];
            string command = "gzip -cd "+inFile;
            stream =
            #if defined(_WIN32)
                      _popen(command.c_str(), "r");
            #else
                    popen(command.c_str(), "r");
            #endif
            char *line_buffer;
            while ((line_buffer = getLine(buffer, MAX_BUFFER, stream, line)) != NULL) {
              if (line_buffer[0] == '#') continue;
              count+=processSNPLine(ifVCF,line_buffer,myChr,index,previousPos);
            }
            #if defined(_WIN32)
                    _pclose(stream);
            #else
                    pclose(stream);
            #endif
        } else {
            while (std::getline(fileSNP,line)) {
                if ((!line.length()) || (line[0] == '#')) continue;
                count+=processSNPLine(ifVCF,(char*)line.c_str(),myChr,index,previousPos);
            }
        }
    } /*else
    {
        while (std::getline(fileSNP,line)) {
            if (! line.length()) continue;
//chr1	11273	C/G	rs72481019
            unsigned int strs_cnt = split((char*)line.c_str(), '\t', strs);
            if (line[0] == '#')
                {
                strs_cnt = 0;
                }
            if (strs_cnt>=3) {
                string chr = strs[0];
                processChrName(chr);
                if (chr.compare(myChr)!=0) {
                    index++;
                    SNP_atChr_->push_back(SNPatChr(chr));
                    myChr = chr;

                }
                int position = atoi(strs[1]);
                if (position!=previousPos) {
                    (*SNP_atChr_)[index].push_SNP(SNPposition(position,strs[4],"+",strs[3]));
                    count++;
                    previousPos = position;
                }
            }
        }
    }*/

    cout << "..read " << count <<" SNP positions\n";
	fileSNP.close();

#ifdef PROFILE_TRACE
	std::cout << "PROFILING [tid=" << pthread_self() << "]: " << inFile << " read in " << (time(NULL)-t0) << " seconds [readSNPs]\n" << std::flush;
#endif

}

SNPinGenome::~SNPinGenome()
{
    //dtor
}

void SNPinGenome::readMateFile(std::string const& mateFile, std::string const& inputFormat, int minimalTotalLetterCountPerPosition, int minimalQualityPerPosition) {
  assignValues(mateFile, inputFormat, minimalTotalLetterCountPerPosition,minimalQualityPerPosition);
  pileup_read = true;
}

void SNPinGenome::readMateFile(std::string const& mateFile, std::string const& inputFormat, int minimalTotalLetterCountPerPosition, int minimalQualityPerPosition, GenomeCopyNumber& genomeCopyNumber, std::string const& chrLenFileName, int windowSize, int step,  std::string targetBed) {
  // must perform partly GenomeCopyNumber::readCopyNumber line #114, all but fillMyHash
    genomeCopyNumber.initCopyNumber(chrLenFileName, windowSize, step, targetBed);
    assignValues(mateFile, inputFormat, minimalTotalLetterCountPerPosition,minimalQualityPerPosition, &genomeCopyNumber);

  pileup_read = true;
}

void SNPinGenome::setBinAt(int indexSNP,int index,int bin) {
    (*SNP_atChr_)[indexSNP].setBinAt(index, bin);
}

void SNPinGenome::setWESanalysis(bool WESanalysis) {
    WESanalysis_=WESanalysis;
}
void SNPinGenome::setCopyNumberFromPileup(bool CopyNumberFromPileup) {
    CopyNumberFromPileup_=CopyNumberFromPileup;
}

void SNPinGenome::perform(std::string const& mateFile, const std::string& inputFormat, int minimalTotalLetterCountPerPosition,
 int minimalQualityPerPosition, bool noisyData,bool CompleteGenomicsData,GenomeCopyNumber& genomeCopyNumber, double breakPointThreshold,
 int breakPointType, int minCNAlength, const char* what)
{
    if (!pileup_read) {
	  assignValues(mateFile, inputFormat, minimalTotalLetterCountPerPosition,minimalQualityPerPosition);
	}

	if (genomeCopyNumber.ifHasRatio()) {
        cout << "..Adding BAF info to the " << what << " dataset" << std::endl;
        genomeCopyNumber.addBAFinfo(*this);
        cout << "..Recalculate breakpoints using BAF profiles" << std::endl;
        genomeCopyNumber.calculateBAFBreakpoints(breakPointThreshold,breakPointType);
        cout << "..Recalculate median values" << std::endl;
        genomeCopyNumber.calculateCopyNumberMedians(minCNAlength, noisyData, CompleteGenomicsData);

        cout << "..Reannotate copy numbers" << std::endl;
        if (genomeCopyNumber.getWESanalysis() == false)
            {genomeCopyNumber.calculateCopyNumberProbs_and_genomeLength(breakPointType);}
        else
            {genomeCopyNumber.calculateCopyNumberProbs_and_exomeLength(breakPointType);}

    }

	cout << "..Done" << std::endl;
}

long SNPinGenome::processPileUPLine(int & positionCount, char* line, string & oldChr, int & sNPpositionToProceed,int minimalTotalLetterCountPerPosition,int & index, int minimalQualityPerPosition, GenomeCopyNumber* p_genomeCopyNumber) {

    if (*line == 0) return 0;
    if (line[0] == '#') return 0;
    if ( line[0] == '@') return 0;

    //we will count letters after . and , in corresponding positions:
    //std::vector<char*> strs(12);
    char* strs[64];
	unsigned int strs_cnt = split(line, '\t', strs);
    if (strs_cnt <= 4) {
        return 0;
    }
    //string chr = strs[0];
    //if (chr.compare(oldChr)!=0) {
	if (strcmp(strs[0], oldChr.c_str()) != 0) {
	    string chr = strs[0];
	    oldChr = chr;
	    oldChr = strs[0];
		processChrName(chr);
		index = findIndex(chr);
		positionCount = 0;
		if ( index == NA) {
		  cout << "will skip chr"<<chr<<"\n";
		  return 0;
		}
		sNPpositionToProceed = (*SNP_atChr_)[index].getPositionAt(positionCount);
    }

    if ( index == NA) {
            return 0;
    }

    int currentPosition = atoi(strs[1]);
    long valueToReturn = 0;
    if (CopyNumberFromPileup_ == true)     {
        if (WESanalysis_ == false)     {
            if (p_genomeCopyNumber) {
                string chr = strs[0];
                processChrName(chr);
                int lindex = p_genomeCopyNumber->findIndex(chr);
                if (lindex != NA) {
                    if (valueToReturn = strccnt(strs[4], '^')) {
                        ChrCopyNumber& chrCopyNumber = p_genomeCopyNumber->getChrCopyNumberAt(lindex);
                        int step = p_genomeCopyNumber->getStep();
                        for (int i=0; i<valueToReturn; i++)
                            chrCopyNumber.mappedPlusOneAtI(currentPosition, step);
                        }
                }
            }
        } else {
            if (p_genomeCopyNumber) {
                string chr = strs[0];
                 processChrName(chr);
                int lindex = p_genomeCopyNumber->findIndex(chr);
                if (lindex != NA) {
                    if (valueToReturn = strccnt(strs[4], '^')) {
                        ChrCopyNumber& chrCopyNumber = p_genomeCopyNumber->getChrCopyNumberAt(lindex);
                        int l = 0;
                        bool leftIsInTheWindow = false;
                        if ((currentPosition - 1 < chrCopyNumber.getEndAtBin(l) && (currentPosition > (chrCopyNumber.getCoordinateAtBin(l)  ))))   {
                            leftIsInTheWindow = true;
                        }  else    {
                            while ((!((currentPosition - 1 < chrCopyNumber.getEndAtBin(l)) &&
                                (currentPosition > (chrCopyNumber.getCoordinateAtBin(l) ))))  &&
                                (l <  chrCopyNumber.getExons_Countchr()))    {
                                leftIsInTheWindow = false;
                                l++;
                            }   if (l < chrCopyNumber.getExons_Countchr() -1 )    {
                                leftIsInTheWindow = true;
                            }
                        }
                        if (leftIsInTheWindow == true)  {
                            int step = p_genomeCopyNumber->getStep();
                            for (int i=0; i<valueToReturn; i++)
                                chrCopyNumber.mappedPlusOneAtI(currentPosition,step, l);
                        }
                    }
                }
            }
        }
    }

    if (currentPosition==sNPpositionToProceed) {
        try {
            float localBAF=addInfoFromAPileUp(atoi(strs[3]),minimalTotalLetterCountPerPosition,(*SNP_atChr_)[index].getNucleotideAt(positionCount),
                                           index,positionCount,sNPpositionToProceed,strs[4], minimalQualityPerPosition,strs[5]);
//                    if (localBAF>=0) {
//                        heterozygousBAFs.push_back(localBAF);
//                        heterozygousBAFs05.push_back(fabs(localBAF-0.5));
//                        heterozygousBAFposition.push_back(sNPpositionToProceed);
//                    }
        } catch (const char * error){
            //will simply forget about this SNP position..
        }
    }
    if (currentPosition>sNPpositionToProceed && sNPpositionToProceed != NA) {
	    while (currentPosition>sNPpositionToProceed && sNPpositionToProceed != NA) {
                (*SNP_atChr_)[index].setValueAt(positionCount,NA);
                positionCount++;
                if (positionCount<(*SNP_atChr_)[index].getSize())
                    sNPpositionToProceed = (*SNP_atChr_)[index].getPositionAt(positionCount);
                else
                    sNPpositionToProceed=NA;
        }
        if (currentPosition==sNPpositionToProceed) {
                 try {
                    float localBAF=addInfoFromAPileUp(atoi(strs[3]),minimalTotalLetterCountPerPosition,(*SNP_atChr_)[index].getNucleotideAt(positionCount),
                                           index,positionCount,sNPpositionToProceed,strs[4], minimalQualityPerPosition,strs[5]);
//                        if (localBAF>=0) {
//                            heterozygousBAFs.push_back(localBAF);
//                            heterozygousBAFs05.push_back(fabs(localBAF-0.5));
//                            heterozygousBAFposition.push_back(sNPpositionToProceed);
//
//                        }
                } catch (const char * error){
                    //will simply forget about this SNP position..
                }
        }
    }
	return valueToReturn;
}

void SNPinGenome::readPileUP(FILE* stream, int minimalTotalLetterCountPerPosition, int minimalQualityPerPosition, GenomeCopyNumber* p_genomeCopyNumber) {
    string oldChr = "azeaze";
    int sNPpositionToProceed;
    int positionCount = 0;
    string line;
    int index=NA;
	char* line_buffer = NULL;
	char buffer[MAX_BUFFER];
	long normalCount = 0;
	long count = 0;

	while ((line_buffer = getLine(buffer, MAX_BUFFER, stream, line)) != NULL) {
	  normalCount += processPileUPLine(positionCount, line_buffer, oldChr, sNPpositionToProceed,minimalTotalLetterCountPerPosition,index,minimalQualityPerPosition, p_genomeCopyNumber);
	  count++;
	}

	if (p_genomeCopyNumber) {
	  p_genomeCopyNumber->finishCopyNumber(normalCount);
	}
	std::cout << count << " lines read\n";
}

void SNPinGenome::assignValues(std::string const& inFile, string inputFormat, int minimalTotalLetterCountPerPosition, int minimalQualityPerPosition, GenomeCopyNumber* p_genomeCopyNumber) {
//    vector <float> heterozygousBAFs;
//    vector <float> heterozygousBAFs05;
//    vector <int> heterozygousBAFposition;

    if ((inputFormat.compare("pileup")==0 || inputFormat.compare("SAMtools pileup")==0)){
        //cout << "..will not proceed with BAF profiles: wrong inputFormat \n";
        cout << "..use \"pileup\" format of reads to calculate BAF profiles\n";

        cout << "..Starting reading "<< inFile << " to calculate BAF profiles" << std::endl;
#ifdef PROFILE_TRACE
		time_t t0 = time(NULL);
#endif

        if (inFile.substr(inFile.size()-3,3).compare(".gz")==0) {
		  string command = "gzip -cd "+inFile;
		  FILE* stream =
            #if defined(_WIN32)
				_popen(command.c_str(), "r");
			#else
				popen(command.c_str(), "r");
			#endif
				readPileUP(stream, minimalTotalLetterCountPerPosition, minimalQualityPerPosition, p_genomeCopyNumber);

            #if defined(_WIN32)
                    _pclose(stream);
            #else
                    pclose(stream);
            #endif

        } else {
		    FILE *stream = fopen(inFile.c_str(), "r");
			if (!stream) {
			  cerr << "Error: unable to open " + inFile + "\n" ;
			  exit(-1);
			}
			readPileUP(stream, minimalTotalLetterCountPerPosition, minimalQualityPerPosition,  p_genomeCopyNumber);
			fclose(stream);
        }

#ifdef PROFILE_TRACE
		std::cout << "PROFILING [tid=" << pthread_self() << "]: " << inFile << " read in " << (time(NULL)-t0) << " seconds [assignValues]\n" << std::flush;
#endif


        //BEGIN TMP:
//        std::vector<int> bpfinal;
//        int tmp = calculateBreakpoints_general(-0.0001, heterozygousBAFs05.size(), heterozygousBAFs05, bpfinal, 0, 4);
//        bpfinal.push_back(heterozygousBAFs05.size());
//        std::vector<float> medians;
//        std::vector<float> medianBAFs;
//
//        int start = 0;
//        for (int i = 0; i<bpfinal.size(); i++) {
//            float localMedian = get_median(heterozygousBAFs05, start, bpfinal[i]);
//            medians.push_back(localMedian);
//            for (int j = start; j<bpfinal[i];j++) {
//                medianBAFs.push_back(localMedian+0.5);
//            }
//            start = bpfinal[i];
//        }
//        std::ofstream file ("/bioinfo/users/vboeva/Desktop/Neuroblastome/ANALYSES_CNG/Valentina_analysis/freec/rowBAFS.txt");
//        for (int i = 0; i< heterozygousBAFs.size(); i++)
//            file << heterozygousBAFs[i]<<"\t"<<medianBAFs[i] <<"\t"<<heterozygousBAFs05[i]<<"\t"<<heterozygousBAFposition[i]<<"\n";
//        file.close();
//        heterozygousBAFs.clear();
//        medianBAFs.clear();
//        heterozygousBAFs05.clear();
//        heterozygousBAFposition.clear();
         //END TMP:
    } else {//if (inputFormat.compare("eland")==0 || inputFormat.compare("Eland")==0) {
        cerr << "\nWrong data format: "<< inputFormat << "\n Currently, Control-FREEC only accepts SAM pileup files to calculate BAF profiles\n";
        exit(1);
    }


}

int SNPinGenome::findIndex (string chromosome) const {
    for (int unsigned i = 0; i < SNP_atChr_->size(); i++) {
        if (chromosome.compare((*SNP_atChr_)[i].getChromosome())==0) {
            return i;
        }
    }
    return NA;
}

const SNPatChr& SNPinGenome::SNP_atChr(int index) const {
    return (*SNP_atChr_)[index];
}

SNPatChr& SNPinGenome::SNP_atChr(int index) {
    return (*SNP_atChr_)[index];
}

float SNPinGenome::addInfoFromAPileUp (int totalLetterCount, int minimalTotalLetterCountPerPosition,char whatToLook, int index,int &positionCount, int &sNPpositionToProceed, const char * pileup,int minimalQualityPerPosition, const char * quality){
    float value=NA;
    float status = NA;
    string pileupShort = "";
    if (totalLetterCount>=minimalTotalLetterCountPerPosition) {
        //check that there is not deletions/insertions at this position:
        if (strccnt(pileup, int('+'))+strccnt(pileup, int('-'))==0) {


            if (minimalQualityPerPosition>0) {

                pileupShort=pileup;
          //  cout << quality << "\t"<< pileupShort <<"\t";
                string qualityS = quality;
                chomp(qualityS);
                deleteChar(pileupShort, '^', 1);
                strkeepOnly(pileupShort, ".,ACGTNacgtn*");
         //   cout << pileupShort <<"\n";
                if (pileupShort.length()==qualityS.length()) {
                    filterWithQualities(pileupShort,qualityS,minimalQualityPerPosition);
                    deleteChar(pileupShort, '*');
                    pileup = pileupShort.c_str();
                    totalLetterCount = pileupShort.length();
                } else {
                    cout << "position to skip: " << pileupShort << " " << qualityS << "\n";
                    //skip this position..
                    totalLetterCount = -1;
                }
            }

            if (totalLetterCount>=minimalTotalLetterCountPerPosition) {

                //we will recalculate totalLetterCount = refCount+countForOverLetter
                int refCount =strccnt(pileup, '.')+strccnt(pileup, ',');
                int countForOverLetter=strccnt(pileup, whatToLook)+strccnt(pileup, tolower(whatToLook));

                //value = countForOverLetter*1./totalLetterCount;


                if (refCount+countForOverLetter>=minimalTotalLetterCountPerPosition) {
                    if (countForOverLetter+refCount>totalLetterCount) {
                            countForOverLetter=totalLetterCount-refCount;
                            if (countForOverLetter<0)
                                countForOverLetter=0;
                    }
                    value = countForOverLetter*1./(refCount+countForOverLetter);

                    //value = fabs(countForOverLetter*1./totalLetterCount-0.5);
                  //  cout << "pileup: "<< pileup << " ; whatToLook: "<< whatToLook<< "\n";
                   // cout << "countForOverLetter-1: "<< countForOverLetter-1 << "; totalLetterCount: "<< totalLetterCount << "; ERROR_PER_POS: "<< ERROR_PER_POS << "\n";
                    double prob = binomialcdistribution(countForOverLetter-1, totalLetterCount, ERROR_PER_POS); //probability that there > k Beta Allel bases if the SNP is HomoZ.
                    status = 0;
                    if (prob<0.05) { //0.01

                        if (binomialcdistribution(totalLetterCount-countForOverLetter-1, totalLetterCount, ERROR_PER_POS)<0.01)
                            //decide that SNP is not homoZ.
                            status = 0.5; //at this point, to indicate possible AB
                    }
                }
            }
        }
    }
    (*SNP_atChr_)[index].setValueAt(positionCount,value);
    (*SNP_atChr_)[index].setStatusAt(positionCount,status);

    positionCount++;
    if (positionCount<(*SNP_atChr_)[index].getSize())
        sNPpositionToProceed = (*SNP_atChr_)[index].getPositionAt(positionCount);
    else
        sNPpositionToProceed=NA;
    if (status>0)
        return value;
    return NA;
}

void* SNPinGenome_perform_wrapper(void *arg)
{
  SNPinGenomePerformArgWrapper* warg = (SNPinGenomePerformArgWrapper*)arg;
  warg->snpingenome.perform(warg->mateFile, warg->inputFormat, warg->minimalTotalLetterCountPerPosition, warg->minimalQualityPerPosition, warg->noisyData,warg->CompleteGenomicsData,warg->genomeCopyNumber, warg->breakPointThreshold, warg->breakPointType, warg->minCNAlength, warg->what);
  return NULL;
}


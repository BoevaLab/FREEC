#ifndef HEADER_3CEFE72622D637A0
#define HEADER_3CEFE72622D637A0

#pragma once
#ifndef SVFINDER_H
#define SVFINDER_H


//#include <errno.h>
//#include <stdlib.h>

//#include <time.h>
//#include <math.h>
#include <iostream>
#include <stdio.h>
#include <string.h>

#ifdef _WIN32
//x32 Windows definitions
#include "io.h"   // For access() in Windows.
#endif

#include <unistd.h>  // for access() in eveaLinux
#include <sys/types.h>  // For stat().
#include <sys/stat.h>   // For stat().


#include "myFunc.h"



//using namespace std ;

#define MAX_ALPHABET 4
#define NA -1
#define FileNameLength 281 //maximal length of file name
#define MaxWordLength 281 // LLIB
#define Infinity 1000000000000
#define GOODSD 0.5
#define TELO_CENTRO_FLANCS 50000
//#define GOODMAPPABILITY 0.85

#define MAXUncertainty 0.5

#define false 0
#define true 1


typedef int SYMB;
typedef long double ldouble;
extern int verbose;
extern char *ULetters[90];
extern double minMappabilityPerWindow;
extern bool uniqueMatch;

//extern double PI = 3.141592;

//#include "GetDensity.h"
#include "ConfigFile.h"
#include "Chameleon.h"
#include "Help.h"
#include "GenomeDensity.h"
#include "myFunc.h"
#include "SNPinGenome.h"
//#include "GenomeCopyNumber.h"


int runWithDefinedPloidy(int ploidy, GenomeCopyNumber & sampleCopyNumber, GenomeCopyNumber & controlCopyNumber, bool isControlIsPresent, int forceGC,
        bool has_BAF,bool ifTargeted,bool WESanalysis,
        int degree,int intercept,bool logLogNorm,float minExpectedGC,float maxExpectedGC,float knownContamination,float breakPointThreshold,int breakPointType,int minCNAlength,
        int teloCentroFlanks, std::vector<double> & RSS, std::vector<double> &percentage_GenExpl,bool contaminationAdjustment,std::vector<double> &contamination,
        ThreadPool * thrPool,ThreadPoolManager * thrPoolManager,std::string makePileup, float seekSubclones,
         std::string myName,std::vector<int> &unexplainedChromosomes, bool CompleteGenomicsData,bool normalization) ;


#endif //SVFINDER_H

#endif // header guard

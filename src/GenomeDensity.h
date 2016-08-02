#ifndef HEADER_569927DB13EA570C
#define HEADER_569927DB13EA570C

#pragma once
#ifndef __GENOME_DENSITY_H__
#define __GENOME_DENSITY_H__


#include <string>
#include <sstream>
#include <iostream>
#include <map>
#include <fstream>
#include <vector>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>

#include "KernelVector.h"
#include "ChrDensity.h"
#include "myFunc.h"


class GenomeDensity
{
public:
	GenomeDensity();
	GenomeDensity(std::string const& mateFile ,std::string const& inputFormat, std::string const& matesOrientation, std::string const& chrLenFile );
	~GenomeDensity(void);
	void smooth();
	void printDensity(std::string key,int i,int j);
	void calculateLogRatio(GenomeDensity controlGD, std::string type);
	void writeToFile(std::string outFileName,std::string type, bool binary);
	int getTotalNumberOfPairs();
private:
	std::map<std::string,ChrDensity> gd_;
	double mu_;
	double sigma_;
	long totalNumberOfPairs_;
};


#endif

#endif // header guard 

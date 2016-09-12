#ifndef HEADER_EE7515834790E569
#define HEADER_EE7515834790E569

#pragma once
#ifndef _ENTRY_CNV_H
#define _ENTRY_CNV_H

#include <iostream>
#include <string>
#include <sstream>

#include "SVfinder.h"

class EntryCNV
{
public:
	EntryCNV(std::string chr, int start, int end, int startCoord, int endCoord, int copyNumber);
    EntryCNV(std::string chr, int start, int end, int startCoord, int endCoord, int copyNumber, float estimatedBAFuncertainty, std::string medianBAFSymbol, bool hasBAF);

	~EntryCNV(void);
	std::string getChr();
	int getStartCoord();
	int getEndCoord();
    int getStart();
	int getEnd();
	int getCopyNumber();
	std::string printEntryCNV(float plody);
	void setType (std::string type);
	void setGermlinePercent (float germlinePercent) ;
    float compare(EntryCNV c2, int overlapPrecision, int & left, int & right, float ploidy, float controlPloidy); //return overlap


private:
	std::string chr_;
	int start_;
	int end_;
	int startCoord_;
	int endCoord_;
	int copyNumber_;
	float estimatedBAFuncertainty_;
	std::string medianBAFSymbol_;
	std::string type_; //somatic or germline
	float germlinePercent_; //percent of germline CNA or LOH
	bool isBAFassessed_;
};
#endif

#endif // header guard

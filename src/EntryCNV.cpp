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


#include "EntryCNV.h"

using namespace std ;

EntryCNV::EntryCNV(string chr, int start, int end, int startCoord, int endCoord, int copyNumber) {
	chr_ = chr;
	start_ = start;
	end_ = end;
	startCoord_ = startCoord;
	endCoord_ = endCoord;
	copyNumber_ = copyNumber;
	estimatedBAFuncertainty_ = NA;
	medianBAFSymbol_ = "";
	type_="";
	isBAFassessed_=0;
}

EntryCNV::EntryCNV(string chr, int start, int end, int startCoord, int endCoord, int copyNumber,float estimatedBAFuncertainty, std::string medianBAFSymbol, bool hasBAF) {
	chr_ = chr;
	start_ = start;
	end_ = end;
	startCoord_ = startCoord;
	endCoord_ = endCoord;
	copyNumber_ = copyNumber;
	estimatedBAFuncertainty_=estimatedBAFuncertainty;
	medianBAFSymbol_=medianBAFSymbol;
	type_="";
	if (hasBAF && copyNumber==0) {
        estimatedBAFuncertainty_=-1;
        medianBAFSymbol_="-";
	}
    isBAFassessed_=1;
}

std::string EntryCNV::printEntryCNV(float ploidy){

    std::stringstream ss;

    ss << chr_ << "\t"<< startCoord_<< "\t"<< endCoord_<<"\t"<<copyNumber_<<"\t";
    if (copyNumber_>ploidy)
        ss << "gain";
    else if (copyNumber_<ploidy)
        ss << "loss";
    else
        ss << "neutral";

    if (isBAFassessed_!=0)
             ss << "\t"<< medianBAFSymbol_<< "\t"<< estimatedBAFuncertainty_;

    if (type_ != "" && type_.compare("normal")!=0)
            ss << "\t"<< type_<< "\t"<< germlinePercent_;
    if (type_.compare("normal")==0)
            ss << "\t"<< "neutral"<< "\t"<< germlinePercent_;
	return ss.str();
}

std::string EntryCNV::getChr(){
	return chr_;
}
int EntryCNV::getStartCoord(){
	return startCoord_;
}
int EntryCNV::getEndCoord(){
	return endCoord_;
}
int EntryCNV::getStart(){
	return start_;
}
int EntryCNV::getEnd(){
	return end_;
}
int EntryCNV::getCopyNumber(){
	return copyNumber_;
}

EntryCNV::~EntryCNV(void)
{
}

void EntryCNV::setGermlinePercent (float germlinePercent) {
    germlinePercent_=germlinePercent;
}

void EntryCNV::setType (string type) {
    type_=type;
}

float EntryCNV::compare(EntryCNV c2, int overlapPrecision, int & left, int & right, float ploidy, float controlPloidy) { //return overlap
    if (chr_.compare(c2.getChr())!=0)
        return 0;
    if ((copyNumber_>ploidy && c2.getCopyNumber()>controlPloidy) || (copyNumber_<ploidy && c2.getCopyNumber()<controlPloidy) || (copyNumber_==ploidy && c2.getCopyNumber()==controlPloidy)) {
        int startControl = c2.getStart()-overlapPrecision; //to get larger window in the control..
        int endControl = c2.getEnd()+overlapPrecision;

        if (endControl<start_ || startControl>end_)
            return 0;
        left = max(start_,startControl);
        right = min(end_,endControl);
        return (right-left+1)/(end_-start_+1);
    } else {
        return 0;
    }

}


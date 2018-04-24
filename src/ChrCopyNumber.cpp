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


#include "ChrCopyNumber.h"

using namespace std ;

ChrCopyNumber::ChrCopyNumber(void)
{
	isMedianCalculated_ = false;
	isSmoothed_ = false;
	normalContamination_=0;
}

ChrCopyNumber::ChrCopyNumber(std::string const& chrName) {
	chromosome_ = chrName;
	isMedianCalculated_ = false;
	isSmoothed_ = false;
	ploidy_=NA;
	normalContamination_=0;
}

ChrCopyNumber::ChrCopyNumber(int windowSize, int chrLength, std::string const& chrName) {
	windowSize_ = windowSize;
	chrLength_ = chrLength;
	chromosome_ = chrName;
	step_=windowSize;
	isMedianCalculated_ = false;
	isSmoothed_ = false;
	ploidy_=NA;
	if (windowSize ==0) {
        cerr << "Error: windowSize is set to Zero\n";
        exit(-1);
	}
	length_ = chrLength/windowSize+1;
	coordinates_ = vector<int>(length_);
	readCount_ = vector<float>(length_,0);
	for (int i = 0; i<length_; i++) {
		coordinates_[i] = i*windowSize;
	}
	normalContamination_=0;
}

ChrCopyNumber::ChrCopyNumber(int windowSize, int chrLength, std::string const& chrName, int step, std::string targetBed) {
	windowSize_ = windowSize;
	step_=step;
	chrLength_ = chrLength;
	chromosome_ = chrName;
	isMedianCalculated_ = false;
	isSmoothed_ = false;
	ploidy_=NA;
	float meanTargetRegionLength=0;
	if (targetBed == "")     {
        if (windowSize ==0) {
            cerr << "Error: windowSize is set to Zero\n";
            exit(-1);
        }
        length_ = chrLength/step+1;
        coordinates_ = vector<int>(length_);
        readCount_ = vector<float>(length_,0);
        if (step<windowSize) {
            ends_ = vector<int>(length_,0);
            for (int i = 0; i<length_; i++) {
                ends_[i] = i*step+windowSize_-1;
            }
        }
        for (int i = 0; i<length_; i++) {
            coordinates_[i] = i*step;
            }
    }   else   {
        std::string const& captureFile = targetBed ;
        ifstream file (captureFile.c_str());
        if (file.is_open())     {
            cout << "..Reading "<< captureFile << "\n";
            cout << "..Your file must be in .BED format, and it must be sorted\n";
            std::string line;

			// length_ should end up equal to the number of exons matching this chromosome.
			length_ = 0;
			size_t tried = 0;


            while (std::getline(file,line))   {

					// Avoid catching a carriage return if the file uses windows-style line endings
					if(line.size() != 0 && line[line.size() - 1] == '\r')
						line.erase(line.size() - 1);

					size_t offset = 0;
					// Treat 'chr1' and '1' equivalently
					if(line.substr(0, 3) == "chr")
						offset += 3;
					size_t chrstart = offset;
					advance_to(line, offset, '\t');

					if(offset == line.size())
						continue;

					if(line.substr(chrstart, offset - chrstart) == chromosome_) {

						++tried;

						std::string start, endstr;
						++offset;
						size_t startoff = offset;
						advance_to(line, offset, '\t');
						start = line.substr(startoff, offset - startoff);
						if(start.size() == 0)
							continue;

						coordinates_.push_back(atoi(start.c_str()));

						++offset;
						size_t endoff = offset;
						advance_to(line, offset, '\t');
						endstr = line.substr(endoff, offset - endoff);
						if(endstr.size() == 0)
							continue;

						ends_.push_back(atoi(endstr.c_str()));

						advance_to(line, offset, ':'); ++offset;
						advance_to(line, offset, ':'); ++offset;
						size_t geneoff = offset;
						advance_to(line, offset, '\t');

						if(offset != geneoff)
							genes_names.push_back(line.substr(geneoff, offset - geneoff));
						else {
							// Print chr:start-end if a name isn't supplied.
							std::ostringstream oss;
							oss << chromosome_ << ":" << start << "-" << endstr;
							oss.flush();
							genes_names.push_back(oss.str());
						}
                        meanTargetRegionLength+=atoi(endstr.c_str())-atoi(start.c_str());
						++length_;

					}

                }

            if(int(tried) != length_) {
                std::cerr << "Warning: skipped " << (tried - length_) << " lines due to formatting problems\n";
            }

            exons_Countchr_ = length_;
            meanTargetRegionLength/=length_;
			readCount_ = vector<float>(exons_Countchr_,0);

			cout << "Number of exons analysed in chromosome "<< chromosome_ << " : " << exons_Countchr_ << "\n";
            cout << "Average exon length in chromosome "<< chromosome_ << " : " << meanTargetRegionLength << "\n";
            if (meanTargetRegionLength <30) {
                cerr << "WARNING: check your file with targeted regions: the average length of targeted regions is unexpectedly short\n";
            }

        }else {
			std::cerr << "Failed to open " << captureFile << "\n";
			exit(1);
		}

    }


//            while (std::getline(file,line) && l < exons_Countchr)
//                {
//                if (chr_namestmp[j] == chromosome_)
//                    {
//                    int i = 0;
//                    while (line[i] != '\t')
//                        {
//                        chr_names[l] += line[i];
//                        i++;
//                        }
//                        i++;
//
//                    while (line[i] != '\t')
//                        {
//                        coordinatesTmp_[l] += line[i];
//                        i++;
//                        }
//                        i++;
//                    while (line[i] != '\t')
//                        {
//                        endsTmp_[l] += line[i];
//                        i++;
//                        }
//                        i++;
//                    while (line[i] != ':')
//                        {
//                        i++;
//                        }
//                        i++;
//                    while (line[i] != ':')
//                        {
//                        i++;
//                        }
//                        i++;
//                    while (line[i] != '\t')
//                        {
//                        genes_names[l] += line[i];
//                        i++;
//                        }
//                        i++;
//                    l++;
//                }
//                    j++;
//                }
//            }
//            else
//            {
//            cerr << "Error: Unable to open file "+captureFile+"\n";
//            exit(-1);
//            }
//        readCount_ = vector<float>(exons_Countchr,0);
//        coordinates_ = vector<int>(exons_Countchr);
//        ends_ = vector<int>(exons_Countchr,0);
//        copy_number_subc = vector<int>(exons_Countchr,0);
//        population_subc = vector<float>(exons_Countchr,0);
//            for (int i = 0; i<exons_Countchr; i++)
//            {
//            ends_[i] = atoi(endsTmp_[i].c_str()); //Chaque fin d'exons
//            }
//            for (int i = 0; i<exons_Countchr; i++)
//            {
//            coordinates_[i] = atoi(coordinatesTmp_[i].c_str()); //Chaque début d'exons
//            }
//        cout << "Number of exons analysed in chromosome "<< chromosome_ << " : " << exons_Countchr << "\n";
//        }
	normalContamination_=0;
}


void ChrCopyNumber::mappedPlusOneAtI(int i, int step, int l) {
  if (l == -1)  {
    int pos = i/step;
    if ((int)readCount_.size()<=pos) {
		//should not normally happen unless we are at the very end of file
		cout << "Reaching end of file for chr "<<chromosome_ <<", position " << i <<"\n";
		//readCount_.resize(pos+1);
		//length_ = pos+1;
	} else {
		readCount_[pos]++;
        while (step*(pos-1)+windowSize_>i && pos>=1) {
            readCount_[--pos]++;
        }
	}
  }  else  {
    int pos = l;
    if ((int)readCount_.size()<=pos)   {
        //should not normally happen unless we are at the very end of file
        cout << "Reaching end of file for chr "<<chromosome_ <<", position " << i <<"\n";
        //readCount_.resize(pos+1);
        //length_ = pos+1;
    }    else    {
        readCount_[pos]++;
        //while (step*(pos-1)+windowSize_>i && pos>=1)
          //{
          //readCount_[--pos]++;
          //}
    }
  }
}

void ChrCopyNumber::setValueAt(int i, float val) {
	readCount_[i] = val;
}
void ChrCopyNumber::setCN_subc(int i, int CN_subc)
{
    //cerr << copy_number_subc_[i]; //WHAT IS THIS OUTPUT, CARINO?
    copy_number_subc_[i] = CN_subc; //THIS VECTOR IS EMPTY!
}

int ChrCopyNumber::getCN_subc(int i)
{
    return copy_number_subc_[i];
}


void ChrCopyNumber::setPopulation_subc(int i, float pop_subc)
{
    population_subc_[i] = pop_subc;
}

float ChrCopyNumber::getPopulation_subc(int i)
{
    return population_subc_[i] ;
}


std::string ChrCopyNumber::getGeneNameAtBin(int i)
{
    if (i<int(genes_names.size()))
        return genes_names[i];
    return "";
}

float ChrCopyNumber::getValueAt(int i) {
	return readCount_[i];
}

int ChrCopyNumber::getCoordinateAtBin(int i) {
	return coordinates_[i];
}
int ChrCopyNumber::getEndAtBin(int i) {
	if ((int)ends_.size()>i)
        return ends_[i];
    else
        return coordinates_[i]+windowSize_-1;
}

int ChrCopyNumber::getExons_Countchr() {
	return exons_Countchr_;
}

void ChrCopyNumber::setNotNprofileAt(int i, float value) {
    if (i<int(notNprofile_.size()))
        notNprofile_[i] = value;
    else
        cout <<"out of boundaries!!\n";
}

void ChrCopyNumber::setMappabilityProfileAt(int i, float value) {
    mappabilityProfile_[i] = value;
}

void ChrCopyNumber::setBAFat(int i, float value) {
    BAF_[i] = value;
}

void ChrCopyNumber::addToReadCount(float f) {
	readCount_.push_back(f);
}
void ChrCopyNumber::addToCoordinates(int i) {
	coordinates_.push_back(i);
}

void ChrCopyNumber::addToEnds(int i) {
	ends_.push_back(i);
}


void ChrCopyNumber::setWindowSize(int windowSize) {
	windowSize_ = windowSize;
}

void ChrCopyNumber::setPloidy(int ploidy) {
    ploidy_=ploidy;
}

void ChrCopyNumber::setNormalContamination(float normalContamination) {
    normalContamination_=normalContamination;
}


void ChrCopyNumber::setStep(int step) {
	step_ = step;
}

void ChrCopyNumber::setVectorLength(int length){
	length_ = length;
}

void ChrCopyNumber::setCN_subcLength(int len)
{
    copy_number_subc_.clear();
    copy_number_subc_ = vector<int>(len,0);
}
void ChrCopyNumber::setpop_subcLength(int len)
{
    population_subc_.clear();
    population_subc_ = vector<float>(len,0);
}

void ChrCopyNumber::setChrLength(int chrLength) {
	chrLength_ = chrLength;
}

std::vector <float> ChrCopyNumber::getValues() {
	return readCount_;
}


void ChrCopyNumber::removeLowReadCountWindows(ChrCopyNumber control,const int RCThresh) {
    if (length_!=control.getLength()) {
        cerr << "Warning: control length is not equal to the sample length for chromosome " << chromosome_ << "\n";
        cerr << "Sample: " << length_ << " windows; control: "<< control.getLength()<<" windows\n";
        return;
    }
    for (int i = 0; i<length_; i++) {
		if (control.getValueAt(i) < RCThresh){
			readCount_[i]=NA;
			control.setValueAt(i,0);
		}
    }
}

void ChrCopyNumber::removeLowReadCountWindows(const int RCThresh) {
    for (int i = 0; i<length_; i++) {
		if (readCount_[i] < RCThresh){
			readCount_[i]=0;
		}
    }
}

void ChrCopyNumber::fillInRatio(bool islog) {
    if ((int)ratio_.size()!=length_)
		ratio_.resize(length_);

    for (int i = 0; i<length_; i++) {
        if (readCount_[i]>=0 && !(mappabilityProfile_.size() > 0 && mappabilityProfile_[i] <= minMappabilityPerWindow)) {
            if (islog) {
                ratio_[i] = log(readCount_[i]+1)/log(2.0);
            }else {
                ratio_[i] = readCount_[i];
            }
        } else {
            ratio_[i] = NA;
        }
    }
}

void ChrCopyNumber::calculateRatioLog(ChrCopyNumber control, const double * a, const int degree){
	if ((int)ratio_.size()!=length_)
		ratio_.resize(length_);
	for (int i = 0; i<length_; i++) {
		if ((control.getLength()>i)&&(control.getValueAt(i) != 0)){
            if (mappabilityProfile_.size() == 0 || mappabilityProfile_[i] > minMappabilityPerWindow) {
                ratio_[i] = readCount_[i]/polynomial(control.getValueAt(i),a,1,degree);
                if (ratio_[i]<0)
                    ratio_[i] = NA;
            } else
				ratio_[i] = NA;
		} else
			if (readCount_[i]==0)
				//ratio_[i] = 1;
				ratio_[i] = NA;
			else
				ratio_[i] = NA;
		//cout << readCount_[i] << "\t" << control.getValueAt(i) << "\t" << ratio_[i] << "\n";
	}
}

void ChrCopyNumber::calculateRatio(ChrCopyNumber control, const double * a, const int degree){
	if ((int)ratio_.size()!=length_)
		ratio_.resize(length_);
	for (int i = 0; i<length_; i++) {
		if ((control.getLength()>i)&&(control.getValueAt(i) != 0)){

			if (mappabilityProfile_.size() == 0 || mappabilityProfile_[i] > minMappabilityPerWindow) {
                ratio_[i] = readCount_[i]/polynomial(control.getValueAt(i),a,1,degree);
				if (ratio_[i]<0)
					ratio_[i] = NA;
			} else
                ratio_[i] = NA;

		} else
			if (readCount_[i]==0)
				//ratio_[i] = 1;
				ratio_[i] = NA;
			else
				ratio_[i] = NA;
		//cout << readCount_[i] << "\t" << control.getValueAt(i) << "\t" << ratio_[i] << "\n";
	}
}

void ChrCopyNumber::recalculateRatio(ChrCopyNumber control){
	for (int i = 0; i<length_; i++) {
	    float controlRatio = control.getRatioAtBin(i);
		if ((control.getLength()>i)&&(controlRatio > 0)){
			if (mappabilityProfile_.size() == 0 || mappabilityProfile_[i] > minMappabilityPerWindow) {
				ratio_[i] = ratio_[i]/controlRatio;
				if (ratio_[i]<0)
					ratio_[i] = NA;
			}	else
				ratio_[i] = NA;
		} else
			ratio_[i] = NA;
		//cout << readCount_[i] << "\t" << control.getValueAt(i) << "\t" << ratio_[i] << "\n";
	}
}

void ChrCopyNumber::calculateRatio(ChrCopyNumber control, double a0, double a1) {
	if ((int)ratio_.size()!=length_)
		ratio_.resize(length_);
	for (int i = 0; i<length_; i++) {
		if ((control.getLength()>i)&&(control.getValueAt(i) != 0) && (readCount_[i] > 0) && (control.getValueAt(i) > 0))
			if (mappabilityProfile_.size() == 0 || mappabilityProfile_[i] > minMappabilityPerWindow)
				ratio_[i] = float(readCount_[i]/(control.getValueAt(i)*a0+a1));
			else
				ratio_[i] = NA;
        else if ((readCount_[i] < 0) || (control.getValueAt(i) < 0))
            {
            ratio_[i] = NA;
            }
		else
			if (readCount_[i]==0)
				//ratio_[i] = 1;
				ratio_[i] = NA;
			else
				ratio_[i] = NA;
		//cout << readCount_[i] << "\t" << control.getValueAt(i) << "\t" << ratio_[i] << "\n";
	}
}

void ChrCopyNumber::calculateRatio(ChrCopyNumber control, float normalizationConst) {
	ratio_ = vector<float>(length_,0);
	for (int i = 0; i<length_; i++) {
		if ((control.getLength()>i)&&(control.getValueAt(i) != 0))
			if (mappabilityProfile_.size() == 0 || mappabilityProfile_[i] > minMappabilityPerWindow)
				ratio_[i] = readCount_[i]/control.getValueAt(i)/normalizationConst;
			else
				ratio_[i] = NA;
		else
			if (readCount_[i]==0)
				//ratio_[i] = 1;
				ratio_[i] = NA;
			else
				ratio_[i] = NA;
		//cout << readCount_[i] << "\t" << control.getValueAt(i) << "\t" << ratio_[i] << "\n";
	}
}



void ChrCopyNumber::setAllNormal() {
	if ((int)ratio_.size()!=length_)
		ratio_.resize(length_);
	for (int i = 0; i<length_; i++) {
		if (readCount_[i] != NA && i<int(notNprofile_.size()) &&notNprofile_[i]>0) {
            ratio_[i] = 1;
		} else {
			ratio_[i] = NA;
		}
	}
}

void ChrCopyNumber::recalculateRatio(double *a, int degree) {
	float x;
	for (int i = 0; i<length_; i++) {
		if (ratio_[i] != NA) {
			x = GCprofile_[i];

                //use a threshold, but correct using notN profile
                if (mappabilityProfile_.size()>0) {
                    if ((x>0)&&(mappabilityProfile_[i]>minMappabilityPerWindow)) //if ((x>0)&&(notNprofile_[i]!=0))
                        ratio_[i] = ratio_[i]/polynomial(x,a,1,degree);
                    else
                        ratio_[i] = NA;
                } else {
                    if ((x>0)&&(notNprofile_[i]>minMappabilityPerWindow)) //if ((x>0)&&(notNprofile_[i]!=0))
                        ratio_[i] = ratio_[i]/polynomial(x,a,1,degree);
                    else
                        ratio_[i] = NA;

                }


		}
		if ((ratio_[i] != NA)&&(ratio_[i] < 0))
			ratio_[i] = 0; //this happens if  polynomial(x,a,b,c,d) < 0
	}
}


void ChrCopyNumber::calculateRatio(double *a, int degree) {
	if ((int)ratio_.size()!=length_)
		ratio_.resize(length_);
	float x;
	for (int i = 0; i<length_; i++) {
		if (readCount_[i] != NA) {
			x = GCprofile_[i];

            // if uniqueMatch, do correction to mappability
			if (uniqueMatch) {
                //use mappabilityProfile_ and correct
                if ((x>0)&&(mappabilityProfile_[i]>minMappabilityPerWindow)) //if ((x>0)&&(notNprofile_[i]!=0))
                    ratio_[i] = readCount_[i]/polynomial(x,a,1,degree)/mappabilityProfile_[i];
                else
                    ratio_[i] = NA;

			} else {
                //use a threshold, but correct using notN profile
                if (mappabilityProfile_.size()>0) {
                    if ((x>0)&&(mappabilityProfile_[i]>minMappabilityPerWindow)) //if ((x>0)&&(notNprofile_[i]!=0))
                        ratio_[i] = readCount_[i]/polynomial(x,a,1,degree)/notNprofile_[i];
                    else
                        ratio_[i] = NA;
                } else {
                    if ((x>0)&&(notNprofile_[i]>minMappabilityPerWindow)) //if ((x>0)&&(notNprofile_[i]!=0))
                        ratio_[i] = readCount_[i]/polynomial(x,a,1,degree)/notNprofile_[i];
                    else
                        ratio_[i] = NA;

                }
			}

		} else {
			ratio_[i] = NA;
		}
		if ((ratio_[i] != NA)&&(ratio_[i] < 0))
			ratio_[i] = 0; //this happens if  polynomial(x,a,b,c,d) < 0
	}
}

//void ChrCopyNumber::calculateRatio(double a,double b, double c,double d) {
//	if ((int)ratio_.size()!=length_)
//		ratio_.resize(length_);
//	float x;
//	for (int i = 0; i<length_; i++) {
//		if (readCount_[i] != NA) {
//			x = GCprofile_[i];
//			if ((x>0)&&(notNprofile_[i]!=0))
//				ratio_[i] = readCount_[i]/polynomial(x,a,b,c,d)/notNprofile_[i];
//			else
//				ratio_[i] = NA;
//		} else {
//			ratio_[i] = NA;
//		}
//		if ((ratio_[i] != NA)&&(ratio_[i] < 0))
//			ratio_[i] = 0; //this happens if  polynomial(x,a,b,c,d) < 0
//	}
//}

bool ChrCopyNumber::ifHasRatio() {
    if (ratio_.size()>0) return 1;
    return 0;
}

void ChrCopyNumber::recalculateRatio (float constant) {
	for (int i = 0; i<length_; i++)
		if (ratio_[i] != NA)
			ratio_[i] /= constant;
}

void ChrCopyNumber::recalculateLogRatio (float constant) {
	for (int i = 0; i<length_; i++)
		if (ratio_[i] != NA)
			ratio_[i] -= constant;
}

void ChrCopyNumber::recalculateRatioWithContam (float contamination, float normGenytype, bool isLogged) { //normGenytype==1 if AB, normGenytype==0.5 if A
	if (!isLogged) {
        for (int i = 0; i<length_; i++)
            if (ratio_[i] != NA) {
			//ratio_[i] = (ratio_[i]-contamination*normGenytype)/(1-contamination); //correct only for ploidy 2
                ratio_[i] = (ratio_[i]*(1-contamination+2*contamination/ploidy_) -contamination*normGenytype/ploidy_*2)/(1-contamination);

			if (ratio_[i]<0)
				ratio_[i] = 0;
		}
    } else {
        for (int i = 0; i<length_; i++)
            if (ratio_[i] != NA) {
                float realCopy = pow(2,ratio_[i]);
                ratio_[i] = (realCopy*(1-contamination+2*contamination/ploidy_) -contamination*normGenytype/ploidy_*2)/(1-contamination);
                if (ratio_[i]<0)
                    ratio_[i] = NA;
                else {
                    ratio_[i]=log2(ratio_[i]);
                }

		}

    }
}


int ChrCopyNumber::calculateBreakpoints(double threshold, int normalChrLength, int breakPointType) {
    int chrLen = calculateBreakpoints_general(threshold,length_,ratio_,bpfinal_,normalChrLength,breakPointType, getChromosome());
 	return chrLen;
}

int ChrCopyNumber::calculateBAFBreakpoints(double threshold, int normalChrLength, int breakPointType) {
    // bpfinal_ should already contain copy number breakpoints
    std::vector <int> bpBAF;
 	//find breakpoints in the BAF profile
 	int chrLen = calculateBreakpoints_general(threshold,length_,BAF_,bpBAF,normalChrLength,breakPointType, getChromosome());
 	//add detected breakpoints to the breakpoints detected using copy number profiles
 	bpfinal_ = merge_no_dups(bpfinal_, bpBAF);
 	sort (bpfinal_.begin(), bpfinal_.end());
 	return chrLen;
}

int ChrCopyNumber::getCoveredPart(int breakPointStart, int breakPointEnd) { //for exome-seq: get length of the genome covered by the targeted region (from breakPointStart to breakPointEnd)
    int lengthCovered = 0;
    for (int i = breakPointStart; i<=breakPointEnd; i++) {
        lengthCovered+=ends_[i]-coordinates_[i]+1;
    }
    return lengthCovered;
}

void ChrCopyNumber::calculateCopyNumberMedian(int ploidy, int minCNAlength, bool noisyData,  bool CompleteGenomicsData, bool isLogged){ //create median profiles using 'bpfinal_' and store them in medianProfile_, info about medians themselves is stored in medianValues_ and about SD in sd_, lengths of fragments in bpLengths_
	if (ploidy!=ploidy_) {
        cerr << "..Warning: in calculateCopyNumberMedian() class's ploidy is different from "<<ploidy<<"\n";
        ploidy_=ploidy;
	}
	int breakPointStart = 0;
	int breakPointEnd;
	float median;

	bpfinal_.push_back(length_-1); //add last point

	vector <int> seg_ends;
	vector <int> seg_starts;

    if (int(medianProfile_.size())!=length_) {
        	medianProfile_ = vector <float> (length_,NA);
    }

	//for BAF:
	bool isBAFpresent = false;
	//float medianBAF;
	float estimatedBAF;
	float fittedBAF;
	float uncertainty;
	string medianBAFSym;
	if (BAF_.size()>0) {
        isBAFpresent = true;
        //medianBAFProfile_ = vector <float> (length_,NA);
        estimatedBAFProfile_ = vector <float> (length_,NA);
        fittedBAFProfile_ = vector <float> (length_,NA);
        medianBAFSymbol_ = vector <std::string> (length_,"-");
        estimatedBAFuncertainty_ = vector <float> (length_,NA);
        cout << "..Control: BAF profile is present\n";
	}

    //clear existing values:
    sd_.clear();
    medianValues_.clear();
    fragmentNotNA_lengths_.clear();
	fragment_lengths_.clear();
	BAFsymbPerFrag_.clear();
	estBAFuncertaintyPerFrag_.clear();

	for (int i  = 0; i < (int)bpfinal_.size();i++) {
		breakPointEnd = bpfinal_[i];
		//int ndatapoints = breakPointEnd-breakPointStart+1;
		vector<float> data;
        //vector<float> dataBAF;
		int notNA = 0;
		string BAFValuesInTheSegment = ""; //here we merge all SNP Values of this segment
		for (int j = breakPointStart; j <= breakPointEnd; j++)
		 {
			if (ratio_[j] != NA) {
				data.push_back(ratio_[j]);
				notNA++;
				if (isBAFpresent && BAF_[j]!=NA && BAFvalues_[j] != ""){
                    //dataBAF.push_back(BAF_[j]);
                    if (BAFValuesInTheSegment != ""){
                        if (BAFvalues_[j].compare(BAFvalues_[j-1])!=0) {
                            string suffix = findSmallestSuffix(BAFvalues_[j-1],BAFvalues_[j]);
                            if (suffix != "")
                                BAFValuesInTheSegment += ";"+suffix;
                                //BAFValuesInTheSegment += ";"+BAFvalues_[j];
                        }
                    }
                    else
                        BAFValuesInTheSegment = BAFvalues_[j];
                }
			}
         }
		int totalCount = breakPointEnd-breakPointStart+1;

        bool ifHomoz = false;
        float locMedian=NA;
        if (int(data.size())>=minCNAlength && data.size()>0) {
            locMedian = get_median(data); //including the last point
            if (isLogged)
                locMedian=pow(2, locMedian);
        }
        if (isBAFpresent && notNA > 100 && locMedian < 1 && noisyData ) {
            vector<string>heteroValuesPerWindowStrings = split(BAFValuesInTheSegment, ';');
            int numberofBAFpoints =heteroValuesPerWindowStrings.size();
            heteroValuesPerWindowStrings.clear();
            double threshold;
            if (step_ != 0) {
                threshold = 0.0001*step_*notNA;
            } else {
                threshold = 0.00001*getCoveredPart(breakPointStart, breakPointEnd);
            }

            if (numberofBAFpoints < threshold)
                ifHomoz = true;
            cout <<"chr"<<chromosome_ <<":"<<breakPointStart<<"-"<<breakPointEnd<<" totalCount=" << totalCount << " notNA="<<notNA<<" numberofBAFpoints="<<numberofBAFpoints<<" ratio="<< 1.*numberofBAFpoints/threshold<<" hom:"<<ifHomoz<<"\n";
        }
		//check is totalCount is notNA is greater than minCNAlength
		if (totalCount<minCNAlength && i<(int)bpfinal_.size()-1) { //merge this fragment with next one
			//do nothing, do to the next cycle, keep the same position of "breakPointStart"
		} else {
			//fragment_lengths_.push_back(totalCount);
			//if (notNA <= (totalCount+0.5)/2) { // no more than 1/2 of NA windows
			if (notNA==0 || notNA<sqrt(float(totalCount))){ //instead of (notNA == 1 && totalCount>2) meaning that only one point is not sufficient if there are more than 2 windows to set the median
				//median = get_median(data);
				median = NA;
				if (isBAFpresent) {
                    //medianBAF = NA;
                    estimatedBAF = NA;
                    fittedBAF=NA;
                    medianBAFSym = "-";
                    uncertainty = NA;
				}
			} else {
				median = get_median(data); //including the last point
				if (isLogged)
                    median=pow(2, median);
				if (isBAFpresent) {
//                    if (dataBAF.size()>0)
//                        medianBAF = get_median(dataBAF);
//                    else
//                        medianBAF = NA;
                    //medianBAF = NA;
                    getBAFinfo(BAFValuesInTheSegment,median*ploidy,estimatedBAF, fittedBAF, medianBAFSym,uncertainty, normalContamination_,ploidy_,noisyData,ifHomoz, CompleteGenomicsData);
				}
			}


			//check the previous median: may be one can merge the two regions?
			//if(totalCount==1 && median == NA && breakPointStart != 0){
			//	median = previousMedian;
			//	//join two segments
			//	medianValues_.pop_back();
			//	notNA += fragmentNotNA_lengths_.back();
			//	fragmentNotNA_lengths_.pop_back();
			//	//don't change the SDev
			//	breakPointStart = seg_starts.back();
			//	seg_ends.pop_back();
			//}
			//else if (round_by_ploidy(median,ploidy)==round_by_ploidy(previousMedian,ploidy)) {

			//	//join two segments
			//	medianValues_.pop_back();
			//	median = (median+previousMedian)/2;			//approximatelly true...
			//	notNA += fragmentNotNA_lengths_.back();
			//	fragmentNotNA_lengths_.pop_back();
			//	float sdLocal = (sd_.back()+sd(data,median))/2;  //approximatelly true...
			//	sd_.pop_back();
			//	sd_.push_back(sdLocal);
			//	breakPointStart = seg_starts.back();
			//	seg_ends.pop_back();
			//} else {
            sd_.push_back(sd(data,median));
            seg_starts.push_back(breakPointStart);
			//}


			seg_ends.push_back(breakPointEnd);
			for (int j = breakPointStart; j<= breakPointEnd; j++) {
                medianProfile_[j] = median;
                if (isBAFpresent) {
                        //medianBAFProfile_[j] = medianBAF;
                    medianBAFSymbol_[j] = medianBAFSym;
                    estimatedBAFProfile_[j]=estimatedBAF;
                    fittedBAFProfile_[j]=fittedBAF;
                    estimatedBAFuncertainty_[j]=uncertainty;
                }
			}
			medianValues_.push_back(median);
			fragmentNotNA_lengths_.push_back(notNA);
			if (isBAFpresent) {
                estBAFuncertaintyPerFrag_.push_back(uncertainty);
                BAFsymbPerFrag_.push_back(medianBAFSym);
               // cout << "..Control: adding "<<medianBAFSym<< " to a fragment with median*ploidy="<<median*ploidy<< "\n";

			}

			breakPointStart = breakPointEnd+1;
			data.clear();

			if (isBAFpresent) {
//                dataBAF.clear();
            }
			//previousMedian = median;

		}
	}

	bpfinal_.clear();
	bpfinal_ = seg_ends;
	seg_ends.clear();
	seg_starts.clear();
//recalc fragment_lengths_
	fragment_lengths_.push_back(bpfinal_[0]+1);
	for (int i = 1; i < int(bpfinal_.size()); i++)
		fragment_lengths_.push_back(bpfinal_[i]-bpfinal_[i-1]);

	bpfinal_.pop_back(); //delete last point which is (length_-1)
	isMedianCalculated_ = true;
}


//old version, only with ploidy
//void ChrCopyNumber::calculateCopyNumberMedian(int ploidy){ //create median profiles using 'bpfinal_' and store them in medianProfile_, info about medians themselves is stored in medianValues_ and about SD in sd_, lengths of fragments in bpLengths_
//	int breakPointStart = 0;
//	int breakPointEnd;
//	float median;
//	float previousMedian = -22; //something strange
//
//	bpfinal_.push_back(length_-1); //add last point
//
//	vector <int> seg_ends;
//	vector <int> seg_starts;
//
//	medianProfile_ = vector <float> (length_);
//	for (int i  = 0; i < (int)bpfinal_.size();i++) {
//		breakPointEnd = bpfinal_[i];
//		//int ndatapoints = breakPointEnd-breakPointStart+1;
//		vector<float> data;
//		int notNA = 0;
//		for (int j = breakPointStart; j <= breakPointEnd; j++)
//			if (ratio_[j] != NA) {
//				data.push_back(ratio_[j]);
//				notNA++;
//			}
//		int totalCount = breakPointEnd-breakPointStart+1;
//		//fragment_lengths_.push_back(totalCount);
//		if (notNA <= (totalCount+0.5)/3) { // no more than 2/3 of NA windows
//			median = NA;
//		} else {
//			median = get_median(data); //including the last point!!!! Ask Kevin or check if it is correct!!
//		}
//
//		//check the previous median: may be one can merge the two regions?
//
//		if(totalCount==1 && median == NA && breakPointStart != 0){
//			median = previousMedian;
//			//join two segments
//			medianValues_.pop_back();
//			notNA += fragmentNotNA_lengths_.back();
//			fragmentNotNA_lengths_.pop_back();
//			//don't change the SDev
//			breakPointStart = seg_starts.back();
//			seg_ends.pop_back();
//		}
//		else if (round_by_ploidy(median,ploidy)==round_by_ploidy(previousMedian,ploidy)) {
//
//			//join two segments
//			medianValues_.pop_back();
//			median = (median+previousMedian)/2;			//approximatelly true...
//			notNA += fragmentNotNA_lengths_.back();
//			fragmentNotNA_lengths_.pop_back();
//			float sdLocal = (sd_.back()+sd(data,median))/2;  //approximatelly true...
//			sd_.pop_back();
//			sd_.push_back(sdLocal);
//			breakPointStart = seg_starts.back();
//			seg_ends.pop_back();
//		} else {
//			sd_.push_back(sd(data,median));
//			seg_starts.push_back(breakPointStart);
//		}
//
//		seg_ends.push_back(breakPointEnd);
//		for (int j = breakPointStart; j<= breakPointEnd; j++)
//			medianProfile_[j] = median;
//		medianValues_.push_back(median);
//		fragmentNotNA_lengths_.push_back(notNA);
//
//		breakPointStart = breakPointEnd+1;
//		data.clear();
//		previousMedian = median;
//	}
//
//	bpfinal_.clear();
//	bpfinal_ = seg_ends;
//	seg_starts.clear();
////recalc fragment_lengths_
//	//fragment_lengths_.clear();
//	fragment_lengths_.push_back(bpfinal_[0]+1);
//	for (int i = 1; i < int(bpfinal_.size()); i++)
//		fragment_lengths_.push_back(bpfinal_[i]-bpfinal_[i-1]);
//
//	bpfinal_.pop_back(); //delete last point which is (length_-1)
//	isMedianCalculated_ = true;
//}

//old version
void ChrCopyNumber::calculateCopyNumberMedian(){ //create median profiles using 'bpfinal_' and store them in medianProfile_, info about medians themselves is stored in medianValues_ and about SD in sd_, lengths of fragments in bpLengths_

	int breakPointStart = 0;
	int breakPointEnd;
	float median;
	bpfinal_.push_back(length_-1); //add last point
	medianProfile_ = vector <float> (length_);
	for (int i  = 0; i < (int)bpfinal_.size();i++) {
		breakPointEnd = bpfinal_[i];
		//int ndatapoints = breakPointEnd-breakPointStart+1;
		vector<float> data;
		int notNA = 0;
		for (int j = breakPointStart; j <= breakPointEnd; j++)
			if (ratio_[j] != NA) {
				data.push_back(ratio_[j]);
				notNA++;
			}
		int totalCount = breakPointEnd-breakPointStart+1;
		fragment_lengths_.push_back(totalCount);
		if (notNA <= (totalCount+0.5)/2) {
			median = NA;
		} else {
			median = get_median(data); //including the last point!!!! Ask Kevin or check if it is correct!!
		}
		medianValues_.push_back(median);
		fragmentNotNA_lengths_.push_back(notNA);
		sd_.push_back(sd(data,median));
		for (int j = breakPointStart; j<= breakPointEnd; j++)
			medianProfile_[j] = median;

		breakPointStart = breakPointEnd+1;
		data.clear();
	}
	bpfinal_.pop_back(); //delete last point which is (length_-1)
	isMedianCalculated_ = true;
}

float ChrCopyNumber::getMedianProfileAtI (int i) {
	return medianProfile_[i];
}
float ChrCopyNumber::getMedianValuesAt (int i) {
	return medianValues_[i];
}

float ChrCopyNumber::getCGprofileAt(int i) {
	return GCprofile_[i];
}
float ChrCopyNumber::getNotNprofileAt(int i) {
	return notNprofile_[i];
}

float ChrCopyNumber::getMappabilityProfileAt(int i) {
	return mappabilityProfile_[i];
}

std::vector <float>	ChrCopyNumber::getRatio() {
	return ratio_;
}
bool ChrCopyNumber::isMedianCalculated() {
	return isMedianCalculated_;
}
bool ChrCopyNumber::isSmoothed(){
	return isSmoothed_;
}

void ChrCopyNumber::setIsSmoothed(bool value) {
	isSmoothed_ = value;
	smoothedProfile_.clear();
}
//void ChrCopyNumber::printLog2Ratio(std::ofstream const& file) {
//	/*for (int i = 0; i<length_; i++) {
//		file << ratio_[i] << "\n";
//	}*/
//}

int ChrCopyNumber::getLength() {
	return length_;
}
int ChrCopyNumber::getChrLength() {
	return chrLength_;
}


int ChrCopyNumber::getMappabilityLength() {
	return mappabilityProfile_.size();
}


std::string ChrCopyNumber::getChromosome() {
	return chromosome_;
}

float ChrCopyNumber::getRatioAtBin(int i) {
	return ratio_[i];
}

float ChrCopyNumber::getEstimatedBAFuncertaintyAtBin(int i) {
    return estBAFuncertaintyPerFrag_[i];
}

std::string ChrCopyNumber::getBAFsymbPerFrg (int i)  {
    return BAFsymbPerFrag_[i];
}

std::vector <float> ChrCopyNumber::getMedianValues () {
	return medianValues_;
}

std::vector <float> ChrCopyNumber::getSDs () {
	return sd_;
}
std::vector <int> ChrCopyNumber::getFragmentLengths_notNA () {
	return fragmentNotNA_lengths_;
}

std::vector <int> ChrCopyNumber::getFragmentLengths () {
	return fragment_lengths_;
}


int ChrCopyNumber::getFragmentLengthsAt (int i) {
	return fragment_lengths_[i];
}

int ChrCopyNumber::getFragmentLengths_notNA_At(int i) {
	return fragmentNotNA_lengths_[i];
}

std::vector <int> ChrCopyNumber::getBreakPoints() {
	return bpfinal_;
}

int ChrCopyNumber::getNumberOfFragments() {
	return bpfinal_.size()+1;
}

int ChrCopyNumber::nextNoNAIndex(int i1, int ploidy, int min_fragment) {
	for (int i = i1+1; i<(int)medianValues_.size(); i++)
		//if (round_by_ploidy(medianValues_[i], ploidy)>=0) // !!! >= and not > 0.
        if (getLevelAt(i, ploidy)>=0) // !!! >= and not > 0.
			if (fragment_lengths_[i] > min_fragment)
				return i;
	return NA;
}

float ChrCopyNumber::nextNoNAMedian(int i1, int ploidy) {

	for (int i = i1+1; i<(int)medianValues_.size(); i++)
		if (getLevelAt(i, ploidy)>0)
			return medianValues_[i];
	return NA;
}

int ChrCopyNumber::nextNoNALength(int i1, int ploidy) {

	for (int i = i1+1; i<(int)medianValues_.size(); i++)
		if (getLevelAt(i, ploidy)>0)
			return fragment_lengths_[i];
	return 0;
}

int ChrCopyNumber::getNumberOfGoodFragments() {
	int count = 0;
	for (int i = 0; i < (int)fragmentNotNA_lengths_.size(); i++)
		if (fragmentNotNA_lengths_[i] != 0)
			count++;
	return count;
}


double ChrCopyNumber::getXiSum(int ploidy, float minSD) {
	double sum = 0;
	const double PI = 3.141592;

	for (int i  = 0; i <= (int)bpfinal_.size();i++) {

		double theta_i = getLevelAt(i,ploidy);
		//bpLengths_[i];
		//sd_[i];
		double d = medianValues_[i]-theta_i;
		if (!(d == 0 && fragmentNotNA_lengths_[i] == 0)) {
			double sqrtVar = sqrt(float(PI/2./fragmentNotNA_lengths_[i]))*max(minSD,sd_[i]);
			double term = pow(d/sqrtVar,2);
			sum += term;
		}
	}


	return sum;
}

double ChrCopyNumber::calculateXiSum(int ploidy, map <float,float> &sds, map <float,float> &meds) {
	double sum = 0;
	const double PI = 3.141592;

	for (int i  = 0; i <= (int)bpfinal_.size();i++) {

		float level = getLevelAt(i,ploidy);
		if (level != NA) {
			double theta_i = meds.find(level)->second;
			//bpLengths_[i];
			//sd_[i];
			double d = medianValues_[i]-theta_i;
			if (fragmentNotNA_lengths_[i] != 0) {
				double sqrtVar = sqrt(float(PI/2./fragmentNotNA_lengths_[i]))*max(sds.find(level)->second,sd_[i]);
				double term = pow(d/sqrtVar,2);
				sum += term;
			}
		}
	}
		/*ofstream myfile;
		 myfile.open ("example.txt");
		 myfile <<"ratio\tmedianValues\n";
		for (int i  = 0; i < (int)ratio_.size();i++) {
			myfile << ratio_[i] << "\t" <<medianProfile_[i] <<"\n";
		}
		myfile.close();*/
	return sum;
}

double ChrCopyNumber::calculateXiSum(int ploidy, map <float,float> &sds) {
	double sum = 0;
	const double PI = 3.141592;

	for (int i  = 0; i <= (int)bpfinal_.size();i++) {

		float level = getLevelAt(i,ploidy);
		if (level != NA) {
			//bpLengths_[i];
			//sd_[i];
			double d = medianValues_[i]-level;
			if ((fragmentNotNA_lengths_[i] != 0)&&(sd_[i] != 0)) {
				float localSV = sd_[i];
				//float globalSV = sds.find(level)->second;
				//double sqrtVar = sqrt(float(PI/2./fragmentNotNA_lengths_[i]))*max(sds.find(level)->second,sd_[i]);
				double sqrtVar = sqrt(float(PI/2./fragmentNotNA_lengths_[i]))*localSV;
				double term = pow(d/sqrtVar,2);
				sum += term;
			}
		}
	}
	return sum;
}

void ChrCopyNumber::deleteFlanks(int telo_centromeric_flanks) {
	int maxRegionLengthToDelete = int(telo_centromeric_flanks/windowSize_);
	for (int i = 0; i < (int)medianValues_.size(); i++) {
		if (medianValues_[i]==NA) {
			if ((i-1>=0)&&(fragmentNotNA_lengths_[i-1]<=maxRegionLengthToDelete)) {
				deleteFragment(i-1);
				for (int j=i-2; j>0; j--) {
					if (fragmentNotNA_lengths_[j]<=maxRegionLengthToDelete)
						deleteFragment(j);
					else
						break;
				}
			}
			if ((i+1<(int)medianValues_.size())&&(fragmentNotNA_lengths_[i+1]<=maxRegionLengthToDelete))
				deleteFragment(i+1);
		}
	}
	//detele telomeric flanks:
	if (medianValues_[0]!=NA)
		if (fragmentNotNA_lengths_[0]<=maxRegionLengthToDelete)
			deleteFragment(0);

	if ((medianValues_[medianValues_.size()-1]!=NA)&&(fragmentNotNA_lengths_[medianValues_.size()-1]<=maxRegionLengthToDelete))
				deleteFragment(medianValues_.size()-1);

}


int ChrCopyNumber::removeLargeExons(float threshold) {
    int howManyRemoved = 0;
    for (int i =0; i< length_; i++) {
        if (ends_[i]-coordinates_[i]>threshold) {
            howManyRemoved++;
            readCount_[i]=NA;
        }
    }
    return howManyRemoved;
}

void ChrCopyNumber::recalcFlanks(int telo_centromeric_flanks, int minNumberOfWindows) {
    int maxRegionLengthToDelete = int(telo_centromeric_flanks/step_);
	for (int i = 0; i < (int)medianValues_.size(); i++) {
		if (medianValues_[i]==NA && fragment_lengths_[i]>=minNumberOfWindows) {
		    int left = i;
		    int right = i;
			if ((i-1>=0)&&(fragmentNotNA_lengths_[i-1]<=maxRegionLengthToDelete && medianValues_[i-1]!=NA )) {
				left = i-1;
				for (int j=i-2; j>0; j--) {
					if (fragmentNotNA_lengths_[j]<=maxRegionLengthToDelete && medianValues_[j]!=NA )
						left = j;
					else
						break;
				}
			}
			if ((i+1<(int)medianValues_.size())&&(fragmentNotNA_lengths_[i+1]<=maxRegionLengthToDelete)&& medianValues_[i+1]!=NA) {
				right = i+1;
				for (int j=i+2; j<(int)fragment_lengths_.size()-1; j++) {
					if (fragmentNotNA_lengths_[j]<=maxRegionLengthToDelete && medianValues_[j]!=NA )
						right = j;
					else
						break;
				}
            }
            if (left < i-1)
                recalcFlanksForIndeces (left,i-1);
            if (i+1 < right)
                recalcFlanksForIndeces (i+1, right);
		}
	}
	//recalculate telomeric flanks:
	if (medianValues_[0]!=NA)
		if (fragmentNotNA_lengths_[0]<=maxRegionLengthToDelete) {
            int i = 0;
            int right = 0;
            if ((i+1<(int)medianValues_.size())&&(fragmentNotNA_lengths_[i+1]<=maxRegionLengthToDelete)) {
				right = i+1;
				for (int j=i+2; j<(int)fragment_lengths_.size()-1; j++) {
					if (fragmentNotNA_lengths_[j]<=maxRegionLengthToDelete)
						right = j;
					else
						break;
				}
            }
            recalcFlanksForIndeces (i, right);
		}

	if ((medianValues_[medianValues_.size()-1]!=NA)&&(fragmentNotNA_lengths_[medianValues_.size()-1]<=maxRegionLengthToDelete)) {
        int i = medianValues_.size()-1;
        int left = medianValues_.size()-1;
        if ((i-1>=0)&&(fragmentNotNA_lengths_[i-1]<=maxRegionLengthToDelete)) {
            left = i-1;
			for (int j=i-2; j>0; j--) {
				if (fragmentNotNA_lengths_[j]<=maxRegionLengthToDelete)
					left = j;
                else
					break;
            }
            recalcFlanksForIndeces (left,i);
        }
	}
}

void ChrCopyNumber::recalcFlanksForIndeces (int i_start, int i_end) {
        vector <float> data;
        int notNA = 0;
        int totalCount = 0;
        int breakPointStart = 0;
        int breakPointEnd = 0;
        for (int i = i_start; i<=i_end; i++) { //collect data points:
            breakPointStart = 0;
            if (i>0)
                breakPointStart = bpfinal_[i-1]+1;
            breakPointEnd = bpfinal_[i];
            for (int j = breakPointStart; j <= breakPointEnd; j++)
                if (ratio_[j] != NA) {
                    data.push_back(ratio_[j]);
                    notNA++;
                }
            totalCount += breakPointEnd-breakPointStart+1;
        }
        float median;
        if (notNA==0 ||(notNA == 1 && totalCount>2)) {
           // median = get_median(data);
            median = NA;
        } else {
            median = get_median(data); //including the last point!!!!
        }
        float sd_local = sd(data,median);

        for (int i = i_start; i<=i_end; i++) {
            medianValues_[i] = median;
            sd_[i] = sd_local;
        }

        breakPointStart = 0;
        if (i_start>0)
                breakPointStart = bpfinal_[i_start-1]+1;
        for (int j = breakPointStart; j<= breakPointEnd; j++)
            medianProfile_[j] = median;
}

void ChrCopyNumber::clearCGcontent () {
	GCprofile_.clear ();
}
void ChrCopyNumber::clearNonNpercent () {
	notNprofile_.clear ();
}

void ChrCopyNumber::clearMappabilityProfile () {
	mappabilityProfile_.clear ();
}

void ChrCopyNumber::addToCGcontent (float valueToAdd) {
	GCprofile_.push_back(valueToAdd);
}
void ChrCopyNumber::addToNonNpercent (float valueToAdd) {
	notNprofile_.push_back(valueToAdd);
}

void ChrCopyNumber::addToMappabilityProfile (float valueToAdd) {
	mappabilityProfile_.push_back(valueToAdd);
}

void ChrCopyNumber::addToGenes_name(string i)
{
    genes_names.push_back(i);
}

void ChrCopyNumber::createMappabilityProfile() {
    for (int i=0; i<length_; i++)
        mappabilityProfile_.push_back(0);
}

void ChrCopyNumber::checkOrCreateNotNprofileWithZeros() {
    int sizeOfNonN = notNprofile_.size();
    if (sizeOfNonN<length_) {
        notNprofile_ = vector <float> (length_);
    }
    for (int i = 0; i< length_; i++)
        notNprofile_[i]=0;
}

void ChrCopyNumber::fillCGprofile(std::string const& chrFolder) {
	GCprofile_ = vector <float> (length_);
	notNprofile_ = vector <float> (length_);
	ifstream file;
	string filename = chromosome_;
	string possibleFilenames[] = {filename,filename+".fa",filename+".fasta","chr"+filename+".fa","chr"+filename+".fasta"};
	for (int i = 0; i < 5; i++) {

		string myFilename = possibleFilenames[i];
		string myFullPath = pathAppend(chrFolder,myFilename);
		file.open(myFullPath.c_str());

		if(!file.is_open())
			file.clear();
		else
			i = 6;
	}


    if (!file.is_open())	 {
		//	throw ("Unable to open fasta file for chr "+chromosome_+" in "+chrFolder+"\n");
        cerr << "Unable to open fasta file for chr "+chromosome_+" in folder "+chrFolder+"\n\nPlease note, "<< chrFolder << " should be a folder, not a file!\n\n";
        exit (-1);
	}
	//string myString;
	char letter; // here we read the chromosome letter by letter. Can do it better!
	file >> letter;
	int count = 0;
	int countCG = 0;
	int countN = 0;
	string line;
    string text = "";

	if (letter == '>')
		getline (file,line);
	else {
		count = 1;
		countCG = isCG(letter);
		countN = isN(letter);
		text.push_back(letter);
	}

	if (ends_.size()==0) { //all windows have equal length => can use the same windowsize for all windows
		for (int i = 0; i<length_; i++) {
			if (file.eof()) {
				GCprofile_[i] = NA;
				//cout << "End-of-file reached.." << endl;
			}
			while((!file.eof()) && (count < windowSize_)) {
				file>>letter;
				countCG += isCG(letter);
				countN += isN(letter);
				count ++;
			}
			notNprofile_[i] = float(count-countN)/count;
			if (count == countN)
				GCprofile_[i] = NA;
			else
				GCprofile_[i] = float(countCG)/(count-countN);
			//reset
			countCG = 0;
			countN = 0;
			count = 0;
		}
	} else {
		int start, end;
		for (int i = 0; i<length_; i++) {
			if (file.eof()) {
				GCprofile_[i] = NA;
				//cout << "End-of-file reached.." << endl;
			}
			start = coordinates_[i];
			end = ends_[i];
			while((!file.eof()) && (count < start)) {
				file>>letter;
				count ++;
			}
			while((!file.eof()) && (count <= end)) {
				file>>letter;
				text.push_back(letter);
				count ++;
			}
			notNprofile_[i] = 1;
			//notNprofile_[i] = float(end-start+1-countN)/(end-start+1);
			/*if (end-start+1 == countN)
				GCprofile_[i] = NA;
			else */

			countCG =0;
			countN = 0;
            for(int j = 0; j < (int)text.length(); j++) //++j????
                if (text[j] == 'C' || text[j] == 'G' || text[j] == 'c' || text[j] == 'g')
                    countCG++;
                else if (text[j] == 'N')
                    countN++;
            if (end-start+1-countN>0)
                GCprofile_[i] = float(countCG)/(end-start+1-countN);
            else
                GCprofile_[i] = NA;
			//reset
			if (i+1<length_) {
			    int nextStart = coordinates_[i+1];
			    if (nextStart<=end) {
                    //count = end;
                    //and delete prefix in text;
                    int howMuchToDelete = nextStart - start;
                    if (howMuchToDelete<0) {
                        cerr << "Error: your BED file with coordinates of targeted regions does not seem to be sorted\nCheck chromosome "<<chromosome_<<"\n";
                        cout << "Exit Control-FREEC: before reruning, please, sort the BED file with coordinates of the targeted regions\n";
                        exit(-1);
                    }
                    if (howMuchToDelete==0) {
                        cerr << "Error: your BED file with coordinates of targeted regions may contain duplicates\nCheck chromosome "<<chromosome_<<"\n";
                        cout << "Exit Control-FREEC: before reruning sort the BED file with coordinates of the targeted regions and REMOVE DUPLICATED REGIONS\n";
                        exit(-1);
                    }
                    text = text.substr(howMuchToDelete);
			    } else {
                    text = "";
			    }
			}

		}
	}
	file.close();

}

void ChrCopyNumber::deleteFragment(int i) {
	fragmentNotNA_lengths_[i] = 0;
	medianValues_[i] = NA;
	sd_[i] = NA;
}

float ChrCopyNumber::getBAFat (int i){
	return BAF_[i];
}

float ChrCopyNumber::getBAFProfileAt (int i){
	return estimatedBAFProfile_[i];
}

float ChrCopyNumber::getFittedBAFProfileAt (int i){
	return fittedBAFProfile_[i];
}

std::string ChrCopyNumber::getBAFsymbolAt (int i) {
	return medianBAFSymbol_[i];
}

float ChrCopyNumber::getEstimatedBAFuncertaintyAtI(int i) {
    if (i>0 &&i<int(estimatedBAFuncertainty_.size()))
        return estimatedBAFuncertainty_[i];
    return NA;
}



float ChrCopyNumber::getSmoothedProfileAtI(int i) {
	return smoothedProfile_[i];
}

float ChrCopyNumber::getSmoothedForInterval(int start , int end) {
   return get_median (smoothedProfile_,start,end);
}


void ChrCopyNumber::pushSmoothedProfile(float value) {
	smoothedProfile_.push_back(value);
}

int ChrCopyNumber::getEndsSize() {
	return ends_.size();
}

void ChrCopyNumber::setLookingForSubclones(bool value) {
    isLookingForSubclones_=value;
    if (value) {
        if (coordinates_.size()==0) {cerr << "Warning: you should intialize the ChrCopyNumber object before calling this function!!!\n";}
        if (copy_number_subc_.size()==0) {
            copy_number_subc_=vector <int> (coordinates_.size(),0);
        }
        if (population_subc_.size()==0) {
            population_subc_=vector <float> (coordinates_.size(),0.0);
        }
    }
}



ChrCopyNumber::~ChrCopyNumber(void)
{
	coordinates_.clear();
	readCount_.clear();
	smoothedProfile_.clear();
	fragmentNotNA_lengths_.clear(); //TODO all other vectors
	length_ = 0;
	copy_number_subc_.clear();
	population_subc_.clear();
}

void ChrCopyNumber::createBAF(float value) {
    for (int i=0; i<length_; i++)
        BAF_.push_back(value);
}
void ChrCopyNumber::createBAFvalues() {
    for (int i=0; i<length_; i++)
        BAFvalues_.push_back("");
}

void ChrCopyNumber::addBAFinfo(SNPinGenome & snpingenome,int indexSNP) {

    SNPatChr SNPsatChr = snpingenome.SNP_atChr(indexSNP);

    //create a vector with BAF
    createBAF(NA);
    createBAFvalues();
//BAFvalues_
    int totalSNPnumber = SNPsatChr.getSize() ;
    cout << "..Total Number of SNPs: "<< totalSNPnumber <<"\n";
    int SNPcount = 0;
    float currentBAF = SNPsatChr.getValueAt(SNPcount);
    float currentBAFstatus = SNPsatChr.getStatusAt(SNPcount);
    int getSNPpos = SNPsatChr.getPositionAt(SNPcount);
    float minBAF;


    for (int i = 0; i<length_; i++) {
        int left = coordinates_[i];
        int right = getEndAtBin(i);
        if (getSNPpos>=left && getSNPpos <=right) {
          //  snpingenome.SNP_atChr(indexSNP).setBinAt(SNPcount,i);
            snpingenome.setBinAt(indexSNP,SNPcount,i);

            if (currentBAFstatus !=0 && currentBAF != NA) { //there are values that indicate that this SNP can be heterozygios
                if (BAFvalues_[i] != "")
                    BAFvalues_[i] += ";";
                stringstream ss (stringstream::in | stringstream::out);
                ss << currentBAF;
                BAFvalues_[i] += ss.str();
            }
            minBAF = BAF_[i];
            if (minBAF==NA) {
                BAF_[i]=currentBAF;
            } else {
                if (fabs(minBAF-0.5) > fabs(currentBAF-0.5)){
                        BAF_[i]=currentBAF;
                //BAF_[i]=min(minBAF,currentBAF);

                }


            }
        } else if (getSNPpos<left) {
		     if (SNPcount+1 < totalSNPnumber) { // EV bug fix: 2013-01-17
                SNPcount++;
                getSNPpos = SNPsatChr.getPositionAt(SNPcount);
                currentBAF = SNPsatChr.getValueAt(SNPcount);
                currentBAFstatus = SNPsatChr.getStatusAt(SNPcount);
                if (windowSize_ == 0 && step_ == 0)
                    {
                    i=max(i-2,0); // assuming that the positions of SNPs are sorted
                    }
                else
                    {
                    i=max(-1,i-windowSize_/step_-2); //changed -1 to -2 in version 10.9
                    }
                    //cout << SNPcount << " out of "<< totalSNPnumber<<"\n";
            }
        }
    }
    if (ratio_.size()==0) {
        cerr << "Warning: Normalized read counts (ratio_) has not been initialized; check your parameters\n";
    }
    for (int i = 0; i<length_; i++) {
        if (ratio_.size()>i && ratio_[i]==NA && BAF_[i]!=NA)     //set BAF=NA in windows with ratio==NA to remove the noise from windows with low mappability
            BAF_[i]=NA;

        if (BAF_[i]==0)     //remove windows with 100% AA counts
            BAF_[i]=NA;

        if (BAF_[i]!=NA && BAF_[i]!=0 && BAF_[i]!=1) {
            //recalculate if using BAFvalues_[i]
            vector<string>heteroValuesPerWindowStrings = split(BAFvalues_[i], ';');
            if (heteroValuesPerWindowStrings.size()>0) {
                vector<float>heteroValuesPerWindow;
                for (int unsigned j = 0; j < heteroValuesPerWindowStrings.size(); j++) {
                    stringstream ss;
                    float f;
                    ss << heteroValuesPerWindowStrings[j];
                    ss >> f;
                    heteroValuesPerWindow.push_back(fabs(f-0.5));
                }
                float median = get_median(heteroValuesPerWindow)+0.5;
                BAF_[i] = median;
            } else {
//                if (BAF_[i]>0 && BAF_[i] <0.5)  //put the noise on top
//                    BAF_[i] = 1-BAF_[i];
                    BAF_[i] = NA; //delete all homoz.!
            }
        }
    }
}


//find closest value to k/ploidy for medianValue
float ChrCopyNumber::getLevelAt(int unsigned i, int ploidy) {
    float value = medianValues_[i];
	float valueToReturn;
	if (value >0)
		valueToReturn = float(int(value*ploidy + 0.5))/ploidy;
    else
        valueToReturn = float(int(value*ploidy))/ploidy;

    //use BAF for ambigious cases:
    if (valueToReturn>0 && estBAFuncertaintyPerFrag_.size()>i && BAFsymbPerFrag_[i].compare("-")!=0 &&BAFsymbPerFrag_[i]!="" && BAFsymbPerFrag_[i].length()!=valueToReturn*ploidy && (estBAFuncertaintyPerFrag_[i] < MAXUncertainty)) {
        //change Level value if there is no incertainty:
        valueToReturn=BAFsymbPerFrag_[i].length()*1./ploidy;
    }

    return valueToReturn;

}


void ChrCopyNumber::setRCountToZeroForNNNN() {
    for (unsigned int i = 0; i< notNprofile_.size(); i++) {
        if (notNprofile_[i]==0) {
            readCount_[i]=0;
        }
    }
    for (unsigned int i = notNprofile_.size(); i< readCount_.size(); i++) {
            readCount_[i]=0;
    }
}

void* ChrCopyNumber_calculateBreakpoint_wrapper(void *arg)
{
  ChrCopyNumberCalculateBreakpointArgWrapper* warg = (ChrCopyNumberCalculateBreakpointArgWrapper*)arg;
  int result = warg->chrCopyNumber.calculateBreakpoints(warg->breakPointThreshold, 0, warg->breakPointType);
  if (result <= 0) {
	cerr << "..failed to run segmentation on chr" << warg->chrCopyNumber.getChromosome() << "\n";
  }
  return NULL;
}

void* ChrCopyNumber_calculateBAFBreakpoint_wrapper(void *arg)
{
  ChrCopyNumberCalculateBreakpointArgWrapper* warg = (ChrCopyNumberCalculateBreakpointArgWrapper*)arg;
  int result = warg->chrCopyNumber.calculateBAFBreakpoints(warg->breakPointThreshold, 0, warg->breakPointType);
  if (result == 0) {
	cerr << "..failed to run BAF segmentation on chr" << warg->chrCopyNumber.getChromosome() << "\n";
  }
  return NULL;
}

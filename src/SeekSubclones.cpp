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


#include "SeekSubclones.h"

using namespace std;

SeekSubclones::SeekSubclones(void)
{
}

SeekSubclones::~SeekSubclones(void)
{
}


SeekSubclones::SeekSubclones(GenomeCopyNumber & samplecopynumber, int ploidy, std::string myName, float minimal_pop) {
    ploidy_ = ploidy;
    minimal_pop_  = minimal_pop/100;
    getSegmentsInfo(samplecopynumber, myName);
}


void SeekSubclones::getSegmentsInfo(GenomeCopyNumber & samplecopynumber, std::string myName)
{
    string::size_type pos = 0;
	map<string,int>::iterator it;
    ofstream myfile;
    std::string Newfile = myName + "_subclones" +  ".txt";
    double bonfer_correction = 0;
    int thresholdOnChrLengthForSubcloneDetection=20;
    int numberOfChromosomes = samplecopynumber.getNumberOfChromosomes();
    for (int index=0; index< numberOfChromosomes; index++) {
        int length = samplecopynumber.getChrCopyNumberAt(index).getLength();
        if (length<thresholdOnChrLengthForSubcloneDetection) {
            continue;
        }
        string chrNumber = samplecopynumber.getChrCopyNumberAt(index).getChromosome();
        if ( ( pos = chrNumber.find("X")) != string::npos ) 		//exclude X and Y from the analysis
            continue;
        if ( ( pos = chrNumber.find("Y")) != string::npos )
            continue;
        bonfer_correction += samplecopynumber.getChrCopyNumberAt(index).getBreakPoints().size();
    }
    if (bonfer_correction == 0)  {  bonfer_correction = 1;  }
    myfile.open(Newfile.c_str());

    for (int index=0; index< numberOfChromosomes; index++) {
        int length = samplecopynumber.getChrCopyNumberAt(index).getLength();

        if (length<thresholdOnChrLengthForSubcloneDetection) {
            continue;
        }
        string chrNumber = samplecopynumber.getChrCopyNumberAt(index).getChromosome();
        processChrName(chrNumber);

        if ( ( pos = chrNumber.find("X")) != string::npos ) 		//exclude X and Y from the analysis
            continue;
        if ( ( pos = chrNumber.find("Y")) != string::npos )
            continue;

        vector <int> breakpoints = samplecopynumber.getChrCopyNumberAt(index).getBreakPoints() ;
        int bpstart=0;
        int bpend=0;
        if (breakpoints.back()<length) {breakpoints.push_back(length);}
        if (breakpoints.size()!=samplecopynumber.getChrCopyNumberAt(index).getMedianValues().size()) {
            cerr << "Problem with the number of breakpoints:" <<breakpoints.size()<< " != " << samplecopynumber.getChrCopyNumberAt(index).getMedianValues().size()<<"\n";
            exit(-1);
        }
        for (unsigned int fragmentCount=0; fragmentCount<breakpoints.size(); fragmentCount++) {
            bpend=breakpoints[fragmentCount];
            float fragmentMedian=samplecopynumber.getChrCopyNumberAt(index).getMedianValuesAt(fragmentCount);

            float expected =round_by_ploidy(fragmentMedian,ploidy_);
            if (samplecopynumber.getChrCopyNumberAt(index).isSmoothed()) {
                expected=samplecopynumber.getChrCopyNumberAt(index).getSmoothedForInterval(bpstart,bpend);
            }

            vector<float>::const_iterator first = samplecopynumber.getChrCopyNumberAt(index).getRatio().begin() + bpstart;
            vector<float>::const_iterator last = samplecopynumber.getChrCopyNumberAt(index).getRatio().begin() + bpend;
            vector<float> data(first, last);

            if (data.size()> 0 && expected >= 0) {
                bool subclonedetected = SignTest(data, expected, bonfer_correction);
                //bool subclonedetected = PercentageTest(data, expected);
                if (subclonedetected == true) {
                    EstimateSubclonalPopulation(data, expected, ploidy_); //fills in population_ and copynumber_
                    if (copynumber_.size() > 0) {
                        myfile << "Possible subclones for fragment chr" << chrNumber << ":" << samplecopynumber.getChrCopyNumberAt(index).getCoordinateAtBin(bpstart) << "-" << samplecopynumber.getChrCopyNumberAt(index).getEndAtBin(bpend-1) << "\n";
                        myfile << "Major clone is suggest to have " << expected*ploidy_ << " copies\n";
                        myfile << "\t Copy number in Subclone (different possibilities) \t Subclonal population \n";
                        for (unsigned int k = 0; k < copynumber_.size(); k++) {
                            myfile << "\t" << copynumber_[k] << "\t" << population_[k]*100 << "% \n";
                        }
                        myfile << "\n";
                        for (int k = bpstart; k < bpend; k++) { // CARINO, WHY WAS IT: k = bpstart-2; k < i-1; ???
                            samplecopynumber.getChrCopyNumberAt(index).setCN_subc(k, copynumber_[0]);
                            samplecopynumber.getChrCopyNumberAt(index).setPopulation_subc(k, population_[0]);
                        }
                    }
                    copynumber_.clear();
                    population_.clear();
                }
            }
            data.clear();
            bpstart=bpend;
        }
        breakpoints.clear();
	}
    myfile.close();
}


bool SeekSubclones::SignTest(std::vector <float>& data, float& threshold, int bonfer_correction)
{
    int upvalues = 0;
    int downvalues = 0;
    bool subclone = false;
    for (unsigned int i = 0; i < data.size(); i++) {
        if (data[i]!=NA){
            if (data[i] < threshold)  {
                downvalues++;
            }
            if (data[i] > threshold)  {
                upvalues++;
            }
        }
    }
    int max_count = upvalues;
    if (upvalues < downvalues)
        {
        max_count = downvalues;
        }
    double result = 2*binomialcdistribution(max_count-1, upvalues+downvalues, 0.5); //CARINO! YOU CALULATE A WRONG P-VALUE HERE!! SHOULD BE max_count-1; I CORRECTED IT.
    if (result < 0.001/bonfer_correction && result != 0)
        {
        subclone = true;
        }
    return subclone;
}

bool SeekSubclones::PercentageTest(std::vector <float>& data, float& threshold)
{
    int upvalues = 0;
    int downvalues = 0;
    bool subclone = false;
    for (unsigned int i = 0; i < data.size(); i++)
        {
        if (data[i] < threshold)
            {
            downvalues++;
            }
        if (data[i] > threshold)
            {
            upvalues++;
            }
        }
    int max_count = upvalues;
    int min_count = downvalues;
    if (upvalues < downvalues)
        {
        max_count = downvalues;
        min_count = upvalues;
        }
    double percentage = 1;
    if (max_count > 100)
        {
        percentage = (float)min_count/(float)max_count;
        }
    else {percentage = 1;}
    if (percentage < 0.20 && (min_count + max_count > 100))
        {
        subclone = true;
        }
    return subclone;
}


void SeekSubclones::EstimateSubclonalPopulation(vector <float> data, float threshold, int ploidy_)
{
    float sumtmp = 0;
    threshold = threshold*ploidy_;
    int countData = 0;
    for (unsigned int i = 0; i < data.size(); i++) {
        if (data[i]!=NA) {
            sumtmp += data[i];
            countData++;
        }

    }
    if (countData==0) {return;}
    float mean = ploidy_*sumtmp/countData;
    if (mean > threshold)
        {
        int i =1;
        float pop = 1;
        int iter = 0;
        while (pop > minimal_pop_ && iter < 100)
            {
            if ((((mean - threshold)/(i)) > minimal_pop_) && (((mean - threshold)/(i)) < 1) ) // && (threshold + i != ploidy_)
                {
                copynumber_.push_back(threshold + i);
                population_.push_back((mean - threshold)/(i));
                }
            pop = (mean - threshold)/(i);
            i++;
            iter++;
            }
        }
    else if (mean < threshold)
        {
        int i =1;
        float pop = 1;
        int iter = 0;
        while (pop > minimal_pop_ && iter < 100 && (threshold - i >= 0))
            {
            if ((((-mean + threshold)/(i)) > minimal_pop_) && (((-mean + threshold)/(i)) < 1) ) //GOT IT EXCEPT for "&& (threshold - i != ploidy_)" SO I REMOVED IT
                {
                copynumber_.push_back(threshold - i);
                population_.push_back((-mean + threshold)/(i));
                }
            pop = (-mean + threshold)/(i);
            i++;
            iter++;
            }
        }
}

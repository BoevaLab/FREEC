#pragma once
#ifndef SEEKSUBCLONES_H
#define SEEKSUBCLONES_H

#include<cstdlib>
#include<ctime>
#include "GenomeCopyNumber.h"
#include "binomialdistr.h"
#include "myFunc.h"
#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <map>
#include <fstream>

class SeekSubclones
{
    public:
        SeekSubclones(GenomeCopyNumber & samplecopynumber, int ploidy, std::string myName, float minimal_pop);
        SeekSubclones();
        ~SeekSubclones(void);
        void getSegmentsInfo(GenomeCopyNumber & samplecopynumber, std::string myName);
        bool SignTest(std::vector <float>& data, float& threshold, int bonfer_correction);
        void EstimateSubclonalPopulation(std::vector <float> data, float threshold, int ploidy_);
        bool PercentageTest(std::vector <float>& data, float& threshold);
    protected:
        int ploidy_;
        float minimal_pop_;
        std::vector <int> copynumber_;
        std::vector <float> population_;
    private:
};

#endif // SEEKSUBCLONES_H

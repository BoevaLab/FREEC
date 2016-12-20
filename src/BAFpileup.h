#pragma once
#ifndef BAFPILEUP_H
#define BAFPILEUP_H

#include <string>
#include <vector>
#include "GenomeCopyNumber.h"

class BAFpileup
{
    public:
        BAFpileup();
        void makepileup(GenomeCopyNumber & sampleCopyNumber, GenomeCopyNumber & controlCopyNumber,
            std::string sample_MateFile,std::string control_MateFile,
            std::string outputDir,std::string makepileup, std::string const& mateFileName,
            std::string const& inputFormat, std::string const& matesOrientation,
            std::string pathToSamtools, std::string chrLen, std::string controlName,
            std::string pathToSambamba, std::string SambambaThreads, std::string targetBed = "",
            std::string pathToBedtools = "", std::string fastaFile="", int minQualPerPos=0);
        float calculateFlankLength(std::string const& mateFileName ,std::string const& inputFormat, std::string const& matesOrientation, std::string pathToSamtools,std::string pathToSambamba, std::string SambambaThreads);
        void calculateNewBoundaries(std::string targetBed, int flanks, std::string bedFileWithRegionsOfInterest);
        std::string intersectWithBedtools(std::string makeminipileup, std::string outputDir,  std::string bedFileWithRegionsOfInterest, std::string chrLen);
        void createBedFileWithChromosomeLengths (std::string bedFileWithRegionsOfInterest, std::string chrLenFile, bool doesNeedChrPrefix);
        std::string createPileUpFile(std::string  outputDir, std::string samtools_path, std::string pathToSambamba, std::string SambambaThreads, std::string control_tumor, std::string intersected, std::string fastaFile,int minQualPerPos);
        std::vector < std::vector<float> >computeBAF(GenomeCopyNumber & sampleorcontrol, std::string minipileup, std::string outputDir, std::string filename);
        std::vector <int> coordinates_;
        std::vector <int> ends_;
        std::vector <int> bpfinalBAF_;
        std::vector < std::vector<std::string> >snp_pos_pileup;
        int exons_Count;
        int exons_Counttmp;
        std::vector <int> findBreakpoints(GenomeCopyNumber & sampleCopyNumber, double threshold, int breakPointType);
        std::vector < std::vector<float> >BAFtumor;
        std::vector < std::vector<std::string> > chr;
        std::vector < std::vector<int> > A_nb;
        std::vector < std::vector<int> > B_nb;
        std::vector < std::vector<std::string> > AB_nb;

    private:

        int length_;
        std::vector<std::string> snp_pos;
        std::vector<std::string> ref_base;
        std::vector<std::string> alt_base;
        std::string chromosome_;
        std::vector <std::string> coordinatesTmp_;
        std::vector <std::string> endsTmp_;
        std::vector <std::string> chr_names;
        std::vector <int> length_with_flanks;
        std::vector <std::string> length_with_flanksTmp;
        std::vector <std::string> strand;
        std::vector <std::string> ref_name;
        std::string pathToBedtools_;
};

#endif // BAFPILEUP_H

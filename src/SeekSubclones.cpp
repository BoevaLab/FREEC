#include "SeekSubclones.h"

using namespace std;

SeekSubclones::SeekSubclones(void)
{
}

SeekSubclones::~SeekSubclones(void)
{
}


SeekSubclones::SeekSubclones(GenomeCopyNumber & samplecopynumber, int ploidy, std::string outputDir, float minimal_pop)
{
ploidy_ = ploidy;
minimal_pop_  = minimal_pop;
getSegmentsInfo(samplecopynumber, outputDir);
}


void SeekSubclones::getSegmentsInfo(GenomeCopyNumber & samplecopynumber, std::string outputDir)
{
    string::size_type pos = 0;
	map<string,int>::iterator it;
    ofstream myfile;
    std::string Newfile = outputDir + "Subclones" +  ".txt";
    double bonfer_correction = 0;
    for ( it=samplecopynumber.chromosomesInd_.begin() ; it != samplecopynumber.chromosomesInd_.end(); it++ )
        {
        string chrNumber = (*it).first;
		if ( ( pos = chrNumber.find("chr", pos)) != string::npos )
			chrNumber.replace( pos, 3, "" );
        if ( ( pos = chrNumber.find("X", pos)) != string::npos ) 		//exclude X and Y from the analysis
            continue;
        if ( ( pos = chrNumber.find("Y", pos)) != string::npos )
            continue;
		int index = samplecopynumber.findIndex(chrNumber);
		if (index > 22)
            {
            continue;
            }
        bonfer_correction += samplecopynumber.chrCopyNumber_[index].getBreakPoints().size();
        }
    if (bonfer_correction == 0)
        {
        bonfer_correction = 1;
        }
    myfile.open(Newfile.c_str());
	for ( it=samplecopynumber.chromosomesInd_.begin() ; it != samplecopynumber.chromosomesInd_.end(); it++ ) {
		string chrNumber = (*it).first;
		if ( ( pos = chrNumber.find("chr", pos)) != string::npos )
			chrNumber.replace( pos, 3, "" );
        if ( ( pos = chrNumber.find("X", pos)) != string::npos ) 		//exclude X and Y from the analysis
            continue;
        if ( ( pos = chrNumber.find("Y", pos)) != string::npos )
            continue;
		int index = samplecopynumber.findIndex(chrNumber);
		if (index > 22)
            {
            continue;
            }
		int length = samplecopynumber.chrCopyNumber_[index].getLength();
		int i = 0;
		while(i < length) {
            vector <float> data;
			float expected = round_by_ploidy(samplecopynumber.chrCopyNumber_[index].getMedianProfileAtI(i),ploidy_);
            int bpstart = i;
            if (expected < 0)
                {
                while (round_by_ploidy(samplecopynumber.chrCopyNumber_[index].getMedianProfileAtI(i),ploidy_) < 0 && i < length)
                    {
                    i++;
                    }
                }
            float threshold = expected;
            while (expected == threshold && i < length)
                {
                float observed = samplecopynumber.chrCopyNumber_[index].getRatioAtBin(i);
                data.push_back(observed);
                expected = round_by_ploidy(samplecopynumber.chrCopyNumber_[index].getMedianProfileAtI(i),ploidy_);
                i++;
                }
                data.pop_back();
                if (data.size()> 0 && threshold >= 0)
                {
                bool subclonedetected = SignTest(data, threshold, bonfer_correction);
                //bool subclonedetected = PercentageTest(data, threshold);
                if (subclonedetected == true)
                {
                    EstimateSubclonalPopulation(data, threshold, ploidy_);
                    if (copynumber.size() > 0)
                        {myfile << "Possible subclones for fragment chr" << chrNumber << ":" << samplecopynumber.chrCopyNumber_[index].getCoordinateAtBin(bpstart) << "-" << samplecopynumber.chrCopyNumber_[index].getEndAtBin(i) << "\n";
                        myfile << "Considering only one clonal population, it would have a copy number of : " << threshold*ploidy_ << "\n";
                        myfile << "\t Copynumber in Subclone \t Subclonal population \n";
                    for (int k = 0; k < copynumber.size(); k++)
                        {
                        myfile << "\t" << copynumber[k] << "\t" << population[k]*100 << "% \n";
                        }
                    for (int k = bpstart-2; k < i-1; k++)
                        {
                        if (copynumber.size() > 0 && k >= 0 && index != NA)
                            {
                            samplecopynumber.chrCopyNumber_[index].setCN_subc(k, copynumber[0]);
                            samplecopynumber.chrCopyNumber_[index].setPopulation_subc(k, population[0]);
                            }
                        else if (copynumber.size() == 0 && k >= 0)
                            {
                            samplecopynumber.chrCopyNumber_[index].setCN_subc(k, -1);
                            samplecopynumber.chrCopyNumber_[index].setPopulation_subc(k, -1);
                            }
                        }
                    }
                copynumber.clear();
                population.clear();
                }
        data.clear();
		i++;
                }
		}
	}
    myfile.close();
}


bool SeekSubclones::SignTest(std::vector <float>& data, float& threshold, int bonfer_correction)
{
    int upvalues = 0;
    int downvalues = 0;
    bool subclone = false;
    for (int i = 0; i < data.size(); i++)
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
    if (upvalues < downvalues)
        {
        max_count = downvalues;
        }
    double result = 2*binomialcdistribution(max_count, upvalues+downvalues, 0.5);
    if (result < 0.01/bonfer_correction && result != 0)
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
    for (int i = 0; i < data.size(); i++)
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
    float meantmp = 0;
    threshold = threshold*ploidy_;
    for (int i = 0; i < data.size(); i++)
        {
        data[i] = ploidy_*data[i];
        meantmp += data[i];
        }
    float mean = meantmp/data.size();
    if (mean > threshold)
        {
        int i =1;
        float pop = 1;
        int iter = 0;
        while (pop > minimal_pop_ && iter < 100)
            {
            if ((((mean - threshold)/(i)) > minimal_pop_) && (((mean - threshold)/(i)) < 1) && (threshold + i != ploidy_))
                {
                copynumber.push_back(threshold + i);
                population.push_back((mean - threshold)/(i));
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
        while (pop > minimal_pop_ && iter < 100)
            {
            if ((((-mean + threshold)/(i)) > minimal_pop_) && (threshold - i > 0) && (((-mean + threshold)/(i)) < 1) && (threshold - i != ploidy_))
                {
                copynumber.push_back(threshold - i);
                population.push_back((-mean + threshold)/(i));
                }
            pop = (-mean + threshold)/(i);
            i++;
            iter++;
            }
        }
}

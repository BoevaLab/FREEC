#include "RSSerror.h"

using namespace std;

RSSerror::RSSerror()
{
}

long double calculateRSS(GenomeCopyNumber & samplecopynumber, int ploidy)
{
    string::size_type pos = 0;
    vector<float> observedvalues;
    vector<float> expectedvalues;
	map<string,int>::iterator it;
	for ( it=samplecopynumber.chromosomesInd_.begin() ; it != samplecopynumber.chromosomesInd_.end(); it++ ) {
		string chrNumber = (*it).first;
		if ( ( pos = chrNumber.find("chr", pos)) != string::npos )
			chrNumber.replace( pos, 3, "" );
        if ( ( pos = chrNumber.find("X", pos)) != string::npos ) 		//exclude X and Y from the analysis
            continue;
        if ( ( pos = chrNumber.find("Y", pos)) != string::npos )
            continue;
		int index = samplecopynumber.findIndex(chrNumber);
		int length = samplecopynumber.chrCopyNumber_[index].getLength();
		for (int i = 0; i< length; i++) {
            float observed = 0;
            float expected = 0;
                {
                observed = samplecopynumber.chrCopyNumber_[index].getRatioAtBin(i);
                expected = observed;
                if (samplecopynumber.chrCopyNumber_[index].isMedianCalculated()) {
                    expected = samplecopynumber.chrCopyNumber_[index].getMedianProfileAtI(i);
                    if (samplecopynumber.chrCopyNumber_[index].isSmoothed())
                        expected = samplecopynumber.chrCopyNumber_[index].getSmoothedProfileAtI(i);
                    }
                }
			observedvalues.push_back(observed);
			expectedvalues.push_back(expected);
		}
	}

    long double RSS = 0;
    for (int i = 0; i < (int)observedvalues.size(); i++)
        {
        if ((observedvalues[i]!=NA) && (expectedvalues[i]!=NA))
            {
            long double diff = (long double)observedvalues[i] - (long double)round(ploidy*expectedvalues[i])/ploidy;
            RSS = RSS + (long double)pow(diff,2);
            }
        }
    return RSS;
}

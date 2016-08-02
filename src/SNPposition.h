#ifndef SNPPOSITION_H
#define SNPPOSITION_H

#include <string>
#include <vector>

#include "myFunc.h"

class SNPposition
{
    public:
	    SNPposition(int position, char* letters, const char* strand, const char* ref); //from TXT
        SNPposition(int position, char* alt); //for VCF
        virtual ~SNPposition();
        int getPosition();
        char getNucleotide();
        float getValue();
        float getStatus();
        void setFrequency(float freq);
        void setStatus(float status);
        void setBin(int i);
        int getBin();
    protected:
    private:
        unsigned int position_;
        char nucleotide_;
        float freq_;
        float status_; //0 - AA/BB, 0.5 - AB
        int bin_; //bin in the copy number array
};

#endif // SNPPOSITION_H

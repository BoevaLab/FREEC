#ifndef SNPATCHR_H
#define SNPATCHR_H

#include <string>
#include <vector>

#include "SNPposition.h"


class SNPatChr
{
    public:
        SNPatChr(const std::string&);
        SNPatChr();
        void setChromosome(const std::string& chromosome);
        void push_SNP(const SNPposition& snp);
        const std::string& getChromosome ();
        int getPositionAt (int index);
        char getNucleotideAt(int index);
        float getValueAt(int index);
        float getStatusAt(int index);
        void setValueAt(int positionCount,float value);
        void setBinAt(int index,int bin);
        int getBinAt(int index); // return NA(-1) if there is not correspondance

        void setStatusAt(int index,float value) ;
        int getSize();
        virtual ~SNPatChr();
    protected:
    private:
    std::vector <SNPposition> SNPpositionArray_;
    std::string chromosome_;
};

#endif // SNPATCHR_H

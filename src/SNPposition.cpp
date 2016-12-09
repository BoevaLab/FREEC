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

#include "SNPposition.h"
#include <assert.h>

SNPposition::SNPposition(int position, char* alt) //for a VCF line
{
    position_=position;
    if (strlen(alt) == 1)   {
        nucleotide_ = alt[0];
    }  else  {
        char* strs[4];
        unsigned strs_cnt = split(alt, ',', strs);
        nucleotide_ = strs[0][0];
    }
    freq_ = 0; // EV: must be initialized
    status_ = 0; // EV: must be initialized
    bin_=NA; //before initialization
}


SNPposition::SNPposition(int position, char* letters, const char* strand, const char* ref) //for .txt line
{
    position_=position;
    char* strs[4];
    bool reverse = strcmp(strand, "-") == 0;
    if (strlen(letters) == 1)    { //should not get here
        nucleotide_ = letters[0];
    }   else    {
        unsigned strs_cnt = split(letters, '/', strs);
        char c_ref;
        if (reverse) {
            c_ref = complement(ref[0]);
        } else {
            c_ref = ref[0];
        }
        if (strs[0][0]==c_ref) {
            nucleotide_ = strs[1][0];
        } else {
            nucleotide_=strs[0][0];
        }
    }
    if (reverse) {
        nucleotide_ = complement(nucleotide_);
    }
    freq_ = 0; // EV: must be initialized
    status_ = 0; // EV: must be initialized
    bin_=NA; //before initialization
}


SNPposition::~SNPposition()
{
    //dtor
}

int SNPposition::getPosition() {
    return position_;
}

char SNPposition::getNucleotide() {
    return nucleotide_;
}

void SNPposition::setFrequency(float freq) {
        freq_ = freq;
}

void SNPposition::setStatus(float status) {
        status_ = status;
}

float SNPposition::getValue() {
    return freq_ ;
}

float SNPposition::getStatus() {
    return status_ ;
}

void SNPposition::setBin(int i) {
    bin_=i;
}
int SNPposition::getBin() {
    return bin_;
}

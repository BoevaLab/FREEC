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

#include "SNPatChr.h"

SNPatChr::SNPatChr(const std::string& chromosome)
{
    chromosome_ = chromosome;
}


SNPatChr::SNPatChr()
{
    //chromosome_ = NULL;
}

SNPatChr::~SNPatChr()
{
    //dtor
}

void SNPatChr::push_SNP(const SNPposition& snp) {
   SNPpositionArray_.push_back(snp);
}

void SNPatChr::setChromosome(const std::string& chromosome) {
    chromosome_ = chromosome;
}

const std::string& SNPatChr::getChromosome () {
    return chromosome_;
}

int SNPatChr::getPositionAt (int index) {
    return SNPpositionArray_[index].getPosition();
}

char SNPatChr::getNucleotideAt(int index) {
    return SNPpositionArray_[index].getNucleotide();
}

void SNPatChr::setValueAt(int index,float value) {
    SNPpositionArray_[index].setFrequency(value);
}

void SNPatChr::setStatusAt(int index,float value) {
    SNPpositionArray_[index].setStatus(value);
}

float SNPatChr::getValueAt(int index) {
    return SNPpositionArray_[index].getValue();
}

float SNPatChr::getStatusAt(int index) {
    return SNPpositionArray_[index].getStatus();
}

int SNPatChr::getSize() {
    return SNPpositionArray_.size();
}

void SNPatChr::setBinAt(int index,int bin) {
    SNPpositionArray_[index].setBin(bin);
}
int SNPatChr::getBinAt(int index) {
    return SNPpositionArray_[index].getBin();
}

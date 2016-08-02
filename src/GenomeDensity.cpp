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


#include "GenomeDensity.h"
using namespace std ;

GenomeDensity::GenomeDensity() {
}
GenomeDensity::GenomeDensity(std::string const& mateFileName ,std::string const& inputFormat, std::string const& matesOrientation, std::string const& chrLenFileName )
{
	//read fileWithChrLengths and iniciate genome density array gd_

	std::vector<std::string> names;
	std::vector<int> lengths;
	readFileWithGenomeInfo(chrLenFileName, names, lengths);
	/*for (int i = 0; i < (int)names.size(); i++) {
		ChrDensity chrDensity (lengths[i]+1);
		gd_.insert(pair<string,ChrDensity> (names[i],chrDensity));
	}	*/

	//read mateFileName and calculate densities, mu and sigma
	int count = 0;
	cout << "..Starting reading "<< mateFileName << "\n";
	std::ifstream fileMates (mateFileName.c_str());
	vector<float> insertSizeVector;
	string line;
	if ((inputFormat.compare("eland")==0 || inputFormat.compare("Eland")==0)&&(matesOrientation.compare("0")!=0)) {
		
		while (std::getline(fileMates,line)) {

			if (! line.length()) continue;
			//if (line[0] == '#') continue;	
		
			std::vector<std::string> strs = split(line, '\t');
			if (strs.size()==14) {
				
				if (strs[11].compare(strs[12])!=0) continue;
				//check for normal orientation
				if ((matesOrientation.compare(strs[8]+strs[10])==0 && atoi(strs[13].c_str())>0) || (matesOrientation.compare(strs[10]+strs[8])==0 && atoi(strs[13].c_str())<0)) {
					//cout << line <<"\n";
					string chr = strs[11];
					//if (gd_.find(chr)==gd_.end()) continue;
					int left = atoi(strs[7].c_str());
					int right = atoi(strs[9].c_str());
					int length = strs[1].length();					
					//gd_.find(chr)->second.coveragePlusOneAtI(left+length,right);
					float insertSize = (float)abs(right - left); //- length; 
					insertSizeVector.push_back (insertSize);
					//sumInsertLengths += insertSize;
					//sumInsertLengthsX2 += insertSize*insertSize;
					count++;
				}
			}
			strs.clear();
		}
	} else if ((inputFormat.compare("bowtie")==0 || inputFormat.compare("Bowtie")==0)&&(matesOrientation.compare("0")==0)){
		while (std::getline(fileMates,line)) {

			if (! line.length()) continue;
			//if (line[0] == '#') continue;	
		
			std::vector<std::string> strs = split(line, '\t');
			if (strs.size() > 4) {
				string chr = strs[2];
				string::size_type pos = 0;
				if ( ( pos = chr.find("chr", pos)) != string::npos ) 
					chr.replace( pos, 3, "" );	
				else {
					if (chr.find("NC_")== string::npos) {
						int $len = chr.length();
						chr = chr.substr($len-2); 
						if (chr.compare("23"))
							chr = "X";
						else if (chr.compare("24"))
							chr = "Y";
						else if (chr.compare("25"))
							chr = "M";
						else if (chr[0]=='0')
							chr = chr[1];						
					}
				} 
				if (gd_.find(chr)==gd_.end()) continue;
				int left = atoi(strs[3].c_str());
				
				int length = strs[1].length();
				
				for (int i = left; i<left+length; i++)
					gd_.find(chr)->second.coveragePlusOneAtI(i);				
				//sumInsertLengths += insertSize;
				//sumInsertLengthsX2 += insertSize*insertSize;
				
			}
			strs.clear();
		}

	}

	fileMates.close();

	//float sumInsertLengths = 0;
	//float sumInsertLengthsX2 = 0;


	if (count>0) {  //calculation of median and variation for insert size for mate pair data
		int numberOfPointsToEvaluateIQR = 100000;
		vector <float> selectedValues(numberOfPointsToEvaluateIQR) ;	
		int step = count/numberOfPointsToEvaluateIQR;

		for (int i = 0; i<numberOfPointsToEvaluateIQR; i++) 
			selectedValues[i] = insertSizeVector[i*step];
		
		sort(selectedValues.begin(),selectedValues.end());

		float median = get_median(selectedValues);
		float iqr = get_iqr(selectedValues);
		mu_ = median;
		sigma_ = iqr/1.349;
		cout << "\tmu "<< mu_ << "\n";
		cout << "\tsigma "<< sigma_ << "\n";
		cout << "\tcount "<< count << "\n";
		selectedValues.clear();
	}
	totalNumberOfPairs_ = count;
	insertSizeVector.clear();
}

int GenomeDensity::getTotalNumberOfPairs() {
	return totalNumberOfPairs_;
}

void GenomeDensity::smooth() {
	//double bindwidth = sigma_;
	int r = int (sigma_);
	cout << "\t..creating a kernel vector.. with radius"<< r <<" \n";
	KernelVector kernelVector(r); 
	map<string,ChrDensity>::iterator it;	
	for ( it=gd_.begin() ; it != gd_.end(); it++ ) {
		cout << "\t..smoothing for "<<(*it).first << "\n";
		(*it).second.smooth(kernelVector);	
	}
}

void GenomeDensity::writeToFile(std::string outFileName,std::string type, bool ifBinary) {
	std::ofstream file;
	if (ifBinary)
		file.open(outFileName.c_str(), ios::in | ios::binary);
	else 
		file.open(outFileName.c_str());
	
	class IntPoint {
	public:
		string chr;
		int pos;
		int value;
	};
	class FloatPoint {
	public:
		string chr;
		int pos;
		float value;
	};
	const string sep = "\t";
	map<string,ChrDensity>::iterator it;	
	if (type.compare("norm")==0) {
		IntPoint x;	
		for ( it=gd_.begin() ; it != gd_.end(); it++ ) {
			x.chr = (*it).first;
			for (int i = 0; i<(*it).second.getLength();i++) {
				x.pos = i;				
				x.value = (*it).second.getCoverageAtI(i);	
				if (ifBinary)
					file.write((char*)&x, sizeof (IntPoint));
				else {					
					file << x.chr+sep;
					file << x.pos;					
					file << sep;
					file << x.value;
					file << "\n";
				}
			}			
		}
	}
	else {
		FloatPoint x;
		for ( it=gd_.begin() ; it != gd_.end(); it++ ) {
			x.chr = (*it).first;
			for (int i = 0; i<(*it).second.getLength();i++) {
				x.pos = i;				
				x.value = (*it).second.getDensityAtI(i);	
				if (ifBinary)
					file.write((char*)&x, sizeof (FloatPoint));
				else {									
					file << x.chr+sep;
					file << x.pos;					
					file << sep;
					file << x.value;
					file << "\n";
				}
			}			
		}
	}
	file.close();
}

void GenomeDensity::printDensity(std::string key,int i,int j) {
	if (gd_.find(key)==gd_.end())
		return ;
	
	for (int k = i;k<j;k++) 
		cout << gd_.find(key)->second.getDensityAtI(k) << "\n";

}

void GenomeDensity::calculateLogRatio(GenomeDensity controlGD, std::string type) {
	map<string,ChrDensity>::iterator it;	
	for ( it=gd_.begin() ; it != gd_.end(); it++ ) {		
		cout << "..calculating density ratio for " << (*it).first << "\n";
		(*it).second.calculateLogRatio(controlGD.gd_.find((*it).first)->second,type);
	}
}

GenomeDensity::~GenomeDensity(void)
{
	gd_.clear();
}

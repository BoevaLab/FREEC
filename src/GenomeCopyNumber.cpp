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


#include "GenomeCopyNumber.h"

using namespace std ;

GenomeCopyNumber::GenomeCopyNumber(void)
{
	ploidy_ = NA;
	step_=NA;
	hasBAF_=0;
	ifUsedControl_ = false;
	normalContamination_=0;
	sex_="";
	SeekingSubc_ = false;
	isMappUsed_=false;
    totalNumberOfPairs_=0;
	normalNumberOfPairs_=0;
}

bool GenomeCopyNumber::isMappUsed() {return isMappUsed_;}

void GenomeCopyNumber::readCopyNumber(std::string const& mateFileName ,std::string const& inputFormat, std::string const& matesOrientation, std::string const& chrLenFileName, float coefficientOfVariation ) {
	//first get the number of reads and the genome size
	//long genomeSize = 0;
	long readNumber = 0;
	int windowSize;
	//reading the file with genome information
	std::vector<std::string> names;
	std::vector<int> lengths;
	readFileWithGenomeInfo(chrLenFileName, names, lengths);
	refGenomeSize_ = sum(lengths);
	cout << "\t total genome size:\t" << refGenomeSize_ << "\n";
	if ((inputFormat.compare("pileup")==0 || inputFormat.compare("SAMtools pileup")==0)) {
        	readNumber = getReadNumberFromPileup(mateFileName);
	} else {
	        readNumber = getLineNumber(mateFileName, pathToSamtools_, pathToSambamba_, SambambaThreads_);
	}
	cout << "\t read number:\t" << readNumber << "\n";
	cout << "\t coefficientOfVariation:\t" << coefficientOfVariation << "\n";
	windowSize = round_f(float(1./(coefficientOfVariation*coefficientOfVariation)/readNumber*refGenomeSize_));
	cout << "\t evaluated window size:\t" << windowSize << "\n";
	for (int i = 0; i < (int) names.size(); i++) {
		ChrCopyNumber chrCopyNumber(windowSize, lengths[i],names[i]);
		chromosomesInd_.insert(pair<string, int> (names[i],i));
		chrCopyNumber_.push_back(chrCopyNumber);
	}
	//cout << "data structure created";
	//read mateFileName and calculate copyNumber

	GenomeCopyNumber::fillMyHash(mateFileName ,inputFormat, matesOrientation, windowSize, windowSize);

}

long GenomeCopyNumber::getNormalNumberOfPairs(void) {
	return normalNumberOfPairs_;
}

int GenomeCopyNumber::getWindowSize(void) {
	return windowSize_;
}

void GenomeCopyNumber::setStep(int step) {
    step_=step;
}

void GenomeCopyNumber::setNormalContamination(float normalContamination) {
    normalContamination_=normalContamination;
    vector<ChrCopyNumber>::iterator it;
	for ( it=chrCopyNumber_.begin() ; it != chrCopyNumber_.end(); it++ )
            it->setNormalContamination(normalContamination);
}

bool GenomeCopyNumber::ifHasBAF() {
    return hasBAF_;
}

void GenomeCopyNumber::setSex(std::string sex) {
    sex_=sex;
    //remove chrY if female:
    if (sex.compare("XX")==0) {
        vector<ChrCopyNumber>::iterator it;
        for ( it=chrCopyNumber_.begin() ; it != chrCopyNumber_.end(); it++ ) {
            string::size_type pos = 0;
            if ( ( pos = it->getChromosome().find("Y", pos)) != string::npos ) {
                chrCopyNumber_.erase(it);
                cout << "..Will not consider chrY..\n";

                int success = chromosomesInd_.erase ("Y");
                if (success) {
                    cout << "..Erased chrY from the list of chromosomes\n";
                } else {
                    cout << "..Was not able to erase chrY from the list of chromosomes\n";
                }
                return;
            }
        }

    }
}

void GenomeCopyNumber::setBAFtrue() {
    hasBAF_=1;
}

void GenomeCopyNumber::readCopyNumber(std::string const& mateFileName ,std::string const& inputFormat, std::string const& matesOrientation, std::string const& chrLenFileName,
                    int windowSize , int step, std::string targetBed ) {
    if (WESanalysis == false)  {
        if (step == NA)
            step = windowSize;
            //check if step value is correct:
        if (step <= 0 || step > windowSize) {
            cerr << "step  should be a positive interger value less than or equal to the window size\n";
            exit(-1);
        }
        step_=step;
    }
    else {
        step_ = step;
    }

	//reading the file with genome information
	std::vector<std::string> names;
	std::vector<int> lengths;
	readFileWithGenomeInfo(chrLenFileName, names, lengths);
	refGenomeSize_ = sum(lengths);
	for (int i = 0; i < (int) names.size(); i++) {
		ChrCopyNumber chrCopyNumber(windowSize, lengths[i],names[i], step, targetBed);
		chromosomesInd_.insert(pair<string, int> (names[i],i));
		chrCopyNumber_.push_back(chrCopyNumber);
	}
	//read mateFileName and calculate copyNumber


	GenomeCopyNumber::fillMyHash(mateFileName ,inputFormat, matesOrientation, windowSize, step, targetBed);


    //remove chromosomes that do to have reads or exons:


    if (targetBed != "") {
        vector <int> toErase;
        for (int i = 0; i < (int) names.size(); i++) {
            if (chrCopyNumber_[i].getExons_Countchr()==0) {
                toErase.push_back(i);
            }
        }
        for (int i=toErase.size()-1;i>=0;i--) {
            chrCopyNumber_.erase(chrCopyNumber_.begin()+toErase[i]);
            chromosomesInd_.erase(names[toErase[i]]);
        }
    }

}

void GenomeCopyNumber::initCopyNumber(std::string const& chrLenFileName, int windowSize , int step, std::string targetBed) {

  ThreadPoolManager::getInstance()->lock();
  if (WESanalysis == false)  {
    if (step == NA) {
	  step = windowSize;
	}
    //check if steo value is correct:
    if (step <= 0 || step > windowSize) {
        cerr << "step  should be a positive interger value less than or equal to the window size\n";
        exit(-1);
    }
    step_=step;
	windowSize_ = windowSize;
	//reading the file with genome information
	std::vector<std::string> names;
	std::vector<int> lengths;
	readFileWithGenomeInfo(chrLenFileName, names, lengths);
	refGenomeSize_ = sum(lengths);
	for (int i = 0; i < (int) names.size(); i++) {
		ChrCopyNumber chrCopyNumber(windowSize, lengths[i],names[i], step);
		chromosomesInd_.insert(pair<string, int> (names[i],i));
		chrCopyNumber_.push_back(chrCopyNumber);
	}
   }   else   {
    std::vector<std::string> names;
	std::vector<int> lengths;
	readFileWithGenomeInfo(chrLenFileName, names, lengths);
	refGenomeSize_ = sum(lengths);
	for (int i = 0; i < (int) names.size(); i++) {
		ChrCopyNumber chrCopyNumber(windowSize, lengths[i],names[i], step, targetBed);
		chromosomesInd_.insert(pair<string, int> (names[i],i));
		chrCopyNumber_.push_back(chrCopyNumber);
    }
   }
   ThreadPoolManager::getInstance()->unlock();
}

void GenomeCopyNumber::finishCopyNumber(long normalCount) {
	totalNumberOfPairs_ = normalCount;
	normalNumberOfPairs_ = normalCount;
	cout << normalNumberOfPairs_<< " reads used to compute copy number profile\n";
}

long GenomeCopyNumber::getTotalNumberOfPairs() {
	return totalNumberOfPairs_;
}

void GenomeCopyNumber::calculateBreakpoints(double breakPointThreshold, int breakPointType) {
	cout << "..Calculating breakpoints, breakPointThreshold = " <<breakPointThreshold<<"\n";
	ThreadPool* thrPool = ThreadPoolManager::getInstance()->newThreadPool("ChrCopyNumber_calculateBreakpoint");
	for (vector<ChrCopyNumber>::iterator it = chrCopyNumber_.begin(); it != chrCopyNumber_.end(); it++) {
	  ChrCopyNumber& chrCopyNumber = *it;
	  ChrCopyNumberCalculateBreakpointArgWrapper* bkpArg = new ChrCopyNumberCalculateBreakpointArgWrapper(chrCopyNumber, breakPointThreshold, breakPointType);
	  thrPool->addThread(ChrCopyNumber_calculateBreakpoint_wrapper, bkpArg);
	}

	thrPool->run();
	delete thrPool;
}

void GenomeCopyNumber::calculateBAFBreakpoints(double breakPointThreshold, int breakPointType) {
	cout << "..Calculating breakpoints for BAF, breakPointThreshold = " <<breakPointThreshold<<"\n";
#if 1
	ThreadPool* thrPool = ThreadPoolManager::getInstance()->newThreadPool("ChrCopyNumber_calculateBAFBreakpoint");
	for (vector<ChrCopyNumber>::iterator it = chrCopyNumber_.begin(); it != chrCopyNumber_.end(); it++) {
	  ChrCopyNumber& chrCopyNumber = *it;
	  ChrCopyNumberCalculateBreakpointArgWrapper* bkpArg = new ChrCopyNumberCalculateBreakpointArgWrapper(chrCopyNumber, breakPointThreshold, breakPointType);
	  thrPool->addThread(ChrCopyNumber_calculateBAFBreakpoint_wrapper, bkpArg);
	}

	thrPool->run();
	delete thrPool;
#else
	vector<ChrCopyNumber>::iterator it=chrCopyNumber_.begin();
	//calculate breakpoints for the first chromosome and get its length to normalize other graphs
	int firstChromLen = 0; //we don't know it
    if (it != chrCopyNumber_.end()) {
        cout << "..processing chromosome " <<it->getChromosome()<<"\n";
		firstChromLen = it->calculateBAFBreakpoints(breakPointThreshold, firstChromLen,breakPointType);
		if (firstChromLen == 0) {
            cerr << "..failed to run segmentation on chr" << it->getChromosome() << "\n..no chromosome length normalization will be performed\n";
		}
    }
    it++;
	for (  ; it != chrCopyNumber_.end(); it++ ) {
		cout << "..processing chromosome " <<it->getChromosome()<<"\n";
		int result = it->calculateBAFBreakpoints(breakPointThreshold, firstChromLen,breakPointType);
		if (result == 0) {
            cerr << "..failed to run segmentation on chr" << it->getChromosome() << "\n";
		}
	}
#endif
}

void GenomeCopyNumber::fillMyHash(std::string const& mateFileName ,std::string const& inputFormat_str, std::string const& matesOrientation_str, int windowSize , int step, std::string targetBed) {
	//read mateFileName and calculate copyNumber
	long count = 0;
	long normalCount = 0;
	int bin = 0;

	cout << "..Starting reading "<< mateFileName << "\n";
	//if (mateFileName.substr(mateFileName.size()-1,3).compare(".gz")==0) {
	//	igzstream in( mateFileName );
	//	if ( ! in.good()) {
	//		std::cerr << "ERROR: Opening file `" << mateFileName << "' failed.\n";
	//		return EXIT_FAILURE;
	//	}
	//	in.close();
	//}
	std::ifstream fileMates (mateFileName.c_str());
	vector<float> insertSizeVector;
	windowSize_ = windowSize;
	step_=step;
	string line;
	if (!fileMates.is_open()) {
	    cerr << "Error: unable to open "+mateFileName+"\n" ;
	    exit(-1);
	}

#ifdef PROFILE_TRACE
	time_t t0 = time(NULL);
#endif

	MateOrientation matesOrientation = getMateOrientation(matesOrientation_str);

	InputFormat inputFormat;
	char* line_buffer;
    if (inputFormat_str.compare("bam")==0 || inputFormat_str.compare("BAM")==0 || mateFileName.substr(mateFileName.size()-4,4).compare(".bam")==0 || mateFileName.substr(mateFileName.size()-3,3).compare(".gz")==0) {

        fileMates.close();

        FILE *stream;
        char buffer[MAX_BUFFER];
        string command;

        string myInputFormat=inputFormat_str;
        if (mateFileName.substr(mateFileName.size()-3,3).compare(".gz")!=0) {
                if (pathToSambamba_ != "")
                    {
                    command = pathToSambamba_ + " view -t " + SambambaThreads_ + " " + mateFileName;
                    myInputFormat="sam";       //will try to use existing sambamba
                    cout << "..sambamba should be installed to be able to read BAM files; will use the following command for sambamba: "<<pathToSambamba_ + " view -t " + SambambaThreads_ + " " + mateFileName << "\n";
                    }
                else
                    {
                    command = pathToSamtools_ + " view "+mateFileName;
                    myInputFormat="sam";       //will try to use existing samtools
                    cout << "..samtools should be installed to be able to read BAM files; will use the following command for samtools: "<<pathToSamtools_ + " view "+mateFileName<<"\n";
                    }
        }
        else {
            command = "gzip -cd "+mateFileName;
		}

		inputFormat = getInputFormat(myInputFormat);
        stream =
#if defined(_WIN32)
		  _popen(command.c_str(), "r");
#else
		popen(command.c_str(), "r");
#endif

		while ((line_buffer = getLine(buffer, MAX_BUFFER, stream, line)) != NULL) {
		  count++;
		  normalCount+=processRead(inputFormat,matesOrientation,line_buffer, bin,targetBed, mateFileName);
		}
        #if defined(_WIN32)
				_pclose(stream);
		#else
				pclose(stream);
		#endif
		cout << "..finished reading "<<mateFileName<<endl;
    } else {
 	        inputFormat = getInputFormat(inputFormat_str);
    		while (std::getline(fileMates,line)) {
                count++;
                if ((inputFormat_str.compare("bowtie")==0 || inputFormat_str.compare("Bowtie")==0)&&(matesOrientation_str.compare("0")!=0)){
                    string line2;
                    std::getline(fileMates,line2);
                    count++;
                    normalCount+=processReadWithBowtie(inputFormat_str,matesOrientation_str,line,line2);
                } else {
				    normalCount+=processRead(inputFormat,matesOrientation,line.c_str(), bin,targetBed, mateFileName);
				}
    		}
            fileMates.close();
    }

#ifdef PROFILE_TRACE
	std::cout << "PROFILING [tid=" << pthread_self() << "]: " << mateFileName << " read in " << (time(NULL)-t0) << " seconds [fillMyHash]\n" << std::flush;
#endif

	totalNumberOfPairs_ = normalCount;
	normalNumberOfPairs_ = normalCount;
	cout << count<< " lines read..\n";
	cout << normalNumberOfPairs_<< " reads used to compute copy number profile\n";
	if (normalNumberOfPairs_==0) {
        cerr << "\nError: FREEC was not able to extract reads from " << mateFileName;
        cerr << "\n\nCheck your parameters: inputFormat and mateOrientation\n";
        cerr << "Use \"matesOrientation=0\" if you have single end reads\n";
        cerr << "Check the list of possible input formats at http://bioinfo-out.curie.fr/projects/freec/tutorial.html#CONFIG\n\n";

        if ((inputFormat_str.compare("SAM")==0)||(inputFormat_str.compare("sam")==0)||(inputFormat_str.compare("BAM")==0)||(inputFormat_str.compare("bam")==0)) {
            if (matesOrientation_str.compare("0")!=0) {
                cerr << "If you use sorted SAM or BAM, please set \"mateOrientation=0\"; then FREEC will not try to detect pairs with normal orientation and insert size. Instead, it will keep all pairs from the input file\n\n";
            }
        }


        exit(-1);
	}
}

int GenomeCopyNumber::findIndex (std::string const& chr) {
	if (chromosomesInd_.find(chr) == chromosomesInd_.end()) {return NA;}
	return chromosomesInd_.find(chr)->second;
}

void GenomeCopyNumber::recalculateRatio (float contamination) {
	if (contamination >= 1 || contamination <= 0) {
        cerr << "contamination should be between 0 and 1\n";
        exit(-1);
	}
	vector<ChrCopyNumber>::iterator it;

    for ( it=chrCopyNumber_.begin() ; it != chrCopyNumber_.end(); it++ ) {
        if (sex_.compare("XY")==0 && (it->getChromosome().find("X")!=string::npos || it->getChromosome().find("Y")!=string::npos)) {
            //should take into account that normally one has only one copy of X and Y..
            it->recalculateRatioWithContam(contamination,0.5, isRatioLogged_);
	    } else
            it->recalculateRatioWithContam(contamination,1, isRatioLogged_);
	}

}

void GenomeCopyNumber::recalculateRatioUsingCG (int degree, bool intercept, float minExpectedGC, float maxExpectedGC) {

	if (degree > MAXDEGREE) {
        cerr << "polynomial degree should be < 10\n";
        exit (-1);
	}
	int maximalNumberOfCopies = ploidy_*2;
	float interval = float (0.01) ;

	//first guess about parameters
	const int npoints = degree+2;


	double around [MAXDEGREE+2];
	for (int i = 0; i<npoints; i++) {
		around[i] = minExpectedGC + (maxExpectedGC-minExpectedGC)/(npoints-1)*i; //0.55-0.35
	}  // for degree = 3 one will get : { 0.35, 0.40, 0.45, 0.5, 0.55 }; for 0.35 and 0.55


	double yValues [MAXDEGREE+2];
	for (int i = 0; i <npoints; i++)
		yValues[i] = calculateMedianRatioAround(interval, float(around[i]));
	int nvars = degree; //fit by cubic polynomial
	ap::real_2d_array xy;
	xy.setlength(npoints,nvars+1);
	for (int i = 0; i <npoints; i++) {
		xy(i,degree) = yValues[i];
		xy(i,degree-1) = around[i];
		for (int j = degree-2; j>=0; j--) {
			xy(i,j)=xy(i,j+1)*around[i];
		}
		/* this is equal to
		xy(i,0) = around[i]*around[i]*around[i];
		xy(i,1) = around[i]*around[i];
		xy(i,2) = around[i];
		xy(i,3) = yValues[i]; */
	}
	linearmodel lm;
	int info;
	lrreport ar;

    if  (intercept)
        lrbuild(xy,npoints,nvars,info,lm,ar);
	else
		lrbuildz(xy,npoints,nvars,info,lm,ar);
	if (info != 1) {
		cerr << "Error in the first linear regression (the first guess about parameters), code: " << info <<"\n";
	}
	ap::real_1d_array v;
    v.setlength(nvars+int(intercept));
	lrunpack(lm,v,nvars);
	//cout << v(0) << "\t" << v(1)<< "\t" <<v(2)<< "\t" <<v(3)<< "\n";
	double a[MAXDEGREE+1];
	for (int i = 0; i <degree; i++) {
			a[i] = v(i);
	} /* this is equal to
	double a0 = v(0);
	double a1 = v(1);
	double a2 = v(2);
	double a3 = v(3); */

    if  (intercept)
		a[degree] = v(degree);
	else
		a[degree] = 0;

	vector <float> y; //y ~ ax^2+bx+c
	vector <float> x;

	//fill x and y:
	vector<ChrCopyNumber>::iterator it;
	for ( it=chrCopyNumber_.begin() ; it != chrCopyNumber_.end(); it++ ) {
		if (! (it->getChromosome().find("X")!=string::npos || it->getChromosome().find("Y")!=string::npos)) {
             //use a threshold, but correct using notN profile
            if (it->getMappabilityLength()>0) {
                for (int i = 0; i< it->getLength(); i++) {
                    if ((it->getRatioAtBin(i)>0)&&(it->getMappabilityProfileAt(i)>minMappabilityPerWindow)) {
                        x.push_back(it->getCGprofileAt(i));
                        y.push_back(it->getRatioAtBin(i));
                    }
                }
            } else {
                for (int i = 0; i< it->getLength(); i++) {
                    if ((it->getRatioAtBin(i)>0)&&(it->getNotNprofileAt(i)>minMappabilityPerWindow)) {
                        x.push_back(it->getCGprofileAt(i));
                        y.push_back(it->getRatioAtBin(i));
                    }
                }
            }

        }
	}
	int maximalNumberOfIterations = 100;
    	float rmserror = runEM(x,y,a,degree,maximalNumberOfIterations,ploidy_,maximalNumberOfCopies, intercept, normalContamination_);
	if (rmserror == -1) {
		cerr << "Error in EM => unable to calculate normalized profile\n";
		return;
	}
	cout << "root mean square error = " << rmserror << "\n";
	if (degree == 3) {
		cout << "Y = " << a[0] << "*x*x*x+" << a[1] << "*x*x+" << a[2] << "*x+" << a[3] <<"\n";
	}else if (degree == 2) {
		cout << "Y = " << a[0] << "*x*x+" << a[1] << "*x+" << a[2]  <<"\n";
	}else if (degree == 1) {
		cout << "Y = " << a[0] << "*x+" << a[1] << "\n";
	}
	if (degree > 3) {
        cout << "Y = ";
        for (int i=0; i<degree;i++) {
            cout << a[i]<<"*x^" <<  degree-i <<"+" ;
        }
        cout << a[degree] <<"\n";
	}
	for ( it=chrCopyNumber_.begin() ; it != chrCopyNumber_.end(); it++ ) {
		it->recalculateRatio(a,degree);
	}

}
void GenomeCopyNumber::setAllNormal () {
	vector<ChrCopyNumber>::iterator it;

	for ( it=chrCopyNumber_.begin() ; it != chrCopyNumber_.end(); it++ ) {
		it->setAllNormal();
	}
}

int GenomeCopyNumber::calculateRatioUsingCG (bool intercept, float minExpectedGC, float maxExpectedGC) { //returns 1 if the number of iteration is less than max; otherwise 0

     //try degree 3 and 4, SELECT THE ONE WITH LESS ITTERATIONS:
    int minDegreeToTest = 3;
    int maxDegreeToTest = 4;
    int selectedDegree=minDegreeToTest;
    int maximalNumberOfIterations = 300;
    int bestNumberOfIterations=maximalNumberOfIterations;
   // int bestSqareError=INFINITY;
    int valueToReturn = 0;


    vector <float> y; //y ~ ax^2+bx+c
    vector <float> x;

    //check whether there are chromosomes without X and Y:
    vector<ChrCopyNumber>::iterator it;
    int countAutosomes = 0;
    for ( it=chrCopyNumber_.begin() ; it != chrCopyNumber_.end(); it++ ) {
        if (! (it->getChromosome().find("X")!=string::npos || it->getChromosome().find("Y")!=string::npos)) {
            countAutosomes++;
        }
    }
    //fill x and y:
    for ( it=chrCopyNumber_.begin() ; it != chrCopyNumber_.end(); it++ ) {

            if (countAutosomes>0 && !(it->getChromosome().find("X")!=string::npos || it->getChromosome().find("Y")!=string::npos) || countAutosomes==0) {
                // if uniqueMatch, do correction to mappability
                if (uniqueMatch) {
                    //use mappabilityProfile_ and correct
                    for (int i = 0; i< it->getLength(); i++) {
                        if ((it->getValueAt(i)>0)&&(it->getMappabilityProfileAt(i)>minMappabilityPerWindow)) {
                            x.push_back(it->getCGprofileAt(i));
                            y.push_back(it->getValueAt(i)/it->getMappabilityProfileAt(i));
                        }
                    }
                } else {
                    //use a threshold, but correct using notN profile
                    if (it->getMappabilityLength()>0) {
                        for (int i = 0; i< it->getLength(); i++) {
                            if ((it->getValueAt(i)>0)&&(it->getMappabilityProfileAt(i)>minMappabilityPerWindow)) {
                                x.push_back(it->getCGprofileAt(i));
                                y.push_back(it->getValueAt(i)/it->getNotNprofileAt(i));
                            }
                        }
                    } else {
                        for (int i = 0; i< it->getLength(); i++) {
                            if ((it->getValueAt(i)>0)&&(it->getNotNprofileAt(i)>minMappabilityPerWindow)) {
                                x.push_back(it->getCGprofileAt(i));
                                y.push_back(it->getValueAt(i)/it->getNotNprofileAt(i));
                            }
                        }
                    }
                }
            }
    }

    for (int degree=minDegreeToTest;degree<=maxDegreeToTest; degree++) {
        if (degree > MAXDEGREE) {
           break;
        }
        int maximalNumberOfCopies = ploidy_*2;
        float interval = float (0.01) ;

        //first guess about parameters
        const int npoints = degree+2;


        double around [MAXDEGREE+2];
        for (int i = 0; i<npoints; i++) {
            around[i] = minExpectedGC + (maxExpectedGC-minExpectedGC)/(npoints-1)*i; //0.55-0.35
        }  // for degree = 3 one will get : { 0.35, 0.40, 0.45, 0.5, 0.55 }; for 0.35 and 0.55


        double yValues [MAXDEGREE+2];
        for (int i = 0; i <npoints; i++)
            yValues[i] = calculateMedianAround(interval, float(around[i]));
        int nvars = degree; //fit by cubic polynomial
        ap::real_2d_array xy;
        xy.setlength(npoints,nvars+1);
        for (int i = 0; i <npoints; i++) {
            xy(i,degree) = yValues[i];
            xy(i,degree-1) = around[i];
            for (int j = degree-2; j>=0; j--) {
                xy(i,j)=xy(i,j+1)*around[i];
            }
            /* this is equal to
            xy(i,0) = around[i]*around[i]*around[i];
            xy(i,1) = around[i]*around[i];
            xy(i,2) = around[i];
            xy(i,3) = yValues[i]; */
        }
        linearmodel lm;
        int info;
        lrreport ar;

        if  (intercept)
            lrbuild(xy,npoints,nvars,info,lm,ar);
        else
            lrbuildz(xy,npoints,nvars,info,lm,ar);
        if (info != 1) {
            cerr << "Error in the first linear regression (the first guess about parameters), code: " << info <<"\n";
        }
        ap::real_1d_array v;
        v.setlength(nvars+int(intercept));
        lrunpack(lm,v,nvars);
        //cout << v(0) << "\t" << v(1)<< "\t" <<v(2)<< "\t" <<v(3)<< "\n";
        double a[MAXDEGREE+1];
        for (int i = 0; i <degree; i++) {
                a[i] = v(i);
        } /* this is equal to
        double a0 = v(0);
        double a1 = v(1);
        double a2 = v(2);
        double a3 = v(3); */

        if  (intercept)
            a[degree] = v(degree);
        else
            a[degree] = 0;


        int realNumberOfIterations = maximalNumberOfIterations;
        float rmserror = runEM(x,y,a,degree,realNumberOfIterations,ploidy_,maximalNumberOfCopies, intercept, normalContamination_);
        if (rmserror == -1) {
            cerr << "Error in EM => unable to calculate normalized profile\n";
            return 0;
        }
        cout << "root mean square error = " << rmserror << "\n";
        if (realNumberOfIterations<bestNumberOfIterations) {
            selectedDegree = degree;
            bestNumberOfIterations=realNumberOfIterations;
        }

    }

    int degree = selectedDegree;
	if (degree > MAXDEGREE) {
        cerr << "polynomial degree should be < 10\n";
        exit (-1);
	}
	int maximalNumberOfCopies = ploidy_*2;
	float interval = float (0.01) ;

	//first guess about parameters
	const int npoints = degree+2;


	double around [MAXDEGREE+2];
	for (int i = 0; i<npoints; i++) {
		around[i] = minExpectedGC + (maxExpectedGC-minExpectedGC)/(npoints-1)*i; //0.55-0.35
	}  // for degree = 3 one will get : { 0.35, 0.40, 0.45, 0.5, 0.55 }; for 0.35 and 0.55


	double yValues [MAXDEGREE+2];
	for (int i = 0; i <npoints; i++)
		yValues[i] = calculateMedianAround(interval, float(around[i]));
	int nvars = degree; //fit by cubic polynomial
	ap::real_2d_array xy;
	xy.setlength(npoints,nvars+1);
	for (int i = 0; i <npoints; i++) {
		xy(i,degree) = yValues[i];
		xy(i,degree-1) = around[i];
		for (int j = degree-2; j>=0; j--) {
			xy(i,j)=xy(i,j+1)*around[i];
		}
		/* this is equal to
		xy(i,0) = around[i]*around[i]*around[i];
		xy(i,1) = around[i]*around[i];
		xy(i,2) = around[i];
		xy(i,3) = yValues[i]; */
	}
	linearmodel lm;
	int info;
	lrreport ar;

    if  (intercept)
        lrbuild(xy,npoints,nvars,info,lm,ar);
	else
		lrbuildz(xy,npoints,nvars,info,lm,ar);
	if (info != 1) {
		cerr << "Error in the first linear regression (the first guess about parameters), code: " << info <<"\n";
	}
	ap::real_1d_array v;
    v.setlength(nvars+int(intercept));
	lrunpack(lm,v,nvars);
	//cout << v(0) << "\t" << v(1)<< "\t" <<v(2)<< "\t" <<v(3)<< "\n";
	double a[MAXDEGREE+1];
	for (int i = 0; i <degree; i++) {
			a[i] = v(i);
	} /* this is equal to
	double a0 = v(0);
	double a1 = v(1);
	double a2 = v(2);
	double a3 = v(3); */

    if  (intercept)
		a[degree] = v(degree);
	else
		a[degree] = 0;

    int realNumberOfIterations = maximalNumberOfIterations;

    float rmserror = runEM(x,y,a,degree,realNumberOfIterations,ploidy_,maximalNumberOfCopies, intercept, normalContamination_);
	if (rmserror == -1) {
		cerr << "Error in EM => unable to calculate normalized profile\n";
		return 0;
	}
	cout << "root mean square error = " << rmserror << "\n";
	if (degree == 3) {
		cout << "Y = " << a[0] << "*x*x*x+" << a[1] << "*x*x+" << a[2] << "*x+" << a[3] <<"\n";
	}else if (degree == 2) {
		cout << "Y = " << a[0] << "*x*x+" << a[1] << "*x+" << a[2]  <<"\n";
	}else if (degree == 1) {
		cout << "Y = " << a[0] << "*x+" << a[1] << "\n";
	}
    if (degree > 3) {
        cout << "Y = ";
        for (int i=0; i<degree;i++) {
            cout << a[i]<<"*x^" <<  degree-i <<"+" ;
        }
        cout << a[degree] <<"\n";
	}


	for ( it=chrCopyNumber_.begin() ; it != chrCopyNumber_.end(); it++ ) {
		it->calculateRatio(a,degree);
	}
    x.clear();
    y.clear();

    if (realNumberOfIterations<maximalNumberOfIterations)
        valueToReturn = 1;
    return valueToReturn;
}



int GenomeCopyNumber::calculateRatioUsingCG (int degree, bool intercept, float minExpectedGC, float maxExpectedGC) {

	if (degree > MAXDEGREE) {
        cerr << "polynomial degree should be < 10\n";
        exit (-1);
	}
	int maximalNumberOfCopies = ploidy_*2;
	float interval = float (0.01) ;

	//first guess about parameters
	const int npoints = degree+2;

    int valueToReturn = 0;

	double around [MAXDEGREE+2];
	for (int i = 0; i<npoints; i++) {
		around[i] = minExpectedGC + (maxExpectedGC-minExpectedGC)/(npoints-1)*i; //0.55-0.35
	}  // for degree = 3 one will get : { 0.35, 0.40, 0.45, 0.5, 0.55 }; for 0.35 and 0.55


	double yValues [MAXDEGREE+2];
	for (int i = 0; i <npoints; i++)
		yValues[i] = calculateMedianAround(interval, float(around[i]));
	int nvars = degree; //fit by cubic polynomial
	ap::real_2d_array xy;
	xy.setlength(npoints,nvars+1);
	for (int i = 0; i <npoints; i++) {
		xy(i,degree) = yValues[i];
		xy(i,degree-1) = around[i];
		for (int j = degree-2; j>=0; j--) {
			xy(i,j)=xy(i,j+1)*around[i];
		}
		/* this is equal to
		xy(i,0) = around[i]*around[i]*around[i];
		xy(i,1) = around[i]*around[i];
		xy(i,2) = around[i];
		xy(i,3) = yValues[i]; */
	}
	linearmodel lm;
	int info;
	lrreport ar;

    if  (intercept)
        lrbuild(xy,npoints,nvars,info,lm,ar);
	else
		lrbuildz(xy,npoints,nvars,info,lm,ar);
	if (info != 1) {
		cerr << "Error in the first linear regression (the first guess about parameters), code: " << info <<"\n";
	}
	ap::real_1d_array v;
    v.setlength(nvars+int(intercept));
	lrunpack(lm,v,nvars);
	//cout << v(0) << "\t" << v(1)<< "\t" <<v(2)<< "\t" <<v(3)<< "\n";
	double a[MAXDEGREE+1];
	for (int i = 0; i <degree; i++) {
			a[i] = v(i);
	} /* this is equal to
	double a0 = v(0);
	double a1 = v(1);
	double a2 = v(2);
	double a3 = v(3); */

    if  (intercept)
		a[degree] = v(degree);
	else
		a[degree] = 0;

	vector <float> y; //y ~ ax^2+bx+c
	vector <float> x;

	//fill x and y:
	vector<ChrCopyNumber>::iterator it;
	for ( it=chrCopyNumber_.begin() ; it != chrCopyNumber_.end(); it++ ) {

		if (! (it->getChromosome().find("X")!=string::npos || it->getChromosome().find("Y")!=string::npos)) {
            // if uniqueMatch, do correction to mappability
			if (uniqueMatch) {
                //use mappabilityProfile_ and correct
                for (int i = 0; i< it->getLength(); i++) {
                    if ((it->getValueAt(i)>0)&&(it->getMappabilityProfileAt(i)>minMappabilityPerWindow)) {
                        x.push_back(it->getCGprofileAt(i));
                        y.push_back(it->getValueAt(i)/it->getMappabilityProfileAt(i));
                    }
                }
			} else {
                //use a threshold, but correct using notN profile
                if (it->getMappabilityLength()>0) {
                    for (int i = 0; i< it->getLength(); i++) {
                        if ((it->getValueAt(i)>0)&&(it->getMappabilityProfileAt(i)>minMappabilityPerWindow)) {
                            x.push_back(it->getCGprofileAt(i));
                            y.push_back(it->getValueAt(i)/it->getNotNprofileAt(i));
                        }
                    }
                } else {
                    for (int i = 0; i< it->getLength(); i++) {
                        if ((it->getValueAt(i)>0)&&(it->getNotNprofileAt(i)>minMappabilityPerWindow)) {
                            x.push_back(it->getCGprofileAt(i));
                            y.push_back(it->getValueAt(i)/it->getNotNprofileAt(i));
                        }
                    }
                }
			}
        }
	}
	int maximalNumberOfIterations = 300;
	int realNumberOfIterations = maximalNumberOfIterations;
    float rmserror = runEM(x,y,a,degree,realNumberOfIterations,ploidy_,maximalNumberOfCopies, intercept, normalContamination_);
	if (rmserror == -1) {
		cerr << "Error in EM => unable to calculate normalized profile\n";
		return 0;
	}

	if (realNumberOfIterations<maximalNumberOfIterations)
        valueToReturn = 1;

	cout << "root mean square error = " << rmserror << "\n";
	if (degree == 3) {
		cout << "Y = " << a[0] << "*x*x*x+" << a[1] << "*x*x+" << a[2] << "*x+" << a[3] <<"\n";
	}else if (degree == 2) {
		cout << "Y = " << a[0] << "*x*x+" << a[1] << "*x+" << a[2]  <<"\n";
	}else if (degree == 1) {
		cout << "Y = " << a[0] << "*x+" << a[1] << "\n";
	}
    if (degree > 3) {
        cout << "Y = ";
        for (int i=0; i<degree;i++) {
            cout << a[i]<<"*x^" <<  degree-i <<"+" ;
        }
        cout << a[degree] <<"\n";
	}


	for ( it=chrCopyNumber_.begin() ; it != chrCopyNumber_.end(); it++ ) {
		it->calculateRatio(a,degree);
	}
    x.clear();
    y.clear();
    return valueToReturn;
}

void GenomeCopyNumber::setPloidy(int ploidy) {
	ploidy_ = ploidy;
	vector<ChrCopyNumber>::iterator it;
	for ( it=chrCopyNumber_.begin() ; it != chrCopyNumber_.end(); it++ )
            it->setPloidy(ploidy);

}


double GenomeCopyNumber::calculateMedianRatioAround (float interval, float around) {

	float maxCG = around+interval;
	float minCG = around-interval;

	vector <float> myValuesAround;
	vector<ChrCopyNumber>::iterator it;
	for ( it=chrCopyNumber_.begin() ; it != chrCopyNumber_.end(); it++ ) {
		if (! (it->getChromosome().find("X")!=string::npos || it->getChromosome().find("Y")!=string::npos))
			for (int i = 0; i< it->getLength(); i++) {
				if ((it->getCGprofileAt(i)<=maxCG)&&(it->getCGprofileAt(i)>=minCG))
					if (it->getRatioAtBin(i) >0) //non-zero values
						myValuesAround.push_back(it->getRatioAtBin(i));
			}
	}
    if (myValuesAround.size()==0) {
        cerr << "Error: zero reads in windows with the GC-content around "<<around<<"\n";
        cerr << "Unable to proceed..\n";
        cerr << "Try to rerun the program with higher number of reads\n";
        exit(-1);
    }
	float median = get_median(myValuesAround);
	myValuesAround.clear();
	return median;
}

int GenomeCopyNumber::calculateSDReadCountPerWindow(int mean) {
    vector <float> myValues;
	vector<ChrCopyNumber>::iterator it;
	for ( it=chrCopyNumber_.begin() ; it != chrCopyNumber_.end(); it++ ) {
		if (! (it->getChromosome().find("X")!=string::npos || it->getChromosome().find("Y")!=string::npos))
			for (int i = 0; i< it->getLength(); i++) {
				if (it->getValueAt(i) >0) //non-zero values
                    myValues.push_back(it->getValueAt(i));
			}
	}

    if (myValues.size()==0) {
            cerr << "Ersqrt(139000)ror: no windows with reads\n";
            cerr << "Unable to proceed..\n";
            exit(-1);
    }

	float sd = get_sd(myValues,mean);
	myValues.clear();
	return floor(sd);
}


int GenomeCopyNumber::calculateMedianReadCountPerWindow() {
    vector <float> myValues;
	vector<ChrCopyNumber>::iterator it;
	for ( it=chrCopyNumber_.begin() ; it != chrCopyNumber_.end(); it++ ) {
		if (! (it->getChromosome().find("X")!=string::npos || it->getChromosome().find("Y")!=string::npos))
			for (int i = 0; i< it->getLength(); i++) {
				if (it->getValueAt(i) >0) //non-zero values
                    myValues.push_back(it->getValueAt(i));
			}
	}

    if (myValues.size()==0) {
            cerr << "Error: no windows with reads\n";
            cerr << "Unable to proceed..\n";
            exit(-1);
    }

	float median = get_median(myValues);
	myValues.clear();
	return floor(median);
}

double GenomeCopyNumber::calculateMedianAround (float interval, float around) {

	float maxCG = around+interval;
	float minCG = around-interval;

	vector <float> myValuesAround;
	vector<ChrCopyNumber>::iterator it;


    int countAutosomes = 0;
    for ( it=chrCopyNumber_.begin() ; it != chrCopyNumber_.end(); it++ ) {
        if (! (it->getChromosome().find("X")!=string::npos || it->getChromosome().find("Y")!=string::npos)) {
            countAutosomes++;
        }
    }


	for ( it=chrCopyNumber_.begin() ; it != chrCopyNumber_.end(); it++ ) {
		if (countAutosomes>0 && ! (it->getChromosome().find("X")!=string::npos || it->getChromosome().find("Y")!=string::npos) || countAutosomes==0)
			for (int i = 0; i< it->getLength(); i++) {
				if ((it->getCGprofileAt(i)<=maxCG)&&(it->getCGprofileAt(i)>=minCG))
					if (it->getValueAt(i) >0) //non-zero values
						myValuesAround.push_back(it->getValueAt(i));
			}
	}
    if (myValuesAround.size()==0) {
        cerr << "Error: zero reads in windows with the GC-content around "<<around<< " with interval "<< interval<<", will try again with "<< interval*4<<"\n";
        interval=interval*4;
        maxCG = around+interval;
        minCG = around-interval;

        for ( it=chrCopyNumber_.begin() ; it != chrCopyNumber_.end(); it++ ) {
            if (! (it->getChromosome().find("X")!=string::npos || it->getChromosome().find("Y")!=string::npos))
                for (int i = 0; i< it->getLength(); i++) {
                    if ((it->getCGprofileAt(i)<=maxCG)&&(it->getCGprofileAt(i)>=minCG))
                        if (it->getValueAt(i) >0) //non-zero values
                            myValuesAround.push_back(it->getValueAt(i));
                }
        }
        if (myValuesAround.size()==0) {
            cerr << "Error: zero reads in windows with the GC-content around "<<around<< " with interval "<< interval<<"\n";
            cerr << "Unable to proceed..\n";
            cerr << "Try to rerun the program with higher number of reads\n";
            exit(-1);
        }
    }
	float median = get_median(myValuesAround);
	myValuesAround.clear();
	return median;
}

double GenomeCopyNumber::calculateMedianAround (GenomeCopyNumber & controlCopyNumber, float interval, float around) {

	double maxVal = around+interval;
	double minVal = around-interval;
    float median;

	vector <float> myValuesAround;
	vector<ChrCopyNumber>::iterator it;
	for ( it=chrCopyNumber_.begin() ; it != chrCopyNumber_.end(); it++ ) {
		if (! (it->getChromosome().find("X")!=string::npos || it->getChromosome().find("Y")!=string::npos)) {
			vector <float> controlcounts = controlCopyNumber.getChrCopyNumber(it->getChromosome()).getValues() ;
            //check that everything is all right:
            if (int(controlcounts.size())!=it->getLength()) {
                cerr << "Possible Error: calculateMedianAround ()\n";
            }
			for (int i = 0; i< it->getLength(); i++) {
				float controlValue = controlcounts[i];
				if ((controlValue<=maxVal)&&(controlValue>=minVal))
					myValuesAround.push_back(it->getValueAt(i));
			}
			controlcounts.clear();
        }
	}
    if (myValuesAround.size()==0) {
        cerr << "Warning: zero reads in windows with the Read Count around "<<around<<" in the Contol dataset\n";
        cerr << "May be your window is too small? which value of the window size do you obtain when you set coefficientOfVariation=0.05? \n";
        //exit(-1);
        return NA;
    } else {
        median = get_median(myValuesAround);
        myValuesAround.clear();
    }
	return median;

}

void GenomeCopyNumber::removeLowReadCountWindows(GenomeCopyNumber & controlCopyNumber,int RCThresh) {
    cout << "..will remove all windows with read count in the control less than "<<RCThresh<<"\n";
	vector<ChrCopyNumber>::iterator it;
    for ( it=chrCopyNumber_.begin() ; it != chrCopyNumber_.end(); it++ ) {
        it->removeLowReadCountWindows(controlCopyNumber.getChrCopyNumber(it->getChromosome()), RCThresh);
	}
}

void GenomeCopyNumber::removeLowReadCountWindowsFromControl(int RCThresh) {
    cout << "..will process the control file as well: removing all windows with read count in the control less than "<<RCThresh<<"\n";
	vector<ChrCopyNumber>::iterator it;
    for ( it=chrCopyNumber_.begin() ; it != chrCopyNumber_.end(); it++ ) {
        it->removeLowReadCountWindows(RCThresh);
	}
}

int GenomeCopyNumber::getNumberOfChromosomes() {
    return chrCopyNumber_.size();
}

bool GenomeCopyNumber::ifHasRatio() {
    if (chrCopyNumber_[0].ifHasRatio()) return 1;
    return 0;
}

long double GenomeCopyNumber::calculateRSS(int ploidy)
{
    string::size_type pos = 0;
    vector<float> observedvalues;
    vector<float> expectedvalues;
	map<string,int>::iterator it;
	for (it=chromosomesInd_.begin() ; it != chromosomesInd_.end(); it++ ) {
		string chrNumber = (*it).first;
		if ( ( pos = chrNumber.find("chr")) != string::npos )
			chrNumber.replace( pos, 3, "" );
        if ( ( pos = chrNumber.find("X")) != string::npos ) 		//exclude X and Y from the analysis
            continue;
        if ( ( pos = chrNumber.find("Y")) != string::npos )
            continue;
		int index = findIndex(chrNumber);
		int length = chrCopyNumber_[index].getLength();
		for (int i = 0; i< length; i++) {
            float observed = chrCopyNumber_[index].getRatioAtBin(i);
            if (observed!=NA) {
                if (isRatioLogged_) {
                    observed=pow(2,observed);
                }
                float expected = observed;
                if (chrCopyNumber_[index].isMedianCalculated()) {
                    expected = chrCopyNumber_[index].getMedianProfileAtI(i);
                    if (chrCopyNumber_[index].isSmoothed())
                        expected = chrCopyNumber_[index].getSmoothedProfileAtI(i);

                }
                observedvalues.push_back(observed);
                expectedvalues.push_back(expected);
            }
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
    if (observedvalues.size()==0) {
        return 0;
    }
    double normRSS = (RSS/observedvalues.size());
    observedvalues.clear();expectedvalues.clear();
    return normRSS;
}


void GenomeCopyNumber::calculateRatioUsingCG( GenomeCopyNumber & controlCopyNumber) {
    //since the raio should be already normalized, just devide sample/control
	vector<ChrCopyNumber>::iterator it;

	for ( it=chrCopyNumber_.begin() ; it != chrCopyNumber_.end(); it++ ) {
        if (sex_.compare("XY")==0 && (it->getChromosome().find("X")!=string::npos || it->getChromosome().find("Y")!=string::npos)) {
            //should take into account that normally one has only one copy of X and Y..
            it->recalculateRatio(controlCopyNumber.getChrCopyNumber(it->getChromosome()));
            it->recalculateRatio(2);
	    } else
            it->recalculateRatio(controlCopyNumber.getChrCopyNumber(it->getChromosome()));
	}
}

void GenomeCopyNumber::calculateRatioUsingCG_Regression( GenomeCopyNumber & controlCopyNumber) {

    int degree = 1;
    bool intercept=0;
    cout << "..polynomial degree set to 1\n";
    cout << "Initial guess for polynomial:\n";

    //first guess about parameters

    double a[MAXDEGREE+1];

    a[degree] = 0;
    a[0]=1;

    cout << "Y = " << a[0] << "*x+" << a[1] << "\n";


    int maximalNumberOfIterations = 300;
    int maximalNumberOfCopies = ploidy_*2;

    vector <float> y; //y ~ a0x+a1
    vector <float> x;
    vector <int> positions;

    //fill x and y:
    vector<ChrCopyNumber>::iterator it;
    for ( it=chrCopyNumber_.begin() ; it != chrCopyNumber_.end(); it++ ) {
        if (! (it->getChromosome().find("X")!=string::npos || it->getChromosome().find("Y")!=string::npos)) {
            vector <float> controlcounts = controlCopyNumber.getChrCopyNumber(it->getChromosome()).getRatio() ;
            //check that everything is all right:
            if (int(controlcounts.size())!=it->getLength()) {
                cerr << "Possible Error: calculateMedianAround ()\n";
            }
            for (int i = 0; i< it->getLength(); i++) {
                if (it->getRatioAtBin(i)>0 && controlcounts[i]>0) {
                    x.push_back(controlcounts[i]);
                    y.push_back(it->getRatioAtBin(i));
                    positions.push_back(it->getCoordinateAtBin(i));
                }
            }
            controlcounts.clear();
        }
    }
    if (degree==1) {
        const char * nametmp = "/home/vboeva/NAS_public/data/projects/FREEC/TEST_FREEC_CLBMA_WES/output_WES/xy.txt";
        std::ofstream file;
        file.open(nametmp);
        for ( unsigned int i=0 ;i<x.size(); i++ ) {
             file << positions[i]<<"\t"<<x[i] <<"\t" << y[i] <<"\n" ;
        }
        file.close();
    }
    float rmserror = runEM(x,y,a,degree,maximalNumberOfIterations,ploidy_,maximalNumberOfCopies, intercept, normalContamination_);
    if (rmserror == -1) {
        cerr << "Error in EM => unable to calculate normalized profile\n";
        return;
    }
    cout << "root mean square error = " << rmserror << "\n";

    cout << "Y = " << a[0] << "*x+" << a[1] << "\n";


    //devide sample/control with the identified constant that should not be too far away from 1:

	for ( it=chrCopyNumber_.begin() ; it != chrCopyNumber_.end(); it++ ) {
            it->recalculateRatio(controlCopyNumber.getChrCopyNumber(it->getChromosome()));
            it->recalculateRatio(a[0]);
	}

	for ( it=chrCopyNumber_.begin() ; it != chrCopyNumber_.end(); it++ ) {
        if (sex_.compare("XY")==0 && (it->getChromosome().find("X")!=string::npos || it->getChromosome().find("Y")!=string::npos)) {
            //should take into account that normally one has only one copy of X and Y..
            it->recalculateRatio(2);
	    }
	}

	//XXX
}

int GenomeCopyNumber::fillInRatio() {
    vector <float> countValues;
    vector<ChrCopyNumber>::iterator it;

    for ( it=chrCopyNumber_.begin() ; it != chrCopyNumber_.end(); it++ ) {
        if ((sex_.compare("XY")==0) && (it->getChromosome().find("X")!=string::npos || it->getChromosome().find("Y")!=string::npos)) {
                //should take into account that normally one has only one copy of X and Y..
                it->fillInRatio(isRatioLogged_);
        } else {
                it->fillInRatio(isRatioLogged_);
                for (int i = 0; i< it->getLength(); i++) {
                    if (it->getValueAt(i)>0) {
                        countValues.push_back(it->getRatioAtBin(i));
                    }
                }
        }
    }
    float median=get_medianNotNA(countValues);
    for ( it=chrCopyNumber_.begin() ; it != chrCopyNumber_.end(); it++ ) {
        if (!isRatioLogged_) {
                it->recalculateRatio(median);
        } else {
            it->recalculateLogRatio(median);
        }
    }
    return 1;
}

int GenomeCopyNumber::calculateRatio( GenomeCopyNumber & controlCopyNumber, int degree, bool intercept) {


    int maximalNumberOfIterations = 300;
    int maximalNumberOfCopies = ploidy_*2;
    int realNumberOfIterations = maximalNumberOfIterations;

    int successfulFit = 0;

	if(isRatioLogged_) {
        intercept=1; degree=1;//because it is loglogscale
        vector <float> y; //y ~ a0x+a1
        vector <float> x;

        //fill x and y:
        vector<ChrCopyNumber>::iterator it;
        for ( it=chrCopyNumber_.begin() ; it != chrCopyNumber_.end(); it++ ) {
            if (! (it->getChromosome().find("X")!=string::npos || it->getChromosome().find("Y")!=string::npos)) {
                vector <float> controlcounts = controlCopyNumber.getChrCopyNumber(it->getChromosome()).getValues() ;
                //check that everything is all right:
                if (int(controlcounts.size())!=it->getLength()) {
                    cerr << "Possible Error: calculateMedianAround ()\n";
                }
                for (int i = 0; i< it->getLength(); i++) {
                    if (it->getValueAt(i)>0 && controlcounts[i]>0) {
                        x.push_back(log(controlcounts[i]));
                        y.push_back(log(it->getValueAt(i)));
                    }
                }
                controlcounts.clear();
            }
        }
        cout << "Initial guess for polynomial:\n";

        int nvars = degree; //1 if the fit is linear
		ap::real_2d_array xy;
		int npoints = x.size();
		xy.setlength(npoints,nvars+1);
		int pos = 0;
		for (int i = 0; i <(int)x.size(); i++) {
				xy(pos,degree) = y[i];
				xy(pos,degree-1) = x[i];;
				for (int j = degree-2; j>=0; j--) {
					xy(pos,j)=xy(pos,j+1)*x[i];
				}
				pos++;
		}
		linearmodel lm;
		int info;
		lrreport ar;
		if  (intercept)
			lrbuild(xy,npoints,nvars,info,lm,ar);
		else
			lrbuildz(xy,npoints,nvars,info,lm,ar);

        if (info != 1) {
			cerr << "Error in linear regression, code: " << info <<"\n";
		}
		ap::real_1d_array v;
		v.setlength(nvars+1);
		lrunpack(lm,v,nvars);
		double a[MAXDEGREE+1];

		for (int i = 0; i <degree; i++) {
			a[i] = v(i);
		}
		if (intercept)
			a[degree] = v(degree); //intercept


		if (degree == 1) {
            cout << "log(Y) = " << a[0] << "*log(x)+" << a[1] << "\n";
        }
        if (degree > 1) {
            cout << "log(Y) = ";
            for (int i=0; i<degree;i++) {
                cout << a[i]<<"*log(x)^" <<  degree-i <<"+" ;
            }
            cout << a[degree] <<"\n";
        }
        float rmserror = runEMlog(x,y,a,degree,realNumberOfIterations,maximalNumberOfIterations,ploidy_,maximalNumberOfCopies, intercept, normalContamination_);
        if (rmserror == -1) {
            cerr << "Error in EM => unable to calculate normalized profile\n";
            return 0;
        }
        cout << "root mean square error = " << rmserror << "\n";

        if (realNumberOfIterations<maximalNumberOfIterations) {
            successfulFit = 1;
        }

		if (degree == 1) {
            cout << "log(Y) = " << a[0] << "*log(x)+" << a[1] << "\n";
        }
        if (degree > 1) {
            cout << "log(Y) = ";
            for (int i=0; i<degree;i++) {
                cout << a[i]<<"*log(x)^" <<  degree-i <<"+" ;
            }
            cout << a[degree] <<"\n";
        }

        for ( it=chrCopyNumber_.begin() ; it != chrCopyNumber_.end(); it++ ) {
            if ((sex_.compare("XY")==0) && (it->getChromosome().find("X")!=string::npos || it->getChromosome().find("Y")!=string::npos)) {
                //should take into account that normally one has only one copy of X and Y..
                it->calculateRatioLog(controlCopyNumber.getChrCopyNumber(it->getChromosome()),a,degree);
                it->recalculateRatio(2);
            } else {
                it->calculateRatioLog(controlCopyNumber.getChrCopyNumber(it->getChromosome()),a,degree);
            }
        }


	} else {

        if (WESanalysis == true){
            degree = 1;
            cout << "..polynomial degree set to 1\n";
        }
        cout << "Initial guess for polynomial:\n";

        //first guess about parameters
        const int npoints = degree+2;

        int medianReadCountPerWindowForControl = controlCopyNumber.calculateMedianReadCountPerWindow();
        int sdReadCountPerWindowForControl = controlCopyNumber.calculateSDReadCountPerWindow(medianReadCountPerWindowForControl);

        double around [MAXDEGREE+2] ;

		int lquar, uquar;
		controlCopyNumber.calculateReadCountQuartiles(lquar, uquar);
		double minVal = lquar;
		double maxVal = uquar;

        for (int i = 0; i<npoints; i++) {
            around[i] = minVal + (maxVal-minVal)/(npoints-1)*i;
        }
        // instead of { meanReadCountPerWindowForControl,meanReadCountPerWindowForControl-2*sqrt(double(meanReadCountPerWindowForControl)),meanReadCountPerWindowForControl+2*sqrt(double(meanReadCountPerWindowForControl)) };

        //float interval = float (5) ;
        float interval = (maxVal-minVal)/(npoints-1)/2 ;

        //check if maxVal and minVal are chosen all right:
        bool valuesAreAllRight = 1;
        for (int i = 0; i <npoints; i=i+npoints-1) {
            double value = calculateMedianAround(controlCopyNumber, interval, float(around[i]));
            if (value==NA)
                valuesAreAllRight = 0;
        }
        if (!valuesAreAllRight){
        //we need to set other values of max and min..
            minVal = max(medianReadCountPerWindowForControl-sdReadCountPerWindowForControl,(int)floor(medianReadCountPerWindowForControl-sqrt(medianReadCountPerWindowForControl)));
            maxVal = min(medianReadCountPerWindowForControl+sdReadCountPerWindowForControl,(int)floor(medianReadCountPerWindowForControl+sqrt(medianReadCountPerWindowForControl)));
            for (int i = 0; i<npoints; i++) {
                around[i] = minVal + (maxVal-minVal)/(npoints-1)*i;
            }
            interval = (maxVal-minVal)/(npoints-1)/2 ;
            valuesAreAllRight = 1;
            for (int i = 0; i <npoints; i=i+npoints-1) {
                double value = calculateMedianAround(controlCopyNumber, interval, float(around[i]));
                if (value==NA)
                    valuesAreAllRight = 0;
            }
            if (!valuesAreAllRight){
                cerr <<"Error: variation in read count per window is too small.\n";
                cerr <<"Unable to proceed..\n";
                exit(-1);

            }

        }


        double yValues [MAXDEGREE+2];
        for (int i = 0; i <npoints; i++)
            yValues[i] = calculateMedianAround(controlCopyNumber, interval, float(around[i]));
        int nvars = degree;
        ap::real_2d_array xy;
        xy.setlength(npoints,nvars+1);
        for (int i = 0; i <npoints; i++) {
            xy(i,degree-1) = around[i];
            xy(i,degree) = yValues[i];
            for (int j = degree-2; j>=0; j--) {
                xy(i,j)=xy(i,j+1)*around[i];
            }
        }
        linearmodel lm;
        int info;
        lrreport ar;
        if  (intercept)
            lrbuild(xy,npoints,nvars,info,lm,ar);
        else
            lrbuildz(xy,npoints,nvars,info,lm,ar);
        if (info != 1) {
            cerr << "Error in the first linear regression (the first guess about parameters), code: " << info <<"\n";
        }
        ap::real_1d_array v;
        v.setlength(nvars+int(intercept));
        lrunpack(lm,v,nvars);
        //cout << v(0) << "\t" << v(1)<< "\t" <<v(2)<< "\t" <<v(3)<< "\n";
        double a[MAXDEGREE+1];
        for (int i = 0; i <degree; i++) {
                a[i] = v(i);
        }
        if  (intercept)
            a[degree] = v(degree);
        else
            a[degree] = 0;
        vector <float> y; //y ~ a0x+a1
        vector <float> x;


        if (degree == 3) {
            cout << "Y = " << a[0] << "*x*x*x+" << a[1] << "*x*x+" << a[2] << "*x+" << a[3] <<"\n";
        }else if (degree == 2) {
            cout << "Y = " << a[0] << "*x*x+" << a[1] << "*x+" << a[2]  <<"\n";
        }else if (degree == 1) {
            cout << "Y = " << a[0] << "*x+" << a[1] << "\n";
        }
        if (degree > 3) {
            cout << "Y = ";
            for (int i=0; i<degree;i++) {
                cout << a[i]<<"*x^" <<  degree-i <<"+" ;
            }
            cout << a[degree] <<"\n";
        }

        //fill x and y:
        vector<ChrCopyNumber>::iterator it;
        for ( it=chrCopyNumber_.begin() ; it != chrCopyNumber_.end(); it++ ) {
            if (! (it->getChromosome().find("X")!=string::npos || it->getChromosome().find("Y")!=string::npos)) {
                vector <float> controlcounts = controlCopyNumber.getChrCopyNumber(it->getChromosome()).getValues() ;
                //check that everything is all right:
                if (int(controlcounts.size())!=it->getLength()) {
                    cerr << "Possible Error: calculateMedianAround ()\n";
                }
                for (int i = 0; i< it->getLength(); i++) {
                    if (it->getValueAt(i)>0) {
                        x.push_back(controlcounts[i]);
                        y.push_back(it->getValueAt(i));
                    }
                }
                controlcounts.clear();
            }
        }
	//const char * nametmp = "/bioinfo/users/vboeva/Desktop/TMP/Lena/patientT/xy.txt";
	//std::ofstream file;
	//file.open(nametmp);
	//for ( int i=0 ;i<x.size(); i++ ) {
	//	 file << x[i] <<"\t" << y[i] <<"\n" ;
	//}
	//file.close(); exit(0);

        float rmserror = runEM(x,y,a,degree,realNumberOfIterations,ploidy_,maximalNumberOfCopies, intercept, normalContamination_);
        if (rmserror == -1) {
            cerr << "Error in EM => unable to calculate normalized profile\n";
            return 0;
        }
        cout << "root mean square error = " << rmserror << "\n";

        if (realNumberOfIterations<maximalNumberOfIterations) {
            successfulFit = 1;
        }

        if (degree == 3) {
            cout << "Y = " << a[0] << "*x*x*x+" << a[1] << "*x*x+" << a[2] << "*x+" << a[3] <<"\n";
        }else if (degree == 2) {
            cout << "Y = " << a[0] << "*x*x+" << a[1] << "*x+" << a[2]  <<"\n";
        }else if (degree == 1) {
            cout << "Y = " << a[0] << "*x+" << a[1] << "\n";
        }

        if (degree > 3) {
            cout << "Y = ";
            for (int i=0; i<degree;i++) {
                cout << a[i]<<"*x^" <<  degree-i <<"+" ;
            }
            cout << a[degree] <<"\n";
        }
        for ( it=chrCopyNumber_.begin() ; it != chrCopyNumber_.end(); it++ ) {
            if ((sex_.compare("XY")==0) && (it->getChromosome().find("X")!=string::npos || it->getChromosome().find("Y")!=string::npos)) {
                //should take into account that normally one has only one copy of X and Y..
                it->calculateRatio(controlCopyNumber.getChrCopyNumber(it->getChromosome()),a, degree);
                it->recalculateRatio(2);
            } else {
                it->calculateRatio(controlCopyNumber.getChrCopyNumber(it->getChromosome()),a, degree);
            }
        }
    }
    return successfulFit;
}

void GenomeCopyNumber::calculateReadCountQuartiles(int& lower, int& upper) {

	vector <float> myValues;
	vector<ChrCopyNumber>::iterator it;
	for ( it=chrCopyNumber_.begin() ; it != chrCopyNumber_.end(); it++ ) {
		if (! (it->getChromosome().find("X")!=string::npos || it->getChromosome().find("Y")!=string::npos))
			for (int i = 0; i< it->getLength(); i++) {
				if (it->getValueAt(i) >0) //non-zero values
					myValues.push_back(it->getValueAt(i));
			}
	}

	if (myValues.size()==0) {
		cerr << "Error: no windows with reads\n";
		cerr << "Unable to proceed..\n";
		exit(-1);
	}

	sort(myValues.begin(), myValues.end());
	lower = myValues[myValues.size() / 4];
	upper = myValues[(3 * myValues.size()) / 4];
}

float GenomeCopyNumber::getMedianRatio() {
	vector<float>selectedValues;
	vector<ChrCopyNumber>::iterator it;
	for ( it=chrCopyNumber_.begin() ; it != chrCopyNumber_.end(); it++ ) {
		vector<float> chr_values = it->getRatio();
		for (int i = 0; i<(int)chr_values.size(); i++)
			if (chr_values[i]!=NA)
				selectedValues.push_back(chr_values[i]);
			else
				cout << chr_values[i] <<"\n";
	}
	float median = get_median(selectedValues);
	selectedValues.clear();
	return median;
}

float GenomeCopyNumber::calculateNormalizationConstant(GenomeCopyNumber & controlCopyNumber) {
	float median,medianControl;
	median = getMedianCopyNumber();
	medianControl = controlCopyNumber.getMedianCopyNumber();
	return median/medianControl;;
}

float GenomeCopyNumber::getMedianCopyNumber() {
	vector<float>selectedValues;
	vector<ChrCopyNumber>::iterator it;
	for ( it=chrCopyNumber_.begin() ; it != chrCopyNumber_.end(); it++ ) {
		vector<float> chr_values = it->getValues();
		for (int i = 0; i<(int)chr_values.size(); i++)
			//if (chr_values[i]!=0)
				selectedValues.push_back(chr_values[i]);
	}
	float median = get_median(selectedValues);
	//float mean = get_mean(selectedValues);
	selectedValues.clear();
    return median;
}

ChrCopyNumber GenomeCopyNumber::getChrCopyNumber(std::string const& chr) {
	int index = findIndex(chr);
	if (index==NA) throw "No such chromosome";
	return chrCopyNumber_[index];
}
void GenomeCopyNumber::printRatio(std::string const& outFile, bool ifBedGraphOutPut, bool printNA) {
    std::ofstream file (outFile.c_str());
    map<string,int>::iterator it;

    if (ifBedGraphOutPut==false) {
        file << "Chromosome\tStart\tRatio\tMedianRatio\tCopyNumber";


        if (hasBAF_==1){
            file << "\tBAF\testimatedBAF\tGenotype\tUncertaintyOfGT";
        }
        if (WESanalysis == true)
            {
            file << "\tGene";
            }
        if (SeekingSubc_ == true)
            {
            file << "\tSubclone_CN\tSubclone_Population";
            }
        file << "\n";
        for ( it=chromosomesInd_.begin() ; it != chromosomesInd_.end(); it++ ) {
            printRatio((*it).first,file,printNA);
        }

    } else {
        //print a track for LOH if hasBAF_
        if (hasBAF_) {
            file << "track type=bedGraph name=\"copy neutral LOH from "<<outFile<<"\" description=\"Copy Neutral LOH\" color=100,100,100\n";
            for ( it=chromosomesInd_.begin() ; it != chromosomesInd_.end(); it++ ) {
                printRatioBedGraph((*it).first,file, "LOH");
            }
        }
        file << "track type=bedGraph name=\"gains from "<<outFile<<"\" description=\"Gains\" color=255,40,20\n";
        for ( it=chromosomesInd_.begin() ; it != chromosomesInd_.end(); it++ ) {
            printRatioBedGraph((*it).first,file, "gain");
        }
        file << "track type=bedGraph name=\"losses from "<<outFile<<"\" description=\"Losses\" color=20,40,255\n";
        for ( it=chromosomesInd_.begin() ; it != chromosomesInd_.end(); it++ ) {
            printRatioBedGraph((*it).first,file, "loss");
        }
        file << "track type=bedGraph name=\"neutral copy number from "<<outFile<<"\" description=\"Neutral Copy Number\" color=20,255,40\n";
        for ( it=chromosomesInd_.begin() ; it != chromosomesInd_.end(); it++ ) {
            printRatioBedGraph((*it).first,file, "normal");
        }
    }
    file.close();
}

void GenomeCopyNumber::printBAF(std::string const& outFile, SNPinGenome& snpingenome, std::string const& matefile) {
    string myName = outFile;
    string outfile = myName + "_BAF.txt";
	std::ofstream file (outfile.c_str());
	map<string,int>::iterator it;
	file << "Chromosome\tPosition\tBAF\tFittedA\tFittedB\tA\tB\tuncertainty";
	file << "\n";
	for ( it=chromosomesInd_.begin() ; it != chromosomesInd_.end(); it++ ) {

		string chrNumber = (*it).first;

        int indexSNP = snpingenome.findIndex(chrNumber);
        if (indexSNP == NA) {
		    cerr << "An error occurred in GenomeCopyNumber::printBAF: could not find an SNP index for "<<chrNumber<<"\n";
		    exit(-1);
        }
        printBAF(chrNumber,file,snpingenome.SNP_atChr(indexSNP), myName, matefile);


	}
	file.close();
}


void GenomeCopyNumber::deleteFlanks(int telo_centromeric_flanks) {
	telo_centromeric_flanks_ = telo_centromeric_flanks;
	vector<ChrCopyNumber>::iterator it;
	for ( it=chrCopyNumber_.begin() ; it != chrCopyNumber_.end(); it++ ) {
		it->deleteFlanks(telo_centromeric_flanks);
	}
}

void GenomeCopyNumber::recalcFlanks(int telo_centromeric_flanks, int minNumberOfWindows) {
	telo_centromeric_flanks_ = telo_centromeric_flanks;
	vector<ChrCopyNumber>::iterator it;
	for ( it=chrCopyNumber_.begin() ; it != chrCopyNumber_.end(); it++ ) {
		it->recalcFlanks(telo_centromeric_flanks, minNumberOfWindows);
	}
}


void GenomeCopyNumber::calculateCopyNumberProbs_and_exomeLength(int breakPointType)
    {
	//estimationOfExomeSize_ = 0;
	vector<ChrCopyNumber>::iterator it;
	unsigned long long count = 0;
    int endsSize=0;
	string NormalBAF,NormalBAF_XY;
    float normalXYploidy=ploidy_*0.5; //will only use it when the genome is male
    estimationOfGenomeSize_=0;
    if (hasBAF_)
        {
            NormalBAF = getNormalBAFforPloidy(ploidy_);
            NormalBAF_XY=getXYBAFforPloidy(ploidy_);//will only use it when the genome is male
        }

    CNVs_.clear();
    copyNumberProbs_.clear();

	for ( it=chrCopyNumber_.begin() ; it != chrCopyNumber_.end(); it++ )
        {
		float previousLevel = NA;
		int lengthOfPreviousLevel = 0;
		float nextLevel = NA;
		int lengthOfNextLevel = 0;
		it->setIsSmoothed(true);
		string chr = it->getChromosome();
		cout << "..Annotation of CNVs for "<< chr <<"\n";
		//will use these three variables to collect CNVs:
		int start = 0;
		int end = NA;
		int cnumber = NA;
		string BAFSym = "-";
        string BAFprev = "-";
        string BAFnext = "-";
        string lBAF = "-";

        float lUncertainty = NA;
        float BAFUncertainty = NA;

        float normalLevel=1;
        if (sex_.compare("XY")==0 && (it->getChromosome().find("X")!=string::npos || it->getChromosome().find("Y")!=string::npos)) {
            normalLevel = 0.5;
        }
		for (int i = 0; i<(int)it->getMedianValues().size(); i++) {
			float level = it->getLevelAt(i, ploidy_);
			int fragmentLength =  it->getFragmentLengths()[i];
            if (hasBAF_) {
                BAFSym = it->getBAFsymbPerFrg(i);
                BAFUncertainty = it->getEstimatedBAFuncertaintyAtBin(i);
//                cout << "..Control: read "<< BAFSym << "\n";
            }

			if (level == NA)
                {
				int nextIndex = it->nextNoNAIndex(i, ploidy_,0);
				BAFUncertainty = NA;
				if (nextIndex == NA) {
					nextLevel = NA;
					lengthOfNextLevel = 0;
					BAFnext = "-";
				} else {
					nextLevel = it->getLevelAt(nextIndex, ploidy_);
					lengthOfNextLevel = it->getFragmentLengthsAt(nextIndex);
                    if (hasBAF_)
                        BAFnext = it->getBAFsymbPerFrg(nextIndex);
				}
				//nextLevel = round_by_ploidy(it->nextNoNAMedian(i, ploidy_), ploidy_);
				//lengthOfNextLevel = it->nextNoNALength(i, ploidy_);

                if (breakPointType==NOCALL) {
                    if ((previousLevel != NA) && (nextLevel != NA) && previousLevel==nextLevel && previousLevel==normalLevel) {
                        if (BAFprev.compare(BAFnext)==0) {
                            level = previousLevel;
                        }
                    }
                } else {
                    if (breakPointType==LARGECLOSE || breakPointType==SIMPLERIGHT || breakPointType==HALFLENGTH) {
                        if ((previousLevel != NA) && (nextLevel != NA) ) {
                            if (lengthOfPreviousLevel>lengthOfNextLevel) {
                                level = previousLevel;
                                if (hasBAF_) {
                                    BAFSym = BAFprev;
                                }
                            } else {
                                level = nextLevel;
                                if (hasBAF_)
                                    BAFSym = BAFnext;
                            }
                            //if (previousLevel == 1) //priority to "normal" copy number
                            //	level = 1;
                            //if (nextLevel == 1)
                            //	level = 1;
                        } else
                            if (nextLevel == NA) {
                                level = previousLevel;
                                if (hasBAF_)
                                    BAFSym = BAFprev;
                            } else {
                                level = nextLevel;
                                if (hasBAF_)
                                    BAFSym = BAFnext;
                            }
                    } else if (breakPointType==NORMALLEVEL) {
                        if ((previousLevel == normalLevel) || (nextLevel == normalLevel) ) {
                            level = normalLevel;
                            if (BAFprev.compare(BAFnext)==0)
                                BAFSym = BAFnext;
                            else
                                BAFSym = "-";

                        } else {
                            if ((previousLevel != NA) && (nextLevel != NA) ) {
                                if (lengthOfPreviousLevel>lengthOfNextLevel) {
                                    level = previousLevel;
                                    if (hasBAF_)
                                        BAFSym = BAFprev;
                                } else {
                                    level = nextLevel;
                                    if (hasBAF_)
                                        BAFSym = BAFnext;
                                }
                            } else
                                if (nextLevel == NA) {
                                    level = previousLevel;
                                    if (hasBAF_)
                                        BAFSym = BAFprev;
                                } else {
                                    level = nextLevel;
                                     if (hasBAF_)
                                        BAFSym = BAFnext;
                                }
                        }
                    }

                    if (level == NA)
                        level = normalLevel; //should never happen
                }
			}
			int copyNumber = round_f(level*ploidy_);
			//if (copyNumber < 0) //should happen only with "NOCALL"
           //     copyNumber = NA;
			if (copyNumberProbs_.count(copyNumber) == 0) {
				copyNumberProbs_.insert ( pair<int,double>(copyNumber,0) );
			}
			if (it->getEndsSize()==0) {
				copyNumberProbs_.find(copyNumber)->second += fragmentLength;
				if (copyNumber > NA) {
                    estimationOfGenomeSize_ += copyNumber*fragmentLength;
                    count+=fragmentLength;
				}

				if (cnumber == NA) {
					end = fragmentLength-1;
					cnumber = round_f(level*ploidy_);
					lBAF = BAFSym;
                    lUncertainty = BAFUncertainty ;
                    //if (cnumber < 0) //should happen only with "NOCALL"
                    //    cnumber = NA;
				} else {
					if (round_f(level*ploidy_) != cnumber || lBAF.compare(BAFSym)!=0 || lUncertainty != BAFUncertainty) {
					    int realEndOfTheCNV=(end+1)*windowSize_; //check that CNV is not larger than chr size
					    if (realEndOfTheCNV > it->getChrLength())
                            realEndOfTheCNV=it->getChrLength();

                        if (hasBAF_ && lBAF!="" && lBAF.compare("-")!=0)
                            cnumber = lBAF.length();

                        if (sex_.compare("XY")==0 && (it->getChromosome().find("X")!=string::npos || it->getChromosome().find("Y")!=string::npos)) {
                            if (cnumber != normalXYploidy && cnumber != NA) //push previous entry
                                if (hasBAF_)
                                    CNVs_.push_back(EntryCNV(it->getChromosome(),start,end,start*windowSize_,realEndOfTheCNV,cnumber,lUncertainty, lBAF,hasBAF_)); //save previous CNV
                                else
                                    CNVs_.push_back(EntryCNV(it->getChromosome(),start,end,start*windowSize_,realEndOfTheCNV,cnumber)); //save previous CNV
                            else if (hasBAF_ && lBAF.compare(NormalBAF_XY)!=0 && cnumber == normalXYploidy && lBAF!=""&& lBAF.compare("-")!=0) //abnormal BAF
                                    CNVs_.push_back(EntryCNV(it->getChromosome(),start,end,start*windowSize_,realEndOfTheCNV,cnumber,lUncertainty, lBAF,hasBAF_)); //save previous CNV


						} else {

                            if (cnumber != ploidy_ && cnumber != NA) //push previous entry
                                if (hasBAF_)
                                    CNVs_.push_back(EntryCNV(it->getChromosome(),start,end,start*windowSize_,realEndOfTheCNV,cnumber,lUncertainty, lBAF,hasBAF_)); //save previous CNV
                                else
                                    CNVs_.push_back(EntryCNV(it->getChromosome(),start,end,start*windowSize_,realEndOfTheCNV,cnumber)); //save previous CNV
                            else if (hasBAF_ && lBAF.compare(NormalBAF)!=0 && cnumber == ploidy_ && lBAF!=""&& lBAF.compare("-")!=0) //abnormal BAF
                                    CNVs_.push_back(EntryCNV(it->getChromosome(),start,end,start*windowSize_,realEndOfTheCNV,cnumber,lUncertainty, lBAF,hasBAF_)); //save previous CNV
						}

						//fill corresponding smoothed profile:

						for (int j = start; j<= end; j++)
							it->pushSmoothedProfile(float(cnumber)/ploidy_);

						start = end+1;
						end = start+fragmentLength-1;
						cnumber = round_f(level*ploidy_);
                        lBAF = BAFSym;
                        lUncertainty = BAFUncertainty ;
                        //if (cnumber < 0) //should happen only with "NOCALL"
                        //    cnumber = NA;
					} else {
						end = end+fragmentLength;  //should almost never happen if  (hasBAF_)
					}
				}
			} else {
				if (cnumber == NA) {
					end = fragmentLength-1; //should only once, in the beginning
					cnumber = round_f(level*ploidy_);
                    lBAF = BAFSym;
                    lUncertainty = BAFUncertainty ;
					//if (cnumber < 0) //should happen only with "NOCALL"
                    //        cnumber = NA;
				} else {
					if (round_f(level*ploidy_) != cnumber || lBAF.compare(BAFSym)!=0 || lUncertainty != BAFUncertainty) { //save previous value:
						int realLength = it->getEndAtBin(end)-it->getCoordinateAtBin(start)+1;
						copyNumberProbs_.find(cnumber)->second += realLength;
						if (cnumber>=0) {
                            estimationOfGenomeSize_ += cnumber*realLength;
                            count+=realLength;
                        }

                        int realEndOfTheCNV=it->getEndAtBin(end); //check that CNV is not larger than chr size
					    if (realEndOfTheCNV > it->getChrLength())
                            realEndOfTheCNV=it->getChrLength();


                        if (hasBAF_ && lBAF!="" && lBAF.compare("-")!=0)
                            cnumber = lBAF.length();


                        if (sex_.compare("XY")==0 && (it->getChromosome().find("X")!=string::npos || it->getChromosome().find("Y")!=string::npos)) {

                           if (cnumber != normalXYploidy && cnumber >= 0) //push previous entry
                                if (hasBAF_)
                                    CNVs_.push_back(EntryCNV(it->getChromosome(),start,end,it->getCoordinateAtBin(start),realEndOfTheCNV,cnumber,lUncertainty, lBAF,hasBAF_)); //save previous CNV
                                else
                                    CNVs_.push_back(EntryCNV(it->getChromosome(),start,end,it->getCoordinateAtBin(start),realEndOfTheCNV,cnumber)); //save previous CNV
                            else if (hasBAF_ && lBAF.compare(NormalBAF_XY)!=0 && cnumber == normalXYploidy && lBAF!=""&& lBAF.compare("-")!=0) //abnormal BAF
                                    CNVs_.push_back(EntryCNV(it->getChromosome(),start,end,it->getCoordinateAtBin(start),realEndOfTheCNV,cnumber,lUncertainty, lBAF,hasBAF_)); //save previous CNV


						} else {

                            if (cnumber != ploidy_ && cnumber >= 0) //push previous entry
                                if (hasBAF_)
                                    CNVs_.push_back(EntryCNV(it->getChromosome(),start,end,it->getCoordinateAtBin(start),realEndOfTheCNV,cnumber,lUncertainty, lBAF,hasBAF_)); //save previous CNV
                                else
                                    CNVs_.push_back(EntryCNV(it->getChromosome(),start,end,it->getCoordinateAtBin(start),realEndOfTheCNV,cnumber)); //save previous CNV
                            else if (hasBAF_ && lBAF.compare(NormalBAF)!=0 && cnumber == ploidy_ && lBAF!=""&& lBAF.compare("-")!=0) //abnormal BAF
                                    CNVs_.push_back(EntryCNV(it->getChromosome(),start,end,it->getCoordinateAtBin(start),realEndOfTheCNV,cnumber,lUncertainty, lBAF,hasBAF_)); //save previous CNV
                        }



						//fill corresponding smoothed profile:
						for (int j = start; j<= end; j++)
							it->pushSmoothedProfile(float(cnumber)/ploidy_);

						start = end+1;
						end = start+fragmentLength-1;
						cnumber = round_f(level*ploidy_);
                        lBAF = BAFSym;
                        lUncertainty = BAFUncertainty ;
						//if (cnumber < 0) //should happen only with "NOCALL"
                        //    cnumber = NA;
					} else {
						end = end+fragmentLength;
					}
				}

			}
		}
		//save the last CNV for this chromosomerealLength

        if (hasBAF_ && lBAF!="" && lBAF.compare("-")!=0)
            cnumber = lBAF.length();

		if (it->getEndsSize()==0) {
            int realEndOfTheCNV=(end+1)*windowSize_; //check that CNV is not larger than chr size
            if (realEndOfTheCNV > it->getChrLength())
                realEndOfTheCNV=it->getChrLength();

            if (sex_.compare("XY")==0 && (it->getChromosome().find("X")!=string::npos || it->getChromosome().find("Y")!=string::npos)) {
                if ((cnumber != normalXYploidy)&&(cnumber != NA))
                    CNVs_.push_back(EntryCNV(it->getChromosome(),start,end,start*windowSize_,realEndOfTheCNV,cnumber)); //save previous CNV

                if ((cnumber != normalXYploidy)&&(cnumber != NA))
                    if (hasBAF_)
                        CNVs_.push_back(EntryCNV(it->getChromosome(),start,end,start*windowSize_,realEndOfTheCNV,cnumber,lUncertainty, lBAF,hasBAF_)); //save previous CNV
                    else
                        CNVs_.push_back(EntryCNV(it->getChromosome(),start,end,start*windowSize_,realEndOfTheCNV,cnumber)); //save previous CNV
                else if (hasBAF_ && lBAF.compare(NormalBAF_XY)!=0 && cnumber == normalXYploidy && lBAF!=""&& lBAF.compare("-")!=0) //abnormal BAF
                    CNVs_.push_back(EntryCNV(it->getChromosome(),start,end,start*windowSize_,realEndOfTheCNV,cnumber,lUncertainty, lBAF,hasBAF_)); //save previous CNV

            } else {

                if ((cnumber != ploidy_)&&(cnumber != NA))
                    if (hasBAF_)
                        CNVs_.push_back(EntryCNV(it->getChromosome(),start,end,start*windowSize_,realEndOfTheCNV,cnumber,lUncertainty, lBAF,hasBAF_)); //save previous CNV
                    else
                        CNVs_.push_back(EntryCNV(it->getChromosome(),start,end,start*windowSize_,realEndOfTheCNV,cnumber)); //save previous CNV
                else if (hasBAF_ && lBAF.compare(NormalBAF)!=0 && cnumber == ploidy_ && lBAF!=""&& lBAF.compare("-")!=0) //abnormal BAF
                    CNVs_.push_back(EntryCNV(it->getChromosome(),start,end,start*windowSize_,realEndOfTheCNV,cnumber,lUncertainty, lBAF,hasBAF_)); //save previous CNV

            }


			for (int j = start; j<= end; j++)
				it->pushSmoothedProfile(float(cnumber)/ploidy_);
		} else {
		    endsSize = it->getEndsSize();
			int realLength = it->getEndAtBin(end)-it->getCoordinateAtBin(start)+1;
			copyNumberProbs_.find(cnumber)->second += realLength;
			if (cnumber>0) {
                estimationOfGenomeSize_ += cnumber*realLength;
                count+=realLength;
            }
            int realEndOfTheCNV=it->getEndAtBin(end); //check that CNV is not larger than chr size
            if (realEndOfTheCNV > it->getChrLength())
                    realEndOfTheCNV=it->getChrLength();

            if (sex_.compare("XY")==0 && (it->getChromosome().find("X")!=string::npos || it->getChromosome().find("Y")!=string::npos)) {

                if ((cnumber != normalXYploidy)&&(cnumber != NA))
                    if (hasBAF_)
                        CNVs_.push_back(EntryCNV(it->getChromosome(),start,end,it->getCoordinateAtBin(start),realEndOfTheCNV,cnumber,lUncertainty, lBAF,hasBAF_)); //save previous CNV
                    else
                        CNVs_.push_back(EntryCNV(it->getChromosome(),start,end,it->getCoordinateAtBin(start),realEndOfTheCNV,cnumber)); //save previous CNV
                else if (hasBAF_ && lBAF.compare(NormalBAF_XY)!=0 && cnumber == normalXYploidy && lBAF!=""&& lBAF.compare("-")!=0) //abnormal BAF
                    CNVs_.push_back(EntryCNV(it->getChromosome(),start,end,it->getCoordinateAtBin(start),realEndOfTheCNV,cnumber,lUncertainty, lBAF,hasBAF_)); //save previous CNV

            } else {
                if ((cnumber != ploidy_)&&(cnumber != NA))
                    if (hasBAF_)
                        CNVs_.push_back(EntryCNV(it->getChromosome(),start,end,it->getCoordinateAtBin(start),realEndOfTheCNV,cnumber,lUncertainty, lBAF,hasBAF_)); //save previous CNV
                    else
                        CNVs_.push_back(EntryCNV(it->getChromosome(),start,end,it->getCoordinateAtBin(start),realEndOfTheCNV,cnumber)); //save previous CNV
                else if (hasBAF_ && lBAF.compare(NormalBAF)!=0 && cnumber == ploidy_ && lBAF!=""&& lBAF.compare("-")!=0) //abnormal BAF
                    CNVs_.push_back(EntryCNV(it->getChromosome(),start,end,it->getCoordinateAtBin(start),realEndOfTheCNV,cnumber,lUncertainty, lBAF,hasBAF_)); //save previous CNV

            }


			for (int j = start; j<= end; j++)
				it->pushSmoothedProfile(float(cnumber)/ploidy_);

		}

	}
	if (endsSize==0) {
		estimationOfGenomeSize_ *= windowSize_;
	} else {
		//do nothing
	}
	map<int,double>::iterator it2;
	for ( it2=copyNumberProbs_.begin() ; it2 != copyNumberProbs_.end(); it2++ ) {
		if ((*it2).first>0) {
            cout << (*it2).first << "\t"<< (*it2).second/count<< "\n";
		}

	}
}

void GenomeCopyNumber::calculateCopyNumberProbs_and_genomeLength(int breakPointType) {
	estimationOfGenomeSize_ = 0;

	vector<ChrCopyNumber>::iterator it;
	unsigned long long count = 0;
    int min_fragment= telo_centromeric_flanks_/step_;
	int endsSize=0;
	string NormalBAF,NormalBAF_XY;
    float normalXYploidy=ploidy_*0.5; //will only use it when the genome is male

    if (hasBAF_) {
            NormalBAF = getNormalBAFforPloidy(ploidy_);
            NormalBAF_XY=getXYBAFforPloidy(ploidy_);//will only use it when the genome is male
    }

    CNVs_.clear();
    copyNumberProbs_.clear();

    //check it they are empty:
    cout << "copyNumberProbs_ contains:\n";
    map<int, double>::iterator itProb;
    for ( itProb=copyNumberProbs_.begin() ; itProb != copyNumberProbs_.end(); itProb++ )
    cout << (*itProb).first << " => " << (*itProb).second << "\n";

	for ( it=chrCopyNumber_.begin() ; it != chrCopyNumber_.end(); it++ ) {
		float previousLevel = NA;
		int lengthOfPreviousLevel = 0;
		float nextLevel = NA;
		int lengthOfNextLevel = 0;
		it->setIsSmoothed(true);
		string chr = it->getChromosome();
		cout << "..Annotation of CNVs for "<< chr <<"\n";
		//will use these three variables to collect CNVs:
		int start = 0;
		int end = NA;
		int cnumber = NA;
		string BAFSym = "-";
        string BAFprev = "-";
        string BAFnext = "-";
        string lBAF = "-";

        float lUncertainty = NA;
        float BAFUncertainty = NA;

        float normalLevel=1;
        if (sex_.compare("XY")==0 && (it->getChromosome().find("X")!=string::npos || it->getChromosome().find("Y")!=string::npos)) {
            normalLevel = 0.5;
        }

		for (int i = 0; i<(int)it->getMedianValues().size(); i++) {
			float level = it->getLevelAt(i, ploidy_);
			int fragmentLength =  it->getFragmentLengths()[i];
            if (hasBAF_) {
                BAFSym = it->getBAFsymbPerFrg(i);
                BAFUncertainty = it->getEstimatedBAFuncertaintyAtBin(i);
//                cout << "..Control: read "<< BAFSym << "\n";
            }

			if (breakPointType==HALFLENGTH && i<(int)it->getMedianValues().size()-1) {
                min_fragment = (it->getFragmentLengths()[i+1]+1)/2;
			}
			if (level != NA) {
				if (fragmentLength>min_fragment) {
					previousLevel = level;
					lengthOfPreviousLevel = fragmentLength;
					BAFprev = BAFSym;
				}
			} else {
                if (breakPointType==HALFLENGTH) {
                    min_fragment = (it->getFragmentLengths()[i]+1)/2;
                }
				int nextIndex = it->nextNoNAIndex(i, ploidy_,min_fragment);
				BAFUncertainty = NA;
				if (nextIndex == NA) {
					nextLevel = NA;
					lengthOfNextLevel = 0;
					BAFnext = "-";
				} else {
					nextLevel = it->getLevelAt(nextIndex, ploidy_);
					lengthOfNextLevel = it->getFragmentLengthsAt(nextIndex);
                    if (hasBAF_)
                        BAFnext = it->getBAFsymbPerFrg(nextIndex);
				}
				//nextLevel = round_by_ploidy(it->nextNoNAMedian(i, ploidy_), ploidy_);
				//lengthOfNextLevel = it->nextNoNALength(i, ploidy_);

                if (breakPointType==NOCALL) {
                    if ((previousLevel != NA) && (nextLevel != NA) && previousLevel==nextLevel && previousLevel==normalLevel) {
                        if (BAFprev.compare(BAFnext)==0) {
                            level = previousLevel;
                        }
                    }
                } else {
                    if (breakPointType==LARGECLOSE || breakPointType==SIMPLERIGHT || breakPointType==HALFLENGTH) {
                        if ((previousLevel != NA) && (nextLevel != NA) ) {
                            if (lengthOfPreviousLevel>lengthOfNextLevel) {
                                level = previousLevel;
                                if (hasBAF_) {
                                    BAFSym = BAFprev;
                                }
                            } else {
                                level = nextLevel;
                                if (hasBAF_)
                                    BAFSym = BAFnext;
                            }
                            //if (previousLevel == 1) //priority to "normal" copy number
                            //	level = 1;
                            //if (nextLevel == 1)
                            //	level = 1;
                        } else
                            if (nextLevel == NA) {
                                level = previousLevel;
                                if (hasBAF_)
                                    BAFSym = BAFprev;
                            } else {
                                level = nextLevel;
                                if (hasBAF_)
                                    BAFSym = BAFnext;
                            }
                    } else if (breakPointType==NORMALLEVEL) {
                        if ((previousLevel == normalLevel) || (nextLevel == normalLevel) ) {
                            level = normalLevel;
                            if (BAFprev.compare(BAFnext)==0)
                                BAFSym = BAFnext;
                            else
                                BAFSym = "-";

                        } else {
                            if ((previousLevel != NA) && (nextLevel != NA) ) {
                                if (lengthOfPreviousLevel>lengthOfNextLevel) {
                                    level = previousLevel;
                                    if (hasBAF_)
                                        BAFSym = BAFprev;
                                } else {
                                    level = nextLevel;
                                    if (hasBAF_)
                                        BAFSym = BAFnext;
                                }
                            } else
                                if (nextLevel == NA) {
                                    level = previousLevel;
                                    if (hasBAF_)
                                        BAFSym = BAFprev;
                                } else {
                                    level = nextLevel;
                                     if (hasBAF_)
                                        BAFSym = BAFnext;
                                }
                        }
                    }

                    if (level == NA)
                        level = normalLevel; //should never happen
                }
			}
			int copyNumber = round_f(level*ploidy_);
			//if (copyNumber < 0) //should happen only with "NOCALL"
            //    copyNumber = NA;
			if (copyNumberProbs_.count(copyNumber) == 0) {
				copyNumberProbs_.insert ( pair<int,double>(copyNumber,0) );
			}
			if (it->getEndsSize()==0) {
				copyNumberProbs_.find(copyNumber)->second += fragmentLength;
				if (copyNumber > NA) {
                    estimationOfGenomeSize_ += copyNumber*fragmentLength;
                    count+=fragmentLength;
				}

				if (cnumber == NA) {
					end = fragmentLength-1;
					cnumber = round_f(level*ploidy_);
					lBAF = BAFSym;
                    lUncertainty = BAFUncertainty ;
                    //if (cnumber < 0) //should happen only with "NOCALL"
                    //    cnumber = NA;
				} else {
					if (round_f(level*ploidy_) != cnumber || lBAF.compare(BAFSym)!=0 || lUncertainty != BAFUncertainty) {
					    int realEndOfTheCNV=(end+1)*windowSize_; //check that CNV is not larger than chr size
					    if (realEndOfTheCNV > it->getChrLength())
                            realEndOfTheCNV=it->getChrLength();

                        if (hasBAF_ && lBAF!="" && lBAF.compare("-")!=0)
                            cnumber = lBAF.length();

                        if (sex_.compare("XY")==0 && (it->getChromosome().find("X")!=string::npos || it->getChromosome().find("Y")!=string::npos)) {
                            if (cnumber != normalXYploidy && cnumber != NA) //push previous entry
                                if (hasBAF_)
                                    CNVs_.push_back(EntryCNV(it->getChromosome(),start,end,start*windowSize_,realEndOfTheCNV,cnumber,lUncertainty, lBAF,hasBAF_)); //save previous CNV
                                else
                                    CNVs_.push_back(EntryCNV(it->getChromosome(),start,end,start*windowSize_,realEndOfTheCNV,cnumber)); //save previous CNV
                            else if (hasBAF_ && lBAF.compare(NormalBAF_XY)!=0 && cnumber == normalXYploidy && lBAF!=""&& lBAF.compare("-")!=0) //abnormal BAF
                                    CNVs_.push_back(EntryCNV(it->getChromosome(),start,end,start*windowSize_,realEndOfTheCNV,cnumber,lUncertainty, lBAF,hasBAF_)); //save previous CNV


						} else {

                            if (cnumber != ploidy_ && cnumber != NA) //push previous entry
                                if (hasBAF_)
                                    CNVs_.push_back(EntryCNV(it->getChromosome(),start,end,start*windowSize_,realEndOfTheCNV,cnumber,lUncertainty, lBAF,hasBAF_)); //save previous CNV
                                else
                                    CNVs_.push_back(EntryCNV(it->getChromosome(),start,end,start*windowSize_,realEndOfTheCNV,cnumber)); //save previous CNV
                            else if (hasBAF_ && lBAF.compare(NormalBAF)!=0 && cnumber == ploidy_ && lBAF!=""&& lBAF.compare("-")!=0) //abnormal BAF
                                    CNVs_.push_back(EntryCNV(it->getChromosome(),start,end,start*windowSize_,realEndOfTheCNV,cnumber,lUncertainty, lBAF,hasBAF_)); //save previous CNV
						}

						//fill corresponding smoothed profile:

						for (int j = start; j<= end; j++)
							it->pushSmoothedProfile(float(cnumber)/ploidy_);

						start = end+1;
						end = start+fragmentLength-1;
						cnumber = round_f(level*ploidy_);
                        lBAF = BAFSym;
                        lUncertainty = BAFUncertainty ;
                        //if (cnumber < 0) //should happen only with "NOCALL"
                        //    cnumber = NA;
					} else {
						end = end+fragmentLength;  //should almost never happen if  (hasBAF_)
					}
				}
			} else {
				if (cnumber == NA) {
					end = fragmentLength-1; //should only once, in the beginning
					cnumber = round_f(level*ploidy_);
                    lBAF = BAFSym;
                    lUncertainty = BAFUncertainty ;
					//if (cnumber < 0) //should happen only with "NOCALL"
                    //        cnumber = NA;
				} else {
					if (round_f(level*ploidy_) != cnumber || lBAF.compare(BAFSym)!=0 || lUncertainty != BAFUncertainty) { //save previous value:
						int realLength = it->getEndAtBin(end)-it->getCoordinateAtBin(start)+1;
						copyNumberProbs_.find(cnumber)->second += realLength;
						if (cnumber>=0) {
                            estimationOfGenomeSize_ += cnumber*realLength;
                            count+=realLength;
                        }

                        int realEndOfTheCNV=it->getEndAtBin(end); //check that CNV is not larger than chr size
					    if (realEndOfTheCNV > it->getChrLength())
                            realEndOfTheCNV=it->getChrLength();


                        if (hasBAF_ && lBAF!="" && lBAF.compare("-")!=0)
                            cnumber = lBAF.length();


                        if (sex_.compare("XY")==0 && (it->getChromosome().find("X")!=string::npos || it->getChromosome().find("Y")!=string::npos)) {

                           if (cnumber != normalXYploidy && cnumber >= 0) //push previous entry
                                if (hasBAF_)
                                    CNVs_.push_back(EntryCNV(it->getChromosome(),start,end,it->getCoordinateAtBin(start),realEndOfTheCNV,cnumber,lUncertainty, lBAF,hasBAF_)); //save previous CNV
                                else
                                    CNVs_.push_back(EntryCNV(it->getChromosome(),start,end,it->getCoordinateAtBin(start),realEndOfTheCNV,cnumber)); //save previous CNV
                            else if (hasBAF_ && lBAF.compare(NormalBAF_XY)!=0 && cnumber == normalXYploidy && lBAF!=""&& lBAF.compare("-")!=0) //abnormal BAF
                                    CNVs_.push_back(EntryCNV(it->getChromosome(),start,end,it->getCoordinateAtBin(start),realEndOfTheCNV,cnumber,lUncertainty, lBAF,hasBAF_)); //save previous CNV


						} else {

                            if (cnumber != ploidy_ && cnumber >= 0) //push previous entry
                                if (hasBAF_)
                                    CNVs_.push_back(EntryCNV(it->getChromosome(),start,end,it->getCoordinateAtBin(start),realEndOfTheCNV,cnumber,lUncertainty, lBAF,hasBAF_)); //save previous CNV
                                else
                                    CNVs_.push_back(EntryCNV(it->getChromosome(),start,end,it->getCoordinateAtBin(start),realEndOfTheCNV,cnumber)); //save previous CNV
                            else if (hasBAF_ && lBAF.compare(NormalBAF)!=0 && cnumber == ploidy_ && lBAF!=""&& lBAF.compare("-")!=0) //abnormal BAF
                                    CNVs_.push_back(EntryCNV(it->getChromosome(),start,end,it->getCoordinateAtBin(start),realEndOfTheCNV,cnumber,lUncertainty, lBAF,hasBAF_)); //save previous CNV
                        }



						//fill corresponding smoothed profile:
						for (int j = start; j<= end; j++)
							it->pushSmoothedProfile(float(cnumber)/ploidy_);

						start = end+1;
						end = start+fragmentLength-1;
						cnumber = round_f(level*ploidy_);
                        lBAF = BAFSym;
                        lUncertainty = BAFUncertainty ;
						//if (cnumber < 0) //should happen only with "NOCALL"
                        //    cnumber = NA;
					} else {
						end = end+fragmentLength;
					}
				}

			}
		}
		//save the last CNV for this chromosomerealLength

        if (hasBAF_ && lBAF!="" && lBAF.compare("-")!=0)
            cnumber = lBAF.length();

		if (it->getEndsSize()==0) {
            int realEndOfTheCNV=(end+1)*windowSize_; //check that CNV is not larger than chr size
            if (realEndOfTheCNV > it->getChrLength())
                realEndOfTheCNV=it->getChrLength();

            if (sex_.compare("XY")==0 && (it->getChromosome().find("X")!=string::npos || it->getChromosome().find("Y")!=string::npos)) {
                if ((cnumber != normalXYploidy)&&(cnumber != NA))
                    CNVs_.push_back(EntryCNV(it->getChromosome(),start,end,start*windowSize_,realEndOfTheCNV,cnumber)); //save previous CNV

                if ((cnumber != normalXYploidy)&&(cnumber != NA))
                    if (hasBAF_)
                        CNVs_.push_back(EntryCNV(it->getChromosome(),start,end,start*windowSize_,realEndOfTheCNV,cnumber,lUncertainty, lBAF,hasBAF_)); //save previous CNV
                    else
                        CNVs_.push_back(EntryCNV(it->getChromosome(),start,end,start*windowSize_,realEndOfTheCNV,cnumber)); //save previous CNV
                else if (hasBAF_ && lBAF.compare(NormalBAF_XY)!=0 && cnumber == normalXYploidy && lBAF!=""&& lBAF.compare("-")!=0) //abnormal BAF
                    CNVs_.push_back(EntryCNV(it->getChromosome(),start,end,start*windowSize_,realEndOfTheCNV,cnumber,lUncertainty, lBAF,hasBAF_)); //save previous CNV

            } else {

                if ((cnumber != ploidy_)&&(cnumber != NA))
                    if (hasBAF_)
                        CNVs_.push_back(EntryCNV(it->getChromosome(),start,end,start*windowSize_,realEndOfTheCNV,cnumber,lUncertainty, lBAF,hasBAF_)); //save previous CNV
                    else
                        CNVs_.push_back(EntryCNV(it->getChromosome(),start,end,start*windowSize_,realEndOfTheCNV,cnumber)); //save previous CNV
                else if (hasBAF_ && lBAF.compare(NormalBAF)!=0 && cnumber == ploidy_ && lBAF!=""&& lBAF.compare("-")!=0) //abnormal BAF
                    CNVs_.push_back(EntryCNV(it->getChromosome(),start,end,start*windowSize_,realEndOfTheCNV,cnumber,lUncertainty, lBAF,hasBAF_)); //save previous CNV

            }


			for (int j = start; j<= end; j++)
				it->pushSmoothedProfile(float(cnumber)/ploidy_);
		} else {
		    endsSize = it->getEndsSize();
			int realLength = it->getEndAtBin(end)-it->getCoordinateAtBin(start)+1;
			copyNumberProbs_.find(cnumber)->second += realLength;
			if (cnumber>0) {
                estimationOfGenomeSize_ += cnumber*realLength;
                count+=realLength;
            }
            int realEndOfTheCNV=it->getEndAtBin(end); //check that CNV is not larger than chr size
            if (realEndOfTheCNV > it->getChrLength())
                    realEndOfTheCNV=it->getChrLength();

            if (sex_.compare("XY")==0 && (it->getChromosome().find("X")!=string::npos || it->getChromosome().find("Y")!=string::npos)) {

                if ((cnumber != normalXYploidy)&&(cnumber != NA))
                    if (hasBAF_)
                        CNVs_.push_back(EntryCNV(it->getChromosome(),start,end,it->getCoordinateAtBin(start),realEndOfTheCNV,cnumber,lUncertainty, lBAF,hasBAF_)); //save previous CNV
                    else
                        CNVs_.push_back(EntryCNV(it->getChromosome(),start,end,it->getCoordinateAtBin(start),realEndOfTheCNV,cnumber)); //save previous CNV
                else if (hasBAF_ && lBAF.compare(NormalBAF_XY)!=0 && cnumber == normalXYploidy && lBAF!=""&& lBAF.compare("-")!=0) //abnormal BAF
                    CNVs_.push_back(EntryCNV(it->getChromosome(),start,end,it->getCoordinateAtBin(start),realEndOfTheCNV,cnumber,lUncertainty, lBAF,hasBAF_)); //save previous CNV

            } else {
                if ((cnumber != ploidy_)&&(cnumber != NA))
                    if (hasBAF_)
                        CNVs_.push_back(EntryCNV(it->getChromosome(),start,end,it->getCoordinateAtBin(start),realEndOfTheCNV,cnumber,lUncertainty, lBAF,hasBAF_)); //save previous CNV
                    else
                        CNVs_.push_back(EntryCNV(it->getChromosome(),start,end,it->getCoordinateAtBin(start),realEndOfTheCNV,cnumber)); //save previous CNV
                else if (hasBAF_ && lBAF.compare(NormalBAF)!=0 && cnumber == ploidy_ && lBAF!=""&& lBAF.compare("-")!=0) //abnormal BAF
                    CNVs_.push_back(EntryCNV(it->getChromosome(),start,end,it->getCoordinateAtBin(start),realEndOfTheCNV,cnumber,lUncertainty, lBAF,hasBAF_)); //save previous CNV

            }


			for (int j = start; j<= end; j++)
				it->pushSmoothedProfile(float(cnumber)/ploidy_);

		}

	}
	if (endsSize==0) {
		estimationOfGenomeSize_ *= windowSize_;
	} else {
		//do nothing
	}
	map<int,double>::iterator it2;
	for ( it2=copyNumberProbs_.begin() ; it2 != copyNumberProbs_.end(); it2++ ) {
		if ((*it2).first>0) {
            cout << (*it2).first << "\t"<< (*it2).second/count<< "\n";
		}

	}

}

void GenomeCopyNumber::readGemMappabilityFile(std::string const& inFile) {
	ifstream file (inFile.c_str());
	string line;
    string::size_type pos = 0;
    string currentChr = "";
    int index = NA;
    int count = 0;
    int uniqueCount = 0;
    int localWindowSize = windowSize_;
    float ratio;
    int positionInd = 0;
    int startPos = 0;
    int endPos = 0;
    int lastEnd = -1;
    string text = "";
	if (file.is_open())	{
	    cout << "..Reading "<< inFile << "\n";
	    int countChromosomes = 0;
	    int countChromosomesOutOfIndex = 0;
	    isMappUsed_=true;
		while (! file.eof() )	{
			getline (file,line);
			if (! line.length()) continue;
			pos = 0;
			if ( (line.substr(0,1) == "~") && (line.substr(0,2) != "~~")){
			//if ( ( pos = line.find("~chr", pos)) != string::npos ){
                //save the last info:
                if (count > endPos && index != NA) { // endPos ?
                    int howMuchToDelete = startPos-lastEnd-1;
                    text = text.substr(howMuchToDelete);
                    uniqueCount =0;
                    for(int i = 0; i < (int)text.length(); ++i)
                        if (text[i] == '!')
                            uniqueCount++;

                    ratio = float(uniqueCount)/localWindowSize;
                    chrCopyNumber_[index].setMappabilityProfileAt(positionInd, ratio);
                }
                //restore all variables for a new chromosome
                if ( ( pos = line.find("~chr", pos)) != string::npos ){
                    currentChr = line.substr(4);
                } else {
                    currentChr = line.substr(1);
                }
                cout << "..Reading mappability for chromosome " << currentChr << "\n";
                index = findIndex(currentChr);

                positionInd = 0;
				count = 0;
				uniqueCount = 0;
                text = "";
                lastEnd = -1;

                //cout <<  "..Index for chromosome " << currentChr << ": "<< index << "\n";

				if (index == NA) {
				    cout <<  "skipping chromosome " << currentChr << "\n";
				    countChromosomesOutOfIndex ++;
				    //return;
				    // do not return! they can be other "good" chromosomes afterwords!
				    //do nothing!!! wait for the next chromosome!

                } else {
                    countChromosomes++;
                    startPos = chrCopyNumber_[index].getCoordinateAtBin(positionInd);
                    endPos = chrCopyNumber_[index].getEndAtBin(positionInd);
                    localWindowSize = endPos-startPos+1;
                    int maxInd = chrCopyNumber_[index].getLength()-1;
                    cout << "..Control: Last window: from " << chrCopyNumber_[index].getCoordinateAtBin(maxInd) << " to " << chrCopyNumber_[index].getEndAtBin(maxInd) <<"\n";
                    chrCopyNumber_[index].createMappabilityProfile();
                }

			} else if (index != NA) {
                count += line.length();
                text.append(line);

                if (count > endPos) {
                    int howMuchToDelete = startPos-lastEnd-1;
                    if (howMuchToDelete<0)
                        howMuchToDelete=0; //this should never happen
                    text = text.substr(howMuchToDelete);
                    lastEnd = endPos;
                    string substr_ = text.substr(0,localWindowSize);
                    uniqueCount = 0;
                    for(int i = 0; i < (int)substr_.length(); ++i)
                        if (substr_[i] == '!')
                            uniqueCount++;


                    if (positionInd+1<chrCopyNumber_[index].getLength()) {
                        int nextStart = chrCopyNumber_[index].getCoordinateAtBin(positionInd+1);
                        if (nextStart<=endPos) {
                            //count = end;
                            //and delete prefix in text;
                            int howMuchToDelete = nextStart - startPos;
                            text = text.substr(howMuchToDelete);
                        } else {
                            text = text.substr(localWindowSize); //actually this should be the same as "text.substr(howMuchToDelete);"
                        }

                    }

                    //count -= localWindowSize;

                    ratio = float(uniqueCount)/localWindowSize;
                    chrCopyNumber_[index].setMappabilityProfileAt(positionInd, ratio);
                    positionInd ++;

                    if (positionInd<chrCopyNumber_[index].getLength()) {
                        startPos = chrCopyNumber_[index].getCoordinateAtBin(positionInd);
                        endPos = chrCopyNumber_[index].getEndAtBin(positionInd);
                        localWindowSize = endPos-startPos+1;
                    } else {
                        //should not read this chromosome any more
                        index = NA;
                    }
                }
			}
		}
		file.close();
		//save the very last info (int text variable)
        if (count >= startPos && index != NA) { //endPos
            uniqueCount =0;
               for(int i = 0; i < (int)text.length(); ++i)
                   if (text[i] == '!')
                        uniqueCount++;
            ratio = float(uniqueCount)/localWindowSize; //  /max((int)text.length(),(int)localWindowSize);
            chrCopyNumber_[index].setMappabilityProfileAt(positionInd, ratio);
        }

		cout << "..file " << inFile << " is read\n";
		cout << "..Mappability profile was taken into account for "<< countChromosomes <<" chromosomes\n";
		cout << ".."<< countChromosomesOutOfIndex <<" chromosomes were skipped\n";

//		if (chrCopyNumber_[1].getLength() == chrCopyNumber_[1].getMappabilityLength())
//            cout << "Mappability profile has been set correctly";
	} else{
        cerr << "Error: Unable to open file "+inFile+"\n";
        exit(-1);
	}
}

int GenomeCopyNumber::readCGprofile(std::string const& inFile) {
	ifstream file (inFile.c_str());
	string line;
	int count = 0;
	int observedStep = 0;
	if (file.is_open())	{
		while (! file.eof() )	{
			getline (file,line);
			if (! line.length()) continue;
			std::vector<std::string> strs = split(line, '\t');
			if (strs.size()>=4) {
				string currentChr = strs[0];
				if (observedStep==0)
                    observedStep = ceil(strtod(strs[1].c_str(), NULL));
				float CGperc =(float)strtod(strs[2].c_str(), NULL);
				float nonNperc =(float)strtod(strs[3].c_str(), NULL);
				string::size_type pos = 0;
				if ( ( pos = currentChr.find("chr", pos)) != string::npos )
					currentChr.replace( pos, 3, "" );
				int index = findIndex(currentChr);
				if (index != NA) {
				    chrCopyNumber_[index].addToCGcontent(CGperc);
                    chrCopyNumber_[index].addToNonNpercent(nonNperc);
                    count++;
                    if (strs.size()==5) { //means that there are also mappability values in the 5th colomn
                        float MappPerc =(float)strtod(strs[4].c_str(), NULL);
                        chrCopyNumber_[index].addToMappabilityProfile(MappPerc);
                    }
                    if (observedStep!=0 && observedStep!=step_ && step_!=NA && step_!=0) {
                        file.close();
                        chrCopyNumber_[index].clearCGcontent();
                        chrCopyNumber_[index].clearNonNpercent();
                        chrCopyNumber_[index].clearMappabilityProfile();
                        return observedStep;

                    }
				}
			}

			strs.clear();
		}
		file.close();
		cout << "file " << inFile << " is read\n";
		if (count==0){
            cerr << "Your GC-content file "<<inFile<< " is empty or is in a wrong format\n\nPlease use chomosome sequences (option \"chrFiles\") to recreate it!\n\n";
            exit(-1);
		}
	} else {
	    cerr << "Unable to open file "+inFile+"\n";
	    exit (-1);
    }

    return observedStep;
}

void GenomeCopyNumber::readCopyNumber(std::string const& inFile) {
  if (WESanalysis == false)
    {
	totalNumberOfPairs_ = 0;
	normalNumberOfPairs_ = 0;
	windowSize_ = 0;
	ifstream file (inFile.c_str());
	string line;
	string currentChr = "";
	int chrCount = -1;
	refGenomeSize_ = 0;
	string::size_type pos = 0;
	int count = 0;
	if (file.is_open())	{
		while (! file.eof() )	{
			getline (file,line);

			if (! line.length()) continue;

			std::vector<std::string> strs = split(line, '\t');
			if (strs.size()==3) {
				if (strs[0].compare(currentChr)!=0) {
					chrCount ++;
					currentChr = strs[0];
					chromosomesInd_.insert(pair<string, int> (currentChr,chrCount));
					chrCopyNumber_.push_back(ChrCopyNumber(currentChr));
					windowSize_ = atoi(strs[1].c_str());
				}
				float cp = (float)strtod(strs[2].c_str(), NULL);
				chrCopyNumber_[chrCount].addToReadCount(cp);
				chrCopyNumber_[chrCount].addToCoordinates(atoi(strs[1].c_str()));
				normalNumberOfPairs_ += (int)cp;
				if (windowSize_ == 0)
					windowSize_ = atoi(strs[1].c_str());
                count++;
			}
			else if (strs.size()==4) {
				string chrNumber = strs[0];
				if ( ( pos = chrNumber.find("chr", pos)) != string::npos )
					chrNumber.replace( pos, 3, "" );

				if (chrNumber.compare(currentChr)!=0) {
					chrCount ++;
					currentChr = chrNumber;
					chromosomesInd_.insert(pair<string, int> (currentChr,chrCount));
					chrCopyNumber_.push_back(ChrCopyNumber(currentChr));
					windowSize_ = atoi(strs[2].c_str())-atoi(strs[1].c_str())+1;
				}
				float cp = (float)strtod(strs[3].c_str(), NULL);
				chrCopyNumber_[chrCount].addToReadCount(cp);
				chrCopyNumber_[chrCount].addToCoordinates(atoi(strs[1].c_str()));
				chrCopyNumber_[chrCount].addToEnds(atoi(strs[2].c_str()));
				normalNumberOfPairs_ += (int)cp;
				if (count>0 && step_==NA) {
				  // int first = chrCopyNumber_[chrCount].getCoordinateAtBin(0);
				//	int second = chrCopyNumber_[chrCount].getCoordinateAtBin(1);
					step_ = chrCopyNumber_[chrCount].getCoordinateAtBin(1)-chrCopyNumber_[chrCount].getCoordinateAtBin(0);
				}

				if (windowSize_ == 0)
					windowSize_ = atoi(strs[1].c_str());

                count++;
			}

			strs.clear();
		}
		file.close();
		cout << "file " << inFile << " read\n";
		if (step_==NA)
            step_=windowSize_;
		//fill other values
		map<string,int>::iterator it;
		for ( it=chromosomesInd_.begin() ; it != chromosomesInd_.end(); it++ ) {
			chrCopyNumber_[(*it).second].setWindowSize(windowSize_);
			int length = chrCopyNumber_[(*it).second].getValues().size();
			chrCopyNumber_[(*it).second].setVectorLength(length);
			int chromosomeLength=length*step_;
			if (chrCopyNumber_[(*it).second].getEndsSize()>0)
                chromosomeLength=chrCopyNumber_[(*it).second].getEndAtBin(chrCopyNumber_[(*it).second].getEndsSize()-1);

            chrCopyNumber_[(*it).second].setChrLength(chromosomeLength);

			chrCopyNumber_[(*it).second].setStep(step_);
			refGenomeSize_ += chromosomeLength;
			//cout <<(*it).first <<"\t" <<(*it).second << "\t" <<windowSize_ << "\t" << length << "\t" <<length*windowSize_ << "\n";
		}

	} else {
	    cerr << "Unable to open file "+inFile+"\n";
	    exit(-1);
	    //throw "Unable to open file "+inFile+"\n";
    }
	totalNumberOfPairs_ = normalNumberOfPairs_;
	if (step_==NA)
            step_=windowSize_;
	if (WESanalysis == false)   {cout << "\t evaluated genome size:\t" << refGenomeSize_ << "\n";}
} else {
    totalNumberOfPairs_ = 0;
	normalNumberOfPairs_ = 0;
	windowSize_ = 0;
	ifstream file (inFile.c_str());
	string line;
	string currentChr = "";
	int chrCount = -1;
	refGenomeSize_ = 0;
	string::size_type pos = 0;
	int count = 0;
	if (file.is_open())	{
		while (! file.eof() )	{
			getline (file,line);

			if (! line.length()) continue;

			std::vector<std::string> strs = split(line, '\t');
				string chrNumber = strs[0];
				if ( ( pos = chrNumber.find("chr", pos)) != string::npos )
					chrNumber.replace( pos, 3, "" );

				if (chrNumber.compare(currentChr)!=0) {
					chrCount ++;
					currentChr = chrNumber;
					chromosomesInd_.insert(pair<string, int> (currentChr,chrCount));
					chrCopyNumber_.push_back(ChrCopyNumber(currentChr));
					windowSize_ = 0;
				}
				float cp = (float)strtod(strs[3].c_str(), NULL);
				chrCopyNumber_[chrCount].addToReadCount(cp);
				chrCopyNumber_[chrCount].addToCoordinates(atoi(strs[1].c_str()));
				chrCopyNumber_[chrCount].addToEnds(atoi(strs[2].c_str()));
				if (strs.size()>=4)
                    chrCopyNumber_[chrCount].addToGenes_name(strs[4].c_str());
				normalNumberOfPairs_ += (int)cp;
				if (count>0 && step_==NA) {
					step_ = 0;
				}
                count++;

			strs.clear();
		}
		file.close();
		cout << "file " << inFile << " read\n";
		if (step_==NA)
            step_=windowSize_;
		//fill other values
		map<string,int>::iterator it;
		for ( it=chromosomesInd_.begin() ; it != chromosomesInd_.end(); it++ ) {
			chrCopyNumber_[(*it).second].setWindowSize(windowSize_);
			int length = chrCopyNumber_[(*it).second].getValues().size();
			chrCopyNumber_[(*it).second].setVectorLength(length);
			if (SeekingSubc_ == true)
                {
                chrCopyNumber_[(*it).second].setCN_subcLength(length+3);
                chrCopyNumber_[(*it).second].setpop_subcLength(length+3);
                }

            int chromosomeLength=length*step_;
			if (chrCopyNumber_[(*it).second].getEndsSize()>=0)
                chromosomeLength=chrCopyNumber_[(*it).second].getEndAtBin(chrCopyNumber_[(*it).second].getEndsSize()-1);

            chrCopyNumber_[(*it).second].setChrLength(chromosomeLength);

			chrCopyNumber_[(*it).second].setStep(0);
			refGenomeSize_ += chrCopyNumber_[(*it).second].getLength();
			//cout <<(*it).first <<"\t" <<(*it).second << "\t" <<windowSize_ << "\t" << length << "\t" <<length*windowSize_ << "\n";
		}

	} else {
	    cerr << "Unable to open file "+inFile+"\n";
	    exit(-1);
	    //throw "Unable to open file "+inFile+"\n";
    }
	totalNumberOfPairs_ = normalNumberOfPairs_;
	if (step_==NA)
            step_=windowSize_;
    if (WESanalysis == false)
       {cout << "\t evaluated genome size:\t" << refGenomeSize_ << "\n";}
	}
}


void GenomeCopyNumber::printCopyNumber(std::string const& outFile) {
	const char * name = outFile.c_str();
	std::ofstream file;
	file.open(name);
	map<string,int>::iterator it;
	for ( it=chromosomesInd_.begin() ; it != chromosomesInd_.end(); it++ ) {
		printCopyNumber((*it).first,file);
	}
	file.close();
	cout << "printing counts into "<<outFile <<"\n";
}

void GenomeCopyNumber::printCGprofile(std::string const& outFile) {
	const char * name = outFile.c_str();
	std::ofstream file;
	file.open(name);
	map<string,int>::iterator it;
	for ( it=chromosomesInd_.begin() ; it != chromosomesInd_.end(); it++ ) {
		printCGprofile((*it).first,file);
	}
	file.close();
	cout << "CG-content printed into "<<outFile <<"\n";
}

void GenomeCopyNumber::printRatio(std::string const& chr, std::string const& outFile, bool printNA) {
	cout << "printing ratio for "<<chr<<" into "<<outFile <<"\n";
	std::ofstream file (outFile.c_str());
	printRatio(chr,file,printNA);
	file.close();
}

void GenomeCopyNumber::printCopyNumber(std::string const& chr, std::string const& outFile) {
	std::ofstream file (outFile.c_str());
	printCopyNumber(chr,file);
	file.close();
}

void GenomeCopyNumber::printInfo(std::ofstream & file) {
    file << "Output_Ploidy\t" <<ploidy_<< endl;
    file << "Sample_Purity\t" <<1-normalContamination_<< endl;

}

void GenomeCopyNumber::printCNVs (std::string const& outFile) {
	std::ofstream file (outFile.c_str());
	for (int i = 0; i < (int)CNVs_.size(); i++) {
	    if (CNVs_[i].getCopyNumber()>-1) {
            if (sex_.compare("XY")==0 && (CNVs_[i].getChr().find("X")!=string::npos || CNVs_[i].getChr().find("Y")!=string::npos)) {
                file << CNVs_[i].printEntryCNV(ploidy_*0.5);
	        } else {
                file << CNVs_[i].printEntryCNV(ploidy_);
	        }

            file << "\n";
	    }
	}
	file.close();
}

vector <EntryCNV> GenomeCopyNumber::getCNVs () {
    return CNVs_;
}

int GenomeCopyNumber::getPloidy() {
    return ploidy_;
}


void GenomeCopyNumber::calculateSomaticCNVs (std::vector <EntryCNV> controlCNVs, int controlPloidy) {
    ifUsedControl_=true;

    int overlapPrecision = 3; //+-3 windows => the same CNA

	std::cout << "..Calculate somatic CNVs" << std::endl;
    for (int i = 0; i < (int)CNVs_.size(); i++) {
        string type = "somatic";
        vector <int>lefts;
        vector <int> rights;
        for (int unsigned j = 0; j < controlCNVs.size(); j++) {
            int left, right;
            float overlap;
            if (sex_.compare("XY")==0 && (CNVs_[i].getChr().find("X")!=string::npos || CNVs_[i].getChr().find("Y")!=string::npos)) {
                overlap=CNVs_[i].compare(controlCNVs[j], overlapPrecision, left, right,ploidy_*0.5,controlPloidy*0.5);
            } else {
                overlap=CNVs_[i].compare(controlCNVs[j], overlapPrecision, left, right,ploidy_,controlPloidy);
            }
            if (overlap==1) {
                type = "germline";
                CNVs_[i].setGermlinePercent(overlap);
                break;
            } else if (overlap>0) {
                lefts.push_back(left);
                rights.push_back(right);
            }
        }
        CNVs_[i].setType(type);
        if (type == "somatic" && lefts.size()>0) {
            //we will devide this CNV into peaces..
            int germlineLength = calculateTotalLength(lefts,rights);
            CNVs_[i].setGermlinePercent(germlineLength*1./(CNVs_[i].getEnd()-CNVs_[i].getStart()+1));
            lefts.clear();
            rights.clear();
        } else  if (type == "somatic") {
            CNVs_[i].setGermlinePercent(0);
        }
	}

}
void GenomeCopyNumber::printRatioBedGraph(std::string const& chr, std::ofstream & file, std::string const& typeCNA) {
    string::size_type pos = 0;
    float value;
    int position;
    string myType;
	string chrNumber = chr;
	if ( ( pos = chrNumber.find("chr", pos)) != string::npos )
		chrNumber.replace( pos, 3, "" );
	map<string,ChrCopyNumber>::iterator it;
	int index = findIndex(chrNumber);
	if (index == NA) {return;}
	if (!chrCopyNumber_[index].isMedianCalculated()) {cerr << "Error: unable to write a BedGraph with without caclulation of median values..\n"; return;}
	int length = chrCopyNumber_[index].getLength();

    if (typeCNA.compare("LOH")==0) {
        if (!hasBAF_) {cerr << "Error: unable to write LOH in a BedGraph with without caclulation of genotypes..\n"; return;}
    }

    for (int i = 0; i< length; i++) {
            value=chrCopyNumber_[index].getRatioAtBin(i);
            if (isRatioLogged_ && value!=NA)
                value=pow(2,value);
            position=chrCopyNumber_[index].getCoordinateAtBin(i);
            float valueToPrint;
			if (chrCopyNumber_[index].isSmoothed())
				valueToPrint = chrCopyNumber_[index].getSmoothedProfileAtI(i)*ploidy_;
			else
				valueToPrint = round_by_ploidy(chrCopyNumber_[index].getMedianProfileAtI(i),ploidy_)*ploidy_;

            if (valueToPrint<0)
                valueToPrint = NA;
            if (hasBAF_) {
                 if (valueToPrint>0 && chrCopyNumber_[index].getBAFsymbolAt(i).compare("-")!=0 &&chrCopyNumber_[index].getBAFsymbolAt(i)!="" && chrCopyNumber_[index].getBAFsymbolAt(i).length()!=valueToPrint ) { //&& (chrCopyNumber_[index].getEstimatedBAFuncertaintyAtI(i)< MAXUncertainty)
                    valueToPrint=chrCopyNumber_[index].getBAFsymbolAt(i).length();
                 }
            }
            myType="NA";
            if (valueToPrint>ploidy_ && valueToPrint != NA) {
                myType = "gain";
            } else if (valueToPrint<ploidy_&& valueToPrint != NA) {
                myType = "loss";
            } else if (valueToPrint==ploidy_) {
                myType = "normal";
            }
            if (hasBAF_ && valueToPrint==ploidy_ && chrCopyNumber_[index].getBAFProfileAt(i)==0) {
                myType="LOH";
            }
            if (myType.compare(typeCNA)==0 && value != NA)
                file << "chr" <<chrNumber << "\t"<<position<< "\t"<< position + windowSize_ << "\t"<<value*ploidy_<<"\n";
	}
}

void GenomeCopyNumber::printRatio(std::string const& chr, std::ofstream & file, bool printNA) {
	string::size_type pos = 0;
	string chrNumber = chr;
	if ( ( pos = chrNumber.find("chr", pos)) != string::npos )
		chrNumber.replace( pos, 3, "" );
	map<string,ChrCopyNumber>::iterator it;
	int index = findIndex(chrNumber);
	if (index == NA) {return;}
	//cout << "..index found "<<index<<"\n";
	int length = chrCopyNumber_[index].getLength();
	//cout <<length<<" == "<<chrCopyNumber_[index].getValues().size() <<"\n";
	for (int i = 0; i< length; i++) {
        float ratioToPrint = chrCopyNumber_[index].getRatioAtBin(i);
        if (isRatioLogged_ && ratioToPrint!= NA) {
            ratioToPrint=pow(2,ratioToPrint);
        }
        if (printNA || ratioToPrint!=NA) {//process this this window
                file << chrNumber <<"\t"<<chrCopyNumber_[index].getCoordinateAtBin(i)+1<<"\t"<<ratioToPrint ;
                if (chrCopyNumber_[index].isMedianCalculated()) {
                    file << "\t"<<chrCopyNumber_[index].getMedianProfileAtI(i) ;
                    float valueToPrint;
                    if (chrCopyNumber_[index].isSmoothed())
                        valueToPrint = chrCopyNumber_[index].getSmoothedProfileAtI(i)*ploidy_;
                    else
                        valueToPrint = round_by_ploidy(chrCopyNumber_[index].getMedianProfileAtI(i),ploidy_)*ploidy_;
                    if (valueToPrint<0)
                        valueToPrint = NA;
                    if (hasBAF_) {
                         if (valueToPrint>0 && chrCopyNumber_[index].getBAFsymbolAt(i).compare("-")!=0 &&chrCopyNumber_[index].getBAFsymbolAt(i)!="" && chrCopyNumber_[index].getBAFsymbolAt(i).length()!=valueToPrint ) { //&& (chrCopyNumber_[index].getEstimatedBAFuncertaintyAtI(i)< MAXUncertainty)
                    //change Level value if there is no incertainty:
                            valueToPrint=chrCopyNumber_[index].getBAFsymbolAt(i).length();
                         }
                    }
                    file <<"\t"<< valueToPrint;
                }
            if (hasBAF_) {
                file <<"\t"<< chrCopyNumber_[index].getBAFat(i);
                file <<"\t"<< (1-chrCopyNumber_[index].getBAFProfileAt(i));

                if (chrCopyNumber_[index].getBAFsymbolAt(i)!="")
                    file <<"\t"<< chrCopyNumber_[index].getBAFsymbolAt(i);
                else
                    file <<"\t0";
                file <<"\t"<< chrCopyNumber_[index].getEstimatedBAFuncertaintyAtI(i);
            }

            if (WESanalysis == true && chrCopyNumber_[index].getGeneNameAtBin(i)!= "")   {
                file <<  "\t" << chrCopyNumber_[index].getGeneNameAtBin(i);
            }
            if (SeekingSubc_ == true)  {
                file << "\t" << chrCopyNumber_[index].getCN_subc(i) << "\t" << chrCopyNumber_[index].getPopulation_subc(i); //check that it is still there
            }

            file << "\n";
        }
	}
	cout <<"..writing ratio for " <<chr <<"->done!\n";
}

void GenomeCopyNumber::printBAF(std::string const& chr, std::ofstream & file, SNPatChr& snpAtChrom, std::string myName, std::string const& matefile) {
	string::size_type pos = 0;
	string chrNumber = chr;
	if ( ( pos = chrNumber.find("chr", pos)) != string::npos )
		chrNumber.replace( pos, 3, "" );
	map<string,ChrCopyNumber>::iterator it;
	int index = findIndex(chrNumber);
	if (index == NA) {return;}
	int length = snpAtChrom.getSize();
    for (int i = 0; i< length; i++) {

        if (snpAtChrom.getStatusAt(i)!=0 && snpAtChrom.getValueAt(i)>=0) {
            int position = snpAtChrom.getPositionAt(i);

            file << chrNumber <<"\t"<<position<<"\t"<<snpAtChrom.getValueAt(i)<<"\t" ;
            int WinNumber;
            if (WESanalysis == false)  {
                WinNumber = position/step_;
            }  else  {
                WinNumber = snpAtChrom.getBinAt(i);
               // findWinNumber(position, myName, matefile);
            }
            float Avalue = NA;
            if (WinNumber!=NA)
                Avalue = chrCopyNumber_[index].getBAFProfileAt(WinNumber);
            float Bvalue = 1-Avalue;
            if (Avalue==NA) {Bvalue=NA;}

            float Afitted = NA;
            if (WinNumber!=NA)
               Afitted=chrCopyNumber_[index].getFittedBAFProfileAt(WinNumber);
            float Bfitted = 1-Afitted;
            if (Afitted==NA) {Bfitted=NA;}

            file << Afitted <<"\t"<<Bfitted<<"\t"<< Avalue <<"\t"<<Bvalue<<"\t"<<chrCopyNumber_[index].getEstimatedBAFuncertaintyAtI(WinNumber)<<"\n";

        }
    }
	cout <<"..writing BAF profiles for " <<chr <<"->done!\n";
}
/*
int GenomeCopyNumber::findWinNumber(int position, std::string myName, std::string const& matefile)
{
    string line;
    string pileup;
    pileup = myName + "_sample.cpn";
    ifstream file (pileup.c_str());
    int WinNumber = 0;
    bool found = false;
    while (std::getline(file,line) && found == false)
    {
        int i = 0;
        while (line[i] != '\t')
            {
            i++;
            }
        i++;
        string snp_pos;
        while (line[i] != '\t')
            {
            snp_pos += line[i];
            i++;
            }
        if (atoi(snp_pos.c_str()) >= position)
            {
            found  = true;
            }
        WinNumber++;
    }
    file.close();
    return WinNumber;
}
*/
void GenomeCopyNumber::printCGprofile(std::string const& chr, std::ofstream & file) {
	string::size_type pos = 0;
	string chrNumber = chr;
	if ( ( pos = chrNumber.find("chr", pos)) != string::npos )
		chrNumber.replace( pos, 3, "" );
	map<string,ChrCopyNumber>::iterator it;
	int index = findIndex(chrNumber);
	if (index == NA) {return;}
	int length = chrCopyNumber_[index].getLength();
	if (chrCopyNumber_[index].getMappabilityLength()>0)
        for (int i = 0; i< length; i++)
            file << chrNumber <<"\t"<< chrCopyNumber_[index].getCoordinateAtBin(i) <<"\t"<< chrCopyNumber_[index].getCGprofileAt(i) <<"\t"<<chrCopyNumber_[index].getNotNprofileAt(i) <<"\t"<< chrCopyNumber_[index].getMappabilityProfileAt(i) << "\n";
    else
        for (int i = 0; i< length; i++)
            file << chrNumber <<"\t"<< chrCopyNumber_[index].getCoordinateAtBin(i) <<"\t"<< chrCopyNumber_[index].getCGprofileAt(i) <<"\t"<<chrCopyNumber_[index].getNotNprofileAt(i) << "\n";
}

void GenomeCopyNumber::printCopyNumber(std::string const& chr, std::ofstream & file) {
	string::size_type pos = 0;
	string chrNumber = chr;
	if ( ( pos = chrNumber.find("chr", pos)) != string::npos )
		chrNumber.replace( pos, 3, "" );
	map<string,ChrCopyNumber>::iterator it;
	int index = findIndex(chrNumber);
	if (index == NA) {return;}
	int length = chrCopyNumber_[index].getLength();
	if (WESanalysis != false)
        {
        for (int i = 0; i< length; i++)
            file << chrNumber <<"\t"<< chrCopyNumber_[index].getCoordinateAtBin(i) <<"\t"<< chrCopyNumber_[index].getEndAtBin(i) <<"\t"<< chrCopyNumber_[index].getValueAt(i) << "\t" << chrCopyNumber_[index].getGeneNameAtBin(i) << "\n";
        }
    else {
	if (windowSize_ == step_) {
		for (int i = 0; i< length; i++)
            file << chrNumber <<"\t"<< chrCopyNumber_[index].getCoordinateAtBin(i) <<"\t"<< chrCopyNumber_[index].getValueAt(i) << "\n";
	} else {
        for (int i = 0; i< length; i++)
            file << chrNumber <<"\t"<< chrCopyNumber_[index].getCoordinateAtBin(i) <<"\t"<< chrCopyNumber_[index].getEndAtBin(i) <<"\t"<< chrCopyNumber_[index].getValueAt(i) << "\n";
	}
	}

}

void GenomeCopyNumber::calculateCopyNumberMedians (int minCNAlength, bool noisyData, bool CompleteGenomicsData) {

	cout << "..calculating medians for copy numbers\n";
	vector<ChrCopyNumber>::iterator it;
	for ( it=chrCopyNumber_.begin() ; it != chrCopyNumber_.end(); it++ ) {
        cout << "..calculating medians for " << it->getChromosome()<< "\n";
		it->calculateCopyNumberMedian(ploidy_, minCNAlength, noisyData, CompleteGenomicsData, isRatioLogged_);
	}
}

void GenomeCopyNumber::shiftNeutalRatioTo1() {

}

void GenomeCopyNumber::calculatePloidy(int minCNAlength) {
	//check if median for copy numbers is calculated
	if (chrCopyNumber_[0].getMedianValues().size()==0)
		calculateCopyNumberMedians(minCNAlength,0,0);
	cout << "..Getting ploidy information\n";
	int ploidy[] = {1,2,3,4,5,6,7,8,9,10};
	vector <double> scores;
	for (int i = 0; i < 6 ; i++) {
		scores.push_back(calculateXiSum(ploidy[i]));
	}
	int i = argmin(scores);
	ploidy_ = ploidy[i];
	cout << "..selected ploidy: " <<ploidy_<<"\n";
}

double GenomeCopyNumber::calculateXiSum(int ploidy) {
	int totalNumberOfBP = 0;
	double totalSum = 0;
	map <float,float> sds;
	//map <float,float> meds;
	//float varZero = calculateVarianceForNormalCopy(ploidy);
	//calculateSDAndMed(ploidy,sds,meds);
	calculateSDs(ploidy,sds);
	vector<ChrCopyNumber>::iterator it;
	cout << "..Calculating sum of squares for ploidy = "<<ploidy<<"\n";
	for ( it=chrCopyNumber_.begin() ; it != chrCopyNumber_.end(); it++ ) {
		//cout << "..processing chromosome " <<it->getChromosome()<<"\n";
		totalNumberOfBP += it->getNumberOfGoodFragments();
		//totalSum += it->getXiSum(ploidy,sqrt(varZero));
		//totalSum += it->calculateXiSum(ploidy,sds,meds);
		totalSum += it->calculateXiSum(ploidy,sds);
	}
	sds.clear();
	//meds.clear();
	cout << "ploidy\t" << ploidy << "\ttest-value\t"<<totalSum<<"\n";
	cout << "ploidy\t" << totalNumberOfBP << "\ttest-value\t"<<totalNumberOfBP<<"\n";
	return totalSum/totalNumberOfBP+15*ploidy;
}

void GenomeCopyNumber::calculateSDAndMed(int ploidy, map <float,float> &sds,map <float,float> &meds) {

	float value,med;
	vector<ChrCopyNumber>::iterator it;

	//fill the map
	map <float,vector <float> > mymap;
	for ( it=chrCopyNumber_.begin() ; it != chrCopyNumber_.end(); it++ ) {
		for (int i = 0; i<it->getLength(); i++) {
			med = it->getMedianProfileAtI(i);
			float level = round_by_ploidy(med, ploidy);
			value = it->getRatioAtBin(i);
			if (isRatioLogged_&& value!=NA)
                value=pow(2,value);
			if (value != NA) {
				if (mymap.count(level) == 0) {
					vector <float> a;
					mymap.insert ( pair<float,vector <float> >(level,a) );
				}
				mymap.find(level)->second.push_back(value);
			}
		}
	}

	//get median and variance for each level
	map<float,vector <float> >::iterator it2;
	for ( it2=mymap.begin() ; it2 != mymap.end(); it2++ ) {
		float level = (*it2).first;
		float median = get_median((*it2).second);
	//	float mean = get_mean((*it2).second);
		float sdev = sd((*it2).second,median);
		meds.insert(pair<float,float>(level,median));
		sds.insert(pair<float,float>(level,sdev));
		(*it2).second.clear();
	}
	mymap.clear();

}

void GenomeCopyNumber::calculateSDs(int ploidy, map <float,float> &sds) {

	float value,med;
	vector<ChrCopyNumber>::iterator it;

	//fill the map
	map <float,vector <float> > mymap;
	for ( it=chrCopyNumber_.begin() ; it != chrCopyNumber_.end(); it++ ) {
		for (int i = 0; i<it->getLength(); i++) {
			med = it->getMedianProfileAtI(i);
			float level = round_by_ploidy(med, ploidy);
			value = it->getRatioAtBin(i);
			if (isRatioLogged_&& value!=NA)
                value=pow(2,value);
			if (value != NA) {
				if (mymap.count(level) == 0) {
					vector <float> a;
					mymap.insert ( pair<float,vector <float> >(level,a) );
				}
				mymap.find(level)->second.push_back(value);
			}
		}
	}

	//get median and variance for each level
	map<float,vector <float> >::iterator it2;
	for ( it2=mymap.begin() ; it2 != mymap.end(); it2++ ) {
		float level = (*it2).first;
	//	float mean = get_mean((*it2).second);
		float sdev = sd((*it2).second,level);
		sds.insert(pair<float,float>(level,sdev));
		(*it2).second.clear();
	}
	mymap.clear();

}


float GenomeCopyNumber::calculateVarianceForNormalCopy(int ploidy) { //geting the total variance for points with median around 1 (+-0.5/ploidy) ---  OLD
	float variance = 0;
	int count = 0;
	double lowBoundary = 1-0.5/ploidy;
	double highBoundary = 1+0.5/ploidy;
	float value,med;
	vector<ChrCopyNumber>::iterator it;

	ofstream myfile;

	myfile.open ("example.txt");


	cout << "..Calculating variance for neutral copies\n";
	for ( it=chrCopyNumber_.begin() ; it != chrCopyNumber_.end(); it++ ) {
		//cout << "..processing chromosome " <<it->getChromosome()<<"\n";
		for (int i = 0; i<it->getLength(); i++) {
			med = it->getMedianProfileAtI(i);
			if ((med>lowBoundary)&&(med < highBoundary)) {
				value = it->getRatioAtBin(i);
				if (isRatioLogged_&& value!=NA)
                    value=pow(2,value);
				if (value != NA) {
					myfile << value-1 << "\n";
					variance += (value-1)*(value-1);
					count++;
				}
			}
		}
	}
	myfile.close();
	return float(variance/(count-1));
}

void GenomeCopyNumber::printPloidy(std::string const& outFile) {
	std::ofstream file (outFile.c_str());
	printPloidy(file);
	file.close();
}

void GenomeCopyNumber::printPloidy(std::ofstream & file) {
	file << "Sample ploidy is " << ploidy_ << " with p-value " << ploidy_pvalue_ << std::endl;
}

void GenomeCopyNumber::fillCGprofile(std::string const& chrFolder) {
	//reading the file with genome information
	vector<ChrCopyNumber>::iterator it;
	for ( it=chrCopyNumber_.begin() ; it != chrCopyNumber_.end(); it++ ) {
		it->fillCGprofile(chrFolder);
	}
}

float GenomeCopyNumber::evaluateContamination () {
	float contam = 0;
	string::size_type pos = 0;
	vector <float> values;
    vector <float> weights;
	map<string,int>::iterator it;
	for ( it=chromosomesInd_.begin() ; it != chromosomesInd_.end(); it++ ) {
		string chrNumber = (*it).first;
		if ( ( pos = chrNumber.find("chr")) != string::npos )
			chrNumber.replace( pos, 3, "" );
        if ( ( pos = chrNumber.find("X")) != string::npos ) 		//exclude X and Y from the analysis
            continue;
        if ( ( pos = chrNumber.find("Y")) != string::npos )
            continue;
		int index = findIndex(chrNumber);
		if (index == NA) {
		    cerr << "An error occurred in GenomeCopyNumber::evaluateContamination: could not find an index for "<<chrNumber<<"\n";
		    return 0;
        }
		int length = chrCopyNumber_[index].getLength();
		for (int i = 0; i< length; i++) {
			float observed = chrCopyNumber_[index].getRatioAtBin(i);
			if (observed!=NA) {
                if(isRatioLogged_)
                    observed=pow(2,observed);
                float expected = observed;
                if (chrCopyNumber_[index].isMedianCalculated()) {
                    expected = chrCopyNumber_[index].getMedianProfileAtI(i) ;
                    if (chrCopyNumber_[index].isSmoothed() && WESanalysis == false)
                        expected = chrCopyNumber_[index].getSmoothedProfileAtI(i);
                }
                if (!(expected == 1 || expected <= 0 || expected >= 1+2.0/ploidy_ || observed > 3 || observed <= 0)
                    && (((1>observed)&&(1>expected))||((1<observed)&&(1<expected)))) {// should it be something related to ploidy_ and not 2
                    float p = (observed-expected)/(observed-expected+2/ploidy_*(1-observed));
                    if (p>-0.5 && p<1.5) {
                        values.push_back(p);
                        weights.push_back(chrCopyNumber_[index].getFragmentLengths_notNA_At(i));
                    }
                }
			}
		}
	}
	contam = get_median(values);
	contam = get_mean(values);
    contam = get_weighted_mean(values,weights);

    if (contam<0) {
        cout << "\t..Evaluation of contamination using contamination = (observedLevel-expectedLevel)/(1-expectedLevel) produced a negative value: "<<contam<<"\n";
        cout << "\t..This probably means that there is not contamination at all or your \"ploidy\" value is not correct\n Will continue with contamination = 0\n";
        contam=0;
    }


	values.clear();
	return contam;
}

float GenomeCopyNumber::evaluateContaminationwithLR () {
	string::size_type pos = 0;
	vector <float> observed_values;
    vector <float> expected_values;
	map<string,int>::iterator it;

	for ( it=chromosomesInd_.begin() ; it != chromosomesInd_.end(); it++ ) {
		string chrNumber = (*it).first;
		if ( ( pos = chrNumber.find("chr")) != string::npos )
			chrNumber.replace( pos, 3, "" );
        if ( ( pos = chrNumber.find("X")) != string::npos ) 		//exclude X and Y from the analysis
            continue;
        if ( ( pos = chrNumber.find("Y")) != string::npos )
            continue;
		int index = findIndex(chrNumber);
		if (index == NA) {
		    cerr << "An error occurred in GenomeCopyNumber::evaluateContamination: could not find an index for "<<chrNumber<<"\n";
		    return 0;
        }
		int length = chrCopyNumber_[index].getLength();
		for (int i = 0; i< length; i++) {
			float observed = chrCopyNumber_[index].getRatioAtBin(i);
			if (observed!=NA) {
                if (isRatioLogged_)
                    observed=pow(2,observed);
                float expected = observed;
                if (chrCopyNumber_[index].isMedianCalculated()) {
                    expected = round_by_ploidy(chrCopyNumber_[index].getMedianProfileAtI(i),ploidy_) ;
                }
                if (expected!=NA && expected<1+2.0/ploidy_){
                    observed_values.push_back(observed-1);
                    expected_values.push_back(expected-1);
                }
			}
		}
	}

    vector <float> x = expected_values;
    vector <float> y = observed_values;

    int n = x.size();
    float sum_x = 0;
    float sum_y = 0;
    float sum_x2 = 0;
    float sum_y2 = 0;
    float sum_xy = 0;
    for (int i = 0; i < n; i++)     {
        sum_x += x[i];
        sum_y += y[i];
        sum_x2 += x[i] * x[i];
        sum_y2 += y[i] * y[i];
        sum_xy += x[i] * y[i];
    }

    float a = (n * sum_xy - sum_x * sum_y) / (n * sum_x2 - sum_x * sum_x); //CARINO
    a =sum_xy/sum_x2; //VALENTINA
      //  float b = (sum_x2 * sum_y - sum_x * sum_xy) / (n * sum_x2 - sum_x * sum_x);
    float contam_a =  (1-a)/(1-a+2*a/(float)ploidy_);
      //  float contam_b = -b/(-b+(2*b/ploidy_)-(2/ploidy_));
        /*cerr << a << "\t" << b  << "\n";
        cerr << contam_a << "\t" << contam_b;*/
        // contam = 2*(1 - a); //contam_a/2 + contam_b/2; CARINO
    if (contam_a<0) {
        cout << "\t..Evaluation of contamination produced a negative value: "<<contam_a<<"\n";
        cout << "\t..This probably means that there is not contamination at all or your \"ploidy\" value is not correct\n Will continue with contamination = 0\n";
        contam_a=0;
    }
    if (contam_a > 1) {
        cout << "\t..Evaluation of contamination produced a value over 1: "<<contam_a<<"\n";
        cout << "\t..This probably means that there is not contamination at all or your \"ploidy\" value is not correct\n Will continue with contamination = 0\n";
        contam_a=0;
    }
	return contam_a;
}


GenomeCopyNumber::~GenomeCopyNumber(void)
{
	chrCopyNumber_.clear();
	CNVs_.clear();
	chromosomesInd_.clear();
	copyNumberProbs_.clear();
}


void GenomeCopyNumber::addBAFinfo(SNPinGenome & snpingenome) {
	string::size_type pos = 0;
	map<string,int>::iterator it;

    hasBAF_=true;

	for ( it=chromosomesInd_.begin() ; it != chromosomesInd_.end(); it++ ) {
		string chrNumber = (*it).first;
		if ( ( pos = chrNumber.find("chr", pos)) != string::npos )
			chrNumber.replace( pos, 3, "" );
		int index = findIndex(chrNumber);
		if (index == NA) {
		    cerr << "An error occurred in GenomeCopyNumber::addBAFinfo: could not find an index for "<<chrNumber<<"\n";
		    exit(-1);
        }
        int indexSNP = snpingenome.findIndex(chrNumber);
        if (indexSNP == NA) {
		    cerr << "An error occurred in GenomeCopyNumber::addBAFinfo: could not find an SNP index for "<<chrNumber<<"\n";
		    exit(-1);
        }
        cout << "..Calculate BAF per window for chr"<< chrNumber << "\n";

        chrCopyNumber_[index].addBAFinfo(snpingenome,indexSNP);

		//int length = chrCopyNumber_[index].getLength();

		//create a vector with BAF
//		chrCopyNumber_[index].createBAF(NA);
//        int totalSNPnumber = snpingenome.SNP_atChr(indexSNP).getSize() ;
//		int SNPcount = 0;
//		float currentBAF = snpingenome.SNP_atChr(indexSNP).getValueAt(SNPcount);
//		int getSNPpos = snpingenome.SNP_atChr(indexSNP).getPositionAt(SNPcount);
//        float minBAF;
//		for (int i = 0; i<length; i++) {
//		    int left = chrCopyNumber_[index].getCoordinateAtBin(i);
//		    int right = chrCopyNumber_[index].getEndAtBin(i);
//		    if (getSNPpos>=left && getSNPpos <=right) {
//		        minBAF = chrCopyNumber_[index].getBAFat(i);
//                if (minBAF==NA) {
//                    chrCopyNumber_[index].setBAFat(i,currentBAF);
//                } else {
//                    chrCopyNumber_[index].setBAFat(i,min(minBAF,currentBAF));
//                }
//		    } else if (getSNPpos<left) {
//                if (SNPcount < totalSNPnumber) {
//                    SNPcount++;
//                    getSNPpos = snpingenome.SNP_atChr(indexSNP).getPositionAt(SNPcount);
//                    currentBAF = snpingenome.SNP_atChr(indexSNP).getValueAt(SNPcount);
//                    i=max(-1,i-windowSize_/step_-1);
//                    //cout << SNPcount << " out of "<< totalSNPnumber<<"\n";
//                }
//            }
//		}
	}
}

/*
int GenomeCopyNumber::processRead(std::string const& inputFormat, std::string const& matesOrientation,std::string const line) {
    if (! line.length()) return 0;
    int valueToReturn = 0;

	if ((inputFormat.compare("sam")==0 || inputFormat.compare("SAM")==0 )&&(matesOrientation.compare("0")!=0)){

        if ( line[0] == '@')
            return 0;
        string chr1, chr2;
		char orient1, orient2;
        int left,right ;
        if (getSAMinfo(line.c_str(),chr1,chr2,orient1,orient2,left,right)) {
		  char orient1_2[] = {orient1, orient2, 0};
		  char orient2_1[] = {orient2, orient1, 0};
		  if ((!strcmp(matesOrientation.c_str(), orient1_2) && (right-left>0)) || (!strcmp(matesOrientation.c_str(), orient2_1) && (right-left<0))) {
                    int index = findIndex(chr1);
                    if (index!=NA) {
                                chrCopyNumber_[index].mappedPlusOneAtI(min(left,right),step_);
                                valueToReturn=1;
                    }
            }
        }
        return valueToReturn;
	}

    string chr1, chr2, orient1, orient2;
    int left, right,insertSize;

	if ((inputFormat.compare("eland")==0 || inputFormat.compare("Eland")==0)&&(matesOrientation.compare("0")!=0)) {
            if (getELANDinfo(line,chr1,chr2,orient1,orient2,left,right,insertSize)) {
                if ((matesOrientation.compare(orient1+orient2)==0 && insertSize>0) || (matesOrientation.compare(orient1+orient2)==0 && insertSize<0)) {
					int ind = findIndex(chr1);
					if (ind!=NA) {
                        chrCopyNumber_[ind].mappedPlusOneAtI(min(left,right), step_);
                        valueToReturn= 1;
                    }
				}
            }
            return valueToReturn;
	}

	if ((inputFormat.compare("bowtie")==0 || inputFormat.compare("Bowtie")==0)&&(matesOrientation.compare("0")==0)){

			std::vector<std::string> strs = split(line, '\t');
			if (strs.size() > 4) {
				string chr = strs[2];
				processChrName(chr);
				int index = findIndex(chr);
				if (index!=NA) {
                    int left = atoi(strs[3].c_str());
                    chrCopyNumber_[index].mappedPlusOneAtI(left, step_);
                    valueToReturn= 1;
                }
			}
			strs.clear();
			return valueToReturn;

	}
	if ((inputFormat.compare("psl")==0 || inputFormat.compare("BLAT")==0)&&(matesOrientation.compare("0")==0)){

			std::vector<std::string> strs = split(line, '\t');
			if (strs.size() > 17) {
				string chr = strs[13];
				processChrName(chr);
				int index = findIndex(chr);
				if (index!=NA) {
                    int left = atoi(strs[15].c_str());
                    chrCopyNumber_[index].mappedPlusOneAtI(left,step_);
                    valueToReturn= 1;
				}
			}
			strs.clear();
			return valueToReturn;

	}
	if ((inputFormat.compare("arachne")==0 || inputFormat.compare("BED")==0 || inputFormat.compare("bed")==0 || inputFormat.compare("ARACHNE")==0)&&(matesOrientation.compare("0")==0)){
			std::vector<std::string> strs = split(line, '\t');
			if (strs.size() > 1) {
				string chr = strs[0];
				processChrName(chr);
				int index = findIndex(chr);
				if (index!=NA){
                    int left = atoi(strs[1].c_str());
                    chrCopyNumber_[index].mappedPlusOneAtI(left,step_);
                    valueToReturn= 1;
				}

			} else {
                strs.clear();
                strs = split(line, ' ');
                if (strs.size() > 1) {
                    string chr = strs[0];
                    processChrName(chr);
                    int index = findIndex(chr);
                    if (index!=NA) {
                        int left = atoi(strs[1].c_str());
                        chrCopyNumber_[index].mappedPlusOneAtI(left,step_);
                        valueToReturn= 1;
                    }
                }
			}
			strs.clear();
			return valueToReturn;

	}
	if ((inputFormat.compare("SOAP")==0 || inputFormat.compare("soap")==0 || inputFormat.compare("Soap")==0)&&(matesOrientation.compare("0")==0)){

			std::vector<std::string> strs = split(line, '\t');
			if (strs.size() > 8) {
				string chr = strs[7];
				processChrName(chr);
				int index = findIndex(chr);
				if (index!=NA) {
                    int left = atoi(strs[8].c_str());
                    chrCopyNumber_[index].mappedPlusOneAtI(left,step_);
                    valueToReturn= 1;
				}
			}
			strs.clear();
            return valueToReturn;

	}
	if ((inputFormat.compare("SOAP")==0 || inputFormat.compare("soap")==0 || inputFormat.compare("Soap")==0)&&(matesOrientation.compare("0")!=0)){

			std::vector<std::string> strs = split(line, '\t');
			if (strs.size() > 8) {
				string chr = strs[7];
				processChrName(chr);
				int index = findIndex(chr);
				if (index!=NA) {
                    int left = atoi(strs[8].c_str());
                    chrCopyNumber_[index].mappedPlusOneAtI(left,step_);
                    valueToReturn= 1;
				}
			}
			strs.clear();
            return valueToReturn;

	}


	if ((inputFormat.compare("sam")==0 || inputFormat.compare("SAM")==0 )&&(matesOrientation.compare("0")==0)){

        if ( line[0] == '@')
            return 0;

        string chr1,chr2,orient1,orient2;
        std::vector<std::string> strs = split(line, '\t');
        if (strs.size() > 3) {
                    string chr = strs[2];
                    processChrName(chr);
                    if (chr.compare("*")!=0 ) {
                        int left = atoi(strs[3].c_str());
                        int index = findIndex(chr);
                        if (index!=NA) {
                            chrCopyNumber_[index].mappedPlusOneAtI(left,step_);
                            valueToReturn=1;
                        }
                    }
        }
        strs.clear();
        return valueToReturn;
	}

	if ((inputFormat.compare("pileup")==0 || inputFormat.compare("SAMtools pileup")==0)){
			if (line[0] == '#') return 0;
			if ( line[0] == '@') return 0;

//chr1	10147	c	2	,.	`!
//chr1	10148	c	3	,$.^~,	b;9
//chr1	10149	c	2	.,	PI
//chr1	10150	c	2	.$,	H1
// we will simply count "^"

			std::vector<std::string> strs = split(line, '\t');
			if (strs.size() > 4) {
			    if (valueToReturn=strccnt(strs[4].c_str(), '^')) {
                    string chr = strs[0];
                    processChrName(chr);
                    int left = atoi(strs[1].c_str());
                    int index = findIndex(chr);
                    if (index==NA)
                        valueToReturn=0;
                    else {
                        for (int i=0; i<valueToReturn; i++)
                            chrCopyNumber_[index].mappedPlusOneAtI(left,step_);
                    }
			    }
			}
			strs.clear();
            return valueToReturn;

	}

    cerr << "Format of the input file ( "<<inputFormat<<" )was not recognized..\n";
    exit (-1);
}
*/

int GenomeCopyNumber::processRead(InputFormat inputFormat, MateOrientation matesOrientation, const char* line_buffer, int & prevInd,std::string targetBed, std::string mateFileName)
{

    int read_Size =150; // in case it is not initialized (e.g. for Pileup files)

    if (!*line_buffer) {
        return 0;
    }

    int valueToReturn = 0;

    if (inputFormat == SAM_INPUT_FORMAT && matesOrientation != SINGLE_END_SORTED_SAM) {
        if (line_buffer[0] == '@')
            return 0;
        if (WESanalysis == true)  {
            string chr1, chr2;
            char orient1, orient2;
            int left,right, insert_size;

            if (getSAMinfo(line_buffer,chr1,chr2,orient1,orient2,left,right, insert_size)) {
                char orient1_2_str[] = {orient1, orient2, 0};
                char orient2_1_str[] = {orient2, orient1, 0};
                MateOrientation orient1_2 = getMateOrientation(orient1_2_str);
                MateOrientation orient2_1 = getMateOrientation(orient2_1_str);
                if ((matesOrientation == orient1_2 && right-left>0) || (matesOrientation == orient2_1 && right-left<0)) {
                    left = min(left, right);
                    right  = left + abs(insert_size);
                    int index = findIndex(chr1);
                    if (index!=NA && chrCopyNumber_[index].getExons_Countchr() != 0) {
                        bool leftIsInTheWindow = false;
                        int l = 0;
                       // try around this value first
                        for (l=prevInd-1; l<=prevInd+30 && l<chrCopyNumber_[index].getExons_Countchr();l++ ) {
                            if ((left - 1 < chrCopyNumber_[index].getEndAtBin(l) && (left > (chrCopyNumber_[index].getCoordinateAtBin(l) - insert_size)))) {
                                leftIsInTheWindow = true;
                                prevInd=l;
                                break;
                            }
                        }
                        if (!leftIsInTheWindow) {
                             //start from the beginning of the chromosome
                            for (l=0; left>=chrCopyNumber_[index].getCoordinateAtBin(l)- 30*insert_size && l<chrCopyNumber_[index].getExons_Countchr();l++ ) {
                                if ((left - 1 < chrCopyNumber_[index].getEndAtBin(l) && (left > (chrCopyNumber_[index].getCoordinateAtBin(l) - insert_size)))) {
                                    leftIsInTheWindow = true;
                                    prevInd=l;
                                    break;
                                }
                            }
                        }
                        if (leftIsInTheWindow == false)   {
                            valueToReturn = 0;
                            return valueToReturn;
                        }
                        if ((right >  chrCopyNumber_[index].getCoordinateAtBin(l)) && (leftIsInTheWindow == true))   {
                            step_ = chrCopyNumber_[index].getEndAtBin(l) - chrCopyNumber_[index].getCoordinateAtBin(l) + insert_size +1;
                            chrCopyNumber_[index].mappedPlusOneAtI(left,step_, l);
                            valueToReturn = 1;
                            step_ = 0;
                        }
                        if (right < chrCopyNumber_[index].getCoordinateAtBin(l))  {
                            valueToReturn = 0;
                            return valueToReturn;
                        }
                    }
                }
            }
        } else  {
            string chr1, chr2;
            char orient1, orient2;
            int left,right ;
            if (getSAMinfo(line_buffer,chr1,chr2,orient1,orient2,left,right))  {
                char orient1_2_str[] = {orient1, orient2, 0};
                char orient2_1_str[] = {orient2, orient1, 0};
                MateOrientation orient1_2 = getMateOrientation(orient1_2_str);
                MateOrientation orient2_1 = getMateOrientation(orient2_1_str);
                if ((matesOrientation == orient1_2 && right-left>=0) || (matesOrientation == orient2_1 && right-left<=0)) {
                    int index = findIndex(chr1);
                    if (index!=NA)  {
                        chrCopyNumber_[index].mappedPlusOneAtI(min(left,right),step_);
                        valueToReturn=1;
                    }
                }
            }
        }
        return valueToReturn;
    }

    if (inputFormat == SAM_INPUT_FORMAT && matesOrientation == SINGLE_END_SORTED_SAM)  {

        if (line_buffer[0] == '@')
            return 0;

        // EV: to be optimized
        string chr1,chr2;
        char* strs[32];
        unsigned int strs_cnt = split((char*)line_buffer, '\t', strs);
        if (strs_cnt > 3) {
            if (WESanalysis == false) {
                string chr = strs[2];
                processChrName(chr);
                if (chr.compare("*")!=0 ) {
                    int left = atoi(strs[3]);
                    int index = findIndex(chr);
                    if (index!=NA) {
                        chrCopyNumber_[index].mappedPlusOneAtI(left,step_);
                        valueToReturn=1;
                    }
                }
            } else {
                string sequence = strs[9];
                read_Size = sequence.size();
                string chr = strs[2];
                processChrName(chr);
                int index = findIndex(chr);
                //chr = "chr" + chr; //I do not understand why it was written here. Valentina
                int l = 0;
                if (chr.compare("*")!=0 ) {
                    int left = atoi(strs[3]);
                    if (index!=NA) {
                        int right = left + read_Size;
                        bool leftIsInTheWindow = false;
                        if ((left - 1 < chrCopyNumber_[index].getEndAtBin(l) && (left > (chrCopyNumber_[index].getCoordinateAtBin(l) - read_Size))))   {
                            leftIsInTheWindow = true;
                        }   else   {
                            while ((!((left - 1 < chrCopyNumber_[index].getEndAtBin(l)) && (left > (chrCopyNumber_[index].getCoordinateAtBin(l) - read_Size))))  && (l <  chrCopyNumber_[index].getExons_Countchr()))
                            {
                                leftIsInTheWindow = false;
                                l++;
                            }
                            if (l < chrCopyNumber_[index].getExons_Countchr())  {
                                leftIsInTheWindow = true;
                            } else if (l == (chrCopyNumber_[index].getExons_Countchr() - 1)) {
                                leftIsInTheWindow = false;
                            }
                        }
                        if (leftIsInTheWindow == false) {
                            valueToReturn = 0;
                            return valueToReturn;
                        }
                        if ((right >  chrCopyNumber_[index].getCoordinateAtBin(l)) && (leftIsInTheWindow == true)) {
                            step_ = chrCopyNumber_[index].getEndAtBin(l) - chrCopyNumber_[index].getCoordinateAtBin(l) + read_Size +1;
                            chrCopyNumber_[index].mappedPlusOneAtI(left,step_, l);
                            valueToReturn = 1;
                            step_ = 0;
                        }
                        if (right < chrCopyNumber_[index].getCoordinateAtBin(l)) {
                            valueToReturn = 0;
                            return valueToReturn;
                        }
                    }
                }
            }
        }
        return valueToReturn;
    }


    if (inputFormat == SAM_PILEUP_INPUT_FORMAT) {
        if (line_buffer[0] == '#') return 0;
        if (line_buffer[0] == '@') return 0;

        //chr1	10147	c	2	,.	`!
        //chr1	10148	c	3	,$.^~,	b;9
        //chr1	10149	c	2	.,	PI
        //chr1	10150	c	2	.$,	H1
        // we will simply count "^"

        char* strs[32];
        unsigned int strs_cnt = split((char*)line_buffer, '\t', strs);
        if (strs_cnt > 4) {
            valueToReturn=strccnt(strs[4], '^');
            if (valueToReturn) {
                string chr = string(strs[0]);
                processChrName(chr);
                int left = atoi(strs[1]);
                if (WESanalysis == false)   {
                    int index = findIndex(chr);
                    if (index==NA)
                        valueToReturn=0;
                    else {
                        for (int i=0; i<valueToReturn; i++)    {
                            chrCopyNumber_[index].mappedPlusOneAtI(left,step_);
                        }
                    }
                } else {
                    int l = 0;
                    int index = findIndex(chr);
                    if (index==NA)   {
                        valueToReturn = 0;
                        return valueToReturn;
                    }
                    bool leftIsInTheWindow = false;
                    if ((left - 1 < chrCopyNumber_[index].getEndAtBin(l) && (left > (chrCopyNumber_[index].getCoordinateAtBin(l) - read_Size))))   {
                        leftIsInTheWindow = true;
                    }  else    {
                        while ((!((left - 1 < chrCopyNumber_[index].getEndAtBin(l)) &&
                            (left > (chrCopyNumber_[index].getCoordinateAtBin(l) - read_Size))))  &&
                            (l <  chrCopyNumber_[index].getExons_Countchr()))    {
                            leftIsInTheWindow = false;
                            l++;
                        }   if (l < chrCopyNumber_[index].getExons_Countchr() -1 )    {
                            leftIsInTheWindow = true;
                        }
                    }
                    if (leftIsInTheWindow == false)  {
                        valueToReturn = 0;
                        return valueToReturn;
                    }   else    {
                        step_ = chrCopyNumber_[index].getEndAtBin(l) - chrCopyNumber_[index].getCoordinateAtBin(l) + read_Size +1;
                        for (int i=0; i<valueToReturn; i++)    {
                            chrCopyNumber_[index].mappedPlusOneAtI(left,step_, l);
                        }
                        step_ = 0;
                        return valueToReturn;
                    }
                }
            }
        }
        return valueToReturn;
    }

    if (inputFormat == ELAND_INPUT_FORMAT && matesOrientation != SINGLE_END_SORTED_SAM) {
        string chr1, chr2, orient1, orient2;
        int left, right,insertSize;

        if (getELANDinfo(line_buffer,chr1,chr2,orient1,orient2,left,right,insertSize)) {
          MateOrientation orient1_2 = getMateOrientation(orient1+orient2);
          MateOrientation orient2_1 = getMateOrientation(orient2+orient1);
          if ((matesOrientation == orient1_2 && insertSize>0) || (matesOrientation == orient2_1 && insertSize<0)) {
            // EV: warning orient2_1 not used: port "as is" from 5.9 version, bug ? - it was a bug! Valentina
            int ind = findIndex(chr1);
            if (ind!=NA) {
              chrCopyNumber_[ind].mappedPlusOneAtI(min(left,right), step_);
              valueToReturn= 1;
            }
          }
        }
        return valueToReturn;
    }

    if (inputFormat == BOWTIE_INPUT_FORMAT && matesOrientation == SINGLE_END_SORTED_SAM) {
        std::vector<std::string> strs = split(line_buffer, '\t');
        if (strs.size() > 4) {
          string chr = strs[2];
          processChrName(chr);
          int index = findIndex(chr);
          if (index!=NA) {
            int left = atoi(strs[3].c_str());
            chrCopyNumber_[index].mappedPlusOneAtI(left, step_);
            valueToReturn= 1;
          }
        }
        strs.clear();
        return valueToReturn;
    }

    if (inputFormat == PSL_INPUT_FORMAT && matesOrientation == SINGLE_END_SORTED_SAM) {
        std::vector<std::string> strs = split(line_buffer, '\t');
        if (strs.size() > 17) {
          string chr = strs[13];
          processChrName(chr);
          int index = findIndex(chr);
          if (index!=NA) {
            int left = atoi(strs[15].c_str());
            chrCopyNumber_[index].mappedPlusOneAtI(left,step_);
            valueToReturn= 1;
          }
        }
        strs.clear();
        return valueToReturn;
    }

    if (inputFormat == ARACHNE_BED_INPUT_FORMAT && matesOrientation == SINGLE_END_SORTED_SAM) {
        std::vector<std::string> strs = split(line_buffer, '\t');
        if (strs.size() > 1) {
          string chr = strs[0];
          processChrName(chr);
          int index = findIndex(chr);
          if (index!=NA){
            int left = atoi(strs[1].c_str());
            chrCopyNumber_[index].mappedPlusOneAtI(left,step_);
            valueToReturn= 1;
          }
        } else {
          strs.clear();
          strs = split(line_buffer, ' ');
          if (strs.size() > 1) {
            string chr = strs[0];
            processChrName(chr);
            int index = findIndex(chr);
            if (index!=NA) {
              int left = atoi(strs[1].c_str());
              chrCopyNumber_[index].mappedPlusOneAtI(left,step_);
              valueToReturn= 1;
            }
          }
        }
        strs.clear();
        return valueToReturn;
    }

    if (inputFormat == SOAP_INPUT_FORMAT && matesOrientation == SINGLE_END_SORTED_SAM) {
        std::vector<std::string> strs = split(line_buffer, '\t');
        if (strs.size() > 8) {
          string chr = strs[7];
          processChrName(chr);
          int index = findIndex(chr);
          if (index!=NA) {
            int left = atoi(strs[8].c_str());
            chrCopyNumber_[index].mappedPlusOneAtI(left,step_);
            valueToReturn= 1;
          }
        }
        strs.clear();
        return valueToReturn;
    }

    if (inputFormat == SOAP_INPUT_FORMAT && matesOrientation != SINGLE_END_SORTED_SAM) {
        std::vector<std::string> strs = split(line_buffer, '\t');
        if (strs.size() > 8) {
          string chr = strs[7];
          processChrName(chr);
          int index = findIndex(chr);
          if (index!=NA) {
            int left = atoi(strs[8].c_str());
            chrCopyNumber_[index].mappedPlusOneAtI(left,step_);
            valueToReturn= 1;
          }
        }
        strs.clear();
        return valueToReturn;
    }

    cerr << "Format of the input file ( "<<inputFormat<<" )was not recognized..\n";
    exit (-1);

}

int GenomeCopyNumber::processReadWithBowtie(std::string const& inputFormat, std::string const& matesOrientation,std::string const line,std::string const line2) {

    int valueToReturn = 0;

    string chr1, chr2, orient1, orient2;

	if ((inputFormat.compare("bowtie")==0 || inputFormat.compare("Bowtie")==0)&&(matesOrientation.compare("0")!=0)){

			if (! line2.length()) return 0;
			std::vector<std::string> strs = split(line, '\t');
			std::vector<std::string> strs2 = split(line2, '\t');
			if ((strs.size() > 4)&&(strs2.size() > 4)) {
				string chr1 = strs[2];
				string chr2 = strs2[2];

				processChrName(chr1);
				processChrName(chr2);
				if (chr1.compare(chr2)==0) {

                    string orient1 = strs[1];
                    string orient2 = strs2[1];
                    if (orient1.compare("+")==0)
                        orient1 = "F";
                    else
                        orient1 = "R";

                    if (orient2.compare("+")==0)
                        orient2 = "F";
                    else
                        orient2 = "R";
                    int left = atoi(strs[3].c_str());
                    int right = atoi(strs2[3].c_str());

                    if ((matesOrientation.compare(orient1+orient2)==0 && (right-left>0)) || (matesOrientation.compare(orient2+orient1)==0 && (right-left<0))) {
                        int index = findIndex(chr1);
                        if (index!=NA) {
                            chrCopyNumber_[index].mappedPlusOneAtI(min(left,right),step_);
                            valueToReturn=1;
                        }
                    }
                }
			}
			strs.clear();
			strs2.clear();
            return valueToReturn;
	}
	return 0;
}
float GenomeCopyNumber::removeLargeExons(float iqrToKeep) {
    float maxLength = 0;
	vector<ChrCopyNumber>::iterator it;
	int totalNumberExons=0;
	float numberOfRemovedExons = 0;
	vector <float> exonLengths;

    for ( it=chrCopyNumber_.begin() ; it != chrCopyNumber_.end(); it++ ) {
        for (int i=0; i < it->getLength(); i++) {
            exonLengths.push_back(it->getEndAtBin(i)-it->getCoordinateAtBin(i));
            totalNumberExons++;
        }
	}
    maxLength=get_iqr(exonLengths)/2*iqrToKeep+get_median(exonLengths);

    for ( it=chrCopyNumber_.begin() ; it != chrCopyNumber_.end(); it++ ) {
        numberOfRemovedExons+=it->removeLargeExons(maxLength);
	}
	return numberOfRemovedExons/totalNumberExons;
}
int GenomeCopyNumber::focusOnCapture (std::string const& captureFile) {
	ifstream file (captureFile.c_str());
	string line;
    //string::size_type pos = 0;
    string currentChr = "";
    int index = NA;
    int count = 0;
    int mapCount;
    float ratio;
    int endShift = - windowSize_/step_ + 1;
    unsigned long int minRegion = refGenomeSize_;
    refGenomeSize_=0; //will be recalculated!!!

    map<int,int> chrRead;
    map<int,int>::iterator itMapChr;

    int averageReadLength=400; averageReadLength=150;//from version 6.6

	if (file.is_open())	{
	    cout << "..Reading "<< captureFile << "\n";
        cout << "..Your file must be in .BED format, and it must be sorted\n";

		while (! file.eof() )	{
			getline (file,line);
			if (! line.length()) continue;
            std::vector<std::string> strs = split(line, '\t');
            if (strs.size()<3) continue;

			string chr = strs[0];
            processChrName(chr);
            int positionS = atoi(strs[1].c_str())-averageReadLength; if (positionS<=0) positionS=1;
            int positionE = atoi(strs[2].c_str())-1; //because it is 1-based

            if (positionS >= positionE || positionS<0) continue;

			if ( currentChr.compare(chr)!=0 ){ //start reading a new chromosome
                cout << "..Reading capture for chromosome " << chr << "\n";
                index = findIndex(chr);
                if (index == NA) {
				    cout <<  "skipping chromosome " << chr << "\n";
				    continue;
                }
                //restore all variables for a new chromosome
                currentChr = chr;
                //check if the file is sorted:
                itMapChr=chrRead.find(index);
                if (itMapChr != chrRead.end()  ) {
                    cerr << "Your .bed file should be sorted. \n..Unable to proceed..\n";
                    exit(-1);
                }
                chrRead[index]=1;

                //check is notNprofile_ exists, and create it if it does not exists:
                chrCopyNumber_[index].checkOrCreateNotNprofileWithZeros();
                //cout <<  "..Index for chromosome " << currentChr << ": "<< index << "\n";
			}
            count ++;

            int leftWindow = positionS/step_;
            int rightWindow = positionE/step_ + endShift;
            if (rightWindow<leftWindow)
                rightWindow=leftWindow;
            unsigned long int RegLength = positionE - positionS;
            if (RegLength<minRegion)
                minRegion=RegLength;
            refGenomeSize_ += RegLength;

            for (int i=leftWindow; i <= rightWindow; i++) {
                mapCount=windowSize_;
                if (positionS>i*step_)
                    mapCount -= (positionS-i*step_);
                if (positionE<i*step_+windowSize_)
                    mapCount += (positionE-(i*step_+windowSize_));
                ratio = float(mapCount)/windowSize_;
                chrCopyNumber_[index].setNotNprofileAt(i, ratio);
            }
		}
		file.close();
		cout << "file " << captureFile << " is read\n";
	} else{
        cerr << "Error: Unable to open file "+captureFile+"\n";
        exit(-1);
	}
    //EXCLUDE all windows that are not in capture
    cout << "..Setting read counts to Zero for all windows outside of capture\n";
    vector<ChrCopyNumber>::iterator it=chrCopyNumber_.begin();

    cout << "..Total size of captured regions "<< refGenomeSize_<<"bp\n";

	for (  ; it != chrCopyNumber_.end(); it++ ) {
		cout << "..processing chromosome " <<it->getChromosome()<<"\n";
		it->setRCountToZeroForNNNN();
	}
    return int(minRegion);
}

void GenomeCopyNumber::setSamtools(std::string const& pathToSamtools) {
    pathToSamtools_=pathToSamtools;
}

void GenomeCopyNumber::setSambamba(std::string const& pathToSambamba, std::string const& SambambaThreads) {
    pathToSambamba_=pathToSambamba;
    SambambaThreads_ = SambambaThreads;
}

bool GenomeCopyNumber::getWESanalysis()
{
    return WESanalysis;
}

void GenomeCopyNumber::setWESanalysis(bool WESgiven)
{
WESanalysis = WESgiven;
}

void GenomeCopyNumber::setIfLogged(bool isRatioLogged) {
    isRatioLogged_=isRatioLogged;
}

void GenomeCopyNumber::setmakingPileup(bool makingPileup_given)
{
makingPileup = makingPileup_given;
}


void GenomeCopyNumber::setSeekSubclones(bool seekSubclones)
{
    SeekingSubc_ = seekSubclones;
    vector<ChrCopyNumber>::iterator it;
	for ( it=chrCopyNumber_.begin() ; it != chrCopyNumber_.end(); it++ )
            it->setLookingForSubclones(seekSubclones);
}

void* GenomeCopyNumber_readMateFile_wrapper(void *arg)
{
  GenomeCopyNumberReadMateFileArgWrapper* warg = (GenomeCopyNumberReadMateFileArgWrapper*)arg;
  if (warg->p_genomeCopyNumber) {
	warg->snpInGenome.readMateFile(warg->mateFile, warg->inputFormat, warg->minimalTotalLetterCountPerPosition, warg->minimalQualityPerPosition, *warg->p_genomeCopyNumber, warg->chrLenFileName, warg->windowSize, warg->step, warg->targetBed);
  } else {
	warg->snpInGenome.readMateFile(warg->mateFile, warg->inputFormat, warg->minimalTotalLetterCountPerPosition, warg->minimalQualityPerPosition);
  }
  return NULL;
}

void* GenomeCopyNumber_calculateBreakpoint_wrapper(void *arg)
{
  GenomeCopyNumberCalculateBreakpointArgWrapper* warg = (GenomeCopyNumberCalculateBreakpointArgWrapper*)arg;
  warg->genomeCopyNumber.calculateBreakpoints(warg->breakPointThreshold, warg->breakPointType);
  return NULL;
}

double GenomeCopyNumber::Percentage_GenomeExplained(int & unexplainedChromosomes)
{
	vector<ChrCopyNumber>::iterator it;
	long double fragment_median2 ;
	long double numberOfPoints=0;
	double sum_frags = 0;
	unexplainedChromosomes=0;
    bool unexplained=0;
    int threshold = 10;

	for ( it=chrCopyNumber_.begin() ; it != chrCopyNumber_.end(); it++ ) {
		if (! (it->getChromosome().find("X")!=string::npos || it->getChromosome().find("Y")!=string::npos)) {
			int fragmentLength=0;
			if (unexplained)unexplainedChromosomes++;
            unexplained=0;
            for (int i = 0; i< it->getNumberOfFragments(); i++) { // for each fragment:
                fragment_median2 = it->getMedianValuesAt(i);
                fragmentLength=it->getFragmentLengthsAt(i);
                if (fragmentLength>threshold) {
                    numberOfPoints+=fragmentLength;
                    if (fragment_median2!=NA && abs(fragment_median2-round_by_ploidy(fragment_median2,ploidy_)) >= 1.0/3/ploidy_) {
                        //cout <<  "Unexplained segment: "<<fragment_median2 << "\t";
                            unexplained=1;
                            sum_frags+=fragmentLength;
                    }
                }
            }
        }
	}
	cout << "Total proportion of unexplained regions: "<<sum_frags << " out of " << numberOfPoints<<" = " << sum_frags/numberOfPoints<< endl ;
	return (1-sum_frags/numberOfPoints);
}

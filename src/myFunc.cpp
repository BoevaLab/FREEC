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


#include "myFunc.h"
#include <assert.h>
#include <pthread.h>

#if defined(_WIN32) || (defined(__APPLE__) && defined(__MACH__))
//x32 Windows definitions
#include <limits.h>
#else
//other platforms
#include <values.h>
#endif


#include "ThreadPool.h"

using namespace std ;

myFunc::myFunc(void)
{
}

myFunc::~myFunc(void)
{
}

unsigned int split(char* str_ori, char delim, char* elems[])
{
  const char* str = str_ori;
  unsigned int last_jj = 0;
  unsigned int jj = 0;
  unsigned int elem_cnt = 0;
  char c;
  for (; c = *str++; ++jj) {
	if (c == delim) {
	  str_ori[jj] = 0;
	  elems[elem_cnt++] = &str_ori[last_jj];
	  last_jj = jj+1;
	}
  }
  elems[elem_cnt++] = &str_ori[last_jj];
  return elem_cnt;
}

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while(std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

bool checkChrLen(const std::string &chrLenFile,const std::string &targetBed) {
    std::vector<std::string> names;
	std::vector<int> lengths;
    std::vector<std::string> names_bed;
	readFileWithGenomeInfo(chrLenFile, names, lengths);
    readChrNamesInBed(targetBed, names_bed);
    bool toReturn = true;

    if(names.empty()){
        cerr << "Error:Cound not read "<< chrLenFile<<"\n";
        exit(1);
    }
    if(names_bed.empty()){
        cerr << "Error:Cound not read "<< targetBed<<"\n";
        exit(1);
    }

    for (int i=0; i<names.size();i++) {
        if(std::find(names_bed.begin(), names_bed.end(), names[i]) != names_bed.end()) {
                /* names_bed contains names[i]; everything is OK */
        } else {
                /* names_bed does not contain names[i] */
                toReturn = false;
                cerr << "Error: chromosome "<< names[i]<< " present in your "<<chrLenFile << " file was not detected in your file with capture regions " <<targetBed<<"\n";
                cerr << "Please solve this issue and rerun Control-FREEC\n";
                cerr << "For example, you can remove chromosome "<< names[i]<<" from your "<<chrLenFile<<"\n";
        }
    }
    return toReturn;
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    return split(s, delim, elems);
}

// Calculate sum
// ---------------------------------------------------------------------------
float get_sum(const std::vector<float>& data) {

	double sum = 0;
	int arrayLength = data.size();
	for (int i=0; i< arrayLength; i++)
			sum+=data[i];
	return (float)sum;
}

// Calculate median across individuals
// ---------------------------------------------------------------------------
float get_median(const std::vector<float> & myvector) {
  vector<float> data (myvector);
  int ndatapoints = data.size();
  if (ndatapoints==0) {
    cerr << "Error: zero values to calculate medians..\n";
    exit(-1);
  }
  sort(data.begin(),data.end());
  // Get median
  float median_value = ndatapoints % 2 == 1 ? data[(ndatapoints-1)/2] : (data[ndatapoints/2 - 1] + data[ndatapoints/2])/float(2.0);
  data.clear();
  return median_value;
}

float get_medianNotNA(const std::vector<float> & myvector) {
  vector<float> data (myvector);

  data.erase(std::remove(data.begin(), data.end(), NA), data.end());

  int ndatapoints = data.size();
  if (ndatapoints==0) {
    cerr << "Error: zero values to calculate medians..\n";
    exit(-1);
  }
  sort(data.begin(),data.end());
  // Get median
  float median_value = ndatapoints % 2 == 1 ? data[(ndatapoints-1)/2] : (data[ndatapoints/2 - 1] + data[ndatapoints/2])/float(2.0);
  data.clear();
  return median_value;
}

// Calculate sd across individuals around given mean
// ---------------------------------------------------------------------------
float get_sd (const std::vector<float>& data, float mean) {
    double sum = 0;
	int arrayLength = data.size();
	for (int i=0; i< arrayLength; i++)
			sum+=(data[i]-mean)*(data[i]-mean);
    return (float)sqrt(sum/arrayLength);
}

// Calculate argmin for a vector for the first smallest value in a range
// ---------------------------------------------------------------------------
int argmin(const std::vector<double> & myvector) {
	int argmin = 0;
	for (int i=1; i< (int)myvector.size(); i++)
			if (myvector[i]<myvector[argmin])
				argmin = i;

	return argmin;
}

// Calculate median across individuals for a subset
// ---------------------------------------------------------------------------
float get_median(const std::vector<float> & myvector, int start, int end) {
	if ((start>=0)&&(start<end)&&(end <= (int)myvector.size())) {
		int ndatapoints = end-start;
		vector<float> data (ndatapoints);
		for (int i = start; i < end; i++)
			data[i-start] = myvector[i];
        sort(data.begin(),data.end());
		// Get median
		float median_value = ndatapoints % 2 == 1 ? data[(ndatapoints-1)/2] : (data[ndatapoints/2 - 1] + data[ndatapoints/2])/float(2.0);
		data.clear();
		return median_value;
	}
	else {
		cout << "Wrong vector boundaries for median calculations\n";
		return 0;
	}
}

// Calculate mean across individuals
// ---------------------------------------------------------------------------
float get_mean(const std::vector<float>& data) {

	double sum = 0;
	int arrayLength = data.size();
	for (int i=0; i< arrayLength; i++)
			sum+=data[i];
	return (float)sum/arrayLength;
}

// Calculate weighted mean across individuals
// ---------------------------------------------------------------------------
float get_weighted_mean(const std::vector<float>& data, const std::vector<float>& weights) {
	double sum = 0;
	double totalWeights = 0;
	int unsigned arrayLength = data.size();
	if (arrayLength != weights.size()) {
		cerr << "Warning: arrays should have the same size in \"get_weighted_mean\"..\n";
		return NA;
	}
	for (int unsigned i=0; i< arrayLength; i++) {
        sum+=data[i]*weights[i];
        totalWeights += weights[i];
	}
	return (float)sum/totalWeights;
}


// Calculate inter quartile range across individuals
// ---------------------------------------------------------------------------
float get_iqr(const std::vector<float>& myvector) {
  vector<float> data (myvector);
  int ndatapoints = data.size();
  if (ndatapoints==0) {
    cerr << "Error: zero values to calculate medians..\n";
    exit(-1);
  }
  sort(data.begin(),data.end());
  // Get IQR
  float lower_quartile, upper_quartile;

  if( ndatapoints % 2 == 1){
    float fl = (ndatapoints-1)/float(4.0);
    float fu = 3*(ndatapoints-1)/float(4.0);

    if( fmod(fl,float(1.0)) == 0 ){
      lower_quartile = data[ (int) fl ];
      upper_quartile = data[ (int) fu ];
    }
    else{
      lower_quartile = (data[ (int) fl ] + data[ (int) fl + 1 ])/float(2.0);
      upper_quartile = (data[ (int) fu ] + data[ (int) fu + 1 ])/float(2.0);
    }
  }
  else{
    double fl = (ndatapoints)/float(4.0);
    double fu = 3*(ndatapoints)/float(4.0);

    if( fmod(fl,1.0) == 0 ){
      lower_quartile = data[ (int) fl ];
      upper_quartile = data[ (int) fu ];
    }
    else{
      lower_quartile = (data[ (int) fl ] + data[ (int) fl + 1 ])/float(2);
      upper_quartile = (data[ (int) fu ] + data[ (int) fu + 1 ])/float(2);
    }
  }

  return upper_quartile - lower_quartile;
}

void readChrNamesInBed(const std::string &targetBed, std::vector<std::string>&names_bed){
    ifstream file(targetBed.c_str());
	if (!file.is_open()) {
        cerr << "Error: unable to open "+targetBed+"\n" ;
        exit(-1);
	}
	string line;
	string name;
	bool isFai=0;
	while (std::getline(file,line)) {
		if (! line.length()) continue;
		if (line[0] == '#') continue;

		std::vector<std::string> strs = split(line, '\t');
		if (strs.size()<2) {
		    continue;
		}
        name  = strs[0];
		strs.clear();
		myReplace(name, " ", "");

		//delete "Chr"
		string::size_type pos = 0;
		if ( (pos = name.find("chr", pos)) != string::npos )
			name.replace( pos, 3, "" );

        names_bed.push_back(name);
	}
	file.close();
}


void readFileWithGenomeInfo(const std::string &chrLenFileName, std::vector<std::string>& names, std::vector<int>& lengths) {
	//reading the file with genome information
	ifstream file(chrLenFileName.c_str());
	if (!file.is_open()) {
        cerr << "Error: unable to open "+chrLenFileName+"\n" ;
        exit(-1);

	}
	string line;
	string name;
	int value = 0;
	bool isFai=0;
	if (chrLenFileName.substr(chrLenFileName.length()-3,3).compare("fai")==0) {isFai=1;}
	while (std::getline(file,line)) {

		if (! line.length()) continue;
		if (line[0] == '#') continue;

		std::vector<std::string> strs = split(line, '\t');
		if (strs.size()<2) {
		    cerr << "uncorrect file with chromosomes "<< chrLenFileName <<"\nUse tab-delimited format:\n1\tchr1\t249250621\nor\nchr1\t249250621\n";
		}
		if (strs.size()==2 || strs[0].substr(0,3)=="chr" || isFai) {
			name  = strs[0];
			value = atoi(strs[1].c_str());
		}
		if (strs.size()>=3 && strs[0].substr(0,3)!="chr" && !isFai) {
			name  = strs[1];
			value = atoi(strs[2].c_str());
		}
		strs.clear();
		myReplace(name, " ", "");

		//delete "Chr"
		string::size_type pos = 0;
		if ( (pos = name.find("chr", pos)) != string::npos )
			name.replace( pos, 3, "" );
		if (value>0) {
		    names.push_back(name);
            lengths.push_back(value);
            //cout << name << "\t" << value << "\n";
		}

	}
	file.close();
	cout << "..File "<<chrLenFileName<<" was read\n";
}
unsigned long sum(const std::vector<int>& data) {
	unsigned long sum = 0;
	for (int i=0; i< (int) data.size(); i++) {
		sum+=data[i];
		//cout << data[i] << "\n";
	}
	return sum;
}

long getLineNumber(std::string const& fileName, const std::string& pathToSamtools, const std::string& pathToSambamba, const std::string& SambambaThreads) {
	string line ;
	long count = 0;
	FILE *stream;
    char buffer[MAX_BUFFER];

    if (fileName.substr(fileName.size()-3,3).compare(".gz")==0) {
        string command = "gzip -cd "+fileName;
        stream =
        #if defined(_WIN32)
            _popen(command.c_str(), "r");
        #else
            popen(command.c_str(), "r");
        #endif
        while ( fgets(buffer, MAX_BUFFER, stream) != NULL ) {
            count++;
        }
        #if defined(_WIN32)
				_pclose(stream);
		#else
				pclose(stream);
		#endif
    } else if (fileName.substr(fileName.size()-4,4).compare(".bam")==0) {
        string command = "";
        if (pathToSambamba != "")
            {
            command = pathToSambamba + " view -t " + SambambaThreads + " " +fileName;
            //myInputFormat="sam";       //will try to use existing samtools
            cout << "..sambamba should be installed to be able to read BAM files\n";
            }
        else
            {
            command = pathToSamtools + " view "+fileName;
            //myInputFormat="sam";       //will try to use existing samtools
            cout << "..samtools should be installed to be able to read BAM files\n";
            }
        stream =
            #if defined(_WIN32)
				_popen(command.c_str(), "r");
			#else
				popen(command.c_str(), "r");
			#endif

        while ( fgets(buffer, MAX_BUFFER, stream) != NULL ) {
            count++;
        }
        #if defined(_WIN32)
				_pclose(stream);
		#else
				pclose(stream);
		#endif
    }else {
    	ifstream file(fileName.c_str()) ;
        while( getline( file, line ) ) count++ ;
        file.close();
    }


	return count;
}

void advance_to(const std::string& haystack, size_t& offset, char needle) {

	// I think std::string should behave correctly in this case, but just to be sure...
	if(offset >= haystack.size())
		return;
	offset = haystack.find(needle, offset);
	if(offset == std::string::npos)
		offset = haystack.size();

}

long getReadNumberFromPileup(std::string const& fileName) {
	string line ;
	long count = 0;
	std::vector<std::string> strs;
	// we will simply count "^"
#ifdef PROFILE_TRACE
	time_t t0 = time(NULL);
#endif
    if (fileName.substr(fileName.size()-3,3).compare(".gz")==0) {
        FILE *stream;
        char buffer[MAX_BUFFER];
        string command = "gzip -cd "+fileName;
        stream =
            #if defined(_WIN32)
				_popen(command.c_str(), "r");
			#else
				popen(command.c_str(), "r");
			#endif
        while ( fgets(buffer, MAX_BUFFER, stream) != NULL ) {
            line = buffer;
            if (! line.length()) continue;
            if (line[0] == '#') continue;
            if ( line[0] == '@') continue;

            strs = split(line, '\t');
            if (strs.size() > 4) {
                int toadd;
                if (toadd=strccnt(strs[4].c_str(), '^')) {
                        count += toadd;
                }
            }
            strs.clear();
        }
        #if defined(_WIN32)
				_pclose(stream);
		#else
				pclose(stream);
		#endif

    } else {
        ifstream fileMates(fileName.c_str()) ;
        while (std::getline(fileMates,line)) {

            if (! line.length()) continue;
            if (line[0] == '#') continue;
            if ( line[0] == '@') continue;

            strs = split(line, '\t');
            if (strs.size() > 4) {
                int toadd;
                if (toadd=strccnt(strs[4].c_str(), '^')) {
                        count += toadd;
                }
            }
            strs.clear();
        }
        fileMates.close();
    }

#ifdef PROFILE_TRACE
	std::cout << "PROFILING [tid=" << pthread_self() << "]: " << fileName << " read in " << (time(NULL)-t0) << " seconds [getReadNumberFromPileup]\n";
#endif
    return count;
}


int factorial (int num)
{
 if (num==1)
  return 1;
 return factorial(num-1)*num; // recursive call
}

//This function returns the index of the maximum value in the vector 'data'. When there
//are several equal maximum elements then the lowest index is returned.
int get_max_index(const std::vector<float>& data) {
	int ind = 0;
	float max = data[0];
	for (int i=1; i< (int) data.size(); i++) {
		if (max<data[i]) {
			max = data[i];
			ind = i;
		}
	}
	return ind;
}

//This function returns the index of the minimum value in the vector v. When there
//are several equal minimum elements then the lowest index is returned.
inline int get_min_index(const std::vector<float>& data, int size) {
	int ind = 0;
	float min = data[ind];
	for (int ii = 1; ii < size; ii++) {
		if (min > data[ii]) {
		    min = data[ii];
			ind = ii;
		}
	}
	return ind;
}

inline int get_min_index(const float* data, int size) {
	int ind = 0;
	float min = data[ind];
	for (int ii = 1; ii < size; ii++) {
		if (min > data[ii]) {
		    min = data[ii];
			ind = ii;
		}
	}
	return ind;
}

int get_min_index(const std::vector<float>& data) {
    return get_min_index(data, data.size());
}

//This function subtracts the elements of vector b from the elements of vector a, a0i =
//ai ?�� bi. The two vectors must have the same length.
void vector_sub (std::vector<float>& a, const std::vector<float>& b) {
	if (a.size() != b.size()) {
		cout << "Error: The two vectors must have the same length\n";
		return;
	}
	for (int i = 0; i < (int)a.size(); i++) {
		a[i] -= b[i];
	}
}
//This function multiplies the elements of vector a by the constant factor x, a0i = xai.
void vector_scale(std::vector<float>& a, float const b) {
	for (int i = 0; i < (int)a.size(); i++) {
		a[i] *= b;
	}
}

//This function adds the constant value x to the elements of the vector a, a0i = ai + x.
void vector_add_constant(std::vector<float>& a, float const b) {
	for (int i = 0; i < (int)a.size(); i++) {
		a[i] += b;
	}
}
//This function evaluate standard diviation supposing that E(a) = b.
float sd(std::vector<float>& a, float const b) {
	if (a.size() == 0)
		return NA;
	if (a.size() == 1)
		return b;
	double sum_sq = 0;
	for (int i=0; i< (int) a.size(); i++) {
		sum_sq+=(a[i]-b)*(a[i]-b);
		//cout << data[i] << "\n";
	}
	return (float)sqrt(sum_sq/(a.size()-1));

}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//function that finds the next breakpoint in the segment x:
void findnextbreakpoint( const std::vector<float>& x, int & besti, int & bestsign,float & bestlambda, int segLength, float a,float b) {

	if (segLength > 2) {
		vector <float> v (segLength);

		float x1 = x[0];
		float xn = x[segLength-1];
		for (int i=0; i<segLength; i++)
			v[i] = x1+i*(xn-x1)/(segLength-1);
		int jj;
		vector_sub(v,x);
		if (a==b) {
			if (a==1) {
				jj = get_max_index(v);
				bestlambda = v[jj]/2;
				bestsign = -1;
			} else {
				jj = get_min_index(v);
				bestlambda = -v[jj]/2;
				bestsign = 1;
			}
			besti = jj;
		} else {
			vector <float> f2 (segLength-2);
			for (int i=0;i<segLength-2;i++)
				f2[i] = (float)i+1;
			float myconst = (b-a)/(segLength-1);
			vector_scale(f2,myconst);
			vector_add_constant(f2,a);
			//I recoded the above three lines in another way:

			int counta = 0;
			//a strange way to initialize bestmax: it's just so that in the
			//following FOR loop, we have for sure that a1 > bestmax.
			float bestmax = v[1]/(f2[0]+1) - 1;
			int bestindex = 0; //but it will be initialized anyway
			for (int i=1; i<segLength-1; i++)  {
				counta++;
				float a1 = v[i]/(f2[i-1]+1);
				float a2 = v[i]/(f2[i-1]-1);
				if (a1 > bestmax || a2 > bestmax) {
					bestindex = counta;
					bestmax = a1*(a1>=a2) + a2*(a1 < a2);
					bestsign = -(a1 >= a2) + (a1 < a2);
				}
			}
			bestlambda = bestmax;
			besti = bestindex;
			f2.clear();
		}
		v.clear();
	} else {
		// no more breakpoint
		besti = -1;
		bestlambda = -1;
		bestsign = 0;
	}
}

static void calculateBreakpoints_perform(float** V, int** jump, float** J, int from, int to, int kk)
{
    //std::cout << "calculateBreakpoints_perform from=" << from << " to=" << to << " kk=" << kk << std::endl;
	float* mk = new float[kk];
	for (int ii = from; ii <= to; ii++) {
		float* pV = V[ii+1];
		int* pJ = jump[ii];
		for (int jj = ii+1; jj <= kk; jj++) {

		    int ind = 0;
			float min = MAXFLOAT;
			float* rr = mk;
			float* pp = &V[ii][ii];
			float** qq = &J[1+ii];
			for (int gg = ii; gg < jj; gg++) {
			  float rr_val = *pp++ + (*qq++)[jj];
			  if (rr_val < min) {
				min = rr_val;
				ind = gg-ii;
			  }
			  *rr++ = rr_val;
			}
			ind++;
			pV[jj] = mk[ind-1];
			pJ[jj] = ind+ii;
		}
	}
	delete [] mk;
}

int calculateBreakpoints_general(double threshold, int length, const vector<float>& ratio,vector<int>& bpfinal, int normal_length_dummy, int breakPointType, const std::string& chr) {

#ifdef PROFILE_TRACE
	time_t t0 = time(NULL);
#endif
    if (ratio.size()==0) {
        cout << "..You have zero windows with reads. Will try to continue anyway..\n" ;
        ThreadPoolManager::getInstance()->unlock();
        return -1;
    }
	const int maxValue = 3;
	float absoluteMax = *max_element(ratio.begin(),ratio.end());
//	const float minL = 0.5;

	//define a vector to stock centered ratio:
	vector <float> Y;
	//vector to define shift
	vector <int> shift (length);
	vector <int> shift_bp;
	shift[0] = 0;

	bool isLogRat = 0;

//	const float miniC = float(.001);
//	const float maxiC = float(0.001); //.05
	for (int i = 0; i < length; i++) {
		if (ratio[i]!=NA) {
			if (i!=0)
				shift[i]=shift[i-1];

			if (isLogRat) {
                if (ratio[i] < maxValue)
                    Y.push_back(log(ratio[i]));
                else
                    Y.push_back(log((ratio[i]-maxValue)/(absoluteMax-maxValue)+maxValue));
			} else {
                if (ratio[i] < maxValue)
                    Y.push_back(ratio[i]);
                else
                    Y.push_back((ratio[i]-maxValue)/(absoluteMax-maxValue)+maxValue);
			}
			shift_bp.push_back(shift[i]);
		} else {
			if (i!=0)
				shift[i]=shift[i-1]+1;
			else
				shift[i]=1;
		}
	}
	int n = Y.size();
	shift.clear();
//test - this block should be deleted.
	//for (int i = 2700; i < 3000; i++) {
	//		Y[i] = 1;
	//}

	//k = maximum number of breakpoints to find: (default:
	//20 percent of points for now) //CHANGED...
	int k	= (n/5)-2;
	if (k>2500) k=2500;
	//chromosome is too small to search breakpoints
	if (k<=0)
		return -1;

	//define a vector to stock the ordered list of found breakpoints:
	vector <int> bp (k);

    if (n==1) {
    	ThreadPoolManager::getInstance()->unlock();
        cerr << "You have 1 window with reads. It is not normal. Please check you parameters\nIf it happens for chrY, maybe your sample is female? Then use sex=XX\nWill try to continue anyway..\n" ;
        return -1;
    }
    if (n<5) {
    	ThreadPoolManager::getInstance()->unlock();
        cerr << "You have "<< n <<" windows with reads. It is not normal. Please check you parameters\nIf it happens for chrY, maybe your sample is female? Then use sex=XX\nWill try to continue anyway..\n" ;
        return -1;
    }
#if 0
    double adjust = 1;
    if (normal_length!=0) {
        adjust = n*1.0/normal_length;
       // cout << "here we will adjust by "<<adjust<<"\n";
    }
#endif
	//find lacking breakpoints
	vector <int> bp_vale;
	vector <float> diffs;

	float coeff = 0.8;
	if (n*(1-coeff)>10000) {
        coeff = 1- 10000.0/n;
        coeff = max(float(0.8),coeff);
	}
	for (int i = 4;i < n-4; i++)
		diffs.push_back(fabs(Y[i+1]-Y[i-1]));

	if (diffs.size()>10) {
        sort (diffs.begin(), diffs.end());
        float maxD = diffs[int(diffs.size()*coeff)-1]; // *0.95 for 5%
        if (maxD>0)
            for (int i = 10;i < n-10; i++)
                if (fabs(Y[i+1]-Y[i-1]) >= maxD && fabs(Y[i+2]-Y[i-2]) >= maxD && fabs(Y[i+3]-Y[i-3]) >= maxD&& fabs(Y[i+4]-Y[i-4]) >= maxD)
                    bp_vale.push_back(i);

	}
	diffs.clear();


   // if (maxD==0) {
    //    cout << "..All values seem to be the same..\n";
    //}


	// center the Y vector:
	float meanY = get_median(Y); //changed to median
	bool ifAllTheSame = true;

	for (int i = 0;i < n; i++) {
		Y[i] = Y[i]-meanY;
		if (Y[i]!=0)
            ifAllTheSame=false;
    }
    if (ifAllTheSame) {
        cout << "..all values are the same\n";
        ThreadPoolManager::getInstance()->unlock();
        return -1;
    }
	/* define vector to stock the cumsum: */
	vector <float> u (n);
	u[0] = Y[0];
	for (int i = 1;i < n; i++)
		u[i] = u[i-1]+Y[i];



	vector <float> c (n);
	vector <float> absC (n);
	c[0] = 0; absC[0] = 0;
	for (int i = 1; i < n; i++)  {
		c[i] = i*u[n-1]/n - u[i-1];
		absC[i] = (float) fabs((float)c[i]);
	}

	u.clear();

	//initialise the stack matrix:
	float ** stack;
	int row = 7;
	int col = 2*k+1;
	stack = new float * [row];
	for (int i = 0; i<row; i++)
		stack[i] = new float[col];


	//initialize two vectors that we will use to order the lambdas:
	vector <int> arrows (2*k+1);
	vector <float> lambdas (2*k+1);

	//find the first breakpoint:
	int nb=1;
	//int last=1; KEVIN, THIS VARIABLE IS NOT USED
	int l = get_max_index(absC)+1;

	bp[0]=l;
	stack[0][0] = 1;
	stack[1][0] = 0;
	stack[2][0] = (float)n;
	stack[3][0] = 0;
	stack[4][0] = (float)l;
	stack[5][0] = (float)1 - 2*(c[l-1] < 0);
	stack[6][0] = absC[l-1];

	absC.clear(); //don't need it any more

	//the main loop:

    int stackcounter = 0;
	int currentstackcolumn = 0;
	float bestlambda;
	int besti,bestsign;
	while (nb < k)  {

		// process the LEFT interval of the current breakpoint:
		int mystart = (int)stack[0][currentstackcolumn];
		int mystop  = (int)stack[4][currentstackcolumn];
		int segLength = mystop-mystart+1;
		vector <float> xx (segLength);

		int countz = 0;
		for (int i = mystart; i <= mystop; i++) {
			xx[countz] = c[i-1];
			countz++;
		}

		float a = stack[1][currentstackcolumn];
		float b = stack[5][currentstackcolumn];

		findnextbreakpoint( xx, besti,bestsign,bestlambda,segLength,a,b);

		//if the interval is non-empty:
		if (besti > 0)	{
			//tell the program to move across one row of stack:
			stackcounter++;
			stack[0][stackcounter] = stack[0][currentstackcolumn];
			stack[1][stackcounter] = stack[1][currentstackcolumn];
			stack[2][stackcounter] = stack[4][currentstackcolumn];
			stack[3][stackcounter] = stack[5][currentstackcolumn];
			stack[4][stackcounter] = besti + stack[0][currentstackcolumn];
			stack[5][stackcounter] = (float) bestsign;
			stack[6][stackcounter] = bestlambda;
			//what follows in the "if" loop is an inefficient way to put the first entry into the
			//vectors "arrows" and "lambdas":
			if (stackcounter == 1)	{
				arrows[0] = 1;
				lambdas[1] = bestlambda;
			} else {
				int currentarrow = 0;
				int previousarrow = 0;
				int ac = arrows[currentarrow];
				while (bestlambda <= lambdas [ac] && lambdas[ac] != 0 && arrows[ac] != 0 ) 	{
					previousarrow = ac;
					ac = arrows[ac];
					currentarrow = ac;
				}
				int a0 = arrows[0];
				// update the arrows vector to give the new total order of the lambdas:
				//... if it came in first:
				if (bestlambda > lambdas[a0])	{
					arrows[stackcounter] = a0;
					arrows[0] = stackcounter;
				}
				//...if it came in last:
				else if ( arrows[ac] == 0 )
					arrows[ac] = stackcounter;
				else {
				//.. or if it came in in the middle:
					arrows[previousarrow] = stackcounter;
					arrows[stackcounter] = currentarrow;
				}
				// update the lambdas vector:
				lambdas[stackcounter] = bestlambda;
			}
      	}
		xx.clear();
		//process the RIGHT interval of the current breakpoint:
		mystart = (int)stack[4][currentstackcolumn];
		mystop = (int)stack[2][currentstackcolumn];
		segLength = mystop-mystart+1;
		vector <float>  xxx (segLength);
		countz = 0;
		for (int i = mystart; i <= mystop; i++) {
			xxx[countz] = c[i-1];
			countz++;
		}
		a = stack[5][currentstackcolumn];
		b = stack[3][currentstackcolumn];

		findnextbreakpoint(xxx, besti,bestsign,bestlambda,segLength,a,b);
		//if the interval is non-empty:
        if (besti > 0) 	{
			//tell the program to move across one row of stack:
			stackcounter++;
			stack[0][stackcounter] = stack[4][currentstackcolumn];
			stack[1][stackcounter] = stack[5][currentstackcolumn];
			stack[2][stackcounter] = stack[2][currentstackcolumn];
			stack[3][stackcounter] = stack[3][currentstackcolumn];
			stack[4][stackcounter] = besti + stack[4][currentstackcolumn];
			stack[5][stackcounter] = (float) bestsign;
			stack[6][stackcounter] = bestlambda;
			//what follows in the "if" loop is an inefficient way to put the first entry into the
			//vectors "arrows" and "lambdas":
			if (stackcounter == 1) {
				arrows[0]= 1;
				lambdas[1] = bestlambda;
			}
			else    {
				int currentarrow = 0;
				int previousarrow = 0;
				int ac3 = arrows[currentarrow];
				while (bestlambda <= lambdas[ac3] && lambdas[ac3] != 0 && arrows[ac3] != 0 ) {
					previousarrow = ac3;
					ac3 = arrows[ac3];
					currentarrow = ac3;
				}
				int a03 = arrows[0];
				// update the arrows vector to give the new total order of the lambdas:

				//... if it came in first:
				if (bestlambda > lambdas[a03] ) {
					arrows[stackcounter] = a03;
					arrows[0] = stackcounter;
				}
				//...if it came in last:
				else if ( arrows[ac3] == 0 )
					arrows[ac3] = stackcounter;
				else {
				//.. or if it came in in the middle:
					arrows[previousarrow] = stackcounter;
					arrows[stackcounter] = currentarrow;
				}
				// update the lambdas vector:
				lambdas[stackcounter] = bestlambda;
			}
      	}
		xxx.clear();
		currentstackcolumn = arrows[0];
		//fix up arrows and lambdas:
		arrows[0] = arrows[currentstackcolumn];
		//print the next breakpoint:
		//cout << stack[4][currentstackcolumn] << "\n";
		bp[nb] = stack[4][currentstackcolumn];   //'='?�: conversion de 'float' en 'int', perte possible de donn?�es
		nb++;
    }

	arrows.clear();
	lambdas.clear();
	c.clear();

	for (int i = 0; i<row; i++)
		delete [] stack[i];
	delete [] stack;

	//%%%%%%%% dynamic programming step %%%%%%%%
	//REMARK: to try and sidestep the nightmare of the difference between matlab
	//and "1" indexing and C with its "0" indexing, I'm here going to pretend
	//i'm in matlab by just defining matrices with one extra row and column,
	//thus allowing retention of matlab indexing.



	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//sort the elements of bp:

	//bp_vale
	for (int i = 0; i<(int)bp_vale.size(); i++)
		bp.push_back(bp_vale[i]);

	//cout << bp.size() << " vs k= " <<k <<"\n";
	sort (bp.begin(), bp.end());
	for (int i = bp.size()-1; i>0; i--)
		if (bp[i] == bp[i-1])
			bp.erase(bp.begin()+i);

	//fix "k" just so that it matches up with the "k" in dp.m:
	k = bp.size()+1;

//	printf("%d\n",k);

	vector <int> bb (k+2);
	bb[1] = 1;
	bb[k+1] = n+1;

	for (int i = 2; i<= k; i++) {
		bb[i] = bp[i-2];
		//cout << bb[i] << "\n";
	}

	float** J ;
	row = k+1;
	col = k+1;
	J = new float * [row];
	for (int i = 0; i<row; i++)
		J[i] = new float[col];


	vector <float> s (n+2);
	s[1] = 0;

	for (int i=2; i <n+2; i++)  {
		s[i] = s[i-1]+Y[i-2];
		//printf("%f\n",gsl_vector_get(s,i));
	}

	vector <float> v (n+2);
	v[1] = 0;

	for (int i=2; i<n+2; i++)
		v[i] = v[i-1] + Y[i-2]*Y[i-2];

	for (int i=1; i <= k; i++) {
		for (int j=i; j<=k; j++) {
			int Istart =  bb[i];
			int Iend = bb[j+1] - 1;

			float st = s[Iend+1] - s[Istart];
			st = st*st;

			int su = Iend - Istart + 1;
			st = st/su;
			st = v[Iend+1] - v[Istart] - st;
			J[i][j] = st;
		}
	}

	bb.clear();
	v.clear();
	s.clear();

	int ** jump; //k+1,k+1
	jump = new int * [row];
	for (int i = 0; i<row; i++)
		jump[i] = new int[col];

	float ** V;
	row = k+2; col = k+1;
	V = new float * [row];
	for (int i = 0; i<row; i++)
		V[i] = new float[col];

	for (int i=1; i<=k; i++)
		V[1][i] = J[1][i];

#ifdef PROFILE_TRACE
	std::cout << "PROFILING [tid=" << pthread_self() << "]: breakpoint prologue computed in " << (time(NULL)-t0) << " seconds\n" << std::flush;
#endif

#ifdef PROFILE_TRACE
	t0 = time(NULL);
#endif


	calculateBreakpoints_perform(V, jump, J, 1, k-1, k);

#ifdef PROFILE_TRACE
	std::cout << "PROFILING [tid=" << pthread_self() << "]: main breakpoint loop computed in " << (time(NULL)-t0) << " seconds\n" << std::flush;
#endif

	for (int i = 0; i<k+1; i++)
		delete[] J[i];
	delete[] J;

#ifdef PROFILE_TRACE
	t0 = time(NULL);
#endif
	vector <float> Ltemp (k+1);
	vector <float> L (k+1);
	for (int i=1; i<=k; i++) {
		Ltemp[i] = V[i][k];
		//printf("%f\n",Ltemp[i]);
	}

	for (int i = 0; i<k+2; i++)
		delete[] V[i];
	delete[] V;



	// start point of modifications 20/09/2011 Kevin/Valentina

	//normalise:
	//first take the log:
	vector <float> logL (k+1);
	//cout << "log(L) values :";
	for (int i=1; i<=k; i++) {
		logL[i] = log(Ltemp[i]);
		if (i < 1) {
		//cout << logL[i] << "\n";
		}
		}
	float L0 = logL[k] - logL[1];
	for (int i=1; i<=k; i++) {
		L[i] = ((logL[k] - logL[i])/L0) * (k) + 1;
		}


	//normalise:
	//float L0 = Ltemp[1]; //*adjust to keep less breakpoints for short chromosomes
	//for (int i=1; i<=k; i++) {
	//	L[i] = Ltemp[i]/L0*adjust;
		//printf("%f\n",L[i]);
	//}


	//find the discrete slope of the slope:
	vector <float> gradL (k);
	vector <float> gradgradL (k-1);
	//calculate gradL:
	for (int i=2; i<=k; i++) {
		gradL[i-1] = L[i] - L[i-1];
		}
	//calculate gradgradL:
	//cout << "2nd derivative values :";
	for (int i=2; i<=k-1; i++) {
		gradgradL[i-1] = gradL[i] - gradL[i-1];
		if (i < 300) {
		//cout << gradgradL[i-1] << "\n";
		}
		}

	//find the largest index at which the gradgradL is larger than some pre-defined cutoff:
	float imax = 0;
	for (int i=1; i<=k-2; i++) {
		if (gradgradL[i] > threshold)
		{
		imax = i;
		}
		}
	//index correction:
	imax = imax + 1;

	//output
	//cout << L[imax] << "\n";

	// keep rk from old code, to keep it simple:
	int rk = imax;

	//int rk = 1;
	//find the estimated slope at the first point:
	//float myGrad = L[rk+1] - L[rk];
	//printf("%f\n",myGrad);

	//keep going until either we get to the end of L, or we get a slope greater than
	// the cutoff "threshold":
	//while (myGrad <= threshold && rk < k) {
	//	rk++;
	//	myGrad = L[rk+1] - L[rk-1];
	//	myGrad = myGrad/2;
		//if (L[rk-1] > minL) //continue
			//myGrad = threshold-1;
		//printf("%f %f\n",myGrad, L[rk-1]);
	//}
	//rk = rk-1;


	//cout << L[rk] << "\n";

	//printf("%d\n",rk);

	// end point of modifications 20/09/2011 Kevin/Valentina


	// a vector for final breakpoints
	bpfinal.clear();
	bpfinal  = vector <int>(rk);

	vector <int> rjumps (rk);

	if (rk) {
		rjumps[rk-1] = jump[rk][k];

		for (int i = rk-1; i>0; i--)
			rjumps[i-1] = jump[i][rjumps[i]-1];


	}

	//for (int i=0;i<rk;i++)
	//{
	//   //printf("%f\n",rjumps[i]);
	//}

	//convert back to the exact, final, chosen breakpoints:
	for (int i=0; i<rk; i++) {
	   //check if the temporary -2 and -2 in this line is actually correct!
	   bpfinal[i] = bp[rjumps[i]-2]-2;
	}

	for (int i = 0; i<k+1; i++)
		delete[] jump[i];
	delete[] jump;

	Y.clear();
	Ltemp.clear();
	L.clear();
	rjumps.clear();
	//end !!!

	ThreadPoolManager::getInstance()->lock();
	cout << "Chromosome: " << chr << "\n";
	cout << "Total windows:" << length << "\n";
	cout << "Not NA windows:" << n << "\n";
	cout << "Final breakpoints:\n";

	//recalculate breakpoints: with shift
	if (breakPointType==SIMPLERIGHT ) {
        for (int i = 0; i<(int)bpfinal.size(); i++) {
            bpfinal[i] += shift_bp[bpfinal[i]];
            cout << bpfinal[i] << "\t";
        }
    }
    else if (breakPointType==LARGECLOSE || breakPointType==NORMALLEVEL || breakPointType==NOCALL || breakPointType==HALFLENGTH) { //we should split a breakpoint into two if there is a difference in "shift_bp"
        vector <int> bpWithAddedPoints;
        //if there are '-1' in the beginning, create a fragment
        if ((int)bpfinal.size()>0 && shift_bp[0]>0) {
            bpWithAddedPoints.push_back(shift_bp[0]-1);
            cout << bpWithAddedPoints.back() << "\t";
        }
        for (int i = 0; i<(int)bpfinal.size(); i++) {
            if (bpfinal[i]!=-1){
                bpWithAddedPoints.push_back(shift_bp[bpfinal[i]]+bpfinal[i]);
                cout << bpWithAddedPoints.back() << "\t";
                if (shift_bp[bpfinal[i]]!=shift_bp[bpfinal[i]+1]) { //add a second breakpoint
                    bpWithAddedPoints.push_back(shift_bp[bpfinal[i]+1]+bpfinal[i]);
                    cout << bpWithAddedPoints.back() << "\t";
                }
            }
        }
        //if there are '-1' in the end, create a fragment
        if ((int)bpfinal.size()>0 && shift_bp[n-1]+n!=length) {
            bpWithAddedPoints.push_back(shift_bp[n-1]+n-1);
            cout << bpWithAddedPoints.back() << "\t";
        }
        //sort(bpWithAddedPoints.begin(), bpWithAddedPoints.end()); //they should be already sorted...
        bpfinal=bpWithAddedPoints;
        bpWithAddedPoints.clear();
    }

	cout << "\n";

#ifdef PROFILE_TRACE
	std::cout << "PROFILING [tid=" << pthread_self() << "]: breakpoint epilogue computed in " << (time(NULL)-t0) << " seconds\n" << std::flush;
#endif
	ThreadPoolManager::getInstance()->unlock();

	shift_bp.clear();
    return n;
}

//find closest value to k/ploidy for medianValue
float round_by_ploidy(float value, int ploidy) {
	if (value >0)
		return float(int(value*ploidy + 0.5))/ploidy;
	return float(int(value*ploidy))/ploidy;
}

int round_f(float r) {
    return (r > 0.0) ? floor(r + 0.5) : ceil(r - 0.5);
}
//int round(double r) {
//    return (r > 0.0) ? floor(r + 0.5) : ceil(r - 0.5);
//}


string pathAppend(const string& p1, const string& p2) {

   char sep = '/';
   string tmp = p1;
   char sep2=sep;
#if defined(_WIN32)
  sep2 = '\\';
#endif
  char lastSymb=  p1[p1.length( )-1];
  if (lastSymb != sep && lastSymb != sep2) { // Need to add a
     tmp += sep;                // path separator
     return(tmp + p2);
  }
  else
     return(p1 + p2);
}

int isCG (const char & a) {
	if ((a == 'C')||(a == 'c')||(a == 'G')||(a == 'g'))
		return 1;
	/*if ((a == 'N')||(a == 'n'))
		return NA;*/
	return 0;
}
int isN (const char & a) {
	if ((a == 'N')||(a == 'n'))
		return 1;
	return 0;
}
int isSpaceCharacter (const char & a)  {
	if ((a == '\n')||(a == ' ')||(a == '\t'))
		return 1;
	return 0;
}
float polynomial(const float x, const double a, const double b, const double c){ //ax^2+bx+c
	return float(a*x*x+b*x+c);
}
float polynomial(const float x, const double * a, double ratio, int degree){
	double y = a[degree];
	for (int i = 0; i < degree; i++) {
		y += pow(x,degree-i)*a[i];
	}
	return float(y*ratio);
}

float polynomial(const float x, const double a0, const double a1, const double a2, const double a3){ //a0x^3+a1x^2+a2x+a3
	return float(a0*x*x*x+a1*x*x+a2*x+a3);
}
float runEM_linear (const vector<float>& x,const vector<float>& y,double & a0,double & a1, int maximalNumberOfIterations,int ploidy, int maximalNumberOfCopies) {

	float rmserror = -1;

	if (x.size() != y.size()) {
		cerr << "Error: object size is different";
		return 0;
	}
	vector <int> cluster (x.size());
	vector <float> res (maximalNumberOfCopies+1);

	int count = 0;

	for (count = 0; count < maximalNumberOfIterations; count++) {

		for (int i = 0; i <(int)x.size(); i++) {
			for (int j = 0; j <= maximalNumberOfCopies; j++)
				res[j] = fabs(j*a0/ploidy*x[i]+ j*a1/ploidy -y[i]);
			cluster[i] = get_min_index(res);
		}
		int npoints = 0;
		for (int i = 0; i <(int)x.size(); i++)
			if (cluster[i] == ploidy)
				npoints++;
		int nvars = 1; //fit by linear function
		ap::real_2d_array xy;
		xy.setlength(npoints,nvars+1);
		int pos = 0;
		for (int i = 0; i <(int)x.size(); i++) {
			if (cluster[i] == ploidy){
				xy(pos,0) = x[i];
				xy(pos,1) = y[i];
				pos++;
			}
		}

		linearmodel lm;
		int info;
		lrreport ar;
		lrbuildz(xy,npoints,nvars,info,lm,ar);
		if (info != 1) {
			cerr << "Error in linear regression, code: " << info <<"\n";
			break;
		}
		rmserror = float(ar.rmserror);
		ap::real_1d_array v;
		v.setlength(nvars+1);
		lrunpack(lm,v,nvars);
		//cout << v(0);  //x
		//cout << v(1);  //x^3
		//cout << v(2);  //x^2
		//cout << v(3);  //intercept

		cout << a0 - v(0) << "\t"; //x
		//cout << a1 - v(1) << "\n"; //intercept
		//if (fabs(a0 - v(0)) + fabs(a1 - v(1)) == 0) {
		if (fabs(a0 - v(0)) == 0) {
			break;
		}
		a0 = v(0);
		//a1 = v(1);
	}
	cluster.clear();
	res.clear();
	cout << "Number of EM iterations :" << count << "\n";
	return rmserror;
}

float runEM (const vector<float>& x,const vector<float>& y,double & a0,double & a1, double & a2, double & a3, int maximalNumberOfIterations,int ploidy, int maximalNumberOfCopies) {

	float rmserror = -1;

	if (x.size() != y.size()) {
		cerr << "Error: object size is different";
		return 0;
	}
	vector <int> cluster (x.size());
	vector <float> res (maximalNumberOfCopies+1);

	int count = 0;

	for (count = 0; count < maximalNumberOfIterations; count++) {

		for (int i = 0; i <(int)x.size(); i++) {
			for (int j = 0; j <= maximalNumberOfCopies; j++)
				res[j] = fabs(polynomial(x[i],j*a0/ploidy,j*a1/ploidy,j*a2/ploidy,j*a3/ploidy)-y[i]);
			cluster[i] = get_min_index(res);
		}
		int npoints = 0;
		for (int i = 0; i <(int)x.size(); i++)
			if (cluster[i] == ploidy)
				npoints++;
		int nvars = 3; //fit by cubic polynomial
		ap::real_2d_array xy;
		xy.setlength(npoints,nvars+1);
		int pos = 0;
		for (int i = 0; i <(int)x.size(); i++) {
			if (cluster[i] == ploidy){
				xy(pos,2) = x[i];
				xy(pos,0) = x[i]*x[i]*x[i];
				xy(pos,1) = x[i]*x[i];
				xy(pos,3) = y[i];
				pos++;
			}
		}

		linearmodel lm;
		int info;
		lrreport ar;
		lrbuild(xy,npoints,nvars,info,lm,ar);
		if (info != 1) {
			cerr << "Error in linear regression, code: " << info <<"\n";
			break;
		}
		rmserror = float(ar.rmserror);
		ap::real_1d_array v;
		v.setlength(nvars+1);
		lrunpack(lm,v,nvars);
		//cout << v(0);  //x
		//cout << v(1);  //x^3
		//cout << v(2);  //x^2
		//cout << v(3);  //intercept

		cout << a0 - v(0) << "\t"; //x^3
		cout << a1 - v(1) << "\t"; //x^2
		cout << a2 - v(2) << "\t"; //x
		cout << a3 - v(3) << "\n"; //intercept
		if (fabs(a0 - v(0)) + fabs(a1 - v(1)) + fabs(a2 - v(2)) +fabs(a3 - v(3)) == 0) {
			break;
		}
		a0 = v(0);
		a1 = v(1);
		a2 = v(2);
		a3 = v(3);
	}
	cluster.clear();
	res.clear();
	cout << "Number of EM iterations :" << count << "\n";
	return rmserror;
}
std::string stringFromBool (bool value) {
    if (value) return "True";
    return "False";
}

float runEM (const vector<float>& x,const vector<float>& y,double * a, int degree, int & NumberOfIterations,int ploidy, int maximalNumberOfCopies, bool intercept, float contamination) {
    if (contamination==0) contamination=0.3; //starting from v11.2 - to improve the fit
	float rmserror = -1;

	if (x.size() != y.size()) {
		cerr << "Error: object size is different\n";
		return 0;
	}
	vector <int> cluster (x.size());
	vector <float> res (maximalNumberOfCopies+1);

	int count = 0;

	for (count = 0; count < NumberOfIterations; count++) {

		for (int i = 0; i <(int)x.size(); i++) {
			for (int j = 0; j <= maximalNumberOfCopies; j++)
				res[j] = fabs(polynomial(x[i],a, ( float(j)*(1-contamination)+2* contamination ) /  (ploidy*(1-contamination)+2*contamination), degree)-y[i]);
			cluster[i] = get_min_index(res);
		}
		int npoints = 0;
		for (int i = 0; i <(int)x.size(); i++)
			if (cluster[i] == ploidy)
				npoints++;
		if (npoints==0) {
            cerr << "ERROR: there was a problem in the initial guess of the polynomial. Please contact the support team of change your input parameters. Exit.\n";
            exit(-1);
		}
		int nvars = degree; //3 if fit by cubic polynomial
		ap::real_2d_array xy;
		xy.setlength(npoints,nvars+1);
		int pos = 0;
		for (int i = 0; i <(int)x.size(); i++) {
			if (cluster[i] == ploidy){
				xy(pos,degree) = y[i];
				xy(pos,degree-1) = x[i];;
				for (int j = degree-2; j>=0; j--) {
					xy(pos,j)=xy(pos,j+1)*x[i];
				}
				//xy(pos,0) = x[i]*x[i]*x[i];
				//xy(pos,1) = x[i]*x[i];
				//xy(pos,2) = x[i];
				//xy(pos,3) = y[i];
				pos++;
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
			cerr << "Error in linear regression, code: " << info <<"\n";
			break;
		}
		rmserror = float(ar.rmserror);
		ap::real_1d_array v;
		v.setlength(nvars+1);
		lrunpack(lm,v,nvars);
		double changeInValue = 0;
		//cout << v(0);  //x
		//cout << v(1);  //x^3
		//cout << v(2);  //x^2
		//cout << v(3);  //intercept
		for (int i = 0; i <degree; i++) {
			cout << a[i] - v(i) << "\t";
			changeInValue += fabs(a[i] - v(i));
		}
		if (intercept) {
			cout << a[degree] - v(degree) << "\n"; //intercept
			changeInValue += fabs(a[degree] - v(degree));
		}
		else
			cout << "\n";

		if (changeInValue == 0) {
			break;
		}
		for (int i = 0; i <degree; i++) {
			a[i] = v(i);
		}
		if (intercept)
			a[degree] = v(degree); //intercept
	}
	cluster.clear();
	res.clear();
	cout << "Number of EM iterations :" << count << "\n";
	NumberOfIterations=count;
	return rmserror;
}


float runEMlog (const vector<float>& x,const vector<float>& y,double * a, int degree, int & realNumberOfIterations, int maximalNumberOfIterations,int ploidy, int maximalNumberOfCopies, bool intercept, float contamination) {

	float rmserror = -1;

	if (x.size() != y.size()) {
		cerr << "Error: object size is different";
		return 0;
	}
	vector <int> cluster (x.size());
	vector <float> res (maximalNumberOfCopies+1);

	int count = 0;

	for (count = 0; count < maximalNumberOfIterations; count++) {

		for (int i = 0; i <(int)x.size(); i++) {
			for (int j = 0; j <= maximalNumberOfCopies; j++) {
                float conLevel = log(1./10);
                if (j>=0) conLevel = log(float(j)/ploidy);
                res[j] = fabs(polynomial(x[i],a,1,degree)-y[i]+ conLevel);
			}
			cluster[i] = get_min_index(res);
		}
		int npoints = 0;
		for (int i = 0; i <(int)x.size(); i++)
			if (cluster[i] == ploidy)
				npoints++;
		int nvars = degree; //3 if fit by cubic polynomial
		ap::real_2d_array xy;
		xy.setlength(npoints,nvars+1);
		int pos = 0;
		for (int i = 0; i <(int)x.size(); i++) {
			if (cluster[i] == ploidy){
				xy(pos,degree) = y[i];
				xy(pos,degree-1) = x[i];;
				for (int j = degree-2; j>=0; j--) {
					xy(pos,j)=xy(pos,j+1)*x[i];
				}
				//xy(pos,0) = x[i]*x[i]*x[i];
				//xy(pos,1) = x[i]*x[i];
				//xy(pos,2) = x[i];
				//xy(pos,3) = y[i];
				pos++;
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
			cerr << "Error in linear regression, code: " << info <<"\n";
			break;
		}
		rmserror = float(ar.rmserror);
		ap::real_1d_array v;
		v.setlength(nvars+1);
		lrunpack(lm,v,nvars);
		double changeInValue = 0;
		//cout << v(0);  //x
		//cout << v(1);  //x^3
		//cout << v(2);  //x^2
		//cout << v(3);  //intercept
		for (int i = 0; i <degree; i++) {
			cout << a[i] - v(i) << "\t";
			changeInValue += fabs(a[i] - v(i));
		}
		if (intercept) {
			cout << a[degree] - v(degree) << "\n"; //intercept
			changeInValue += fabs(a[degree] - v(degree));
		}
		else
			cout << "\n";

		if (changeInValue == 0) {
			break;
		}
		for (int i = 0; i <degree; i++) {
			a[i] = v(i);
		}
		if (intercept)
			a[degree] = v(degree); //intercept
	}
	cluster.clear();
	res.clear();
	cout << "Number of EM iterations :" << count << "\n";
	realNumberOfIterations=count;
	return rmserror;
}

void processChrName(string & chr) {
	string::size_type pos = 0;
	if ( ( pos = chr.find("chr", pos)) != string::npos ) {
		chr.replace( pos, 3, "" );
		return ;
	} else {
		pos = 0;
		if ((pos = chr.find("NC_000", pos))!= string::npos) {
			//do nothing
		} else
		if ((pos = chr.find("NC_", pos))!= string::npos) {
			int $len = chr.length();
			chr = chr.substr($len-2);
			if (chr.compare("23")==0)
				chr = "X";
			else if (chr.compare("24")==0)
				chr = "Y";
			else if (chr.compare("25")==0)
				chr = "M";
			else if (chr[0]=='0' || chr[0]=='.')
				chr = chr[1];
		}
	}
}


size_t strccnt(const char *s, int c)
{
const unsigned char *us = (const unsigned char *) s;
const unsigned char uc = c;
size_t n = 0;
if (!uc) return 1;
while (*us) if (*us++ == uc) n++;
return n;
}


void strkeepOnly(char *s, const char *c) {
    string newString;
    for (unsigned int i = 0; i< strlen(s);i++) {
        for (unsigned int j = 0; j< strlen(c);j++){
            if (s[i]==c[j])
                newString+=s[i];
                break;
        }
    }
    strcpy(s, newString.c_str());
    return;
}
void strkeepOnly(string & s, const char *c) {
    string newString;
    for (unsigned int i = 0; i< s.length(); i++) {
        for (unsigned int j = 0; j< strlen(c);j++){
            if (s[i]==c[j]) {
                newString+=s[i];
                break;
            }
        }
    }
    s=newString;
    return;
}

std::string int2string (int a) {
    stringstream ss;
    ss << a;
    string str = ss.str();
    return str;
}

void deleteChar(std::string & s, char c, int moreLettersToDelete) {
    string newString;
    for (unsigned int i = 0; i< s.length(); i++) {
        if (s[i]!=c) {
                newString+=s[i];
        } else {
            i+=moreLettersToDelete;
        }
    }
    s=newString;
    return;
}


void deleteChar(std::string & s, char c) {
    string newString;
    for (unsigned int i = 0; i< s.length(); i++) {
        if (s[i]!=c) {
                newString+=s[i];
        }
    }
    s=newString;
    return;
}


void chomp (char* s) {
    int end = strlen(s) - 1;
    if (end >= 0 && s[end] == '\n')
      s[end] = '\0';
}

void chomp (string & s) {
    if (s.at(s.length() - 1) == '\n')
        s.erase (s.end()-1);
}

void filterWithQualities(string & pileupShort,string & qualityS, int minimalQualityPerPosition) {
    string newString;
    for (unsigned int i = 0; i< pileupShort.length(); i++) {
        int q = int (qualityS[i]);
        if (q>=minimalQualityPerPosition)
            newString+=pileupShort[i];
    }
    pileupShort=newString;
    return;

}

char complement(char n) {
  switch(n)  {
    case 'A': return 'T';
    case 'T': return 'A';
    case 'C': return 'G';
    case 'G': return 'C';
    case 'a': return 't';
    case 't': return 'a';
    case 'c': return 'g';
    case 'g': return 'c';
    case 'N': return 'N';

  }
  return 'N';
}

vector<int> merge_no_dups(const vector<int>& v1, const vector<int>& v2);


void getBAFinfo(std::string BAFValuesInTheSegment,float copyNumber,float &estimatedBAF,float &fittedBAF,std::string &medianBAFSym,
float & uncertainty, float normalContamination,int ploidy, bool noisyData, bool ifHomoz, bool CompleteGenomicsData) {

    bool fixedMu = true;
    vector <int> copyNumbers;

    if (copyNumber==NA) {
        estimatedBAF = NA;
        fittedBAF=NA;
        medianBAFSym = "-";
        uncertainty = NA;
        return;
    }

    //getCopyNumbers(copyNumber,copyNumbers);
    getCopyNumbers(copyNumber,copyNumbers,ploidy, noisyData);

    int preferedCN=round_f(copyNumber);
    int winningCN=NA;

    if (ifHomoz) {
        //align Homozygous status to this fragment

        cout << "..too few points to fit the data, average copy number="<<copyNumber<<"\n";
        cout << "suggest homozygous region\n";
        uncertainty = NA;
        estimatedBAF = 0;
        fittedBAF=NA;
        medianBAFSym="";
        if (round_f(copyNumber)>=1) {
            medianBAFSym="A";
            for (int i = 1; i<round_f(copyNumber); i++) {
                    medianBAFSym+="A";
            }
        }

        copyNumbers.clear();
        return;
    }

    double maxLogLikelyHood = -INFINITY;

    vector <float> BAFs;
    vector<string>heteroValuesPerWindowStrings = split(BAFValuesInTheSegment, ';');
    if (heteroValuesPerWindowStrings.size()>0) {
            for (int unsigned j = 0; j < heteroValuesPerWindowStrings.size(); j++) {
                stringstream ss;
                float f;
                ss << heteroValuesPerWindowStrings[j];
                ss >> f;
                BAFs.push_back(f);
            }
    } else {
        estimatedBAF = NA;
        fittedBAF = NA;
        uncertainty = NA;
        medianBAFSym="";
        if (round_f(copyNumber)==1){medianBAFSym="A";}
        if (round_f(copyNumber)>1) {
            medianBAFSym="-";
        }
        return;
    }
    heteroValuesPerWindowStrings.clear();

    uncertainty = NA;
    vector <double> LogLikelyHoods;
    vector <int> testedCN;
    vector <string> medianBAFSyms;
    vector <float>estimatedBAFs;
    vector <float>fittedBAFs;


    for (int unsigned i=0; i<copyNumbers.size(); i++) {
        int myCopyNumber = copyNumbers.at(i);
        if(myCopyNumber==0) {
            medianBAFSym = "";
            estimatedBAF = NA;
            fittedBAF = NA;
            uncertainty = NA;
        } else if (myCopyNumber==1) {
            medianBAFSym = "A";
            estimatedBAF = 0;
            fittedBAF = NA;
            uncertainty = NA;
        }
        else {
            int Bcount = 0;
            int Acount = myCopyNumber;
            while (Bcount<=Acount) {
                vector <float> mu;
                float mu_=HOMOZYG_MEAN;
                mu.push_back(mu_);
                mu.push_back(1-mu_);
                if (normalContamination && Bcount==0) {
                    mu_=max(mu_,normalContamination/(myCopyNumber*(1-normalContamination)+2*normalContamination));
                    if (mu_!=mu[0]) {
                        mu.push_back(mu_);
                        mu.push_back(1-mu_);
                    }
                }

                if (Bcount==Acount) {
                    mu_=0.5;
                    mu.push_back(mu_);
                } else if (Bcount!=0){
                    mu_ =(Bcount*(1.-normalContamination)+normalContamination) /(myCopyNumber*(1-normalContamination)+2*normalContamination);
                    mu.push_back(mu_);
                    mu.push_back(1-mu_);
                }
                double LogLikelyHood = calculateLogLikelyHoodNormalMixtureForBAFs(BAFs,mu,middleComponentMinWeight, fixedMu,CompleteGenomicsData);
                LogLikelyHoods.push_back(LogLikelyHood);
                testedCN.push_back(myCopyNumber);

                string medianBAFSymX = "A";
                for (int i = 1; i<Acount; i++) {
                    medianBAFSymX+="A";
                }
                for (int i = 0; i<Bcount; i++) {
                    medianBAFSymX+="B";
                }
                medianBAFSyms.push_back(medianBAFSymX);

                estimatedBAFs.push_back(Bcount*1./myCopyNumber);
                fittedBAFs.push_back(mu_);

                if (maxLogLikelyHood<LogLikelyHood) {
                    maxLogLikelyHood=LogLikelyHood;
                    medianBAFSym = medianBAFSymX;
                    estimatedBAF = Bcount*1./myCopyNumber;
                    fittedBAF=mu_;
                    winningCN=myCopyNumber;
                }



                mu.clear();
                Bcount++;
                Acount--;
            }
        }
    }
    if (LogLikelyHoods.size()>1) {
        vector <double> LogLikelyHoods_sec=LogLikelyHoods;
        sort (LogLikelyHoods.begin(), LogLikelyHoods.end());
        uncertainty = 1./(LogLikelyHoods[LogLikelyHoods.size()-1]-LogLikelyHoods[LogLikelyHoods.size()-2]);
   //     uncertainty = exp(LogLikelyHoods[LogLikelyHoods.size()-2]-LogLikelyHoods[LogLikelyHoods.size()-1]);//since v9.4 : p(secondBest)/p(best)
        if (copyNumbers.size()>=1) {
            if (copyNumbers[0]==1 && copyNumbers[1]==2 && copyNumber<1.5 && (medianBAFSym.compare("AA")==0 || (medianBAFSym.compare("AB")==0 && uncertainty >0.1))) {
                medianBAFSym="A";
                estimatedBAF = 0;
                fittedBAF=NA;
            }
            if (copyNumbers[0]>1 && LogLikelyHoods.size()>=3) { //should be always >=3, but just in case
                if (medianBAFSym[medianBAFSym.size()-1] == 'A') { //means AAAAA
                    //one should choose between AAAA and AAA
                    medianBAFSym = "A";
                    for (int i = 1; i<copyNumbers[0]; i++) {
                        medianBAFSym+="A";
                    }
                    if (copyNumber-copyNumbers[0]>=0.5)
                        medianBAFSym+="A";
                    //end recalculate uncertainty:
                   // cout << LogLikelyHoods[LogLikelyHoods.size()-1]<< " " << LogLikelyHoods[LogLikelyHoods.size()-2] <<"\n";
                    uncertainty = 1./(LogLikelyHoods[LogLikelyHoods.size()-1]-LogLikelyHoods[LogLikelyHoods.size()-3]);
                }
            }
            if (copyNumbers[0]>1 && winningCN!=preferedCN) { //check uncertainty for ambigous cases
                float preferedCNLogLH=-INFINITY;
                int bestInd=NA;
                for(unsigned int jj=0; jj<testedCN.size();jj++) {
                    if (testedCN[jj]==preferedCN && preferedCNLogLH<LogLikelyHoods_sec[jj]) {
                        preferedCNLogLH=LogLikelyHoods_sec[jj];
                        bestInd=jj;
                    }
                }
                float uncertWinVsPref=1./(LogLikelyHoods[LogLikelyHoods.size()-1]-preferedCNLogLH); //since v9.4
                if (uncertWinVsPref>0.1) {
                    uncertainty=uncertWinVsPref;
                    medianBAFSym=medianBAFSyms[bestInd];
                    estimatedBAF=estimatedBAFs[bestInd];
                    fittedBAF=fittedBAFs[bestInd];
                }
            }
        }
        LogLikelyHoods_sec.clear();
    }
    if (uncertainty>1)
        uncertainty = 100;
    else if (uncertainty>0)
        uncertainty *= 100;

    LogLikelyHoods.clear();
    BAFs.clear();
    copyNumbers.clear();
    testedCN.clear();
    medianBAFSyms.clear();
    estimatedBAFs.clear();
    fittedBAFs.clear();
}


void getCopyNumbers (float copyNumber, vector <int> & copyNumbers) {
    int lowCopy = int(floor(copyNumber));
    int highCopy = int(ceil(copyNumber));
    float minDiff = 0.45; minDiff = 0.23; //TODO : external parameter!!
    if (highCopy==1) {
        copyNumbers.push_back(round_f(copyNumber));
        return;
    }
    if (copyNumber-lowCopy<minDiff)
        copyNumbers.push_back(lowCopy);
    else if (highCopy-copyNumber<minDiff)
        copyNumbers.push_back(highCopy);
    else {
        copyNumbers.push_back(lowCopy);
        copyNumbers.push_back(highCopy);
    }
}

void getCopyNumbers (float copyNumber, vector <int> & copyNumbers,int ploidy, bool noisyData) {
    int lowCopy = int(floor(copyNumber));
    int highCopy = int(ceil(copyNumber));

    if (highCopy==1) {
            copyNumbers.push_back(round_f(copyNumber));
            return;

    }
    if(highCopy==lowCopy){ //not possible!!
            copyNumbers.push_back(round_f(copyNumber));
            return;

    }
    if (!noisyData) {
        float minDiff = 0.45;
        if (copyNumber-lowCopy<minDiff)
            copyNumbers.push_back(lowCopy);
        else if (highCopy-copyNumber<minDiff)
            copyNumbers.push_back(highCopy);
        else {
            copyNumbers.push_back(lowCopy);
            copyNumbers.push_back(highCopy);
        }
    } else { //prioretize normal copy number status
        float minDiffNoNorm = 0.35;
        float minDiff = 0.15;

        if (lowCopy == ploidy) {

            if (copyNumber-lowCopy<minDiffNoNorm)
                copyNumbers.push_back(lowCopy);
            else if (highCopy-copyNumber<minDiff)
                copyNumbers.push_back(highCopy);
            else {
                copyNumbers.push_back(lowCopy);
                copyNumbers.push_back(highCopy);
            }
        } else if (highCopy == ploidy){
            if (copyNumber-lowCopy<minDiff)
                copyNumbers.push_back(lowCopy);
            else if (highCopy-copyNumber<minDiffNoNorm)
                copyNumbers.push_back(highCopy);
            else {
                copyNumbers.push_back(lowCopy);
                copyNumbers.push_back(highCopy);
            }
        } else {

            if (copyNumber-lowCopy<minDiff)
                copyNumbers.push_back(lowCopy);
            else if (highCopy-copyNumber<minDiff)
                copyNumbers.push_back(highCopy);
            else {
                copyNumbers.push_back(lowCopy);
                copyNumbers.push_back(highCopy);
            }
        }
    }


}

double calculateLogLikelyHoodNormalMixtureForBAFs(vector <float> X,vector <float> mu,float middleComponentMinW, bool isMuFixed, bool CompleteGenomicsData) {
    double LogLikelyHood_i = -INFINITY ;
    int numberOfStates = mu.size();
    int N = X.size();
    int iterationCount = 0;
    int maxIterationCount = 10000;

    float maxSigma = 0.07; //used only if isMuFixed==TRUE
    float minOmega = 0.6;  //used only if isMuFixed==TRUE; for the sum of the two components in case of CN>2
    float minMinOmega = 0.15;  //used only if isMuFixed==TRUE; for one component

    //START TMP
//    std::ofstream file ("/bioinfo/users/vboeva/Desktop/Neuroblastome/ANALYSES_CNG/Valentina_analysis/freec/rowBAFS_cn16.txt");
//    for (int t = 0; t<N; t++) {
//        file << X[t]<<"\n";
//    }
//    file.close();
    //END TMP

    if (CompleteGenomicsData) { //correct mu:
        for (int i=0;i<numberOfStates;i++) {
            if (mu[i]!=HOMOZYG_MEAN)
                mu[i]*=0.82;
        }
    }

    vector <float> Omega;
    //set Omega_0:
    if (numberOfStates==2) {
        Omega.push_back(0.9);
        Omega.push_back(0.1);
    } else if (numberOfStates==3) {
        Omega.push_back(0.18);
        Omega.push_back(0.02);
        Omega.push_back(0.8);
    } else if (numberOfStates==4) {
        Omega.push_back(0.18);
        Omega.push_back(0.02);
        Omega.push_back(0.4);
        Omega.push_back(0.4);
    } else {
        cerr << "Warning: in fitting by a mixture of normals, the maximum number of compontents is 4. Using "<<numberOfStates <<" components\n";
        return NA;
    }

    vector <float> sigma;
    //set Sigma_0:
    if (numberOfStates==2) {
        sigma.push_back(0.2);
        sigma.push_back(0.2);
    } else if (numberOfStates==3) {
        sigma.push_back(0.2);
        sigma.push_back(0.2);
        sigma.push_back(maxSigma);
    } else if (numberOfStates==4) {
        sigma.push_back(0.2);
        sigma.push_back(0.2);
        sigma.push_back(maxSigma);
        sigma.push_back(maxSigma);
    }

	float ** h ;
	h = new float * [numberOfStates];
	for (int s = 0; s<numberOfStates; s++) {
        h[s] = new float[N];
	}


    for (int t = 0; t<N; t++) {
        float sum_t = 0;
        for (int i = 0; i<numberOfStates; i++) {
            sum_t+=Omega.at(i)*NormalDistributionDensity(X[t],mu[i],sigma[i]);
        }
        for (int s = 0; s<numberOfStates; s++) {
            h[s][t] = Omega.at(s)*NormalDistributionDensity(X[t],mu[s],sigma[s])/sum_t;
            if (t==0 || t == 1 || t == 376 || t == 377) {
//                   cout << s << "," << t << ", X[t]=" << X[t] << ", h[s][t] =" <<h[s][t] << ", Omega.at(s)=" <<Omega.at(s)<< "\n";
            }
        }
    }



    double LogLikelyHood_iPlusOne = 0;
//    for (int s = 0; s<numberOfStates; s++) {
//        for (int t = 0; t<N; t++) {
//            LogLikelyHood_iPlusOne += (log(Omega[s])+log(NormalDistributionDensity(X[t],mu[s],sigma[s])))*h[s][t];
//        }
//    }
//    LogLikelyHood_iPlusOne = 0;
    for (int t = 0; t<N; t++) {
        double sum_t = 0;
        for (int s = 0; s<numberOfStates; s++) {
            sum_t += Omega[s]*NormalDistributionDensity(X[t],mu[s],sigma[s]);
        }
        LogLikelyHood_iPlusOne += log(sum_t);
    }


    while (fabs(LogLikelyHood_i-LogLikelyHood_iPlusOne)>=0.0001 && iterationCount<maxIterationCount) {
        LogLikelyHood_i = LogLikelyHood_iPlusOne;
        iterationCount ++;
        //new interation:
        for (int s = 0; s<numberOfStates; s++) {
            float hSum = 0;
            float hXSum = 0;
            float varSum = 0;
            for (int t = 0; t<N; t++) {
                hSum += h[s][t];
                hXSum += h[s][t]*X[t];//DELETE THIS?????
                //varSum +=h[s][t]*X[t]*X[t];
                varSum +=h[s][t]*(X[t]-mu[s])*(X[t]-mu[s]);
            }
            Omega[s] = hSum/N;
            if (!isMuFixed)
                mu[s] = hXSum/hSum ;
            //sigma[s] = sqrt(varSum/hSum-mu[s]*mu[s]);

            if (hSum<=0) {
                cerr << "No enough values to fit the data.. !!!\n";
                return LogLikelyHood_i;
            }

            sigma[s] = sqrt(varSum/hSum);

            if (varSum/hSum<=0) {
                cerr << "No enough values to fit the data.. !!!\n";
                return LogLikelyHood_i;
            }

//            cout << s<< ": mu="<<mu[s] << ", w=" <<Omega[s]<< ", s=" <<sigma[s]<< "\n";

        }


        if (isMuFixed && numberOfStates==3 ) {
            if (Omega[2]<minOmega) {
                float prop = Omega[1]/Omega[0];
                Omega[2]=minOmega;
                Omega[0] = (1-Omega[2])/(1+prop);
                Omega[1] = 1-Omega[2]-Omega[0];

            }
            if (sigma[2]>maxSigma)
                sigma[2]=maxSigma;
        }
        if (isMuFixed && numberOfStates==4) {

            if (Omega[2]+Omega[3]<minOmega){
                float osum = Omega[2]+Omega[3];
                Omega[2]=(minOmega*Omega[2])/osum;
                Omega[3]=(minOmega*Omega[3])/osum;
                float prop = Omega[1]/Omega[0];
                Omega[0] = (1-minOmega)/(1+prop);
                Omega[1] = 1-minOmega-Omega[0];

            }
            if (Omega[2]<minMinOmega||Omega[3]<minMinOmega) {
                if (Omega[2]<minMinOmega&&Omega[3]<minMinOmega) { //it is not possible
                    Omega[2]=max(minOmega/2,minMinOmega);
                    Omega[3]=Omega[2];
                    float osum = Omega[2]+Omega[3];
                    float prop = Omega[1]/Omega[0];
                    Omega[0] = (1-osum)/(1+prop);
                    Omega[1] = 1-osum-Omega[0];
                } else if (Omega[2]<minMinOmega) {
                    Omega[2]=minMinOmega;
                    float osum = min(Omega[2]+Omega[3],1.f);
                    float prop = Omega[1]/Omega[0];
                    Omega[0] = (1-osum)/(1+prop);
                    Omega[1] = 1-osum-Omega[0];
                } else {
                    Omega[3]=minMinOmega;
                    float osum = min(Omega[2]+Omega[3],1.f);
                    float prop = Omega[1]/Omega[0];
                    Omega[0] = (1-osum)/(1+prop);
                    Omega[1] = 1-osum-Omega[0];
                }
            }

            if (sigma[2]>maxSigma)
                sigma[2]=maxSigma;
            if (sigma[3]>maxSigma)
                sigma[3]=maxSigma;
        }
//        for (int s = 0; s<numberOfStates; s++) {
//            cout << s<< ": mu="<<mu[s] << ", w=" <<Omega[s]<< ", s=" <<sigma[s]<< "\n";
//        }
        //recalculate h:

        for (int t = 0; t<N; t++) {
            float sum_t = 0;
            for (int i = 0; i<numberOfStates; i++) {
                sum_t+=Omega.at(i)*NormalDistributionDensity(X[t],mu[i],sigma[i]);
            }
            for (int s = 0; s<numberOfStates; s++) {
                h[s][t] = Omega.at(s)*NormalDistributionDensity(X[t],mu[s],sigma[s])/sum_t;
            }
        }


//        LogLikelyHood_iPlusOne = 0;
//        for (int s = 0; s<numberOfStates; s++) {
//            for (int t = 0; t<N; t++) {
//                LogLikelyHood_iPlusOne += (log(Omega[s])+log(NormalDistributionDensity(X[t],mu[s],sigma[s])))*h[s][t];
//            }
//        }

        LogLikelyHood_iPlusOne = 0;
        for (int t = 0; t<N; t++) {
            double sum_t = 0;
            for (int s = 0; s<numberOfStates; s++) {
                sum_t += Omega[s]*NormalDistributionDensity(X[t],mu[s],sigma[s]);
            }
            LogLikelyHood_iPlusOne += log(sum_t);
        }


    }

    if (iterationCount==maxIterationCount)
        cerr <<"..Warning: done in maximum number of iterations: " <<iterationCount << "\n";


    for (int i = 0; i<numberOfStates; i++)
		delete (h[i]);
	delete (h);
	sigma.clear();
    Omega.clear();

    return LogLikelyHood_iPlusOne ;
}

double NormalDistributionDensity (double x, double mu, double sigma) {
    double result;
    if (sigma)
        result = 0.398942449/sigma*exp(-(x-mu)*(x-mu)/2/sigma/sigma);
    else {
        if (x==mu) {
            result = 0.2419708*INFINITY;
        } else
            return 0;
    }
    return result;
}

string findSmallestSuffix (string s1, string s2) {

    if (s1 == "")
        return  s2;

   // cout <<"s1: "<<s1<< "\ns2: "<<s2 << "\n";

    string suffix = "";
    vector<string>values1 = split(s1, ';');
    vector<string>values2 = split(s2, ';');

    int size1 = values1.size();
    int size2 = values2.size();

    //cout <<"s1 size: "<<size1<< "; s2 size: "<<size2 << "\n";

    int startSuffix = 0;

    for (int i=max(0,size1-size2); i<size1; i++) {

        if (values1[i].compare(values2[0])==0) {
        //found possible overlap
            int realOverlap = 1;

            for (int j=i+1; j<size1; j++) {
                if (values1[j].compare(values2[j-i])!=0) {
                    j=size1;
                    realOverlap=0;
                }
            }
            if (realOverlap) {
                startSuffix = size1 - i;
                i=size1;
            }
        }
    }
    if (startSuffix == 0)
        return s2;

    if (startSuffix >= size2) //should never happen
        return "";

    suffix = values2[startSuffix];
    for (int i=startSuffix+1; i<size2; i++) {
        suffix += ";"+ values2[i];
    }

    return suffix;
}

std::string getNormalBAFforPloidy(int ploidy) {
    if (ploidy>1) {
        int Bcount = ploidy/2;
        int Acount = ploidy-Bcount;
        string BAF = "A";
        for (int i = 1; i<Acount;i++)
            BAF += "A";
        for (int i = 0; i< Bcount;i++)
            BAF += "B";
        return BAF;
    } else {
        if (ploidy==1)
            return "A";
    }
    return "";
}

std::string getXYBAFforPloidy(int ploidy) {
    if (ploidy>1) {
        int Acount = ploidy/2;
        string BAF = "A";
        for (int i = 1; i<Acount;i++)
            BAF += "A";
        return BAF;
    } else {
        if (ploidy==1)
            return "A";
    }
    return "";
}

bool getELANDinfo(std::string line,std::string &chr1,std::string &chr2,std::string &orient1,std::string &orient2,int &left,int &right, int &insertSize) {
    if (! line.length())
        return false;

    std::vector<std::string> strs = split(line, '\t');
    if (strs.size()>=14) {
        chr1=strs[11];
        chr2=strs[12];
        processChrName(chr1);
        processChrName(chr2);

        orient1 = strs[8];
        orient2 =strs[10];

        left = atoi(strs[7].c_str());
		right = atoi(strs[9].c_str());

		insertSize=atoi(strs[13].c_str());

        if (chr1.compare(chr2)!=0) {//Can be avoided, TODO {
            insertSize=0;
            return false;
        }
    } else {
        return false;
    }
    strs.clear();
    return true;
}

bool getSAMinfo(const char* line, std::string &chr1, std::string &chr2, char& orient1, char& orient2, int &left,int &right, int &insert_size) {
    if (!*line) {
	  return false;
	}
	if (line[0] == '@') {
	  return false;
	}

	char* strs[32];
	unsigned int strs_cnt = split((char*)line, '\t', strs);
	if (strs_cnt < 7)
        return false;

    chr1 = strs[2];
	chr2 = strs[6];
    processChrName(chr1);
    processChrName(chr2);
    if (chr1.compare("*")==0 || chr2.compare("*")==0) {
        return false;
	}

	if (chr2.compare("=")==0) {
	  chr2 = chr1;
	}

	unsigned int mask = atoi(strs[1]);
    if (mask & 0x0010) {
	  orient1 = 'R';
	} else {
	  orient1 = 'F';
	}

    if (mask & 0x0020) {
        orient2 = 'R';
	} else {
        orient2 = 'F';
	}

    left = atoi(strs[3]);
    right = atoi(strs[7]);

    insert_size = atoi(strs[8]);

    if (chr1.compare(chr2)!=0) {
        return false;//Can be avoided, TODO
	}
    return true;

}

bool getSAMinfo(const std::string& line,std::string &chr1,std::string &chr2,std::string &orient1,std::string &orient2,int &left,int &right) {
    if (! line.length())
        return false;

    if ( line[0] == '@')
        return false;

    std::vector<std::string> strs = split(line, '\t');
    if (strs.size() < 7)
        return false;

    chr1 = strs[2];
    chr2 = strs[6];
    processChrName(chr1);
    processChrName(chr2);
    if (chr1.compare("*")==0 || chr2.compare("*")==0)
        return false;

    if(chr2.compare("=")==0)
        chr2=chr1;



    if (atoi(strs[1].c_str()) & 0x0010)
        orient1 = "R";
    else
        orient1 = "F";

    if (atoi(strs[1].c_str()) & 0x0020)
        orient2 = "R";
    else
        orient2 = "F";

    left = atoi(strs[3].c_str());
    right = atoi(strs[7].c_str());

    strs.clear();
    if (chr1.compare(chr2)!=0)
        return false;//Can be avoided, TODO
    return true;

}

int calculateTotalLength(std::vector <int> lefts,std::vector <int> rights) {
    if (lefts.size()==0 || rights.size()==0 || lefts.size()!=rights.size())
        return 0;
    int length = rights[0]-lefts[0]+1;
        for (int unsigned i=1; i< rights.size(); i++) {
            if (lefts[i]<lefts[i-1]) {
                cerr << "..unsorted list in calculateTotalLength()..";
                cerr << "..somatic/germline attribution can be incorrect..";
                return 0;
            }
            if (lefts[i]>rights[i-1])
                length += rights[i]-lefts[i]+1;
            else if (rights[i]>rights[i-1]) {
                length += rights[i]-rights[i-1];
            }
        }
    return length;
}

void myReplace(std::string& str, const std::string& oldStr, const std::string& newStr) {
  size_t pos = 0;
  while((pos = str.find(oldStr, pos)) != std::string::npos) {
     str.replace(pos, oldStr.length(), newStr);
     pos += newStr.length();
  }
}


#if defined(_WIN32) || (defined(__APPLE__) && defined(__MACH__))
double expm1(double x) {
        if (fabs(x) < 1e-5)
                return x + 0.5*x*x;
        else
                return exp(x) - 1.0;
}
double log1p(double x) {
        if (fabs(x) < 1e-5)
                return (x-x*x/2+x*x*x/3);
        else
                return log(x+1);
}
#endif

MateOrientation getMateOrientation(std::string const& matesOrientation)
{
  if (matesOrientation.compare("0") == 0) {
	return SINGLE_END_SORTED_SAM;
  }

  if (matesOrientation.compare("RF") == 0) {
	return ILLUMINA_MATE_PAIRS;
  }

  if (matesOrientation.compare("FR") == 0) {
	return ILLUMNINA_PAIRED_END;
  }

  if (matesOrientation.compare("FF") == 0 || matesOrientation.compare("RR") == 0 || matesOrientation.compare("FF/RR") == 0 || matesOrientation.compare("RR/FF") == 0) {
	return SOLID_MATE_PAIRS;
  }
  cerr << "Error: you have set an unknown read orientation: \""<<matesOrientation<<"\" ; it does not match 0, RF, FR or FF; please correct your config file.\n";
  exit(1);
  return SINGLE_END_SORTED_SAM; // instead of UNKNOWN_MATE_ORIENTATION;
}

InputFormat getInputFormat(std::string const& inputFormat)
{
  if (inputFormat.compare("sam")==0 || inputFormat.compare("SAM")==0) {
	return SAM_INPUT_FORMAT;
  }

    if (inputFormat.compare("bam")==0 || inputFormat.compare("BAM")==0) {
	return SAM_INPUT_FORMAT;
  }

  if (inputFormat.compare("eland")==0 || inputFormat.compare("Eland")==0) {
	return ELAND_INPUT_FORMAT;
  }

  if (inputFormat.compare("bowtie")==0 || inputFormat.compare("Bowtie")==0) {
	return BOWTIE_INPUT_FORMAT;
  }

  if (inputFormat.compare("psl")==0 || inputFormat.compare("BLAT")==0) {
	return PSL_INPUT_FORMAT;
  }

  if (inputFormat.compare("arachne")==0 || inputFormat.compare("BED")==0 || inputFormat.compare("bed")==0 || inputFormat.compare("ARACHNE")==0) {
	return ARACHNE_BED_INPUT_FORMAT;
  }

  if (inputFormat.compare("SOAP")==0 || inputFormat.compare("soap")==0 || inputFormat.compare("Soap")==0) {
	return SOAP_INPUT_FORMAT;
  }

  if (inputFormat.compare("pileup")==0 || inputFormat.compare("SAMtools pileup")==0) {
	return SAM_PILEUP_INPUT_FORMAT;
  }

  return UNKNOWN_INPUT_FORMAT;
}

char* getLine(char* buffer, int buffer_size, FILE* stream, std::string& line)
{
  while (fgets(buffer, buffer_size, stream) != NULL) {
	size_t len = strlen(buffer);
	if (buffer[len-1] == '\n') {
	  return buffer;
	}
	line = buffer;
	for (;;) {
	  if (line.at(line.length() - 1) == '\n') {
		break;
	  }
	  if (fgets(buffer, MAX_BUFFER, stream) != NULL) {
		line.append(buffer);
	  } else {
		break;
	  }
	}
	return (char*)line.c_str();
  }
  return NULL;
}


 std::vector<float> get_quartiles(vector<float> vect)
{
    float q1 = 0, q2 = 0, q3 = 0;
    if (vect.size() > 10)
    {
    std::sort(vect.begin(), vect.end());
    q2 = get_median(vect, 0, vect.size());
    q1 = get_median(vect, 0, round(vect.size()/2));
    q3 = get_median(vect, round(vect.size()/2), vect.size());
    }
    vector<float> quartiles;
    quartiles.push_back(q1);
    quartiles.push_back(q2);
    quartiles.push_back(q3);
    return quartiles;
}


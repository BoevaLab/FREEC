/*************************************************************************
Copyright (c) 2010-2016, Valentina BOEVA.

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


#include "BAFpileup.h"

using namespace std;

BAFpileup::BAFpileup()
{
    //ctor
}

void BAFpileup::makepileup(GenomeCopyNumber & sampleCopyNumber, GenomeCopyNumber & controlCopyNumber,
        std::string sample_MateFile, std::string control_Matefile, std::string outputDir, std::string makeminipileup,
        std::string const& mateFileName ,std::string const& inputFormat, std::string const& matesOrientation,
        std::string pathToSamtools, std::string pathToSambamba, std::string SambambaThreads, std::string chrLenFileName, std::string controlName, std::string targetBed,  std::string pathToBedtools,
        std::string fastaFile, int minQualPerPos)
{
    //create a .bed file with regions of interest to create a minipileup: targeted + flanks for WES or all chromosomes for WGS:
    std::string bedFileWithRegionsOfInterest = outputDir + "_NewCaptureRegions" +  ".bed";

    if (targetBed != "")
        {
        int flanks = calculateFlankLength(mateFileName, inputFormat, matesOrientation, pathToSamtools, pathToSambamba, SambambaThreads);
        calculateNewBoundaries(targetBed, flanks, bedFileWithRegionsOfInterest);
        }
    else
        {
        createBedFileWithChromosomeLengths(bedFileWithRegionsOfInterest,chrLenFileName, true);
        }
    pathToBedtools_=pathToBedtools; // /*
    string intersected = intersectWithBedtools(makeminipileup, outputDir, bedFileWithRegionsOfInterest, chrLenFileName);
    string sampleOutFileName = createPileUpFile( outputDir, pathToSamtools , pathToSambamba, SambambaThreads, sample_MateFile, intersected, fastaFile,minQualPerPos);

    //BAFtumor = computeBAF(sampleCopyNumber, _sample, outputDir, "_sample");

    if (controlName.compare("")!=0) {
        string controlOutFileName = createPileUpFile( controlName, pathToSamtools, pathToSambamba, SambambaThreads, control_Matefile, intersected, fastaFile,minQualPerPos);
    }
    //computeBAF(controlCopyNumber, _control, outputDir, "_control");
    remove(intersected.c_str()); // */
}

float BAFpileup::calculateFlankLength(std::string const& mateFileName, std::string const& inputFormat_str, std::string const& matesOrientation_str, std::string pathToSamtools_, std::string pathToSambamba, std::string SambambaThreads)
{
        if (matesOrientation_str=="0") return 0; // do not add anything in case of single end data
        if (getInputFormat(inputFormat_str)!=SAM_INPUT_FORMAT)  return 0;
        std::ifstream fileMates (mateFileName.c_str());
        vector<float> insertSizeVector;
        string line;
        if (!fileMates.is_open())        {
            cerr << "Error: unable to open "+mateFileName+"\n" ;
            exit(-1);
        }
        fileMates.close();

        #ifdef PROFILE_TRACE
            time_t t0 = time(NULL);
        #endif

       // MateOrientation matesOrientation = getMateOrientation(matesOrientation_str); // will not consider read orientation to increase speed
        char* line_buffer;
        FILE *stream;
        char buffer[MAX_BUFFER];
        bool zgOrbam=0;

        int numberOfReadsToCheck=3000000;
        int j = 0;
        int fragmentLength=0;
        if(mateFileName.substr(mateFileName.size()-3,3).compare("bam")==0 || mateFileName.substr(mateFileName.size()-3,3).compare(".gz")==0) {
            string command;
            if (mateFileName.substr(mateFileName.size()-3,3).compare("bam")==0) {
                if (pathToSambamba != "")     {
                     command = pathToSambamba + " view -t " + SambambaThreads + " " + mateFileName;
                }    else    {
                     command = pathToSamtools_ + " view "+mateFileName;
                }
            }
            if (mateFileName.substr(mateFileName.size()-3,3).compare(".gz")==0) {
                          command = "gzip -c -d "+mateFileName;
            }
            stream =
            #if defined(_WIN32)
                _popen(command.c_str(), "r");
            #else
                popen(command.c_str(), "r");
            #endif
            zgOrbam=1;

            while (((line_buffer = getLine(buffer, MAX_BUFFER, stream, line)) != NULL) && j<numberOfReadsToCheck) {
                if (line_buffer[0] == '@')
                    continue;

                char* strs[32];
                unsigned int strs_cnt = split((char*)line_buffer, '\t', strs);

                if (strs_cnt > 8 && strs[6][0]=='=') {
                    int frLen=atoi(strs[8]);
                    if(frLen>0 && frLen< 10000) {
                        fragmentLength+=frLen;
                        j++;
                    }
                }
            }

        } else {
            fileMates.open(mateFileName.c_str());
            while( getline( fileMates, line ) ) {
                if (line[0] == '@')
                    continue;
                vector <string> strs;
                strs = split(line, '\t');
                if (strs.size() > 8 && strs[6][0]=='=') {
                    int frLen=atoi(strs[8].c_str());
                    if(frLen>0 && frLen< 10000) {
                        fragmentLength+=frLen;
                        j++;
                    }
                }
            }
            fileMates.close();
        }


        fragmentLength = fragmentLength/j;
        float flanks = fragmentLength/2;
        if (zgOrbam) {
            #if defined(_WIN32)
				_pclose(stream);
            #else
                    pclose(stream);
            #endif

        }

        cout << "..will increase flanking regions by "<< flanks << " bp"<<endl;
        return flanks;
}

void BAFpileup::calculateNewBoundaries(std::string targetBed, int flanks, std::string bedFileWithRegionsOfInterest)
{
        std::string const& captureFile = targetBed ;
        ifstream file (captureFile.c_str());

        ofstream myfile;
        myfile.open(bedFileWithRegionsOfInterest.c_str());


        if (file.is_open())  {
            std::string line;
            exons_Count = 0;
            exons_Counttmp = 0;
            while (std::getline(file,line)) {
                vector <string> strs = split(line, '\t');                // Print new capture regions
                myfile << strs[0] << "\t" << atoi(strs[1].c_str())-flanks <<
                 "\t" << atoi(strs[2].c_str())+flanks<< "\t" << strs[0]<<":"
                 <<atoi(strs[1].c_str())-flanks << "-" << atoi(strs[2].c_str())+flanks<< "\t"<<
                 atoi(strs[2].c_str())- atoi(strs[1].c_str())+2*flanks<< "\t" << "+" << "\n";
            }
            myfile.close();
        } else       {
                cerr << "Error: Unable to open file "+captureFile+"\n";
                myfile.close();
                exit(-1);
        }
}

void BAFpileup::createBedFileWithChromosomeLengths(std::string bedFileWithRegionsOfInterest, std::string chrLenFileName, bool doesNeedChrPrefix) {

    //reading the file with the chromosome length information
    std::vector<std::string> chr_names;
	std::vector<int> lengths;
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
		    chr_names.push_back(name);
            lengths.push_back(value);
            //cout << name << "\t" << value << "\n";
		}
	}
	file.close();
	cout << "..File "<<chrLenFileName<<" was read to create a miniPileup\n";

    ofstream myfile;
    myfile.open(bedFileWithRegionsOfInterest.c_str());
    for (unsigned int i = 0; i < chr_names.size(); i++)  {
        if(doesNeedChrPrefix)
            myfile << "chr"<< chr_names[i] << "\t" << "1" << "\t" << lengths[i] << "\n";
        else
            myfile << chr_names[i] << "\t" << "1" << "\t" << lengths[i] << "\n";
    }
    myfile.close();
}

std::string BAFpileup::intersectWithBedtools(std::string makeminipileup, std::string outputDir, std::string bedFileWithRegionsOfInterest, std::string chrLen)
{
    FILE *stream;
    std::string intersected = outputDir + "_SNPinNewCaptureRegions.bed";

    if (makeminipileup.substr(makeminipileup.size()-3,3).compare("vcf")==0 || makeminipileup.substr(makeminipileup.size()-6,6).compare("vcf.gz")==0) {
        intersected = outputDir + "_SNPinNewCaptureRegions.vcf";
    }

    string command = pathToBedtools_ +" intersect -a " + makeminipileup + " -b " + bedFileWithRegionsOfInterest + " > " + intersected;

    stream =
    #if defined(_WIN32)
        _popen(command.c_str(), "w");
    #else
        popen(command.c_str(), "w");
    #endif

    #if defined(_WIN32)
    _pclose(stream);
    #else
    pclose(stream);
    #endif

    remove(bedFileWithRegionsOfInterest.c_str());

    if (makeminipileup.substr(makeminipileup.size()-3,3).compare("bed")==0) {
        return intersected;
    }

    //else: VCF input format; need to transform into BED
    ifstream file (intersected.c_str());
    ofstream myfile;
    std::string intersectedBed = outputDir + "_SNPinNewCaptureRegions.bed";
    myfile.open(intersectedBed.c_str());
    std:: string line;
    int nb_snp = 0;
    while (std::getline(file,line))  {
        nb_snp++;
        std::vector<std::string> strs = split(line, '\t');
        myfile << strs[0] << "\t"<< atoi(strs[1].c_str()) -1 << "\t" <<  atoi(strs[1].c_str()) << "\t" << strs[3] << "\t" << strs[4]<<"\n";
    }
    cout << "..will use "<<nb_snp<<" SNPs for calculation of the BAF profile\n";
    myfile.close();
    file.close();
    remove(intersected.c_str());
    return intersectedBed;
}

std::string BAFpileup::createPileUpFile(std::string outputDir, std::string samtools_path, std::string pathToSambamba,std::string SambambaThreads , std::string control_MateFile, std::string intersected, std::string fastaFile, int minQualPerPos)
{
    string minipileup = outputDir + "_minipileup" +".pileup";
    FILE *stream;

    string command;
    if (pathToSambamba != "")    {
        string samtools_arg = "--samtools -f " +fastaFile+ " -d 8000 -Q "+int2string(minQualPerPos)+ " -q 1 -l " + intersected;
        command  = pathToSambamba + " mpileup -t " + SambambaThreads + " -o " + minipileup + " " + control_MateFile + " " + samtools_arg ;
    }  else   {
         command = samtools_path + " mpileup -f "+fastaFile+" -d 8000 -Q "+int2string(minQualPerPos)+" -q 1 -l " + intersected + " " + control_MateFile + " > " + minipileup; //discard reads wit 0 mapping quality
    }

    stream =
    #if defined(_WIN32)
        _popen(command.c_str(), "w");
    #else
        popen(command.c_str(), "w");
    #endif

    #if defined(_WIN32)
        _pclose(stream);
    #else
        pclose(stream);
    #endif
    cout << "..If you have got an error at this step and a mini-pileup file is empty, check that you are using samtools v1.1 or later and provide a corresponding path in your config file"<<endl;
    return minipileup;
}

/*
std::vector < std::vector<float> >BAFpileup::computeBAF(GenomeCopyNumber & sampleorcontrol, std::string minipileup, std::string outputDir,  std::string filename)
{

    int numberofObjects = 0;

    vector<ChrCopyNumber>::iterator it;
    for ( it=sampleorcontrol.chrCopyNumber_.begin() ; it != sampleorcontrol.chrCopyNumber_.end(); it++ )
        {
            numberofObjects++;
        }

    string line;
    ifstream file (minipileup.c_str());
    int nb_snp = 0;
    while (std::getline(file,line))
        {
        nb_snp++;
        }

    file.clear();
    file.seekg(0);

    chr = vector < vector<string> >(numberofObjects);
    snp_pos_pileup = vector < vector<string> >(numberofObjects);
    A_nb = vector < vector<int> >(numberofObjects);
    B_nb = vector < vector<int> >(numberofObjects);
    AB_nb = vector < vector<string> >(numberofObjects);


    int k = 0;
    for ( it=sampleorcontrol.chrCopyNumber_.begin() ; it != sampleorcontrol.chrCopyNumber_.end(); it++ )
        {
        file.clear();
        file.seekg(0);
        int i = 0;
        int l = 0;
        while (std::getline(file,line) && (k < numberofObjects))
            {
            i = 0;
            string chrtmp;
            while (line[i] != '\t')
                {
                chrtmp += line[i];
                i++;
                }
                i++;
            if (chrtmp == "chr" + it->getChromosome())
            {

                chr[k].push_back(chrtmp);
                string snp_pos_pileuptmp;
                while (line[i] != '\t')
                    {
                    snp_pos_pileuptmp += line[i];
                    i++;
                    }
                    snp_pos_pileup[k].push_back(snp_pos_pileuptmp);
                    i++;

                int j = 0;
                while ((snp_pos_pileup[k][l] != snp_pos[j]) && (j < nb_snp))
                    {
                    j++;
                    }

                while (line[i] != '\t')
                    {
                    i++;
                    }
                    i++;
                string AB_nbtmp;
                while (line[i] != '\t')
                    {
                    AB_nbtmp += line[i];
                    i++;
                    }
                    AB_nb[k].push_back(AB_nbtmp);
                    i++;
                char* alternate_base = const_cast<char*>(alt_base[j].c_str());
                char alternate_basetmp = *alternate_base;
                char* reference_base = const_cast<char*>(ref_base[j].c_str());
                char reference_basetmp = *reference_base;
                B_nb[k].push_back(0);
                A_nb[k].push_back(0);
                while (line[i] != '\t')
                    {
                    if ((char)line[i] == reference_basetmp || (char)line[i] == (reference_basetmp+32) || ((char)line[i] == (reference_basetmp-32)))
                        {
                        B_nb[k][l]++;
                        }
                    if ((char)line[i] == alternate_basetmp || (char)line[i] == (alternate_basetmp+32) || ((char)line[i] == (alternate_basetmp-32)))
                        {
                        A_nb[k][l]++;
                        }
                    i++;
                    }
                    i++;
                l++;
                }
            }
        k++;
        }

    string Newfile = outputDir + "_BAF" + ".txt";
    ofstream myfile;
    std::vector < std::vector<float> > BAF = vector < vector<float> > (numberofObjects);

    int j = 0;
    for ( it=sampleorcontrol.chrCopyNumber_.begin() ; it != sampleorcontrol.chrCopyNumber_.end(); it++ )
        {
        for (int i = 0; i < chr[j].size(); i++)
        {
            BAF[j].push_back(((float)B_nb[j][i]) / ((float)A_nb[j][i] + (float)B_nb[j][i]));
        }
        j++;
    }

    j = 0;
    myfile.open(Newfile.c_str());
    for ( it=sampleorcontrol.chrCopyNumber_.begin() ; it != sampleorcontrol.chrCopyNumber_.end(); it++ )
        {
        for (int i = 0; i < chr[j].size(); i++)
        {
            myfile << chr[j][i] << "\t" << snp_pos_pileup[j][i] << "\t" << AB_nb[j][i] << "\t" << B_nb[j][i] << "\t" << A_nb[j][i] << "\t" << BAF[j][i]<<"\n";
        }
        j++;
    }
    myfile.close();
    return BAF;
}
*/

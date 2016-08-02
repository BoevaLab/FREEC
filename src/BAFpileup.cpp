#include "BAFpileup.h"

using namespace std;

BAFpileup::BAFpileup()
{
    //ctor
}

void BAFpileup::makepileup(GenomeCopyNumber & sampleCopyNumber, GenomeCopyNumber & controlCopyNumber,
        std::string sample_MateFile, std::string control_Matefile, std::string outputDir, std::string makeminipileup,
        std::string const& mateFileName ,std::string const& inputFormat, std::string const& matesOrientation,
        std::string pathToSamtools, std::string chrLen, std::string controlName, std::string targetBed,  std::string pathToBedtools,
        std::string fastaFile, int minQualPerPos)
{
    if (targetBed != "")
        {
        int flanks = calculateFlankLength(mateFileName, inputFormat, matesOrientation, pathToSamtools);
        calculateNewBoundaries(targetBed, flanks, outputDir);
        }
    else
        {
        printfile(outputDir,chrLen);
        }
    pathToBedtools_=pathToBedtools; // /*
    string intersected = intersectWithBedtools(makeminipileup, outputDir, targetBed, chrLen);
    string _sample = createPileUpFile( outputDir, pathToSamtools, sample_MateFile, intersected, fastaFile,minQualPerPos);
    //BAFtumor = computeBAF(sampleCopyNumber, _sample, outputDir, "_sample");
    string _control = createPileUpFile( controlName, pathToSamtools, control_Matefile, intersected, fastaFile,minQualPerPos);
    //computeBAF(controlCopyNumber, _control, outputDir, "_control");
    remove(intersected.c_str()); // */
}

float BAFpileup::calculateFlankLength(std::string const& mateFileName, std::string const& inputFormat_str, std::string const& matesOrientation_str, std::string pathToSamtools_)
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

        MateOrientation matesOrientation = getMateOrientation(matesOrientation_str);
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
                      command = pathToSamtools_ + " view "+mateFileName;
            }
            if (mateFileName.substr(mateFileName.size()-3,3).compare(".gz")==0) {
                          command = "gzip -c -d "+mateFileName;
            }
            stream =
            #if defined(_WIN32) || (defined(__APPLE__) && defined(__MACH__))
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
            #if defined(_WIN32) || (defined(__APPLE__) && defined(__MACH__))
				_pclose(stream);
            #else
                    pclose(stream);
            #endif

        }

        cout << "..will increase flanking regions by "<< flanks << " bp"<<endl;
        return flanks;
}

void BAFpileup::calculateNewBoundaries(std::string targetBed, int flanks, std::string outputDir)
{
        std::string const& captureFile = targetBed ;
        ifstream file (captureFile.c_str());

        ofstream myfile;
        std::string Newfile = outputDir + "_NewCaptureRegions" +  ".bed";
        myfile.open(Newfile.c_str());


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

void BAFpileup::printfile(std::string outputDir, std::string chrLen)
{
    std::string const& captureFile = chrLen ;
    ifstream file (captureFile.c_str());
    std::string line;

    int compt = 0;
    while (std::getline(file,line))
    {
    compt++;
    }

    int l=0;
    file.clear();
    file.seekg(0);

    chr_names = vector<string>(compt);
    ends_ = vector<int>(compt,0);

    int i = 0;
    if (file.is_open() & l < compt)
        {
        while (std::getline(file,line))
            {
            i = 0;
            while (line[i] != '\t')
                {
                i++;
                }
            i++;
            while (line[i] != '\t')
                {
                chr_names[l] += line[i];
                i++;
                }
            i++;
            string endtmp;
            while (line[i] != '\t')
                {
                endtmp += line[i];
                i++;
                }
            ends_[l] = atoi(endtmp.c_str());
            endtmp.clear();
            l++;
            }
        }

    ofstream myfile;
    std::string Newfile = outputDir + "_NewCaptureRegions" +  ".bed";
    myfile.open(Newfile.c_str());
    for (int i = 0; i < chr_names.size(); i++)
        {
        myfile << chr_names[i] << "\t" << "0" << "\t" << ends_[i] << "\n";
        }
    myfile.close();
}

std::string BAFpileup::intersectWithBedtools(std::string makeminipileup, std::string outputDir, std::string targetBed, std::string chrLen)
{
    FILE *stream;
    string NewCaptRegions = outputDir + "_NewCaptureRegions" +  ".bed";
    string command = pathToBedtools_ +" intersect -a " + makeminipileup + " -b " + NewCaptRegions + " > " + outputDir + "_SNPinNewCaptureRegions.bed";

    stream =
    #if defined(_WIN32) || (defined(__APPLE__) && defined(__MACH__))
        _popen(command.c_str(), "w");
    #else
        popen(command.c_str(), "w");
    #endif

    #if defined(_WIN32) || (defined(__APPLE__) && defined(__MACH__))
    _pclose(stream);
    #else
    pclose(stream);
    #endif

    remove(NewCaptRegions.c_str());
    std::string intersected = outputDir + "_SNPinNewCaptureRegions.bed";

    if (makeminipileup.substr(makeminipileup.size()-3,3).compare("bed")==0) {
        return intersected;
    }
    //else: VCF input format; need to transform into BED
    ifstream file (intersected.c_str());
    std:: string line;

    int nb_snp = 0;
    while (std::getline(file,line))
                {
                    nb_snp++;
                }
    int l = 0;
    int i = 0;

    file.clear();
    file.seekg(0);

    vector<string> chromosome = vector<string>(nb_snp);
    snp_pos = vector<string>(nb_snp);
    ref_base = vector<string>(nb_snp);
    alt_base = vector<string>(nb_snp);
    while (std::getline(file,line) && l < nb_snp)
                {
                    i = 0;
                    while (line[i] != '\t')
                        {
                        chromosome[l] += line[i];
                        i++;
                        }
                        i++;
                    while (line[i] != '\t')
                        {
                        snp_pos[l] += line[i];
                        i++;
                        }
                        i++;
                    while (line[i] != '\t')
                        {
                        i++;
                        }
                        i++;
                    while (line[i] != '\t')
                        {
                        ref_base[l] += line[i];
                        i++;
                        }
                        i++;
                    while (line[i] != '\t')
                        {
                        alt_base[l] += line[i];
                        i++;
                        }
                        i++;
                    while (line[i] != '\n')
                        {
                        i++;
                        }
                    l++;
                }
    ofstream myfile;
    myfile.open(intersected.c_str());
    for (i = 0; i < nb_snp; i++)
    {
    myfile << chromosome[i] << "\t"<< atoi(snp_pos[i].c_str()) -1 << "\t" <<  atoi(snp_pos[i].c_str()) << "\t" << ref_base[i] << "\t" << alt_base[i] <<"\n";
    }
    myfile.close();
    return intersected;
}


std::string BAFpileup::createPileUpFile(std::string outputDir, std::string samtools_path,std::string control_MateFile, std::string intersected, std::string fastaFile, int minQualPerPos)
{
    string minipileup = outputDir + "_minipileup" +".pileup";
    FILE *stream;
    string command = samtools_path + " mpileup -f "+fastaFile+" -d8000 -Q "+int2string(minQualPerPos)+" -q 1 -l " + intersected + " " + control_MateFile + " > " + minipileup; //discard reads wit 0 mapping quality


     stream =
    #if defined(_WIN32) || (defined(__APPLE__) && defined(__MACH__))
        _popen(command.c_str(), "w");
    #else
        popen(command.c_str(), "w");
    #endif

    #if defined(_WIN32) || (defined(__APPLE__) && defined(__MACH__))
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

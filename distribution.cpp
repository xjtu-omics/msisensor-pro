/*
 * distribution.cpp for MSIsensor
 * Copyright (c) 2013 Beifang Niu && Kai Ye WUGSC All Rights Reserved.
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

// System header files
#include <iomanip>
#include <iostream>
#include <fstream>
#include <cmath>
#include <assert.h>
#include <unistd.h>
#include <cstdio>
#include <string>
#include <vector>
#include <algorithm>
#include<sys/stat.h>
#include <sys/types.h>
#include <map>
//#include <io.h>

// Static function declaration
#include "param.h"
#include "polyscan.h"
#include "distribution.h"
#include "utilities.h"
#include "sample.h"


// branch
#include "cmds.h"

Param paramd;
PolyScan polyscan;
Sample sample;
//Train train;
std::map <std::string,int > SitesSupport;


std::string homoFile;
std::string bamListFile;
std::string normalBam;
std::string TrainBamConfig;//Yelab
std::vector<std::string> TrainName;//Yelab
std::vector<std::string> TrainBam;//Yelab
std::string tumorBam;
std::string bedFile;
std::string disFile;



std::ifstream finH;
std::ifstream finM;
std::ifstream finB;
std::ofstream foutD;
std::ofstream foutO;

std::string one_region;

int loadFilepathFromConfig(std::string TrainBamConfig, std::vector<std::string>& TrainBam) {//Yelab
	std::ifstream fin;
	std::string oneLine = "";
	fin.open(TrainBamConfig.c_str());
	if(!fin)
	{
		std::cerr << "Open configure file failed, please provide valid configure file ! \n"; exit(0);
	}
	while (getline(fin, oneLine)) {
//		getline(fin, oneLine);
		std::string name="";
		std::string bam="";

		std::stringstream linestream(oneLine);
		linestream >> name;
		linestream >> bam;
		const char *pathname=bam.c_str();
		if ( access(pathname, R_OK)==-1){
			std::cerr << "Open bam file failed( "<<bam<<" ), please provide valid bam file with its index file ! \n"; exit(0);;
		}
		else{
            std::cout<< "load bam:"<<bam<<" OK!\n";
		}
		if (!name.empty()){
			TrainName.push_back(name);
			TrainBam.push_back(bam);
		}
	}
	fin.close();
}
void DisUsage(void) {
    std::cerr<<"\nUsage:  msisensor-pro msi [options] \n\n"
        <<"       -d   <string>   homopolymers and microsatellites file\n"
        <<"       -n   <string>   normal bam file with index\n"
        <<"       -t   <string>   tumor  bam file with index\n"
        <<"       -o   <string>   output prefix\n\n"

        <<"       -e   <string>   bed file, optional\n"
        <<"       -f   <double>   FDR threshold for somatic sites detection, default="<<paramd.fdrThreshold<<"\n"
//        <<"       -i   <double>   minimal comentropy threshold for somatic sites detection (just for tumor only data), default="<<paramd.comentropyThreshold<<"\n"
        <<"       -c   <int>      coverage threshold for msi analysis, WXS: 20; WGS: 15, default="<<paramd.covCutoff<<"\n"
        <<"       -z   <int>      coverage normalization for paired tumor and normal data, 0: no; 1: yes, default="<<paramd.Normalization<<"\n"
        <<"       -r   <string>   choose one region, format: 1:10000000-20000000\n"
//        <<"       -l   <int>      minimal homopolymer size, default="<<paramd.MininalHomoSize<<"\n"
        <<"       -p   <int>      minimal homopolymer size for distribution analysis, default="<<paramd.MininalHomoForDis<<"\n"
        <<"       -m   <int>      maximal homopolymer size for distribution analysis, default="<<paramd.MaxHomoSize<<"\n"

        <<"       -q   <int>      minimal microsatellite size, default="<<paramd.MinMicrosate<<"\n"
        <<"       -s   <int>      minimal microsatellite size for distribution analysis, default="<<paramd.MinMicrosateForDis<<"\n"
        <<"       -w   <int>      maximal microsatellite size for distribution analysis, default="<<paramd.MaxMicrosateForDis<<"\n"

        <<"       -u   <int>      span size around window for extracting reads, default="<<paramd.DisSpan<<"\n"
        <<"       -b   <int>      threads number for parallel computing, default="<<paramd.numberThreads<<"\n"
        <<"       -x   <int>      output homopolymer only, 0: no; 1: yes, default="<<paramd.HomoOnly<<"\n"
        <<"       -y   <int>      output microsatellites only, 0: no; 1: yes, default="<<paramd.MicrosateOnly<<"\n"
		<<"       -0   <int>      output site have no read coverage, 1: no; 0: yes, default=" << paramd.outputzeroDis << "\n"

        <<"       \n"
        <<"       -h   help\n\n"
        << std::endl;
    exit(1);
}


// add by YeLab
void ProUsage(void) {


	std::cerr << "\nUsage:  msisensor pro [options] \n\n"
		<< "       -d   <string>   homopolymer and microsates file\n"
		<< "       -t   <string>   tumor bam file\n"
		<< "       -o   <string>   output prefix\n\n"

		<< "       -e   <string>   bed file, optional\n"
//		<< "       -f   <double>   FDR threshold for somatic sites detection, default=" << paramd.fdrThreshold << "\n"
		<< "       -i   <double>   minimal msihunter threshold for instable sites detection (just for tumor only data), default=" << paramd.hunterThreshold << "\n"
		<< "       -c   <int>      coverage threshold for msi analysis, WXS: 20; WGS: 15, default=" << paramd.covCutoff << "\n"

//      << "       -z   <int>      coverage normalization for paired tumor and normal data, 0: no; 1: yes, default=" << paramd.Normalization << "\n"
		<< "       -r   <string>   choose one region, format: 1:10000000-20000000\n"

//		<< "       -l   <int>      minimal homopolymer size, default=" << paramd.MininalHomoSize << "\n"
		<< "       -p   <int>      minimal homopolymer size for distribution analysis, default=" << paramd.MininalHomoForDis << "\n"
		<< "       -m   <int>      maximal homopolymer size for distribution analysis, default=" << paramd.MaxHomoSize << "\n"

//		<< "       -q   <int>      minimal microsates size, default=" << paramd.MinMicrosate << "\n"
		<< "       -s   <int>      minimal microsatellite size for distribution analysis, default=" << paramd.MinMicrosateForDis << "\n"
		<< "       -w   <int>      maximal microsatellite size for distribution analysis, default=" << paramd.MaxMicrosateForDis << "\n"

		<< "       -u   <int>      span size around window for extracting reads, default=" << paramd.DisSpan << "\n"
		<< "       -b   <int>      threads number for parallel computing, default=" << paramd.numberThreads << "\n"
		<< "       -x   <int>      output homopolymer only, 0: no; 1: yes, default=" << paramd.HomoOnly << "\n"
		<< "       -y   <int>      output microsatellite only, 0: no; 1: yes, default=" << paramd.MicrosateOnly << "\n"
		<< "       -0   <int>      output site have no read coverage, 1: no; 0: yes, default=" << paramd.outputzeroDis << "\n"
		<< "       \n"
		<< "       -h   help\n\n"
		<< std::endl;
	exit(1);
}
void EntropyDisUsage(void) {
    std::cerr<<"\nUsage:  msisensor entropy [options] \n\n"
        <<"       -d   <string>   homopolymer and microsatellite file\n"
        <<"       -t   <string>   tumor  bam file\n"
        <<"       -o   <string>   output prefix\n\n"

        <<"       -e   <string>   bed file, optional\n"
        <<"       -i   <double>   minimal comentropy threshold for somatic sites detection (just for tumor only data), default="<<paramd.comentropyThreshold<<"\n"
        <<"       -c   <int>      coverage threshold for msi analysis, WXS: 20; WGS: 15, default="<<paramd.covCutoff<<"\n"
        <<"       -r   <string>   choose one region, format: 1:10000000-20000000\n"
        <<"       -l   <int>      minimal homopolymer size, default="<<paramd.MininalHomoSize<<"\n"
        <<"       -p   <int>      minimal homopolymer size for distribution analysis, default="<<paramd.MininalHomoForDis<<"\n"
        <<"       -m   <int>      maximal homopolymer size for distribution analysis, default="<<paramd.MaxHomoSize<<"\n"

        <<"       -q   <int>      minimal microsatellite size, default="<<paramd.MinMicrosate<<"\n"
        <<"       -s   <int>      minimal microsatellite size for distribution analysis, default="<<paramd.MinMicrosateForDis<<"\n"
        <<"       -w   <int>      maximal microsatellite size for distribution analysis, default="<<paramd.MaxMicrosateForDis<<"\n"

        <<"       -u   <int>      span size around window for extracting reads, default="<<paramd.DisSpan<<"\n"
        <<"       -b   <int>      threads number for parallel computing, default="<<paramd.numberThreads<<"\n"
        <<"       -x   <int>      output homopolymer only, 0: no; 1: yes, default="<<paramd.HomoOnly<<"\n"
        <<"       -y   <int>      output microsatellite only, 0: no; 1: yes, default="<<paramd.MicrosateOnly<<"\n"
        <<"       \n"
        <<"       -h   help\n\n"
        << std::endl;
    exit(1);
}


//add by yelab
void TrainUsage(void){
	std::cerr << "\nUsage:  msisensor baseline [options] \n\n"
		<< "       -d   <string>   homopolymer and microsatellite file\n"
		<< "       -i   <string>   configure files for building baseline (text file) \n"
		<< "            e.g.\n"
		<< "              case1\t/path/to/case1_sorted.bam\n"
		<< "              case2\t/path/to/case1_sorted.bam\n"
		<< "              case2\t/path/to/case1-sorted.bam\n"
		<< "       -o   <string>   output directory\n\n"

//		<< "       -e   <string>   bed file, optional\n"
//		<< "       -i   <double>   minimal comentropy threshold for somatic sites detection (just for tumor only data), default=" << paramd.comentropyThreshold << "\n"
		<< "       -c   <int>      coverage threshold for msi analysis, WXS: 20; WGS: 15, default=" << paramd.covCutoff << "\n"
//		<< "       -r   <string>   choose one region, format: 1:10000000-20000000\n"
//		<< "       -l   <int>      minimal homopolymer size, default=" << paramd.MininalHomoSize << "\n"
//		<< "       -p   <int>      minimal homopolymer size for distribution analysis, default=" << paramd.MininalHomoForDis << "\n"
//		<< "       -m   <int>      maximal homopolymer size for distribution analysis, default=" << paramd.MaxHomoSize << "\n"
//		<< "       -z   <string>   config file of training bam set, default=" << paramd.MaxHomoSize << "\n"
//		<< "       -a   <double>   cut off of pass possibility, default=" << paramd.MaxHomoSize << "\n"
//
//		<< "       -q   <int>      minimal microsatellite size, default=" << paramd.MinMicrosate << "\n"
//		<< "       -s   <int>      minimal microsatellite size for distribution analysis, default=" << paramd.MinMicrosateForDis << "\n"
//		<< "       -w   <int>      maximal microsatellite size for distribution analysis, default=" << paramd.MaxMicrosateForDis << "\n"
//
//		<< "       -u   <int>      span size around window for extracting reads, default=" << paramd.DisSpan << "\n"
//		<< "       -b   <int>      threads number for parallel computing, default=" << paramd.numberThreads << "\n"
//		<< "       -x   <int>      output homopolymer only, 0: no; 1: yes, default=" << paramd.HomoOnly << "\n"
//		<< "       -y   <int>      output microsatellite only, 0: no; 1: yes, default=" << paramd.MicrosateOnly << "\n"
		<< "       -0   <int>      output site have no read coverage, 1: no; 0: yes, default=" << paramd.outputzeroDis << "\n"
		<< "       \n"
		<< "       -h   help\n\n"
		<< std::endl;
	exit(1);
}


int dGetOptions(int rgc, char *rgv[]) {
    int i;
    for (i=1; i<rgc; i++) {
        if (rgv[i][0] != '-') return i;
        switch(rgv[i][1]) {
            case 'd': homoFile = rgv[++i]; break;
            case 'n': normalBam  = rgv[++i]; break;
            case 't': tumorBam  = rgv[++i]; break;
            case 'o': disFile  = rgv[++i]; break;
            case 'e': bedFile  = rgv[++i]; break;
            case 'r': one_region = rgv[++i]; break;
            case 'f': paramd.fdrThreshold  = atof(rgv[++i]); break;
//            case 'i': paramd.comentropyThreshold = atof(rgv[++i]); break;
            case 'c': paramd.covCutoff = atoi(rgv[++i]); break;
            case 'z': paramd.Normalization = atoi(rgv[++i]); break;
            case 'l': paramd.MininalHomoSize = atoi(rgv[++i]); break;
            case 'p': paramd.MininalHomoForDis = atoi(rgv[++i]); break;
            case 'u': paramd.DisSpan = atoi(rgv[++i]); break;
            case 'm': paramd.MaxHomoSize = atoi(rgv[++i]); break;
            case 'q': paramd.MinMicrosate = atoi(rgv[++i]); break;
            case 's': paramd.MinMicrosateForDis = atoi(rgv[++i]); break;
            case 'w': paramd.MaxMicrosateForDis = atoi(rgv[++i]); break;
            case 'b': paramd.numberThreads = atoi(rgv[++i]); break;
            case 'x': paramd.HomoOnly= atoi(rgv[++i]); break;
            case 'y': paramd.MicrosateOnly = atoi(rgv[++i]); break;
            case '0': paramd.outputzeroDis = atoi(rgv[++i]); break;
            break;
            case 'h':DisUsage();
            case '?':DisUsage();
        }
    }
    return i;
}


//add by YeLab
int dGetProOptions(int rgc, char *rgv[]) {
	int i;
	for (i = 1; i<rgc; i++) {
		if (rgv[i][0] != '-') return i;
		switch (rgv[i][1]) {
		case 'd': homoFile = rgv[++i]; break;
		//case 'n': normalBam = rgv[++i]; break;
		case 't': tumorBam = rgv[++i]; break;
		case 'o': disFile = rgv[++i]; break;
		case 'e': bedFile = rgv[++i]; break;
		case 'r': one_region = rgv[++i]; break;
	//	case 'f': paramd.fdrThreshold = atof(rgv[++i]); break;
		case 'i': paramd.hunterThreshold = atof(rgv[++i]); break;
		case 'c': paramd.covCutoff = atoi(rgv[++i]); break;
	//	case 'z': paramd.Normalization = atoi(rgv[++i]); break;
	//	case 'l': paramd.MininalHomoSize = atoi(rgv[++i]); break;
		case 'p': paramd.MininalHomoForDis = atoi(rgv[++i]); break;
		case 'u': paramd.DisSpan = atoi(rgv[++i]); break;
		case 'm': paramd.MaxHomoSize = atoi(rgv[++i]); break;
		case 'q': paramd.MinMicrosate = atoi(rgv[++i]); break;
		case 's': paramd.MinMicrosateForDis = atoi(rgv[++i]); break;
		case 'w': paramd.MaxMicrosateForDis = atoi(rgv[++i]); break;
		case 'b': paramd.numberThreads = atoi(rgv[++i]); break;
		case 'x': paramd.HomoOnly = atoi(rgv[++i]); break;
		case 'y': paramd.MicrosateOnly = atoi(rgv[++i]); break;
		case '0': paramd.outputzeroDis = atoi(rgv[++i]); break;
			break;
		case 'h':ProUsage();
		case '?':ProUsage();
		}
	}
	return i;
}
//add by yelab
int tGetOptions(int rgc, char *rgv[]) {
    int i;
    for (i=1; i<rgc; i++) {
        if (rgv[i][0] != '-') return i;
        switch(rgv[i][1]) {
            case 'd': homoFile = rgv[++i]; break;
            case 'i':TrainBamConfig = rgv[++i]; break;
            case 'o': disFile  = rgv[++i]; break;

            case 'e': bedFile  = rgv[++i]; break;
            case 'r': one_region = rgv[++i]; break;
			case '0': paramd.outputzeroDis=atoi(rgv[++i]); break;
			case 'a': paramd.NormalcovCutoff = atof(rgv[++i]); break;//!
            //case 'i': paramd.comentropyThreshold = atof(rgv[++i]); break;
            case 'c': paramd.covCutoff = atoi(rgv[++i]); break;
            case 'l': paramd.MininalHomoSize = atoi(rgv[++i]); break;
            case 'p': paramd.MininalHomoForDis = atoi(rgv[++i]); break;
            case 'u': paramd.DisSpan = atoi(rgv[++i]); break;
            case 'm': paramd.MaxHomoSize = atoi(rgv[++i]); break;
            case 'q': paramd.MinMicrosate = atoi(rgv[++i]); break;
            case 's': paramd.MinMicrosateForDis = atoi(rgv[++i]); break;
            case 'w': paramd.MaxMicrosateForDis = atoi(rgv[++i]); break;
            case 'b': paramd.numberThreads = atoi(rgv[++i]); break;
            case 'x': paramd.HomoOnly= atoi(rgv[++i]); break;
            case 'y': paramd.MicrosateOnly = atoi(rgv[++i]); break;
            break;
            case 'h':TrainUsage();
            case '?':TrainUsage();
        }
    }
    paramd.homoFile=homoFile;
    //add code for directory



    return i;
}




int dGetEntropyOptions(int rgc, char *rgv[]) {
    int i;
    for (i=1; i<rgc; i++) {
        if (rgv[i][0] != '-') return i;
        switch(rgv[i][1]) {
            case 'd': homoFile = rgv[++i]; break;
            case 't': tumorBam  = rgv[++i]; break;
            case 'o': disFile  = rgv[++i]; break;
            case 'e': bedFile  = rgv[++i]; break;
            case 'r': one_region = rgv[++i]; break;
            case 'i': paramd.comentropyThreshold = atof(rgv[++i]); break;
            case 'c': paramd.covCutoff = atoi(rgv[++i]); break;
            case 'l': paramd.MininalHomoSize = atoi(rgv[++i]); break;
            case 'p': paramd.MininalHomoForDis = atoi(rgv[++i]); break;
            case 'u': paramd.DisSpan = atoi(rgv[++i]); break;
            case 'm': paramd.MaxHomoSize = atoi(rgv[++i]); break;
            case 'q': paramd.MinMicrosate = atoi(rgv[++i]); break;
            case 's': paramd.MinMicrosateForDis = atoi(rgv[++i]); break;
            case 'w': paramd.MaxMicrosateForDis = atoi(rgv[++i]); break;
            case 'b': paramd.numberThreads = atoi(rgv[++i]); break;
            case 'x': paramd.HomoOnly= atoi(rgv[++i]); break;
            case 'y': paramd.MicrosateOnly = atoi(rgv[++i]); break;
            break;
            case 'h':EntropyDisUsage();
            case '?':EntropyDisUsage();
        }
    }
    return i;
}
//add by yelab : train normal data and find the p values for normal control
int TrainMsiP(int argc, char *argv[]) {
	paramd.train=true;
	if (argc == 1) TrainUsage();
	for (int i = 0; i<argc; i++) {
		std::cout << argv[i] << ' ';
	}
	Initial_Time();
	std::cout << "Start at:  " << Curr_Time() << std::endl;

	int noptions = tGetOptions(argc, argv);
	// process user defined region
	if (!one_region.empty()) {
		if (!polyscan.ParseOneRegion(one_region)) {
			std::cerr << "fatal error: Please give correct defined region format (-r) \n";
			exit(1);
		}
		polyscan.ifUserDefinedRegion = true;
	}
	else {
		polyscan.ifUserDefinedRegion = false;
	}
	// reading bed file if is exist
	finB.open(bedFile.c_str());
	if (finB) {
		std::cout << "loading bed regions ..." << std::endl;
		polyscan.LoadBeds(finB);
		polyscan.BedFilterorNot();
	}
	finB.close();

	loadFilepathFromConfig(TrainBamConfig, TrainBam);
	for (int j = 0; j < TrainBam.size(); j++) {
		polyscan.LoadBamn(TrainBam[j],TrainName[j]);
	}
	finH.open(homoFile.c_str());
	if (!finH) {
		std::cerr << "fatal error: failed to open homopolymer and microsatellites file\n";
		exit(1);
	}
	std::cout << "loading homopolymer and microsatellite sites ..." << std::endl;
	polyscan.LoadHomosAndMicrosates(finH);
	finH.close();
//	for (int j = 0; j < TrainBam.size(); j++) {
//		polyscan.SplitWindows();
//			//polyscan.TestWindows();
//			std::cout << "\nTotal loading windows:  " << polyscan.totalWindowsNum << " \n\n";
//			std::cout << "\nTotal loading homopolymer and microsatellites:  " << polyscan.totalHomosites << " \n\n";
//			;
//	}

	//polyscan.TestHomos();
	polyscan.SplitWindows();
	//polyscan.TestWindows();
	std::cout << "\nTotal loading windows:  " << polyscan.totalWindowsNum << " \n\n";
	std::cout << "\nTotal loading homopolymer and microsatellites:  " << polyscan.totalHomosites << " \n\n";

	// change code to one sample
	//polyscan.GetHomoDistribution(sample, disFile);
	// control distribution for tumor only input
	//    if (!normalBam.empty() && !tumorBam.empty()) {
	//        polyscan.GetHomoDistribution(sample, disFile);
	//    }
	const char* per=disFile.c_str();
	int h=mkdir(per,00755);
	const char* per1=(disFile+"/detail").c_str();
	int isCreate=mkdir(per1,00755);
//	const char* per=(prefix+"/"+totalBamNormals[j].sName).c_str();
	//		int isCreate = mkdir(per,00755);


//	 if( !isCreate )
//	   printf("create path:%s\n",disFile);
//	   else


	polyscan.GetNormalDistrubution(sample, disFile);
//	if (normalBam.empty() && !tumorBam.empty()) {
//		polyscan.GetNormalDistrubution(sample, disFile);
//	}

	std::cout << "\nTotal time consumed:  " << Cal_AllTime() << " secs\n\n";

	return 0;
}


int HomoAndMicrosateDisMsi(int argc, char *argv[]) {
    if (argc == 1) DisUsage();
    for (int i=0; i<argc; i++) {
        std::cout <<argv[i]<<' ';
    }
    Initial_Time();
    std::cout <<"Start at:  "<<Curr_Time() << std::endl;

    int noptions = dGetOptions(argc, argv);
    // process user defined region
    if (!one_region.empty()) {
        if (!polyscan.ParseOneRegion(one_region)) {
            std::cerr<<"fatal error: Please give correct defined region format (-r) \n";
            exit(1);
        }
        polyscan.ifUserDefinedRegion = true;
    } else {
        polyscan.ifUserDefinedRegion = false;
    }
    // reading bed file if is exist
    finB.open(bedFile.c_str());
    if (finB) {
        std::cout << "loading bed regions ..." << std::endl;
        polyscan.LoadBeds(finB);
        polyscan.BedFilterorNot();
    }
    finB.close();

    // load bam files
    //polyscan.LoadBams( normalBam, tumorBam );
//    if (!normalBam.empty() && !tumorBam.empty()) {
//        polyscan.LoadBams( normalBam, tumorBam );
//    }
//    // just for tumor only data
//    if (normalBam.empty() && !tumorBam.empty()) {
//        polyscan.LoadBam(tumorBam);
//    }
    polyscan.LoadBams( normalBam, tumorBam );
    // check homo/microsate file
    finH.open(homoFile.c_str());
    if (!finH) {
        std::cerr<<"fatal error: failed to open homopolymer and microsatellites file\n";
        exit(1);
    }
    std::cout << "loading homopolymer and microsatellite sites ..." << std::endl;
    polyscan.LoadHomosAndMicrosates(finH);
    finH.close();
    //polyscan.TestHomos();
    polyscan.SplitWindows();
    //polyscan.TestWindows();
    std::cout << "\nTotal loading windows:  " << polyscan.totalWindowsNum << " \n\n";
    std::cout << "\nTotal loading homopolymer and microsatellites:  " << polyscan.totalHomosites << " \n\n";

    // change code to one sample
    //polyscan.GetHomoDistribution(sample, disFile);
    // control distribution for tumor only input
    if (!normalBam.empty() && !tumorBam.empty()) {
    polyscan.GetHomoDistribution(sample, disFile);
    }
  //  if (normalBam.empty() && !tumorBam.empty()) {
    //    polyscan.GetHomoTumorDistribution(sample, disFile);
    //}

    std::cout << "\nTotal time consumed:  " << Cal_AllTime() << " secs\n\n";

    return 0;
}

int HomoAndMicrosateDisMsiEntropy(int argc, char *argv[]) {
    if (argc == 1) EntropyDisUsage();
    for (int i=0; i<argc; i++) {
        std::cout <<argv[i]<<' ';
    }
    Initial_Time();
    std::cout <<"Start at:  "<<Curr_Time() << std::endl;

    int noptions = dGetEntropyOptions(argc, argv);
    // process user defined region
    if (!one_region.empty()) {
        if (!polyscan.ParseOneRegion(one_region)) {
            std::cerr<<"fatal error: Please give correct defined region format (-r) \n";
            exit(1);
        }
        polyscan.ifUserDefinedRegion = true;
    } else {
        polyscan.ifUserDefinedRegion = false;
     }

    finB.open(bedFile.c_str());
    if (finB) {
        std::cout << "loading bed regions ..." << std::endl;
        polyscan.LoadBeds(finB);
        polyscan.BedFilterorNot();
    }
    finB.close();

    // load bam files
    // just for tumor only data
    polyscan.LoadBam(tumorBam);

    // check homo/microsate file
    finH.open(homoFile.c_str());
    if (!finH) {
        std::cerr<<"fatal error: failed to open homopolymer and microsatellites file\n";
        exit(1);
    }
    std::cout << "loading homopolymer and microsatellite sites ..." << std::endl;
    polyscan.LoadHomosAndMicrosates(finH);
    finH.close();
    //polyscan.TestHomos();
    polyscan.SplitWindows();
    //polyscan.TestWindows();
    std::cout << "\nTotal loading windows:  " << polyscan.totalWindowsNum << " \n\n";
    std::cout << "\nTotal loading homopolymer and microsatellites:  " << polyscan.totalHomosites << " \n\n";

    // change code to one sample
    polyscan.GetHomoTumorDistribution(sample, disFile);
    std::cout << "\nTotal time consumed:  " << Cal_AllTime() << " secs\n\n";
    return 0;
 }






// add by YeLab
int InitHunter(){
	//change some default value of parameter
	paramd.MininalHomoForDis=10;
	paramd.MaxHomoSize=20;
	paramd.MaxMicrosateForDis=20;
	return 0;
}

// add by Yelab for hunter
int HomoAndMicrosateDisMsiPro(int argc, char *argv[]) {
	InitHunter();
	if (argc == 1) ProUsage();
	for (int i = 0; i<argc; i++) {
		std::cout << argv[i] << ' ';
	}
	Initial_Time();
	std::cout << "Start at:  " << Curr_Time() << std::endl;

	int noptions = dGetProOptions(argc, argv);
	// process user defined region
	if (!one_region.empty()) {
		if (!polyscan.ParseOneRegion(one_region)) {
			std::cerr << "fatal error: Please give correct defined region format (-r) \n";
			exit(1);
		}
		polyscan.ifUserDefinedRegion = true;
	}
	else {
		polyscan.ifUserDefinedRegion = false;
	}
	// reading bed file if is exist
	finB.open(bedFile.c_str());
	if (finB) {
		std::cout << "loading bed regions ..." << std::endl;
		polyscan.LoadBeds(finB);
		polyscan.BedFilterorNot();
	}

	// load bam files
	//polyscan.LoadBams( normalBam, tumorBam );
	/*if (!normalBam.empty() && !tumorBam.empty()) {
		polyscan.LoadBams(normalBam, tumorBam);
	}*/
	// just for tumor only data
	if (normalBam.empty() && !tumorBam.empty()) {
		polyscan.LoadBam(tumorBam);
	}
	// check homo/microsate file
	finH.open(homoFile.c_str());
	if (!finH) {
		std::cerr << "fatal error: failed to open homopolymer and microsatellites file\n";
		exit(1);
	}
	std::cout << "loading homopolymer and microsatellite sites ..." << std::endl;
	polyscan.LoadHomosAndMicrosates(finH);
	finH.close();
	//polyscan.TestHomos();
	polyscan.SplitWindows();
	//polyscan.TestWindows();
	std::cout << "\nTotal loading windows:  " << polyscan.totalWindowsNum << " \n\n";
	std::cout << "\nTotal loading homopolymer and microsatellites:  " << polyscan.totalHomosites << " \n\n";

	// change code to one sample
	//polyscan.GetHomoDistribution(sample, disFile);
	// control distribution for tumor only input
	/*if (!normalBam.empty() && !tumorBam.empty()) {
		polyscan.GetHomoDistribution(sample, disFile);
	}*/
//	std::cout<< "normalBam.empty()" <<normalBam.empty()<<"\n";
//	std::cout<< "tumorBam.empty()" <<tumorBam.empty()<<"\n";

	if (normalBam.empty() && !tumorBam.empty()) {
		polyscan.GetHunterTumorDistribution(sample, disFile);

	}

	std::cout << "\nTotal time consumed:  " << Cal_AllTime() << " secs\n\n";

	return 0;
}



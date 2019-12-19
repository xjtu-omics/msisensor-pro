
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
#include <omp.h>

#include "param.h"
#include "refseq.h"
#include "scan.h"
#include "utilities.h"
#include "cmds.h"

Param param;
RefSeq ref;

std::string ref_file;
std::string out_file;
std::ifstream fin_d;
std::ofstream fout;

void ScanUsage(void) {
    std::cerr<<"\nUsage:  msisensor-pro scan [options] \n\n"
        <<"       -d   <string>   reference genome sequences file, *.fasta or *.fa format\n"
        <<"       -o   <string>   output homopolymers and microsatellites file\n\n"

        <<"       -l   <int>      minimal homopolymer(repeat unit length = 1) size, default="<<param.MininalHomoSize<<"\n"
        <<"       -c   <int>      context length, default="<<param.ContextLength<<"\n"
        <<"       -m   <int>      maximal homopolymer size, default="<<param.MaxHomoSize<<"\n"
        <<"       -s   <int>      maximal length of microsatellite, default="<<param.MaxMicrosate<<"\n"
        <<"       -r   <int>      minimal repeat times of microsatellite(repeat unit length>=2), default="<<param.Repeats<<"\n"
        <<"       -p   <int>      output homopolymer only, 0: no; 1: yes, default="<<param.HomoOnly<<"\n"
        <<"       \n"
        <<"       -h   help\n\n"
		<<"Function: \n"
		<<"   This module scan the reference genome to get microsatellites information. You need to input (-d) a reference file (*.fa or *.fasta), and you will get a microsatellites file (-o) for following analysis. If you use GRCh38.d1.vd1 , you can download the file on out github directly. \n\n"
		<<"Example:\n"
		<<"   msisensor-pro scan -d /path/to/reference.fa -o /path/to/reference.list\n\n"
		<<"Note:\n"
		<<"   This module inherits from msisensor.If you use it for your work, please cite:\n"
		<<"   Beifang Niu*, Kai Ye*, Qunyuan Zhang, Charles Lu, Mingchao Xie, Michael D. McLellan, Michael C. Wendl and Li Ding#.MSIsensor: microsatellite instability detection using paired tumor-normal sequence data. Bioinformatics 30, 1015–1016 (2014).\n \n"

        << std::endl;
    exit(1);
}

int mGetOptions(int rgc, char *rgv[]) {
    int i;
    for (i=1; i<rgc; i++) {
        if (rgv[i][0] != '-') return i;
        switch(rgv[i][1]) {
            case 'd': ref_file = rgv[++i]; break;
            case 'o': out_file = rgv[++i]; break;
            case 'l': param.MininalHomoSize = atoi(rgv[++i]); break;
            case 'c': param.ContextLength = atoi(rgv[++i]); break;
            case 'm': param.MaxHomoSize = atoi(rgv[++i]); break;
            case 's': param.MaxMicrosate = atoi(rgv[++i]); break;
            case 'r': param.Repeats = atoi(rgv[++i]); break;
            case 'p': param.HomoOnly = atoi(rgv[++i]); break;
            break;
            case 'h':ScanUsage();
            case '?':ScanUsage();    
        }
    }
    return i;
}

int HomoAndMicrosateScan(int argc, char *argv[]) {
    if (argc == 1) ScanUsage();
    for (int i=0; i<argc; i++) {
        std::cout <<argv[i]<<' ';
    }
    Initial_Time();
    std::cout <<"Start at:  "<<Curr_Time()<< std::endl;
    int noptions = mGetOptions(argc, argv);
    // check refseq file
    fin_d.open(ref_file.c_str());
    if (!fin_d) {
        std::cerr<<"fatal error: failed to open ref file\n";
        exit(1);
    }
    // output calling results
    fout.open(out_file.c_str());
    if (!fout) {
        std::cerr <<"failed to open file: "<<out_file<< std::endl;
        exit(1);
    }
    ref.PouroutHeader(fout);
    // reading refseq and count homo sites
    ref.ScanHomoAndMicrosate(fin_d);
    ref.PouroutBuffer(fout);
    std::cout << "\nTotal time consumed:  " << Cal_AllTime() << " secs\n\n";
    fout.close();
    fin_d.close();

    return 0;
}


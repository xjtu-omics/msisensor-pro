/*
 * cmds.cpp for MSIsensor
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
#include <iostream>
#include <string.h>
#include <stdio.h>
#include <omp.h>

#include "cmds.h"

#ifndef VERSION
#define VERSION "v0.1"
#endif

int usage(void) {
    std::cerr<<"\n\n"
    << "Program: msisensor-pro: Microsatellite Instability (MSI) detection using high-throughput sequencing data. \n"
	<< "         (Support tumor-normal paired samples and tumor-only samples) \n\n"
    << "Version: "<<VERSION<<"\n\n"
	 << "Usage:   msisensor-pro <command> [options]\n\n"
	    << "Key commands:\n\n"
	    << "\t scan\n"
		<< "\t   scan the reference genome to get microsatellites information\n\n"
		<< "\t baseline\n"
		<< "\t   build baseline for tumor only detection\n\n"
	    << "\t msi\n"
		<< "\t   evaluate MSI using paired tumor-normal sequencing data\n\n"
	//	<< " entropy\n"
	//	<< "    evaluate msi using msisensor-entropy (tumor only)\n\n"  // add by Niu lab
		<< "\t pro\n"
		<< "\t   evaluate MSI using single (tumor) sample sequencing data \n"//add by Peng Jia (2018.10.9)
	    << "\n\n"

	<< "Best Practices:\n"
	<< "\t(a) For tumor only samples:\n"
	<< "\t\t1. scan\n"
	<< "\t\t      >> msisensor-pro scan -d /path/to/reference.fa -o /path/to/reference.site\n"

	<< "\t\t2. baseline \n"
	<< "\t\t      >> msisensor-pro baseline -d /path/to/reference.site -i /path/to/configure.txt -o /path/to/baseline/directory\n"

	<< "\t\t3. pro \n"
	<< "\t\t      >> msisensor-pro pro -d /path/to/baseline/directory/reference_baseline.site -t /path/to/tumor/case1_sorted.bam -o /path/to/output\n\n"
	<< "\t(b) For tumor-normal paired samples:\n"
	<< "\t\t1. scan\n"
	<< "\t\t      >> msisensor-pro scan -d /path/to/reference.fa -o /path/to/reference.site\n"
	<< "\t\t2. msi \n"
	<< "\t\t      >> msisensor-pro msi -d /path/to/reference.site -n /path/to/case1_normal_sorted.bam -t /path/to/case1_tumor_sorted.bam -o /path/to/output\n\n"

	<< "Notes:We offer the scan result of GRCh38.d1.vd1 on our github \n\n"
    << "\n\n";
    return 1; 
}

int main(int argc, char **argv) {
    try {
        if (argc < 2) {
            return usage();
        }
        if (strcmp(argv[1], "scan") == 0) {
            // scan homopolymer and microsate
            HomoAndMicrosateScan(argc-1, argv+1);
            return 0;
        }
        else if (strcmp(argv[1], "baseline") == 0) {
            // distribution && msi scoring analysis
            TrainMsiP(argc-1, argv+1);
            return 0;

        }
        else if (strcmp(argv[1], "msi") == 0) {
            // distribution && msi scoring analysis 
            HomoAndMicrosateDisMsi(argc-1, argv+1);
            return 0;
		
        }
//        else if (strcmp(argv[1], "entropy") == 0) {
//              // distribution && msi scoring analysis for tumor only
//             HomoAndMicrosateDisMsiEntropy(argc-1, argv+1);
//             return 0;
//        }
		else if (strcmp(argv[1], "pro") == 0) { // this condition was added by Peng Jia (2018.10.9)
												   // distribution && msi detection using MSIsensor-pro

			HomoAndMicrosateDisMsiPro(argc - 1, argv + 1);
			return 0;
		}
		else {
            std::cerr<<"ERROR: unrecognized command "<<argv[1]<<"\n";
            return 2;
        }
    } catch (const char *e) {
        std::cerr << e << std::endl;
        return 3;
    }

    return 0;
}


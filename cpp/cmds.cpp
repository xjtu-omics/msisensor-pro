
// System header files
#include <iostream>
#include <string.h>
#include <stdio.h>
#include <omp.h>

#include "cmds.h"

#ifndef VERSION
//#define VERSION "v0.1.0"
//#define VERSION "v1.0.a"
#define VERSION "v1.0.2"  
#endif

int usage(void) {
	std::cerr << "\n\n"
			<< "Program: msisensor-pro: Microsatellite Instability (MSI) detection using high-throughput sequencing data. \n"
			<< "         (Support tumor-normal paired samples and tumor-only samples) \n\n"
			<< "Version: " << VERSION << "\n\n"
			<< "Usage:   msisensor-pro <command> [options]\n\n"
			<< "Key Commands:\n\n" << "\t scan\n"
			<< "\t   scan the reference genome to get microsatellites information\n\n"
			<< "\t baseline\n"
			<< "\t   build baseline for tumor only detection\n\n" << "\t msi\n"
			<< "\t   evaluate MSI using paired tumor-normal sequencing data\n\n"
			<< "\t pro\n"
			<< "\t   evaluate MSI using single (tumor) sample sequencing data \n"
			//add by Peng Jia (2018.10.9)
			<< "\n\n"

			<< "Best Practices:\n" << "\t(a) For tumor only samples:\n"
			<< "\t\t1. scan\n"
			<< "\t\t      >> msisensor-pro scan -d /path/to/reference.fa -o /path/to/reference.list\n"

			<< "\t\t2. baseline \n"
			<< "\t\t      >> msisensor-pro baseline -d /path/to/reference.list -i /path/to/configure.txt -o /path/to/baseline/directory\n"

			<< "\t\t3. pro \n"
			<< "\t\t      >> msisensor-pro pro -d /path/to/baseline/directory/reference.list_baseline -t /path/to/case1_tumor_sorted.bam -o /path/to/case1_output\n\n"
			<< "\t(b) For tumor-normal paired samples:\n" << "\t\t1. scan\n"
			<< "\t\t      >> msisensor-pro scan -d /path/to/reference.fa -o /path/to/reference.list\n"
			<< "\t\t2. msi \n"
			<< "\t\t      >> msisensor-pro msi -d /path/to/reference.list -n /path/to/case1_normal_sorted.bam -t /path/to/case1_tumor_sorted.bam -o /path/to/case1_output\n\n"

			<< "Notes:\n"
			<< "    1. If your analysis are based on reference GRCh38.d1.vd1, you can ignore the scan step by downloading the microsatellites information on our github directly. \n\n"
			<< "    2. If you don't have normal samples to build baseline(baseline step for tumor only sample detection), you can download the microsatellites information with baseline on our github or use -i option in pro module to set a hard cutoff directly.\n\n"

			<< "    If you have any questions, please contact with Peng Jia (pengjia@stu.xjtu.edu.cn) or Kai Ye (kaiye@xjtu.edu.cn) .\n"
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
			HomoAndMicrosateScan(argc - 1, argv + 1);
			return 0;
		} else if (strcmp(argv[1], "baseline") == 0) {
			// distribution && msi scoring analysis
			TrainMsiP(argc - 1, argv + 1);
			return 0;

		} else if (strcmp(argv[1], "msi") == 0) {
			// distribution && msi scoring analysis
			HomoAndMicrosateDisMsi(argc - 1, argv + 1);
			return 0;

		} else if ((strcmp(argv[1], "-V") == 0) || (strcmp(argv[1], "-v") == 0)
				|| (strcmp(argv[1], "--version") == 0)) {
			// distribution && msi scoring analysis
			std::cerr << "\n" << "Version: " << VERSION << "\n\n";
			return 0;
		}
		else if (strcmp(argv[1], "pro") == 0) { // this condition was added by Peng Jia (2018.10.9)
												// distribution && msi detection using MSIsensor-pro

			HomoAndMicrosateDisMsiPro(argc - 1, argv + 1);
			return 0;
		} else if (strcmp(argv[1], "-h") || strcmp(argv[1], "--help")) {
			usage();
			return 2;
		} else {
			std::cerr << "ERROR: unrecognized command " << argv[1] << "\n";
			return 2;
		}
	} catch (const char *e) {
		std::cerr << e << std::endl;
		return 3;
	}

	return 0;
}


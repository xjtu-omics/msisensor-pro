
#ifndef _STRUCTS_H_
#define _STRUCTS_H_

#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <map>

// user defined region
struct UserDefinedRegion {
	UserDefinedRegion() :
			chr(""), start(-1), end(-1) {
		//xxxx
	}
	std::string chr;
	int start;
	int end;
};

// bed region 
struct BedRegion {
	BedRegion() :
			start(0), end(0) {
		//xxxx
	}
	int start;
	int end;
};

// bed regions located on one chromosome
struct BedChr {
	BedChr() :
			chr("") {
		//xxx
	}
	std::string chr;
	std::vector<BedRegion> regions_list;
};

// genotype by Kai
struct Genotype {
	Genotype() {
		//xxx
		GT[0] = -1;
		GT[1] = -1;
		WithGenotype = false;
	}
	short GT[2];
	bool WithGenotype;
};

// Bam file pairs
struct BamPairs {
	BamPairs() :
			sName(""), normal_bam(""), tumor_bam("") {
		//xxx
	}
	std::string sName;
	// bam files
	std::string normal_bam;
	std::string tumor_bam;

};

// Bam file tumors
struct BamTumors {
	BamTumors() :
			sName(""), tumor_bam("") {
		//xxx
	}
	std::string sName;
	// bam files
	std::string tumor_bam;

};
// Bam file normals
struct BamNormals {
	BamNormals() :
			sName(""), normal_bam("") {
		//xxx
	}
	std::string sName;
	// bam files
	std::string normal_bam;

};

#endif //_STRUCTS_H_


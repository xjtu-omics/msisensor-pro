
#ifndef _POLYSCAN_H_
#define _POLYSCAN_H_

#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>

#include "param.h"
#include "structs.h"

#include "homo.h"
#include "window.h"
#include "sample.h"
//#include "train.h"

class PolyScan {
public:

	PolyScan();
	~PolyScan();

	// user defined region
	UserDefinedRegion region_one;
	bool ifUserDefinedRegion;
	void eliminate(const char ch, std::string & str);
	bool ParseOneRegion(const std::string &region);

	// read bed regions
	bool ifUserDefinedBed;
	std::map<std::string, bit16_t> chrMaptoIndex;
	std::vector<BedChr> beds;
	void LoadBeds(std::ifstream &fin);
	void BedFilterorNot();

	// load bam list file
	std::vector<BamPairs> totalBamPairs;
	std::vector<BamTumors> totalBamTumors;
	std::vector<BamNormals> totalBamNormals; //add by yelab

	//void LoadBams(std::ifstream &fin);
	void LoadBams(const std::string &nBam, const std::string &tBam);
	void LoadBam(const std::string &tBam);
	void LoadBamn(const std::string &bam, const std::string &Name); //add by yelab
	unsigned int totalBamPairsNum;
	unsigned int totalBamTumorsNum;
	unsigned int totalBamNormalsNum; //add by yelab
	// load homos and microsatellites
	unsigned long totalHomosites;
	//std::vector< HomoSite * > totalSites;
	std::vector<HomoSite> totalSites;
	bool LoadHomosAndMicrosates(std::ifstream &fin);
	void TestHomos();

	std::vector<HomoSite> homosBuffer;

	// windows
	std::vector<Window> totalWindows;
	void SplitWindows();
	void TestWindows();
	unsigned long totalWindowsNum;

	// distribution
	void InithializeDistributions();
	void outputDistributions();
	void releaseDistributions();
	void GetNormalDistrubution(Sample &oneSample, const std::string &prefix);
	void GetHomoDistribution(Sample &oneSample, const std::string &prefix);
	void GetHomoTumorDistribution(Sample &oneSample, const std::string &prefix);
	void GetHunterTumorDistribution(Sample &oneSample,
			const std::string &prefix);

protected:

	// xxxxxx
	// xxxx

};

#endif //_POLYSCAN_H_


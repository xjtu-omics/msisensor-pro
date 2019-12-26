
#ifndef _SAMPLE_H_
#define _SAMPLE_H_

#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>

#include "somatic.h"

// sample
class Sample {
public:
	Sample();
	~Sample();

	std::string outputPrefix;

	std::ofstream output;
	//std::ofstream outputPSomatic;
	std::ofstream outputSomatic;
	std::ofstream outputAll;
	std::ofstream outputGermline;
	std::ofstream outputDistribution;
	std::ofstream outputTrain;

	unsigned numberOfSites;

	unsigned precisionNumS;
	unsigned precisionNumL;

	unsigned numberOfDataPoints;
	unsigned numberOfMsiDataPoints;
	unsigned numberOftotalSites;

	// container for FDR
	std::vector<SomaticSite> totalSomaticSites;

	void iniOutput(const std::string &gavePrefix);
	void iniTumorDisOutput(const std::string &gavePrefix);
	void hunterIniTumorDisOutput(const std::string &gavePrefix); //YeLab
	void trainIniNormalDisOutput(const std::string &gavePrefix); //YeLab
	void pourOutMsiScore();
	void closeOutStream();
	void closeOutStreamTrain();
	void calculateFDR();
	void pourOutSomaticFDR();
	void VerboseInfo();

protected:
	// xxx
};

#endif //_SAMPLE_H_


#ifndef _HOMO_H_
#define _HOMO_H_

#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>

#include "param.h"
#include "sample.h"
#include "chi.h"

// homopolymer site
class HomoSite {
public:
	HomoSite();
	~HomoSite();
	// length of homo
	bit8_t typeLen;
	// homo or microsate
	// A/C/G/T/AC/AGC
	bit64_t homoType;
	bit16_t length;
	bit64_t frontKmer;
	bit64_t endKmer;

	std::string chr;
	std::string site;
	// readable
	std::string transfer;

	std::string bases;
	std::string fbases;
	std::string ebases;

	// location
	int location;
	// added low and high cutoff
	// for filtering the reads
	// without the window
	//
	int lowcut;
	int highcut;
	double thres;

	// distribution
	unsigned short **normalDis;
	unsigned short **tumorDis;

	////// genotyping //////
	unsigned normalCov;
	unsigned tumorCov;
	bool withSufCov;
	bool normalWithSufCov;
	bool somatic;
	bool withGenotype;
	double dif;
	double pValue;
	double comentropy;
	double hunterValueU;
	double hunterValueV;
	int genotype[2];
	////////////////////////

	inline void InitType() {
		genotype[0] = genotype[1] = -2;
	}
	;

	void TransferString();
	void InitialDis();
	void InitialTumorDis();
	void InitialNormalDis();    //add by yelab
	void OutputDis();
	void OutputTumorDis();
	void ReleaseMemory();
	void ReleaseTumorMemory();
	void ReleaseNormalMemory();
	//void PouroutDis(std::ofstream &fout);
	void PouroutDis(Sample &sample);
	void PourTumoroutDis(Sample &sample);
	int outDislabel(Sample &sample);
	int outDislabelTumorOnly(Sample &sample);
	void DisGenotyping(Sample &sample);
	void DisTumorSomatic(Sample &sample);
	void HunterDisTumorSomatic(Sample &sample);
	void HunterTrain(Sample &sample);    //add by yelab
	//// genotyping ///
	void BoolsInitial();
	double DistanceBetweenTwo(unsigned short * FirstOriginal,
			unsigned short * SecondOriginal);
	double Comentropy(unsigned short * tumorDis, unsigned int dispots);
	std::vector<double> Hunterp(unsigned short * tumorDis, unsigned int dispots,
			unsigned int reflen);
	void ComputeGenotype(unsigned short * NormalReadCount);

protected:
	// xxx
};

#endif //_HOMO_H_


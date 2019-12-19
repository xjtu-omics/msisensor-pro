

#ifndef _WINDOW_H_
#define _WINDOW_H_

#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>

#include "param.h"
#include "homo.h"
#include "bamreader.h"
#include "sample.h"

// homopolymer site
class Window {
public:
    Window();
    ~Window();

    int _start;
    int _end;
    unsigned short _siteCount;
    std::string _chr;
    HomoSite *_startSite;
    HomoSite *_endSite;

    void InitialDisW(); 
    void InitialTumorDisW();
	void InitialNormalDisW();//add by yelab
    void OutputDisW();
    void OutputTumorDisW();
    void ClearDis();
    void ClearTumorDis();
	void ClearNormalDis();//add by yelab
    void ChangeStart(); 
    void GetDistribution(std::vector <SPLIT_READ> &readsInWindow);
    void GetTumorDistribution(std::vector <SPLIT_READ> &readsInWindow);
	void GetNormalDistribution(std::vector <SPLIT_READ> &readsInWindow);//add by yelab
    void LoadReads(std::vector <SPLIT_READ> &readsInWindow, const std::string bam);
    void ScanReads(const std::vector <SPLIT_READ> &readsInWindow, unsigned short bamIndex, bool isTumor);
    void ReverseComplement(std::string &theWord);
    unsigned short DoOneRead(const std::string &oneRead, const HomoSite *p);
    void PouroutDisW(Sample &oneSample);
    void PourTumoroutDisW(Sample &oneSample);
	void PouroutNormalH(Sample &oneSample);
    void DisGenotypingW(Sample &oneSample);
    void PouroutTumorSomatic(Sample &oneSample);
	void PouroutTumorSomaticH(Sample &oneSample);

protected:
    //xxxxxx
};

#endif //_WINDOW_H_


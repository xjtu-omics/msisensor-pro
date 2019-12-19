
#ifndef _SOMATIC_H_
#define _SOMATIC_H_

#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>

// homopolymer site
class SomaticSite {
public:
    SomaticSite();
    ~SomaticSite();

    // homo or microsate 
    // A/C/G/T/AC/AGC
    //
    std::string chr;

    // location
    int location;

    // repeat times 
    unsigned short length;

    // homo or microsate
    // content 
    std::string bases;
    // front kmer
    std::string fbases;
    // tail kmer
    std::string ebases;
    // difference between
    // normal and tumor 
    double diff;
    double pValue;
    double FDR;
    unsigned short rank;
    bool somatic;
    
    // output content
    void PourOut(); 
    // sorting based on p-value
    // 
    bool operator < (const SomaticSite& rhs) const { return pValue < rhs.pValue; }

    protected:
        // xxx
};

#endif //_SOMATIC_H_


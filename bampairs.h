

#ifndef _BAMPAIRS_H_
#define _BAMPAIRS_H_

#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>

#include "param.h"

// Bam file pairs
class BamPairs {
public:
    BamPairs();
    ~BamPairs();
    // total pairs
    unsigned int totalParis;
    // bam files
    std::string normal_bam;
    std::string tumor_bam;

protected:
    //xxx
};

#endif //_BAMPAIRS_H_


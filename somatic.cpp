
#include <iostream>
#include <sstream>
#include <bitset>
#include <omp.h>
#include <cmath>
#include "somatic.h"

SomaticSite::SomaticSite()
    : chr("")
    , length(0)
    , location(0)
    , bases("")
    , fbases("")
    , ebases("")
    , diff( 0.0 )
    , pValue( 1.0 )
    , somatic( false )
    , FDR( 1.0 )
    , rank( 1 )
{
    //xxxxxxxx
};


SomaticSite::~SomaticSite() {
    // xxxxx
};

// PourOut values
void SomaticSite::PourOut() {
    std::cerr << chr << "\t"
              << location << "\t"
              << bases << "\t"
              << length << "\t"
              << fbases << "\t"
              << bases << "\t"
              << ebases << "\t"
              << diff << "\t"
              << pValue << "\n";

};



#include <iostream>
#include <sstream>
#include <bitset>
#include <omp.h>
#include <iomanip>

#include "sample.h"
#include "param.h"

extern Param paramd;

Sample::Sample() :
        outputPrefix("test"), output(NULL)
        //, outputPSomatic( NULL )
        , outputSomatic(NULL), outputAll(NULL) //YeLab
        , outputGermline(NULL), outputDistribution(NULL), outputTrain(
        NULL), numberOfSites(0), precisionNumS(2), precisionNumL(5), numberOfDataPoints(
        0), numberOfMsiDataPoints(0), numberOftotalSites(0) {
    //xxxx
    output.precision(precisionNumS);
    // outputPSomatic.precision( precisionNumL );
    outputSomatic.precision(precisionNumL);
    outputGermline.precision(precisionNumL);
    outputDistribution.precision(precisionNumL);

};

Sample::~Sample() {
    // xxxxx
};

void Sample::iniOutput(const std::string &gavePrefix) {
    if (!gavePrefix.empty()) {
        outputPrefix = gavePrefix;
    }
    // init pour out result files
    output.open(outputPrefix.c_str());
    // outputPSomatic.open( (outputPrefix + "_p_somatic").c_str() );
    outputSomatic.open((outputPrefix + "_unstable").c_str());
    outputGermline.open((outputPrefix + "_all").c_str());
    outputDistribution.open((outputPrefix + "_dis").c_str());

    //if ( !output || !outputPSomatic || !outputSomatic || !outputGermline || !outputDistribution ) {
    if (!output || !outputSomatic || !outputGermline || !outputDistribution) {
        std::cerr << "failed to open output files to write !" << std::endl;
        exit(1);
    }
    outputSomatic << "chromosome" << "\t" << "location" << "\t"
                  << "left_flank_bases" << "\t" << "repeat_times" << "\t"
                  << "repeat_unit_bases" << "\t" << "right_flank_bases" << "\t"
                  << "difference" << "\t" << "P_value" << "\t" << "FDR" << "\t"
                  << "rank" << "\n";
//	outputGermline << "chromosome" << "\t" << "location" << "\t"
//			<< "left_flank_bases" << "\t" << "repeat_times" << "\t"
//			<< "repeat_unit_bases" << "\t" << "right_flank_bases" << "\t"
//			<< "genotype" << "\n";
    outputGermline << "chromosome" << "\t" << "location" << "\t"
                   << "left_flank_bases" << "\t" << "repeat_times" << "\t"
                   << "repeat_unit_bases" << "\t" << "right_flank_bases" << "\t"
                   << "\t" << "difference" << "\t" << "P_value" << "\t" << "FDR" << "\n";

}

void Sample::iniTumorDisOutput(const std::string &gavePrefix) {
    if (!gavePrefix.empty()) {
        outputPrefix = gavePrefix;
    }
    // init pour out result files
    output.open(outputPrefix.c_str());
    outputSomatic.open((outputPrefix + "_somatic").c_str());
    outputDistribution.open((outputPrefix + "_dis").c_str());

    //if ( !output || !outputPSomatic || !outputSomatic || !outputGermline || !outputDistribution ) {
    if (!output || !outputDistribution || !outputSomatic) {
        std::cerr << "failed to open output files to write !" << std::endl;
        exit(1);
    }
    outputSomatic << "chromosome" << "\t" << "location" << "\t"
                  << "left_flank_bases" << "\t" << "repeat_times" << "\t"
                  << "repeat_unit_bases" << "\t" << "right_flank_bases" << "\t"
                  << "entropy" << "\n";

}

// add by YeLab, for hunter
void Sample::hunterIniTumorDisOutput(const std::string &gavePrefix) {
    if (!gavePrefix.empty()) {
        outputPrefix = gavePrefix;
    }
    // init pour out result files
    output.open(outputPrefix.c_str());
    outputSomatic.open((outputPrefix + "_unstable").c_str());
    outputDistribution.open((outputPrefix + "_dis").c_str());
    outputAll.open((outputPrefix + "_all").c_str());

    //if ( !output || !outputPSomatic || !outputSomatic || !outputGermline || !outputDistribution ) {
    if (!output || !outputDistribution || !outputSomatic || !outputAll) {
        std::cerr << "failed to open output files to write !" << std::endl;
        exit(1);
    }
    outputSomatic << "chromosome" << "\t" << "location" << "\t"
                  << "left_flank_bases" << "\t" << "repeat_times" << "\t"
                  << "repeat_unit_bases" << "\t" << "right_flank_bases" << "\t"
                  << "pro_p" << "\t" << "pro_q" << "\t" << "CovReads" << "\t"
                  << "threshold" << "\n";
    outputAll << "chromosome" << "\t" << "location" << "\t"
              << "left_flank_bases" << "\t" << "repeat_times" << "\t"
              << "repeat_unit_bases" << "\t" << "right_flank_bases" << "\t"
              << "pro_p" << "\t" << "pro_q" << "\t" << "CovReads" << "\t"
              << "threshold" << "\n";
} //YeLab

void Sample::trainIniNormalDisOutput(const std::string &gavePrefix) { //Yelab
    if (!gavePrefix.empty()) {
        outputPrefix = gavePrefix;
    }

    // init pour out result files
//		output.open(outputPrefix.c_str());
//		outputSomatic.open((outputPrefix + "_somatic").c_str());
//		outputDistribution.open((outputPrefix + "_dis").c_str());
    outputAll.open((outputPrefix + "_all").c_str());

    //if ( !output || !outputPSomatic || !outputSomatic || !outputGermline || !outputDistribution ) {
    if (!outputAll) {
        std::cerr << "failed to open output files to write !" << std::endl;
        exit(1);
    }
    outputAll << "chromosome" << "\t" << "location" << "\t"
              << "left_flank_bases" << "\t" << "repeat_times" << "\t"
              << "repeat_unit_bases" << "\t" << "right_flank_bases" << "\t"
              << "covReads" << "\t" << "pro_p" << "\t" << "pro_q" << "\n";

} //YeLab

void Sample::pourOutMsiScore() {
    output << "Total_Number_of_Sites\tNumber_of_Unstable_Sites\t%" << std::endl;
    if (numberOfDataPoints != 0) {
        output << numberOfDataPoints << "\t" << numberOfMsiDataPoints << "\t"
               << std::fixed
               << (numberOfMsiDataPoints / (double) numberOfDataPoints) * 100.0
               << std::endl;
    } else {
        output << numberOfDataPoints << "\t" << numberOfMsiDataPoints << "\t"
               << std::setprecision(2) << std::fixed << 0.00 << std::endl;
    }
}

void Sample::closeOutStream() {
    output.close();
    //outputPSomatic.close();
    outputSomatic.close();
    outputGermline.close();
    outputDistribution.close();
}

void Sample::closeOutStreamTrain() {
//    output.close();
    //outputPSomatic.close();
//    outputSomatic.close();
    outputAll.close();
//    outputDistribution.close();
}

// FDR determination
void Sample::calculateFDR() {
    // sorting by p_value
    sort(totalSomaticSites.begin(), totalSomaticSites.end());
    unsigned short rank = 1;
    // FDR calculation
    for (std::vector<SomaticSite>::iterator _it = totalSomaticSites.begin();
         _it != totalSomaticSites.end(); ++_it) {
        //_it->PourOut();
        _it->FDR = _it->pValue * numberOfDataPoints / rank;
        if (_it->FDR > paramd.fdrThreshold) {
            rank++;
            continue;
        } else {
            _it->rank = rank;
            _it->somatic = true;
            numberOfMsiDataPoints++;
            rank++;
        }
    }
}

// report somatics && FDR
void Sample::pourOutSomaticFDR() {
    for (std::vector<SomaticSite>::iterator _it = totalSomaticSites.begin();
         _it != totalSomaticSites.end(); ++_it) {
        outputGermline << _it->chr << "\t" << _it->location << "\t"
                       << _it->fbases << "\t" << _it->length << "\t" << _it->bases
                       << "\t" << _it->ebases      //  << "\t" << _it->genotype1 << "|" << _it->genotype2
                       << "\t" << _it->diff << "\t"
                       << _it->pValue << "\t" << _it->FDR << "\n";
        if (!_it->somatic)
            continue;
        outputSomatic << _it->chr << "\t" << _it->location << "\t"
                      << _it->fbases << "\t" << _it->length << "\t" << _it->bases
                      << "\t" << _it->ebases << "\t" << _it->diff << "\t"
                      << _it->pValue << "\t" << _it->FDR << "\t" << _it->rank << "\n";
    }
}

// verbose 
void Sample::VerboseInfo() {
    std::cerr << "\n*** Summary information ***\n\n"
              << "Number of total sites: " << numberOftotalSites << "\n"
              << "Number of sites with enough coverage: " << numberOfDataPoints
              << "\n" << "Number of MSI sites: " << numberOfMsiDataPoints << "\n";
}


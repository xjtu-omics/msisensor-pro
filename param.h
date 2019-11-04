

#ifndef _PARAM_H_
#define _PARAM_H_

#include <string>

// maximal length of 
// microsatellite
#define MAX_MICROSATE_LEN 8
#define MAX_FLANK_REGION 8
#define MAX_TRANSFER_LINE_LENGTH 100
#define MAX_WINDOW 1000000
#define MAX_SPAN_SIZE 1000
#define MAX_READ_LENGTH 110

typedef unsigned char bit8_t;
typedef unsigned short bit16_t;
typedef unsigned bit32_t;
typedef unsigned long long bit64_t;

typedef bit32_t ref_id_t;
typedef bit32_t ref_loc_t;

class Param {
public:
    Param();
    ~Param();

    int max_dbseq_size; 
    int append_dbseq_size; 
    int bufSize;
    // Homo sites
    int MininalHomoSize;
    int ContextLength;
    // Microsate
    unsigned int MaxMicrosate;
    unsigned int Repeats;

    unsigned int MinMicrosate;
    unsigned int MinMicrosateForDis;
    unsigned int MaxMicrosateForDis;
    std::string homoFile;

    // filtering
    int HomoOnly;
    int MicrosateOnly;
    int ncpu;
    int chains;
    // Homo sites
    int MininalHomoForDis;
    int MaxHomoSize;
    int SpanSize;
    int DisSpan;
    // coverage normalization

    int Normalization;

    // output 0 dis , add by Yelab
    int outputzeroDis;

    // Thread number
    unsigned int numberThreads;
    // statistic var
    unsigned int s_dispots; 
    unsigned int PercentPairs;
    unsigned int PercentPairsNumber;
    unsigned int HomoCoverage;
    // window size
    unsigned int windowSize;

    // genotyping 
    unsigned int covCutoff;
    double fdrThreshold;
    double comentropyThreshold;

	double hunterThreshold;//add by YeLab
	double NormalcovCutoff;//add by YeLab
	double sampleRatio; //add by YeLab
	bool train;
	bool pro;
	bool hard;
};

#endif //_PARAM_H_


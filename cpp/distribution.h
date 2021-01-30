
#ifndef _DISTRIBUTION_H_
#define _DISTRIBUTION_H_

#include<vector>
#include<string>
#include<iostream>
#include<fstream>
#include "param.h"

// Samtools header files
#include "htslib/khash.h"
#include "htslib/sam.h"

//using namespace std;

#ifdef __cplusplus
extern "C" {
std::string abs_path(std::string path);

}
#endif

#endif //_DISTRIBUTION_H_


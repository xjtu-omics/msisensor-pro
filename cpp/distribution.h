


#ifndef _DISTRIBUTION_H_
#define _DISTRIBUTION_H_

#include<vector>
#include<string>
#include<iostream>
#include<fstream>
#include "param.h"

// Samtools header files
#include "khash.h"
#include "sam.h"

//using namespace std;

#ifdef __cplusplus
extern "C" 
{
    int32_t bam_get_tid(const bam_header_t * header, const char *seq_name);
    int32_t bam_aux2i(const uint8_t * s);
    void bam_init_header_hash(bam_header_t * header);
    std::string abs_path(std::string path);

}
#endif

#endif //_DISTRIBUTION_H_


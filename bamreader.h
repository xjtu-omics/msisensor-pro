
#ifndef BAMREADER_H
#define	BAMREADER_H

#include "map"
#include "vector"

// Samtools header files
#include "bam.h"
#include "sam.h"
#include "kstring.h"
#include "kseq.h"
#include "khash.h"
#include "ksort.h"

const unsigned g_SpacerBeforeAfter = 10000000;

struct flags_hit {
    flags_hit() {
        mapped = false;
        unique = false;
        sw = false;
        edits = 0;
        suboptimal = false;
    }
    bool mapped;
    bool unique;
    bool sw;
    int edits;
    bool suboptimal;
};

/* initial from kai
// struct of grab reads
struct SPLIT_READ
{
    SPLIT_READ() 
    {        
        FragName = "";
        Name = "";
	Mapped = true;
        ReadSeq = "";
        MatchedD = 0;
        MatchedRelPos = 0;
        MS = 0;
        Tag = "";
    }

    std::string FragName;
    std::string Name;
    bool Mapped;
    char MatchedD; // rename AnchorStrand?
    unsigned int MatchedRelPos;
    short MS; // rename MappingQuality ?
    std::string Tag; // rename SampleName ?
    friend std::ostream& operator<<(std::ostream& os, const SPLIT_READ& splitRead);
    std::string ReadSeq;
};
*/

struct SPLIT_READ {
    SPLIT_READ() {        
        ReadSeq = "";
        Mapped = true;
        MatchedRelPos = 0;
        
    }
    std::string ReadSeq;
    bool Mapped;
    unsigned int MatchedRelPos;
};

struct SupportPerSample {
    int NumPlus;
    int NumMinus;
    int NumUPlus;
    int NumUMinus;
    SupportPerSample() { 
        NumPlus = 0; 
        NumMinus = 0; 
        NumUPlus = 0; 
        NumUMinus = 0; 
    }
};

// read processing
struct HomoSiteforBam {
    HomoSiteforBam(){
        length = 0;
        type = '\0';
        front_kmer = "";
        end_kmer = "";
    }
    unsigned int length;
    std::string front_kmer;
    std::string end_kmer;
    char type;
    std::vector< unsigned int > distribution;
};

void build_record( const bam1_t * mapped_read, 
                   void *data, 
                   const flags_hit *flag_current_read);

bool ReadInBamReads( const char *bam_path, 
                     const std::string & FragName, 
                     unsigned start, 
                     unsigned end, 
                     std::vector < SPLIT_READ > & AllReads, 
                     std::string Tag);

#endif /* READER_H */


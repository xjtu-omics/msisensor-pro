
#include <iostream>
#include <sstream>
#include <bitset>
#include <omp.h>

extern Param paramd;

BamTumors::BamTumors() :
		_start(0), _end(0), _chr(""), _startSite(NULL), _endSite(NULL) {
	//xxxx

}
;

BamTumors::~BamTumors() {
	//xxx
}
;


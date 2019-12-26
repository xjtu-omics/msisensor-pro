#include <iostream>
#include <sstream>
#include <bitset>
#include <omp.h>

#include "bampairs.h"

extern Param paramd;

Bampairs::Bampairs() :
		_start(0), _end(0), _chr(""), _startSite(NULL), _endSite(NULL) {
	//xxxx

}
;

Bampairs::~Bampairs() {
	//xxx
}
;


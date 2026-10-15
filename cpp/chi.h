#ifndef _CHI_H_
#define _CHI_H_

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "param.h"

double get_chisqr_p (unsigned short *first, unsigned short *second);
double X2BetweenTwo (unsigned short *FirstOriginal,
                     unsigned short *SecondOriginal, unsigned int dispots);

#endif //_CHI_H_

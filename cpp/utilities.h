

#ifndef _UTILITIES_H_
#define _UTILITIES_H_

#include<time.h>

time_t Initial_Time();
//time used during the past step
time_t Cal_StepTime();
//total time exhaust
time_t Cal_AllTime();
//current time on string format
char * Curr_Time();

#endif //_UTILITIES_H_


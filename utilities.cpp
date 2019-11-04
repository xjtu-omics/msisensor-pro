

#include "utilities.h"

static time_t timeBegin;
static time_t timeLast;

time_t Initial_Time() {    
    timeBegin = time(NULL);
    timeLast = timeBegin;

    return timeBegin;
};

//time used during the past step
time_t Cal_StepTime() {
    time_t tused = time(NULL)-timeLast;
    timeLast = time(NULL);

    return tused;
};

//total time exhaust
time_t Cal_AllTime() { 
    return time(NULL) - timeBegin;
};

//current time on string format
char * Curr_Time() { 
    time_t t=time(NULL); 
    return ctime(&t);
}


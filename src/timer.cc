
#include "timer.h"

Timer::Timer()
{
    gettimeofday(&start_time, NULL);
}

double
Timer::elapsed() const
{
    struct timeval current_time;
    gettimeofday(&current_time, NULL);

    // the 1000's are for units (sec to ms and us to ms)
    double e = (current_time.tv_sec - start_time.tv_sec) * 1000.0;
    e += (current_time.tv_usec - start_time.tv_usec) / 1000.0;
    return e;
}

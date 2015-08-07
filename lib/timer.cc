/*
2014, the dachshund authors.
*/

#include <sys/time.h>

#include "timer.h"

Timer::Timer()
{
    // Take the first time.
    reset();
}

void
Timer::reset()
{
    t0 = get_time();
}

tv
Timer::get_time()
{
    tv t;
    gettimeofday(&t, NULL);
    return t;
}

double
Timer::interval(const tv t0, const tv t1)
{
    return (t1.tv_sec - t0.tv_sec) + 1.0e-6 * (t1.tv_usec - t0.tv_usec);
}

double
Timer::elapsed() const
{
    tv t1 = get_time();
    return interval(t0, t1);
}

Timer total_timer;

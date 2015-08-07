/*
2014, the dachshund authors.
*/

#ifndef __DS_TIMER_H__
#define __DS_TIMER_H__

#include <cstddef>
#include <sys/time.h>

typedef struct timeval tv;

class Timer {
  public:
    Timer();

    static tv get_time();
    static double interval(const tv t0, const tv t1);

    void reset();
    double elapsed() const;
  private:
    tv t0;
};

extern Timer total_timer;

#endif

/*
2014, the dachshund authors.
*/

#ifndef __DS_TIMER_H__
#define __DS_TIMER_H__

#include <cstddef>
#include <sys/time.h>

class Timer {
  public:
    Timer();

    void stop();
    void restart();
    double elapsed() const;

  private:
    struct timeval start_time;
};

#endif

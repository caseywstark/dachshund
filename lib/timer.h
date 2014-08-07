/*
2014, the dachshund authors.
*/

#ifndef __DS_TIMER_H__
#define __DS_TIMER_H__

#include <string>
#include <vector>

#include <cstddef>
#include <sys/time.h>

typedef struct timeval tv;

class Timer {
  private:
    std::vector<std::string> categories;
    std::vector<std::vector<tv> > times;
  public:
    Timer();
    int record(const int icat);
    int add_category(const std::string category_name);
    double sec_interval(const int icat, const int i0, const int i1) const;
    std::string report();
};

// the global timer
extern Timer g_timer;

#endif

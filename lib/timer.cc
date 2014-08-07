/*
2014, the dachshund authors.
*/

#include <iostream>
#include <sstream>

#include "timer.h"

// the global timer
Timer g_timer;

Timer::Timer()
{
    // Add general cat
    int i = add_category("general");
    // take first time.
    int i0 = record(i);
}

int
Timer::record(const int icat)
{
    int it = times[icat].size();
    tv t;
    gettimeofday(&t, NULL);
    times[icat].push_back(t);
    return it;
}

int
Timer::add_category(std::string cat_name)
{
    categories.push_back(cat_name);

    const size_t ncat = times.size();
    times.resize(ncat + 1);
    return ncat;
}

double
Timer::sec_interval(const int icat, const int i0, const int i1) const
{
    tv t0 = times[icat][i0];
    tv t1 = times[icat][i1];
    double dt = (t1.tv_sec - t0.tv_sec) + 1.0e-6 * (t1.tv_usec - t0.tv_usec);
    return dt;
}

std::string
Timer::report()
{
    // Add final time.
    int i_final = record(0);
    double total_time = sec_interval(0, 0, i_final);

    std::ostringstream s;
    s << "Total timing" << std::endl;
    for (int ic = 0; ic < (int)times.size(); ++ic) {
        int num_times = times[ic].size();
        double total_cat_time = sec_interval(ic, 0, num_times-1);

        s << categories[ic] << ": time " << total_cat_time << " s, fraction "
          << total_cat_time / total_time << std::endl;
        for (int i = 1; i < num_times; ++i) {
            double dt = sec_interval(ic, i-1, i);
            s << "    " << i << ": " << dt << " s, " << dt / total_cat_time
              << std::endl;
        }
    }

    std::string r = s.str();
    return r;
}

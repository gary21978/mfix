#ifndef MFIX_RANDOM_H_
#define MFIX_RANDOM_H_

#include <time.h>

namespace aux {

int calc_random_int () noexcept {

  time_t timer;
  struct tm y2k;

  y2k.tm_hour = 0;   y2k.tm_min = 0; y2k.tm_sec = 0;
  y2k.tm_year = 100; y2k.tm_mon = 0; y2k.tm_mday = 1;

  time(&timer);  /* get current time; same as: timer = time(NULL)  */
  return static_cast<int>(difftime(timer,mktime(&y2k)));
}

} // namespace

#endif

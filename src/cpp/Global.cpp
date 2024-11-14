#include "Global.hpp"



double tempo::cpu_time(void) {
  struct rusage ru;
  getrusage(RUSAGE_SELF, &ru);
  return (double)ru.ru_utime.tv_sec + (double)ru.ru_utime.tv_usec / 1000000;
}

static unsigned long x = 123456789, y = 362436069, z = 521288629;

void tempo::setSeed(const unsigned long x_, const unsigned long y_, const unsigned long z_) {
    x = x_;
    y = y_;
    z = z_;
}


void tempo::seed(const unsigned long s) {
  x = s;
  y = 362436069;
  z = 521288629;
}

unsigned long tempo::random(void) { // period 2^96-1
  unsigned long t;
  x ^= x << 16;
  x ^= x >> 5;
  x ^= x << 1;

  t = x;
  x = y;
  y = z;
  z = t ^ x ^ y;

  return z;
}

/**
* @author Tim Luchterhand
* @date 14.11.24
* @brief
*/

#include "util/random.hpp"

namespace tempo {

    static unsigned long x = 123456789, y = 362436069, z = 521288629;

    void seed(unsigned long x_, unsigned long y_, unsigned long z_) noexcept {
        x = x_;
        y = y_;
        z = z_;
    }

    unsigned long random() noexcept {
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

    auto RNG::operator()() const noexcept -> unsigned long { return random(); }
}
#pragma once

#include <stdint.h>

namespace dyft {

// From http://xoroshiro.di.unimi.it/splitmix64.c
class splitmix64 {
  public:
    splitmix64(uint64_t seed) : x(seed){};

    uint64_t next() {
        uint64_t z = (x += uint64_t(0x9E3779B97F4A7C15));
        z = (z ^ (z >> 30)) * uint64_t(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)) * uint64_t(0x94D049BB133111EB);
        return z ^ (z >> 31);
    }

  private:
    uint64_t x;
};

}  // namespace dyft

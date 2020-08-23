#pragma once

#include <stdint.h>
#include <cassert>

#ifdef __SSE4_2__
#include <xmmintrin.h>
#endif

namespace dyft {

struct bit_tools {
    static const uint8_t POPCNT_TABLE[256];

    static int popcnt(uint8_t x) {
        return POPCNT_TABLE[x];
    }
    static int popcnt(uint16_t x) {
        return POPCNT_TABLE[x & UINT8_MAX] + POPCNT_TABLE[x >> 8];
    }
    static int popcnt(uint32_t x) {
#ifdef __SSE4_2__
        return static_cast<int>(__builtin_popcount(x));
#else
        x = x - ((x >> 1) & 0x55555555);
        x = (x & 0x33333333) + ((x >> 2) & 0x33333333);
        return (0x10101010 * x >> 28) + (0x01010101 * x >> 28);
#endif
    }
    static int popcnt(uint64_t x) {
#ifdef __SSE4_2__
        return static_cast<int>(__builtin_popcountll(x));
#else
        x = x - ((x >> 1) & 0x5555555555555555ull);
        x = (x & 0x3333333333333333ull) + ((x >> 2) & 0x3333333333333333ull);
        x = (x + (x >> 4)) & 0x0f0f0f0f0f0f0f0full;
        return (0x0101010101010101ull * x >> 56);
#endif
    }

    static void set_bit(uint8_t& x, int i, bool bit = true) {
        assert(0 <= i and i <= 8);
        if (bit) {
            x |= (1U << i);
        } else {
            x &= ~(1U << i);
        }
    }
    static void set_bit(uint16_t& x, int i, bool bit = true) {
        assert(0 <= i and i <= 16);
        if (bit) {
            x |= (1U << i);
        } else {
            x &= ~(1U << i);
        }
    }
    static void set_bit(uint32_t& x, int i, bool bit = true) {
        assert(0 <= i and i <= 32);
        if (bit) {
            x |= (1U << i);
        } else {
            x &= ~(1U << i);
        }
    }
    static void set_bit(uint64_t& x, int i, bool bit = true) {
        assert(0 <= i and i <= 64);
        if (bit) {
            x |= (1ULL << i);
        } else {
            x &= ~(1ULL << i);
        }
    }

    static bool get_bit(uint8_t x, int i) {
        assert(0 <= i and i <= 8);
        return (x & (1U << i)) != 0;
    }
    static bool get_bit(uint16_t x, int i) {
        assert(0 <= i and i <= 16);
        return (x & (1U << i)) != 0;
    }
    static bool get_bit(uint32_t x, int i) {
        assert(0 <= i and i <= 32);
        return (x & (1U << i)) != 0;
    }
    static bool get_bit(uint64_t x, int i) {
        assert(0 <= i and i <= 64);
        return (x & (1ULL << i)) != 0;
    }
};

const uint8_t bit_tools::POPCNT_TABLE[256] = {
    0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 1, 2, 2, 3, 2,
    3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3,
    3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5,
    6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4,
    3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4,
    5, 5, 6, 5, 6, 6, 7, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6,
    6, 7, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8};

}  // namespace dyft

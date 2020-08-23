#pragma once

#include "bit_tools.hpp"

namespace dyft {

template <int N>
struct vcode_traits;

template <>
struct vcode_traits<8> {
    using vint_type = uint8_t;
};
template <>
struct vcode_traits<16> {
    using vint_type = uint16_t;
};
template <>
struct vcode_traits<32> {
    using vint_type = uint32_t;
};
template <>
struct vcode_traits<64> {
    using vint_type = uint64_t;
};

template <int N>
struct vcode_tools {
    using vint_type = typename vcode_traits<N>::vint_type;

    static int get_hamdist(const vint_type* x, const vint_type* y, int bits) {
        if (bits == 1) {
            return bit_tools::popcnt(x[0] ^ y[0]);
        } else {
            vint_type diff = 0;
            for (int j = 0; j < bits; ++j) {
                diff |= (x[j] ^ y[j]);
            }
            return bit_tools::popcnt(diff);
        }
    }

    static int get_hamdist(const vint_type* x, const vint_type* y, int bits, int radius) {
        if (bits == 1) {
            return bit_tools::popcnt(x[0] ^ y[0]);
        } else {
            int dist = 0;
            vint_type diff = 0;
            for (int j = 0; j < bits; ++j) {
                diff |= (x[j] ^ y[j]);
                dist = bit_tools::popcnt(diff);
                if (dist > radius) {
                    return dist;
                }
            }
            return dist;
        }
    }

    static vint_type to_vint(const uint8_t* in, int j) {
        vint_type v = vint_type(0);
        for (int i = 0; i < N; ++i) {
            vint_type b = (in[i] >> j) & vint_type(1);
            v |= (b << i);
        }
        return v;
    }

    // static const vint_type* to_vints(const uint8_t* in, int bits) {
    //     static vint_type out[64];
    //     for (int j = 0; j < bits; ++j) {
    //         out[j] = to_vint(in, j);
    //     }
    //     return out;
    // }

    static void to_vints(const uint8_t* in, vint_type* out, int bits) {
        for (int j = 0; j < bits; ++j) {
            out[j] = to_vint(in, j);
        }
    }

    static uint8_t get_int(const vint_type* x, int i, int bits) {
        DEBUG_ABORT_IF_LE(N, i);

        if (bits == 1) {
            return (x[0] >> i) & uint8_t(1);
        } else {
            uint8_t c = uint8_t(0);
            for (int j = 0; j < bits; ++j) {
                c |= (((x[j] >> i) & uint8_t(1)) << j);
            }
            return c;
        }
    }
};

}  // namespace dyft

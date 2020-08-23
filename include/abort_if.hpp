#pragma once

#include <iostream>
#include <string>

#define ABORT_IF_(x, y, z)                                                                                           \
    if ((x)) {                                                                                                       \
        std::cerr << "\033[0;31mERROR in " << __FILE__ << " at line " << std::to_string(__LINE__) << ", if " << (#x) \
                  << " for " << std::to_string(y) << " vs " << std::to_string(z) << "\033[0;0m" << std::endl;        \
        abort();                                                                                                     \
    }

#define ABORT_IF(x)                                                                                                  \
    if ((x)) {                                                                                                       \
        std::cerr << "\033[0;31mERROR in " << __FILE__ << " at line " << std::to_string(__LINE__) << ", if " << (#x) \
                  << "\033[0;0m" << std::endl;                                                                       \
        abort();                                                                                                     \
    }

#define ABORT_IF_EQ(x, y) ABORT_IF_((x) == (y), x, y)
#define ABORT_IF_NE(x, y) ABORT_IF_((x) != (y), x, y)
#define ABORT_IF_LE(x, y) ABORT_IF_((x) <= (y), x, y)
#define ABORT_IF_LT(x, y) ABORT_IF_((x) < (y), x, y)
#define ABORT_IF_GE(x, y) ABORT_IF_((x) >= (y), x, y)
#define ABORT_IF_GT(x, y) ABORT_IF_((x) > (y), x, y)
#define ABORT_IF_OUT(n, x, y) ABORT_IF_((n) < (x), n, x) ABORT_IF_((y) < (n), y, n)

#ifdef NDEBUG
#define DEBUG_ABORT_IF(x)
#define DEBUG_ABORT_IF_EQ(x, y)
#define DEBUG_ABORT_IF_NE(x, y)
#define DEBUG_ABORT_IF_LE(x, y)
#define DEBUG_ABORT_IF_LT(x, y)
#define DEBUG_ABORT_IF_GE(x, y)
#define DEBUG_ABORT_IF_GT(x, y)
#define DEBUG_ABORT_IF_OUT(n, x, y)
#else
#define DEBUG_ABORT_IF(x) ABORT_IF((x))
#define DEBUG_ABORT_IF_EQ(x, y) ABORT_IF_((x) == (y), x, y)
#define DEBUG_ABORT_IF_NE(x, y) ABORT_IF_((x) != (y), x, y)
#define DEBUG_ABORT_IF_LE(x, y) ABORT_IF_((x) <= (y), x, y)
#define DEBUG_ABORT_IF_LT(x, y) ABORT_IF_((x) < (y), x, y)
#define DEBUG_ABORT_IF_GE(x, y) ABORT_IF_((x) >= (y), x, y)
#define DEBUG_ABORT_IF_GT(x, y) ABORT_IF_((x) > (y), x, y)
#define DEBUG_ABORT_IF_OUT(n, x, y) ABORT_IF_((n) < (x), n, x) ABORT_IF_((y) < (n), y, n)
#endif

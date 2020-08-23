#pragma once

#include <boost/timer/timer.hpp>

#include "tinyformat/tinyformat.h"

#define START_TIMER(___name) boost::timer::cpu_timer ___name

#define STOP_TIMER(___name) ___name.stop()

#define STOP_TIMER_V(___name)                                         \
    do {                                                              \
        ___name.stop();                                               \
        tfm::reportf("[%s] %s", #___name, ___name.format(4).c_str()); \
    } while (false);

#define GET_TIMER_SEC(___name) ___name.elapsed().wall / 1000000000.0

#define GET_TIMER_MILLISEC(___name) ___name.elapsed().wall / 1000000.0

#pragma once

#ifdef __APPLE__
#include <mach/mach.h>
#endif

#include <cxxabi.h>
#include <stdint.h>
#include <array>
#include <cassert>
#include <chrono>
#include <exception>
#include <fstream>
#include <iostream>
#include <numeric>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

#include "bit_tools.hpp"

namespace dyft {

template <class T>
inline float get_average(const std::vector<T>& vec) {
    return std::accumulate(vec.begin(), vec.end(), 0.0) / vec.size();
}

inline std::string normalize_filepath(std::string filepath) {
    std::replace(filepath.begin(), filepath.end(), '/', '_');
    std::replace(filepath.begin(), filepath.end(), '.', '_');
    std::replace(filepath.begin(), filepath.end(), ':', '_');
    return filepath;
}

inline int get_hamdist(uint64_t x, uint64_t y) {
    return static_cast<int>(bit_tools::popcnt(x ^ y));
}

// ceil(a / b), cf. https://nariagari-igakusei.com/cpp-division-round-up/
inline constexpr int ceil_div(int a, int b) {
    return (a + (b - 1)) / b;
}

// From Cedar (http://www.tkl.iis.u-tokyo.ac.jp/~ynaga/cedar/)
inline size_t get_process_size_in_bytes() {
#ifdef __APPLE__
    struct task_basic_info t_info;
    mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;
    task_info(current_task(), TASK_BASIC_INFO, reinterpret_cast<task_info_t>(&t_info), &t_info_count);
    return t_info.resident_size;
#else
    FILE* fp = std::fopen("/proc/self/statm", "r");
    size_t dummy(0), vm(0);
    std::fscanf(fp, "%ld %ld ", &dummy, &vm);  // get resident (see procfs)
    std::fclose(fp);
    return vm * ::getpagesize();
#endif
}

inline double to_MiB(size_t bytes) {
    return bytes / (1024.0 * 1024.0);
}

inline double to_GiB(size_t bytes) {
    return bytes / (1024.0 * 1024.0 * 1024.0);
}

}  // namespace dyft

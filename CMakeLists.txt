cmake_minimum_required(VERSION 3.0)
project(dyft)

if(NOT CMAKE_BUILD_TYPE)
  # set(CMAKE_BUILD_TYPE "Debug")
  set(CMAKE_BUILD_TYPE "Release")
endif(NOT CMAKE_BUILD_TYPE)

# C++17 compiler check
if (CMAKE_CXX_COMPILER MATCHES ".*clang.*" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  set(CMAKE_COMPILER_IS_CLANGXX 1)
endif ()
if (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  set(CMAKE_COMPILER_IS_INTEL 1)
endif ()
if ((CMAKE_COMPILER_IS_GNUCXX AND ${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS 7.0) OR (CMAKE_COMPILER_IS_CLANGXX AND ${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS 4.0))
  message(FATAL_ERROR "Your C++ compiler does not support C++17. Please install g++ 7.0 (or greater) or clang 4.0 (or greater)")
else ()
  message(STATUS "Compiler is recent enough to support C++17.")
endif ()

# Boost
set(Boost_USE_MULTITHREADED ON) # https://gitlab.kitware.com/cmake/cmake/issues/19714
set(Boost_USE_STATIC_LIBS ON)
find_package(Boost REQUIRED COMPONENTS system filesystem iostreams timer date_time)
include_directories(${Boost_INCLUDE_DIRS})

# SSE Popcnt
set(BUILTIN_POPCNT 0)

if (DISABLE_SSE4_2)
  message(STATUS "sse4.2 disabled")
elseif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
  # Use /proc/cpuinfo to get the information
  file(STRINGS "/proc/cpuinfo" _cpuinfo)
  if(_cpuinfo MATCHES "(sse4_2)|(sse4a)")
    set(BUILTIN_POPCNT 1)
  endif()
elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
  execute_process(COMMAND sysctl -n machdep.cpu.features OUTPUT_VARIABLE _cpuinfo OUTPUT_STRIP_TRAILING_WHITESPACE)
  # message(STATUS "_cpuinfo is ${_cpuinfo}")
  if(_cpuinfo MATCHES "SSE4.2")
    set(BUILTIN_POPCNT 1)
  endif()
endif()

if(BUILTIN_POPCNT)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4.2")
endif()

# AVX2
set(SIMD_COMPARISON 0)

if (DISABLE_AVX)
  message(STATUS "AVX disabled")
elseif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
  # # Use /proc/cpuinfo to get the information
  file(STRINGS "/proc/cpuinfo" _cpuinfo)
  if(_cpuinfo MATCHES "avx2")
    set(SIMD_COMPARISON 2)
  elseif(_cpuinfo MATCHES "avx")
    set(SIMD_COMPARISON 1)
  endif()
elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
  execute_process(COMMAND sysctl -n machdep.cpu.leaf7_features OUTPUT_VARIABLE _cpuinfo OUTPUT_STRIP_TRAILING_WHITESPACE)
  # message(STATUS "_cpuinfo is ${_cpuinfo}")
  if(_cpuinfo MATCHES "AVX2")
    set(SIMD_COMPARISON 2)
  else()
    execute_process(COMMAND sysctl -n machdep.cpu.features OUTPUT_VARIABLE _cpuinfo OUTPUT_STRIP_TRAILING_WHITESPACE)
    # message(STATUS "_cpuinfo is ${_cpuinfo}")
    if(_cpuinfo MATCHES "AVX1.0")
      set(SIMD_COMPARISON 1)
    endif()
  endif()
endif()

if(SIMD_COMPARISON EQUAL 1)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx")
elseif(SIMD_COMPARISON EQUAL 2)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2")
endif()


# -fvisibility=hidden is enabled to ignore "ld: warning: direct access in function for boost::filesystem."
# https://stackoverflow.com/questions/36567072/why-do-i-get-ld-warning-direct-access-in-main-to-global-weak-symbol-in-this
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1z -pthread -fvisibility=hidden -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DNDEBUG -march=native -O3")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address -fno-omit-frame-pointer -O0 -g -DDEBUG")

message(STATUS "BUILD_TYPE is ${CMAKE_BUILD_TYPE}")
message(STATUS "CXX_FLAGS are ${CMAKE_CXX_FLAGS}")
message(STATUS "CXX_FLAGS_DEBUG are ${CMAKE_CXX_FLAGS_DEBUG}")
message(STATUS "CXX_FLAGS_RELEASE are ${CMAKE_CXX_FLAGS_RELEASE}")

include_directories(include)
add_subdirectory(src)

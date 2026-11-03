cmake_minimum_required(VERSION 3.15)
project(pydicomrt LANGUAGES NONE)

# The compiled core of rigid_registration's default algorithm: a plain shared library with a C ABI,
# loaded through ctypes by pydicomrt/reg/method/_rigid_core.py. No Python headers, so one wheel per
# platform serves every Python version.
#
# The C++ compiler is optional. Without one -- an install from the sdist on a machine with no
# toolchain -- the package installs without the library and the algorithm runs on its numpy
# implementation, about ten times slower. Wheel builds set PYDICOMRT_REQUIRE_CORE=ON so a missing
# compiler fails there instead of shipping a slow wheel.
option(PYDICOMRT_REQUIRE_CORE "Fail when no C++ compiler is available" OFF)

include(CheckLanguage)
check_language(CXX)
if(NOT CMAKE_CXX_COMPILER)
  if(PYDICOMRT_REQUIRE_CORE)
    message(FATAL_ERROR "pydicomrt: no C++ compiler found and PYDICOMRT_REQUIRE_CORE is ON")
  endif()
  message(WARNING "pydicomrt: no C++ compiler found; installing without the compiled rigid core "
                  "(rigid_registration falls back to numpy, about ten times slower)")
  return()
endif()
enable_language(CXX)

add_library(rigid_core SHARED csrc/rigid_core.cpp)
target_compile_features(rigid_core PRIVATE cxx_std_17)
set_target_properties(rigid_core PROPERTIES CXX_VISIBILITY_PRESET hidden)
if(MSVC)
  target_compile_options(rigid_core PRIVATE /O2)
else()
  target_compile_options(rigid_core PRIVATE -O3)
endif()

# std::thread rather than OpenMP: blocking rather than spinning workers, the same threads on every
# platform (Apple clang ships no OpenMP), and no runtime library to bundle into the wheels.
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
target_link_libraries(rigid_core PRIVATE Threads::Threads)

install(TARGETS rigid_core
        LIBRARY DESTINATION pydicomrt/reg/method
        RUNTIME DESTINATION pydicomrt/reg/method)

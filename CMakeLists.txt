cmake_minimum_required (VERSION 3.5)



project(Coupling)


if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt. ")
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED true)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Locate  FFTW3
find_package(FFTW)

# add in Eigen3
find_package (Eigen3 3.3 REQUIRED NO_MODULE)

# add in Eigen3
include(FetchContent)
FetchContent_Declare(
  Eigen3
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG master
)
FetchContent_MakeAvailable(Eigen3)

# find_package (FFTWpp)
include(FetchContent)
FetchContent_Declare(
  FFTWpp
  GIT_REPOSITORY https://github.com/da380/FFTWpp
  GIT_TAG main
)
FetchContent_MakeAvailable(FFTWpp)



add_subdirectory(src)

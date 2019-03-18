cmake_minimum_required (VERSION 2.6)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

project (quasistatic_pic)

find_package(HDF5 REQUIRED COMPONENTS C CXX)
include_directories(${HDF5_INCLUDE_DIRS})

add_executable(quasistatic_pic_2d quasistatic_pic_2d.cpp)
target_link_libraries(quasistatic_pic_2d ${HDF5_LIBRARIES} fftw3)

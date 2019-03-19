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
add_executable(quasistatic_pic_3d quasistatic_pic_3d.cpp)
add_executable(test_2d_fourier test_2d_fourier.cpp)
add_executable(test_1d_fourier test_1d_fourier.cpp)
target_link_libraries(quasistatic_pic_2d ${HDF5_LIBRARIES} fftw3)
target_link_libraries(quasistatic_pic_3d ${HDF5_LIBRARIES} fftw3)
target_link_libraries(test_2d_fourier ${HDF5_LIBRARIES} fftw3)
target_link_libraries(test_1d_fourier ${HDF5_LIBRARIES} fftw3)
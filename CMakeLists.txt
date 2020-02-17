cmake_minimum_required (VERSION 3.10)

project (quasistatic_pic)

find_package(HDF5 REQUIRED COMPONENTS C CXX)
find_package(fmt REQUIRED)

add_executable(quasistatic_pic_2d src/quasistatic_pic_2d.cpp)
add_executable(quasistatic_pic_3d src/quasistatic_pic_3d.cpp src/System_3d.cpp)
target_include_directories(quasistatic_pic_2d PRIVATE include ${HDF5_INCLUDE_DIRS})
target_include_directories(quasistatic_pic_3d PRIVATE include ${HDF5_INCLUDE_DIRS})

target_link_libraries(quasistatic_pic_2d ${HDF5_LIBRARIES} fftw3)
target_link_libraries(quasistatic_pic_3d ${HDF5_LIBRARIES} fftw3 fmt::fmt)
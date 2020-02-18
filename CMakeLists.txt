cmake_minimum_required (VERSION 3.10)

project (quasistatic_pic)

find_package(HDF5 REQUIRED COMPONENTS C CXX)
find_package(fmt REQUIRED)
find_package(cpptoml CONFIG REQUIRED)

if(cpptoml_FOUND)
    message(STATUS "cpptoml config found in ${cpptoml_DIR}")
else()
    message(FATAL_ERROR "cpptoml is not found; provide the path to cpptomlConfig.cmake in CMAKE_PREFIX_PATH")
endif()

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -DNDEBUG")

add_executable(quasistatic_pic_2d src/quasistatic_pic_2d.cpp)
add_executable(quasistatic_pic_3d src/quasistatic_pic_3d.cpp src/System_3d.cpp src/Fourier2d.cpp src/output.cpp src/Config_reader.cpp src/profiles.cpp)
target_include_directories(quasistatic_pic_2d PRIVATE include ${HDF5_INCLUDE_DIRS})
target_include_directories(quasistatic_pic_3d PRIVATE include ${HDF5_INCLUDE_DIRS})

target_link_libraries(quasistatic_pic_2d ${HDF5_LIBRARIES} fftw3)
target_link_libraries(quasistatic_pic_3d ${HDF5_LIBRARIES} fftw3 fmt::fmt cpptoml)

target_compile_options(quasistatic_pic_2d PRIVATE -Wall -Wextra)
target_compile_options(quasistatic_pic_3d PRIVATE -Wall -Wextra)
cmake_minimum_required (VERSION 3.10)

if(POLICY CMP0074)
    cmake_policy(SET CMP0074 NEW)
endif()


project (quasistatic_pic)

find_package(HDF5 REQUIRED COMPONENTS C CXX)
find_package(fmt REQUIRED)
find_package(cpptoml CONFIG REQUIRED)

if(cpptoml_FOUND)
    message(STATUS "cpptoml config found in ${cpptoml_DIR}")
else()
    message(FATAL_ERROR "cpptoml is not found; provide the path to cpptomlConfig.cmake in CMAKE_PREFIX_PATH")
endif()

find_package(OpenMP REQUIRED)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -DNDEBUG")

add_executable(quasistatic_pic src/quasistatic_pic.cpp src/System_3d.cpp src/Fourier2d.cpp src/output.cpp 
               src/Config_reader.cpp src/profiles.cpp src/array_utils.cpp src/array1d.cpp src/array2d.cpp)
target_include_directories(quasistatic_pic PRIVATE ${HDF5_INCLUDE_DIRS})
target_include_directories(quasistatic_pic PRIVATE ${HDF5_INCLUDE_DIRS})

target_link_libraries(quasistatic_pic ${HDF5_LIBRARIES} fftw3 fftw3_omp fmt::fmt cpptoml OpenMP::OpenMP_CXX)

target_compile_options(quasistatic_pic PRIVATE -Wall -Wextra)

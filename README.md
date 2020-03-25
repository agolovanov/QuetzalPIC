QuetzalPIC is a 2D and 3D quasistatic particle-in-cell (PIC) code.

# Build

QuetzalPIC is written in C++14 and uses CMake for building the project. It has been tested to compile and run under Linux.

To build the project, perform
```
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

## Dependencies

The following dependencies are required for QuetzalPIC:

* [HDF5](https://www.hdfgroup.org/solutions/hdf5/) with C++ support.
* [fmt](https://github.com/fmtlib/fmt).
* [cpptoml](https://github.com/skystrife/cpptoml).
* [fftw3](http://www.fftw.org/) with OpenMP support.

# Run

To run the code, run the `quasitatic_pic_3d` executable with the input file as the first argument, e.g.
```
./quasistatic_pic_3d ../inputs/test/laser_test.toml
```
Input files are written in the TOML language. Examples can be found in the inputs folder.
The `OMP_NUM_THREADS` environment variable might be used to set the number of threads.
Output data is written to the current folder in the shell.

Alternatively, the `run.sh` shell script can be used:
```
./run.sh %path_to_toml% %number_of_threads%
```
In this case, output is written the the directory with the same path and name as the `.toml` file.

# Output

All output is written in the HDF5 format.
Any tool capable of reading HDF5 files can be used to analyze the data.

A Python helper located in `utils/reader.py` provides some useful methods.

# Acknowledgements

The development was supported by the Russian Foundation for Basic Research (Grant 18-32-00943).

#pragma once

#include <H5Cpp.h>
#include "array1d.h"
#include "array2d.h"
#include "array3d.h"
#include "System_parameters.h"

class Output_writer {
public:
    Output_writer(Output_parameters output_parameters);

    void initialize_slice_array(ivector3d size, vector3d steps, const array2d & array, const std::string name);
    void write_array(array3d & array, std::string name);
    void write_slice(array2d & array, std::string name, int slice_index);
private:
    Output_parameters output_parameters;
    H5::H5File fields_file;
    H5::H5File fields_xy_file;
};


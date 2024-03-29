#pragma once

#include <H5Cpp.h>
#include "array1d.h"
#include "array2d.h"
#include "array3d.h"
#include "System_parameters.h"
#include "containers.h"
#include <vector>

class Output_writer {
public:
    Output_writer(Output_parameters output_parameters, int count);

    void initialize_slice_array(ivector3d size, vector3d steps, const array2d & array, const std::string name);
    void write_array(array3d & array, std::string name);
    void write_slice(array2d & array, std::string name, int slice_index);
    void write_bunch_parameters(const std::vector<bunch_particle_3d> & particles, const double weight_norm);
    void write_photon_parameters(const std::vector<bunch_particle_3d> & photons, const double weight_norm);
private:
    Output_parameters output_parameters;
    H5::H5File fields_file;
    H5::H5File fields_xy_file;
    H5::H5File bunch_parameters_file;
    H5::H5File photons_file;
};


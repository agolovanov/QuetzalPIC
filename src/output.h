#pragma once

#include <H5Cpp.h>
#include <string>
#include "containers_3d.h"

void write_array(const array3d & array, const std::string name, H5::H5File file);
void write_array(const array2d & array, const std::string name, H5::H5File file);
void initialize_slice_array(ivector3d size, vector3d steps, vector3d origin, const std::string name, H5::H5File file);
void initialize_slice_array(ivector2d size, vector2d steps, vector2d origin, Plane plane, double plane_coordinate,
                            const std::string name, H5::H5File file);
void write_slice(const array2d & slice, int slice_index, const std::string name, H5::H5File file);
void write_slice(const array1d & slice, int slice_index, const std::string name, H5::H5File file);
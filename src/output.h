#pragma once

#include <H5Cpp.h>
#include <string>
#include "containers_3d.h"

void write_array(const array3d & array, const std::string name, H5::H5File file);
void write_array(const array2d & array, const std::string name, H5::H5File file);
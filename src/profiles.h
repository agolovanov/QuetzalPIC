#pragma once

#include <functional>
#include "containers_3d.h"

std::function<double(double, double, double)> gaussian(double amplitude, vector3d width, vector3d r0);
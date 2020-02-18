#pragma once

#include <functional>
#include "containers_3d.h"

std::function<double(double, double, double)> gaussian(double amplitude, dvector3d width, dvector3d r0);
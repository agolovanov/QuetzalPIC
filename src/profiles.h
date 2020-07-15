#pragma once

#include <functional>
#include "containers.h"

std::function<double(double, double, double)> gaussian3d(double amplitude, vector3d width, vector3d r0);
std::function<double(double, double, double)> constant3d(double amplitude = 1.0);
std::function<double(double, double)> constant2d(double amplitude = 1.0);
std::function<double(double, double)> powerlaw2d(double power, double radius, double x0, double y0, double min_factor = 0.0, double max_factor = 1e10, double amplitude = 1.0);
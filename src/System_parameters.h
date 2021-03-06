#pragma once

#include <functional>
#include "containers.h"
#include "profiles.h"

struct Output_parameters {
    bool output3d = false;
    bool output_xy = false;
    double z0;
};

struct System_parameters {
    vector3d l;
    vector3d d;
    int ppcy = 1;
    int ppcz = 1;
    int magnetic_field_iterations = 1;

    std::function<double(double, double, double)> a_sqr = constant3d(0.0);
    std::function<double(double, double, double)> rho = constant3d(0.0);
    std::function<double(double, double)> plasma_profile = constant2d(-1.0);

    Output_parameters output_parameters;
};
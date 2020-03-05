#pragma once

#include <functional>
#include "containers_3d.h"
#include "profiles.h"

struct Output_parameters {
    bool output3d = false;
    bool output_xy = false;
};

struct System_parameters {
    vector3d l;
    vector3d d;
    int ppcy = 1;
    int ppcz = 1;
    int magnetic_field_iterations = 1;

    std::function<double(double, double, double)> a_sqr = constant(0.0);
    std::function<double(double, double, double)> rho = constant(0.0);

    Output_parameters output_parameters;
};
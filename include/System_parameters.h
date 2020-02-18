#pragma once

#include <functional>
#include "containers_3d.h"

struct System_parameters {
    dvector3d l;
    dvector3d d;
    int ppcy = 1;
    int ppcz = 1;
    int magnetic_field_iterations = 1;

    std::function<double(double, double, double)> a_sqr = [] (double x, double y, double z) -> double {return 0;};
    std::function<double(double, double, double)> rho = [] (double x, double y, double z) -> double {return 0;};
};
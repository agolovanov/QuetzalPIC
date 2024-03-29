#pragma once

#include <functional>
#include "containers.h"
#include "profiles.h"
#include "Species_register.h"

struct Output_parameters {
    bool output3d = false;
    bool output_xy = false;
    bool output_bunch = false;
    bool output_photons = true;
    double z0;
};

struct Bunch_parameters {
    std::function<double(double, double, double)> rho = constant3d(0.0);
    ivector3d ppc = {1, 1, 1};
    double gamma = 1.0;
    int species_id;
};

struct System_parameters {
    vector3d l;
    vector3d d;
    double dt;
    double t_end;
    int ppcy = 1;
    int ppcz = 1;
    int magnetic_field_iterations = 1;

    std::function<double(double, double, double)> a_sqr = constant3d(0.0);
    std::function<double(double, double)> plasma_profile = constant2d(-1.0);

    Output_parameters output_parameters;
    std::vector<Bunch_parameters> bunch_parameters_array;

    Species_register species;

    double base_frequency_SI = 1;
    bool qed = false;
};
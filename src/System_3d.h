#pragma once

#include <H5Cpp.h>
#include <string>
#include <functional>
#include <iostream>
#include "System_parameters.h"
#include "containers.h"
#include "array2d.h"
#include "array3d.h"
#include "Fourier2d.h"

class System_3d {
public:
    System_3d(System_parameters & params, std::ostream & out);
    void solve_wakefield();
private:
    void normalize_coordinates(double & y, double & z);
    void solve_poisson_equation(double D=0.0);
    void init_particles(int ppcy, int ppcz, std::function<double(double, double)> plasma_profile);
    void init_a_sqr(std::function<double(double, double, double)> func);
    void increase_minimum(array3d & array, int slice, double value) const;
    void increase_minimum(array2d & array, double value) const;
    void initialize_output(H5::H5File & fields_file, H5::H5File & fields_xy_file) const;
    void output_step(H5::H5File & fields_file, H5::H5File & fields_xy_file, int slice_index);
    void output_full(H5::H5File & fields_file, H5::H5File & fields_xy_file) const;

    ivector3d n;
    vector3d l;
    vector3d d;

    int magnetic_field_iterations;
    double psi_threshold = 1e-4;

    double magnetic_field_D = 10.0;

    std::function<double(double, double, double)> rhobunch;

    Output_parameters output_parameters;
    
    std::ostream & out;

    std::vector<particle> particles;
    array2d psi_middle;
    array2d djy_dxi;
    array2d djz_dxi;
    array2d rho_ion;

    array3d a_sqr;
    array3d susceptibility;

    array2d psi_prev;
    array2d psi;
    array2d jx;
    array2d jy;
    array2d jz;
    array2d jy_next;
    array2d jz_next;
    array2d rho;
    array2d ex;
    array2d ey;
    array2d ez;
    array2d by;
    array2d bz;

    Fourier2d fourier;
};


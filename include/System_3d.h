#pragma once

#include <H5Cpp.h>
#include <string>
#include <functional>
#include "System_parameters.h"
#include "containers_3d.h"
#include "Fourier2d.h"

class System_3d {
public:
    System_3d(System_parameters & params);
    void solve_wakefield();
    void output() const;
private:
    void deposit(double y, double z, double value, array3d & array, int slice, double yshift=0.0, double zshift=0.0);
    void deposit(double y, double z, double value, double * array);
    double array_to_particle(double y, double z, const array3d & array, int slice, double yshift=0.0, double zshift=0.0) const;
    double array_to_particle(double y, double z, const array2d & array) const;
    double array_yder_to_particle(double y, double z, const array3d & array, int slice) const;
    double array_zder_to_particle(double y, double z, const array3d & array, int slice) const;
    void normalize_coordinates(double & y, double & z);
    void solve_poisson_equation(double D=0.0);
    void init_particles(int ppcy, int ppcz);
    void init_a_sqr(std::function<double(double, double, double)> func);

    ivector3d n;
    dvector3d l;
    dvector3d d;

    int magnetic_field_iterations;

    std::vector<particle> particles;
    array2d psi_middle;
    array2d djy_dxi;
    array2d djz_dxi;

    array3d psi;
    array3d psi_source;
    array3d a_sqr;
    array3d jx;
    array3d jy;
    array3d jz;
    array3d rho;
    array3d ex;
    array3d ey;
    array3d ez;
    array3d by;
    array3d bz;

    Fourier2d fourier;

    std::function<double(double, double, double)> rhobunch;
};


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
#include "Output_writer.h"
#include "Species_register.h"

template <class T>
struct Output_reference {
    Output_reference(std::string name, T * ptr) : name(name), ptr(ptr) {}

    std::string name;
    T* ptr;
};

class System_3d {
public:
    System_3d(System_parameters & params, std::ostream & out);
    void run();
private:
    void solve_wakefield(int iteration);
    void normalize_coordinates(double & y, double & z);
    void solve_poisson_equation(double D=0.0);
    size_t count_wake_particles(int ppcy, int ppcz, std::function<double(double, double)> plasma_profile) const;
    void init_wake_particles(int ppcy, int ppcz, std::function<double(double, double)> plasma_profile);
    size_t count_bunch_particles(ivector3d ppc, std::function<double(double, double, double)> rho) const;
    void init_bunch_particles(size_t index, Bunch_parameters bunch);
    void init_a_sqr(std::function<double(double, double, double)> func);
    void increase_minimum(array3d & array, int slice, double value) const;
    void increase_minimum(array2d & array, double value) const;
    void output_step(Output_writer & output_writer, const std::vector<Output_reference<array2d>> & output_arrays_2d,
                     int slice_index);
    void calculate_em_density();

    ivector3d n;
    vector3d l;
    vector3d d;
    double dt;
    double t_end;
    int time_iterations;

    double ppcy;
    double ppcz;
    std::function<double(double, double)> plasma_profile;

    int magnetic_field_iterations;
    double psi_threshold = 1e-4;

    double magnetic_field_D = 10.0;

    double base_frequency_SI = 1.0;

    Output_parameters output_parameters;
    
    std::ostream & out;

    std::vector<wake_particle_2d> wake_particles;
    std::vector<bunch_particle_3d> bunch_particles;

    array2d psi_middle;
    array2d djy_dxi;
    array2d djz_dxi;
    array2d rho_ion;

    array3d a_sqr;
    array3d susceptibility;
    array3d rho_bunch;
    array3d jx_bunch;
    array3d ex;
    array3d ey;
    array3d ez;
    array3d by;
    array3d bz;
    array3d em_energy_density;
    array3d em_sx;

    array2d psi_prev;
    array2d psi;
    array2d jx;
    array2d jy;
    array2d jz;
    array2d jy_next;
    array2d jz_next;
    array2d rho;

    // energy density and Poynting vectors
    array2d particle_energy_density;
    array2d particle_sx;
    
    Fourier2d fourier;

    Species_register species;

    double field_schwinger;
};


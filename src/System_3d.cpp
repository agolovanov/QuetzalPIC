#include <cmath>
#include <iostream>
#include <stdexcept>
#include <H5Cpp.h>
#include <cassert>
#include <fmt/format.h>
#include <vector>
#include <string>

#include "System_3d.h"
#include "containers.h"
#include "array2d.h"
#include "array3d.h"
#include "array_utils.h"
#include "Output_writer.h"

std::string memory_formatter(long bytes) {
    double res = bytes;
    const std::vector<std::string> prefixes {"b", "KiB", "MiB", "GiB", "TiB", "PiB"};
    int index = 0;
    while (res > 1024) {
        res /= 1024;
        index++;
    }
    return fmt::format("{:5.1f} {}", res, prefixes[index]);
}

const std::string ASQR = "aSqr";
const std::string PSI = "psi";
const std::string EX = "ex";
const std::string EY = "ey";
const std::string EZ = "ez";
const std::string BY = "by";
const std::string BZ = "bz";
const std::string RHO = "rho";
const std::string JX = "jx";
const std::string JY = "jy";
const std::string JZ = "jz";
const std::string SUSCEPTIBILITY = "susceptibility";

System_3d::System_3d(System_parameters & params, std::ostream & out) : 
    l(params.l),
    magnetic_field_iterations(params.magnetic_field_iterations),
    rhobunch(params.rho),
    output_parameters(params.output_parameters),
    out(out)
{
    if (magnetic_field_iterations < 0) {
        throw std::invalid_argument("Magnetic field iterations should be non-negative");
    }
    if (params.ppcy <= 0) {
        throw std::invalid_argument("ppcy should be positive");
    }
    if (params.ppcz <= 0) {
        throw std::invalid_argument("ppcz should be positive");
    }

    n.x = static_cast<int>(l.x / params.d.x);
    n.y = static_cast<int>(l.y / params.d.y);
    n.z = static_cast<int>(l.z / params.d.z);

    out << fmt::format("Spatial dimensions:  [{}, {}, {}]", n.x, n.y, n.z) << std::endl;

    d.x = l.x / n.x;
    d.y = l.y / n.y;
    d.z = l.z / n.z;

    out << fmt::format("Steps:               [{}, {}, {}]", d.x, d.y, d.z) << std::endl;
    out << fmt::format("Simulation box size: [{}, {}, {}]", l.x, l.y, l.z) << std::endl;

    const long fourier_memory = sizeof(double) * n.y * n.z + 2 * sizeof(double) * n.y * (n.z / 2 + 1);
    const long array2d_memory = 17l * sizeof(double) * n.y * n.z;
    const long array3d_memory = 2l * sizeof(double) * n.x * n.y * n.z;
    const long particle_memory = static_cast<long>(sizeof(particle)) * params.ppcy * params.ppcz * n.y * n.z;
    const long total_memory = array2d_memory + array3d_memory + particle_memory + fourier_memory;

    out << "Expected RAM usage:\n";
    out << fmt::format("3D arrays: {}\n", memory_formatter(array3d_memory))
        << fmt::format("2D arrays: {}\n", memory_formatter(array2d_memory))
        << fmt::format("Fourier:   {}\n", memory_formatter(fourier_memory))
        << fmt::format("Particles: {}\n", memory_formatter(particle_memory))
        << fmt::format("Total:     {}", memory_formatter(total_memory)) << std::endl;

    std::cout << "----------------------------------------" << std::endl;
    
    init_particles(params.ppcy, params.ppcz, params.plasma_profile);

    fourier = Fourier2d(n.y, n.z);

    const ivector2d size_yz{n.y, n.z};
    const vector2d d_yz{d.y, d.z};
    
    psi_middle = array2d(size_yz, d_yz);
    djy_dxi = array2d(size_yz, d_yz, {0.5 * d.y, 0});
    djz_dxi = array2d(size_yz, d_yz, {0, 0.5 * d.z});
    rho_ion = array2d(size_yz, d_yz);

    a_sqr = array3d(n, d);
    susceptibility = array3d(n, d);

    psi = array2d(size_yz, d, {0, 0, 0}, Plane::YZ);
    psi_prev = array2d(size_yz, d, {0, 0, 0}, Plane::YZ);
    rho = array2d(size_yz, d, {0, 0, 0}, Plane::YZ);
    jx = array2d(size_yz, d, {0, 0, 0}, Plane::YZ);
    jy = array2d(size_yz, d, {-0.5 * d.x, 0.5 * d.y, 0}, Plane::YZ);
    jz = array2d(size_yz, d, {-0.5 * d.x, 0, 0.5 * d.z}, Plane::YZ);
    jy_next = array2d(size_yz, d, {-0.5 * d.x, 0.5 * d.y, 0}, Plane::YZ);
    jz_next = array2d(size_yz, d, {-0.5 * d.x, 0, 0.5 * d.z}, Plane::YZ);
    by = array2d(size_yz, d, {0, 0, 0.5 * d.z}, Plane::YZ);
    bz = array2d(size_yz, d, {0, 0.5 * d.y, 0}, Plane::YZ);
    ex = array2d(size_yz, d, {-0.5 * d.x, 0, 0}, Plane::YZ);
    ey = array2d(size_yz, d, {0, 0.5 * d.y, 0}, Plane::YZ);
    ez = array2d(size_yz, d, {0, 0, 0.5 * d.z}, Plane::YZ);

    init_a_sqr(params.a_sqr);
}

void System_3d::solve_wakefield() {
    auto output_writer = Output_writer(output_parameters);

    std::vector<Output_reference<array3d>> output_arrays_3d;
    output_arrays_3d.push_back(Output_reference<array3d>(ASQR, &a_sqr));
    output_arrays_3d.push_back(Output_reference<array3d>(SUSCEPTIBILITY, &susceptibility));

    std::vector<Output_reference<array2d>> output_arrays_2d;
    output_arrays_2d.push_back(Output_reference<array2d>(PSI, &psi));
    output_arrays_2d.push_back(Output_reference<array2d>(RHO, &rho));
    output_arrays_2d.push_back(Output_reference<array2d>(JX, &jx));
    output_arrays_2d.push_back(Output_reference<array2d>(JY, &jy));
    output_arrays_2d.push_back(Output_reference<array2d>(JZ, &jz));
    output_arrays_2d.push_back(Output_reference<array2d>(EX, &ex));
    output_arrays_2d.push_back(Output_reference<array2d>(EY, &ey));
    output_arrays_2d.push_back(Output_reference<array2d>(EZ, &ez));
    output_arrays_2d.push_back(Output_reference<array2d>(BY, &by));
    output_arrays_2d.push_back(Output_reference<array2d>(BZ, &bz));

    for (auto & output_arr : output_arrays_2d) {
        output_writer.initialize_slice_array(n, d, *(output_arr.ptr), output_arr.name);
    }

    int particle_number = particles.size();

    // rho_ion deposition
    #pragma omp parallel for
    for (int pi = 0; pi < particle_number; pi++) {
        auto & p = particles[pi];
        deposit(p.y, p.z, -p.n, rho_ion);
    }


    bool stop_flag = false;
    const double THRESHOLD_B = 100;

    for (int i = 0; i < n.x; i++) {

        out << "Slice " << i << std::endl;

        #pragma omp parallel for
        for (int j = 0; j < n.y; j++) {
            for (int k = 0; k < n.z; k++) {
                psi_prev(j, k) = psi(j, k);
            }
        }

        #pragma omp parallel for
        for (int j = 0; j < n.y; j++) {
            for (int k = 0; k < n.z; k++) {
                fourier.in[n.z * j + k] = 0.0;
            }
        }

        // psi source deposition
        #pragma omp parallel for
        for (int pi = 0; pi < particle_number; pi++) {
            auto & p = particles[pi];
            deposit(p.y, p.z, -p.n, fourier.in, {d.y, d.z}, {n.y, n.z});
        }

        // cacluate psi
        #pragma omp parallel for
        for (int j = 0; j < n.y; j++) {
            for (int k = 0; k < n.z; k++) {
                fourier.in[n.z * j + k] -= rho_ion(j, k);
            }
        }
        solve_poisson_equation();

        #pragma omp parallel for
        for (int j = 0; j < n.y; j++) {
            for (int k = 0; k < n.z; k++) {
                psi(j, k) = (fourier.in[n.z * j + k]) / n.y / n.z;
            }
        }
        increase_minimum(psi, psi_threshold - 1);

        // calculate initial gamma, px
        #pragma omp parallel for
        for (int pi = 0; pi < particle_number; pi++) {
            auto & p = particles[pi];
            double a_particle = array_to_particle(p.y, p.z, a_sqr, i);
            double psi_particle = array_to_particle(p.y, p.z, psi);
            p.gamma = 0.5 * (1 + p.py * p.py + p.pz * p.pz + a_particle + (1 + psi_particle) * (1 + psi_particle))
                    / (1 + psi_particle);
            p.px = 0.5 * (1 + p.py * p.py + p.pz * p.pz + a_particle - (1 + psi_particle) * (1 + psi_particle))
                    / (1 + psi_particle);
            assert(p.gamma >= 1.0);
        }

        // initial jx, rho deposition
        #pragma omp parallel for
        for (int j = 0; j < n.y; j++) {
            for (int k = 0; k < n.z; k++) {
                jx(j, k) = rhobunch(i * d.x, j * d.y, k * d.z);
                rho(j, k) = rhobunch(i * d.x, j * d.y, k * d.z);
                susceptibility(i, j, k) = 0.0;
            }
        }

        #pragma omp parallel for
        for (int pi = 0; pi < particle_number; pi++) {
            auto & p = particles[pi];
            double vx = p.px / p.gamma;
            deposit(p.y, p.z, p.n * vx / (1 - vx), jx);
            deposit(p.y, p.z, p.n / (1 - vx), rho);
            deposit(p.y, p.z, -p.n / (1 - vx) / p.gamma, susceptibility, i);
        }

        if (i == n.x - 1) {
            break;
        }

        for (int iteration = 0; iteration < magnetic_field_iterations; iteration++) {

            // advance momenta
            #pragma omp parallel for
            for (int pi = 0; pi < particle_number; pi++) {
                auto & p = particles[pi];
                double psi_particle = array_to_particle(p.y, p.z, psi);
                double da_dy_particle = array_yder_to_particle(p.y, p.z, a_sqr, i);
                double da_dz_particle = array_zder_to_particle(p.y, p.z, a_sqr, i);
                double dpsi_dy_particle = array_yder_to_particle(p.y, p.z, psi);
                double dpsi_dz_particle = array_zder_to_particle(p.y, p.z, psi);
                double by_particle = array_to_particle(p.y, p.z, by);
                double bz_particle = array_to_particle(p.y, p.z, bz);

                p.py_next = p.py - d.x * 0.5 * da_dy_particle / (1 + psi_particle);
                p.py_next += d.x * p.gamma * dpsi_dy_particle / (1 + psi_particle);
                p.py_next -= d.x * bz_particle;
                p.pz_next = p.pz - d.x * 0.5 * da_dz_particle / (1 + psi_particle);
                p.pz_next += d.x * p.gamma * dpsi_dz_particle / (1 + psi_particle);
                p.pz_next += d.x * by_particle;
            }

            // advance half coordinate
            #pragma omp parallel for
            for (int pi = 0; pi < particle_number; pi++) {
                auto & p = particles[pi];
                double psi_particle = array_to_particle(p.y, p.z, psi);
                p.y_middle = p.y + 0.5 * d.x * p.py_next / (1 + psi_particle);
                p.z_middle = p.z + 0.5 * d.x * p.pz_next / (1 + psi_particle);
                normalize_coordinates(p.y_middle, p.z_middle);
            }

            // psi_source middle deposition
            #pragma omp parallel for
            for (int j = 0; j < n.y; j++) {
                for (int k = 0; k < n.z; k++) {
                    fourier.in[n.z * j + k] = - rho_ion(j, k);
                }
            }

            #pragma omp parallel for
            for (int pi = 0; pi < particle_number; pi++) {
                auto & p = particles[pi];
                deposit(p.y_middle, p.z_middle, -p.n, fourier.in, {d.y, d.z}, {n.y, n.z});
            }

            // calculate psi_middle
            solve_poisson_equation();

            #pragma omp parallel for
            for (int j = 0; j < n.y; j++) {
                for (int k = 0; k < n.z; k++) {
                    psi_middle(j, k) = (fourier.in[n.z * j + k]) / n.y / n.z;
                }
            }
            increase_minimum(psi_middle, psi_threshold - 1);

            // deposit jy, jz
            #pragma omp parallel for
            for (int j = 0; j < n.y; j++) {
                for (int k = 0; k < n.z; k++){
                    jy_next(j, k) = 0.0;
                    jz_next(j, k) = 0.0;
                }
            }

            #pragma omp parallel for
            for (int pi = 0; pi < particle_number; pi++) {
                auto & p = particles[pi];
                double psi_particle = array_to_particle(p.y_middle, p.z_middle, psi_middle);
                deposit(p.y_middle, p.z_middle, p.n * p.py_next / (1 + psi_particle), jy_next);
                deposit(p.y_middle, p.z_middle, p.n * p.pz_next / (1 + psi_particle), jz_next);
            }

            // calculate new djy_dxi, djz_dxi
            #pragma omp parallel for
            for (int j = 0; j < n.y; j++) {
                for (int k = 0; k < n.z; k++) {
                    djy_dxi(j, k) = (jy(j, k) - jy_next(j, k)) / d.x;
                    djz_dxi(j, k) = (jz(j, k) - jz_next(j, k)) / d.x;
                }
            }

            // calculate new gamma, px
            #pragma omp parallel for
            for (int pi = 0; pi < particle_number; pi++) {
                auto & p = particles[pi];
                double a_particle = array_to_particle(p.y, p.z, a_sqr, i);
                double psi_particle = array_to_particle(p.y, p.z, psi);
                double p_squared = 0.25 * (p.py + p.py_next) * (p.py + p.py_next) + 0.25 * (p.pz + p.pz_next) * (p.pz + p.pz_next);
                p.gamma = 0.5 * (1 + p_squared + a_particle +
                        (1 + psi_particle) * (1 + psi_particle)) / (1 + psi_particle);
                p.px = 0.5 * (1 + p_squared + a_particle - (1 + psi_particle) * (1 + psi_particle)) / (1 + psi_particle);
                assert(p.gamma >= 1.0);
            }

            // new jx, rho deposition
            #pragma omp parallel for
            for (int j = 0; j < n.y; j++) {
                for (int k = 0; k < n.z; k++) {
                    jx(j, k) = rhobunch(i * d.x, j * d.y, k * d.z);
                    rho(j, k) = rhobunch(i * d.x, j * d.y, k * d.z);
                    susceptibility(i, j, k) = 0.0;
                }
            }

            #pragma omp parallel for
            for (int pi = 0; pi < particle_number; pi++) {
                auto & p = particles[pi];
                double vx = p.px / p.gamma;
                deposit(p.y, p.z, p.n * vx / (1 - vx), jx);
                deposit(p.y, p.z, p.n / (1 - vx), rho);
                deposit(p.y, p.z, -p.n / (1 - vx) / p.gamma, susceptibility, i);
            }

            // new source for B_y
            #pragma omp parallel for
            for (int j = 0; j < n.y; j++) {
                for (int k = 0; k < n.z-1; k++) {
                    fourier.in[n.z * j + k] = - magnetic_field_D * by(j, k) - (jx(j, k+1) - jx(j, k)) / d.z + djz_dxi(j, k);
                }
                fourier.in[n.z * j + (n.z-1)] = - magnetic_field_D * by(j, n.z-1) - (jx(j, 0) - jx(j, n.z-1)) / d.z + djz_dxi(j, n.z-1);
            }

            // new guess for B_y
            solve_poisson_equation(magnetic_field_D);

            #pragma omp parallel for shared(stop_flag)
            for (int j = 0; j < n.y; j++) {
                for (int k = 0; k < n.z; k++) {
                    by(j, k) = fourier.in[n.z * j + k] / n.y / n.z;
                    if (fabs(by(j, k)) > THRESHOLD_B) {
                        stop_flag = true;
                    }
                    if (std::isnan(by(j, k))) {
                        stop_flag = true;
                    }
                }
            }

            // new source for B_z
            #pragma omp parallel for
            for (int j = 0; j < n.y-1; j++) {
                for (int k = 0; k < n.z; k++) {
                    fourier.in[n.z * j + k] = - magnetic_field_D * bz(j, k) + (jx(j+1, k) - jx(j, k)) / d.y - djy_dxi(j, k);
                }
            }
            #pragma omp parallel for
            for (int k = 0; k < n.z; k++) {
                fourier.in[n.z * (n.y-1) + k] = - magnetic_field_D * bz(0, k) + (jx(0, k) - jx(n.y-1, k)) / d.y - djy_dxi(n.y-1, k);
            }

            // new guess for B_z
            solve_poisson_equation(magnetic_field_D);

            #pragma omp parallel for shared(stop_flag)
            for (int j = 0; j < n.y; j++) {
                for (int k = 0; k < n.z; k++) {
                    bz(j, k) = fourier.in[n.z * j + k] / n.y / n.z;
                    if (fabs(bz(j, k)) > THRESHOLD_B) {
                        stop_flag = true;
                    }
                    if (std::isnan(bz(j, k))) {
                        stop_flag = true;
                    }
                }
            }

            if (stop_flag) {
                break;
            }

        }

        if (stop_flag) {
            out << "Stop condition encountered, stopping the simulation" << std::endl;
            break;
        }

        #pragma omp parallel for
        for (int pi = 0; pi < particle_number; pi++) {
            auto & p = particles[pi];
            p.py = p.py_next;
            p.pz = p.pz_next;
        }

        // advance coordinate
        #pragma omp parallel for
        for (int pi = 0; pi < particle_number; pi++) {
            auto & p = particles[pi];
            double psi_particle = array_to_particle(p.y, p.z, psi_middle);
            p.y += d.x * p.py / (1 + psi_particle);
            p.z += d.x * p.pz / (1 + psi_particle);
            normalize_coordinates(p.y, p.z);
        }

        output_step(output_writer, output_arrays_2d, i);

        #pragma omp parallel for
        for (int j = 0; j < n.y; j++) {
            for (int k = 0; k < n.z; k++) {
                jy(j, k) = jy_next(j, k);
                jz(j, k) = jz_next(j, k);
            }
        }
    }

    output_step(output_writer, output_arrays_2d, n.x - 1);

    for (auto & output_arr : output_arrays_3d) {
        output_writer.write_array(*(output_arr.ptr), output_arr.name);
    }
}

void System_3d::output_step(Output_writer & output_writer, 
                            const std::vector<Output_reference<array2d>> & output_arrays_2d, int slice_index) {
    const int i = slice_index;
    
    // calculate ex
    if (i > 0) {
        #pragma omp parallel for
        for (int j = 0; j < n.y; j++) {
            for (int k = 0; k < n.z; k++) {
                ex(j, k) = (psi(j, k) - psi_prev(j, k)) / d.x;
            }
        }
    }

    // calculate ey
    #pragma omp parallel for
    for (int j = 0; j < n.y-1; j++) {
        for (int k = 0; k < n.z; k++) {
            ey(j, k) = -(psi(j+1, k) - psi(j, k)) / d.y + bz(j, k);
        }
    }
    for (int k = 0; k < n.z; k++) {
        ey(n.y-1, k) = -(psi(0, k) - psi(n.y-1, k)) / d.y + bz(n.y-1, k);
    }
    
    // calculate ez
    #pragma omp parallel for
    for (int j = 0; j < n.y; j++) {
        for (int k = 0; k < n.z-1; k++) {
            ez(j, k) = -(psi(j, k+1) - psi(j, k)) / d.z - by(j, k);
        }
        ez(j, n.z-1) = -(psi(j, 0) - psi(j, n.z-1)) / d.z - by(j, n.z-1);
    }
    
    for (auto & output_arr : output_arrays_2d) {
        output_writer.write_slice(*(output_arr.ptr), output_arr.name, slice_index);
    }
}

void System_3d::normalize_coordinates(double & y, double & z) {
    y = fmod(y, l.y);
    if (y < 0) {
        y += l.y;
    }

    assert(y >= 0);
    assert(y < l.y);

    z = fmod(z, l.z);
    if (z < 0) {
        z += l.z;
    }
    
    assert(z >= 0);
    assert(z < l.z);
}

void System_3d::solve_poisson_equation(double D) {
    fourier.forward_transform();
    fourier.out[0][0] = 0;
    fourier.out[0][1] = 0;
    const int nz_size = (n.z / 2) + 1;

    #pragma omp parallel for
    for (int j = 0; j < n.y; j++) {
        for (int k = 0; k < nz_size; k++){
            if ((k == 0) && (j == 0)) {
                continue;
            }

            double multiplier = -D + (2 * cos(2 * M_PI * j / n.y) - 2) / d.y / d.y + (2 * cos(2 * M_PI * k / n.z) - 2) / d.z / d.z;

            fourier.out[nz_size * j + k][0] /= multiplier;
            fourier.out[nz_size * j + k][1] /= multiplier;
        }
    }
    fourier.backward_transform();
}

void System_3d::init_particles(int ppcy, int ppcz, std::function<double(double, double)> plasma_profile) {
    assert(ppcy > 0);
    assert(ppcz > 0);
    particles = std::vector<particle>(n.y * n.z * ppcy * ppcz);

    for (int i = 0; i < ppcy * n.y; i++) {
        for (int j = 0; j < ppcz * n.z; j++) {
            int index = ppcz * n.z * i + j;
            const double y = (i + 0.5) * d.y / ppcy;
            const double z = (j + 0.5) * d.z / ppcz;
            double value = plasma_profile(y, z);
            particles[index].y = y;
            particles[index].z = z;
            particles[index].n = value / ppcy / ppcz;
        }
    }
}

void System_3d::init_a_sqr(std::function<double(double, double, double)> func) {
    for (int i = 0; i < n.x; i++) {
        for (int j = 0; j < n.y; j++) {
            for (int k = 0; k < n.z; k++) {
                double x = i * d.x;
                double y = j * d.y;
                double z = k * d.z;
                a_sqr(i, j, k) = func(x, y, z);
            }
        }
    }
}

void System_3d::increase_minimum(array3d & array, int slice, double value) const {
    auto n = array.get_dimensions();
    #pragma omp parallel for
    for (int j = 0; j < n.y; j++) {
        for (int k = 0; k < n.z; k++) {
            if (array(slice, j, k) < value) {
                array(slice, j, k) = value;
            }
        }
    }
}

void System_3d::increase_minimum(array2d & array, double value) const {
    auto n1 = array.get_n1();
    auto n2 = array.get_n2();
    #pragma omp parallel for
    for (int j = 0; j < n1; j++) {
        for (int k = 0; k < n2; k++) {
            if (array(j, k) < value) {
                array(j, k) = value;
            }
        }
    }
}

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <H5Cpp.h>
#include <cassert>
#include <fmt/format.h>
#include <vector>
#include <string>

#include "System_3d.h"
#include "containers_3d.h"
#include "output.h"

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
    const long array2d_memory = 3l * sizeof(double) * n.y * n.z;
    const long array3d_memory = 11l * sizeof(double) * n.x * n.y * n.z;
    const long particle_memory = static_cast<long>(sizeof(particle)) * params.ppcy * params.ppcz * n.y * n.z;
    const long total_memory = array2d_memory + array3d_memory + particle_memory + fourier_memory;

    out << "Expected RAM usage:\n";
    out << fmt::format("3D arrays: {}\n", memory_formatter(array3d_memory))
        << fmt::format("2D arrays: {}\n", memory_formatter(array2d_memory))
        << fmt::format("Fourier:   {}\n", memory_formatter(fourier_memory))
        << fmt::format("Particles: {}\n", memory_formatter(particle_memory))
        << fmt::format("Total:     {}", memory_formatter(total_memory)) << std::endl;

    std::cout << "----------------------------------------" << std::endl;
    
    init_particles(params.ppcy, params.ppcz);

    fourier = Fourier2d(n.y, n.z);

    psi_middle = array2d({n.y, n.z});
    djy_dxi = array2d({n.y, n.z});
    djz_dxi = array2d({n.y, n.z});

    psi = array3d(n, d);
    a_sqr = array3d(n, d);
    jx = array3d(n, d);
    jy = array3d(n, d, {-0.5 * d.x, 0.5 * d.y, 0});
    jz = array3d(n, d, {-0.5 * d.x, 0, 0.5 * d.z});
    rho = array3d(n, d);
    ex = array3d(n, d, {-0.5 * d.x, 0, 0});
    ey = array3d(n, d, {0, 0.5 * d.y, 0});
    ez = array3d(n, d, {0, 0, 0.5 * d.z});
    by = array3d(n, d, {0, 0, 0.5 * d.z});
    bz = array3d(n, d, {0, 0.5 * d.y, 0});

    init_a_sqr(params.a_sqr);
}

void System_3d::solve_wakefield() {
    int particle_number = particles.size();

    bool stop_flag = false;
    const double THRESHOLD_B = 100;

    for (int i = 0; i < n.x; i++) {

        out << "Slice " << i << std::endl;

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
            deposit(p.y, p.z, -p.n, fourier.in);
        }

        // cacluate psi
        #pragma omp parallel for
        for (int j = 0; j < n.y; j++) {
            for (int k = 0; k < n.z; k++) {
                fourier.in[n.z * j + k] -= 1.0;
            }
        }
        solve_poisson_equation();

        #pragma omp parallel for
        for (int j = 0; j < n.y; j++) {
            for (int k = 0; k < n.z; k++) {
                psi(i, j, k) = (fourier.in[n.z * j + k]) / n.y / n.z;
            }
        }
        increase_minimum(psi, i, psi_threshold - 1);

        // calculate initial gamma, px
        #pragma omp parallel for
        for (int pi = 0; pi < particle_number; pi++) {
            auto & p = particles[pi];
            double a_particle = array_to_particle(p.y, p.z, a_sqr, i);
            double psi_particle = array_to_particle(p.y, p.z, psi, i);
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
                jx(i, j, k) = rhobunch(i * d.x, j * d.y, k * d.z);
                rho(i, j, k) = rhobunch(i * d.x, j * d.y, k * d.z);
            }
        }

        #pragma omp parallel for
        for (int pi = 0; pi < particle_number; pi++) {
            auto & p = particles[pi];
            double vx = p.px / p.gamma;
            deposit(p.y, p.z, p.n * vx / (1 - vx), jx, i);
            deposit(p.y, p.z, p.n / (1 - vx), rho, i);
        }

        if (i > 0) {
            #pragma omp parallel for
            for (int j = 0; j < n.y; j++) {
                for (int k = 0; k < n.z; k++) {
                    by(i, j, k) = by(i-1, j, k);
                    bz(i, j, k) = bz(i-1, j, k);
                }
            }
        }

        if (i == n.x - 1) {
            break;
        }

        for (int iteration = 0; iteration < magnetic_field_iterations; iteration++) {

            // advance momenta
            #pragma omp parallel for
            for (int pi = 0; pi < particle_number; pi++) {
                auto & p = particles[pi];
                double psi_particle = array_to_particle(p.y, p.z, psi, i);
                double da_dy_particle = array_yder_to_particle(p.y, p.z, a_sqr, i);
                double da_dz_particle = array_zder_to_particle(p.y, p.z, a_sqr, i);
                double dpsi_dy_particle = array_yder_to_particle(p.y, p.z, psi, i);
                double dpsi_dz_particle = array_zder_to_particle(p.y, p.z, psi, i);
                double by_particle = array_to_particle(p.y, p.z, by, i, 0.0, 0.5);
                double bz_particle = array_to_particle(p.y, p.z, bz, i, 0.5, 0.0);

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
                double psi_particle = array_to_particle(p.y, p.z, psi, i);
                p.y_middle = p.y + 0.5 * d.x * p.py_next / (1 + psi_particle);
                p.z_middle = p.z + 0.5 * d.x * p.pz_next / (1 + psi_particle);
                normalize_coordinates(p.y_middle, p.z_middle);
            }

            // psi_source middle deposition
            #pragma omp parallel for
            for (int j = 0; j < n.y; j++) {
                for (int k = 0; k < n.z; k++) {
                    fourier.in[n.z * j + k] = -1.0; // rho_ion
                }
            }

            #pragma omp parallel for
            for (int pi = 0; pi < particle_number; pi++) {
                auto & p = particles[pi];
                deposit(p.y_middle, p.z_middle, -p.n, fourier.in);
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
                    jy(i+1, j, k) = 0.0;
                    jz(i+1, j, k) = 0.0;
                }
            }

            #pragma omp parallel for
            for (int pi = 0; pi < particle_number; pi++) {
                auto & p = particles[pi];
                double psi_particle = array_to_particle(p.y_middle, p.z_middle, psi_middle);
                deposit(p.y_middle, p.z_middle, p.n * p.py_next / (1 + psi_particle), jy, i+1, 0.5, 0.0);
                deposit(p.y_middle, p.z_middle, p.n * p.pz_next / (1 + psi_particle), jz, i+1, 0.0, 0.5);
            }

            // calculate new djy_dxi, djz_dxi
            #pragma omp parallel for
            for (int j = 0; j < n.y; j++) {
                for (int k = 0; k < n.z; k++) {
                    djy_dxi(j, k) = (jy(i, j, k) - jy(i+1, j, k)) / d.x;
                    djz_dxi(j, k) = (jz(i, j, k) - jz(i+1, j, k)) / d.x;
                }
            }

            // calculate new gamma, px
            #pragma omp parallel for
            for (int pi = 0; pi < particle_number; pi++) {
                auto & p = particles[pi];
                double a_particle = array_to_particle(p.y, p.z, a_sqr, i);
                double psi_particle = array_to_particle(p.y, p.z, psi, i);
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
                    jx(i, j, k) = rhobunch(i * d.x, j * d.y, k * d.z);
                    rho(i, j, k) = rhobunch(i * d.x, j * d.y, k * d.z);
                }
            }

            #pragma omp parallel for
            for (int pi = 0; pi < particle_number; pi++) {
                auto & p = particles[pi];
                double vx = p.px / p.gamma;
                deposit(p.y, p.z, p.n * vx / (1 - vx), jx, i);
                deposit(p.y, p.z, p.n / (1 - vx), rho, i);
            }

            // new source for B_y
            #pragma omp parallel for
            for (int j = 0; j < n.y; j++) {
                for (int k = 0; k < n.z-1; k++) {
                    fourier.in[n.z * j + k] = - magnetic_field_D * by(i, j, k) - (jx(i, j, k+1) - jx(i, j, k)) / d.z + djz_dxi(j, k);
                }
                fourier.in[n.z * j + (n.z-1)] = - magnetic_field_D * by(i, j, n.z-1) - (jx(i, j, 0) - jx(i, j, n.z-1)) / d.z + djz_dxi(j, n.z-1);
            }

            // new guess for B_y
            solve_poisson_equation(magnetic_field_D);

            #pragma omp parallel for shared(stop_flag)
            for (int j = 0; j < n.y; j++) {
                for (int k = 0; k < n.z; k++) {
                    by(i, j, k) = fourier.in[n.z * j + k] / n.y / n.z;
                    if (fabs(by(i, j, k)) > THRESHOLD_B) {
                        stop_flag = true;
                    }
                    if (std::isnan(by(i, j, k))) {
                        stop_flag = true;
                    }
                }
            }

            // new source for B_z
            #pragma omp parallel for
            for (int j = 0; j < n.y-1; j++) {
                for (int k = 0; k < n.z; k++) {
                    fourier.in[n.z * j + k] = - magnetic_field_D * bz(i, j, k) + (jx(i, j+1, k) - jx(i, j, k)) / d.y - djy_dxi(j, k);
                }
            }
            #pragma omp parallel for
            for (int k = 0; k < n.z; k++) {
                fourier.in[n.z * (n.y-1) + k] = - magnetic_field_D * bz(i, 0, k) + (jx(i, 0, k) - jx(i, n.y-1, k)) / d.y - djy_dxi(n.y-1, k);
            }

            // new guess for B_z
            solve_poisson_equation(magnetic_field_D);

            #pragma omp parallel for shared(stop_flag)
            for (int j = 0; j < n.y; j++) {
                for (int k = 0; k < n.z; k++) {
                    bz(i, j, k) = fourier.in[n.z * j + k] / n.y / n.z;
                    if (fabs(bz(i, j, k)) > THRESHOLD_B) {
                        stop_flag = true;
                    }
                    if (std::isnan(bz(i, j, k))) {
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
    }

    // calculate ex
    #pragma omp parallel for
    for (int i = 1; i < n.x; i++) {
        for (int j = 0; j < n.y; j++) {
            for (int k = 0; k < n.z; k++) {
                ex(i, j, k) = (psi(i, j, k) - psi(i-1, j, k)) / d.x;
            }
        }
    }

    // calculate ey
    #pragma omp parallel for
    for (int i = 1; i < n.x; i++) {
        for (int j = 0; j < n.y-1; j++) {
            for (int k = 0; k < n.z; k++) {
                ey(i, j, k) = -(psi(i, j+1, k) - psi(i, j, k)) / d.y + bz(i, j, k);
            }
        }
        for (int k = 0; k < n.z; k++) {
            ey(i, n.y-1, k) = -(psi(i, 0, k) - psi(i, n.y-1, k)) / d.y + bz(i, n.y-1, k);
        }
    }

    // calculate ez
    #pragma omp parallel for
    for (int i = 1; i < n.x; i++) {
        for (int j = 0; j < n.y; j++) {
            for (int k = 0; k < n.z-1; k++) {
                ez(i, j, k) = -(psi(i, j, k+1) - psi(i, j, k)) / d.z - by(i, j, k);
            }
            ez(i, j, n.z-1) = -(psi(i, j, 0) - psi(i, j, n.z-1)) / d.z - by(i, j, n.z-1);
        }
    }
}

void System_3d::output() const {
    if (output_parameters.output3d) {
        H5::H5File fields_file("Fields.h5", H5F_ACC_TRUNC);

        write_array(a_sqr, "aSqr", fields_file);
        write_array(rho, "rho", fields_file);
        write_array(jx, "jx", fields_file);
        write_array(jy, "jy", fields_file);
        write_array(jz, "jz", fields_file);
        write_array(psi, "psi", fields_file);
        write_array(ex, "ex", fields_file);
        write_array(ey, "ey", fields_file);
        write_array(ez, "ez", fields_file);
        write_array(by, "by", fields_file);
        write_array(bz, "bz", fields_file);
    }
    if (output_parameters.output_xy) {
        H5::H5File fields_xy_file("Fields_xy.h5", H5F_ACC_TRUNC);

        const double z0 = 0.5 * l.z;
        write_array(calculate_xy_slice(a_sqr, z0), "aSqr", fields_xy_file);
        write_array(calculate_xy_slice(rho, z0), "rho", fields_xy_file);
        write_array(calculate_xy_slice(jx, z0), "jx", fields_xy_file);
        write_array(calculate_xy_slice(jy, z0), "jy", fields_xy_file);
        write_array(calculate_xy_slice(jz, z0), "jz", fields_xy_file);
        write_array(calculate_xy_slice(psi, z0), "psi", fields_xy_file);
        write_array(calculate_xy_slice(ex, z0), "ex", fields_xy_file);
        write_array(calculate_xy_slice(ey, z0), "ey", fields_xy_file);
        write_array(calculate_xy_slice(ez, z0), "ez", fields_xy_file);

        write_array(calculate_xy_slice(by, z0), "by", fields_xy_file);
        write_array(calculate_xy_slice(bz, z0), "bz", fields_xy_file);
    }
}

void System_3d::deposit(double y, double z, double value, array3d & array, int slice, double yshift, double zshift) {
    int j1 = (int) floor(y / d.y - yshift);
    int j2 = (j1 + 1) % n.y;

    double y_frac = (y / d.y - yshift) - j1;
    if (j1 < 0) {
        j1 += n.y;
    }

    int k1 = (int) floor(z / d.z - zshift);
    int k2 = (k1 + 1) % n.z;
    double z_frac = (z / d.z - zshift) - k1;
    if (k1 < 0) {
        k1 += n.z;
    }

    assert((j1 >= 0) && (j1 < n.y));
    assert((j2 >= 0) && (j2 < n.y));
    assert((k1 >= 0) && (k1 < n.z));
    assert((k2 >= 0) && (k2 < n.z));

    #pragma omp atomic update
    array(slice, j1, k1) += value * (1 - y_frac) * (1 - z_frac);
    #pragma omp atomic update
    array(slice, j2, k1) += value * y_frac * (1 - z_frac);
    #pragma omp atomic update
    array(slice, j1, k2) += value * (1 - y_frac) * z_frac;
    #pragma omp atomic update
    array(slice, j2, k2) += value * y_frac * z_frac;
}

void System_3d::deposit(double y, double z, double value, double * array) {
    int j1 = (int) (y / d.y);
    int j2 = (j1 + 1) % n.y;
    double y_frac = (y / d.y) - j1;

    int k1 = (int) (z / d.z);
    int k2 = (k1 + 1) % n.z;
    double z_frac = (z / d.z) - k1;

    assert((j1 >= 0) && (j1 < n.y));
    assert((j2 >= 0) && (j2 < n.y));
    assert((k1 >= 0) && (k1 < n.z));
    assert((k2 >= 0) && (k2 < n.z));

    #pragma omp atomic update
    array[n.z * j1 + k1] += value * (1 - y_frac) * (1 - z_frac);
    #pragma omp atomic update
    array[n.z * j2 + k1] += value * y_frac * (1 - z_frac);
    #pragma omp atomic update
    array[n.z * j1 + k2] += value * (1 - y_frac) * z_frac;
    #pragma omp atomic update
    array[n.z * j2 + k2] += value * y_frac * z_frac;
}

double System_3d::array_to_particle(double y, double z, const array3d & array, int slice, double yshift, double zshift) const {
    int j1 = (int) floor(y / d.y - yshift);
    int j2 = (j1 + 1) % n.y;
    double y_frac = (y / d.y - yshift) - j1;
    if (j1 < 0) {
        j1 += n.y;
    }

    int k1 = (int) floor(z / d.z - zshift);
    int k2 = (k1 + 1) % n.z;
    double z_frac = (z / d.z - zshift) - k1;
    if (k1 < 0) {
        k1 += n.z;
    }

    return array(slice, j1, k1) * (1 - y_frac) * (1 - z_frac) + array(slice, j2, k1) * y_frac * (1 - z_frac) +
           array(slice, j1, k2) * (1 - y_frac) * z_frac + array(slice, j2, k2) * y_frac * z_frac;
}

double System_3d::array_to_particle(double y, double z, const array2d & array) const {
    int j1 = (int) (y / d.y);
    int j2 = (j1 + 1) % n.y;
    double y_frac = (y / d.y) - j1;

    int k1 = (int) (z / d.z);
    int k2 = (k1 + 1) % n.z;
    double z_frac = (z / d.z) - k1;

    return array(j1, k1) * (1 - y_frac) * (1 - z_frac) + array(j2, k1) * y_frac * (1 - z_frac) +
           array(j1, k2) * (1 - y_frac) * z_frac + array(j2, k2) * y_frac * z_frac;
}

double System_3d::array_yder_to_particle(double y, double z, const array3d & array, int slice) const {
    int j1 = (int) floor(y / d.y - 0.5);
    double y_frac = (y / d.y - 0.5) - j1;
    if (j1 < 0) {
        j1 += n.y;
    }
    int j2 = (j1 + 1) % n.y;
    int j3 = (j1 + 2) % n.y;

    int k1 = (int) (z / d.z);
    int k2 = (k1 + 1) % n.z;
    double z_frac = (z / d.z) - k1;

    double der_y1_z1 = (array(slice, j2, k1) - array(slice, j1, k1)) / d.y;
    double der_y2_z1 = (array(slice, j3, k1) - array(slice, j2, k1)) / d.y;
    double der_y1_z2 = (array(slice, j2, k2) - array(slice, j1, k2)) / d.y;
    double der_y2_z2 = (array(slice, j3, k2) - array(slice, j2, k2)) / d.y;

    return der_y1_z1 * (1 - y_frac) * (1 - z_frac) + der_y2_z1 * y_frac * (1 - z_frac) +
           der_y1_z2 * (1 - y_frac) * z_frac + der_y2_z2 * y_frac * z_frac;
}

double System_3d::array_zder_to_particle(double y, double z, const array3d & array, int slice) const {
    int j1 = (int) (y / d.y);
    int j2 = (j1 + 1) % n.y;
    double y_frac = (y / d.y) - j1;

    int k1 = (int) floor(z / d.z - 0.5);
    double z_frac = (z / d.z - 0.5) - k1;
    if (k1 < 0) {
        k1 += n.z;
    }
    int k2 = (k1 + 1) % n.z;
    int k3 = (k1 + 2) % n.z;

    double der_y1_z1 = (array(slice, j1, k2) - array(slice, j1, k1)) / d.z;
    double der_y2_z1 = (array(slice, j2, k2) - array(slice, j2, k1)) / d.z;
    double der_y1_z2 = (array(slice, j1, k3) - array(slice, j1, k2)) / d.z;
    double der_y2_z2 = (array(slice, j2, k3) - array(slice, j2, k2)) / d.z;

    return der_y1_z1 * (1 - y_frac) * (1 - z_frac) + der_y2_z1 * y_frac * (1 - z_frac) +
           der_y1_z2 * (1 - y_frac) * z_frac + der_y2_z2 * y_frac * z_frac;
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

void System_3d::init_particles(int ppcy, int ppcz) {
    assert(ppcy > 0);
    assert(ppcz > 0);
    particles = std::vector<particle>(n.y * n.z * ppcy * ppcz);

    for (int i = 0; i < ppcy * n.y; i++) {
        for (int j = 0; j < ppcz * n.z; j++) {
            int index = ppcz * n.z * i + j;
            particles[index].y = (i + 0.5) * d.y / ppcy;
            particles[index].z = (j + 0.5) * d.z / ppcz;
            particles[index].n = -1.0 / ppcy / ppcz;
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

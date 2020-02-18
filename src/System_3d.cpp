#include <cmath>
#include <iostream>
#include <stdexcept>
#include <H5Cpp.h>
#include <cassert>

#include "System_3d.h"
#include "containers_3d.h"
#include "output.h"

System_3d::System_3d(System_parameters & params, std::ostream & out) : 
    l(params.l),
    d(params.d),
    magnetic_field_iterations(params.magnetic_field_iterations),
    rhobunch(params.rho),
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

    n.x = static_cast<int>(l.x / d.x);
    n.y = static_cast<int>(l.y / d.y);
    n.z = static_cast<int>(l.z / d.z);

    init_particles(params.ppcy, params.ppcz);

    fourier = Fourier2d(n.y, n.z);

    psi_middle = array2d(n.y, n.z);
    djy_dxi = array2d(n.y, n.z);
    djz_dxi = array2d(n.y, n.z);

    psi = array3d(n);
    psi_source = array3d(n);
    a_sqr = array3d(n);
    jx = array3d(n);
    jy = array3d(n);
    jz = array3d(n);
    rho = array3d(n);
    ex = array3d(n);
    ey = array3d(n);
    ez = array3d(n);
    by = array3d(n);
    bz = array3d(n);

    init_a_sqr(params.a_sqr);
}

void System_3d::solve_wakefield() {
    for (int i = 0; i < n.x; i++) {
        std::cout << "Slice " << i << std::endl;

        // psi_source deposition
        for (auto & p : particles) {
            deposit(p.y, p.z, -p.n, psi_source, i);
        }

        // cacluate psi
        for (int j = 0; j < n.y; j++) {
            for (int k = 0; k < n.z; k++) {
                fourier.in[n.z * j + k] = psi_source(i, j, k) - 1.0;
            }
        }
        solve_poisson_equation();
        for (int j = 0; j < n.y; j++) {
            for (int k = 0; k < n.z; k++) {
                psi(i, j, k) = (fourier.in[n.z * j + k]) / n.y / n.z;
            }
        }

        // calculate initial gamma, px
        for (auto & p : particles) {
            double a_particle = array_to_particle(p.y, p.z, a_sqr, i);
            double psi_particle = array_to_particle(p.y, p.z, psi, i);
            p.gamma = 0.5 * (1 + p.py * p.py + p.pz * p.pz + a_particle + (1 + psi_particle) * (1 + psi_particle))
                    / (1 + psi_particle);
            p.px = 0.5 * (1 + p.py * p.py + p.pz * p.pz + a_particle - (1 + psi_particle) * (1 + psi_particle))
                    / (1 + psi_particle);
        }

        // initial jx, rho deposition
        for (int j = 0; j < n.y; j++) {
            for (int k = 0; k < n.z; k++) {
                jx(i, j, k) = rhobunch(i * d.x, j * d.y, k * d.z);
                rho(i, j, k) = rhobunch(i * d.x, j * d.y, k * d.z);
            }
        }

        for (auto & p : particles) {
            double vx = p.px / p.gamma;
            deposit(p.y, p.z, p.n * vx / (1 - vx), jx, i);
            deposit(p.y, p.z, p.n / (1 - vx), rho, i);
        }

        if (i > 0) {
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
            for (auto & p : particles) {
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
            for (auto & p : particles) {
                double psi_particle = array_to_particle(p.y, p.z, psi, i);
                p.y_middle = p.y + 0.5 * d.x * p.py_next / (1 + psi_particle);
                p.z_middle = p.z + 0.5 * d.x * p.pz_next / (1 + psi_particle);
                normalize_coordinates(p.y_middle, p.z_middle);
            }

            // psi_source middle deposition
            for (int j = 0; j < n.y; j++) {
                for (int k = 0; k < n.z; k++) {
                    fourier.in[n.z * j + k] = -1.0; // rho_ion
                }
            }

            for (auto & p : particles) {
                deposit(p.y_middle, p.z_middle, -p.n, fourier.in);
            }

            // calculate psi_middle
            solve_poisson_equation();
            for (int j = 0; j < n.y; j++) {
                for (int k = 0; k < n.z; k++) {
                    psi_middle(j, k) = (fourier.in[n.z * j + k]) / n.y / n.z;
                }
            }

            // deposit jy, jz
            for (int j = 0; j < n.y; j++) {
                for (int k = 0; k < n.z; k++){
                    jy(i+1, j, k) = 0.0;
                    jz(i+1, j, k) = 0.0;
                }
            }

            for (auto & p : particles) {
                double psi_particle = array_to_particle(p.y_middle, p.z_middle, psi_middle);
                deposit(p.y_middle, p.z_middle, p.n * p.py_next / (1 + psi_particle), jy, i+1, 0.5, 0.0);
                deposit(p.y_middle, p.z_middle, p.n * p.pz_next / (1 + psi_particle), jz, i+1, 0.0, 0.5);
            }

            // calculate new djy_dxi, djz_dxi
            for (int j = 0; j < n.y; j++) {
                for (int k = 0; k < n.z; k++) {
                    djy_dxi(j, k) = (jy(i, j, k) - jy(i+1, j, k)) / d.x;
                    djz_dxi(j, k) = (jz(i, j, k) - jz(i+1, j, k)) / d.x;
                }
            }

            // calculate new gamma, px
            for (auto & p : particles) {
                double a_particle = array_to_particle(p.y, p.z, a_sqr, i);
                double psi_particle = array_to_particle(p.y, p.z, psi, i);
                double p_squared = 0.25 * (p.py + p.py_next) * (p.py + p.py_next) + 0.25 * (p.pz + p.pz_next) * (p.pz + p.pz_next);
                p.gamma = 0.5 * (1 + p_squared + a_particle +
                        (1 + psi_particle) * (1 + psi_particle)) / (1 + psi_particle);
                p.px = 0.5 * (1 + p_squared + a_particle - (1 + psi_particle) * (1 + psi_particle)) / (1 + psi_particle);
            }

            // new jx, rho deposition
            for (int j = 0; j < n.y; j++) {
                for (int k = 0; k < n.z; k++) {
                    jx(i, j, k) = rhobunch(i * d.x, j * d.y, k * d.z);
                    rho(i, j, k) = rhobunch(i * d.x, j * d.y, k * d.z);
                }
            }

            for (auto & p : particles) {
                double vx = p.px / p.gamma;
                deposit(p.y, p.z, p.n * vx / (1 - vx), jx, i);
                deposit(p.y, p.z, p.n / (1 - vx), rho, i);
            }

            // new source for B_y
            for (int j = 0; j < n.y; j++) {
                for (int k = 0; k < n.z-1; k++) {
                    fourier.in[n.z * j + k] = -by(i, j, k) - (jx(i, j, k+1) - jx(i, j, k)) / d.z + djz_dxi(j, k);
                }
                fourier.in[n.z * j + (n.z-1)] = -by(i, j, n.z-1) - (jx(i, j, 0) - jx(i, j, n.z-1)) / d.z + djz_dxi(j, n.z-1);
            }

            // new guess for B_y
            solve_poisson_equation(1.0);

            for (int j = 0; j < n.y; j++) {
                for (int k = 0; k < n.z; k++) {
                    by(i, j, k) = fourier.in[n.z * j + k] / n.y / n.z;
                }
            }

            // new source for B_z
            for (int j = 0; j < n.y-1; j++) {
                for (int k = 0; k < n.z; k++) {
                    fourier.in[n.z * j + k] = -bz(i, j, k) + (jx(i, j+1, k) - jx(i, j, k)) / d.y - djy_dxi(j, k);
                }
            }
            for (int k = 0; k < n.z; k++) {
                fourier.in[n.z * (n.y-1) + k] = -bz(i, 0, k) + (jx(i, 0, k) - jx(i, n.y-1, k)) / d.y - djy_dxi(n.y-1, k);
            }

            // new guess for B_z
            solve_poisson_equation(1.0);

            for (int j = 0; j < n.y; j++) {
                for (int k = 0; k < n.z; k++) {
                    bz(i, j, k) = fourier.in[n.z * j + k] / n.y / n.z;
                }
            }
        }

        for (auto & p : particles) {
            p.py = p.py_next;
            p.pz = p.pz_next;
        }

        // advance coordinate
        for (auto & p : particles) {
            double psi_particle = array_to_particle(p.y, p.z, psi_middle);
            p.y += d.x * p.py / (1 + psi_particle);
            p.z += d.x * p.pz / (1 + psi_particle);
            normalize_coordinates(p.y, p.z);
        }
    }

    // calculate ex
    for (int i = 1; i < n.x; i++) {
        for (int j = 0; j < n.y; j++) {
            for (int k = 0; k < n.z; k++) {
                ex(i, j, k) = (psi(i, j, k) - psi(i-1, j, k)) / d.x;
            }
        }
    }

    // calculate ey
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
    H5::H5File fields_file("Fields.h5", H5F_ACC_TRUNC);

    write_array(a_sqr, "aSqr", fields_file);
    write_array(psi_source, "jx_minus_rho", fields_file);
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

    array(slice, j1, k1) += value * (1 - y_frac) * (1 - z_frac);
    array(slice, j2, k1) += value * y_frac * (1 - z_frac);
    array(slice, j1, k2) += value * (1 - y_frac) * z_frac;
    array(slice, j2, k2) += value * y_frac * z_frac;
}

void System_3d::deposit(double y, double z, double value, double * array) {
    int j1 = (int) (y / d.y);
    int j2 = (j1 + 1) % n.y;
    double y_frac = (y / d.y) - j1;

    int k1 = (int) (z / d.z);
    int k2 = (k1 + 1) % n.z;
    double z_frac = (z / d.z) - k1;

    array[n.z * j1 + k1] += value * (1 - y_frac) * (1 - z_frac);
    array[n.z * j2 + k1] += value * y_frac * (1 - z_frac);
    array[n.z * j1 + k2] += value * (1 - y_frac) * z_frac;
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
    if (y < 0) {
        y += l.y;
    } else if (y >= l.y) {
        y -= l.y;
    }

    if (z < 0) {
        z += l.z;
    } else if (z >= l.z) {
        z -= l.z;
    }
}

void System_3d::solve_poisson_equation(double D) {
    fourier.forward_transform();
    fourier.out[0][0] = 0;
    fourier.out[0][1] = 0;
    const int nz_size = (n.z / 2) + 1;

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
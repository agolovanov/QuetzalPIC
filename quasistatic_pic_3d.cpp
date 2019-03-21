#include <vector>
#include <cmath>
#include <H5Cpp.h>
#include <iostream>
#include <fftw3.h>

#include "containers_3d.h"

const double PI = 3.14159265358979323846;

int nx, ny, nz;
double lx, ly, lz;
double dx, dy, dz;

fftw_plan fftw_forward;
fftw_plan fftw_backward;
double * fftw_in;
fftw_complex * fftw_out;


void deposit(double y, double z, double value, array3d<double> & array, int slice, double yshift=0.0, double zshift=0.0) {
    int j1 = (int) (y / dy - yshift);
    int j2 = (j1 + 1) % ny;
    double y_frac = (y / dy - yshift) - j1;
    if (j1 < 0) {
        j1 += ny;
    }

    int k1 = (int) (z / dz - zshift);
    int k2 = (k1 + 1) % nz;
    double z_frac = (z / dz - zshift) - k1;
    if (k1 < 0) {
        k1 += nz;
    }

    array(slice, j1, k1) += value * (1 - y_frac) * (1 - z_frac);
    array(slice, j2, k1) += value * y_frac * (1 - z_frac);
    array(slice, j1, k2) += value * (1 - y_frac) * z_frac;
    array(slice, j2, k2) += value * y_frac * z_frac;
}

void deposit(double y, double z, double value, double * array) {
    int j1 = (int) (y / dy);
    int j2 = (j1 + 1) % ny;
    double y_frac = (y / dy) - j1;

    int k1 = (int) (z / dz);
    int k2 = (k1 + 1) % nz;
    double z_frac = (z / dz) - k1;

    array[nz * j1 + k1] += value * (1 - y_frac) * (1 - z_frac);
    array[nz * j2 + k1] += value * y_frac * (1 - z_frac);
    array[nz * j1 + k2] += value * (1 - y_frac) * z_frac;
    array[nz * j2 + k2] += value * y_frac * z_frac;
}

double array_to_particle(double y, double z, array3d<double> & array, int slice, double yshift=0.0, double zshift=0.0) {
    int j1 = (int) (y / dy - yshift);
    int j2 = (j1 + 1) % ny;
    double y_frac = (y / dy - yshift) - j1;
    if (j1 < 0) {
        j1 += ny;
    }

    int k1 = (int) (z / dz - zshift);
    int k2 = (k1 + 1) % nz;
    double z_frac = (z / dz - zshift) - k1;
    if (k1 < 0) {
        k1 += nz;
    }

    return array(slice, j1, k1) * (1 - y_frac) * (1 - z_frac) + array(slice, j2, k1) * y_frac * (1 - z_frac) +
           array(slice, j1, k2) * (1 - y_frac) * z_frac + array(slice, j2, k2) * y_frac * z_frac;
}

double array_to_particle(double y, double z, array2d<double> & array) {
    int j1 = (int) (y / dy);
    int j2 = (j1 + 1) % ny;
    double y_frac = (y / dy) - j1;

    int k1 = (int) (z / dz);
    int k2 = (k1 + 1) % nz;
    double z_frac = (z / dz) - k1;

    return array(j1, k1) * (1 - y_frac) * (1 - z_frac) + array(j2, k1) * y_frac * (1 - z_frac) +
           array(j1, k2) * (1 - y_frac) * z_frac + array(j2, k2) * y_frac * z_frac;
}


double array_yder_to_particle(double y, double z, array3d<double> & array, int slice) {
    int j1 = (int) (y / dy - 0.5);
    double y_frac = (y / dy - 0.5) - j1;
    if (j1 < 0) {
        j1 += ny;
    }
    int j2 = (j1 + 1) % ny;
    int j3 = (j1 + 2) % ny;

    int k1 = (int) (z / dz);
    int k2 = (k1 + 1) % nz;
    double z_frac = (z / dz) - k1;

    double der_y1_z1 = (array(slice, j2, k1) - array(slice, j1, k1)) / dy;
    double der_y2_z1 = (array(slice, j3, k1) - array(slice, j2, k1)) / dy;
    double der_y1_z2 = (array(slice, j2, k2) - array(slice, j1, k2)) / dy;
    double der_y2_z2 = (array(slice, j3, k2) - array(slice, j2, k2)) / dy;

    return der_y1_z1 * (1 - y_frac) * (1 - z_frac) + der_y2_z1 * y_frac * (1 - z_frac) +
           der_y1_z2 * (1 - y_frac) * z_frac + der_y2_z2 * y_frac * z_frac;
}

double array_zder_to_particle(double y, double z, array3d<double> & array, int slice) {
    int j1 = (int) (y / dy);
    int j2 = (j1 + 1) % ny;
    double y_frac = (y / dy) - j1;

    int k1 = (int) (z / dz - 0.5);
    double z_frac = (z / dz - 0.5) - k1;
    if (k1 < 0) {
        k1 += nz;
    }
    int k2 = (k1 + 1) % nz;
    int k3 = (k1 + 2) % nz;

    double der_y1_z1 = (array(slice, j1, k2) - array(slice, j1, k1)) / dz;
    double der_y2_z1 = (array(slice, j2, k2) - array(slice, j2, k1)) / dz;
    double der_y1_z2 = (array(slice, j1, k3) - array(slice, j1, k2)) / dz;
    double der_y2_z2 = (array(slice, j2, k3) - array(slice, j2, k2)) / dz;

    return der_y1_z1 * (1 - y_frac) * (1 - z_frac) + der_y2_z1 * y_frac * (1 - z_frac) +
           der_y1_z2 * (1 - y_frac) * z_frac + der_y2_z2 * y_frac * z_frac;
}

void normalize_coordinates(double & y, double & z) {
    if (y < 0) {
        y += ly;
    } else if (y >= ly) {
        y -= ly;
    }

    if (z < 0) {
        z += lz;
    } else if (z >= lz) {
        z -= lz;
    }
}

void solve_poisson_equation(double D=0.0) {
    fftw_execute(fftw_forward);
    fftw_out[0][0] = 0;
    fftw_out[0][1] = 0;
    const int nz_size = (nz / 2) + 1;

    for (int j = 0; j < ny; j++) {
        for (int k = 0; k < nz_size; k++){
            if ((k == 0) && (j == 0)) {
                continue;
            }

            double multiplier = -D + (2 * cos(2 * PI * j / ny) - 2) / dy / dy + (2 * cos(2 * PI * k / nz) - 2) / dz / dz;

            fftw_out[nz_size * j + k][0] /= multiplier;
            fftw_out[nz_size * j + k][1] /= multiplier;
        }
    }
    fftw_execute(fftw_backward);
}

void write_array(array3d<double> & field, std::string name, H5::H5File file) {
    const hsize_t dims[3] {nx, ny, nz};
    H5::DataSpace dataspace(3, dims);

    H5::DataSet dataset = file.createDataSet(name, H5::PredType::NATIVE_DOUBLE, dataspace);
    dataset.write(&(field(0,0,0)), H5::PredType::NATIVE_DOUBLE);
}

double rhobunch(double xi, double y, double z) {
    const double x0 = 4;
    const double y0 = 10;
    const double z0 = 10;

    double xwidth = 2;
    double ywidth = 1.5;
    double zwidth = 1.5;

    double x_prof = (fabs(xi - x0) < xwidth) ? cos(0.5 * PI * (xi - x0) / xwidth) : 0.0;
    x_prof *= x_prof;
    double y_prof = (fabs(y - y0) < ywidth) ? cos(0.5 * PI * (y - y0) / ywidth) : 0.0;
    y_prof *= y_prof;
    double z_prof = (fabs(z - z0) < zwidth) ? cos(0.5 * PI * (z - z0) / zwidth) : 0.0;
    z_prof *= z_prof;

    return -3.0 * x_prof * y_prof * z_prof;
}

int main() {

    nx = 300;
    ny = 300;
    nz = 300;

    fftw_in = fftw_alloc_real(ny * nz);
    fftw_out = fftw_alloc_complex(ny * (nz / 2 + 1));
    fftw_forward = fftw_plan_dft_r2c_2d(ny, nz, fftw_in, fftw_out, FFTW_MEASURE);
    fftw_backward = fftw_plan_dft_c2r_2d(ny, nz, fftw_out, fftw_in, FFTW_MEASURE);

    lx = 20;
    ly = 20;
    lz = 20;

    dx = lx / nx;
    dy = ly / ny;
    dz = lz / nz;

    int ppcy = 2;
    int ppcz = 2;

    const int iterations = 4;

    std::vector<particle> particles(ny * nz * ppcy * ppcz);
    array2d<double> psi_middle(ny, nz);
    array2d<double> djy_dxi(ny, nz);
    array2d<double> djz_dxi(ny, nz);

    array3d<double> psi(nx, ny, nz);
    array3d<double> psi_source(nx, ny, nz);
    array3d<double> dpsi_dy(nx, ny, nz);
    array3d<double> a(nx, ny, nz);
    array3d<double> jx(nx, ny, nz);
    array3d<double> jy(nx, ny, nz);
    array3d<double> jz(nx, ny, nz);
    array3d<double> rho(nx, ny, nz);
    array3d<double> ex(nx, ny, nz);
    array3d<double> ey(nx, ny, nz);
    array3d<double> ez(nx, ny, nz);
    array3d<double> by(nx, ny, nz);
    array3d<double> bz(nx, ny, nz);

    for (int i = 0; i < ppcy * ny; i++) {
        for (int j = 0; j < ppcz * nz; j++) {
            int index = ppcz * nz * i + j;
            particles[index].y = (i + 0.5) * dy / ppcy;
            particles[index].z = (j + 0.5) * dz / ppcz;
            particles[index].n = -1.0 / ppcy / ppcz;
        }
    }

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                double x = i * dx;
                double y = j * dy;
                double z = k * dz;
                const double x0 = 5;
                const double y0 = ly / 2;
                const double z0 = lz / 2;
                const double xsigma = 2;
                const double ysigma = 2;
                const double zsigma = 2;
                a(i, j, k) = 0.0 * exp(- (x-x0) * (x-x0) / xsigma / xsigma
                                       - (y-y0) * (y-y0) / ysigma / ysigma
                                       - (z-z0) * (z-z0) / zsigma / zsigma);
            }

        }
    }

    // main loop
    for (int i = 0; i < nx; i++) {
        std::cout << "Iteration " << i << std::endl;
        // psi_source deposition
        for (auto & p : particles) {
            deposit(p.y, p.z, -p.n, psi_source, i);
        }

        // cacluate psi
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                fftw_in[nz * j + k] = psi_source(i, j, k) - 1.0;
                //fftw_in[nz * j + k] = exp(- (k - 0.5 * nz) * (k - 0.5 * nz) / 20.0 / 20.0) *
                //        exp(- (j - 0.5 * ny) * (j - 0.5 * ny) / 20.0 / 20.0);
            }
        }
        solve_poisson_equation();
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                psi(i, j, k) = (fftw_in[nz * j + k]) / ny / nz;
            }
        }

        // calculate initial gamma, px
        for (auto & p : particles) {
            double a_particle = array_to_particle(p.y, p.z, a, i);
            double psi_particle = array_to_particle(p.y, p.z, psi, i);
            p.gamma = 0.5 * (1 + p.py * p.py + p.pz * p.pz + a_particle + (1 + psi_particle) * (1 + psi_particle))
                    / (1 + psi_particle);
            p.px = 0.5 * (1 + p.py * p.py + p.pz * p.pz + a_particle - (1 + psi_particle) * (1 + psi_particle))
                    / (1 + psi_particle);
        }

        // initial jx, rho deposition
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                jx(i, j, k) = rhobunch(i * dx, j * dy, k * dz);
                rho(i, j, k) = rhobunch(i * dx, j * dy, k * dz);
            }
        }

        for (auto & p : particles) {
            double vx = p.px / p.gamma;
            deposit(p.y, p.z, p.n * vx / (1 - vx), jx, i);
            deposit(p.y, p.z, p.n / (1 - vx), rho, i);
        }

        if (i > 0) {
            for (int j = 0; j < ny; j++) {
                for (int k = 0; k < nz; k++) {
                    bz(i, j, k) = bz(i-1, j, k);
                }
            }
        }

        if (i == nx - 1) {
            break;
        }

        for (int n = 0; n < iterations; n++) {
            // advance momenta
            for (auto & p : particles) {
                double psi_particle = array_to_particle(p.y, p.z, psi, i);
                double da_dy_particle = array_yder_to_particle(p.y, p.z, a, i);
                double da_dz_particle = array_zder_to_particle(p.y, p.z, a, i);
                double dpsi_dy_particle = array_yder_to_particle(p.y, p.z, psi, i);
                double dpsi_dz_particle = array_zder_to_particle(p.y, p.z, psi, i);
                double by_particle = array_to_particle(p.y, p.z, by, i, 0.5, 0.0);
                double bz_particle = array_to_particle(p.y, p.z, bz, i, 0.0, 0.5);

                p.py_next = p.py - dx * 0.5 * da_dy_particle / (1 + psi_particle);
                p.py_next += dx * p.gamma * dpsi_dy_particle / (1 + psi_particle);
                p.py_next -= dx * bz_particle;
                p.pz_next = p.pz - dx * 0.5 * da_dz_particle / (1 + psi_particle);
                p.pz_next += dx * p.gamma * dpsi_dz_particle / (1 + psi_particle);
                p.pz_next += dx * by_particle;
            }

            // advance half coordinate
            for (auto & p : particles) {
                double psi_particle = array_to_particle(p.y, p.z, psi, i);
                p.y_middle = p.y + 0.5 * dx * p.py_next / (1 + psi_particle);
                p.z_middle = p.z + 0.5 * dx * p.pz_next / (1 + psi_particle);
                normalize_coordinates(p.y_middle, p.z_middle);
            }

            // psi_source middle deposition
            for (int j = 0; j < ny; j++) {
                for (int k = 0; k < nz; k++) {
                    fftw_in[nz * j + k] = -1.0; // rho_ion
                }
            }

            for (auto & p : particles) {
                deposit(p.y, p.z, -p.n, fftw_in);
            }

            // calculate psi_middle
            solve_poisson_equation();
            for (int j = 0; j < ny; j++) {
                for (int k = 0; k < ny; k++) {
                    psi_middle(j, k) = (fftw_in[nz * j + k]) / ny / nz;
                }
            }

            // deposit jy, jz
            for (int j = 0; j < ny; j++) {
                for (int k = 0; k < nz; k++){
                    jy(i+1, j, k) = 0.0;
                    jz(i+1, j, k) = 0.0;
                }
            }

            for (auto & p : particles) {
                int j1 = (int) (p.y / dy);
                int j2 = (j1 + 1) % ny;
                double frac = (p.y / dy) - j1;
                double psi_particle = array_to_particle(p.y, p.z, psi_middle);

                deposit(p.y, p.z, p.n * p.py_next / (1 + psi_particle), jy, i+1, 0.5, 0.0);
                deposit(p.y, p.z, p.n * p.pz_next / (1 + psi_particle), jz, i+1, 0.0, 0.5);
            }

            // calculate new djy_dxi, djz_dxi
            for (int j = 0; j < ny; j++) {
                for (int k = 0; k < nz; k++) {
                    djy_dxi(j, k) = (jy(i, j, k) - jy(i+1, j, k)) / dx;
                    djz_dxi(j, k) = (jz(i, j, k) - jz(i+1, j, k)) / dx;
                }
            }

            // calculate new gamma, px
            for (auto & p : particles) {
                double a_particle = array_to_particle(p.y, p.z, a, i);
                double psi_particle = array_to_particle(p.y, p.z, psi, i);
                double p_squared = 0.25 * (p.py + p.py_next) * (p.py + p.py_next) + 0.25 * (p.pz + p.pz_next) * (p.pz + p.pz_next);
                p.gamma = 0.5 * (1 + p_squared + a_particle +
                        (1 + psi_particle) * (1 + psi_particle)) / (1 + psi_particle);
                p.px = 0.5 * (1 + p_squared + a_particle - (1 + psi_particle) * (1 + psi_particle)) / (1 + psi_particle);
            }

            // new jx, rho deposition
            for (int j = 0; j < ny; j++) {
                for (int k = 0; k < nz; k++) {
                    jx(i, j, k) = rhobunch(i * dx, j * dy, k * dz);
                    rho(i, j, k) = rhobunch(i * dx, j * dy, k * dz);
                }
            }

            for (auto & p : particles) {
                double vx = p.px / p.gamma;
                deposit(p.y, p.z, p.n * vx / (1 - vx), jx, i);
                deposit(p.y, p.z, p.n / (1 - vx), rho, i);
            }

            // new source for B_y
            for (int j = 0; j < ny; j++) {
                for (int k = 0; k < nz-1; k++) {
                    fftw_in[nz * j + k] = -by(i, j, k) - (jx(i, j, k+1) - jx(i, j, k)) / dz + djz_dxi(j, k);
                }
                fftw_in[nz * j + (nz-1)] = -by(i, j, nz-1) - (jx(i, j, 0) - jx(i, j, nz-1)) / dz + djz_dxi(j, nz-1);
            }

            // new guess for B_y
            solve_poisson_equation(1.0);

            for (int j = 0; j < ny; j++) {
                for (int k = 0; k < nz; k++) {
                    by(i, j, k) = fftw_in[nz * j + k] / ny / nz;
                }
            }

            // new source for B_z
            for (int j = 0; j < ny-1; j++) {
                for (int k = 0; k < nz; k++) {
                    fftw_in[nz * j + k] = -bz(i, j, k) + (jx(i, j+1, k) - jx(i, j, k)) / dy - djy_dxi(j, k);
                }
            }
            for (int k = 0; k < nz; k++) {
                fftw_in[nz * (ny-1) + k] = -bz(i, 0, k) + (jx(i, 0, k) - jx(i, ny-1, k)) / dy - djy_dxi(ny-1, k);
            }

            // new guess for B_z
            solve_poisson_equation(1.0);

            for (int j = 0; j < ny; j++) {
                for (int k = 0; k < nz; k++) {
                    bz(i, j, k) = fftw_in[nz * j + k] / ny / nz;
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
            p.y += dx * p.py / (1 + psi_particle);
            p.z += dx * p.pz / (1 + psi_particle);
            normalize_coordinates(p.y, p.z);
        }
    }

    // calculate ex
    for (int i = 1; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                ex(i, j, k) = (psi(i, j, k) - psi(i-1, j, k)) / dx;
            }
        }
    }

    // calculate ey
    for (int i = 1; i < nx; i++) {
        for (int j = 0; j < ny-1; j++) {
            for (int k = 0; k < nz; k++) {
                ey(i, j, k) = -(psi(i, j+1, k) - psi(i, j, k)) / dy + bz(i, j, k);
            }
        }
        for (int k = 0; k < nz; k++) {
            ey(i, ny-1, k) = -(psi(i, 0, k) - psi(i, ny-1, k)) / dy + bz(i, ny-1, k);
        }
    }

    // calculate ez
    for (int i = 1; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz-1; k++) {
                ez(i, j, k) = -(psi(i, j, k+1) - psi(i, j, k)) / dz - by(i, j, k);
            }
            ez(i, j, nz-1) = -(psi(i, j, 0) - psi(i, j, nz-1)) / dz - by(i, j, nz-1);
        }
    }

    H5::H5File fields_file("Fields.h5", H5F_ACC_TRUNC);

    write_array(a, "aSqr", fields_file);
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

    fftw_destroy_plan(fftw_forward);
    fftw_destroy_plan(fftw_backward);

    fftw_free(fftw_in);
    fftw_free(fftw_out);
}

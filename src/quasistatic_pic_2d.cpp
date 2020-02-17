#include <vector>
#include <cmath>
#include <H5Cpp.h>
#include <iostream>
#include <fftw3.h>

#include "containers_2d.h"

const double PI = 3.14159265358979323846;

int main() {

    int ny = 400;
    int nx = 400;

    double * fftw_in = fftw_alloc_real(ny);
    fftw_complex * fftw_out = fftw_alloc_complex(ny / 2 + 1);
    fftw_plan fftw_forward = fftw_plan_dft_r2c_1d(ny, fftw_in, fftw_out, FFTW_MEASURE);
    fftw_plan fftw_backward = fftw_plan_dft_c2r_1d(ny, fftw_out, fftw_in, FFTW_MEASURE);

    double ly = 20;
    double lx = 40;

    double dx = lx / nx;
    double dy = ly / ny;

    int ppc = 1;

    const int iterations = 100;

    std::vector<particle> particles(ny * ppc);
    std::vector<double> psi_middle(ny);
    std::vector<double> djy_dxi(ny);

    array2d<double> djy_dxi_arr(nx, ny);
    array2d<double> djx_dy(nx, ny);

    array2d<double> psi(nx, ny);
    array2d<double> psi_source(nx, ny);
    array2d<double> dpsi_dy(nx, ny);
    array2d<double> a(nx, ny);
    array2d<double> jx(nx, ny);
    array2d<double> jy(nx, ny);
    array2d<double> rho(nx, ny);
    array2d<double> ex(nx, ny);
    array2d<double> ey(nx, ny);
    array2d<double> bz(nx, ny);

    for (int i = 0; i < particles.size(); i++) {
        particles[i].y = (i + 0.5) * dy / ppc;
        particles[i].n = -1.0 / ppc;
    }

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            double x = i * dx;
            double y = j * dy;
            const double x0 = 5;
            const double y0 = ly / 2;
            const double xsigma = 2;
            const double ysigma = 2;
            a(i, j) = 1.3 * exp(- (x-x0) * (x-x0) / xsigma / xsigma - (y-y0)*(y-y0) / ysigma / ysigma);
        }
    }

    // main loop
    for (int i = 0; i < nx; i++) {
        std::cout << "Iteration " << i << std::endl;
        // psi_source deposition
        for (auto & p : particles) {
            int j1 = (int) (p.y / dy);
            int j2 = (j1 + 1) % ny;
            double y_frac = (p.y / dy) - j1;
            psi_source(i, j1) -= (1 - y_frac) * p.n;
            psi_source(i, j2) -= y_frac * p.n;
        }

        // cacluate psi
        for (int j = 0; j < ny; j++) {
            fftw_in[j] = psi_source(i, j) - 1.0;
        }
        fftw_execute(fftw_forward);
        fftw_out[0][0] = 0;
        fftw_out[0][1] = 0;
        for (int j = 1; j < (ny / 2) + 1; j++) {
            double k = 2 * PI * j / ly;
            fftw_out[j][0] /= - k * k;
            fftw_out[j][1] /= - k * k;
        }
        fftw_execute(fftw_backward);
        for (int j = 0; j < ny; j++) {
            psi(i, j) = (fftw_in[j] - fftw_in[0]) / ny;
        }

        // calculate initial gamma, px
        for (auto & p : particles) {
            int j1 = (int) (p.y / dy);
            int j2 = (j1 + 1) % ny;
            double y_frac = (p.y / dy) - j1;
            double a_particle = a(i, j1) * (1 - y_frac) + a(i, j2) * y_frac;
            double psi_particle = (1 - y_frac) * psi(i, j1) + y_frac * psi(i, j2);
            p.gamma = 0.5 * (1 + p.py * p.py + a_particle + (1 + psi_particle) * (1 + psi_particle)) / (1 + psi_particle);
            p.px = 0.5 * (1 + p.py * p.py + a_particle - (1 + psi_particle) * (1 + psi_particle)) / (1 + psi_particle);
        }

        // initial jx, rho deposition
        for (auto & p : particles) {
            int j1 = (int) (p.y / dy);
            int j2 = (j1 + 1) % ny;
            double y_frac = (p.y / dy) - j1;
            double vx = p.px / p.gamma;

            jx(i, j1) += (1 - y_frac) * p.n * vx / (1 - vx);
            jx(i, j2) += y_frac * p.n * vx / (1 - vx);

            rho(i, j1) += (1 - y_frac) * p.n / (1 - vx);
            rho(i, j2) += y_frac * p.n / (1 - vx);
        }

        /*
        // initial source for B_z
        for (int j = 0; j < ny-1; j++) {
            fftw_in[j] = (jx(i, j+1) - jx(i, j)) / dy - djy_dxi[j];
        }
        fftw_in[ny-1] = (jx(i, 0) - jx(i, ny-1)) / dy - djy_dxi[ny-1];

        // initial guess for B_z
        fftw_execute(fftw_forward);
        fftw_out[0][0] = 0;
        fftw_out[0][1] = 0;
        for (int j = 1; j < (ny / 2) + 1; j++) {
            double k = 2 * PI * j / ly;
            fftw_out[j][0] /= - k * k;
            fftw_out[j][1] /= - k * k;
        }
        fftw_execute(fftw_backward);
        for (int j = 0; j < ny; j++) {
            bz(i, j) = (fftw_in[j] - fftw_in[0]) / ny;
        }
        */
        if (i > 0) {
            for (int j = 0; j < ny; j++) {
                bz(i, j) = bz(i-1, j);
            }
        }

        if (i == nx - 1) {
            break;
        }

        for (int n = 0; n < iterations; n++) {
            // advance momenta
            for (auto & p : particles) {
                int j1 = (int) (p.y / dy);
                int j2 = (j1 + 1) % ny;
                double frac = (p.y / dy) - j1;
                double psi_particle = (1 - frac) * psi(i, j1) + frac * psi(i, j2);

                j1 = (int) (p.y / dy - 0.5);
                frac = (p.y / dy - 0.5) - j1;
                if (j1 < 0) {
                    j1 += ny;
                }
                j2 = (j1 + 1) % ny;
                int j3 = (j1 + 2) % ny;
                double da_dy_particle = (1 - frac) * (a(i, j2) - a(i, j1)) / dy + frac * (a(i, j3) - a(i, j2)) / dy;
                double dpsi_dy_particle = (1 - frac) * (psi(i, j2) - psi(i, j1)) / dy + frac * (psi(i, j3) - psi(i, j2)) / dy;
                double bz_particle = (1 - frac) * bz(i, j1) + frac * bz(i, j2);
                p.py_next = p.py - dx * 0.5 * da_dy_particle / (1 + psi_particle);
                p.py_next += dx * p.gamma * dpsi_dy_particle / (1 + psi_particle);
                p.py_next -= dx * bz_particle;
            }

            // advance half coordinate
            for (auto & p : particles) {
                int j1 = (int) (p.y / dy);
                int j2 = (j1 + 1) % ny;
                double frac = (p.y / dy) - j1;
                double psi_particle = (1 - frac) * psi(i, j1) + frac * psi(i, j2);
                p.y_middle = p.y + 0.5 * dx * p.py_next / (1 + psi_particle);
                if (p.y_middle < 0) {
                    p.y_middle += ly;
                } else if (p.y_middle >= ly) {
                    p.y_middle -= ly;
                }
            }

            // psi_source middle deposition
            for (int j = 0; j < ny; j++) {
                fftw_in[j] = -1.0; // rho_ion
            }

            for (auto & p : particles) {
                int j1 = (int) (p.y_middle / dy);
                int j2 = (j1 + 1) % ny;
                double y_frac = (p.y_middle / dy) - j1;
                fftw_in[j1] -= (1 - y_frac) * p.n;
                fftw_in[j2] -= y_frac * p.n;
            }

            // calculate psi_middle
            fftw_execute(fftw_forward);
            fftw_out[0][0] = 0;
            fftw_out[0][1] = 0;
            for (int j = 1; j < (ny / 2) + 1; j++) {
                double k = 2 * PI * j / ly;
                fftw_out[j][0] /= - k * k;
                fftw_out[j][1] /= - k * k;
            }
            fftw_execute(fftw_backward);
            for (int j = 0; j < ny; j++) {
                psi_middle[j] = (fftw_in[j] - fftw_in[0]) / ny;
            }

            // deposit jy
            for (int j = 0; j < ny; j++) {
                jy(i+1, j) = 0.0;
            }

            for (auto & p : particles) {
                int j1 = (int) (p.y / dy);
                int j2 = (j1 + 1) % ny;
                double frac = (p.y / dy) - j1;
                double psi_particle = (1 - frac) * psi_middle[j1] + frac * psi_middle[j2];

                j1 = (int) (p.y / dy - 0.5);
                j2 = (j1 + 1) % ny;
                double y_frac = (p.y / dy - 0.5) - j1;

                jy(i+1, j1) += (1 - y_frac) * p.n * p.py_next / (1 + psi_particle);
                jy(i+1, j2) += y_frac * p.n * p.py_next / (1 + psi_particle);
            }

            // calculate new djy_dxi
            for (int j = 0; j < ny; j++) {
                djy_dxi[j] = (jy(i, j) - jy(i+1, j)) / dx;
                djy_dxi_arr(i, j) = djy_dxi[j];
            }

            // calculate new gamma, px
            for (auto & p : particles) {
                int j1 = (int) (p.y / dy);
                int j2 = (j1 + 1) % ny;
                double y_frac = (p.y / dy) - j1;
                double a_particle = a(i, j1) * (1 - y_frac) + a(i, j2) * y_frac;
                double psi_particle = (1 - y_frac) * psi(i, j1) + y_frac * psi(i, j2);
                p.gamma = 0.5 * (1 + 0.25 * (p.py + p.py_next) * (p.py + p.py_next) + a_particle +
                        (1 + psi_particle) * (1 + psi_particle)) / (1 + psi_particle);
                p.px = 0.5 * (1 + 0.25 * (p.py + p.py_next) * (p.py + p.py_next) + a_particle - (1 + psi_particle) * (1 + psi_particle)) / (1 + psi_particle);
            }

            // new jx, rho deposition
            for (int j = 0; j < ny; j++) {
                jx(i, j) = 0.0;
                rho(i, j) = 0.0;
            }

            for (auto & p : particles) {
                int j1 = (int) (p.y / dy);
                int j2 = (j1 + 1) % ny;
                double y_frac = (p.y / dy) - j1;
                double vx = p.px / p.gamma;

                jx(i, j1) += (1 - y_frac) * p.n * vx / (1 - vx);
                jx(i, j2) += y_frac * p.n * vx / (1 - vx);

                rho(i, j1) += (1 - y_frac) * p.n / (1 - vx);
                rho(i, j2) += y_frac * p.n / (1 - vx);
            }

            // new source for B_z
            for (int j = 0; j < ny-1; j++) {
                fftw_in[j] = -bz(i, j) + (jx(i, j+1) - jx(i, j)) / dy - djy_dxi[j];
                djx_dy(i, j) = (jx(i, j+1) - jx(i, j)) / dy;
            }
            fftw_in[ny-1] = -bz(i,0) + (jx(i, 0) - jx(i, ny-1)) / dy - djy_dxi[ny-1];
            djx_dy(i, ny-1) = (jx(i, 0) - jx(i, ny-1)) / dy;

            // new guess for B_z
            fftw_execute(fftw_forward);
            for (int j = 0; j < (ny / 2) + 1; j++) {
                double k = 2 * PI * j / ly;
                fftw_out[j][0] /= (-1 - k * k);
                fftw_out[j][1] /= (-1 - k * k);
            }
            fftw_execute(fftw_backward);
            for (int j = 0; j < ny; j++) {
                bz(i, j) = (fftw_in[j] - fftw_in[0]) / ny;
            }
        }

        for (auto & p : particles) {
            p.py = p.py_next;
        }

        // advance coordinate
        for (auto & p : particles) {
            int j1 = (int) (p.y / dy);
            int j2 = (j1 + 1) % ny;
            double frac = (p.y / dy) - j1;
            double psi_particle = (1 - frac) * psi_middle[j1] + frac * psi_middle[j2];
            p.y += dx * p.py / (1 + psi_particle);
            if (p.y < 0) {
                p.y += ly;
            } else if (p.y >= ly) {
                p.y -= ly;
            }
        }
    }

    // calculate ex
    for (int i = 1; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            ex(i, j) = (psi(i, j) - psi(i-1, j)) / dx;
        }
    }

    // calculate ey
    for (int i = 1; i < nx; i++) {
        for (int j = 0; j < ny-1; j++) {
            ey(i, j) = -(psi(i, j+1) - psi(i, j)) / dy + bz(i, j);
        }
        ey(i, ny-1) = -(psi(i, 0) - psi(i, ny-1)) / dy + bz(i, ny-1);
    }

    H5::H5File fields_file("Fields.h5", H5F_ACC_TRUNC);
    const hsize_t dims[2] {nx, ny};
    H5::DataSpace dataspace(2, dims);

    H5::DataSet dataset = fields_file.createDataSet("aSqr", H5::PredType::NATIVE_DOUBLE, dataspace);
    dataset.write(&(a(0,0)), H5::PredType::NATIVE_DOUBLE);

    dataset = fields_file.createDataSet("jx_minus_rho", H5::PredType::NATIVE_DOUBLE, dataspace);
    dataset.write(&(psi_source(0,0)), H5::PredType::NATIVE_DOUBLE);

    dataset = fields_file.createDataSet("jx", H5::PredType::NATIVE_DOUBLE, dataspace);
    dataset.write(&(jx(0,0)), H5::PredType::NATIVE_DOUBLE);

    dataset = fields_file.createDataSet("rho", H5::PredType::NATIVE_DOUBLE, dataspace);
    dataset.write(&(rho(0,0)), H5::PredType::NATIVE_DOUBLE);

    dataset = fields_file.createDataSet("jy", H5::PredType::NATIVE_DOUBLE, dataspace);
    dataset.write(&(jy(0,0)), H5::PredType::NATIVE_DOUBLE);

    dataset = fields_file.createDataSet("psi", H5::PredType::NATIVE_DOUBLE, dataspace);
    dataset.write(&(psi(0,0)), H5::PredType::NATIVE_DOUBLE);

    dataset = fields_file.createDataSet("ex", H5::PredType::NATIVE_DOUBLE, dataspace);
    dataset.write(&(ex(0,0)), H5::PredType::NATIVE_DOUBLE);

    dataset = fields_file.createDataSet("ey", H5::PredType::NATIVE_DOUBLE, dataspace);
    dataset.write(&(ey(0,0)), H5::PredType::NATIVE_DOUBLE);

    dataset = fields_file.createDataSet("bz", H5::PredType::NATIVE_DOUBLE, dataspace);
    dataset.write(&(bz(0,0)), H5::PredType::NATIVE_DOUBLE);

    dataset = fields_file.createDataSet("djx_dy", H5::PredType::NATIVE_DOUBLE, dataspace);
    dataset.write(&(djx_dy(0,0)), H5::PredType::NATIVE_DOUBLE);

    dataset = fields_file.createDataSet("djy_dxi", H5::PredType::NATIVE_DOUBLE, dataspace);
    dataset.write(&(djy_dxi_arr(0,0)), H5::PredType::NATIVE_DOUBLE);

    fftw_destroy_plan(fftw_forward);
    fftw_destroy_plan(fftw_backward);

    fftw_free(fftw_in);
    fftw_free(fftw_out);
}

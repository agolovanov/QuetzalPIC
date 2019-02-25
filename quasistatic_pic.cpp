#include <vector>
#include <cmath>
#include <H5Cpp.h>
#include <iostream>
#include <fftw3.h>

#include "containers.h"

const double PI = 3.14159265358979323846;

int main() {

    int ny = 800;
    int nx = 800;

    double * fftw_in = fftw_alloc_real(ny);
    fftw_complex * fftw_out = fftw_alloc_complex(ny / 2 + 1);
    fftw_plan fftw_forward = fftw_plan_dft_r2c_1d(ny, fftw_in, fftw_out, FFTW_MEASURE);
    fftw_plan fftw_backward = fftw_plan_dft_c2r_1d(ny, fftw_out, fftw_in, FFTW_MEASURE);

    double ly = 30;
    double lx = 20;

    double dx = lx / nx;
    double dy = ly / ny;

    int ppc = 4;

    std::vector<particle> particles(ny * ppc);

    array2d<double> psi(nx, ny);
    array2d<double> a(nx, ny);
    array2d<double> psi_source(nx, ny);
    array2d<double> jx(nx, ny);
    array2d<double> rho(nx, ny);
    array2d<double> ex(nx, ny);
    array2d<double> ey(nx, ny);

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
            a(i, j) = 3 * exp(- (x-x0) * (x-x0) / xsigma / xsigma - (y-y0)*(y-y0) / ysigma / ysigma);
        }
    }

    // main loop
    for (int i = 0; i < nx; i++) {
        // calculate gamma, px
        for (auto & p : particles) {
            int j1 = (int) (p.y / dy);
            int j2 = (j1 + 1) % ny;
            double y_frac = (p.y / dy) - j1;
            double a_particle = a(i, j1) * (1 - y_frac) + a(i, j2) * y_frac;
            double psi_particle = (1 - y_frac) * psi(i, j1) + y_frac * psi(i, j2);
            p.gamma = 0.5 * (1 + p.py * p.py + a_particle + (1 + psi_particle) * (1 + psi_particle)) / (1 + psi_particle);
            p.px = 0.5 * (1 + p.py * p.py + a_particle - (1 + psi_particle) * (1 + psi_particle)) / (1 + psi_particle);
        }

        // psi_source, jx, rho deposition
        for (auto & p : particles) {
            int j1 = (int) (p.y / dy);
            int j2 = (j1 + 1) % ny;
            double y_frac = (p.y / dy) - j1;
            double vx = p.px / p.gamma;
            psi_source(i, j1) -= (1 - y_frac) * p.n;
            psi_source(i, j2) -= y_frac * p.n;

            jx(i, j1) += (1 - y_frac) * p.n * vx / (1 - vx);
            jx(i, j2) += y_frac * p.n * vx / (1 - vx);

            rho(i, j1) += (1 - y_frac) * p.n / (1 - vx);
            rho(i, j2) += y_frac * p.n / (1 - vx);
        }

        // calculate psi
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
            psi(i, j) = fftw_in[j] / ny;
        }

        // advance momenta
        for (auto & p : particles) {
            int j1 = (int) (p.y / dy);
            int j2 = (j1 + 1) % ny;
            double frac = (p.y / dy) - j1;
            double psi_particle = (1 - frac) * psi(i, j1) + frac * psi(i, j2);

            j1 = (int) (p.y / dy - 0.5);
            j2 = (j1 + 1) % ny;
            int j3 = (j1 + 2) % ny;
            frac = (p.y / dy - 0.5) - j1;
            double da_dy_particle = (1 - frac) * (a(i, j2) - a(i, j1)) / dy + frac * (a(i, j3) - a(i, j2)) / dy;
            double dpsi_dy_particle = (1 - frac) * (psi(i, j2) - psi(i, j1)) / dy + frac * (psi(i, j3) - psi(i, j2)) / dy;
            p.py -= dx * 0.5 * da_dy_particle / (1 + psi_particle);
            p.py += dx * p.gamma * dpsi_dy_particle / (1 + psi_particle);
        }

        // advance coordinates
        for (auto & p : particles) {
            int j1 = (int) (p.y / dy);
            int j2 = (j1 + 1) % ny;
            double frac = (p.y / dy) - j1;
            double psi_particle = (1 - frac) * psi(i, j1) + frac * psi(i, j2);
            p.y += dx * p.py / (1 + psi_particle);
            if (p.y < 0) {
                std::cout << i << "\n";
                std::cout << p.y << "\t";
                p.y += ly;
                std::cout << p.y << "\n";
            } else if (p.y >= ly) {
                std::cout << i << "\n";
                std::cout << p.y << "\t";
                p.y -= ly;
                std::cout << p.y << "\n";
            }
        }
    }

    // calculate ex
    for (int i = 1; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            ex(i, j) = (psi(i, j) - psi(i-1, j)) / dx;
        }
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

    dataset = fields_file.createDataSet("psi", H5::PredType::NATIVE_DOUBLE, dataspace);
    dataset.write(&(psi(0,0)), H5::PredType::NATIVE_DOUBLE);

    dataset = fields_file.createDataSet("ex", H5::PredType::NATIVE_DOUBLE, dataspace);
    dataset.write(&(ex(0,0)), H5::PredType::NATIVE_DOUBLE);

    fftw_destroy_plan(fftw_forward);
    fftw_destroy_plan(fftw_backward);

    fftw_free(fftw_in);
    fftw_free(fftw_out);
}

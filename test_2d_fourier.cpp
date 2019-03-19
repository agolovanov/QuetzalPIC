#include <fftw3.h>
#include <math.h>
#include <iostream>

fftw_plan fftw1;
fftw_plan fftw2;
double * fftw_in;
double * fftw_in_2;
fftw_complex * fftw_out;
fftw_complex * fftw_out_2;

const double PI = 3.14159265358979323846;

int main() {
    double ly=25, lz=10;
    int ny=20, nz=15;

    double dy = ly / ny;
    double dz = ly / nz;

    double ysigma = 0.06 * ly;
    double zsigma = 0.06 * lz;

    fftw_in = fftw_alloc_real(ny * nz);
    fftw_in_2 = fftw_alloc_real(ny * nz);
    fftw_out = fftw_alloc_complex(ny * ((nz / 2) + 1));
    fftw_out_2 = fftw_alloc_complex(ny * ((nz / 2) + 1));

    fftw1 = fftw_plan_dft_r2c_2d(ny, nz, fftw_in, fftw_out, FFTW_MEASURE);
    fftw2 = fftw_plan_dft_r2c_2d(ny, nz, fftw_in_2, fftw_out_2, FFTW_MEASURE);
    //fftw1 = fftw_plan_dft_2d(ny, nz, fftw_in, fftw_out, FFTW_FORWARD, FFTW_MEASURE);
    //fftw2 = fftw_plan_dft_2d(ny, nz, fftw_in_2, fftw_out_2, FFTW_FORWARD, FFTW_MEASURE);

    for (int j = 0; j < ny; j++) {
        for (int k = 0; k < nz; k++) {
            double yrel = dy * j - 0.5 * ly;
            double zrel = dz * k - 0.5 * lz;
            fftw_in[j * nz + k] = exp(-yrel * yrel / ysigma / ysigma - zrel * zrel / zsigma * zsigma);
        }
    }
    fftw_execute(fftw1);

    for (int j = 0; j < ny; j++) {
        for (int k = 0; k < nz; k++) {
            double yrel = dy * j - 0.5 * ly;
            double zrel = dz * k - 0.5 * lz;

            int j_prev = j > 0 ? j - 1 : ny - 1;
            int j_next = j < ny-1 ? j+1 : 0;

            int k_prev = k > 0 ? k - 1 : nz - 1;
            int k_next = k < nz-1 ? k+1 : 0;

            int index = j * nz + k;
            int index_left = j_prev * nz + k;
            int index_right = j_next * nz + k;
            int index_top = j * nz + k_prev;
            int index_bottom = j * nz + k_next;

            fftw_in_2[index] = (fftw_in[index_left] + fftw_in[index_right] - 2 * fftw_in[index]) / dy / dy +
                               (fftw_in[index_top] + fftw_in[index_bottom] - 2 * fftw_in[index]) / dz / dz;
        }
    }
    fftw_execute(fftw2);

    const int nz_size = (nz / 2) + 1;
    for (int j = 0; j < ny; j++) {
        for (int k = 0; k < nz_size; k++) {
            int index = nz_size * j + k;
            double y_multiplier = -(2 * cos(2 * PI * j / ny) - 2) / dy / dy;
            double z_multiplier = -(2 * cos(2 * PI * k / nz) - 2) / dz / dz;
            double multiplier = y_multiplier + z_multiplier;

            double re1 = fftw_out[index][0];
            double im1 = fftw_out[index][1];

            double re2 = fftw_out_2[index][0];
            double im2 = fftw_out_2[index][1];

            double abs1 = sqrt(re1 * re1 + im1 * im1);
            double abs2 = sqrt(re2 * re2 + im2 * im2);

            std::cout << j << ", " << k << ":\t" << re2 / re1 << "\t" << im2 / im1 << "\t expected: " << -multiplier <<  "\n";
        }
    }


    fftw_destroy_plan(fftw1);
    fftw_destroy_plan(fftw2);

    fftw_free(fftw_in);
    fftw_free(fftw_in_2);
    fftw_free(fftw_out);
    fftw_free(fftw_out_2);
}

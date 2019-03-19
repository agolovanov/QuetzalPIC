#include <fftw3.h>
#include <math.h>
#include <iostream>

fftw_plan fftw1;
fftw_plan fftw2;
//double * fftw_in;
fftw_complex * fftw_in;
fftw_complex * fftw_in_2;
fftw_complex * fftw_out;
fftw_complex * fftw_out_2;

const double PI = 3.14159265358979323846;

int main() {
    double l = 5;
    int n = 5;

    double dx = l / n;

    double sigma = 0.03 * l;

    //fftw_in = fftw_alloc_real(ny * nz);
    //fftw_out = fftw_alloc_complex(ny * (nz / 2 + 1));
    //fftw_out_2 = fftw_alloc_complex(ny * (nz / 2 + 1));

    fftw_in = fftw_alloc_complex(n);
    fftw_in_2 = fftw_alloc_complex(n);
    fftw_out = fftw_alloc_complex(n);
    fftw_out_2 = fftw_alloc_complex(n);

    //fftw1 = fftw_plan_dft_r2c_2d(ny, nz, fftw_in, fftw_out, FFTW_MEASURE);
    //fftw2 = fftw_plan_dft_r2c_2d(ny, nz, fftw_in, fftw_out_2, FFTW_MEASURE);
    fftw1 = fftw_plan_dft_1d(n, fftw_in, fftw_out, FFTW_FORWARD, FFTW_MEASURE);
    fftw2 = fftw_plan_dft_1d(n, fftw_in_2, fftw_out_2, FFTW_FORWARD, FFTW_MEASURE);

    for (int j = 0; j < n; j++) {
        double xrel = dx * j - 0.5 * l;
        fftw_in[j][0] = exp(-xrel * xrel / sigma / sigma);
        fftw_in[j][1] = 0;
    }
    fftw_execute(fftw1);

    fftw_in_2[0][0] = (fftw_in[1][0] - 2 * fftw_in[0][0] + fftw_in[n][0]) / dx / dx;
    fftw_in_2[0][1] = 0;
    for (int j = 1; j < n-1; j++) {
        fftw_in_2[j][0] = (fftw_in[j+1][0] - 2 * fftw_in[j][0] + fftw_in[j-1][0]) / dx / dx;
        fftw_in_2[j][1] = 0;
    }
    fftw_in_2[n-1][0] = (fftw_in[0][0] - 2 * fftw_in[n-1][0] + fftw_in[n-2][0]) / dx / dx;
    fftw_in_2[n-1][1] = 0;
    fftw_execute(fftw2);

    const int size = (n / 2) + 1;
    for (int j = 0; j < n; j++) {
        double multiplier = -(2 * cos(2 * PI * j / n) - 2) / dx / dx;

        double re1 = fftw_out[j][0];
        double im1 = fftw_out[j][1];

        double re2 = fftw_out_2[j][0];
        double im2 = fftw_out_2[j][1];

        double abs1 = sqrt(re1 * re1 + im1 * im1);
        double abs2 = sqrt(re2 * re2 + im2 * im2);

        std::cout << j << ":\t" << abs2 / abs1 << "\t expected: " << multiplier <<  "\n";
    }


    fftw_destroy_plan(fftw1);
    fftw_destroy_plan(fftw2);

    fftw_free(fftw_in);
    fftw_free(fftw_out);
    fftw_free(fftw_out_2);
}

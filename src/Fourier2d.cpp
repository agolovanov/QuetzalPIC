#include "Fourier2d.h"

Fourier2d::Fourier2d(int n1, int n2) : n1(n1), n2(n2) {
    in = fftw_alloc_real(n1 * n2);
    out = fftw_alloc_complex(n1 * (n2 / 2 + 1));
    update_plans();
}

Fourier2d::~Fourier2d() {
    fftw_destroy_plan(forward_plan);
    fftw_destroy_plan(backward_plan);

    if (in) {
        fftw_free(in);
    }
    if (out) {
        fftw_free(out);
    }
}

Fourier2d & Fourier2d::operator=(Fourier2d && other) {
    n1 = other.n1;
    n2 = other.n2;
    in = other.in;
    out = other.out;
    other.in = nullptr;
    other.out = nullptr;
    update_plans();
    return *this;
}

void Fourier2d::update_plans() {
    forward_plan = fftw_plan_dft_r2c_2d(n1, n2, in, out, FFTW_MEASURE);
    backward_plan = fftw_plan_dft_c2r_2d(n1, n2, out, in, FFTW_MEASURE);
}

void Fourier2d::forward_transform() {
    fftw_execute(forward_plan);
}

void Fourier2d::backward_transform() {
    fftw_execute(backward_plan);
}
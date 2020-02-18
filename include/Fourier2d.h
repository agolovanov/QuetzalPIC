#pragma once
#include <fftw3.h>

class Fourier2d {
public:
    double * in = nullptr;
    fftw_complex * out = nullptr;
    Fourier2d() = default;
    Fourier2d(int n1, int n2);
    ~Fourier2d();
    void forward_transform();
    void backward_transform();
    Fourier2d & operator=(Fourier2d && other);
private:
    int n1, n2;
    fftw_plan forward_plan;
    fftw_plan backward_plan;

    void update_plans();
};
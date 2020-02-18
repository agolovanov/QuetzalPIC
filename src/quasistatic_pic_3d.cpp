#include <chrono>
#include <fmt/format.h>
#include "containers_3d.h"
#include "System_parameters.h"
#include "System_3d.h"
#include <cmath>

int main() {
    auto t_begin = std::chrono::high_resolution_clock::now();

    System_parameters params;

    params.l = {20, 20, 20};
    params.d = {0.1, 0.1, 0.1};
    params.ppcy = 1;
    params.ppcz = 1;
    params.magnetic_field_iterations = 5;

    const double x0 = 4;
    const double y0 = params.l.y / 2;
    const double z0 = params.l.z / 2;
    const double xsigma = 2;
    const double ysigma = 2;
    const double zsigma = 2;
    const double a0 = sqrt(2.0);

    params.a_sqr = [x0, y0, z0, xsigma, ysigma, zsigma, a0] (double x, double y, double z) -> double {
        return 0.5 * a0 * a0 * exp(- (x-x0) * (x-x0) / xsigma / xsigma
                                   - (y-y0) * (y-y0) / ysigma / ysigma
                                   - (z-z0) * (z-z0) / zsigma / zsigma);
    };

    System_3d system{params};

    system.solve_wakefield();

    system.output();

    auto t_end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(t_end - t_begin).count();

    fmt::print("Finished in {} s\n", duration);
}
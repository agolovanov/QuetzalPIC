#include <chrono>
#include <fmt/format.h>
#include "containers.h"
#include "System_parameters.h"
#include "System_3d.h"
#include "Config_reader.h"
#include <cmath>
#include <omp.h>

double rhobunch(double xi, double y, double z) {
    const double x0 = 4;
    const double y0 = 10;
    const double z0 = 10;

    double xwidth = 2;
    double ywidth = 1.5;
    double zwidth = 1.5;

    double x_prof = (fabs(xi - x0) < xwidth) ? cos(0.5 * M_PI * (xi - x0) / xwidth) : 0.0;
    x_prof *= x_prof;
    double y_prof = (fabs(y - y0) < ywidth) ? cos(0.5 * M_PI * (y - y0) / ywidth) : 0.0;
    y_prof *= y_prof;
    double z_prof = (fabs(z - z0) < zwidth) ? cos(0.5 * M_PI * (z - z0) / zwidth) : 0.0;
    z_prof *= z_prof;

    return 0.0 * x_prof * y_prof * z_prof;
}

int main(int argc, char const *argv[]) {
    std::string config_filename;
    
    if (argc < 2) {
        std::cout << "Provide a path to the input .toml file as the first argument. Aborting." << std::endl;
        return -1;
    } else {
        if (argc > 2) {
            std::cout << "WARNING: more than one command-line argument are not supported. Additional arguments will be ignored" << std::endl;
        }
        config_filename = argv[1];
    }

    std::cout << "----------------------------------------" << std::endl;

    System_parameters params = Config_reader(config_filename, std::cout).get_parameters();

    std::cout << "----------------------------------------" << std::endl;

    std::cout << "OpenMP threads: " << omp_get_max_threads() << std::endl;

    fftw_init_threads();
    fftw_plan_with_nthreads(omp_get_max_threads());

    std::cout << "----------------------------------------" << std::endl;

    auto t_begin = std::chrono::high_resolution_clock::now();

    System_3d system{params, std::cout};

    system.solve_wakefield();

    auto t_end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(t_end - t_begin).count();

    std::cout << "----------------------------------------" << std::endl;

    fmt::print("Finished in {} s\n", duration);

    return 0;
}
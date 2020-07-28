#include <chrono>
#include <fmt/format.h>
#include "System_parameters.h"
#include "System_3d.h"
#include "Config_reader.h"
#include <omp.h>

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

    system.run();

    auto t_end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(t_end - t_begin).count();

    std::cout << "----------------------------------------" << std::endl;

    fmt::print("Finished in {} s\n", duration);

    return 0;
}
#include <chrono>
#include <fmt/format.h>
#include "containers_3d.h"
#include "System_parameters.h"
#include "System_3d.h"

int main() {
    auto t_begin = std::chrono::high_resolution_clock::now();

    System_parameters params;

    params.l = {20, 20, 20};
    params.d = {0.1, 0.1, 0.1};
    params.ppcy = 1;
    params.ppcz = 1;
    params.magnetic_field_iterations = 5;

    System_3d system{params};

    system.solve_wakefield();

    system.output();

    auto t_end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(t_end - t_begin).count();

    fmt::print("Finished in {} s\n", duration);
}
#include "containers_3d.h"
#include "System_parameters.h"
#include "System_3d.h"


int main() {
    System_parameters params;

    params.l = {20, 20, 20};
    params.d = {0.1, 0.1, 0.1};
    params.ppcy = 1;
    params.ppcz = 1;
    params.magnetic_field_iterations = 5;

    System_3d system{params};

    system.solve_wakefield();

    system.output();    
}
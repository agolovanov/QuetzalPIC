#include <iostream>

#include "containers.h"
#include <vector>

int main(int argc, char **argv) {

    int ny = 100;
    int nx = 100;

    double ly = 10;
    double lx = 10;

    std::vector<particle> particles(ny);

    array2d<double> psi(nx, ny);

    for (int i = 0; i < particles.size(); i++) {
        particles[i].y =
    }

}

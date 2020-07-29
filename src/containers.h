#pragma once

#include <vector>
#include <array>

template <class T>
struct vector3d_t {
    T x, y, z;
    vector3d_t<T> & operator/=(T value) {
        x /= value;
        y /= value;
        z /= value;
        return *this;
    }
};

using vector3d = vector3d_t<double>;
using ivector3d = vector3d_t<int>;

using vector2d = std::array<double, 2>;
using ivector2d = std::array<int, 2>;

enum class Plane {
    NONE=0, XY=1, XZ=2, YZ=3
};

enum class Axis {
    NONE=0, X=1, Y=1, Z=1
};

struct wake_particle_2d {
    double y = 0;
    double y_middle = 0;
    double z = 0;
    double z_middle = 0;
    double py = 0;
    double pz = 0;
    double py_next = 0;
    double pz_next = 0;
    double n = 0;
    double gamma = 1;
    double px = 0;
};

struct bunch_particle_3d {
    double x = 0;
    double y = 0;
    double z = 0;
    double px = 0;
    double py = 0;
    double pz = 0;
    double gamma = 1;
    double n = 0; 
};
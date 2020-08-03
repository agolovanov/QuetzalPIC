#include "array_utils.h"

#include <cmath>
#include <stdexcept>
#include <cassert>

template <class T>
array2d_t<T> calculate_xy_slice(const array3d_t<T> & array, double z0) {
    const auto n = array.get_dimensions();
    const auto origin = array.get_origin();
    const auto steps = array.get_steps();
    array2d_t<T> result({n.x, n.y}, vector2d{steps.x, steps.y}, vector2d{origin.x, origin.y}, Plane::XY, z0);

    int k0 = static_cast<int>(floor((z0 - origin.z) / steps.z));
    double z_frac = (z0 - origin.z) / steps.z - k0;
    #pragma omp parallel for
    for (int i = 0; i < n.x; i++) {
        for (int j = 0; j < n.y; j++) {
            result(i, j) = (1 - z_frac) * array(i, j, k0) + z_frac * array(i, j, k0+1);
        }
    }

    return result;
}

template <class T>
array1d_t<T> calculate_y_slice(const array2d_t<T> & array, double z0) {
    if (array.get_plane() != Plane::YZ) {
        throw std::invalid_argument("Plane should be YZ");
    }

    const auto n = array.get_dimensions();
    const auto origin = array.get_origin_3d();
    const auto steps = array.get_steps();
    array1d_t<T> result(n[0], steps[0], origin.y, Axis::Y, vector2d{origin.x, z0});

    int k0 = static_cast<int>(floor((z0 - origin.z) / steps[1]));
    double z_frac = (z0 - origin.z) / steps[1] - k0;
    #pragma omp parallel for
    for (int j = 0; j < n[0]; j++) {
        result(j) = (1 - z_frac) * array(j, k0) + z_frac * array(j, k0+1);
    }

    return result;
}

template array2d_t<double> calculate_xy_slice(const array3d_t<double> & array, double z0);
template array1d_t<double> calculate_y_slice(const array2d_t<double> & array, double z0);

void deposit(double y, double z, double value, array3d & array, int slice) {
    const auto d = array.get_steps();
    const auto n = array.get_dimensions();
    const auto origin = array.get_origin();

    int j1 = (int) floor((y - origin.y) / d.y);
    int j2 = (j1 + 1) % n.y;

    double y_frac = (y - origin.y) / d.y - j1;
    if (j1 < 0) {
        j1 += n.y;
    }

    int k1 = (int) floor((z - origin.z) / d.z);
    int k2 = (k1 + 1) % n.z;
    double z_frac = (z - origin.z) / d.z - k1;
    if (k1 < 0) {
        k1 += n.z;
    }

    assert((j1 >= 0) && (j1 < n.y));
    assert((j2 >= 0) && (j2 < n.y));
    assert((k1 >= 0) && (k1 < n.z));
    assert((k2 >= 0) && (k2 < n.z));

    #pragma omp atomic update
    array(slice, j1, k1) += value * (1 - y_frac) * (1 - z_frac);
    #pragma omp atomic update
    array(slice, j2, k1) += value * y_frac * (1 - z_frac);
    #pragma omp atomic update
    array(slice, j1, k2) += value * (1 - y_frac) * z_frac;
    #pragma omp atomic update
    array(slice, j2, k2) += value * y_frac * z_frac;
}

void deposit(double y, double z, double value, array2d & array) {
    const auto d = array.get_steps();
    const auto n = array.get_dimensions();
    const auto origin = array.get_origin_2d();
    
    int j1 = (int) floor((y - origin[0]) / d[0]);
    int j2 = (j1 + 1) % n[0];

    double y_frac = (y - origin[0]) / d[0] - j1;
    if (j1 < 0) {
        j1 += n[0];
    }

    int k1 = (int) floor((z - origin[1]) / d[1]);
    int k2 = (k1 + 1) % n[1];
    double z_frac = (z - origin[1]) / d[1] - k1;
    if (k1 < 0) {
        k1 += n[1];
    }

    assert((j1 >= 0) && (j1 < n[0]));
    assert((j2 >= 0) && (j2 < n[0]));
    assert((k1 >= 0) && (k1 < n[1]));
    assert((k2 >= 0) && (k2 < n[1]));

    #pragma omp atomic update
    array(j1, k1) += value * (1 - y_frac) * (1 - z_frac);
    #pragma omp atomic update
    array(j2, k1) += value * y_frac * (1 - z_frac);
    #pragma omp atomic update
    array(j1, k2) += value * (1 - y_frac) * z_frac;
    #pragma omp atomic update
    array(j2, k2) += value * y_frac * z_frac;
}

void deposit(double y, double z, double value, double * array, vector2d d, ivector2d n) {
    int j1 = (int) (y / d[0]);
    int j2 = (j1 + 1) % n[0];
    double y_frac = (y / d[0]) - j1;

    int k1 = (int) (z / d[1]);
    int k2 = (k1 + 1) % n[1];
    double z_frac = (z / d[1]) - k1;

    assert((j1 >= 0) && (j1 < n[0]));
    assert((j2 >= 0) && (j2 < n[0]));
    assert((k1 >= 0) && (k1 < n[1]));
    assert((k2 >= 0) && (k2 < n[1]));

    #pragma omp atomic update
    array[n[1] * j1 + k1] += value * (1 - y_frac) * (1 - z_frac);
    #pragma omp atomic update
    array[n[1] * j2 + k1] += value * y_frac * (1 - z_frac);
    #pragma omp atomic update
    array[n[1] * j1 + k2] += value * (1 - y_frac) * z_frac;
    #pragma omp atomic update
    array[n[1] * j2 + k2] += value * y_frac * z_frac;
}

double array_to_particle(double y, double z, const array3d & array, int slice) {
    const auto d = array.get_steps();
    const auto n = array.get_dimensions();
    const auto origin = array.get_origin();
    
    int j1 = (int) floor((y - origin.y) / d.y);
    int j2 = (j1 + 1) % n.y;
    double y_frac = (y - origin.y) / d.y - j1;
    if (j1 < 0) {
        j1 += n.y;
    }

    int k1 = (int) floor((z - origin.z) / d.z);
    int k2 = (k1 + 1) % n.z;
    double z_frac = (z - origin.z) / d.z - k1;
    if (k1 < 0) {
        k1 += n.z;
    }

    return array(slice, j1, k1) * (1 - y_frac) * (1 - z_frac) + array(slice, j2, k1) * y_frac * (1 - z_frac) +
           array(slice, j1, k2) * (1 - y_frac) * z_frac + array(slice, j2, k2) * y_frac * z_frac;
}

double array_to_particle(double y, double z, const array2d & array) {
    const auto d = array.get_steps();
    const auto n = array.get_dimensions();
    const auto origin = array.get_origin_2d();

    int j1 = (int) floor((y - origin[0]) / d[0]);
    int j2 = (j1 + 1) % n[0];
    double y_frac = (y - origin[0]) / d[0] - j1;
    if (j1 < 0) {
        j1 += n[0];
    }

    int k1 = (int) floor((z - origin[1]) / d[1]);
    int k2 = (k1 + 1) % n[1];
    double z_frac = (z - origin[1]) / d[1] - k1;
    if (k1 < 0) {
        k1 += n[1];
    }

    return array(j1, k1) * (1 - y_frac) * (1 - z_frac) + array(j2, k1) * y_frac * (1 - z_frac) +
           array(j1, k2) * (1 - y_frac) * z_frac + array(j2, k2) * y_frac * z_frac;
}

double array_yder_to_particle(double y, double z, const array3d & array, int slice) {
    const auto d = array.get_steps();
    const auto n = array.get_dimensions();

    int j1 = (int) floor(y / d.y - 0.5);
    double y_frac = (y / d.y - 0.5) - j1;
    if (j1 < 0) {
        j1 += n.y;
    }
    int j2 = (j1 + 1) % n.y;
    int j3 = (j1 + 2) % n.y;

    int k1 = (int) (z / d.z);
    int k2 = (k1 + 1) % n.z;
    double z_frac = (z / d.z) - k1;

    double der_y1_z1 = (array(slice, j2, k1) - array(slice, j1, k1)) / d.y;
    double der_y2_z1 = (array(slice, j3, k1) - array(slice, j2, k1)) / d.y;
    double der_y1_z2 = (array(slice, j2, k2) - array(slice, j1, k2)) / d.y;
    double der_y2_z2 = (array(slice, j3, k2) - array(slice, j2, k2)) / d.y;

    return der_y1_z1 * (1 - y_frac) * (1 - z_frac) + der_y2_z1 * y_frac * (1 - z_frac) +
           der_y1_z2 * (1 - y_frac) * z_frac + der_y2_z2 * y_frac * z_frac;
}

double array_zder_to_particle(double y, double z, const array3d & array, int slice) {
    const auto d = array.get_steps();
    const auto n = array.get_dimensions();

    int j1 = (int) (y / d.y);
    int j2 = (j1 + 1) % n.y;
    double y_frac = (y / d.y) - j1;

    int k1 = (int) floor(z / d.z - 0.5);
    double z_frac = (z / d.z - 0.5) - k1;
    if (k1 < 0) {
        k1 += n.z;
    }
    int k2 = (k1 + 1) % n.z;
    int k3 = (k1 + 2) % n.z;

    double der_y1_z1 = (array(slice, j1, k2) - array(slice, j1, k1)) / d.z;
    double der_y2_z1 = (array(slice, j2, k2) - array(slice, j2, k1)) / d.z;
    double der_y1_z2 = (array(slice, j1, k3) - array(slice, j1, k2)) / d.z;
    double der_y2_z2 = (array(slice, j2, k3) - array(slice, j2, k2)) / d.z;

    return der_y1_z1 * (1 - y_frac) * (1 - z_frac) + der_y2_z1 * y_frac * (1 - z_frac) +
           der_y1_z2 * (1 - y_frac) * z_frac + der_y2_z2 * y_frac * z_frac;
}

double array_yder_to_particle(double y, double z, const array2d & array) {
    const auto d = array.get_steps();
    const auto n = array.get_dimensions();

    int j1 = (int) floor(y / d[0] - 0.5);
    double y_frac = (y / d[0] - 0.5) - j1;
    if (j1 < 0) {
        j1 += n[0];
    }
    int j2 = (j1 + 1) % n[0];
    int j3 = (j1 + 2) % n[0];

    int k1 = (int) (z / d[1]);
    int k2 = (k1 + 1) % n[1];
    double z_frac = (z / d[1]) - k1;

    double der_y1_z1 = (array(j2, k1) - array(j1, k1)) / d[0];
    double der_y2_z1 = (array(j3, k1) - array(j2, k1)) / d[0];
    double der_y1_z2 = (array(j2, k2) - array(j1, k2)) / d[0];
    double der_y2_z2 = (array(j3, k2) - array(j2, k2)) / d[0];

    return der_y1_z1 * (1 - y_frac) * (1 - z_frac) + der_y2_z1 * y_frac * (1 - z_frac) +
           der_y1_z2 * (1 - y_frac) * z_frac + der_y2_z2 * y_frac * z_frac;
}

double array_zder_to_particle(double y, double z, const array2d & array) {
    const auto d = array.get_steps();
    const auto n = array.get_dimensions();

    int j1 = (int) (y / d[0]);
    int j2 = (j1 + 1) % n[0];
    double y_frac = (y / d[0]) - j1;

    int k1 = (int) floor(z / d[1] - 0.5);
    double z_frac = (z / d[1] - 0.5) - k1;
    if (k1 < 0) {
        k1 += n[1];
    }
    int k2 = (k1 + 1) % n[1];
    int k3 = (k1 + 2) % n[1];

    double der_y1_z1 = (array(j1, k2) - array(j1, k1)) / d[1];
    double der_y2_z1 = (array(j2, k2) - array(j2, k1)) / d[1];
    double der_y1_z2 = (array(j1, k3) - array(j1, k2)) / d[1];
    double der_y2_z2 = (array(j2, k3) - array(j2, k2)) / d[1];

    return der_y1_z1 * (1 - y_frac) * (1 - z_frac) + der_y2_z1 * y_frac * (1 - z_frac) +
           der_y1_z2 * (1 - y_frac) * z_frac + der_y2_z2 * y_frac * z_frac;
}

void deposit(double x, double y, double z, double value, array3d & array) {
    const auto d = array.get_steps();
    const auto n = array.get_dimensions();
    const auto origin = array.get_origin();

    int i1 = (int) floor((x - origin.x) / d.x);
    int i2 = (i1 + 1);
    double x_frac = (x - origin.x) / d.x - i1;

    if ((i1 < 0) || (i2 >= n.x)) {
        return;
    }

    int j1 = (int) floor((y - origin.y) / d.y);
    int j2 = (j1 + 1) % n.y;
    double y_frac = (y - origin.y) / d.y - j1;
    if (j1 < 0) {
        j1 += n.y;
    }

    int k1 = (int) floor((z - origin.z) / d.z);
    int k2 = (k1 + 1) % n.z;
    double z_frac = (z - origin.z) / d.z - k1;
    if (k1 < 0) {
        k1 += n.z;
    }

    assert((j1 >= 0) && (j1 < n.y));
    assert((j2 >= 0) && (j2 < n.y));
    assert((k1 >= 0) && (k1 < n.z));
    assert((k2 >= 0) && (k2 < n.z));

    #pragma omp atomic update
    array(i1, j1, k1) += value * (1 - x_frac) * (1 - y_frac) * (1 - z_frac);
    #pragma omp atomic update
    array(i1, j2, k1) += value * (1 - x_frac) * y_frac * (1 - z_frac);
    #pragma omp atomic update
    array(i1, j1, k2) += value * (1 - x_frac) * (1 - y_frac) * z_frac;
    #pragma omp atomic update
    array(i1, j2, k2) += value * (1 - x_frac) * y_frac * z_frac;
    #pragma omp atomic update
    array(i2, j1, k1) += value * x_frac * (1 - y_frac) * (1 - z_frac);
    #pragma omp atomic update
    array(i2, j2, k1) += value * x_frac * y_frac * (1 - z_frac);
    #pragma omp atomic update
    array(i2, j1, k2) += value * x_frac * (1 - y_frac) * z_frac;
    #pragma omp atomic update
    array(i2, j2, k2) += value * x_frac * y_frac * z_frac;
}

double array_to_particle(double x, double y, double z, const array3d & array) {
    const auto d = array.get_steps();
    const auto n = array.get_dimensions();
    const auto origin = array.get_origin();

    int i1 = (int) floor((x - origin.x) / d.x);
    int i2 = (i1 + 1);
    double x_frac = (x - origin.x) / d.x - i1;

    if ((i1 < 0) || (i2 >= n.x)) {
        return 0.0;
    }
    
    int j1 = (int) floor((y - origin.y) / d.y);
    int j2 = (j1 + 1) % n.y;
    double y_frac = (y - origin.y) / d.y - j1;
    if (j1 < 0) {
        j1 += n.y;
    }

    int k1 = (int) floor((z - origin.z) / d.z);
    int k2 = (k1 + 1) % n.z;
    double z_frac = (z - origin.z) / d.z - k1;
    if (k1 < 0) {
        k1 += n.z;
    }

    return array(i1, j1, k1) * (1 - x_frac) * (1 - y_frac) * (1 - z_frac) + 
           array(i1, j2, k1) * (1 - x_frac) * y_frac * (1 - z_frac) +
           array(i1, j1, k2) * (1 - x_frac) * (1 - y_frac) * z_frac + 
           array(i1, j2, k2) * (1 - x_frac) * y_frac * z_frac +
           array(i2, j1, k1) * x_frac * (1 - y_frac) * (1 - z_frac) + 
           array(i2, j2, k1) * x_frac * y_frac * (1 - z_frac) +
           array(i2, j1, k2) * x_frac * (1 - y_frac) * z_frac + 
           array(i2, j2, k2) * x_frac * y_frac * z_frac;
}
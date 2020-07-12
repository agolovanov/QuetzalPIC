#include "containers_3d.h"
#include <cmath>
#include <stdexcept>

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
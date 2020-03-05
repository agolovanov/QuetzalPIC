#include "containers_3d.h"
#include <cmath>

template <class T>
array2d_t<T> calculate_xy_slice(const array3d_t<T> & array, double z0) {
    const auto n = array.get_dimensions();
    const auto origin = array.get_origin();
    const auto steps = array.get_steps();
    array2d_t<T> result({n.x, n.y}, {steps.x, steps.y}, {origin.x, origin.y}, Plane::XY, z0);

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

template array2d_t<double> calculate_xy_slice(const array3d_t<double> & array, double z0);
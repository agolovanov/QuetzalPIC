#include "profiles.h"
#include <cmath>

std::function<double(double, double, double)> gaussian(double amplitude, dvector3d width, dvector3d r0) {
    auto func = [amplitude, width, r0] (double x, double y, double z) -> double {
        return amplitude * exp(- (x-r0.x) * (x-r0.x) / width.x / width.x
                               - (y-r0.y) * (y-r0.y) / width.y / width.y
                               - (z-r0.z) * (z-r0.z) / width.z / width.z);
    };
    return func;
}
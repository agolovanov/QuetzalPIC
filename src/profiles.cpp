#include "profiles.h"
#include <cmath>

std::function<double(double, double, double)> gaussian3d(double amplitude, vector3d width, vector3d r0) {
    auto func = [amplitude, width, r0] (double x, double y, double z) -> double {
        return amplitude * exp(- (x-r0.x) * (x-r0.x) / width.x / width.x
                               - (y-r0.y) * (y-r0.y) / width.y / width.y
                               - (z-r0.z) * (z-r0.z) / width.z / width.z);
    };
    return func;
}

std::function<double(double, double, double)> constant3d(double amplitude) {
    // Cast to void added to suppress warnings
    auto func = [amplitude] (double x, double y, double z) -> double { (void)x; (void)y; (void)z; return amplitude; };
    return func;
}

std::function<double(double, double, double)> parabolic3d(double amplitude, double xsize, double rsize, vector3d r0) {
    auto func = [amplitude, xsize, rsize, r0] (double x, double y, double z) -> double {
        const double r_squared = (y - r0.y) * (y - r0.y) + (z - r0.z) * (z - r0.z);
        const double x_rel_squared = (x - r0.x) * (x - r0.x);
        if ((r_squared <= rsize * rsize) && (x_rel_squared < xsize * xsize)) {
            return amplitude * (1 - x_rel_squared / xsize / xsize) * (1 - r_squared / rsize / rsize);
        } else {
            return 0.0;
        }
    };
    return func;
}

std::function<double(double, double)> constant2d(double amplitude) {
    // Cast to void added to suppress warnings
    auto func = [amplitude] (double x, double y) -> double { (void)x; (void)y; return amplitude; };
    return func;
}

std::function<double(double, double)> powerlaw2d(double power, double radius, double x0, double y0, double min_factor, double max_factor, double amplitude) {
    auto func = [amplitude, radius, power, x0, y0, min_factor, max_factor] (double x, double y) -> double { 
        double r = sqrt((x - x0) * (x - x0) + (y - y0) * (y - y0));
        double value = pow(r / radius, power);
        if (value > max_factor) {
            value = max_factor;
        }
        if (value < min_factor) {
            value = min_factor;
        }
        return amplitude * value;
    };
    return func;
}
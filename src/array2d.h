#pragma once

#include "containers.h"

template <class T>
class array2d_t {
public:
    array2d_t() = default;

    array2d_t(const ivector2d n, const vector2d steps = {0, 0}, const vector2d origin = {0, 0}, 
              Plane plane = Plane::NONE, double plane_coordinate = 0);

    array2d_t(const ivector2d n, const vector3d steps, const vector3d origin, Plane plane);

    inline T& operator()(const int i, const int j) {
        return data[n[1] * i + j];
    }

    inline const T& operator()(const int i, const int j) const {
        return data[n[1] * i + j];
    }

    inline auto get_n1() const {
        return n[0];
    }

    inline auto get_n2() const {
        return n[1];
    }

    inline auto get_dimensions() const {
        return n;
    }

    inline auto get_origin_3d() const {
        return origin;
    }

    inline vector2d get_steps() const {
        return steps;
    }

    inline vector2d get_origin_2d() const {
        return origin_2d;
    }

    inline auto get_plane() const { 
        return plane;
    }

    double get_plane_coordinate() const;

private:
    ivector2d n;
    std::vector<T> data;
    vector2d steps;
    vector3d origin;
    vector2d origin_2d;
    Plane plane;
};

using array2d = array2d_t<double>;

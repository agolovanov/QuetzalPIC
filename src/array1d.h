#pragma once

#include "containers.h"

template <class T>
class array1d_t {
public:
    array1d_t() = default;

    array1d_t(const int n, const double step = 0, const double origin = 0, 
              Axis axis = Axis::NONE, vector2d axis_coordinate = {0, 0});

    array1d_t(const int n, const vector3d steps, const vector3d origin, Axis axis);

    inline T& operator()(const int i) {
        return data[i];
    }

    inline const T& operator()(const int i) const {
        return data[i];
    }

    inline auto get_size() const {
        return n;
    }

    inline auto get_step() const {
        return step;
    }

    inline auto get_origin_3d() const {
        return origin;
    }

    double get_origin_1d() const;

    auto get_axis() const { 
        return axis;
    }

    vector2d get_axis_coordinate() const;

private:
    int n;
    std::vector<T> data;
    double step;
    vector3d origin;
    Axis axis;
};

using array1d = array1d_t<double>;
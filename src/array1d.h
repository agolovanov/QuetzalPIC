#pragma once

#include "containers.h"

template <class T>
class array1d_t {
public:
    array1d_t() = default;

    array1d_t(const int n, const double step = 0, const double origin = 0, 
              Axis axis = Axis::NONE, vector2d axis_coordinate = {0, 0}) : 
        n(n), data(n), step(step), axis(axis) { 
        if ((this->axis == Axis::NONE) || (this->axis == Axis::X)) {
            this->origin = {origin, axis_coordinate[0], axis_coordinate[1]};
        } else if (this->axis == Axis::Y) {
            this->origin = {axis_coordinate[0], origin, axis_coordinate[1]};
        } else {
            this->origin = {axis_coordinate[0], axis_coordinate[1], origin};
        }
    }

    array1d_t(const int n, const vector3d steps, const vector3d origin, Axis axis) : 
        n(n), data(n), origin(origin), axis(axis) {
        
        if ((this->axis == Axis::NONE) || (this->axis == Axis::X)) {
            this->step = steps.x;
        } else if (this->axis == Axis::Y) {
            this->step = steps.y;
        } else {
            this->step = steps.z;
        }
    }

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

    inline auto get_origin_1d() const {
        if ((axis == Axis::NONE) || (axis == Axis::X)) {
            return origin.x;
        } else if (axis == Axis::Y) {
            return origin.y;
        } else {
            return origin.z;
        }
    }

    inline auto get_axis() const { 
        return axis;
    }

    inline vector2d get_axis_coordinate() const {
        if ((axis == Axis::NONE) || (axis == Axis::X)) {
            return {origin.y, origin.z};
        } else if (axis == Axis::Y) {
            return {origin.x, origin.z};
        } else {
            return {origin.y, origin.z};
        }
    }

private:
    int n;
    std::vector<T> data;
    double step;
    vector3d origin;
    Axis axis;
};

using array1d = array1d_t<double>;
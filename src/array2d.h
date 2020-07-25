#pragma once

#include "containers.h"

template <class T>
class array2d_t {
public:
    array2d_t() = default;

    array2d_t(const ivector2d n, const vector2d steps = {0, 0}, const vector2d origin = {0, 0}, 
              Plane plane = Plane::NONE, double plane_coordinate = 0) : 
        n(n), data(n[0] * n[1]), steps(steps), origin_2d(origin), plane(plane) { 
        if ((this->plane == Plane::NONE) || (this->plane == Plane::XY)) {
            this->origin = {origin[0], origin[1], plane_coordinate};
        } else if (this->plane == Plane::XZ) {
            this->origin = {origin[0], plane_coordinate, origin[1]};
        } else {
            this->origin = {plane_coordinate, origin[0], origin[1]};
        }
    }

    array2d_t(const ivector2d n, const vector3d steps, const vector3d origin, Plane plane) : 
        n(n), data(n[0] * n[1]), origin(origin), plane(plane) { 

        if ((this->plane == Plane::NONE) || (this->plane == Plane::XY)) {
            this->steps = {steps.x, steps.y};
            this->origin_2d = {origin.x, origin.y};
        } else if (this->plane == Plane::XZ) {
            this->steps = {steps.x, steps.z};
            this->origin_2d = {origin.x, origin.z};
        } else {
            this->steps = {steps.y, steps.z};
            this->origin_2d = {origin.y, origin.z};
        }
    }

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

    inline auto get_plane_coordinate() const {
        if ((plane == Plane::NONE) || (plane == Plane::XY)) {
            return origin.z;
        } else if (plane == Plane::XZ) {
            return origin.y;
        } else {
            return origin.x;
        }
    }

private:
    ivector2d n;
    std::vector<T> data;
    vector2d steps;
    vector3d origin;
    vector2d origin_2d;
    Plane plane;
};

using array2d = array2d_t<double>;

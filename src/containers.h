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


template <class T>
class array3d_t {
public:
    array3d_t() = default;

    array3d_t(const ivector3d n, const vector3d d = {0, 0, 0}, const vector3d origin = {0, 0, 0}) :
        n(n), data(n.x * n.y * n.z), d(d), origin(origin) {}

    inline T& operator()(const size_t i, const size_t j, const size_t k) {
        return data[n.y * n.z * i + n.z * j + k];
    }

    inline const T& operator()(const size_t i, const size_t j, const size_t k) const {
        return data[n.y * n.z * i + n.z * j + k];
    }

    inline auto get_dimensions() const {
        return n;
    }

    inline auto get_n1() const {
        return n.x;
    }

    inline auto get_n2() const {
        return n.y;
    }

    inline auto get_n3() const {
        return n.z;
    }

    inline auto get_origin() const {
        return origin;
    }

    inline auto get_steps() const {
        return d;
    }

private:
    ivector3d n;
    std::vector<T> data;
    vector3d d;
    vector3d origin;
};

using array3d = array3d_t<double>;

using vector2d = std::array<double, 2>;
using ivector2d = std::array<int, 2>;

enum class Plane {
    NONE=0, XY=1, XZ=2, YZ=3
};

template <class T>
class array2d_t {
public:
    array2d_t() = default;

    array2d_t(const ivector2d n, const vector2d steps = {0, 0}, const vector2d origin = {0, 0}, 
              Plane plane = Plane::NONE, double plane_coordinate = 0) : 
        n(n), data(n[0] * n[1]), steps(steps), plane(plane) { 
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
        } else if (this->plane == Plane::XZ) {
            this->steps = {steps.x, steps.z};
        } else {
            this->steps = {steps.y, steps.z};
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
        if ((plane == Plane::NONE) || (plane == Plane::XY)) {
            return {origin.x, origin.y};
        } else if (plane == Plane::XZ) {
            return {origin.x, origin.z};
        } else {
            return {origin.y, origin.z};
        }
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
    Plane plane;
};

using array2d = array2d_t<double>;

enum class Axis {
    NONE=0, X=1, Y=1, Z=1
};

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


struct particle {
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

template <class T>
array2d_t<T> calculate_xy_slice(const array3d_t<T> & array, double z0);

template <class T>
array1d_t<T> calculate_y_slice(const array2d_t<T> & array, double z0);
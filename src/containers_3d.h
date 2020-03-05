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
        n(n), data(n[0] * n[1]), steps(steps), origin(origin), plane(plane), plane_coordinate(plane_coordinate) { }

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

    inline auto get_steps() const {
        return steps;
    }

    inline auto get_origin() const {
        return origin;
    }

    inline auto get_plane() const { 
        return plane;
    }

    inline auto get_plane_coordinate() const {
        return plane_coordinate;
    }

private:
    ivector2d n;
    std::vector<T> data;
    vector2d steps;
    vector2d origin;
    Plane plane;
    double plane_coordinate;
};

using array2d = array2d_t<double>;

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
#pragma once

#include "containers.h"

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
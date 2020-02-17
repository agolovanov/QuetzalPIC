#pragma once

#include <vector>

template <class T>
struct vector3d {
    T x, y, z;
};

using dvector3d = vector3d<double>;
using ivector3d = vector3d<int>;


template <class T>
class array3d_t {
public:
    array3d_t() = default;

    array3d_t(const size_t n1, const size_t n2, const size_t n3) :
            n1(n1), n2(n2), n3(n3), data(n1 * n2 * n3) {
    }

    array3d_t(const ivector3d n) : array3d_t(n.x, n.y, n.z) {}

    inline T& operator()(const size_t i, const size_t j, const size_t k) {
        return data[n2 * n3 * i + n3 * j + k];
    }

    inline const T& operator()(const size_t i, const size_t j, const size_t k) const {
        return data[n2 * n3 * i + n3 * j + k];
    }

    inline auto get_n1() const {
        return n1;
    }

    inline auto get_n2() const {
        return n2;
    }

    inline auto get_n3() const {
        return n3;
    }

private:
    size_t n1, n2, n3;
    std::vector<T> data;
};

using array3d = array3d_t<double>;

template <class T>
class array2d_t {
public:
    array2d_t() = default;

    array2d_t(const size_t n1, const size_t n2) : n1(n1), n2(n2), rows(n1), data(n1 * n2) {
        for (int i = 0; i < n1; i++) {
            rows[i] = &data[n2 * i];
        }
    }

    inline T& operator()(const size_t i, const size_t j) {
        return rows[i][j];
    }

    inline const T& operator()(const size_t i, const size_t j) const {
        return rows[i][j];
    }

    inline auto get_n1() const {
        return n1;
    }

    inline auto get_n2() const {
        return n2;
    }

private:
    size_t n1, n2;
    std::vector<T*> rows;
    std::vector<T> data;
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
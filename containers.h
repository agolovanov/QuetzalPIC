#pragma once

#include <vector>

template <class T>
class array2d {
public:
    array2d(const size_t n1, const size_t n2) : n1(n1), n2(n2), rows(n1), data(n1 * n2) {
        for (int i = 0; i < n1; i++) {
            rows[i] = &data[n2 * i];
        }
    }

    inline T& operator()(const size_t i, const size_t j) {
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


struct particle {
    double y = 0;
    double y_next = 0;
    double py = 0;
    double py_next = 0;
    double n = 0;
    double gamma = 1;
    double px = 0;
};

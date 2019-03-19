#pragma once

#include <vector>

template <class T>
class array3d {
public:
    array3d(const size_t n1, const size_t n2, const size_t n3) :
            n1(n1), n2(n2), n3(n3), ppointers(n1), pointers(n1 * n2), data(n1 * n2 * n3) {
        /*for (int i = 0; i < n1; i++) {
            ppointers[i] = &pointers[n2 * i];
            for (int j = 0; j < n2; j++) {
                int index = n2 * i + j;
                pointers[index] = &data[n3 * index];
            }
        }
        */
    }

    inline T& operator()(const size_t i, const size_t j, const size_t k) {
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
    std::vector<T**> ppointers;
    std::vector<T*> pointers;
    std::vector<T> data;
};

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

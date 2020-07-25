#include "array1d.h"

template <class T>
array1d_t<T>::array1d_t(const int n, const double step, const double origin, Axis axis, vector2d axis_coordinate) : 
        n(n), data(n), step(step), axis(axis) { 
    if ((this->axis == Axis::NONE) || (this->axis == Axis::X)) {
        this->origin = {origin, axis_coordinate[0], axis_coordinate[1]};
    } else if (this->axis == Axis::Y) {
        this->origin = {axis_coordinate[0], origin, axis_coordinate[1]};
    } else {
        this->origin = {axis_coordinate[0], axis_coordinate[1], origin};
    }
}

template <class T>
array1d_t<T>::array1d_t(const int n, const vector3d steps, const vector3d origin, Axis axis) : 
    n(n), data(n), origin(origin), axis(axis) {
    
    if ((this->axis == Axis::NONE) || (this->axis == Axis::X)) {
        this->step = steps.x;
    } else if (this->axis == Axis::Y) {
        this->step = steps.y;
    } else {
        this->step = steps.z;
    }
}

template <class T>
double array1d_t<T>::get_origin_1d() const {
    if ((axis == Axis::NONE) || (axis == Axis::X)) {
        return origin.x;
    } else if (axis == Axis::Y) {
        return origin.y;
    } else {
        return origin.z;
    }
}

template <class T>
vector2d array1d_t<T>::get_axis_coordinate() const {
    if ((axis == Axis::NONE) || (axis == Axis::X)) {
        return {origin.y, origin.z};
    } else if (axis == Axis::Y) {
        return {origin.x, origin.z};
    } else {
        return {origin.y, origin.z};
    }
}

template class array1d_t<double>;
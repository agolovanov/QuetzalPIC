#include "array2d.h"

template <class T>
array2d_t<T>::array2d_t(const ivector2d n, const vector2d steps, const vector2d origin, Plane plane, 
                        double plane_coordinate) : 
        n(n), data(n[0] * n[1]), steps(steps), origin_2d(origin), plane(plane) { 
    if ((this->plane == Plane::NONE) || (this->plane == Plane::XY)) {
        this->origin = {origin[0], origin[1], plane_coordinate};
    } else if (this->plane == Plane::XZ) {
        this->origin = {origin[0], plane_coordinate, origin[1]};
    } else {
        this->origin = {plane_coordinate, origin[0], origin[1]};
    }
}

template <class T>
array2d_t<T>::array2d_t(const ivector2d n, const vector3d steps, const vector3d origin, Plane plane) : 
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

template <class T>
double array2d_t<T>::get_plane_coordinate() const {
    if ((plane == Plane::NONE) || (plane == Plane::XY)) {
        return origin.z;
    } else if (plane == Plane::XZ) {
        return origin.y;
    } else {
        return origin.x;
    }
}

template class array2d_t<double>;
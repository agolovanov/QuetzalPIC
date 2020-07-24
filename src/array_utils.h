#pragma once

#include "array1d.h"
#include "array2d.h"
#include "array3d.h"

template <class T>
array1d_t<T> calculate_y_slice(const array2d_t<T> & array, double z0);

template <class T>
array2d_t<T> calculate_xy_slice(const array3d_t<T> & array, double z0);
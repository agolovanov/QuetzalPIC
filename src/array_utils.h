#pragma once

#include "array1d.h"
#include "array2d.h"
#include "array3d.h"

template <class T>
array1d_t<T> calculate_y_slice(const array2d_t<T> & array, double z0);

template <class T>
array2d_t<T> calculate_xy_slice(const array3d_t<T> & array, double z0);

void deposit(double y, double z, double value, array3d & array, int slice);
void deposit(double y, double z, double value, array2d & array);
void deposit(double y, double z, double value, double * array, vector2d d, ivector2d n);
double array_to_particle(double y, double z, const array3d & array, int slice);
double array_to_particle(double y, double z, const array2d & array);
double array_yder_to_particle(double y, double z, const array3d & array, int slice);
double array_zder_to_particle(double y, double z, const array3d & array, int slice);
double array_yder_to_particle(double y, double z, const array2d & array);
double array_zder_to_particle(double y, double z, const array2d & array);
void deposit(double x, double y, double z, double value, array3d & array);
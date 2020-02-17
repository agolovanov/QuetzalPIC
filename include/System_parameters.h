#pragma once

#include "containers_3d.h"

struct System_parameters {
    dvector3d l;
    dvector3d d;
    int ppcy = 1;
    int ppcz = 1;
    int magnetic_field_iterations = 1;
};
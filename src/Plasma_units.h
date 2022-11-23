#pragma once

#include "constants.h"
#include <cmath>
#include <iostream>

constexpr double frequency_to_density_CGS(const double frequency) {
    return ELECTRON_MASS_CGS * frequency * frequency / (4 * PI * ELEMENTARY_CHARGE_CGS * ELEMENTARY_CHARGE_CGS);
}

constexpr double density_to_frequency(const double density_CGS) {
    return sqrt(4 * PI * ELEMENTARY_CHARGE_CGS * ELEMENTARY_CHARGE_CGS * density_CGS / ELECTRON_MASS_CGS);
}

class Plasma_units {
public:
    const double frequency;
    const double plasma_wavenumber;
    const double wavelength;
    const double density;
    const double spatial_charge_norm;
    const double field_schwinger;
    Plasma_units(double base_frequency_SI) :
        frequency(base_frequency_SI),
        plasma_wavenumber(frequency / SPEED_OF_LIGHT_CGS),
        wavelength(2 * PI / plasma_wavenumber),
        density(frequency_to_density_CGS(frequency)),
        spatial_charge_norm(density / plasma_wavenumber / plasma_wavenumber / plasma_wavenumber),
        field_schwinger(ELECTRON_MASS_CGS * SPEED_OF_LIGHT_CGS * SPEED_OF_LIGHT_CGS / PLANCK_CONST_BAR_CGS / frequency) { }

    Plasma_units() : Plasma_units(1.0) {};
};
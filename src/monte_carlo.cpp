#include "monte_carlo.h"

#include <cstdlib>
#include <cmath>
#include "constants.h"

double rand_double() {
    return static_cast<double>(rand()) / RAND_MAX;
}

double W_new(double g, double r, double chi, double omega) {
    if ((1.0 - r) * g < 1.0 || chi == 0.0 || r == 0.0) {
        return 0.0;
    }
    else if (chi > 0.13333) { // 0.13(3) corresponds to a cut off at 5 omega_c
        //constexpr double coef = 9.1e-28 * 3e10 * 3e10 / (1.05e-27 * omega); // mc^2 / hbar omega
        const double coef = ELECTRON_MASS_CGS * SPEED_OF_LIGHT_CGS * SPEED_OF_LIGHT_CGS / (PLANCK_CONST_BAR_CGS * omega);
        const double y = r / ((1.0 - r) * chi); // hbar * omega / (chi * (gamma * mc^2 - hbar * omega))
        const double x = powf(y, 2.0 / 3.0);
        const double exp_v = exp(-2.0 * y / 3.0);
        const double iAiryapp = exp_v / (2.0 * sqrt(PI)) * powf(x + 0.80049, -0.75);
        const double dAiryapp = - exp_v / (2.0 * sqrt(PI)) * powf((x + 0.70861 * (1.0 - 0.65 * x / (1.0 + x * x))), 0.25);
        return - FINE_STRUCTURE_CONSTANT * coef / (g * g) * (iAiryapp + (2.0 / x + chi * r * sqrt(x)) * dAiryapp);
    }
    else {
        //let rm = 1.0 / (1.0 + 2.0 / (15.0 * chi));
        const double rm = chi / 0.13333;
        //let rm = 1.0;
        r = r * rm;
        const double coef = ELECTRON_MASS_CGS * SPEED_OF_LIGHT_CGS * SPEED_OF_LIGHT_CGS / (PLANCK_CONST_BAR_CGS * omega);
        const double y = r / ((1.0 - r) * chi); // hbar * omega / (chi * (gamma * mc^2 - hbar * omega))
        const double x = powf(y, 2.0 / 3.0);
        const double exp_v = exp(-2.0 * y / 3.0);
        const double iAiryapp = exp_v / (2.0 * sqrt(PI)) * powf(x + 0.80049, -0.75);
        const double dAiryapp = - exp_v / (2.0 * sqrt(PI)) * powf((x + 0.70861 * (1.0 - 0.65 * x / (1.0 + x * x))), 0.25);
        return - rm * FINE_STRUCTURE_CONSTANT * coef / (g * g) * (iAiryapp + (2.0 / x + chi * r * sqrt(x)) * dAiryapp);
    }
}
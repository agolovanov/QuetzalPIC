#include "Species_register.h"

#include <stdexcept>
#include <fmt/format.h>

int Species_register::register_species(std::string name, double charge, double mass) {
    if (is_registered(name)) {
        throw std::invalid_argument(fmt::format("Species [{}] already registered", name));
    }
    species.push_back(Species(name, charge, mass));
    return species.size() - 1;
}

bool Species_register::is_registered(std::string name) const {
    for (const auto & s : species) {
        if (s.name == name) {
            return true;
        }
    }
    return false;
}

int Species_register::find_species(std::string name) const {
    const auto size = species.size();
    for (size_t i = 0; i < size; i++) {
        if (species[i].name == name) {
            return i;
        }
    }
    return -1; // TODO That's C style. Return optional instead? 
}
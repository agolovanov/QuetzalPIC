#pragma once

#include <string>
#include <vector>

struct Species {
    double mass;
    double charge;
    double charge_to_mass_ratio;
    std::string name;
    Species(std::string name, double charge, double mass) : mass(mass), charge(charge), charge_to_mass_ratio(charge / mass), name(name) { };
};

class Species_register {
public:
    int register_species(std::string name, double charge, double mass);
    bool is_registered(std::string name) const;
    int find_species(std::string name) const;
    inline const Species & operator[](int i) const { return species[i]; }
private:
    std::vector<Species> species;  
};
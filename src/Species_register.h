#pragma once

#include <string>
#include <vector>

struct Species {
    double mass;
    double charge;
    double charge_to_mass_ratio;
    bool is_photon;
    std::string name;
    Species(std::string name, double charge, double mass, bool is_photon=false) : mass(mass), charge(charge), 
            charge_to_mass_ratio((charge / mass) ? mass != 0.0 : 0.0), is_photon(is_photon), name(name) { };
};

class Species_register {
public:
    int register_species(std::string name, double charge, double mass, bool is_photon = false);
    bool is_registered(std::string name) const;
    int find_species(std::string name) const;
    inline const Species & operator[](int i) const { return species[i]; }
    inline const std::vector<Species> & get_species() const { return species; }
private:
    std::vector<Species> species;
};
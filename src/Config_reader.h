#pragma once

#include <string>
#include <iostream>
#include <memory>
#include <cpptoml.h>
#include "System_parameters.h"

class Config_reader {
public:
    Config_reader(const std::string & filename, std::ostream & out);
    System_parameters get_parameters() {
        return params;
    }
private:
    System_parameters params;
    std::shared_ptr<cpptoml::table> config;
    std::ostream & out;
    
    template <class T>
    T read_value(const std::string & name);
    
    template <class T>
    T read_value(const std::string & name, std::shared_ptr<cpptoml::table> table);
    
    template <class T>
    T read_value(const std::string & name, T default_value);
    
    template <class T>
    T read_value(const std::string & name, T default_value, std::shared_ptr<cpptoml::table> table);
    void init_laser();
    Bunch_parameters parse_bunch_table(std::shared_ptr<cpptoml::table> bunch_table);
    void init_bunch();
    void init_output_parameters();
    void init_plasma_profile();
};
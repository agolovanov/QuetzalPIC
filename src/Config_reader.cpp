#include <cpptoml.h>

#include <exception>
#include <fmt/format.h>
#include "Config_reader.h"
#include "System_parameters.h"
#include "profiles.h"

class Config_exception : public std::exception {
public:
    Config_exception(std::string message) : message(message) {}
    const char* what() const noexcept { return message.c_str(); }
private:
    std::string message;
};

Config_reader::Config_reader(const std::string & filename, std::ostream & out) : out(out) {
    out << "Reading config " << filename << "\n" << std::endl;
    
    config = cpptoml::parse_file(filename);

    params.l.x = read_value<double>("lx");
    params.l.y = read_value<double>("ly");
    params.l.z = read_value<double>("lz");

    params.d.x = read_value<double>("dx");
    params.d.y = read_value<double>("dy");
    params.d.z = read_value<double>("dz");

    params.ppcy = read_value<int>("ppcy", 1);
    params.ppcz = read_value<int>("ppcz", 1);

    params.magnetic_field_iterations = read_value<int>("magnetic_field_iterations", 0);

    out << std::endl;

    init_laser();

    out << std::endl;

    init_output_parameters();
}

template <class T>
T Config_reader::read_value(const std::string & name) {
    return read_value<T>(name, config);
}

template <class T>
T Config_reader::read_value(const std::string & name, std::shared_ptr<cpptoml::table> table) {
    auto value = table->get_as<T>(name);
    if (value) { 
        out << fmt::format("{} = {}", name, *value) << std::endl;
        return *value;
    } else {
        throw Config_exception(fmt::format("Value [{}] doesn't exist or is of wrong type", name));
    }
}

template <class T>
T Config_reader::read_value(const std::string & name, T default_value) {
    return read_value<T>(name, default_value, config);
}

template <class T>
T Config_reader::read_value(const std::string & name, T default_value, std::shared_ptr<cpptoml::table> table) {
    auto value = table->get_as<T>(name);
    if (value) { 
        out << fmt::format("{} = {}", name, *value) << std::endl;
        return *value;
    } else {
        out << fmt::format("{} = {} (default)", name, default_value) << std::endl;
        return default_value;
    }
}

void Config_reader::init_laser() {
    if (config->contains("laser")) {
        out << "Parsing [laser] ..." << std::endl;
        auto laser_table = config->get_table("laser");

        const auto shape = read_value<std::string>("shape", laser_table);
        if (shape == "gaussian") {
            const double a0 = read_value<double>("a0", laser_table);
            const double xsigma = read_value<double>("xsigma", laser_table);
            const double ysigma = read_value<double>("ysigma", laser_table);
            const double zsigma = read_value<double>("zsigma", laser_table);
            const double x0 = read_value<double>("x0", xsigma, laser_table);
            const double y0 = read_value<double>("y0", 0.5 * params.l.y, laser_table);
            const double z0 = read_value<double>("z0", 0.5 * params.l.z, laser_table);
            
            vector3d width{xsigma, ysigma, zsigma};
            // conversion from full width at 1/e^2 to Gauss paramters
            width /= sqrt(8.0);
            
            params.a_sqr = gaussian(0.5 * a0 * a0, width, {x0, y0, z0});
        } else {
            throw Config_exception(fmt::format("Laser shape \"{}\" is not supported, use \"gaussian\".", shape));
        }

    } else {
        out << "No [laser] block in config" << std::endl;
    }
}

void Config_reader::init_output_parameters() {
    Output_parameters output_parameters;

    if (config->contains("output")) {
        out << "Parsing [output] ..." << std::endl;
        auto output_table = config->get_table("output");

        output_parameters.output3d = read_value<bool>("fields_3d", false, output_table);
        output_parameters.output_xy = read_value<bool>("fields_xy", false, output_table);
    } else {
        out << "No [output] block in config; no output files will be created" << std::endl;
    }

    params.output_parameters = output_parameters;
}
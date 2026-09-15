#pragma once
#include <array>
#include <string>
#include <toml.hpp>
namespace eks {
struct Config {
    std::string mode = "frozen";
    std::string boundary = "analytic";
    double mu = 0.0, fa = 1.0, amplitude = 1.0e-3, width = 2.0;
    double ko_sigma = 0.0;
    std::array<double,3> center{0.0,0.0,0.0};
    std::array<double,3> wavevector{0.5,0.25,0.0};
    bool allow_constraint_violating_id = false;
};
extern Config config;
void read_config(const toml::value& root);
toml::value config_as_toml();
void validate_runtime(unsigned int ts_mode);
void initial_data_physical(double x, double y, double z, double* u);
}
